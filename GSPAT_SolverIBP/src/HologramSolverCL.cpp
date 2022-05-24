#include <GSPAT_SolverIBP.h>
#include <src/HologramSolverCL.h>
#include <Helper/HelperMethods.h>
#include <Helper/OpenCLUtilityFunctions.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <clBLAS.h>
#include <glm/glm.hpp>						//GLM (Maths)
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <OpenCLSolverImpl_Interoperability.h>


char consoleLine[512];
static void printConsole(const char* str) {
	printf("%s", str);

}

static void printError(const char* str) {
	perror(str);
}

HologramSolverCL::HologramSolverCL(int numTransducers)
	: numTransducers(numTransducers), configured(false)
{ 
	//Registers default functions (in case the client did not register anything. It will be ignored it functions already set).
	GSPAT_IBP::RegisterPrintFuncs(printConsole, printConsole, printError);
		
	// Initializes the OpenCL platform
	initializeOpenCL();
	//3. Fill in initial guess for phase/amplitude...(create a method/class specifically for this?)
	for (int i = 0; i < HologramSolution::MAX_POINTS; i++) {
		initialGuess[2 * i] = 1.0f*i / HologramSolution::MAX_POINTS;
		initialGuess[2 * i + 1] = sqrtf(1 - initialGuess[2 * i] * initialGuess[2 * i]);
	}
	//Initialize global resources (shared by all solutions)
	createDirectivityTexture();
	createSolutionPool();
}

void* HologramSolverCL::getSolverContext() {
	OpenCL_ExecutionContext* context = new  OpenCL_ExecutionContext;
	context->context = this->context;
	context->device = this->deviceUsed;
	context->queue = queue;
	return (void*)context;
}

HologramSolverCL::~HologramSolverCL() {
	clFinish(queue);
	////1. Delete Solution pool
	delete solutionPool;
	GSPAT_IBP::printMessage_GSPAT("GSPAT: Thread pool succesfully destroyed");
	//2. Delete buffers
	cl_uint count;
	if (configured) {
		clGetMemObjectInfo(transducerMappings, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
		clReleaseMemObject(transducerMappings);
		clGetMemObjectInfo(phaseCorrections, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
		clReleaseMemObject(phaseCorrections);
		GSPAT_IBP::printMessage_GSPAT("GSPAT: Board config (mappings and phase delays) succesfully destroyed");
		configured = false;
	}
	clGetMemObjectInfo(directivityTexture, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(directivityTexture);
	clGetMemObjectInfo(transducerPositions, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(transducerPositions);
	//3. Delete kernels
	clReleaseKernel(fillAllPointHologramsKernel);
	clReleaseKernel(solveIBPKernel);
	clReleaseKernel(discretiseKernel);
	GSPAT_IBP::printMessage_GSPAT("GSPAT: Kernels succesfully destroyed");
	//4. Delete programs
	clReleaseProgram(program);
	GSPAT_IBP::printMessage_GSPAT("GSPAT: Program succesfully destroyed");
	//5. Delete  queue and context
	clReleaseCommandQueue(queue);
	clReleaseContext(context);
	GSPAT_IBP::printMessage_GSPAT("GSPAT: Context destroyed");
	
}

void HologramSolverCL::setBoardConfig(float* transducerPositions, int* transducerToPINMap, int* phaseAdjust, float* amplitudeAdjust, int numDiscreteLevels ) {
	//0. Unload previous configuraiton
	if (configured) {
		clReleaseMemObject(this->transducerPositions);//cl_mem
		clReleaseMemObject(transducerMappings);
		clReleaseMemObject(phaseCorrections);
		configured = false;
		GSPAT_IBP::printMessage_GSPAT("GSPAT: Discarding previous configuration (transducer IDs, phase adjustments)\n");

	}
	this->numDiscreteLevels = numDiscreteLevels;
	//1. create buffer with transducer positions (bit complex-> encapsulated as function):
	createTransducerPositionBuffer(transducerPositions);	
	//2. Configure PIN mappings (t), phase delays (p),...
	unsigned char* t=new unsigned char[numTransducers];//The mapping should be stored in bytes (not integers)
	float* p=new float[numTransducers];//We store phase adjusts in radians (not degrees)
	for (int i = 0; i < numTransducers; i++) {
		t[i] = (unsigned char)(transducerToPINMap[i]);
		p[i] = (float)( phaseAdjust[i]* M_PI / 180.0f);
	}
	//Now in correct format... let's copy them
	transducerMappings= clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, numTransducers* sizeof(unsigned char), t, NULL);	
	phaseCorrections= clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, numTransducers* sizeof(float), p, NULL);	
	amplitudeCorrections= clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, numTransducers* sizeof(float), amplitudeAdjust, NULL);	
	this->numDiscreteLevels = numDiscreteLevels;
	GSPAT_IBP::printMessage_GSPAT("GSPAT: Board configured (transducer positions, transducerIDs, phase Adjustments)\n");
	//Delete auxiliary memory
	delete t; delete p;	
	configured = true;

	//DEBUG: Check parameters uploaded (Transducer mappings, phase corrections):
	/*cl_int err;
	unsigned char PIN_mapping[512];
	float phaseAdjustGPU[512];
	err = clEnqueueReadBuffer(queue, transducerMappings, CL_TRUE, 0, 512*sizeof(unsigned char), &(PIN_mapping[0]),0, NULL, NULL);
	err = clEnqueueReadBuffer(queue, phaseCorrections, CL_TRUE, 0, 512*sizeof(float), &(phaseAdjustGPU[0]),0, NULL, NULL);
	for (int i = 0; i < 512; i++)phaseAdjustGPU[i] *= (180/M_PI);
	*///END DEBUG
}

void HologramSolverCL::initializeOpenCL(){
	cl_uint num_platforms;
	clGetPlatformIDs(0, NULL, &num_platforms);
	this->platforms = (cl_platform_id*)malloc(num_platforms*sizeof(cl_platform_id));
	clGetPlatformIDs(num_platforms, platforms, NULL);
	//this->platformUsed = platforms[num_platforms-1];
	this->platformUsed = platforms[0];

	// Initializes the OpenCL device
	cl_uint num_devices;
	clGetDeviceIDs( this->platformUsed, CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices);
	this->devices = (cl_device_id*)malloc(num_devices*sizeof(cl_device_id));
	clGetDeviceIDs( this->platformUsed, CL_DEVICE_TYPE_GPU, num_devices, this->devices, NULL);
	this->deviceUsed = this->devices[num_devices-1];

	// Creates the OpenCL context, queue, event and program.
	this->context = clCreateContext(NULL, 1, &(this->deviceUsed), NULL, NULL, NULL);
	this->queue = clCreateCommandQueue(this->context,  this->deviceUsed, 0, NULL);
	//Create kernels:
	char compilerOptions[256];
	sprintf(compilerOptions,"-cl-mad-enable -cl-finite-math-only -cl-fast-relaxed-math -cl-std=CL2.0 -D NUM_TRANSDUCERS=%d", this->numTransducers);
	//sprintf(compilerOptions,"-D NUM_TRANSDUCERS=%d", this->numTransducers);
	OpenCLUtilityFunctions::createProgramFromFile("hologramSolver_IBP.cl", context, deviceUsed, &program, compilerOptions);
	cl_int err;
	fillAllPointHologramsKernel = clCreateKernel(program, "computeFandB", &err);
    if(err < 0) {
		GSPAT_IBP::printError_GSPAT("GSPAT: Couldn't create kernel\n");
    }
	solveIBPKernel = clCreateKernel(program, "solvePhases_IBP", &err);
	if(err < 0) {
		sprintf(consoleLine,"GSPAT: Couldn't create a kernel: %d", err);
		GSPAT_IBP::printError_GSPAT(consoleLine);
	}
	discretiseKernel=clCreateKernel(program, "discretise", &err);
	if (err < 0) {
		sprintf(consoleLine,"GSPAT: Couldn't create a kernel: %d", err);
		GSPAT_IBP::printError_GSPAT(consoleLine);
	};
}

void HologramSolverCL::createDirectivityTexture()
{
	cl_int error;
	//Create 1D texture for transducer's directivity: 
	{
		size_t imageWidth = 500, imageHeight=1;
		float directivityBuffer[500], P_ref=8.02f;
		for (size_t p = 0; p < imageWidth-1; p++)
			directivityBuffer[p] = computeDirectivityForCosAlpha((1.0f*p) / (imageWidth-1));
		directivityBuffer[imageWidth - 1] = P_ref;
		cl_image_format textureFormat;
		textureFormat.image_channel_data_type = CL_FLOAT;
		textureFormat.image_channel_order = CL_LUMINANCE;
		directivityTexture= clCreateImage2D(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, &textureFormat, imageWidth, imageHeight, imageWidth * sizeof(float),  directivityBuffer, &error);			
	}
	if (error == NULL ) {
		GSPAT_IBP::printMessage_GSPAT("GSPAT: Directivity bufffers created successfully\n");
	}
	else 
		GSPAT_IBP::printError_GSPAT("GSPAT: Could not create directivity buffers\n");
}

void HologramSolverCL:: createTransducerPositionBuffer(float* transducerPositions) {
		size_t numElementsInBuffer = this->numTransducers* 4;
		float* positions = new float[numElementsInBuffer];
		//1. Transform positions to homogenenous coordinates.
		for (int t = 0; t < numTransducers;t++) {
				//fill in (x,y,z)
				memcpy(&(positions[4 * t]), &(transducerPositions[3 * t]), 3 * sizeof(float));
				//fill in homogenenous coordinate;
				positions[4 * t + 3] = 1;
			}
		//2. Create an OpenCL buffer and store positions:
		cl_int err;
		this->transducerPositions = clCreateBuffer(context, CL_MEM_READ_WRITE, numElementsInBuffer * sizeof(float), NULL, &err);
		if (err < 0) { GSPAT_IBP::printError_GSPAT("GSPAT: Couldn't create 'transducerPositions' buffer"); return ; }
		err = clEnqueueWriteBuffer(queue, this->transducerPositions, CL_TRUE, 0, numElementsInBuffer* sizeof(float), &(positions[0]), 0, NULL, NULL);
		if (err < 0) { GSPAT_IBP::printError_GSPAT("GSPAT: Couldn't write 'transducerPositions' buffer"); return ; }
		delete positions;
		//BEGIN DEBUG
		//static float result [32 * 16 * 4 ];//Check it was written...
		// err = clEnqueueReadBuffer(queue, transducerPositions, CL_TRUE, 0, 32*16*4*sizeof(float), &(result[0]),0, NULL, NULL);
		//if (err < 0) { perror("Couldn't read hologram"); return NULL; }
		//END DEBUG
	}

void HologramSolverCL::createSolutionPool() {
	solutionPool = new SolutionPool(this, numTransducers, context, queue);
}

Solution* HologramSolverCL::createSolution(int numPoints, int numGeometries, bool phaseOnly, float* positions, float* amplitudes, float* matStarts, float* matEnds, GSPAT::MatrixAlignment a) {
	HologramSolution* nextSolution; 
	//1. Access next available solution from the pool (and remove it from the list of "availables")
	nextSolution = solutionPool->getNextSolution();	
	
	//2. Adjust input data (alignment layout).
	float mStart[16 * HologramSolution::MAX_POINTS], mEnd[16 * HologramSolution::MAX_POINTS];	
	if (a == GSPAT::ColumnMajorAlignment) {
		//glm matrix are column major order... we need to transpose them
		for (int i = 0; i < numPoints; i++) {
			unsigned int offset = 16 * i;
			mStart[offset + 0] = matStarts[offset + 4 * 0 + 0];	mStart[offset + 1] = matStarts[offset + 4 * 1 + 0]; mStart[offset + 2] = matStarts[offset + 4 * 2 + 0]; mStart[offset + 3] = matStarts[offset + 4 * 3 + 0];
			mStart[offset + 4] = matStarts[offset + 4 * 0 + 1]; mStart[offset + 5] = matStarts[offset + 4 * 1 + 1]; mStart[offset + 6] = matStarts[offset + 4 * 2 + 1]; mStart[offset + 7] = matStarts[offset + 4 * 3 + 1];
			mStart[offset + 8] = matStarts[offset + 4 * 0 + 2]; mStart[offset + 9] = matStarts[offset + 4 * 1 + 2]; mStart[offset + 10]= matStarts[offset + 4 * 2 + 2]; mStart[offset + 11]= matStarts[offset + 4 * 3 + 2];
			mStart[offset +12] = matStarts[offset + 4 * 0 + 3]; mStart[offset + 13]= matStarts[offset + 4 * 1 + 3]; mStart[offset + 14]= matStarts[offset + 4 * 2 + 3]; mStart[offset + 15]= matStarts[offset + 4 * 3 + 3];

			mEnd[offset + 0] = matEnds[offset + 4 * 0 + 0];	mEnd[offset + 1] = matEnds[offset + 4 * 1 + 0]; mEnd[offset + 2] = matEnds[offset + 4 * 2 + 0]; mEnd[offset + 3] = matEnds[offset + 4 * 3 + 0];
			mEnd[offset + 4] = matEnds[offset + 4 * 0 + 1]; mEnd[offset + 5] = matEnds[offset + 4 * 1 + 1]; mEnd[offset + 6] = matEnds[offset + 4 * 2 + 1]; mEnd[offset + 7] = matEnds[offset + 4 * 3 + 1];
			mEnd[offset + 8] = matEnds[offset + 4 * 0 + 2]; mEnd[offset + 9] = matEnds[offset + 4 * 1 + 2]; mEnd[offset + 10]= matEnds[offset + 4 * 2 + 2]; mEnd[offset + 11]= matEnds[offset + 4 * 3 + 2];
			mEnd[offset + 12]= matEnds[offset + 4 * 0 + 3]; mEnd[offset + 13]= matEnds[offset + 4 * 1 + 3]; mEnd[offset + 14]= matEnds[offset + 4 * 2 + 3]; mEnd[offset + 15]= matEnds[offset + 4 * 3 + 3];
	
		}
	}
	else {
		memcpy(mStart, matStarts, 16 * numPoints * sizeof(float));
		memcpy(mEnd, matEnds, 16 * numPoints * sizeof(float));
	}
	//3. Update input data to use and return it.
	nextSolution->configureSolution(numPoints, numGeometries, phaseOnly, positions, amplitudes,  mStart, mEnd);
	return nextSolution ;	
}

Solution* HologramSolverCL::createSolutionExternalData(int numPoints, int numGeometries, bool phaseOnly, GSPAT::MatrixAlignment a) {
	HologramSolution* nextSolution; 
	//1. Access next available solution from the pool (and remove it from the list of "availables")
	nextSolution = solutionPool->getNextSolution();	
	nextSolution->configureSolutionExternalData(numPoints, numGeometries, phaseOnly);
	return nextSolution ;	
}

Solution* HologramSolverCL::createTransducersOffSolution() {
	HologramSolution* nextSolution; 
	//1. Access next available solution from the pool (and remove it from the list of "availables")
	nextSolution = solutionPool->getNextSolution();	
	nextSolution->configureSolutionTransducersOff(this->numTransducers);
	return nextSolution ;	
}

Solution* HologramSolverCL::createNewDividerSolution(unsigned char FPSDivider) {
	if (FPSDivider == 0)
		GSPAT_IBP::printWarning_GSPAT("GS_PAT: Incorrect divider value (zero). Command ignored (previous dividir still in use).");
	HologramSolution* nextSolution;
	//1. Access next available solution from the pool (and remove it from the list of "availables")
	nextSolution = solutionPool->getNextSolution();
	nextSolution->configureSolutionNewDivider(this->numTransducers, FPSDivider);
	return nextSolution;
}

void HologramSolverCL::releaseSolution(Solution* solution) {
	((HologramSolution*)solution)->releaseEvents();
	solutionPool->returnSolution((HologramSolution*)solution);
}


void HologramSolverCL::compute(Solution* solution) {
	HologramSolution* hs = ((HologramSolution*)solution);
	updateCLBuffers(hs);
	computeFandB(hs);
	solvePhases_IBP(hs);
	discretise(hs);
	clFlush(queue);
}

void HologramSolverCL::updateCLBuffers(HologramSolution * hs)
{
	cl_int err;
	hs->lockEvents();	
	if (hs->manualData) {
		err = clEnqueueWriteBuffer(queue, hs->positions_CLBuffer, CL_FALSE, NULL, 4 * hs->numPoints *hs->numGeometries * sizeof(float), hs->positions, 0, NULL, &(hs->events[GSPAT_Event::POSITIONS_UPLOADED]));
		err = clEnqueueWriteBuffer(queue, hs->targetAmplitudes_CLBuffer, CL_FALSE, NULL, hs->numPoints* hs->numGeometries * sizeof(float), hs->amplitudes, 0, NULL, &(hs->events[GSPAT_Event::AMPLITUDES_UPLOADED]));
		err = clEnqueueWriteBuffer(queue, hs->pointsReIm_CLBuffer, CL_FALSE, 0, 2 * hs->points() * sizeof(float), initialGuess, 0, NULL, &(hs->events[GSPAT_Event::INITIAL_GUESS_UPLOADED]));
		err = clEnqueueWriteBuffer(queue, hs->matrixStart_CLBuffer, CL_FALSE, NULL, 16 * (hs->numPoints) * sizeof(float), hs->matrixStart, 0, NULL, &(hs->events[GSPAT_Event::MATRIX_0_UPLOADED]));
		err = clEnqueueWriteBuffer(queue, hs->matrixEnd_CLBuffer, CL_FALSE, NULL, 16 * (hs->numPoints) * sizeof(float), hs->matrixEnd, 0, NULL, &(hs->events[GSPAT_Event::MATRIX_G_UPLOADED]));
		if (err < 0) { sprintf(consoleLine, "GSPAT::updateCLBuffers::Couldn't copy input parameters to the Solution"); GSPAT_IBP::printWarning_GSPAT(consoleLine); return; }
	}
}

void HologramSolverCL::computeFandB(HologramSolution* solution) {
	cl_int err;
	//B. Run the kernel
	size_t global_size[3] = {(size_t) numTransducers, 1,(size_t) solution->numPoints*solution->numGeometries};
	size_t local_size[3] =  {(size_t) numTransducers, 1, 1};
	//0. Setup inputs for the kernell (do once)
	err = clSetKernelArg(fillAllPointHologramsKernel, 0, sizeof(cl_mem), &(transducerPositions));
	if (err < 0) { GSPAT_IBP::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 0"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 1, sizeof(cl_mem), &(solution->positions_CLBuffer));
	if (err < 0) { GSPAT_IBP::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 1"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 2, sizeof(cl_mem), &(solution->matrixStart_CLBuffer));
	if (err < 0) { GSPAT_IBP::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 2"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 3, sizeof(cl_mem), &(solution->matrixEnd_CLBuffer));
	if (err < 0) { GSPAT_IBP::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 3"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 4, sizeof(int), &(solution->numPoints));
	if (err < 0) { GSPAT_IBP::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 4"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 5, sizeof(int), &(solution->numGeometries));
	if (err < 0) { GSPAT_IBP::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 5"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 6, sizeof(cl_mem), &(this->directivityTexture));
	if (err < 0) { GSPAT_IBP::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 6"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 7, sizeof(cl_mem), &(solution->singlePointField_CLBuffer));
	if (err < 0) { GSPAT_IBP::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 7"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 8, sizeof(cl_mem), &(solution->singlePointFieldNormalised_CLBuffer));
	if (err < 0) { GSPAT_IBP::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 8"); return ; }
		
	//1.3. Trigger the kernell
	cl_event dataUploaded[] = {solution->events[GSPAT_Event::POSITIONS_UPLOADED], solution->events[GSPAT_Event::MATRIX_0_UPLOADED], solution->events[GSPAT_Event::MATRIX_G_UPLOADED]};
	cl_uint numEvents = solution->manualData ? 3 : 1;
	err |= clEnqueueNDRangeKernel(queue, fillAllPointHologramsKernel, 3, NULL, global_size, local_size, numEvents, dataUploaded, &(solution->events[GSPAT_Event::F_AND_B_READY]));
	if (err < 0) { 
		GSPAT_IBP::printError_GSPAT("GSPAT: computeFandB:: Couldn't enqueue the kernel"); return ; 
	}	
	//{//DEBUG: Print input data: 
	//	float *initialGuess = new float[solution->numPoints * 2];
	//	float* positions = new float[solution->numPoints*solution->numGeometries * 4];
	//	float* amplitudes = new float[solution->numPoints*solution->numGeometries];
	//	float* mStarts = new float[solution->numPoints*16];
	//	float* mEnds = new float[solution->numPoints*16];
	//	err = clEnqueueReadBuffer(queue, solution->pointsReIm_CLBuffer, CL_TRUE, 0, solution->numPoints*2* sizeof(float), initialGuess, 1, &(solution->events[GSPAT_Event::POSITIONS_UPLOADED]), NULL);
	//	err = clEnqueueReadBuffer(queue, solution->positions_CLBuffer, CL_TRUE, 0, solution->numPoints*solution->numGeometries * 4* sizeof(float), positions, 1, &(solution->events[GSPAT_Event::POSITIONS_UPLOADED]), NULL);
	//	err = clEnqueueReadBuffer(queue, solution->targetAmplitudes_CLBuffer, CL_TRUE, 0, solution->numPoints*solution->numGeometries *  sizeof(float), amplitudes, 1, &(solution->events[GSPAT_Event::POSITIONS_UPLOADED]), NULL);
	//	err = clEnqueueReadBuffer(queue, solution->matrixStart_CLBuffer, CL_TRUE, 0, solution->numPoints*16* sizeof(float), mStarts, 1, &(solution->events[GSPAT_Event::POSITIONS_UPLOADED]), NULL);
	//	err = clEnqueueReadBuffer(queue, solution->matrixEnd_CLBuffer, CL_TRUE, 0, solution->numPoints*16*sizeof(float), mEnds, 1, &(solution->events[GSPAT_Event::POSITIONS_UPLOADED]), NULL);
	//	
	//		std::ostringstream msg;
	//		msg << "GS-PAT input(Z=" << solution->numPoints << ", G=" << solution->numGeometries << "):\n";
	//		msg << "   initialGuess :";
	//		for (int i = 0; i < solution->numPoints; i++)msg << "(" << initialGuess[2 * i] << "," << initialGuess[2 * i + 1] << "), ";
	//		msg << "\n   positions (";
	//		for (int i = 0; i < solution->numPoints*solution->numGeometries * 4; i++)
	//			msg << positions[i] << ";"<<(i%4==3?"    ":"")<< (i%4*solution->numPoints==4*solution->numPoints-1?"\n":""); 

	//		msg << ")\n   amplitudes (";
	//		for (int i = 0; i < solution->numPoints*solution->numGeometries * 1; i++)
	//			msg << amplitudes[i] << ";";// << (i%solution->numPoints == solution->numPoints - 1 ? "\n" : ""); msg << ")\n";
	//		msg << "   mStarts:\n";
	//		for (int i = 0; i < solution->numPoints; i++) {
	//			msg << "\t M" << i << "[";
	//			for (int m = 0; m < 16; m++) msg << mStarts[16 * i + m] << ", ";
	//			msg << "]\n";
	//		}
	//		msg << "   mEnds:\n";
	//		for (int i = 0; i < solution->numPoints; i++) {
	//			msg << "\t M" << i << "[";
	//			for (int m = 0; m < 16; m++) msg << mEnds[16 * i + m] << ", ";
	//			msg << "]\n";
	//		}
	//		delete positions; delete amplitudes; delete mStarts; delete mEnds;
	//		GSPAT::printMessage_GSPAT(msg.str().c_str());
	//	
	//
	//}
	//DEBUG:
	//float checkPointHologram[1024 * 2 *2 ];
	//err = clEnqueueReadBuffer(queue, solution->singlePointField_CLBuffer, CL_TRUE, 0, 1024*2*2* sizeof(float), checkPointHologram, 1, &(solution->events[GSPAT_Event::F_AND_B_READY]), NULL);
	//printf("Check point hologram.");
	//END DEBUG
}

void HologramSolverCL::solvePhases_IBP(HologramSolution* solution){
	cl_int err;
	size_t global_size[2] = { (size_t)numTransducers, (size_t)solution->numGeometries }, local_size[2] = { (size_t)numTransducers, 1 };

	//1. Set arguments
	err = clSetKernelArg(solveIBPKernel, 0, sizeof(int), &(solution->numPoints));
	if (err < 0) { sprintf(consoleLine,"GSPAT: solvePhases_IBP::Couldn't set kernel argument 0");GSPAT_IBP::printWarning_GSPAT(consoleLine); return; }
	err = clSetKernelArg(solveIBPKernel, 1, sizeof(cl_mem), &solution->pointsReIm_CLBuffer);
	if (err < 0) { sprintf(consoleLine,"GSPAT: solvePhases_IBP::Couldn't set kernel argument 1"); GSPAT_IBP::printWarning_GSPAT(consoleLine); return; }
	err = clSetKernelArg(solveIBPKernel, 2, sizeof(cl_mem), &solution->targetAmplitudes_CLBuffer);
	if (err < 0) { sprintf(consoleLine,"GSPAT: solvePhases_IBP::Couldn't set kernel argument 2"); GSPAT_IBP::printWarning_GSPAT(consoleLine); return; }
	//err = clSetKernelArg(solveIBPKernel, 3, sizeof(cl_mem), &solution->singlePointField_CLBuffer);
	err = clSetKernelArg(solveIBPKernel, 3, sizeof(cl_mem), &solution->singlePointFieldNormalised_CLBuffer);
	if (err < 0) { sprintf(consoleLine,"GSPAT: solvePhases_IBP::Couldn't set kernel argument 3"); GSPAT_IBP::printWarning_GSPAT(consoleLine); return; }
	err = clSetKernelArg(solveIBPKernel, 4, sizeof(cl_mem), &solution->finalHologramPhases_CLBuffer);
	if (err < 0) { sprintf(consoleLine,"GSPAT: solvePhases_IBP::Couldn't set kernel argument 4"); GSPAT_IBP::printWarning_GSPAT(consoleLine);return; }
	err = clSetKernelArg(solveIBPKernel, 5, sizeof(cl_mem), &solution->finalHologramAmplitudes_CLBuffer);
	if (err < 0) { sprintf(consoleLine,"GSPAT: solvePhases_IBP::Couldn't set kernel argument 5"); GSPAT_IBP::printWarning_GSPAT(consoleLine);return; }
	err = clSetKernelArg(solveIBPKernel, 6, sizeof(cl_mem), &solution->finalHologramReIm_CLBuffer);
	if (err < 0) { sprintf(consoleLine,"GSPAT: solvePhases_IBP::Couldn't set kernel argument 6"); GSPAT_IBP::printWarning_GSPAT(consoleLine);return; }

	//2. Group events to wait for (initialGuessUploaded + amplitudesUploaded+ propagationMatrixReady [numGeometries])
	
	//3. Run kernel with arguments and events:
	err |= clEnqueueNDRangeKernel(queue, solveIBPKernel, 2, NULL, global_size, local_size, 1, &(solution->events[GSPAT_Event::F_AND_B_READY]), &(solution->events[GSPAT_Event::POINT_PHASES_READY]));
	if (err < 0) { 
		GSPAT_IBP::printError_GSPAT("GSPAT: solvePhases_GS:: Couldn't enqueue the kernel"); return; }	
	//DEBUG
	//float pointsReIm_read[32 * 32 * 2];
	//float corrections[32];
	//float estimatedAmplitudes[32 * 32];
	//clEnqueueReadBuffer(queue, solution->pointsReIm_CLBuffer, CL_TRUE, 0, solution->numGeometries*solution->numPoints * 2 * sizeof(float), pointsReIm_read, 0, NULL, NULL);
	//clEnqueueReadBuffer(queue, solution->correction, CL_TRUE, 0, solution->numGeometries * sizeof(float), corrections, 0, NULL, NULL);
	//clEnqueueReadBuffer(queue, solution->amplitudesPerPoint, CL_TRUE, 0, solution->numGeometries*solution->numPoints * sizeof(float), estimatedAmplitudes, 0, NULL, NULL);
	//sprintf(consoleLine,"Hi!\n");
	float phases[512 * 32], amplitudes[512 * 32];
	err = clEnqueueReadBuffer(queue, solution->finalHologramPhases_CLBuffer, CL_TRUE, 0, numTransducers * solution->numGeometries * sizeof(float), phases, 1, &(solution->events[GSPAT_Event::POINT_PHASES_READY]), NULL);
	err = clEnqueueReadBuffer(queue, solution->finalHologramAmplitudes_CLBuffer, CL_TRUE, 0, numTransducers * solution->numGeometries * sizeof(float), amplitudes, 1, &(solution->events[GSPAT_Event::POINT_PHASES_READY]), NULL);
	
	if (err < 0) {
		GSPAT_IBP::printError_GSPAT("GSPAT: discretise::Couldn't read message queue"); return;
	}		
	//END DEBUG

}

void HologramSolverCL::discretise(HologramSolution* solution) {
	cl_int err;
	{
		size_t global_size[3] = { (size_t)solution->numTransducers, 1,(size_t)solution->numGeometries }, local_size[3] = { (size_t)numTransducers,1,1 };
		float solvePhaseOnly = solution->phaseOnly ? 1.0f : 0.0f;
		err = clSetKernelArg(discretiseKernel, 0, sizeof(int), &numDiscreteLevels);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 0"); GSPAT_IBP::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(discretiseKernel, 1, sizeof(cl_mem), &solution->finalHologramPhases_CLBuffer);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 1"); GSPAT_IBP::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(discretiseKernel, 2, sizeof(cl_mem), &solution->finalHologramAmplitudes_CLBuffer);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 2"); GSPAT_IBP::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(discretiseKernel, 3, sizeof(float), &solvePhaseOnly);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 3"); GSPAT_IBP::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(discretiseKernel, 4, sizeof(cl_mem), &phaseCorrections);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 4"); GSPAT_IBP::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(discretiseKernel, 5, sizeof(cl_mem), &amplitudeCorrections);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 5"); GSPAT_IBP::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(discretiseKernel, 6, sizeof(cl_mem), &transducerMappings);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 6"); GSPAT_IBP::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(discretiseKernel, 7, sizeof(cl_mem), &solution->messagesArray_CLBuffer);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 7"); GSPAT_IBP::printMessage_GSPAT(consoleLine);return; }		
		
		err = clEnqueueNDRangeKernel(queue, discretiseKernel, 3, NULL, global_size, local_size, 1, &(solution->events[GSPAT_Event::POINT_PHASES_READY]), &(solution->events[GSPAT_Event::HOLOGRAM_DISCRETISED]));

		if (err < 0) {
			GSPAT_IBP::printError_GSPAT("GSPAT: discretise::Couldn't enqueue the kernel"); return;
		}
		//Read top and bottom messages:
		err = clEnqueueReadBuffer(queue, solution->messagesArray_CLBuffer, CL_FALSE, 0, 2*numTransducers * solution->numGeometries * sizeof(unsigned char), &(solution->messagesArray[0]), 1, &(solution->events[GSPAT_Event::HOLOGRAM_DISCRETISED]), &(solution->events[GSPAT_Event::MESSAGES_READY]));
		if (err < 0) {
			GSPAT_IBP::printError_GSPAT("GSPAT: discretise::Couldn't read message queue"); return;
		}		
	}
}





