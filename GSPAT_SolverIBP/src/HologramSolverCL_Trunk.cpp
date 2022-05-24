#include <GSPAT_Solver.h>
#include <src/HologramSolverCL.h>
#include <../Helper/HelperMethods.h>
#include <src/OpenCLUtilityFunctions.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <clBLAS.h>
#include <glm/glm.hpp>						//GLM (Maths)
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <include\OpenCLSolverImpl_Interoperability.h>

char consoleLine[512];
static void printConsole(const char* str) {
	printf("%s", str);

}

static void printError(const char* str) {
	perror(str);
}

HologramSolverCL::HologramSolverCL(int boardWidth, int boardHeight, float pitch)
	: pitch(pitch), configured(false)
{ 
	//Registers default functions (in case the client did not register anything. It will be ignored it functions already set).
	GSPAT::RegisterPrintFuncs(printConsole, printConsole, printError);

	this->boardSize[0] = boardWidth;
	this->boardSize[1] = boardHeight;
	// Initializes the OpenCL platform
	initializeOpenCL();
	//3. Fill in initial guess for phase/amplitude...(create a method/class specifically for this?)
	for (int i = 0; i < HologramSolution::MAX_POINTS; i++) {
		initialGuess[2 * i] = 1.0f*i / HologramSolution::MAX_POINTS;
		initialGuess[2 * i + 1] = sqrtf(1 - initialGuess[2 * i] * initialGuess[2 * i]);
	}
	//Initialize global resources (shared by all solutions)
	createDirectivityTexture();
	createTransducerPositionBuffer(pitch, boardSize);	//1. create buffer with transducer positions:
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
	GSPAT::printMessage_GSPAT("GSPAT: Thread pool succesfully destroyed");
	//2. Delete buffers
	cl_uint count;
	if (configured) {
		clGetMemObjectInfo(transducerMappings, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
		clReleaseMemObject(transducerMappings);
		clGetMemObjectInfo(phaseCorrections, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
		clReleaseMemObject(phaseCorrections);
		GSPAT::printMessage_GSPAT("GSPAT: Board config (mappings and phase delays) succesfully destroyed");
		configured = false;
	}
	clGetMemObjectInfo(directivityTexture, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(directivityTexture);
	clGetMemObjectInfo(transducerPositions, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(transducerPositions);
	//3. Delete kernels
	clReleaseKernel(fillAllPointHologramsKernel);
	clReleaseKernel(powerMethodKernel);
	clReleaseKernel(addHologramsKernel);
	clReleaseKernel(discretiseKernel);
	GSPAT::printMessage_GSPAT("GSPAT: Kernels succesfully destroyed");
	//4. Delete programs
	clReleaseProgram(program);
	GSPAT::printMessage_GSPAT("GSPAT: Program succesfully destroyed");
	//5. Delete clBlas, queue and context
	clblasTeardown();
	GSPAT::printMessage_GSPAT("GSPAT: clBlas succesfully destroyed");
	//Delete global context
	clReleaseCommandQueue(queue);
	clReleaseContext(context);
	GSPAT::printMessage_GSPAT("GSPAT: Context destroyed");
	
}

void HologramSolverCL::setBoardConfig(int transducerIds[512], int phaseAdjust[512], int numLevels ) {
	//Unload previous configuraiton
	if (configured) {
		clReleaseMemObject(transducerMappings);
		clReleaseMemObject(phaseCorrections);
		configured = false;
		GSPAT::printMessage_GSPAT("GSPAT: Discarding previous configuration (transducer IDs, phase adjustments)\n");

	}
	unsigned char t[512];//The mapping should be stored in bytes (not integers)
	float p[512];//We store phase adjusts in radians
	for (int i = 0; i < 512; i++) {
		t[i] = (unsigned char)(transducerIds[i] - (i>256? 256 : 0));
		p[i] = (float)( phaseAdjust[i]* M_PI / 180.0f);
	}
	transducerMappings= clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, 512* sizeof(unsigned char), t, NULL);	
	phaseCorrections= clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, 512* sizeof(float), p, NULL);	
	numDiscreteLevels = numLevels;
	GSPAT::printMessage_GSPAT("GSPAT: Board configured (transducerIDs, phase Adjustments)\n");
	configured = true;
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
	OpenCLUtilityFunctions::createProgramFromFile("hologramSolver2.cl", context, deviceUsed, &program);
	cl_int err;
	fillAllPointHologramsKernel = clCreateKernel(program, "computeFandB", &err);
    if(err < 0) {
		GSPAT::printError_GSPAT("GSPAT: Couldn't create kernel\n");
    }
	powerMethodKernel = clCreateKernel(program, "solvePhases_GS", &err);
	if(err < 0) {
		sprintf(consoleLine,"GSPAT: Couldn't create a kernel: %d", err);
		GSPAT::printError_GSPAT(consoleLine);
	}
	addHologramsKernel = clCreateKernel(program, "computeActivation", &err);
	if(err < 0) {
		sprintf(consoleLine,"GSPAT: Couldn't create a kernel: %d", err);
		GSPAT::printError_GSPAT(consoleLine);
	}; 
	discretiseKernel=clCreateKernel(program, "discretise", &err);
	if (err < 0) {
		sprintf(consoleLine,"GSPAT: Couldn't create a kernel: %d", err);
		GSPAT::printError_GSPAT(consoleLine);
	};
	//Initialise clBLAS
	err = clblasSetup();
	if (err != CL_SUCCESS) {
		sprintf(consoleLine,"GSPAT: clblasSetup() failed with %d\n", err);
		GSPAT::printError_GSPAT(consoleLine);
	}
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
		GSPAT::printMessage_GSPAT("GSPAT: Directivity bufffers created successfully\n");
	}
	else 
		GSPAT::printError_GSPAT("GSPAT: Could not create directivity buffers\n");
}

void HologramSolverCL:: createTransducerPositionBuffer(float pitch, int boardSize[2]) {
		size_t numElementsInBuffer = boardSize[0] * boardSize[1] * 4;
		float* positions = new float[numElementsInBuffer];//float positions [32*16*4];
		//1. Compute the positions in main memory
		int t_index[2] = { 0,0 };
		for (t_index[1] = 0; t_index[1] < boardSize[1]; t_index[1]++)
			for(t_index[0]=0;t_index[0]<boardSize[0];t_index[0]++) {
				size_t offset = (t_index[0] + boardSize[0] * t_index[1])*4;
				//fill in (x,y,z)
				//computeTransducerPos(t_index, boardSize, pitch, (float*)(positions+offset));
				computeTransducerPos_SideBySide(t_index, pitch, (float*)(positions+offset));
				//fill in homogenenous coordinate;
				positions[offset + 3] = 1;
			}
		//2. Create an OpenCL buffer and store positions:
		cl_int err;
		transducerPositions = clCreateBuffer(context, CL_MEM_READ_WRITE, numElementsInBuffer * sizeof(float), NULL, &err);
		if (err < 0) { GSPAT::printError_GSPAT("GSPAT: Couldn't create 'transducerPositions' buffer"); return ; }
		err = clEnqueueWriteBuffer(queue, transducerPositions, CL_TRUE, 0, 32 * 16 * 4 * sizeof(float), &(positions[0]), 0, NULL, NULL);
		if (err < 0) { GSPAT::printError_GSPAT("GSPAT: Couldn't write 'transducerPositions' buffer"); return ; }
		delete positions;
		//BEGIN DEBUG
		//static float result [32 * 16 * 4 ];//Check it was written...
		// err = clEnqueueReadBuffer(queue, transducerPositions, CL_TRUE, 0, 32*16*4*sizeof(float), &(result[0]),0, NULL, NULL);
		//if (err < 0) { perror("Couldn't read hologram"); return NULL; }
		//END DEBUG
	}

void HologramSolverCL::createSolutionPool() {
	solutionPool = new SolutionPool(this, boardSize[0], boardSize[1], pitch, context, queue);
}

Solution* HologramSolverCL::createSolution(int numPoints, int numGeometries, bool phaseOnly, float* positions, float* amplitudes, float* matStarts, float* matEnds, GSPAT::MatrixAlignment a) {
	HologramSolution* nextSolution; 
	//1. Access next available solution from the pool (and remove it from the list of "availables")
	nextSolution = solutionPool->getNextSolution();	
	
	//2. Adjust input data (alignment layout).
	float mStart[16 * 32], mEnd[16 * 32];	
	if (a == GSPAT::ColumnMajorAlignment) {
		//glm matrix are column major order... we need to transpose them
		for (int i = 0; i < numPoints; i++) {
			unsigned int offset = 16 * i;
			mStart[offset + 0] = matStarts[offset+4*0+0];		mStart[offset + 1] = matStarts[offset + 4 * 1 + 0]; mStart[offset + 2] = matStarts[offset + 4 * 2 + 0]; mStart[offset + 3] = matStarts[offset + 4 * 3 + 0];
			mStart[offset + 4] = matStarts[offset + 4 * 0 + 1]; mStart[offset + 5] = matStarts[offset + 4 * 1 + 1]; mStart[offset + 6] = matStarts[offset + 4 * 2 + 1]; mStart[offset + 7] = matStarts[offset + 4 * 3 + 1];
			mStart[offset + 8] = matStarts[offset + 4 * 0 + 2]; mStart[offset + 9] = matStarts[offset + 4 * 1 + 2]; mStart[offset + 10]= matStarts[offset + 4 * 2 + 2]; mStart[offset + 11]= matStarts[offset + 4 * 3 + 2];
			mStart[offset +12] = matStarts[offset + 4 * 0 + 3]; mStart[offset + 13]= matStarts[offset + 4 * 1 + 3]; mStart[offset + 14]= matStarts[offset + 4 * 2 + 3]; mStart[offset + 15]= matStarts[offset + 4 * 3 + 3];

			mEnd[offset + 0] = matEnds[offset + 4 * 0 + 0];	mEnd[offset + 1] = matEnds[offset + 4 * 1 + 0]; mEnd[offset + 2] = matEnds[offset + 4 * 2 + 0]; mEnd[offset + 3] = matEnds[offset + 4 * 3 + 0];
			mEnd[offset + 4] = matEnds[offset + 4 * 0 + 1]; mEnd[offset + 5] = matEnds[offset + 4 * 1 + 1]; mEnd[offset + 6] = matEnds[offset + 4 * 2 + 1]; mEnd[offset + 7] = matEnds[offset + 4 * 3 + 1];
			mEnd[offset + 8] = matEnds[offset + 4 * 0 + 2]; mEnd[offset + 9] = matEnds[offset + 4 * 1 + 2]; mEnd[offset + 10]= matEnds[offset + 4 * 2 + 2]; mEnd[offset + 11]= matEnds[offset + 4 * 3 + 2];
			mEnd[offset + 12]= matEnds[offset + 4 * 0 + 3]; mEnd[offset + 13]= matEnds[offset + 4 * 1 + 3]; mEnd[offset + 14]= matEnds[offset + 4 * 2 + 3]; mEnd[offset + 15]= matEnds[offset + 4 * 3 + 3];

			/*mEnd[offset + 0] = matEnds[i][0][0]; mEnd[offset + 1] = matEnds[i][1][0]; mEnd[offset + 2] = matEnds[i][2][0]; mEnd[offset + 3] = matEnds[i][3][0];
			mEnd[offset + 4] = matEnds[i][0][1]; mEnd[offset + 5] = matEnds[i][1][1]; mEnd[offset + 6] = matEnds[i][2][1]; mEnd[offset + 7] = matEnds[i][3][1];
			mEnd[offset + 8] = matEnds[i][0][2]; mEnd[offset + 9] = matEnds[i][1][2]; mEnd[offset + 10] = matEnds[i][2][2]; mEnd[offset + 11] = matEnds[i][3][2];
			mEnd[offset + 12] = matEnds[i][0][3]; mEnd[offset + 13] = matEnds[i][1][3]; mEnd[offset + 14] = matEnds[i][2][3]; mEnd[offset + 15] = matEnds[i][3][3];*/
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


void HologramSolverCL::releaseSolution(Solution* solution) {
	((HologramSolution*)solution)->releaseEvents();
	solutionPool->returnSolution((HologramSolution*)solution);

}


void HologramSolverCL::compute(Solution* solution) {
	HologramSolution* hs = ((HologramSolution*)solution);
	updateCLBuffers(hs);
	computeFandB(hs);
	computeR(hs,true);
	solvePhases_GS(hs,true);
	computeActivation(hs);
	discretise(hs);
	clFlush(queue);
}

void HologramSolverCL::updateCLBuffers(HologramSolution * hs)
{
	cl_int err;
	cl_event positions_copied, mStart_copied, mEnd_copied;
	hs->lockEvents();	
	if (hs->manualData) {
		err = clEnqueueWriteBuffer(queue, hs->positions_CLBuffer, CL_FALSE, NULL, 4 * hs->numPoints *hs->numGeometries * sizeof(float), hs->positions, 0, NULL, &(hs->events[GSPAT_Event::POSITIONS_UPLOADED]));
		err = clEnqueueWriteBuffer(queue, hs->targetAmplitudes_CLBuffer, CL_FALSE, NULL, hs->numPoints* hs->numGeometries * sizeof(float), hs->amplitudes, 0, NULL, &(hs->events[GSPAT_Event::AMPLITUDES_UPLOADED]));
		err = clEnqueueWriteBuffer(queue, hs->pointsReIm_CLBuffer, CL_FALSE, 0, 2 * hs->points() * sizeof(float), initialGuess, 0, NULL, &(hs->events[GSPAT_Event::INITIAL_GUESS_UPLOADED]));
		err = clEnqueueWriteBuffer(queue, hs->matrixStart_CLBuffer, CL_FALSE, NULL, 16 * (hs->numPoints) * sizeof(float), hs->matrixStart, 0, NULL, &(hs->events[GSPAT_Event::MATRIX_0_UPLOADED]));
		err = clEnqueueWriteBuffer(queue, hs->matrixEnd_CLBuffer, CL_FALSE, NULL, 16 * (hs->numPoints) * sizeof(float), hs->matrixEnd, 0, NULL, &(hs->events[GSPAT_Event::MATRIX_G_UPLOADED]));
		if (err < 0) { sprintf(consoleLine, "GSPAT::updateCLBuffers::Couldn't copy input parameters to the Solution"); GSPAT::printWarning_GSPAT(consoleLine); return; }
	}
}

/*void  HologramSolverCL::computeWithExternalPropagators(float* pointHologram, float* normalizedPointHologram, Solution* s) {
	HologramSolution* solution = (HologramSolution*)s;
	//0.Copy input parameters: 
	//1. Copy progagators: 
	clEnqueueWriteBuffer(queue, solution->singlePointField_CLBuffer, CL_TRUE, 0, solution->numGeometries*solution->numPoints*512 * 2 * sizeof(float),pointHologram,0, NULL, &(solution->pointHologramsReady[0]));
	clEnqueueWriteBuffer(queue, solution->singlePointFieldNormalised_CLBuffer, CL_TRUE, 0, solution->numGeometries*solution->numPoints*512 * 2 * sizeof(float),normalizedPointHologram,0, NULL, &(solution->pointHologramsReady[1]));
	clEnqueueWriteBuffer(queue, solution->singlePointField_CLBuffer, CL_TRUE, 0, solution->numGeometries*solution->numPoints*512 * 2 * sizeof(float),pointHologram,0, NULL, &(solution->pointHologramsReady[2]));
	clEnqueueWriteBuffer(queue, solution->singlePointFieldNormalised_CLBuffer, CL_TRUE, 0, solution->numGeometries*solution->numPoints*512 * 2 * sizeof(float),normalizedPointHologram,0, NULL, &(solution->pointHologramsReady[3]));
	
	computePropagationMatrix_OpenCL(solution,true);
	//profiler.recordEndPhaseData();
	solvePhasesPower(solution,true);
	//profiler.recordEndPhaseSolved();
	computeSolutionAddition(solution);
	//profiler.recordEndCompute();
	discretiseSolution(solution);
}*/


void HologramSolverCL::computeFandB(HologramSolution* solution) {
	cl_int err;
	//B. Run the kernel
	size_t global_size[3] = { boardSize[0], boardSize[1], solution->numPoints*solution->numGeometries};
	size_t local_size[3] =  { boardSize[0], boardSize[1], 1};
	//0. Setup inputs for the kernell (do once)
	err = clSetKernelArg(fillAllPointHologramsKernel, 0, sizeof(cl_mem), &(transducerPositions));
	if (err < 0) { GSPAT::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 0"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 1, sizeof(cl_mem), &(solution->positions_CLBuffer));
	if (err < 0) { GSPAT::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 1"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 2, sizeof(cl_mem), &(solution->matrixStart_CLBuffer));
	if (err < 0) { GSPAT::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 2"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 3, sizeof(cl_mem), &(solution->matrixEnd_CLBuffer));
	if (err < 0) { GSPAT::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 3"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 4, sizeof(int), &(solution->numPoints));
	if (err < 0) { GSPAT::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 4"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 5, sizeof(int), &(solution->numGeometries));
	if (err < 0) { GSPAT::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 5"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 6, sizeof(cl_mem), &(this->directivityTexture));
	if (err < 0) { GSPAT::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 6"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 7, sizeof(cl_mem), &(solution->singlePointField_CLBuffer));
	if (err < 0) { GSPAT::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 7"); return ; }
	err = clSetKernelArg(fillAllPointHologramsKernel, 8, sizeof(cl_mem), &(solution->singlePointFieldNormalised_CLBuffer));
	if (err < 0) { GSPAT::printWarning_GSPAT("GSPAT: computeFandB::Couldn't set kernel argument 8"); return ; }
		
	//1.3. Trigger the kernell
	cl_event dataUploaded[] = {solution->events[GSPAT_Event::POSITIONS_UPLOADED], solution->events[GSPAT_Event::MATRIX_0_UPLOADED], solution->events[GSPAT_Event::MATRIX_G_UPLOADED]};
	cl_uint numEvents = solution->manualData ? 3 : 1;
	err |= clEnqueueNDRangeKernel(queue, fillAllPointHologramsKernel, 3, NULL, global_size, local_size, numEvents, dataUploaded, &(solution->events[GSPAT_Event::F_AND_B_READY]));
	if (err < 0) { 
		GSPAT::printError_GSPAT("GSPAT: computeFandB:: Couldn't enqueue the kernel"); return ; 
	}	

	{//DEBUG: Print input data: 
	/*	float *initialGuess = new float[solution->numPoints * 2];
		float* positions = new float[solution->numPoints*solution->numGeometries * 4];
		float* amplitudes = new float[solution->numPoints*solution->numGeometries];
		float* mStarts = new float[solution->numPoints*16];
		float* mEnds = new float[solution->numPoints*16];
		err = clEnqueueReadBuffer(queue, solution->pointsReIm_CLBuffer, CL_TRUE, 0, solution->numPoints*2* sizeof(float), initialGuess, 1, &(solution->events[GSPAT_Event::POSITIONS_UPLOADED]), NULL);
		err = clEnqueueReadBuffer(queue, solution->positions_CLBuffer, CL_TRUE, 0, solution->numPoints*solution->numGeometries * 4* sizeof(float), positions, 1, &(solution->events[GSPAT_Event::POSITIONS_UPLOADED]), NULL);
		err = clEnqueueReadBuffer(queue, solution->targetAmplitudes_CLBuffer, CL_TRUE, 0, solution->numPoints*solution->numGeometries *  sizeof(float), amplitudes, 1, &(solution->events[GSPAT_Event::POSITIONS_UPLOADED]), NULL);
		err = clEnqueueReadBuffer(queue, solution->matrixStart_CLBuffer, CL_TRUE, 0, solution->numPoints*16* sizeof(float), mStarts, 1, &(solution->events[GSPAT_Event::POSITIONS_UPLOADED]), NULL);
		err = clEnqueueReadBuffer(queue, solution->matrixEnd_CLBuffer, CL_TRUE, 0, solution->numPoints*16*sizeof(float), mEnds, 1, &(solution->events[GSPAT_Event::POSITIONS_UPLOADED]), NULL);
		
			std::ostringstream msg;
			msg << "GS-PAT input(Z=" << solution->numPoints << ", G=" << solution->numGeometries << "):\n";
			msg << "   initialGuess :";
			for (int i = 0; i < solution->numPoints; i++)msg << "(" << initialGuess[2 * i] << "," << initialGuess[2 * i + 1] << "), ";
			msg << "\n   positions (";
			for (int i = 0; i < solution->numPoints*solution->numGeometries * 4; i++)
				msg << positions[i] << ", "; msg << ")\n";
			msg << "   amplitudes (";
			for (int i = 0; i < solution->numPoints*solution->numGeometries * 1; i++)
				msg << amplitudes[i] << ", "; msg << ")\n";
			msg << "   mStarts:\n";
			for (int i = 0; i < solution->numPoints; i++) {
				msg << "\t M" << i << "[";
				for (int m = 0; m < 16; m++) msg << mStarts[16 * i + m] << ", ";
				msg << "]\n";
			}
			msg << "   mEnds:\n";
			for (int i = 0; i < solution->numPoints; i++) {
				msg << "\t M" << i << "[";
				for (int m = 0; m < 16; m++) msg << mEnds[16 * i + m] << ", ";
				msg << "]\n";
			}
			delete positions; delete amplitudes; delete mStarts; delete mEnds;
			GSPAT::printMessage_GSPAT(msg.str().c_str());
		
	*/
	}
}

void HologramSolverCL::computeR(HologramSolution* solution,bool useDiego) {
	//1. Select sources to use (point holograms or normalised point holograms)
	cl_mem cl_hologramsToUse;	
	if (useDiego)
		cl_hologramsToUse = solution->singlePointFieldNormalised_CLBuffer;
	else
		cl_hologramsToUse = solution->singlePointField_CLBuffer;
	
	//3. compute
	cl_event event = NULL;
	cl_float2 alpha; alpha.x = 1; alpha.y = 0;
	cl_float2 beta; beta.x = beta.y = 0;
	int solutionSize = solution->boardSize[0] * solution->boardSize[1] * solution->numPoints
		, R_matrix_size = solution->numPoints * solution->numPoints;
	for (int g = 0; g < solution->numGeometries; g++) {
		int retCode = clblasCgemm(clblasRowMajor,
			clblasNoTrans, clblasConjTrans,
			solution->numPoints, solution->numPoints, solution->boardSize[0] * solution->boardSize[1],
			alpha,
			solution->singlePointField_CLBuffer, solutionSize*g, solution->boardSize[0] * solution->boardSize[1],
			cl_hologramsToUse, solutionSize*g, solution->boardSize[0] * solution->boardSize[1],
			beta,
			solution->R_CLBuffer, R_matrix_size*g, solution->numPoints,
			1, &queue, 1, &(solution->events[GSPAT_Event::F_AND_B_READY]), &(solution->events[GSPAT_Event::R0_READY + g]));	// this line is different from the previous one
		
		if (retCode != clblasStatus::clblasSuccess) {
			sprintf(consoleLine, "GSPAT: computeR::CL_BLAS error: %d", retCode); GSPAT::printError_GSPAT(consoleLine);
		}
	}
}

void HologramSolverCL::solvePhases_GS(HologramSolution* solution,bool useNorm){
	cl_int err;
	size_t global_size[2] = { solution->numPoints, solution->numGeometries }, local_size[2] = { solution->numPoints, 1 };

	//1. Set arguments
	err = clSetKernelArg(powerMethodKernel, 0, sizeof(cl_mem), &solution->pointsReIm_CLBuffer);
	if (err < 0) { sprintf(consoleLine,"GSPAT: solvePhases_GS::Couldn't set kernel argument 0"); GSPAT::printWarning_GSPAT(consoleLine); return; }
	err = clSetKernelArg(powerMethodKernel, 1, sizeof(cl_mem), &solution->R_CLBuffer);
	if (err < 0) { sprintf(consoleLine,"GSPAT: solvePhases_GS::Couldn't set kernel argument 1"); GSPAT::printWarning_GSPAT(consoleLine); return; }
	err = clSetKernelArg(powerMethodKernel, 2, sizeof(cl_mem), &solution->targetAmplitudes_CLBuffer);
	if (err < 0) { sprintf(consoleLine,"GSPAT: solvePhases_GS::Couldn't set kernel argument 2"); GSPAT::printWarning_GSPAT(consoleLine); return; }
	err = clSetKernelArg(powerMethodKernel, 3, sizeof(cl_mem), &solution->amplitudesPerPoint);
	if (err < 0) { sprintf(consoleLine,"GSPAT: solvePhases_GS::Couldn't set kernel argument 3"); GSPAT::printWarning_GSPAT(consoleLine); return; }
	err = clSetKernelArg(powerMethodKernel, 4, sizeof(cl_mem), &solution->correction);
	if (err < 0) { sprintf(consoleLine,"GSPAT: solvePhases_GS::Couldn't set kernel argument 4"); GSPAT::printWarning_GSPAT(consoleLine); return; }
	//2. Group events to wait for (initialGuessUploaded + amplitudesUploaded+ propagationMatrixReady [numGeometries])
	cl_event* powerEvents = new cl_event[solution->numGeometries + 2];
	memcpy(powerEvents, &(solution->events[GSPAT_Event::R0_READY]), solution->numGeometries * sizeof(cl_event));
	powerEvents[solution->numGeometries] = solution->events[GSPAT_Event::INITIAL_GUESS_UPLOADED];
	powerEvents[solution->numGeometries+1] = solution->events[GSPAT_Event::AMPLITUDES_UPLOADED];
	//3. Run kernel with arguments and events:
	err |= clEnqueueNDRangeKernel(queue, powerMethodKernel, 2, NULL, global_size, local_size, solution->numGeometries+2*solution->manualData, powerEvents, &(solution->events[GSPAT_Event::POINT_PHASES_READY]));
	delete powerEvents;
	if (err < 0) { 
		GSPAT::printError_GSPAT("GSPAT: solvePhases_GS:: Couldn't enqueue the kernel"); return; }	
	//DEBUG
	//float pointsReIm_read[32 * 32 * 2];
	//float corrections[32];
	//float estimatedAmplitudes[32 * 32];
	//clEnqueueReadBuffer(queue, solution->pointsReIm_CLBuffer, CL_TRUE, 0, solution->numGeometries*solution->numPoints * 2 * sizeof(float), pointsReIm_read, 0, NULL, NULL);
	//clEnqueueReadBuffer(queue, solution->correction, CL_TRUE, 0, solution->numGeometries * sizeof(float), corrections, 0, NULL, NULL);
	//clEnqueueReadBuffer(queue, solution->amplitudesPerPoint, CL_TRUE, 0, solution->numGeometries*solution->numPoints * sizeof(float), estimatedAmplitudes, 0, NULL, NULL);
	//sprintf(consoleLine,"Hi!\n");
	//END DEBUG

}

void HologramSolverCL::computeActivation(HologramSolution* solution){
	cl_int err;
	//B. Build hologram by adding point holograms (and applying phase to each point)
	{
		size_t global_size[3] = { solution->boardSize[0], solution->boardSize[1],solution->numGeometries }, local_size[3] = { solution->boardSize[0],solution->boardSize[1],1 };
		size_t N = solution->boardSize[0] * solution->boardSize[1]*solution->numGeometries;
		err = clSetKernelArg(addHologramsKernel, 0, sizeof(int), &(solution->numPoints));
		if (err < 0) { sprintf(consoleLine,"GSPAT: computeActivation::Couldn't set kernel argument 0");GSPAT::printWarning_GSPAT(consoleLine); return; }
		err = clSetKernelArg(addHologramsKernel, 1, sizeof(cl_mem), &solution->singlePointFieldNormalised_CLBuffer);
		if (err < 0) { sprintf(consoleLine,"GSPAT: computeActivation::Couldn't set kernel argument 1"); GSPAT::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(addHologramsKernel, 2, sizeof(cl_mem), &solution->pointsReIm_CLBuffer);
		if (err < 0) { sprintf(consoleLine,"GSPAT: computeActivation::Couldn't set kernel argument 2"); GSPAT::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(addHologramsKernel, 3, sizeof(cl_mem), &solution->correction);
		if (err < 0) { sprintf(consoleLine,"GSPAT: computeActivation::Couldn't set kernel argument 3"); GSPAT::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(addHologramsKernel, 4, sizeof(cl_mem), &solution->finalHologramPhases_CLBuffer);
		if (err < 0) { sprintf(consoleLine,"GSPAT: computeActivation::Couldn't set kernel argument 4"); GSPAT::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(addHologramsKernel, 5, sizeof(cl_mem), &solution->finalHologramAmplitudes_CLBuffer);
		if (err < 0) { sprintf(consoleLine,"GSPAT: computeActivation::Couldn't set kernel argument 5"); GSPAT::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(addHologramsKernel, 6, sizeof(cl_mem), &solution->finalHologramReIm_CLBuffer);
		if (err < 0) { sprintf(consoleLine,"GSPAT: computeActivation::Couldn't set kernel argument 6"); GSPAT::printWarning_GSPAT(consoleLine);return; }


		err |= clEnqueueNDRangeKernel(queue, addHologramsKernel, 3, NULL, global_size, local_size, 1, &(solution->events[GSPAT_Event::POINT_PHASES_READY]), &(solution->events[GSPAT_Event::HOLOGRAM_READY]));

		if (err < 0) {
			GSPAT::printError_GSPAT("GSPAT: computeActivation::Couldn't enqueue the kernel"); return;
		}
		
		//DEBUG: Enable only if you need the hologram with the original Re and Im parts. NOTE: You also need to modify ~HologramSolution, to release the hologramRead event
		/*err |= clEnqueueReadBuffer(queue, solution->finalHologramReIm_CLBuffer, CL_FALSE, 0, N * 2 * sizeof(float), &(solution->transducerReIm[0]), 1, &(solution->finalHologramComputed), &(solution->hologramRead));
		if (err < 0) {
			GSPAT::printError_GSPAT("Couldn't read final hologram"); return;
		}
		err |= clEnqueueReadBuffer(queue, solution->finalHologramPhases_CLBuffer, CL_FALSE, 0, N *  sizeof(float), &(solution->transducerPhases[0]), 1, &(solution->finalHologramComputed), &(solution->phasesAndAmplitudesRead[0]));
		//err |= clEnqueueReadBuffer(queue, solution->finalHologramPhases_CLBuffer, CL_TRUE, 0, N *  sizeof(float), &(solution->transducerPhases[0]), 1, &(solution->finalHologramComputed), NULL);
		if (err < 0) {
			GSPAT::printError_GSPAT("Couldn't read final hologram (phases)"); return;
		}
		err |= clEnqueueReadBuffer(queue, solution->finalHologramAmplitudes_CLBuffer, CL_FALSE, 0, N *  sizeof(float), &(solution->transducerAmplitudes[0]), 1, &(solution->finalHologramComputed), &(solution->phasesAndAmplitudesRead[1]));
		if (err < 0) {
			GSPAT::printError_GSPAT("Couldn't read final hologram (amplitudes)"); return;
		}*/
	}
}

void HologramSolverCL::discretise(HologramSolution* solution) {
	cl_int err;
	{
		size_t global_size[3] = { solution->boardSize[0], solution->boardSize[1],solution->numGeometries }, local_size[3] = { solution->boardSize[0],solution->boardSize[1],1 };
		int discreteLevels = 128;
		float solvePhaseOnly = solution->phaseOnly ? 1.0f : 0.0f;
		err = clSetKernelArg(discretiseKernel, 0, sizeof(int), &discreteLevels);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 0"); GSPAT::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(discretiseKernel, 1, sizeof(cl_mem), &solution->finalHologramPhases_CLBuffer);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 1"); GSPAT::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(discretiseKernel, 2, sizeof(cl_mem), &solution->finalHologramAmplitudes_CLBuffer);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 2"); GSPAT::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(discretiseKernel, 3, sizeof(float), &solvePhaseOnly);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 3"); GSPAT::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(discretiseKernel, 4, sizeof(cl_mem), &phaseCorrections);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 4"); GSPAT::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(discretiseKernel, 5, sizeof(cl_mem), &transducerMappings);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 5"); GSPAT::printWarning_GSPAT(consoleLine);return; }
		err = clSetKernelArg(discretiseKernel, 6, sizeof(cl_mem), &solution->messagesTopArray_CLBuffer);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 6"); GSPAT::printMessage_GSPAT(consoleLine);return; }		
		err = clSetKernelArg(discretiseKernel, 7, sizeof(cl_mem), &solution->messagesBottomArray_CLBuffer);
		if (err < 0) { sprintf(consoleLine,"GSPAT: discretise::Couldn't set kernel argument 7"); GSPAT::printWarning_GSPAT(consoleLine);return; }


		err = clEnqueueNDRangeKernel(queue, discretiseKernel, 3, NULL, global_size, local_size, 1, &(solution->events[GSPAT_Event::HOLOGRAM_READY]), &(solution->events[GSPAT_Event::HOLOGRAM_DISCRETISED]));

		if (err < 0) {
			GSPAT::printError_GSPAT("GSPAT: discretise::Couldn't enqueue the kernel"); return;
		}
		//Read top and bottom messages:
		err = clEnqueueReadBuffer(queue, solution->messagesBottomArray_CLBuffer, CL_FALSE, 0, 512 * solution->numGeometries * sizeof(unsigned char), &(solution->messagesBottomArray[0]), 1, &(solution->events[GSPAT_Event::HOLOGRAM_DISCRETISED]), &(solution->events[GSPAT_Event::MESSAGE_BOTTOM_READY]));
		if (err < 0) {
			GSPAT::printError_GSPAT("GSPAT: discretise::Couldn't read message queue (top array)"); return;
		}
		err = clEnqueueReadBuffer(queue, solution->messagesTopArray_CLBuffer, CL_FALSE, 0, 512 * solution->numGeometries * sizeof(unsigned char), &(solution->messagesTopArray[0]), 1, &(solution->events[GSPAT_Event::HOLOGRAM_DISCRETISED]), &(solution->events[GSPAT_Event::MESSAGE_TOP_READY]));
		if (err < 0) {
			GSPAT::printError_GSPAT("GSPAT: discretise::Couldn't read message queue (top array)"); return;
		}
	}
}





