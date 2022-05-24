#include <src/HologramSolution.h>
#include <GSPAT_SolverV2.h>
#include <iostream>

HologramSolution::HologramSolution(int numTransducers,  cl_context context, cl_command_queue& queue)
	:  numTransducers(numTransducers)
	, pointsReIm_CLBuffer(NULL)
	, positions_CLBuffer(NULL)
	, targetAmplitudes_CLBuffer(NULL)
	, matrixStart_CLBuffer(NULL)
	, matrixEnd_CLBuffer(NULL)
	, singlePointField_CLBuffer(NULL)
	, singlePointFieldNormalised_CLBuffer(NULL)
	, R_CLBuffer(NULL)
	, finalHologramPhases_CLBuffer(NULL)
	, finalHologramAmplitudes_CLBuffer(NULL)
	, finalHologramReIm_CLBuffer(NULL)
	, amplitudesPerPoint(NULL)
	, metric(NULL)
	, correction(NULL)
	, queue(queue)
	, eventsReleased(true)
	, manualData(true)
	, transducersOffSolution(false)
{
	//Initialise buffers:
	//1. Allocate memory dynamically (TODO: Required for transducer setups >512 transducers) :
	transducerReIm=new float[MAX_GEOMETRIES*numTransducers*2];									
	transducerPhases=new float[MAX_GEOMETRIES*numTransducers];									
	transducerAmplitudes=new float [MAX_GEOMETRIES*numTransducers];								
	messagesArray= new unsigned char[MAX_GEOMETRIES*numTransducers*2];
	//2. Prepare solution buffers to hold the focusing hologram
	positions_CLBuffer= clCreateBuffer(context, CL_MEM_READ_ONLY, MAX_GEOMETRIES* MAX_POINTS*4* sizeof(float), NULL, NULL);	
	targetAmplitudes_CLBuffer= clCreateBuffer(context, CL_MEM_READ_ONLY, MAX_GEOMETRIES* MAX_POINTS* sizeof(float), NULL, NULL);	
	matrixStart_CLBuffer= clCreateBuffer(context, CL_MEM_READ_ONLY, 16* MAX_POINTS * sizeof(float), NULL, NULL);
	matrixEnd_CLBuffer= clCreateBuffer(context, CL_MEM_READ_ONLY, 16* MAX_POINTS *sizeof(float), NULL, NULL);
	singlePointField_CLBuffer= clCreateBuffer(context, CL_MEM_READ_WRITE,  numTransducers * MAX_POINTS*2 *MAX_GEOMETRIES * sizeof(float), NULL, NULL);
	singlePointFieldNormalised_CLBuffer= clCreateBuffer(context, CL_MEM_READ_WRITE, numTransducers * MAX_POINTS*2 *MAX_GEOMETRIES* sizeof(float), NULL, NULL);	
	pointsReIm_CLBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, 2 * MAX_GEOMETRIES* MAX_POINTS * sizeof(float), NULL, NULL);
	R_CLBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE, 2 * MAX_POINTS*MAX_POINTS *MAX_GEOMETRIES* sizeof(float), NULL, NULL);
	finalHologramPhases_CLBuffer= clCreateBuffer(context, CL_MEM_READ_WRITE, MAX_GEOMETRIES*numTransducers * sizeof(float), NULL, NULL);
	finalHologramAmplitudes_CLBuffer= clCreateBuffer(context, CL_MEM_READ_WRITE, MAX_GEOMETRIES*numTransducers * sizeof(float), NULL, NULL);
	finalHologramReIm_CLBuffer =clCreateBuffer(context, CL_MEM_READ_WRITE, MAX_GEOMETRIES*numTransducers*2*sizeof(float), NULL, NULL);
	amplitudesPerPoint= clCreateBuffer(context, CL_MEM_READ_WRITE, MAX_POINTS*MAX_GEOMETRIES* sizeof(float), NULL, NULL);
	metric= clCreateBuffer(context, CL_MEM_READ_WRITE, 2*MAX_GEOMETRIES* sizeof(float), NULL, NULL);
	correction= clCreateBuffer(context, CL_MEM_READ_WRITE, MAX_GEOMETRIES* sizeof(float), NULL, NULL);
	messagesArray_CLBuffer =clCreateBuffer(context, CL_MEM_READ_WRITE, MAX_GEOMETRIES*numTransducers*2*sizeof(unsigned char), NULL, NULL);
	//messagesBottomArray_CLBuffer =clCreateBuffer(context, CL_MEM_READ_WRITE, MAX_GEOMETRIES*numTransducers*sizeof(unsigned char), NULL, NULL);
	
}

void HologramSolution::configureSolution(int numPoints, int numGeometries, bool phaseOnly, float* positions,float* amplitudes, float* mStart, float* mEnd) {
	this->numPoints = numPoints;
	this->numGeometries = numGeometries;
	this->phaseOnly = phaseOnly;
	memcpy(this->positions, positions, 4 * numPoints *numGeometries* sizeof(float));
	memcpy(this->amplitudes, amplitudes, numPoints *numGeometries* sizeof(float));
	memcpy(this->matrixStart, mStart, 16 * numPoints * sizeof(float));
	memcpy(this->matrixEnd, mEnd, 16 * numPoints * sizeof(float));
	manualData = true;
	transducersOffSolution = false;
}

void HologramSolution::configureSolutionExternalData(int numPoints, int numGeometries, bool phaseOnly) {
	this->numPoints = numPoints;
	this->numGeometries = numGeometries;
	this->phaseOnly = phaseOnly;
	manualData = false;
	transducersOffSolution = false;
}

void HologramSolution::configureSolutionTransducersOff(int numTransducers) {
	transducersOffSolution = true;
	eventsReleased = true;		//We did not allocate any events, so nothing to release when it returns to Solution pool
	numGeometries = 1;
	//Copy data for a 'transducersOff' message (see AsierInho_V2.cpp::turnTransducersOff)
	memset(messagesArray, 0, numTransducers*2* sizeof(unsigned char));
	for (int b = 0; b < numTransducers; b+=256) {
		messagesArray[0 + 2*b] = 0 + 128;	//Set flag indicating new update (>128)
		messagesArray[25 + 2*b] = 128;	//Set flag indicating to change the LED colour (turn off in this case)
	}
}

void HologramSolution::configureSolutionNewDivider(int numTransducers, unsigned char FPSdivider) {
	transducersOffSolution = true;
	eventsReleased = true;		//We did not allocate any events, so nothing to release when it returns to Solution pool
	//Create custom messages: One message to set the new divider and 15 "turn off"
	numGeometries = 16;
	/////////////////////////////////////////////////////////////////////////////////
	//MESSAGE 0: New dividder:
	//1. Break divider into bits:
	unsigned char dividerBits[] = {
		(FPSdivider / 128) % 2,
		(FPSdivider / 64) % 2,
		(FPSdivider / 32) % 2,
		(FPSdivider / 16) % 2,
		(FPSdivider / 8) % 2,
		(FPSdivider / 4) % 2,
		(FPSdivider / 2) % 2,
		(FPSdivider / 1) % 2,
	};
	//2. Copy data for a 'transducersOff' message (see AsierInho_V2.cpp::turnTransducersOff)
	memset(messagesArray, 0, numTransducers * 2 * sizeof(unsigned char));
	for (int b = 0; b < numTransducers; b += 256) {
		messagesArray[0 + 2 * b] = 0 + 128;	//Set flag indicating new update (>128)
		messagesArray[25 + 2 * b] = 128;	//Set flag indicating to change the LED colour (turn off in this case)
		//Encode new divider: 
		messagesArray[26 + 2 * b] = dividerBits[7] * 128;//lowest bit
		messagesArray[27 + 2 * b] = dividerBits[6] * 128;
		messagesArray[28 + 2 * b] = dividerBits[5] * 128;
		messagesArray[29 + 2 * b] = dividerBits[4] * 128;
		messagesArray[30 + 2 * b] = dividerBits[3] * 128;
		messagesArray[31 + 2 * b] = dividerBits[2] * 128;
		messagesArray[32 + 2 * b] = dividerBits[1] * 128;
		messagesArray[33 + 2 * b] = dividerBits[0] * 128;//highest bit
		messagesArray[34 + 2 * b] =  128;				 //Signal that we are updating the divider.
	}
	//////////////////////////////////////////////////////////////////////////////////
	//MESSAGE 1-15: Turn transducers off:
	for (int msg = 1; msg < 16; msg++) {
		memset(&(messagesArray[numTransducers * 2 * msg]), 0, numTransducers * 2 * sizeof(unsigned char));
		for (int b = 0; b < numTransducers; b += 256) {
			messagesArray[numTransducers * 2 * msg + 2 * b] = 0 + 128;	//Set flag indicating new update (>128)
			messagesArray[numTransducers * 2 * msg + 25 + 2 * b] = 128;	//Set flag indicating to change the LED colour (turn off in this case)
		}
	}

}

float* HologramSolution::finalArrayPhases() {
	cl_int err;	
	if (transducersOffSolution) {
		GSPAT_V2::printWarning_GSPAT("Trying to call 'finalArrayPhases' on a dummy (transducersOFF) solution. Ignored.");
		return NULL;
	}
	err = clEnqueueReadBuffer(queue, this->finalHologramPhases_CLBuffer, CL_TRUE, 0, numTransducers*numGeometries  *  sizeof(float), &(this->transducerPhases[0]), 1, &(events[GSPAT_Event::HOLOGRAM_READY]), NULL);
	if (err < 0) {
		GSPAT_V2::printError_GSPAT("Couldn't read final hologram (phases)\n"); return NULL;
	}

	return &(transducerPhases[0]);
}

float* HologramSolution::finalArrayAmplitudes() {
	if (transducersOffSolution) {
		GSPAT_V2::printWarning_GSPAT("Trying to call 'finalArrayAmplitudes' on a dummy (transducersOFF) solution. Ignored.");
		return NULL;
	}
	cl_int err = clEnqueueReadBuffer(queue, this->finalHologramAmplitudes_CLBuffer, CL_TRUE, 0, numTransducers  * numGeometries * sizeof(float), &(this->transducerAmplitudes[0]), 1, &(events[GSPAT_Event::HOLOGRAM_READY]), NULL);
	if (err < 0) {
		GSPAT_V2::printError_GSPAT("Couldn't read final hologram (amplitudes)"); return NULL;
	}	
		
	return &(transducerAmplitudes[0]);
}

void HologramSolution::finalMessages(unsigned char** out) {
	if (!transducersOffSolution) {
		cl_event messageBuffersRead[] = { events[MESSAGES_READY] };
		cl_int err = clWaitForEvents(1, messageBuffersRead);
		if (err < 0) {
			GSPAT_V2::printError_GSPAT("Couldn't read message queues"); ;
		}
	}
	*out = &(this->messagesArray[0]);	
}

float* HologramSolution::finalHologramReIm() {
	if (transducersOffSolution) {
		GSPAT_V2::printWarning_GSPAT("Trying to call 'finalHologramReIm' on a dummy (transducersOFF) solution. Ignored.");
		return NULL;
	}	
	cl_int err = clEnqueueReadBuffer(queue, this->finalHologramReIm_CLBuffer, CL_TRUE, 0, numTransducers*numGeometries * 2 * sizeof(float), &(this->transducerReIm[0]), 1, &(events[GSPAT_Event::HOLOGRAM_READY]), NULL);
	if (err < 0) {
		GSPAT_V2::printError_GSPAT("Couldn't read final hologram"); return NULL;
	}
	
	return &(transducerReIm[0]);
}

void HologramSolution::releaseEvents() {
	if (eventsReleased)
		return;
	//2. Release all OpenCL events
	//2.1. Events related to input data (depends on whether data was provided manually or externally)
	if(!manualData)
		clReleaseEvent(events[POSITIONS_UPLOADED]);//Data uploaded externally. Notified by a single event.
	else {
		clReleaseEvent(events[POSITIONS_UPLOADED]);
		clReleaseEvent(events[AMPLITUDES_UPLOADED]);
		clReleaseEvent(events[INITIAL_GUESS_UPLOADED]);
		clReleaseEvent(events[MATRIX_0_UPLOADED]);
		clReleaseEvent(events[MATRIX_G_UPLOADED]);
	}
	//2.2. Pipeline events: 
	int lastEvent = R0_READY + geometries();
	for (int e = F_AND_B_READY; e < lastEvent; e++)
		clReleaseEvent(events[e]);
	eventsReleased = true;
}

int HologramSolution::dataBuffers(void* vBuffers) {
	if (vBuffers == NULL)
		return 5;
	cl_mem* pBuffers = (cl_mem*)vBuffers;
	pBuffers[0] = positions_CLBuffer;
	pBuffers[1] = targetAmplitudes_CLBuffer;
	pBuffers[2] = pointsReIm_CLBuffer;
	pBuffers[3] = matrixStart_CLBuffer;
	pBuffers[4] = matrixEnd_CLBuffer;
	return 5;
}

void HologramSolution::dataExternallyUploadedEvent(void* pExternalEvents, int numEvents) {
	if (numEvents != 1) {
		GSPAT_V2::printError_GSPAT("HologramSolution:dataExternallyUploadedEvent-> Only one event expected.");
		return; 
	}
	events[POSITIONS_UPLOADED] = *((cl_event*)pExternalEvents);
}



void HologramSolution::readPropagators(float* singlePointFields, float* singlePointFieldsNormalised, int numGeometries) {
	if (transducersOffSolution) {
		GSPAT_V2::printWarning_GSPAT("Trying to call 'readPropagators' on a dummy (transducersOFF) solution. Ignored.");
		return ;
	}
	clEnqueueReadBuffer(queue, singlePointField_CLBuffer, CL_TRUE,  0, numGeometries*numPoints*numTransducers * 2 * sizeof(float), singlePointFields, 0, NULL, NULL);
	clEnqueueReadBuffer(queue, singlePointFieldNormalised_CLBuffer, CL_TRUE, 0, numGeometries*numPoints*numTransducers * 2 * sizeof(float), singlePointFieldsNormalised, 0, NULL, NULL);
}

void HologramSolution::readMatrixR(float* R, int numGeometries) {
	if (transducersOffSolution) {
		GSPAT_V2::printWarning_GSPAT("Trying to call 'readMatrixR' on a dummy (transducersOFF) solution. Ignored.");
		return ;
	}
	clEnqueueReadBuffer(queue, R_CLBuffer, CL_TRUE, 0, numGeometries*numPoints*numPoints * 2 * sizeof(float), R, 0, NULL, NULL);
}
void HologramSolution::readTargetPointsReIm(float* targetPoints, int numGeometries) {
	if (transducersOffSolution) {
		GSPAT_V2::printWarning_GSPAT("Trying to call 'readTargetPointsReIm' on a dummy (transducersOFF) solution. Ignored.");
		return ;
	}
	clEnqueueReadBuffer(queue, pointsReIm_CLBuffer, CL_TRUE, 0, numGeometries*numPoints* 2 * sizeof(float), targetPoints, 0, NULL, NULL);
}

HologramSolution::~HologramSolution() {
 	releaseEvents();//Just in case the client forgot to do this (Naughty client!)
	//1. Release all OpenCL buffers
	//cl_uint count;
	//clGetMemObjectInfo(positions_CLBuffer, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(positions_CLBuffer);
	//clGetMemObjectInfo(targetAmplitudes_CLBuffer, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(targetAmplitudes_CLBuffer);
	//clGetMemObjectInfo(matrixStart_CLBuffer, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(matrixStart_CLBuffer);
	//clGetMemObjectInfo(matrixEnd_CLBuffer, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(matrixEnd_CLBuffer);
	//clGetMemObjectInfo(singlePointField_CLBuffer, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(singlePointField_CLBuffer);
	//clGetMemObjectInfo(singlePointFieldNormalised_CLBuffer, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(singlePointFieldNormalised_CLBuffer);
	//clGetMemObjectInfo(R_CLBuffer, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(R_CLBuffer);
	//clGetMemObjectInfo(pointsReIm_CLBuffer, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(pointsReIm_CLBuffer);
	//clGetMemObjectInfo(finalHologramPhases_CLBuffer, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(finalHologramPhases_CLBuffer);
	//clGetMemObjectInfo(finalHologramAmplitudes_CLBuffer, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(finalHologramAmplitudes_CLBuffer);
	//clGetMemObjectInfo(finalHologramReIm_CLBuffer, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(finalHologramReIm_CLBuffer);
	//clGetMemObjectInfo(amplitudesPerPoint, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(amplitudesPerPoint);
	//clGetMemObjectInfo(metric, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(metric);
	//clGetMemObjectInfo(correction, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(correction);
	//clGetMemObjectInfo(messagesTopArray_CLBuffer, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	clReleaseMemObject(messagesArray_CLBuffer);
	//clGetMemObjectInfo(messagesBottomArray_CLBuffer, CL_MEM_REFERENCE_COUNT, sizeof(cl_uint), &count, NULL);
	//clReleaseMemObject(messagesBottomArray_CLBuffer);	
	delete transducerReIm;
	delete transducerPhases;
	delete transducerAmplitudes;		
	delete messagesArray;
}