#pragma once
#include <GSPAT_Solver_Prerequisites.h>
#include <GSPAT_Solver.h>
#include <src/HologramSolution.h>
#include <src/SolutionPool.h>
#include <CL/opencl.h>
#include <../Helper/HelperMethods.h>
#include <pthread.h>
#include <deque>
#include <vector>

using namespace GSPAT;

class HologramSolverCL: public Solver {
	/**
		SolutionPool needs to destroy solutions (using Solver's protected method destroySolution), so it is made friend.
	*/
	friend class SolutionPool;
public:
	HologramSolverCL(int numTransducers);
	~HologramSolverCL();
	virtual void* getSolverContext();
	virtual void setBoardConfig( float* transducerPositions, float* transducerNormals, int* transducerToPINMap, int* phaseAdjust, float* amplitudeAdjust, int numDiscreteLevels );
	//virtual void setBoardConfig(char t_Ids[],	char t_phaseAdjust[], int numDiscreteLevels );
	virtual Solution* createSolution(int numPoints, int numGeometries, bool phaseOnly, float* positions, float*amplitudes, float* matStarts, float* matEnds , GSPAT::MatrixAlignment a);
	virtual Solution* createSolutionExternalData(int numPoints, int numGeometries, bool phaseOnly, GSPAT::MatrixAlignment a);
	virtual Solution* createTransducersOffSolution();
	virtual Solution* createNewDividerSolution(unsigned char dividerFPS);
	virtual void compute(Solution* solution);
	//virtual void computeWithExternalPropagators(float* pointHologram, float* normalizedPointHologram, Solution* solution);
	virtual void releaseSolution(Solution* solution);

protected:
	void initializeOpenCL();
	void createDirectivityTexture();
	void createTransducerPositionBuffer(float* transducerPositions);
	void createTransducerNormalBuffer(float* transducerNormals);
	void createSolutionPool();
	/**
		This method computes the directivity of our transducer as a function of the cosine 
		of the angle between the transducer's normal, and the point.
	*/
	float  computeDirectivityForCosAlpha(float cos_alpha, float P_ref=8.02f) {
		return (float)(2*P_ref*_j1(K()*r0()*sqrtf(1-cos_alpha*cos_alpha)) / (K()*r0()*sqrtf(1-cos_alpha*cos_alpha)));
	}

	//Pipeline methods (called every time a solutiuon is computed) :
	void updateCLBuffers(HologramSolution* solution);
	void computeFandB(HologramSolution* solution);
	void computeR(HologramSolution* solution,bool useDiego=false);
	void solvePhases_GS(HologramSolution* solution,bool useNorm);
	void computeActivation(HologramSolution* solution);
	void discretise(HologramSolution* solution);
private: 
	int boardSize[2];
	int numTransducers;
	int numLoopsInKernel; // Added by Ryuji... This is larger than 1 when numTransducers > maxGroupSize of your GPU...
	int numDiscreteLevels;
	float initialGuess[2*HologramSolution::MAX_POINTS*HologramSolution::MAX_GEOMETRIES];	
	bool configured;
	//OpenCL variables:
	const size_t platform_id = 1;
	const size_t device_id = 0;
	cl_platform_id* platforms, platformUsed;
	cl_device_id *devices, deviceUsed;
	cl_context context;
	cl_command_queue queue;
	cl_program program;
	cl_kernel fillAllPointHologramsKernel, powerMethodKernel,  addHologramsKernel, discretiseKernel;
	cl_kernel fillAllPointHologramsLoopKernel;
	//Global buffers (shared among all solutions):
	/** transducerMappings: Describes the pin ID of each transducer in its respective board. 
		Such data is used during the discretization stage to build AsierInho messages.*/
	cl_mem transducerMappings;
	/** phaseCorrections: Describes the phase correction required by each transducer. 
		Such data is used during the discretization stage to build AsierInho messages.*/
	cl_mem phaseCorrections;				//Shared by all solutions
	/** amplitudeCorrections: Describes the phase correction required by each transducer. 
		Such data is used during the discretization stage to build AsierInho messages.*/
	cl_mem amplitudeCorrections;				//Shared by all solutions

	/**
		directivityTexture: Our solver uses a piston model to define the directivity of our transducer sources. Such model describes 
		the scalar response of the transducer at a point in space (in Pascals), as a function of the angle between 
		the transducer’s normal and the point. Computing this angle in the OpenCL shader is (slightly) costly, and 
		it is more convenient to compute the cosine of such angle (i.e. dot product). Thus, the directivity 
		function is precomputed as a function of the cosine, and encoded into a 1D texture (directivityTexture), 
		allowing the shader to access it as if it was a continuous variable.   The texture is created during object 
		construction (see createDirectivityTexture) and used in stage computeFandB. 
	*/
	cl_mem directivityTexture;//rel								//Texture containing the directivity of each transducer, according to cos(alpha)
	
	/**
		transducerPositions: Buffer describing the position of each transducer in homogeneous coordinates (4 floats per position). 
		Positions are automatically computed during object creation in createTransducerPositionBuffer.
		The current implementation is hardcoded for top-bottom arrangements of 512 transducers.
	*/
	cl_mem transducerPositions;	//rel							
	
	/**
		Added by Ryuji
		This allows you to arrange transducers in an arbitrary direction.
	*/
	cl_mem transducerNormals;

	/**
		solutionPool: Contains the pool of solutions to be reused by solver’s clients. 
	*/
	SolutionPool* solutionPool;
};
