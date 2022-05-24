#ifndef _GSPAT_SOLUTIONIMPL
#define _GSPAT_SOLUTIONIMPL
#include <GSPAT_Solution.h>
#include <CL/opencl.h>

using namespace GSPAT;

enum GSPAT_Event {
	POSITIONS_UPLOADED=0, AMPLITUDES_UPLOADED,INITIAL_GUESS_UPLOADED, MATRIX_0_UPLOADED, MATRIX_G_UPLOADED, 
	F_AND_B_READY, POINT_PHASES_READY, HOLOGRAM_READY, HOLOGRAM_DISCRETISED, MESSAGES_READY,
	R0_READY, R1_READY, R2_READY, R3_READY, R4_READY, R5_READY, R6_READY, R7_READY,R8_READY, R9_READY,
	R10_READY, R11_READY, R12_READY, R13_READY, R14_READY, R15_READY, R16_READY, R17_READY,R18_READY, R19_READY,
	R20_READY, R21_READY, R22_READY, R23_READY, R24_READY, R25_READY, R26_READY, R27_READY,R28_READY, R29_READY,
	R30_READY, R31_READY, R32_READY, R33_READY, R34_READY, R35_READY, R36_READY, R37_READY,R38_READY, R39_READY,
	EVENTS_TOTAL
};

class HologramSolution : public Solution {
		friend class HologramSolverCL;
	protected:
		static const int MAX_POINTS=32;
		static const int MAX_GEOMETRIES = 32;
		//CPU Storage for parameters.
		bool manualData; //Indicates if data has been manually provided by client or indirectly uploaded to Solution's buffers.
		bool transducersOffSolution; //Indicatesif this is a "dummy" solution to turn levitator off
		int numPoints;
		int numGeometries;
		bool phaseOnly;
		int numTransducers;
		float positions[MAX_POINTS*MAX_GEOMETRIES*4];//del
		float amplitudes[MAX_POINTS*MAX_GEOMETRIES];//del
		float* transducerReIm;										// Final transducers' state (ReIm, just focussing).
		float* transducerPhases;									//Final transducer phases to send to the device.
		float* transducerAmplitudes;								//Final transducer amplitudes to send to the device.
		unsigned char* messagesArray;								//Message to send to the Asierinho board (top-bottom)
		float matrixStart[16* MAX_POINTS], matrixEnd[16* MAX_POINTS];				//Transformation to apply to first and last geometries (intermediate ones will interpolate -> WARNING! Only OK if transformations are small, which they must be in our case)
		
		//GPU Storage and events required to coordinate execution along the pipeline. 
		cl_command_queue queue;

		//Events
		bool eventsReleased;
		cl_event events[GSPAT_Event::EVENTS_TOTAL];
		//Buffers
		cl_mem positions_CLBuffer;//rel
		cl_mem targetAmplitudes_CLBuffer;//rel			//CL Buffer containing the desired amplitudes per point					
		cl_mem matrixStart_CLBuffer, matrixEnd_CLBuffer;//rel
		cl_mem singlePointField_CLBuffer;//rel			//CL buffer storing single point holograms.
		cl_mem singlePointFieldNormalised_CLBuffer;//rel
		cl_mem R_CLBuffer;//rel							//CL buffer storing propagation matrix R
		cl_mem pointsReIm_CLBuffer;//rel					//CL buffer storing the state of the points.
		cl_mem finalHologramPhases_CLBuffer;//rel		//CL buffer storing the final phases to send to the array (with lev signature, phase only, A=1)
		cl_mem finalHologramAmplitudes_CLBuffer;//rel	//CL buffer storing the final phases to send to the array (with lev signature, phase only, A=1)
		cl_mem finalHologramReIm_CLBuffer;//rel			//CL Buffer storing the "focussing hologram" (Re and Im parts, no lev signature)
		cl_mem amplitudesPerPoint;//rel					//CL Buffer containing estimated amplitudes per point (from matrix R)
		cl_mem metric;//rel								//CL_Buffer containing AVG(A) and STEDEV(A) of the solution.
		cl_mem correction;	//rel						//CL_Buffer containing the overall amplitude corection to apply to the solution to minimize reconstruction error.
		cl_mem messagesArray_CLBuffer;//rel

	public:
		/**
			Returns the number of geometries (i.e. arrangements of points) to be computed in parallel by the solver.
		*/
		virtual int geometries() {
			return numGeometries;
		}
		/**
			Returns the number of points in each geometry to be computed by the solver (all geometries must have the same number of points).
		*/
		virtual int points() {
			return numPoints;
		}
		/**
			Indicates that the solution is being processed by the pipeline. Solution-related events will be automatically filled 
			by the solvre, and they must be released (only) at the end of computation.
		*/
		virtual void lockEvents(){
			eventsReleased = false;
		}
		//Life cycle methods
		HologramSolution(int numTransducers, cl_context context, cl_command_queue& queue);
		void configureSolution(int numPoints, int numGeometries, bool phaseOnly, float* positions, float* amplitudes, float* mStart, float* mEnd);
		void configureSolutionExternalData(int numPoints, int numGeometries, bool phaseOnly);
		void configureSolutionTransducersOff(int numTransducers);
		void configureSolutionNewDivider(int numTransducers, unsigned char  FPSdivider);
		virtual void releaseEvents();
		//Life cycle methods for data computed externally (e.g. Rendering Engine)
		/** Returns the 5 OpenCL data buffers (cl_mem) required to store that data related to the current solution.
			These are stored in the following order: 
				[0] positions_CLBuffer
				[1] amplitudes_CLBuffer
				[2] pointsReIm_CLBuffer
				[3] matrixStart_CLBuffer
				[4] matrixEnd_CLBuffer		*/
		virtual int dataBuffers(void* pBuffers) ;
		/** The external client must provide a single event notifying that all data is ready to be used.*/
		virtual void dataExternallyUploadedEvent(void* pExternalEvents, int numEvents);

		//Query methods (post computation)
		virtual float* finalArrayPhases();
		virtual float* finalArrayAmplitudes();
		virtual float* finalHologramReIm();
		virtual void finalMessages(unsigned char** messages);
		//Methods to read results from intermediate steps (mostly used for debugging)
		virtual void readPropagators(float* singlePointFields, float* singlePointFieldsNormalised, int numGeometries );
		virtual void readMatrixR(float* R, int numGeometries);
		virtual void readTargetPointsReIm(float* targetPoints, int numGeometries);
protected:
	/**
		Solution destructors must remain protected to avoid clients from (accidentally) destroying them.
		Deleting solutions is still available using GSPAT::Solver::destroySolution, but this method is only 
		available to Solvers (i.e. DLL clients cannot use it).	*/
	virtual ~HologramSolution();

	};

#endif
