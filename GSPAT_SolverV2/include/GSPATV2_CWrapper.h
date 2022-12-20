#ifndef _GSPATV2_CWRAPPER
#define _GSPATV2_CWRAPPER
#include <GSPAT_Solver_Prerequisites.h>

extern "C" {
	static const int GSPAT_ColumnMajorAlignment = GSPAT::ColumnMajorAlignment;
	static const int GSPAT_RowMajorAlignment = GSPAT::RowMajorAlignment;

		//General methods for DLL wrapper:
		_GSPAT_Export bool GSPATV2_CWrapper_Initialize();
		_GSPAT_Export bool GSPATV2_CWrapper_Release();
		_GSPAT_Export void GSPATV2_CWrapper_RegisterPrintFuncs(void(*p_Message)(const char*), void(*p_Warning)(const char*), void(*p_Error)(const char*));

		//Wrapper Methods for GSPAT_Solver
		_GSPAT_Export GSPAT_Solver_Handler GSPATV2_CWrapper_createSolver(int numTransducers);
		_GSPAT_Export void GSPATV2_CWrapper_destroySolver(GSPAT_Solver_Handler);
		_GSPAT_Export void GSPATV2_CWrapper_setBoardConfig(GSPAT_Solver_Handler solver, float* transducerPositions, float* transducerNormals, int* transducerToPINMap, int* phaseAdjust, float* amplitudeAdjust, int numDiscreteLevels = 128) ;
		_GSPAT_Export GSPAT_Solution_Handler GSPATV2_CWrapper_createSolution(GSPAT_Solver_Handler solver, int numPoints, int numGeometries, bool phaseOnly, float* positions, float*amplitudes, float* matStartPerPoint, float* matEndPerPoint, int matrixAlignment=GSPAT_ColumnMajorAlignment) ;
		_GSPAT_Export void GSPATV2_CWrapper_compute(GSPAT_Solver_Handler solver, GSPAT_Solution_Handler solution) ;
		_GSPAT_Export void GSPATV2_CWrapper_releaseSolution(GSPAT_Solver_Handler solver, GSPAT_Solution_Handler solution) ;
		//_GSPAT_Export void GSPAT_CWrapper_computeWithExternalPropagators(GSPAT_Solver_Handler solver ,float* pointHologram, float* normalizedPointHologram, GSPAT_Solution_Handler solution) ;
		
		//Wrapper Methods for GSPAT_Solution
		_GSPAT_Export int GSPATV2_CWrapper_getNumGeometries(GSPAT_Solution_Handler s);
		_GSPAT_Export int GSPATV2_CWrapper_getNumPoints(GSPAT_Solution_Handler s);
		_GSPAT_Export void GSPATV2_CWrapper_finalMessages(GSPAT_Solution_Handler s, unsigned char** out);
		_GSPAT_Export float* GSPATV2_CWrapper_finalArrayPhases(GSPAT_Solution_Handler s);
		_GSPAT_Export float* GSPATV2_CWrapper_finalArrayAmplitudes(GSPAT_Solution_Handler s);
		_GSPAT_Export float* GSPATV2_CWrapper_finalHologramReIm(GSPAT_Solution_Handler s);
		_GSPAT_Export void GSPATV2_CWrapper_readPropagators(GSPAT_Solution_Handler s, float* singlePointFields, float* singlePointFieldsNormalised, int numGeometries = 1);
		_GSPAT_Export void GSPATV2_CWrapper_readMatrixR(GSPAT_Solution_Handler s, float* R, int numGeometries = 1);
		_GSPAT_Export void GSPATV2_CWrapper_readTargetPointsReIm(GSPAT_Solution_Handler s, float* targetPoints, int numGeometries = 1);
};
#endif