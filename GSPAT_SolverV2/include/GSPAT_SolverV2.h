#pragma once
#include <GSPAT_Solver.h>

namespace GSPAT_V2 {
	/**
		Create a Solver. Width and Height specify the number of transducers, but currently only 32x16 transducers is supported.

	*/
	_GSPAT_Export GSPAT::Solver* createSolver(int numTransducers);

	void printMessage_GSPAT(const char*);
	void printError_GSPAT(const char*);
	void printWarning_GSPAT(const char*);

	_GSPAT_Export void RegisterPrintFuncs(void(*p_Message)(const char*), void(*p_Warning)(const char*), void(*p_Error)(const char*));
};