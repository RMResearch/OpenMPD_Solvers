#include <GSPAT_SolverIBP.h>
#include <src/HologramSolverCL.h>

GSPAT::Solver* GSPAT_IBP::createSolver(int numTransducers) {
	return new HologramSolverCL(numTransducers);
}

void (*_printMessage_GSPAT) (const char*) = NULL;
void (*_printError_GSPAT)(const char*)= NULL;
void (*_printWarning_GSPAT)(const char*)= NULL;

void GSPAT_IBP::printMessage_GSPAT(const char* str) {
	if (_printMessage_GSPAT) _printMessage_GSPAT(str);
}
void GSPAT_IBP::printError_GSPAT(const char* str) {
	if (_printError_GSPAT) _printError_GSPAT(str);
}
void GSPAT_IBP::printWarning_GSPAT(const char* str) {
	if (_printWarning_GSPAT) _printWarning_GSPAT(str);
}

void GSPAT_IBP::RegisterPrintFuncs(void(*p_Message)(const char*), void(*p_Warning)(const char*), void(*p_Error)(const char*)) {
	if(!_printMessage_GSPAT )_printMessage_GSPAT = p_Message;
	if(!_printError_GSPAT )_printError_GSPAT = p_Error;
	if(!_printWarning_GSPAT )_printWarning_GSPAT = p_Warning;
	printMessage_GSPAT("Got message functions!");
}

