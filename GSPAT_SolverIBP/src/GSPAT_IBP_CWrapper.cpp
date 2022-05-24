#include <GSPAT_IBP_CWrapper.h>
#include <GSPAT_SolverIBP.h>
#include <string>       // std::string
#include <iostream>     // std::cout, std::ostream, std::hex
#include <sstream>
#include <list>

using namespace GSPAT_IBP;
static std::list<GSPAT::Solver*> handlers;

bool GSPAT_IBP_CWrapper_Initialize() {
	return true;
}

bool GSPAT_IBP_CWrapper_Release() {
	//1. Delete all handlers in our lists
	std::list<GSPAT::Solver*>::iterator itHandlers;
	for (itHandlers = handlers.begin(); itHandlers != handlers.end(); itHandlers++) {
		delete(*itHandlers);
		printMessage_GSPAT("Solver Destroyed\n");
	}
	//2.Reset the list to empty
	handlers.clear();	
	return true;
}

void GSPAT_IBP_CWrapper_RegisterPrintFuncs(void(*p_Message)(const char*), void(*p_Warning)(const char*), void(*p_Error)(const char*)) {
	RegisterPrintFuncs(p_Message, p_Warning, p_Error);
	//printMessage_GSPAT("Got print functions!");
}


GSPAT_Solver_Handler GSPAT_IBP_CWrapper_createSolver(int numTransducers) {
	std::stringbuf buffer;             // empty buffer
	std::ostream os (&buffer);      // associate stream buffer to strea
	os << "GSPAT_CWrapper_createSolver(" << numTransducers << ");\n";
	printMessage_GSPAT(buffer.str().c_str());
	GSPAT::Solver* result = GSPAT_IBP::createSolver(numTransducers);
	handlers.push_back(result);
	return (GSPAT_Solver_Handler)result;
}
void GSPAT_IBP_CWrapper_destroySolver(GSPAT_Solver_Handler solver) {
	GSPAT::Solver* target = (GSPAT::Solver*)solver;
	handlers.remove(target);
	printMessage_GSPAT("Solver Destroyed\n");
	delete target;		
}
//Wrapper Methods for GSPAT_Solver
void GSPAT_IBP_CWrapper_setBoardConfig(GSPAT_Solver_Handler solver, float* transducerPositions, int* transducerToPINMap, int* phaseAdjust, float* amplitudeAdjust, int numDiscreteLevels) {
	((GSPAT::Solver*)solver)->setBoardConfig(transducerPositions,transducerToPINMap, phaseAdjust, amplitudeAdjust, numDiscreteLevels);
}
GSPAT_Solution_Handler GSPAT_IBP_CWrapper_createSolution(GSPAT_Solver_Handler solver, int numPoints, int numGeometries, bool phaseOnly, float* positions, float*amplitudes, float* matStartPerPoint, float* matEndPerPoint, int matrixAlignment) {
	//2. Call method:	
	return (GSPAT_Solution_Handler)((GSPAT::Solver*)solver)->createSolution(numPoints, numGeometries, phaseOnly, positions, amplitudes, matStartPerPoint, matEndPerPoint, (GSPAT::MatrixAlignment)matrixAlignment);
}
void GSPAT_IBP_CWrapper_compute(GSPAT_Solver_Handler solver, GSPAT_Solution_Handler solution) {
	((GSPAT::Solver*)solver)->compute((GSPAT::Solution*)solution);
}
void GSPAT_IBP_CWrapper_releaseSolution(GSPAT_Solver_Handler solver, GSPAT_Solution_Handler solution) {
	((GSPAT::Solver*)solver)->releaseSolution((GSPAT::Solution*)solution);
}
/*void GSPAT_IBP_CWrapper_computeWithExternalPropagators(GSPAT_Solver_Handler solver, float* pointHologram, float* normalizedPointHologram, GSPAT_Solution_Handler solution) {
	((Solver*)solver)->computeWithExternalPropagators(pointHologram, normalizedPointHologram, (Solution*)solution);
}*/
//Wrapper Methods for GSPAT_Solution
int GSPAT_IBP_CWrapper_getNumGeometries(GSPAT_Solution_Handler s) {
	return ((GSPAT::Solution*)s)->geometries();
}

int GSPAT_IBP_CWrapper_getNumPoints(GSPAT_Solution_Handler s) {
	return ((GSPAT::Solution*)s)->points();
}

void GSPAT_IBP_CWrapper_finalMessages(GSPAT_Solution_Handler s, unsigned char** out) {
	((GSPAT::Solution*)s)->finalMessages(out);
}
float* GSPAT_IBP_CWrapper_finalArrayPhases(GSPAT_Solution_Handler s) {
	return ((GSPAT::Solution*)s)->finalArrayPhases();
}
float* GSPAT_IBP_CWrapper_finalArrayAmplitudes(GSPAT_Solution_Handler s) {
	return ((GSPAT::Solution*)s)->finalArrayAmplitudes();
}
float* GSPAT_IBP_CWrapper_finalHologramReIm(GSPAT_Solution_Handler s) {
	return ((GSPAT::Solution*)s)->finalHologramReIm();
}
void GSPAT_IBP_CWrapper_readPropagators(GSPAT_Solution_Handler s, float* singlePointFields, float* singlePointFieldsNormalised, int numGeometries) {
	((GSPAT::Solution*)s)->readPropagators(singlePointFields, singlePointFieldsNormalised, numGeometries);
}
void GSPAT_IBP_CWrapper_readMatrixR(GSPAT_Solution_Handler s, float* R, int numGeometries ) {
	((GSPAT::Solution*)s)->readMatrixR(R, numGeometries);
}
void GSPAT_IBP_CWrapper_readTargetPointsReIm(GSPAT_Solution_Handler s, float* targetPoints, int numGeometries) {
	((GSPAT::Solution*)s)->readTargetPointsReIm(targetPoints, numGeometries);
}