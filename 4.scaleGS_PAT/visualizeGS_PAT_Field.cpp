#include <GSPAT_SolverV2.h>
#include <GSPAT_SolverV3.h>
#include <GSPAT_SolverV4.h>
#include <AsierInho_V2.h>
#include <stdio.h>
#include <conio.h>
#include <../Helper/VisualizePlane.h>
void print(const char* str) {
	printf("%s\n", str);
}



void application(void) {
	const int numBoards = 2;
	const int numTransducersPerBoard = 256;
	const size_t numGeometries = 1;
	const size_t numPoints = 1;
	//Create driver and connect to it
	AsierInho_V2::RegisterPrintFuncs(print, print, print);	
	AsierInho_V2::AsierInhoBoard_V2* driver= AsierInho_V2::createAsierInho();
	int boardIDs[] = { 1, 2,3,4 };
	float matBoardToWorld[64] = {/*bottom-left*/
								1, 0, 0, 0.0f,
								0, 1, 0, 0,
								0, 0, 1, 0,
								0, 0, 0, 1,	
								/*top-left*/
								-1, 0, 0, 0.0f,
								0, 1, 0, 0,
								0, 0,-1, 0.2388f,
								0, 0, 0, 1,	
								
	};
	if(!driver->connect(numBoards,boardIDs, matBoardToWorld))	//Device IDs to connect to
		printf("Failed to connect to board.");
	
	//Create solver:
	GSPAT_V2::RegisterPrintFuncs(print, print, print);
	GSPAT::Solver* solver = GSPAT_V2::createSolver(numBoards*numTransducersPerBoard);//Number of transducers used (two boards of 16x16)
	//Configura the solver
	float transducerPositions[numBoards*numTransducersPerBoard * 3], amplitudeAdjust[numBoards*numTransducersPerBoard];
	int mappings[numBoards*numTransducersPerBoard], phaseDelays[numBoards*numTransducersPerBoard], numDiscreteLevels;
	driver->readParameters(transducerPositions, mappings, phaseDelays, amplitudeAdjust, &numDiscreteLevels);
	solver->setBoardConfig(transducerPositions, mappings, phaseDelays, amplitudeAdjust, numDiscreteLevels);
	//Program: Create a trap and visualize it
	float curPos[4 * numPoints*numGeometries];
	float amplitude[numPoints*numGeometries] ;
	for (int g = 0; g < numGeometries; g++) {
		curPos[4*numPoints* g] = 0.00f-0.0f;
		curPos[4*numPoints* g+1] = 0.00f+0;
		curPos[4*numPoints* g+2] =  0.00f +(0.2388f / 2);
		curPos[4*numPoints* g+3] = 1;
		amplitude[numPoints*g] = 17000;		
		}
	unsigned char* msg;
	float m1[] = { 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1
				  ,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
	//a. create a solution
	bool phaseOnly = false;
	GSPAT::Solution* solution = solver->createSolution(numPoints,numGeometries,phaseOnly, curPos, amplitude,m1,m1);
	solver->compute(solution); 
	//b. Send to levitator (if connected)
	solution->finalMessages(&msg);
	driver->updateMessage(msg);
	//c. Retrieve solution field, make phase-only and apply signature
	float* field=solution->finalHologramReIm();
	//float propagatorsNormalised[2 * 512], propagators[2 * 512];
	//solution->readPropagators(propagators, propagatorsNormalised, 1);
	//float field[512 * 2]; memcpy(field, field_GSPAT, 2 * 512 * sizeof(float));*/
	for(int t=0;t<2*numTransducersPerBoard; t++)
	{
		float magnitude;
		if(phaseOnly)magnitude= sqrtf(field[2 * t] * field[2 * t] + field[2 * t + 1] * field[2 * t + 1]);
		else magnitude = 1;
		field[2*t] = (t>=numTransducersPerBoard?-1:1)*field[2*t]/magnitude;		//Real +signature
		field[2*t+1] = (t>=numTransducersPerBoard?-1:1)*field[2*t+1]/magnitude;	//Imag + signature 			
	}
	//c. Visualize field
	{
		float A[] = {curPos[0] -0.015f, curPos[1]+0,curPos[2] + 0.015f }, B[] = {curPos[0]+ 0.015f,curPos[1]+ 0, curPos[2]+ 0.015f }, C[] = {curPos[0] -0.015f, curPos[1]+0, curPos[2]- 0.015f };
		int imgRes[] = { 64,64 };
		VisualizePlane::visualize(A,B,C, imgRes , field , numBoards*numTransducersPerBoard, transducerPositions);
	}

	char c;
	scanf("%c", &c);
	//d. Destroy solution
	solver->releaseSolution(solution);
		
	//Deallocate the AsierInho controller: 
	driver->turnTransducersOff();
	driver->disconnect();
	delete driver;
	delete solver;
}

void main(){
	while (true) {
		printf("Press any key to start GSPAT.\n Press 'X' for exit\n");
		char option = getch();
		if (option == 'x' || option == 'X')
			return;
		//Run application:
		application();
	}
}