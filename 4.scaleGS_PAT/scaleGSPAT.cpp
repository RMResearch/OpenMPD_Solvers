#include <GSPAT_SolverV2.h>
#include <AsierInho_V2.h>
#include <stdio.h>
#include <conio.h>
#include <../Helper/VisualizePlane.h>
void print(const char* str) {
	printf("%s\n", str);
}



void application(void) {
	const int numBoards = 4;
	const int numTransducersPerBoard = 256;
	const size_t numGeometries = 32;
	const size_t numPoints = 1;
	//Create driver and connect to it
	AsierInho_V2::RegisterPrintFuncs(print, print, print);	
	AsierInho_V2::AsierInhoBoard_V2* driver= AsierInho_V2::createAsierInho();
	int boardIDs[] = { 1, 2,3,4,5,6,7,8 };
	float matBoardToWorld[128] = {/*bottom-left*/
								1, 0, 0, -0.08f,
								0, 1, 0, 0,
								0, 0, 1, 0,
								0, 0, 0, 1,	
								/*top-left*/
								-1, 0, 0, -0.08f,
								0, 1, 0, 0,
								0, 0,-1, 0.2388f,
								0, 0, 0, 1,	
								/*bottom-right*/
								1, 0, 0, 0.08f,
								0, 1, 0, 0,
								0, 0, 1, 0,
								0, 0, 0, 1,	
								/*top-right*/
								-1, 0, 0, 0.08f,
								0, 1, 0, 0,
								0, 0,-1, 0.2388f,
								0, 0, 0, 1,	
								/*bottom-left*/
								1, 0, 0, -0.16f,
								0, 1, 0, 0,
								0, 0, 1, 0,
								0, 0, 0, 1,	
								/*top-left*/
								-1, 0, 0, -0.16f,
								0, 1, 0, 0,
								0, 0,-1, 0.2388f,
								0, 0, 0, 1,	
								/*bottom-right*/
								1, 0, 0, 0.16f,
								0, 1, 0, 0,
								0, 0, 1, 0,
								0, 0, 0, 1,	
								/*top-right*/
								-1, 0, 0, 0.16f,
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
		curPos[4*numPoints* g] = -0.0f;
		curPos[4*numPoints* g+1] = 0;
		curPos[4*numPoints* g+2] = 0.2388f / 2;
		curPos[4*numPoints* g+3] = 1;
		amplitude[numPoints*g] = 10000;
		/*curPos[4*numPoints* g+4] = -0.05f;
		curPos[4*numPoints* g+5] = 0;
		curPos[4*numPoints* g+6] = 0.15f;
		curPos[4*numPoints* g+7] = 1;
		amplitude[2*g+1] = 10000;*/
		}
	unsigned char* msg;
	float m1[] = { 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1
				  ,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
	//a. create a solution
	GSPAT::Solution* solution = solver->createSolution(numPoints,numGeometries,false, curPos, amplitude,m1,m1);
	solver->compute(solution); 
	//b. Send to levitator (if connected)
	solution->finalMessages(&msg);
	driver->updateMessage(msg);
	//c. Visualize field
	{
		float A[] = { -0.24f, 0, 0.2388f / 2 + 0.12f }, B[] = { 0.24f, 0, 0.2388f / 2 + 0.12f }, C[] = { -0.24f, 0, 0.2388f / 2 - 0.12f };
		int imgRes[] = { 2*128,128 };
		VisualizePlane::visualize(A,B,C, imgRes , solution->finalHologramReIm(), numBoards*numTransducersPerBoard, transducerPositions);
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