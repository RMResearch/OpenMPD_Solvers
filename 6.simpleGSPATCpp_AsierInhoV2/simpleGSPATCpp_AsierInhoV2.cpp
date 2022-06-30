#include <GSPAT_SolverV2.h>
#include <GSPAT_SolverIBP.h>
#include <GSPAT_SolverNaive.h>
#include <AsierInho_V2.h>
#include <stdio.h>
#include <conio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <../Helper/VisualizePlane.h>
void print(const char* str) {
	printf("%s\n", str);
}

void application(void) {
	const size_t numPoints = 1; //Change this if you want (but make sure to change also each point's position )
	const size_t numGeometries = 32;
	//Create driver and connect to it
	AsierInho_V2::RegisterPrintFuncs(print, print, print);	
	AsierInho_V2::AsierInhoBoard_V2* driver= AsierInho_V2::createAsierInho();
	//Create solver:
	GSPAT_V2::RegisterPrintFuncs(print, print, print);
	GSPAT::Solver* solver = GSPAT_V2::createSolver(512);//Number of transducers used (two boards of 16x16)
	//Connect ot the driver (configured for top-bottom manually) - See 7.simpleGSPATC_AsierInhoV2 to see how to do this automatically.
	int boardIDs[] = {  61 , 60 };
	float matBoardToWorld[32] = {/*bottom*/
								1, 0, 0, 0,
								0, 1, 0, 0,
								0, 0, 1, 0,
								0, 0, 0, 1,	
								/*top*/
								-1, 0, 0, 0,
								0, 1, 0, 0,
								0, 0,-1, 0.2388f,
								0, 0, 0, 1,	
	};
	
	if(!driver->connect(2, boardIDs, matBoardToWorld))	//Device IDs to connect to
		printf("Failed to connect to board.");
	float transducerPositions[512 * 3], amplitudeAdjust[512];
	int mappings[512], phaseDelays[512], numDiscreteLevels;
	//driver->readAdjustments(mappings, phaseDelays);
	driver->readParameters(transducerPositions, mappings, phaseDelays, amplitudeAdjust, &numDiscreteLevels);
	solver->setBoardConfig(transducerPositions, mappings, phaseDelays, amplitudeAdjust, numDiscreteLevels);
	//Program: Create a trap and move it with the keyboard
	float curPos[4 * numGeometries*numPoints];
	float amplitude[numGeometries*numPoints] ;
	for (int g = 0; g < numGeometries; g++) {
		for (int p = 0; p < numPoints; p++) {
			float angle = 2 * M_PI / numPoints;
			curPos[4 * g*numPoints + 4*p + 0] = 0.045f * cos(p*angle);
			curPos[4 * g*numPoints + 4*p + 1] = 0.045f * sin(p*angle);
			curPos[4 * g*numPoints + 4*p + 2] = 0.1f;
			curPos[4 * g*numPoints + 4*p + 3] = 1;
			amplitude[g*numPoints + p] = 1;
		}
	}
	unsigned char* msg;
	float mat[] = { 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1 }, m[16*numPoints];
	for (int p = 0; p < numPoints; p++)
		memcpy(&(m[16 * p]), mat, 16 * sizeof(float));
	//a. create a solution
	GSPAT::Solution* solution = solver->createSolution(numPoints,numGeometries,true, curPos, amplitude,m,m);
	solver->compute(solution); 
	solution->finalMessages(&msg);
	driver->updateMessages(msg,numGeometries);	
	solver->releaseSolution(solution);

	//b. Main loop (finished when space bar is pressed):
	printf("\n Place a bead at (%f,%f,%f)\n ", curPos[0], curPos[1], curPos[2]);
	printf("Use keys A-D , W-S and Q-E to move the bead\n");
	printf("Press 'X' to destroy the solver.\n");

	bool finished = false;
	while (!finished) {
		printf("\n Points at (%f,%f,%f)\n ", curPos[0], curPos[1], curPos[2]);
		//Update 3D position from keyboard presses
		switch (getch()) {
			case 'a':
				for (int g = 0; g < numGeometries; g++)
					for (int p = 0; p < numPoints; p++) 
						curPos[4 * g*numPoints + 4*p + 0] += 0.0005f; break;
			case 'd':
				for (int g = 0; g < numGeometries; g++)
					for (int p = 0; p < numPoints; p++) 
						curPos[4 * g*numPoints + 4*p + 0] -= 0.0005f; break;
			case 'w':
				for (int g = 0; g < numGeometries; g++)
					for (int p = 0; p < numPoints; p++) 
						curPos[4 * g*numPoints + 4*p +1] += 0.0005f; break;
			case 's':
				for (int g = 0; g < numGeometries; g++)
					for (int p = 0; p < numPoints; p++) 
						curPos[4 * g*numPoints + 4*p +1] -= 0.0005f; break;
			case 'q':
				for (int g = 0; g < numGeometries; g++)
					for (int p = 0; p < numPoints; p++) 
						curPos[4 * g*numPoints + 4*p +2] += 0.00025f; break;
			case 'e':
				for (int g = 0; g < numGeometries; g++)
					for (int p = 0; p < numPoints; p++) 
						curPos[4 * g*numPoints + 4*p + 2] -= 0.00025f; break;
			case 'x': 
			case 'X':	
				finished = true;

		}
		//Create the trap and send to the board:
		GSPAT::Solution* solution = solver->createSolution(numPoints,numGeometries,true, curPos, amplitude,m,m);
		solver->compute( solution);
		solution->finalMessages(&msg);
		driver->updateMessages(msg,numGeometries);		
		solver->releaseSolution(solution);
	}
	//Deallocate the AsierInho controller: 
	/*solution=solver->createTransducersOffSolution();
	solution->finalMessages(&msg);
	driver->updateMessages(msg,numGeometries);		
	solver->releaseSolution(solution);*/
	//driver->turnTransducersOff();
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