#include <GSPAT_SolverV4.h>
#include <GSPAT_SolverV3.h>
#include <GSPAT_SolverV2.h>
#include <GSPAT_SolverIBP.h>
#include <GSPAT_SolverNaive.h>
#include <AsierInho.h>
#include <stdio.h>
#include <conio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <../Helper/VisualizePlane.h>
#include <Windows.h>
void print(const char* str) {
	printf("%s\n", str);
}

void application(void) {
	const size_t numPoints = 1; //Change this if you want (but make sure to change also each point's position )
	const size_t numGeometries = 1;//Please, do not change this. Multiple update Messages require AsierInhoV2 (see 6.simpleGSPATCpp_AsierInhoV2)
	
								   //Create driver and connect to it
	AsierInho::RegisterPrintFuncs(print, print, print);	
	AsierInho::AsierInhoBoard* driver= AsierInho::createAsierInho();
	//Create solver:
	GSPAT_IBP::RegisterPrintFuncs(print, print, print);
	GSPAT::Solver* solver = GSPAT_IBP::createSolver(512);//Number of transducers used (two boards of 16x16)
	if(!driver->connect(AsierInho::BensDesign, 32, 30))	//Device IDs to connect to
		printf("Failed to connect to board.");
	float transducerPositions[512 * 3], amplitudeAdjust[512];
	int mappings[512], phaseDelays[512], numDiscreteLevels;
	//driver->readAdjustments(mappings, phaseDelays);
	driver->readParameters(transducerPositions, mappings, phaseDelays, amplitudeAdjust, &numDiscreteLevels);
	//DEBUG: Test removing noisy transducers (see "Examples/AsierInho/0.testTransducers").
	/*for (int i = 0; i < 512; i++)
		amplitudeAdjust[i] = 1.0f;*/
	//END DEBUG
	solver->setBoardConfig(transducerPositions, mappings, phaseDelays, amplitudeAdjust, numDiscreteLevels);
	//Program: Create a trap and move it with the keyboard
	static float radius = 0;// 0.02f;
	float curPos[4 * numGeometries*numPoints];
	float amplitude[numGeometries*numPoints] ;
	for (int g = 0; g < numGeometries; g++) {
		for (int p = 0; p < numPoints; p++) {
			float angle = 2 * M_PI / numPoints;
			curPos[4 * g*numPoints + 4*p + 0] = radius * cos(p*angle);
			curPos[4 * g*numPoints + 4*p + 1] = radius * sin(p*angle);
			curPos[4 * g*numPoints + 4*p + 2] = 0.12f;
			curPos[4 * g*numPoints + 4*p + 3] = 1;
			amplitude[g*numPoints + p] = 10000;
		}
	}
	unsigned char* msg;
	float m1[] = { 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1 };
	//a. create a solution
	GSPAT::Solution* solution = solver->createSolution(numPoints,numGeometries,false, curPos, amplitude,m1,m1);
	solver->compute(solution); 
	solution->finalMessages(&msg);
	for (int s = 0; s < 16; s++)//Fill FPGA buffers so update is processed directly
		driver->updateMessage(msg);
	
	/*{
		float A[] = { -0.1f, 0, 0.2388f / 2 + 0.1f }, B[] = { 0.1f, 0, 0.2388f / 2 + 0.1f }, C[] = { -0.1f, 0, 0.2388f / 2 - 0.1f };
		int imgRes[] = { 128,128 };
		VisualizePlane::visualize(A,B,C, imgRes , solution->finalHologramReIm(), 512, transducerPositions);
		VisualizePlane::visualizeFromPhases(A,B,C, imgRes , solution->finalArrayPhases(), 512, transducerPositions);

	}*/
	solver->releaseSolution(solution);

	//b. Main loop (finished when space bar is pressed):
	printf("\n Place a bead at (%f,%f,%f)\n ", curPos[0], curPos[1], curPos[2]);
	printf("Use keys A-D , W-S and Q-E to move the bead\n");
	printf("Press 'X' to destroy the solver.\n");

	bool finished = false;
	while (!finished) {
		//Update 3D position from keyboard presses
		switch (getch()) {
			case 'a':
				for (int g = 0; g < numGeometries; g++)
					for (int p = 0; p < numPoints; p++) 
						curPos[4 * g*numPoints + 4*p + 0] += 0.0005f; break;
			case 'd':
				for (int g = 0; g < numGeometries; g++)
					for (int p = 0; p < numPoints; p++) 
						curPos[4 * g*numPoints + 4*p +0] -= 0.0005f; break;
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
		GSPAT::Solution* solution = solver->createSolution(numPoints,numGeometries,true, curPos, amplitude,m1,m1);
		solver->compute( solution);
		solution->finalMessages(&msg);
		for (int s = 0; s < 16; s++)//Fill FPGA buffers so update is processed directly
			driver->updateMessage(msg);		
		solver->releaseSolution(solution);
		printf("\n Bead location (%f,%f,%f)\n ", curPos[0], curPos[1], curPos[2]);
	
	}
	//Deallocate the AsierInho controller: 
	driver->turnTransducersOff();
	Sleep(100);
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