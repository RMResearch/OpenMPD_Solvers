#include <GSPATV2_CWrapper.h>
#include <AsierInho_V2_CWrapper.h>
#include <stdio.h>
#include <conio.h>

void print(const char* str) {
	printf("%s\n", str);
}

void application(void) {
	//Create handler and connect to it
	AsierInho_V2_Handler handler = AsierInho_CWrapper_createHandler();
	GSPATV2_CWrapper_Initialize();
	AsierInho_CWrapper_RegisterPrintFuncs(print, print, print);
	GSPATV2_CWrapper_RegisterPrintFuncs(print, print, print);
	GSPAT_Solver_Handler solver = GSPATV2_CWrapper_createSolver(512);
	while(!AsierInho_CWrapper_connectTopBottom(handler, 2, 4))
		printf("Failed to connect to board.");
	float transducerPositions[512 * 3], amplitudeAdjust[512];
	int mappings[512], phaseDelays[512], numDiscreteLevels;
	//AsierInho_CWrapper_readAdjustments(handler, mappings, phaseDelays);
	AsierInho_CWrapper_readParameters(handler, transducerPositions, mappings, phaseDelays, amplitudeAdjust, &numDiscreteLevels);
	GSPATV2_CWrapper_setBoardConfig(solver, transducerPositions, mappings, phaseDelays, amplitudeAdjust, numDiscreteLevels);
	//Program: Create a trap and move it with the keyboard
	float curPos[] = { 0,0,0.1f , 1}, amplitude =1;
	unsigned char* msg;
	float m1[] = { 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1 };
	//a. create a solution
	GSPAT_Solution_Handler solution = GSPATV2_CWrapper_createSolution(solver, 1, 1,true, curPos, &amplitude,m1,m1);
	GSPATV2_CWrapper_compute(solver, solution); 
	GSPATV2_CWrapper_finalMessages(solution, &msg);
	AsierInho_CWrapper_updateMessage(handler, msg);
	GSPATV2_CWrapper_releaseSolution(solver, solution);
	//char c;  scanf("%c", &c);//Wait so that we can put the bead in place

	//b. Main loop (finished when space bar is pressed):
	printf("\n Place a bead at (%f,%f,%f)\n ", curPos[0], curPos[1], curPos[2]);
	printf("Use keys A-D , W-S and Q-E to move the bead\n");
	printf("Press 'X' to destroy the solver.\n");

	bool finished = false;
	while (!finished) {
		//Update 3D position from keyboard presses
		switch (getch()) {
			case 'a':curPos[0] += 0.0005f; break;
			case 'd':curPos[0] -= 0.0005f; break;
			case 'w':curPos[1] += 0.0005f; break;
			case 's':curPos[1] -= 0.0005f; break;
			case 'q':curPos[2] += 0.00025f; break;
			case 'e':curPos[2] -= 0.00025f; break;
			case 'x': 
			case 'X':	
				finished = true;

		}
		//Create the trap and send to the board:
		GSPAT_Solution_Handler solution = GSPATV2_CWrapper_createSolution(solver, 1, 1,true, curPos, &amplitude,m1,m1);
		GSPATV2_CWrapper_compute(solver, solution);
		GSPATV2_CWrapper_finalMessages(solution, &msg);
		AsierInho_CWrapper_updateMessage(handler, msg);
		GSPATV2_CWrapper_releaseSolution(solver, solution);
	}
	//Deallocate the AsierInho controller: 
	AsierInho_CWrapper_turnTransducersOff(handler);
	AsierInho_CWrapper_disconnect(handler);
	AsierInho_CWrapper_destroyHandler(handler);
	AsierInho_CWrapper_Release();
	GSPATV2_CWrapper_destroySolver(solver);
	GSPATV2_CWrapper_Release();
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