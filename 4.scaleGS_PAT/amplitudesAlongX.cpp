#include <GSPAT_SolverV2.h>
#include <AsierInho_V2.h>
#include <../Helper/GorKovComputation.h>
#include <../Helper/VisualizePlane.h>
#include "CImg\CImg.h"

void print(const char* str) {
	printf("%s\n", str);
}

using namespace cimg_library;


void generateField(float position[3], GSPAT::Solver* solver, float* transducerPositions, lapack_complex_float transducerState[512]) {
	float amplitudes[] = { 17000 };
	bool phaseOnly = true;
	const int numTransducersPerBoard = 256;
	float m1[] = { 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1 };
	GSPAT::Solution* solution = solver->createSolution(1,1,phaseOnly, position, amplitudes,m1,m1);
	solver->compute(solution); 
	//c. Retrieve solution field, make phase-only and apply signature
	float* field=solution->finalHologramReIm();
	for(int t=0;t<2*numTransducersPerBoard; t++)
	{
		float magnitude;
		if(phaseOnly)magnitude= sqrtf(field[2 * t] * field[2 * t] + field[2 * t + 1] * field[2 * t + 1]);
		else magnitude = 1;
		//Include Signature:
		//transducerState[t].real= (t>=numTransducersPerBoard?-1:1)*field[2*t]/magnitude;		//Real +signature
		//transducerState[t].imag= (t>=numTransducersPerBoard?-1:1)*field[2*t+1]/magnitude;	//Imag + signature 			
		//Avoid signature:
		transducerState[t].real= field[2*t]/magnitude;		//Real +signature
		transducerState[t].imag= field[2*t+1]/magnitude;	//Imag + signature 		
		
	}
	solver->releaseSolution(solution);
}


float testPoint(float position[], GSPAT::Solver* solver, float* transducerPositions, float sliceSize=0.03f, float numSlices=2) {
	//0. Generate field.
	lapack_complex_float transducerState[512];
	generateField(position, solver, transducerPositions, transducerState);
	//1. Visualize:
	/*{
		float A[] = {-0.075f, 0, 0.12f + 0.075f }, B[] = {0.075f, 0, 0.12f + 0.075f }, C[] = {-0.075f, 0, 0.12f - 0.075f};
		int imgRes[] = { 128,128 };
		VisualizePlane::visualize(A,B,C, imgRes , (float*)transducerState , 512, transducerPositions);
	}*/
	//2.Compute Amplitude:
	lapack_complex_float field_P = propagateFieldToPoint(position, transducerState, transducerPositions, 512,8.0f);
	return sqrtf(field_P.real*field_P.real + field_P.imag*field_P.imag);
}


void main(void) {
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
	//Configure the solver
	float transducerPositions[numBoards*numTransducersPerBoard * 3], amplitudeAdjust[numBoards*numTransducersPerBoard];
	int mappings[numBoards*numTransducersPerBoard], phaseDelays[numBoards*numTransducersPerBoard], numDiscreteLevels;
	driver->readParameters(transducerPositions, mappings, phaseDelays, amplitudeAdjust, &numDiscreteLevels);
	solver->setBoardConfig(transducerPositions, mappings, phaseDelays, amplitudeAdjust, numDiscreteLevels);
	//Generate data
	FILE* dataFile;  fopen_s(&dataFile, "TrapAmplitudes.csv", "w");
	fprintf(dataFile, "X, Y, Z, Amplitude(Pa)\n");
	for (float x = 0.0f; x <= 0.075f; x += 0.001f) {
		float position[] = { x,0,0.12f };
		float amplitude = testPoint(position, solver, transducerPositions);
		fprintf(dataFile,"%f, %f, %f, %f\n", position[0], position[1], position[2], amplitude);
	}	
	fclose(dataFile);
}


void saveForcesAroundTrap(char* fileName, float*transducerState, int boardSize[2], float pitch, float A[3], float B[3], float C[3], int imageSize[2]) {
	//0. Create output files (data file, force buffers and RGB images).
	FILE* dataFile;  fopen_s(&dataFile, fileName, "w");
	float* forcesX = new float[imageSize[1] * imageSize[0]];
	float* forcesY = new float[imageSize[1] * imageSize[0]];
	float* forcesZ = new float[imageSize[1] * imageSize[0]];
	float maxAbsForceX = 0, maxAbsForceY = 0, maxAbsForceZ = 0;
	CImg<unsigned char> imgX = CImg<unsigned char>(imageSize[0], imageSize[1], 1, 1);
	CImg<unsigned char> imgY = CImg<unsigned char>(imageSize[0], imageSize[1], 1, 1);
	CImg<unsigned char> imgZ = CImg<unsigned char>(imageSize[0], imageSize[1], 1, 1);


	//1.Initialize iterator vectors.
	float step_AB[3], step_AC[3];
	step_AB[0] = (B[0] - A[0]) / imageSize[0]; //Step vector in direction AB
	step_AB[1] = (B[1] - A[1]) / imageSize[0];
	step_AB[2] = (B[2] - A[2]) / imageSize[0];
	step_AC[0] = (C[0] - A[0]) / imageSize[1]; //Step vector in direction AC
	step_AC[1] = (C[1] - A[1]) / imageSize[1];
	step_AC[2] = (C[2] - A[2]) / imageSize[1];
	//2. Compute forces for each pixel in the image
	for(int py=0;py<imageSize[1]; py++)
		for (int px = 0; px < imageSize[0]; px++) {
			//2.a Compute 3D coordinates of the pixel
			float pixelCoords[3];
			pixelCoords[0] = A[0] + px*step_AB[0] + py*step_AC[0];
			pixelCoords[1] = A[1] + px*step_AB[1] + py*step_AC[1];
			pixelCoords[2] = A[2] + px*step_AB[2] + py*step_AC[2];
			//2.b Compute forces:
			float force[3];
			GorKovComputation::computeAcousticForceAtPoint(force, pixelCoords, (_lapack_complex_float*)transducerState, boardSize, pitch);
			//2.c. Store raw output:
			fprintf_s(dataFile,"%f,%f,%f,%f,%f,%f,\n", pixelCoords[0], pixelCoords[1], pixelCoords[2], force[0], force[1], force[2]);
			//2.d. Store in force buffers, and get preprocessing info (max forces will be mapped to maz colour intensity)
			forcesX[py*imageSize[0] + px] = force[0];
			forcesY[py*imageSize[0] + px] = force[1];
			forcesZ[py*imageSize[0] + px] = force[2];			
			if (fabs(force[0]) > maxAbsForceX)maxAbsForceX = fabs(force[0]);
			if (fabs(force[1]) > maxAbsForceY)maxAbsForceY = fabs(force[1]);
			if (fabs(force[2]) > maxAbsForceZ)maxAbsForceZ = fabs(force[2]);		
		}
	fclose(dataFile);
		//3. Save auxiliary images: 
		for(int py=0;py<imageSize[1]; py++)
			for (int px = 0; px < imageSize[0]; px++) {
				if(maxAbsForceX>0)
					imgX(px, py, 0) = (unsigned char)(128 + 127 * forcesX[py*imageSize[0] + px] / maxAbsForceX);
				if(maxAbsForceY>0)
					imgY(px, py, 0) = (unsigned char)(128 + 127 * forcesY[py*imageSize[0] + px] / maxAbsForceY);
				if(maxAbsForceZ>0)
					imgZ(px, py, 0) = (unsigned char)(128 + 127 * forcesZ[py*imageSize[0] + px] / maxAbsForceZ);
			
			}
		//4. Show auxiliary images: 
		CImgDisplay x(imgX, "X forces", false);
		CImgDisplay y(imgY, "Y forces", false);
		CImgDisplay z(imgZ, "Z forces", false);
		imgX.save("Xforces_normalised.bmp");
		imgY.save("Yforces_normalised.bmp");
		imgZ.save("Zforces_normalised.bmp");
		char c;  scanf("%c", &c);

}