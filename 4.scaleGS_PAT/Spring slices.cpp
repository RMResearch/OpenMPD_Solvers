#include <GSPAT_SolverV2.h>
#include <AsierInho_V2.h>
#include <../Helper/GorKovComputation.h>
#include <../Helper/VisualizePlane.h>
#include "TrapForceModel\FileNameGenerators.h"
#include "TrapForceModel\FindTrapCentre.h"
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
		transducerState[t].real= (t>=numTransducersPerBoard?-1:1)*field[2*t]/magnitude;		//Real +signature
		transducerState[t].imag= (t>=numTransducersPerBoard?-1:1)*field[2*t+1]/magnitude;	//Imag + signature 			
	}
	solver->releaseSolution(solution);
}

void generateData(float trapPosition[], float angle, bool corrected,  float* transducerPositions, lapack_complex_float transducerState[512], float A_w[], float B_w[], float C_w[], float X_w[], float Y_w[], float Z_w[]) {
	static const size_t resolution = 200;
	int imageSize[] = { resolution,resolution };
	VisualizePlane::visualize(A_w, B_w, C_w, imageSize, (float*)transducerState, 512, transducerPositions);
	float* forcesX = new float[resolution*resolution];
	float* forcesZ = new float[resolution*resolution];
	float maxAbsForceX = 0, maxAbsForceZ=0; 
	//1.Initialize iterator vectors.
	float step_AB[3], step_AC[3];
	step_AB[0] = (B_w[0] - A_w[0]) / resolution; //Step vector in direction AB
	step_AB[1] = (B_w[1] - A_w[1]) / resolution;
	step_AB[2] = (B_w[2] - A_w[2]) / resolution;
	step_AC[0] = (C_w[0] - A_w[0]) / resolution; //Step vector in direction AC
	step_AC[1] = (C_w[1] - A_w[1]) / resolution;
	step_AC[2] = (C_w[2] - A_w[2]) / resolution;
	//2. Compute forces for each position in the slice (in local and world coords)
	for (int py = 0; py < resolution; py++) {
		for (int px = 0; px < resolution; px++) {
			//Compute 3D coordinates of the pixel
			float pixelCoords[3];
			pixelCoords[0] = A_w[0] + (px + 0.5f)*step_AB[0] + (py + 0.5f)*step_AC[0];
			pixelCoords[1] = A_w[1] + (px + 0.5f)*step_AB[1] + (py + 0.5f)*step_AC[1];
			pixelCoords[2] = A_w[2] + (px + 0.5f)*step_AB[2] + (py + 0.5f)*step_AC[2];
			//Compute force along X:
			float force[3];
			static float Ax = 1, Az = 1, Fx = -467.91877075f, Fz = 1307.830115f;
			force[0] = -(px - resolution / 2.0f) / resolution;
			force[2] = (py - resolution/2.0f)/resolution;

			//Store force (in world coords) in buffers
			forcesX[py*resolution + px] = force[0];
			forcesZ[py*resolution + px] = force[2];
			if (fabs(force[0]) > maxAbsForceX)
				maxAbsForceX = fabs(force[0]);
			if (fabs(force[2]) > maxAbsForceZ)
				maxAbsForceZ = fabs(force[2]);
		}
	}
	//3. Generate images
	CImg<unsigned char> imgX = CImg<unsigned char>(resolution  , resolution  , 1, 3);
	CImg<unsigned char> imgZ = CImg<unsigned char>(resolution  , resolution  , 1, 3);
	for (int py = 0; py < resolution; py++) {
		for (int px = 0; px < resolution; px++) {
			//Plot X force:
			imgX(px, py, 0) = imgX(px, py, 1) = imgX(px, py, 2) = 0;
			if(forcesX[py*resolution + px]>0)
				imgX(px, py, 1) = (unsigned char)(255*forcesX[py*resolution + px] / maxAbsForceX);
			else
				imgX(px, py, 0) = (unsigned char)(-255*forcesX[py*resolution + px] / maxAbsForceX);
			//Plot Z force:
			imgZ(px, py, 0) = imgZ(px, py, 1) = imgZ(px, py, 2) = 0;
			if(forcesZ[py*resolution + px]>0)
				imgZ(px, py, 1) = (unsigned char)(255*forcesZ[py*resolution + px] / maxAbsForceZ);
			else
				imgZ(px, py, 0) = (unsigned char)(-255*forcesZ[py*resolution + px] / maxAbsForceZ);	
		}
	}
	imgX.save("SpringSlice_forcesX.bmp");
	imgZ.save("SpringSlice_forcesZ.bmp");
}

void testPoint(float position[], GSPAT::Solver* solver, float* transducerPositions, float sliceSize=0.02f) {
	//0. Generate field.
	lapack_complex_float transducerState[512];
	generateField(position, solver, transducerPositions, transducerState);
	//1.Generate slices: 
	//a. local coordinates for the slices
	float A_l[] = { -sliceSize / 2 , 0, sliceSize / 2 };
	float B_l[] = {  sliceSize / 2 , 0, sliceSize / 2 };
	float C_l[] = { -sliceSize / 2 , 0,-sliceSize / 2 };
	//b. generate each slice, applying a rotation (around Z axis):
	int n = 0;
		//Rotation angle (cos and sin) 
		float angle = 0;
		float cos_angle = cosf(angle*M_PI / 180);
		float sin_angle = sinf(angle*M_PI / 180);
		//Corner coordinates in World space (after applying the rotation)
		float A_w[] = { cos_angle*A_l[0] - sin_angle*A_l[1] + position[0]
					  , sin_angle*A_l[0] + cos_angle*A_l[1] + position[1]
					  , A_l[2]								+ position[2]};
		float B_w[] = { cos_angle*B_l[0] - sin_angle*B_l[1] + position[0]
					  , sin_angle*B_l[0] + cos_angle*B_l[1] + position[1]
					  , B_l[2]								+ position[2]};
		float C_w[] = { cos_angle*C_l[0] - sin_angle*C_l[1] + position[0]
					  , sin_angle*C_l[0] + cos_angle*C_l[1] + position[1]
					  , C_l[2]								+ position[2]};
		//Local axis in world coordinates:
		float X_w[] = { cos_angle,  sin_angle, 0 };
		float Y_w[] = {-sin_angle,  cos_angle, 0 };
		float Z_w[] = { 0,0,1 };
		generateData(position, angle, false, transducerPositions, transducerState, A_w, B_w, C_w, X_w, Y_w, Z_w);
		//generateData(NULL, angle, false, NULL, pitch, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
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
	//Configura the solver
	float transducerPositions[numBoards*numTransducersPerBoard * 3], amplitudeAdjust[numBoards*numTransducersPerBoard];
	int mappings[numBoards*numTransducersPerBoard], phaseDelays[numBoards*numTransducersPerBoard], numDiscreteLevels;
	driver->readParameters(transducerPositions, mappings, phaseDelays, amplitudeAdjust, &numDiscreteLevels);
	solver->setBoardConfig(transducerPositions, mappings, phaseDelays, amplitudeAdjust, numDiscreteLevels);
	
	float position[] = { 0,0,0+0.2388f / 2 };
	testPoint(position, solver, transducerPositions);
		
}
