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

struct trapStats {
		float trapPos[3];
		float trapCentrePos[3];
		//Maximum force values needed to normalise.
		float maxAbsForceX_w = 0, maxAbsForceY_w = 0, maxAbsForceZ_w = 0;
		float maxAbsForceX_l = 0, maxAbsForceY_l = 0, maxAbsForceZ_l = 0;
		//Max and min forces (max positive and min negative)
		float maxX_l = 0, minX_l = 0, maxY_l = 0, minY_l = 0, maxZ_l = 0, minZ_l = 0;
		float maxX_w = 0, minX_w = 0, maxY_w = 0, minY_w = 0, maxZ_w = 0, minZ_w = 0;
		//Positions for these max/min forces
		float maxX_l_pos[3], maxY_l_pos[3], maxZ_l_pos[3], minX_l_pos[3], minY_l_pos[3], minZ_l_pos[3];
		float maxX_w_pos[3], maxY_w_pos[3], maxZ_w_pos[3], minX_w_pos[3], minY_w_pos[3], minZ_w_pos[3];
		//Distances between max/min points (compute at the end)
		float distMaxMinX_l = 0, distMaxMinY_l = 0, distMaxMinZ_l = 0, distMaxMinX_w = 0, distMaxMinY_w = 0, distMaxMinZ_w = 0;
};

void saveForcesAroundTrap(char* fileName, float*transducerState, int boardSize[2], float pitch, float A[3], float B[3], float C[3], int resolution[2]);
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

void updateStats(struct trapStats& stats, float pos[3], float fX_l, float fY_l, float fZ_l, float fX_w, float fY_w, float fZ_w) {
	//Update absolute forces in world
	if (fabs(fX_w) > stats.maxAbsForceX_w)stats.maxAbsForceX_w = fabs(fX_w);
	if (fabs(fY_w) > stats.maxAbsForceY_w)stats.maxAbsForceY_w = fabs(fY_w);
	if (fabs(fZ_w) > stats.maxAbsForceZ_w)stats.maxAbsForceZ_w = fabs(fZ_w);
	//Update absolute forces in local			
	if (fabs(fX_l) > stats.maxAbsForceX_l)stats.maxAbsForceX_l = fabs(fX_l);
	if (fabs(fY_l) > stats.maxAbsForceY_l)stats.maxAbsForceY_l = fabs(fY_l);
	if (fabs(fZ_l) > stats.maxAbsForceZ_l)stats.maxAbsForceZ_l = fabs(fZ_l);
	//Update min forces in world
	if (fX_w < stats.minX_w) {stats.minX_w = fX_w; memcpy(stats.minX_w_pos, pos, 3 * sizeof(float));}
	if (fY_w < stats.minY_w) {stats.minY_w = fY_w; memcpy(stats.minY_w_pos, pos, 3 * sizeof(float));}
	if (fZ_w < stats.minZ_w) {stats.minZ_w = fZ_w; memcpy(stats.minZ_w_pos, pos, 3 * sizeof(float));}
	//Update max forces in world
	if (fX_w > stats.maxX_w) {stats.maxX_w = fX_w; memcpy(stats.maxX_w_pos, pos, 3 * sizeof(float));}
	if (fY_w > stats.maxY_w) {stats.maxY_w = fY_w; memcpy(stats.maxY_w_pos, pos, 3 * sizeof(float));}
	if (fZ_w > stats.maxZ_w) {stats.maxZ_w = fZ_w; memcpy(stats.maxZ_w_pos, pos, 3 * sizeof(float));}
	//Update min forces in local
	if (fX_l < stats.minX_l) {stats.minX_l = fX_l; memcpy(stats.minX_l_pos, pos, 3 * sizeof(float));}
	if (fY_l < stats.minY_l) {stats.minY_l = fY_l; memcpy(stats.minY_l_pos, pos, 3 * sizeof(float));}
	if (fZ_l < stats.minZ_l) {stats.minZ_l = fZ_l; memcpy(stats.minZ_l_pos, pos, 3 * sizeof(float));}
	//Update max forces in local
	if (fX_l > stats.maxX_l) {stats.maxX_l = fX_l; memcpy(stats.maxX_l_pos, pos, 3 * sizeof(float));}
	if (fY_l > stats.maxY_l) {stats.maxY_l = fY_l; memcpy(stats.maxY_l_pos, pos, 3 * sizeof(float));}
	if (fZ_l > stats.maxZ_l) {stats.maxZ_l = fZ_l; memcpy(stats.maxZ_l_pos, pos, 3 * sizeof(float));}
	

}


void updateStatsDistances(struct trapStats& stats) {
	stats.distMaxMinX_l = sqrtf( (stats.maxX_l_pos[0]-stats.minX_l_pos[0]) * (stats.maxX_l_pos[0]-stats.minX_l_pos[0])
								/*+(stats.maxX_l_pos[1]-stats.minX_l_pos[1]) * (stats.maxX_l_pos[1]-stats.minX_l_pos[1])
								+(stats.maxX_l_pos[2]-stats.minX_l_pos[2]) * (stats.maxX_l_pos[2]-stats.minX_l_pos[2])*/);
	stats.distMaxMinY_l = sqrtf( /*(stats.maxY_l_pos[0]-stats.minY_l_pos[0]) * (stats.maxY_l_pos[0]-stats.minY_l_pos[0])
								+*/(stats.maxY_l_pos[1]-stats.minY_l_pos[1]) * (stats.maxY_l_pos[1]-stats.minY_l_pos[1])
								/*+(stats.maxY_l_pos[2]-stats.minY_l_pos[2]) * (stats.maxY_l_pos[2]-stats.minY_l_pos[2])*/);
	stats.distMaxMinZ_l = sqrtf( /*(stats.maxZ_l_pos[0]-stats.minZ_l_pos[0]) * (stats.maxZ_l_pos[0]-stats.minZ_l_pos[0])
								+(stats.maxZ_l_pos[1]-stats.minZ_l_pos[1]) * (stats.maxZ_l_pos[1]-stats.minZ_l_pos[1])
								+*/(stats.maxZ_l_pos[2]-stats.minZ_l_pos[2]) * (stats.maxZ_l_pos[2]-stats.minZ_l_pos[2]));
	stats.distMaxMinX_w = sqrtf( (stats.maxX_w_pos[0]-stats.minX_w_pos[0]) * (stats.maxX_w_pos[0]-stats.minX_w_pos[0])
								/*+(stats.maxX_w_pos[1]-stats.minX_w_pos[1]) * (stats.maxX_w_pos[1]-stats.minX_w_pos[1])
								+(stats.maxX_w_pos[2]-stats.minX_w_pos[2]) * (stats.maxX_w_pos[2]-stats.minX_w_pos[2])*/);
	stats.distMaxMinY_w = sqrtf( /*(stats.maxY_w_pos[0]-stats.minY_w_pos[0]) * (stats.maxY_w_pos[0]-stats.minY_w_pos[0])
								+*/(stats.maxY_w_pos[1]-stats.minY_w_pos[1]) * (stats.maxY_w_pos[1]-stats.minY_w_pos[1])
								/*+(stats.maxY_w_pos[2]-stats.minY_w_pos[2]) * (stats.maxY_w_pos[2]-stats.minY_w_pos[2])*/);
	stats.distMaxMinZ_w = sqrtf( /*(stats.maxZ_w_pos[0]-stats.minZ_w_pos[0]) * (stats.maxZ_w_pos[0]-stats.minZ_w_pos[0])
								+(stats.maxZ_w_pos[1]-stats.minZ_w_pos[1]) * (stats.maxZ_w_pos[1]-stats.minZ_w_pos[1])
								+*/(stats.maxZ_w_pos[2]-stats.minZ_w_pos[2]) * (stats.maxZ_w_pos[2]-stats.minZ_w_pos[2]));

}


void generateData(float trapPosition[], float angle, bool corrected,  float* transducerPositions, lapack_complex_float transducerState[512], float A_w[], float B_w[], float C_w[], float X_w[], float Y_w[], float Z_w[]) {
	static bool firstTime = true;
	static const size_t resolutionX = 64;
	static const size_t resolutionY = 64;
	int imageSize[] = { resolutionX,resolutionY };
	VisualizePlane::visualize(A_w, B_w, C_w, imageSize, (float*)transducerState, 512, transducerPositions);
	//Buffers storing forces per axis (local and world axis... only local should be relevant)
	float* positions=new float[3 * resolutionX*resolutionY];
	float* forcesX_l=new float[resolutionX  * resolutionY ];
	float* forcesY_l=new float[resolutionX  * resolutionY ];
	float* forcesZ_l=new float[resolutionX  * resolutionY ];
	float* forcesX_w=new float[resolutionX  * resolutionY ];
	float* forcesY_w=new float[resolutionX  * resolutionY ];
	float* forcesZ_w=new float[resolutionX  * resolutionY ];
	//Stats to obtain:
	struct trapStats _stats;
	memcpy(_stats.trapPos, trapPosition, 3 * sizeof(float));
	FindTrapCentre(trapPosition, _stats.trapCentrePos, transducerState, transducerPositions, 512);
	//1.Initialize iterator vectors.
	float step_AB[3], step_AC[3];
	step_AB[0] = (B_w[0] - A_w[0]) / resolutionX; //Step vector in direction AB
	step_AB[1] = (B_w[1] - A_w[1]) / resolutionX;
	step_AB[2] = (B_w[2] - A_w[2]) / resolutionX;
	step_AC[0] = (C_w[0] - A_w[0]) / resolutionY; //Step vector in direction AC
	step_AC[1] = (C_w[1] - A_w[1]) / resolutionY;
	step_AC[2] = (C_w[2] - A_w[2]) / resolutionY;
	//2. Compute forces for each position in the slice (in local and world coords)
	for(int py=0;py<resolutionY; py++)
		for (int px = 0; px < resolutionX; px++) {
			//2.a Compute 3D coordinates of the pixel
			float pixelCoords[3];
			//USED FOR ORIGINAL DATA
			/*pixelCoords[0] = A_w[0] + (px)*step_AB[0] + (py)*step_AC[0];
			pixelCoords[1] = A_w[1] + (px)*step_AB[1] + (py)*step_AC[1];
			pixelCoords[2] = A_w[2] + (px)*step_AB[2] + (py)*step_AC[2];*/

			//CORRECT WAY (CENTRE)
			pixelCoords[0] = A_w[0] + (px+0.5f)*step_AB[0] + (py+0.5f)*step_AC[0];
			pixelCoords[1] = A_w[1] + (px+0.5f)*step_AB[1] + (py+0.5f)*step_AC[1];
			pixelCoords[2] = A_w[2] + (px+0.5f)*step_AB[2] + (py+0.5f)*step_AC[2];
			memcpy(&(positions[3 * (py*resolutionX + px)]), pixelCoords, 3 * sizeof(float));
			//2.b Compute forces:
			float force[3];
			GorKovComputation::computeAcousticForceAtPoint(force, pixelCoords, (_lapack_complex_float*)transducerState, transducerPositions,512);
			//2.d. Store force (in world coords) in buffers, and get preprocessing info (max forces will be mapped to maz colour intensity)
			forcesX_w[py*resolutionX + px] = force[0];
			forcesY_w[py*resolutionX + px] = force[1];
			forcesZ_w[py*resolutionX + px] = force[2];			
			//2.e. Store force (in local coords) in buffers, and get preprocessing info (max forces will be mapped to maz colour intensity)
			forcesX_l[py*resolutionX + px] = force[0]*X_w[0] + force[1]*X_w[1] + force[2]*X_w[2];
			forcesY_l[py*resolutionX + px] = force[0]*Y_w[0] + force[1]*Y_w[1] + force[2]*Y_w[2];
			forcesZ_l[py*resolutionX + px] = force[0]*Z_w[0] + force[1]*Z_w[1] + force[2]*Z_w[2];			
			//2.f. Update stats
			updateStats(_stats, pixelCoords, forcesX_l[py*resolutionX + px], forcesY_l[py*resolutionX + px], forcesZ_l[py*resolutionX + px], forcesX_w[py*resolutionX + px], forcesY_w[py*resolutionX + px], forcesZ_w[py*resolutionX + px]);			
		}

	//3. Store data in files:
	//3.a. Generate filenames
	char dataFileName[512], forcesFileName[512], imgX_l_FileName[512], imgY_l_FileName[512], imgZ_l_FileName[512], imgX_w_FileName[512], imgY_w_FileName[512], imgZ_w_FileName[512];
	FILE* dataFile, *forcesFile;  
	FileNameGenerator::generateTrapDataFileName(dataFileName, corrected, trapPosition[0], trapPosition[1], trapPosition[2], angle);
	FileNameGenerator::generateTrapForcesFileName(forcesFileName, corrected, trapPosition[0], trapPosition[1], trapPosition[2], angle);
	fopen_s(&dataFile, dataFileName, (firstTime? "w":"a"));//We clear the summary file first time we run this
	fopen_s(&forcesFile, forcesFileName, "w");
	FileNameGenerator::generateTrapImageFileName(imgX_w_FileName, corrected, trapPosition[0], trapPosition[1], trapPosition[2], angle, 'X');
	FileNameGenerator::generateTrapImageFileName(imgY_w_FileName, corrected, trapPosition[0], trapPosition[1], trapPosition[2], angle, 'Y');
	FileNameGenerator::generateTrapImageFileName(imgZ_w_FileName, corrected, trapPosition[0], trapPosition[1], trapPosition[2], angle, 'Z');
	FileNameGenerator::generateTrapImageFileName(imgX_l_FileName, corrected, trapPosition[0], trapPosition[1], trapPosition[2], angle, 'U');
	FileNameGenerator::generateTrapImageFileName(imgY_l_FileName, corrected, trapPosition[0], trapPosition[1], trapPosition[2], angle, 'V');
	FileNameGenerator::generateTrapImageFileName(imgZ_l_FileName, corrected, trapPosition[0], trapPosition[1], trapPosition[2], angle, 'W');
	//3.b. images (in both local and world axis). 
	CImg<unsigned char> imgX_l = CImg<unsigned char>(resolutionX  , resolutionY  , 1, 1);
	CImg<unsigned char> imgY_l = CImg<unsigned char>(resolutionX  , resolutionY  , 1, 1);
	CImg<unsigned char> imgZ_l = CImg<unsigned char>(resolutionX  , resolutionY  , 1, 1);
	CImg<unsigned char> imgX_w = CImg<unsigned char>(resolutionX  , resolutionY  , 1, 1);
	CImg<unsigned char> imgY_w = CImg<unsigned char>(resolutionX  , resolutionY  , 1, 1);
	CImg<unsigned char> imgZ_w = CImg<unsigned char>(resolutionX  , resolutionY  , 1, 1);
	//3.c. Store data (per point)
	for(int py=0;py<resolutionY; py++)
		for (int px = 0; px < resolutionX; px++) {
			if(_stats.maxAbsForceX_w>0)
				imgX_w(px, py, 0) = (unsigned char)(128 + 127 * forcesX_w[py*resolutionX + px] / _stats.maxAbsForceX_w);
			if(_stats.maxAbsForceY_w>0)
				imgY_w(px, py, 0) = (unsigned char)(128 + 127 * forcesY_w[py*resolutionX + px] / _stats.maxAbsForceY_w);
			if(_stats.maxAbsForceZ_w>0)
				imgZ_w(px, py, 0) = (unsigned char)(128 + 127 * forcesZ_w[py*resolutionX + px] / _stats.maxAbsForceZ_w);
			if(_stats.maxAbsForceX_l>0)
				imgX_l(px, py, 0) = (unsigned char)(128 + 127 * forcesX_l[py*resolutionX + px] / _stats.maxAbsForceX_l);
			if(_stats.maxAbsForceY_l>0)
				imgY_l(px, py, 0) = (unsigned char)(128 + 127 * forcesY_l[py*resolutionX+ px] / _stats.maxAbsForceY_l);
			if(_stats.maxAbsForceZ_l>0)
				imgZ_l(px, py, 0) = (unsigned char)(128 + 127 * forcesZ_l[py*resolutionX+ px] / _stats.maxAbsForceZ_l);
			fprintf_s(forcesFile,"%f,%f,%f,%f,%f,%f,%f,%f,%f,\n", positions[3*(py*resolutionX + px)], positions[3*(py*resolutionX + px)+1], positions[3*(py*resolutionX + px)+2], forcesX_l[py*resolutionX + px], forcesY_l[py*resolutionX + px], forcesZ_l[py*resolutionX + px], forcesX_w[py*resolutionX + px], forcesY_w[py*resolutionX + px], forcesZ_w[py*resolutionX + px]);
			
		}
	//3.d. Store stats summary: 
	updateStatsDistances(_stats);
	if (firstTime) {//Add the header (only) the first time we write to the file
		fprintf_s(dataFile, "trap.X,trap.Y,trap.Z,centre.X,centre.Y,centre.Z,minfX_l, maxfX_l,minfY_l, maxfYl, minfZ_l, maxfZ_l,min_fX_w, maxfX_w,min_fY_w, maxfY_w, min_fZ_w, maxfZ_w, distMinMaxX_l, distMinMaxY_l, distMinMaxZ_l,distMinMaxX_w, distMinMaxY_w, distMinMaxZ_w,\n");
		firstTime = false;
	}
	fprintf_s(dataFile,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,\n", _stats.trapPos[0],_stats.trapPos[1],_stats.trapPos[2],_stats.trapCentrePos[0],_stats.trapCentrePos[1],_stats.trapCentrePos[2],_stats.minX_l, _stats.maxX_l,_stats.minY_l, _stats.maxY_l, _stats.minZ_l, _stats.maxZ_l,_stats.minX_w,_stats.maxX_w,_stats.minY_w, _stats.maxY_w, _stats.minZ_w, _stats.maxZ_w, _stats.distMaxMinX_l, _stats.distMaxMinY_l, _stats.distMaxMinZ_l,_stats.distMaxMinX_w, _stats.distMaxMinY_w, _stats.distMaxMinZ_w);
	fclose(dataFile);
	fclose(forcesFile);
	imgX_l.save(imgX_l_FileName);
	imgY_l.save(imgY_l_FileName);
	imgZ_l.save(imgZ_l_FileName);
	imgX_w.save(imgX_w_FileName);
	imgY_w.save(imgY_w_FileName);
	imgZ_w.save(imgZ_w_FileName);

	delete positions;
	delete forcesX_l;
	delete  forcesY_l;
	delete forcesZ_l;
	delete forcesX_w;
	delete forcesY_w;
	delete forcesZ_w;
}

void testPoint(float position[], GSPAT::Solver* solver, float* transducerPositions, float sliceSize=0.03f, float numSlices=2) {
	//0. Generate field.
	lapack_complex_float transducerState[512];
	generateField(position, solver, transducerPositions, transducerState);
	//1.Generate slices: 
	//a. local coordinates for the slices
	float A_l[] = { -sliceSize / 2 , 0, sliceSize / 2 };
	float B_l[] = {  sliceSize / 2 , 0, sliceSize / 2 };
	float C_l[] = { -sliceSize / 2 , 0,-sliceSize / 2 };
	//b. generate each slice, applying a rotation (around Z axis):
	for (int n = 0; n < numSlices; n++) {
		//Rotation angle (cos and sin) 
		float angle = 180.0f*n / numSlices;
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
	
	
	for(float z=0.0f; z<=0.04f; z+=0.01f)
		for(float y=0.0f; y<=0.04f; y+=0.01f)
			for (float x = 0.0f; x <= 0.04f; x += 0.01f) {
				float position[] = { x,y,z+0.2388f / 2 };
				//float position[] = { 0,0,0+0.2388f / 2 };
				testPoint(position, solver, transducerPositions);
			}	
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