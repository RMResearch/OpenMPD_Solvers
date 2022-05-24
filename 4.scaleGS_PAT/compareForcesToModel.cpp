#include "HologramSolver\HologramSolverOurs.h"
#include "HologramSolver/HologramSolverAsier.h"
#include "HologramSolver/HologramSolverNaive.h"
#include "Helper/GorKovComputation.h"
#include "TrapForceModel\FileNameGenerators.h"
#include "TrapForceModel\FindTrapCentre.h"

void computeForcesFromModel(float trapPosition[3], float particlePosition[3], float force_w[3], float &dZ, float &dR) {
	dZ = trapPosition[2] - particlePosition[2];
	dR = sqrtf((trapPosition[0] - particlePosition[0])*(trapPosition[0] - particlePosition[0])+(trapPosition[1] - particlePosition[1])*(trapPosition[1] - particlePosition[1]));
	float theta = atan2f((trapPosition[1] - particlePosition[1]), (trapPosition[0] - particlePosition[0]));
	static const float A_r = -0.000729;
	static const float A_z = -0.004405;
	//Build quadratics: 
	float Qx = 1 - 625 * (1 - 0.790123f)*trapPosition[0]*trapPosition[0];
	float Qy = 1 - 625 * (1 - 0.75583f)*trapPosition[1]*trapPosition[1];
	float Qz = 1 - 625 * (1 - 0.979424f)*(trapPosition[2]-0.1194)*(trapPosition[2]-0.1194);
	float Rx = 1 - 625 * (1 - 0.819781f)*trapPosition[0]*trapPosition[0];
	float Ry = 1 - 625 * (1 - 0.819781f)*trapPosition[1]*trapPosition[1];
	float Rz = 1 - 625 * (1 - 0.906087f)*(trapPosition[2]-0.1194)*(trapPosition[2]-0.1194);
	//Compute forces (vertical and radial)
	//float fZ = A_z*Rx*Ry*Rz*sinf(-2 * M_PI*dZ / 0.00492f)*cos(2 * M_PI*dR / 0.02527f);
	float fZ = A_z*Rx*Ry*Rz*sinf(-2 * M_PI*dZ / 0.00475f)*cos(2 * M_PI*dR / 0.021f);
	float fR = A_r*Qx*Qy*Qz*sinf(-2 * M_PI*dR / 0.01325f)*cos(2 * M_PI*dZ / 0.00492f);
	//Transform forces to euclidean coordinates:
	force_w[0] = fR * cosf(theta);
	force_w[1] = fR * sinf(theta);
	force_w[2] = fZ;



}

void compareData(FILE* compareFile, float trapPosition[], float angle, bool corrected, float X_w[], float Y_w[], float Z_w[]) {
	static bool firstTime = true;
	//1. Open file for the current trapPosition and angle
	char forcesFileName[512];
	FILE* forcesFile;  
	FileNameGenerator::generateTrapForcesFileName(forcesFileName, corrected, trapPosition[0], trapPosition[1], trapPosition[2], angle);
	fopen_s(&forcesFile, forcesFileName, "r");
	//2. Read each entry and compare it to the model
	while (!feof(forcesFile)) {
		//2.1. Reading data
		float position_w[3], forceGorkov_w[3], forceGorkov_l[3];
		fscanf(forcesFile, "%f,%f,%f,%f,%f,%f,%f,%f,%f,\n", &(position_w[0]), &(position_w[1]), &(position_w[2])
															 , &(forceGorkov_l[0]), &(forceGorkov_l[1]), &(forceGorkov_l[2]) 	
															 , &(forceGorkov_w[0]), &(forceGorkov_w[1]), &(forceGorkov_w[2]));
		//2.2. Compute forces from model: 
		float forceModel_w[3], forceModel_l[3], distZ, distR;
		computeForcesFromModel(trapPosition, position_w, forceModel_w, distZ, distR);
		forceModel_l[0] = forceModel_w[0]*X_w[0] + forceModel_w[1]*X_w[1] + forceModel_w[2]*X_w[2];
		forceModel_l[1] = forceModel_w[0]*Y_w[0] + forceModel_w[1]*Y_w[1] + forceModel_w[2]*Y_w[2];
		forceModel_l[2] = forceModel_w[0]*Z_w[0] + forceModel_w[1]*Z_w[1] + forceModel_w[2]*Z_w[2];

		//2.3. Compare and store: 
		if(distZ<0.002f && distZ>-0.002f && distR<0.004)
			fprintf(compareFile,"%f,%f,%f,  %f,%f,%f,  %f,%f,  %f,%f,%f,  %f,%f,%f,  %f,%f,%f,  %f,%f,%f,\n",
																trapPosition[0],trapPosition[1],trapPosition[2],
																position_w[0], position_w[1], position_w[2],
																distZ, distR,
																forceGorkov_w[0],forceGorkov_w[1], forceGorkov_w[2],
																forceGorkov_l[0],forceGorkov_l[1], forceGorkov_l[2],
																forceModel_w[0],forceModel_w[1], forceModel_w[2],
																forceModel_l[0],forceModel_l[1], forceModel_l[2]																
																);	
	}
	fclose(forcesFile);
}

void testPoint(FILE* compareFile, float position[], float numSlices=2) {
	//b. generate each slice, applying a rotation (around Z axis):
	for (int n = 0; n < numSlices; n++) {
		//Rotation angle (cos and sin) 
		float angle = 180.0f*n / numSlices;
		float cos_angle = cosf(angle*M_PI / 180);
		float sin_angle = sinf(angle*M_PI / 180);		
		//Local axis in world coordinates:
		float X_w[] = { cos_angle,  sin_angle, 0 };
		float Y_w[] = {-sin_angle,  cos_angle, 0 };
		float Z_w[] = { 0,0,1 };
		compareData(compareFile, position, angle, false, X_w, Y_w, Z_w);
	}
}


void main(void) {
	FILE* compareFile;
	fopen_s(&compareFile,"compareModelToGorkov.csv", "w");
	fprintf(compareFile, "x,y,z, (x+dx), (y+dy), (z+dz), dZ, dR, GKx_w,GKy_w,GKz_w, GKx_l,GKy_l,GKz_l, Mx_w,My_w,Mz_w, Mx_l,My_l,Mz_l,\n");
	
	/*for(float z=-0.04f; z<=0.04f; z+=0.01f)
		for(float y=-0.04f; y<=0.04f; y+=0.01f)
			for (float x = -0.04f; x <= 0.04f; x += 0.01f) {
	*/for(float z=	0.0f; z<=0.04f; z+=0.01f)
		for(float y=0.0f; y<=0.04f; y+=0.01f)
			for (float x = 0.0f; x <= 0.04f; x += 0.01f) {
				float position[] = { x,y,z+0.2388f / 2 };				
				testPoint(compareFile, position);
			}	
	fclose(compareFile);
}

