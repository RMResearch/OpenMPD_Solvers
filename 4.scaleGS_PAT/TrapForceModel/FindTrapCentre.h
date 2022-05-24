#pragma once
#include "Helper/GorKovComputation.h"
#include <cstring>

void FindTrapCentre(float trapPosition[], float curTrapCentre[], float curMinGorkov, float curGradDirection[], float curStepSize, int numStepsRemaining, _lapack_complex_float* transducerState, int boardSize[2], float pitch ) {
	//0. Normalise direction (should be externally initialised to direction of Gor'kov gadient)
	float normalisedDir[3];
	float amplitude = sqrtf(curGradDirection[0]*curGradDirection[0]+curGradDirection[1]*curGradDirection[1]+curGradDirection[2]*curGradDirection[2]);
	normalisedDir[0] = curGradDirection[0]/amplitude;
	normalisedDir[1] = curGradDirection[1]/amplitude;
	normalisedDir[2] = curGradDirection[2]/amplitude;
	//1. Compute next position, based on direction and step size 
	//   Please note that the force already represents the negative gradient ("downhill" direction)).
	float nextPosition[3];
	nextPosition[0] = curTrapCentre[0] + curStepSize*normalisedDir[0];
	nextPosition[1] = curTrapCentre[1] + curStepSize*normalisedDir[1];
	nextPosition[2] = curTrapCentre[2] + curStepSize*normalisedDir[2];

	//2. Assess new position
	float newGorkovPotential=GorKovComputation::computeGorKovAtPoint(nextPosition, transducerState, boardSize, pitch);

	//2.1. If Gorkov potential decreases, take it and call recursively from the new position
	if (newGorkovPotential < curMinGorkov) {
		printf("Step taken (%f)\n", newGorkovPotential);
		memcpy(curTrapCentre, nextPosition, 3 * sizeof(float));
		float newGradDirection[3]; //Please note that the force already represents the negative gradient ("downhill" direction).
		GorKovComputation::computeAcousticForceAtPoint(newGradDirection, nextPosition, transducerState, boardSize, pitch);
		FindTrapCentre(trapPosition, curTrapCentre, newGorkovPotential, newGradDirection, curStepSize, numStepsRemaining, transducerState, boardSize, pitch);
		return; 
	}
	//2.2. If it does not, reduce step size and N.
	numStepsRemaining--;
	curStepSize /= 2;

	//2.2.a. If N>0, try again (recursively)
	if (numStepsRemaining > 0) {
		return FindTrapCentre(trapPosition, curTrapCentre, curMinGorkov, curGradDirection, curStepSize, numStepsRemaining, transducerState, boardSize, pitch);		
	}
	//2.2.b Otherwise, return best solution found so far (this happens automatically through arguments). 
	return;
}

void FindTrapCentre(float trapPosition[], float trapCentre[], _lapack_complex_float* transducerState, int boardSize[2], float pitch) {
	//0. Initialize current position from the trap centre position.
	memcpy(trapCentre, trapPosition, 3 * sizeof(float));
	//1. Compute initialization parameters:
	float curMinGorkov = GorKovComputation::computeGorKovAtPoint(trapCentre, transducerState, boardSize, pitch);
	float gorkovGradient[3];
	GorKovComputation::computeAcousticForceAtPoint(gorkovGradient, trapCentre, transducerState, boardSize, pitch);
	float stepSize = lambda() / 4;
	int maxSteps = 6;
	//2. Launch implementation method.
	FindTrapCentre(trapPosition, trapCentre, curMinGorkov, gorkovGradient, stepSize, maxSteps, transducerState, boardSize, pitch);
}


void FindTrapCentre(float trapPosition[], float curTrapCentre[], float curMinGorkov, float curGradDirection[], float curStepSize, int numStepsRemaining, _lapack_complex_float* transducerState, float* transducerPositions, int numTransducers) {
	//0. Normalise direction (should be externally initialised to direction of Gor'kov gadient)
	float normalisedDir[3];
	float amplitude = sqrtf(curGradDirection[0]*curGradDirection[0]+curGradDirection[1]*curGradDirection[1]+curGradDirection[2]*curGradDirection[2]);
	normalisedDir[0] = curGradDirection[0]/amplitude;
	normalisedDir[1] = curGradDirection[1]/amplitude;
	normalisedDir[2] = curGradDirection[2]/amplitude;
	//1. Compute next position, based on direction and step size 
	//   Please note that the force already represents the negative gradient ("downhill" direction)).
	float nextPosition[3];
	nextPosition[0] = curTrapCentre[0] + curStepSize*normalisedDir[0];
	nextPosition[1] = curTrapCentre[1] + curStepSize*normalisedDir[1];
	nextPosition[2] = curTrapCentre[2] + curStepSize*normalisedDir[2];

	//2. Assess new position
	float newGorkovPotential=GorKovComputation::computeGorKovAtPoint(nextPosition, transducerState, transducerPositions, numTransducers);

	//2.1. If Gorkov potential decreases, take it and call recursively from the new position
	if (newGorkovPotential < curMinGorkov) {
		printf("Step taken (%f)\n", newGorkovPotential);
		memcpy(curTrapCentre, nextPosition, 3 * sizeof(float));
		float newGradDirection[3]; //Please note that the force already represents the negative gradient ("downhill" direction).
		GorKovComputation::computeAcousticForceAtPoint(newGradDirection, nextPosition, transducerState, transducerPositions, numTransducers);
		FindTrapCentre(trapPosition, curTrapCentre, newGorkovPotential, newGradDirection, curStepSize, numStepsRemaining, transducerState, transducerPositions, numTransducers);
		return; 
	}
	//2.2. If it does not, reduce step size and N.
	numStepsRemaining--;
	curStepSize /= 2;

	//2.2.a. If N>0, try again (recursively)
	if (numStepsRemaining > 0) {
		return FindTrapCentre(trapPosition, curTrapCentre, curMinGorkov, curGradDirection, curStepSize, numStepsRemaining, transducerState, transducerPositions, numTransducers);		
	}
	//2.2.b Otherwise, return best solution found so far (this happens automatically through arguments). 
	return;
}


void FindTrapCentre(float trapPosition[], float trapCentre[], _lapack_complex_float* transducerState, float* transducerPositions, int numTransducers) {
	//0. Initialize current position from the trap centre position.
	memcpy(trapCentre, trapPosition, 3 * sizeof(float));
	//1. Compute initialization parameters:
	float curMinGorkov = GorKovComputation::computeGorKovAtPoint(trapCentre, transducerState, transducerPositions, numTransducers);
	float gorkovGradient[3];
	GorKovComputation::computeAcousticForceAtPoint(gorkovGradient, trapCentre, transducerState, transducerPositions, numTransducers);
	float stepSize = lambda() / 4;
	int maxSteps = 6;
	//2. Launch implementation method.
	FindTrapCentre(trapPosition, trapCentre, curMinGorkov, gorkovGradient, stepSize, maxSteps, transducerState, transducerPositions, numTransducers);
}
