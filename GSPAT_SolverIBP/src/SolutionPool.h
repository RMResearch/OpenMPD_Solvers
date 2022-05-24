#pragma once
#include "HologramSolution.h"
#include <pthread.h>
#include <deque>
#include <vector>

class HologramSolverCL;

class SolutionPool {
	static const int MAX_SOLUTIONS = 4;
	HologramSolverCL* solver;
	HologramSolution* solutionPool[MAX_SOLUTIONS];
	std::deque<Solution*> availableSolutions;
	pthread_mutex_t mutexSolutionPool;
	pthread_mutex_t newSolutionsAvailable;

public: 
	SolutionPool(HologramSolverCL* solver, int numTransducers, cl_context context, cl_command_queue queue);
	~SolutionPool();
	HologramSolution* getNextSolution();	
	void returnSolution(HologramSolution* solution);
};
