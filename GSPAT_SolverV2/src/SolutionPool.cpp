#include <src/SolutionPool.h>
#include <src/HologramSolverCL.h>

SolutionPool::SolutionPool(HologramSolverCL* solver, int numTransducers, cl_context context, cl_command_queue queue) {
	this->solver = solver;
	//Create repository of solutions:
	pthread_mutex_init(&mutexSolutionPool, NULL);
	pthread_mutex_init(&newSolutionsAvailable, NULL);
	pthread_mutex_lock(&newSolutionsAvailable);
	for (int i = 0; i < MAX_SOLUTIONS; i++) {
		solutionPool [i]= new HologramSolution(numTransducers, context, queue);
		availableSolutions.push_back(solutionPool[i]);
	}	
}

SolutionPool::~SolutionPool() {	
	pthread_mutex_lock(&mutexSolutionPool);
	for (int s = 0; s < MAX_SOLUTIONS; s++)
		solver->destroySolution(solutionPool[s]);
	availableSolutions.clear();
	pthread_mutex_destroy(&mutexSolutionPool);
	pthread_mutex_destroy(&newSolutionsAvailable);
}

HologramSolution* SolutionPool::getNextSolution() {
	pthread_mutex_lock(&mutexSolutionPool);
	while (availableSolutions.size() == 0) {
		pthread_mutex_unlock(&mutexSolutionPool); //let others access the solution pool
		//sprintf(consoleLine,"Waiting for empty slot to write!\n");
		pthread_mutex_lock(&newSolutionsAvailable);//while we wait for new solutions
		pthread_mutex_lock(&mutexSolutionPool);		//...and check again		
	}
	//...a solution is available. Let's get it!
	HologramSolution* nextSolution = (HologramSolution*) (availableSolutions.front());
	availableSolutions.pop_front();
	pthread_mutex_unlock(&mutexSolutionPool);
	return nextSolution;
}

void SolutionPool::returnSolution(HologramSolution* solution) {
	pthread_mutex_lock(&mutexSolutionPool);
	availableSolutions.push_back(solution);
	pthread_mutex_unlock(&newSolutionsAvailable);//Notify new solutions are available (in case the pool was empty). 
	pthread_mutex_unlock(&mutexSolutionPool);
}