#pragma once
#include <../Helper/TimeFunctions.h>

/**
	Helper class to time how much the different stages of the solver take. 
	It records the times and provides time elapsed for each stage. 
	Check HologramSolver to see how it is used:
		- times are recorded during computation (see HologramSolver::compute)
		- time elapsed is retrieved during HologramSolver::printStats 
*/
class HologramProfiler {
protected:	
	struct timeval cycleStarted, sceneUpdated, sceneRendered, hologramComputed, hologramDiscretized, boardUpdated;
	int numSamples;
	float _sceneUpdateTimeAVG, _renderTimeAVG,  _hologramComputationAVG, _hologramDiscretizationAVG, _boardUpdateTimeAVG, _totalTimeAVG;
public: 
	HologramProfiler()
		: numSamples(0)
		, _sceneUpdateTimeAVG(0)
		, _renderTimeAVG(0)		
		, _hologramComputationAVG(0)
		, _hologramDiscretizationAVG(0)
		, _boardUpdateTimeAVG(0)
		, _totalTimeAVG(0)
	{
		;
	}
	//Record events:
	void recordStartCycle();
	void recordSceneUpdated();
	void recordSceneRendered();
	void recordHologramComputed();
	void recordHologramDiscretized();
	void recordBoardUpdated();
	//Get times for each stage in seconds (only last stage):
	float getSceneUpdateTime();
	float getRenderTime();
	float getHologramComputationTime();
	float getHologramDiscretizationTime();
	float getBoardUpdateTime();
	float getTotalTime();
	//Get average times (over the whole application lifetime)
	inline float getAVGSceneUpdateTime() {
		return _sceneUpdateTimeAVG;
	}
	inline float getAVGRenderTime() {
		return _renderTimeAVG;
	}	
	inline float getAVGHologramComputationTime() {
		return _hologramComputationAVG;
	}
	inline float getAVGHologramDiscretizationTime() {
		return _hologramDiscretizationAVG;
	}
	inline float getAVGBoardUpdateTime() {
		return _boardUpdateTimeAVG;
	}
	inline float getAVGTotalTime() {
		return _totalTimeAVG;
	}
};