#include <src/HologramProfiler.h>
#include <../Helper/TimeFunctions.h>
#include <stdlib.h>
#include <stdio.h>
#include <Windows.h>
#include <stdint.h>

void HologramProfiler::recordStartCycle() {
	gettimeofday(&cycleStarted,0x0);
}
void HologramProfiler::recordSceneUpdated(){
	gettimeofday(&sceneUpdated,0x0);
}
void HologramProfiler::recordSceneRendered(){
	gettimeofday(&sceneRendered,0x0);
}
void HologramProfiler::recordHologramComputed(){
	gettimeofday(&hologramComputed,0x0);
}
void HologramProfiler::recordHologramDiscretized(){
	gettimeofday(&hologramDiscretized,0x0);
}
void HologramProfiler::recordBoardUpdated() {
	gettimeofday(&boardUpdated, 0x0);
	//Update Average values: 
	numSamples++;
	_sceneUpdateTimeAVG = (((float)numSamples - 1)*_sceneUpdateTimeAVG + getSceneUpdateTime()) / numSamples;
	_renderTimeAVG = (((float)numSamples - 1)*_renderTimeAVG + getRenderTime()) / numSamples;
	_hologramComputationAVG = (((float)numSamples - 1)*_hologramComputationAVG + getHologramComputationTime()) / numSamples;
	_hologramDiscretizationAVG = (((float)numSamples - 1)*_hologramDiscretizationAVG + getHologramDiscretizationTime()) / numSamples;
	_boardUpdateTimeAVG = (((float)numSamples - 1)*_boardUpdateTimeAVG + getBoardUpdateTime()) / numSamples;
	_totalTimeAVG = (((float)numSamples - 1)*_totalTimeAVG + getTotalTime()) / numSamples;
}

float HologramProfiler::getSceneUpdateTime() {
	return computeTimeElapsed(cycleStarted, sceneUpdated);
}

float HologramProfiler::getRenderTime () {
	return computeTimeElapsed(sceneUpdated, sceneRendered);
}
float HologramProfiler::getHologramComputationTime() {
	return computeTimeElapsed(sceneRendered, hologramComputed);
}
float HologramProfiler::getHologramDiscretizationTime() {
	return computeTimeElapsed(hologramComputed, hologramDiscretized);
}
float HologramProfiler::getBoardUpdateTime() {
	return computeTimeElapsed(hologramDiscretized, boardUpdated);
}
float HologramProfiler::getTotalTime() {
	return computeTimeElapsed(cycleStarted, boardUpdated);
}