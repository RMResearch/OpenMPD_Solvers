#pragma once
#include <stdio.h>

namespace FileNameGenerator {
	void generateTrapDataFileName(char outputBuffer[512], bool corrected, float posX, float posY, float posZ, float rotYDegrees) {
		sprintf(outputBuffer, "./ForceModelData/%s/TrapStatsSummary.csv", (corrected ? "Corrected" : "Uncorrected"));
	}
	void generateTrapForcesFileName(char outputBuffer[512], bool corrected, float posX, float posY, float posZ, float rotYDegrees) {
		sprintf(outputBuffer, "./ForceModelData/%s/TrapForcesAt(%.2f,%.2f,%.2f)_rot%.1f.csv", (corrected ? "Corrected" : "Uncorrected"), posX, posY, posZ, rotYDegrees);
	}
	void generateTrapImageFileName(char outputBuffer[512], bool corrected, float posX, float posY, float posZ, float rotYDegrees, char axis) {
		sprintf(outputBuffer, "./ForceModelData/%s/TrapForcesAt(%.2f,%.2f,%.2f)_rot%.1f_%c.bmp", (corrected ? "Corrected" : "Uncorrected"), posX, posY, posZ, rotYDegrees, axis);
	}
};


