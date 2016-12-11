/*This function is responsible for the output field variables.*/

#ifndef FILESWRITING_H
#define FILESWRITING_H

void WriteResults(char* OutputFile, float* CoordX, float* CoordY, float* U, float* V, float* Rho, int* Fluid, int* n, int* m, int* postproc_prog);

#endif
