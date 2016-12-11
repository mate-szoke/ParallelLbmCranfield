/*This function compute the residuals.*/

#ifndef ComputeResiduals_H
#define ComputeResiduals_H

#include "include/ShellFunctions.h"

void ComputeResiduals(struct CellProps *Cells, double* Residuals,
                    double* sumVel0, double* sumVel1, double* sumRho0,
                    double* sumRho1, int ComputeDragLift, int* iter, int* Iterations);

//void ComputeResiduals(struct CellProps *Cells, double* Residuals, double* sumVel0, double* sumVel1, double* sumRho0, double* sumRho1, int ComputeDragLift, int* iter, int* Iterations);

#endif
