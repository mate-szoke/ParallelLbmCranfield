/*This function compute the residuals.*/

#ifndef ComputeResiduals_H
#define ComputeResiduals_H

void ComputeResiduals(struct CellProps_var *Cells_var_h, struct CellProps_const *Cells_const_h,
	                  float* Residuals, int* m, int* n, float* sumVel0, float* sumVel1,
	                  float* sumRho0, float* sumRho1, int ComputeDragLift);

#endif
