/*This function compute the residuals.*/

#ifndef ComputeResiduals_H
#define ComputeResiduals_H

void ComputeResiduals(int* BoundaryID, float* F, float* METAF, float* DragF, float* LiftF, float* Residuals, int* m, int* n, int ComputeDragLift);

#endif
