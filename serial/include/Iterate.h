/* This function is the SOLVER itself */

#ifndef ITERATE_H
#define ITERATE_H

#include "ShellFunctions.h" // convenience

void CollisionStep(struct CellProps *Cells, int* m, int* n, MyReal* w,
               		 int* cx, int* cy, int* opp, MyReal Omega, MyReal OmegaA,
               		 MyReal **tm, MyReal **stmiv, int CollisionModel);

// void StreamingStep(struct CellProps *Cells, int* m, int* n,  int* cx, int* cy);
void StreamingStep(struct CellProps *Cells, int* m, int* n,  int* c);

// void HandleBoundariesStep(struct CellProps *Cells, int* m, int* n,  int* cx, int* cy, int* opp, int OutletProfile, int CurvedBoundaries);
void HandleBoundariesStep(struct CellProps *Cells, int* cx, int* cy, int* c, int* opp, int OutletProfile, int CurvedBoundaries, int* n, int* m);
void UpdateMacroscopicStep(struct CellProps *Cells, int* m, int* n,  int* cx, int* cy, int CalculateDragLift);                        

void CalculateDragLiftForcesStep(struct CellProps *Cells,  int* m, int* n, int CalculateDragLift);

void Iteration(char* NodeDataFile, char* BCconnectorDataFile,
               float Uavg,        float Vavg,        float rho_ini, float Viscosity,
               int InletProfile,   int CollisionModel, int CurvedBoundaries,
               int OutletProfile,  int Iterations,     int AutosaveAfter,
               int AutosaveEvery,  int postproc_prog,  int CalculateDragLift,
               float ConvergenceCritVeloc,            float ConvergenceCritRho,
               int ReadFormerData);


#endif
