/* This function is the SOLVER itself */

#ifndef ITERATE_H
#define ITERATE_H

#include <upc.h>

#include "ShellFunctions.h" // convenience

void CollisionStep(double* w, int* cx, int* cy, int* opp,
                   double Omega, double OmegaA, float **tm,
                   float **stmiv, int CollisionModel, int lm, int ln,
                   int MyAffinity[NN*MM]);

void StreamingStep(int* c, int lm, int ln, int MyAffinity[NN*MM]);

void HandleBoundariesStep(int* cx, int* cy, int* opp, int OutletProfile, int CurvedBoundaries, int lm, int ln, int MyAffinity[NN*MM]);

void UpdateMacroscopicStep(int* cx, int* cy, int CalculateDragLift, int lm, int ln, int MyAffinity[NN*MM]);

//void CalculateDragLiftForcesStep(struct CellProps *Cells, shared []  int* m, shared [] int* n, int CalculateDragLift);

void Iteration(char* NodeDataFile, char* BCconnectorDataFile,
               float Uavg, float Vavg, float rho_ini, float Viscosity,
               int InletProfile,  int CollisionModel, int CurvedBoundaries,
               int OutletProfile, int Iterations, int AutosaveAfter,
               int AutosaveEvery, int postproc_prog, int CalculateDragLift);


#endif
