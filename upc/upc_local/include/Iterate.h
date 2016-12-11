/* This function is the SOLVER itself */

#ifndef ITERATE_H
#define ITERATE_H

#include "ShellFunctions.h" // convenience

void putCellsToShared(struct CellProps *Cells);

void getSharedToCells(struct CellProps *Cells);

void putCellsToWCells(struct CellProps *Cells);

void CollisionStep(struct CellProps *Cells, double* w, int* cx, int* cy, int* opp,
                   double Omega, double OmegaA, double **tm,
                   double **stmiv, int CollisionModel);

void StreamingStep(struct CellProps *Cells, int* c);

void HandleBoundariesStep(struct CellProps *Cells, int* cx, int* cy, int* c, int* opp, int OutletProfile, int CurvedBoundaries);

void UpdateMacroscopicStep(struct CellProps *Cells, int* cx, int* cy, int CalculateDragLift);

void CalculateDragLiftForcesStep(struct CellProps *Cells, shared []  int* m, shared [] int* n, int CalculateDragLift);

void Iteration(char* NodeDataFile, char* BCconnectorDataFile,
               float Uavg,         float Vavg,         float rho_ini, float Viscosity,
               int InletProfile,   int CollisionModel, int CurvedBoundaries,
               int OutletProfile,  int Iterations,     int AutosaveAfter,
               int AutosaveEvery,  int postproc_prog,  int CalculateDragLift,
               float ConvergenceCritVeloc, float ConvergenceCritRho);


#endif
