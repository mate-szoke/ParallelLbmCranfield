/* This function is the SOLVER itself */

#ifndef ITERATE_H
#define ITERATE_H

void Iteration(char* NodeDataFile, char* BCconnectorDataFile,
               float Uavg, float Vavg, float rho_ini, float Viscosity,
               int InletProfile,  int CollisionModel, int CurvedBoundaries,
               int OutletProfile, int Iterations, int AutosaveAfter,
               int AutosaveEvery, int postproc_prog, int CalculateDragLift);

#endif
