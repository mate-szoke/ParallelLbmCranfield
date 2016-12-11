/*The first function reads the main solver parameters. The further ones
read the two input files of the mesh and compute some basic parameter from them.*/

#ifndef FILESREADING_H
#define FILESREADING_H

float **ReadNodes(char* NodeDataFile);

float **ReadBCconn(char* BCconnectorDataFile);

void CompDataNode(float **Nodes);

void CompDataConn(float** BCconn);

void ReadIniData(char* IniFileName, float* Uavg, float* Vavg, float* rho_ini,
                 float* Viscosity, int* InletProfile, int* CollisionModel,
                 int* CurvedBoundaries, int* OutletProfile, int* Iterations,
                 int* AutosaveEvery, int* AutosaveAfter, int* PostprocProg,
                 int* CalculateDragLift, float* ConvergenceCritVeloc,
                 float* ConvergenceCritRho);

#endif
