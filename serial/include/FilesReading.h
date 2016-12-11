/*The first function reads the main solver parameters. The further ones
read the two input files of the mesh and compute some basic parameter from them.*/

#ifndef FILESREADING_H
#define FILESREADING_H

float **ReadNodes(char* NodeDataFile, int* NumOfLines);

float **ReadBCconn(char* BCconnectorDataFile, int* NumOfLines);

void CompDataNode(MyReal* Delta, int* m, int* n, float **Nodes, int* NumNodes);

void CompDataConn(int* NumInletNodes, MyReal* MaxInletCoordY, MyReal* MinInletCoordY,
				  float** BCconn, int* NumConn, MyReal* Delta);

void ReadIniData(char* IniFileName, float* Uavg, float* Vavg, float* rho_ini,
                 float* Viscosity, int* InletProfile, int* CollisionModel,
                 int* CurvedBoundaries, int* OutletProfile, int* Iterations,
                 int* AutosaveEvery, int* AutosaveAfter, int* PostprocProg,
                 int* CalculateDragLift, float* ConvergenceCritVeloc,
                 float* ConvergenceCritRho);

#endif
