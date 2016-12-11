/*The first function reads the main solver parameters. The further ones
read the two input files of the mesh and compute some basic parameter from them.*/

#ifndef FILESREADING_H
#define FILESREADING_H

void ReadNodLines(char* NodeDataFile, int* NumOfLines);

void ReadNodes(char* NodeDataFile, int* NumOfLines, int *Nodes0, int *Nodes1, float *Nodes2, float *Nodes3, int *Nodes4);

void ReadBCconLines(char* BCconnectorDataFile, int* NumOfLines);

void ReadBCconn(char* BCconnectorDataFile, int* NumOfLines, int* BCconn0,
                int* BCconn1, int* BCconn2, int* BCconn3, float* BCconn4,
                float* BCconn5, int* BCconn6);

void CompDataNode(float* Delta, int* m,  int* n, int *Nodes0, int *Nodes1,
                  float *Nodes2, float *Nodes3, int *Nodes4,  int* NumNodes);

void CompDataConn(int* NumInletNodes, float* MaxInletCoordY,
	float* MinInletCoordY, int* BCconn0, int* BCconn1, int* BCconn2,
  int* BCconn3, float* BCconn4, float* BCconn5, int* BCconn6, int* NumConn, float* Delta);

void ReadIniData(char* IniFileName, float* Uavg, float* Vavg, float* rho_ini,
                 float* Viscosity, int* InletProfile, int* CollisionModel,
                 int* CurvedBoundaries, int* OutletProfile, int* Iterations,
                 int* AutosaveEvery, int* AutosaveAfter, int* PostprocProg,
                 int* CalculateDragLift);

#endif
