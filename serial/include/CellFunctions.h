/* These functionc are responsible for the initialization and includes the collision model.
The microscopic and macroscopic variables are refreshed based on the last two functions.*/

#ifndef CELLFUNCTIONS_H
#define CELLFUNCTIONS_H

//void D2Q9Vars(MyReal* w, int* cx, int* cy, int* opp );
void D2Q9Vars(MyReal* w, int* cx, int* cy, int* opp, int* c, int* n);

void MRTInitializer(MyReal** tm, MyReal** stmiv, MyReal Omega);

void CellIni(float **Nod, float **Con, int* NumN, int* NumC, int A, int B,  
             MyReal* MinInletCoordY, MyReal* MaxInletCoordY, MyReal* Delta,
             float Uavg, float Vavg, int InletProfile, int CollisionModel,
             int* opp, struct CellProps *Cells, float rho_ini, int* n);

void BGKW(struct CellProps *Cells, int i, MyReal* w, int* cx, int* cy, MyReal Omega);

void TRT(struct CellProps *Cells, int i, MyReal* w, int* cx, int* cy, int* opp, MyReal Omega, MyReal OmegaA);

void MRT(struct CellProps *Cells, int i, MyReal** tm, MyReal** stmiv);

void UpdateF(struct CellProps *Cells, int i);


void UpdateMacroscopic(struct CellProps *Cells, int j, int i, int* cx, int* cy, int CalculateDragLift);

void CalculateDragLiftForces(struct CellProps *Cells,  int i,int j, int CalculateDragLift);

#endif
