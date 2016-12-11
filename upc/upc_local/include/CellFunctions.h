/* These functionc are responsible for the initialization and includes the collision model.
The microscopic and macroscopic variables are refreshed based on the last two functions.*/

#ifndef CELLFUNCTIONS_H
#define CELLFUNCTIONS_H

void D2Q9Vars(double* w, int* cx, int* cy, int* opp, int* c);

void MRTInitializer(double** tm, double** stmiv, double Omega);

void CellIni(struct CellProps *Cells,
             float **Nod,             
             float **Con,
             int A, // i              
             int B, // j
             float Uavg,
             float Vavg,              
             int InletProfile,
             int CollisionModel,      
             int* opp,
             float rho_ini,
             int shift);

void BGKW(struct CellProps *Cells, int i, double* w, int* cx, int* cy, double Omega);

void TRT(struct CellProps *Cells, int i, double* w, int* cx, int* cy, int* opp, double Omega, double OmegaA);

void MRT(struct CellProps *Cells, int i, double** tm, double** stmiv);

void UpdateF(struct CellProps *Cells, int i);

void UpdateMacroscopic(struct CellProps *Cells, int i, int* cx, int* cy, int CalculateDragLift);

void CalculateDragLiftForces(struct CellProps *Cells, int j, int i, int CalculateDragLift, shared [] int* n, shared [] int* m);

#endif
