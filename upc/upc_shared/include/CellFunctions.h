/* These functionc are responsible for the initialization and includes the collision model.
The microscopic and macroscopic variables are refreshed based on the last two functions.*/

#ifndef CELLFUNCTIONS_H
#define CELLFUNCTIONS_H

void D2Q9Vars(double* w, int* cx, int* cy, int* opp, int* c);

void MRTInitializer(float** tm, float** stmiv, double Omega);

void CellIni(float **Nod,             
             float **Con,
             int A,                   
             int B, // DiscreteX and DiscreteY 
             float Uavg,
             float Vavg,              
             int InletProfile,
             int CollisionModel,      
             int* opp,
             float rho_ini,
             int lm,
             int ln) ;

void BGKW(int i, double* w, int* cx, int* cy, double Omega);

void TRT(int i, double* w, int* cx, int* cy, int* opp, double Omega, double OmegaA);

void MRT(int i, float** tm, float** stmiv);

void UpdateF(int i);

void UpdateMacroscopic(int i, int* cx, int* cy, int CalculateDragLift);

//void CalculateDragLiftForces(struct CellProps *Cells, int j, int i, int CalculateDragLift, shared [] int* n, shared [] int* m);

#endif
