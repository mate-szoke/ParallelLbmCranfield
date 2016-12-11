/* These functionc are responsible for the different boundary conditions.*/

#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

void InletBC(int j, int i, int lm, int ln);

void OutletBoundaries(int j, int i, int lm, int ln);

void WallBC(int j, int i, int* opp, int lm, int ln);

void CurvedWallBoundaries(int j, int i, int* opp, int lm, int ln);

#endif
