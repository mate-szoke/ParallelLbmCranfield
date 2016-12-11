/* These functionc are responsible for the different boundary conditions.*/

#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include "include/ShellFunctions.h" // convenience

void InletBC(struct CellProps *Cells, int j, int i);

void OutletBoundaries(struct CellProps *Cells, int j, int i);

void WallBC(struct CellProps *Cells, int j, int i, int* opp);

void CurvedWallBoundaries(struct CellProps *Cells, int j, int i, int* opp);

#endif
