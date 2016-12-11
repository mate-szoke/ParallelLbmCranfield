/* These functionc are responsible for the different boundary conditions.*/

#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

//void InletBC(struct CellProps *Cells, int j, int i);
void InletBC(struct CellProps *Cells, int j, int i, int* n, int* m);

//void OutletBoundaries(struct CellProps *Cells, int j, int i);
void OutletBoundaries(struct CellProps *Cells, int j, int i, int* n, int* m);

//void WallBC(struct CellProps *Cells, int j, int i, int* opp);
void WallBC(struct CellProps *Cells, int j, int i, int* opp, int* n, int* m);

//void CurvedWallBoundaries(struct CellProps *Cells, int j, int i, int* opp);
void CurvedWallBoundaries(struct CellProps *Cells, int j, int i, int* opp, int* n, int* m);

#endif
