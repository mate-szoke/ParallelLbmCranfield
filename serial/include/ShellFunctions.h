/* The Cell structure is defined here. The functions are responsible for
scalar-vector-matrix allocation and this header incluzdes the multifunctional
"min" and "max" functions too (comparison of two number). */

#include <stdbool.h>  // bool variables
#include <string.h>   // string handling

#ifndef SHELLFUNCTIONS_H
#define SHELLFUNCTIONS_H

typedef double MyReal;

/////////////////// STRUCT FOR CELLS ////////////////////
struct CellProps
{
  int    Fluid;    
  int    Corner;          
  int    StreamLattice[9];
  int    ID;             
  int    Boundary;    
  int    BoundaryID;  
  int    BC_ID[9];         
  MyReal Q[9];           
  MyReal CoordX;   
  MyReal CoordY;
  MyReal U;
  MyReal V;
  MyReal Rho;
  MyReal Uo;    
  MyReal Vo;    
  MyReal DragF;
  MyReal LiftF;
  MyReal F[9];
  MyReal Feq[9];
  MyReal METAF[9];
  MyReal Fneighbours[9];	
};

int    *Create1DArrayInt(int length);
MyReal *Create1DArrayMyReal(int length);

int    **Create2DArrayInt(int width, int height);
MyReal **Create2DArrayMyReal(int width, int height);

int    ***Create3DArrayInt(int width, int height, int depth);
MyReal ***Create3DArrayMyReal(int width, int height, int depth);
bool   ***Create3DArrayBool(int width, int height, int depth);

void CreateDirectory(char* MainWorkDir);

void StringAddition(char* first, char* second, char* result);


#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })


#endif
