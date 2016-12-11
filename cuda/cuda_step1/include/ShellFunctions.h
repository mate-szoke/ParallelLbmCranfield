/* The Cell structure is defined here. The functions are responsible for
scalara-vector-matrix allocation and this header includes the multifunctional
"min" and "max" functions too (comparison of two number). */

#include <stdbool.h>  // bool variables
#include <string.h>   // string handling

#ifndef SHELLFUNCTIONS_H
#define SHELLFUNCTIONS_H


/////////////////// STRUCT FOR CELLS ////////////////////
struct CellProps_const
{
  int    Fluid;    
  int    Corner;          
  int    ID;             
  int    Boundary;    
  int    BoundaryID;  
  float  CoordX;   
  float  CoordY;
};

struct CellProps_const_9d
{
  int    StreamLattice;
  int    BC_ID;         
  float  Q;           
};

struct CellProps_var
{
  float U;
  float V;
  float Rho;
  float Uo;    
  float Vo;    
  float DragF;
  float LiftF;
};


struct CellProps_var_9d
{
  float F;
  float Feq;
  float METAF;
  float Fneighbours;
  float fmom;
  float fmeq;
  float sumb;
};

int    *Create1DArrayInt(int length);
float  *Create1DArrayFloat(int length);
double *Create1DArrayDouble(int length);

int    **Create2DArrayInt(int width, int height);
float  **Create2DArrayFloat(int width, int height);
double **Create2DArrayDouble(int width, int height);

int    ***Create3DArrayInt(int width, int height, int depth);
float  ***Create3DArrayFloat(int width, int height, int depth);
double ***Create3DArrayDouble(int width, int height, int depth);
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
