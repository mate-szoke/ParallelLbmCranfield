/* The Cell structure is defined here. The functions are responsible for
scalar-vector-matrix allocation and this header includes the multifunctional
"min" and "max" functions too (comparison of two number). */

#include <stdbool.h>  // bool variables
#include <string.h>   // string handling
#include <upc.h>    

#ifndef SHELLFUNCTIONS_H
#define SHELLFUNCTIONS_H


/////////////////////////////////////////////////////////
/////////////////// STRUCT FOR CELLS ////////////////////
/////////////////////////////////////////////////////////

typedef double MyReal;

struct CellProps
{
  int    Fluid;    
  int    Corner;          
  int    ID;             
  int    Boundary;    
  int    BoundaryID;  
  int    ThreadNumber;         
  int    BC_ID[9];
  int    StreamLattice[9];

  MyReal CoordX;   
  MyReal CoordY;
  MyReal U;
  MyReal V;
  MyReal Rho;
  MyReal Uo;    
  MyReal Vo;    
  MyReal DragF;
  MyReal LiftF;

  MyReal Q[9];           
  MyReal F[9];
  MyReal Feq[9];
  MyReal METAF[9];
  MyReal Fneighbours[9];  
};


////////////////////////////////////////////////////
///////////////////// DEFINE ///////////////////////
////////////////////////////////////////////////////
#include "BlockSizeDefiner.h"

shared [BLOCKSIZE]  int    Fluid[THREADS*BLOCKSIZE];    
shared [BLOCKSIZE]  int    Corner[THREADS*BLOCKSIZE];    
shared [BLOCKSIZE]  int    ID[THREADS*BLOCKSIZE];    
shared [BLOCKSIZE]  int    Boundary[THREADS*BLOCKSIZE];    
shared [BLOCKSIZE]  int    BoundaryID[THREADS*BLOCKSIZE];    
shared [BLOCKSIZE]  int    ThreadNumber[THREADS*BLOCKSIZE];      

shared [BLOCKSIZE]  int    BC_ID[9][THREADS*BLOCKSIZE];    
shared [BLOCKSIZE]  int    StreamLattice[9][THREADS*BLOCKSIZE];    


shared [BLOCKSIZE]  MyReal CoordX[THREADS*BLOCKSIZE];
shared [BLOCKSIZE]  MyReal CoordY[THREADS*BLOCKSIZE];
shared [BLOCKSIZE]  MyReal U[THREADS*BLOCKSIZE];
shared [BLOCKSIZE]  MyReal V[THREADS*BLOCKSIZE];
shared [BLOCKSIZE]  MyReal Rho[THREADS*BLOCKSIZE];
shared [BLOCKSIZE]  MyReal Uo[THREADS*BLOCKSIZE];
shared [BLOCKSIZE]  MyReal Vo[THREADS*BLOCKSIZE];
shared [BLOCKSIZE]  MyReal DragF[THREADS*BLOCKSIZE];
shared [BLOCKSIZE]  MyReal LiftF[THREADS*BLOCKSIZE];


shared [BLOCKSIZE]  MyReal Q[9][THREADS*BLOCKSIZE];
shared [BLOCKSIZE]  MyReal F[9][THREADS*BLOCKSIZE];
shared [BLOCKSIZE]  MyReal Feq[9][THREADS*BLOCKSIZE];
shared [BLOCKSIZE]  MyReal METAF[9][THREADS*BLOCKSIZE];
shared [BLOCKSIZE]  MyReal Fneighbours[9][THREADS*BLOCKSIZE]; 

////////////////////////////////////////////////////
////////////// SHARED DECLARATIONS /////////////////
////////////////////////////////////////////////////

//shared [BLOCKSIZE] struct CellProps *Cells;
shared [BLOCKSIZE] struct CellProps Cells[THREADS*BLOCKSIZE];
//Cells = (shared [BLOCKSIZE] struct CellProps*)upc_all_alloc(THREADS, BLOCKSIZE*sizeof(struct CellProps));

shared [] int    *NumNodes;       // This will store the number of lines of the read files
shared [] int    *NumConn;        // This will store the number of lines of the read files
shared [] int    *n;                 
shared [] int    *m;              // number of nodes in the x and y directions
shared [] int    *NumInletNodes;  // number of inlet nodes
shared [] float  *Delta;          // grid spacing
shared [] double *MaxInletCoordY; // maximum inlet coordinate in y
shared [] double *MinInletCoordY; // minimum inlet coordinate in y

shared int loop;


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
