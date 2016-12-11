/* The Cell structure is defined here. The functions are responsible for
scalar-vector-matrix allocation and this header includes the multifunctional
"min" and "max" functions too (comparison of two number). */

#include <stdbool.h>  // bool variables
#include <string.h>   // string handling

#ifndef SHELLFUNCTIONS_H
#define SHELLFUNCTIONS_H

////////////////////////////////////////////////////
///////////////////// DEFINE ///////////////////////
////////////////////////////////////////////////////

//#define UPC_MAX_BLOCK_SIZE 145000000

#include "BlockSizeDefiner.h"

//#define BLOCKSIZE 802*((int)(83/THREADS+1))
//#define NN 802
//#define MM 83

////////////////////////////////////////////////////
////////////// SHARED DECLARATIONS /////////////////
////////////////////////////////////////////////////

//shared [BLOCKSIZE+2*NN] struct CellProps  *SCells;
shared [BLOCKSIZE]      struct CellProps  *WCells; // Writing Cells: cells to write data
shared [2*NN]           struct CellProps  *BCells; // Boundary cells
shared [4]              double sResiduals[4*THREADS]; // variable to store residuals

shared [] int    *NumNodes;       // This will store the number of lines of the read files
shared [] int    *NumConn;        // This will store the number of lines of the read files
shared [] int    *n;              // number of nodes in the x direction   
shared [] int    *m;              // number of nodes in the y direction
shared [] int    *NumInletNodes;  // number of inlet nodes
shared [] double *Delta;          // grid spacing
shared [] double *MaxInletCoordY; // maximum inlet coordinate in y
shared [] double *MinInletCoordY; // minimum inlet coordinate in y

shared int loop;

typedef double MyReal;

/////////////////////////////////////////////////////////
/////////////////// STRUCT FOR CELLS ////////////////////
/////////////////////////////////////////////////////////

struct CellProps
{
  int    Fluid;    
  int    Corner;          
  int    StreamLattice[9];
  int    ID;             
  int    Boundary;    
  int    BoundaryID;  
  int    BC_ID[9];
  int    ThreadNumber;         
  MyReal  Q[9];           
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
