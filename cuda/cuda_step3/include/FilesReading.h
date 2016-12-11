#ifndef FILESREADING_H
#define FILESREADING_H

#include "FloatType.h"
#include "Arguments.h"

/**
 * Get the number of lines in a file
 * @param[in] filename filename
 * @return number of lines
 */
int getNumberOfLines(const char *filename);

/**
 * Read init file (usually SetUpData.ini)
 * @param[in] filename filename
 * @param[out] args input parameters
 */
void readInitFile(const char* filename, Arguments *args);

/**
 * Read node file (usually D2node.dat)
 * Arrays supplied are allocated inside, freeing these arrays is the caller's responsibility
 *
 * @param[in] filename filename
 * @param[out] ni,nj node ID (x,y)
 * @param[out] nx,ny node coordinate (x,y)
 * @param[out] nf node type (0: solid, 1:fluid)
 * @return number of lines read
 */
int readNodeFile(const char *filename, int **ni, int **nj, FLOAT_TYPE **nx, FLOAT_TYPE **ny, int **nf);

/**
 * Read boundary conditions file (usually BCconnectors.dat)
 * Arrays supplied are allocated inside, freeing these arrays is the caller's responsibility
 *
 * @verbatim
   6   2   5
     \ | /
   3 - 0 - 1
     / | \
   7   4   8 @endverbatim
 *
 * @param[in] filename filename
 * @param[out] ni,nj node ID (x,y)
 * @param[out] dir direction
 * @param[out] bc BC type (1: wall, 2: inlet, 3: outlet)
 * @param[out] bcx,bcy BC coordinate (x,y)
 * @param[out] id boundary ID (for drag/lift computation)
 * @return number of lines read
 */
int readConnFile(const char *filename, int **ni, int **nj, int **dir, int **bc,
                 FLOAT_TYPE **bcx, FLOAT_TYPE **bcy, int **id);

/**
 * Read result file (usually FinalData.csv)
 * @param[in] filename filename
 * @param[out] results coordinates (x,y), velocity (u,v), vel_mag, rho, pressure
 * @param[out] fluid node type (0: solid, 1:fluid)
 * @return number of lines read
 */
int readResultFile(const char *filename, FLOAT_TYPE ***results, int **fluid);

/**
 * @brief Compare result files
 * @note Use it to compare FinalData.csv with previous result.
 *
 * @param f1,f2 files to compare
 * @return 0 for identical files and 1 for different files
 */
__host__ int compareFiles(const char* f1, const char* f2);

/**
 * Gets the last value of the array
 * Use it with nodeIdx to get number rows, and with nodeIdy to get number of columns
 * @param arr input array
 * @param n array size
 * @return last element plus one
 */
int getLastValue(int *arr, int n);

/**
 * Returns the grid spacing of the mesh
 * @param ni,nj node ID (x,y)
 * @param nx node coordinate x
 * @param n number of nodes
 * @return grid spacing
 */
FLOAT_TYPE getGridSpacing(int *ni, int *nj, FLOAT_TYPE *nx, int n);

/**
 * Get number of inlet nodes
 * @param bc BC type (1: wall, 2: inlet, 3: outlet)
 * @param dir direction
 * @param n number of BC nodes
 * @return number of inlet nodes
 */
int getNumInletNodes(int *bc, int *dir, int n);

/**
 * Get the maximum coordinate of the inlet in the Y direction
 * @param bc BC type (1: wall, 2: inlet, 3: outlet)
 * @param dir direction
 * @param bcy BC coordinate y
 * @param delta grid spacing
 * @param n number of BC nodes
 * @return maximum coordinate of the inlet in the Y direction
 */
FLOAT_TYPE getMaxInletCoordY(int *bc, int *dir, FLOAT_TYPE *bcy, FLOAT_TYPE delta, int n);

/**
 * Get the minimum coordinate of the inlet in the Y direction
 * @param bc BC type (1: wall, 2: inlet, 3: outlet)
 * @param dir direction
 * @param bcy BC coordinate y
 * @param delta grid spacing
 * @param n number of BC nodes
 * @return minimum coordinate of the inlet in the Y direction
 */
FLOAT_TYPE getMinInletCoordY(int *bc, int *dir, FLOAT_TYPE *bcy, FLOAT_TYPE delta, int n);

/**
 * @brief Compute values from node data
 * @deprecated use #getLastValue and #getGridSpacing instead
 *
 * @param[out] Delta grid spacing
 * @param[out] m number of columns
 * @param[out] n number of rows
 * @param[in] Nodes0 node id x
 * @param[in] Nodes1 node id y
 * @param[in] Nodes2 node coordinate x
 * @param[in] Nodes3 node coordinate y
 * @param[in] Nodes4 node type
 * @param[in] NumNodes number of nodes
 */
void CompDataNode(FLOAT_TYPE* Delta, int* m,  int* n, int *Nodes0, int *Nodes1,
                  FLOAT_TYPE *Nodes2, FLOAT_TYPE *Nodes3, int *Nodes4,  int* NumNodes);

/**
 * @brief Compute values from BC data
 * @deprecated use #getNumInletNodes, #getMaxInletCoordY and #getMinInletCoordY instead
 *
 * @param[out] NumInletNodes number of inlet nodes
 * @param[out] MaxInletCoordY maximum inlet coordinate y
 * @param[out] MinInletCoordY minimum inlet coordinate y
 * @param[in] BCconn0 node id x
 * @param[in] BCconn1 node id y
 * @param[in] BCconn2 lattice id (direction)
 * @param[in] BCconn3 boundary type
 * @param[in] BCconn4 BC coordinate x
 * @param[in] BCconn5 BC coordinate y
 * @param[in] BCconn6 BC ID for drag/lift
 * @param[in] NumConn number of BC
 * @param[in] Delta grid spacing
 */
void CompDataConn(int* NumInletNodes, FLOAT_TYPE* MaxInletCoordY,
	FLOAT_TYPE* MinInletCoordY, int* BCconn0, int* BCconn1, int* BCconn2,
  int* BCconn3, FLOAT_TYPE* BCconn4, FLOAT_TYPE* BCconn5, int* BCconn6, int* NumConn, FLOAT_TYPE* Delta);

#endif
