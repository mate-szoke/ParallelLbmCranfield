/**
 * Utility functions for the unittests
 * @file TestUtils.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#ifndef TEST_UTILS
#define TEST_UTILS

#include "FloatType.h"

static int M,N;
static FLOAT_TYPE *AFH, *BFH, *CFH, *DFH, *EFH, *FFH, *GFH, *HFH;
static FLOAT_TYPE *AFD, *BFD, *CFD, *DFD, *EFD, *FFD, *GFD, *HFD;
static int *AIH, *BIH, *CIH, *DIH, *EIH, *FIH;
static int *AID, *BID, *CID, *DID, *EID, *FID;

/**
 * @brief Add two integer vectors on the GPU
 *
 * @param[in,out] array1_d vector (result stored here)
 * @param[in] array2_d vector
 * @param size vector size
 */
__global__ void gpuArrayAddInt(int *array1_d, int *array2_d, int size);

/**
 * @brief Add two floating point vectors on the GPU
 *
 * @param[in,out] array1_d vector (result stored here)
 * @param[in] array2_d vector
 * @param size vector size
 */
__global__ void gpuArrayAddFlt(FLOAT_TYPE *array1_d, FLOAT_TYPE *array2_d, int size);

/**
 * @brief Add two integer vectors on the host
 *
 * @param[in,out] array1_h vector (result stored here)
 * @param[in] array2_h vector
 * @param size vector size
 */
void hostArrayAddInt(int *array1_h, int *array2_h, int size);

/**
 * @brief Add two floating point vectors on the host
 *
 * @param[in,out] array1_h vector (result stored here)
 * @param[in] array2_h vector
 * @param size vector size
 */
void hostArrayAddFlt(FLOAT_TYPE *array1_h, FLOAT_TYPE *array2_h, int size);

/**
 * @brief Summaries integer vector on the host
 *
 * @param array_h vector
 * @param size vector size
 * @return sum of vector
 */
int sumHostInt(int *array_h, int size);

/**
 * @brief Summaries floating point vector on the host
 *
 * @param array_h vector
 * @param size vector size
 * @return sum of vector
 */
FLOAT_TYPE sumHostFlt(FLOAT_TYPE *array_h, int size);

/**
 * @brief Create fluid conditions for the lid driven cavity
 *
 * @param[out] fluid fluid array
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 */
void createLidBcFluid(int *fluid, int m, int n);

/**
 * @brief Create node boundary conditions for the lid driven cavity
 *
 * @param[out] boundary node boundary array
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 */
void createLidBcBoundary(int *boundary, int m, int n);

/**
 * @brief Create lattice boundary conditions for the lid driven cavity
 *
 * @param[out] bcid lattice boundary array
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 */
void createLidBcBcId(int *bcid, int m, int n);

/**
 * @brief Create corner conditions for the lid driven cavity
 *
 * @param[out] corner corner boundary array
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 */
void createLidBcCorner(int *corner, int m, int n);

/**
 * @brief Boundary condition bitmask size for the lid driven cavity
 *
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 * @return number of boundary conditions
 */
int getLidBcMaskSize(int m, int n);

/**
 * @brief Create boundary condition indices for the lid driven cavity
 *
 * @param[out] index boundary condition indices array
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 */
void createLidBcIdx(int *index, int m, int n);

/**
 * @brief Create boundary conditions bitmask for the lid driven cavity (collapsed)
 *
 * @param[out] mask boundary condition bitmask array
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 */
void createLidBcMask(int *mask, int m, int n);


/**
 * @brief Create boundary conditions bitmask for the lid driven cavity (full size)
 *
 * @param[out] mask boundary condition bitmask array
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 */
void createLidBcMaskFull(int *mask, int m, int n);

/**
 * @brief Create boundary IDs for the lid driven cavity
 *
 * @param bid boundary id array
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 */
void createLidBoundaryId(int *bid, int m, int n);

/**
 * @brief Create coordinates for the lid driven cavity
 *
 * @param[out] x,y coordinates
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 */
void fillLidCoordinates(FLOAT_TYPE *x, FLOAT_TYPE *y, int m, int n);

/**
 * @brief Fill distribution function in only one direction
 *
 * @param[out] f distribution fucntion
 * @param[in]  fu value to fill
 * @param[in]  dir direction
 * @param[in]  size array size
 */
void fillFDir(FLOAT_TYPE *f, FLOAT_TYPE fu, int dir, int size);

/**
 * @brief Create boundary conditions bitmask for the channel
 *
 * @param[out] mask boundary condition bitmask array
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 */
void createChannelBcMask(int *mask, int m, int n);

/**
 * @brief Create lattice boundary conditions for the channel
 *
 * @param[out] bcid lattice boundary array
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 */
void createChannelBcBcId(int *bcid, int m, int n);

/**
 * @brief Create corner conditions for the channel
 *
 * @param[out] corner corner boundary array
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 */
void createChannelBcCorner(int *corner, int m, int n);

/**
 * @brief Compare integer arrays by element
 *
 * @param array1_h,array2_h arrays to compare
 * @param size array size
 * @return number of different elements
 */
int compareArraysInt(int *array1_h, int *array2_h, int size);

/**
 * @brief Compare floating point arrays by element
 *
 * @param array1_h,array2_h arrays to compare
 * @param size array size
 * @return number of different elements
 */
int compareArraysFlt(FLOAT_TYPE *array1_h, FLOAT_TYPE *array2_h, int size);

/**
 * @brief Compare GPU int array on CPU
 * @details Compares two array value by value
 *
 * @param a_d,b_d arrays to compare
 * @param size array size
 * @return number of differences
 */
__host__ int compareGpuArrayInt(int *a_d, int *b_d, int size);

/**
 * @brief Compare GPU float array on CPU
 * @details Compares two array value by value
 *
 * @param a_d,b_d arrays to compare
 * @param size array size
 * @return number of differences
 */
__host__ int compareGpuArrayFlt(FLOAT_TYPE *a_d, FLOAT_TYPE *b_d, int size);

/**
 * @brief Compare two vectors on the GPU
 * @details Compare two vectors on the GPU by summing up the squared differences
 * (Frobenius norm: \f$\sqrt{\sum_{i=1}^{m}\sum_{j=1}^{n}|a_{ij}-b_{ij}|^2}\f$)
 *
 * @param F1_d,F2_d vectors to compare
 * @param temp9a_d,temp9b_d temporary vectors
 * @param size vector size
 * @return Frobenius norm of the difference of the vectors
 */
__host__ FLOAT_TYPE compareGpuArraySumFlt(FLOAT_TYPE *F1_d, FLOAT_TYPE *F2_d, FLOAT_TYPE *temp9a_d, FLOAT_TYPE *temp9b_d, int size);

/**
 * @brief Compare two vectors on the GPU
 * @details Compare two vectors on the GPU by summing up the squared differences
 * (Frobenius norm: \f$\sqrt{\sum_{i=1}^{m}\sum_{j=1}^{n}|a_{ij}-b_{ij}|^2}\f$)
 *
 * @param F1_d,F2_d vectors to compare
 * @param temp9a_d,temp9b_d temporary vectors
 * @param size vector size
 * @return Frobenius norm of the difference of the vectors
 */
__host__ FLOAT_TYPE compareGpuArraySumInt(int *F1_d, int *F2_d, FLOAT_TYPE *temp9a_d, FLOAT_TYPE *temp9b_d, int size);

#define MAX_LEN 80 ///< max length for banner

/**
 * @brief Print a formated banner for the test case
 *
 * @param testname test name
 */
void printBanner(const char *testname);

#endif