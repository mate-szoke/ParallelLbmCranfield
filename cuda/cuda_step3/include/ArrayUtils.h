/**
 * Memory allocation wrappers for host and gpu arrays
 * @file ArrayUtils.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */

#ifndef ARRAY_UTILS_H
#define ARRAY_UTILS_H

#include <stdbool.h>  // bool variables
#include "FloatType.h"

#define THREADS 256 ///< number of threads used for kernels

/**
 * wraps sizeof(int)
 * @param a number of ints
 */
#define SIZEINT(a) ((a)*sizeof(int))
/**
 * wraps sizeof(FLOAT_TYPE)
 * @param a number of floats/doubles
 */
#define SIZEFLT(a) ((a)*sizeof(FLOAT_TYPE))

/**
 * Get a uniform random number
 * @param seed seed for the random number (initialise with time(NULL))
 * @return uniform random number
 */
FLOAT_TYPE getRandom(unsigned long *seed);

/// Array creation options
typedef enum
{
    ARRAY_NONE, ///< Simple allocation
    ARRAY_ZERO, ///< Allocation and memset to zero
    ARRAY_FILL, ///< Allocation and fill with given value
    ARRAY_COPY, ///< Allocation and copy from given array
    ARRAY_CPYD, ///< Allocation and copy from given GPU array
    ARRAY_RAND  ///< Allocation and fill with random values
} ArrayOption;

/**
 * Create 1D int array on the host
 * @param length length of the array
 * @param op options (default: none)
 * @param fill value to fill array (default: 0)
 * @param copy array to copy into the new array (default: NULL)
 * @return new array
 */
int *createHostArrayInt(int length, ArrayOption op=ARRAY_NONE, int fill=0, int *copy=NULL);

/**
 * Create 1D int array on the device
 * @param length length of the array
 * @param op options (default: none)
 * @param fill value to fill array (default: 0)
 * @param copy host array to copy into the new array (default: NULL)
 * @return new array
 */
int *createGpuArrayInt(int length, ArrayOption op=ARRAY_NONE, int fill=0, int *copy=NULL);

/**
 * Create 1D float/double array on the device
 * @param length length of the array
 * @param op options (default: none)
 * @param fill value to fill array (default: 0)
 * @param copy array to copy into the new array (default: NULL)
 * @return new array
 */
FLOAT_TYPE *createHostArrayFlt(int length, ArrayOption op=ARRAY_NONE, FLOAT_TYPE fill=0, FLOAT_TYPE *copy=NULL);

/**
 * Create 1D float/double array on the device
 * @param length length of the array
 * @param op options (default: none)
 * @param fill value to fill array (default: 0)
 * @param copy host array to copy into the new array (default: NULL)
 * @return new array
 */
FLOAT_TYPE *createGpuArrayFlt(int length, ArrayOption op=ARRAY_NONE, FLOAT_TYPE fill=0, FLOAT_TYPE *copy=NULL);

/**
 *  @brief Allocate 2D integer array
 *  @note not used
 *
 *  @param [in] width array width
 *  @param [in] height array height
 *  @return pointer to the allocated array
 */
int    **create2DHostArrayInt(int width, int height);

/**
 *  @brief Allocate 2D floating point array
 *  @note not used
 *
 *  @param [in] width array width
 *  @param [in] height array height
 *  @return pointer to the allocated array
 */
FLOAT_TYPE  **create2DHostArrayFlt(int width, int height);

/**
 *  @brief Allocate 3D integer array
 *  @note not used
 *
 *  @param [in] width array width
 *  @param [in] height array height
 *  @param [in] depth array depth
 *  @return pointer to the allocated array
 */
int    ***create3DHostArrayInt(int width, int height, int depth);

/**
 *  @brief Allocate 3D floating point array
 *  @note not used
 *
 *  @param [in] width array width
 *  @param [in] height array height
 *  @param [in] depth array depth
 *  @return pointer to the allocated array
 */
FLOAT_TYPE  ***create3DHostArrayFlt(int width, int height, int depth);

/**
 *  @brief Allocate 3D boolean array
 *  @note not used
 *
 *  @param [in] width array width
 *  @param [in] height array height
 *  @param [in] depth array depth
 *  @return pointer to the allocated array
 */
bool   ***create3DHostArrayBool(int width, int height, int depth);

/**
 * Fill int array with given value
 * @param array_h array on device
 * @param fill value to fill
 * @param size size of the array
 */
void hostArrayFillInt(int *array_h, int fill, int size);

/**
 * GPU kernel to fill int array with given value
 * @param array_d array on device
 * @param size size of the array
 * @warning fill value should be set into variable fill_id in constant memory
 */
__global__ void gpuArrayFillInt(int *array_d, int size);

/**
 * Fill float/double array with given value
 * @param array_h array on device
 * @param fill value to fill
 * @param size size of the array
 */
void hostArrayFillFlt(FLOAT_TYPE *array_h, FLOAT_TYPE fill, int size);

/**
 * GPU kernel to fill float/double array with given value
 * @param array_d array on device
 * @param size size of the array
 * @warning fill value should be set into variable fill_fd in constant memory
 */
__global__ void gpuArrayFillFlt(FLOAT_TYPE *array_d, int size);

/**
 * Fill float/double array with random value
 * @param array_h array on device
 * @param size size of the array
 * @param r range of random number
 */
void hostArrayFillRandom(FLOAT_TYPE *array_h, int size, FLOAT_TYPE r=1.0);

/**
 * GPU kernel to fill float/double array with random value
 * @param array_d array on device
 * @param seed seed for random generator (spiked with global id)
 * @param size size of the array
 */
 __global__ void gpuArrayFillRandom(FLOAT_TYPE *array_d, unsigned long seed, int size);

 /**
  * Free all host arrays
  * @param as array of arrays
  * @param n number of arrays
  */
void freeAllHost(void **as, int n);

 /**
  * Free all GPU arrays
  * @param as array of arrays
  * @param n number of arrays
  */
void freeAllGpu(void **as, int n);

#endif