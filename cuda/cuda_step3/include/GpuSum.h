/**
 * Functions for sum on GPU
 * @file GpuSum.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#ifndef GPU_SUM_H
#define GPU_SUM_H

#include "FloatType.h"
#include <cuda.h>

/**
 * @brief Compute square of the difference of two vectors: \f$(A-B)^2\f$
 *
 * @param[in] A,B input vector
 * @param[out] C result vector
 * @param[in] size vector size
 */
__global__ void gpu_sqsub(FLOAT_TYPE *A, FLOAT_TYPE *B, FLOAT_TYPE *C, int size);

/**
 * @brief Compute square of the difference of two vectors: \f$(A-B)^2\f$
 *
 * @param[in] A,B input vector
 * @param[out] C result vector
 * @param[in] size vector size
 */
__global__ void gpu_sqsubi(int *A, int *B, FLOAT_TYPE *C, int size);

/**
 * @brief Sum of a vector
 *
 * @param[in] A input vector
 * @param[out] B sum of vector (in the first element)
 * @param[in] size vector size
 */
__global__ void gpu_sum(FLOAT_TYPE *A, FLOAT_TYPE *B, int size);

/**
 * @brief Sum of a vector (for 256 threads only)
 * @note faster than #gpu_sum but works only with 256 threads
 *
 * @param[in] A input vector
 * @param[out] B sum of vector (in the first element)
 * @param[in] size vector size
 */
__global__ void gpu_sum256(FLOAT_TYPE *A, FLOAT_TYPE *B, int size);

/**
 * @brief Conditional vector copy
 * @deprecated use #gpu_cond_copy_mask
 *
 * @param[in]  A input vector
 * @param[out] B output vector
 * @param[in]  cond array with the conditional values
 * @param[in]  value value to compare the conditionals
 * @param[in]  size vector size
 */
__global__ void gpu_cond_copy(FLOAT_TYPE *A, FLOAT_TYPE *B, int *cond, int value, int size);

/**
 * @brief Conditional vector copy for BC bitmask
 *
 * @param[in] A input vector
 * @param[out] B output vector
 * @param[in] mask array of conditional values
 * @param[in] value value to compare to
 * @param[in] size vector size
 */
__global__ void gpu_cond_copy_mask(FLOAT_TYPE *A, FLOAT_TYPE *B, int *mask, int value, int size);

/**
 * @brief Host function for GPU vector sum (calls #gpu_sum256)
 *
 * @param C input vector
 * @param D temporary vector
 * @param size vector size
 * @return sum of the vector
 */
__host__ FLOAT_TYPE gpu_sum_h(FLOAT_TYPE *C, FLOAT_TYPE *D, int size);

#endif