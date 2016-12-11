/**
 * Functions for residual computations and result comparing
 * @file ComputeResiduals.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */

#ifndef ComputeResiduals_H
#define ComputeResiduals_H

#include "FloatType.h"

/**
 * @brief This function compute the residuals
 * @warning Exits the whole program if the iteration diverges
 *
 * @param[in] boundaryId boundary types
 * @param[in] f distribution function
 * @param[in] fColl distribution function (previous step)
 * @param[in] drag drag
 * @param[in] lift lift
 * @param[out] residuals array for residuals and lift/drag
 * @param[in] m number of columns
 * @param[in] n number of rows
 * @param[in] computeDragLift boundary ID for drag/lift computation
 * @deprecated use #GpuComputeResidMask instead
 */
void ComputeResiduals(int* boundaryId, FLOAT_TYPE* f, FLOAT_TYPE* fColl, FLOAT_TYPE* drag,
                      FLOAT_TYPE* lift, FLOAT_TYPE* residuals, int* m, int* n, int computeDragLift);

/**
 * @brief Compute the residuals on the GPU
 * @warning Exits the whole program if the iteration diverges
 * @deprecated use #GpuComputeResidMask
 *
 * @param[in] boundaryId_d boundary types
 * @param[in] f_d distribution function
 * @param[in] fColl_d distribution function (previous step)
 * @param[in] drag_d drag
 * @param[in] lift_d lift
 * @param[in] temp9a_d,temp9b_d temporary array for summation (size: 9xMxN)
 * @param[in] tempA_d, tempB_d temporary array for summation (size: MxN)
 * @param[out] residuals array for residuals and lift/drag
 * @param[in] m number of columns
 * @param[in] n number of rows
 * @param[in] computeDragLift boundary ID for drag/lift computation
 */
__host__ void GpuComputeResiduals(int* boundaryId_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d,
                                  FLOAT_TYPE* drag_d, FLOAT_TYPE* lift_d,
                                  FLOAT_TYPE* temp9a_d, FLOAT_TYPE* temp9b_d, FLOAT_TYPE* tempA_d,
                                  FLOAT_TYPE* tempB_d,
                                  FLOAT_TYPE* residuals, int* m, int* n, int computeDragLift);

/**
 * @brief Compute the residual on the GPU for BC bitmask
 * @warning Exits the whole program if the iteration diverges
 *
 * @param bcMask_d BC bitmask
 * @param f_d distribution function
 * @param fColl_d distribution function (previous step)
 * @param drag_d drag
 * @param lift_d lift
 * @param temp9a_d,temp9b_d temporary array for summation (size: 9xMxN)
 * @param tempA_d, tempB_d temporary array for summation (size: MxN)
 * @param residuals array for residuals and lift/drag
 * @param m number of columns
 * @param n number of rows
 * @param computeDragLift boundary ID for drag/lift computation
 */
__host__ void GpuComputeResidMask(int* bcMask_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d,
                                  FLOAT_TYPE* drag_d, FLOAT_TYPE* lift_d,
                                  FLOAT_TYPE* temp9a_d, FLOAT_TYPE* temp9b_d, FLOAT_TYPE* tempA_d,
                                  FLOAT_TYPE* tempB_d,
                                  FLOAT_TYPE* residuals, int* m, int* n, int computeDragLift);

/**
 * @brief Compute residual on GPU
 * @details compute the norm of the difference of f and fColl
 *  
 * @param f_d               distribution function
 * @param fColl_d           distribution function collision step
 * @param temp9a_d,temp9b_d temporary array for summation (size: 9xMxN)
 * @param m                 number of columns
 * @param n                 number of rows
 * @return norm of the difference
 */
__host__ FLOAT_TYPE computeResidual(FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d,
                                    FLOAT_TYPE *temp9a_d, FLOAT_TYPE *temp9b_d,
                                    int m, int n);

/**
 * @brief Compute drag or lift on GPU
 * 
 * @param bcMask_d         BC bitmask
 * @param dl_d             drag or lift
 * @param tempA_d, tempB_d temporary array for summation (size: MxN)
 * @param m                number of columns
 * @param n                number of rows
 * @param boundaryId       boundary to calculate drag/lift on
 * @return drag/lift
 */
__host__ FLOAT_TYPE computeDragLift(int *bcMask_d, FLOAT_TYPE *dl_d,
                                    FLOAT_TYPE *tempA_d, FLOAT_TYPE *tempB_d,
                                    int m, int n, int boundaryId);

#endif
