#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "ComputeResiduals.h"
#include "GpuFunctions.h"
#include "FilesReading.h"
#include "ShellFunctions.h"
#include "ArrayUtils.h"
#include "BcMacros.h"
#include "GpuSum.h"

void ComputeResiduals(int* boundaryId, FLOAT_TYPE* f, FLOAT_TYPE* fColl, FLOAT_TYPE* drag,
                      FLOAT_TYPE* lift, FLOAT_TYPE* residuals, int* m, int* n, int computeDragLift)
{
  // Create variables for residuals
  FLOAT_TYPE ResDrag=0.0;
  FLOAT_TYPE ResLift=0.0;
  // Loop variables
  int i, k;

  FLOAT_TYPE L2n   = 0.0;  // L2 norm
  FLOAT_TYPE L2n_w = 0.0;  // L2 norm weighted

  // sum up velocity and density
  for (i=0; i<((*m)*(*n)); i++)
  {
    for (k=0; k<9; k++)
    { L2n = L2n + pow((f[i+k*((*m)*(*n))]-fColl[i+k*((*m)*(*n))]),2); }

    if (boundaryId[i]==computeDragLift)  // BE AWARE!! check in the STAR-CD files the ID of the surface where you want to check the drag/lift
    {
      ResDrag+=drag[i];
      ResLift+=lift[i];
    }
  }
  // Calculate residuals
  L2n_w = sqrt(L2n/((*m)*(*n)));
  L2n = sqrt(L2n);
  // Write them to vector
  residuals[0] = L2n;
  residuals[1] = L2n_w;
  // Lift and drag
  residuals[2] = ResDrag;
  residuals[3] = ResLift;

  if(L2n!=L2n) // if density residuals are NaN
  {
    printf("\nDIVERGENCE!\n");
    exit(1); // ERROR!
  }
}

__host__ void GpuComputeResiduals(int* boundaryId_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d,
                                  FLOAT_TYPE* drag_d, FLOAT_TYPE* lift_d,
                                  FLOAT_TYPE* temp9a_d, FLOAT_TYPE* temp9b_d, FLOAT_TYPE* tempA_d,
                                  FLOAT_TYPE* tempB_d,
                                  FLOAT_TYPE* residuals, int* m, int* n, int computeDragLift)
{
  // Create variables for residuals
  FLOAT_TYPE ResDrag=0.0;
  FLOAT_TYPE ResLift=0.0;

  FLOAT_TYPE L2n   = 0.0;  // L2 norm
  FLOAT_TYPE L2n_w = 0.0;  // L2 norm weighted

  dim3 grid_dim((int)(9*(*m)*(*n)/THREADS)+1);
  gpu_sqsub <<<grid_dim, THREADS>>>(f_d, fColl_d, temp9a_d, 9*(*m)*(*n));
  L2n = gpu_sum_h(temp9a_d, temp9b_d, 9*(*m)*(*n));

  gpu_cond_copy <<<grid_dim, THREADS>>>(tempA_d, drag_d, boundaryId_d, computeDragLift, (*m)*(*n));
  ResDrag = gpu_sum_h(tempA_d, tempB_d, (*m)*(*n));
  gpu_cond_copy <<<grid_dim, THREADS>>>(tempA_d, lift_d, boundaryId_d, computeDragLift, (*m)*(*n));
  ResLift = gpu_sum_h(tempA_d, tempB_d, (*m)*(*n));

  // Calculate residuals
  L2n_w = sqrt(L2n/((*m)*(*n)));
  L2n = sqrt(L2n);
  // Write them to vector
  residuals[0] = L2n;
  residuals[1] = L2n_w;

  // Lift and drag
  residuals[2] = ResDrag;
  residuals[3] = ResLift;


  if(L2n!=L2n) // if density residuals are NaN
  {
    printf("\nDIVERGENCE! (%g)\n", L2n);
    exit(1); // ERROR!
  }
}

__host__ void GpuComputeResidMask(int* bcMask_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d,
                                  FLOAT_TYPE* drag_d, FLOAT_TYPE* lift_d,
                                  FLOAT_TYPE* temp9a_d, FLOAT_TYPE* temp9b_d, FLOAT_TYPE* tempA_d,
                                  FLOAT_TYPE* tempB_d,
                                  FLOAT_TYPE* residuals, int* m, int* n, int computeDragLift)
{
  FLOAT_TYPE ResDrag=0.0;
  FLOAT_TYPE ResLift=0.0;

  FLOAT_TYPE L2n   = 0.0;  // L2 norm
  FLOAT_TYPE L2n_w = 0.0;  // L2 norm weighted

  dim3 grid_dim((int)(9*(*m)*(*n)/THREADS)+1);
  gpu_sqsub <<<grid_dim, THREADS>>>(f_d, fColl_d, temp9a_d, 9*(*m)*(*n));
  L2n = gpu_sum_h(temp9a_d, temp9b_d, 9*(*m)*(*n));

  gpu_cond_copy_mask <<<grid_dim, THREADS>>>(tempA_d, drag_d, bcMask_d, computeDragLift, (*m)*(*n));
  ResDrag = gpu_sum_h(tempA_d, tempB_d, (*m)*(*n));
  gpu_cond_copy_mask <<<grid_dim, THREADS>>>(tempA_d, lift_d, bcMask_d, computeDragLift, (*m)*(*n));
  ResLift = gpu_sum_h(tempA_d, tempB_d, (*m)*(*n));

  L2n_w = sqrt(L2n/((*m)*(*n)));
  L2n = sqrt(L2n);

  residuals[0] = L2n;
  residuals[1] = L2n_w;
  residuals[2] = ResDrag;
  residuals[3] = ResLift;

  if(L2n!=L2n) // if density residuals are NaN
  {
    printf("\nDIVERGENCE! (%g)\n", L2n);
    exit(1); // ERROR!
  }
}

__host__ FLOAT_TYPE computeResidual(FLOAT_TYPE *f_d, FLOAT_TYPE *fTemp_d,
                                    FLOAT_TYPE *temp9a_d, FLOAT_TYPE *temp9b_d,
                                    int m, int n)
{
  dim3 bpg9((9*m*n-1)/THREADS+1);
  gpu_sqsub <<<bpg9,THREADS>>>(f_d, fTemp_d, temp9a_d, 9*m*n);
  return sqrt(gpu_sum_h(temp9a_d, temp9b_d, 9*m*n));
}

__host__ FLOAT_TYPE computeDragLift(int *bcMask_d, FLOAT_TYPE *dl_d,
                                    FLOAT_TYPE *tempA_d, FLOAT_TYPE *tempB_d,
                                    int m, int n, int boundaryId)
{
  dim3 bpg((m*n-1)/THREADS+1);
  gpu_cond_copy_mask <<<bpg,THREADS>>>(tempA_d, dl_d, bcMask_d, boundaryId, m*n);
  return gpu_sum_h(tempA_d, tempB_d, m*n);
}
