#include "GpuFunctions.h"
#include "BcMacros.h"
#include "GpuConstants.h"

__global__ void gpu_update_new(int *fluid_d, FLOAT_TYPE* rho_d, FLOAT_TYPE* u_d, FLOAT_TYPE* v_d,
                               int *bcMask_d, FLOAT_TYPE* drag_d, FLOAT_TYPE* lift_d,
                               FLOAT_TYPE* coordX_d, FLOAT_TYPE* coordY_d, FLOAT_TYPE* f_d)
{
    int k;

    int bidx=blockIdx.x;
    int tidx=threadIdx.x;

    int ind = tidx + bidx*blockDim.x;

    int ms = width_d*height_d;

    if (ind<(width_d*height_d))
    {
        if (fluid_d[ind]==1)
        {
            // Update macroscopic
            rho_d[ind]=0;
            u_d[ind]=0;
            v_d[ind]=0;
            for (k=0; k<9; k++)
            {
                rho_d[ind] = rho_d[ind] + f_d[ind+k*ms];
                u_d[ind] = u_d[ind] + f_d[ind+k*ms]*cx_d[k];
                v_d[ind]= v_d[ind] + f_d[ind+k*ms]*cy_d[k];
            }

            u_d[ind] = u_d[ind] / rho_d[ind];
            v_d[ind] = v_d[ind] / rho_d[ind];

            if ((bcMask_d[ind] & BC_OUTL_E) == BC_OUTL_E) // for outlet on the right
            {
                ///@todo code: probably should handle outlet on other sides
                v_d[ind]=0.0;
            }

            //   DRAG/LIFT FORCE
            if (dlBoundaryId_d != 0 && (bcMask_d[ind] & BND_ID_ALL) == BOUND_ID(dlBoundaryId_d))
            {
                drag_d[ind] = 0.33333333*rho_d[ind]*(20-coordX_d[ind])*0.2;
                lift_d[ind] = 0.33333333*rho_d[ind]*(20-coordY_d[ind])*0.2;
            }
        }
    }
}

__global__ void gpuUpdateMacro(int *fluid_d, FLOAT_TYPE* rho_d, FLOAT_TYPE* u_d, FLOAT_TYPE* v_d,
                               int *bcMask_d, FLOAT_TYPE* drag_d, FLOAT_TYPE* lift_d,
                               FLOAT_TYPE* coordX_d, FLOAT_TYPE* coordY_d, FLOAT_TYPE* f_d)
{
    int ind = threadIdx.x + blockIdx.x*blockDim.x;
    int ms = width_d*height_d;

    FLOAT_TYPE r,u,v;

    if (ind < ms)
    {
        if (fluid_d[ind]==1)
        {
            r = u = v = 0.0;
            r = f_d[ind     ] + f_d[ind+  ms] + f_d[ind+2*ms] + f_d[ind+3*ms] + f_d[ind+4*ms] + f_d[ind+5*ms] +
                f_d[ind+6*ms] + f_d[ind+7*ms] + f_d[ind+8*ms];
            u = f_d[ind+  ms] - f_d[ind+3*ms] + f_d[ind+5*ms] - f_d[ind+6*ms] - f_d[ind+7*ms] + f_d[ind+8*ms];
            v = f_d[ind+2*ms] - f_d[ind+4*ms] + f_d[ind+5*ms] + f_d[ind+6*ms] - f_d[ind+7*ms] - f_d[ind+8*ms];

            rho_d[ind] = r;
            u_d[ind] = u / r;
            ///@todo code: probably should handle outlet on other sides
            v_d[ind] = ((bcMask_d[ind] & BC_OUTL_E) == BC_OUTL_E) ? 0.0 : v / r;

            //   DRAG/LIFT FORCE
            if (dlBoundaryId_d != 0 && (bcMask_d[ind] & BND_ID_ALL) == BOUND_ID(dlBoundaryId_d))
            {
                // printf("draglift: %d\n",ind);
                drag_d[ind] = 0.33333333*r*(20-coordX_d[ind])*0.2;
                lift_d[ind] = 0.33333333*r*(20-coordY_d[ind])*0.2;
            }
        }
    }
}

__global__ void gpu_update_macro(int* fluid_d, FLOAT_TYPE* rho_d, FLOAT_TYPE* u_d, FLOAT_TYPE* v_d,
                                int* bcId_d, int* boundaryId_d, FLOAT_TYPE* drag_d,
                                FLOAT_TYPE* lift_d, FLOAT_TYPE* coordX_d, FLOAT_TYPE* coordY_d, FLOAT_TYPE* f_d)
{
    int k;

    int bidx=blockIdx.x;
    int tidx=threadIdx.x;

    int ind = tidx + bidx*blockDim.x;

    int ms = width_d*height_d;

    if (ind<(width_d*height_d))
    {
        if (fluid_d[ind]==1)
        {
            // Update macroscopic
            rho_d[ind]=0;
            u_d[ind]=0;
            v_d[ind]=0;
            for (k=0; k<9; k++)
            {
                rho_d[ind] = rho_d[ind] + f_d[ind+k*ms];
                u_d[ind] = u_d[ind] + f_d[ind+k*ms]*cx_d[k];
                v_d[ind]= v_d[ind] + f_d[ind+k*ms]*cy_d[k];
            }

            u_d[ind] = u_d[ind] / rho_d[ind];
            v_d[ind] = v_d[ind] / rho_d[ind];

            if (bcId_d[ind+ms]==3) // for outlet on the right
            {
                ///@todo code: probably should handle outlet on other sides
                v_d[ind]=0.0;
            }

            //   DRAG/LIFT FORCE
            if (dlBoundaryId_d != 0 && boundaryId_d[ind]==dlBoundaryId_d)
            {
                drag_d[ind] = 0.33333333*rho_d[ind]*(20-coordX_d[ind])*0.2;
                lift_d[ind] = 0.33333333*rho_d[ind]*(20-coordY_d[ind])*0.2;
            }
        }
    }
}