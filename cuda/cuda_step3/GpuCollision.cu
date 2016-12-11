/**
 * Collision model
 * @file GpuCollision.cu
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */

#include "GpuFunctions.h"
#include "GpuConstants.h"

__global__ void gpu_bgk(int* fluid_d, FLOAT_TYPE* fEq_d, FLOAT_TYPE* rho_d, FLOAT_TYPE* u_d,
                        FLOAT_TYPE* v_d, FLOAT_TYPE* fColl_d, FLOAT_TYPE* f_d)
{
    int bidx=blockIdx.x;
    int tidx=threadIdx.x;

    int ind = tidx + bidx*blockDim.x;

    int ind_s = ind - (width_d*height_d) *( (int)(ind/(width_d*height_d)) );

    int ind_c = (int)(ind/(width_d*height_d));  /// latticeId

    if (ind<(9*width_d*height_d))
    {
        if (fluid_d[ind_s]==1) //if fluid
        {
            // Collision
            ///NOTE rho * w * ( 1 + 3(u*ex + v*ey) + 4.5(u*ex + v*ey)^2 - 1.5(u^2 + v^2) )
            fEq_d[ind]  = rho_d[ind_s]*w_d[ind_c]*( 1.0 + 3.0*(u_d[ind_s]*cx_d[ind_c] + v_d[ind_s]*cy_d[ind_c])
                        + 4.5*(u_d[ind_s]*cx_d[ind_c] + v_d[ind_s]*cy_d[ind_c])
                        * (u_d[ind_s]*cx_d[ind_c] + v_d[ind_s]*cy_d[ind_c])
                        - 1.5*(u_d[ind_s] * u_d[ind_s] + v_d[ind_s]* v_d[ind_s]) );
            ///NOTE collision_freq*Feq + ( 1-coll_freq)*F
            fColl_d[ind] = omega_d*fEq_d[ind]+(1.0-omega_d)*f_d[ind];
            // printf("CF_%d: %.2f\n", ind, fColl_d[ind]);
        }
    }
}

/**
 * @brief Compute the equilibrum distribution without the collision frequency
 * @details ...
 * @todo doc: insert name or reference, explain method
 *
 * @param u,v velocity
 * @param uc,vc velocity component, see: #cx_d #cy_d
 * @param rho density
 * @param w lattice weight, see: #w_d
 * @return equlibrum distribution function
 */
__device__ FLOAT_TYPE feqc(FLOAT_TYPE u, FLOAT_TYPE uc, FLOAT_TYPE v, FLOAT_TYPE vc, FLOAT_TYPE rho, FLOAT_TYPE w)
{
    FLOAT_TYPE vec = u*uc + v*vc;
    return rho * w * (1. + 3.*vec + 4.5*vec*vec - 1.5*(u*u + v*v));
}

__global__ void gpuCollBgkw(int *fluid_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *u_d,
                            FLOAT_TYPE *v_d, FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d)
{
    int ind = blockIdx.x * blockDim.x + threadIdx.x;
    int ms = width_d*height_d;
    FLOAT_TYPE r, u, v;
    if (ind < ms && fluid_d[ind] == 1)
    {
        u = u_d[ind];
        v = v_d[ind];
        r = rho_d[ind];
        fColl_d[ind]      = omega_d * feqc(u,  0, v,  0, r, 4./9.)  + (1.0-omega_d) * f_d[ind];
        fColl_d[ind+1*ms] = omega_d * feqc(u,  1, v,  0, r, 1./9.)  + (1.0-omega_d) * f_d[ind+1*ms];
        fColl_d[ind+2*ms] = omega_d * feqc(u,  0, v,  1, r, 1./9.)  + (1.0-omega_d) * f_d[ind+2*ms];
        fColl_d[ind+3*ms] = omega_d * feqc(u, -1, v,  0, r, 1./9.)  + (1.0-omega_d) * f_d[ind+3*ms];
        fColl_d[ind+4*ms] = omega_d * feqc(u,  0, v, -1, r, 1./9.)  + (1.0-omega_d) * f_d[ind+4*ms];
        fColl_d[ind+5*ms] = omega_d * feqc(u,  1, v,  1, r, 1./36.) + (1.0-omega_d) * f_d[ind+5*ms];
        fColl_d[ind+6*ms] = omega_d * feqc(u, -1, v,  1, r, 1./36.) + (1.0-omega_d) * f_d[ind+6*ms];
        fColl_d[ind+7*ms] = omega_d * feqc(u, -1, v, -1, r, 1./36.) + (1.0-omega_d) * f_d[ind+7*ms];
        fColl_d[ind+8*ms] = omega_d * feqc(u,  1, v, -1, r, 1./36.) + (1.0-omega_d) * f_d[ind+8*ms];
    }
    // if (ind == 0)
    // {
    //     printf("GPU omega: %f\n", omega_d);
    // }
}

__global__ void gpuCollTrt(int *fluid_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *u_d,
                        FLOAT_TYPE *v_d, FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d)
{
    int ind = blockIdx.x * blockDim.x + threadIdx.x;
    int ms = width_d*height_d;
    FLOAT_TYPE r, u, v;
    if (ind < ms && fluid_d[ind] == 1)
    {
        u = u_d[ind];
        v = v_d[ind];
        r = rho_d[ind];

        FLOAT_TYPE feq0 = feqc(u,  0, v,  0, r, 4./9.);
        FLOAT_TYPE feq1 = feqc(u,  1, v,  0, r, 1./9.);
        FLOAT_TYPE feq2 = feqc(u,  0, v,  1, r, 1./9.);
        FLOAT_TYPE feq3 = feqc(u, -1, v,  0, r, 1./9.);
        FLOAT_TYPE feq4 = feqc(u,  0, v, -1, r, 1./9.);
        FLOAT_TYPE feq5 = feqc(u,  1, v,  1, r, 1./36.);
        FLOAT_TYPE feq6 = feqc(u, -1, v,  1, r, 1./36.);
        FLOAT_TYPE feq7 = feqc(u, -1, v, -1, r, 1./36.);
        FLOAT_TYPE feq8 = feqc(u,  1, v, -1, r, 1./36.);

        FLOAT_TYPE f0 = f_d[ind];
        FLOAT_TYPE f1 = f_d[ind+1*ms];
        FLOAT_TYPE f2 = f_d[ind+2*ms];
        FLOAT_TYPE f3 = f_d[ind+3*ms];
        FLOAT_TYPE f4 = f_d[ind+4*ms];
        FLOAT_TYPE f5 = f_d[ind+5*ms];
        FLOAT_TYPE f6 = f_d[ind+6*ms];
        FLOAT_TYPE f7 = f_d[ind+7*ms];
        FLOAT_TYPE f8 = f_d[ind+8*ms];

        fColl_d[ind]      = f0 - 0.5 * omega_d * (f0+f0 - feq0-feq0) - 0.5 * omegaA_d * (f0-f0 - feq0+feq0);
        fColl_d[ind+1*ms] = f1 - 0.5 * omega_d * (f1+f3 - feq1-feq3) - 0.5 * omegaA_d * (f1-f3 - feq1+feq3);
        fColl_d[ind+2*ms] = f2 - 0.5 * omega_d * (f2+f4 - feq2-feq4) - 0.5 * omegaA_d * (f2-f4 - feq2+feq4);
        fColl_d[ind+3*ms] = f3 - 0.5 * omega_d * (f3+f1 - feq3-feq1) - 0.5 * omegaA_d * (f3-f1 - feq3+feq1);
        fColl_d[ind+4*ms] = f4 - 0.5 * omega_d * (f4+f2 - feq4-feq2) - 0.5 * omegaA_d * (f4-f2 - feq4+feq2);
        fColl_d[ind+5*ms] = f5 - 0.5 * omega_d * (f5+f7 - feq5-feq7) - 0.5 * omegaA_d * (f5-f7 - feq5+feq7);
        fColl_d[ind+6*ms] = f6 - 0.5 * omega_d * (f6+f8 - feq6-feq8) - 0.5 * omegaA_d * (f6-f8 - feq6+feq8);
        fColl_d[ind+7*ms] = f7 - 0.5 * omega_d * (f7+f5 - feq7-feq5) - 0.5 * omegaA_d * (f7-f5 - feq7+feq5);
        fColl_d[ind+8*ms] = f8 - 0.5 * omega_d * (f8+f6 - feq8-feq6) - 0.5 * omegaA_d * (f8-f6 - feq8+feq6);

        // fColl_d[ind]      = f0 - omega_d * (0.5*(f0+f0) - 0.5*(feq0+feq0)) - omegaA_d * (0.5*(f0-f0) - 0.5*(feq0-feq0));
        // fColl_d[ind+1*ms] = f1 - omega_d * (0.5*(f1+f3) - 0.5*(feq1+feq3)) - omegaA_d * (0.5*(f1-f3) - 0.5*(feq1-feq3));
        // fColl_d[ind+2*ms] = f2 - omega_d * (0.5*(f2+f4) - 0.5*(feq2+feq4)) - omegaA_d * (0.5*(f2-f4) - 0.5*(feq2-feq4));
        // fColl_d[ind+3*ms] = f3 - omega_d * (0.5*(f3+f1) - 0.5*(feq3+feq1)) - omegaA_d * (0.5*(f3-f1) - 0.5*(feq3-feq1));
        // fColl_d[ind+4*ms] = f4 - omega_d * (0.5*(f4+f2) - 0.5*(feq4+feq2)) - omegaA_d * (0.5*(f4-f2) - 0.5*(feq4-feq2));
        // fColl_d[ind+5*ms] = f5 - omega_d * (0.5*(f5+f7) - 0.5*(feq5+feq7)) - omegaA_d * (0.5*(f5-f7) - 0.5*(feq5-feq7));
        // fColl_d[ind+6*ms] = f6 - omega_d * (0.5*(f6+f8) - 0.5*(feq6+feq8)) - omegaA_d * (0.5*(f6-f8) - 0.5*(feq6-feq8));
        // fColl_d[ind+7*ms] = f7 - omega_d * (0.5*(f7+f5) - 0.5*(feq7+feq5)) - omegaA_d * (0.5*(f7-f5) - 0.5*(feq7-feq5));
        // fColl_d[ind+8*ms] = f8 - omega_d * (0.5*(f8+f6) - 0.5*(feq8+feq6)) - omegaA_d * (0.5*(f8-f6) - 0.5*(feq8-feq6));
    }
}

__global__ void gpu_trt1(int* fluid_d, FLOAT_TYPE* fEq_d, FLOAT_TYPE* rho_d, FLOAT_TYPE* u_d, FLOAT_TYPE* v_d)
{
    int bidx=blockIdx.x;
    int tidx=threadIdx.x;

    int ind = tidx + bidx*blockDim.x;

    int ind_s = ind - (width_d*height_d) *( (int)(ind/(width_d*height_d)) );

    int ind_c = (int)(ind/(width_d*height_d));

    if (ind<(9*width_d*height_d))
    {
        if (fluid_d[ind_s]==1)
        {
            fEq_d[ind]  = rho_d[ind_s]*w_d[ind_c] * (1.0+3.0*(u_d[ind_s]*cx_d[ind_c]+v_d[ind_s]*cy_d[ind_c])
                        + 4.5 *(u_d[ind_s]*cx_d[ind_c]+v_d[ind_s]*cy_d[ind_c]) * (u_d[ind_s]*cx_d[ind_c]+v_d[ind_s]*cy_d[ind_c])
                        - 1.5*(u_d[ind_s]*u_d[ind_s]+v_d[ind_s]*v_d[ind_s]));
        }
    }
}

__global__ void gpu_trt2(int* fluid_d, FLOAT_TYPE* fEq_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d)
{
    int bidx=blockIdx.x;
    int tidx=threadIdx.x;

    int ind = tidx + bidx*blockDim.x;

    int ind_s = ind - (width_d*height_d) *( (int)(ind/(width_d*height_d)) );

    int ind_c = (int)(ind/(width_d*height_d));

    if (ind<(9*width_d*height_d))
    {
        if (fluid_d[ind_s]==1)
        {
            fColl_d[ind]    = f_d[ind]
                            - ((0.5*(f_d[ind]   + f_d[ind_s+opp_d[ind_c]]))
                            - (0.5* (fEq_d[ind] + fEq_d[ind_s+opp_d[ind_c]])))*omega_d
                            - ((0.5*(f_d[ind]   - f_d[ind_s+opp_d[ind_c]]))
                            - (0.5* (fEq_d[ind] - fEq_d[ind_s+opp_d[ind_c]])))*omegaA_d;
        }
    }
}

__global__ void gpuCollMrt(int* fluid_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *u_d, FLOAT_TYPE *v_d, FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d)
{
    int ind = blockIdx.x * blockDim.x + threadIdx.x;
    int ms = width_d*height_d;
    FLOAT_TYPE mEq[9], m[9], collision[9], f[9];

    FLOAT_TYPE r,u,v;

    if (ind < ms && fluid_d[ind] == 1)
    {
        r = rho_d[ind];
        u = u_d[ind];
        v = v_d[ind];

        // m^eq = (rho, e, epsilon, jx, qx, jy, qy, pxx, pxy)^T
        mEq[0] = r;
        mEq[1] = r * (-2. + 3. * r * (u*u + v*v));
        mEq[2] = r * (1.  - 3. * r * (u*u + v*v));
        mEq[3] =  r * u;
        mEq[4] = -r * u;
        mEq[5] =  r * v;
        mEq[6] = -r * v;
        mEq[7] = r * (u*u - v*v);
        mEq[8] = r * u * v;

        f[0] = f_d[ind];
        f[1] = f_d[ind+ms];
        f[2] = f_d[ind+2*ms];
        f[3] = f_d[ind+3*ms];
        f[4] = f_d[ind+4*ms];
        f[5] = f_d[ind+5*ms];
        f[6] = f_d[ind+6*ms];
        f[7] = f_d[ind+7*ms];
        f[8] = f_d[ind+8*ms];

        // m = Mf
        for (int i=0; i<9; ++i)
        {
            m[i] = 0;
            for (int j=0; j<9; ++j)
            {
                m[i] += velMomMap_d[i*9+j] * f[j];
            }
        }

        // Diff = M^-1 * S * (m - m^eq)
        for (int i=0; i<9; ++i)
        {
            collision[i] = 0;
            for (int j=0; j<9; ++j)
            {
                collision[i] += momCollMtx_d[i*9+j] * (m[j] - mEq[j]);
            }
        }

        // fColl = f - M^-1 * S * (m - m^eq) >>> MRT equation
        fColl_d[ind     ] = f[0] - collision[0];
        fColl_d[ind+  ms] = f[1] - collision[1];
        fColl_d[ind+2*ms] = f[2] - collision[2];
        fColl_d[ind+3*ms] = f[3] - collision[3];
        fColl_d[ind+4*ms] = f[4] - collision[4];
        fColl_d[ind+5*ms] = f[5] - collision[5];
        fColl_d[ind+6*ms] = f[6] - collision[6];
        fColl_d[ind+7*ms] = f[7] - collision[7];
        fColl_d[ind+8*ms] = f[8] - collision[8];
    }
}

__global__ void gpu_mrt1(int* fluid_d, FLOAT_TYPE* rho_d, FLOAT_TYPE* u_d, FLOAT_TYPE* v_d, FLOAT_TYPE* f_d, FLOAT_TYPE* mEq_d, FLOAT_TYPE* m_d)
{

    int l;

    int bidx=blockIdx.x;
    int tidx=threadIdx.x;

    int ind = tidx + bidx*blockDim.x;

    int ms=width_d*height_d;
    int ind_s = ind - (width_d*height_d) *( (int)(ind/(width_d*height_d)) );

    int ind_c = (int)(ind/(width_d*height_d));

    if (ind<(9*width_d*height_d))
    {
        if (fluid_d[ind_s]==1)
        {
            mEq_d[ind_s]      = rho_d[ind_s];
            mEq_d[ind_s+1*ms] = rho_d[ind_s]*(-2.0 + 3.0*rho_d[ind_s] * (u_d[ind_s]*u_d[ind_s] + v_d[ind_s]*v_d[ind_s]));
            mEq_d[ind_s+2*ms] = rho_d[ind_s]*(1.0 -  3.0*rho_d[ind_s] * (u_d[ind_s]*u_d[ind_s] + v_d[ind_s]*v_d[ind_s]));
            mEq_d[ind_s+3*ms] = rho_d[ind_s]*u_d[ind_s];
            mEq_d[ind_s+4*ms] =-rho_d[ind_s]*u_d[ind_s];
            mEq_d[ind_s+5*ms] = rho_d[ind_s]*v_d[ind_s];
            mEq_d[ind_s+6*ms] =-rho_d[ind_s]*v_d[ind_s];
            mEq_d[ind_s+7*ms] = rho_d[ind_s] * (u_d[ind_s]*u_d[ind_s] - v_d[ind_s]*v_d[ind_s]);
            mEq_d[ind_s+8*ms] = rho_d[ind_s]*u_d[ind_s]*v_d[ind_s];


            m_d[ind]=0;
            for (l=0; l<9;l++)
                m_d[ind]=m_d[ind] + velMomMap_d[ind_c*9+l]*f_d[ind_s+l*ms];
        }
    }
}

__global__ void gpu_mrt2(int* fluid_d, FLOAT_TYPE* collision_d, FLOAT_TYPE* m_d, FLOAT_TYPE* mEq_d, FLOAT_TYPE* fColl_d, FLOAT_TYPE* f_d)
{
    int l;

    int bidx=blockIdx.x;
    int tidx=threadIdx.x;

    int ind = tidx + bidx*blockDim.x;

    int ms=width_d*height_d;
    int ind_s = ind - (width_d*height_d) *( (int)(ind/(width_d*height_d)) );

    int ind_c = (int)(ind/(width_d*height_d));

    if (ind<(9*width_d*height_d))
    {
        if (fluid_d[ind_s]==1)
        {
            collision_d[ind] = 0;

            for (l=0; l<9;l++)
                collision_d[ind] = collision_d[ind] + (momCollMtx_d[ind_c*9+l]*(m_d[ind_s+l*ms]-mEq_d[ind_s+l*ms]));

            fColl_d[ind] = f_d[ind] - collision_d[ind];
        }
    }
}
