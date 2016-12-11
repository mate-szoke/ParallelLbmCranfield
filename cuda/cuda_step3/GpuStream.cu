#include <stdio.h>

#include "GpuFunctions.h"
#include "BcMacros.h"
#include "GpuConstants.h"

__global__ void gpu_streaming(int* fluid_d, int* stream_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d)
{
    int bidx=blockIdx.x;
    int tidx=threadIdx.x;

    int ind = tidx + bidx*blockDim.x;

    int ind_s = ind - (width_d*height_d) *( (int)(ind/(width_d*height_d)) );

    int ind_c = (int)(ind/(width_d*height_d));

    int ms = width_d*height_d;

    if (ind<(9*width_d*height_d))
    {
        if (fluid_d[ind_s]==1)
        {
            // STREAMING
            if ( (ind > ms && stream_d[ind-ms]) == 1)
            {
                f_d[ind] = fColl_d[ind+c_d[ind_c]];
            }
            else
            {
                f_d[ind] = fColl_d[ind];
            }
        }
    }
}

__global__ void gpuStreaming(int* fluid_d, int* stream_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d)
{
    int ind = blockIdx.x * blockDim.x + threadIdx.x;
    int ms = width_d*height_d;
    FLOAT_TYPE *f, *mf;
    int n = height_d;
    if (ind < ms && fluid_d[ind] == 1)
    {
        f_d[ind] = fColl_d[ind];
        f = f_d + ms;
        mf = fColl_d + ms;
        f[ind]      = (stream_d[ind]      == 1) ? mf[ind-1]        : mf[ind];
        f[ind+ms]   = (stream_d[ind+ms]   == 1) ? mf[ind+ms-n]     : mf[ind+ms];
        f[ind+2*ms] = (stream_d[ind+2*ms] == 1) ? mf[ind+2*ms+1]   : mf[ind+2*ms];
        f[ind+3*ms] = (stream_d[ind+3*ms] == 1) ? mf[ind+3*ms+n]   : mf[ind+3*ms];
        f[ind+4*ms] = (stream_d[ind+4*ms] == 1) ? mf[ind+4*ms-n-1] : mf[ind+4*ms];
        f[ind+5*ms] = (stream_d[ind+5*ms] == 1) ? mf[ind+5*ms-n+1] : mf[ind+5*ms];
        f[ind+6*ms] = (stream_d[ind+6*ms] == 1) ? mf[ind+6*ms+n+1] : mf[ind+6*ms];
        f[ind+7*ms] = (stream_d[ind+7*ms] == 1 && ind < ms-n+1) ? mf[ind+7*ms+n-1] : mf[ind+7*ms];
    }
}
