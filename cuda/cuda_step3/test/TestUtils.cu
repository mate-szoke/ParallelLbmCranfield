#include "TestUtils.h"
#include "GpuFunctions.h"
#include "BcMacros.h"
#include "ArrayUtils.h"
#include "GpuSum.h"
#include <cuda.h>
#include <stdio.h>

__global__ void gpuArrayAddInt(int *array1_d, int *array2_d, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i<size) array1_d[i] += array2_d[i];
}

__global__ void gpuArrayAddFlt(FLOAT_TYPE *array1_d, FLOAT_TYPE *array2_d, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i<size) array1_d[i] += array2_d[i];
}

void hostArrayAddInt(int *array1_h, int *array2_h, int size)
{
    int i;
    for (i=0; i<size; ++i) array1_h[i] += array2_h[i];
}

void hostArrayAddFlt(FLOAT_TYPE *array1_h, FLOAT_TYPE *array2_h, int size)
{
    int i;
    for (i=0; i<size; ++i) array1_h[i] += array2_h[i];
}

int sumHostInt(int *array_h, int size)
{
    int i,s=0;
    for (i=0; i<size; ++i) s += array_h[i];
    return s;
}

FLOAT_TYPE sumHostFlt(FLOAT_TYPE *array_h, int size)
{
    int i;
    FLOAT_TYPE s=0.;
    for (i=0; i<size; ++i) s += array_h[i];
    return s;
}

void createLidBcFluid(int *fluid, int m, int n)
{
    int i;
    for(i=0; i<m*n; ++i) fluid[i] = 1;
}

void createLidBcBoundary(int *boundary, int m, int n)
{
    int i;
    for(i=0; i<m*n; ++i)
    {
        if (i%m == 0 || i%m == m-1 || i<m || i>=m*(n-1)) boundary[i] = 1;
        else boundary[i] = 0;
        if (i >= (n-1)*m) boundary[i] = 2; //top
    }
}

void createLidBcBcId(int *bcid, int m, int n)
{
    int i,j;
    for(i=0; i<m*n; ++i)
    {
        for (j=0; j<9; ++j)
        {
            bcid[i+j*n*m] = 0;
        }
        if (i%m == 0)
        {
            bcid[i+3*n*m] = bcid[i+6*n*m] = bcid[i+7*n*m] = 1;
        }
        if (i%m == m-1)
        {
            bcid[i+1*n*m] = bcid[i+5*n*m] = bcid[i+8*n*m] = 1;
        }
        if (i < m)
        {
            bcid[i+4*n*m] = bcid[i+7*n*m] = bcid[i+8*n*m] = 1;
        }
        if (i >= (n-1)*m) //top
        {
            bcid[i+2*n*m] = bcid[i+5*n*m] = bcid[i+6*n*m] = 2;
        }
    }
    bcid[(n-1)*m+6*n*m] = bcid[n*m-1+5*n*m] = 1;
}

void createLidBcCorner(int *corner, int m, int n)
{
    int i;
    for(i=0; i<m*n; ++i) corner[i] = 0;

    corner[(n-1)*m] = corner[m*n-1] = 1;
}

int getLidBcMaskSize(int m, int n)
{
    return 2*m+2*n-4;
}

void createLidBcIdx(int *index, int m, int n)
{
    int i,j=0;
    for (i=0; i<m*n; ++i)
    {
        if (i%m == 0 || i%m == m-1 || i < m || i >= (n-1)*m)
        {
            index[j++] = i;
        }
    }
    if (j != getLidBcMaskSize(m,n))
    {
        fprintf(stderr, "INDEX Size mismatch %d != %d\n", j, getLidBcMaskSize(m,n));
    }
}

void createLidBcMask(int *mask, int m, int n)
{
    int i,j=0;
    for (i=0; i<m*n; ++i)
    {
        if (i%m == 0 || i%m == m-1 || i < m || i >= (n-1)*m)
        {
            if (i%m == 0) //west
            {
                mask[j] = BC_FLUID | BC_WALL_B | BC_WALL_W | BC_WALL_NW | BC_WALL_SW | BOUND_ID(3);
            }
            if (i%m == m-1) //east
            {
                mask[j] = BC_FLUID | BC_WALL_B | BC_WALL_E | BC_WALL_NE | BC_WALL_SE | BOUND_ID(1);
            }
            if (i < m) //south
            {
                mask[j] = BC_FLUID | BC_WALL_B | BC_WALL_S | BC_WALL_SE | BC_WALL_SW | BOUND_ID(4);
            }
            if (i >= (n-1)*m) //north
            {
                mask[j] = BC_FLUID | BC_INLT_B | BC_INLT_N | BC_INLT_NE | BC_INLT_NW | BOUND_ID(2);
            }
            ++j;
        }
    }
    if (j != getLidBcMaskSize(m,n))
    {
        fprintf(stderr, "MASK Size mismatch! %d != %d\n", j, getLidBcMaskSize(m,n));
    }
    //nw corner
    mask[j-m] = BC_FLUID | BC_CORNER | BC_WALL_B | BC_WALL_NW | BC_WALL_W | BC_WALL_SW | BC_INLT_N | BC_INLT_NE | BOUND_ID(2);
    //ne corner
    mask[j-1] = BC_FLUID | BC_CORNER | BC_WALL_B | BC_WALL_NE | BC_WALL_E | BC_WALL_SE | BC_INLT_N | BC_INLT_NW | BOUND_ID(2);
}

void createLidBcMaskFull(int *mask, int m, int n)
{
    int i;
    for (i=0; i<m*n; ++i)
    {
        if (i%m == 0 || i%m == m-1 || i < m || i >= (n-1)*m)
        {
            if (i%m == 0) //west
            {
                mask[i] = BC_FLUID | BC_WALL_B | BC_WALL_W | BC_WALL_NW | BC_WALL_SW | BOUND_ID(3);
            }
            if (i%m == m-1) //east
            {
                mask[i] = BC_FLUID | BC_WALL_B | BC_WALL_E | BC_WALL_NE | BC_WALL_SE | BOUND_ID(1);
            }
            if (i < m) //south
            {
                mask[i] = BC_FLUID | BC_WALL_B | BC_WALL_S | BC_WALL_SE | BC_WALL_SW | BOUND_ID(4);
            }
            if (i >= (n-1)*m) //north
            {
                mask[i] = BC_FLUID | BC_INLT_B | BC_INLT_N | BC_INLT_NE | BC_INLT_NW | BOUND_ID(2);
            }
        }
    }
    //nw corner
    mask[m*n-m] = BC_FLUID | BC_CORNER | BC_WALL_B | BC_WALL_NW | BC_WALL_W | BC_WALL_SW | BC_INLT_N | BC_INLT_NE | BOUND_ID(2);
    //ne corner
    mask[m*n-1] = BC_FLUID | BC_CORNER | BC_WALL_B | BC_WALL_NE | BC_WALL_E | BC_WALL_SE | BC_INLT_N | BC_INLT_NW | BOUND_ID(2);
}

void createLidBoundaryId(int *bid, int m, int n)
{
    int i;
    for (i=0; i<m*n; ++i)
    {
        if (i%m == 0) //west
            {
                bid[i] = 3;
            }
            if (i%m == m-1) //east
            {
                bid[i] = 1;
            }
            if (i < m) //south
            {
                bid[i] = 4;
            }
            if (i >= (n-1)*m) //north
            {
                bid[i] = 2;
            }
    }
}

void fillLidCoordinates(FLOAT_TYPE *x, FLOAT_TYPE *y, int m, int n)
{
    int i,j;
    for(i=0; i<m; ++i)
    {
        for(j=0; j<n; ++j)
        {
            x[i*m+j] = (FLOAT_TYPE)i / (FLOAT_TYPE)m;
            y[i*m+j] = (FLOAT_TYPE)j / (FLOAT_TYPE)n;
        }
    }
}

void fillFDir(FLOAT_TYPE *f, FLOAT_TYPE fu, int dir, int size)
{
    int i;
    for (i=0; i<size; ++i) f[i+dir*size] = fu;
}

void createChannelBcMask(int *mask, int m, int n)
{
    int i,j=0;
    for (i=0; i<m*n; ++i)
    {
        if (i%m == 0 || i%m == m-1 || i < m || i >= (n-1)*m)
        {
            if (i%m == 0) //west
            {
                mask[j] |= BC_FLUID | BC_INLT_B | BC_INLT_W | BC_INLT_NW | BC_INLT_SW;
            }
            if (i%m == m-1) //east
            {
                mask[j] |= BC_FLUID | BC_OUTL_B | BC_OUTL_E | BC_OUTL_NE | BC_OUTL_SE;
            }
            if (i < m) //south
            {
                mask[j] |= BC_FLUID | BC_WALL_B | BC_WALL_S | BC_WALL_SE | BC_WALL_SW;
            }
            if (i >= (n-1)*m) //north
            {
                mask[j] |= BC_FLUID | BC_WALL_B | BC_WALL_N | BC_WALL_NE | BC_WALL_NW;
            }
            ++j;
        }
    }
    if (j != getLidBcMaskSize(m,n))
    {
        fprintf(stderr, "MASK Size mismatch! %d != %d\n", j, getLidBcMaskSize(m,n));
    }
    //nw corner
    mask[j-m] = BC_FLUID | BC_CORNER | BC_WALL_B | BC_WALL_NW | BC_INLT_W | BC_INLT_SW | BC_WALL_N | BC_WALL_NE;
    //ne corner
    mask[j-1] = BC_FLUID | BC_CORNER | BC_WALL_B | BC_WALL_NE | BC_OUTL_E | BC_OUTL_SE | BC_WALL_N | BC_WALL_NW;
    //sw corner
    mask[0]   = BC_FLUID | BC_CORNER | BC_WALL_B | BC_WALL_SW | BC_WALL_S | BC_WALL_SE | BC_INLT_W | BC_INLT_NW;
    //se corner
    mask[m-1] = BC_FLUID | BC_CORNER | BC_WALL_B | BC_WALL_SE | BC_WALL_S | BC_WALL_SW | BC_OUTL_E | BC_OUTL_NE;
}



void createChannelBcBcId(int *bcid, int m, int n)
{
    int i,j;
    for(i=0; i<m*n; ++i)
    {
        for (j=0; j<9; ++j)
        {
            bcid[i+j*n*m] = 0;
        }
        if (i%m == 0) //west
        {
            bcid[i+3*n*m] = bcid[i+6*n*m] = bcid[i+7*n*m] = 2;
        }
        if (i%m == m-1) //east
        {
            bcid[i+1*n*m] = bcid[i+5*n*m] = bcid[i+8*n*m] = 3;
        }
        if (i < m) //south
        {
            bcid[i+4*n*m] = bcid[i+7*n*m] = bcid[i+8*n*m] = 1;
        }
        if (i >= (n-1)*m) //north
        {
            bcid[i+2*n*m] = bcid[i+5*n*m] = bcid[i+6*n*m] = 1;
        }
    }
    bcid[(n-1)*m+6*n*m] = bcid[n*m-1+5*n*m] = bcid[7*m*n] = bcid[m-1+8*m*n] = 1;
}

void createChannelBcCorner(int *corner, int m, int n)
{
    int i;
    for(i=0; i<m*n; ++i) corner[i] = 0;

    corner[(n-1)*m] = corner[m*n-1] = corner[0] = corner[m-1] = 1;
}

int compareArraysInt(int *array1_h, int *array2_h, int size)
{
    int i,n=0;
    for (i=0; i<size; ++i)
    {
        if (array1_h[i] != array2_h[i])
        {
            ++n;
            printf("%d: %d != %d\n", i, array1_h[i], array2_h[i]);
        }
    }
    return n;
}

int compareArraysFlt(FLOAT_TYPE *array1_h, FLOAT_TYPE *array2_h, int size)
{
    int i,n=0;
    for (i=0; i<size; ++i)
    {
        if (array1_h[i] != array2_h[i])
        {
            ++n;
            printf("%d: %g != %g\n", i, array1_h[i], array2_h[i]);
        }
    }
    return n;
}

__host__ int compareGpuArrayInt(int *a_d, int *b_d, int size)
{
  int *a = (int*)malloc(size*sizeof(int));
  int *b = (int*)malloc(size*sizeof(int));
  cudaMemcpy(a, a_d, size*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(b, b_d, size*sizeof(int), cudaMemcpyDeviceToHost);
  int n = compareArraysInt(a, b, size);
  free(a);
  free(b);
  return n;
}

__host__ int compareGpuArrayFlt(FLOAT_TYPE *a_d, FLOAT_TYPE *b_d, int size)
{
  FLOAT_TYPE *a = (FLOAT_TYPE*)malloc(size*sizeof(FLOAT_TYPE));
  FLOAT_TYPE *b = (FLOAT_TYPE*)malloc(size*sizeof(FLOAT_TYPE));
  cudaMemcpy(a, a_d, size*sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);
  cudaMemcpy(b, b_d, size*sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);
  int n = compareArraysFlt(a, b, size);
  free(a);
  free(b);
  return n;
}

__host__ FLOAT_TYPE compareGpuArraySumFlt(FLOAT_TYPE *F1_d, FLOAT_TYPE *F2_d, FLOAT_TYPE *temp9a_d, FLOAT_TYPE *temp9b_d, int size)
{
  dim3 gridDim(size/THREADS+1);
  gpu_sqsub<<<gridDim, THREADS>>>(F1_d, F2_d, temp9a_d, size);
  return sqrt(gpu_sum_h(temp9a_d, temp9b_d, size));
}

__host__ FLOAT_TYPE compareGpuArraySumInt(int *F1_d, int *F2_d, FLOAT_TYPE *temp9a_d, FLOAT_TYPE *temp9b_d, int size)
{
  dim3 gridDim(size/THREADS+1);
  gpu_sqsubi<<<gridDim, THREADS>>>(F1_d, F2_d, temp9a_d, size);
  return gpu_sum_h(temp9a_d, temp9b_d, size);
}

void printBanner(const char *testname)
{
    int l = strlen(testname);
    if (l >= MAX_LEN-4)
    {
        printf("\n%s\n", testname);
    }
    else
    {
        printf("\n");
        int i,sl = (MAX_LEN - l) / 2;
        for (i=0; i<sl-1; ++i)
        {
            printf("=");
        }
        printf(" %s ", testname);
        sl = MAX_LEN - sl - 1 - l;
        for (i=0; i<sl-1; ++i)
        {
            printf("=");
        }
        printf("\n");
    }
}
