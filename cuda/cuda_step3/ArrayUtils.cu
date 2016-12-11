/**
 * Memory allocation wrappers for host and gpu arrays
 * @file ArrayUtils.cu
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <time.h>
#include <assert.h>    // ensure successfull allocation
#include "ArrayUtils.h"

__constant__ FLOAT_TYPE fill_fd; ///< value to fill into floating point array
__constant__ int        fill_id; ///< value to fill into integer array

FLOAT_TYPE getRandom(unsigned long *seed) {
    *seed = (*seed * 279470273u) % 4294967291u;
    return (FLOAT_TYPE)*seed / 4294967291.;
}

/**
 * Get a uniform random number on GPU
 * @param seed seed for the random number (initialise with time(NULL))
 * @return uniform random number
 */
__device__ FLOAT_TYPE getRandomDev(unsigned long seed) {
    seed = (seed * 279470273u) % 4294967291u;
    return (FLOAT_TYPE)seed / 4294967291.;
}

int *createGpuArrayInt(int length, ArrayOption op, int fill, int *copy)
{
    int *array_d;
    cudaError err = cudaMalloc((void**)&array_d, SIZEINT(length));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Error allocating on GPU: %d\n", err);
        return NULL;
    }
    dim3 bpg((length-1)/THREADS+1);
    cudaMemcpyToSymbol(fill_id, &fill, SIZEINT(1));
    switch (op)
    {
        case ARRAY_ZERO:
            cudaMemset(array_d, 0, SIZEINT(length));
        break;
        case ARRAY_FILL:
            gpuArrayFillInt<<<bpg,THREADS>>>(array_d, length);
        break;
        case ARRAY_COPY:
            cudaMemcpy(array_d, copy, SIZEINT(length), cudaMemcpyHostToDevice);
        break;
        case ARRAY_CPYD:
            cudaMemcpy(array_d, copy, SIZEINT(length), cudaMemcpyDeviceToDevice);
        break;
    }
#ifdef DEBUG
    printf("cm - %p (%ldB)\n", array_d, SIZEINT(length));
#endif
    return array_d;
}

FLOAT_TYPE *createGpuArrayFlt(int length, ArrayOption op, FLOAT_TYPE fill, FLOAT_TYPE *copy)
{
    FLOAT_TYPE *array_d;
    cudaError err = cudaMalloc((void**)&array_d, SIZEFLT(length));
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Error allocating on GPU: %d\n", err);
        return NULL;
    }
    unsigned long seed = time(NULL);
    dim3 bpg((length-1)/THREADS+1);
    cudaMemcpyToSymbol(fill_fd, &fill, SIZEFLT(1));
    switch (op)
    {
        case ARRAY_ZERO:
            cudaMemset(array_d, 0, SIZEFLT(length));
        break;
        case ARRAY_FILL:
            gpuArrayFillFlt<<<bpg,THREADS>>>(array_d, length);
        break;
        case ARRAY_COPY:
            cudaMemcpy(array_d, copy, SIZEFLT(length), cudaMemcpyHostToDevice);
        break;
        case ARRAY_CPYD:
            cudaMemcpy(array_d, copy, SIZEFLT(length), cudaMemcpyDeviceToDevice);
        break;
        case ARRAY_RAND:
            gpuArrayFillRandom<<<bpg,THREADS>>>(array_d, seed, length);
        break;
    }
#ifdef DEBUG
    printf("cm - %p (%ldB)\n", array_d, SIZEFLT(length));
#endif
    return array_d;
}


int *createHostArrayInt(int length, ArrayOption op, int fill, int *copy)
{
    int *array_h;
    array_h = (int*)malloc(SIZEINT(length));
    switch (op)
    {
        case ARRAY_ZERO:
            memset(array_h, 0, SIZEINT(length));
        break;
        case ARRAY_FILL:
            hostArrayFillInt(array_h, fill, length);
        break;
        case ARRAY_COPY:
            memcpy(array_h, copy, SIZEINT(length));
        break;
    }
    return array_h;
}

FLOAT_TYPE *createHostArrayFlt(int length, ArrayOption op, FLOAT_TYPE fill, FLOAT_TYPE *copy)
{
    FLOAT_TYPE *array_h;
    array_h = (FLOAT_TYPE*)malloc(SIZEFLT(length));
    switch (op)
    {
        case ARRAY_ZERO:
            memset(array_h, 0, SIZEFLT(length));
        break;
        case ARRAY_FILL:
            hostArrayFillFlt(array_h, fill, length);
            // printf("fill: %f\n", array_h[0]);
        break;
        case ARRAY_COPY:
            memcpy(array_h, copy, SIZEFLT(length));
        break;
        case ARRAY_RAND:
            hostArrayFillRandom(array_h, length, (fill)?fill:1.0);
        break;
    }
    return array_h;
}

int **create2DHostArrayInt(int width, int height)
{
    int **MyMatrix;
    int i;
    MyMatrix = (int **)calloc(height,sizeof(int*));
    assert(MyMatrix != NULL);
    for (i = 0; i < height; i++)
        MyMatrix[i] = (int *)calloc(width,sizeof(int));
    assert(MyMatrix != NULL);
    return MyMatrix;
}

FLOAT_TYPE **create2DHostArrayFlt(int width, int height)
{
    FLOAT_TYPE **MyMatrix;
    int i;
    MyMatrix = (FLOAT_TYPE **)calloc(height,sizeof(FLOAT_TYPE*));
    assert(MyMatrix != NULL);
    for (i = 0; i < height; i++)
        MyMatrix[i] = (FLOAT_TYPE *)calloc(width,sizeof(FLOAT_TYPE));
    assert(MyMatrix != NULL);
    return MyMatrix;
}

int ***create3DHostArrayInt(int width, int height, int depth)
{
    int ***MyMatrix;
    int i, j, k;

    MyMatrix = (int ***)calloc(height,sizeof(int**));
    assert(MyMatrix != NULL);
    for (i = 0; i < height; i++)
    {
        MyMatrix[i] = (int **)calloc(width,sizeof(int*));
        assert(MyMatrix != NULL);
        for (j = 0; j < width; j++)
        {
            MyMatrix[i][j] = (int *)calloc(depth,sizeof(int));
            assert(MyMatrix != NULL);
            for (k = 0; k < depth; k++)
                MyMatrix[i][j][k] = 0;
        }
    }
    return MyMatrix;
}


FLOAT_TYPE ***create3DHostArrayFlt(int width, int height, int depth)
{
    FLOAT_TYPE ***MyMatrix;
    int i, j, k;

    MyMatrix = (FLOAT_TYPE ***)calloc(height,sizeof(FLOAT_TYPE**));
    assert(MyMatrix != NULL);
    for (i = 0; i < height; i++)
    {
        MyMatrix[i] = (FLOAT_TYPE **)calloc(width,sizeof(FLOAT_TYPE*));
        assert(MyMatrix != NULL);
        for (j = 0; j < width; j++)
        {
            MyMatrix[i][j] = (FLOAT_TYPE *)calloc(depth,sizeof(FLOAT_TYPE));
            assert(MyMatrix != NULL);
            for (k = 0; k < depth; k++)
                MyMatrix[i][j][k] = 0;
        }
    }
    return MyMatrix;
}

bool ***create3DHostArrayBool(int width, int height, int depth)
{
    bool ***MyMatrix;
    int i, j, k;

    MyMatrix = (bool ***)calloc(height,sizeof(bool**));
    assert(MyMatrix != NULL);
    for (i = 0; i < height; i++)
    {
        MyMatrix[i] = (bool **)calloc(width,sizeof(bool*));
        assert(MyMatrix != NULL);
        for (j = 0; j < width; j++)
        {
            MyMatrix[i][j] = (bool *)calloc(depth,sizeof(bool));
            assert(MyMatrix != NULL);
            for (k = 0; k < depth; k++)
                MyMatrix[i][j][k] = 0;
        }
    }
    return MyMatrix;
}

__global__ void gpuArrayFillInt(int *array_d, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i<size) array_d[i] = fill_id;
}

__global__ void gpuArrayFillFlt(FLOAT_TYPE *array_d, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i<size) array_d[i] = fill_fd;
}

__global__ void gpuArrayFillRandom(FLOAT_TYPE *array_d, unsigned long seed, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i<size) array_d[i] = getRandomDev(seed+i);
}

void hostArrayFillInt(int *array_h, int fill, int size)
{
    int i;
    for (i=0; i<size; ++i) array_h[i] = fill;
}

void hostArrayFillFlt(FLOAT_TYPE *array_h, FLOAT_TYPE fill, int size)
{
    int i;
    // printf("fill_in:%f -> %f\n", fill, array_h[0]);
    for (i=0; i<size; ++i) array_h[i] = fill;
}

void hostArrayFillRandom(FLOAT_TYPE *array_h, int size, FLOAT_TYPE r)
{
    int i;
    unsigned long seed = time(NULL);
    for (i=0; i<size; ++i) array_h[i] = r*getRandom(&seed);
}

void freeAllHost(void **as, int n)
{
    int i;
    for (i=0; i<n; ++i)
    {
        free(as[i]);
        as[i] = NULL;
    }
    as = NULL;
}

void freeAllGpu(void **as, int n)
{
    int i;
    for (i=0; i<n; ++i)
    {
        cudaFree(as[i]);
        as[i] = NULL;
    }
    as = NULL;
}
