/**
 * Unittests for GPU sum functions
 * @file TestGpuSum.cu
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#include "FloatType.h"
#include "CuTest.h"
#include "GpuSum.h"
#include "ArrayUtils.h"
#include "TestUtils.h"

/**
 * @brief Unittest for #gpu_cond_copy
 *
 * @param tc test case
 * @test
 *   - Allocate arrays and fill them with random values
 *   - call function
 *   - check the sum of the arrays against predefined values
 */
void testGpuCondCopy(CuTest *tc)
{
    printBanner("Test gpu_cond_copy");
    unsigned long seed = time(NULL);
    N = 265;
    M = 222;
    int cond = 1;

    AFH = createHostArrayFlt(M*N);
    BFH = createHostArrayFlt(M*N);
    AIH = createHostArrayInt(M*N);

    int i;
    for (i=0; i<N*M; ++i)
    {
        AFH[i] = getRandom(&seed);
        AIH[i] = (i<N) ? cond : 0;
    }

    AFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(M*N);
    AID = createGpuArrayInt(M*N, ARRAY_COPY, 0, AIH);

    gpu_cond_copy <<< (N*M-1)/THREADS+1, THREADS >>> (BFD, AFD, AID, 1, N*M);

    cudaMemcpy(BFH, BFD, SIZEFLT(M*N), cudaMemcpyDeviceToHost);

    int b = 1;
    for (i=0; i<N*M; ++i)
    {
        b &= (BFH[i] == AFH[i] && i<N) || (BFH[i] == 0 && i>=N);
    }

    CuAssertIntEquals(tc, cond, b);
}

///Clean up after test case
void cleanupTestGpuCondCopy()
{
    cudaFree(AFD); cudaFree(BFD); cudaFree(AID);
    free(AFH); free(BFH); free(AIH);
}

/**
 * @brief Unittest for #gpu_sqsub
 *
 * @param tc test case
 * @test
 *   - Allocate arrays and fill them with fixed values
 *   - call function
 *   - check the sum of the arrays against predefined values
 */
void testGpuSquareSubstract(CuTest *tc)
{
    printBanner("Test gpu_sqsub");
    N = 211;
    M = 259;

    AFH = createHostArrayFlt(M*N);
    BFH = createHostArrayFlt(M*N);
    CFH = createHostArrayFlt(M*N);

    int i;
    for (i=0; i<N*M; ++i)
    {
        AFH[i] = 5.0;
        BFH[i] = 2.0;
    }

    AFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, BFH);
    CFD = createGpuArrayFlt(M*N);

    gpu_sqsub <<< (N*M-1)/THREADS+1, THREADS >>> (AFD, BFD, CFD, N*M);

    cudaMemcpy(CFH, CFD, SIZEFLT(M*N), cudaMemcpyDeviceToHost);

    int b = 1;
    for (i=0; i<N*M; ++i)
    {
        b &= CFH[i] == 9.0;
    }

    CuAssertIntEquals(tc, 1, b);
}

///Clean up after test case
void cleanupTestGpuSquareSubsctract()
{
    cudaFree(AFD); cudaFree(BFD); cudaFree(CFD);
    free(AFH); free(BFH); free(CFH);
}

/**
 * @brief Unittest for #gpu_sum
 *
 * @param tc test case
 * @test
 *   - Allocate arrays and fill them with fixed values
 *   - call function
 *   - check the sum of the arrays against predefined values
 */
void testGpuSum(CuTest *tc)
{
    printBanner("Test gpu_sum");
    dim3 grid_dim;
    FLOAT_TYPE result = 0.0;
    FLOAT_TYPE num = 1.0;
    N = 2307;
    M = 255;

    AFH = createHostArrayFlt(M*N);
    BFH = createHostArrayFlt(M*N);

    int i;
    for (i=0; i<N*M; ++i)
    {
        AFH[i] = num;
    }
    AFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(M*N);

    int remaining = N*M;
    int shared_size = THREADS * sizeof(FLOAT_TYPE);
    int req_blocks = 0;
    while (remaining > 1)
    {
        req_blocks = (remaining - 1) / THREADS / 2 + 1;
        grid_dim.x = static_cast<int>(ceil(sqrt(req_blocks)));
        grid_dim.y = (req_blocks - 1) / grid_dim.x + 1;
        gpu_sum <<<grid_dim, THREADS, shared_size >>>(AFD, BFD, remaining);

        //swap
        FLOAT_TYPE *temp = AFD;
        AFD = BFD;
        BFD = temp;

        remaining = req_blocks;
        cudaMemcpy(&result, AFD, sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);
    }

    cudaMemcpy(&result, AFD, sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);

    CuAssertDblEquals(tc, N*M*num, result, 0.0001);
}

/**
 * @brief Unittest for #gpu_sum256
 *
 * @param tc test case
 * @test
 *   - Allocate arrays and fill them with fixed values
 *   - call function
 *   - check the sum of the arrays against predefined values
 */
void testGpuSum256(CuTest *tc)
{
    printBanner("Test gpu_sum256");
    dim3 grid_dim;
    FLOAT_TYPE result = 0.0;
    FLOAT_TYPE num = 1.0;
    N = 2307;
    M = 255;

    AFH = createHostArrayFlt(M*N);
    BFH = createHostArrayFlt(M*N);

    int i;
    for (i=0; i<N*M; ++i)
    {
        AFH[i] = num;
    }
    AFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(M*N);

    int remaining = N*M;
    int shared_size = THREADS * sizeof(FLOAT_TYPE);
    int req_blocks = 0;
    while (remaining > 1)
    {
        req_blocks = (remaining - 1) / THREADS / 2 + 1;
        grid_dim.x = static_cast<int>(ceil(sqrt(req_blocks)));
        grid_dim.y = (req_blocks - 1) / grid_dim.x + 1;
        gpu_sum256 <<<grid_dim, THREADS, shared_size >>>(AFD, BFD, remaining);

        //swap
        FLOAT_TYPE *temp = AFD;
        AFD = BFD;
        BFD = temp;

        remaining = req_blocks;
        cudaMemcpy(&result, AFD, sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);
    }

    cudaMemcpy(&result, AFD, sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);

    CuAssertDblEquals(tc, N*M*num, result, 0.0001);
}

/**
 * @brief Unittest for #gpu_sum_h
 *
 * @param tc test case
 * @test
 *   - Allocate arrays and fill them with fixed values
 *   - call function
 *   - check the sum of the arrays against predefined values
 */
void testGpusumHost(CuTest *tc)
{
    printBanner("Test gpu_sum_h");
    FLOAT_TYPE result = 0.0;
    FLOAT_TYPE num = 1.0;
    N = 2907;
    M = 242;

    AFH = createHostArrayFlt(M*N);
    BFH = createHostArrayFlt(M*N);

    int i;
    for (i=0; i<N*M; ++i)
    {
        AFH[i] = num;
    }
    AFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(M*N);

    result = gpu_sum_h(AFD, BFD, N*M);

    CuAssertDblEquals(tc, N*M*num, result, 0.0001);
}

///Clean up after test case
void cleanTestGpuSum()
{

    cudaFree(AFD); cudaFree(BFD);
    free(AFH); free(BFH);
}

CuSuite* gpuSumGetSuite()
{
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TCLN(suite, testGpuCondCopy,         cleanupTestGpuCondCopy);
    SUITE_ADD_TCLN(suite, testGpuSquareSubstract,  cleanupTestGpuSquareSubsctract);
    SUITE_ADD_TCLN(suite, testGpuSum,              cleanTestGpuSum);
    SUITE_ADD_TCLN(suite, testGpuSum256,           cleanTestGpuSum);
    SUITE_ADD_TCLN(suite, testGpusumHost,          cleanTestGpuSum);
    return suite;
}