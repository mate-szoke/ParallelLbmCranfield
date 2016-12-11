/**
 * Unit tests for the residual computations
 * @file TestComputeResiduals.cu
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#include "CuTest.h"
#include "GpuFunctions.h"
#include "ComputeResiduals.h"
#include "FloatType.h"
#include "ArrayUtils.h"
#include "TestUtils.h"
#include "GpuSum.h"

/**
 * @brief Unittest for #GpuComputeResiduals
 *
 * @param tc test case
 * @test
 *   - Prepare arrays
 *   - call function
 *   - check the sum of the arrays against predefined values
 */
void testGpuComputeResiduals(CuTest *tc)
{
    printBanner("Test GpuComputeResiduals");
    N = 211;
    M = 259;
    int cond = 1;
    FLOAT_TYPE a=5.0, b=2.0, c=1.0, d=3.0;

    AIH = createHostArrayInt(  M*N); //boundary id

    int i;
    for (i=0; i<N*M; ++i)
    {
        AIH[i] = (i<N) ? cond : 0;
    }

    AFD = createGpuArrayFlt(9*M*N, ARRAY_FILL, a);
    BFD = createGpuArrayFlt(9*M*N, ARRAY_FILL, b);
    CFD = createGpuArrayFlt(  M*N, ARRAY_FILL, c);
    DFD = createGpuArrayFlt(  M*N, ARRAY_FILL, d);
    EFD = createGpuArrayFlt(9*M*N, ARRAY_ZERO);
    FFD = createGpuArrayFlt(9*M*N, ARRAY_ZERO);
    GFD = createGpuArrayFlt(  M*N, ARRAY_ZERO);
    HFD = createGpuArrayFlt(  M*N, ARRAY_ZERO);
    AID = createGpuArrayInt(  M*N, ARRAY_COPY, 0, AIH);

    FLOAT_TYPE result[4];
    for (i=0; i<4; ++i)
    {
        result[i] = 0.0;
    }

    GpuComputeResiduals(AID, AFD, BFD, CFD, DFD, EFD, FFD, GFD, HFD, result, &M, &N, cond);

    CuAssertDblEquals(tc, sqrt((a-b)*(a-b)*N*M*9), result[0], 0.0001);
    CuAssertDblEquals(tc, sqrt((a-b)*(a-b)*N*M*9/M/N), result[1], 0.0001);
    CuAssertDblEquals(tc, c*N, result[2], 0.0001);
    CuAssertDblEquals(tc, d*N, result[3], 0.0001);
}

///Clean up after test case
void cleanupTestGpuComputeResiduals()
{

    cudaFree(AFD); cudaFree(BFD); cudaFree(CFD); cudaFree(DFD); cudaFree(EFD); cudaFree(FFD);
    cudaFree(GFD); cudaFree(HFD); cudaFree(AID);
    free(AIH);
}

/**
 * @brief Unittest for #GpuComputeResidMask
 *
 * @param tc test case
 * @test
 *   - Prepare arrays
 *   - call function
 *   - check the sum of the arrays against predefined values
 */
void testGpuComputeResidMask(CuTest *tc)
{
    printBanner("Test gpuComputeResidMask");
    N = 257;
    M = 215;
    FLOAT_TYPE a=5.0, b=2.0, c=1.0, d=3.0;
    int boundId = 2;

    AIH = createHostArrayInt(M*N);               //bcmask
    createLidBcMaskFull(AIH, M, N);
    AID = createGpuArrayInt(M*N, ARRAY_COPY, 0, AIH);

    AFD = createGpuArrayFlt(9*M*N, ARRAY_FILL, a); //f
    BFD = createGpuArrayFlt(9*M*N, ARRAY_FILL, b); //fTemp
    CFD = createGpuArrayFlt(  M*N, ARRAY_FILL, c); //drag
    DFD = createGpuArrayFlt(  M*N, ARRAY_FILL, d); //lift
    EFD = createGpuArrayFlt(9*M*N, ARRAY_ZERO);    //temp
    FFD = createGpuArrayFlt(9*M*N, ARRAY_ZERO);    //temp
    GFD = createGpuArrayFlt(  M*N, ARRAY_ZERO);    //temp
    HFD = createGpuArrayFlt(  M*N, ARRAY_ZERO);    //temp

    FLOAT_TYPE result[4];
    int i;
    for (i=0; i<4; ++i)
    {
        result[i] = 0.0;
    }

    GpuComputeResidMask(AID, AFD, BFD, CFD, DFD, EFD, FFD, GFD, HFD, result, &M, &N, boundId);

    CuAssertDblEquals(tc, sqrt((a-b)*(a-b)*N*M*9), result[0], 0.0001);
    CuAssertDblEquals(tc, sqrt((a-b)*(a-b)*N*M*9/M/N), result[1], 0.0001);
    CuAssertDblEquals(tc, c*M, result[2], 0.0001);
    CuAssertDblEquals(tc, d*M, result[3], 0.0001);
}

///Clean up after test case
void cleanupTestGpuComputeResidMask()
{
    cudaFree(AFD); cudaFree(BFD); cudaFree(CFD); cudaFree(DFD); cudaFree(EFD); cudaFree(FFD);
    cudaFree(GFD); cudaFree(HFD); cudaFree(AID);
    free(AIH);
}

/**
 * @brief Unittest for #computeResidual
 *
 * @param tc test case
 * @test
 *   - Prepare arrays
 *   - call function
 *   - check the sum of the arrays against predefined values
 */
void testComputeResidual(CuTest *tc)
{
    printBanner("Test computeResidual");
    N = 457;
    M = 391;
    FLOAT_TYPE a=5.0, b=2.0;

    AFD = createGpuArrayFlt(9*M*N, ARRAY_FILL, a);
    BFD = createGpuArrayFlt(9*M*N, ARRAY_FILL, b);
    EFD = createGpuArrayFlt(9*M*N);
    FFD = createGpuArrayFlt(9*M*N);

    FLOAT_TYPE result = computeResidual(AFD, BFD, EFD, FFD, M, N);
    CuAssertDblEquals(tc, sqrt((a-b)*(a-b)*N*M*9), result, 0.0001);
}

///Clean up after test case
void cleanupTestComputeResidual()
{
    cudaFree(AFD); cudaFree(BFD); cudaFree(EFD); cudaFree(FFD);
}

/**
 * @brief Unittest for #computeDragLift
 *
 * @param tc test case
 * @test
 *   - Prepare arrays
 *   - call function
 *   - check the sum of the arrays against predefined values
 */
void testComputeDragLift(CuTest *tc)
{
    printBanner("Test computeDragLift");
    N = 3333;
    M = 3234;
    FLOAT_TYPE a=5.0;
    int boundId = 2;

    AFD = createGpuArrayFlt(M*N, ARRAY_FILL, a); //drag or lift
    BFD = createGpuArrayFlt(M*N);                //temp
    CFD = createGpuArrayFlt(M*N);                //temp

    AIH = createHostArrayInt(M*N);               //bcmask
    createLidBcMaskFull(AIH, M, N);
    AID = createGpuArrayInt(M*N, ARRAY_COPY, 0, AIH);

    FLOAT_TYPE result = computeDragLift(AID, AFD, BFD, CFD, M, N, boundId);
    CuAssertDblEquals(tc, a*M, result, 0.0001);
}

///Clean up after test case
void cleanupTestComputeDragLift()
{
    cudaFree(AFD); cudaFree(BFD); cudaFree(CFD); cudaFree(AID);
    free(AIH);
}

CuSuite* computeResidualsGetSuite()
{
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TCLN(suite, testGpuComputeResiduals, cleanupTestGpuComputeResiduals);
    SUITE_ADD_TCLN(suite, testGpuComputeResidMask, cleanupTestGpuComputeResidMask);
    SUITE_ADD_TCLN(suite, testComputeResidual,     cleanupTestComputeResidual);
    SUITE_ADD_TCLN(suite, testComputeDragLift,     cleanupTestComputeDragLift);
    return suite;
}
