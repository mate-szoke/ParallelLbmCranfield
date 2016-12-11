/**
 * Unittests for the boundary conditions
 * @file TestGpuBoundaries.cu
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#include "CuTest.h"
#include "GpuFunctions.h"
#include "ShellFunctions.h"
#include "ComputeResiduals.h"
#include "FloatType.h"
#include "TestUtils.h"
#include "ArrayUtils.h"
#include "GpuConstants.h"

static FLOAT_TYPE fSum[2];

/**
 * @brief Initialise conditions for BC testcases
 *
 * @param af value to fill the distribution function
 * @param bf value to fill the u0 array
 */
void BoundaryTestInit(FLOAT_TYPE af, FLOAT_TYPE bf)
{
    AIH = createHostArrayInt(M*N, ARRAY_FILL, 0); //fluid
    BIH = createHostArrayInt(M*N, ARRAY_FILL, 0); //boundary
    CIH = createHostArrayInt(9*M*N, ARRAY_FILL, 0); //bcid
    DIH = createHostArrayInt(M*N, ARRAY_FILL, 0); //corner

    createLidBcFluid(AIH, M, N);
    createLidBcBoundary(BIH, M, N);
    createLidBcBcId(CIH, M, N);
    createLidBcCorner(DIH, M, N);

    AID = createGpuArrayInt(M*N, ARRAY_COPY, 0, AIH);
    BID = createGpuArrayInt(M*N, ARRAY_COPY, 0, BIH);
    CID = createGpuArrayInt(9*M*N, ARRAY_COPY, 0, CIH);
    DID = createGpuArrayInt(M*N, ARRAY_COPY, 0, DIH);

    AFH = createHostArrayFlt(9*M*N, ARRAY_FILL, af); //f
    BFH = createHostArrayFlt(M*N, ARRAY_FILL, bf); //uo
    CFH = createHostArrayFlt(M*N, ARRAY_FILL, 0.); //vo

    AFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, BFH);
    CFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, CFH);
}

/**
 * @brief Test #gpu_boundaries1
 *
 * @param tc test case
 * @test
 *  - Prepare boundary conditions for the lid driven cavity
 *  - Run the function
 *  - Check the sum of all the arrays against predefined values
 */
void testGpuBoundaries1(CuTest *tc)
{
    printBanner("Test gpu_boundaries1");
    M = 64;
    N = 64;
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    dim3 tpb(THREADS);
    dim3 bpg9((int)(9*M*N/THREADS)+1);

    FLOAT_TYPE af = 0.002;
    FLOAT_TYPE bf = 0.3;

    BoundaryTestInit(af, bf);

    int sumAI = sumHostInt(AIH, M*N);
    int sumBI = sumHostInt(BIH, M*N);
    int sumCI = sumHostInt(CIH, 9*M*N);
    int sumDI = sumHostInt(DIH, M*N);
    FLOAT_TYPE sumAF = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBF = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCF = sumHostFlt(CFH, M*N);

    CuAssertIntEquals(tc, N*M, sumAI);
    CuAssertIntEquals(tc, 3*M+2*N-4, sumBI);
    CuAssertIntEquals(tc, N*6+M*9-6, sumCI);
    CuAssertIntEquals(tc, 2, sumDI);
    FLOAT_TYPE tmpaf = sumAF;
    FLOAT_TYPE tmpbf = sumBF;
    CuAssertDblEquals(tc, 0.0, sumCF, 0.01);

    gpu_boundaries1<<<bpg9,tpb>>>(AID, BID, CID, AFD, BFD, CFD, DID);

    cudaMemcpy(AIH, AID, SIZEINT(N*M),   cudaMemcpyDeviceToHost);
    cudaMemcpy(BIH, BID, SIZEINT(N*M),   cudaMemcpyDeviceToHost);
    cudaMemcpy(CIH, CID, SIZEINT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(DIH, DID, SIZEINT(N*M),   cudaMemcpyDeviceToHost);
    cudaMemcpy(AFH, AFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, SIZEFLT(N*M),   cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, SIZEFLT(N*M),   cudaMemcpyDeviceToHost);

    sumAI = sumHostInt(AIH, M*N);
    sumBI = sumHostInt(BIH, M*N);
    sumCI = sumHostInt(CIH, 9*M*N);
    sumDI = sumHostInt(DIH, M*N);
    sumAF = sumHostFlt(AFH, 9*M*N);
    sumBF = sumHostFlt(BFH, M*N);
    sumCF = sumHostFlt(CFH, M*N);

    CuAssertIntEquals(tc, N*M, sumAI);
    CuAssertIntEquals(tc, 3*M+2*N-4, sumBI);
    CuAssertIntEquals(tc, N*6+M*9-6, sumCI);
    CuAssertIntEquals(tc, 2, sumDI);
    CuAssertDblEquals(tc, tmpaf, sumAF, 0.00001);
    CuAssertDblEquals(tc, tmpbf, sumBF, 0.00001);
    CuAssertDblEquals(tc, 0.0, sumCF, 0.01);
}

/**
 * @brief Test #gpu_boundaries2
 *
 * @param tc test case
 * @param boundaryType straight or curved
 * @test
 *  - Prepare boundary conditions for the lid driven cavity
 *  - Run the function
 *  - Check the sum of all the arrays against predefined values
 */
void runTestGpuBoundaries2(CuTest *tc, BoundaryType boundaryType)
{
    M = 64;
    N = 64;
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    cudaMemcpyToSymbol(boundaryType_d, &boundaryType, sizeof(int));
    dim3 tpb(THREADS);
    dim3 bpg9((int)(9*M*N/THREADS)+1);

    FLOAT_TYPE bf = 0.002;
    FLOAT_TYPE cf = 0.5;
    FLOAT_TYPE df = 0.003;

    AIH = createHostArrayInt(M*N, ARRAY_FILL, 0);
    BIH = createHostArrayInt(9*M*N, ARRAY_FILL, 0);

    createLidBcFluid(AIH, M, N);
    createLidBcBcId(BIH, M, N);

    AID = createGpuArrayInt(M*N, ARRAY_COPY, 0, AIH);
    BID = createGpuArrayInt(9*M*N, ARRAY_COPY, 0, BIH);

    AFH = createHostArrayFlt(9*M*N, ARRAY_FILL, 0.);
    BFH = createHostArrayFlt(9*M*N, ARRAY_FILL, bf);
    CFH = createHostArrayFlt(9*M*N, ARRAY_FILL, cf);
    DFH = createHostArrayFlt(9*M*N, ARRAY_FILL, df);

    AFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, BFH);
    CFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, CFH);
    DFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, DFH);

    int sumAIa = sumHostInt(AIH, M*N);
    int sumBIa = sumHostInt(BIH, 9*M*N);
    FLOAT_TYPE sumAFa = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFa = sumHostFlt(BFH, 9*M*N);
    FLOAT_TYPE sumCFa = sumHostFlt(CFH, 9*M*N);
    FLOAT_TYPE sumDFa = sumHostFlt(DFH, 9*M*N);

    gpu_boundaries2<<<bpg9,tpb>>>(AID, AFD, BFD, BID, CFD, DFD);

    cudaMemcpy(AIH, AID, SIZEINT(M*N),   cudaMemcpyDeviceToHost);
    cudaMemcpy(BIH, BID, SIZEINT(9*M*N), cudaMemcpyDeviceToHost);
    cudaMemcpy(AFH, AFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(DFH, DFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);

    int sumAIb = sumHostInt(AIH, M*N);
    int sumBIb = sumHostInt(BIH, 9*M*N);
    FLOAT_TYPE sumAFb = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFb = sumHostFlt(BFH, 9*M*N);
    FLOAT_TYPE sumCFb = sumHostFlt(CFH, 9*M*N);
    FLOAT_TYPE sumDFb = sumHostFlt(DFH, 9*M*N);

    CuAssertIntEquals(tc, sumAIa, sumAIb);
    CuAssertIntEquals(tc, sumBIa, sumBIb);
    CuAssertDblEquals(tc, sumBFa, sumBFb, 0.00001);
    CuAssertDblEquals(tc, sumCFa, sumCFb, 0.00001);
    if (boundaryType == CURVED)
    {
        // printf(" f was: %f is: %f\n", sumDFa, sumDFb);
        fSum[0] = sumDFb;
        CuAssertTrue(tc, sumDFa > sumDFb);
        CuAssertTrue(tc, sumAFa < sumAFb);
    }
    else if (boundaryType == 2)
    {
        CuAssertDblEquals(tc, sumAFa, sumAFb, 0.00001);
        CuAssertDblEquals(tc, sumDFa, sumDFb, 0.00001);
    }
}

///Clean up after test case
void cleanupTestGpuBoundaries2()
{
    cudaFree(AID); cudaFree(BID); cudaFree(AFD); cudaFree(BFD); cudaFree(CFD); cudaFree(DFD);
    free(AIH); free(BIH); free(AFH); free(BFH); free(CFH); free(DFH);
}

///Test case for curved wall @param tc test case
void testGpuBoundaries2_wcb(CuTest *tc)
{
    printBanner("Test gpu_boundaries2 curved");
    runTestGpuBoundaries2(tc, CURVED);
}

///Test case for straight wall @param tc test case
void testGpuBoundaries2_wocb(CuTest *tc)
{
    printBanner("Test gpu_boundaries2 straight");
    runTestGpuBoundaries2(tc, STRAIGHT);
}

/**
 * @brief Test #gpu_boundaries3
 *
 * @param tc test case
 * @test
 *  - Prepare boundary conditions for the lid driven cavity
 *  - Run the function
 *  - Check the sum of all the arrays against predefined values
 */
void testGpuBoundaries3(CuTest *tc)
{
    printBanner("Test gpu_boundaries3");
    M = 64;
    N = 64;
    OutletProfile outl = OUTLET_SECOND;
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    cudaMemcpyToSymbol(outletProfile_d, &outl, sizeof(int));
    dim3 tpb(THREADS);
    dim3 bpg9((int)(9*M*N/THREADS)+1);

    FLOAT_TYPE af = 0.002;
    FLOAT_TYPE bf = 0.3;

    BoundaryTestInit(af, bf);

    int sumAIa = sumHostInt(AIH, M*N);   //fluid
    int sumCIa = sumHostInt(CIH, 9*M*N); //bcid
    int sumDIa = sumHostInt(DIH, M*N);   //corner
    FLOAT_TYPE sumAFa = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFa = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCFa = sumHostFlt(CFH, M*N);

    gpu_boundaries3<<<bpg9,tpb>>>(AID, CID, AFD, BFD, CFD, DID);

    cudaMemcpy(AIH, AID, SIZEINT(M*N),   cudaMemcpyDeviceToHost);
    cudaMemcpy(CIH, CID, SIZEINT(9*M*N), cudaMemcpyDeviceToHost);
    cudaMemcpy(AFH, AFD, SIZEFLT(9*M*N), cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, SIZEFLT(M*N),   cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, SIZEFLT(M*N),   cudaMemcpyDeviceToHost);

    int sumAIb = sumHostInt(AIH, M*N);
    int sumCIb = sumHostInt(CIH, 9*M*N);
    int sumDIb = sumHostInt(DIH, M*N);
    FLOAT_TYPE sumAFb = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFb = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCFb = sumHostFlt(CFH, M*N);

    CuAssertIntEquals(tc, sumAIa, sumAIb);
    CuAssertIntEquals(tc, sumCIa, sumCIb);
    CuAssertIntEquals(tc, sumDIa, sumDIb);
    CuAssertDblEquals(tc, sumAFa, sumAFb, 0.00001);
    CuAssertDblEquals(tc, sumBFa, sumBFb, 0.00001);
    CuAssertDblEquals(tc, sumCFa, sumCFb, 0.00001);
}

///Clean up after test case
void cleanupTestGpuBoundaries()
{
    cudaFree(AID); cudaFree(BID); cudaFree(CID); cudaFree(DID); cudaFree(AFD); cudaFree(BFD); cudaFree(CFD);
    free(AIH); free(BIH); free(CIH); free(DIH); free(AFH); free(BFH); free(CFH);
}

/**
 * @brief Test #gpuBcInlet
 *
 * @param tc test case
 * @test
 *  - Prepare boundary conditions for the lid driven cavity
 *  - Run the function
 *  - Check the sum of all the arrays against predefined values
 */
void testGpuBcInlet(CuTest *tc)
{
    printBanner("Test gpuBcInlet");
    M = 64;
    N = 64;
    int size = getLidBcMaskSize(M,N);
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    dim3 tpb(THREADS);
    dim3 bpgB((int)(size/THREADS)+1);

    FLOAT_TYPE af = 0.002;
    FLOAT_TYPE bf = 0.3;

    AIH = createHostArrayInt(size, ARRAY_FILL, 0);
    BIH = createHostArrayInt(size, ARRAY_FILL, 0);

    createLidBcIdx(AIH, M, N);
    createLidBcMask(BIH, M, N);

    AID = createGpuArrayInt(size, ARRAY_COPY, 0, AIH);
    BID = createGpuArrayInt(size, ARRAY_COPY, 0, BIH);

    AFH = createHostArrayFlt(9*M*N, ARRAY_FILL, af);
    BFH = createHostArrayFlt(M*N, ARRAY_FILL, bf);
    CFH = createHostArrayFlt(M*N, ARRAY_FILL, 0.);

    AFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, BFH);
    CFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, CFH);

    int sumAIa = sumHostInt(AIH, size);
    int sumBIa = sumHostInt(BIH, size);
    FLOAT_TYPE sumAFa = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFa = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCFa = sumHostFlt(CFH, M*N);

    gpuBcInlet<<<bpgB,tpb>>>(AID, BID, AFD, BFD, CFD, size);

    cudaMemcpy(AIH, AID, SIZEINT(size),  cudaMemcpyDeviceToHost);
    cudaMemcpy(BIH, BID, SIZEINT(size),  cudaMemcpyDeviceToHost);
    cudaMemcpy(AFH, AFD, SIZEFLT(9*M*N), cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, SIZEFLT(M*N),   cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, SIZEFLT(M*N),   cudaMemcpyDeviceToHost);

    int sumAIb = sumHostInt(AIH, size);
    int sumBIb = sumHostInt(BIH, size);
    FLOAT_TYPE sumAFb = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFb = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCFb = sumHostFlt(CFH, M*N);

    CuAssertIntEquals(tc, sumAIa, sumAIb);
    CuAssertIntEquals(tc, sumBIa, sumBIb);
    CuAssertDblEquals(tc, sumAFa, sumAFb, 0.000001);
    CuAssertDblEquals(tc, sumBFa, sumBFb, 0.000001);
    CuAssertDblEquals(tc, sumCFa, sumCFb, 0.000001);
}

///Clean up after test case
void cleanupTestGpuBcInlet()
{
    cudaFree(AID); cudaFree(BID); cudaFree(AFD); cudaFree(BFD); cudaFree(CFD);
    free(AIH); free(BIH); free(AFH); free(BFH); free(CFH);
}


/**
 * @brief Test #gpuBcWall
 *
 * @param tc test case
 * @param boundaryType straight or curved
 * @test
 *  - Prepare boundary conditions for the lid driven cavity
 *  - Run the function
 *  - Check the sum of all the arrays against predefined values
 */
void runTestGpuBcWall(CuTest *tc, BoundaryType boundaryType)
{
    M = 64;
    N = 64;
    int size = getLidBcMaskSize(M, N);
    cudaMemcpyToSymbol(boundaryType_d, &boundaryType, sizeof(int));
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    dim3 tpb(THREADS);
    dim3 bpgB((int)(size/THREADS)+1);

    FLOAT_TYPE af = 0.003;
    FLOAT_TYPE bf = 0.002;

    AIH = createHostArrayInt(size, ARRAY_FILL, 0);
    BIH = createHostArrayInt(size, ARRAY_FILL, 0);

    createLidBcIdx(AIH, M, N);
    createLidBcMask(BIH, M, N);

    AID = createGpuArrayInt(size, ARRAY_COPY, 0, AIH);
    BID = createGpuArrayInt(size, ARRAY_COPY, 0, BIH);

    AFH = createHostArrayFlt(9*M*N, ARRAY_FILL, af);
    BFH = createHostArrayFlt(9*M*N, ARRAY_FILL, bf);
    CFH = createHostArrayFlt(8*size, ARRAY_FILL, 0.5);

    AFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, BFH);
    CFD = createGpuArrayFlt(8*size, ARRAY_COPY, 0, CFH);

    int sumAIa = sumHostInt(AIH, size);
    int sumBIa = sumHostInt(BIH, size);
    FLOAT_TYPE sumAFa = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFa = sumHostFlt(BFH, 9*M*N);
    FLOAT_TYPE sumCFa = sumHostFlt(CFH, 8*size);

    gpuBcWall<<<bpgB,tpb>>>(AID, BID, AFD, BFD, CFD, size);

    cudaMemcpy(AIH, AID, SIZEINT(size),   cudaMemcpyDeviceToHost);
    cudaMemcpy(BIH, BID, SIZEINT(size),   cudaMemcpyDeviceToHost);
    cudaMemcpy(AFH, AFD, SIZEFLT(9*M*N),  cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, SIZEFLT(9*M*N),    cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, SIZEFLT(8*size), cudaMemcpyDeviceToHost);

    int sumAIb = sumHostInt(AIH, size);
    int sumBIb = sumHostInt(BIH, size);
    FLOAT_TYPE sumAFb = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFb = sumHostFlt(BFH, 9*M*N);
    FLOAT_TYPE sumCFb = sumHostFlt(CFH, 8*size);

    CuAssertIntEquals(tc, sumAIa, sumAIb);
    CuAssertIntEquals(tc, sumBIa, sumBIb);
    if (boundaryType == CURVED)
    {
        // printf(" f was: %f is: %f\n", sumAFa, sumAFb);
        fSum[1] = sumAFb;
        CuAssertTrue(tc, sumAFa > sumAFb);
    }
    else if (boundaryType == 2)
    {
        CuAssertDblEquals(tc, sumAFa, sumAFb, 0.00001);
    }
    CuAssertDblEquals(tc, sumBFa, sumBFb, 0.000001);
    CuAssertDblEquals(tc, sumCFa, sumCFb, 0.000001);
}

///Clean up after test case
void cleanupTestGpuBcWall()
{
    cudaFree(AID); cudaFree(BID); cudaFree(AFD); cudaFree(BFD); cudaFree(CFD);
    free(AIH); free(BIH); free(AFH); free(BFH); free(CFH);
}

///Test case for curved wall @param tc test case
void testGpuBcWall_wcb(CuTest *tc)
{
    printBanner("Test gpuBcWall curved");
    runTestGpuBcWall(tc, CURVED);
}

///Test case for straight wall @param tc test case
void testGpuBcWall_wocb(CuTest *tc)
{
    printBanner("Test gpuBcWall straight");
    runTestGpuBcWall(tc, STRAIGHT);
}

/**
 * @brief Test #gpuBcOutlet
 *
 * @param tc test case
 * @test
 *  - Prepare boundary conditions for the lid driven cavity
 *  - Run the function
 *  - Check the sum of all the arrays against predefined values
 */

void testGpuBcOutlet(CuTest *tc)
{
    printBanner("Test gpuBcOutlet");
    M = 64;
    N = 64;
    int size = getLidBcMaskSize(M,N);
    OutletProfile outl = OUTLET_SECOND;
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    cudaMemcpyToSymbol(outletProfile_d, &outl, sizeof(int));
    dim3 tpb(THREADS);
    dim3 bpgB((int)(size/THREADS)+1);

    FLOAT_TYPE af = 0.002;
    FLOAT_TYPE bf = 0.3;

    AIH = createHostArrayInt(size, ARRAY_FILL, 0);
    BIH = createHostArrayInt(size, ARRAY_FILL, 0);

    createLidBcIdx(AIH, M, N);
    createLidBcMask(BIH, M, N);

    AID = createGpuArrayInt(size, ARRAY_COPY, 0, AIH);
    BID = createGpuArrayInt(size, ARRAY_COPY, 0, BIH);

    AFH = createHostArrayFlt(9*M*N, ARRAY_FILL, af);
    BFH = createHostArrayFlt(M*N, ARRAY_FILL, bf);
    CFH = createHostArrayFlt(M*N, ARRAY_FILL, 0.);

    AFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, BFH);
    CFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, CFH);

    int sumAIa = sumHostInt(AIH, size);
    int sumBIa = sumHostInt(BIH, size);
    FLOAT_TYPE sumAFa = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFa = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCFa = sumHostFlt(CFH, M*N);

    gpuBcOutlet<<<bpgB,tpb>>>(AID, BID, AFD, BFD, CFD, size);

    cudaMemcpy(AIH, AID, SIZEINT(size),  cudaMemcpyDeviceToHost);
    cudaMemcpy(BIH, BID, SIZEINT(size),  cudaMemcpyDeviceToHost);
    cudaMemcpy(AFH, AFD, SIZEFLT(9*M*N), cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, SIZEFLT(M*N),   cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, SIZEFLT(M*N),   cudaMemcpyDeviceToHost);

    int sumAIb = sumHostInt(AIH, size);
    int sumBIb = sumHostInt(BIH, size);
    FLOAT_TYPE sumAFb = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFb = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCFb = sumHostFlt(CFH, M*N);

    CuAssertIntEquals(tc, sumAIa, sumAIb);
    CuAssertIntEquals(tc, sumBIa, sumBIb);
    CuAssertDblEquals(tc, sumAFa, sumAFb, 0.000001);
    CuAssertDblEquals(tc, sumBFa, sumBFb, 0.000001);
    CuAssertDblEquals(tc, sumCFa, sumCFb, 0.000001);
}

///Clean up after test case
void cleanupTestGpuBcOutlet()
{
    cudaFree(AID); cudaFree(BID); cudaFree(AFD); cudaFree(BFD); cudaFree(CFD);
    free(AIH); free(BIH); free(AFH); free(BFH); free(CFH);
}

/**
 * @brief Test to compare results from #gpuBcWall and #gpu_boundaries2
 *
 * @param tc test case
 * @test
 *  - Prepare boundary conditions for the lid driven cavity
 *  - Run the function
 *  - Check the sum of all the arrays against predefined values
 */
void testCompareWall(CuTest *tc)
{
    printBanner("Test compare wall");
    CuAssertDblEquals(tc, fSum[0], fSum[1], 1e-12);
}

/**
 * @brief Test to compare results from #gpuBcOutlet and #gpu_boundaries3
 *
 * @param tc test case
 * @param op outlet profile
 * @test
 *  - Prepare boundary conditions for the lid driven cavity
 *  - Run the function
 *  - Check the sum of all the arrays between the two algorithms
 */
void runTestCompareOutlet(CuTest *tc, OutletProfile op)
{
    M = 64;
    N = 64;
    int size = getLidBcMaskSize(M,N);
    OutletProfile outl = op;
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    cudaMemcpyToSymbol(outletProfile_d, &outl, sizeof(int));
    dim3 tpb(THREADS);
    dim3 bpgB((int)(size/THREADS)+1);
    dim3 bpg9((int)(9*M*N/THREADS)+1);
    FLOAT_TYPE bf = 0.3;

    AIH = createHostArrayInt(size, ARRAY_FILL, 0);
    BIH = createHostArrayInt(size, ARRAY_FILL, 0);

    createLidBcIdx(AIH, M, N);
    createChannelBcMask(BIH, M, N);

    AID = createGpuArrayInt(size, ARRAY_COPY, 0, AIH); //bcidx
    BID = createGpuArrayInt(size, ARRAY_COPY, 0, BIH); //bcmask

    AFH = createHostArrayFlt(9*M*N, ARRAY_RAND);   //f
    BFH = createHostArrayFlt(M*N, ARRAY_FILL, bf); //fTemp
    CFH = createHostArrayFlt(M*N, ARRAY_FILL, 0.); //q

    AFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, BFH);
    CFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, CFH);

    int sumAIa = sumHostInt(AIH, size);
    int sumBIa = sumHostInt(BIH, size);
    FLOAT_TYPE sumAFa = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFa = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCFa = sumHostFlt(CFH, M*N);

    //old
    CIH = createHostArrayInt(M*N, ARRAY_FILL, 0);   //fluid
    DIH = createHostArrayInt(9*M*N, ARRAY_FILL, 0); //bcid
    EIH = createHostArrayInt(M*N, ARRAY_FILL, 0);   //corner

    createLidBcFluid(CIH, M, N);
    createChannelBcBcId(DIH, M, N);
    createChannelBcCorner(EIH, M, N);

    CID = createGpuArrayInt(M*N, ARRAY_COPY, 0, CIH);
    DID = createGpuArrayInt(9*M*N, ARRAY_COPY, 0, DIH);
    EID = createGpuArrayInt(M*N, ARRAY_COPY, 0, EIH);

    DFH = createHostArrayFlt(9*M*N, ARRAY_COPY, 0, AFH); //f
    DFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, DFH);

    FLOAT_TYPE sumDFa = sumHostFlt(DFH, 9*M*N);

    gpuBcOutlet<<<bpgB,tpb>>>(AID, BID, AFD, BFD, CFD, size);
    gpu_boundaries3<<<bpg9,tpb>>>(CID, DID, DFD, BFD, CFD, EID);

    cudaMemcpy(AIH, AID, SIZEINT(size),  cudaMemcpyDeviceToHost);
    cudaMemcpy(BIH, BID, SIZEINT(size),  cudaMemcpyDeviceToHost);
    cudaMemcpy(AFH, AFD, SIZEFLT(9*M*N), cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, SIZEFLT(M*N),   cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, SIZEFLT(M*N),   cudaMemcpyDeviceToHost);

    cudaMemcpy(DFH, DFD, SIZEFLT(9*M*N), cudaMemcpyDeviceToHost);

    int sumAIb = sumHostInt(AIH, size);
    int sumBIb = sumHostInt(BIH, size);
    FLOAT_TYPE sumAFb = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFb = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCFb = sumHostFlt(CFH, M*N);

    FLOAT_TYPE sumDFb = sumHostFlt(DFH, 9*M*N);

    CuAssertIntEquals(tc, sumAIa, sumAIb);
    CuAssertIntEquals(tc, sumBIa, sumBIb);
    CuAssertDblEquals(tc, sumBFa, sumBFb, 0.000001);
    CuAssertDblEquals(tc, sumCFa, sumCFb, 0.000001);

    CuAssertDblEquals_Msg(tc, "f", sumAFb, sumDFb, 1e-12);
    //CuAssertIntEquals_Msg(tc, "f-diff", 0, compareArraysFlt(AFH, DFH, 9*M*N));

    EFD = createGpuArrayFlt(9*M*N);
    FFD = createGpuArrayFlt(9*M*N);
    //CuAssertIntEquals_Msg(tc, "f", 0, compareArraysFlt(AFH, DFH, 9*M*N));
    FLOAT_TYPE val = computeResidual(AFD, DFD, EFD, FFD, M, N);
    CuAssertDblEquals_Msg(tc, "f-l2", 0, val, 1e-6);
}

///compare outlet conditions with Zhu-He profile @param tc test case
void testCompareOutletProfile1(CuTest *tc)
{
    printBanner("Test compare outlet (zhu he)");
    runTestCompareOutlet(tc, OUTLET);
}

///compare outlet conditions with 2nd order open boundary @param tc test case
void testCompareOutletProfile2(CuTest *tc)
{
    printBanner("Test compare outlet (second)");
    runTestCompareOutlet(tc, OUTLET_SECOND);
}

///compare outlet conditions with 1st order open boundary @param tc test case
void testCompareOutletProfile3(CuTest *tc)
{
    printBanner("Test compare outlet (first)");
    runTestCompareOutlet(tc, OUTLET_FIRST);
}

///Clean up after test case
void cleanupTestCompareOutlet()
{
    cudaFree(AID); cudaFree(BID); cudaFree(CID); cudaFree(DID); cudaFree(EID); cudaFree(AFD); cudaFree(BFD); cudaFree(CFD); cudaFree(DFD); cudaFree(EFD); cudaFree(FFD);
    free(AIH); free(BIH); free(CIH); free(DIH); free(EIH); free(AFH); free(BFH); free(CFH); free(DFH);
}

/**
 * @brief Test to compare results from #gpuBcInlet and #gpu_boundaries1
 *
 * @param tc test case
 * @test
 *  - Prepare boundary conditions for the lid driven cavity
 *  - Run the function
 *  - Check the sum of all the arrays between the two algorithms
 */
void testCompareInlet(CuTest *tc)
{
    printBanner("Test compare inlet");
    M = 64;
    N = 64;
    int size = getLidBcMaskSize(M,N);
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    dim3 tpb(THREADS);
    dim3 bpgB((int)(size/THREADS)+1);
    dim3 bpg9((int)(9*M*N/THREADS)+1);

    FLOAT_TYPE af = 0.002;
    FLOAT_TYPE bf = 0.3;

    AIH = createHostArrayInt(size, ARRAY_FILL, 0); //bcidx
    BIH = createHostArrayInt(size, ARRAY_FILL, 0); //bcmask

    createLidBcIdx(AIH, M, N);
    createLidBcMask(BIH, M, N);

    CIH = createHostArrayInt(M*N, ARRAY_FILL, 0); //boundary
    DIH = createHostArrayInt(9*M*N, ARRAY_FILL, 0); //bcid
    EIH = createHostArrayInt(M*N, ARRAY_FILL, 0); //corner
    FIH = createHostArrayInt(M*N, ARRAY_FILL, 0); //fluid

    createLidBcBoundary(CIH, M, N);
    createLidBcBcId(DIH, M, N);
    createLidBcCorner(EIH, M, N);
    createLidBcFluid(FIH, M, N);

    AID = createGpuArrayInt(size, ARRAY_COPY, 0, AIH);
    BID = createGpuArrayInt(size, ARRAY_COPY, 0, BIH);
    CID = createGpuArrayInt(M*N, ARRAY_COPY, 0, CIH);
    DID = createGpuArrayInt(9*M*N, ARRAY_COPY, 0, DIH);
    EID = createGpuArrayInt(M*N, ARRAY_COPY, 0, EIH);
    FID = createGpuArrayInt(M*N, ARRAY_COPY, 0, FIH);

    AFH = createHostArrayFlt(9*M*N, ARRAY_RAND, af); //f1
    BFH = createHostArrayFlt(M*N, ARRAY_FILL, bf); //u0
    CFH = createHostArrayFlt(M*N, ARRAY_FILL, 0.); //v0
    DFH = createHostArrayFlt(9*M*N, ARRAY_COPY, 0, AFH); //f2

    AFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, BFH);
    CFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, CFH);
    DFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, DFH);

    int sumAIa = sumHostInt(AIH, size);
    int sumBIa = sumHostInt(BIH, size);
    int sumCIa = sumHostInt(CIH, M*N);
    int sumDIa = sumHostInt(DIH, 9*M*N);
    int sumEIa = sumHostInt(EIH, M*N);
    int sumFIa = sumHostInt(FIH, M*N);
    FLOAT_TYPE sumAFa = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFa = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCFa = sumHostFlt(CFH, M*N);
    FLOAT_TYPE sumDFa = sumHostFlt(DFH, 9*M*N);

    gpuBcInlet<<<bpgB,tpb>>>(AID, BID, AFD, BFD, CFD, size);
    gpu_boundaries1<<<bpg9,tpb>>>(FID, CID, DID, DFD, BFD, CFD, EID);

    cudaMemcpy(AIH, AID, SIZEINT(size),  cudaMemcpyDeviceToHost);
    cudaMemcpy(BIH, BID, SIZEINT(size),  cudaMemcpyDeviceToHost);
    cudaMemcpy(CIH, CID, SIZEINT(M*N),  cudaMemcpyDeviceToHost);
    cudaMemcpy(DIH, DID, SIZEINT(9*M*N),  cudaMemcpyDeviceToHost);
    cudaMemcpy(EIH, EID, SIZEINT(M*N),  cudaMemcpyDeviceToHost);
    cudaMemcpy(FIH, FID, SIZEINT(M*N),  cudaMemcpyDeviceToHost);
    cudaMemcpy(AFH, AFD, SIZEFLT(9*M*N), cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, SIZEFLT(M*N),   cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, SIZEFLT(M*N),   cudaMemcpyDeviceToHost);
    cudaMemcpy(DFH, DFD, SIZEFLT(9*M*N),   cudaMemcpyDeviceToHost);

    int sumAIb = sumHostInt(AIH, size);
    int sumBIb = sumHostInt(BIH, size);
    int sumCIb = sumHostInt(CIH, M*N);
    int sumDIb = sumHostInt(DIH, 9*M*N);
    int sumEIb = sumHostInt(EIH, M*N);
    int sumFIb = sumHostInt(FIH, M*N);
    FLOAT_TYPE sumAFb = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFb = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCFb = sumHostFlt(CFH, M*N);
    FLOAT_TYPE sumDFb = sumHostFlt(DFH, 9*M*N);

    CuAssertIntEquals(tc, sumAIa, sumAIb);
    CuAssertIntEquals(tc, sumBIa, sumBIb);
    CuAssertIntEquals(tc, sumCIa, sumCIb);
    CuAssertIntEquals(tc, sumDIa, sumDIb);
    CuAssertIntEquals(tc, sumEIa, sumEIb);
    CuAssertIntEquals(tc, sumFIa, sumFIb);
    // CuAssertDblEquals(tc, sumAFa, sumAFb, 0.000001);
    CuAssertDblEquals(tc, sumBFa, sumBFb, 0.000001);
    CuAssertDblEquals(tc, sumCFa, sumCFb, 0.000001);
    // CuAssertDblEquals(tc, sumDFa, sumDFb, 0.000001);
    CuAssertDblEquals(tc, sumAFb, sumDFb, 1e-12);
}

///Clean up after test case
void cleanupTestCompareInlet()
{
    cudaFree(AID); cudaFree(BID); cudaFree(CID); cudaFree(DID); cudaFree(EID); cudaFree(FID);
    cudaFree(AFD); cudaFree(BFD); cudaFree(CFD); cudaFree(DFD);
    free(AIH); free(BIH); free(CIH); free(DIH); free(EIH); free(FIH); free(AFH); free(BFH); free(CFH); free(DFH);
}

CuSuite* gpuBoundariesGetSuite()
{
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TCLN(suite, testGpuBoundaries1,        cleanupTestGpuBoundaries);
    SUITE_ADD_TCLN(suite, testGpuBoundaries2_wcb,    cleanupTestGpuBoundaries2);
    SUITE_ADD_TCLN(suite, testGpuBoundaries2_wocb,   cleanupTestGpuBoundaries2);
    SUITE_ADD_TCLN(suite, testGpuBoundaries3,        cleanupTestGpuBoundaries);
    SUITE_ADD_TCLN(suite, testGpuBcInlet,            cleanupTestGpuBcInlet);
    SUITE_ADD_TCLN(suite, testGpuBcWall_wcb,         cleanupTestGpuBcWall);
    SUITE_ADD_TCLN(suite, testGpuBcWall_wocb,        cleanupTestGpuBcWall);
    SUITE_ADD_TCLN(suite, testGpuBcOutlet,           cleanupTestGpuBcOutlet);
    SUITE_ADD_TCLN(suite, testCompareOutletProfile1, cleanupTestCompareOutlet);
    SUITE_ADD_TCLN(suite, testCompareOutletProfile2, cleanupTestCompareOutlet);
    SUITE_ADD_TCLN(suite, testCompareOutletProfile3, cleanupTestCompareOutlet);
    SUITE_ADD_TEST(suite, testCompareWall);
    SUITE_ADD_TCLN(suite, testCompareInlet,          cleanupTestCompareInlet);
    return suite;
}