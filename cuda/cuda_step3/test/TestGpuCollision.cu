/**
 * Unittests for the boundary conditions
 * @file TestGpuCollision.cu
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
#include "CellFunctions.h"

/**
 * @brief Test to compare results from #gpuCollMrt and #gpu_mrt1 #gpu_mrt2
 *
 * @param tc test case
 * @test
 *  - Prepare boundary conditions for the lid driven cavity
 *  - Run the function
 *  - Check the sum of all the arrays between the two algorithms
 */
void testCompareMrt(CuTest *tc)
{
    M = 1281;
    N = 1282;
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    dim3 tpb(THREADS);
    dim3 bpg1((int)(M*N/THREADS)+1);
    dim3 bpg9((int)(9*M*N/THREADS)+1);

    FLOAT_TYPE omega = 0.02;
    FLOAT_TYPE *velMomMap  = createHostArrayFlt(81);
    FLOAT_TYPE *momCollMtx = createHostArrayFlt(81);
    MRTInitializer(velMomMap, momCollMtx, omega);

    cudaMemcpyToSymbol(velMomMap_d, velMomMap, 81*sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(momCollMtx_d, momCollMtx, 81*sizeof(FLOAT_TYPE));

    printBanner("Test compare mrt models");
    AIH = createHostArrayInt(M*N); //fluid
    AFH = createHostArrayFlt(M*N, ARRAY_RAND, 0.1); //rho
    BFH = createHostArrayFlt(M*N, ARRAY_RAND, 0.05); //u
    CFH = createHostArrayFlt(M*N, ARRAY_ZERO); //v
    DFH = createHostArrayFlt(9*M*N, ARRAY_FILL, 0.003); //f
    EFH = createHostArrayFlt(9*M*N, ARRAY_ZERO); //fColl
    FFH = createHostArrayFlt(9*M*N, ARRAY_ZERO); //mEq
    GFH = createHostArrayFlt(9*M*N, ARRAY_ZERO); //m
    HFH = createHostArrayFlt(9*M*N, ARRAY_ZERO); //coll

    createLidBcFluid(AIH, M, N);

    AID = createGpuArrayInt(M*N, ARRAY_COPY, 0, AIH);
    AFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, BFH);
    CFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, CFH);
    DFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, DFH);
    EFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, EFH);
    FFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, FFH);
    GFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, GFH);
    HFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, HFH);

    FLOAT_TYPE *EFH2 = createHostArrayFlt(9*M*N, ARRAY_COPY, 0, EFH);
    FLOAT_TYPE *EFD2 = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, EFH);

    int sumAIa = sumHostInt(AIH, M*N);
    FLOAT_TYPE sumAFa = sumHostFlt(AFH, M*N);
    FLOAT_TYPE sumBFa = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCFa = sumHostFlt(CFH, M*N);
    FLOAT_TYPE sumDFa = sumHostFlt(DFH, 9*M*N);
    FLOAT_TYPE sumEFa = sumHostFlt(EFH, 9*M*N);
    FLOAT_TYPE sumEF2a = sumHostFlt(EFH2, 9*M*N);
    FLOAT_TYPE sumFFa = sumHostFlt(FFH, 9*M*N);
    FLOAT_TYPE sumGFa = sumHostFlt(GFH, 9*M*N);
    FLOAT_TYPE sumHFa = sumHostFlt(HFH, 9*M*N);

    gpuCollMrt<<<bpg1,tpb>>>(AID, AFD, BFD, CFD, DFD, EFD);
    gpu_mrt1<<<bpg9,tpb>>>(AID, AFD, BFD, CFD, DFD, FFD, GFD);
    gpu_mrt2<<<bpg9,tpb>>>(AID, HFD, GFD, FFD, EFD2, DFD);

    cudaMemcpy(AIH, AID, SIZEINT(N*M),   cudaMemcpyDeviceToHost);
    cudaMemcpy(AFH, AFD, SIZEFLT(N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, SIZEFLT(N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, SIZEFLT(N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(DFH, DFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(EFH, EFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(EFH2, EFD2, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(FFH, FFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(GFH, GFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(HFH, HFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);

    int sumAIb = sumHostInt(AIH, M*N);
    FLOAT_TYPE sumAFb = sumHostFlt(AFH, M*N);
    FLOAT_TYPE sumBFb = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCFb = sumHostFlt(CFH, M*N);
    FLOAT_TYPE sumDFb = sumHostFlt(DFH, 9*M*N);
    FLOAT_TYPE sumEFb = sumHostFlt(EFH, 9*M*N);
    FLOAT_TYPE sumEF2b = sumHostFlt(EFH2, 9*M*N);
    FLOAT_TYPE sumFFb = sumHostFlt(FFH, 9*M*N);
    FLOAT_TYPE sumGFb = sumHostFlt(GFH, 9*M*N);
    FLOAT_TYPE sumHFb = sumHostFlt(HFH, 9*M*N);

    CuAssertIntEquals_Msg(tc, "fluid", sumAIa, sumAIb);
    CuAssertDblEquals_Msg(tc, "rho", sumAFa, sumAFb, 0.00001);
    CuAssertDblEquals_Msg(tc, "u", sumBFa, sumBFb, 0.00001);
    CuAssertDblEquals_Msg(tc, "v", sumCFa, sumCFb, 0.00001);
    CuAssertDblEquals_Msg(tc, "f", sumDFa, sumDFb, 0.00001);
    CuAssertDblEquals_Msg(tc, "fColl", sumEFb, sumEF2b, 0.00001);
    FLOAT_TYPE val = computeResidual(EFD, EFD2, FFD, GFD, M, N);
    CuAssertDblEquals_Msg(tc, "fCollNorm", val, 0, 0.00001);

    free(EFH2); cudaFree(EFD2);
}

///Clean up after test case
void cleanupTestCompareMrt()
{
    free(AIH); free(AFH);free(BFH);free(CFH);free(DFH);free(EFH);free(FFH);free(GFH);free(HFH);
    cudaFree(AID); cudaFree(AFD);cudaFree(BFD);cudaFree(CFD);cudaFree(DFD);cudaFree(EFD);cudaFree(FFD);cudaFree(GFD);cudaFree(HFD);
}

/**
 * @brief Test to compare results from #gpuCollTrt and #gpu_trt1 #gpu_trt2
 *
 * @param tc test case
 * @test
 *  - Prepare boundary conditions for the lid driven cavity
 *  - Run the function
 *  - Check the sum of all the arrays between the two algorithms
 */
void testCompareTrt(CuTest *tc)
{
    printBanner("Test compare trt models");
    M = 1281;
    N = 1282;
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    dim3 tpb(THREADS);
    dim3 bpg1((int)(M*N/THREADS)+1);
    dim3 bpg9((int)(9*M*N/THREADS)+1);

    FLOAT_TYPE omega = 0.02;
    FLOAT_TYPE *velMomMap  = createHostArrayFlt(81);
    FLOAT_TYPE *momCollMtx = createHostArrayFlt(81);
    MRTInitializer(velMomMap, momCollMtx, omega);

    cudaMemcpyToSymbol(velMomMap_d, velMomMap, 81*sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(momCollMtx_d, momCollMtx, 81*sizeof(FLOAT_TYPE));

    AIH = createHostArrayInt(M*N); //fluid
    AFH = createHostArrayFlt(M*N, ARRAY_RAND, 0.1); //rho
    BFH = createHostArrayFlt(M*N, ARRAY_RAND, 0.05); //u
    CFH = createHostArrayFlt(M*N, ARRAY_ZERO); //v
    DFH = createHostArrayFlt(9*M*N, ARRAY_FILL, 0.000003); //f
    EFH = createHostArrayFlt(9*M*N, ARRAY_ZERO); //fColl
    FFH = createHostArrayFlt(9*M*N, ARRAY_ZERO); //mEq
    GFH = createHostArrayFlt(9*M*N, ARRAY_ZERO); //m
    HFH = createHostArrayFlt(9*M*N, ARRAY_ZERO); //coll

    createLidBcFluid(AIH, M, N);

    AID = createGpuArrayInt(M*N, ARRAY_COPY, 0, AIH);
    AFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, BFH);
    CFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, CFH);
    DFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, DFH);
    EFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, EFH);
    FFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, FFH);
    GFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, GFH);
    HFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, HFH);

    FLOAT_TYPE *EFH2 = createHostArrayFlt(9*M*N, ARRAY_COPY, 0, EFH);
    FLOAT_TYPE *EFD2 = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, EFH);

    int sumAIa = sumHostInt(AIH, M*N);
    FLOAT_TYPE sumAFa = sumHostFlt(AFH, M*N);
    FLOAT_TYPE sumBFa = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCFa = sumHostFlt(CFH, M*N);
    FLOAT_TYPE sumDFa = sumHostFlt(DFH, 9*M*N);
    FLOAT_TYPE sumEFa = sumHostFlt(EFH, 9*M*N);
    FLOAT_TYPE sumEF2a = sumHostFlt(EFH2, 9*M*N);
    FLOAT_TYPE sumFFa = sumHostFlt(FFH, 9*M*N);
    FLOAT_TYPE sumGFa = sumHostFlt(GFH, 9*M*N);
    FLOAT_TYPE sumHFa = sumHostFlt(HFH, 9*M*N);

    gpuCollTrt<<<bpg1,tpb>>>(AID, AFD, BFD, CFD, DFD, EFD);
    gpu_trt1<<<bpg9,tpb>>>(AID, FFD, BFD, CFD, DFD);
    gpu_trt2<<<bpg9,tpb>>>(AID, FFD, DFD, EFD2);

    cudaMemcpy(AIH, AID, SIZEINT(N*M),   cudaMemcpyDeviceToHost);
    cudaMemcpy(AFH, AFD, SIZEFLT(N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, SIZEFLT(N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, SIZEFLT(N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(DFH, DFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(EFH, EFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(EFH2, EFD2, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(FFH, FFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(GFH, GFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(HFH, HFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);

    int sumAIb = sumHostInt(AIH, M*N);
    FLOAT_TYPE sumAFb = sumHostFlt(AFH, M*N);
    FLOAT_TYPE sumBFb = sumHostFlt(BFH, M*N);
    FLOAT_TYPE sumCFb = sumHostFlt(CFH, M*N);
    FLOAT_TYPE sumDFb = sumHostFlt(DFH, 9*M*N);
    FLOAT_TYPE sumEFb = sumHostFlt(EFH, 9*M*N);
    FLOAT_TYPE sumEF2b = sumHostFlt(EFH2, 9*M*N);
    FLOAT_TYPE sumFFb = sumHostFlt(FFH, 9*M*N);
    FLOAT_TYPE sumGFb = sumHostFlt(GFH, 9*M*N);
    FLOAT_TYPE sumHFb = sumHostFlt(HFH, 9*M*N);

    CuAssertIntEquals_Msg(tc, "fluid", sumAIa, sumAIb);
    CuAssertDblEquals_Msg(tc, "rho", sumAFa, sumAFb, 0.00001);
    CuAssertDblEquals_Msg(tc, "u", sumBFa, sumBFb, 0.00001);
    CuAssertDblEquals_Msg(tc, "v", sumCFa, sumCFb, 0.00001);
    CuAssertDblEquals_Msg(tc, "f", sumDFa, sumDFb, 0.00001);
    CuAssertDblEquals_Msg(tc, "fColl", sumEFb, sumEF2b, 0.00001);
    FLOAT_TYPE val = computeResidual(EFD, EFD2, FFD, GFD, M, N);
    CuAssertDblEquals_Msg(tc, "fCollNorm", val, 0, 0.00001);

    free(EFH2); cudaFree(EFD2);
}

///Clean up after test case
void cleanupTestCompareTrt()
{
    free(AIH); free(AFH);free(BFH);free(CFH);free(DFH);free(EFH);free(FFH);free(GFH);free(HFH);
    cudaFree(AID); cudaFree(AFD);cudaFree(BFD);cudaFree(CFD);cudaFree(DFD);cudaFree(EFD);cudaFree(FFD);cudaFree(GFD);cudaFree(HFD);
}

CuSuite* gpuCollisionGetSuite()
{
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TCLN(suite, testCompareMrt, cleanupTestCompareMrt);
    SUITE_ADD_TCLN(suite, testCompareTrt, cleanupTestCompareTrt);
    return suite;
}