/**
 * Unittests for macroscopic value update
 * @file TestGpuUpdateMacro.cu
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#include "CuTest.h"
#include "GpuFunctions.h"
#include "ShellFunctions.h"
#include "TestUtils.h"
#include "ArrayUtils.h"
#include "GpuConstants.h"

static const int cx[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
static const int cy[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
static const FLOAT_TYPE hf1 = 0.3;
static const FLOAT_TYPE hf2 = 0.4;
static const int boundaryId = 2;
static FLOAT_TYPE rhoSum[3];
static FLOAT_TYPE uSum[3];

/**
 * @brief Unittest for #gpu_update_macro
 *
 * @param tc test case
 * @test
 *   - Prepare arrays for lid driven cavity
 *   - call function
 *   - check the sum of the arrays against predefined values
 */
void testGpuUpdateMacro_(CuTest *tc)
{
    printBanner("Test gpu_update_macro");
    M = 256;
    N = 256;
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    cudaMemcpyToSymbol(cx_d,  cx,  9*sizeof(int));
    cudaMemcpyToSymbol(cy_d,  cy,  9*sizeof(int));
    dim3 tpb(THREADS);
    dim3 bpg1((int)((M*N-1)/THREADS)+1);

    AIH = createHostArrayInt(  M*N, ARRAY_FILL, 0);
    BIH = createHostArrayInt(9*M*N, ARRAY_FILL, 0);
    CIH = createHostArrayInt(  M*N, ARRAY_FILL, 0);

    createLidBcFluid(AIH, M, N);
    createLidBcBcId(BIH, M, N);
    createLidBcBoundary(CIH, M, N);

    AID = createGpuArrayInt(  M*N, ARRAY_COPY, 0, AIH);  //fluid
    BID = createGpuArrayInt(9*M*N, ARRAY_COPY, 0, BIH);  //bcid
    CID = createGpuArrayInt(  M*N, ARRAY_COPY, 0, CIH);  //boun

    AFH = createHostArrayFlt(  M*N, ARRAY_FILL, 0);   //rho
    BFH = createHostArrayFlt(  M*N, ARRAY_FILL, 0);   //u
    CFH = createHostArrayFlt(  M*N, ARRAY_FILL, 0);   //v
    DFH = createHostArrayFlt(  M*N, ARRAY_FILL, 0);   //drag
    EFH = createHostArrayFlt(  M*N, ARRAY_FILL, 0);   //lift
    FFH = createHostArrayFlt(  M*N, ARRAY_FILL, 0);   //x
    GFH = createHostArrayFlt(  M*N, ARRAY_FILL, 0);   //y
    HFH = createHostArrayFlt(9*M*N, ARRAY_FILL, hf1); //f

    fillLidCoordinates(FFH, GFH, M, N);
    fillFDir(HFH, hf2, 1, M*N);

    AFD = createGpuArrayFlt(  M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(  M*N, ARRAY_COPY, 0, BFH);
    CFD = createGpuArrayFlt(  M*N, ARRAY_COPY, 0, CFH);
    DFD = createGpuArrayFlt(  M*N, ARRAY_COPY, 0, DFH);
    EFD = createGpuArrayFlt(  M*N, ARRAY_COPY, 0, EFH);
    FFD = createGpuArrayFlt(  M*N, ARRAY_COPY, 0, FFH);
    GFD = createGpuArrayFlt(  M*N, ARRAY_COPY, 0, GFH);
    HFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, HFH);

    int sumAIa = sumHostInt(AIH,   M*N);
    int sumBIa = sumHostInt(BIH, 9*M*N);
    int sumCIa = sumHostInt(CIH,   M*N);
    FLOAT_TYPE sumAFa = sumHostFlt(AFH,   M*N);
    FLOAT_TYPE sumBFa = sumHostFlt(BFH,   M*N);
    FLOAT_TYPE sumCFa = sumHostFlt(CFH,   M*N);
    FLOAT_TYPE sumDFa = sumHostFlt(DFH,   M*N);
    FLOAT_TYPE sumEFa = sumHostFlt(EFH,   M*N);
    FLOAT_TYPE sumFFa = sumHostFlt(FFH,   M*N);
    FLOAT_TYPE sumGFa = sumHostFlt(GFH,   M*N);
    FLOAT_TYPE sumHFa = sumHostFlt(HFH, 9*M*N);

    gpu_update_macro<<<bpg1,tpb>>>(AID, AFD, BFD, CFD, BID, CID, DFD, EFD, FFD, GFD, HFD);

    cudaMemcpy(AIH, AID, SIZEINT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(BIH, BID, SIZEINT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(CIH, CID, SIZEINT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(AFH, AFD, SIZEFLT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, SIZEFLT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, SIZEFLT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(DFH, DFD, SIZEFLT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(EFH, EFD, SIZEFLT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(FFH, FFD, SIZEFLT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(GFH, GFD, SIZEFLT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(HFH, HFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);

    int sumAIb = sumHostInt(AIH,   M*N);
    int sumBIb = sumHostInt(BIH, 9*M*N);
    int sumCIb = sumHostInt(CIH,   M*N);
    FLOAT_TYPE sumAFb = sumHostFlt(AFH,   M*N);
    FLOAT_TYPE sumBFb = sumHostFlt(BFH,   M*N);
    FLOAT_TYPE sumCFb = sumHostFlt(CFH,   M*N);
    FLOAT_TYPE sumDFb = sumHostFlt(DFH,   M*N);
    FLOAT_TYPE sumEFb = sumHostFlt(EFH,   M*N);
    FLOAT_TYPE sumFFb = sumHostFlt(FFH,   M*N);
    FLOAT_TYPE sumGFb = sumHostFlt(GFH,   M*N);
    FLOAT_TYPE sumHFb = sumHostFlt(HFH, 9*M*N);

    FLOAT_TYPE rhoB = 8*M*N*hf1+M*N*hf2;
    FLOAT_TYPE uB   = (hf2/(8*hf1+hf2))*M*N-(hf1/(8*hf1+hf2))*M*N;

    CuAssertIntEquals_Msg(tc, "fluid", sumAIa, sumAIb);
    CuAssertIntEquals_Msg(tc, "bcid",  sumBIa, sumBIb);
    CuAssertIntEquals_Msg(tc, "bndid", sumCIa, sumCIb);
    CuAssertDblEquals_Msg(tc, "rho",   rhoB,   sumAFb, 0.002*9*M*N); //TODO that's a huge numerical error
    CuAssertDblEquals_Msg(tc, "u",     uB,     sumBFb, 2);
    CuAssertDblEquals_Msg(tc, "v",     sumCFa, sumCFb, 0.000001);
    CuAssertDblEquals_Msg(tc, "drag",  sumDFa, sumDFb, 0.000001);
    CuAssertDblEquals_Msg(tc, "lift",  sumEFa, sumEFb, 0.000001);
    CuAssertDblEquals_Msg(tc, "x",     sumFFa, sumFFb, 0.000001);
    CuAssertDblEquals_Msg(tc, "y",     sumGFa, sumGFb, 0.000001);
    CuAssertDblEquals_Msg(tc, "f",     sumHFa, sumHFb, 0.000001);
    //printf(" rho_sum: %f\tf_sum: %f\n", sumAFb, sumHFb);
    rhoSum[0] = sumAFb;
    uSum[0]   = sumBFb;
}

///Clean up after test case
void cleanupTestGpuUpdateMacro()
{
    cudaFree(AID); cudaFree(BID); cudaFree(CID);
    cudaFree(AFD); cudaFree(BFD); cudaFree(CFD); cudaFree(DFD); cudaFree(EFD); cudaFree(FFD); cudaFree(GFD); cudaFree(HFD);
    free(AIH); free(BIH); free(CIH);
    free(AFH); free(BFH); free(CFH); free(DFH); free(EFH); free(FFH); free(GFH); free(HFH);
}

/**
 * @brief Prepare host array for lid driven cavity
 *
 * @param hf value to fill the distribution function
 */
void runTestGpuUpdatePart1(FLOAT_TYPE hf)
{
    AIH = createHostArrayInt(M*N, ARRAY_FILL, 0);     //fluid
    BIH = createHostArrayInt(M*N, ARRAY_FILL, 0);     //bcmask

    createLidBcFluid(AIH, M, N);
    createLidBcMaskFull(BIH, M, N);

    AFH = createHostArrayFlt(  M*N, ARRAY_FILL, 0);   //rho
    BFH = createHostArrayFlt(  M*N, ARRAY_FILL, 0);   //u
    CFH = createHostArrayFlt(  M*N, ARRAY_FILL, 0);   //v
    DFH = createHostArrayFlt(  M*N, ARRAY_FILL, 0);   //drag
    EFH = createHostArrayFlt(  M*N, ARRAY_FILL, 0);   //lift
    FFH = createHostArrayFlt(  M*N, ARRAY_FILL, 0);   //x
    GFH = createHostArrayFlt(  M*N, ARRAY_FILL, 0);   //y
    HFH = createHostArrayFlt(9*M*N, ARRAY_FILL, hf);  //f

    fillLidCoordinates(FFH, GFH, M, N);
}

/**
 * @brief Compute sums and copy host arrays to GPU arrays
 *
 * @param IntSums integer arrays sums
 * @param FltSums float arrays sums
 */
void runTestGpuUpdatePart2(int *IntSums, FLOAT_TYPE *FltSums)
{

    IntSums[0] = sumHostInt(AIH, M*N);
    IntSums[1] = sumHostInt(BIH, M*N);

    FltSums[0] = sumHostFlt(AFH, M*N);
    FltSums[1] = sumHostFlt(BFH, M*N);
    FltSums[2] = sumHostFlt(CFH, M*N);
    FltSums[3] = sumHostFlt(DFH, M*N);
    FltSums[4] = sumHostFlt(EFH, M*N);
    FltSums[5] = sumHostFlt(FFH, M*N);
    FltSums[6] = sumHostFlt(GFH, M*N);
    FltSums[7] = sumHostFlt(HFH, 9*M*N);

    AID = createGpuArrayInt(  M*N, ARRAY_COPY, 0, AIH);  //fluid
    BID = createGpuArrayInt(  M*N, ARRAY_COPY, 0, BIH);  //bcmask

    AFD = createGpuArrayFlt(  M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(  M*N, ARRAY_COPY, 0, BFH);
    CFD = createGpuArrayFlt(  M*N, ARRAY_COPY, 0, CFH);
    DFD = createGpuArrayFlt(  M*N, ARRAY_COPY, 0, DFH);
    EFD = createGpuArrayFlt(  M*N, ARRAY_COPY, 0, EFH);
    FFD = createGpuArrayFlt(  M*N, ARRAY_COPY, 0, FFH);
    GFD = createGpuArrayFlt(  M*N, ARRAY_COPY, 0, GFH);
    HFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, HFH);
}

/**
 * @brief Copy GPU array to host arrays and compute their sums
 *
 * @param IntSums integer arrays sums
 * @param FltSums float arrays sums
 */
void runTestGpuUpdatePart3(int *IntSums, FLOAT_TYPE *FltSums)
{
    cudaMemcpy(AIH, AID, SIZEINT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(BIH, BID, SIZEINT(  N*M), cudaMemcpyDeviceToHost);

    cudaMemcpy(AFH, AFD, SIZEFLT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, SIZEFLT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, SIZEFLT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(DFH, DFD, SIZEFLT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(EFH, EFD, SIZEFLT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(FFH, FFD, SIZEFLT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(GFH, GFD, SIZEFLT(  N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(HFH, HFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);

    IntSums[0] = sumHostInt(AIH,   M*N);
    IntSums[1] = sumHostInt(BIH,   M*N);

    FltSums[0] = sumHostFlt(AFH,   M*N);
    FltSums[1] = sumHostFlt(BFH,   M*N);
    FltSums[2] = sumHostFlt(CFH,   M*N);
    FltSums[3] = sumHostFlt(DFH,   M*N);
    FltSums[4] = sumHostFlt(EFH,   M*N);
    FltSums[5] = sumHostFlt(FFH,   M*N);
    FltSums[6] = sumHostFlt(GFH,   M*N);
    FltSums[7] = sumHostFlt(HFH, 9*M*N);
}

/**
 * @brief Unittest for #gpu_update_new
 *
 * @param tc test case
 * @test
 *   - Prepare arrays for lid driven cavity
 *   - call function
 *   - check the sum of the arrays against predefined values
 */
void testGpuUpdateNew(CuTest *tc)
{
    printBanner("Test gpu_update_new");
    M = 256;
    N = 256;
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    cudaMemcpyToSymbol(dlBoundaryId_d, &boundaryId, sizeof(int));
    cudaMemcpyToSymbol(cx_d,  cx,  9*sizeof(int));
    cudaMemcpyToSymbol(cy_d,  cy,  9*sizeof(int));
    dim3 tpb(THREADS);
    dim3 bpg1((int)((M*N-1)/THREADS)+1);

    int IntSumsA[2], IntSumsB[2];
    FLOAT_TYPE FltSumsA[8], FltSumsB[8];

    runTestGpuUpdatePart1(hf1);
    fillFDir(HFH, hf2, 1, M*N);
    runTestGpuUpdatePart2(IntSumsA, FltSumsA);
    gpu_update_new<<<bpg1,tpb>>>(AID, AFD, BFD, CFD, BID, DFD, EFD, FFD, GFD, HFD);
    runTestGpuUpdatePart3(IntSumsB, FltSumsB);

    FLOAT_TYPE rhoB = 8*M*N*hf1+M*N*hf2;
    FLOAT_TYPE uB   = (hf2/(8*hf1+hf2))*M*N-(hf1/(8*hf1+hf2))*M*N;

    CuAssertIntEquals_Msg(tc, "fluid", IntSumsA[0], IntSumsB[0]);
    CuAssertIntEquals_Msg(tc, "mask ", IntSumsA[1], IntSumsB[1]);

    CuAssertDblEquals_Msg(tc, "rho  ", rhoB,        FltSumsB[0], 0.002*9*M*N); //TODO that's a huge numerical error +9*M*N*hf
    CuAssertDblEquals_Msg(tc, "u    ", uB,          FltSumsB[1], 2);
    CuAssertDblEquals_Msg(tc, "v    ", FltSumsA[2], FltSumsB[2], 0.000001);
    CuAssert             (tc, "drag ", FltSumsA[3] < FltSumsB[3]);
    CuAssert             (tc, "lift ", FltSumsA[4] < FltSumsB[4]);
    CuAssertDblEquals_Msg(tc, "coorx", FltSumsA[5], FltSumsB[5], 0.000001);
    CuAssertDblEquals_Msg(tc, "coory", FltSumsA[6], FltSumsB[6], 0.000001);
    CuAssertDblEquals_Msg(tc, "f    ", FltSumsA[7], FltSumsB[7], 0.000001);
    // printf(" rho_sum: %f\tf_sum: %f\n", FltSumsB[0], FltSumsB[7]);
    rhoSum[1] = FltSumsB[0];
    uSum[1]   = FltSumsB[1];
}

/**
 * @brief Unittest for #gpuUpdateMacro
 *
 * @param tc test case
 * @test
 *   - Prepare arrays for lid driven cavity
 *   - call function
 *   - check the sum of the arrays against predefined values
 */
void testGpuUpdateMacro(CuTest *tc)
{
    printBanner("Test gpuUpdateMacro");
    M = 256;
    N = 256;
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    cudaMemcpyToSymbol(dlBoundaryId_d, &boundaryId, sizeof(int));
    dim3 tpb(THREADS);
    dim3 bpg1((int)((M*N-1)/THREADS)+1);

    int IntSumsA[2], IntSumsB[2];
    FLOAT_TYPE FltSumsA[8], FltSumsB[8];

    runTestGpuUpdatePart1(hf1);
    fillFDir(HFH, hf2, 1, M*N);
    runTestGpuUpdatePart2(IntSumsA, FltSumsA);
    gpuUpdateMacro<<<bpg1,tpb>>>(AID, AFD, BFD, CFD, BID, DFD, EFD, FFD, GFD, HFD);
    runTestGpuUpdatePart3(IntSumsB, FltSumsB);

    FLOAT_TYPE rhoB = 8*M*N*hf1+M*N*hf2;
    FLOAT_TYPE uB   = (hf2/(8*hf1+hf2))*M*N-(hf1/(8*hf1+hf2))*M*N;

    CuAssertIntEquals_Msg(tc, "fluid", IntSumsA[0], IntSumsB[0]);
    CuAssertIntEquals_Msg(tc, "mask ", IntSumsA[1], IntSumsB[1]);

    CuAssertDblEquals_Msg(tc, "rho  ", rhoB,        FltSumsB[0], 0.002*9*M*N); //TODO that's a huge numerical error +9*M*N*hf
    CuAssertDblEquals_Msg(tc, "u    ", uB,          FltSumsB[1], 2);
    CuAssertDblEquals_Msg(tc, "v    ", FltSumsA[2], FltSumsB[2], 0.000001);
    CuAssert             (tc, "drag ", FltSumsA[3] < FltSumsB[3]);
    CuAssert             (tc, "lift ", FltSumsA[4] < FltSumsB[4]);
    CuAssertDblEquals_Msg(tc, "coorx", FltSumsA[5], FltSumsB[5], 0.000001);
    CuAssertDblEquals_Msg(tc, "coory", FltSumsA[6], FltSumsB[6], 0.000001);
    CuAssertDblEquals_Msg(tc, "f    ", FltSumsA[7], FltSumsB[7], 0.000001);
    // printf(" rho_sum: %f\tf_sum: %f\n", FltSumsB[0], FltSumsB[7]);
    rhoSum[2] = FltSumsB[0];
    uSum[2]   = FltSumsB[1];
}

///Clean up after test case
void cleanupTestGpuUpdate()
{
    cudaFree(AID); cudaFree(BID);
    cudaFree(AFD); cudaFree(BFD); cudaFree(CFD); cudaFree(DFD); cudaFree(EFD); cudaFree(FFD); cudaFree(GFD); cudaFree(HFD);
    free(AIH); free(BIH);
    free(AFH); free(BFH); free(CFH); free(DFH); free(EFH); free(FFH); free(GFH); free(HFH);
}

/**
 * @brief Unittest to compare results of #gpu_update_macro, #gpu_update_new, #gpuUpdateMacro
 *
 * @param tc test case
 * @test
 *   - compare results from previous runs
 */
void testCompareUpdateMacro(CuTest *tc)
{
    printBanner("Test compare gpu_update_* functions");
    CuAssertDblEquals_Msg(tc, "rho_macro != rho_new", rhoSum[0], rhoSum[1], 0.000001);
    CuAssertDblEquals_Msg(tc, "rho_macro != rho_new", rhoSum[0], rhoSum[2], 0.000001);
    CuAssertDblEquals_Msg(tc, "u_macro != u_new",     uSum[0],   uSum[1],   0.000001);
    CuAssertDblEquals_Msg(tc, "u_macro != u_new",     uSum[0],   uSum[2],   0.000001);
}

CuSuite *gpuMacroGetSuite()
{
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TCLN(suite, testGpuUpdateMacro_, cleanupTestGpuUpdateMacro);
    SUITE_ADD_TCLN(suite, testGpuUpdateNew,    cleanupTestGpuUpdate);
    SUITE_ADD_TCLN(suite, testGpuUpdateMacro,  cleanupTestGpuUpdate);
    SUITE_ADD_TEST(suite, testCompareUpdateMacro);
    return suite;
}