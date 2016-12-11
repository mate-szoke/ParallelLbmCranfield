/**
 * Unittests for the boundary conditions
 * @file TestGpuStream.cu
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
 * @brief Test to compare results from #gpuStreaming and #gpu_streaming
 *
 * @param tc test case
 * @test
 *  - Prepare boundary conditions for the lid driven cavity
 *  - Run the function
 *  - Check the sum of all the arrays between the two algorithms
 */
void testCompareGpuStream(CuTest *tc)
{
    printBanner("Test compare streaming");
    M = 1281;
    N = 1282;
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    dim3 tpb(THREADS);
    dim3 bpg1((int)(M*N/THREADS)+1);
    dim3 bpg9((int)(9*M*N/THREADS)+1);

    AIH = createHostArrayInt(M*N); //fluid
    createLidBcFluid(AIH, M, N);

    BIH = createHostArrayInt(9*M*N, ARRAY_FILL, 1); //stream9
    CIH = createHostArrayInt(8*M*N, ARRAY_FILL, 1); //stream8

    AFH = createHostArrayFlt(9*M*N, ARRAY_ZERO); //f1
    BFH = createHostArrayFlt(9*M*N, ARRAY_ZERO); //f2
    CFH = createHostArrayFlt(9*M*N, ARRAY_FILL, 0.0002); //fColl2

    AID = createGpuArrayInt(M*N, ARRAY_COPY, 0, AIH);
    BID = createGpuArrayInt(9*M*N, ARRAY_COPY, 0, BIH);
    CID = createGpuArrayInt(8*M*N, ARRAY_COPY, 0, CIH);
    AFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, BFH);
    CFD = createGpuArrayFlt(9*M*N, ARRAY_COPY, 0, CFH);
    DFD = createGpuArrayFlt(9*M*N);

    int sumAIa = sumHostInt(AIH, M*N);
    int sumBIa = sumHostInt(BIH, 9*M*N);
    int sumCIa = sumHostInt(CIH, 8*M*N);
    FLOAT_TYPE sumAFa = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFa = sumHostFlt(BFH, 9*M*N);
    FLOAT_TYPE sumCFa = sumHostFlt(CFH, 9*M*N);

    gpu_streaming<<<bpg9,tpb>>>(AID, BID, AFD, CFD);
    gpuStreaming<<<bpg1,tpb>>>(AID, CID, BFD, CFD);

    cudaMemcpy(AIH, AID, SIZEINT(N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(BIH, BID, SIZEINT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(CIH, CID, SIZEINT(8*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(AFH, AFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, SIZEFLT(9*N*M), cudaMemcpyDeviceToHost);

    int sumAIb = sumHostInt(AIH, M*N);
    int sumBIb = sumHostInt(BIH, 9*M*N);
    int sumCIb = sumHostInt(CIH, 8*M*N);
    FLOAT_TYPE sumAFb = sumHostFlt(AFH, 9*M*N);
    FLOAT_TYPE sumBFb = sumHostFlt(BFH, 9*M*N);
    FLOAT_TYPE sumCFb = sumHostFlt(CFH, 9*M*N);

    CuAssertIntEquals_Msg(tc, "fluid", sumAIa, sumAIb);
    CuAssertIntEquals_Msg(tc, "stream9", sumBIa, sumBIb);
    CuAssertIntEquals_Msg(tc, "stream8", sumCIa, sumCIb);
    CuAssertDblEquals_Msg(tc, "f", sumAFb, sumBFb, 0.00001);
    FLOAT_TYPE val = computeResidual(AFD, BFD, CFD, DFD, M, N);
    CuAssertDblEquals_Msg(tc, "fNorm", val, 0, 0.00001);
    CuAssertDblEquals_Msg(tc, "fColl", sumCFa, sumCFb, 0.00001);
}

///Clean up after test case
void cleanupTestCompareGpuStream()
{
    free(AIH);free(BIH);free(CIH);free(AFH);free(BFH);free(CFH);
    cudaFree(AID);cudaFree(BID);cudaFree(CID);cudaFree(AFD);cudaFree(BFD);cudaFree(CFD);cudaFree(DFD);
}

CuSuite* gpuStreamGetSuite()
{
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TCLN(suite, testCompareGpuStream, cleanupTestCompareGpuStream);
    return suite;
}