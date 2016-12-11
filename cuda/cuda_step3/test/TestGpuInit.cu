/**
 * Unit tests for the init functions
 * @file TestGpuInit.cu
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#include "CuTest.h"
#include "GpuFunctions.h"
#include "ArrayUtils.h"
#include "GpuConstants.h"
#include "TestUtils.h"
#include "Arguments.h"

static const FLOAT_TYPE rho = 5.0;
static const FLOAT_TYPE u = 0.05;
static const FLOAT_TYPE v = 0.0;
static const FLOAT_TYPE ySum = 2016; //it's a magic number, I am sorry (not sorry)
static const FLOAT_TYPE maxY = 0.9;
static const FLOAT_TYPE minY = 0.1;

/**
 * @brief Unittest for #gpu_init_1
 *
 * @param inlt inlet profile
 * @param u0,v0 input velocity
 * @test Unittest for #gpu_init_1
 *   - Prepare arrays
 *   - call function
 *   - check the sum of the arrays against predefined values
 *   - test for all inlet profiles
 */
void runTestGpuInit1(InletProfile inlt, FLOAT_TYPE u0, FLOAT_TYPE v0)
{
    M = 64;
    N = 64;
    cudaMemcpyToSymbol(width_d, &M, sizeof(int));
    cudaMemcpyToSymbol(height_d, &N, sizeof(int));
    cudaMemcpyToSymbol(inletProfile_d, &inlt, sizeof(int));

    cudaMemcpyToSymbol(rhoIn_d, &rho, sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(uIn_d, &u, sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(vIn_d, &v, sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(maxInletCoordY_d, &maxY, sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(minInletCoordY_d, &minY, sizeof(FLOAT_TYPE));

    dim3 tpb(THREADS); // THREADS/block
    dim3 bpg1((int)(M*N/THREADS)+1); // blocks/grid

    AFH = createHostArrayFlt(M*N, ARRAY_FILL, rho); //Rho
    BFH = createHostArrayFlt(M*N, ARRAY_FILL, u0); //Uo
    CFH = createHostArrayFlt(M*N, ARRAY_FILL, v0); //Vo
    DFH = createHostArrayFlt(M*N, ARRAY_ZERO); //U
    EFH = createHostArrayFlt(M*N, ARRAY_ZERO); //V
    FFH = createHostArrayFlt(M*N, ARRAY_ZERO); //Y
    GFH = createHostArrayFlt(M*N, ARRAY_ZERO); //X - temp

    fillLidCoordinates(GFH, FFH, M, N);
    free(GFH);

    AFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, AFH);
    BFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, BFH);
    CFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, CFH);
    DFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, DFH);
    EFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, EFH);
    FFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, FFH);
    cudaFree(DFD); cudaFree(EFD);

    gpu_init_1<<<bpg1,tpb>>>(AFD, BFD, CFD, DFD, EFD, FFD, M*N);
    DFD = createGpuArrayFlt(M*N, ARRAY_CPYD, 0, BFD);
    EFD = createGpuArrayFlt(M*N, ARRAY_CPYD, 0, CFD);

    cudaMemcpy(AFH, AFD, M*N*sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, M*N*sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, M*N*sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);
    cudaMemcpy(DFH, DFD, M*N*sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);
    cudaMemcpy(EFH, EFD, M*N*sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);
    cudaMemcpy(FFH, FFD, M*N*sizeof(FLOAT_TYPE), cudaMemcpyDeviceToHost);
}

///Test case for inlet @param tc test case
void testGpuInit1_inlet(CuTest *tc)
{
    printBanner("Test gpu_init_1 inlet");

    runTestGpuInit1(INLET, 0, 0);

    FLOAT_TYPE uSum = 0.0;
    FLOAT_TYPE inletLenghth2 = (maxY - minY) * (maxY - minY);
    int i;
    for (i=0; i<M*N; ++i)
    {
        uSum += 6 * u * (FFH[i] - minY) * (maxY - FFH[i]) / inletLenghth2;
    }

    CuAssertDblEquals_Msg(tc, "rho", M*N*rho, sumHostFlt(AFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "u0",  uSum,    sumHostFlt(BFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "v0",  M*N*v,   sumHostFlt(CFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "u",   uSum,    sumHostFlt(DFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "v",   M*N*v,   sumHostFlt(EFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "y",   ySum,    sumHostFlt(FFH, M*N), 0.01);
}

///Test case for no inlet @param tc test case
void testGpuInit1_noInlet(CuTest *tc)
{
    printBanner("Test gpu_init_1 no inlet");

    runTestGpuInit1(NO_INLET, u, v);

    CuAssertDblEquals_Msg(tc, "rho", M*N*rho, sumHostFlt(AFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "u0",  M*N*u,   sumHostFlt(BFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "v0",  M*N*v,   sumHostFlt(CFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "u",   M*N*u,   sumHostFlt(DFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "v",   M*N*v,   sumHostFlt(EFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "y",   ySum,    sumHostFlt(FFH, M*N), 0.01);
}

///Test case for pulsatile @param tc test case
void testGpuInit1_pulsatile(CuTest *tc)
{
    printBanner("Test gpu_init_1 pulsatile");

    runTestGpuInit1(PULSATILE_INLET, 0, 0);

    CuAssertDblEquals_Msg(tc, "rho", M*N*rho, sumHostFlt(AFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "u0",  0,       sumHostFlt(BFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "v0",  0,       sumHostFlt(CFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "u",   0,       sumHostFlt(DFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "v",   0,       sumHostFlt(EFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "y",   ySum,    sumHostFlt(FFH, M*N), 0.01);
}

///Clean up after test case
void cleanupTestGpuInit1()
{
    cudaFree(AFD); AFD = NULL;
    cudaFree(BFD); BFD = NULL;
    cudaFree(CFD); CFD = NULL;
    cudaFree(DFD); DFD = NULL;
    cudaFree(EFD); EFD = NULL;
    cudaFree(FFD); FFD = NULL;
    free(AFH); AFH = NULL;
    free(BFH); BFH = NULL;
    free(CFH); CFH = NULL;
    free(DFH); DFH = NULL;
    free(EFH); EFH = NULL;
    free(FFH); FFH = NULL;
}

/**
 * @brief Unittest for #gpuInitInletProfile
 *
 * @param tc test case
 * @test
 *   - Prepare arrays for inlet profile
 *   - call function
 *   - check the sum of the arrays against predefined values
 */
void testGpuInitInletProfile(CuTest *tc)
{
    printBanner("Test gpuInitInletProfile");
    M = 64;
    N = 64;
    cudaMemcpyToSymbol(uIn_d, &u, sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(vIn_d, &v, sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(maxInletCoordY_d, &maxY, sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(minInletCoordY_d, &minY, sizeof(FLOAT_TYPE));

    dim3 tpb(THREADS); // THREADS/block
    dim3 bpg1((int)(M*N/THREADS)+1); // blocks/grid

    AFH = createHostArrayFlt(M*N); //u0
    BFH = createHostArrayFlt(M*N); //v0
    CFH = createHostArrayFlt(M*N); //y
    DFH = createHostArrayFlt(M*N); //x not needed

    fillLidCoordinates(DFH, CFH, M, N);
    free(DFH);

    AFD = createGpuArrayFlt(M*N);
    BFD = createGpuArrayFlt(M*N);
    CFD = createGpuArrayFlt(M*N, ARRAY_COPY, 0, CFH);

    gpuInitInletProfile<<<bpg1,tpb>>>(AFD, BFD, CFD, M*N);

    FLOAT_TYPE uSum = 0.0;
    FLOAT_TYPE inletLenghth2 = (maxY - minY) * (maxY - minY);
    int i;
    for (i=0; i<M*N; ++i)
    {
        uSum += 6 * u * (CFH[i] - minY) * (maxY - CFH[i]) / inletLenghth2;
    }
    cudaMemcpy(AFH, AFD, SIZEFLT(M*N), cudaMemcpyDeviceToHost);
    cudaMemcpy(BFH, BFD, SIZEFLT(M*N), cudaMemcpyDeviceToHost);
    cudaMemcpy(CFH, CFD, SIZEFLT(M*N), cudaMemcpyDeviceToHost);

    CuAssertDblEquals_Msg(tc, "u0",  uSum,    sumHostFlt(AFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "v0",  M*N*v,   sumHostFlt(BFH, M*N), 0.01);
    CuAssertDblEquals_Msg(tc, "y",   ySum,    sumHostFlt(CFH, M*N), 0.01);
}

///Clean up after test case
void cleanupTestGpuInitInletProfile()
{
    cudaFree(AFD); cudaFree(BFD); cudaFree(CFD);
    free(AFH); free(BFH); free(CFH);
}

CuSuite *gpuInitGetSuite()
{
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TCLN(suite, testGpuInit1_inlet,      cleanupTestGpuInit1);
    SUITE_ADD_TCLN(suite, testGpuInit1_noInlet,    cleanupTestGpuInit1);
    SUITE_ADD_TCLN(suite, testGpuInit1_pulsatile,  cleanupTestGpuInit1);
    SUITE_ADD_TCLN(suite, testGpuInitInletProfile, cleanupTestGpuInitInletProfile);
    return suite;
}