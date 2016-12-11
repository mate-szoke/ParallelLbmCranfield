/**
 * Unittests for the solver
 * @file TestIterate.cu
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 * @warning This test needs actual meshes, so put them into the test directory or give their path in
 * the code.
 * @todo create a command-line argument for test meshes
 */
#include "GpuFunctions.h"       // GPU kernels
#include "ShellFunctions.h"     // For convenience
#include "FilesReading.h"       // For reading files
#include "CellFunctions.h"      // For cell modifications
#include "ComputeResiduals.h"   // residuals
#include "LogWriter.h"
#include "Iterate.h"
#include "ArrayUtils.h"
#include "Check.h"
#include "TestUtils.h"
#include "CuTest.h"


static int  *nodeIdX, *nodeIdY, *nodeType, *bcNodeIdX, *bcNodeIdY, *latticeId, *bcType, *bcBoundId;
static FLOAT_TYPE *nodeX, *nodeY,*bcX, *bcY;
static int *fluid_d, *bcNodeIdX_d, *bcNodeIdY_d, *latticeId_d, *bcType_d, *bcBoundId_d;
static FLOAT_TYPE *coordX_d, *coordY_d, *bcX_d, *bcY_d;
static FLOAT_TYPE *normi, *dragSum, *liftSum;
static FLOAT_TYPE *u, *v, *rho, *u_d, *v_d, *rho_d, *u0_d, *v0_d;
static FLOAT_TYPE *drag_d, *lift_d, *f_d, *fColl_d, *temp9a_d, *temp9b_d, *tempA_d, *tempB_d;
static int *mask, *bcMask, *bcIdx, *stream, *stream_d, *bcIdxCollapsed_d, *bcMaskCollapsed_d, *bcMask_d, *bcIdx_d;
static FLOAT_TYPE *q, *q_d, *qCollapsed_d;
static FLOAT_TYPE *qOld, *qOld_d, *rhoOld_d, *uOld_d, *vOld_d, *u0Old_d, *v0Old_d, *liftOld_d, *dragOld_d;
static int *bcId, *bcId_d, *boundary_d, *boundaryId, *boundaryId_d, *corner_d, *streamOld_d, *streamOld;
static FLOAT_TYPE *m_d, *mEq_d, *coll_d, *fEq_d, *f2_d, *fc2_d, *fN_d;
static int *cornerCheck_d, *fluidCheck_d, *bcIdCheck_d, *boundaryCheck_d;

/**
 * @brief Test the solver by having all of its steps done the old and the new way
 *
 * @param args input arguments
 * @param inFn input files
 * @param tc test case
 */
void runTestIteration(Arguments *args, InputFilenames *inFn, CuTest *tc)
{
  FLOAT_TYPE delta;          // grid spacing
  int n,m;                   // number of nodes in the x and y directions
  FLOAT_TYPE maxInletCoordY; // maximum inlet coordinate in y
  FLOAT_TYPE minInletCoordY; // minimum inlet coordinate in y
  int numInletNodes;         // number of inlet nodes

  int numNodes = readNodeFile(inFn->node, &nodeIdX, &nodeIdY, &nodeX, &nodeY, &nodeType);

  fluid_d  = createGpuArrayInt(numNodes, ARRAY_COPY, 0, nodeType);
  coordX_d = createGpuArrayFlt(numNodes, ARRAY_COPY, 0., nodeX);
  coordY_d = createGpuArrayFlt(numNodes, ARRAY_COPY, 0., nodeY);

  int numConns = readConnFile(inFn->bc, &bcNodeIdX, &bcNodeIdY, &latticeId, &bcType, &bcX, &bcY, &bcBoundId);

  bcNodeIdX_d  = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcNodeIdX);
  bcNodeIdY_d  = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcNodeIdX);
  latticeId_d  = createGpuArrayInt(numConns, ARRAY_COPY, 0, latticeId);
  bcType_d     = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcType);
  bcBoundId_d  = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcBoundId);
  bcX_d        = createGpuArrayFlt(numConns, ARRAY_COPY, 0., bcX);
  bcY_d        = createGpuArrayFlt(numConns, ARRAY_COPY, 0., bcY);

  m = getLastValue(nodeIdY, numNodes);
  n = getLastValue(nodeIdX, numNodes);
  delta = getGridSpacing(nodeIdX, nodeIdY, nodeX, numNodes);
  numInletNodes  = getNumInletNodes(bcType, latticeId, numConns);
  maxInletCoordY = getMaxInletCoordY(bcType, latticeId, bcY, delta, numConns);
  minInletCoordY = getMinInletCoordY(bcType, latticeId, bcY, delta, numConns);

  initConstants(args, maxInletCoordY, minInletCoordY, delta, m, n);
  dim3 tpb(THREADS); // THREADS/block
  dim3 bpg1((int)(  m*n/THREADS)+1);     // blocks/grid  MxN
  dim3 bpg8((int)(8*m*n/THREADS)+1);     // blocks/grid 8MxN
  dim3 bpg9((int)(9*m*n/THREADS)+1);     // blocks/grid 9MxN
  dim3 bpgBC((int)(numConns/THREADS)+1); // blocks/grid N_BC

  // residuals
  normi    = createHostArrayFlt(args->iterations, ARRAY_ZERO);
  dragSum = createHostArrayFlt(args->iterations, ARRAY_ZERO);
  liftSum = createHostArrayFlt(args->iterations, ARRAY_ZERO);

  u   = createHostArrayFlt(m*n, ARRAY_ZERO);
  v   = createHostArrayFlt(m*n, ARRAY_ZERO);
  rho = createHostArrayFlt(m*n, ARRAY_ZERO);

  rho_d  = createGpuArrayFlt(m*n, ARRAY_FILL, args->rho);

  if (args->inletProfile == NO_INLET)
  {
    u0_d = createGpuArrayFlt(m*n, ARRAY_FILL, args->u);
    v0_d = createGpuArrayFlt(m*n, ARRAY_FILL, args->v);
  }
  else
  {
    u0_d = createGpuArrayFlt(m*n, ARRAY_ZERO);
    v0_d = createGpuArrayFlt(m*n, ARRAY_ZERO);
  }
  if (args->inletProfile == INLET)
  {
    gpuInitInletProfile<<<bpg1,tpb>>>(u0_d, v0_d, coordY_d, m*n);
  }

  drag_d = createGpuArrayFlt(m*n, ARRAY_ZERO);
  lift_d = createGpuArrayFlt(m*n, ARRAY_ZERO);

  f_d     = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  fColl_d = createGpuArrayFlt(9*m*n, ARRAY_ZERO);

  temp9a_d = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  temp9b_d = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  tempA_d  = createGpuArrayFlt(  m*n, ARRAY_ZERO);
  tempB_d  = createGpuArrayFlt(  m*n, ARRAY_ZERO);

  mask   = createHostArrayInt(m*n, ARRAY_ZERO);
  bcMask = createHostArrayInt(m*n, ARRAY_ZERO);
  bcIdx  = createHostArrayInt(m*n, ARRAY_ZERO);

  u_d = createGpuArrayFlt(m*n, ARRAY_CPYD, 0, u0_d);
  v_d = createGpuArrayFlt(m*n, ARRAY_CPYD, 0, v0_d);

  stream = createHostArrayInt(8*m*n, ARRAY_FILL, 1);
  q      = createHostArrayFlt(8*m*n, ARRAY_FILL, 0.5);

  int bcCount = initBoundaryConditions(bcNodeIdX, bcNodeIdY, q, bcBoundId, nodeType, bcX, bcY, nodeX, nodeY,
                                       latticeId, stream, bcType, bcMask, bcIdx, mask,
                                       delta, m, n, numConns);

  bcIdxCollapsed_d    = createGpuArrayInt(bcCount,   ARRAY_ZERO);
  bcMaskCollapsed_d   = createGpuArrayInt(bcCount,   ARRAY_ZERO);
  qCollapsed_d        = createGpuArrayFlt(8*bcCount, ARRAY_ZERO);

  dim3 bpgB((int)(bcCount/THREADS)+1); // blocks/grid
  bcMask_d = createGpuArrayInt(m*n, ARRAY_COPY, 0, bcMask);
  bcIdx_d  = createGpuArrayInt(m*n, ARRAY_COPY, 0, bcIdx);

  collapseBc(bcIdx, bcIdxCollapsed_d, bcMask, bcMaskCollapsed_d, q, qCollapsed_d, mask, m, n, bcCount);

  stream_d   = createGpuArrayInt(8*m*n, ARRAY_COPY, 0, stream);
  q_d = createGpuArrayFlt(8*m*n, ARRAY_COPY, 0, q);

  /************************************************************************************************/
  //old init
  FLOAT_TYPE qLat[9]={0.,1.,1.,1.,1.,sqrt(2),sqrt(2),sqrt(2),sqrt(2)};
  qOld       = createHostArrayFlt(9*m*n);
  bcId       = createHostArrayInt(9*m*n, ARRAY_ZERO);
  streamOld  = createHostArrayInt(9*m*n);
  boundaryId = createHostArrayInt(  m*n);
  int ind, ind9, i;
  for(i=0;i<numConns;i++)
  {
    ind = bcNodeIdX[i] + bcNodeIdY[i] * n;
    ind9 = ind + m*n * latticeId[i];
    bcId[ind9] = bcType[i];
    boundaryId[ind] = bcBoundId[i];
    qOld[ind9] = sqrt(pow( bcX[i]-nodeX[ind],2 ) + pow( bcY[i]-nodeY[ind],2) ) / (delta * qLat[latticeId[i]]);
  }
  bcId_d       = createGpuArrayInt(9*m*n, ARRAY_COPY, 0, bcId);
  boundaryId_d = createGpuArrayInt(  m*n, ARRAY_COPY, 0, boundaryId);
  qOld_d       = createGpuArrayFlt(9*m*n, ARRAY_COPY, 0, qOld);

  corner_d     = createGpuArrayInt(  m*n, ARRAY_ZERO);
  rhoOld_d     = createGpuArrayFlt(  m*n);
  boundary_d   = createGpuArrayInt(  m*n, ARRAY_ZERO);
  streamOld_d  = createGpuArrayInt(9*m*n);
  uOld_d       = createGpuArrayFlt(  m*n);
  vOld_d       = createGpuArrayFlt(  m*n);
  u0Old_d      = createGpuArrayFlt(  m*n);
  v0Old_d      = createGpuArrayFlt(  m*n);
  liftOld_d    = createGpuArrayFlt(  m*n);
  dragOld_d    = createGpuArrayFlt(  m*n);

  gpu_init<<<bpg9,tpb>>>(corner_d, rhoOld_d, qOld_d, boundary_d, coordY_d, streamOld_d, bcId_d, uOld_d, vOld_d, u0Old_d, v0Old_d);

  CHECK(cudaMemcpy(streamOld, streamOld_d, SIZEINT(9*m*n), cudaMemcpyDeviceToHost));
  CuAssertIntEquals_Msg(tc, "stream", 0, compareArraysInt(streamOld+m*n, stream, 8*m*n));
  CuAssertDblEquals_Msg(tc, "u0",     0, compareGpuArraySumFlt(rho_d, rhoOld_d, tempA_d, tempB_d, m*n), 1e-12);
  CuAssertDblEquals_Msg(tc, "u0",     0, compareGpuArraySumFlt(u0_d,  u0Old_d,  tempA_d, tempB_d, m*n), 1e-12);
  CuAssertDblEquals_Msg(tc, "u0",     0, compareGpuArraySumFlt(v0_d,  v0Old_d,  tempA_d, tempB_d, m*n), 1e-12);

  //check bc
  cornerCheck_d   = createGpuArrayInt(  m*n, ARRAY_ZERO);
  boundaryCheck_d = createGpuArrayInt(  m*n, ARRAY_ZERO);
  fluidCheck_d    = createGpuArrayInt(  m*n, ARRAY_CPYD, 0, fluid_d);
  bcIdCheck_d     = createGpuArrayInt(9*m*n, ARRAY_ZERO);
  gpu_convert<<<bpgB,tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, fluidCheck_d, boundaryCheck_d, cornerCheck_d, bcIdCheck_d, bcCount);
  CuAssertIntEquals_Msg(tc, "corner1",  0, compareGpuArrayInt(corner_d, cornerCheck_d, m*n));
  CuAssertDblEquals_Msg(tc, "corner2",  0, compareGpuArraySumInt(corner_d,   cornerCheck_d, tempA_d,  tempB_d,    m*n), 1e-12);
  CuAssertIntEquals_Msg(tc, "boundary1", 0, compareGpuArrayInt(boundary_d, boundaryCheck_d, m*n));
  CuAssertDblEquals_Msg(tc, "boundary2", 0, compareGpuArraySumInt(boundary_d, boundaryCheck_d, tempA_d,  tempB_d,    m*n), 1e-12);
  CuAssertIntEquals_Msg(tc, "fluid1",    0, compareGpuArrayInt(fluid_d, fluidCheck_d, m*n));
  CuAssertDblEquals_Msg(tc, "fluid2",    0, compareGpuArraySumInt(fluid_d,    fluidCheck_d,    tempA_d,  tempB_d,    m*n), 1e-12);
  CuAssertDblEquals_Msg(tc, "bcid",      0, compareGpuArraySumInt(bcId_d,     bcIdCheck_d,     temp9a_d, temp9b_d, 9*m*n), 1e-12);

  m_d    = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  mEq_d  = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  coll_d = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  f2_d   = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  fc2_d  = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  fEq_d  = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  fN_d   = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  /************************************************************************************************/

  FLOAT_TYPE val;

  int iter = 0;
  while (iter < args->iterations)
  {
    printf("\niter%d\n", iter);
    switch(args->collisionModel)
    {

      case BGKW:
        CHECK(cudaMemcpy(f2_d, f_d, SIZEFLT(9*m*n), cudaMemcpyDeviceToDevice));
        gpuCollBgkw<<<bpg1,tpb>>>(fluid_d, rho_d, u_d, v_d, f_d, fColl_d);
        gpu_bgk<<<bpg9,tpb>>>(fluid_d, fEq_d, rhoOld_d, uOld_d, vOld_d, fc2_d, f2_d);
        val = compareGpuArraySumFlt(fColl_d, fc2_d, temp9a_d, temp9b_d, 9*m*n);
        CuAssertDblEquals_Msg(tc, "fCollBgkw", 0, val, 1e-12);
      break;

      case TRT:
        CHECK(cudaMemcpy(f2_d, f_d, SIZEFLT(9*m*n), cudaMemcpyDeviceToDevice));
        gpuCollTrt<<<bpg1,tpb>>>(fluid_d, rho_d, u_d, v_d, f_d, fColl_d);
        gpu_trt1<<<bpg9,tpb>>>(fluid_d, fEq_d, rhoOld_d, uOld_d, vOld_d);
        gpu_trt2<<<bpg9,tpb>>>(fluid_d, fEq_d, f2_d, fc2_d);
        val = compareGpuArraySumFlt(fColl_d, fc2_d, temp9a_d, temp9b_d, 9*m*n);
        CuAssertDblEquals_Msg(tc, "fCllTrt", 0, val, 1e-12);
      break;

      case MRT:
        CHECK(cudaMemcpy(f2_d, f_d, SIZEFLT(9*m*n), cudaMemcpyDeviceToDevice));
        gpuCollMrt<<<bpg1,tpb>>>(fluid_d, rho_d, u_d, v_d, f_d, fColl_d);
        gpu_mrt1<<<bpg9,tpb>>>(fluid_d, rho_d, u_d, v_d, f2_d, mEq_d, m_d);
        gpu_mrt2<<<bpg9,tpb>>>(fluid_d, coll_d, m_d, mEq_d, fc2_d, f2_d);
        val = compareGpuArraySumFlt(fColl_d, fc2_d, temp9a_d, temp9b_d, 9*m*n);
        CuAssertDblEquals_Msg(tc, "fCollMrt", 0, val, 1e-8);
      break;
    }
    gpuStreaming<<<bpg1,tpb>>>(fluid_d, stream_d, f_d, fColl_d);
    gpu_streaming<<<bpg9,tpb>>>(fluid_d, streamOld_d+m*n, f2_d, fc2_d);
    val = compareGpuArraySumFlt(f_d, f2_d, temp9a_d, temp9b_d, 9*m*n);
    CuAssertDblEquals_Msg(tc, "fstream", 0, val, 1e-12);

    gpuBcInlet <<<bpgB,tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d, u0_d, v0_d, bcCount);
    gpu_boundaries1<<<bpg9,tpb>>>(fluid_d, boundary_d, bcId_d, f2_d, u0Old_d, v0Old_d, corner_d);
    val = compareGpuArraySumFlt(f_d, f2_d, temp9a_d, temp9b_d, 9*m*n);
    CuAssertDblEquals_Msg(tc, "finlet", 0, val, 1e-12);

    gpuBcWall  <<<bpgB,tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d, fColl_d, qCollapsed_d, bcCount);
    gpu_boundaries2<<<bpg9,tpb>>>(fluid_d, fN_d, fc2_d, bcId_d, qOld_d, f2_d);
    CuAssertIntEquals_Msg(tc, "fwall", 0, compareGpuArrayFlt(f_d, f2_d, 9*m*n));
    val = compareGpuArraySumFlt(f_d, f2_d, temp9a_d, temp9b_d, 9*m*n);
    CuAssertDblEquals_Msg(tc, "fwall", 0, val, 1e-12);

    gpuBcOutlet<<<bpgB,tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d, u0_d, v0_d, bcCount);
    gpu_boundaries3<<<bpg9,tpb>>>(fluid_d, bcId_d, f2_d, u0Old_d, v0Old_d, corner_d);
    val = compareGpuArraySumFlt(f_d, f2_d, temp9a_d, temp9b_d, 9*m*n);
    CuAssertDblEquals_Msg(tc, "foutlet", 0, val, 1e-6);

    gpuUpdateMacro<<<bpg1,tpb>>>(fluid_d, rho_d, u_d, v_d, bcMask_d, drag_d, lift_d, coordX_d, coordY_d, f_d);
    gpu_update_macro<<<bpg9,tpb>>>(fluid_d, rhoOld_d, uOld_d, vOld_d, bcId_d, boundaryId_d, dragOld_d, liftOld_d, coordX_d, coordY_d, f2_d);
    val = compareGpuArraySumFlt(rho_d, rhoOld_d, tempA_d, tempB_d, m*n);
    CuAssertDblEquals_Msg(tc, "rhomacro", 0, val, 1e-12);
    val = compareGpuArraySumFlt(u_d, uOld_d, tempA_d, tempB_d, m*n);
    CuAssertDblEquals_Msg(tc, "umacro", 0, val, 1e-7);
    val = compareGpuArraySumFlt(v_d, vOld_d, tempA_d, tempB_d, m*n);
    CuAssertDblEquals_Msg(tc, "vmacro", 0, val, 1e-12);

    FLOAT_TYPE r = computeResidual(f_d, fColl_d, temp9a_d, temp9b_d, m, n);
    if (r != r)
    {
      CuFail_Line(tc, __FILE__, __LINE__, NULL, "divergence");
    }
    normi[iter] = r;
    if (args->boundaryId > 0)
    {
      dragSum[iter] = computeDragLift(bcMask_d, drag_d, tempA_d, tempB_d, m, n, args->boundaryId);
      liftSum[iter] = computeDragLift(bcMask_d, lift_d, tempA_d, tempB_d, m, n, args->boundaryId);
    }

    iter++; // update loop variable
  }
}

///test for lid-driven cavity, using MRT model @param tc test case
void testIterationLidMrt(CuTest *tc)
{
  printBanner("Test iteration with MRT on lid-driven cavity");
  Arguments args;
  args.u         = 0.05;
  args.v         = 0.;
  args.rho       = 5;
  args.viscosity = 0.01;
  args.inletProfile   = NO_INLET;
  args.outletProfile  = OUTLET;
  args.collisionModel = MRT;
  args.boundaryType   = STRAIGHT;
  args.outputFormat   = PARAVIEW;
  args.iterations    = 2;
  args.autosaveEvery = 0;
  args.autosaveAfter = 0;
  args.boundaryId    = 0;

  InputFilenames inFn;
  strcpy(inFn.node  , "test/lidnode.256");
  strcpy(inFn.bc    , "test/lidconn.256");
  strcpy(inFn.init  , "SetUpData.ini");
  strcpy(inFn.result, "Results/");
  strcpy(inFn.final , "Results/FinalData.csv");

  runTestIteration(&args, &inFn, tc);
}

///test for expansion, using BGKW model @param tc test case
void testIterationExpBgkw(CuTest *tc)
{
  printBanner("Test iteration with BGKW on expansion");
  Arguments args;
  args.u         = 0.05;
  args.v         = 0.;
  args.rho       = 5;
  args.viscosity = 0.01;
  args.inletProfile   = NO_INLET;
  args.outletProfile  = OUTLET_FIRST;
  args.collisionModel = BGKW;
  args.boundaryType   = STRAIGHT;
  args.outputFormat   = PARAVIEW;
  args.iterations    = 2;
  args.autosaveEvery = 0;
  args.autosaveAfter = 0;
  args.boundaryId    = 0;

  InputFilenames inFn;
  strcpy(inFn.node  , "test/expnode.256");
  strcpy(inFn.bc    , "test/expconn.256");
  strcpy(inFn.init  , "SetUpData.ini");
  strcpy(inFn.result, "Results/");
  strcpy(inFn.final , "Results/FinalData.csv");

  runTestIteration(&args, &inFn, tc);
}

///test for expansion, using MRT model @param tc test case
void testIterationExpMrt(CuTest *tc)
{
  printBanner("Test iteration with MRT on expansion");
  Arguments args;
  args.u         = 0.05;
  args.v         = 0.;
  args.rho       = 5;
  args.viscosity = 0.01;
  args.inletProfile   = NO_INLET;
  args.outletProfile  = OUTLET_FIRST;
  args.collisionModel = MRT;
  args.boundaryType   = STRAIGHT;
  args.outputFormat   = PARAVIEW;
  args.iterations    = 2;
  args.autosaveEvery = 0;
  args.autosaveAfter = 0;
  args.boundaryId    = 0;

  InputFilenames inFn;
  strcpy(inFn.node  , "test/expnode.256");
  strcpy(inFn.bc    , "test/expconn.256");
  strcpy(inFn.init  , "SetUpData.ini");
  strcpy(inFn.result, "Results/");
  strcpy(inFn.final , "Results/FinalData.csv");

  runTestIteration(&args, &inFn, tc);
}

///Clean up after test case
void cleanupTestIteration()
{
  void *hostArrays[] = {nodeIdX, nodeIdY, nodeX, nodeY, nodeType, bcNodeIdX, bcNodeIdY, latticeId,
                        bcType, bcX, bcY, bcBoundId, u, v, rho, mask, bcMask, bcIdx, stream, q,
                        normi, dragSum, liftSum, qOld, bcId, boundaryId, streamOld};
  void *gpuArrays[] = { coordX_d, coordY_d, fluid_d, bcNodeIdX_d, bcNodeIdY_d, latticeId_d,
                        bcType_d, bcX_d, bcY_d, bcBoundId_d, u_d, v_d, rho_d, u0_d, v0_d, drag_d, lift_d,
                        f_d, fColl_d, temp9a_d, temp9b_d, tempA_d, tempB_d,
                        bcMask_d, bcMaskCollapsed_d, bcIdx_d, bcIdxCollapsed_d, stream_d, q_d,
                        qCollapsed_d, m_d, mEq_d, coll_d, f2_d, fc2_d,
                        bcId_d, boundaryId_d, qOld_d, corner_d, rhoOld_d, boundary_d, streamOld_d,
                        uOld_d, vOld_d, u0Old_d, v0Old_d, liftOld_d, dragOld_d, fEq_d, fN_d,
                        cornerCheck_d, fluidCheck_d, bcIdCheck_d, boundaryCheck_d
                      };
  freeAllHost(hostArrays, sizeof(hostArrays)/sizeof(hostArrays[0]));
  freeAllGpu (gpuArrays,  sizeof(gpuArrays) /sizeof(gpuArrays[0]) );
}

CuSuite* iterateGetSuite()
{
  CuSuite *suite = CuSuiteNew();
  SUITE_ADD_TCLN(suite, testIterationLidMrt,  cleanupTestIteration);
  SUITE_ADD_TCLN(suite, testIterationExpBgkw, cleanupTestIteration);
  SUITE_ADD_TCLN(suite, testIterationExpMrt,  cleanupTestIteration);
  return suite;
}
