#include <stdio.h>                      // printf();
#include <math.h>                       // need to compile with -lm
#include <stdlib.h>                     // for calloc();
#include <stdbool.h>                    // Include for bool type variables!
#include <string.h>                     // String operations
#include <time.h>                       // time functions
#include <errno.h>
#include "GpuFunctions.h"       // GPU kernels
#include "ShellFunctions.h"     // For convenience
#include "FilesReading.h"       // For reading files
#include "FilesWriting.h"       // For writing files e.g. tecplot
#include "CellFunctions.h"      // For cell modifications
#include "ComputeResiduals.h"   // residuals
#include "LogWriter.h"
#include "Iterate.h"
#include "ArrayUtils.h"
#include "Check.h"
#include "TestUtils.h"

int Iteration(InputFilenames *inFn, Arguments *args)
{
  // Time measurement: declaration, begin
  clock_t tStart = clock();

  FILE* logFile;               // file for log
  char autosaveFilename[768];  // autosave filename
  char outputFilename[768];    // initial data will be written to this file
  char finalFilename[768];     // final data will be written to this file
  char logFilename[768];       // path of the .log file
  char residualsFilename[768]; // path of the residuals file
  char timeFilename[768];      // path of time measurement file

  logFilename[0]       = '\0';
  residualsFilename[0] = '\0';
  timeFilename[0]      = '\0';

  if (strlen(inFn->result))
  {
    strcat(logFilename, inFn->result);
    strcat(residualsFilename, inFn->result);
    strcat(timeFilename, inFn->result);
  }
  strcat(logFilename,       "lbmsolver.log");
  strcat(residualsFilename, "residuals.dat");
  strcat(timeFilename,      "runtimes.dat");

  int autosaveIt = 1;        // autosave i variable, will be incremented after every autosave
  int numNodes,numConns;     // This will store the number of lines of the read files
  FLOAT_TYPE delta;          // grid spacing
  int n,m;                   // number of nodes in the x and y directions
  FLOAT_TYPE maxInletCoordY; // maximum inlet coordinate in y
  FLOAT_TYPE minInletCoordY; // minimum inlet coordinate in y
  int numInletNodes;         // number of inlet nodes

  int *nodeIdX, *nodeIdY, *nodeType, *bcNodeIdX, *bcNodeIdY, *latticeId, *bcType, *bcBoundId;
  FLOAT_TYPE *nodeX, *nodeY,*bcX, *bcY;

  FLOAT_TYPE taskTime[9];
  int i;
  for (i=0; i<9; ++i)
  {
    taskTime[i] = 0.0;
  }

  clock_t tInstant1, tInstant2; // Time measurement points, universal
  clock_t tIterStart, tIterEnd; // Time measurement points: main loop

  // cuda time measurement variables
  cudaEvent_t start, stop;
  float cudatime;
  CHECK(cudaEventCreate(&start));
  CHECK(cudaEventCreate(&stop));

  numNodes = readNodeFile(inFn->node, &nodeIdX, &nodeIdY, &nodeX, &nodeY, &nodeType);
  if (numNodes == 0)
  {
    return 2;
  }

  int *fluid_d = createGpuArrayInt(numNodes, ARRAY_COPY, 0, nodeType);
  FLOAT_TYPE *coordX_d = createGpuArrayFlt(numNodes, ARRAY_COPY, 0., nodeX);
  FLOAT_TYPE *coordY_d = createGpuArrayFlt(numNodes, ARRAY_COPY, 0., nodeY);

  numConns = readConnFile(inFn->bc, &bcNodeIdX, &bcNodeIdY, &latticeId, &bcType, &bcX, &bcY, &bcBoundId);
  if (numConns == 0)
  {
    return 2;
  }

  int *bcNodeIdX_d  = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcNodeIdX);
  int *bcNodeIdY_d  = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcNodeIdX);
  int *latticeId_d  = createGpuArrayInt(numConns, ARRAY_COPY, 0, latticeId);
  int *bcType_d     = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcType);
  int *bcBoundId_d  = createGpuArrayInt(numConns, ARRAY_COPY, 0, bcBoundId);
  FLOAT_TYPE *bcX_d = createGpuArrayFlt(numConns, ARRAY_COPY, 0., bcX);
  FLOAT_TYPE *bcY_d = createGpuArrayFlt(numConns, ARRAY_COPY, 0., bcY);

  m = getLastValue(nodeIdY, numNodes);
  n = getLastValue(nodeIdX, numNodes);
  delta = getGridSpacing(nodeIdX, nodeIdY, nodeX, numNodes);
  numInletNodes  = getNumInletNodes(bcType, latticeId, numConns);
  maxInletCoordY = getMaxInletCoordY(bcType, latticeId, bcY, delta, numConns);
  minInletCoordY = getMinInletCoordY(bcType, latticeId, bcY, delta, numConns);

  writeInitLog(logFilename, args, delta, m, n, numInletNodes, maxInletCoordY, minInletCoordY);
  logFile = fopen(logFilename, "a");

  // In case of no autosave
  sprintf(autosaveFilename, "NOWHERE!");

  initConstants(args, maxInletCoordY, minInletCoordY, delta, m, n);

  dim3 tpb(THREADS); // THREADS/block
  dim3 bpg1((int)(  m*n/THREADS)+1);     // blocks/grid  MxN
  dim3 bpg8((int)(8*m*n/THREADS)+1);     // blocks/grid 8MxN
  dim3 bpg9((int)(9*m*n/THREADS)+1);     // blocks/grid 9MxN
  dim3 bpgBC((int)(numConns/THREADS)+1); // blocks/grid N_BC

  // residuals
  FLOAT_TYPE *norm    = createHostArrayFlt(args->iterations, ARRAY_ZERO);
  FLOAT_TYPE *dragSum = createHostArrayFlt(args->iterations, ARRAY_ZERO);
  FLOAT_TYPE *liftSum = createHostArrayFlt(args->iterations, ARRAY_ZERO);

  fprintf(logFile,"\n:::: Initializing ::::\n");
  printf("\n:::: Initializing ::::\n");
  CHECK(cudaEventRecord(start, 0));

  FLOAT_TYPE *u   = createHostArrayFlt(m*n, ARRAY_ZERO);
  FLOAT_TYPE *v   = createHostArrayFlt(m*n, ARRAY_ZERO);
  FLOAT_TYPE *rho = createHostArrayFlt(m*n, ARRAY_ZERO);

  FLOAT_TYPE *rho_d  = createGpuArrayFlt(m*n, ARRAY_FILL, args->rho);

  FLOAT_TYPE *u0_d, *v0_d;
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

  FLOAT_TYPE *drag_d = createGpuArrayFlt(m*n, ARRAY_ZERO);
  FLOAT_TYPE *lift_d = createGpuArrayFlt(m*n, ARRAY_ZERO);

  FLOAT_TYPE *f_d     = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  FLOAT_TYPE *fColl_d = createGpuArrayFlt(9*m*n, ARRAY_ZERO);

  FLOAT_TYPE *temp9a_d = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  FLOAT_TYPE *temp9b_d = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  FLOAT_TYPE *tempA_d  = createGpuArrayFlt(  m*n, ARRAY_ZERO);
  FLOAT_TYPE *tempB_d  = createGpuArrayFlt(  m*n, ARRAY_ZERO);

  int *mask   = createHostArrayInt(m*n, ARRAY_ZERO);
  int *bcMask = createHostArrayInt(m*n, ARRAY_ZERO);
  int *bcIdx  = createHostArrayInt(m*n, ARRAY_ZERO);

  FLOAT_TYPE *u_d = createGpuArrayFlt(m*n, ARRAY_CPYD, 0, u0_d);
  FLOAT_TYPE *v_d = createGpuArrayFlt(m*n, ARRAY_CPYD, 0, v0_d);

  int *stream     = createHostArrayInt(8*m*n, ARRAY_FILL, 1);
  FLOAT_TYPE *q   = createHostArrayFlt(8*m*n, ARRAY_FILL, 0.5);

  int bcCount = initBoundaryConditions(bcNodeIdX, bcNodeIdY, q, bcBoundId, nodeType, bcX, bcY, nodeX, nodeY,
                                       latticeId, stream, bcType, bcMask, bcIdx, mask,
                                       delta, m, n, numConns);

  int *bcIdxCollapsed_d    = createGpuArrayInt(bcCount,   ARRAY_ZERO);
  int *bcMaskCollapsed_d   = createGpuArrayInt(bcCount,   ARRAY_ZERO);
  FLOAT_TYPE *qCollapsed_d = createGpuArrayFlt(8*bcCount, ARRAY_ZERO);

  dim3 bpgB((int)(bcCount/THREADS)+1); // blocks/grid
  int *bcMask_d = createGpuArrayInt(m*n, ARRAY_COPY, 0, bcMask);
  int *bcIdx_d  = createGpuArrayInt(m*n, ARRAY_COPY, 0, bcIdx);

  collapseBc(bcIdx, bcIdxCollapsed_d, bcMask, bcMaskCollapsed_d, q, qCollapsed_d, mask, m, n, bcCount);

  int *stream_d   = createGpuArrayInt(8*m*n, ARRAY_COPY, 0, stream);
  FLOAT_TYPE *q_d = createGpuArrayFlt(8*m*n, ARRAY_COPY, 0, q);

  CHECK(cudaMemcpy(u,   u_d,   SIZEFLT(m*n), cudaMemcpyDeviceToHost));
  CHECK(cudaMemcpy(v,   v_d,   SIZEFLT(m*n), cudaMemcpyDeviceToHost));
  CHECK(cudaMemcpy(rho, rho_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));

  CHECK(cudaEventRecord(stop, 0));
  CHECK(cudaEventSynchronize(stop));
  CHECK(cudaEventElapsedTime(&cudatime, start, stop));
  taskTime[T_INIT] += cudatime/1000;

  fclose(logFile);
  writeNodeNumbers(logFilename, numNodes, numConns, bcCount);
  logFile = fopen(logFilename, "a");

  void *hostArrays[] = {nodeIdX, nodeIdY, nodeX, nodeY, nodeType, bcNodeIdX, bcNodeIdY, latticeId,
                        bcType, bcX, bcY, bcBoundId, u, v, rho, mask, bcMask, bcIdx, stream, q,
                        norm, dragSum, liftSum};
  void *gpuArrays[]  = {coordX_d, coordY_d, fluid_d, bcNodeIdX_d, bcNodeIdY_d, latticeId_d,
                        bcType_d, bcX_d, bcY_d, bcBoundId_d, u_d, v_d, rho_d, u0_d, v0_d, drag_d, lift_d,
                        f_d, fColl_d, temp9a_d, temp9b_d, tempA_d, tempB_d,
                        bcMask_d, bcMaskCollapsed_d, bcIdx_d, bcIdxCollapsed_d, stream_d, q_d,
                        qCollapsed_d};

  fprintf(logFile,"\n:::: Initialization done! ::::\n");

  printf("Initialization took %f seconds\n",         taskTime[T_INIT]);


  // Write Initialized data
  switch(args->outputFormat)
  {
    case PARAVIEW: sprintf(outputFilename, "%sInitialData.csv", inFn->result); break;
    case TECPLOT : sprintf(outputFilename, "%sInitialData.dat", inFn->result); break;
  }

  tInstant1 = clock(); // Start measuring time
  WriteResults(outputFilename, nodeX, nodeY, u, v, rho, nodeType, n, m, args->outputFormat);
  tInstant2 = clock();
  taskTime[T_WRIT] += (FLOAT_TYPE)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

  printf("\nInitialized data was written to %s\n", outputFilename);

  ////////////////// ITERATION ///////////////////////

  fprintf(logFile,"\n:::: Start Iterations ::::\n");
  printf("\n:::: Start Iterations ::::\n");

  printf("%d is the number of iterations \n", args->iterations);

  //for old mrt
  // FLOAT_TYPE *m_d    = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  // FLOAT_TYPE *mEq_d  = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  // FLOAT_TYPE *coll_d = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  // FLOAT_TYPE *f2_d   = createGpuArrayFlt(9*m*n, ARRAY_ZERO);
  // FLOAT_TYPE *fc2_d  = createGpuArrayFlt(9*m*n, ARRAY_ZERO);

  tIterStart = clock(); // Start measuring time of main loop

  int iter = 0;
  while (iter < args->iterations)
  {
    CHECK(cudaThreadSynchronize());
    CHECK(cudaEventRecord(start, 0)); // Start measuring time
    switch(args->collisionModel)
    {

      case BGKW:
        gpuCollBgkw<<<bpg1,tpb>>>(fluid_d, rho_d, u_d, v_d, f_d, fColl_d);
      break;

      case TRT:
        gpuCollTrt<<<bpg1,tpb>>>(fluid_d, rho_d, u_d, v_d, f_d, fColl_d);
      break;

      case MRT:
        // CHECK(cudaMemcpy(f2_d, f_d, SIZEFLT(9*m*n), cudaMemcpyDeviceToDevice));
        gpuCollMrt<<<bpg1,tpb>>>(fluid_d, rho_d, u_d, v_d, f_d, fColl_d);
        // gpu_mrt1<<<bpg9,tpb>>>(fluid_d, rho_d, u_d, v_d, f2_d, mEq_d, m_d);
        // gpu_mrt2<<<bpg9,tpb>>>(fluid_d, coll_d, m_d, mEq_d, fc2_d, f2_d);
        // printf("diff: %g   ", compareGpuArraySumFlt(fColl_d, fc2_d, temp9a_d, temp9b_d, 9*m*n));
      break;
    }


    CHECK(cudaEventRecord(stop, 0));
    CHECK(cudaEventSynchronize(stop));
    CHECK(cudaEventElapsedTime(&cudatime, start, stop));
    taskTime[T_COLL] += cudatime;

    ////////////// STREAMING ///////////////
    CHECK(cudaThreadSynchronize());
    CHECK(cudaEventRecord(start, 0));

    gpuStreaming<<<bpg1,tpb>>>(fluid_d, stream_d, f_d, fColl_d);

    CHECK(cudaEventRecord(stop, 0));
    CHECK(cudaEventSynchronize(stop));
    CHECK(cudaEventElapsedTime(&cudatime, start, stop));
    taskTime[T_STRM] += cudatime;

    // make the host block until the device is finished with foo
    CHECK(cudaThreadSynchronize());

    // check for error
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess)
    {
      // print the CUDA error message and exit
      printf("CUDA error: %s\n", cudaGetErrorString(error));
      exit(-1);
    }

    ////////////// BOUNDARIES ///////////////
    CHECK(cudaThreadSynchronize());
    CHECK(cudaEventRecord(start, 0));

    gpuBcInlet <<<bpgB,tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d, u0_d, v0_d, bcCount);
    gpuBcWall  <<<bpgB,tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d, fColl_d, qCollapsed_d, bcCount);
    gpuBcOutlet<<<bpgB,tpb>>>(bcIdxCollapsed_d, bcMaskCollapsed_d, f_d, u0_d, v0_d, bcCount);

    CHECK(cudaEventRecord(stop, 0));
    CHECK(cudaEventSynchronize(stop));
    CHECK(cudaEventElapsedTime(&cudatime, start, stop));
    taskTime[T_BNDC] += cudatime;


    // UPDATE VELOCITY AND DENSITY
    CHECK(cudaThreadSynchronize());
    CHECK(cudaEventRecord(start, 0));

    gpuUpdateMacro<<<bpg1,tpb>>>(fluid_d, rho_d, u_d, v_d, bcMask_d, drag_d, lift_d, coordX_d, coordY_d, f_d);

    tInstant2 = clock();
    CHECK(cudaEventRecord(stop, 0));
    CHECK(cudaEventSynchronize(stop));
    CHECK(cudaEventElapsedTime(&cudatime, start, stop));
    taskTime[T_MACR] += cudatime;

    // COMPUTE RESIDUALS
    CHECK(cudaThreadSynchronize());
    CHECK(cudaEventRecord(start, 0));

    /*FLOAT_TYPE r = computeResidual(f_d, fColl_d, temp9a_d, temp9b_d, m, n);
    if (r != r)
    {
      fprintf(stderr, "\nDIVERGENCE!\n");

      writeResiduals(residualsFilename, norm, dragSum, liftSum, m*n, iter+1);
      cudaEventDestroy(start);
      cudaEventDestroy(stop);

      freeAllHost(hostArrays, sizeof(hostArrays)/sizeof(hostArrays[0]));
      freeAllGpu (gpuArrays,  sizeof(gpuArrays) /sizeof(gpuArrays[0]) );

      return 1; // ERROR!
    }
    norm[iter] = r;
    if (args->boundaryId > 0)
    {
      dragSum[iter] = computeDragLift(bcMask_d, drag_d, tempA_d, tempB_d, m, n, args->boundaryId);
      liftSum[iter] = computeDragLift(bcMask_d, lift_d, tempA_d, tempB_d, m, n, args->boundaryId);
    }
    */
    CHECK(cudaEventRecord(stop, 0));
    CHECK(cudaEventSynchronize(stop));
    CHECK(cudaEventElapsedTime(&cudatime, start, stop));
    taskTime[T_RESI] += cudatime;

    printf("Iterating... %d/%d (%3.1f %%)\r", iter+1, args->iterations, (FLOAT_TYPE)(iter+1)*100/(FLOAT_TYPE)(args->iterations));

    iter++; // update loop variable

    ////////////// Autosave ///////////////

    if(iter == (args->autosaveEvery*autosaveIt))
    {
      autosaveIt++;
      if(iter > args->autosaveAfter)
      {
        printf("autosave\n\n");
        //////////// COPY VARIABLES TO HOST ////////////////
        CHECK(cudaMemcpy(u,   u_d,   SIZEFLT(m*n), cudaMemcpyDeviceToHost));
        CHECK(cudaMemcpy(v,   v_d,   SIZEFLT(m*n), cudaMemcpyDeviceToHost));
        CHECK(cudaMemcpy(rho, rho_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));

        switch(args->outputFormat)
        {
          case PARAVIEW: sprintf(autosaveFilename, "%sautosave_iter%05d.csv", inFn->result, iter); break;
          case TECPLOT : sprintf(autosaveFilename, "%sautosave_iter%05d.dat", inFn->result, iter); break;
        }

      tInstant1 = clock(); // Start measuring time
      WriteResults(autosaveFilename, nodeX, nodeY, u, v, rho, nodeType, n, m, args->outputFormat);
      tInstant2 = clock();
      taskTime[T_WRIT] += (FLOAT_TYPE)(tInstant2-tInstant1) / CLOCKS_PER_SEC;
      }
    }
  }     ////////////// END OF MAIN WHILE CYCLE! ///////////////

  tIterEnd = clock(); // End measuring time of main loop
  taskTime[T_ITER] = (FLOAT_TYPE)(tIterEnd - tIterStart ) / CLOCKS_PER_SEC;

  clock_t tEnd = clock();
  taskTime[T_OALL] = (FLOAT_TYPE)(tEnd - tStart) / CLOCKS_PER_SEC; // Calculate elapsed time
  taskTime[T_COLL] /= 1000;
  taskTime[T_STRM] /= 1000;
  taskTime[T_BNDC] /= 1000;
  taskTime[T_MACR] /= 1000;
  taskTime[T_RESI] /= 1000;

  fclose(logFile);
  writeEndLog(logFilename, taskTime);
  writeTimerLog(timeFilename, taskTime);
  //writeResiduals(residualsFilename, norm, dragSum, liftSum, m*n, args->iterations);

  // Write final data
  CHECK(cudaMemcpy(u,   u_d,   SIZEFLT(m*n), cudaMemcpyDeviceToHost));
  CHECK(cudaMemcpy(v,   v_d,   SIZEFLT(m*n), cudaMemcpyDeviceToHost));
  CHECK(cudaMemcpy(rho, rho_d, SIZEFLT(m*n), cudaMemcpyDeviceToHost));
  switch(args->outputFormat)
  {
    case PARAVIEW: sprintf(finalFilename, "%sFinalData.csv", inFn->result); break;
    case TECPLOT : sprintf(finalFilename, "%sFinalData.dat", inFn->result); break;
  }
  WriteResults(finalFilename, nodeX, nodeY, u, v, rho, nodeType, n, m, args->outputFormat);

  // Write information for user
  printf("\n\nLog was written to %s\n",               logFilename);
  printf("Last autosave result can be found at %s\n", autosaveFilename);
  //printf("residuals were written to %s\n",            residualsFilename);
  printf("Profiling results were written to %s\n",    timeFilename);
  printf("Final results were written to %s\n",        finalFilename);

  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  freeAllHost(hostArrays, sizeof(hostArrays)/sizeof(hostArrays[0]));
  freeAllGpu (gpuArrays,  sizeof(gpuArrays) /sizeof(gpuArrays[0]) );

  return 0;
}
