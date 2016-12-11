/**
 * Main file
 * @file main.cu
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#include <string.h>              // String operations
#include <stdio.h>

#include "Iterate.h"             // Iteration takes place
#include "ShellFunctions.h"      // For convenience
#include "FilesReading.h"        // For reading files
#include "FilesWriting.h"        // For writing files e.g. tecplot
#include "CellFunctions.h"       // For cell modifications
#include "ComputeResiduals.h"
#include "Arguments.h"
#include "AllTests.h"

/**
 * Run all unittests
 * @return 0: success, n: number of failures
 */
int runAllTests(void)
{
  CuString *output = CuStringNew();
  CuSuite* suite = CuSuiteNew();

  CuSuiteAddSuite(suite, computeResidualsGetSuite());
  CuSuiteAddSuite(suite, gpuSumGetSuite());
  CuSuiteAddSuite(suite, gpuBoundariesGetSuite());
  CuSuiteAddSuite(suite, gpuInitGetSuite());
  CuSuiteAddSuite(suite, gpuMacroGetSuite());
  CuSuiteAddSuite(suite, gpuCollisionGetSuite());
  CuSuiteAddSuite(suite, gpuStreamGetSuite());
  CuSuiteAddSuite(suite, iterateGetSuite());

  CuSuiteRun(suite);
  CuSuiteSummary(suite, output);
  CuSuiteDetails(suite, output);
  printf("\nSUMMARY\n%s\n", output->buffer);

  cudaDeviceReset();

  return suite->failCount;
}

/**
 * @brief Main function
 * @details
 *  - handles input parameters
 *  - creates directory for results
 *  - runs the solver
 *
 * @param argc number of command line parameters
 * @param argv command line parameters
 *
 * @return 0: success n: error
 */
int main(int argc, char* argv[])
{
  Arguments args;
  args.u         = 0.;
  args.v         = 0.;
  args.rho       = 0.;
  args.viscosity = 0.;
  args.inletProfile   = NO_INLET;
  args.outletProfile  = OUTLET_SECOND;
  args.collisionModel = BGKW;
  args.boundaryType   = STRAIGHT;
  args.outputFormat   = PARAVIEW;
  args.iterations    = 5000;
  args.autosaveEvery = 0;
  args.autosaveAfter = 0;
  args.boundaryId    = 0;

  InputFilenames inFn;
  strcpy(inFn.node  , "Mesh/D2node.dat");
  strcpy(inFn.bc    , "Mesh/BCconnectors.dat");
  strcpy(inFn.init  , "SetUpData.ini");
  strcpy(inFn.result, "Results/");
  strcpy(inFn.final , "Results/FinalData.csv");

  if (argc > 1)
  {
    switch (handleArguments(argc, argv, &inFn, &args))
    {
      case INIT:
        readInitFile(inFn.init, &args);
      break;
      case HELP:
        return 0;
      case TEST:
        return runAllTests();
      case COMPARE:
        return compareFiles(inFn.comp, inFn.final);
      case ERROR:
        return 1;
    }
  }
  else
  {
    readInitFile(inFn.init, &args);
  }

  CreateDirectory(inFn.result);

  Iteration(&inFn, &args);

  cudaDeviceReset(); //needed for cuda-gdb and cuda-memcheck

  return 0;
}

