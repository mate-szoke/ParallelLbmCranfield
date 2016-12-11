/**
 * All unittests need to be declared in this header
 * @file AllTests.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#ifndef TEST_HEADERS
#define TEST_HEADERS

#include "CuTest.h"

/// Unittests for residual computations
CuSuite* computeResidualsGetSuite();
/// Unittests for sum on GPU
CuSuite* gpuSumGetSuite();
/// Unittests for boundary conditions
CuSuite* gpuBoundariesGetSuite();
/// Unittests for initialisation
CuSuite* gpuInitGetSuite();
/// Unittests for macroscopic values
CuSuite* gpuMacroGetSuite();
/// Unittests for collision models
CuSuite* gpuCollisionGetSuite();
/// Unittests for streaming
CuSuite* gpuStreamGetSuite();
/// Unittests for the solver
CuSuite* iterateGetSuite();

#endif
