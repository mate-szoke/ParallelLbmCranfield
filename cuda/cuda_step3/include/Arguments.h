/**
 * This file contains the acceptable arguments for the solver.
 * @file Arguments.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */

#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include "FloatType.h"

/// Inlet profile options
typedef enum
{
    INLET=1,        ///< inlet profile
    NO_INLET,       ///< no inlet profile
    PULSATILE_INLET ///< pulsatile inlet profile @warning not implemented
} InletProfile;

/// Outlet profile options
typedef enum
{
    OUTLET=1,      ///< outlet profile
    OUTLET_SECOND, ///< open boundary second order
    OUTLET_FIRST   ///< open boundary first order
} OutletProfile;

/// Collision models
typedef enum
{
    BGKW=1, ///< BGKW method
    TRT,    ///< TRT method
    MRT     ///< MRT method
} CollisionModel;

/// Boundary types
typedef enum
{
    CURVED=1, ///< curved boundaries
    STRAIGHT  ///< straight boundaries
} BoundaryType;

/// Output formats
typedef enum
{
    PARAVIEW=1, ///< paraview format (.csv)
    TECPLOT     ///< tecplot format (.dat)
} OutputFormat;

/// Result from command line argument parsing
typedef enum
{
    NORMAL=0, ///< normal execution, parameters read
    HELP,     ///< print help
    INIT,     ///< read values from init file
    TEST,     ///< run unit tests
    COMPARE,  ///< compare results
    ERROR     ///< error happened
} ArgResult;

/// Input parameters
typedef struct arg
{
    FLOAT_TYPE u;                  ///< velocity x component
    FLOAT_TYPE v;                  ///< velocity y component
    FLOAT_TYPE rho;                ///< density
    FLOAT_TYPE viscosity;          ///< viscosity
    InletProfile inletProfile;     ///< inlet profile
    OutletProfile outletProfile;   ///< outlet profile
    CollisionModel collisionModel; ///< collision model
    BoundaryType boundaryType;     ///< boundary type
    OutputFormat outputFormat;     ///< output format
    int iterations;                ///< number of iterations
    int autosaveEvery;             ///< autosave every n iteration
    int autosaveAfter;             ///< autosave after nth iteration
    int boundaryId;                ///< boundary ID for drag/lift calculation
} Arguments;

/// Input filenames
typedef struct inputf
{
    char init[512];   ///< init file name
    char node[512];   ///< node file name
    char bc[512];     ///< bc file name
    char result[512]; ///< results directory name
    char comp[512];   ///< filename to compare
    char final[512];  ///< filename to compare to
} InputFilenames;

/**
 * Parse command line arguments
 * @param argc number of command line arguments
 * @param argv command line arguments
 * @param inFn input filenames
 * @param args input parameters
 * @return task that need to be done, or error
 *
 * @test MANUAL:
 * In order to test argument handling try some of these cases
 *   - run with -h, should display help
 *   - run without arguments, should work as written in SetUpData.ini
 *   - run with -f \<ini\> and any other argument, should work as written in SetUpData.ini
 *   - run with -o \<folder\> should output results to given folder
 *   - run with -i \<n\> should do iterations n times
 *   - run with -n \<file\> with non-existent file, should print out error message
 *   - try using long argument options as well, should work as the short version
 *   - run in compare mode with one file given, should compare Results/FinalData.csv to given file
 *   - run in compare mode with two files, should compare the 2 given files
 */
ArgResult handleArguments(int argc, char **argv, InputFilenames *inFn, Arguments *args);

#endif
