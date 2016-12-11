// In a Linux System compile and run as (e.g. Ubuntu):
// gcc main.c Iterate.c ShellFunctions.c FilesReading.c FilesWriting.c CellFunctions.c BoundaryConditions.c ComputeResiduals.c -lm -o LBMSolver && ./LBMSolver
// Intel C compiler:
// icc *.c -lm -o LBMSolver
// GNU C compiler:
// gcc *.c -lm -o LBMSolver

#include <string.h>                      // String operations
#include <stdio.h>                       // Standard input output commands (e.g. printf)
#include "include/Iterate.h"             // Iteration takes place
#include "include/ShellFunctions.h"      // For convenience
#include "include/FilesIO.h"             // For reading and writing file
//#include "include/FilesReading.h"        // For reading files
//#include "include/FilesWriting.h"        // For writing files e.g. tecplot
#include "include/CellFunctions.h"       // For cell modifications
#include "include/BoundaryConditions.h"  // boundary conditions



int main(int argc, char* argv[])
{

   printf("\n###############################################\n");
     printf("#### This is a 2D lattice-Boltzmann solver ####\n");
     printf("#### Written in C, running on single core  ####\n");
     printf("###############################################\n");

  /////////////////////////////////////////////////////////////////////////
  //////////////////////////// INITIAL DATA ///////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  // Name of folder to which we write the files
  char MainWorkDir[] = "Results";  

  // Create the working directory, continue if it already exist!
  CreateDirectory(MainWorkDir);
  
  ///////////////////// Declare Simulation Variables //////////////////////
  
  // inlet parameters  
  float Uavg, Vavg, rho_ini, Viscosity, ConvergenceCritVeloc, ConvergenceCritRho;
  // integer (logical) inlet parameters
  
  int InletProfile, CollisionModel, CurvedBoundaries, OutletProfile,
      CalculateDragLift, ReadFormerData; 
  // numbers regarding the iterations
  int Iterations, AutosaveEvery, AutosaveAfter, PostprocProg;       

  // Import Simulation Variables
  char IniFileName[] = "SetUpData.ini";
  ReadIniData(IniFileName,    &Uavg, &Vavg, &rho_ini, &Viscosity,
              &InletProfile,  &CollisionModel, &CurvedBoundaries,
              &OutletProfile, &Iterations, &AutosaveEvery,
              &AutosaveAfter, &PostprocProg, &CalculateDragLift,
              &ConvergenceCritVeloc, &ConvergenceCritRho, &ReadFormerData);

  // Mesh files
  char NodeDataFile[]="Mesh/D2node.dat";
  char BCconnectorDataFile[]="Mesh/BCconnectors.dat";

  /////////                                                       /////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  //////////////////////
  // START THE SOLVER //
  //////////////////////
 
  Iteration(NodeDataFile,          // data from mesher            
            BCconnectorDataFile,   // data from mesher            
            Uavg,                  // mean x velocity 
            Vavg,                  // mean y velocity 
            rho_ini,               // initial density
            Viscosity,             // viscosity of fluid 
            InletProfile,          // do we have inlet profile? 
            CollisionModel,        // which collision model to use
            CurvedBoundaries,      // do we have curved boundaries?
            OutletProfile,         // do we have outlet profile?
            Iterations,            // how many iterations to perform
            AutosaveAfter,         // autosave after this many iteration
            AutosaveEvery,         // autosave every #th iteration
            PostprocProg,          // postproc with Tecplot or ParaView
            CalculateDragLift,     // 0: no calculation, 1: calc on BC_ID (1), 2: calc on BC_ID (2), etc
            ConvergenceCritVeloc,  // Convergence criterion for velocity
            ConvergenceCritRho,    // Convergence criterion for density
            ReadFormerData);       // Do we want to read former data?

  ///////////////////////
  // END OF THE SOLVER //
  ///////////////////////

  return 0;
}

