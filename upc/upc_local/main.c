// In a Linux System compile and run as (e.g. Ubuntu):
// upcc -T=4 main.c Iterate.c ShellFunctions.c FilesReading.c FilesWriting.c CellFunctions.c BoundaryConditions.c ComputeResiduals.c -lm -o LBMSolver && upcrun LBMSolver

#include <string.h>                      // String operations
#include <upc_relaxed.h>                 // Required for UPC 

#include "include/Iterate.h"             // Iteration takes place
#include "include/ShellFunctions.h"      // For convenience
#include "include/FilesReading.h"        // For reading files
#include "include/FilesWriting.h"        // For writing files e.g. tecplot
#include "include/CellFunctions.h"       // For cell modifications
#include "include/BoundaryConditions.h"  // boundary conditions

int main(int argc, char* argv[])
{

  /////////////////////////////////////////////////////////////////////////
  //////////////////////////// INITIAL DATA ///////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  // Name of folder to which we write the files
  char MainWorkDir[] = "Results";  

  if(MYTHREAD==0)
  {
     printf("\n###############################################\n");
       printf("#### This is a 2D lattice-Boltzmann solver ####\n");
       printf("#### Written in UPC, running on %02d threads ####\n", THREADS);
       printf("###############################################\n");
    // Create the working directory, continue if it already exist!
    CreateDirectory(MainWorkDir);
  }
  ///////////////////// Declare Simulation Variables //////////////////////
  
  // inlet parameters  
    float Uavg, Vavg, rho_ini, Viscosity;
  
  // integer (logical) inlet parameters
    int InletProfile, CollisionModel, CurvedBoundaries, OutletProfile, CalculateDragLift; 
  
  // numbers regarding the iterations
  int Iterations, AutosaveEvery, AutosaveAfter, PostprocProg;       
  float ConvergenceCritVeloc, ConvergenceCritRho;
  
  // Import Simulation Variables
  char IniFileName[] = "SetUpData.ini";
  ReadIniData(IniFileName,    &Uavg, &Vavg, &rho_ini, &Viscosity,
              &InletProfile,  &CollisionModel, &CurvedBoundaries,
              &OutletProfile, &Iterations, &AutosaveEvery,
              &AutosaveAfter, &PostprocProg, &CalculateDragLift,
              &ConvergenceCritVeloc,  &ConvergenceCritRho);
    
  // Mesh files
  char NodeDataFile[]="Mesh/D2node.dat";
  char BCconnectorDataFile[]="Mesh/BCconnectors.dat";
  
  /////////////////////////////////////////////////////////////////////////
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
            AutosaveAfter,         // autosave after this many iterations
            AutosaveEvery,         // autosave every #th iteration
            PostprocProg,          // postproc with Tecplot or ParaView
            CalculateDragLift,     // 0: no calculation, 1: calc on BC_ID (1), 2: calc on BC_ID (2), etc
            ConvergenceCritVeloc,  // convergence criterion for velocity 
            ConvergenceCritRho);   // convergence criterion for density

  ///////////////////////
  // END OF THE SOLVER //
  ///////////////////////

  return 0;
}

