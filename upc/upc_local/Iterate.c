#include <stdio.h>                      // printf();
#include <math.h>                       // need to compile with -lm
#include <stdlib.h>                     // for calloc();
#include <stdbool.h>                    // Include for bool type variables!
#include <string.h>                     // String operations
#include <time.h>                       // time functions
#include <upc_relaxed.h>                // Required for UPC 

////////////////////////////////////////////////////
////////////////// OWN HEADERS /////////////////////
////////////////////////////////////////////////////

#include "include/ShellFunctions.h"     // For convenience
#include "include/Iterate.h"            // Iteration takes place
#include "include/FilesReading.h"       // For reading files
#include "include/FilesWriting.h"       // For writing files e.g. tecplot
#include "include/CellFunctions.h"      // For cell modifications
#include "include/BoundaryConditions.h" // boundary conditions
#include "include/ComputeResiduals.h"   // Residuals

////////////////////////////////////////////////////
/////////////////// ITERATION //////////////////////
////////////////////////////////////////////////////

void Iteration(char* NodeDataFile, char* BCconnectorDataFile,
               float Uavg,         float Vavg,         float rho_ini, float Viscosity,
               int InletProfile,   int CollisionModel, int CurvedBoundaries,
               int OutletProfile,  int Iterations,     int AutosaveAfter,
               int AutosaveEvery,  int postproc_prog,  int CalculateDragLift,
               float ConvergenceCritVeloc, float ConvergenceCritRho)
{

  clock_t tStart;     // Time measurement: declaration, begin
  if(MYTHREAD==0)
    tStart = clock(); // BEGIN OF OVERALL TIME MEASUREMENT

  ////////////////////////////////////////////////////
  ///////////////////// Declare //////////////////////
  ////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////
  /////////// LOCAL VARIABLES FOR EACH CPU ///////////
  ////////////////////////////////////////////////////

  int i, j, k, iter = 0;                  // variables for loops
  FILE* resid_file;                       // file for residuals
  FILE* log_file;                         // file for log
  FILE* TimeMeasurementFile;              // file for time measurement results
  char IterOutputFile[50];                // write results to this file after the iterations
  char AutosaveOutputFile[50];            // autosave filename
  char OutputFile[50];                    // initial data will be written to this file
  char FinalOutputFile[50];               // final data will be written to this file
  char logFile[] = "Results/logFile.log"; // path of the .log file
  int  AutosaveI = 1;                     // autosave i variable, will be incremented after every autosave
  int* ppp;                               // pointer of the postproc_prog variable
  
  struct CellProps *Cells;                // Struct for Cells

  // Time measurement variables
  float tInitialization  = 0.0; // Time measurement of Initialization
  float tIteration       = 0.0; // Time measurement of Iteration
  float tCollision       = 0.0; // Time measurement of Collision
  float tUpdateF         = 0.0; // Time measurement of UpdateF
  float tStreaming       = 0.0; // Time measurement of Streaming
  float tBoundaries      = 0.0; // Time measurement of Boundaries
  float tUpdateMacro     = 0.0; // Time measurement of Update Macroscopic vars
  float tResiduals       = 0.0; // Time measurement of calculating residuals
  float tWriting         = 0.0; // Time measurement of writing data
  float tBCells          = 0.0; // Time measurement of handling boundaries
  clock_t tInstant1, tInstant2; // Time measurement points: universal
  clock_t tIterStart, tIterEnd; // Time measurement points: main loop
  
  // Variables for residuals
  double *Residuals;
  double *sumVel0;
  double *sumVel1;
  double *sumRho0;
  double *sumRho1; 
 
  float **Nodes;          // matrices for the nodes
  float **BCconn;         // matrices for connections

  double Omega;
  double OmegaA;          // collision frequency from the viscosity
  double **tm;            // variable for the MRT collision model
  double **stmiv;         // variable for the MRT collision model
  
  // D2Q9 Variables of the lattice
  double* w;              // weight values for the directions
  int*   cx;              // x coordinate of the discrete lattice directions
  int*   cy;              // y coordinate of the discrete lattice directions
  int*  opp;              // opposite vector
  int*    c;              // shift of lattice directions written in vector form
  
  ////////////////////////////////////////////////////
  //////////////////// ALLOCATE //////////////////////
  ////////////////////////////////////////////////////
  
  // allocate residuals
  sumVel0   = Create1DArrayDouble(1);
  sumVel1   = Create1DArrayDouble(1);
  sumRho0   = Create1DArrayDouble(1);
  sumRho1   = Create1DArrayDouble(1);
  Residuals = Create1DArrayDouble(4); 
  
  ////////////////////////////
  // THESE ARE SHARED STUFF //
  ////////////////////////////

  // allocate mesh properties  :: DEFINED IN ShellFunctions.h
  Delta          = (shared [] double*)upc_alloc(1*sizeof(double));
  m              = (shared []    int*)upc_alloc(1*sizeof(int));
  n              = (shared []    int*)upc_alloc(1*sizeof(int));
  NumNodes       = (shared []    int*)upc_alloc(1*sizeof(int));
  NumConn        = (shared []    int*)upc_alloc(1*sizeof(int));
  MaxInletCoordY = (shared [] double*)upc_alloc(1*sizeof(double));
  MinInletCoordY = (shared [] double*)upc_alloc(1*sizeof(double));
  NumInletNodes  = (shared []    int*)upc_alloc(1*sizeof(int));

  ///////////////////////////////////////////////////////////////////////
  //// WE HAVE BETTER TO STORE THIS ON EACH CPU, MAKE THINGS EASIER  ////
  ///////////////////////////////////////////////////////////////////////

  // D2Q9 Variables of the lattice
    w = Create1DArrayDouble(9); // weight values for the directions
    c = Create1DArrayInt(9);    // 
   cx = Create1DArrayInt(9);    // x coordinate of the discrete lattice directions
   cy = Create1DArrayInt(9);    // y coordinate of the discrete lattice directions
  opp = Create1DArrayInt(9);    // opposite vector

  ////////////////////////////////////////////////////
  ///////////////////// Read data ////////////////////
  ////////////////////////////////////////////////////
  
  Nodes   = ReadNodes(NodeDataFile);          // Read Node data
  BCconn  = ReadBCconn(BCconnectorDataFile);  // Read BCconn data
  CompDataNode(Nodes);
  CompDataConn(BCconn);

  ////////////////////////////////////////////////////
  /////////////// Print info to log //////////////////
  ////////////////////////////////////////////////////

  if(MYTHREAD==0) // Print information to log file
  {
    // Check whether we got back what we wanted :), write to log file!
    log_file = fopen(logFile, "w");  // open log file
    ppp      = &postproc_prog;       // for convenience ppp points to postproc_prog
    fprintf(log_file,"This is the 2D lattice Boltzmann *.log file\n\n");
    fprintf(log_file,"\n:::: Imported variables from the *.ini file :::: \n");
    fprintf(log_file,">>> Uavg              : %3.6f\n", Uavg);
    fprintf(log_file,">>> Vavg              : %3.6f\n", Vavg);
    fprintf(log_file,">>> Initial density   : %2.1f\n", rho_ini);
    fprintf(log_file,">>> Viscosity         : %3.8f\n", Viscosity);
    fprintf(log_file,">>> # of iterations   : %1.1d\n", Iterations);
    fprintf(log_file,">>> Autosave after    : %1.1d\n", AutosaveAfter);
    fprintf(log_file,">>> Autosave every    : %1.1d\n", AutosaveEvery);
    fprintf(log_file,">>> Convergence Veloc : %3.8f\n", ConvergenceCritVeloc);
    fprintf(log_file,">>> Convergence Rho   : %3.8f\n", ConvergenceCritRho);
    switch(CollisionModel)         // 1: BGKW, 2: TRT, 3: MRT
    {
      case 1: fprintf(log_file,">>> CollisionModel    : BGKW\n"); break;
      case 2: fprintf(log_file,">>> CollisionModel    : TRT\n" ); break;
      case 3: fprintf(log_file,">>> CollisionModel    : MRT\n" ); break;
    }
    switch(InletProfile)                      // 1:ON, 2:OFF
    {
      case 1: fprintf(log_file,">>> InletProfile      : ON\n" ); break;
      case 2: fprintf(log_file,">>> InletProfile      : OFF\n"); break;
    }
    switch(OutletProfile)                     // 1:ON, 2:OFF
    {
      case 1: fprintf(log_file,">>> OutletProfile     : ON\n" ); break;
      case 2: fprintf(log_file,">>> OutletProfile     : OFF\n"); break;
    }
    switch(CurvedBoundaries)                  // 1:ON, 2:OFF
    {
      case 1: fprintf(log_file,">>> CurvedBoundaries  : ON\n" ); break;
      case 2: fprintf(log_file,">>> CurvedBoundaries  : OFF\n"); break;
    }
    switch(postproc_prog)   // 1->Paraview (*.csv)     2->Tecplot 
    {
      case 1: fprintf(log_file,">>> Results format    : Paraview (*.csv)\n" ); break;
      case 2: fprintf(log_file,">>> Results format    : Tecplot (*.dat)\n"); break;
    }
    if (CalculateDragLift != 0)
              fprintf(log_file,">>> Drag, lift @ BC   : %d\n", CalculateDragLift);
    else 
              fprintf(log_file,">>> Drag, lift was not calculated\n");

    fprintf(log_file,"\n:::: Calculated variables from mesh :::: \n");
    fprintf(log_file,">>> Grid spacing        = %f\n", *Delta);
    fprintf(log_file,">>> # of nodes in x (n) = %d\n", *n);
    fprintf(log_file,">>> # of nodes in y (m) = %d\n", *m);
    fprintf(log_file,">>> NumInletNodes       = %d\n", *NumInletNodes);
    fprintf(log_file,">>> MaxInletCoordY      = %f\n", *MaxInletCoordY);
    fprintf(log_file,">>> MinInletCoordY      = %f\n", *MinInletCoordY);

    fprintf(log_file,"\n:::: Parallel properties :::: \n");
    fprintf(log_file,">>> # of threads        = %d\n", THREADS);
    fprintf(log_file,">>> BlockSize           = %d\n", BLOCKSIZE);

    // In case of no autosave
    sprintf(AutosaveOutputFile, "NOWHERE!");

  } // END OF THREAD ZERO

  ////////////////////////////////////////////////////
  ///////////////// INITIALIZE ///////////////////////
  ////////////////////////////////////////////////////

  // Fill up D2Q9 variables
  D2Q9Vars(w, cx, cy, opp, c);

  // Calculate collision frequency
  Omega  = 1.0/(3.*Viscosity+0.5);
  OmegaA = 8*(2-Omega)/(8-Omega);


  // Initialize variables for MRT Collision model, if used
  tm    = Create2DArrayDouble(9, 9);
  stmiv = Create2DArrayDouble(9, 9);

  if (CollisionModel == 3)
    MRTInitializer(tm, stmiv, Omega);


  //////////////////////////////////////////////////////
  // Allocate structure for the cell properties (see ShellFunctions.h)
  WCells = (shared [BLOCKSIZE]   struct CellProps*)upc_all_alloc(THREADS, BLOCKSIZE*sizeof(struct CellProps));
  BCells = (shared [2*NN]        struct CellProps*)upc_all_alloc(THREADS,    (*n)*2*sizeof(struct CellProps));
  //sResiduals =                 (shared [4] double*)upc_all_alloc(THREADS,                   4*sizeof(double));
  Cells = calloc(BLOCKSIZE+2*(*n),sizeof(struct CellProps));
  //////////////////////////////////////////////////////

  // PUT THREAD INFO TO BCELLS (boundary cells)
  upc_forall(i = 0; i < 2*(*n)*THREADS; i++; &BCells[i])
  { (BCells+i)->ThreadNumber = MYTHREAD; }

  if(MYTHREAD==0)
  {
    fprintf(log_file,"\n:::: Initializing ::::\n");
    printf("\n:::: Initializing ::::\n");
  } // END OF THREAD ZERO

  ////////////////////////////////////////////////////
  ///////////////// INITIALIZE CELLS /////////////////
  ////////////////////////////////////////////////////
  
  upc_barrier;         // Synchronise
  loop = 0;            // This will measure that how much is done of initialization
  tInstant1 = clock(); // Measure time of initialization

  for(j = 1;  j < MLIM+1;  j++)
  {
    for(i=0; i<*n; i++)
    {
      CellIni(Cells,
              Nodes,            // Nodes
							BCconn,           // BCconn
							i,                // index
							j,                // index
							Uavg,             // INPUT PARAMETER
							Vavg,             // INPUT PARAMETER
							InletProfile,     // INPUT PARAMETER
							CollisionModel,   // INPUT PARAMETER
							opp,              // Opposite direction
							rho_ini,          // Initial density
              -(2*MYTHREAD+2)+(MYTHREAD*(MLIM+2)+1) );  // shift cells index
      //printf("Cells[%d].ID = %d \n", j*(*n)+i, Cells[j*(*n)+i].ID );
      loop++;
    }
    if (MYTHREAD==0)
      printf("Initializing... %2.1f %%\r",(float)(loop)*100/(float)( (*n)*(*m) ));
  }
  tInstant2 = clock(); // Measure time of initialization 
  tInitialization = (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

  upc_barrier;         // Synchronise

  // PUT INITIALIZED DATA TO BOUNDARIES
  putCellsToShared(Cells);
  
  // COPY CELLS TO WCELLS (writing cells) TO WRITE DATA
  putCellsToWCells(Cells);
  
  // Write boundary cells to Results to see how mesh was distributed
  char  fnMemCopyRes[50]; 
  if(MYTHREAD==0)
  {
    switch(postproc_prog)
    { case 1: sprintf(fnMemCopyRes, "Results/MyBoundaryCells.csv");   break;
      case 2: sprintf(fnMemCopyRes, "Results/MyBoundaryCells.dat");   break; }

    WriteBCells(fnMemCopyRes, ppp);
  }

  // We dont need these matrices anymore
  free(Nodes);
  free(BCconn);
  upc_barrier;         // Synchronise
  
  ////////////////////////////////////////////////////
  /////////////// SAVE INITIALIZATION ////////////////
  ////////////////////////////////////////////////////

  if(MYTHREAD==0)
  {
    fprintf(log_file,"\n:::: Initialization done! ::::\n");
    printf("\n:::: Initialization done! ::::\n");

    // Write Initialized data 
    switch(postproc_prog)
      { case 1: sprintf(OutputFile, "Results/InitialData.csv");   break;
        case 2: sprintf(OutputFile, "Results/InitialData.dat");   break; }

    WriteResults(OutputFile, ppp);
    printf("\nInitialized data was written to %s\n", OutputFile);

    // Open residuals file
    resid_file = fopen("Results/residuals.dat", "w");
    fprintf(resid_file,"Iter Time Vel_res Rho_res Drag Lift\n");

    fprintf(log_file,"\n:::: Start Iterations ::::\n");
    printf("\n:::: Start Iterations ::::\n");
  } // END OF THREAD ZERO
  
  ////////////////////////////////////////////////////
  /////////////////// ITERATE ////////////////////////
  ////////////////////////////////////////////////////

  Residuals[0] = 100;
  Residuals[1] = 100;

  upc_barrier;         // Synchronise
  
  tIterStart = clock(); // Start measuring time of main loop  
  while (Residuals[0] > ConvergenceCritVeloc && Residuals[1] > ConvergenceCritRho && iter<Iterations)
  {

    //////////////// COLLISION ////////////////
		tInstant1 = clock();
    CollisionStep(Cells, w, cx, cy,  opp,  Omega,  OmegaA, tm, stmiv, CollisionModel); ////////////////////// !!!!!!!!!!!!!!!!! CX CY!
    tInstant2 = clock();
    tCollision = tCollision + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;
    
    ////////////// UPDATE DISTR ///////////////
		tInstant1 = clock(); // Start measuring time
    for(i = (*n);  i < (MLIM+1)*(*n);  i++)
      { UpdateF(Cells, i); }
    tInstant2 = clock();
    tUpdateF = tUpdateF + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

    // PUT THREAD-BOUNDARY CELLS TO SHARED
    tInstant1 = clock(); // Start measuring time
    putCellsToShared(Cells); 
    tInstant2 = clock(); // End of time measuremet
    tBCells = tBCells + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

    upc_barrier;         // Synchronise

    //////////////// COPY SHARED BCELLS TO CELLS ////////////////    
    tInstant1 = clock();
    getSharedToCells(Cells);
    tInstant2 = clock();
    tBCells = tBCells + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

    ////////////// STREAMING ///////////////
    tInstant1 = clock(); // Start measuring time
    StreamingStep(Cells, c);
    tInstant2 = clock();
    tStreaming = tStreaming + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;
    
    ////////////// BOUNDARIES ///////////////
    tInstant1 = clock(); // Start measuring time
    HandleBoundariesStep(Cells, cx, cy, c, opp, OutletProfile, CurvedBoundaries);
    tInstant2 = clock();
    tBoundaries = tBoundaries + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;
    
    // UPDATE VELOCITY AND DENSITY
    tInstant1 = clock(); // Start measuring time
    UpdateMacroscopicStep(Cells, cx, cy, CalculateDragLift);
    tInstant2 = clock();
    tUpdateMacro = tUpdateMacro + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

    ////////////// Residuals ///////////////
    tInstant1 = clock(); // Start measuring time
  	  ComputeResiduals(Cells, Residuals, sumVel0, sumVel1, sumRho0, sumRho1, CalculateDragLift, &iter, &Iterations);
    //ComputeResiduals(Cells, Residuals, sumVel0, sumVel1, sumRho0, sumRho1, CalculateDragLift, &iter, &Iterations);
    tInstant2 = clock(); // End of time measuremet
    tResiduals = tResiduals + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

    if(MYTHREAD==0)
      fprintf(resid_file,"%d %5.4e %5.4e %5.4e %f %f\n", iter, (iter+1.0)*(*Delta), Residuals[0], Residuals[1], Residuals[2], Residuals[3]);

    iter++; // update loop variable 

    if(iter%100==0 && MYTHREAD==0){
      printf("Iterations: %05d/%05d || ", iter, Iterations);
      printf("Residuals: l2norm  %e; L2_norm_weighted  %e\n", Residuals[0], Residuals[1]);
    }
    
    ////////////// Autosave ///////////////
    if(iter == (AutosaveEvery*AutosaveI))
    {
      AutosaveI++;
      if(iter>AutosaveAfter)
      {
        tInstant1 = clock(); // Start measuring time
        switch(postproc_prog) {
          case 1: sprintf(AutosaveOutputFile, "Results/autosave_iter%05d.csv", iter); break;
          case 2: sprintf(AutosaveOutputFile, "Results/autosave_iter%05d.dat", iter); break; }
        putCellsToWCells(Cells); // Put information to WCells and write (Write Cells)
        if (MYTHREAD==0) // AUTOSAVE
          WriteResults(AutosaveOutputFile, ppp);
        tInstant2 = clock(); // End of time measurement
        tWriting = tWriting + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;
      }
    }
        //////////////////////////////////////////////////////
  }     ////////////// END OF MAIN WHILE CYCLE ///////////////
        //////////////////////////////////////////////////////

  tIterEnd = clock(); // End measuring time of main loop
  tIteration = (float)(tIterEnd - tIterStart ) / CLOCKS_PER_SEC;

  upc_barrier;         // Synchronise

  if(MYTHREAD==0) // EXPORT DATA, TIME MEASUREMENT RESULTS
  {
    // Close residuals file
    fclose(resid_file);

    clock_t tEnd = clock();
    float tOverall = (float)(tEnd - tStart) / CLOCKS_PER_SEC; // Calculate elapsed time

    fprintf(log_file,"\nOverall calculations took %f seconds\n", tOverall);
    fprintf(log_file,"Main while loop took %f seconds\n",        tIteration);
    fprintf(log_file,"Initialization took %f seconds\n",         tInitialization);
    fprintf(log_file,"Collision took %f seconds\n",              tCollision);
    fprintf(log_file,"UpdateF took %f seconds\n",                tUpdateF);
    fprintf(log_file,"Streaming took %f seconds\n",              tStreaming);
    fprintf(log_file,"Calculating Boundaries took %f seconds\n", tBoundaries);
    fprintf(log_file,"Update Macroscopic took %f seconds\n",     tUpdateMacro);
    fprintf(log_file,"Calculating Residuals took %f seconds\n",  tResiduals);
    fprintf(log_file,"Writing results took %f seconds\n",        tWriting);
    fprintf(log_file,"Copying boundary cells took %f seconds\n", tBCells);

  
    // end time measurement, close log file
    fprintf(log_file,"\n:::: Iterations done! ::::\n");
    fclose(log_file);

    TimeMeasurementFile = fopen("Results/ParallelTimeMeasuerment.dat","w");
    fprintf(TimeMeasurementFile,"tOverall %f\n",        tOverall);
    fprintf(TimeMeasurementFile,"tIteration %f\n",      tIteration);
    fprintf(TimeMeasurementFile,"tInitialization %f\n", tInitialization);
    fprintf(TimeMeasurementFile,"tCollision %f\n",      tCollision);
    fprintf(TimeMeasurementFile,"tUpdateF %f\n",        tUpdateF);
    fprintf(TimeMeasurementFile,"tStreaming %f\n",      tStreaming);
    fprintf(TimeMeasurementFile,"tBoundaries %f\n",     tBoundaries);
    fprintf(TimeMeasurementFile,"tUpdateMacro %f\n",    tUpdateMacro);
    fprintf(TimeMeasurementFile,"tResiduals %f\n",      tResiduals);
    fprintf(TimeMeasurementFile,"tWriting %f\n",        tWriting);
    fprintf(TimeMeasurementFile,"tBCells %f\n",         tBCells);
    fprintf(TimeMeasurementFile,"THREADS %d\n",         THREADS);
     fclose(TimeMeasurementFile);
  
    // Write final data
    switch(postproc_prog){
      case 1: sprintf(FinalOutputFile, "Results/FinalData.csv"); break;
      case 2: sprintf(FinalOutputFile, "Results/FinalData.dat"); break;
    }
    WriteResults(FinalOutputFile,  ppp);

    // Write information for user
    printf("\n\nLog was written to %s\n", logFile);
    printf("Last autosave result can be found at %s\n", AutosaveOutputFile);
    printf("Residuals were written to Results/residuals.dat\n");
    printf("Profiling results were written to Results/ParallelTimeMeasuerment.dat\n");
    printf("Final results were written to %s\n", FinalOutputFile);

    WriteBCells(fnMemCopyRes, ppp);
    puts("BCells were written!");
  } // END OF THREAD ZERO

  
  ////////////////////////////////////////////////////
  ///////////////// End of line //////////////////////
  ////////////////////////////////////////////////////

  upc_barrier;         // Synchronise

  // FREE POINTERS

  upc_free(Delta);
  upc_free(m);
  upc_free(n);
  upc_free(MaxInletCoordY);
  upc_free(MinInletCoordY);
  upc_free(NumInletNodes);
  upc_free(NumNodes);
  upc_free(NumConn);
  
  free(Cells);
  free(w);
  free(cx);
  free(cy);
  free(c);
  free(opp);
  free(Residuals);

} // END OF ITERATIVE PROCESS

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////


/////////////// PUT DATA TO BOUNDARIES///////////////
void putCellsToShared(struct CellProps *Cells)
{
  upc_memput( &BCells[(*n)*(2*MYTHREAD)  ], &Cells[ (*n)       ], (*n)*sizeof(struct CellProps) ); // FIRST LINE
  upc_memput( &BCells[(*n)*(2*MYTHREAD+1)], &Cells[ (MLIM)*(*n)], (*n)*sizeof(struct CellProps) ); // LAST LINE
}

/////////////// GET DATA FROM BOUNDARIES ///////////////
void getSharedToCells(struct CellProps *Cells)
{
  
  if (MYTHREAD == THREADS-1)  // LAST THREAD ::: copy only the first line
  { 
    upc_memget( &Cells[   0            ], &BCells[(*n)*(2*MYTHREAD-1)], (*n)*sizeof(struct CellProps) ); // FIRST LINE

  } else if (MYTHREAD == 0){ // FIRST THREAD ::: copy only the last line

    upc_memget( &Cells[ (MLIM+1)*(*n)  ], &BCells[(*n)*(2*MYTHREAD+2)], (*n)*sizeof(struct CellProps) ); // LAST  LINE 

  } else {

    upc_memget( &Cells[   0            ], &BCells[(*n)*(2*MYTHREAD-1)], (*n)*sizeof(struct CellProps) ); // FIRST LINE
    upc_memget( &Cells[ (MLIM+1)*(*n)  ], &BCells[(*n)*(2*MYTHREAD+2)], (*n)*sizeof(struct CellProps) ); // LAST  LINE

  }

}

/////////////// PUT DATA TO WCELLS TO WRITE ///////////////
void putCellsToWCells(struct CellProps *Cells)
{
  upc_memput( &WCells[(1 -(2*MYTHREAD+2) +( MYTHREAD*(MLIM+2)+1 ))*(*n)], &Cells[1*(*n)], MLIM*(*n)*sizeof(struct CellProps) );
}


////////////////// COLLISION STEP ///////////////////
void CollisionStep(struct CellProps *Cells,
                   double* w, int* cx, int* cy, int* opp,
                   double Omega, double OmegaA, double **tm,
                   double **stmiv, int CollisionModel)
{

  int i, j;

  switch(CollisionModel)
  {
    // BGKW
      case 1:
        for(i = (*n);  i < (MLIM+1)*(*n);  i++)
        {
          if( (Cells+i)->Fluid == 1 )
            BGKW(Cells, i, w, cx, cy, Omega);
        }
      break;

      // TRT 
      case 2:
        for(i = (*n);  i < (MLIM+1)*(*n);  i++)
        {
          if( (Cells+i)->Fluid == 1)
              TRT (Cells, i, w, cx, cy, opp, Omega, OmegaA);
        }
      break;

      // MRT 
      case 3:
        for(i = (*n);  i < (MLIM+1)*(*n);  i++)
        {
          if( (Cells+i)->Fluid == 1)
            MRT(Cells, i, tm, stmiv); 
        }
      break;
  }

  /* OLD CODE, simpler but switches at every cell
    for (i=0; i<*m;i++)
    {
      for (j=0; j<*n;j++)
      {
        if (Cells[j][i].Fluid == 1)
        {
          switch(CollisionModel)
          {
            case 1:  BGKW(Cells, j, i, w, cx, cy, Omega);               break;
            case 2:  TRT (Cells, j, i, w, cx, cy, opp, Omega, OmegaA);  break;
            case 3:  MRT (Cells, j, i, tm, stmiv);                      break;
          }
        }
      }
    }
  */

} // END OF COLLISION STEP

////////////////// STREAMING STEP ///////////////////
void StreamingStep(struct CellProps *Cells, int* c)
{
  int i, j, k;
  //upc_forall(i=0; i<((*m)*(*n)); i++; &Cells[i])
  for(i = (*n);  i < (MLIM+1)*(*n);  i++)
  {
    if ( (Cells+i)->Fluid == 1 )
    {
      for(k=0; k<9; k++)
      {
        if ( ((Cells+i)->StreamLattice[k]) == 1)
        {
          (Cells +i)->F[k] = ( Cells+i+c[k] )-> METAF[k];

          // Cells[j][i].F[k] = Cells [j-cx[k]] [i-cy[k]].METAF[k]; // SERIAL CODE // 
          // C++ code
          // Cells[j][i].setF( Cells[j-cy[k]][i-cx[k]].METAF[k], k );
          // setF(newF, k){ F[k] = newF; };
        }
      }
    }
  }
} // END OF STREAMING STEP

///////////////// BOUNDARIES STEP ///////////////////
void HandleBoundariesStep(struct CellProps *Cells, int* cx, int* cy, int* c, int* opp, int OutletProfile, int CurvedBoundaries)
{
  int i, j, k;

  for(j = 1;  j < MLIM+1;  j++)
  {
    for(i=0; i<*n; i++)
    {
      if ((Cells +j*(*n)+i)->Fluid == 1)
      {
        // INLET
        InletBC(Cells, j, i);

        // WALL
        switch(CurvedBoundaries)
        {
          // curved boundaries
          case 1:
            for(k=0; k<9; k++)
              (Cells +j*(*n)+i)->Fneighbours[k] = ( Cells+j*(*n)+i +c[k] )-> METAF[k];

            //(Cells +j*(*n)+i)->Fneighbours[k] = (Cells + (j-cx[k])*(*m) + i-cy[k])->METAF[k]; // serial code, good
            // C++ code
            // Cells[j][i].setFneighbours(Cells[j-INI.getcy()[k]][i-1*INI.getcx()[k]].getMETAF()[k],k);
          
            CurvedWallBoundaries(Cells, j, i, opp);
          break;

          // bounceback boundaries
          case 2:
            WallBC(Cells, j, i, opp);
          break;
        }

        // OUTLET
        switch(OutletProfile)
        {
          // set profile in outlet
          case 1:
            OutletBoundaries(Cells, j, i);
          break;
          
          // OPEN BOUNDARY
          case 2 :
            if ((Cells +j*(*n)+i)->BC_ID[1]==3)
            { 
              
              (Cells +j*(*n)+i)->F[1] = 2*( (Cells+(*n)*(j)+i-1)->F[1] ) - (Cells+(*n)*(j)+i-2)->F[1]; 
              (Cells +j*(*n)+i)->F[5] = 2*( (Cells+(*n)*(j)+i-1)->F[5] ) - (Cells+(*n)*(j)+i-2)->F[5];
              (Cells +j*(*n)+i)->F[8] = 2*( (Cells+(*n)*(j)+i-1)->F[8] ) - (Cells+(*n)*(j)+i-2)->F[8];

              // C++ code
              //Cells[j][i].setF(2*Cells[j][i-1].getF()[1]-Cells[j][i-2].getF()[1],1);
              //Cells[j][i].setF(2*Cells[j][i-1].getF()[5]-Cells[j][i-2].getF()[5],5);
              //Cells[j][i].setF(2*Cells[j][i-1].getF()[8]-Cells[j][i-2].getF()[8],8);

            }
            /* !!!!  THIS IS NOT NECESSARY !!!! 
            if ((Cells +j*(*n)+i)->BC_ID[2]==3)
            {
              // FILL!!
            }
            
            if ((Cells +j*(*n)+i)->BC_ID[3]==3)
            {
              // FILL!!
            }
            
            if ((Cells +j*(*n)+i)->BC_ID[4]==3)
            {  
              // FILL!!
            }
            */
          break;
        }
      } 
    }
  }
}

/////////////// UPDATE MACRO STEP ///////////////////
void UpdateMacroscopicStep(struct CellProps *Cells, int* cx, int* cy, int CalculateDragLift)
{
  int i, j;
  //upc_forall(i=0; i<((*m)*(*n)); i++; &Cells[i])
  for(i = (*n);  i < (MLIM+1)*(*n);  i++)
  {
    if ((Cells+i)->Fluid == 1)
      UpdateMacroscopic(Cells, i, cx, cy, CalculateDragLift);
  }
}
