#include <stdio.h>                      // printf();
#include <math.h>                       // need to compile with -lm
#include <stdlib.h>                     // for calloc();
#include <stdbool.h>                    // Include for bool type variables!
#include <string.h>                     // String operations
#include <time.h>                       // time functions
#include <upc.h>                        // Required for UPC 
//#include <upc_cray.h>                   // Required for UPC 
//#include <mpi.h>

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

double MPI_Wtime(void){
	double mytime=0.0;
	mytime = clock();
  mytime = (double)(mytime/CLOCKS_PER_SEC);
	return mytime;
}


////////////////////////////////////////////////////
/////////////////// ITERATION //////////////////////
////////////////////////////////////////////////////

void Iteration(char* NodeDataFile, char* BCconnectorDataFile,
               float Uavg,         float Vavg,         float rho_ini, float Viscosity,
               int InletProfile,   int CollisionModel, int CurvedBoundaries,
               int OutletProfile,  int Iterations,     int AutosaveAfter,
               int AutosaveEvery,  int postproc_prog,  int CalculateDragLift)
{
  //clock_t tStart;     // Time measurement: declaration, begin
  double tStart;     // Time measurement: declaration, begin
  if(MYTHREAD==0)
    tStart = MPI_Wtime(); // BEGIN OF OVERALL TIME MEASUREMENT

  ////////////////////////////////////////////////////
  ///////////////////// Declare //////////////////////
  ////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////
  /////////// LOCAL VARIABLES FOR EACH CPU ///////////
  ////////////////////////////////////////////////////

  int i, j, k, iter = 0;                  //variables for loops
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
  int lm, ln;                             // local m and n variables

  // Time measurement variables
  double tInitialization = 0.0;  // Time measurement of Initialization
  double tIteration      = 0.0;  // Time measurement of Iteration
  double tCollision      = 0.0;  // Time measurement of Collision
  double tUpdateF        = 0.0;  // Time measurement of UpdateF
  double tStreaming      = 0.0;  // Time measurement of Streaming
  double tBoundaries     = 0.0;  // Time measurement of Boundaries
  double tUpdateMacro    = 0.0;  // Time measurement of Update Macroscopic vars
  double tResiduals      = 0.0;  // Time measurement of calculating residuals
  double tWriting        = 0.0;  // Time measurement of writing data
  double tInstant1, tInstant2; // Time measurement points, universal
  double tIterStart, tIterEnd; // Time measurement points: main loop
  //clock_t tInstant1, tInstant2; // Time measurement points, universal
  //clock_t tIterStart, tIterEnd; // Time measurement points: main loop
  
  // Variables for residuals
  double *Residuals;
  double *sumVel0;
  double *sumVel1;
  double *sumRho0;
  double *sumRho1; 

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////// VARIABLES FOR ALL THREADS :::: SHARED ::: ////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  
  float **Nodes;          // matrices for the nodes
  float **BCconn;         // matrices for connections

  double Omega, OmegaA;   // collision frequency from the viscosity
  float **tm;             // variable for the MRT collision model
  float **stmiv;          // variable for the MRT collision model

  // D2Q9 Variables of the lattice
  double* w;    // weight values for the directions
  int*   cx;    // x coordinate of the discrete lattice directions
  int*   cy;    // y coordinate of the discrete lattice directions
  int*  opp;    // opposite vector
  int*    c;    // shift of lattice directions written in vector form
  
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
  Delta          = (shared [] float*)upc_alloc(1*sizeof(float));
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
  /////////////// Compute constants //////////////////
  ////////////////////////////////////////////////////

  
  if(MYTHREAD==0) // Print information to log file
  {
    // Check whether we got back what we wanted :), write to log file!
    log_file = fopen(logFile, "w");  // open log file
    ppp      = &postproc_prog;       // for convenience ppp points to postproc_prog
    fprintf(log_file,"This is the 2D lattice Boltzmann *.log file\n\n");
    fprintf(log_file,"\n:::: Imported variables from the *.ini file :::: \n");
    fprintf(log_file,">>> Uavg             : %3.6f\n", Uavg);
    fprintf(log_file,">>> Vavg             : %3.6f\n", Vavg);
    fprintf(log_file,">>> Initial density  : %2.1f\n", rho_ini);
    fprintf(log_file,">>> Viscosity        : %3.8f\n", Viscosity);
    fprintf(log_file,">>> # of iterations  : %1.1d\n", Iterations);
    fprintf(log_file,">>> Autosave after   : %1.1d\n", AutosaveAfter);
    fprintf(log_file,">>> Autosave every   : %1.1d\n", AutosaveEvery);
    switch(CollisionModel)         // 1: BGKW, 2: TRT, 3: MRT
    {
      case 1: fprintf(log_file,">>> CollisionModel   : BGKW\n"); break;
      case 2: fprintf(log_file,">>> CollisionModel   : TRT\n" ); break;
      case 3: fprintf(log_file,">>> CollisionModel   : MRT\n" ); break;
    }
    switch(InletProfile)                      // 1:ON, 2:OFF
    {
      case 1: fprintf(log_file,">>> InletProfile     : ON\n" ); break;
      case 2: fprintf(log_file,">>> InletProfile     : OFF\n"); break;
    }
    switch(OutletProfile)                     // 1:ON, 2:OFF
    {
      case 1: fprintf(log_file,">>> OutletProfile    : ON\n" ); break;
      case 2: fprintf(log_file,">>> OutletProfile    : OFF\n"); break;
    }
    switch(CurvedBoundaries)                  // 1:ON, 2:OFF
    {
      case 1: fprintf(log_file,">>> CurvedBoundaries : ON\n" ); break;
      case 2: fprintf(log_file,">>> CurvedBoundaries : OFF\n"); break;
    }
    switch(postproc_prog)   // 1->Paraview (*.csv)     2->Tecplot 
    {
      case 1: fprintf(log_file,">>> Results format   : Paraview (*.csv)\n" ); break;
      case 2: fprintf(log_file,">>> Results format   : Tecplot (*.dat)\n"); break;
    }
    if (CalculateDragLift != 0)
              fprintf(log_file,">>> Drag, lift @ BC  : %d\n", CalculateDragLift);
    else 
              fprintf(log_file,">>> Drag, lift was not calculated\n");

    fprintf(log_file,"\n:::: Calculated variables from mesh :::: \n");
    fprintf(log_file,">>> Grid spacing        = %f\n", *Delta);
    fprintf(log_file,">>> # of nodes in x (n) = %d\n", *n);
    fprintf(log_file,">>> # of nodes in y (m) = %d\n", *m);
    fprintf(log_file,">>> NumInletNodes       = %d\n", *NumInletNodes);
    fprintf(log_file,">>> MaxInletCoordY      = %f\n", *MaxInletCoordY);
    fprintf(log_file,">>> MinInletCoordY      = %f\n", *MinInletCoordY);
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

  // Calculate collision freq
  Omega  = 1.0/(3.*Viscosity+0.5);
  OmegaA = 8*(2-Omega)/(8-Omega);

  // Initialize variables for MRT Collision model, if used
  if (CollisionModel == 3)
  {
    tm    = Create2DArrayFloat(9, 9);
    stmiv = Create2DArrayFloat(9, 9);
    MRTInitializer(tm, stmiv, Omega);
  }  

  //////////////////////////////////////////////////////
  // Allocate structure for the cell properties >>> see ShellFunctions.h
  // OLD ALLOCATION TYPE:
  // Cells = (shared [BLOCKSIZE] struct CellProps*)upc_all_alloc(THREADS, BLOCKSIZE*sizeof(struct CellProps));
  //////////////////////////////////////////////////////

  if(MYTHREAD==0)
  {
    fprintf(log_file,"\n:::: Initializing ::::\n");
    printf("\n:::: Initializing ::::\n");
  } // END OF THREAD ZERO

  ////////////////////////////////////////////////////
  ///////////////// INITIALIZE CELLS /////////////////
  ////////////////////////////////////////////////////

  // LocalN and LocalM variables
  lm = *m;
  ln = *n;
  //printf("ln = %d, lm = %d\n", ln, lm);

  // Find affinity of variables.
  int MyAffinity[MM*NN];
  for(j=0; j<lm; j++)
  {
    for(i=0; i<ln; i++)
    {
      MyAffinity[j*(ln)+i] = upc_threadof(&Cells[j*(ln)+i]);
    }
  }

  upc_barrier;             // Synchronise
  loop = 0;                // This will measure that how much is done of initialization
  tInstant1 = MPI_Wtime(); // Measure time of initialization
  for(j=0; j<lm; j++)
  {
    upc_forall(i=0; i<ln; i++; MyAffinity[j*(ln)+i])
    {
      CellIni(Nodes,            // Nodes
							BCconn,           // BCconn
							i,                // index
							j,                // index
							Uavg,             // INPUT PARAMETER
							Vavg,             // INPUT PARAMETER
							InletProfile,     // INPUT PARAMETER
							CollisionModel,   // INPUT PARAMETER
							opp,              // Opposite direction
							rho_ini,          // Initial density
              lm, ln);          // mesh size
      loop++;
    } //} 
    if (MYTHREAD==0)
     printf("Initializing... %2.1f %%\r",(float)(loop)*100/(float)((ln)*(lm)));
  }
  upc_barrier;
  tInstant2 = MPI_Wtime(); // Measure time of initialization
  tInitialization = (tInstant2-tInstant1) ;

  // We dont need these matrices anymore
  free(Nodes);
  free(BCconn);


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

    WriteResults(OutputFile, ppp, lm, ln);
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

  upc_barrier;
  
  tIterStart = MPI_Wtime(); // Start measuring time of main loop  
  while (iter<Iterations)
  {
    //upc_barrier;
    //////////////// COLLISION ////////////////
		tInstant1 = MPI_Wtime();
    CollisionStep(w, cx, cy,  opp,  Omega,  OmegaA, tm, stmiv, CollisionModel, lm, ln, MyAffinity);
    tInstant2 = MPI_Wtime();
    tCollision = tCollision + (float)(tInstant2-tInstant1) ;
    
    ////////////// UPDATE DISTR ///////////////
		tInstant1 = MPI_Wtime(); // Start measuring time
    upc_forall(i=0; i<((lm)*(ln)); i++; MyAffinity[i])
    {
       for(k=0;k<9;k++)
        F[k][i] = METAF[k][i];
    } 
    tInstant2 = MPI_Wtime();
    tUpdateF = tUpdateF + (float)(tInstant2-tInstant1) ;

    //////////// BARRIER //////////////////////
    tInstant1 = MPI_Wtime(); 
    upc_barrier;
    tInstant2 = MPI_Wtime();
    tResiduals = tResiduals + (float)(tInstant2-tInstant1) ;
    
    ////////////// STREAMING ///////////////
    tInstant1 = MPI_Wtime(); // Start measuring time
    StreamingStep(c, lm, ln, MyAffinity);
    tInstant2 = MPI_Wtime();
    tStreaming = tStreaming + (float)(tInstant2-tInstant1) ;
    
    ////////////// BOUNDARIES ///////////////
    tInstant1 = MPI_Wtime(); // Start measuring time
    HandleBoundariesStep(cx, cy, opp, OutletProfile, CurvedBoundaries, lm, ln, MyAffinity);
    tInstant2 = MPI_Wtime();
    tBoundaries = tBoundaries + (float)(tInstant2-tInstant1) ;
    
    // UPDATE VELOCITY AND DENSITY
    tInstant1 = MPI_Wtime(); // Start measuring time
    UpdateMacroscopicStep(cx, cy, CalculateDragLift, lm, ln, MyAffinity);
    tInstant2 = MPI_Wtime();
    tUpdateMacro = tUpdateMacro + (float)(tInstant2-tInstant1) ;
    
    /* RESIDUAL CALCULATION IS SWITCHED OFF
    if(MYTHREAD==0)
    {
      ////////////// Residuals ///////////////
      tInstant1 = MPI_Wtime(); // Start measuring time
      ComputeResiduals(Residuals, sumVel0, sumVel1, sumRho0, sumRho1, CalculateDragLift, lm, ln);
      fprintf(resid_file,"%d %5.4e %5.4e %5.4e %f %f\n", iter, (iter+1.0)*(*Delta), Residuals[0], Residuals[1], Residuals[2], Residuals[3]);
      tInstant2 = MPI_Wtime();
      tResiduals = tResiduals + (float)(tInstant2-tInstant1) ;
      //printf("Iterating... %d/%d (%3.1f %%)\r", iter+1, Iterations, (float)(iter+1)*100/(float)(Iterations));
    }
    */
    iter++; // update loop variable

    if(iter%100==0 && MYTHREAD==0){
      printf("Iterations: %06d/%06d || ", iter, Iterations);
      printf("Residuals: Velocity  %e; Density  %e\n", Residuals[0], Residuals[1]);
    }
    
    if (MYTHREAD==0) // AUTOSAVE
    { 
      ////////////// Autosave ///////////////
      if(iter == (AutosaveEvery*AutosaveI))
      {
        AutosaveI++;
        if(iter>AutosaveAfter)
        {
          tInstant1 = MPI_Wtime(); // Start measuring time
          switch(postproc_prog) {
            case 1: sprintf(AutosaveOutputFile, "Results/autosave_iter%05d.csv", iter); break;
            case 2: sprintf(AutosaveOutputFile, "Results/autosave_iter%05d.dat", iter); break; }
          WriteResults(AutosaveOutputFile, ppp, lm, ln);
          tInstant2 = MPI_Wtime(); // End of time measurement
          tWriting = tWriting + (float)(tInstant2-tInstant1) ;
        }
      }
    } // END OF THREAD ZERO

  }     ////////////// END OF MAIN WHILE CYCLE ///////////////

  tIterEnd = MPI_Wtime(); // End measuring time of main loop
  tIteration = (float)(tIterEnd - tIterStart ) ;


  if(MYTHREAD==0) // EXPORT DATA, TIME MEASUREMENT RESULTS
  {
    // Close residuals file
    fclose(resid_file);

    double tEnd = MPI_Wtime();
    float tOverall = (float)(tEnd - tStart) ; // Calculate elapsed time

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
    fprintf(TimeMeasurementFile,"THREADS %d\n",         THREADS);
     fclose(TimeMeasurementFile);

    // Write final data
    switch(postproc_prog){
      case 1: sprintf(FinalOutputFile, "Results/FinalData.csv"); break;
      case 2: sprintf(FinalOutputFile, "Results/FinalData.dat"); break;
    }
    WriteResults(FinalOutputFile,  ppp, lm, ln);

    // Write information for user
    printf("\n\nLog was written to %s\n", logFile);
    printf("Last autosave result can be found at %s\n", AutosaveOutputFile);
    printf("Residuals were written to Results/residuals.dat\n");
    printf("Profiling results were written to Results/ParallelTimeMeasuerment.dat\n");
    printf("Final results were written to %s\n", FinalOutputFile);


  } // END OF THREAD ZERO

  ////////////////////////////////////////////////////
  ///////////////// End of line //////////////////////
  ////////////////////////////////////////////////////

  // FREE POINTERS
  upc_free(Delta);
  upc_free(m);
  upc_free(n);
  upc_free(MaxInletCoordY);
  upc_free(MinInletCoordY);
  upc_free(NumInletNodes);
  upc_free(NumNodes);
  upc_free(NumConn);
  //upc_free(Cells);
  
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




////////////////// COLLISION STEP ///////////////////
void CollisionStep(double* w, int* cx, int* cy, int* opp,
                   double Omega, double OmegaA, float **tm,
                   float **stmiv, int CollisionModel, int lm, int ln,
                   int MyAffinity[NN*MM])
{

  int i, j;

  switch(CollisionModel)
  {
    // BGKW
      case 1:
        upc_forall(i=0; i<((lm)*(ln)); i++; MyAffinity[i])
        //for(i=0; i<(ln)*(lm); i++) { if ( upc_threadof( &Cells[i] ) == MYTHREAD )
        {
          if( Fluid[i] == 1 )
            BGKW(i, w, cx, cy, Omega);
        } //}
      break;

      // TRT 
      case 2:
        upc_forall(i=0; i<((lm)*(ln)); i++; MyAffinity[i])
        //for(i=0; i<(ln)*(lm); i++) { if ( upc_threadof( &Cells[i] ) == MYTHREAD )
        {
          if( Fluid[i] == 1)
              TRT (i, w, cx, cy, opp, Omega, OmegaA);
        } //}
      break;

      // MRT 
      case 3:
        upc_forall(i=0; i<((lm)*(ln)); i++; MyAffinity[i])
        //for(i=0; i<(ln)*(lm); i++) { if ( upc_threadof( &Cells[i] ) == MYTHREAD )
        {
          if( Fluid[i] == 1)
            MRT(i, tm, stmiv); 
        } //}
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
void StreamingStep(int* c, int lm, int ln, int MyAffinity[NN*MM])
{
  int i, j, k;
  upc_forall(i=0; i<((lm)*(ln)); i++; MyAffinity[i])
  //for(i=0; i<(*n)*(*m); i++) { if ( upc_threadof( &Cells[i] ) == MYTHREAD )
  {
    if ( Fluid[i] == 1 )
    {
      for(k=0; k<9; k++)
      {
        if ( (StreamLattice[k][i]) == 1)
        {
          F[k][i] = METAF[k] [i + c[k]];

          // Cells[j][i].F[k] = Cells [j-cx[k]] [i-cy[k]].METAF[k]; // SERIAL CODE // 
          // C++ code
          // Cells[j][i].setF( Cells[j-cy[k]][i-cx[k]].METAF[k], k );
          // setF(newF, k){ F[k] = newF; };
        }
      }
    }
  } //}
}  // END OF STREAMING STEP

///////////////// BOUNDARIES STEP ///////////////////
void HandleBoundariesStep(int* cx, int* cy, int* opp, int OutletProfile, int CurvedBoundaries, int lm, int ln, int MyAffinity[NN*MM])
{
  int i, j, k;
  for(j=0; j<lm; j++)
  {
    upc_forall(i=0; i<ln; i++; MyAffinity[j*(ln)+i])
    //for(i=0; i<*n; i++) { if ( upc_threadof( &Cells[j*(*n)+i] ) == MYTHREAD )
    {

      if (Fluid[j*(ln)+i] == 1)
      {
        // INLET
        InletBC(j, i, lm, ln);

        // WALL
        switch(CurvedBoundaries)
        {
          // curved boundaries
          case 1:
            puts("CurvedBoundaries!!!");
            for(k=0; k<9; k++)
              Fneighbours[k][j*(ln)+i] = METAF[k][(j-cx[k])*(ln) + i-cy[k]];
            // C++ code
            // Cells[j][i].setFneighbours(Cells[j-INI.getcy()[k]][i-1*INI.getcx()[k]].getMETAF()[k],k);
          
            CurvedWallBoundaries(j, i, opp, lm, ln);
          break;

          // bounceback boundaries
          case 2:
            WallBC(j, i, opp, lm, ln);
          break;
        }

        // OUTLET
        switch(OutletProfile)
        {
          // set profile in outlet
          case 1:
            OutletBoundaries(j, i, lm, ln);
          break;
          
          // OPEN BOUNDARY
          case 2 :
            if (BC_ID[1][j*(ln)+i]==3)
            { 
              
              F[1][j*(ln)+i] = 2*( F[1][(ln)*(j)+i-1] ) - F[1][(ln)*(j)+i-2]; 
              F[5][j*(ln)+i] = 2*( F[5][(ln)*(j)+i-1] ) - F[5][(ln)*(j)+i-2];
              F[8][j*(ln)+i] = 2*( F[8][(ln)*(j)+i-1] ) - F[8][(ln)*(j)+i-2];

              // C++ code
              //Cells[j][i].setF(2*Cells[j][i-1].getF()[1]-Cells[j][i-2].getF()[1],1);
              //Cells[j][i].setF(2*Cells[j][i-1].getF()[5]-Cells[j][i-2].getF()[5],5);
              //Cells[j][i].setF(2*Cells[j][i-1].getF()[8]-Cells[j][i-2].getF()[8],8);

            }
            if (BC_ID[2][j*(ln)+i]==3)
            {
              // FILL!!
            }
            if (BC_ID[3][j*(ln)+i]==3)
            {
              // FILL!!
            }
            if (BC_ID[4][j*(ln)+i]==3)
            {  
              // FILL!!
            }
          break;
        }
      } 
    } //}
  }
}

/////////////// UPDATE MACRO STEP ///////////////////
void UpdateMacroscopicStep(int* cx, int* cy, int CalculateDragLift, int lm, int ln, int MyAffinity[NN*MM])
{
  int i, j;
  upc_forall(i=0; i<((lm)*(ln)); i++; MyAffinity[i])
  //for(i=0; i<(*n)*(*m); i++) { if ( upc_threadof( &Cells[i] ) == MYTHREAD )
  {
    if (Fluid[i] == 1)
      UpdateMacroscopic(i, cx, cy, CalculateDragLift);
  } //}
}


