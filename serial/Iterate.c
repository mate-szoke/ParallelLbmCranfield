#include <stdio.h>                      // printf();
#include <math.h>                       // need to compile with -lm
#include <stdlib.h>                     // for calloc();
#include <stdbool.h>                    // Include for bool type variables!
#include <string.h>                     // String operations
#include <time.h>                       // time functions

#include "include/Iterate.h"            // Iteration takes place
#include "include/ShellFunctions.h"     // For convenience
#include "include/FilesIO.h"            // For reading & writing files
//#include "include/FilesWriting.h"       // For writing files e.g. tecplot
#include "include/CellFunctions.h"      // For cell modifications
#include "include/BoundaryConditions.h" // boundary conditions
#include "include/ComputeResiduals.h"   // Residuals


void Iteration(char* NodeDataFile, char* BCconnectorDataFile,
               float Uavg,         float Vavg,         float rho_ini, float Viscosity,
               int InletProfile,   int CollisionModel, int CurvedBoundaries,
               int OutletProfile,  int Iterations,     int AutosaveAfter,
               int AutosaveEvery,  int postproc_prog,  int CalculateDragLift,
               float ConvergenceCritVeloc, float ConvergenceCritRho, int ReadFormerData)
{

  ////////////////////////////////////////////////////
  ///////////////////// Declare //////////////////////
  ////////////////////////////////////////////////////
  
  // Time measurement: declaration, begin
  clock_t tStart = clock();

  FILE* resid_file;                       // file for residuals
  FILE* log_file;                         // file for log
  FILE* TimeMeasurementFile;              // file for time measurement results
  char IterOutputFile[50];                // write results to this file after the iterations
  char AutosaveOutputFile[50];            // autosave filename
  char OutputFile[50];                    // initial data will be written to this file
  char FinalOutputFile[50];               // final data will be written to this file
  char logFile[] = "Results/logFile.log"; // path of the .log file
  int i, j, k, iter = 0;                  // variables for loops

  // Variables for residuals
  MyReal *Residuals;

  int AutosaveI = 1;       // autosave i variable, will be incremented after every autosave
  int* ppp;                // pointer of the postproc_prog variable
  int *NumNodes,*NumConn;  // This will store the number of lines of the read files
  MyReal *Delta;           // grid spacing
  int *n,*m;               // number of nodes in the x and y directions
  MyReal *MaxInletCoordY;  // maximum inlet coordinate in y
  MyReal *MinInletCoordY;  // minimum inlet coordinate in y
  int *NumInletNodes;      // number of inlet nodes
  float **Nodes,**BCconn;  // matrices for the nodes and connections
  MyReal Omega, OmegaA;    // collision frequency from the viscosity
  MyReal **tm, **stmiv;    // variables for the MRT collision model

  float tInitialization  = 0.0;  // Time measurement of Initialization
  float tIteration       = 0.0;  // Time measurement of Iteration
  float tCollision       = 0.0;  // Time measurement of Collision
  float tUpdateF         = 0.0;  // Time measurement of UpdateF
  float tStreaming       = 0.0;  // Time measurement of Streaming
  float tBoundaries      = 0.0;  // Time measurement of Boundaries
  float tUpdateMacro     = 0.0;  // Time measurement of Update Macroscopic vars
  float tResiduals       = 0.0;  // Time measurement of calculating residuals
  float tWriting         = 0.0;  // Time measurement of writing data
  clock_t tInstant1, tInstant2;  // Time measurement points, universal
  clock_t tIterStart, tIterEnd;  // Time measurement points: main loop
  
  ////////////////////////////////////////////////////
  //////////////////// Allocate //////////////////////
  ////////////////////////////////////////////////////

  // allocate residuals
  Residuals = Create1DArrayMyReal(4); 

  // allocate mesh properties  
  Delta          = Create1DArrayMyReal(1);
  m              = Create1DArrayInt(1);
  n              = Create1DArrayInt(1);
  NumNodes       = Create1DArrayInt(1);
  NumConn        = Create1DArrayInt(1);
  MaxInletCoordY = Create1DArrayMyReal(1);
  MinInletCoordY = Create1DArrayMyReal(1);
  NumInletNodes  = Create1DArrayInt(1);

  // D2Q9 Variables of the lattice
  MyReal* w = Create1DArrayMyReal(9); // weight values for the directions
  int*    c = Create1DArrayInt(9);   
  int*   cx = Create1DArrayInt(9);    // x coordinate of the discrete lattice directions
  int*   cy = Create1DArrayInt(9);    // y coordinate of the discrete lattice directions
  int*  opp = Create1DArrayInt(9);    // opposite vector

  ////////////////////////////////////////////////////
  ///////////////////// Read data ////////////////////
  ////////////////////////////////////////////////////

  // Read Node data
  Nodes = ReadNodes(NodeDataFile, NumNodes); 
  
  // Read BCconn data
  BCconn  = ReadBCconn(BCconnectorDataFile, NumConn); 

  ////////////////////////////////////////////////////
  /////////////// Compute constants //////////////////
  ////////////////////////////////////////////////////

  CompDataNode(Delta, m,  n, Nodes, NumNodes);

  CompDataConn(NumInletNodes, MaxInletCoordY,
               MinInletCoordY, BCconn, NumConn, Delta);

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
  fprintf(log_file,">>> Former data iter# : %1.1d\n", ReadFormerData);
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
  // End of checking

  ////////////////////////////////////////////////////
  ///////////////// INITIALIZE ///////////////////////
  ////////////////////////////////////////////////////

  // In case of no autosave
  sprintf(AutosaveOutputFile, "NOWHERE!");

  // Fill up D2Q9 variables
  D2Q9Vars(w, cx, cy, opp, c, n);
  
  // Calculate collision freq
  Omega  = 1.0/(3.*Viscosity+0.5);
  OmegaA = 8*(2-Omega)/(8-Omega);

  // Initialize variables for MRT Collision model, if used
  tm    = Create2DArrayMyReal(9, 9);
  stmiv = Create2DArrayMyReal(9, 9);
  if (CollisionModel == 3)
  {
    MRTInitializer(tm, stmiv, Omega);
  }  

  // Create structure for the cell properties (see ShellFunctions.h)
  struct CellProps *Cells; 
  // allocate...
  Cells = calloc((*n)*(*m),sizeof(struct CellProps));


  // initializing matrix of struct Cells
  fprintf(log_file,"\n:::: Initializing ::::\n");
  printf("\n:::: Initializing ::::\n");
  tInstant1 = clock(); // Measure time of initialization

//  printf("Ny %d \n", *n);
//  printf("Nx %d \n", *m);

  for(j=0;j<*m;j++)
  {
    for(i=0;i<*n;i++)
    {
      CellIni(Nodes,            // Nodes
    	BCconn,           // BCconn
    	NumNodes,         // Number of lines in Nodes
	   	NumConn,          // Number of lines in BCconn
		  i,                // A in the functions
		  j,                // B in the functions
		  MinInletCoordY,   // min inlet y coordinate
		  MaxInletCoordY,   // max inlet y coordinate
		  Delta,            // grid spacing 
		  Uavg,             // INPUT PARAMETER
		  Vavg,             // INPUT PARAMETER
		  InletProfile,     // INPUT PARAMETER
		  CollisionModel,   // INPUT PARAMETER
		  opp,              // Opposite
		  Cells,            // Struct of Cells
		  rho_ini,          // initial density    
		  n);		  // number of nodes in the
//    printf("Index %d \n", i);
    printf("Initializing... %2.1f %%\r",(float)j*100/(float)(*m));
    }
  }



  tInstant2 = clock(); // Measure time of initialization
  tInitialization = (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

  fprintf(log_file,"\n:::: Initialization done! ::::\n");
  // Write Initialized data 
  switch(postproc_prog)
    {
      case 1: sprintf(OutputFile, "Results/InitialData.csv");   break;
      case 2: sprintf(OutputFile, "Results/InitialData.dat");   break;
    }

  WriteResults(OutputFile, Cells, n, m, ppp);

  printf("\nInitialized data was written to %s\n", OutputFile);

   // Open residuals file
   resid_file = fopen("Results/residuals.dat", "w");
   fprintf(resid_file,"Iter L2_norm L2_norm_weighted Drag Lift\n");

   ////////////////////////////////////////////////////
   /////////////////// ITERATE ////////////////////////
   ////////////////////////////////////////////////////

   fprintf(log_file,"\n:::: Start Iterations ::::\n");
   printf("\n:::: Start Iterations ::::\n");

   // Set residuals to be big
   Residuals[0] = 100.0;
   Residuals[1] = 100.0;
  
   tIterStart = clock(); // Start measuring time of main loop  
   while (Residuals[0] > ConvergenceCritVeloc && Residuals[1] > ConvergenceCritRho && iter<Iterations)
   {

     ////////////// COLLISION ///////////////
     tInstant1 = clock();
     CollisionStep( Cells, m,  n,  w, cx, cy,  opp,  Omega,  OmegaA, tm, stmiv, CollisionModel);
     tInstant2 = clock();
     tCollision = tCollision + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;


     ////////////// UPDATE DISTR ///////////////
     tInstant1 = clock(); // Start measuring time
     for (i=0; i<(*m)*(*n); i++)
       { UpdateF(Cells, i); }
     tInstant2 = clock();
     tUpdateF = tUpdateF + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;
     
     ////////////// STREAMING ///////////////
     tInstant1 = clock(); // Start measuring time
     StreamingStep(Cells, m, n, c);
     tInstant2 = clock();
     tStreaming = tStreaming + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

     ////////////// BOUNDARIES ///////////////
     tInstant1 = clock(); // Start measuring time
     HandleBoundariesStep(Cells, cx, cy, c, opp, OutletProfile, CurvedBoundaries, n, m);
     tInstant2 = clock();
     tBoundaries = tBoundaries + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

     // UPDATE VELOCITY AND DENSITY
     tInstant1 = clock(); // Start measuring time
     UpdateMacroscopicStep(Cells, m, n, cx, cy, CalculateDragLift);
     tInstant2 = clock();
     tUpdateMacro = tUpdateMacro + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;
    
     ////////////// Residuals ///////////////
     tInstant1 = clock(); // Start measuring time
     ComputeResiduals(Cells, Residuals, m, n, CalculateDragLift);
     fprintf(resid_file,"%d %5.4e %5.4e %f %f\n", iter, Residuals[0], Residuals[1], Residuals[2], Residuals[3]);
     tInstant2 = clock();
     tResiduals = tResiduals + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;


     // update loop variable
     iter++;  

     if(iter%100==0){
       printf("Iterations: %05d/%05d || ", iter, Iterations);
       printf("Residuals: l2norm  %e; L2_norm_weighted  %e\n", Residuals[0], Residuals[1]);
     }
    
     ////////////// Autosave ///////////////
     //if(iter == (AutosaveEvery*AutosaveI))
     if(iter%AutosaveEvery == 0)
     {
       AutosaveI++;
       if(iter>AutosaveAfter)
       {
         switch(postproc_prog){
           case 1: sprintf(AutosaveOutputFile, "Results/autosave_iter%05d.csv", iter);  break;
           case 2: sprintf(AutosaveOutputFile, "Results/autosave_iter%05d.dat", iter);   break;
         }
         tInstant1 = clock(); // Start measuring time
         WriteResults(AutosaveOutputFile, Cells, n, m, ppp);
         tInstant2 = clock();
         tWriting = tWriting + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;
       }
     }

   }     ////////////// END OF MAIN WHILE CYCLE! ///////////////


	 tIterEnd = clock(); // End measuring time of main loop
   tIteration = (float)(tIterEnd - tIterStart ) / CLOCKS_PER_SEC;

   clock_t tEnd = clock();
   float tOverall = (float)(tEnd - tStart) / CLOCKS_PER_SEC; // Calculate overall elapsed time

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
  
   // Close residuals file
   fclose(resid_file);
  
   // Write the time measurements to a separate dat file
   TimeMeasurementFile = fopen("Results/SerialTimeMeasuerment.dat","w");
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
    fclose(TimeMeasurementFile);

   // Write final data
   switch(postproc_prog){
     case 1: sprintf(FinalOutputFile, "Results/FinalData.csv"); break;
     case 2: sprintf(FinalOutputFile, "Results/FinalData.dat"); break;
   }
   WriteFinalResults(FinalOutputFile, Cells, n, m, ppp);


   // Write information for user
   printf("\n\nLog was written to %s\n", logFile);
   printf("Last autosave result can be found at %s\n", AutosaveOutputFile);
   printf("Residuals were written to Results/residuals.dat\n");
   printf("Profiling results were written to Results/SerialTimeMeasuerment.dat\n");
   printf("Final results were written to %s\n", FinalOutputFile);

  ////////////////////////////////////////////////////
  ///////////////// End of line //////////////////////
  ////////////////////////////////////////////////////

  // FREE POINTERS
  free(Delta);
  free(m);
  free(n);
  free(MaxInletCoordY);
  free(MinInletCoordY);
  free(NumInletNodes);
  free(NumNodes);
  free(NumConn);
  free(Nodes);
  free(BCconn);
  free(w);
  free(cx);
  free(cy);
  free(opp);
  free(Cells);
  free(Residuals);

}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////




// ////////////////// COLLISION STEP ///////////////////
 void CollisionStep(struct CellProps *Cells, int* m, int* n, MyReal* w,
                int* cx, int* cy, int* opp, MyReal Omega, MyReal OmegaA,
                MyReal **tm, MyReal **stmiv, int CollisionModel)
 {

   int i;

   switch(CollisionModel)
   {
     // BGKW
     case 1:

       for (i=0; i<(*m)*(*n);i++)
       {
         if ((Cells+i)->Fluid == 1)
         {
           BGKW(Cells, i, w, cx, cy, Omega);
         }
       }
       
       break;

       // TRT 
     case 2:
       for (i=0; i<(*m)*(*n);i++)
       {
         if ((Cells+i)->Fluid == 1)
         {
           TRT (Cells, i, w, cx, cy, opp, Omega, OmegaA);
         }
       }
       break;

       // MRT 
     case 3:
       for (i=0; i<(*m)*(*n);i++)
       {
         if ((Cells+i)->Fluid == 1)
         {
           MRT (Cells, i, tm, stmiv); 
         }
       }
       break;

   }


 }

 ////////////////// STREAMING STEP ///////////////////
 void StreamingStep(struct CellProps *Cells, int* m, int* n,  int* c)
 {
   int i, k;

     for (i=0; i<(*m)*(*n); i++)
     {
       if ( (Cells+i)->Fluid == 1 )
       {
         for(k=0; k<9; k++)
         {
           if ( ((Cells+i)->StreamLattice[k]) == 1)
           {
             //(Cells+i)->F[k] = Cells [j-cx[k]] [i-cy[k]].METAF[k]; OLD CODE, 2D VERSION
             (Cells +i)->F[k] = ( Cells+i+c[k] )-> METAF[k];
           }
         }
       }
     }
 }

 ///////////////// BOUNDARIES STEP ///////////////////
void HandleBoundariesStep(struct CellProps *Cells, int* cx, int* cy, int* c, int* opp, int OutletProfile, int CurvedBoundaries, int* n, int* m)
{
  int i, j, k;

  for(j=0; j<*m; j++)
  {
    for(i=0; i<*n; i++)
    {
      if ((Cells +j*(*n)+i)->Fluid == 1)
      {
        // INLET
        InletBC(Cells, j, i, n, m);

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
          
            CurvedWallBoundaries(Cells, j, i, opp, n, m);
          break;

          // bounceback boundaries
          case 2:
            WallBC(Cells, j, i, opp, n, m);
          break;
        }

        // OUTLET
        switch(OutletProfile)
        {
          // set profile in outlet
          case 1:
            OutletBoundaries(Cells, j, i, n, m);
          break;
          
          // OPEN BOUNDARY
          case 2 :
            if ((Cells +j*(*n)+i)->BC_ID[1]==3)
            { 
              
              (Cells +j*(*n)+i)->F[1] = 2*( (Cells+(*n)*(j)+i-1)->F[1] ) - (Cells+(*n)*(j)+i-2)->F[1]; 
              (Cells +j*(*n)+i)->F[5] = 2*( (Cells+(*n)*(j)+i-1)->F[5] ) - (Cells+(*n)*(j)+i-2)->F[5];
              (Cells +j*(*n)+i)->F[8] = 2*( (Cells+(*n)*(j)+i-1)->F[8] ) - (Cells+(*n)*(j)+i-2)->F[8];


            }

          break;
        }
      } 
    }
  }
}


 /////////////// UPDATE MACRO STEP ///////////////////
 void UpdateMacroscopicStep(struct CellProps *Cells, int* m, int* n,  int* cx,
                           int* cy, int CalculateDragLift)
 {
   int i, j;

     for (i=0; i<(*m)*(*n);i++)
     {
       if ((Cells+i)->Fluid==1)
       {
         UpdateMacroscopic(Cells, j, i, cx, cy, CalculateDragLift);
       }
     }

 }

 ///////////////// LIFTDRAG ONLY /////////////////////  >>>>>>>>> OUT OF USE !!!! 
 void CalculateDragLiftForcesStep(struct CellProps *Cells,  int* m, int* n, int CalculateDragLift)
 {

   int i, j;

   for (i=0; i<(*m)*(*n);i++)
     {
       if ((Cells+i)->Fluid==1)
       {
         CalculateDragLiftForces(Cells, j, i, CalculateDragLift);
       }
     }

 }
