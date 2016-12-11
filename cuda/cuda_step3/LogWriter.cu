#include <stdio.h>
#include <errno.h>
#include "LogWriter.h"

void writeInitLog(const char* filename, Arguments *args, FLOAT_TYPE delta, int m, int n,
                  int numInletNodes, FLOAT_TYPE maxInletCoordY, FLOAT_TYPE minInletCoordY)
{
  FILE* logFile = fopen(filename, "w");  // open log file
  fprintf(logFile,"This is the 2D lattice Boltzmann *.log file\n\n");
  fprintf(logFile,"\n:::: Imported variables from the *.ini file :::: \n");
  fprintf(logFile,">>> Uavg              : %3.6f\n", args->u);
  fprintf(logFile,">>> Vavg              : %3.6f\n", args->v);
  fprintf(logFile,">>> Initial density   : %2.1f\n", args->rho);
  fprintf(logFile,">>> Viscosity         : %3.8f\n", args->viscosity);
  fprintf(logFile,">>> # of iterations   : %1.1d\n", args->iterations);
  fprintf(logFile,">>> Autosave after    : %1.1d\n", args->autosaveAfter);
  fprintf(logFile,">>> Autosave every    : %1.1d\n", args->autosaveEvery);

  switch(args->collisionModel)
  {
    case BGKW: fprintf(logFile,">>> CollisionModel   : BGKW\n"); break;
    case TRT : fprintf(logFile,">>> CollisionModel   : TRT\n" ); break;
    case MRT : fprintf(logFile,">>> CollisionModel   : MRT\n" ); break;
  }
  switch(args->inletProfile)
  {
    case INLET          : fprintf(logFile,">>> InletProfile     : ON\n" ); break;
    case NO_INLET       : fprintf(logFile,">>> InletProfile     : OFF\n"); break;
    case PULSATILE_INLET: fprintf(logFile,">>> InletProfile     : PULSATILE\n"); break;
  }
  switch(args->outletProfile)
  {
    case OUTLET       : fprintf(logFile,">>> OutletProfile    : ON\n" ); break;
    case OUTLET_SECOND: fprintf(logFile,">>> OutletProfile    : OPEN (2nd)\n"); break;
    case OUTLET_FIRST : fprintf(logFile,">>> OutletProfile    : OPEN (1st)\n"); break;
  }
  switch(args->boundaryType)
  {
    case CURVED  : fprintf(logFile,">>> CurvedBoundaries : ON\n" ); break;
    case STRAIGHT: fprintf(logFile,">>> CurvedBoundaries : OFF\n"); break;
  }
  switch(args->outputFormat)
  {
    case PARAVIEW: fprintf(logFile,">>> Results format   : Paraview (*.csv)\n" ); break;
    case TECPLOT : fprintf(logFile,">>> Results format   : Tecplot (*.dat)\n"); break;
  }
  if (args->boundaryId > 0)
            fprintf(logFile,">>> Drag, lift @ BC  : %d\n", args->boundaryId);
  else
            fprintf(logFile,">>> Drag, lift was not calculated\n");

  fprintf(logFile,"\n:::: Calculated variables from mesh :::: \n");
  fprintf(logFile,">>> Grid spacing             = %f\n", delta);
  fprintf(logFile,">>> # of nodes in x (n)      = %d\n", n);
  fprintf(logFile,">>> # of nodes in y (m)      = %d\n", m);
  fprintf(logFile,">>> NumInletNodes            = %d\n", numInletNodes);
  fprintf(logFile,">>> MaxInletCoordY           = %f\n", maxInletCoordY);
  fprintf(logFile,">>> MinInletCoordY           = %f\n", minInletCoordY);
  fclose(logFile);
}

void writeNodeNumbers(const char* filename, int numNodes, int numConns, int bcCount)
{
  FILE* logFile = fopen(filename, "a");  // open log file
  fprintf(logFile,">>> Number of nodes          = %d\n", numNodes);
  fprintf(logFile,">>> Number of conditions     = %d\n", numConns);
  fprintf(logFile,">>> Number of boundary nodes = %d\n", bcCount);
  fclose(logFile);
}

void writeEndLog(const char *filename, FLOAT_TYPE *taskTimes)
{
  FILE* logFile = fopen(filename, "a");
  fprintf(logFile,"\nOverall calculations took %f seconds",   taskTimes[T_OALL]);
  fprintf(logFile,"\nMain while loop took %f seconds\n",      taskTimes[T_ITER]);
  fprintf(logFile,"Initialization took %f seconds\n",         taskTimes[T_INIT]);
  fprintf(logFile,"Collision took %f seconds\n",              taskTimes[T_COLL]);
  fprintf(logFile,"Streaming took %f seconds\n",              taskTimes[T_STRM]);
  fprintf(logFile,"Calculating Boundaries took %f seconds\n", taskTimes[T_BNDC]);
  fprintf(logFile,"Update Macroscopic took %f seconds\n",     taskTimes[T_MACR]);
  fprintf(logFile,"Calculating residuals took %f seconds\n",  taskTimes[T_RESI]);
  fprintf(logFile,"Writing results took %f seconds\n",        taskTimes[T_WRIT]);

  fprintf(logFile,"\n:::: Iterations done! ::::\n");
  fclose(logFile);
}

void writeTimerLog(const char *filename, FLOAT_TYPE *taskTimes)
{
  FILE* timeFile = fopen(filename, "w");
  fprintf(timeFile,"Overall %f\n",        taskTimes[T_OALL]);
  fprintf(timeFile,"Iteration %f\n",      taskTimes[T_ITER]);
  fprintf(timeFile,"Initialisation %f\n", taskTimes[T_INIT]);
  fprintf(timeFile,"Collision %f\n",      taskTimes[T_COLL]);
  fprintf(timeFile,"UpdateF %f\n",        0.000000000000000);
  fprintf(timeFile,"Streaming %f\n",      taskTimes[T_STRM]);
  fprintf(timeFile,"Boundaries %f\n",     taskTimes[T_BNDC]);
  fprintf(timeFile,"Macroscopic %f\n",    taskTimes[T_MACR]);
  fprintf(timeFile,"Residuals %f\n",      taskTimes[T_RESI]);
  fprintf(timeFile,"FileWriting %f\n",   taskTimes[T_WRIT]);
  fclose(timeFile);
}

void writeResiduals(const char *filename, FLOAT_TYPE *norm, FLOAT_TYPE *drag, FLOAT_TYPE *lift, int size, int n)
{
  FILE *residualsFile = fopen(filename, "w");
  if (!residualsFile)
  {
    fprintf(stderr, "Error reading file %s: %s\n", filename, strerror(errno));
    return;
  }
  fprintf(residualsFile, "Iter L2_norm L2_norm_weighted Drag Lift\n");
  int i;
  for (i=0; i<n; ++i)
  {
    fprintf(residualsFile,"%d %5.4e %5.4e %f %f\n", i, norm[i], norm[i]/(FLOAT_TYPE)size, drag[i], lift[i]);
  }
  fclose(residualsFile);
}
