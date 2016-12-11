#include <stdio.h>   // for calloc();
#include "FilesWriting.h"
#include <math.h>

void WriteResults(char* OutputFile, FLOAT_TYPE* CoordX, FLOAT_TYPE* CoordY, FLOAT_TYPE* U,
                  FLOAT_TYPE* V, FLOAT_TYPE* Rho, int* Fluid, int n, int m, OutputFormat outputFormat)
{
  int i;                      // Loop variables
  FILE * fp1;                 // file pointer to output file
  fp1 = fopen(OutputFile, "w"); // open file
  switch(outputFormat)
  {
    case PARAVIEW:

      fprintf(fp1, "x,y,u,v,vel_mag,rho,press,fluid\n");
      for(i=0; i<m*n; ++i)
      {

        fprintf(fp1, "%f, %f, %f, %f, %f, %f, %f, %d\n",
                CoordX[i], // x
                CoordY[i], // y
                U[i],      // u
                V[i],      // v
                sqrt(pow(U[i],2)+pow(V[i],2)), // u magnitude
                Rho[i],    // density
                Rho[i]/3,  // pressure
                Fluid[i]);
      }
      fclose(fp1);
    break;

  	case TECPLOT:
    	fprintf(fp1, "Title = \"LBM results\"\n");
    	fprintf(fp1, "Variables = \"x\",\"y\",\"u\",\"v\",\"vel mag\",\"rho\",\"press\",\"fluid\"\n");
    	fprintf(fp1, "Zone i=%d, j=%d, f=point\n",n,m);

    	for(i=0; i<m*n; ++i)
    	{

        fprintf(fp1, "%f %f %f %f %f %f %f %d\n",
                CoordX[i], // x
                CoordY[i], // y
                U[i],      // u
                V[i],      // v
                sqrt(pow(U[i],2)+pow(V[i],2)), // u magnitude
                Rho[i],    // density
                Rho[i]/3,  // pressure
                Fluid[i]);
     	}

      fclose(fp1);
  	break;
  }


}
