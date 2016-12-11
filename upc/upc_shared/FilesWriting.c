#include <stdio.h>   // for calloc();
#include "math.h"
//#include <upc_cray.h>                 // Required for UPC 

#include "include/ShellFunctions.h"


void WriteResults(char* OutputFile, int* postproc_prog, int lm, int ln)
{
  int i, j=0;                   // Loop variables
  FILE * fp1;                 // file pointer to output file
  fp1=fopen(OutputFile, "w"); // open file
  switch(*postproc_prog)
  {
    case 1: // ParaView
    	fprintf(fp1, "x,y,u,v,vel_mag,rho,press,fluid,ThID\n");
      //for(j=0;j<lm;j++)
    	{
          for(i=0;i<((ln)*(lm));i++)
          {
    	    fprintf(fp1, "%f, %f, %f, %f, %f, %f, %f, %d, %d\n",
                  CoordX[j*(ln)+i], // x
                  CoordY[j*(ln)+i], // y
                  U[j*(ln)+i],      // u
                  V[j*(ln)+i],      // v
                  sqrt(pow(U[j*(ln)+i],2)+pow(V[j*(ln)+i],2)), // u magnitude
                  Rho[j*(ln)+i],    // density
                  (Rho[j*(ln)+i])/3,  // pressure
                  Fluid[j*(ln)+i], // fluid or solid
                  ThreadNumber[j*(ln)+i]/*
                  F[1][j*(ln)+i],
                  F[5][j*(ln)+i],
                  F[6][j*(ln)+i],
                  ID[j*(ln)+i],
                  (Cells+i)->Boundary[j*(ln)+i]*/);
          }
      	}

    	fclose(fp1);
    break;

  	case 2: // TECPLOT
    	fprintf(fp1, "Title = \"LBM results\"\n");
    	fprintf(fp1, "Variables = \"x\",\"y\",\"u\",\"v\",\"u mag\",\"rho\",\"press\",\"fluid\",\"ThID\"\n");
    	fprintf(fp1, "Zone i=%d, j=%d, f=point\n",ln,lm);

    	//for(j=0;j<ln;j++)
    	{
          for(i=0;i<((ln)*(lm));i++)
          {
    	    fprintf(fp1, "%f %f %f %f %f %f %f %d %d\n",
            	    CoordX[j*(ln)+i], // x
                  CoordY[j*(ln)+i], // y
            	    U[j*(ln)+i],      // u
                  V[j*(ln)+i],      // v
            	    sqrt(pow(U[j*(ln)+i],2)+pow(V[j*(ln)+i],2)), // u magnitude
            	    Rho[j*(ln)+i],    // density
                  (Rho[j*(ln)+i])/3,  // pressure
                  Fluid[j*(ln)+i], // fluid or solid
                  ThreadNumber[j*(ln)+i]);
          }
      	}

      fclose(fp1); 
  	break;
  }
	
	
}

