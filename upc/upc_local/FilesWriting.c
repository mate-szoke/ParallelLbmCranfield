#include <stdio.h>   // for calloc();
#include "math.h"
#include <upc_relaxed.h>                 // Required for UPC 

#include "include/ShellFunctions.h"


void WriteResults(char* OutputFile, int* postproc_prog)
{
  int i, j;                   // Loop variables
  FILE * fp1;                 // file pointer to output file
  fp1=fopen(OutputFile, "w"); // open file
  switch(*postproc_prog)
  {
    case 1: // ParaView
    	fprintf(fp1, "x,y,u,v,vel_mag,rho,press,fluid,ThID\n");
    	j = 0;
      //for(j=0;j<*m;j++)
    	{
          for(i=0; i<((*n)*(*m));i++)
          {
    	    fprintf(fp1, "%f, %f, %f, %f, %f, %f, %f, %d, %d\n",
                  (WCells+j*(*n)+i)->CoordX, // x
                  (WCells+j*(*n)+i)->CoordY, // y
                  (WCells+j*(*n)+i)->U,      // u
                  (WCells+j*(*n)+i)->V,      // v
                  sqrt(pow((WCells+j*(*n)+i)->U,2)+pow((WCells+j*(*n)+i)->V,2)), // u magnitude
                  (WCells+j*(*n)+i)->Rho,    // density
                  ((WCells+j*(*n)+i)->Rho)/3,  // pressure
                  (WCells+j*(*n)+i)->Fluid, // fluid or solid
                  (WCells+i)->ThreadNumber/*
                  (Cells+j*(*n)+i)->F[1],
                  (Cells+j*(*n)+i)->F[5],
                  (Cells+j*(*n)+i)->F[6],
                  (Cells+j*(*n)+i)->ID,
                  (Cells+i)->Boundary*/);
          }
      	}

    	fclose(fp1);
    break;

  	case 2: // TECPLOT
    	fprintf(fp1, "Title = \"LBM results\"\n");
    	fprintf(fp1, "Variables = \"x\",\"y\",\"u\",\"v\",\"u mag\",\"rho\",\"press\",\"fluid\"\n");
    	fprintf(fp1, "Zone i=%d, j=%d, f=point\n",*n,*m);

    	for(j=0;j<*n;j++)
    	{
          for(i=0;i<*m;i++)
          {
    	    fprintf(fp1, "%f %f %f %f %f %f %f %d\n",
            	    (WCells+j*(*m)+i)->CoordX, // x
                  (WCells+j*(*m)+i)->CoordY, // y
            	    (WCells+j*(*m)+i)->U,      // u
                  (WCells+j*(*m)+i)->V,      // v
            	    sqrt(pow((WCells+j*(*m)+i)->U,2)+pow((WCells+j*(*m)+i)->V,2)), // u magnitude
            	    (WCells+j*(*m)+i)->Rho,    // density
                  ((WCells+j*(*m)+i)->Rho)/3,  // pressure
                  (WCells+j*(*m)+i)->Fluid); // fluid or solid
          }
      	}

      fclose(fp1); 
  	break;
  }
	
	
}


void WriteBCells(char* OutputFile, int* postproc_prog)
{
  int i, j;                   // Loop variables
  FILE * fp1;                 // file pointer to output file
  fp1=fopen(OutputFile, "w"); // open file
  switch(*postproc_prog)
  {
    case 1: // ParaView
      fprintf(fp1, "x,y,u,v,vel_mag,rho,press,fluid,ThID\n");
          for(i=0;i<(2*THREADS*(*n));i++)
          {
            fprintf(fp1, "%f, %f, %f, %f, %f, %f, %f, %d, %d\n",
                    (BCells+i)->CoordX, // x
                    (BCells+i)->CoordY, // y
                    (BCells+i)->U,      // u
                    (BCells+i)->V,      // v
                    sqrt(pow((BCells+i)->U,2)+pow((BCells+i)->V,2)), // u magnitude
                    (BCells+i)->Rho,    // density
                    ((BCells+i)->Rho)/3,  // pressure
                    (BCells+i)->Fluid, // fluid or solid
                    (BCells+i)->ThreadNumber/*
                    (Cells+j*(*n)+i)->F[1],
                    (Cells+j*(*n)+i)->F[5],
                    (Cells+j*(*n)+i)->F[6],
                    (Cells+j*(*n)+i)->ID,
                    (Cells+i)->Boundary*/);
          }

      fclose(fp1);
    break;

    case 2: // TECPLOT
      fprintf(fp1, "Title = \"LBM results\"\n");
      fprintf(fp1, "Variables = \"x\",\"y\",\"u\",\"v\",\"u mag\",\"rho\",\"press\",\"fluid\"\n");
      fprintf(fp1, "Zone i=%d, j=%d, f=point\n",*n,*m);

      for(j=0;j<*n;j++)
      {
          for(i=0;i<*m;i++)
          {
          fprintf(fp1, "%f %f %f %f %f %f %f %d\n",
                  (BCells+j*(*m)+i)->CoordX, // x
                  (BCells+j*(*m)+i)->CoordY, // y
                  (BCells+j*(*m)+i)->U,      // u
                  (BCells+j*(*m)+i)->V,      // v
                  sqrt(pow((BCells+j*(*m)+i)->U,2)+pow((BCells+j*(*m)+i)->V,2)), // u magnitude
                  (BCells+j*(*m)+i)->Rho,    // density
                  ((BCells+j*(*m)+i)->Rho)/3,  // pressure
                  (BCells+j*(*m)+i)->Fluid); // fluid or solid
          }
        }

      fclose(fp1); 
    break;
  }
  
  
}