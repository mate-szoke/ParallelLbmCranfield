#include <stdio.h>   // for calloc();
#include "include/ShellFunctions.h"
#include "math.h"

void WriteResults(char* OutputFile, struct CellProps_const *Cells_const_h,
                  struct CellProps_var *Cells_var_h, struct CellProps_var_9d *Cells_var_9d_h, int* n, int* m, int* postproc_prog)
{
  int i;                   // Loop variables
  FILE * fp1;                 // file pointer to output file
  fp1=fopen(OutputFile, "w"); // open file
  switch(*postproc_prog)
  {
    case 1: // ParaView

      fprintf(fp1, "x,y,u,v,vel_mag,rho,press,fluid\n");
      for(i=0;i<((*m)*(*n));i++)
      {
        
        fprintf(fp1, "%f, %f, %f, %f, %f, %f, %f, %d\n",
                Cells_const_h[i].CoordX, // x
                Cells_const_h[i].CoordY, // y
                Cells_var_h[i].U,      // u
                Cells_var_h[i].V,      // v
                sqrt(pow(Cells_var_h[i].U,2)+pow(Cells_var_h[i].V,2)), // u magnitude
                Cells_var_h[i].Rho,    // density
                Cells_var_h[i].Rho/3,  // pressure
                Cells_const_h[i].Fluid);
      }
      fclose(fp1);
    break;

  	case 2: // TECPLOT
    	fprintf(fp1, "Title = \"LBM results\"\n");
    	fprintf(fp1, "Variables = \"x\",\"y\",\"u\",\"v\",\"vel mag\",\"rho\",\"press\",\"fluid\"\n");
    	fprintf(fp1, "Zone i=%d, j=%d, f=point\n",*n,*m);

    	for(i=0;i<((*m)*(*n));i++)
    	{

   	    fprintf(fp1, "%f %f %f %f %f %f %f %d\n",
          	    Cells_const_h[i].CoordX, // x
                Cells_const_h[i].CoordY, // y
          	    Cells_var_h[i].U,      // u
                Cells_var_h[i].V,      // v
          	    sqrt(pow(Cells_var_h[i].U,2)+pow(Cells_var_h[i].V,2)), // u magnitude
          	    Cells_var_h[i].Rho,    // density
                Cells_var_h[i].Rho/3,  // pressure
                Cells_const_h[i].Fluid); // fluid or solid
     	}

      fclose(fp1); 
  	break;
  }
	
	
}
