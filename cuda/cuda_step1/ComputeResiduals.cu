#include <stdio.h> 
#include <stdlib.h> 
#include "math.h"
#include "include/ShellFunctions.h"

void ComputeResiduals(struct CellProps_var *Cells_var_h, struct CellProps_const *Cells_const_h, float* Residuals, int* m, int* n,
          float* sumVel0, float* sumVel1, float* sumRho0, float* sumRho1, int ComputeDragLift)
{
  // Create variables for residuals
  float ResVel=0.0;
  float ResRho=0.0;
  float ResDrag=0.0;
  float ResLift=0.0;
  // Loop variables
  int i;

  *sumVel1 = 0;
  *sumRho1 = 0;

  // sum up velocity and density
  for (i=0; i<((*m)*(*n)); i++)
  {
    *sumVel1=*sumVel1+sqrt(pow(Cells_var_h[i].U,2)+pow(Cells_var_h[i].V,2));
    *sumRho1=*sumRho1+Cells_var_h[i].Rho;
    
    if (Cells_const_h[i].BoundaryID==ComputeDragLift)  // BE AWARE!! check in the STAR-CD files the ID of the surface where you want to check the drag/lift
    {        
      ResDrag+=Cells_var_h[i].DragF;
      ResLift+=Cells_var_h[i].LiftF;
    }
  }
  // Lift and drag
  Residuals[2] = ResDrag;
  Residuals[3] = ResLift;

  // Calculate residuals
  ResVel=sqrt(pow(((*sumVel0-*sumVel1)/max(*sumVel0,*sumVel1)),2));
  ResRho=sqrt(pow(((*sumRho0-*sumRho1)/max(*sumRho0,*sumRho1)),2));

  // be safe
  if(ResVel<1)
  {
    Residuals[0]=ResVel;
  }
  else
  {
    Residuals[0]=1.0;
  }

  if(ResRho<1)
  {
    Residuals[1]=ResRho;
  }
  else
  {
    Residuals[1]=1.0;
  }

  if(ResRho!=ResRho) // if density residuals are NaN
  {
    printf("\nDensity divergence!\n");
    exit(1); // ERROR!
  }
}
