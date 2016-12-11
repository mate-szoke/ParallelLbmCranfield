#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>

#include <upc.h>                 // Required for UPC 
//#include <upc_cray.h>                 // Required for UPC 
#include "include/ShellFunctions.h"

void ComputeResiduals(double* Residuals, double* sumVel0, double* sumVel1, double* sumRho0, double* sumRho1, int ComputeDragLift, int lm, int ln)
{

  *sumVel0=*sumVel1;
  *sumRho0=*sumRho1;

  // Create variables for residuals
  double ResVel   = 0.0;
  double ResRho   = 0.0;
  double ResDrag  = 0.0;
  double ResLift  = 0.0;
  double LsumVel1 = 0.0;
  double LsumRho1 = 0.0;
  // Loop variables
  int i, j;

  *sumVel1 = 0;
  *sumRho1 = 0;

  // sum up velocity and density
  for (i=0; i<(lm)*(ln);i++)
  {
    //for(i=0; i<*n; i++)
    //{
      //*sumVel1=*sumVel1+sqrt(pow(U[i],2)+pow(V[i],2));
      //*sumRho1=*sumRho1+Rho[i];
      
      LsumVel1 = LsumVel1 + sqrt( (U[i])*(U[i])  +  (V[i])*(V[i])  );
      LsumRho1 = LsumRho1 + Rho[i];

      if (BoundaryID[i]==ComputeDragLift)  // BE AWARE!! check in the STAR-CD files the ID of the surface where you want to check the drag/lift. In ths case it's 1
      {        
        ResDrag+=DragF[i];
        ResLift+=LiftF[i];
      }
    //}
  }

  *sumVel1 = LsumVel1;
  *sumRho1 = LsumRho1;
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
    upc_global_exit(1);
    puts("you shall not see this!");
    //exit(1); // ERROR!
  }
}

