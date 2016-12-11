#include <stdio.h> 
#include <stdlib.h> 
#include "math.h"
#include "include/ShellFunctions.h"

void ComputeResiduals(int* BoundaryID, float* F, float* METAF, float* DragF, float* LiftF, float* Residuals, int* m, int* n, int ComputeDragLift)
{
  // Create variables for residuals
  float ResDrag=0.0;
  float ResLift=0.0;
  // Loop variables
  // Loop variables
  int i, k;

  float L2n   = 0.0;  // L2 norm
  float L2n_w = 0.0;  // L2 norm weighted

  // sum up velocity and density
  for (i=0; i<((*m)*(*n)); i++)
  {
    for (k=0; k<9; k++)
    { L2n = L2n + pow((F[i+k*((*m)*(*n))]-METAF[i+k*((*m)*(*n))]),2); }
    
    if (BoundaryID[i]==ComputeDragLift)  // BE AWARE!! check in the STAR-CD files the ID of the surface where you want to check the drag/lift
    {        
      ResDrag+=DragF[i];
      ResLift+=LiftF[i];
    }
  }



  // Calculate residuals
  L2n_w = sqrt(L2n/((*m)*(*n)));
  L2n = sqrt(L2n);
  // Write them to vector
  Residuals[0] = L2n;
  Residuals[1] = L2n_w;

  // Lift and drag
  Residuals[2] = ResDrag;
  Residuals[3] = ResLift;


  if(L2n!=L2n) // if density residuals are NaN
  {
    printf("\nDIVERGENCE!\n");
    exit(1); // ERROR!
  }
}
