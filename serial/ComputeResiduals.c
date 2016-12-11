#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include "include/ShellFunctions.h"

void ComputeResiduals(struct CellProps *Cells, MyReal* Residuals, int* m, int* n, int ComputeDragLift)
{

  // Create variables for residuals
  MyReal ResDrag = 0.0;
  MyReal ResLift = 0.0;
  // Loop variables
  int i, j, k;

  MyReal L2n   = 0.0;  // L2 norm
  MyReal L2n_w = 0.0;  // L2 norm weighted
  
  // sum
  for (i=0; i<*m;i++)
  {
    for (j=0; j<*n;j++)
    {

      for (k=0; k<9;k++)
      { L2n = L2n + pow(((Cells+i)->F[k]-(Cells+i)->METAF[k]),2); }

      if ((Cells+i)->BoundaryID==ComputeDragLift)  // BE AWARE!! check in the STAR-CD files the ID of the surface where you want to check the drag/lift. In ths case it's 1
      {        
        ResDrag+=(Cells+i)->DragF;
        ResLift+=(Cells+i)->LiftF;
      }
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