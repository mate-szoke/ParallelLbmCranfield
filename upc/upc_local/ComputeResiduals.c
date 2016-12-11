
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include <upc_relaxed.h>                 // Required for UPC 
#include <time.h>

#include "include/ShellFunctions.h"


void ComputeResiduals(struct CellProps *Cells, double* Residuals,
                    double* sumVel0, double* sumVel1, double* sumRho0,
                    double* sumRho1, int ComputeDragLift, int* iter, int* Iterations)
{
  
  // Update variables
  *sumVel0=*sumVel1;
  *sumRho0=*sumRho1;
  
  // Create variables for residuals
  //double ResVel   = 0.0;
  //double ResRho   = 0.0;
  double ResDrag  = 0.0;
  double ResLift  = 0.0;
  //double LsumVel1 = 0.0;
  //double LsumRho1 = 0.0;
  double PUTtmp[4];
  double GETtmp[4*THREADS];
  clock_t tInstant1, tInstant2; // Time measurement points: universal

  double L2n   = 0.0;  // L2 norm
  double L2n_w = 0.0;  // L2 norm weighted

  // Loop variables
  int i, j, k;

  *sumVel1 = 0;
  //*sumRho1 = 0;
  
  // sum up velocity and density
  for(i = (*n);  i < (MLIM+1)*(*n);  i++)
  {

    //LsumVel1 = LsumVel1 + sqrt( ((Cells+i)->U)*((Cells+i)->U)  +  ((Cells+i)->V)*((Cells+i)->V)  );
    //LsumRho1 = LsumRho1 + (Cells+i)->Rho;

      for (k=0; k<9;k++)
      { L2n = L2n + pow(((Cells+i)->F[k]-(Cells+i)->METAF[k]),2); }


    if ((Cells+i)->BoundaryID == ComputeDragLift)  // BE AWARE!! check in the STAR-CD files the ID of the surface where you want to check the drag/lift. 
    {        
      ResDrag += (Cells+i)->DragF;
      ResLift += (Cells+i)->LiftF;
    }
  }

  PUTtmp[0] = L2n;
  PUTtmp[1] = L2n;
  PUTtmp[2] = ResDrag;
  PUTtmp[3] = ResLift;
  
  *sumVel1 = 0;
  //*sumRho1 = 0;
  ResDrag  = 0;
  ResLift  = 0;

  upc_memput( &sResiduals[0+4*MYTHREAD] , &PUTtmp[0],  4*sizeof(double) );
//  tInstant1 = clock(); // Start measuring time
  upc_barrier;
//  tInstant2 = clock(); // End of time measuremet
  upc_memget( &GETtmp[0] , &sResiduals[0] , (THREADS*4)*sizeof(double) );
//  *tSendRcv =  (tInstant2-tInstant1) ;

  //printf("iter %d TH%d Res Synced\n", iter, MYTHREAD);
  
  if(MYTHREAD==0)
  {
    //printf("Th%d residualcalc\n",MYTHREAD);
    for (i = 0; i < THREADS; i++)
    {
      
      *sumVel1 = *sumVel1 + GETtmp[0+4*i];
      //*sumRho1 = *sumRho1 + GETtmp[1+4*i];
      ResDrag  = ResDrag  + GETtmp[2+4*i];
      ResLift  = ResLift  + GETtmp[3+4*i];
      
    }

  L2n = *sumVel1;
  // Calculate residuals
  L2n_w = sqrt(L2n/((*m)*(*n)));
  L2n = sqrt(L2n);

  Residuals[0] = L2n;
  Residuals[1] = L2n_w;
  Residuals[2] = ResDrag;
  Residuals[3] = ResLift;
 
  }

  if(L2n!=L2n) // if density residuals are NaN
  {
    printf("\nDIVERGENCE!\n");
    exit(1); // ERROR!
    *iter  = *Iterations+1;
    Residuals[0] = 1;
    Residuals[1] = 1;
    Residuals[2] = 0;
    Residuals[3] = 0;
    
    //upc_global_exit(1);
    //exit(1); // ERROR!
  }
}

/*

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include <upc_relaxed.h>                 // Required for UPC 

#include "include/ShellFunctions.h"

void ComputeResiduals(struct CellProps *Cells, double* Residuals, double* sumVel0, double* sumVel1, double* sumRho0, double* sumRho1, int ComputeDragLift, int* iter, int* Iterations)
{
  // Update variables
  *sumVel0=*sumVel1;
  *sumRho0=*sumRho1;
  
  // Create variables for residuals
  double ResVel  = 0.0;
  double ResRho  = 0.0;
  double ResDrag = 0.0;
  double ResLift = 0.0;

  // Loop variables
  int i, j;

  *sumVel1 = 0;
  *sumRho1 = 0;
  
  // sum up velocity and density
  
  for (j = 1; j<MLIM+1; j++)
  {
    for(i=0; i<*n; i++)
    {
      *sumVel1 = *sumVel1 + sqrt( pow( (Cells+j*(*n)+i)->U ,2) + pow( (Cells+j*(*n)+i)->V ,2) );
      *sumRho1 = *sumRho1 + (Cells+j*(*n)+i)->Rho;
      
      if ((Cells+j*(*n)+i)->BoundaryID == ComputeDragLift)  // BE AWARE!! check in the STAR-CD files the ID of the surface where you want to check the drag/lift. 
      {        
        ResDrag += (Cells+j*(*n)+i)->DragF;
        ResLift += (Cells+j*(*n)+i)->LiftF;
      }

    }
  }

  sResiduals[0+4*MYTHREAD] = *sumVel1;
  sResiduals[1+4*MYTHREAD] = *sumRho1;
  sResiduals[2+4*MYTHREAD] = ResDrag;
  sResiduals[3+4*MYTHREAD] = ResLift;

  *sumVel1 = 0;
  *sumRho1 = 0;
  ResDrag  = 0;
  ResLift  = 0;

  upc_barrier;

  //printf("iter %d TH%d Res Synced\n", iter, MYTHREAD);

  if(MYTHREAD==0)
  {
    //printf("Th%d residualcalc\n",MYTHREAD);
    for (i = 0; i < THREADS; i++)
    {
      *sumVel1 = *sumVel1 + sResiduals[0+4*i];
      *sumRho1 = *sumRho1 + sResiduals[1+4*i];
      ResDrag  = ResDrag  + sResiduals[2+4*i];
      ResLift  = ResLift  + sResiduals[3+4*i];
    }

  // Calculate residuals
  ResVel = sqrt( pow( ((*sumVel0-*sumVel1)/max(*sumVel0,*sumVel1)) ,2) );
  ResRho = sqrt( pow( ((*sumRho0-*sumRho1)/max(*sumRho0,*sumRho1)) ,2) );

  Residuals[0] = ResVel;
  Residuals[1] = ResRho;
  Residuals[2] = ResDrag;
  Residuals[3] = ResLift;
  
  }


  if(ResRho!=ResRho) // if density residuals are NaN
  {
    if(MYTHREAD==0){
      printf("ResRho = %f\n", ResRho);
      printf("\nDensity divergence!\n");
    }
    *iter  = *Iterations+1;
    Residuals[0] = 1;
    Residuals[1] = 1;
    Residuals[2] = 0;
    Residuals[3] = 0;
    
    //upc_global_exit(1);
    //exit(1); // ERROR!
  }
}
*/
