#include <stdio.h>                  // printf();
#include <stdlib.h>                 // for malloc();
#include <stdbool.h>                // for bool type variables!
#include <math.h>                   // for sin,cos,pow... compile with -lm
#include <upc_relaxed.h>            // Required for UPC 

#include "include/ShellFunctions.h" // convenience


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  ////////////////////// Inlet boundary treatment /////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

void InletBC(struct CellProps *Cells, int j, int i)
{
  double RhoN = 0.0, RhoW = 0.0;
  if ((Cells+j*(*n)+i)->Boundary == 2) // if inlet boundary
  {
    if ((Cells+j*(*n)+i)->BC_ID[1]==2 && (Cells+j*(*n)+i)->Corner!=1)
    {
      // corner treatment has to be worked out 
    }

    // inlet on the top (to be completed with the y velocity)
    if ((Cells+j*(*n)+i)->BC_ID[2]==2 && (Cells+j*(*n)+i)->Corner!=1)
    {

      RhoN             = (Cells+j*(*n)+i)->F[0] + (Cells+j*(*n)+i)->F[1] + (Cells+j*(*n)+i)->F[3]+
                      2*((Cells+j*(*n)+i)->F[2] + (Cells+j*(*n)+i)->F[6] + (Cells+j*(*n)+i)->F[5]);

      (Cells+j*(*n)+i)->F[4] = (Cells+j*(*n)+i)->F[2];

      (Cells+j*(*n)+i)->F[8] = (Cells+j*(*n)+i)->F[6]+RhoN*((Cells+j*(*n)+i)->Uo)/6.0;

      (Cells+j*(*n)+i)->F[7] = (Cells+j*(*n)+i)->F[5]-RhoN*((Cells+j*(*n)+i)->Uo)/6.0;
    }

    // inlet on the left (to be completed with the x velocity)
    if ((Cells+j*(*n)+i)->BC_ID[3]==2 && (Cells+j*(*n)+i)->Corner!=1)
    {
        RhoW             = ((Cells+j*(*n)+i)->F[0]+(Cells+j*(*n)+i)->F[2]+(Cells+j*(*n)+i)->F[4]+
                           2.0*((Cells+j*(*n)+i)->F[3]+(Cells+j*(*n)+i)->F[6]+
                           (Cells+j*(*n)+i)->F[7]))/(1.0-(Cells+j*(*n)+i)->Uo);

        (Cells+j*(*n)+i)->F[1] = (Cells+j*(*n)+i)->F[3]+2*RhoW*((Cells+j*(*n)+i)->Uo)/3;

        (Cells+j*(*n)+i)->F[5] = (Cells+j*(*n)+i)->F[7]+RhoW*((Cells+j*(*n)+i)->Uo)/6;

        (Cells+j*(*n)+i)->F[8] = (Cells+j*(*n)+i)->F[6]+RhoW*((Cells+j*(*n)+i)->Uo)/6;
    }

    if ((Cells+j*(*n)+i)->BC_ID[4]==2 && (Cells+j*(*n)+i)->Corner!=1)
    {
      // corner treatment has to be worked out 
    }
  }

}


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  ////////////////////// Outlet boundary treatment ////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

void OutletBoundaries(struct CellProps *Cells, int j, int i)
{
  double RhoE;
   if ((Cells+j*(*n)+i)->BC_ID[1]==3) // outlet boundary on the right side of the domain
  {

      RhoE = ((Cells+j*(*n)+i)->F[0]+(Cells+j*(*n)+i)->F[2]+(Cells+j*(*n)+i)->F[4]
           +2.0*((Cells+j*(*n)+i)->F[1]+(Cells+j*(*n)+i)->F[5]+(Cells+j*(*n)+i)->F[8]))/(1-(Cells+j*(*n)+i)->Uo);

      (Cells+j*(*n)+i)->F[3] = (Cells+j*(*n)+i)->F[1]-2*RhoE*((Cells+j*(*n)+i)->Uo)/3.0;

      (Cells+j*(*n)+i)->F[7] = (Cells+j*(*n)+i)->F[5]-RhoE*((Cells+j*(*n)+i)->Uo)/6.0;

      (Cells+j*(*n)+i)->F[6] = (Cells+j*(*n)+i)->F[8]-RhoE*((Cells+j*(*n)+i)->Uo)/6.0;
  }
  if ((Cells+j*(*n)+i)->BC_ID[2]==3)
  {
      // FILL!!
  }
  if ((Cells+j*(*n)+i)->BC_ID[3]==3)
  {
      // FILL!!
  }
  if ((Cells+j*(*n)+i)->BC_ID[4]==3)
  {
      // FILL!!
  }

}


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////// Wall boundary treatment /////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

void WallBC(struct CellProps *Cells, int j, int i, int* opp)
{

// WALL BC (half-way bounceback)
  int k;
  for(k=0;k<9;k++)
  {
    if ((Cells+j*(*n)+i)->BC_ID[k]==1) // if wall boundary
    {
      (Cells+j*(*n)+i)->F[opp[k]] = (Cells+j*(*n)+i)->F[k];
    }
  }

//WALL BC (full-way bounceback)
// if (Boundary==1){
//       for(k=0;k<9;k++){
//                F[opp[k]]=F[k];
//        }
//     }

}


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////// Curved wall boundary treatment //////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

void CurvedWallBoundaries(struct CellProps *Cells, int j, int i, int* opp)
{
  int k=0;
  for(k=0;k<9;k++)
  {
    if ((Cells+j*(*n)+i)->BC_ID[k]==1) // if wall
    {
      if ((Cells+j*(*n)+i)->Q[k]<0.5) // if the distance from the boundary is less than 0.5?
      { 
        (Cells+j*(*n)+i)->F[opp[k]] = 2*(Cells+j*(*n)+i)->Q[k]*(Cells+j*(*n)+i)->METAF[k]
                              +(1-2*(Cells+j*(*n)+i)->Q[k])*(Cells+j*(*n)+i)->Fneighbours[k];
      }
      else
      {
        (Cells+j*(*n)+i)->F[opp[k]] = (Cells+j*(*n)+i)->METAF[k]/2/(Cells+j*(*n)+i)->Q[k]
                              +(2*(Cells+j*(*n)+i)->Q[k]-1)/(2*(Cells+j*(*n)+i)->Q[k])*(Cells+j*(*n)+i)->METAF[opp[k]];
      }
    }
  }
}
