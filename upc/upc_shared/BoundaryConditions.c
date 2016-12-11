#include <stdio.h>                  // printf();
#include <stdlib.h>                 // for malloc();
#include <stdbool.h>                // for bool type variables!
#include <math.h>                   // for sin,cos,pow... compile with -lm
//#include <upc_cray.h>            // Required for UPC 

#include "include/ShellFunctions.h" // convenience


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  ////////////////////// Inlet boundary treatment /////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

void InletBC(int j, int i, int lm, int ln)
{
  double RhoN = 0.0, RhoW = 0.0;
  if (Boundary[j*(ln)+i] == 2) // if inlet boundary
  {
    if (BC_ID[1][j*(ln)+i]==2 && Corner[j*(ln)+i]!=1)
    {
      // corner treatment has to be worked out 
    }

    // inlet on the top (to be completed with the y velocity)
    if (BC_ID[2][j*(ln)+i]==2 && Corner[j*(ln)+i]!=1)
    {

      RhoN             = F[0][j*(ln)+i] + F[1][j*(ln)+i] + F[3][j*(ln)+i]+
                      2*(F[2][j*(ln)+i] + F[6][j*(ln)+i] + F[5][j*(ln)+i]);

      F[4][j*(ln)+i] = F[2][j*(ln)+i];

      F[8][j*(ln)+i] = F[6][j*(ln)+i]+RhoN*(Uo[j*(ln)+i])/6.0;

      F[7][j*(ln)+i] = F[5][j*(ln)+i]-RhoN*(Uo[j*(ln)+i])/6.0;
    }

    // inlet on the left (to be completed with the x velocity)
    if (BC_ID[3][j*(ln)+i]==2 && Corner[j*(ln)+i]!=1)
    {
        RhoW             = (F[0][j*(ln)+i]+F[2][j*(ln)+i]+F[4][j*(ln)+i]+
                           2.0*(F[3][j*(ln)+i]+F[6][j*(ln)+i]+
                           F[7][j*(ln)+i]))/(1.0-Uo[j*(ln)+i]);

        F[1][j*(ln)+i] = F[3][j*(ln)+i]+2*RhoW*(Uo[j*(ln)+i])/3;

        F[5][j*(ln)+i] = F[7][j*(ln)+i]+RhoW*(Uo[j*(ln)+i])/6;

        F[8][j*(ln)+i] = F[6][j*(ln)+i]+RhoW*(Uo[j*(ln)+i])/6;
    }

    if (BC_ID[4][j*(ln)+i]==2 && Corner[j*(ln)+i]!=1)
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

void OutletBoundaries(int j, int i, int lm, int ln)
{
  double RhoE;
   if (BC_ID[1][j*(ln)+i]==3) // outlet boundary on the right side of the domain
  {

      RhoE = (F[0][j*(ln)+i]+F[2][j*(ln)+i]+F[4][j*(ln)+i]
           +2.0*(F[1][j*(ln)+i]+F[5][j*(ln)+i]+F[8][j*(ln)+i]))/(1-Uo[j*(ln)+i]);

      F[3][j*(ln)+i] = F[1][j*(ln)+i]-2*RhoE*(Uo[j*(ln)+i])/3.0;

      F[7][j*(ln)+i] = F[5][j*(ln)+i]-RhoE*(Uo[j*(ln)+i])/6.0;

      F[6][j*(ln)+i] = F[8][j*(ln)+i]-RhoE*(Uo[j*(ln)+i])/6.0;
  }
  if (BC_ID[2][j*(ln)+i]==3)
  {
      // FILL!!
  }
  if (BC_ID[3][j*(ln)+i]==3)
  {
      // FILL!!
  }
  if (BC_ID[4][j*(ln)+i]==3)
  {
      // FILL!!
  }

}


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////// Wall boundary treatment /////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

void WallBC(int j, int i, int* opp, int lm, int ln)
{

// WALL BC (half-way bounceback)
  int k;
  for(k=0;k<9;k++)
  {
    if (BC_ID[k][j*(ln)+i]==1) // if wall boundary
    {
      F[opp[k]][j*(ln)+i] = F[k][j*(ln)+i];
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

void CurvedWallBoundaries(int j, int i, int* opp, int lm, int ln)
{
  int k=0;
  for(k=0;k<9;k++)
  {
    if (BC_ID[k][j*(ln)+i]==1) // if wall
    {
      if (Q[k][j*(ln)+i]<0.5) // if the distance from the boundary is less than 0.5?
      { 
        F[opp[k]][j*(ln)+i] = 2*Q[k][j*(ln)+i] * METAF[k][j*(ln)+i] + (1-2*Q[k][j*(ln)+i])*Fneighbours[k][j*(ln)+i];
      }
      else
      {
        F[opp[k]][j*(ln)+i] = METAF[k][j*(ln)+i]/2/Q[k][j*(ln)+i]
                              +(2*Q[k][j*(ln)+i]-1)/(2*Q[k][j*(ln)+i])*METAF[opp[k]][j*(ln)+i];
      }
    }
  }
}
