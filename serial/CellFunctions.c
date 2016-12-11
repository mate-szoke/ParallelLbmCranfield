#include <stdio.h>                  // printf();
#include <stdlib.h>                 // for malloc();
#include <stdbool.h>                // for bool type variables!
#include <math.h>                   // for sin,cos,pow... compile with -lm
#include "include/ShellFunctions.h" // convenience

/*==================================================
=========Initialization for the MRT model===========
==================================================*/
// This function fills up tm and stimv with variables
void MRTInitializer(MyReal** tm, MyReal** stmiv, MyReal Omega)
{
  // RETURN THESE VALUES:


  ///////////// Declarations ////////////////
  int i, j, l;  // loop variables
  
  // declarations for this collision model
  MyReal sumcc;
  MyReal sm[9];
  MyReal ev[9][9];
  const MyReal a1=1./36.;
  MyReal tminv[9][9]=
  {
      {4.*a1, -4.*a1,  4.*a1,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0},
      {4.*a1,    -a1, -2.*a1,  6.*a1, -6.*a1,    0.0,    0.0,  9.*a1,    0.0},
      {4.*a1,    -a1, -2.*a1,    0.0,    0.0,  6.*a1, -6.*a1, -9.*a1,    0.0},
      {4.*a1,    -a1, -2.*a1, -6.*a1,  6.*a1,    0.0,    0.0,  9.*a1,    0.0},
      {4.*a1,    -a1, -2.*a1,    0.0,    0.0, -6.*a1,  6.*a1, -9.*a1,    0.0},
      {4.*a1,  2.*a1,     a1,  6.*a1,  3.*a1,  6.*a1,  3.*a1,    0.0,  9.*a1},
      {4.*a1,  2.*a1,     a1, -6.*a1, -3.*a1,  6.*a1,  3.*a1,    0.0, -9.*a1},
      {4.*a1,  2.*a1,     a1, -6.*a1, -3.*a1, -6.*a1, -3.*a1,    0.0,  9.*a1},
      {4.*a1,  2.*a1,     a1,  6.*a1,  3.*a1, -6.*a1, -3.*a1,    0.0, -9.*a1}
  };

  MyReal temp[9][9] = {
    {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0},
    {-4.,-1.,-1.,-1.,-1.,2.0,2.0,2.0,2.0},
    {4.0,-2.,-2.,-2.,-2.,1.0,1.0,1.0,1.0},
    {0.0,1.0,0.0,-1.,0.0,1.0,-1.,-1.,1.0},
    {0.0,-2.,0.0,2.0,0.0,1.0,-1.,-1.,1.0},
    {0.0,0.0,1.0,0.0,-1.,1.0,1.0,-1.,-1.},
    {0.0,0.0,-2.,0.0,2.0,1.0,1.0,-1.,-1.},
    {0.0,1.0,-1.,1.0,-1.,0.0,0.0,0.0,0.0},
    {0.0,0.0,0.0,0.0,0.0,1.0,-1.,1.0,-1.} 
  };


  ///////////// Fill up variables ////////////////
  
  // Filling up tm
  for (i = 0; i < 9; i++)
  {
    for (j = 0; j < 9; j++)
      tm[i][j] = temp[i][j];
  }

  // Filling up stimv
  sm[0] = 1.0;
  sm[1] = 1.4;
  sm[2] = 1.4;
  sm[3] = 1.0;
  sm[4] = 1.2;
  sm[5] = 1.0;
  sm[6] = 1.2;
  sm[7] = Omega;
  sm[8] = Omega;

  for(i=0;i<9;i++)
  {
    for(j=0;j<9;j++)
    {
      sumcc=0.0;
      for(l=0;l<9;l++)
      {
        sumcc=sumcc+tminv[i][l]*tm[l][j];
      }

      ev[i][j]=sumcc;
    }
  }

  for(i=0;i<9;i++)
  {
    for(j=0;j<9;j++)
    {
        stmiv[i][j]=tminv[i][j]*sm[j];
    }
  }

} // End of function


/*==================================================
============Initialization of the cells=============
==================================================*/

void CellIni(float **Nod,  float **Con,
             int* NumN,     int* NumC,
             int A, int B, // DiscreteX and DiscreteY 
             MyReal* MinInletCoordY, MyReal* MaxInletCoordY,
             MyReal* Delta, float Uavg, float Vavg,
             int InletProfile, int CollisionModel,
             int* opp, struct CellProps *Cells, float rho_ini, int* n) 
{
  ///////////// DECLARATION /////////////////
  int i, j, l, k;  // loop variables
  MyReal Qlat[9]={0.,1.,1.,1.,1.,sqrt(2),sqrt(2),sqrt(2),sqrt(2)}; // Q lattice, length of the discrete directions

  
  // FIND ID of the actual cell

  /* old version
  for (i=0; i<*NumN; i++ )
  {
      if ( ( Nod[i][0] == A ) && ( Nod[i][1] == B ) )
      {
          Cells[A][B].ID = i;
      }
  }
  */
  i = 0;
  for (i=0; i<*NumN; i++ )
  {
      if ( ( Nod[i][0] == A ) && ( Nod[i][1] == B ) )
      {
        (Cells+B*(*n)+A)->ID = i;
        //Cells[A][B].ID = i;
      }
  }

  // FIND X and Y of the actual cell
  (Cells+B*(*n)+A)->CoordX = Nod[ (Cells+B*(*n)+A)->ID ][2];
  (Cells+B*(*n)+A)->CoordY = Nod[ (Cells+B*(*n)+A)->ID ][3];

  // CHECK FLUID OR NOT
  (Cells+B*(*n)+A)->Fluid  = Nod[ (Cells+B*(*n)+A)->ID ][4];


  // REMEMBER BCCONNECTOR FILE
  //  ___________________________________________________________________________________________________
  // |    0    |   1   |       2       |     3     |          4         |         5        |     6       |
  // |    i    |   j   |  ID lattice   |   BC type |   X intersection   |   Y intersection | Boundary ID |
  // |_________|_______|_______________|___________|____________________|__________________|_____________|

  // INITIALIZE VARIABLEs
  (Cells+B*(*n)+A)->U   = 0;
  (Cells+B*(*n)+A)->V   = 0;
  (Cells+B*(*n)+A)->Rho = rho_ini;

  (Cells+B*(*n)+A)->Boundary   = 0;  // IT IS NOT BOUNDARY NODE
  (Cells+B*(*n)+A)->BoundaryID = 0;  // IT IS NOT BOUNDARY NODE

  for(i=0;i<9;i++)
  {
    (Cells+B*(*n)+A)->BC_ID[i]= 0  ; // IT IS NOT BOUNDARY LATTICE
    (Cells+B*(*n)+A)->Q[i]    = 0.5; 
  }

  //SEARCH FOR BC TYPE, BOUNDARY ID AND DISTANCES IN THE LATTICE
  for(i=0;i<*NumC;i++)
  {
      if ( ( (int)Con[i][0] == A ) && ( (int)Con[i][1] == B) )
      {
          for(j=1; j<9;j++)
          {
              if ( Con[i][2] == j )
              {
               (Cells+B*(*n)+A)->BC_ID[j]   = Con[i][3];
               (Cells+B*(*n)+A)->BoundaryID = Con[i][6];

               // find distance from the boundary
               (Cells+B*(*n)+A)->Q[j] = sqrt(pow( Con[i][4]-((Cells+B*(*n)+A)->CoordX),2 ) + pow( Con[i][5]-((Cells+B*(*n)+A)->CoordY),2) ) / ((*Delta)*Qlat[j]);
              }
          }
      }
  }

  // not in the corner
  (Cells+B*(*n)+A)->Corner=0;

  for(j=1; j<9;j++)
  {
    if ((Cells+B*(*n)+A)->BC_ID[j] != 0)
    {
      if ((Cells+B*(*n)+A)->Boundary == 0)// if Boundary condition in the node is 0 it becomes equal to the BC of the lattice direction
      {
        (Cells+B*(*n)+A)->Boundary=(Cells+B*(*n)+A)->BC_ID[j];
      }else{// if in the same node there are lattice directions with different BC (corners) the BC of the node is WALL (assuming that it's impossibe to findoutlet and inlet together)
        if (((Cells+B*(*n)+A)->Boundary) < ((Cells+B*(*n)+A)->BC_ID[j]) ) 
        {
          (Cells+B*(*n)+A)->Boundary=1;
          (Cells+B*(*n)+A)->Corner=1;
        }
        if ((Cells+B*(*n)+A)->Boundary > (Cells+B*(*n)+A)->BC_ID[j])
        {
          (Cells+B*(*n)+A)->Boundary=1;
          (Cells+B*(*n)+A)->Corner=1;
        }
      }
    }
  }

  
  //Boundary=0 NO BC
  //Boundary=1 WALL
  //Boundary=2 INLET
  //Boundary=3 OUTLET


  // BC ON CORNERS IS WALL!!!! (this operation is useful for wall condition, which checks the single direction)
  if ( ((Cells+B*(*n)+A)->Corner) == 1)
  {
    if ((Cells+B*(*n)+A)->BC_ID[1] != 0 && (Cells+B*(*n)+A)->BC_ID[2] != 0 )
    {
        (Cells+B*(*n)+A)->BC_ID[5]=1;
    }
    if ((Cells+B*(*n)+A)->BC_ID[1] != 0 && (Cells+B*(*n)+A)->BC_ID[4] != 0 )
    {
        (Cells+B*(*n)+A)->BC_ID[8]=1;
    }
    if ((Cells+B*(*n)+A)->BC_ID[2] != 0 && (Cells+B*(*n)+A)->BC_ID[3] != 0 )
    {
        (Cells+B*(*n)+A)->BC_ID[6]=1;
    }
    if ((Cells+B*(*n)+A)->BC_ID[3] != 0 && (Cells+B*(*n)+A)->BC_ID[4] != 0 )
    {
        (Cells+B*(*n)+A)->BC_ID[7]=1;
    }
  }


  // INITIALIZE STREAMING (STREAM EVERYWHERE)
  for(i=0;i<9;i++)  
  {
    (Cells+B*(*n)+A)->StreamLattice[i] = 1;
    //printf("StreamLattice[%d] = %d\n",i, Cells[A][B].StreamLattice[i]);
  }


  // DON'T STREAM FROM OUTSIDE OF THE DOMAIN
  for(k=0;k<9;k++)
  {
    if ((Cells+B*(*n)+A)->BC_ID[k]!=0)
    {
        (Cells+B*(*n)+A)->StreamLattice[opp[k]]= 0 ;
    }
  }


  // INLET VELOCITY // THIS IS CRAPPY, NOT USED!
  switch(InletProfile)
  {
    case 1:
      //Uo=1.5*Uavg*(1-4*(pow((CoordY-MinInletCoordY-0.5*(MaxInletCoordY-MinInletCoordY)),2)));
      //Uo=4*1.5*Uavg*CoordY*(41-CoordY)/(41*41);
      (Cells+B*(*n)+A)->Uo = 4*1.5*Uavg*(((Cells+B*(*n)+A)->CoordY)-(*MinInletCoordY))*(((*MaxInletCoordY)-(*MinInletCoordY))-(((Cells+B*(*n)+A)->CoordY)-(*MinInletCoordY)))/(((*MaxInletCoordY)-(*MinInletCoordY))*((*MaxInletCoordY)-(*MinInletCoordY)));
      (Cells+B*(*n)+A)->Vo = Vavg;
    break;
    case 2:
      (Cells+B*(*n)+A)->Uo = Uavg;
      (Cells+B*(*n)+A)->Vo = Vavg;
    break;

  }
  (Cells+B*(*n)+A)->U = (Cells+B*(*n)+A)->Uo;


} // End of function


  /*==================================================
  ========Creating constant lattice parameters========
  ==================================================*/

void D2Q9Vars(MyReal* w, int* cx, int* cy, int* opp, int* c, int* n)
{
  // Fill up variables with constants
  //  ID lattice
  //        6       2       5
  //          \     |     /
  //            \   |   /
  //              \ | /
  //        3 - - - 0 - - - 1
  //              / | \
  //            /   |   \
  //          /     |     \
  //        7       4      8
  
  int i;
  // D2Q9 properties
    
  w[0]=4./9.;

  for (i=1; i<5; i++ )
      w[i]=1./9.;

  for (i=5; i<9; i++ )
      w[i]=1./36.;
  
  cx[0] =  0;
  cx[1] =  1;
  cx[2] =  0;
  cx[3] = -1;
  cx[4] =  0;
  cx[5] =  1;
  cx[6] = -1;
  cx[7] = -1;
  cx[8] =  1;

  cy[0] =  0;
  cy[1] =  0;
  cy[2] =  1;
  cy[3] =  0;
  cy[4] = -1;
  cy[5] =  1;
  cy[6] =  1;
  cy[7] = -1;
  cy[8] = -1;

  opp[0] = 0;
  opp[1] = 3;
  opp[2] = 4;
  opp[3] = 1;
  opp[4] = 2;
  opp[5] = 7;
  opp[6] = 8;
  opp[7] = 5;
  opp[8] = 6;

  c[0] =         0;
  c[1] =        -1;
  c[2] = -1*(*n)  ;
  c[3] =         1;
  c[4] =    (*n)  ;
  c[5] = -1*(*n)-1;
  c[6] = -1*(*n)+1;
  c[7] =    (*n)+1;
  c[8] =    (*n)-1;

} // End of function



  /*==================================================
  =============Function for the BGKW model============
  ==================================================*/



void BGKW(struct CellProps *Cells, int i, MyReal* w, int* cx, int* cy, MyReal Omega)
{

  MyReal T1 = ((Cells+i)->U) * ((Cells+i)->U) + 
              ((Cells+i)->V) * ((Cells+i)->V);
  MyReal T2 = 0.0;
  int k;
  for (k=0; k<9; k++)
  {
    T2                  = ((Cells+i)->U)*cx[k] + ((Cells+i)->V)*cy[k];
    (Cells+i)->Feq[k]   = ((Cells+i)->Rho)*w[k]*( 1.0 + 3.0*T2 + 4.5*T2*T2 - 1.5*T1 );
    (Cells+i)->METAF[k] = Omega*((Cells+i)->Feq[k])+(1.0-Omega)*((Cells+i)->F[k]);
  }
}

  /*==================================================
  =============Function for the TRT model=============
  ==================================================*/


void TRT(struct CellProps *Cells, int i, MyReal* w, int* cx, int* cy, int* opp, MyReal Omega, MyReal OmegaA)
{
  MyReal T1 = ((Cells +i)->U)*((Cells +i)->U)+((Cells +i)->V)*((Cells +i)->V);
  MyReal T2 = 0.0;
  int k;
  for (k=0; k<9; k++)
  {
      T2                  = ((Cells +i)->U)   *cx[k] + ((Cells +i)->V)*cy[k];
      (Cells +i)->Feq[k]  = ((Cells +i)->Rho) *w[k]*(1.0+3.0*T2+4.5*T2*T2-1.5*T1);
  }
  
  float F_p[9],Feq_p[9],F_m[9],Feq_m[9];
  for (k=0; k<9; k++)
  {
      F_p[k]   = 0.5*( (Cells +i)->F[k]   + (Cells +i)->F[opp[k]]   );
      Feq_p[k] = 0.5*( (Cells +i)->Feq[k] + (Cells +i)->Feq[opp[k]] );
      F_m[k]   = 0.5*( (Cells +i)->F[k]   - (Cells +i)->F[opp[k]]   );
      Feq_m[k] = 0.5*( (Cells +i)->Feq[k] - (Cells +i)->Feq[opp[k]] );
      
      (Cells +i)->METAF[k] = (Cells +i)->F[k] - 
                             (F_p[k]-Feq_p[k])*Omega - 
                             (F_m[k]-Feq_m[k])*OmegaA;
  }
}

  /*==================================================
  =============Function for the MRT model=============
  ==================================================*/

void MRT(struct CellProps *Cells, int i, MyReal** tm, MyReal** stmiv)
{
  int k, l;
  MyReal fmom[9],fmeq[9];
  MyReal U, V;
  MyReal suma,sumb;

  U = (Cells +i)->U;
  V = (Cells +i)->V;

  fmeq[0] =  (Cells +i)->Rho;
  fmeq[1] = ((Cells +i)->Rho)*(-2.0+3.0*((Cells +i)->Rho)*(U*U+V*V));
  fmeq[2] = ((Cells +i)->Rho)*(1.0-3.0*((Cells +i)->Rho)*(U*U+V*V));
  fmeq[3] = ((Cells +i)->Rho)*((Cells +i)->U);
  fmeq[4] =-((Cells +i)->Rho)*((Cells +i)->U);
  fmeq[5] = ((Cells +i)->Rho)*((Cells +i)->V);
  fmeq[6] =-((Cells +i)->Rho)*((Cells +i)->V);
  fmeq[7] = ((Cells +i)->Rho)*(U*U-V*V);
  fmeq[8] = ((Cells +i)->Rho)*U*V;
  
  for (k=0; k<9;k++)
  {
    suma=0.0;
    for (l=0; l<9;l++)
      suma = suma + tm[k][l]*((Cells +i)->F[l]);

    fmom[k]=suma;
  }

  for (k=0; k<9;k++)
  {
    sumb=0.0;
    for (l=0; l<9;l++)
      sumb = sumb + stmiv[k][l]*(fmom[l]-fmeq[l]);
    
    (Cells +i)->METAF[k] = ((Cells +i)->F[k]) - sumb;
  }
}


  /*==================================================
  ======Function to update the distribution fct.======
  ==================================================*/

void UpdateF(struct CellProps *Cells, int i)
{
  int k;
  for(k=0;k<9;k++)
    {(Cells+i)->F[k] = (Cells+i)->METAF[k];}
}


  /*==================================================
  ======Function to update the macroscopic var.=======
  ==================================================*/

void UpdateMacroscopic(struct CellProps *Cells, int j, int i, int* cx, int* cy, int CalculateDragLift)
{
  MyReal Ssum, Usum, Vsum;
  int k;

  if ((Cells+i)->Fluid==1)
  {
    Ssum=0.0;
    for (k=0; k<9; k++)
        Ssum = Ssum+(Cells+i)->F[k];

    (Cells+i)->Rho = Ssum;

    //printf("Rho[%d][%d] = %f\n", j, i, Ssum);

    Usum = 0.0;
    Vsum = 0.0;
    for (k=0; k<9; k++)
    {
      Usum = Usum + (Cells+i)->F[k]*cx[k];
      Vsum = Vsum + (Cells+i)->F[k]*cy[k];
    }
    (Cells+i)->U = Usum/(Cells+i)->Rho;
    (Cells+i)->V = Vsum/(Cells+i)->Rho;
  }

  if ((Cells+i)->BC_ID[1]==3) // for outlet on the right
  {
      (Cells+i)->V=0.0;
  }

  //   DRAG/LIFT FORCE
  if (CalculateDragLift != 0 && (Cells+i)->BoundaryID==CalculateDragLift)
  {
     (Cells+i)->DragF = (Cells+i)->Rho/3*(20-(Cells+i)->CoordX)/5;
     (Cells+i)->LiftF = (Cells+i)->Rho/3*(20-(Cells+i)->CoordY)/5;
  }

}

void CalculateDragLiftForces(struct CellProps *Cells, int j, int i, int CalculateDragLift)
{
  //   DRAG/LIFT FORCE
  if (CalculateDragLift != 0 && (Cells+i)->BoundaryID==CalculateDragLift)
  {
     (Cells+i)->DragF = (Cells+i)->Rho/3*(20-(Cells+i)->CoordX)/5;
     (Cells+i)->LiftF = (Cells+i)->Rho/3*(20-(Cells+i)->CoordY)/5;
  }
}
