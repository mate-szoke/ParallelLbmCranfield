#include <stdio.h>                  // printf();
#include <stdlib.h>                 // for malloc();
#include <stdbool.h>                // for bool type variables!
#include <math.h>                   // for sin,cos,pow... compile with -lm
//#include <upc_cray.h>            // Required for UPC 

#include "include/ShellFunctions.h" // convenience

/*==================================================
=========Initialization for the MRT model===========
==================================================*/
// This function fills up tm and stimv with variables
void MRTInitializer(float** tm, float** stmiv, double Omega)
{
  // RETURN THESE VALUES:
  // float tm[9][9];
  // float stmiv[9][9];

  ///////////// Declarations ////////////////
  int i, j, l;  // loop variables
  
  // declarations for this collision model
  float sumcc;
  float sm[9];
  float ev[9][9];
  /*
  const float a1=1./36.;
  float tminv[9][9]=
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
  };*/

  float tminv[9][9] =
  {
      {4.*(1./36), -4.*(1./36),  4.*(1./36),         0.0,         0.0,         0.0,         0.0,         0.0,         0.0},
      {4.*(1./36),    -(1./36), -2.*(1./36),  6.*(1./36), -6.*(1./36),         0.0,         0.0,  9.*(1./36),         0.0},
      {4.*(1./36),    -(1./36), -2.*(1./36),         0.0,         0.0,  6.*(1./36), -6.*(1./36), -9.*(1./36),         0.0},
      {4.*(1./36),    -(1./36), -2.*(1./36), -6.*(1./36),  6.*(1./36),         0.0,         0.0,  9.*(1./36),         0.0},
      {4.*(1./36),    -(1./36), -2.*(1./36),         0.0,         0.0, -6.*(1./36),  6.*(1./36), -9.*(1./36),         0.0},
      {4.*(1./36),  2.*(1./36),     (1./36),  6.*(1./36),  3.*(1./36),  6.*(1./36),  3.*(1./36),         0.0,  9.*(1./36)},
      {4.*(1./36),  2.*(1./36),     (1./36), -6.*(1./36), -3.*(1./36),  6.*(1./36),  3.*(1./36),         0.0, -9.*(1./36)},
      {4.*(1./36),  2.*(1./36),     (1./36), -6.*(1./36), -3.*(1./36), -6.*(1./36), -3.*(1./36),         0.0,  9.*(1./36)},
      {4.*(1./36),  2.*(1./36),     (1./36),  6.*(1./36),  3.*(1./36), -6.*(1./36), -3.*(1./36),         0.0, -9.*(1./36)}
  };      



  float temp[9][9] = {
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

void CellIni(float **Nod,             
             float **Con,
             int A,                   
             int B, // DiscreteX and DiscreteY 
             float Uavg,
             float Vavg,              
             int InletProfile,
             int CollisionModel,      
             int* opp,
             float rho_ini,
             int lm,
             int ln) 
{
  ///////////// DECLARATION /////////////////
  int i, j, l, k;  // loop variables
  float Qlat[9]={0.,1.,1.,1.,1.,sqrt(2),sqrt(2),sqrt(2),sqrt(2)}; // Q lattice, length of the discrete directions

  
  // FIND ID of the actual cell
  i = 0;
  while(i<*NumNodes)
  //for (i=0; i<*NumNodes; i++ )
  {
      if ( ( (int)Nod[i][0] == A ) && ( (int)Nod[i][1] == B ) )
      {
        ID[B*(ln)+A] = i;  //////////////////////////////////////// ????????????????????????????
        i = 2*(*NumNodes);
        //Cells[A][B].ID = i;
      }
      i++;
  }

  // FIND X and Y of the actual cell
  CoordX[B*(ln)+A] = Nod[ ID[B*(ln)+A] ][2];
  CoordY[B*(ln)+A] = Nod[ ID[B*(ln)+A] ][3];
  
  // CHECK FLUID OR NOT
  Fluid[B*(ln)+A]  = Nod[ ID[B*(ln)+A] ][4];

  // Which thread does it belongs to
  ThreadNumber[B*(ln)+A] = MYTHREAD;
  
  //if( MYTHREAD != 0) {
  //printf("TH%d; A = %03d; B = %03d; ID = %04d; CoordX = %03.1f; CoordY = %03.1f; Density = %03.1f; U = %03.1f; \n",
  //  ThreadNumber, A, B, ID, CoordX, CoordY, Rho,  U);
  //}
  //printf(" B*(*m)+A = %d \n", (B*(*m)+A) );
  
  // REMEMBER BCCONNECTOR FILE
  //  ___________________________________________________________________________________________________
  // |    0    |   1   |       2       |     3     |          4         |         5        |     6       |
  // |    i    |   j   |  ID lattice   |   BC type |   X intersection   |   Y intersection | Boundary ID |
  // |_________|_______|_______________|___________|____________________|__________________|_____________|


  // INITIALIZE VARIABLEs
  U[B*(ln)+A]   = 0;
  V[B*(ln)+A]   = 0;
  Rho[B*(ln)+A] = rho_ini;

  Boundary[B*(ln)+A]   = 0;  // IT IS NOT BOUNDARY NODE
  BoundaryID[B*(ln)+A] = 0;  // IT IS NOT BOUNDARY NODE


  for(i=0;i<9;i++)
  {
    BC_ID[i][B*(ln)+A]= 0  ; // IT IS NOT BOUNDARY LATTICE
    Q[i][B*(ln)+A]    = 0.5; 
  }


  //SEARCH FOR BC TYPE, BOUNDARY ID AND DISTANCES IN THE LATTICE
  for(i=0;i<*NumConn;i++)
  {
      if ( ( (int)Con[i][0] == A ) && ( (int)Con[i][1] == B ) )
      {
          for(j=1; j<9;j++)
          {
              if ( Con[i][2] == j )
              {
               BC_ID[j][B*(ln)+A]   = Con[i][3];
               BoundaryID[B*(ln)+A] = Con[i][6];

               // find distance from the boundary
               Q[j][B*(ln)+A] = sqrt(pow( Con[i][4]-(CoordX[B*(ln)+A]),2 ) + pow( Con[i][5]-(CoordY[B*(ln)+A]),2) ) / ((*Delta)*Qlat[j]);
              }
          }
      }
  }

  // not in the corner
  Corner[B*(ln)+A]=0;

  for(j=1; j<9;j++)
  {
    if (BC_ID[j][B*(ln)+A] != 0)
    {
      if (Boundary[B*(ln)+A] == 0)// if Boundary condition in the node is 0 it becomes equal to the BC of the lattice direction
      {
        Boundary[B*(ln)+A]=BC_ID[j][B*(ln)+A];
      }else{// if in the same node there are lattice directions with different BC (corners) the BC of the node is WALL (assuming that it's impossibe to findoutlet and inlet together)
        if ((Boundary[B*(ln)+A]) < (BC_ID[j][B*(ln)+A]) ) 
        {
          Boundary[B*(ln)+A]=1;
          Corner[B*(ln)+A]=1;
        }
        if (Boundary[B*(ln)+A] > BC_ID[j][B*(ln)+A])
        {
          Boundary[B*(ln)+A]=1;
          Corner[B*(ln)+A]=1;
        }
      }
    }
  }

  
  //Boundary=0 NO BC
  //Boundary=1 WALL
  //Boundary=2 INLET
  //Boundary=3 OUTLET


  // BC ON CORNERS IS WALL!!!! (this operation is useful for wall condition, which checks the single direction)
  if ( (Corner[B*(ln)+A]) == 1)
  {
    if (BC_ID[1][B*(ln)+A] != 0 && BC_ID[2][B*(ln)+A] != 0 )
    {
        BC_ID[5][B*(ln)+A]=1;
    }
    if (BC_ID[1][B*(ln)+A] != 0 && BC_ID[4][B*(ln)+A] != 0 )
    {
        BC_ID[8][B*(ln)+A]=1;
    }
    if (BC_ID[2][B*(ln)+A] != 0 && BC_ID[3][B*(ln)+A] != 0 )
    {
        BC_ID[6][B*(ln)+A]=1;
    }
    if (BC_ID[3][B*(ln)+A] != 0 && BC_ID[4][B*(ln)+A] != 0 )
    {
        BC_ID[7][B*(ln)+A]=1;
    }
  }


  // INITIALIZE STREAMING (STREAM EVERYWHERE)
  for(i=0;i<9;i++)  
  {
    StreamLattice[i][B*(ln)+A] = 1;
    //printf("StreamLattice[%d] = %d\n",i, Cells[A][B].StreamLattice[i]);
  }


  // DON'T STREAM FROM OUTSIDE OF THE DOMAIN
  for(k=0;k<9;k++)
  {
    if (BC_ID[k][B*(ln)+A]!=0)
    {
        StreamLattice[opp[k]][B*(ln)+A]= 0 ;
    }
  }


  // INLET VELOCITY
  switch(InletProfile)
  {
    case 1:
      //Uo=1.5*Uavg*(1-4*(pow((CoordY-MinInletCoordY-0.5*(MaxInletCoordY-MinInletCoordY)),2)));
      //Uo=4*1.5*Uavg*CoordY*(41-CoordY)/(41*41);
      Uo[B*(ln)+A] = 4*1.5*Uavg*((CoordY[B*(ln)+A])-(*MinInletCoordY))*(((*MaxInletCoordY)-(*MinInletCoordY))-((CoordY[B*(ln)+A])-(*MinInletCoordY)))/(((*MaxInletCoordY)-(*MinInletCoordY))*((*MaxInletCoordY)-(*MinInletCoordY)));
      Vo[B*(ln)+A] = Vavg;
      break;
    case 2:
      Uo[B*(ln)+A] = Uavg;
      Vo[B*(ln)+A] = Vavg;
    break;

  }
  U[B*(ln)+A] = Uo[B*(ln)+A];


} // End of function


  /*==================================================
  ========Creating constant lattice parameters========
  ==================================================*/

void D2Q9Vars(double* w, int* cx, int* cy, int* opp, int* c)
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

void BGKW(int i, double* w, int* cx, int* cy, double Omega)
{

  double T1 = (U[i]) * (U[i]) + 
              (V[i]) * (V[i]);
  double T2 = 0.0;
  int k;
  for (k=0; k<9; k++)
  {
    T2                          = (U[i])*cx[k] + (V[i])*cy[k];
    Feq[k][i]   = (Rho[i])*w[k]*( 1.0 + 3.0*T2 + 4.5*T2*T2 - 1.5*T1 );
    METAF[k][i] = Omega*(Feq[k][i])+(1.0-Omega)*(F[k][i]);

	/* Without T1 and T2 ::: longer execution time!
    Cells[j][i].Feq[k]   = Cells[j][i].Rho*w[k]*( 1.0 + 3.0*(Cells[j][i].U*cx[k] + Cells[j][i].V*cy[k])
    					   + 4.5*(Cells[j][i].U*cx[k] + Cells[j][i].V*cy[k])*(Cells[j][i].U*cx[k] + Cells[j][i].V*cy[k])
    					   - 1.5*(Cells[j][i].U * Cells[j][i].U + Cells[j][i].V * Cells[j][i].V) );
    Cells[j][i].METAF[k] = Omega*Cells[j][i].Feq[k]+(1.0-Omega)*Cells[j][i].F[k];*/

  }
}

  /*==================================================
  =============Function for the TRT model=============
  ==================================================*/

void TRT(int i, double* w, int* cx, int* cy, int* opp, double Omega, double OmegaA)
{
  double T1 = (U[i])*(U[i])+(V[i])*(V[i]);
  double T2 = 0.0;
  int k;
  for (k=0; k<9; k++)
  {
      T2                  = (U[i])   *cx[k] + (V[i])*cy[k];
      Feq[k][i]  = (Rho[i]) *w[k]*(1.0+3.0*T2+4.5*T2*T2-1.5*T1);
  }
  
  float F_p[9],Feq_p[9],F_m[9],Feq_m[9];
  for (k=0; k<9; k++)
  {
      F_p[k]   = 0.5*(   F[k][i] +   F[opp[k]][i] );
      Feq_p[k] = 0.5*( Feq[k][i] + Feq[opp[k]][i] );
      F_m[k]   = 0.5*(   F[k][i] -   F[opp[k]][i] );
      Feq_m[k] = 0.5*( Feq[k][i] - Feq[opp[k]][i] );
      
      METAF[k][i] = F[k][i] - 
                             (F_p[k]-Feq_p[k])*Omega - 
                             (F_m[k]-Feq_m[k])*OmegaA;
  }
}


  /*==================================================
  =============Function for the MRT model=============
  ==================================================*/

void MRT(int i, float** tm, float** stmiv)
{
  int k, l;
  float fmom[9],fmeq[9];
  float UU, VV;
  double suma,sumb;

  UU = U[i];
  VV = V[i];

  fmeq[0] =  Rho[i];
  fmeq[1] = (Rho[i])*(-2.0+3.0*(Rho[i])*(UU*UU+VV*VV));
  fmeq[2] = (Rho[i])*(1.0-3.0*(Rho[i])*(UU*UU+VV*VV));
  fmeq[3] = (Rho[i])*(U[i]);
  fmeq[4] =-(Rho[i])*(U[i]);
  fmeq[5] = (Rho[i])*(V[i]);
  fmeq[6] =-(Rho[i])*(V[i]);
  fmeq[7] = (Rho[i])*(UU*UU-VV*VV);
  fmeq[8] = (Rho[i])*UU*VV;
  
  for (k=0; k<9;k++)
  {
    suma=0.0;
    for (l=0; l<9;l++)
      suma = suma + tm[k][l]*(F[l][i]);

    fmom[k]=suma;
  }

  for (k=0; k<9;k++)
  {
    sumb=0.0;
    for (l=0; l<9;l++)
      sumb = sumb + stmiv[k][l]*(fmom[l]-fmeq[l]);
    
    METAF[k][i] = F[k][i] - sumb;
  }
}


  /*==================================================
  ======Function to update the distribution fct.======
  ==================================================*/

void UpdateF(int i)
{
  int k;
  for(k=0;k<9;k++)
    F[k][i] = METAF[k][i];
}


  /*==================================================
  ======Function to update the macroscopic var.=======
  ==================================================*/

void UpdateMacroscopic(int i, int* cx, int* cy, int CalculateDragLift)
{
  double Ssum, Usum, Vsum;
  int k;

  if (Fluid[i]==1)
  {
    Ssum=0.0;
    for (k=0; k<9; k++)
        Ssum = Ssum+F[k][i];

    Rho[i] = Ssum;

    //printf("Rho[%d][%d] = %f\n", j, i, Ssum);

    Usum = 0.0;
    Vsum = 0.0;
    for (k=0; k<9; k++)
    {
      Usum = Usum + (F[k][i])*cx[k];
      Vsum = Vsum + (F[k][i])*cy[k];
    }
    U[i] = Usum/(Rho[i]);
    V[i] = Vsum/(Rho[i]);
  }

  if (BC_ID[1][i]==3) // for outlet on the right
  {
      V[i]=0.0;
  }

  //   DRAG/LIFT FORCE
  if (CalculateDragLift != 0 && BoundaryID[i]==CalculateDragLift)
  {
     DragF[i] = (Rho[i])/3*(20-CoordX[i])/5;
     LiftF[i] = (Rho[i])/3*(20-CoordY[i])/5;
  }

}
/*
void CalculateDragLiftForces(struct CellProps *Cells, int j, int i, int CalculateDragLift, shared [] int* n, shared [] int* m)
{
  //   DRAG/LIFT FORCE
  if (CalculateDragLift != 0 && (Cells +j*(*m)+i)->BoundaryID==CalculateDragLift)
  {
     (Cells +j*(*m)+i)->DragF = ((Cells +j*(*m)+i)->Rho)/3*(20-(Cells +j*(*m)+i)->CoordX)/5;
     (Cells +j*(*m)+i)->LiftF = ((Cells +j*(*m)+i)->Rho)/3*(20-(Cells +j*(*m)+i)->CoordY)/5;
  }
}*/
