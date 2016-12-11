#include <cuda.h>                       // CUDA
#include <stdio.h>                      // printf();
#include <math.h>                       // need to compile with -lm
#include <stdlib.h>                     // for calloc();
#include <stdbool.h>                    // Include for bool type variables!
#include <string.h>                     // String operations
#include <time.h>                       // time functions
#include "include/ShellFunctions.h"     // For convenience
#include "include/FilesReading.h"       // For reading files
#include "include/FilesWriting.h"       // For writing files e.g. tecplot
#include "include/CellFunctions.h"      // For cell modifications
#include "include/ComputeResiduals.h"   // Residuals

#define threads 256

__constant__ int InletProfile_d[1];
__constant__ int NumNodes_d[1];
__constant__ int NumConn_d[1];
__constant__ int CurvedBoundaries_d[1];
__constant__ int OutletProfile_d[1];
__constant__ int CalculateDragLift_d[1];
__constant__ int cx_d[9]; 
__constant__ int cy_d[9];
__constant__ int width_d[1];
__constant__ int height_d[1];
__constant__ int c_d[9];
__constant__ int opp_d[9];
__constant__ float Delta_d[1];
__constant__ float w_d[9];
__constant__ float Qlat_d[9];
__constant__ float omega_d[1];
__constant__ float omegaA_d[1];
__constant__ float rho_ini_d[1];
__constant__ float Uavg_d[1];
__constant__ float Vavg_d[1];
__constant__ float MinInletCoordY_d[1];
__constant__ float MaxInletCoordY_d[1];
__constant__ float tm_d[81];
__constant__ float stmiv_d[81];

  ////////////////////////////////////////////////////
  ///////////////// KERNEL FCTS //////////////////////
  ////////////////////////////////////////////////////

  ////////////////////////////////////////////////////
  /////////////// Initialization1 ////////////////////
  ////////////////////////////////////////////////////

__global__ void gpu_init1(int *Fluid_d, int *ID_d, int *Boundary_d, int *BoundaryID_d, float *CoordX_d, float *CoordY_d,
                          int *BC_ID_d, float *Q_d, float *U_d, float *V_d, float *Rho_d,
                         int *Nodes0_d, int *Nodes1_d, float *Nodes2_d, float *Nodes3_d, int *Nodes4_d,
                         int *BCconn0_d, int *BCconn1_d, int *BCconn2_d, int *BCconn3_d, float *BCconn4_d,
                         float *BCconn5_d, int *BCconn6_d)
{

  int bidx=blockIdx.x;
  int tidx=threadIdx.x;

  int ind = tidx + bidx*blockDim.x;

  if (ind<((*width_d)*(*height_d)))
  {

    ID_d[ind] = ind;

    CoordX_d[ind] = Nodes2_d[ind];
    CoordY_d[ind] = Nodes3_d[ind];
  
    // CHECK FLUID OR NOT
    Fluid_d[ind] = Nodes4_d[ind];

    U_d[ind] = 0;
    V_d[ind] = 0;
    Rho_d[ind] = *rho_ini_d;

    Boundary_d[ind] = 0;
    BoundaryID_d[ind] = 0;
  }

  if (ind<(9*(*width_d)*(*height_d)))
  {
      BC_ID_d[ind]= 0  ; // IT IS NOT BOUNDARY LATTICE
      Q_d[ind]    = 0.5; 
  }
}


  ////////////////////////////////////////////////////
  /////////////// Initialization2 ////////////////////
  ////////////////////////////////////////////////////

__global__ void gpu_init2(int *Corner_d, int *Boundary_d, float *CoordY_d, int *StreamLattice_d, int *BC_ID_d, float *U_d, float *V_d, float *Uo_d, float *Vo_d,
                         int *Nodes0_d, int *Nodes1_d, float *Nodes2_d, float *Nodes3_d, int *Nodes4_d,
                         int *BCconn0_d, int *BCconn1_d, int *BCconn2_d, int *BCconn3_d, float *BCconn4_d,
                         float *BCconn5_d, int *BCconn6_d)
{

  int bidx=blockIdx.x;
  int tidx=threadIdx.x;

  int ind = tidx + bidx*blockDim.x;

  int ind_s = ind - ((*width_d)*(*height_d)) *( (int)(ind/((*width_d)*(*height_d))) );

  int ms=(*width_d)*(*height_d);

  int k;

  if (ind<((*width_d)*(*height_d)))
  {
    // not in the corner
    Corner_d[ind]=0;
  }

  if (ind<(9*(*width_d)*(*height_d)))
  {

    if (BC_ID_d[ind]!=0)
    {
      if (Boundary_d[ind_s]==0)// if Boundary condition in the node is 0 it becomes equal to the BC of the lattice direction
      {
        Boundary_d[ind_s]=BC_ID_d[ind];
      }
      else
      {// if in the same node there are lattice directions with different BC (corners) the BC of the node is WALL (assuming that it's impossibe to findoutlet and inlet together)
        if (Boundary_d[ind_s]<BC_ID_d[ind])
        {
          Boundary_d[ind_s]=1;
          Corner_d[ind_s]=1;
        }
        if (Boundary_d[ind_s]>BC_ID_d[ind])
        {
          Boundary_d[ind_s]=1;
          Corner_d[ind_s]=1;
        }
      }
    }


    // BC ON CORNERS IS WALL!!!! (this operation is useful for wall condition, which checks the single direction)
    if (Corner_d[ind_s]==1)
    {
      if (BC_ID_d[ind_s+(*width_d)*(*height_d)]!=0 && BC_ID_d[ind_s+2*(*width_d)*(*height_d)]!=0)
      {
          BC_ID_d[ind_s+5*(*width_d)*(*height_d)]=1;
      }
      if (BC_ID_d[ind_s+(*width_d)*(*height_d)]!=0 && BC_ID_d[ind_s+4*(*width_d)*(*height_d)]!=0)
      {
          BC_ID_d[ind_s+8*(*width_d)*(*height_d)]=1;
      }
      if (BC_ID_d[ind_s+2*(*width_d)*(*height_d)]!=0 && BC_ID_d[ind_s+3*(*width_d)*(*height_d)]!=0)
      {
          BC_ID_d[ind_s+6*(*width_d)*(*height_d)]=1;
      }
      if (BC_ID_d[ind_s+3*(*width_d)*(*height_d)]!=0 && BC_ID_d[ind_s+4*(*width_d)*(*height_d)]!=0)
      {
          BC_ID_d[ind_s+7*(*width_d)*(*height_d)]=1;
      }
    }

    // INITIALIZE STREAMING (STREAM EVERYWHERE)
    StreamLattice_d[ind] = 1;

  }

  if (ind<(9*(*width_d)*(*height_d)))
  {
    // DON'T STREAM FROM OUTSIDE OF THE DOMAIN
    for(k=0;k<9;k++)
    {
      if (BC_ID_d[ind_s+k*ms]!=0)
      {
        StreamLattice_d[ind_s+opp_d[k]]= 0 ;
      }
    }
  }

  if (ind<((*width_d)*(*height_d)))
  {  
  // INLET VELOCITY
    switch(*InletProfile_d)
    {
      case 1:
        Uo_d[ind_s] = 4*1.5*(*Uavg_d)*(CoordY_d[ind_s]-(*MinInletCoordY_d))*(((*MaxInletCoordY_d)-
                 (*MinInletCoordY_d))-(CoordY_d[ind_s]-(*MinInletCoordY_d)))/(((*MaxInletCoordY_d)-
                 (*MinInletCoordY_d))*((*MaxInletCoordY_d)-(*MinInletCoordY_d)));

        Vo_d[ind_s] = *Vavg_d;
      break;
      case 2:
        Uo_d[ind_s] = *Uavg_d;
        Vo_d[ind_s] = *Vavg_d;
      break;
      case 3:
        Uo_d[ind_s]= 0;
        Vo_d[ind_s] = 0;
      break;

    }
    U_d[ind_s] = Uo_d[ind_s];
    V_d[ind_s] = Vo_d[ind_s];

  } 


}


  ////////////////////////////////////////////////////
  ////////////////// BGKW model //////////////////////
  ////////////////////////////////////////////////////

__global__ void gpu_bgk(int* Fluid_d, float* Feq_d, float* Rho_d, float* U_d, float* V_d, float* METAF_d, float* F_d)
{

  int bidx=blockIdx.x;
  int tidx=threadIdx.x;

  int ind = tidx + bidx*blockDim.x;

  int ind_s = ind - ((*width_d)*(*height_d)) *( (int)(ind/((*width_d)*(*height_d))) );

  int ind_c = (int)(ind/((*width_d)*(*height_d)));

  if (ind<(9*(*width_d)*(*height_d)))
  {
    if (Fluid_d[ind_s]==1)
    {
      // Collision
      Feq_d[ind]  = Rho_d[ind_s]*w_d[ind_c]*( 1.0 + 3.0*(U_d[ind_s]*cx_d[ind_c] + V_d[ind_s]*cy_d[ind_c])
                + 4.5*(U_d[ind_s]*cx_d[ind_c] + V_d[ind_s]*cy_d[ind_c])
                * (U_d[ind_s]*cx_d[ind_c] + V_d[ind_s]*cy_d[ind_c])
                - 1.5*(U_d[ind_s] * U_d[ind_s] + V_d[ind_s]* V_d[ind_s]) );
      METAF_d[ind] = (*omega_d)*Feq_d[ind]+(1.0-(*omega_d))*F_d[ind];
    }
  }
}



  ////////////////////////////////////////////////////
  ////////////////// TRT model1 //////////////////////
  ////////////////////////////////////////////////////

__global__ void gpu_trt1(int* Fluid_d, float* Feq_d, float* Rho_d, float* U_d, float* V_d)
{
  
  int bidx=blockIdx.x;
  int tidx=threadIdx.x;

  int ind = tidx + bidx*blockDim.x;

  int ind_s = ind - ((*width_d)*(*height_d)) *( (int)(ind/((*width_d)*(*height_d))) );

  int ind_c = (int)(ind/((*width_d)*(*height_d)));

  if (ind<(9*(*width_d)*(*height_d)))
  {
    if (Fluid_d[ind_s]==1)
    {
      Feq_d[ind] = Rho_d[ind_s]*w_d[ind_c] * (1.0+3.0*(U_d[ind_s]*cx_d[ind_c]+V_d[ind_s]*cy_d[ind_c])
               + 4.5 *(U_d[ind_s]*cx_d[ind_c]+V_d[ind_s]*cy_d[ind_c]) * (U_d[ind_s]*cx_d[ind_c]+V_d[ind_s]*cy_d[ind_c])
               - 0.5*(U_d[ind_s]*U_d[ind_s]+V_d[ind_s]*V_d[ind_s]));
    }
  }
}

  ////////////////////////////////////////////////////
  ////////////////// TRT model2 //////////////////////
  ////////////////////////////////////////////////////

__global__ void gpu_trt2(int* Fluid_d, float* Feq_d, float* F_d, float* METAF_d)
{
  
  int bidx=blockIdx.x;
  int tidx=threadIdx.x;

  int ind = tidx + bidx*blockDim.x;

  int ind_s = ind - ((*width_d)*(*height_d)) *( (int)(ind/((*width_d)*(*height_d))) );

  int ind_c = (int)(ind/((*width_d)*(*height_d)));

  if (ind<(9*(*width_d)*(*height_d)))
  {
    if (Fluid_d[ind_s]==1)
    {
      METAF_d[ind] = F_d[ind] 
                   - ((0.5*(F_d[ind]   + F_d[ind_s+opp_d[ind_c]]))
                   - (0.5* (Feq_d[ind] + Feq_d[ind_s+opp_d[ind_c]])))*(*omega_d)
                   - ((0.5*(F_d[ind]   - F_d[ind_s+opp_d[ind_c]]))
                   - (0.5* (Feq_d[ind] - Feq_d[ind_s+opp_d[ind_c]])))*(*omegaA_d);

    }
  }
}


  ////////////////////////////////////////////////////
  ////////////////// MRT model1 //////////////////////
  ////////////////////////////////////////////////////

__global__ void gpu_mrt1(int* Fluid_d, float* Rho_d, float* U_d, float* V_d, float* F_d, float* fmeq_d, float* fmom_d)
{
  
  int l;

  int bidx=blockIdx.x;
  int tidx=threadIdx.x;

  int ind = tidx + bidx*blockDim.x;
  
  int ms=(*width_d)*(*height_d);
  int ind_s = ind - ((*width_d)*(*height_d)) *( (int)(ind/((*width_d)*(*height_d))) );

  int ind_c = (int)(ind/((*width_d)*(*height_d)));

  if (ind<(9*(*width_d)*(*height_d)))
  {
    if (Fluid_d[ind_s]==1)
    {

      fmeq_d[ind_s]      = Rho_d[ind_s];
      fmeq_d[ind_s+1*ms] = Rho_d[ind_s]*(-2.0+3.0*Rho_d[ind_s]*
                             (U_d[ind_s]*U_d[ind_s]+V_d[ind_s]*V_d[ind_s]));
      fmeq_d[ind_s+2*ms] = Rho_d[ind_s]*(1.0-3.0*Rho_d[ind_s]*
                             (U_d[ind_s]*U_d[ind_s]+V_d[ind_s]*V_d[ind_s]));
      fmeq_d[ind_s+3*ms] = Rho_d[ind_s]*U_d[ind_s];
      fmeq_d[ind_s+4*ms] =-Rho_d[ind_s]*U_d[ind_s];
      fmeq_d[ind_s+5*ms] = Rho_d[ind_s]*V_d[ind_s];
      fmeq_d[ind_s+6*ms] =-Rho_d[ind_s]*V_d[ind_s];
      fmeq_d[ind_s+7*ms] = Rho_d[ind_s]*(U_d[ind_s]*U_d[ind_s]-
                             V_d[ind_s]*V_d[ind_s]);
      fmeq_d[ind_s+8*ms] = Rho_d[ind_s]*U_d[ind_s]*V_d[ind_s];
  

      fmom_d[ind]=0;
      for (l=0; l<9;l++)
      fmom_d[ind]=fmom_d[ind] + tm_d[ind_c*9+l]*F_d[ind_s+l*ms];

    }
  }
}


  ////////////////////////////////////////////////////
  ////////////////// MRT model2 //////////////////////
  ////////////////////////////////////////////////////

__global__ void gpu_mrt2(int* Fluid_d, float* sumb_d, float* fmom_d, float* fmeq_d, float* METAF_d, float* F_d)
{
  
  int l;

  int bidx=blockIdx.x;
  int tidx=threadIdx.x;

  int ind = tidx + bidx*blockDim.x;
  
  int ms=(*width_d)*(*height_d);
  int ind_s = ind - ((*width_d)*(*height_d)) *( (int)(ind/((*width_d)*(*height_d))) );

  int ind_c = (int)(ind/((*width_d)*(*height_d)));

  if (ind<(9*(*width_d)*(*height_d)))
  {
    if (Fluid_d[ind_s]==1)
    {
      sumb_d[ind] = 0;

        for (l=0; l<9;l++)
          sumb_d[ind] = sumb_d[ind] + (stmiv_d[ind_c*9+l]*(fmom_d[ind_s+l*ms]-fmeq_d[ind_s+l*ms]));

        METAF_d[ind] = F_d[ind] - sumb_d[ind];

    }
  }
}





  ////////////////////////////////////////////////////
  //////////// UPDATE DISTR. FCT. ////////////////////
  ////////////////////////////////////////////////////

__global__ void gpu_update_f(int* Fluid_d, float* F_d, float* METAF_d)
{

  int bidx=blockIdx.x;
  int tidx=threadIdx.x;

  int ind = tidx + bidx*blockDim.x;;

  int ind_s = ind - ((*width_d)*(*height_d)) *( (int)(ind/((*width_d)*(*height_d))) );


  if (ind<(9*(*width_d)*(*height_d)))
  {
    if (Fluid_d[ind_s]==1)
    {
      // Update F
      F_d[ind] = METAF_d[ind];
    }
  }
}


  ////////////////////////////////////////////////////
  ///////////////// STREAMING ////////////////////////
  ////////////////////////////////////////////////////


__global__ void gpu_streaming(int* Fluid_d, int* StreamLattice_d, float* F_d, float* METAF_d)
{

  int bidx=blockIdx.x;
  int tidx=threadIdx.x;

  int ind = tidx + bidx*blockDim.x;

  int ind_s = ind - ((*width_d)*(*height_d)) *( (int)(ind/((*width_d)*(*height_d))) );

  int ind_c = (int)(ind/((*width_d)*(*height_d)));

  if (ind<(9*(*width_d)*(*height_d)))
  {
    if (Fluid_d[ind_s]==1)
    {
      // STREAMING 
      if ( (StreamLattice_d[ind]) == 1 )
      {
        F_d[ind] = METAF_d[ind+c_d[ind_c]];
      }
    }
  }
}



  ////////////////////////////////////////////////////
  //////////////// boundaries1 ///////////////////////
  ////////////////////////////////////////////////////

__global__ void gpu_boundaries1(int* Fluid_d, int* Boundary_d, int* BC_ID_d, float* F_d, float* Uo_d, float* Vo_d, int* Corner_d)
{

  int bidx=blockIdx.x;
  int tidx=threadIdx.x;

  int ind = tidx + bidx*blockDim.x;

  int ind_s = ind - ((*width_d)*(*height_d)) *( (int)(ind/((*width_d)*(*height_d))) );

  int ms=(*width_d)*(*height_d);

  if (ind<(9*(*width_d)*(*height_d)))
  {
    if (Fluid_d[ind_s]==1)
    {
    
      //inlet BC 
      if (Boundary_d[ind_s]==2) // if inlet boundary
      {
        
        /* if (_d[ind_s+ms].BC_ID==2 && _d[ind_s].Corner!=1)
        {
          // corner treatment has to be work out 
        }*/

        // inlet on the top (to be completed with the y velocity)
        if (BC_ID_d[ind_s+2*ms]==2 && Corner_d[ind_s]!=1)
        {

          F_d[ind_s+4*ms] = F_d[ind_s+2*ms];

          F_d[ind_s+8*ms] = F_d[ind_s+6*ms]
                                                    + (F_d[ind_s] + F_d[ind_s+1*ms]
                                                    + F_d[ind_s+3*ms]
                                                    + 2*(F_d[ind_s+2*ms]
                                                    + F_d[ind_s+6*ms]
                                                    + F_d[ind_s+5*ms]))*Uo_d[ind_s]/6.0;

          F_d[ind_s+7*ms] = F_d[ind_s+5*ms]
                                                    - (F_d[ind_s] + F_d[ind_s+1*ms]
                                                    + F_d[ind_s+3*ms]
                                                    + 2*(F_d[ind_s+2*ms]
                                                    + F_d[ind_s+6*ms]
                                                    + F_d[ind_s+5*ms]))*Uo_d[ind_s]/6.0;
        }

        // inlet on the left (to be completed with the x velocity)
        if (BC_ID_d[ind_s+3*ms]==2 && Corner_d[ind_s]!=1)
        {

            F_d[ind_s+1*ms] = F_d[ind_s+3*ms]
                              + 2*((F_d[ind_s]+F_d[ind_s+2*ms]
                              + F_d[ind_s+4*ms]
                              + 2.*(F_d[ind_s+3*ms]
                              + F_d[ind_s+6*ms]
                              + F_d[ind_s+7*ms]))
                              / (1.0-Uo_d[ind_s]))*Uo_d[ind_s]/3;

            F_d[ind_s+5*ms] = F_d[ind_s+7*ms]
                              +  ((F_d[ind_s]+F_d[ind_s+2*ms]
                              + F_d[ind_s+4*ms]
                              + 2.*(F_d[ind_s+3*ms]
                              + F_d[ind_s+6*ms] 
                              + F_d[ind_s+7*ms]))
                              / (1.0-Uo_d[ind_s]))*Uo_d[ind_s]/6;

            F_d[ind_s+8*ms] = F_d[ind_s+6*ms]
                              +  ((F_d[ind_s]+ F_d[ind_s+2*ms]
                              +F_d[ind_s+4*ms] 
                              + 2.*(F_d[ind_s+3*ms]
                              +F_d[ind_s+6*ms]
                              + F_d[ind_s+7*ms]))
                              / (1.0-Uo_d[ind_s]))*Uo_d[ind_s]/6;
        }
    
        /*if (_d[ind_s+4*ms].BC_ID==2 && _d[ind_s].Corner!=1)
        {
          // corner treatment has to be work out 
        }*/
      }
    }
  }  
}



  ////////////////////////////////////////////////////
  //////////////// boundaries2 ///////////////////////
  ////////////////////////////////////////////////////


__global__ void gpu_boundaries2(int* Fluid_d, float* Fneighbours_d, float* METAF_d, int* BC_ID_d, float* Q_d, float* F_d)
{

  int bidx=blockIdx.x;
  int tidx=threadIdx.x;

  int ind = tidx + bidx*blockDim.x;

  int ind_s = ind - ((*width_d)*(*height_d)) *( (int)(ind/((*width_d)*(*height_d))) );

  int ind_c = (int)(ind/((*width_d)*(*height_d)));

  if (ind<(9*(*width_d)*(*height_d)))
  {
    if (Fluid_d[ind_s]==1)
    {


      // WALL
      switch(*CurvedBoundaries_d)
      {
        // curved boundaries
      case 1:
 

          Fneighbours_d[ind]=METAF_d[ind+c_d[ind_c]];

          //CURVED WALL BC 

          if (BC_ID_d[ind]==1) // if wall
          {
            if (Q_d[ind]<0.5) // if the distance from the boundary is less than 0.5?
            { 
              F_d[ind_s+opp_d[ind_c]] = 2*Q_d[ind]*METAF_d[ind]
                                +(1-2*Q_d[ind])*Fneighbours_d[ind];
            }
            else
            {
              F_d[ind_s+opp_d[ind_c]] = METAF_d[ind]/2/Q_d[ind]
                                + (2*Q_d[ind]-1)/(2*Q_d[ind])
                                * METAF_d[ind_s+opp_d[ind_c]];
            }
          }

      break;

        // bounceback boundaries
      case 2:         
        //WALL BC (half-way bounceback)
        if (BC_ID_d[ind]==1) // if wall boundary
        {
          F_d[ind_s+opp_d[ind_c]] = F_d[ind];
        }
      break;
      }
    }
  }
}



  ////////////////////////////////////////////////////
  //////////////// boundaries3 ///////////////////////
  ////////////////////////////////////////////////////

__global__ void gpu_boundaries3(int* Fluid_d, int* BC_ID_d, float* F_d, float* Uo_d, float* Vo_d, int* Corner_d)
{

  int bidx=blockIdx.x;
  int tidx=threadIdx.x;

  int ind = tidx + bidx*blockDim.x;

  int ind_s = ind - ((*width_d)*(*height_d)) *( (int)(ind/((*width_d)*(*height_d))) );

  int ms=(*width_d)*(*height_d);

  if (ind<(9*(*width_d)*(*height_d)))
  {
    if (Fluid_d[ind_s]==1)
    {
      // OUTLET
      switch(*OutletProfile_d)
      {
        // set profile in outlet
        case 1:
        //outlet BC 
              if (BC_ID_d[ind_s+ms]==3) // if outlet boundary on the right side of the domain
          {

            F_d[ind_s+3*ms] = F_d[ind_s+1*ms]-2*((F_d[ind_s]
                                        + F_d[ind_s+2*ms]+F_d[ind_s+4*ms]
                                        + 2.*(F_d[ind_s+1*ms]+F_d[ind_s+5*ms]
                                        + F_d[ind_s+8*ms]))/(1-Uo_d[ind_s]))*Uo_d[ind_s]/3;

            F_d[ind_s+7*ms] = F_d[ind_s+5*ms]-((F_d[ind_s]
                                        + F_d[ind_s+2*ms]+F_d[ind_s+4*ms]
                                        + 2.*(F_d[ind_s+4*ms]+F_d[ind_s+5*ms]
                                        + F_d[ind_s+8*ms]))/(1-Uo_d[ind_s]))*Uo_d[ind_s]/6;

            F_d[ind_s+6*ms] = F_d[ind_s+8*ms]-((F_d[ind_s]
                                        + F_d[ind_s+2*ms]+F_d[ind_s+4*ms]
                                        + 2.*(F_d[ind_s+1*ms]+F_d[ind_s+5*ms]
                                        + F_d[ind_s+8*ms]))/(1-Uo_d[ind_s]))*Uo_d[ind_s]/6;
          }
          if (BC_ID_d[ind_s+2*ms]==3 && Corner_d[ind_s]!=1)
          {
          // FILL!!
          }
          if (BC_ID_d[ind_s+3*ms]==3 && Corner_d[ind_s]!=1)
          {
          // FILL!!
          }
          if (BC_ID_d[ind_s+4*ms]==3 && Corner_d[ind_s]!=1)
          {
          // FILL!!
          }
        break;

         // OPEN BOUNDARY
        case 2 :
          if (BC_ID_d[ind_s+ms]==3)
          {
            F_d[ind_s+ms] = 2*F_d[ind_s+ms-1]-F_d[ind_s+ms-2];
            F_d[ind_s+5*ms] = 2*F_d[ind_s+5*ms-1]-F_d[ind_s+5*ms-2];
            F_d[ind_s+8*ms] = 2*F_d[ind_s+8*ms-1]-F_d[ind_s+8*ms-2];
          }
          if (BC_ID_d[ind+2*ms]==3 && Corner_d[ind_s]!=1)
          {
          // FILL!!
          }
          if (BC_ID_d[ind+3*ms]==3 && Corner_d[ind_s]!=1)
          {
          // FILL!!
          }
          if (BC_ID_d[ind+4*ms]==3 && Corner_d[ind_s]!=1)
          {
          // FILL!!
          }
        break;
        // first order outlet
        case 3 :
          if (BC_ID_d[ind_s+ms]==3)
          {
            F_d[ind_s+ms] = F_d[ind_s-1+ms];
            F_d[ind_s+5*ms] = F_d[ind_s-1+5*ms];
            F_d[ind_s+8*ms] = F_d[ind_s-1+8*ms];
          }
          if (BC_ID_d[ind+2*ms]==3 && Corner_d[ind_s]!=1)
          {
          // FILL!!
          }
          if (BC_ID_d[ind+3*ms]==3 && Corner_d[ind_s]!=1)
          {
          // FILL!!
          }
          if (BC_ID_d[ind+4*ms]==3 && Corner_d[ind_s]!=1)
          {
          // FILL!!
          }
        break;
      }

    }
  }
}




  ////////////////////////////////////////////////////
  //////////// UPDATE MACROSCOPIC ////////////////////
  ////////////////////////////////////////////////////



__global__ void gpu_update_macro(int* Fluid_d, float* Rho_d, float* U_d, float* V_d,
                                int* BC_ID_d, int* BoundaryID_d, float* DragF_d,
                                float* LiftF_d, float* CoordX_d, float* CoordY_d, float* F_d)
{
  
  int k;

  int bidx=blockIdx.x;
  int tidx=threadIdx.x;

  int ind = tidx + bidx*blockDim.x;

  int ms = (*width_d)*(*height_d);

  if (ind<((*width_d)*(*height_d)))
  {
    if (Fluid_d[ind]==1)
    {
      // Update macroscopic
      Rho_d[ind]=0;
      U_d[ind]=0;
      V_d[ind]=0;
      for (k=0; k<9; k++)
      {        
        Rho_d[ind] = Rho_d[ind] + F_d[ind+k*ms];
        U_d[ind] = U_d[ind] + F_d[ind+k*ms]*cx_d[k];
        V_d[ind]= V_d[ind] + F_d[ind+k*ms]*cy_d[k];
      }

      U_d[ind] = U_d[ind] / Rho_d[ind];
      V_d[ind] = V_d[ind] / Rho_d[ind];

      if (BC_ID_d[ind+ms]==3) // for outlet on the right
      {
        V_d[ind]=0.0;
      }

      //   DRAG/LIFT FORCE
      if (*CalculateDragLift_d != 0 && BoundaryID_d[ind]==*CalculateDragLift_d)
      {
         DragF_d[ind] = 0.33333333*Rho_d[ind]*(20-CoordX_d[ind])*0.2;
         LiftF_d[ind] = 0.33333333*Rho_d[ind]*(20-CoordY_d[ind])*0.2;
      }
    }
  }
}




  ////////////////////////////////////////////////////
  ///////////////// SOLVER ///////////////////////////
  ////////////////////////////////////////////////////




void Iteration(char* NodeDataFile, char* BCconnectorDataFile,
               float Uavg,         float Vavg,         float rho_ini, float Viscosity,
               int InletProfile,   int CollisionModel, int CurvedBoundaries,
               int OutletProfile,  int Iterations,     int AutosaveAfter,
               int AutosaveEvery, int postproc_prog,   int CalculateDragLift)
{

  ////////////////////////////////////////////////////
  ///////////////////// Declare //////////////////////
  ////////////////////////////////////////////////////
  
  int bx; // blocks in X direction


  // Time measurement: declaration, begin
  clock_t tStart = clock();

  FILE* resid_file;                       // file for residuals
  FILE* log_file;                         // file for log
  FILE* TimeMeasurementFile;              // file for time measurement results
  // char IterOutputFile[50];                // write results to this file after the iterations
  char AutosaveOutputFile[50];            // autosave filename
  char OutputFile[50];                    // initial data will be written to this file
  char FinalOutputFile[50];               // final data will be written to this file
  char logFile[] = "Results/logFile.log"; // path of the .log file
  int i, j, A, B, ind, ind9, iter = 0;                  //variables for loops

  // Variables for residuals
  float *Residuals;

  float Qlat[9]={0.,1.,1.,1.,1.,sqrt(2),sqrt(2),sqrt(2),sqrt(2)};

  int AutosaveI = 1;      // autosave i variable, will be incremented after every autosave
  int* ppp;               // pointer of the postproc_prog variable
  int *NumNodes,*NumConn; // This will store the number of lines of the read files
  float *Delta;           // grid spacing
  int *n,*m;              // number of nodes in the x and y directions
  float *MaxInletCoordY; // maximum inlet coordinate in y
  float *MinInletCoordY; // minimum inlet coordinate in y
  int *NumInletNodes;     // number of inlet nodes
  int *Nodes0, *Nodes1, *Nodes4, *BCconn0, *BCconn1, *BCconn2, *BCconn3, *BCconn6; // vectors for the nodes and connections
  float *Nodes2, *Nodes3,*BCconn4, *BCconn5; // vectors for the nodes and connections
  float Omega, OmegaA;   // collision frequency from the viscosity


  float tInitialization = 0.0;  // Time measurement of Initialization
  float tIteration      = 0.0;  // Time measurement of Iteration
  float tCollision      = 0.0;  // Time measurement of Collision
  float tUpdateF        = 0.0;  // Time measurement of UpdateF
  float tStreaming      = 0.0;  // Time measurement of Streaming
  float tBoundaries     = 0.0;  // Time measurement of Boundaries
  float tUpdateMacro    = 0.0;  // Time measurement of Update Macroscopic vars
  float tResiduals      = 0.0;  // Time measurement of calculating residuals
  float tWriting        = 0.0;  // Time measurement of writing data

  clock_t tInstant1, tInstant2; // Time measurement points, universal
  clock_t tIterStart, tIterEnd; // Time measurement points: main loop

  // cuda time measurement variables
  cudaEvent_t start, stop;
  float cudatime;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);  
 

  ////////////////////////////////////////////////////
  ////////////////// Declare GPU /////////////////////
  ////////////////////////////////////////////////////

  // scalar on host
  float *omega_h, *omegaA_h;
  int *OutletProfile_h;
  int *CurvedBoundaries_h;
  int *CalculateDragLift_h;

  // vectors on device (variables
  int *Nodes0_d, *Nodes1_d, *Nodes4_d, *BCconn0_d, *BCconn1_d, *BCconn2_d, *BCconn3_d, *BCconn6_d; // vectors for the nodes and connections
  float *Nodes2_d, *Nodes3_d, *BCconn4_d, *BCconn5_d; // vectors for the nodes and connections
  


  ////////////////////////////////////////////////////
  //////////////////// Allocate //////////////////////
  ////////////////////////////////////////////////////

  // allocate residuals
  Residuals = Create1DArrayFloat(4); 

  // allocate mesh properties  
  Delta          = Create1DArrayFloat(1);
  m              = Create1DArrayInt(1);
  n              = Create1DArrayInt(1);
  NumNodes       = Create1DArrayInt(1);
  NumConn        = Create1DArrayInt(1);
  MaxInletCoordY = Create1DArrayFloat(1);
  MinInletCoordY = Create1DArrayFloat(1);
  NumInletNodes  = Create1DArrayInt(1);


  // D2Q9 Variables of the lattice
  float* w = Create1DArrayFloat(9); // weight values for the directions
  int*   c_h = Create1DArrayInt(9);    // x coordinate of the discrete lattice directions
  int*   cx = Create1DArrayInt(9);    // x coordinate of the discrete lattice directions
  int*   cy = Create1DArrayInt(9);    // y coordinate of the discrete lattice directions
  int*  opp = Create1DArrayInt(9);    // opposite vector




  ////////////////////////////////////////////////////
  ///////////////////// Read data ////////////////////
  ////////////////////////////////////////////////////

  ReadNodLines(NodeDataFile, NumNodes);

  // allocate vectors for nodes
  Nodes0 = (int *)calloc(*(NumNodes),sizeof(int));
  Nodes1 = (int *)calloc(*(NumNodes),sizeof(int));
  Nodes2 = (float *)calloc(*(NumNodes),sizeof(float));
  Nodes3 = (float *)calloc(*(NumNodes),sizeof(float));
  Nodes4 = (int *)calloc(*(NumNodes),sizeof(int));

  cudaMemcpyToSymbol(NumNodes_d, NumNodes, sizeof(int));

  cudaMalloc((void**)&Nodes0_d, (*NumNodes) * sizeof(int));
  cudaMalloc((void**)&Nodes1_d, (*NumNodes) * sizeof(int));
  cudaMalloc((void**)&Nodes2_d, (*NumNodes) * sizeof(float));
  cudaMalloc((void**)&Nodes3_d, (*NumNodes) * sizeof(float));
  cudaMalloc((void**)&Nodes4_d, (*NumNodes) * sizeof(int));



  // Read Node data
  ReadNodes(NodeDataFile, NumNodes, Nodes0, Nodes1, Nodes2, Nodes3, Nodes4);

  ReadBCconLines(BCconnectorDataFile, NumConn);

  // allocate vectors for connectors
  BCconn0 = (int *)calloc(*(NumConn),sizeof(int));
  BCconn1 = (int *)calloc(*(NumConn),sizeof(int));
  BCconn2 = (int *)calloc(*(NumConn),sizeof(int));
  BCconn3 = (int *)calloc(*(NumConn),sizeof(int));
  BCconn4 = (float *)calloc(*(NumConn),sizeof(float));
  BCconn5 = (float *)calloc(*(NumConn),sizeof(float));
  BCconn6 = (int *)calloc(*(NumConn),sizeof(int)); 

  cudaMemcpyToSymbol(NumConn_d, NumConn, sizeof(int));

  cudaMalloc((void**)&BCconn0_d, (*NumConn) * sizeof(int));
  cudaMalloc((void**)&BCconn1_d, (*NumConn) * sizeof(int));
  cudaMalloc((void**)&BCconn2_d, (*NumConn) * sizeof(int));
  cudaMalloc((void**)&BCconn3_d, (*NumConn) * sizeof(int));
  cudaMalloc((void**)&BCconn4_d, (*NumConn) * sizeof(float));
  cudaMalloc((void**)&BCconn5_d, (*NumConn) * sizeof(float));
  cudaMalloc((void**)&BCconn6_d, (*NumConn) * sizeof(int));

  // Read BCconn data
  ReadBCconn(BCconnectorDataFile, NumConn, BCconn0, BCconn1, BCconn2, BCconn3, BCconn4, BCconn5, BCconn6); 

  ////////////////////////////////////////////////////
  /////////////// Compute constants //////////////////
  ////////////////////////////////////////////////////

  CompDataNode(Delta, m,  n, Nodes0, Nodes1, Nodes2, Nodes3, Nodes4, NumNodes);

  cudaMemcpyToSymbol(Delta_d, Delta, sizeof(float));

  CompDataConn(NumInletNodes, MaxInletCoordY,
               MinInletCoordY, BCconn0, BCconn1, BCconn2, BCconn3, BCconn4,
               BCconn5, BCconn6, NumConn, Delta);

  cudaMemcpy(Nodes0_d, Nodes0, (*NumNodes)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(Nodes1_d, Nodes1, (*NumNodes)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(Nodes2_d, Nodes2, (*NumNodes)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(Nodes3_d, Nodes3, (*NumNodes)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(Nodes4_d, Nodes4, (*NumNodes)*sizeof(int), cudaMemcpyHostToDevice);

  cudaMemcpy(BCconn0_d, BCconn0, (*NumConn)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(BCconn1_d, BCconn1, (*NumConn)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(BCconn2_d, BCconn2, (*NumConn)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(BCconn3_d, BCconn3, (*NumConn)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(BCconn4_d, BCconn4, (*NumConn)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(BCconn5_d, BCconn5, (*NumConn)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(BCconn6_d, BCconn6, (*NumConn)*sizeof(int), cudaMemcpyHostToDevice);

  cudaMemcpyToSymbol(InletProfile_d, &InletProfile, sizeof(int));
  cudaMemcpyToSymbol(rho_ini_d, &rho_ini, sizeof(float));
  cudaMemcpyToSymbol(Uavg_d, &Uavg, sizeof(float));
  cudaMemcpyToSymbol(Vavg_d, &Vavg, sizeof(float));
  cudaMemcpyToSymbol(MinInletCoordY_d, MinInletCoordY, sizeof(float));
  cudaMemcpyToSymbol(MaxInletCoordY_d, MaxInletCoordY, sizeof(float));

  // Check whether we got back what we wanted :), wtite to log file!
  log_file = fopen(logFile, "w");  // open log file
  ppp      = &postproc_prog;       // for convenience ppp points to postproc_prog
  fprintf(log_file,"This is the 2D lattice Boltzmann *.log file\n\n");
  fprintf(log_file,"\n:::: Imported variables from the *.ini file :::: \n");
  fprintf(log_file,">>> Uavg              : %3.6f\n", Uavg);
  fprintf(log_file,">>> Vavg              : %3.6f\n", Vavg);
  fprintf(log_file,">>> Initial density   : %2.1f\n", rho_ini);
  fprintf(log_file,">>> Viscosity         : %3.8f\n", Viscosity);
  fprintf(log_file,">>> # of iterations   : %1.1d\n", Iterations);
  fprintf(log_file,">>> Autosave after    : %1.1d\n", AutosaveAfter);
  fprintf(log_file,">>> Autosave every    : %1.1d\n", AutosaveEvery);

  fprintf(log_file,">>> Convergence Veloc : UNKNOWN OPTION\n");
  fprintf(log_file,">>> Convergence Rho   : UNKNOWN OPTION\n");
  fprintf(log_file,">>> Former data iter# : UNKNOWN OPTION\n");


  switch(CollisionModel)         // 1: BGKW, 2: TRT, 3: MRT
  {
    case 1: fprintf(log_file,">>> CollisionModel   : BGKW\n"); break;
    case 2: fprintf(log_file,">>> CollisionModel   : TRT\n" ); break;
    case 3: fprintf(log_file,">>> CollisionModel   : MRT\n" ); break;
  }
  switch(InletProfile)                      // 1:ON, 2:OFF
  {
    case 1: fprintf(log_file,">>> InletProfile     : ON\n" ); break;
    case 2: fprintf(log_file,">>> InletProfile     : OFF\n"); break;
  }
  switch(OutletProfile)                     // 1:ON, 2:OFF
  {
    case 1: fprintf(log_file,">>> OutletProfile    : ON\n" ); break;
    case 2: fprintf(log_file,">>> OutletProfile    : OFF\n"); break;
  }
  switch(CurvedBoundaries)                  // 1:ON, 2:OFF
  {
    case 1: fprintf(log_file,">>> CurvedBoundaries : ON\n" ); break;
    case 2: fprintf(log_file,">>> CurvedBoundaries : OFF\n"); break;
  }
  switch(postproc_prog)   // 1->Paraview (*.csv)     2->Tecplot 
  {
    case 1: fprintf(log_file,">>> Results format   : Paraview (*.csv)\n" ); break;
    case 2: fprintf(log_file,">>> Results format   : Tecplot (*.dat)\n"); break;
  }
  if (CalculateDragLift != 0)
            fprintf(log_file,">>> Drag, lift @ BC  : %d\n", CalculateDragLift);
  else 
            fprintf(log_file,">>> Drag, lift was not calculated\n");


  fprintf(log_file,"\n:::: Calculated variables from mesh :::: \n");
  fprintf(log_file,">>> Grid spacing        = %f\n", *Delta);
  fprintf(log_file,">>> # of nodes in x (n) = %d\n", *n);
  fprintf(log_file,">>> # of nodes in y (m) = %d\n", *m);
  fprintf(log_file,">>> NumInletNodes       = %d\n", *NumInletNodes);
  fprintf(log_file,">>> MaxInletCoordY      = %f\n", *MaxInletCoordY);
  fprintf(log_file,">>> MinInletCoordY      = %f\n", *MinInletCoordY);
  // End of checking

  ////////////////////////////////////////////////////
  ///////////////// Allocate GPU /////////////////////
  ////////////////////////////////////////////////////

  // host scalar
  omega_h = Create1DArrayFloat(1);
  omegaA_h = Create1DArrayFloat(1);
  OutletProfile_h = Create1DArrayInt(1);
  CurvedBoundaries_h = Create1DArrayInt(1);
  CalculateDragLift_h = Create1DArrayInt(1);
  

  // declare and allocate 1D constant gpu and cpu vectors
  int *Fluid, *Corner, *ID, *Boundary, *BoundaryID;
  float *CoordX, *CoordY;

  Fluid = Create1DArrayInt((*n)*(*m));
  Corner = Create1DArrayInt((*n)*(*m));          
  ID = Create1DArrayInt((*n)*(*m));             
  Boundary = Create1DArrayInt((*n)*(*m));    
  BoundaryID = Create1DArrayInt((*n)*(*m));  
  CoordX = Create1DArrayFloat((*n)*(*m));   
  CoordY = Create1DArrayFloat((*n)*(*m));

  int *Fluid_d, *Corner_d, *ID_d, *Boundary_d, *BoundaryID_d;
  float *CoordX_d, *CoordY_d;

  cudaMalloc((void**)&Fluid_d, (*m)*(*n) * sizeof(int));
  cudaMalloc((void**)&Corner_d, (*m)*(*n) * sizeof(int));
  cudaMalloc((void**)&ID_d, (*m)*(*n) * sizeof(int));
  cudaMalloc((void**)&Boundary_d, (*m)*(*n) * sizeof(int));
  cudaMalloc((void**)&BoundaryID_d, (*m)*(*n) * sizeof(int));
  cudaMalloc((void**)&CoordX_d, (*m)*(*n) * sizeof(float));
  cudaMalloc((void**)&CoordY_d, (*m)*(*n) * sizeof(float));

// declare and allocate 9D constant gpu and cpu vectors
  int *StreamLattice, *BC_ID;
  float *Q;  

  StreamLattice = Create1DArrayInt(9*(*n)*(*m));
  BC_ID = Create1DArrayInt(9*(*n)*(*m));          
  Q = Create1DArrayFloat(9*(*n)*(*m));             

  int *StreamLattice_d, *BC_ID_d;
  float *Q_d;

  cudaMalloc((void**)&StreamLattice_d, 9*(*m)*(*n) * sizeof(int));
  cudaMalloc((void**)&BC_ID_d, 9*(*m)*(*n) * sizeof(int));
  cudaMalloc((void**)&Q_d, 9*(*m)*(*n) * sizeof(float));


 // declare and allocate 1D variable gpu and cpu vectors
  float *U, *V, *Rho, *Uo, *Vo, *DragF, *LiftF;

  U = Create1DArrayFloat((*n)*(*m));
  V = Create1DArrayFloat((*n)*(*m));
  Rho = Create1DArrayFloat((*n)*(*m));
  Uo = Create1DArrayFloat((*n)*(*m));
  Vo = Create1DArrayFloat((*n)*(*m));
  DragF = Create1DArrayFloat((*n)*(*m));
  LiftF = Create1DArrayFloat((*n)*(*m));


  float *U_d, *V_d, *Rho_d, *Uo_d, *Vo_d, *DragF_d, *LiftF_d;

  cudaMalloc((void**)&U_d, (*m)*(*n) * sizeof(float));
  cudaMalloc((void**)&V_d, (*m)*(*n) * sizeof(float));
  cudaMalloc((void**)&Rho_d, (*m)*(*n) * sizeof(float));
  cudaMalloc((void**)&Uo_d, (*m)*(*n) * sizeof(float));
  cudaMalloc((void**)&Vo_d, (*m)*(*n) * sizeof(float));
  cudaMalloc((void**)&DragF_d, (*m)*(*n) * sizeof(float));
  cudaMalloc((void**)&LiftF_d, (*m)*(*n) * sizeof(float));


// declare and allocate 9D constant gpu and cpu vectors
  float *F, *Feq, *METAF, *Fneighbours, *fmom, *fmeq, *sumb;  
    
  F = Create1DArrayFloat(9*(*n)*(*m));
  Feq = Create1DArrayFloat(9*(*n)*(*m));  
  METAF = Create1DArrayFloat(9*(*n)*(*m));  
  Fneighbours = Create1DArrayFloat(9*(*n)*(*m));  
  fmom = Create1DArrayFloat(9*(*n)*(*m));  
  fmeq = Create1DArrayFloat(9*(*n)*(*m));  
  sumb = Create1DArrayFloat(9*(*n)*(*m));               

  float *F_d, *Feq_d, *METAF_d, *Fneighbours_d, *fmom_d, *fmeq_d, *sumb_d;

  cudaMalloc((void**)&F_d, 9*(*m)*(*n) * sizeof(float));
  cudaMalloc((void**)&Feq_d, 9*(*m)*(*n) * sizeof(float));
  cudaMalloc((void**)&METAF_d, 9*(*m)*(*n) * sizeof(float));
  cudaMalloc((void**)&Fneighbours_d, 9*(*m)*(*n) * sizeof(float));
  cudaMalloc((void**)&fmom_d, 9*(*m)*(*n) * sizeof(float));
  cudaMalloc((void**)&fmeq_d, 9*(*m)*(*n) * sizeof(float));
  cudaMalloc((void**)&sumb_d, 9*(*m)*(*n) * sizeof(float));


  ////////////////////////////////////////////////////
  ///////////////// INITIALIZE ///////////////////////
  ////////////////////////////////////////////////////

  // In case of no autosave
  sprintf(AutosaveOutputFile, "NOWHERE!");


  //CONSTANT LATTICE QUANTITIES
  w[0]=4./9.;

  for (i=1; i<5; i++ )
      w[i]=1./9.;

  for (i=5; i<9; i++ )
      w[i]=1./36.;


  opp[0] = 0;
  opp[1] = 3*(*m)*(*n);
  opp[2] = 4*(*m)*(*n);
  opp[3] = 1*(*m)*(*n);
  opp[4] = 2*(*m)*(*n);
  opp[5] = 7*(*m)*(*n);
  opp[6] = 8*(*m)*(*n);
  opp[7] = 5*(*m)*(*n);
  opp[8] = 6*(*m)*(*n);

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

  c_h[0]=0;
  c_h[1]=-1;
  c_h[2]=-1*(*n);
  c_h[3]=1;
  c_h[4]=(*n);
  c_h[5]=-1*(*n)-1;
  c_h[6]=-1*(*n)+1;
  c_h[7]=(*n)+1;
  c_h[8]=(*n)-1;
  
  // Calculate collision freq
  Omega  = 1.0/(3.*Viscosity+0.5);
  OmegaA = 8*(2-Omega)/(8-Omega);

  omega_h=&Omega;
  omegaA_h=&OmegaA;
  OutletProfile_h=&OutletProfile;
  CalculateDragLift_h=&CalculateDragLift;
  CurvedBoundaries_h=&CurvedBoundaries;

  ////////////////////////////////////////////////////
  ///////////// COPY CONSTANTS TO GPU ////////////////
  ////////////////////////////////////////////////////

  cudaMemcpyToSymbol(cx_d, cx, 9*sizeof(int));
  cudaMemcpyToSymbol(cy_d, cy, 9*sizeof(int));
  cudaMemcpyToSymbol(c_d, c_h, 9*sizeof(int));
  cudaMemcpyToSymbol(opp_d, opp, 9*sizeof(int));

  cudaMemcpyToSymbol(OutletProfile_d, OutletProfile_h, sizeof(int));
  cudaMemcpyToSymbol(CurvedBoundaries_d, CurvedBoundaries_h, sizeof(int));
  cudaMemcpyToSymbol(CalculateDragLift_d, CalculateDragLift_h, sizeof(int));

  cudaMemcpyToSymbol(width_d, m, sizeof(int));
  cudaMemcpyToSymbol(height_d, n, sizeof(int));
  cudaMemcpyToSymbol(Qlat_d, Qlat, 9*sizeof(float));
  cudaMemcpyToSymbol(w_d, w, 9*sizeof(float));
  cudaMemcpyToSymbol(omega_d, omega_h, sizeof(float));
  cudaMemcpyToSymbol(omegaA_d, omegaA_h, sizeof(float));

  // Initialize variables for MRT Collision model, if used
  if (CollisionModel == 3)
  {
    float *tm, *stmiv;    // variables for the MRT collision model
    tm    = Create1DArrayFloat(81);
    stmiv = Create1DArrayFloat(81);
    MRTInitializer(tm, stmiv, Omega);

    cudaMemcpyToSymbol(tm_d, tm, 81*sizeof(float));
    cudaMemcpyToSymbol(stmiv_d, stmiv, 81*sizeof(float));

  }



  // Create structure for the cell properties (see ShellFunctions.h)

  bx=(int)(9*(*m)*(*n)/threads)+1;

  dim3 tpb(threads); // threads/block
  dim3 bpg(bx); // blocks/grid

  // initializing matrix of struct Cells
  fprintf(log_file,"\n:::: Initializing ::::\n");
  printf("\n:::: Initializing ::::\n");
  tInstant1 = clock(); // Measure time of initialization
  

  ////////////////////////////////////////////////////
  ///////////// COPY VARIABLES TO GPU ////////////////
  ////////////////////////////////////////////////////
  
/*
  cudaMemcpy(Fluid_d, Fluid, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(Corner_d, Corner, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(ID_d, ID, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(Boundary_d, Boundary, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(BoundaryID_d, BoundaryID, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(CoordX_d, CoordX, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(CoordY_d, CoordY, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(StreamLattice_d, StreamLattice, 9*(*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(BC_ID_d, BC_ID, 9*(*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(Q_d, Q, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(U_d, U, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(V_d, V, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(Rho_d, Rho, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(Uo_d, Uo, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(Vo_d, Vo, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(DragF_d, DragF, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(LiftF_d, LiftF, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(F_d, F_d, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(Feq_d, Feq, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(METAF_d, METAF, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(Fneighbours_d, Fneighbours, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(fmom_d, fmom, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(fmeq_d, fmeq, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(sumb_d, sumb, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
*/

/*
  cudaMemcpy(Fluid_d, Fluid, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(ID_d, ID, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(Boundary_d, Boundary, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(BoundaryID_d, BoundaryID, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(CoordX_d, CoordX, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(CoordY_d, CoordY, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(BC_ID_d, BC_ID, 9*(*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(Q_d, Q, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(U_d, U, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(V_d, V, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(Rho_d, Rho, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
*/


  cudaMemcpy(Fluid_d, Fluid, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(Corner_d, Corner, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(ID_d, ID, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(Boundary_d, Boundary, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(BoundaryID_d, BoundaryID, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(CoordX_d, CoordX, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(CoordY_d, CoordY, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(StreamLattice_d, StreamLattice, 9*(*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(BC_ID_d, BC_ID, 9*(*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(Q_d, Q, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(U_d, U, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(V_d, V, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(Rho_d, Rho, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(Uo_d, Uo, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(Vo_d, Vo, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(DragF_d, DragF, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(LiftF_d, LiftF, (*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);

  cudaMemcpy(F_d, F, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(Feq_d, Feq, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(METAF_d, METAF, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(Fneighbours_d, Fneighbours, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);

switch(CollisionModel)
    {
      
      case 1:
      break;
      
      case 2:
      break;
      
      case 3:
          cudaMemcpy(fmom_d, fmom, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
          cudaMemcpy(fmeq_d, fmeq, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
          cudaMemcpy(sumb_d, sumb, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);
      break;
      
    }




  gpu_init1<<<bpg,tpb>>>(Fluid_d, ID_d, Boundary_d, BoundaryID_d, CoordX_d, CoordY_d, BC_ID_d, Q_d, U_d, V_d, Rho_d,
                        Nodes0_d, Nodes1_d, Nodes2_d, Nodes3_d, Nodes4_d, BCconn0_d,
                        BCconn1_d, BCconn2_d, BCconn3_d, BCconn4_d, BCconn5_d, BCconn6_d);

  cudaMemcpy(BC_ID, BC_ID_d, 9*(*m)*(*n)*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(BoundaryID, BoundaryID_d, (*m)*(*n)*sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(Q, Q_d, 9*(*m)*(*n)*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(CoordX, CoordX_d, (*m)*(*n)*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(CoordY, CoordY_d, (*m)*(*n)*sizeof(float), cudaMemcpyDeviceToHost);


  for(B=0;B<*m;B++)
  {
  for(A=0;A<*n;A++)
  {
    ind=A+B*(*n);

    for(i=0;i<*NumConn;i++)
    {
      if ( ( *(BCconn0+i) == A ) && ( *(BCconn1+i) == B ) )
      {
          for(j=1; j<9;j++)
          {
              if ( *(BCconn2+i) == j )
              {
                ind9=A+B*(*n)+(*m)*(*n)*j;
                *(BC_ID + ind9)   = *(BCconn3+i);
                *(BoundaryID + ind) = *(BCconn6+i);

                // find distance from the boundary
                *(Q + ind9) = sqrt(pow( *(BCconn4+i)-*(CoordX+ind),2 ) + pow( *(BCconn5+i)-*(CoordY+ind),2) ) / ((*Delta)*Qlat[j]);
              }
          }
      }
    }
  }
  }  

  cudaMemcpy(BC_ID_d, BC_ID, 9*(*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(BoundaryID_d, BoundaryID, (*m)*(*n)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(Q_d, Q, 9*(*m)*(*n)*sizeof(float), cudaMemcpyHostToDevice);


  gpu_init2<<<bpg,tpb>>>(Corner_d, Boundary_d, CoordY_d, StreamLattice_d, BC_ID_d, U_d, V_d, Uo_d, Vo_d,
                        Nodes0_d, Nodes1_d, Nodes2_d, Nodes3_d, Nodes4_d,
                        BCconn0_d, BCconn1_d, BCconn2_d, BCconn3_d, BCconn4_d,
                        BCconn5_d, BCconn6_d );

  cudaMemcpy(U, U_d, (*m)*(*n)*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(V, V_d, (*m)*(*n)*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(Rho, Rho_d, (*m)*(*n)*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(Fluid, Fluid_d, (*m)*(*n)*sizeof(int), cudaMemcpyDeviceToHost);


  tInstant2 = clock(); // Measure time of initialization
  tInitialization = (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;

  fprintf(log_file,"\n:::: Initialization done! ::::\n");

  printf("Initialization took %f seconds\n",         tInitialization);


  // Write Initialized data 
  switch(postproc_prog)
    {
      case 1: sprintf(OutputFile, "Results/InitialData.csv");   break;
      case 2: sprintf(OutputFile, "Results/InitialData.dat");   break;
    }

  //WriteResults(OutputFile, Cells_const_h, Cells_const_9d_h, Cells_var_h, Cells_var_9d_h, n, m, ppp);
  WriteResults(OutputFile, CoordX, CoordY, U, V, Rho, Fluid, n, m, ppp);
  
  printf("\nInitialized data was written to %s\n", OutputFile);

  // Open residuals file
  resid_file = fopen("Results/residuals.dat", "w");
  fprintf(resid_file,"Iter L2_norm L2_norm_weighted Drag Lift\n");


  ////////////////////////////////////////////////////
  /////////////////// ITERATE ////////////////////////
  ////////////////////////////////////////////////////

  fprintf(log_file,"\n:::: Start Iterations ::::\n");
  printf("\n:::: Start Iterations ::::\n");

  printf("%d is the number of iterations \n",         Iterations);


  tIterStart = clock(); // Start measuring time of main loop
  while (iter<Iterations)
  {
    
    ////////////// COLLISION ///////////////
    cudaEventRecord(start, 0); // Start measuring time

    switch(CollisionModel)
    {
      
      case 1:
        gpu_bgk<<<bpg,tpb>>>(Fluid_d, Feq_d, Rho_d, U_d, V_d, METAF_d, F_d);
      break;
      
      case 2:
        gpu_trt1<<<bpg,tpb>>>(Fluid_d, Feq_d, Rho_d, U_d, V_d);
        gpu_trt2<<<bpg,tpb>>>(Fluid_d, Feq_d, F_d, METAF_d);
      break;
      
      case 3:
        gpu_mrt1<<<bpg,tpb>>>(Fluid_d, Rho_d, U_d, V_d, F_d, fmeq_d, fmom_d);
        gpu_mrt2<<<bpg,tpb>>>(Fluid_d, sumb_d, fmom_d, fmeq_d, METAF_d, F_d);
      break;
      
    }

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&cudatime, start, stop);
    tCollision = tCollision + cudatime;


    ////////////// UPDATE DISTR ///////////////
  	cudaEventRecord(start, 0);

    gpu_update_f<<<bpg,tpb>>>(Fluid_d, F_d, METAF_d);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&cudatime, start, stop);
    tUpdateF = tUpdateF + cudatime;

    
    ////////////// STREAMING ///////////////
    cudaEventRecord(start, 0);

    gpu_streaming<<<bpg,tpb>>>(Fluid_d, StreamLattice_d, F_d, METAF_d);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&cudatime, start, stop);
    tStreaming = tStreaming + cudatime;

  // make the host block until the device is finished with foo
  cudaThreadSynchronize();

  // check for error
  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
  {
    // print the CUDA error message and exit
    printf("CUDA error: %s\n", cudaGetErrorString(error));
    exit(-1);
  }

    ////////////// BOUNDARIES ///////////////
    cudaEventRecord(start, 0);

    gpu_boundaries1<<<bpg,tpb>>>(Fluid_d, Boundary_d, BC_ID_d, F_d, Uo_d, Vo_d, Corner_d);

    gpu_boundaries2<<<bpg,tpb>>>(Fluid_d, Fneighbours_d, METAF_d, BC_ID_d, Q_d, F_d);
    gpu_boundaries3<<<bpg,tpb>>>(Fluid_d, BC_ID_d, F_d, Uo_d, Vo_d, Corner_d);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&cudatime, start, stop);
    tBoundaries = tBoundaries + cudatime;


    // UPDATE VELOCITY AND DENSITY
    cudaEventRecord(start, 0);

    gpu_update_macro<<<bpg,tpb>>>(Fluid_d, Rho_d, U_d, V_d, BC_ID_d, BoundaryID_d, DragF_d, LiftF_d, CoordX_d, CoordY_d, F_d);

    tInstant2 = clock();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&cudatime, start, stop);
    tUpdateMacro = tUpdateMacro + cudatime;
    

    ////////////////////////////////////////////////////
    //////////// COPY VARIABLES TO HOST ////////////////
    ////////////////////////////////////////////////////
/*
    cudaMemcpy(F, F_d, (*m)*(*n)*9*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(METAF, METAF_d, (*m)*(*n)*9*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(DragF, DragF_d, (*m)*(*n)*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(LiftF, LiftF_d, (*m)*(*n)*sizeof(float), cudaMemcpyDeviceToHost);
*/
    ////////////// Residuals ///////////////
/*
    tInstant1 = clock(); // Start measuring time
    
    ComputeResiduals(BoundaryID, F, METAF, DragF, LiftF, Residuals, m, n, CalculateDragLift);
    fprintf(resid_file,"%d %5.4e %5.4e %f %f\n", iter, Residuals[0], Residuals[1], Residuals[2], Residuals[3]);
    
    tInstant2 = clock();
    tResiduals = tResiduals + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;
*/
    printf("Iterating... %d/%d (%3.1f %%)\r", iter+1, Iterations, (float)(iter+1)*100/(float)(Iterations));

    // update loop variable
    iter++;  
   
    ////////////// Autosave ///////////////
    
    if(iter == (AutosaveEvery*AutosaveI))
    {
      AutosaveI++;
      if(iter>AutosaveAfter)
      {

        ////////////////////////////////////////////////////
        //////////// COPY VARIABLES TO HOST ////////////////
        ////////////////////////////////////////////////////

        cudaMemcpy(U, U_d, (*m)*(*n)*sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(V, V_d, (*m)*(*n)*sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(Rho, Rho_d, (*m)*(*n)*sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(Fluid, Fluid_d, (*m)*(*n)*sizeof(int), cudaMemcpyDeviceToHost);
  
        switch(postproc_prog)
        {
          case 1: sprintf(AutosaveOutputFile, "Results/autosave_iter%05d.csv", iter);
          break;
          case 2: sprintf(AutosaveOutputFile, "Results/autosave_iter%05d.dat", iter);
          break;
        }

      tInstant1 = clock(); // Start measuring time
      WriteResults(AutosaveOutputFile, CoordX, CoordY, U, V, Rho, Fluid, n, m, ppp);
      tInstant2 = clock();
      tWriting = tWriting + (float)(tInstant2-tInstant1) / CLOCKS_PER_SEC;
      }
    }
  

  }     ////////////// END OF MAIN WHILE CYCLE! ///////////////

  tIterEnd = clock(); // End measuring time of main loop
  tIteration = (float)(tIterEnd - tIterStart ) / CLOCKS_PER_SEC;

  ////////////////////////////////////////////////////
  //////////// COPY VARIABLES TO HOST ////////////////
  ////////////////////////////////////////////////////
  cudaMemcpy(U, U_d, (*m)*(*n)*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(V, V_d, (*m)*(*n)*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(Rho, Rho_d, (*m)*(*n)*sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(Fluid, Fluid_d, (*m)*(*n)*sizeof(int), cudaMemcpyDeviceToHost);

  clock_t tEnd = clock();
  float tOverall = (float)(tEnd - tStart) / CLOCKS_PER_SEC; // Calculate elapsed time

  fprintf(log_file,"\nOverall calculations took %f seconds", tOverall);
  fprintf(log_file,"\nMain while loop took %f seconds\n",      tIteration);
  fprintf(log_file,"Initialization took %f seconds\n",         tInitialization);
  fprintf(log_file,"Collision took %f seconds\n",              tCollision/1000);
  fprintf(log_file,"UpdateF took %f seconds\n",                tUpdateF/1000);
  fprintf(log_file,"Streaming took %f seconds\n",              tStreaming/1000);
  fprintf(log_file,"Calculating Boundaries took %f seconds\n", tBoundaries/1000);
  fprintf(log_file,"Update Macroscopic took %f seconds\n",     tUpdateMacro/1000);
  fprintf(log_file,"Calculating Residuals took %f seconds\n",  tResiduals);
  fprintf(log_file,"Writing results took %f seconds\n",        tWriting);

  // end time measurement, close log file
  fprintf(log_file,"\n:::: Iterations done! ::::\n");
  fclose(log_file);

  // Close residuals file
  fclose(resid_file);
  
  // Write the time measurements to a separate dat file
  TimeMeasurementFile = fopen("Results/SerialTimeMeasuerment.dat","w");
  fprintf(TimeMeasurementFile,"tOverall %f\n",        tOverall);
  fprintf(TimeMeasurementFile,"tIteration %f\n",      tIteration);
  fprintf(TimeMeasurementFile,"tInitialization %f\n", tInitialization);
  fprintf(TimeMeasurementFile,"tCollision %f\n",      tCollision/1000);
  fprintf(TimeMeasurementFile,"tUpdateF %f\n",        tUpdateF/1000);
  fprintf(TimeMeasurementFile,"tStreaming %f\n",      tStreaming/1000);
  fprintf(TimeMeasurementFile,"tBoundaries %f\n",     tBoundaries/1000);
  fprintf(TimeMeasurementFile,"tUpdateMacro %f\n",    tUpdateMacro/1000);
  fprintf(TimeMeasurementFile,"tResiduals %f\n",      tResiduals);
  fprintf(TimeMeasurementFile,"tWriting %f\n",        tWriting);
   fclose(TimeMeasurementFile);
  
  // Write final data
  switch(postproc_prog)
  {
    case 1: sprintf(FinalOutputFile, "Results/FinalData.csv"); break;
    case 2: sprintf(FinalOutputFile, "Results/FinalData.dat"); break;
  }
  WriteResults(FinalOutputFile, CoordX, CoordY, U, V, Rho, Fluid, n, m, ppp);

  // Write information for user
  printf("\n\nLog was written to %s\n", logFile);
  printf("Last autosave result can be found at %s\n", AutosaveOutputFile);
  printf("Residuals were written to Results/residuals.dat\n");
  printf("Profiling results were written to Results/SerialTimeMeasuerment.dat\n");
  printf("Final results were written to %s\n", FinalOutputFile);

  ////////////////////////////////////////////////////
  ///////////////// End of line //////////////////////
  ////////////////////////////////////////////////////

  // FREE POINTERS
  free(Delta);
  free(m);
  free(n);
  free(MaxInletCoordY);
  free(MinInletCoordY);
  free(NumInletNodes);
  free(NumNodes);
  free(NumConn);
  free(w);
  free(cx);
  free(cy);
  free(opp);
  free(Residuals);

  cudaFree (cx_d);
  cudaFree (cy_d);
  cudaFree (w_d);
  cudaFree (omega_d);
  //...

}
