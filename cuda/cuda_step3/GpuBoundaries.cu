#include "GpuFunctions.h"
#include "BcMacros.h"
#include "GpuConstants.h"

__global__ void gpuBcInlet(int *bcIdx_d, int *bcMask_d, FLOAT_TYPE* f_d, FLOAT_TYPE* u0_d, FLOAT_TYPE* v0_d, int size)
{
  int bci = blockIdx.x * blockDim.x + threadIdx.x;
  int ms=width_d*height_d;

  if (bci < size)
  {
    int ind = bcIdx_d[bci];
    if (bcMask_d[bci] & BC_FLUID)
    {
      if (bcMask_d[bci] & BC_INLT_B)
      {
        if (bcMask_d[bci] & BC_INLT_N && !(bcMask_d[bci] & BC_CORNER))
        {
          ///@todo code: simplify code even though it will change accuracy
          f_d[ind+4*ms] = f_d[ind+2*ms];

          f_d[ind+8*ms] = f_d[ind+6*ms]
                          + ( f_d[ind] + f_d[ind+1*ms]
                            + f_d[ind+3*ms]
                            + 2*( f_d[ind+2*ms]
                                + f_d[ind+6*ms]
                                + f_d[ind+5*ms] )
                            ) * u0_d[ind]/6.0;

          f_d[ind+7*ms] = f_d[ind+5*ms]
                          - ( f_d[ind] + f_d[ind+1*ms]
                              + f_d[ind+3*ms]
                              + 2*( f_d[ind+2*ms]
                                  + f_d[ind+6*ms]
                                  + f_d[ind+5*ms])
                            ) * u0_d[ind]/6.0;
          //printf("BCF%d 4:%f, 7:%f, 8:%f\n", ind, f_d[ind+4*ms], f_d[ind+7*ms], f_d[ind+8*ms]);
        }
        if (bcMask_d[bci] & BC_INLT_W && !(bcMask_d[bci] & BC_CORNER))
        {
          f_d[ind+1*ms] = f_d[ind+3*ms]
                        + 2.* ( ( f_d[ind]
                                + f_d[ind+2*ms]
                                + f_d[ind+4*ms]
                                + 2.*(f_d[ind+3*ms]
                                    + f_d[ind+6*ms]
                                    + f_d[ind+7*ms])
                                ) / (1.0-u0_d[ind])
                              )*u0_d[ind]/3;

          f_d[ind+5*ms] = f_d[ind+7*ms]
                        + ( ( f_d[ind]
                            + f_d[ind+2*ms]
                            + f_d[ind+4*ms]
                            + 2.*(f_d[ind+3*ms]
                                + f_d[ind+6*ms]
                                + f_d[ind+7*ms])
                            ) / (1.0-u0_d[ind])
                          )*u0_d[ind]/6;

          f_d[ind+8*ms] = f_d[ind+6*ms]
                        + ( ( f_d[ind]
                            + f_d[ind+2*ms]
                            + f_d[ind+4*ms]
                            + 2.*(f_d[ind+3*ms]
                                + f_d[ind+6*ms]
                                + f_d[ind+7*ms])
                            ) / (1.0-u0_d[ind])
                          )*u0_d[ind]/6;
        }
        ///@todo code: compute inlet on other sides
      }
    }
  }
}

__global__ void gpuBcWall(int *bcIdx_d, int *bcMask_d, FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d, FLOAT_TYPE *q_d, int size)
{
  int bci = blockIdx.x * blockDim.x + threadIdx.x;
  int ms=width_d*height_d;
  int dir;

  if (bci < size)
  {
    int ind = bcIdx_d[bci];
    if (bcMask_d[bci] & BC_FLUID)
    {
      for (dir=1; dir<9; ++dir)
      {
        if ((bcMask_d[bci] & BC_MASK(BC_ALL, dir)) == (bcMask_d[bci] & BC_MASK(BC_WALL, dir)) && (bcMask_d[bci] & BC_MASK(BC_WALL, dir)))
        {
          // printf("%d: %X-%X-%X\n", ind, bcMask_d[bci], bcMask_d[bci] & BC_MASK(BC_ALL, dir), bcMask_d[bci] & BC_MASK(BC_WALL, dir));
          switch(boundaryType_d)
          {
            case CURVED: //curved
              if (q_d[bci*8+dir-1] < 0.5) // if the distance from the boundary is less than 0.5?
              {
                f_d[ind+opp_d[dir]] = 2*q_d[bci*8+dir] * fColl_d[ind+dir*ms] +
                                      (1 - 2*q_d[bci*8+dir-1]) * fColl_d[ind+dir*ms+c_d[dir]];

                // printf("if %d: dir: %d, Q:%f, F:%f; ", ind, dir, q_d[bci*8+dir], f_d[ind+opp_d[dir]]);
              }
              else
              {
                f_d[ind+opp_d[dir]] = fColl_d[ind+dir*ms]/2/q_d[bci*8+dir-1] + (2*q_d[bci*8+dir-1]-1)/(2*q_d[bci*8+dir-1])*fColl_d[ind+opp_d[dir]];

                // printf("else %d: dir: %d, Q:%f,  F: %f; ", ind, dir, q_d[bci*8+dir], f_d[ind+opp_d[dir]]);
              }
            break;
            case STRAIGHT: //half-way bounce back
              f_d[ind+opp_d[dir]] = f_d[ind+dir*ms];
            break;
          }
        }
      }
    }
  }
}

__global__ void gpuBcOutlet(int *bcIdx_d, int *bcMask_d, FLOAT_TYPE *f_d, FLOAT_TYPE *u0_d, FLOAT_TYPE *v0_d, int size)
{
  int bci = blockIdx.x * blockDim.x + threadIdx.x;
  int ms=width_d*height_d;

  if (bci < size)
  {
    int ind = bcIdx_d[bci];
    if (bcMask_d[bci] & BC_FLUID)
    {
      switch (outletProfile_d)
      {
        case OUTLET:
          if ((bcMask_d[bci] & BC_OUTL_E) == BC_OUTL_E)
          {
            ///@todo code: simplify code even though it will change accuracy
            f_d[ind+3*ms] = f_d[ind+1*ms]
                          - 2*( ( f_d[ind]
                                + f_d[ind+2*ms]
                                + f_d[ind+4*ms]
                                + 2*( f_d[ind+1*ms]
                                    + f_d[ind+5*ms]
                                    + f_d[ind+8*ms])
                                )/(1-u0_d[ind])
                              )*u0_d[ind]/3;

            f_d[ind+7*ms] = f_d[ind+5*ms]
                          - ( ( f_d[ind]
                              + f_d[ind+2*ms]
                              + f_d[ind+4*ms]
                              + 2*( f_d[ind+1*ms]
                                  + f_d[ind+5*ms]
                                  + f_d[ind+8*ms])
                              )/(1-u0_d[ind])
                            )*u0_d[ind]/6;

            f_d[ind+6*ms] = f_d[ind+8*ms]
                          - ( ( f_d[ind]
                              + f_d[ind+2*ms]
                              + f_d[ind+4*ms]
                              + 2*( f_d[ind+1*ms]
                                  + f_d[ind+5*ms]
                                  + f_d[ind+8*ms])
                              )/(1-u0_d[ind])
                            )*u0_d[ind]/6;
          }
          if (bcMask_d[bci] & BC_OUTL_N)
          {
            ///@todo code: fill north-side outlet
          }
          if (bcMask_d[bci] & BC_OUTL_W)
          {
            ///@todo code: fill west-side outlet
          }
          if (bcMask_d[bci] & BC_OUTL_S)
          {
            ///@todo code: fill south-side outlet
          }
        break;
        case OUTLET_SECOND: //open boundary
          if ((bcMask_d[bci] & BC_OUTL_E) == BC_OUTL_E)
          {
            f_d[ind+  ms] = 2*f_d[ind+  ms-1] - f_d[ind+  ms-2];
            f_d[ind+5*ms] = 2*f_d[ind+5*ms-1] - f_d[ind+5*ms-2];
            f_d[ind+8*ms] = 2*f_d[ind+8*ms-1] - f_d[ind+8*ms-2];
            // if (ind==514) printf("%d: X:%X 1:%f, 5:%f, 8:%f\n", ind, BC_OUTL_E, f_d[ind+1*ms], f_d[ind+5*ms], f_d[ind+8*ms]);
          }
          if (bcMask_d[bci] & BC_OUTL_N)
          {
            ///@todo code: fill north-side outlet
          }
          if (bcMask_d[bci] & BC_OUTL_W)
          {
            ///@todo code: fill west-side outlet
          }
          if (bcMask_d[bci] & BC_OUTL_S)
          {
            ///@todo code: fill south-side outlet
          }
        break;
        case OUTLET_FIRST: //first order
          if ((bcMask_d[bci] & BC_OUTL_E) == BC_OUTL_E)
          {
            f_d[ind+  ms] = f_d[ind-1+  ms];
            f_d[ind+5*ms] = f_d[ind-1+5*ms];
            f_d[ind+8*ms] = f_d[ind-1+8*ms];
          }
          if (bcMask_d[bci] & BC_OUTL_N)
          {
            ///@todo code: fill north-side outlet
          }
          if (bcMask_d[bci] & BC_OUTL_W)
          {
            ///@todo code: fill west-side outlet
          }
          if (bcMask_d[bci] & BC_OUTL_S)
          {
            ///@todo code: fill south-side outlet
          }
        break;
      }
    }
  }
}

__global__ void gpu_convert(int *bcIdx_d, int *bcMask_d, int *fluid_d, int *boundary_d, int *corner_d, int *bcId_d, int size)
{
  int bci = blockIdx.x * blockDim.x + threadIdx.x;
  int ms=width_d*height_d;
  int dir;

  if (bci < size)
  {
    int ind = bcIdx_d[bci];
    // printf("X-%d ", ind);
    fluid_d[ind] = (bcMask_d[bci] & BC_FLUID) ? 1 : 0;
    boundary_d[ind] = (bcMask_d[bci] & BC_B(BC_ALL))>>16;
    corner_d[ind] = (bcMask_d[bci] & BC_CORNER) ? 1 : 0;
    for (dir=1; dir<9; ++dir)
    {
      // if (bcMask_d[bci] & BC_MASK(BC_ALL, dir)) printf("(%d)",ind);
      bcId_d[ind + dir*ms] = (bcMask_d[bci] & BC_MASK(BC_ALL, dir)) >> ((dir-1)*2);
    }
  }
}

__global__ void gpu_boundaries1(int* fluid_d, int* boundary_d, int* bcId_d, FLOAT_TYPE* f_d,
                                FLOAT_TYPE* u0_d, FLOAT_TYPE* v0_d, int* corner_d)
{
    int bidx=blockIdx.x;
    int tidx=threadIdx.x;

    int ind = tidx + bidx*blockDim.x;

    int ind_s = ind - (width_d*height_d) *( (int)(ind/(width_d*height_d)) );

    int ms=width_d*height_d;

    if (ind<(9*width_d*height_d))
    {
        if (fluid_d[ind_s]==1)
        {
            //inlet BC
            if (boundary_d[ind_s]==2) // if inlet boundary
            {
                // inlet on the top (to be completed with the y velocity)
                ///@todo code: complete with y velocity
                if (bcId_d[ind_s+2*ms]==2 && corner_d[ind_s]!=1)
                {
                    // printf("BCFBEF%d 4:%f, 7:%f, 8:%f\n", ind_s, f_d[ind_s+4*ms], f_d[ind_s+7*ms], f_d[ind_s+8*ms]);
                    f_d[ind_s+4*ms] = f_d[ind_s+2*ms];

                    f_d[ind_s+8*ms] = f_d[ind_s+6*ms]
                                    + ( f_d[ind_s] + f_d[ind_s+1*ms]
                                        + f_d[ind_s+3*ms]
                                        + 2*(f_d[ind_s+2*ms]
                                        + f_d[ind_s+6*ms]
                                        + f_d[ind_s+5*ms]))*u0_d[ind_s]/6.0;

                    f_d[ind_s+7*ms] = f_d[ind_s+5*ms]
                                    - ( f_d[ind_s] + f_d[ind_s+1*ms]
                                        + f_d[ind_s+3*ms]
                                        + 2*(f_d[ind_s+2*ms]
                                        + f_d[ind_s+6*ms]
                                        + f_d[ind_s+5*ms]))*u0_d[ind_s]/6.0;
                    // printf("BCFAFT%d 4:%f, 7:%f, 8:%f\n", ind_s, f_d[ind_s+4*ms], f_d[ind_s+7*ms], f_d[ind_s+8*ms]);
                }

                // inlet on the left (to be completed with the x velocity)
                if (bcId_d[ind_s+3*ms]==2 && corner_d[ind_s]!=1)
                {
                    f_d[ind_s+1*ms] = f_d[ind_s+3*ms]
                                      + 2*((f_d[ind_s]+f_d[ind_s+2*ms]
                                          + f_d[ind_s+4*ms]
                                          + 2.*(f_d[ind_s+3*ms]
                                          + f_d[ind_s+6*ms]
                                          + f_d[ind_s+7*ms]))
                                          / (1.0-u0_d[ind_s]))*u0_d[ind_s]/3;

                    f_d[ind_s+5*ms] = f_d[ind_s+7*ms]
                                      +  (( f_d[ind_s]+f_d[ind_s+2*ms]
                                          + f_d[ind_s+4*ms]
                                          + 2.*(f_d[ind_s+3*ms]
                                          + f_d[ind_s+6*ms]
                                          + f_d[ind_s+7*ms]))
                                          / (1.0-u0_d[ind_s]))*u0_d[ind_s]/6;

                    f_d[ind_s+8*ms] = f_d[ind_s+6*ms]
                                      + ((f_d[ind_s]+ f_d[ind_s+2*ms]
                                          + f_d[ind_s+4*ms]
                                          + 2.*(f_d[ind_s+3*ms]
                                          + f_d[ind_s+6*ms]
                                          + f_d[ind_s+7*ms]))
                                          / (1.0-u0_d[ind_s]))*u0_d[ind_s]/6;
                }
                ///@todo code: compute inlet on other sides
            } //if inlet
        } //if fluid
    } //if index
}

__global__ void gpu_boundaries2(int* fluid_d, FLOAT_TYPE* fNeighbours_d, FLOAT_TYPE* fColl_d, int* bcId_d, FLOAT_TYPE* q_d, FLOAT_TYPE* f_d)
{

    int bidx=blockIdx.x;
    int tidx=threadIdx.x;

    int ind = tidx + bidx*blockDim.x;

    int ind_s = ind - (width_d*height_d) *( (int)(ind/(width_d*height_d)) );

    int ind_c = (int)(ind/(width_d*height_d));

    if (ind<(9*width_d*height_d))
    {
        if (fluid_d[ind_s]==1)
        {
            // WALL
            switch(boundaryType_d)
            {
                // curved boundaries
                case CURVED:
                    fNeighbours_d[ind]=fColl_d[ind+c_d[ind_c]];

                    if (bcId_d[ind]==1) // if wall
                    {
                        if (q_d[ind]<0.5) // if the distance from the boundary is less than 0.5?
                        {
                            f_d[ind_s+opp_d[ind_c]] = 2*q_d[ind]*fColl_d[ind]
                                                    +(1-2*q_d[ind])*fNeighbours_d[ind];
                        }
                        else
                        {
                            f_d[ind_s+opp_d[ind_c]] = fColl_d[ind]/2/q_d[ind]
                                                    + (2*q_d[ind]-1)/(2*q_d[ind])
                                                    * fColl_d[ind_s+opp_d[ind_c]];
                        }
                    }
                break;

                // bounceback boundaries
                case STRAIGHT:
                    //WALL BC (half-way bounceback)
                    if (bcId_d[ind]==1) // if wall boundary
                    {
                        f_d[ind_s+opp_d[ind_c]] = f_d[ind];
                        // printf("ind-%d: %f\n", ind, f_d[ind_s+opp_d[ind_c]]);
                    }
                break;
            } //switch boundary
        } //if fluid
    } //if index
}

__global__ void gpu_boundaries3(int* fluid_d, int* bcId_d, FLOAT_TYPE* f_d, FLOAT_TYPE* u0_d, FLOAT_TYPE* v0_d, int* corner_d)
{

    int bidx=blockIdx.x;
    int tidx=threadIdx.x;

    int ind = tidx + bidx*blockDim.x;

    int ind_s = ind - (width_d*height_d) *( (int)(ind/(width_d*height_d)) );

    int ms=width_d*height_d;

    if (ind<(width_d*height_d))
    {
        if (fluid_d[ind_s]==1)
        {
          // printf("F%d ", ind_s);
            // OUTLET
            switch(outletProfile_d)
            {
                // set profile in outlet
                case OUTLET:
                    //outlet BC
          // printf("O1%d ", ind_s);
                    if (bcId_d[ind_s+ms]==3) // if outlet boundary on the right side of the domain
                    {
                        f_d[ind_s+3*ms] = f_d[ind_s+1*ms]-2*((f_d[ind_s]
                                        + f_d[ind_s+2*ms]+f_d[ind_s+4*ms]
                                        + 2.*(f_d[ind_s+1*ms]+f_d[ind_s+5*ms]
                                        + f_d[ind_s+8*ms]))/(1-u0_d[ind_s]))*u0_d[ind_s]/3;

                        f_d[ind_s+7*ms] = f_d[ind_s+5*ms]-((f_d[ind_s]
                                        + f_d[ind_s+2*ms]+f_d[ind_s+4*ms]
                                        + 2.*(f_d[ind_s+1*ms]+f_d[ind_s+5*ms]
                                        + f_d[ind_s+8*ms]))/(1-u0_d[ind_s]))*u0_d[ind_s]/6;

                        f_d[ind_s+6*ms] = f_d[ind_s+8*ms]-((f_d[ind_s]
                                        + f_d[ind_s+2*ms]+f_d[ind_s+4*ms]
                                        + 2.*(f_d[ind_s+1*ms]+f_d[ind_s+5*ms]
                                        + f_d[ind_s+8*ms]))/(1-u0_d[ind_s]))*u0_d[ind_s]/6;
                        // printf("BCF%d 4:%f, 7:%f, 8:%f\n", ind_s, f_d[ind_s+3*ms], f_d[ind_s+7*ms], f_d[ind_s+6*ms]);
                    }
                    if (bcId_d[ind_s+2*ms]==3 && corner_d[ind_s]!=1)
                    {
                        ///@todo code: fill north-side outlet
                    }
                    if (bcId_d[ind_s+3*ms]==3 && corner_d[ind_s]!=1)
                    {
                        ///@todo code: fill west-side outlet
                    }
                    if (bcId_d[ind_s+4*ms]==3 && corner_d[ind_s]!=1)
                    {
                        ///@todo code: fill south-side outlet
                    }
                break;

                // OPEN BOUNDARY
                case OUTLET_SECOND :
          // printf("O2_%d_%d ", ind_s, bcId_d[ind_s+ms]);
                    if (bcId_d[ind_s+ms]==3)
                    {
                        f_d[ind_s+  ms] = 2*f_d[ind_s+  ms-1]-f_d[ind_s+  ms-2];
                        f_d[ind_s+5*ms] = 2*f_d[ind_s+5*ms-1]-f_d[ind_s+5*ms-2];
                        f_d[ind_s+8*ms] = 2*f_d[ind_s+8*ms-1]-f_d[ind_s+8*ms-2];
                        // if (ind_s==514) printf("BCF%d 1:%f, 5:%f, 8:%f\n", ind_s, f_d[ind_s+1*ms], f_d[ind_s+5*ms], f_d[ind_s+8*ms]);
                    }
                    if (bcId_d[ind+2*ms]==3 && corner_d[ind_s]!=1)
                    {
                        ///@todo code: fill north-side outlet
                    }
                    if (bcId_d[ind+3*ms]==3 && corner_d[ind_s]!=1)
                    {
                        ///@todo code: fill west-side outlet
                    }
                    if (bcId_d[ind+4*ms]==3 && corner_d[ind_s]!=1)
                    {
                        ///@todo code: fill south-side outlet
                    }
                break;
                // first order outlet
                case OUTLET_FIRST :
          // printf("O3%d ", ind_s);
                    if (bcId_d[ind_s+ms]==3)
                    {
                        f_d[ind_s+ms] = f_d[ind_s-1+ms];
                        f_d[ind_s+5*ms] = f_d[ind_s-1+5*ms];
                        f_d[ind_s+8*ms] = f_d[ind_s-1+8*ms];
                    }
                    if (bcId_d[ind+2*ms]==3 && corner_d[ind_s]!=1)
                    {
                        ///@todo code: fill north-side outlet
                    }
                    if (bcId_d[ind+3*ms]==3 && corner_d[ind_s]!=1)
                    {
                        ///@todo code: fill west-side outlet
                    }
                    if (bcId_d[ind+4*ms]==3 && corner_d[ind_s]!=1)
                    {
                        ///@todo code: fill south-side outlet
                    }
                break;
            } //end switch oulet
        } //if fluid
    } //if index
}
