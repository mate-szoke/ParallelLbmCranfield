#include <stdio.h>
#include "GpuFunctions.h"
#include "CellFunctions.h"
#include "ArrayUtils.h"
#include "BcMacros.h"

__constant__ InletProfile inletProfile_d;
__constant__ BoundaryType boundaryType_d;
__constant__ OutletProfile outletProfile_d;
__constant__ int dlBoundaryId_d;
__constant__ int cx_d[9];
__constant__ int cy_d[9];
__constant__ int width_d;
__constant__ int height_d;
__constant__ int c_d[9];
__constant__ int opp_d[9];
__constant__ FLOAT_TYPE delta_d;
__constant__ FLOAT_TYPE w_d[9];
__constant__ FLOAT_TYPE omega_d;
__constant__ FLOAT_TYPE omegaA_d;
__constant__ FLOAT_TYPE rhoIn_d;
__constant__ FLOAT_TYPE uIn_d;
__constant__ FLOAT_TYPE vIn_d;
__constant__ FLOAT_TYPE minInletCoordY_d;
__constant__ FLOAT_TYPE maxInletCoordY_d;
__constant__ FLOAT_TYPE velMomMap_d[81];
__constant__ FLOAT_TYPE momCollMtx_d[81];

__host__ void initConstants(Arguments *args,
                            FLOAT_TYPE maxInletCoordY, FLOAT_TYPE minInletCoordY,
                            FLOAT_TYPE delta, int m, int n)
{
    //CONSTANT LATTICE QUANTITIES
    int s = m*n;
    FLOAT_TYPE w[9] = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
    int opp[9] = { 0, 3*s,  4*s, 1*s, 2*s,    7*s,    8*s, 5*s, 6*s };
    int cx[9]  = { 0,   1,    0,  -1,   0,      1,     -1,  -1,   1 };
    int cy[9]  = { 0,   0,    1,   0,  -1,      1,      1,  -1,  -1 };
    int c[9]   = { 0,  -1, -1*n,   1,   n, -1*n-1, -1*n+1, n+1, n-1 };

    // Calculate collision freq
    FLOAT_TYPE omega  = 1.0/(3.*args->viscosity+0.5);
    FLOAT_TYPE omegaA = 8*(2-omega)/(8-omega);

    cudaMemcpyToSymbol(cx_d,  cx,  9*sizeof(int));
    cudaMemcpyToSymbol(cy_d,  cy,  9*sizeof(int));
    cudaMemcpyToSymbol(c_d,   c,   9*sizeof(int));
    cudaMemcpyToSymbol(opp_d, opp, 9*sizeof(int));

    cudaMemcpyToSymbol(outletProfile_d, &args->outletProfile, sizeof(OutletProfile));
    cudaMemcpyToSymbol(boundaryType_d,  &args->boundaryType,  sizeof(BoundaryType));
    cudaMemcpyToSymbol(dlBoundaryId_d,  &args->boundaryId,    sizeof(int));

    cudaMemcpyToSymbol(width_d,  &m,      sizeof(int));
    cudaMemcpyToSymbol(height_d, &n,      sizeof(int));
    cudaMemcpyToSymbol(w_d,      w,     9*sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(omega_d,  &omega,  sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(omegaA_d, &omegaA, sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(delta_d,  &delta,  sizeof(FLOAT_TYPE));

    cudaMemcpyToSymbol(inletProfile_d,   &args->inletProfile, sizeof(InletProfile));
    cudaMemcpyToSymbol(rhoIn_d,          &args->rho,          sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(uIn_d,            &args->u,            sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(vIn_d,            &args->v,            sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(minInletCoordY_d, &minInletCoordY,     sizeof(FLOAT_TYPE));
    cudaMemcpyToSymbol(maxInletCoordY_d, &maxInletCoordY,     sizeof(FLOAT_TYPE));

    // Initialize variables for MRT Collision model, if used
    if (args->collisionModel == MRT)
    {
        FLOAT_TYPE *velMomMap  = createHostArrayFlt(81);
        FLOAT_TYPE *momCollMtx = createHostArrayFlt(81);
        MRTInitializer(velMomMap, momCollMtx, omega);

        cudaMemcpyToSymbol(velMomMap_d, velMomMap, 81*sizeof(FLOAT_TYPE));
        cudaMemcpyToSymbol(momCollMtx_d, momCollMtx, 81*sizeof(FLOAT_TYPE));
    }
}

__global__ void gpu_init(int *corner_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *q_d, int *boundary_d,
                         FLOAT_TYPE *coordY_d, int *stream_d, int *bcId_d, FLOAT_TYPE *u_d,
                         FLOAT_TYPE *v_d, FLOAT_TYPE *u0_d, FLOAT_TYPE *v0_d)
{
    int bidx=blockIdx.x;
    int tidx=threadIdx.x;

    int ind = tidx + bidx*blockDim.x;

    int ind_s = ind - (width_d*height_d) *( (int)(ind/(width_d*height_d)) );

    int ms=width_d*height_d;

    int k;

    if (ind<(width_d*height_d))
    {
        // not in the corner
        rho_d[ind] = rhoIn_d;
        corner_d[ind]=0;
    }

    if (ind<(9*width_d*height_d))
    {
        if (q_d[ind] == 0)
        {
            q_d[ind] = 0.5;
        }
        if (bcId_d[ind]!=0) //1wall,2inlet,3outlet
        {
            if (boundary_d[ind_s]==0)// if Boundary condition in the node is 0 it becomes equal to the BC of the lattice direction
            { //wallID
                boundary_d[ind_s]=bcId_d[ind];
            }
            else
            {// if in the same node there are lattice directions with different BC (corners) the BC of the node is WALL (assuming that it's impossibe to findoutlet and inlet together)
                if (boundary_d[ind_s]<bcId_d[ind])
                {
                  boundary_d[ind_s]=1;
                  corner_d[ind_s]=1;
                  // printf("cornerD1 %d\n", ind_s);
                }
                if (boundary_d[ind_s]>bcId_d[ind])
                {
                  boundary_d[ind_s]=1;
                  corner_d[ind_s]=1;
                  // printf("cornerD2 %d\n", ind_s);
                }
            }
        }


        // BC ON CORNERS IS WALL!!!! (this operation is useful for wall condition, which checks the single direction)
        if (corner_d[ind_s]==1)
        {
            if (bcId_d[ind_s+width_d*height_d]!=0 && bcId_d[ind_s+2*width_d*height_d]!=0)
            {
                bcId_d[ind_s+5*width_d*height_d]=1;
            }
            if (bcId_d[ind_s+width_d*height_d]!=0 && bcId_d[ind_s+4*width_d*height_d]!=0)
            {
                bcId_d[ind_s+8*width_d*height_d]=1;
            }
            if (bcId_d[ind_s+2*width_d*height_d]!=0 && bcId_d[ind_s+3*width_d*height_d]!=0)
            {
                bcId_d[ind_s+6*width_d*height_d]=1;
            }
            if (bcId_d[ind_s+3*width_d*height_d]!=0 && bcId_d[ind_s+4*width_d*height_d]!=0)
            {
                bcId_d[ind_s+7*width_d*height_d]=1;
            }
        }

        // INITIALIZE STREAMING (STREAM EVERYWHERE)
        stream_d[ind] = 1;
    }

    if (ind<(9*width_d*height_d))
    {
        // DON'T STREAM FROM OUTSIDE OF THE DOMAIN
        for(k=0;k<9;k++)
        {
            if (bcId_d[ind_s+k*ms]!=0)
            {
                //if (ind<width_d*height_d) printf("STR: %d, BCID: %d(%d)\n", ind_s+opp_d[k], ind_s+k*ms, bcId_d[ind_s+k*ms]);
                stream_d[ind_s+opp_d[k]]= 0 ; //no stream on opposite of the wall BC
            }
        }
    }

    if (ind<(width_d*height_d))
    {
        // INLET VELOCITY
        switch(inletProfile_d)
        {
            case INLET:
                u0_d[ind_s] = 4*1.5*uIn_d*(coordY_d[ind_s]-minInletCoordY_d)*((maxInletCoordY_d-
                         minInletCoordY_d)-(coordY_d[ind_s]-minInletCoordY_d))/((maxInletCoordY_d-
                         minInletCoordY_d)*(maxInletCoordY_d-minInletCoordY_d));
                /// 4*1.5*Uavg* distance_from_min * distance_from_max / length^2
                v0_d[ind_s] = vIn_d;
            break;
            case NO_INLET:
                u0_d[ind_s] = uIn_d;
                v0_d[ind_s] = vIn_d;
            break;
            case PULSATILE_INLET:
                u0_d[ind_s] = 0;
                v0_d[ind_s] = 0;
            break;

        }
        u_d[ind_s] = u0_d[ind_s];
        v_d[ind_s] = v0_d[ind_s];
    }
}

__global__ void gpu_init_1(FLOAT_TYPE *rho_d, FLOAT_TYPE *u0_d, FLOAT_TYPE *v0_d,
                           FLOAT_TYPE *u_d, FLOAT_TYPE *v_d, FLOAT_TYPE *coordY_d, int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    FLOAT_TYPE inletLenghth2 = 0.;

    if (idx < size)
    {
        // rho_d[idx] = rhoIn_d;

        switch(inletProfile_d)
        {
            case INLET:
                ///NOTE 4*1.5*Uavg* distance_from_min * distance_from_max / length^2
                inletLenghth2 = (maxInletCoordY_d - minInletCoordY_d) * (maxInletCoordY_d - minInletCoordY_d);
                u0_d[idx] = 4*1.5*uIn_d * (coordY_d[idx] - minInletCoordY_d) * (maxInletCoordY_d - coordY_d[idx]) / inletLenghth2;
                v0_d[idx] = vIn_d;
            break;
            case NO_INLET:
                // u0_d[idx] = uIn_d;
                // v0_d[idx] = vIn_d;
            break;
            case PULSATILE_INLET:
                // u0_d[idx]= 0;
                // v0_d[idx] = 0;
            break;

        }
        // u_d[idx] = u0_d[idx];
        // v_d[idx] = v0_d[idx];
    }
}

__global__ void gpuInitInletProfile(FLOAT_TYPE *u0_d, FLOAT_TYPE *v0_d, FLOAT_TYPE *y_d, int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    FLOAT_TYPE inletLenghth2 = 0.;

    if (idx < size)
    {
        inletLenghth2 = (maxInletCoordY_d - minInletCoordY_d) * (maxInletCoordY_d - minInletCoordY_d);
        u0_d[idx] = 4 * 1.5 * uIn_d * (y_d[idx] - minInletCoordY_d) * (maxInletCoordY_d - y_d[idx]) / inletLenghth2;
        v0_d[idx] = vIn_d;
    }
}

__global__ void gpu_init_8(int *stream_d, FLOAT_TYPE *q_d, int size)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < size)
    {
        stream_d[i] = 1;
        q_d[i] = 0.5;
    }
}

__host__ int initBoundaryConditions(int *bcNodeIdX, int *bcNodeIdY, FLOAT_TYPE *q, int *bcBoundId, int *fluid,
                                    FLOAT_TYPE *bcX, FLOAT_TYPE *bcY, FLOAT_TYPE *nodeX, FLOAT_TYPE *nodeY,
                                    int *latticeId, int *stream, int *bcType, int *bcMask, int *bcIdx,
                                    int *mask, FLOAT_TYPE delta, int m, int n, int size)
{
    int bci; //bound_dim
    int ms = m * n;
    FLOAT_TYPE dx, dy;
    int opp[9] = {0, 3*ms, 4*ms, 1*ms, 2*ms, 7*ms, 8*ms, 5*ms, 6*ms};
    FLOAT_TYPE qLat[9]={0.,1.,1.,1.,1.,sqrt(2),sqrt(2),sqrt(2),sqrt(2)};

    for (bci=0; bci<size; ++bci)
    {
        int ind = bcNodeIdX[bci] + bcNodeIdY[bci] * n; //node_dim
        int dir = latticeId[bci];

        bcIdx[ind] = ind;
        bcMask[ind] |= BC_MASK(bcType[bci], dir); //was BC_ID
        bcMask[ind] |= BC_F(fluid[ind]);

        if ( !(bcMask[ind] & BC_B(BC_ALL)) )
        {
            bcMask[ind] |= BC_B(bcType[bci]);
        }
        else if ( (bcMask[ind] & BC_B(BC_ALL)) != BC_B(bcType[bci]) )
        {
            bcMask[ind] &= ~(BC_B(BC_ALL)); //delete previous
            bcMask[ind] |= BC_WALL_B;
            bcMask[ind] |= BC_CORNER;
            // printf("corner %d\n", ind);
        }

        bcMask[ind] |= BOUND_ID(bcBoundId[bci]); //BoundId must be between 0..255
        dx = bcX[bci]-nodeX[ind];
        dy = bcY[bci]-nodeY[ind];
        q[ind*8 + dir-1] = sqrt( dx*dx + dy*dy ) / (delta * qLat[dir]);
        //q_d[ind+(dir-1)*ms] = sqrt( dx*dx + dy*dy ) / (*delta_d * qLat_d[dir]);
        stream[ind + opp[dir] - ms]= 0; //no stream on opposite of the wall BC
        //printf("ind:%d(%d,%d) str:%d\n", ind, bcNodeIdX[bci], bcNodeIdY[bci], ind + opp[dir] - ms);

        mask[ind] = 1;
    }
    for (bci=0; bci<size; ++bci)
    {
        int ind = bcNodeIdX[bci] + bcNodeIdY[bci] * n; //node_dim
        //corners
        if (bcMask[ind] & BC_CORNER)
        {
            if (bcMask[ind] & BC_E(BC_ALL) && bcMask[ind] & BC_N(BC_ALL))
            {
                bcMask[ind] &= ~(BC_NE(BC_ALL));
                bcMask[ind] |= BC_WALL_NE;
            }
            if (bcMask[ind] & BC_W(BC_ALL) && bcMask[ind] & BC_N(BC_ALL))
            {
                bcMask[ind] &= ~(BC_NW(BC_ALL));
                bcMask[ind] |= BC_WALL_NW;
            }
            if (bcMask[ind] & BC_W(BC_ALL) && bcMask[ind] & BC_S(BC_ALL))
            {
                bcMask[ind] &= ~(BC_SW(BC_ALL));
                bcMask[ind] |= BC_WALL_SW;
            }
            if (bcMask[ind] & BC_E(BC_ALL) && bcMask[ind] & BC_S(BC_ALL))
            {
                bcMask[ind] &= ~(BC_SE(BC_ALL));
                bcMask[ind] |= BC_WALL_SE;
            }
        }
    }
    int i,num=0;
    for (i=0; i<m*n; ++i)
    {
        num += mask[i];
    }
    return num;
}

__host__ void collapseBc(int *bcIdx, int *bcIdxCollapsed_d, int *bcMask, int *bcMaskCollapsed_d,
                         FLOAT_TYPE *q, FLOAT_TYPE *qCollapsed_d, int *mask, int m, int n, int size)
{
    int *bcIdxCollapsed = (int*)malloc(size*sizeof(int));
    int *bcMaskCollapsed = (int*)malloc(size*sizeof(int));
    FLOAT_TYPE *QCollapsed = (FLOAT_TYPE*)malloc(size*8*sizeof(FLOAT_TYPE));

    int flyId = 0;
    int i,j;
    for (i=0; i<n*m; ++i)
    {
        if (mask[i])
        {
            bcIdxCollapsed[flyId] = bcIdx[i];
            bcMaskCollapsed[flyId] = bcMask[i];
            for (j=0; j<8; ++j)
            {
                QCollapsed[flyId*8+j] = q[i*8+j];
            }
            flyId++;
        }
    }

    cudaMemcpy(bcIdxCollapsed_d, bcIdxCollapsed, size*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(bcMaskCollapsed_d, bcMaskCollapsed, size*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(qCollapsed_d, QCollapsed, 8*size*sizeof(FLOAT_TYPE), cudaMemcpyHostToDevice);
    free(bcIdxCollapsed);
    free(bcMaskCollapsed);
    free(QCollapsed);
}
