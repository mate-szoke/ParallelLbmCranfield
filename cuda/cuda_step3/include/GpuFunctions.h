/**
 * Header for all the kernels
 * @file GpuFunctions.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#ifndef GPU_FUNCTIONS_H
#define GPU_FUNCTIONS_H

#include <cuda.h>
#include <stdio.h>
#include "FloatType.h"
#include "Arguments.h"

/**
 * @brief Initialise GPU constants
 * @see GpuConstant.h for the constants
 *
 * @param args input parameters
 * @param maxInletCoordY maximum inlet coordinate y
 * @param minInletCoordY minimum inlet coordinate y
 * @param delta grid spacing
 * @param m number of columns
 * @param n number of rows
 */
__host__ void initConstants(Arguments *args,
                            FLOAT_TYPE maxInletCoordY, FLOAT_TYPE minInletCoordY,
                            FLOAT_TYPE delta, int m, int n);

/**
 * @brief Initialise boundary conditions and macroscopic values
 * @deprecated use #gpu_init_1 and #initBoundaryConditions instead
 *
 * @param corner_d corner array
 * @param rho_d density array
 * @param Q_d grid center ratio array
 * @param boundary_d node boundary type
 * @param coordY_d node coordinates y
 * @param stream_d streaming array
 * @param bcId_d lattice boundary type
 * @param u_d,v_d velocity vectors (x,y)
 * @param u0_d,v0_d input velocity vectors (x,y)
 */
__global__ void gpu_init(int *corner_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *Q_d, int *boundary_d,
                         FLOAT_TYPE *coordY_d, int *stream_d, int *bcId_d, FLOAT_TYPE *u_d,
                         FLOAT_TYPE *v_d, FLOAT_TYPE *u0_d, FLOAT_TYPE *v0_d);

/**
 * @brief Initialise macroscopic values
 * @deprecated use #createGpuArrayFlt and #gpuInitInletProfile
 *
 * @param rho_d density
 * @param u0_d,v0_d input velocity
 * @param u_d,v_d velocity
 * @param coordY_d node coordinates y
 * @param size vector size (MxN)
 */
__global__ void gpu_init_1(FLOAT_TYPE *rho_d, FLOAT_TYPE *u0_d, FLOAT_TYPE *v0_d,
                           FLOAT_TYPE *u_d, FLOAT_TYPE *v_d, FLOAT_TYPE *coordY_d, int size);

/**
 * @brief Initialise initial velocity vector for inlet profile
 * @details This function only sets the x velocity (u) to the following profile
 * \f[ l = (y_{max}^{inlet} - y_{min}^{inlet})^2 \f]
 * \f[ u_0 = \frac{6}{l} * u_{in} * (y - y_{min}^{inlet}) * (y_{max}^{inlet} - y) \f]
 *
 * @param u0_d, v0_d initial velocity
 * @param y_d node coordinates y
 * @param size vector size (MxN)
 */
__global__ void gpuInitInletProfile(FLOAT_TYPE *u0_d, FLOAT_TYPE *v0_d, FLOAT_TYPE *y_d, int size);

/**
 * @brief Initialise stream and q arrays
 * @deprecated use #createGpuArrayFlt instead
 *
 * @param stream_d streaming array
 * @param Q_d grid center ratio array
 * @param size vector size (8xMxN)
 */
__global__ void gpu_init_8(int *stream_d, FLOAT_TYPE *Q_d, int size);

/**
 * @brief Initialise boundary conditions
 * @details Fills the bitmask with appropriate boundary conditions and counts the boundary
 * conditions group by the nodes. Fills the mask array for #collapseBc to collapse the MxN
 * size arrays into Nbc size arrays. Nodes treated as corners if their lattices hold more than
 * one kind of boundaries (e.g. inlet and wall). The corner node and the corner lattice than
 * treated as wall, and the corner flag is set in the bitmask.
 *
 * Streaming is enabled everywhere by default, but it is disabled on the opposite side of the
 * boundary lattices.
 *
 * The q vector contains the ratio between the grid spacing and the lattice distance from the
 * center of the node. If all boundaries are straight lines q vector contains 0.5 everywhere.
 * If the boundary doesn't go through the center of the node (curved conditions), the following
 * interpolation is applied:
 * \f[ \delta_{grid} = |x_{0,0}-x_{1,0}| \f]
 * \f[ q_{lat} = [0;1;1;1;1;\sqrt{2};\sqrt{2};\sqrt{2};\sqrt{2}] \f]
 * \f[ q = \frac{\sqrt{ (x_{BC}-x_{node})^2 + (y_{BC}-y_{node})^2 }}{\delta_{grid} q_{lat}} \f]
 * @see BcMacros.h
 *
 * @param[in]  bcNodeIdX, bcNodeIdY node ID (i,j)
 * @param[out] q grid center ratio array
 * @param[in]  bcBoundId boundary ID for drag/lift computation
 * @param[in]  fluid fluid condition
 * @param[in]  bcX,bcY BC coordinates (x,y)
 * @param[in]  nodeX,nodeY node coordinates (x,y)
 * @param[in]  latticeId lattice id (0..8) array
 * @param[out] stream streaming array
 * @param[in]  bcType boundary type
 * @param[out] bcMask boundary conditions bitmask
 * @param[out] bcIdx boundary condition indices
 * @param[out] mask mask for collapsing (needed by #collapseBc)
 * @param[in]  delta grid spacing
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 * @param[in]  size vector size (MxN)
 * @return number of boundary conditions
 */
__host__ int initBoundaryConditions(int *bcNodeIdX, int *bcNodeIdY, FLOAT_TYPE *q, int *bcBoundId, int *fluid,
                                    FLOAT_TYPE *bcX, FLOAT_TYPE *bcY, FLOAT_TYPE *nodeX, FLOAT_TYPE *nodeY,
                                    int *latticeId, int *stream, int *bcType, int *bcMask, int *bcIdx,
                                    int *mask, FLOAT_TYPE delta, int m, int n, int size);

/**
 * @brief Collapse boundary condition arrays
 *
 * @param[in]  bcIdx boundary condition indices
 * @param[out] bcIdxCollapsed_d boundary condition indeces collapsed
 * @param[in]  bcMask [description]
 * @param[out] bcMaskCollapsed_d [description]
 * @param[in]  q grid center ratio array
 * @param[out] qCollapsed_d grid center ratio array collapsed
 * @param[in]  mask mask for collapsing
 * @param[in]  m number of columns
 * @param[in]  n number of rows
 * @param[in]  size number of boundary conditions
 */
__host__ void collapseBc(int *bcIdx, int *bcIdxCollapsed_d, int *bcMask, int *bcMaskCollapsed_d,
                         FLOAT_TYPE *q, FLOAT_TYPE *qCollapsed_d, int *mask, int m, int n, int size);

/**
 * @brief Inlet boundary conditions
 * @details Computes the effect of the inlet conditions, using ...
 * North inlet conditions (with \f$ v_0 = 0 \f$):
 * \f[ f_4 = f_2 \f]
 * \f[ f_7 = f_5 - u_0(\frac{1}{6}(f_0 + f_1 + f_3) - \frac{1}{3}(f_2 + f_5 + f_6)) \f]
 * \f[ f_8 = f_6 + u_0(\frac{1}{6}(f_0 + f_1 + f_3) + \frac{1}{3}(f_2 + f_5 + f_6)) \f]
 * West inlet conditions:
 * \f[ f_1 = f_3 + \frac{2}{3} \frac{u_0}{1-u_0} (f_0 + f_2 + f_4 + 2(f_3 + f_6 + f_7)) \f]
 * \f[ f_5 = f_7 + \frac{1}{6} \frac{u_0}{1-u_0} (f_0 + f_2 + f_4 + 2(f_3 + f_6 + f_7)) \f]
 * \f[ f_8 = f_6 + \frac{1}{6} \frac{u_0}{1-u_0} (f_0 + f_2 + f_4 + 2(f_3 + f_6 + f_7)) \f]
 * The function goes through the boundary conditions and only act on lattices having inlet
 * conditions that are not corner nodes. In the corners, wall conditions apply.
 *
 * @param[in]  bcIdx_d boundary condition indices
 * @param[in]  bcMask_d boundary conditions bitmask
 * @param[out] f_d distribution function
 * @param[in]  u0_d,v0_d input velocity vectors (x,y)
 * @param[in]  size number of boundary conditions
 * @ingroup solver
 */
__global__ void gpuBcInlet(int *bcIdx_d, int *bcMask_d, FLOAT_TYPE* f_d, FLOAT_TYPE* u0_d, FLOAT_TYPE* v0_d, int size);

/**
 * @brief Wall boundary conditions
 * @details Computes the effect of the wall conditions. For straight wall the computation we use
 * the half-way bounce-back method: \f$ f_{opp(i)} = f_i \f$, e.g. on an east-side wall:
 * \f[ f_3 = f_1; f_6 = f_8; f_7 = f_5 \f]
 * If the boundaries are not straight the following interpolation is done using q vector as defined
 * in #initBoundaryConditions, where \f$ f_{c(i)} \f$ means the same lattice in the neighbouring node in the
 * given direction (more about this in #gpuStreaming):
 * \f[ f_{opp(i)} = 2q_i f_i^{coll} + (1-2q_i)f_{c(i)}^{coll}, q_i<\frac{1}{2} \f]
 * \f[ f_{opp(i)} = \frac{1}{2q_i} f_i^{coll} + \frac{2q_i-1}{2q_i} f_{opp(i)}^{coll}, q_i>=\frac{1}{2} \f]
 *
 * @param[in]  bcIdx_d boundary condition indices
 * @param[in]  bcMask_d boundary conditions bitmask
 * @param[out] f_d distribution function
 * @param[in]  fColl_d distribution function from the collison step
 * @param[in]  Q_d grid center ratio array
 * @param[in]  size number of boundary conditions
 * @ingroup solver
 */
__global__ void gpuBcWall(int *bcIdx_d, int *bcMask_d, FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d, FLOAT_TYPE *Q_d, int size);

/**
 * @brief Outlet boundary conditions
 * @details Computes the effect of the outlet conditions, using ...
 * East outlet conditions:
 * \f[ f_3 = f_1 - \frac{2}{3} \frac{u_0}{1-u_0} (f_0 + f_2 + f_4 + 2(f_1 + f_5 + f_8)) \f]
 * \f[ f_6 = f_8 - \frac{1}{6} \frac{u_0}{1-u_0} (f_0 + f_2 + f_4 + 2(f_1 + f_5 + f_8)) \f]
 * \f[ f_7 = f_5 - \frac{1}{6} \frac{u_0}{1-u_0} (f_0 + f_2 + f_4 + 2(f_1 + f_5 + f_8)) \f]
 * Second order open boundary east outlet conditions
 * \f[ f_1 = 2f_1^{(-1)} - f_1^{(-2)} \f]
 * \f[ f_5 = 2f_5^{(-1)} - f_5^{(-2)} \f]
 * \f[ f_8 = 2f_8^{(-1)} - f_8^{(-2)} \f]
 * First order open boundary east outlet conditions
 * \f[ f_1 = f_1^{(-1)} \f]
 * \f[ f_5 = f_5^{(-1)} \f]
 * \f[ f_8 = f_8^{(-1)} \f]
 *
 * @param[in]  bcIdx_d boundary condition indices
 * @param[in]  bcMask_d boundary conditions bitmask
 * @param[out] f_d distribution function
 * @param[in]  u0_d,v0_d input velocity vectors (x,y)
 * @param[in]  size number of boundary conditions
 * @ingroup solver
 */
__global__ void gpuBcOutlet(int *bcIdx_d, int *bcMask_d, FLOAT_TYPE *f_d, FLOAT_TYPE *u0_d, FLOAT_TYPE *v0_d, int size);

/**
 * @brief Convert new bitmask to the old boundary condition arrays
 * @note used only for validation
 *
 * @param[in]  bcIdx_d boundary condition indices
 * @param[in]  bcMask_d boundary conditions bitmask
 * @param[out] fluid_d fluid condition
 * @param[out] boundary_d node boundary type
 * @param[out] corner_d corner array
 * @param[out] bcId_d lattice boundary type
 * @param[in]  size number of boundary conditions
 */
__global__ void gpu_convert(int *bcIdx_d, int *bcMask_d, int *fluid_d, int *boundary_d, int *corner_d, int *bcId_d, int size);

/**
 * @brief BGKW collision model
 * @deprecated use #gpuCollBgkw
 *
 * @param[in]  fluid_d fluid condition
 * @param[out] fEq_d distribution function without frequency (unused)
 * @param[in]  rho_d density
 * @param[in]  u_d,v_d velocity vectors (x,y)
 * @param[out] fColl_d distribution function
 * @param[in]  f_d distribution function from boundary step
 */
__global__ void gpu_bgk(int* fluid_d, FLOAT_TYPE* fEq_d, FLOAT_TYPE* rho_d, FLOAT_TYPE* u_d, FLOAT_TYPE* v_d, FLOAT_TYPE* fColl_d, FLOAT_TYPE* f_d);

/**
 * @brief BGKW collision model (Bhatnagar–Gross–Krook-Welander)
 * @details Computation for every lattice :
 * \f[ f_i^{eq} = \rho w (1 + 3(u c_{xi} + v c_{yi}) + \frac{9}{2}(u c_{xi} + v c_{yi})^2 - \frac{3}{2}(u^2 + v^2)) \f]
 * \f[ f_i^{coll} = \omega f_i^{eq} + (1-\omega) f_i^{eq} \f]
 *
 * @param[in]  fluid_d fluid condition
 * @param[in]  rho_d density
 * @param[in]  u_d,v_d velocity vectors (x,y)
 * @param[in]  f_d distribution function from boundary step
 * @param[out] fColl_d distribution function
 * @ingroup solver
 */
__global__ void gpuCollBgkw(int* fluid_d, FLOAT_TYPE* rho_d, FLOAT_TYPE* u_d, FLOAT_TYPE* v_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d);

/**
 * @brief TRT collision model (two-relaxation-time)
 * @details Computation for every lattice :
 * \f[ f_i^{eq} = \rho w (1 + 3(u c_{xi} + v c_{yi}) + \frac{9}{2}(u c_{xi} + v c_{yi})^2 - \frac{3}{2}(u^2 + v^2)) \f]
 * \f[ f_i^{coll} = f_i - \frac{1}{2}\omega(f_i+f_i^{-1} - f_i^{eq}-f_i^{eq,-1}) - \frac{1}{2}\omega_a(f_i-f_i^{-1} - f_i^{eq}+f_i^{eq,-1}) \f]
 * \f[ f_i^{-1} = f_{opp(i)} \f]
 *
 * @param[in]  fluid_d fluid condition
 * @param[in]  rho_d   density
 * @param[in]  u_d,v_d velocity vectors (x,y)
 * @param[in]  f_d     distribution function from boundary step
 * @param[out] fColl_d distribution function
 * @ingroup solver
 */
__global__ void gpuCollTrt(int *fluid_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *u_d, FLOAT_TYPE *v_d, FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d);

/**
 * @brief TRT collision model (step 1)
 * @deprecated use #gpuCollTrt
 *
 * @param[in]  fluid_d fluid condition
 * @param[out] fEq_d distribution function without frequency
 * @param[in]  rho_d density
 * @param[in]  u_d,v_d velocity vectors (x,y)
 */
__global__ void gpu_trt1(int* fluid_d, FLOAT_TYPE* fEq_d, FLOAT_TYPE* rho_d, FLOAT_TYPE* u_d, FLOAT_TYPE* v_d);

/**
 * @brief TRT collision model (step 2))
 * @deprecated use #gpuCollTrt
 *
 * @param[in]  fluid_d fluid condition
 * @param[in]  fEq_d distribution function without frequency
 * @param[out] fColl_d distribution function
 * @param[in]  f_d distribution function from boundary step
 */
__global__ void gpu_trt2(int* fluid_d, FLOAT_TYPE* fEq_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d);

/**
 * @brief MRT collision model (multiple-relaxation-time)
 * @details The MRT collision model computes new values for the distribution function in momentum
 * space and uses the following formula:
 * \f[ f_i^{coll} = f_i - \mathbf{M}^{-1}\mathbf{S}(\mathbf{m}-\mathbf{m}^{eq}) = f_i - \mathbf{M}^{-1}\mathbf{S}(\mathbf{M}\mathbf{f}-\mathbf{m}^{eq}) \f]
 * \f[ mathbf{m} = (\rho, e, \epsilon, j_x, q_x, j_y, q_y, p_{xx}, p_{xy})^T \f]
 * \f[ m_0^{eq} = \rho \f]
 * \f[ m_1^{eq} = -2\rho + 3(j_x^2+j_y^2) \f]
 * \f[ m_2^{eq} = \rho - 3(j_x^2+j_y^2) \f]
 * \f[ m_3^{eq} = j_x = \rho u_x = \sum_i f_i^{eq}c_{ix} \f]
 * \f[ m_4^{eq} = -j_x \f]
 * \f[ m_5^{eq} = j_y = \rho u_y = \sum_i f_i^{eq}c_{iy} \f]
 * \f[ m_6^{eq} = -j_y \f]
 * \f[ m_7^{eq} = j_x^2+j_y^2 \f]
 * \f[ m_8^{eq} = j_x j_y \f]
 * For more details check A. A. Mohamad - Lattice Boltzmann Method (book)
 *
 * @param[in]  fluid_d fluid condition
 * @param[in]  rho_d density
 * @param[in]  u_d,v_d velocity vectors (x,y)
 * @param[in]  f_d distribution function from boundary step
 * @param[out] fColl_d distribution function
 * @ingroup solver
 */
__global__ void gpuCollMrt(int* fluid_d, FLOAT_TYPE *rho_d, FLOAT_TYPE *u_d, FLOAT_TYPE *v_d, FLOAT_TYPE *f_d, FLOAT_TYPE *fColl_d);

/**
 * @brief MRT collision model (step 1)
 * @deprecated use #gpuCollMrt
 *
 * @param[in]  fluid_d fluid condition
 * @param[in]  rho_d density
 * @param[in]  u_d,v_d velocity vectors (x,y)
 * @param[in]  f_d distribution function from boundary step
 * @param[out] m_d,mEq_d intermediary arrays
 */
__global__ void gpu_mrt1(int* fluid_d, FLOAT_TYPE* rho_d, FLOAT_TYPE* u_d, FLOAT_TYPE* v_d, FLOAT_TYPE* f_d, FLOAT_TYPE* mEq_d, FLOAT_TYPE* m_d);

/**
 * @brief MRT collision model (step 2)
 * @deprecated use #gpuCollMrt
 *
 * @param[in]  fluid_d fluid condition
 * @param[out] fColl_d distribution function
 * @param[in]  f_d distribution function from boundary step
 * @param[in]  m_d,mEq_d intermediary arrays
 * @param[out] collision_d intermediary arrays (unused)
 */
__global__ void gpu_mrt2(int* fluid_d, FLOAT_TYPE* collision_d, FLOAT_TYPE* m_d, FLOAT_TYPE* mEq_d, FLOAT_TYPE* fColl_d, FLOAT_TYPE* f_d);

/**
 * @brief Streaming
 * @deprecated use #gpuStreaming
 *
 * @param fluid_d fluid condition
 * @param stream_d streaming array
 * @param f_d distribution function
 * @param fColl_d distribution function from the collison step
 */
__global__ void gpu_streaming(int* fluid_d, int* stream_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d);

/**
 * @brief Streaming
 * @details  Basically the distribution function propagates to neighbouring cells the following way
 * |     | i-4 | i-3 | i-2 | i-1 |  i  | i+1 | i+2 | i+3 | i+4 |
 * |:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
 * | j+4 | f6  |     |     |     | f2  |     |     |     | f5  |
 * | j+3 |     |  \\ |     |     | \|  |     |     |  /  |     |
 * | j+2 |     |     |  \\ |     | \|  |     |  /  |     |     |
 * | j+1 |     |     |     | f6  | f2  | f5  |     |     |     |
 * | j   | f3  |  -  |  -  | f3  | f0  | f1  |  -  |  -  | f1  |
 * | j-1 |     |     |     | f7  | f4  | f8  |     |     |     |
 * | j-2 |     |     |  /  |     | \|  |     |  \\ |     |     |
 * | j-3 |     |  /  |     |     | \|  |     |     |  \\ |     |
 * | j-4 | f7  |     |     |     | f4  |     |     |     | f8  |
 *
 * @param fluid_d fluid condition
 * @param stream_d streaming array
 * @param f_d distribution function
 * @param fColl_d distribution function from the collison step
 * @ingroup solver
 */
__global__ void gpuStreaming(int* fluid_d, int* stream_d, FLOAT_TYPE* f_d, FLOAT_TYPE* fColl_d);

/**
 * @brief Inlet boundary conditions
 * @deprecated use #gpuBcInlet
 *
 * @param[in]  fluid_d fluid condition
 * @param[in]  boundary_d node boundary type
 * @param[in]  bcId_d lattice boundary type
 * @param[out] f_d distribution function
 * @param[in]  u0_d,v0_d input velocity vectors (x,y)
 * @param[in]  corner_d corner array
 */
__global__ void gpu_boundaries1(int* fluid_d, int* boundary_d, int* bcId_d, FLOAT_TYPE* f_d, FLOAT_TYPE* u0_d, FLOAT_TYPE* v0_d, int* corner_d);

/**
 * @brief Wall boundary conditions
 * @deprecated use #gpuBcWall
 *
 * @param[in]  fluid_d fluid condition
 * @param[out] Fneighbours_d intermediary array for curved wall (unused)
 * @param[in]  bcId_d lattice boundary type
 * @param[out] f_d distribution function
 * @param[in]  Q_d grid center ratio array
 * @param[in]  fColl_d distribution function from the collison step
 */
__global__ void gpu_boundaries2(int* fluid_d, FLOAT_TYPE* Fneighbours_d, FLOAT_TYPE* fColl_d, int* bcId_d, FLOAT_TYPE* Q_d, FLOAT_TYPE* f_d);

/**
 * @brief Outlet boundary conditions
 * @deprecated use #gpuBcOutlet
 *
 * @param[in]  fluid_d fluid condition
 * @param[in]  bcId_d lattice boundary type
 * @param[out] f_d distribution function
 * @param[in]  u0_d,v0_d input velocity vectors (x,y)
 * @param[in]  corner_d corner array
 */
__global__ void gpu_boundaries3(int* fluid_d, int* bcId_d, FLOAT_TYPE* f_d, FLOAT_TYPE* u0_d, FLOAT_TYPE* v0_d, int* corner_d);

/**
 * @brief Update macroscopic values
 * @deprecated use #gpuUpdateMacro
 *
 * @param[in]  fluid_d fluid condition
 * @param[out] rho_d density
 * @param[out] u_d,v_d velocity vectors (x,y)
 * @param[in]  bcId_d lattice boundary type
 * @param[in]  boundaryId_d boundary ID for lift/drag
 * @param[out] drag_d drag
 * @param[out] lift_d lift
 * @param[in]  coordX_d,coordY_d node coordinates x,y
 * @param[in]  f_d distribution function
 */
__global__ void gpu_update_macro(int* fluid_d, FLOAT_TYPE* rho_d, FLOAT_TYPE* u_d, FLOAT_TYPE* v_d,
                                int* bcId_d, int* boundaryId_d, FLOAT_TYPE* drag_d,
                                FLOAT_TYPE* lift_d, FLOAT_TYPE* coordX_d, FLOAT_TYPE* coordY_d, FLOAT_TYPE* f_d);

/**
 * @brief Update macroscopic values using bitmask
 * @deprecated use #gpuUpdateMacro
 *
 * @param[in]  fluid_d fluid condition
 * @param[out] rho_d density
 * @param[out] u_d,v_d velocity vectors (x,y)
 * @param[in]  bcMask_d boundary conditions bitmask
 * @param[out] drag_d drag
 * @param[out] lift_d lift
 * @param[in]  coordX_d,coordY_d node coordinates x,y
 * @param[in]  f_d distribution function
 */
__global__ void gpu_update_new(int *fluid_d, FLOAT_TYPE* rho_d, FLOAT_TYPE* u_d, FLOAT_TYPE* v_d,
                               int *bcMask_d, FLOAT_TYPE* drag_d, FLOAT_TYPE* lift_d,
                               FLOAT_TYPE* coordX_d, FLOAT_TYPE* coordY_d, FLOAT_TYPE* f_d);

/**
 * @brief Update macroscopic values using bitmask
 * @details Computation of the macroscopic values for a node are the following:
 * \f[ \rho = \sum_{i=0}^8 f_i \f]
 * \f[ \vec{v} = (u,v) = \frac{1}{\rho} \sum_{i=1}^8 f_i \vec{c_i} \f]
 * \f[ F_{drag} = \frac{\rho}{15} (20-x) \f]
 * \f[ F_{lift} = \frac{\rho}{15} (20-y) \f]
 *
 * @see for \f$ \vec{c_i} \f$ see #cx_d and #cy_d in GpuConstants.h
 *
 * @param[in]  fluid_d fluid condition
 * @param[out] rho_d density
 * @param[out] u_d,v_d velocity vectors (x,y)
 * @param[in]  bcMask_d boundary conditions bitmask (see BcMacros.h)
 * @param[out] drag_d drag
 * @param[out] lift_d lift
 * @param[in]  coordX_d,coordY_d node coordinates x,y
 * @param[in]  f_d distribution function
 * @ingroup solver
 */
__global__ void gpuUpdateMacro(int *fluid_d, FLOAT_TYPE* rho_d, FLOAT_TYPE* u_d, FLOAT_TYPE* v_d,
                               int *bcMask_d, FLOAT_TYPE* drag_d, FLOAT_TYPE* lift_d,
                               FLOAT_TYPE* coordX_d, FLOAT_TYPE* coordY_d, FLOAT_TYPE* f_d);

#endif
