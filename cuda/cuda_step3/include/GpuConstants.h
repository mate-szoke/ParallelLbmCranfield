/**
 * Header for variables stored in GPU constant memory
 * @file GpuConstants.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 *
 * @details These constants need to be extern in other files because they can only be declared in
 * one source file. If other object files need these variable, just include this header and don't
 * forget to compile them with -rdc=true
 *
 * Velocity unit vector components (#cx_d, #cy_d)
 * @verbatim
   (-1,1)   (0,1)   (1,1)
          6   2   5
            \ | /
  (-1,0) 3 -(0,0)- 1 (1,0)
            / | \
          7   4   8
   (-1,-1)  (0,-1)  (1,-1) @endverbatim
 * Lattice weights (#w_d)
 * @verbatim
    (1/36)   (1/9)   (1/36)
           6   2   5
             \ | /
    (1/9) 3 -(4/9)- 1 (1/9)
             / | \
           7   4   8
    (1/36)   (1/9)   (1/36) @endverbatim
 * Opposite lattices (#opp_d)
 * @verbatim
        (8)   (4)   (7)
           6   2   5
             \ | /
       (1) 3 -(0)- 1 (3)
             / | \
           7   4   8
        (5)   (2)   (6) @endverbatim
 */
#ifndef GPU_CONSTANTS_H
#define GPU_CONSTANTS_H

#ifndef RELEASE

extern __constant__ InletProfile inletProfile_d;   ///< inlet profile
extern __constant__ BoundaryType boundaryType_d;   ///< boundary type
extern __constant__ OutletProfile outletProfile_d; ///< outlet profile
extern __constant__ int dlBoundaryId_d;            ///< boudary ID
extern __constant__ int cx_d[9];                   ///< velocity x unit vector components
extern __constant__ int cy_d[9];                   ///< velocity y unit vector components
extern __constant__ int width_d;                   ///< number of columns
extern __constant__ int height_d;                  ///< number of rows
extern __constant__ int c_d[9];                    ///< direction offset
extern __constant__ int opp_d[9];                  ///< opposite lattice offset
extern __constant__ FLOAT_TYPE delta_d;            ///< grid spacing
extern __constant__ FLOAT_TYPE w_d[9];             ///< lattice weights
extern __constant__ FLOAT_TYPE omega_d;            ///< collision frequency for D2Q9 \f$ \omega = \frac{1}{3\nu + 0.5} \f$
extern __constant__ FLOAT_TYPE omegaA_d;           ///< assymetric collision frequency \f$ \omega_a = \frac{8(2-\omega)}{8-\omega} \f$
extern __constant__ FLOAT_TYPE rhoIn_d;            ///< input density
extern __constant__ FLOAT_TYPE uIn_d;              ///< input velocity x
extern __constant__ FLOAT_TYPE vIn_d;              ///< input velocity y
extern __constant__ FLOAT_TYPE minInletCoordY_d;   ///< maximum inlet coordinate y
extern __constant__ FLOAT_TYPE maxInletCoordY_d;   ///< minimum inlet coordinate y
extern __constant__ FLOAT_TYPE velMomMap_d[81];    ///< MRT constants: mapping between velocity and momentum space \f$ \mathbf{M} \f$
extern __constant__ FLOAT_TYPE momCollMtx_d[81];   ///< MRT constants: collision matrix in momentum space \f$ \mathbf{M}^{-1}\mathbf{S} \f$

#endif

#endif