/**
 * The solver itself
 * @file Iterate.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 * @ingroup solver
 *
 * @details Most of the time the same variable name was used for the same data. Here is a list for
 * most of them. The suffix _d is used to indicate GPU array. The old convention is the following:
 *
 * name     | size  | description                   | available values   | used in step
 * ---------|-------|-------------------------------|--------------------|--------------------------
 * corner   | MxN   | corner boundary condition     | 0: none, 1: corner | boundaries
 * rho      | MxN   | density                       | any (float)        | collision, macro
 * q        | 9xMxN | lattice length / grid spacing | any (float)        | boundaries
 * stream   | 9xMxN | streaming indicator           | 0: no, 1: stream   | streaming
 * boundary | MxN   | node boundary type            | 0: none, 1: wall, 2: inlet, 3: outlet | boundaries
 * bc_id    | 9xMxN | lattice boundary type         | 0: none, 1: wall, 2: inlet, 3: outlet | boundaries, macro
 * u,v      | MxN   | velocity                      | any (float)        | collision, macro
 * u0,v0    | MxN   | input velocity                | any (float)        | boundaries
 * coord x,y| MxN   | node coordinates              | any (float)        | macro
 * fluid    | MxN   | node fluid condition          | 0: solid, 1: fluid | collision, stream, boundary, macro
 * f        | 9xMxN | distribution function         | any (float)        | collision, stream, boundary, macro
 * meta-f   | 9xMxN | temporary distribution func   | any (float)        | collision, stream, boundary, macro
 *
 * Later the arrays containing boundary conditions are combined into a bitmask. Also the q and
 * stream arrays were reduced in size beacuse there is no point storing zero information for the
 * node itself (0) only the lattices (1..8). Changes in new convention:
 *
 * name     | size  | description                   | available values   | used in step
 * ---------|-------|-------------------------------|--------------------|--------------------------
 * bc-mask  | Nbc   | BC bitmask (BcMacros.h)       | bitmask            | boundaries, macro
 * bc-idx   | Nbc   | BC node ID                    | 0..MxN             | boundaries
 * q        | 8xMxN | lattice length / grid spacing | any (float)        | boundaries
 * stream   | 8xMxN | streaming indicator           | 0: no, 1: stream   | streaming
 *
 * The iteration takes the following steps:
 *  0. initialisation (run once)
 *  1. collision model
 *  2. streaming
 *  3. boundary conditions
 *  4. macroscopic values
 */
/**
 * @defgroup solver Solver
 * @brief CUDA LBM Solver functions
 * @details Steps:
 *  1. collision model
 *  2. streaming
 *  3. boundary conditions
 *  4. macroscopic values
 */
#ifndef ITERATE_H
#define ITERATE_H

#include "Arguments.h"

 /**
 *  @brief LBM solver
 *
 *  @param [in] inFn filenames
 *  @param [in] args input parameters
 *  @return 0 if everything went okay, 1 if iteration diverged, 2 if file read error
 */
int Iteration(InputFilenames *inFn, Arguments *args);

#endif
