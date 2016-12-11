/**
 * This file contains initialisation for the MRT collision model
 * @file CellFunctions.h
 * @author Istvan Tamas Jozsa (jozsait@gmail.com)
 */

#ifndef CELLFUNCTIONS_H
#define CELLFUNCTIONS_H

#include "FloatType.h"

/**
 * This function is responsible for the MRT initialization
 * @details The MRT collision model computes new values for the distribution function in momentum
 * space and uses the following formula:
 * \f[ f_i = f_i - \mathbf{M}^{-1}\mathbf{S}(\mathbf{m}-\mathbf{m}^{eq}) \f]
 * \f[ \mathbf{S} = diag(1.0, 1.4, 1.4, 1.0, 1.2, 1.0, 1.2, \omega, \omega) \f]
 * For more details check A. A. Mohamad - Lattice Boltzmann Method (book)
 *
 * @param[out] velMomMap mapping between velocity and momentum space \f$ \mathbf{M} \f$
 * @param[out] momCollMtx collision matrix in momentum space \f$ \mathbf{M}^{-1}\mathbf{S} \f$
 * @param omega collision frequency
 */
void MRTInitializer(FLOAT_TYPE* velMomMap, FLOAT_TYPE* momCollMtx, FLOAT_TYPE omega);

#endif
