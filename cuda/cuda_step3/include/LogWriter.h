/**
 * Log writing function
 * @file LogWriter.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#ifndef LOG_WRITER_H
#define LOG_WRITER_H

#include "FloatType.h"
#include "Arguments.h"

/// Enumeration of the tasks
enum taskTime
{
    T_INIT=0, ///< initialisation time
    T_ITER,   ///< iteration time
    T_COLL,   ///< collision time
    T_STRM,   ///< streaming time
    T_BNDC,   ///< bc time
    T_MACR,   ///< macroscopic time
    T_RESI,   ///< residuals time
    T_WRIT,   ///< file write time
    T_OALL    ///< overall time
};

/**
 *  @brief Print initialisation log
 *
 *  @param [in] filename       output filename
 *  @param [in] args           input parameters
 *  @param [in] delta          grid spacing
 *  @param [in] m              number of columns
 *  @param [in] n              number of rows
 *  @param [in] numInletNodes  number of inlet nodes
 *  @param [in] maxInletCoordY maximum inlet coordinate y
 *  @param [in] minInletCoordY minimum inlet coordinate y
 */
void writeInitLog(const char* filename, Arguments *args, FLOAT_TYPE delta, int m, int n,
                  int numInletNodes, FLOAT_TYPE maxInletCoordY, FLOAT_TYPE minInletCoordY);

/**
 *  @brief Print node numbers log
 *
 *  @param [in] filename       output filename
 *  @param [in] numNodes       number of nodes
 *  @param [in] numConns       number of conditions
 *  @param [in] bcCount        number of boundary nodes
 */
void writeNodeNumbers(const char* filename, int numNodes, int numConns, int bcCount);

/**
 * @brief Print time results into the log
 *
 * @param filename output filename
 * @param taskTimes time values
 */
void writeEndLog(const char *filename, FLOAT_TYPE *taskTimes);

/**
 * @brief Print time results into the time file
 *
 * @param filename output filename
 * @param taskTimes time values
 */
void writeTimerLog(const char *filename, FLOAT_TYPE *taskTimes);

/**
 * @brief Print residual data
 *
 * @param filename output filename
 * @param norm norm of the distribution function change \f$ f - f_{coll} \f$
 * @param drag,lift drag and lift
 * @param size mesh size (MxN)
 * @param n number of iterations
 */
void writeResiduals(const char *filename, FLOAT_TYPE *norm, FLOAT_TYPE *drag, FLOAT_TYPE *lift, int size, int n);

#endif
