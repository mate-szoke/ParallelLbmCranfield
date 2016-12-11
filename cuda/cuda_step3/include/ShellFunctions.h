/**
 * Allocator and shell functions
 * @file ShellFunctions.h
 * @author Istvan Tamas Jozsa (jozsait@gmail.com)
 */

#ifndef SHELLFUNCTIONS_H
#define SHELLFUNCTIONS_H

#include <string.h>   // string handling
#include "FloatType.h"

/**
 *  @brief Make directory
 *
 *  @param [in] MainWorkDir directory to be created
 */
void CreateDirectory(const char* MainWorkDir);

/**
 *  @brief Concatenate strings
 *  @note not used
 *
 *  @param [in] first  first string
 *  @param [in] second second string
 *  @param [out] result concatenated string
 */
void StringAddition(char* first, char* second, char* result);

/**
 *  @brief Get maximum of two values
 *
 *  @param [in] a,b values to compare
 *  @return maximum of (a,b)
 */
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

/**
 *  @brief Get minimum of two values
 *
 *  @param [in] a,b values to compare
 *  @return minimum of (a,b)
 */
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })


#endif
