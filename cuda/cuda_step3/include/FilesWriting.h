/**
 * Functions for file writing
 * @file FilesWriting.h
 * @author Istvan Tamas Jozsa (jozsait@gmail.com)
 */

#ifndef FILESWRITING_H
#define FILESWRITING_H

#include "FloatType.h"
#include "Arguments.h"

/**
 * @brief Print results to output file
 *
 * @param OutputFile name of the outpu file
 * @param CoordX,CoordY x,y coordinates
 * @param U,V velocity x,y coordinates
 * @param Rho density
 * @param Fluid node type (0: solid, 1: fluid)
 * @param n number of rows
 * @param m number of columns
 * @param outputFormat output file format
 */
void WriteResults(char* OutputFile, FLOAT_TYPE* CoordX, FLOAT_TYPE* CoordY, FLOAT_TYPE* U,
                  FLOAT_TYPE* V, FLOAT_TYPE* Rho, int* Fluid, int n, int m, OutputFormat outputFormat);

#endif
