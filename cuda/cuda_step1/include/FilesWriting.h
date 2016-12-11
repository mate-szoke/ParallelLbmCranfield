/*This function is responsible for the output field variables.*/

#ifndef FILESWRITING_H
#define FILESWRITING_H

void WriteResults(char* OutputFile, struct CellProps_const *Cells_const_h,
                  struct CellProps_var *Cells_var_h, struct CellProps_var_9d *Cells_var_9d_h, int* n, int* m, int* postproc_prog);

#endif
