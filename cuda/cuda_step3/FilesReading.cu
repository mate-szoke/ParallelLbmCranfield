/**
 * File reading functions
 * @file FilesReading.h
 * @author Adam Koleszar (adam.koleszar@gmail.com)
 */
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "FilesReading.h"
#include "ArrayUtils.h"
#include "TestUtils.h"

int getNumberOfLines(const char *filename)
{
  FILE *f = fopen(filename, "r");
  if (!f)
  {
    fprintf(stderr, "Error reading file %s: %s\n", filename, strerror(errno));
    return 0;
  }

  int lines=0;
  while (!feof(f))
  {
    if (fgetc(f) == '\n')
    {
      lines++;
    }
  }
  fclose(f);

  return lines;
}

/**
 * @brief Set enum options from init file
 *
 * @param args input parameters
 * @param ipr inlet profile
 * @param coll collision model
 * @param curved boundary type
 * @param opr outlet profile
 * @param format output format
 */
void setOptions(Arguments *args, int ipr, int coll, int curved, int opr, int format)
{
    args->inletProfile = (InletProfile)ipr;
    args->collisionModel = (CollisionModel)coll;
    args->boundaryType = (BoundaryType)curved;
    args->outletProfile = (OutletProfile)opr;
    args->outputFormat = (OutputFormat)format;
}

void readInitFile(const char* filename, Arguments *args)
{
  int ipr, coll, curved, opr, format;
  FILE *f_init = fopen(filename,"r");
  fscanf(f_init,FLOAT_FORMAT, &(args->u));
  fscanf(f_init,FLOAT_FORMAT, &(args->v));
  fscanf(f_init,FLOAT_FORMAT, &(args->rho));
  fscanf(f_init,FLOAT_FORMAT, &(args->viscosity));
  fscanf(f_init,"%d", &ipr);
  fscanf(f_init,"%d", &coll);
  fscanf(f_init,"%d", &curved);
  fscanf(f_init,"%d", &opr);
  fscanf(f_init,"%d", &(args->iterations));
  fscanf(f_init,"%d", &(args->autosaveEvery));
  fscanf(f_init,"%d", &(args->autosaveAfter));
  fscanf(f_init,"%d", &format);
  fscanf(f_init,"%d", &(args->boundaryId));
  fclose(f_init);

  setOptions(args, ipr, coll, curved, opr, format);
}

int readNodeFile(const char *filename, int **ni, int **nj, FLOAT_TYPE **nx, FLOAT_TYPE **ny, int **nf)
{
  int n = getNumberOfLines(filename);
  if (!n)
  {
    return 0;
  }

  *ni = createHostArrayInt(n);
  *nj = createHostArrayInt(n);
  *nx = createHostArrayFlt(n);
  *ny = createHostArrayFlt(n);
  *nf = createHostArrayInt(n);

  FILE *f = fopen(filename, "r");
  int i;
  for (i=0; i<n; ++i)
  {
    fscanf(f, "%d %d "FLOAT_FORMAT" "FLOAT_FORMAT" %d", (*ni)+i, (*nj)+i, (*nx)+i, (*ny)+i, (*nf)+i);
  }
  fclose(f);
  return n;
}

int readConnFile(const char *filename, int **ni, int **nj, int **dir, int **bc,
                 FLOAT_TYPE **bcx, FLOAT_TYPE **bcy, int **id)
{
  int n = getNumberOfLines(filename);
  if (!n)
  {
    return 0;
  }

  *ni = createHostArrayInt(n);
  *nj = createHostArrayInt(n);
  *dir = createHostArrayInt(n);
  *bc = createHostArrayInt(n);
  *bcx = createHostArrayFlt(n);
  *bcy = createHostArrayFlt(n);
  *id = createHostArrayInt(n);

  FILE *f = fopen(filename, "r");
  int i;
  for (i=0; i<n; ++i)
  {
    fscanf(f, "%d %d %d %d "FLOAT_FORMAT" "FLOAT_FORMAT" %d", (*ni)+i, (*nj)+i, (*dir)+i,
                                                              (*bc)+i, (*bcx)+i, (*bcy)+i,
                                                              (*id)+i);
  }
  fclose(f);
  return n;
}

int readResultFile(const char *filename, FLOAT_TYPE ***results, int **fluid)
{
  int i,n = getNumberOfLines(filename);
  if (!n)
  {
    return 0;
  }

  FILE *f = fopen(filename, "r");
  char firstline[256];
  fscanf(f, "%s", firstline);
  n -= 1;

  *fluid = createHostArrayInt(n);
  *results = create2DHostArrayFlt(n, 7);

  for (i=0; i<n; ++i)
  {
    fscanf(f, FLOAT_FORMAT", "FLOAT_FORMAT", "FLOAT_FORMAT", "FLOAT_FORMAT", "FLOAT_FORMAT", "FLOAT_FORMAT", "FLOAT_FORMAT", %d",
           results[0][0]+i, results[0][1]+i, results[0][2]+i, results[0][3]+i, results[0][4]+i, results[0][5]+i, results[0][6]+i, (*fluid)+i);
  }
  fclose(f);
  return n;
}

__host__ int compareFiles(const char* f1, const char* f2)
{
  printf("Comparing results\n");
  int l1 = getNumberOfLines(f1);
  int l2 = getNumberOfLines(f2);
  if (l1 != l2)
  {
    printf("Line number mismatch\n");
    exit(1);
  }
  int i;

  const char *columnName[] = {"x", "y", "u", "v", "vel_mag", "rho", "pressure"};

  int *fluid1;
  int *fluid2;
  FLOAT_TYPE **res1;
  FLOAT_TYPE **res2;

  printf("Reading files...");
  readResultFile(f1, &res1, &fluid1);
  printf("...");
  readResultFile(f2, &res2, &fluid2);
  printf("...done\nComparing results...\n");

  dim3 gridDim(l1/THREADS + 1);
  FLOAT_TYPE *da = createGpuArrayFlt(l1);
  FLOAT_TYPE *db = createGpuArrayFlt(l1);
  FLOAT_TYPE *dc = createGpuArrayFlt(l1);
  FLOAT_TYPE *dd = createGpuArrayFlt(l1);

  FLOAT_TYPE result[7];

  for (i=0; i<7; ++i)
  {
    cudaMemcpy(da, res1[i], SIZEFLT(l1), cudaMemcpyHostToDevice);
    cudaMemcpy(db, res2[i], SIZEFLT(l1), cudaMemcpyHostToDevice);
    result[i] = compareGpuArraySumFlt(da, db, dc, dd, l1);
  }

  int fluid = 0;
  for (i=0; i<l1; ++i)
  {
    fluid += (fluid1[i]-fluid2[i]);
  }
  printf("     array |            diff |      diff/nodes\n"
          "----------------------------------------------\n");
  printf("    fluids | %15d | %15d\n", fluid, fluid/l1);

  int b = 0;
  for (i=0; i<7; ++i)
  {
    printf("%10s | %15g | %15g\n", columnName[i], result[i], result[i]/l1);
    b |= result[i] > 0.001;
    free(res1[i]);
    free(res2[i]);
  }

  free(fluid1); free(fluid2);
  free(res1); free(res2);

  return b;
}

int getLastValue(int *arr, int n)
{
  return arr[n-1] + 1;
}

FLOAT_TYPE getGridSpacing(int *ni, int *nj, FLOAT_TYPE *nx, int n)
{
  int i;
  FLOAT_TYPE delta1, delta2;
  for (i=0; i<n; ++i)
  {
    if (ni[i] == 0 && nj[i] == 0)
    {
      delta1 = nx[i];
      break;
    }
  }
  for (i=0; i<n; ++i)
  {
    if (ni[i] == 1 && nj[i] == 0)
    {
      delta2 = nx[i];
      break;
    }
  }

  return (FLOAT_TYPE)fabs(delta2 - delta1);
}

int getNumInletNodes(int *bc, int *dir, int n)
{
  int i;
  int nodes = 0;
  for (i=0; i<n; ++i)
  {
    if (bc[i] == 2 && dir[i] >= 1 && dir[i] <= 4)
    {
      ++nodes;
    }
  }

  return nodes;
}

FLOAT_TYPE getMaxInletCoordY(int *bc, int *dir, FLOAT_TYPE *bcy, FLOAT_TYPE delta, int n)
{
  int i=0;
  FLOAT_TYPE maxY;
  while (bc[i] != 2) //inlet
  {
    maxY = bcy[++i];
  }
  for (i=0; i<n; ++i)
  {
    if (bc[i] == 2 && dir[i] >= 1 && dir[i] <= 4)
    {
      maxY = (bcy[i] > maxY) ? bcy[i] : maxY;
    }
  }

  return maxY + delta/2;
}

FLOAT_TYPE getMinInletCoordY(int *bc, int *dir, FLOAT_TYPE *bcy, FLOAT_TYPE delta, int n)
{
  int i=0;
  FLOAT_TYPE minY;
  while (bc[i] != 2) //inlet
  {
    minY = bcy[++i];
  }

  for (i=0; i<n; ++i)
  {
    if (bc[i] == 2 && dir[i] >= 1 && dir[i] <= 4)
    {
      minY = (bcy[i] < minY) ? bcy[i] : minY;
    }
  }

  return minY - delta/2;
}

void CompDataNode(FLOAT_TYPE *Delta, int *m, int *n, int *ni, int *nj, FLOAT_TYPE *nx, int size)
{
  int i; // variable for the loop
  FLOAT_TYPE DeltaP1, DeltaP2; // local grid spacing

  *n = *(ni+size-1)+1; // number of rows
  *m = *(nj+size-1)+1; // number of columns

  for(i=0;i<size;i++)
  {
    if(*(ni+i)==0 && *(nj+i)==0)
    {
      DeltaP1=*(nx+i);
    }
    if(*(ni+i)==1 && *(nj+i)==0)
    {
      DeltaP2=*(nx+i);
    }
  }

  *Delta = (max(DeltaP1,DeltaP2)-min(DeltaP1,DeltaP2)); // grid spacing
}

void CompDataConn(int* NumInletNodes, FLOAT_TYPE* MaxInletCoordY,
	FLOAT_TYPE* MinInletCoordY, int* BCconn0, int* BCconn1, int* BCconn2,
  int* BCconn3, FLOAT_TYPE* BCconn4, FLOAT_TYPE* BCconn5, int* BCconn6, int* NumConn, FLOAT_TYPE* Delta)
{
  int i=0; // counter

  while(*(BCconn3+i)!=2)
  {
      MaxInletCoordY[0] = *(BCconn5+i+1); // maximum Y coordinate of the inlet line
      MinInletCoordY[0] = *(BCconn5+i+1); // minimum Y coordinate of the inlet line
      i++;
  }

  *NumInletNodes = 0; // number of inlet nodes

  for (i=0; i< *NumConn;i++)
  {
      if(*(BCconn3+i)==2){
          if(*(BCconn2+i)==1 || *(BCconn2+i)==2 || *(BCconn2+i)==3 || *(BCconn2+i)==4){
              if(*(BCconn5+i)>*MaxInletCoordY){
                  *MaxInletCoordY = *(BCconn5+i);
              }
              if(*(BCconn5+i)<MinInletCoordY[0]){
                  *MinInletCoordY = *(BCconn5+i);
              }
              *NumInletNodes=*NumInletNodes+1;
          }
      }
  }

  (*MaxInletCoordY) = (*MaxInletCoordY)+(*Delta)/2;
  (*MinInletCoordY) = (*MinInletCoordY)-(*Delta)/2;
}
