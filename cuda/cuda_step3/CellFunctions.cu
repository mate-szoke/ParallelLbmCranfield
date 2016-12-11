#include <stdio.h>                  // printf();
#include <stdlib.h>                 // for malloc();
#include <stdbool.h>                // for bool type variables!
#include <math.h>                   // for sin,cos,pow... compile with -lm
#include "CellFunctions.h" // convenience

void MRTInitializer(FLOAT_TYPE* velMomMap, FLOAT_TYPE* momCollMtx, FLOAT_TYPE omega)
{
  int i, j;  // loop variables

  FLOAT_TYPE sm[9];
  const FLOAT_TYPE a1=1./36.;

  //this is M^-1
  FLOAT_TYPE tminv[81]=
     {4.*a1, -4.*a1,  4.*a1,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
      4.*a1,    -a1, -2.*a1,  6.*a1, -6.*a1,    0.0,    0.0,  9.*a1,    0.0,
      4.*a1,    -a1, -2.*a1,    0.0,    0.0,  6.*a1, -6.*a1, -9.*a1,    0.0,
      4.*a1,    -a1, -2.*a1, -6.*a1,  6.*a1,    0.0,    0.0,  9.*a1,    0.0,
      4.*a1,    -a1, -2.*a1,    0.0,    0.0, -6.*a1,  6.*a1, -9.*a1,    0.0,
      4.*a1,  2.*a1,     a1,  6.*a1,  3.*a1,  6.*a1,  3.*a1,    0.0,  9.*a1,
      4.*a1,  2.*a1,     a1, -6.*a1, -3.*a1,  6.*a1,  3.*a1,    0.0, -9.*a1,
      4.*a1,  2.*a1,     a1, -6.*a1, -3.*a1, -6.*a1, -3.*a1,    0.0,  9.*a1,
      4.*a1,  2.*a1,     a1,  6.*a1,  3.*a1, -6.*a1, -3.*a1,    0.0, -9.*a1};

  //this is M
  FLOAT_TYPE temp[81] =
   {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
    -4.,-1.,-1.,-1.,-1.,2.0,2.0,2.0,2.0,
    4.0,-2.,-2.,-2.,-2.,1.0,1.0,1.0,1.0,
    0.0,1.0,0.0,-1.,0.0,1.0,-1.,-1.,1.0,
    0.0,-2.,0.0,2.0,0.0,1.0,-1.,-1.,1.0,
    0.0,0.0,1.0,0.0,-1.,1.0,1.0,-1.,-1.,
    0.0,0.0,-2.,0.0,2.0,1.0,1.0,-1.,-1.,
    0.0,1.0,-1.,1.0,-1.,0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,1.0,-1.,1.0,-1.};

  for (i = 0; i < 81; i++)
  {
      velMomMap[i] = temp[i];
  }

  //This is diag(S)
  sm[0] = 1.0;
  sm[1] = 1.4;
  sm[2] = 1.4;
  sm[3] = 1.0;
  sm[4] = 1.2;
  sm[5] = 1.0;
  sm[6] = 1.2;
  sm[7] = omega;
  sm[8] = omega;


  for(i=0;i<9;i++)
  {
    for(j=0;j<9;j++)
    {
        momCollMtx[i*9+j]=tminv[i*9+j]*sm[j];
    }
  }

} // End of function
