#include <stdio.h>                  // printf();
#include <stdlib.h>                 // for malloc();
#include <stdbool.h>                // for bool type variables!
#include <math.h>                   // for sin,cos,pow... compile with -lm
#include "include/ShellFunctions.h" // convenience

/*==================================================
=========Initialization for the MRT model===========
==================================================*/
// This function fills up tm and stimv with variables
void MRTInitializer(float* tm, float* stmiv, float Omega)
{

  ///////////// Declarations ////////////////
  int i, j;  // loop variables
  
  // declarations for this collision model
  float sm[9];
  // float ev[9][9];
  const float a1=1./36.;
  float tminv[81]=
     {4.*a1, -4.*a1,  4.*a1,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
      4.*a1,    -a1, -2.*a1,  6.*a1, -6.*a1,    0.0,    0.0,  9.*a1,    0.0,
      4.*a1,    -a1, -2.*a1,    0.0,    0.0,  6.*a1, -6.*a1, -9.*a1,    0.0,
      4.*a1,    -a1, -2.*a1, -6.*a1,  6.*a1,    0.0,    0.0,  9.*a1,    0.0,
      4.*a1,    -a1, -2.*a1,    0.0,    0.0, -6.*a1,  6.*a1, -9.*a1,    0.0,
      4.*a1,  2.*a1,     a1,  6.*a1,  3.*a1,  6.*a1,  3.*a1,    0.0,  9.*a1,
      4.*a1,  2.*a1,     a1, -6.*a1, -3.*a1,  6.*a1,  3.*a1,    0.0, -9.*a1,
      4.*a1,  2.*a1,     a1, -6.*a1, -3.*a1, -6.*a1, -3.*a1,    0.0,  9.*a1,
      4.*a1,  2.*a1,     a1,  6.*a1,  3.*a1, -6.*a1, -3.*a1,    0.0, -9.*a1};
     
  float temp[81] =
   {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
    -4.,-1.,-1.,-1.,-1.,2.0,2.0,2.0,2.0,
    4.0,-2.,-2.,-2.,-2.,1.0,1.0,1.0,1.0,
    0.0,1.0,0.0,-1.,0.0,1.0,-1.,-1.,1.0,
    0.0,-2.,0.0,2.0,0.0,1.0,-1.,-1.,1.0,
    0.0,0.0,1.0,0.0,-1.,1.0,1.0,-1.,-1.,
    0.0,0.0,-2.,0.0,2.0,1.0,1.0,-1.,-1.,
    0.0,1.0,-1.,1.0,-1.,0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,1.0,-1.,1.0,-1.};


  ///////////// Fill up variables ////////////////
  
  // Filling up tm
  for (i = 0; i < 81; i++)
  {
      tm[i] = temp[i];
  }


  // Filling up stimv
  sm[0] = 1.0;
  sm[1] = 1.4;
  sm[2] = 1.4;
  sm[3] = 1.0;
  sm[4] = 1.2;
  sm[5] = 1.0;
  sm[6] = 1.2;
  sm[7] = Omega;
  sm[8] = Omega;


  for(i=0;i<9;i++)
  {
    for(j=0;j<9;j++)
    {
        stmiv[i*9+j]=tminv[i*9+j]*sm[j];
    }
  }

} // End of function
