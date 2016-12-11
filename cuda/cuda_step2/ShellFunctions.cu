#include <stdlib.h>    // for calloc();
#include <assert.h>    // ensure successfull allocation
#include <stdbool.h>   // bool variables
#include <stdio.h>     // printf...
#include <string.h> 

#include <dirent.h>    // Directory management
#include <sys/stat.h>  // system commands ?
#include <sys/types.h> // Extra types of variables, google it
#include <unistd.h>    // standard symbolic constants and types


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////// 1D array allocators /////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// This function creates a 1D Array (vector) :: INT
// It returns a pointer
int *Create1DArrayInt(int length)
{
	int *MyArray;
	MyArray = (int *)calloc(length,sizeof(int));
 	assert(MyArray != NULL);
	return MyArray;
}


// This function creates a 1D Array (vector) :: FLOAT
// It returns a pointer
float *Create1DArrayFloat(int length)
{
	float *MyArray;
	MyArray = (float *)calloc(length,sizeof(float));
	assert(MyArray != NULL);
	return MyArray;
}


// This function creates a 1D Array (vector) :: DOUBLE
// It returns a pointer
double *Create1DArrayDouble(int length)
{
	double *MyArray;
	MyArray = (double *)calloc(length,sizeof(double));
	assert(MyArray != NULL);
	return MyArray;
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////// 2D array allocators /////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// This function creates a 2D Array (matrix) :: FLOAT
// It returns a pointer
int **Create2DArrayInt(int width, int height)
{
	int **MyMatrix;
	int i;
	MyMatrix = (int **)calloc(height,sizeof(int*));
	assert(MyMatrix != NULL);
	for (i = 0; i < height; i++)
		MyMatrix[i] = (int *)calloc(width,sizeof(int));
	assert(MyMatrix != NULL);
	return MyMatrix;
}


// This function creates a 2D Array (matrix) :: FLOAT
// It returns a pointer
float **Create2DArrayFloat(int width, int height)
{
	float **MyMatrix;
	int i;
	MyMatrix = (float **)calloc(height,sizeof(float*));
	assert(MyMatrix != NULL);
	for (i = 0; i < height; i++)
		MyMatrix[i] = (float *)calloc(width,sizeof(float));
	assert(MyMatrix != NULL);
	return MyMatrix;
}


// This function creates a 2D Array (matrix) :: DOUBLE
// It returns a pointer
double **Create2DArrayDouble(int width, int height)
{
	double **MyMatrix;
	int i;
	MyMatrix = (double **)calloc(height,sizeof(double*));
	assert(MyMatrix != NULL);
	for (i = 0; i < height; i++)
		MyMatrix[i] = (double *)calloc(width,sizeof(double));
	assert(MyMatrix != NULL);
	return MyMatrix;
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////// 3D array allocators /////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// This function creates a 3D Array (matrix of vectors) :: INT
// It returns a pointer
int ***Create3DArrayInt(int width, int height, int depth)
{
	int ***MyMatrix;
	int i, j, k;
	
	MyMatrix = (int ***)calloc(height,sizeof(int**));
	assert(MyMatrix != NULL);
	for (i = 0; i < height; i++)
	{
		MyMatrix[i] = (int **)calloc(width,sizeof(int*));
		assert(MyMatrix != NULL);
		for (j = 0; j < width; j++)
		{
			MyMatrix[i][j] = (int *)calloc(depth,sizeof(int));
			assert(MyMatrix != NULL);
			for (k = 0; k < depth; k++)
				MyMatrix[i][j][k] = 0;
		}
	}
	return MyMatrix;
}


// This function creates a 3D Array (matrix of vectors) :: FLOAT
// It returns a pointer
float ***Create3DArrayFloat(int width, int height, int depth)
{
	float ***MyMatrix;
	int i, j, k;
	
	MyMatrix = (float ***)calloc(height,sizeof(float**));
	assert(MyMatrix != NULL);
	for (i = 0; i < height; i++)
	{
		MyMatrix[i] = (float **)calloc(width,sizeof(float*));
		assert(MyMatrix != NULL);
		for (j = 0; j < width; j++)
		{
			MyMatrix[i][j] = (float *)calloc(depth,sizeof(float));
			assert(MyMatrix != NULL);
			for (k = 0; k < depth; k++)
				MyMatrix[i][j][k] = 0;
		}
	}
	return MyMatrix;
}


// This function creates a 3D Array (matrix of vectors) :: DOUBLE
// It returns a pointer
double ***Create3DArrayDouble(int width, int height, int depth)
{
	double ***MyMatrix;
	int i, j, k;
	
	MyMatrix = (double ***)calloc(height,sizeof(double**));
	assert(MyMatrix != NULL);
	for (i = 0; i < height; i++)
	{
		MyMatrix[i] = (double **)calloc(width,sizeof(double*));
		assert(MyMatrix != NULL);
		for (j = 0; j < width; j++)
		{
			MyMatrix[i][j] = (double *)calloc(depth,sizeof(double));
			assert(MyMatrix != NULL);
			for (k = 0; k < depth; k++)
				MyMatrix[i][j][k] = 0;
		}
	}
	return MyMatrix;
}


// This function creates a 3D Array (matrix of vectors) :: DOUBLE
// It returns a pointer
bool ***Create3DArrayBool(int width, int height, int depth)
{
	bool ***MyMatrix;
	int i, j, k;
	
	MyMatrix = (bool ***)calloc(height,sizeof(bool**));
	assert(MyMatrix != NULL);
	for (i = 0; i < height; i++)
	{
		MyMatrix[i] = (bool **)calloc(width,sizeof(bool*));
		assert(MyMatrix != NULL);
		for (j = 0; j < width; j++)
		{
			MyMatrix[i][j] = (bool *)calloc(depth,sizeof(bool));
			assert(MyMatrix != NULL);
			for (k = 0; k < depth; k++)
				MyMatrix[i][j][k] = 0;
		}
	}
	return MyMatrix;
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//////////////////////////// Create directory ///////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void CreateDirectory(char* MainWorkDir)
{
    // Create the working directory, exit if it already exist!
    DIR *MyDirVar;
    
    MyDirVar = opendir (MainWorkDir);   // Open the following directory
    if (MyDirVar != NULL)               // if dp equals to NULL dir does not exist
    {
        //printf("Directory %s already exist, move on.\n", MainWorkDir);
    }   
    else
    {
        mkdir(MainWorkDir, S_IRWXU|S_IRGRP|S_IXGRP); // Create the directiory
        //printf("Directory does not exist yet... creating: %s\n", MainWorkDir);
    }
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////// Add string //////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void StringAddition(char* first, char* second, char* result)
{
  strcat(result, first); 
  strcat(result, second); 
}
