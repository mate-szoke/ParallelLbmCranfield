#include <stdlib.h>    // for calloc();
#include <assert.h>    // ensure successfull allocation
#include <stdbool.h>   // bool variables
#include <stdio.h>     // printf...
#include <string.h> 
#include <dirent.h>    // Directory management
#include <sys/stat.h>  // system commands ?
#include <sys/types.h> // Extra types of variables, google it
#include <unistd.h>    // standard symbolic constants and types

#include "include/ShellFunctions.h" // convenience

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
	MyArray = calloc(length,sizeof(int));
 	assert(MyArray != NULL);
	return MyArray;
}


// This function creates a 1D Array (vector) :: MyReal
// It returns a pointer
MyReal *Create1DArrayMyReal(int length)
{
	MyReal *MyArray;
	MyArray = calloc(length,sizeof(MyReal));
	assert(MyArray != NULL);
	return MyArray;
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////// 2D array allocators /////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// This function creates a 2D Array (matrix) :: MyReal
// It returns a pointer
int **Create2DArrayInt(int width, int height)
{
	int **MyMatrix;
	int i;
	MyMatrix = calloc(height,sizeof(int*));
	assert(MyMatrix != NULL);
	for (i = 0; i < height; i++)
		MyMatrix[i] = calloc(width,sizeof(int));
	assert(MyMatrix != NULL);
	return MyMatrix;
}


// This function creates a 2D Array (matrix) :: MyReal
// It returns a pointer
MyReal **Create2DArrayMyReal(int width, int height)
{
	MyReal **MyMatrix;
	int i;
	MyMatrix = calloc(height,sizeof(MyReal*));
	assert(MyMatrix != NULL);
	for (i = 0; i < height; i++)
		MyMatrix[i] = calloc(width,sizeof(MyReal));
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
	
	MyMatrix = calloc(height,sizeof(int**));
	assert(MyMatrix != NULL);
	for (i = 0; i < height; i++)
	{
		MyMatrix[i] = calloc(width,sizeof(int*));
		assert(MyMatrix != NULL);
		for (j = 0; j < width; j++)
		{
			MyMatrix[i][j] = calloc(depth,sizeof(int));
			assert(MyMatrix != NULL);
			for (k = 0; k < depth; k++)
				MyMatrix[i][j][k] = 0;
		}
	}
	return MyMatrix;
}


// This function creates a 3D Array (matrix of vectors) :: MyReal
// It returns a pointer
MyReal ***Create3DArrayMyReal(int width, int height, int depth)
{
	MyReal ***MyMatrix;
	int i, j, k;
	
	MyMatrix = calloc(height,sizeof(MyReal**));
	assert(MyMatrix != NULL);
	for (i = 0; i < height; i++)
	{
		MyMatrix[i] = calloc(width,sizeof(MyReal*));
		assert(MyMatrix != NULL);
		for (j = 0; j < width; j++)
		{
			MyMatrix[i][j] = calloc(depth,sizeof(MyReal));
			assert(MyMatrix != NULL);
			for (k = 0; k < depth; k++)
				MyMatrix[i][j][k] = 0;
		}
	}
	return MyMatrix;
}


// This function creates a 3D Array (matrix of vectors) :: BOOL
// It returns a pointer
bool ***Create3DArrayBool(int width, int height, int depth)
{
	bool ***MyMatrix;
	int i, j, k;
	
	MyMatrix = calloc(height,sizeof(bool**));
	assert(MyMatrix != NULL);
	for (i = 0; i < height; i++)
	{
		MyMatrix[i] = calloc(width,sizeof(bool*));
		assert(MyMatrix != NULL);
		for (j = 0; j < width; j++)
		{
			MyMatrix[i][j] = calloc(depth,sizeof(bool));
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
