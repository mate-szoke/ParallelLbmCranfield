#include <stdlib.h>    // for calloc();
#include <stdbool.h>   // bool variables
#include <stdio.h>     // printf...
#include <string.h>

#if defined(_WIN32) || defined(_WIN64)
#include <direct.h>
#include <errno.h>
void CreateDirectory(const char *dir)
{
  if (_mkdir(dir) != 0)
  {
    if (errno != EEXIST)
    {
      fprintf(stderr, "Error creating directory: %s", strerror(errno));
      exit (1);
    }
  }
}
#else
#include <dirent.h>    // Directory management
#include <sys/stat.h>  // system commands ?
#include <sys/types.h> // Extra types of variables, google it
#include <unistd.h>    // standard symbolic constants and types

#include "include/ShellFunctions.h"

void CreateDirectory(const char* MainWorkDir)
{
    // Create the working directory, exit if it already exist!
    DIR *MyDirVar;

    MyDirVar = opendir (MainWorkDir);   // Open the following directory
    if (MyDirVar != NULL)               // if dp equals to NULL dir does not exist
    {
        //printf("Directory %s already exist, move on.\n", MainWorkDir);
        closedir(MyDirVar);
    }
    else
    {
        mkdir(MainWorkDir, S_IRWXU|S_IRGRP|S_IXGRP); // Create the directiory
        //printf("Directory does not exist yet... creating: %s\n", MainWorkDir);
    }
}
#endif

void StringAddition(char* first, char* second, char* result)
{
  strcat(result, first);
  strcat(result, second);
}

