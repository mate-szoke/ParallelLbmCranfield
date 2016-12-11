// File importer!

#include <stdio.h>                  // printf();
#include <stdlib.h>                 // for malloc();
#include "include/ShellFunctions.h" // For convenience


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  ///////////////// Read the file including the nodes /////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

void ReadNodLines(char* NodeDataFile, int* NumOfLines)
{
  ////////////////////////////////////////////////////
  ///////////////////// Declare //////////////////////
  ////////////////////////////////////////////////////

  //declare the file pointers of the nodes
  FILE *fp_nodes;              

  // number of lines in the files
  int lines;               

  // number of lines in the files
  int Ir1,Ir2,Ir3;

  // variables to read floats
  float Fr1,Fr2;           

  // D2node.dat includes:
  //  _______________________________________________________________________________
  // |         1          |        2         |    3    |    4    |         5        |
  // |   node index i int | node index j int | x coord | y coord | solid/fluid int  |
  // |____________________|__________________|_________|_________|__________________|
  //  solid/fluid: 0 -> solid; 1 -> fluid

  fp_nodes = fopen(NodeDataFile,"r"); // open the file to count the lines

  if (fp_nodes==NULL) // if the file does not exist
  {
    fprintf(stderr, "Can't open the file %s\n", NodeDataFile);
    exit(0); // ERROR! 
  }
  else // if the file exist the following while goes to the end
  {
    lines=0;
    while (fscanf(fp_nodes, "%d %d %f %f %d",&Ir1,&Ir2,&Fr1,&Fr2,&Ir3) == 5)
    {
      lines++; // counter of the lines
    }

    fclose(fp_nodes); // close the file
   }

  *NumOfLines = lines; // number of lines
}


void ReadNodes(char* NodeDataFile, int* NumOfLines, int *Nodes0, int *Nodes1,
               float *Nodes2, float *Nodes3, int *Nodes4)
{

  // variable for loops
  int j;

  //declare the file pointers of the nodes
  FILE *fp_nodes;  

  fp_nodes = fopen(NodeDataFile,"r"); // and open again to read the data

  for (j=0;j<(*NumOfLines);j++) // read the data from the file to the variables
  {
    fscanf(fp_nodes,"%d %d %f %f %d", Nodes0+j,Nodes1+j,
    Nodes2+j,Nodes3+j,Nodes4+j);
  }

  fclose(fp_nodes); // close the file
}

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////// Read the file including the connections /////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

void ReadBCconLines(char* BCconnectorDataFile, int* NumOfLines)
{
  ////////////////////////////////////////////////////
  ///////////////////// Declare //////////////////////
  ////////////////////////////////////////////////////

  //declare the file pointers of the connectors
  FILE *fp_connect;

  // number of lines in the files
  int lines; 

  // number of lines in the files
  int Ir1,Ir2,Ir3,Ir4,Ir5;

  // variables to read floats
  float Fr1,Fr2;  

  // BCconnectors.dat includes:
  //  __________________________________________________________________________________
  // |        1      |        2     |     3     |    4    |    5     |    6     |  7   |
  // |  node index i | node index j | latticeID | BC type | x coord  | y coord  | ???  |
  // |_______________|______________| _________ |_________|__________|__________|______| 
  //     
  // lattice ID is based on the following speed model and it depends on the BC
  //  ID lattice
  //        6       2       5
  //          \     |     /
  //            \   |   /
  //              \ | /
  //        3 - - - 0 - - - 1
  //              / | \
  //            /   |   \
  //          /     |     \
  //        7       4      8
  // BC types: *1->wall; *2->inlet; *3->outlet
  // WHAT IS IN THE LAST COLUMN???


  fp_connect = fopen(BCconnectorDataFile,"r"); // open the file to count the lines

  if (fp_connect==NULL) // if the file does not exist
  {
    fprintf(stderr, "Can't open the file %s!\n",BCconnectorDataFile);
    exit(1); // ERROR!
  }
  else // if the file exist the following while goes to the end
  {
    lines=0;
    while (fscanf(fp_connect, "%d %d %d %d %f %f %d",&Ir1,&Ir2,&Ir3,
            &Ir4,&Fr1,&Fr2,&Ir5) == 7)
    {
      lines++; // counter of the lines
    } 

    fclose(fp_connect); // close the file
  }

  *NumOfLines = lines; // number of lines
}




void ReadBCconn(char* BCconnectorDataFile, int* NumOfLines, int* BCconn0,
                int* BCconn1, int* BCconn2, int* BCconn3, float* BCconn4,
                float* BCconn5, int* BCconn6)
{

  //declare the file pointers of the connectors
  FILE *fp_connect;

  // variable for loops
  int j;  
  

  fp_connect = fopen(BCconnectorDataFile,"r"); // and open again to read the data

  for (j=0;j<(*NumOfLines);j++) // read the data from the file to the variables
  {
  fscanf(fp_connect,"%d %d %d %d %f %f %d",BCconn0+j,BCconn1+j,
         BCconn2+j,BCconn3+j,BCconn4+j,BCconn5+j,BCconn6+j);
  }
    
  fclose(fp_connect); // close the file
}


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  ///////////////////// Constants from D2node.dat /////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

void CompDataNode(float* Delta, int* m,  int* n, int *Nodes0, int *Nodes1,
                  float *Nodes2, float *Nodes3, int *Nodes4,  int* NumNodes)
{

  ////////////////////////////////////////////////////
  ///////////////////// Declare //////////////////////
  ////////////////////////////////////////////////////

  int i; // variable for the loop
  float DeltaP1, DeltaP2; // local grid spacing

  *n = *(Nodes0+*NumNodes-1)+1; // number of rows
  *m = *(Nodes1+*NumNodes-1)+1; // number of columns

  for(i=0;i<*NumNodes;i++)
  {
    if(*(Nodes0+i)==0 && *(Nodes1+i)==0)
    {
      DeltaP1=*(Nodes2+i);
    }
    if(*(Nodes0+i)==1 && *(Nodes1+i)==0)
    {
      DeltaP2=*(Nodes2+i);
    }
  }

  *Delta = (max(DeltaP1,DeltaP2)-min(DeltaP1,DeltaP2)); // grid spacing 
}


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////// Constants from BCconnectors.dat /////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

void CompDataConn(int* NumInletNodes, float* MaxInletCoordY,
	float* MinInletCoordY, int* BCconn0, int* BCconn1, int* BCconn2,
  int* BCconn3, float* BCconn4, float* BCconn5, int* BCconn6, int* NumConn, float* Delta)
{

  ////////////////////////////////////////////////////
  ///////////////////// Declare //////////////////////
  ////////////////////////////////////////////////////

  int i=0; // counter

  while(*(BCconn3+i)!=2)
  {
      MaxInletCoordY[0] = *(BCconn5+i+1); // maximum Y coordinate of the inlet line
      MinInletCoordY[0] = *(BCconn5+i+1); // minimum Y coordinate of the inlet line
      i++;
  }

  *NumInletNodes = 0; // unmber of inlet nodes

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


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  ///////////////////////// Read the *.ini file ///////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

void ReadIniData(char* IniFileName, float* Uavg, float* Vavg, float* rho_ini,
                 float* Viscosity, int* InletProfile, int* CollisionModel,
                 int* CurvedBoundaries, int* OutletProfile, int* Iterations,
                 int* AutosaveEvery, int* AutosaveAfter, int* PostprocProg,
                 int* CalculateDragLift)
{
  FILE *f_init; // *.ini file pointer
  f_init = fopen(IniFileName,"r"); // open the file
  fscanf(f_init,"%f", Uavg); // U velocity to initialize
  fscanf(f_init,"%f", Vavg); // V velocity to initialize
  fscanf(f_init,"%f", rho_ini); // Rho velocity to initialize
  fscanf(f_init,"%f", Viscosity); // Viscosity
  fscanf(f_init,"%d", InletProfile); // inlet profile (yes/no)
  fscanf(f_init,"%d", CollisionModel); // collision model (BGKW/TRT/MRT)
  fscanf(f_init,"%d", CurvedBoundaries); // curved boundaries (yes/no)
  fscanf(f_init,"%d", OutletProfile); // outlet profile (yes/no)
  fscanf(f_init,"%d", Iterations); // # of iterations
  fscanf(f_init,"%d", AutosaveEvery); // autosave every #th of iteration
  fscanf(f_init,"%d", AutosaveAfter); // autosave after the #th iteration
  fscanf(f_init,"%d", PostprocProg); // program of postp rocessing (Paraview/TECplot)
  fscanf(f_init,"%d", CalculateDragLift);  // calculate drag & lift? if > 0 than on which BC_ID
  fclose(f_init); // close the file
}
