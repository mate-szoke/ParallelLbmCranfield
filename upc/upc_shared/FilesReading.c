// File importer!

#include <stdio.h>                  // printf();
#include <stdlib.h>                 // for malloc();
//#include <upc_cray.h>            // Required for UPC 

#include "include/ShellFunctions.h" // For convenience

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
///////////////// Read the file including the nodes /////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

float **ReadNodes(char* NodeDataFile)
{

  ////////////////////////////////////////////////////
  ///////////////////// Declare //////////////////////
  ////////////////////////////////////////////////////

  //declare the file pointers of the nodes
  FILE *fp_nodes;

  // variable for loops
  int i, j;               

  // number of lines in the files
  int lines;               

  // number of lines in the files
  int Ir1,Ir2,Ir3;

  // variables to read floats
  float Fr1,Fr2;           

  // variable for the Nodal data
  float **Nodes;

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
    return 0; 
  }
  else // if the file exist the following while goes to the end
  {
    lines=0;
    while (fscanf(fp_nodes, "%d %d %f %f %d",&Ir1,&Ir2,&Fr1,&Fr2,&Ir3) == 5)
    { lines++; } // counter of the lines

    fclose(fp_nodes); // close the file

    // allocate matrix for nodes
    Nodes = calloc(lines,sizeof(float*));
    for (i = 0; i < lines; ++i)
    Nodes[i] = calloc(5,sizeof(float));

    fp_nodes = fopen(NodeDataFile,"r"); // and open again to read the data

    for (j=0;j<lines;j++) // read the data from the file to the variables
    {
      fscanf(fp_nodes,"%f %f %f %f %f", &Nodes[j][0],&Nodes[j][1],
      &Nodes[j][2],&Nodes[j][3],&Nodes[j][4]);
    }

    fclose(fp_nodes); // close the file
    
    *NumNodes = lines; // number of lines

    return Nodes;
  }
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////// Read the file including the connections /////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

float **ReadBCconn(char* BCconnectorDataFile)
{

  ////////////////////////////////////////////////////
  ///////////////////// Declare //////////////////////
  ////////////////////////////////////////////////////

  //declare the file pointers of the connectors
  FILE *fp_connect;

  // variable for loops
  int i, j;               

  // number of lines in the files
  int lines; 

  // number of lines in the files
  int Ir1,Ir2,Ir3,Ir4,Ir5;

  // variables to read floats
  float Fr1,Fr2;  

  // variable for the conncection data
  float **BCconn;

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
    return 0;
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

    // allocate matrix for connectors
    BCconn = malloc(lines*sizeof(float*));
    for (i = 0; i < lines; ++i)
    BCconn[i] = malloc(7*sizeof(float));

    fp_connect = fopen(BCconnectorDataFile,"r"); // and open again to read the data

    for (j=0;j<lines;j++) // read the data from the file to the variables
    {
    fscanf(fp_connect,"%f %f %f %f %f %f %f",&BCconn[j][0],&BCconn[j][1],
           &BCconn[j][2], &BCconn[j][3],&BCconn[j][4],&BCconn[j][5],&BCconn[j][6]);
    }
    
    fclose(fp_connect); // close the file

    *NumConn = lines; // number of lines

    return BCconn;
  }
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
///////////////////// Constants from D2node.dat /////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void CompDataNode(float **Nodes)
{

  ////////////////////////////////////////////////////
  ///////////////////// Declare //////////////////////
  ////////////////////////////////////////////////////

  int i; // variable for the loop
  float DeltaP1, DeltaP2; // local grid spacing

  *n = Nodes[*NumNodes-1][0]+1; // number of rows
  *m = Nodes[*NumNodes-1][1]+1; // number of columns

  for(i=0;i<*NumNodes;i++)
  {
    if(Nodes[i][0]==0 && Nodes[i][1]==0)
    {
      DeltaP1=Nodes[i][2];
    }
    if(Nodes[i][0]==1 && Nodes[i][1]==0)
    {
      DeltaP2=Nodes[i][2];
    }
  }

  *Delta = (max(DeltaP1,DeltaP2)-min(DeltaP1,DeltaP2)); // grid spacing 
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////// Constants from BCconnectors.dat /////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void CompDataConn(float** BCconn)
{

  ////////////////////////////////////////////////////
  ///////////////////// Declare //////////////////////
  ////////////////////////////////////////////////////

  int i=0; // counter

  while(BCconn[i][3]!=2)
  {
      MaxInletCoordY[0] = BCconn[i+1][5]; // maximum Y coordinate of the inlet line
      MinInletCoordY[0] = BCconn[i+1][5]; // minimum Y coordinate of the inlet line
      i++;
  }

  *NumInletNodes = 0; // unmber of inlet nodes

  for (i=0; i< *NumConn;i++)
  {
      if(BCconn[i][3]==2){
          if(BCconn[i][2]==1 || BCconn[i][2]==2 || BCconn[i][2]==3 || BCconn[i][2]==4){
              if(BCconn[i][5]>*MaxInletCoordY){
                  *MaxInletCoordY = BCconn[i][5];
              }
              if(BCconn[i][5]<MinInletCoordY[0]){
                  *MinInletCoordY = BCconn[i][5];
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
  f_init = fopen(IniFileName,"r");         // open the file
  fscanf(f_init,"%f", Uavg);               // U velocity to initialize
  fscanf(f_init,"%f", Vavg);               // V velocity to initialize
  fscanf(f_init,"%f", rho_ini);            // Rho velocity to initialize
  fscanf(f_init,"%f", Viscosity);          // Viscosity
  fscanf(f_init,"%d", InletProfile);       // inlet profile (yes/no)
  fscanf(f_init,"%d", CollisionModel);     // collision model (BGKW/TRT/MRT)
  fscanf(f_init,"%d", CurvedBoundaries);   // curved boundaries (yes/no)
  fscanf(f_init,"%d", OutletProfile);      // outlet profile (yes/no)
  fscanf(f_init,"%d", Iterations);         // # of iterations
  fscanf(f_init,"%d", AutosaveEvery);      // autosave every #th of iteration
  fscanf(f_init,"%d", AutosaveAfter);      // autosave after the #th iteration
  fscanf(f_init,"%d", PostprocProg);       // program of postp rocessing (Paraview/TECplot)
  fscanf(f_init,"%d", CalculateDragLift);  // calculate drag & lift? if > 0 than on which BC_ID
  fclose(f_init);                          // close the file
}
