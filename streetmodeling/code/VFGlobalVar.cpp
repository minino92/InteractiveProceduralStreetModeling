#include "VFDataStructure.h"

////VFGlobalVar.cpp
/* Define the global variables may be used in the whole program!!
*/


/*-----------------------------------------------------------------------*/
/////Define global vector field related variables here 06/23/05
////Global variables

int MaxNumSingularElems;                 //Maximum number of singular elements
int MaxNumRegularElems;                  //Maximum number of regular elements
int MaxNumSingularities;                 //Maximum number of being captured singularities
int MaxNumTrajectories;                  //Maximum number of possible trajectories
                                         //(it should be flexible for future pen-and-ink sketch)
int MaxNumLinesegsPerTraj;               //Maximum number of line segments for each trajectory

SingularElement *singularelem;          //Singular elememts' list
int cur_singelem_index;
RegularElement *regularelem;            //regular elememts' list
int cur_regelem_index;
Singularities *singularities;           //being captured singularites' list
int cur_singularity_index;
LineSeg **trajectories;                 //trajectories' list
int cur_traj_index;
int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory


////Define global variables for underneath mesh

Polygon3D Object;

void nullfun()
{
}