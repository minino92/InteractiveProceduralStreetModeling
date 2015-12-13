#include "lib/nr.h"
#include "lib/icVector.h"
#include "lib/icMatrix.h"

/*----------------------------------------------------------------------*/
void AllocateVarforTopologyEdit();
void InitCancellationAndMovement();
void FinalizeTopologyEdit();

void BuildTheSparseLinearSystem2(Vec_DP &tsa, Vec_INT &tija, Mat_DP &bound_v);
void Cancel_RegionSmooth();

void AddToRegionTriangles(int triangleID, int type);

void UpdateBoundary(int type);
void GetRegionNormals(int type);

//void Cancel_GetRegion(int repell_triangle, int attract_triangle, int repellID, int attractID);
//void Cancel_GrowRepellerRegion(int which_triangle, int target_triangle, int singularID, int initsaddlelength);
//void Cancel_GrowAttractorRegion(int which_triangle, int target_triangle, int singularID, int initsaddlelength);
void Cancel_GetIntersectedRegion(int repellID, int attractID);
void Cancel_Growing(int type, int target_triangle);
void Cancel_InitSaddleRegion(int singularID, int type, double percentage);
bool AdaptivePairCancel(int repell_triangle, int attract_triangle, int repellID, int attractID);
void PairCancellation(int, int);
bool IntersectedRegionSingCapture();

void Cancel_GrowAttractorRegion(int which_triangle, int target_triangle, int singularID, int initsaddlelength,
								int MultiCancel, double length_percentage);
void Cancel_GrowRepellerRegion(int which_triangle, int target_triangle, int singularID, int initsaddlelength,
								int MultiCancel, double length_percentage);
bool AdaptiveGetMultRegion_new(int *repellers, int num_repellers, int *attractors, int num_attractors);

/*---------------------------------------------------------------*/
////11/06/05
//void AdaptivePairCancel_adv(int, int, int, int);
//void AdaptivePairCancel_adv2(int repeller, int attractor/*, int saddle1, int saddle2*/); //11/19/05
void AllocateVarforMultiRegion();
void FianlizeMultiRegion();

void InitMultiRegionVars();
void Cancel_GetIntersectedRegion_adv(int repellID, int attractID);
void Cancel_GetRegion_adv(int repell_triangle, int attract_triangle, int repellID, int attractID);
void Cancel_GrowRepellerRegion_adv(int which_triangle, int target_triangle, int singularID,
								   double length_percentage);

void Cancel_GrowAttractorRegion_adv(int which_triangle, int target_triangle, int singularID, 
									double length_percentage);

void CopyRegion(int*, int*, int);
bool IstheMediaNode(int);

void AdaptiveRegionGenerate();   //for adaptive region generation of a pair of repeller and attractor
void PairCancelTwoSameType();    //for pair cancellation of two same type of nodes

//// 11/20/05
////Get region for multiple repellers and attractors
void GetMultRegion(int *repellers, int num_repellers, int *attractors, int num_attractors);
void GetMultRegion2(int *repellers, int num_repellers, int *attractors, int num_attractors);

bool IsRepeated(int *a, int b, int num);
void AddToRegionforMultRegion(int triangleID, int type, int &num);
void InitSaddleGrowforMultRegion(int singularID, int type, 
									double length_percentage,
								     int *sametypesings, int numsings);
void InitSaddleGrowforMultRegion_new(int singularID, int type, 
									 int initsaddlelength, double length_percentage,
								     int *sametypesings, int numsings);

void GetMultIntersection(int *repellers, int num_repellers, int *attractors, int num_attractors);

void MultCancellation(int *repellers, int num_repellers, int *attractors, int num_attractor); //main routine for multiple cancellation
//int IntersectedRegionSingCount();
int IntersectedRegionSingCount(int &totalindex);

int CountNumTriangleonSeparatrix(int index);

/**---------Routines for setting fence for region growing---------**/
////2/16/06

void SetFences(int *repellers, int num_repellers,
			   int *attractors, int num_attractors,
			   int *intervals, int num_intervals);
void SetFenceForASep(int index);
void SetFenceForOneSaddleandOtherSingCancel(int saddleNode, int singNode);
bool SepoftheSaddleIntheList(int *intervals, int num_intervals, int sepindex);
void SetFenceForMultiPairCancel(int *repellers, int num_repellers,
			   int *attractors, int num_attractors,
			   int *intervals, int num_intervals);
void SetConnectedSeps(int *repellers, int num_repellers,
			   int *attractors, int num_attractors,
			   int *intervals, int num_intervals);

void ClearAllFences();

/*---------------------------------------------------------------*/
////11/06/05
void PairCancelThroughGraph(int node1, int node2); //we may perform pair cancellation through selecting nodes in the graph
/*---------------------------------------------------------------*/


void Move_Growing();     ////growing the source(repeller) region for movement
void Move_GetRegion(int triangleID, int singID, int target_triangle, double newx, double newy);
void Move_SetVectorsAtNewPos(int, int);
void Move_GrowRepellerRegion(int triangleID, int target_triangle);
void Move_GrowAttractorRegion(int triangleID, int target_triangle, int singularID, double newx, double newy);
void Move_GetIntersectedRegion(int sourceID, int target_triangle);

////This routine calculate the incoming/outgoing eigen vector
void Move_GetVirtualSaddle(int old_singID, double x, double y, icVector2 &in, icVector2 &out);  
void Move_InitSaddleRegion(int triangleID, double newx, double newy);
void Move_GetInitRegionofVirtualSaddle(int Init_triangle, icVector2 out, double newx, double newy);


void TestingRoutine(int, int);
void TestDisplayRegion(unsigned int mode);

void Move_TestingRoutine(int, int, int, double, double);

/////Routine for smoothing region validation
int CalEulerValue(int *trianglelist, int num);

/////Routines for automatic rotation
void RotateFieldBasedOnSing(int sing_index);
void CompensateRotate(icMatrix3x3 compense_matrix);

///////////////

int *Extend_link(int *edge_link, int Num_edges);


//////////////////////////////////////////////////////////////
////For two connection orbit cases
bool TwoConnectOrbits(int saddle, int singID); //judge whether there are two connections orbits
icVector2 GetSaddleTriangleDirection(int saddle, int traj);
void SetExtraBoundary(int saddle, int singID, icVector2 saddle_dir, int traj);
void SetExtraBoundary2(int saddle, int singID, icVector2 saddle_dir, int traj);

void PairCancelWithTwoConnections(int saddle, int singID);
void InitSaddleRegionForTwoConnections(int singularID, int type, int connect_sing, int initsaddlelength);
void GrowRepellRegionForTwoConnections(int singularID, int type);
void GrowAttractRegionForTwoConnections(int singularID, int type);
void GetIntersectRegionForTwoConnections(int repellID, int attractID);

void SetBoundaryTriangles();
void ExtendIniStrip(int traj, double length);
void SetExtraBoundary3(int saddle, int singID, icVector2 saddle_dir, int traj);

///Temporary routine to get one of the two connection orbits for testing 2/23/06
int GetOneConnection(int saddle, int singID);

bool PeriodicDetect_Growing(int type);
bool PeriodicDetect_Growing(int type, int centersing, int &flag);
