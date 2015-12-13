
#include "lib/nr.h"
#include "lib/icVector.h"



/*----------------------------------------------------------------------*/
////Functions for region smoothing

void InitRegionSmooth();
void FinalizeSmoothing();
void AddToPointList(double, double, int);
void EndPointList();

double GetLength(int ver1_id, int ver2_id);
void RecursiveDijistra(int s_id, int d_id, int cur_v_id, double dis_to_s);

void FindOutInnerVerts();
bool InRegion(int vert_id);           ////judge wether a vertex is inside the region
double GetDirectionalAngBetween2Vec(icVector2, icVector2);

int GetFirstInsideVertex();

void BuildTheSparseLinearSystem(Vec_DP &tsa, Vec_INT &tija, Mat_DP &bound_v);
void BubbleSorting(int *a, int size);
void RegionSmooth();
void ResetSmoothVars();

void Undo();
void Redo();

void SeedSearch(int seed);
void NonRecursiveSearch(int seed);

void AllocateVarforSmoothing();

void SavePostField();
void SavePreviousField();

