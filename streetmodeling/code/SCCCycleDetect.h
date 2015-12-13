/*
SCCCyclceDetect.h
The head file for the strongly connected component based limit cycle detection 
*/

bool IsValidSCC(int);

void GetSpPtsForValidSCC();

void SCCCycleDetect();

void RegularDetect(int, int);

void SpecialPtsBasedDetect(int, int);


bool IsPreviousCycle_SCC(int scc_index, int, int &);

bool IsPreviousCycle_all(int triangle, int &pre_limit_index, int type);

int *GetDisc(double p[3], int triangle, double dsep, double discsize,
			 int *NearbyTriangles, int &num_triangles);

void DetectLimitCycle_SCC();


void ResetEdgeIntersections(int scc_index);


/*******************************************************************************/
//Routines for building the connections between limit cycles and other elements

int TraceInTriangleForConnection(double g[2], int &face_id, int type, int &flag);
void TraceandBuildConnection(double s[2], int triangle, int singID, int type);
void TraceForConnection_nonsaddle(int singID);
void TraceForConnection_saddle(int saddle);
void ConnectionForSingCycle();

void TraceForConnection_cycle(int cycle, int connected_side);
void TraceandBuildConnection_cycle(double s[2], int triangle, int cycleID, int type);
void ConnectionForCyclePair();

void BuildConnectionForCycle();

