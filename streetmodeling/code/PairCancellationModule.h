/*  PairCancellationModule.h */

/* function declaration */
void Gen_PairCancellation(int repeller, int attractor, int type);

void SingularityPairCancel(int, int);
void LimitCycleSingularityCancel(int, int);
void LimitCyclePairCancel(int, int);

bool IsDirectlyConnected(int node1, int node2); //judge whether two nodes are directly connected in C-graph or not


/* functions for different kinds of singularity pair cancellations */
void DirectlyConnectedSingularityPair(int, int);
void DoublyConnectedSingularityPair(int, int);
void SimplyConnectedSingularityPair(int, int);

void IndirectlyConnectedSingularityPair(int, int);
void IndirectlyDoublyConnectedSingularityPair(int, int, int *, int);
void IndirectlySimplyConnectedSingularityPair(int, int, int *, int); //simply indirectly connected pair

void GrowAttractRegionForTwoConnections_test(int singularID, int target_triangle, double percentage);
void GrowRepellRegionForTwoConnections_test(int singularID, int target_triangle, double percentage);
void InitSaddleRegionForTwoConnections_test(int saddle, int type, int connect_sing, 
											int initsaddlelength, double percentage);

/* functions for different kinds of singularity and limit cycle pair cancellation */
bool IndirectlyConnectedSingandCyclePair(int sing, int cycle);
bool Ada_Getregion_IndirectlyConnectedSingandCyclePair(int sing, int cycle, int *Media, int num_media);
void InitSaddle_IndirectlyConnectedSingandCycle(int saddle, int type, int sing, int cycle, double percentage);
void RemoveFenceFromSepToCycle(int cycle);
int GetMinTriangleNumofSeps(int *Media, int num, int sing, int cycle);

bool DirectlyConnectedSingandCyclePair(int sing, int cycle);

/* functions for sink and source cancellation with double connections */
void InitSaddle_IndirectlyDoublyConnectedSingularityPair(int saddle, int sing1, int sing2,
													int type, double percentage);
bool Ada_Getregion_IndirectlyDoublyConnectedSingularityPair(int sing1, int sing2, 
													   int inter_saddle);
void IndirectlyDoublyConnectedSingularityPair(int sing1, int sing2, int *Media, int num_media);

int GetMinTriangleNumofSepsForDoublyConnectedSingularityPair(int saddle, int sing1, int sing2);

void RemoveFenceForIndirectlyDoublyPair(int, int);

void SetFenceForSeps(int except_sa);
void RemoveLimitCycleFence(int limitID);
void RemoveFenceForSaddle(int saddle);
