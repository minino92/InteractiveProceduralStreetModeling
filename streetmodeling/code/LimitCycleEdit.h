
//////LimitCycleEdit.h

void InnerLimitCycleSearch(int seed);
int GetBeginningVert(int limitcycle_index);
bool GetInnerLimitCycleRegion(int limitcycle_index, int *trianglelist, int num);
//void CancelLimitCycle(int index);
//void CancelPairLimitCycles(int inner, int outer);
//void CancelOneLimitCycle(int index);

bool LimitSaddleCancel(int saddle, int limitID);


////routines for limit cycle and singularity pair cancellation
void CancelLimitSingPair(int singID, int limitID);
//void GetLimitSaddlePairRegion(int saddleID, int limitID);
//void GetLimitSaddlePairRegion2(int saddleID, int limitID, int singID);
void GetInnerVerts(int singID, int *Medianodes, int Num_medianodes);

void Cancel_InitLimitSaddleRegion(int singularID, int type, int initsaddlelength);
void Cancel_InitLimitSaddleRegion2(int singularID, int type, int singID);
void InitSaddleRegionforLimitPairCancel(int saddle, int type, int initlength);

bool SaddleSingConnectionJudge(int singID, int saddleID);
//int FindTheConnectedSaddle(int singID);
bool FindTheConnectedSaddle(int singID, int &saddle1, int &saddle2);
//void GrowInitLimitRegion(int SaddleID, int type);
void GrowInitLimitRegion(int SaddleID, int otherSaddleID, int type);
void Region_Growing(int type);
//int GetLimitRegionFromSaddle(int singID, int type);
bool GetLimitRegionFromSaddle(int singID, int type);
void GetLimitSingSmoothingRegion(int singID, int limitID);
void GetLimitSingSmoothingRegion2(int singID, int limitID, int *MediaNodes, int Num_medianodes);

void GetLimitSingIntersectRegion(int singID/*, int saddleID*/);
void CancelEmbededLimitCycles(int inner, int outer);
void GetLimitPairIntersectRegion();

void GrowRegionforALimitCycle(int limitID);
void GrowAttractLimitCycle_Region(int limitID);
void GrowRepellLimitCycle_Region(int limitID);


/////Routines for limit cycle relocation 3/8/06
void GetALimitCycleRegion(int limitID);
double GetTheLengthofWholeCycle(int limitID);
double GetIntervalLength(double whole_len, int num_controlpts);

//void SampleControPts(int num_controlpts, int limitID, ctr_point *control_pts);
//
//void GetOneSamplePoint(int limitID, int cur_lineindex, 
//					   ctr_point prectrlpt, ctr_point &curctrlpt, 
//					   int &new_lineindex, double interval);


void SetRingRegionBoundary(int type);
void BuildSmoothRingRegion(int type);
void GeneNewLimitCycle(int type);


////New routines for assigning vector values on the boundary of the new limit cycle
double GetDistance(double x1, double y1,
				   double x2, double y2,
				   double outx, double outy);
void FindAllDistance();
void FindTheDisForOneVert(int *triangles, int num_triangles, int vertid);
void ResetVectorOnBoundary(int type);


///////Routines for cancellation of limit cycle and saddle with two connecting orbits
////3/11/06

bool IsTwoConnections(int limitID, int saddle); //judge whether there are two connecting orbits between them

void LimitSaddleCancelWithTwoOrbits(int limitID, int saddle);
void InitSaddleRegionForTwoOrbits(int limitID, int saddle, double initsaddlelength);
void GrowSaddleRegionForTwoOrbits(int limitID, int saddle);
void GetIntersectRegionForLimitSaddleWithTwoOrbits(int limitID, int saddle);

void SetFenceForALimitCycle(int limitID);

void SetExtraBoundary(int limitID, int saddle, int traj); //set extra boundary??

///Temporary routine to get one of the two connection orbits for testing 2/23/06
int GetOneConnection(int limitID, int saddle, int type);
void SetFence_LimitCycles(int inner, int outer, int *MediaNodes, int Num_MediaNodes);
void SetFence_LimitCyclesexcept(int limitID);
