// Test_LimitcycleDetect.h

void BinaryFindtheFixedPoint_new(double v1[2], double v2[2], int type, double FixPoint[2], int &triangleID);
void GettheClosedStreamline_new(int repellorattract, double FixPoint[2], int &triangleID);
bool CloseToSingularity(double point[2], int singID);
bool IsASingularity(double alpha[3], int triangle);
void ChooseBeginningPointForRegularSing(int singID, double begin_p[2], int &begin_triangle);
void GetBeginPoint(int singID, int type, double begin_p[2], int &begin_triangle, int Saddleornot, icVector2 sep_vec);
void GetNextLevelBeginPoint(int singID, double begin_p[2], int type);
void UpdateListInSingularity(int singID, int limitcycle);
void BuildHandlerforLimitCycle(int cur_limitcycle_index);
void DetectFromASingularity(int singID, int type, int Saddleornot);


///////////////////////////////////////

bool IntheRegion(int triangle, int type);
void GetTheIntersectRegion();
void GetAShareEdge_new();
void DetectALimitCycle(double begin_p[2], int begin_triangle, int type, int inner);
void GetCenterofATriangle(int triangle, double center[2]);
bool IsAtFieldBoundary(int triangle);
void GetACellCycle_new(int singID, int type);
void DetectLimitCycle_new();

