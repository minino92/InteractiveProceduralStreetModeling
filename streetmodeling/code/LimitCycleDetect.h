////LimitCycleDetect.h

////This module implements the detection of limit cycle using the closed streamline detection method

/*-------------------------------------------------------------------------------------*/
////Routines for global limit cycle detection 05/27/05
////it will be put to limit cycle detector file
////Note all the following routines should be compatible with local tracing
////especial TriangleDetect
void InitLimitCycleDetection();            ////initialize the limit cycle detection
//void InitLimitCycleStructure();            ////initialize the limit cycle data structure
void ResetLimitCycles();

void finalizeLimitCycleDetect();           ////release memory allocate to variables of limit cycle detection
void LimitCycleDetect();                   ////begin from the capture singularities

bool TriangleSearch(int *acycle, int num_cycletriangles, int oneTriangle, int &position);
bool FindCellCycle(int *trilist, int originNum, int *acycle, int &CellNum, int cur_triangleID); //search the whole triangle list
	                                                                //with current triangle ID to find a cycle
bool closedStreamlineTest(int repellorattract);                  //closed streamline testing
bool backwardExitTest(double x, double y, int repellorattract); //test one potential exit
void boundaryBuilding();
void PotentialVertsBuilding();
void AddToPotentialVertList(int cur_v_id);

void GetNextTestBeginPoint(double &bx, double &by);

void GettheClosedStreamline(int repellorattract); ////after we locate a cell cycle, using binary search to find the 
	                                        ////accurate closed streamline

void GetNextLevelBeginPoint(double &next_x, double &next_y, int limitcycle_index, int &flag);

////routine for cell cycle detection using local tracing
void local_CellCycleDetect(double &x, double &y, int type, int &flag);
bool local_TravelCellCycleTest(double x, double y, int repellorattract);  //close streamline testing using local tracing

////
void Double_CellCycleDetect(double &x, double &y, int &begin_triangle, int type, int &flag);
void  First_CellCycleDetect(double &x, double &y, int &begin_triangle, 
							int *acycle, int &num_localcell, int type, int &flag);
bool Second_CellCycleTest(double &x, double &y, int type, int *localcycle, int num_localcell, int &flag);

////
void  First_CellCycleDetect_2(double &x, double &y, int &begin_triangle, 
							int *acycle, int &num_localcell, int type, int &flag);
bool  Second_CellCycleTest_2(double &x, double &y, int type, int *localcycle, 
						    int num_localcell, int &cur_Triangle, int &flag);
void  Double_CellCycleDetect_2(double &x, double &y,  int &begin_triangle, int type, int &flag);


/*--------------------------------*/
////routines for locate the closed streamline accurately
void BinaryFindtheFixedPoint(double v1[2], double v2[2], int type);
//void local_GetNextIntersection(double beginx, double beginy, double &interx, double &intery, int type, int &flag);
void  local_GetNextIntersection(double beginx, double beginy, int thetriangle,
								double &interx, double &intery, 
								int type, int &flag);

void StoreCurrentCellCycle(int *acellcycle, int num);
void GetASharedEdge(double center[2]);

/*--------------------------------*/
////routines for one step Runge Kutta
void OneBackwardRungeKutta(double x, double y, double &nextx, double &nexty, int repellorattract);

////Initialize the variables for limit cycle detection
void AllocVarforLimitCycleDetect();

////Check whether current detected cell cycle is the same of one of previous limit cycles
bool IsPreviousLimitCycle(int &limit_index);

/*---------------------------------------------------------------------*/
////routines for new method of limit cycle detection
void GetLimitCycleRegion(int singID, int type, int inner);

void LimitCycleDetectNew();  ////using region growing to locate limit cycle

/*---------------------------------------------------------------------*/
////Update the corresponding information for the generation of Conley relation graph
void UpdateListInSaddle(int saddleID, int limitcycle);
void UpdateListInLimitCycle(int limitcycle, int saddleID);
void UpdateCycleListInLimitCycle(int limitcycle1, int limitcycle2);
