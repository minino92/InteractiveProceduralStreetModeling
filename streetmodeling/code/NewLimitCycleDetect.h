//NewLimitCycleDetect.h  1/24/06


void DetectLimitCycle_new2();

void GetACellCycle_new2(int singID, int type);

void CellCycleDetect_new(double &x, double &y, int &begin_triangle, int type, int &flag);

void DetectALimitCycle_new(double begin_p[2], int begin_triangle, int type, int inner, int singID, int &flag);
int GetEdgeIndexofTriangle(int triangle1, int triangle2);

int FindTheMostAccessTriangle();

int OuterBoundaryTriangle(int type, int &anintriangle);

//////////////////////////////////////////////////////////
////Using parallel tracing to locate a cell cycle 3/22/06
bool ParallelCellCycleLocate(int out_t, double out_gp[2],
							 int in_t, double in_gp[2],
							 int *cellcycle, int &num_cells, int type);

bool IsSameIntArray(int *a1, int num_a1, int *a2, int num_a2); //judge whether a two integer arrayes are the same (not considering the order)

bool IsSimilarIntArray(int *a1, int num_a1, int *a2, int num_a2, double percent);

void GetALimitCycle(int out_t, double out_gp[2],
					int in_t, double in_gp[2],
					int *cellcycle, int &num_cells, 
					int singID, int type, int inner);


//////////////////////////
void GetInitTriangleStrip(int type);



/////New boundary extraction 4/11/06 these following routines can be put to a library
////Need to be verified!!! 5/22/06
void AllocBoundaryList();
void GetAllBoundaries(int *region, int &num);
void GetACellCycle_3(int singID, int type);  //Using new method to determine the embedded periodic orbits
void ShowAllBoundaries();
int GetAnOuterBoundaryTriangle(int *regiontriangles, int nums);

bool ReachMeshBoundary();




/////Revisit the limit cycle detection 5/22/06

void GetCellCycles(int singID, int type);
void GetBackwardRegion(int *region, int num, int type, int *backwardregion, int num_backward);
bool GetOneRingShapedRegion2(int *ring, int &num, int type); //get the region that grow from one boundary
bool GetOneRingShapedRegion(int *ring, int &num, int type); //get the intersection of all the regions
bool FindALimitCycle(int *trianglestrip, int num_triangles, double bgp[2], int begin_triangle, 
					 int inner, int singID, int type);
bool IsPreviousDetectedLimitCycle(int *cycle, int num, int &pre_limitcycle_id);
bool OneCellCycleDetect(int *trianglestrip, int &num_triangles, double bgp[2], int &begin_triangle, int type);
bool FindFixedPoint(int *cellcycle, int num_cells, int type, double fixedpoint[2], 
					int &fix_triangle, int &OneEndVert);
void BuildHandlerforLimitCycle(int cur_limitcycle_index, double fixedpoint[2], int ver);

//Edge *GetAGoodEdge(int *cellcycle, int num_cells);


////A general integer processing routine
bool RemoveOneElem(int *a, int b, int &num);