////
////LimitCycleCreator.h  

/*
Method1: Use the control points to calculate the tangent vectors on them, then use the 
        convergent or divergent elements to form a limit cycle
Method2:Use curve design and constrained optimization method to design a limit cycle
*/

#include "lib/icVector.h"

typedef struct LocalPts{
	double x[2];
}LocalPts;

typedef struct DesignTriangleCycle{
	int *DesignCellCycle;
	int num_triangles_designcurve;
	icVector2 *Appro_dir;
	LocalPts *bases;
	int MaxNumTrianglesDesignCurve;	
}DesignTriangleCycle;

///////////////////////////////////////////////////
////Limit cycle shape design
void AllocShapeDesignVars();
void FinalizeShapeDesign();
void AddToShapeCtrPtsList(double x, double y);
void CalHermiteCurve();
void CalOpenHermiteCurve();
void GenerateLimitCycle(int type);

void BuildSmoothRegion();

void NormalizeVectorsOnBoundary();
void InitLimitCycleShapeDesign();
void ResetLimitCycleShapeDesign();

void GetNormalsForBoundaryEdges(int type);

void GetCellCycleforDesignCurve();
void BuildBoundaryEdgeList(int *DesignCellCycle, int num_celltriangle);
void BuildBoundaryEdgeList(int *region, int num, int type);
void GetBoundary();
void LimitCycleGeneration(int type);

bool ShareTheSameVertex(int t1, int t2, int &v); ////added on 2/21/06
////Routines for new limit cycle shape design methods
void GetOrientationofBoundaries(int &same_orient, int &clockwise);
void GetVectorsforAllBoundaryVerts(int same_orient, int clockwise, int type);
void GetVectorsOnBoundaries_new(int type);
void GetVectorsOnBoundaries(int type);

void GetDesignCellCycle_new();   ////new method to get design cell cycle
void GetNextTriangle_new(int &face_id, double pre_p[2], double cur_p[2], int &Passvertornot, double alpha[3]);
void GetNextTriangle_new2(int &face_id, double pre_p[2], double cur_p[2], 
						  int *designcellcycle, int &num_cells);
void CrossVert_new(int &face_id, double cur_p[2], double pre_p[2],int &passornot);

void ExtendDesignCellCycle();

void GetDesignCellCycle_new2();  ////use edges to simulate the curve 2/5/06


//////////////////////////////////////////////
////Save and set previous boundary vertices
////2/20/06
void SaveBoundaryVerts();
void SetBoundaryVerts();

void ResetCurBoundary();

int GetResolution();
double GetWholeLengthofDesignCurve();

bool InBoundary(int triangleid); //3/29/06


////Save the design curve for debugging
void SaveDesignCurve();
void LoadDesignCurve();

//#define	_DEBUG
//
//#ifdef _DEBUG
//	fprintf(stderr, "asdfasdf");
//#endif

/////
int GetOneTriangle(int &pos, double globalp[2], int &cur_triangle, int num_curvepts, 
				   int *DesignCurveCellCycle, int &num_triangles_designcurve);
int GetOneTriangle2(int &pos, double globalp[2], int &cur_triangle, int num_curvepts, 
				   int *DesignCurveCellCycle, int &num_triangles_designcurve, int &Passvertornot);

void ConnectTwoTriangles2(int triangle1, int triangle2, int SingleSharingVert, 
						  int *DesignCurveCellCycle, int &num_triangles_designcurve);

icVector2 GlobalToLocal_Vector(icVector2 glob_vec, int triangle);
void GlobalToLocal_Pt(double gp[2], int triangle, double lp[2]);

void GetDesignCellCycle_new4();
void GetBoundaryVerts(int *DesignCellCycle, int num_triangles_designcurve);
void GetVectorsOnBoundaries(int *Boundaryverts, int numverts, int type);

icVector2 GetALocalVec(double loc_p[2], icVector2 loc_vec, double loc_base[2], int type);
void LimitCycleGeneration2(int type);
bool CreateALimitCycle(/*int *trianglelist, int num*/);
void UnsetTriangleList(int *trianglelist, int num);
