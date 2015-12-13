////LocalTracing.h
#include "lib/icVector.h"
#include "lib/nr.h"

#define TRACESTEPS 700



/*--------------------------------------------------------------------------------------*/

////Routines and variables for local tracing 05/24/05
void CalLocalTracing(int face_id, double x, double y, int type);     //Drawing trajectory in local linear field

int TraceInATriangle(int &face_id, double globalp[2], int type, int &flag);

int TraceInATriangle2(int &face_id, double globalp[2], int type, int &flag);

int TraceParticleforOneStep(double &x, double &y, int &face_id, int type);

void Get2DBarycentricFacters(int face_id, double a, double b, double alpha[3]);

icVector2 GetVectorAtPoints(int face_id, double alpha[3]);

//bool StoreToGlobalList(CurvePoints *temp, int num);

bool ToNextPoint(double first[2], double second[2], int &face_id, double alpha[3], int type);
bool ToNextPoint2(double first[2], double second[2], int &face_id, double alpha[3], int type);
bool ToNextPoint3(double first[2], double second[2], int &face_id, double alpha[3], int type); //2/19/06
bool Euler_ToNextPoint(double first[2], double second[2], int &face_id, double alpha[3], int type);
bool get_nextpt_euler(double first[2], double second[2], int &face_id, double alpha[3], int type);
bool get_nextpt_2ndeuler(double first[2], double second[2], int &face_id, double alpha[3], int type);
bool compute_next_pt(double first[2], double second[2], int &face_id, double alpha[3], int type);
double get_shortestedge_tri(int tri);
bool get_nextpt_2ndeuler_2(double first[2], double second[2], int &face_id, double alpha[3], int type);
bool get_nextpt_RK4_adp(double first[2], double second[2], int &face_id, double alpha[3], int type);
bool get_nextpt_RK4(double first[2], double second[2], int &face_id, double alpha[3], int type);


void GetNextTriangle(int &face_id, double pre[2], double cur[2], double param_t[2], int type, 
					 int &PassVertornot, double alpha[3]);

void TriangleThroughVertex(int vert_id, int &theone, int type);

void localderive(const DP t, Vec_I_DP &y, Vec_O_DP &dydx);

void localinverse_derive(const DP t, Vec_I_DP &y, Vec_O_DP &dydx);

int TriangleDetect(double x, double y);

void  CrossBoundary(int prev_id,  double point2D[2], int &which_edge);
void  CrossBoundary3(double pre[2], double cur[2], int face_id, double alpha[3], 
					 int &which_edge, double t[2]); //4/29/06

void CrossVertex(int &face_id, double cur_p[2], double pre_p[2], int type, int &passornot);
void CrossVertex2(int &face_id, double cur_p[2], double pre_p[2], int type, int &passornot); //4/30/06

void PassEdge(int &face_id, int which_edge);

void GetIntersection(double p1[2], double p2[2], double q1[2], double q2[2], double t[2]);
int GetIntersection2(double PointA[2], double PointB[2], double PointC[2], double PointD[2], double t[2]);
extern int cal_majIntersect(double PointA[2], double PointB[2], double PointC[2], 
					 double PointD[2], double t[2]);
int cal_intersect_2(double PointA[2], double PointB[2], double PointC[2], double PointD[2], double t[2]);
int cal_intersect(double PointA[2], double PointB[2], double PointC[2], double PointD[2], double t[2]);

void CalSingleSeparatrix(int Triangle_ID, double x, double y, int inout);  ////for calculate a single separatrix

void CalSeparatrices();

bool LoopDetect(int *triangles, int num, int oneTriangle);

bool FallIntheTriangle(int triangleID, double pos[2]);

void FixSepBeginningPos(int &triangleID, icVector2 sep_vector, double sing_center[2], double newpos[2]);

/*-------------------------------------------------------------------------------------*/


double GetSepLength(int trajID);
