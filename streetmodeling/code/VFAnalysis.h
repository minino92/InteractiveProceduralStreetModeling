
////VFAnalysis.h

/*This file provide routines fulfilling the analysis of current vector field (generally based on the vector field
stored in 2D local frame)

The inputs of this module are underneath mesh, vector field stored at each vertex
The output is the results, such as singularities, limit cycles, separatrices

Most of the calculations are performed here, especially local tracing calculation for trajectory
*/


////
//The followings are some possible useful routines from previous program

////singularities extraction
void CaptureSing(void);
void ComputeSing(void);
int GetSingType(int MarkTriangleID);

int get_singularity_type(int MarkTriangleID);
void get_loc_Jacobian(int tri, double Jacobian[2][2], double &c, double &f);
int get_singularity_type_loc(int tri, double Jacobian[2][2], double x_loc, double y_loc, double alpha[3]);
void compute_fixedpt_info(void);

////functions for calculate trajectory
void CalSingleTrajectory(double x, double y);
void CalSingleSeparatrix(double x, double y, int inout);
void CalSeparatrices();


////limit cycle detection







//////functions for matrix inversion
////These routines should be moved to another library
double *MatrixOpp(double A[],int m,int n); //inverse 
double *MatrixInver(double A[],int m,int n); //zhuanzhi
double Surplus(double A[],int m,int n); //hanglieshi


double cal_determinant2x2(double a[][2]);
bool cal_inverse2x2(double a[][2], double inverse_a[][2]);
//void cal_eigenvectors(icMatrix2x2 mat, double evalues[2], icVector2 ev[2]);
