
#include "lib/icMatrix.h"
#include "lib/icVector.h"

void get_deformation_tensor(icVector2 x[3], icVector2 x_img[3], icMatrix2x2 mat);
void get_deformation_tensor(double vx[3], double vy[3], 
							 double vx_img[3], double vy_img[3], 
							 double mat[2][2]);
int get_type_deformation(double mat[2][2]);
int get_type_deformation(icMatrix2x2 mat);
void get_vf_formula_tri(int tri, double Jacobian[2][2], double &c, double &f);
void cal_eigenvectors(icMatrix2x2 mat, double evalues[2], icVector2 ev[2]);
int cal_eigenval_matrix2x2(icMatrix2x2 mat, double evalues[2]);
 void cal_Jacobian_All_Vers();