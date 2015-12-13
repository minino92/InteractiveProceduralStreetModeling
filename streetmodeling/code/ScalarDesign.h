/*
scalardesign.h
the head file of the scalar field creation
*/
#include "lib/icMatrix.h"
#include "lib/icVector.h"

void cal_eigen_vector_asym(icMatrix2x2 m, double phi, icVector2 ev[2]);
void cal_eigenvecs_onevert_quad_asym(int ver);
void vis_scalarfield();
double cal_phi_at(double x, double y);
void cal_phi_allverts();
void pro_phi_allverts();
void normalize_scalar_phi();
void cal_all_eigenvecs_asym();
void compute_phi_in_quad(int face, double x, double y, double &phi);

void init_scalar_singular_elemlist();
void add_new_scalar_singular_elem(double x, double y, unsigned char type);
void test_scalar_field();