/*
tensordesign.h
the head file of the tensor field creation
*/
#pragma once

void init_ten_designelems();
void init_regular_ten_designelems();


void init_tenelem_EditBox(int index, double x, double y);

void addto(double x, double y, int triangle, int type);

void get_tensor(double x, double y, double t[4]);

void cal_alltensor();
void cal_alltensor_quad();

/*region based tensor field design*/
void cal_tensorvals_quad_inReg(/*int regionid*/);
void cal_tensorvals_quad_inReg(int regionid);

void get_tensor_inReg(double x, double y, double t[4], int regionid, bool);
void cal_eigenvecs_quad_inReg(int regionid);


void update_tenSingularElem_transform(int TransformType, int elem_id);

void update_tenElem_EditBox(int elem_id);

void update_tenElem_Transform(int TransformType, int choose_ID);

void set_ten_regBasis(double x, double y, int type);

void set_ten_regDir(double x, double y);

void update_ten_regDir(int id, double x, double y);

void cal_eigen_vector(icMatrix2x2 m, icVector2 ev[2]);

void cal_sym_eigenvecs_vertices();

void save_tenField_perVer(char *filename);
void load_tenField_perVer(char *filename);

bool save_tenField_elems(char *filename);
bool load_tenField_elems(char *filename);

void save_cur_field();

void init_pre_ten();
void store_cur_ten();

