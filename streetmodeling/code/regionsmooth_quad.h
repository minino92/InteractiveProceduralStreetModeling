/*
regionsmooth_quad.h
*/

void alloc_quad_regionsmooth();
void init_quad_regionsmooth();


int get_firstInsideVertex();
bool is_inregion(int vert_id);
bool is_inregion(double x, double y);
bool is_inregion_appro(double x, double y, double threshold);
void mark_boundVerts();
void search_innerVerts_nonrecursive(int seed);
void find_innerVerts();
void find_boundarycells_oneline(double start[2], int start_cell, double end[2],
								int end_cell, int *celllist, int &ncells);
void find_boundarycells();
void display_innerverts();

void smooth_Jac_quadregion();
