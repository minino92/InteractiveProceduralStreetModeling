
/*
     LocalRegionEdit.h
*/

#pragma once

void find_inner_cells();
void find_inner_intersections();
void recal_intersections_on_RegBoundary();
void convert_to_tensorline();
void search_for_intersection_along_edge(int edgeid,  double intersect[2], int &which_cell);
void remove_street_in_reg();
void mark_inner_verts_with_new_RegIndex();
void update_RegIndex_inner_and_boundary_cells();
void replace_inReg();
void compute_intersections_inReg();
void construct_sub_graph_inReg();
void search_connection_inReg();
void update_street_network();

int compute_contour_in_one_cell(int cellid, double dist, double from[2], int &from_edge,
								 double to[2], int &to_edge);
void get_sorted_contour_cells();
void mark_contour_cells(double dist);
bool is_contour_cell(int id, double dist);
void compute_the_contour(double dist);
void compute_level_set_contour(double dist);

bool cal_intersect_at_graph_edge(int edgeid, double intersect[2], int &which_samp,
								 double A[2], double B[2]);
void recal_intersections_on_RegBoundary_2();
void recal_intersections_on_RegBoundary_3();

bool cal_smallest_dist_from_one_cell(int cellid, double A[2], icVector2 dir,
									 double &smallestdist, int &which_intersect, int except_intersect);
double cal_smallest_dist_to_one_edge(int edgeid, double p[2], int &which_samp);
bool cal_smallest_dist_to_intersects(int cellid, double A[2], double ang, icVector2 dir, 
									 double &dist, int &which_intersect, int except_intersect);
bool cal_smallest_dist_to_edges(int &cellid, double A[2], double ang, icVector2 dir,
									 double &dist, int &which_edge, int &which_samp, int except_edge);
bool cal_smallest_dist_to_one_edge_in_a_cell(int cellid, double A[2], icVector2 dir, double &dist,
											   int &which_edge, int &which_samp, int except_edge);
//void merge_close_intersections_on_RegBoundary(double threshold);
void merge_subgraph_to_streetnet(bool majormin, bool);
void merge_subgraph_to_streetnet_2(bool majormin, bool startorend);
void merge_subgraph_to_streetnet_3(bool majormin, bool startorend);
void merge_subgraph_to_streetnet_4(bool majormin, bool startorend);

void connect_outer_deadends_Reg();

void update_street_network();

bool find_closest_intersection(double p[2], double &dist, int &which_intersect,
							   icVector2 goDir, double cosang);

int trace_comp_nearby_intersect(double startp[2], int start_cell, icVector2 goDir, 
								int except_intersect, int except_edge, 
								int &which_intersect, double intersect_p[2],
								int &which_edge, int &which_samp, int MaxSearchCells);

bool search_nearby_intersection_around_oneCell(int cell, double threshold, double pt[2], int except_intersect,
											   icVector2 goDir, double cosang, int &which_intersect);

//void get_brush_approDir();