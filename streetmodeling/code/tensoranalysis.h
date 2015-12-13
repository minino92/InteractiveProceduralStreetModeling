/*
tensoranalysis.h
the head file of the tensor field analysis
*/

/*for degenerate point detections 09/25/2007*/
void locate_degeneratepts_tranvec(void);
void compute_tensor_at(int tri, double alpha[3], icMatrix2x2 &ten);

void compute_separatrixVec_degpt(int degpt_index);
int solve_ten_cubic(double a, double b, double d, double e, double solutions[3]);
int solve_ten_quadratic(double a, double b, double c, double solutions[2]);
int solve_ten_cubic_2(double a, double b, double c, double d, double solutions[4]);
int solve_ten_cubic_3(double a, double b, double c, double d, double solutions[4]);

/*detect degenerated points using quad-based mesh 09/25/2007*/
void locate_degpts_cells_tranvec_quad(void);
void compute_degpts_pos_tranvec_quad(void);
void compute_onedegpt_pos_tranvec_quad(int cellid, double &x, double &y);
void compute_a_alongx_degptlocate(double &a, icVector2 v0, icVector2 v1, icVector2 v2, icVector2 v3,
								  icVector2 &v0v1, icVector2 &v2v3);
void cal_all_eigenvecs_quad();
void cal_eigenvecs_onevert_quad(int vertid);
void normalized_tensorfield_quad();
void init_degpts();
void init_degpts_tri();

/*for tensor line tracing*/
void ten_trace(int face_id, double x, double y);
int trace_ten_in_one_tri(int &face_id, double globalp[2], int &flag);
bool compute_next_pt_tensor(double first[2], double second[2], int &face_id, double alpha[3]);
void get_next_tri(int &face_id, double pre[2], double cur[2], double param_t[2], 
					 int &PassVertornot, double alpha[3]);
void cross_AVertex(int &face_id, double cur_p[2], double pre_p[2], int &passornot);


/*computing tensor lines under quad-based mesh 09/25/2007*/

void compute_tensor_at_quad(int face, double x, double y, icMatrix2x2 &ten);
void get_x_y_ranges(int id, double &xstart, double &xend, 
					double &ystart, double &yend);
bool is_in_cell(int id, double x, double y);
bool is_in_reg_cell(int id, double x, double y);

double bilinear_interpolate(double a, double b, 
						  double f00, double f01, double f10, double f11);
void ten_trace_quad(int face_id, double x, double y, int type);
int trace_ten_in_one_tri_quad(int &face_id, double globalp[2], int &flag, int type);
bool compute_next_pt_tensor_quad_global(double first[2], double second[2], int &face_id);
//void derive_tensor_quad_global(const DP t, Vec_I_DP &y, Vec_O_DP &dydx);
bool cross_vertex_ten_quad(int &face_id, double cur_p[2], double pre_p[2], int &passornot, int type);
void get_cell_through_ver(int vertid, int &cell, int type);
void get_next_cell(int &face_id, double pre[2], double cur[2], 
					 int &PassVertornot, int type);
void get_next_cell_2(int &face_id, double pre[2], double cur[2], 
					 int &PassVertornot, int type);
bool get_nextpt_2ndeuler_ten_quad(double first[2], double second[2], int &face_id, int type);
bool get_nextpt_RK45_ten_quad(double first[2], double second[2], int &face_id, int type);

void global_test_trace(double x, double y);
bool get_nextpt_2ndeuler_ten_quad_gl(double first[2], double second[2]);
double bilinear_interpolate(double a, double b, 
						  double f00, double f01, double f10, double f11);

void RK23_2d(double pre_p[2], double next_p[2], double &hstep, double &hnext,
			 double eps, double &eps_did, void deriv(double cur_p[2], double vec[2]));
bool get_nextpt_RK23_ten_quad(double first[2], double second[2], int &face_id, int type);
