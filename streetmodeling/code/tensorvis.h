/*
tensorvis.h
the head file of the tensor field visualization
*/

void cal_eigenvector_vertices();

void cal_eigenvector_perver(int ver);

/*************************************************/
void major_vis();
void major_vis_quad();

void minor_vis();
void minor_vis_quad();

void mix_vis();
void mix_vis_quad();

void render_ibfv_tens(bool major_minor, bool x_y);
void render_ibfv_tens_quad(bool major_minor, bool x_y);


void render_alpha_map(bool major_minor);
void render_alpha_map_quad(bool major_minor);

void tensor_init_tex();
void init_texture_objs();

void make_tens_Patterns(void);

void render_majorfield();
void render_majorfield_quad();

void render_minorfield();
void render_minorfield_quad();

void vis_alpha_map();

void init_degpts();

void clear_all_tens();

void reset_inland();

void render_a_map(unsigned char *map);


//void display_quadmesh(GLenum mode);

/*the following code we generate the regular quadractic mesh 09/25/2007
NOTE: the result will be saved into the "Object" variable directly.
*/
typedef struct QuadVer{
	int id;      /*id = i*ydim+j*/
	double x, y;
}QuadVer;

//bool gen_regquad_mesh(int xdim, int ydim, double xstart, double xend, 
//				   double ystart, double yend);
//bool gen_regquad_vertices(int xdim, int ydim, double xstart, double xend, 
//				   double ystart, double yend);
//void finalize_quad_verts();
//void gen_regquad_faces();  
//void copy_to_vertexlist();
