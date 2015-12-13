

/*
SketchDesign.h
*/

void init_brushlist();
void add_to_current_brushPts(double x, double y);


/*Testing codes*/
void mark_vertices(double xstart, double xend, double ystart, double yend,
				   unsigned char *map);
void load_the_test_img(char *filename);
int get_region_id_for_cell(int cellid);
unsigned char get_region_id(double x, double y);
void convert_one_brush_to_atraj(int id);
void init_brush_trajs();

void copy_sketchlines_to_major_level1();
void convert_majRoads_to_sketches();

void init_sketchline_intersectionlists();
void init_sketchnet();
void cal_sketchlines_intersects();
void cal_boundary_intersects();
void init_domain_boundaries();
void release_domain_boundaries();

void search_for_sketchbased_graph();
void compute_sketchlineintersects_in_cell(int cellid);
void update_sketchcurveinfo();
void cal_blocks_sketchnet();


void get_linesegs_extendskecthcurve(int sketchid, double p1[2], int cell1, 
									double p2[2], int cell2, icVector2 dir, /*int lineid,*/
									/*int &nlines, */Trajectory *traj, int type);
void extend_one_sketchcurve(int id);
void extend_sketchcurves();
void init_sketchlineinfo_incells();


void init_sketchblocklist();

void get_linesegs_anytwopts(double p1[2], int cell1, 
						double p2[2], int cell2,
						Trajectory *traj, int type, int MaxIters);

/*  probably don't need it any more   */
void store_sketchline_to_majandmin();


void reset_region_indexes_all_designElems();
void gen_sketch_based_tensorfield(double widtheachbound);

void add_to_edgelist_one_cell(int cellid, int edgeid);

void save_sketch_brushes(char *filename);
bool load_sketch_brushes(char *filename);


