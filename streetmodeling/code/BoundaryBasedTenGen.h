/*
BoundaryBasedTenGen.h
*/

void cal_ten_with_fast_marching_quad(double disthred, int which_bound);
void cal_seeds_basedon_multibounds();
void cal_one_seed_from_one_sample(double start[2], double end[2]);
void place_tensorlines_consider_bounds();
void cal_ten_with_fast_marching_allbounds(double disthred);

void cal_one_boundSeed_inward(double pt[2], icVector2 line_dir,
							  double move_dist);

double cal_weight_for_a_seed(double sx, double sy, icVector2 norm, double disc_radius);
bool is_inland_pixel(double x, double y, double xstart, double xrang, double ystart, double yrang,
				 double dx, double dy,  unsigned char *map, int width);
