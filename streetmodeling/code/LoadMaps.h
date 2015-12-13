
/******************************************************/
/*                                                    */
/*                      LoadMaps.h                    */

void get_mask_map(double xstart, double xend, double ystart, double yend,
				  unsigned char *map, int width, int height);
void get_fitted_map(unsigned char *map, int width, int height,
					unsigned char *fittedmap, int fitwidth, int fitheight);
void get_density_value(double xstart, double xend, double ystart, double yend,
				  unsigned char *map, int width, int height);
