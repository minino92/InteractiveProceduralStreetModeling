/*
	generate the approximate boundary directions based on the 
	boundaries extracted from the water map 
*/
#include "stdafx.h"

#include "VFDataStructure.h"

#include "ImgBoundaryExtract.h"

#include "tensordesign.h"
#include "evenlystreamlines.h"
#include "HermiteCurve.h"
#include "RegionSmoothing.h"
#include "regionsmooth_quad.h"
#include "tensoranalysis.h"

#include "BoundaryBasedTenGen.h"


#include "shareinterfacevars.h"
extern SharedInterfaceVars sharedvars;

extern QuadMesh *quadmesh;

extern MapBoundaryList *mapboundarylist;

extern int *boundarycells;
extern int nboundarycells;
extern int *boundvertlist;
extern int nboundverts;
extern ctr_point *control_pts;        // allocate our control point array
extern int resolution;    // how many points in our output array
extern int num_shapecontrol_pts;
extern int num_curvepts_output;
extern ctr_point *out_pts;

extern TenRegularElem *ten_regularelems;
extern int nten_regelems;

extern MinHeap *narrowband ;

extern EvenStreamlinePlace *major;
extern EvenStreamlinePlace *minor;
extern bool brushinterfaceOn;

extern void CalOpenHermiteCurve();

double BoundRegionWidth = 5.;

extern int get_region_id_for_cell(int cellid);

extern unsigned char *fittedmap1;

extern double DistAwayBoundaries;


////////////////////////////////////////////////////////////////////////////////////

void get_imgbound_approDir()
{
	int i,j;
	MapBoundary *curboundary;

	for(i=0;i<mapboundarylist->nmapboundaries;i++)
	{
		curboundary=&mapboundarylist->mapboundarylist[i];
		for(j=0;j<curboundary->nelems-1;j++)
		{
			if(j%2==0)  /* use fewer elements along the boundaries */
			{
			set_ten_regBasis(curboundary->pts[j]->x, curboundary->pts[j]->y, 0);
			set_ten_regDir(curboundary->pts[j+1]->x, curboundary->pts[j+1]->y);
			}
		}
	}
}

void get_imgbound_approRegElems_onebound(int which_bound)
{
	int i;
	MapBoundary *curboundary=&mapboundarylist->mapboundarylist[which_bound];
	for(i=0;i<curboundary->nelems-1;i++)
	{
		set_ten_regBasis(curboundary->pts[i]->x, curboundary->pts[i]->y, 0);
		set_ten_regDir(curboundary->pts[i+1]->x, curboundary->pts[i+1]->y);
	}
}

/*
*/

void cal_init_multibound_verts_onebound(double disthred, int which_bound)
{
	int i, j, k;
	int curve_pos=0;
	QuadCell *face;
	QuadVertex *v;
	double start[2], end[2];
	double distanceLine, disttoline;
	int pre_pos = 0;

	num_shapecontrol_pts = 0;
	
	MapBoundary *curboundary=&mapboundarylist->mapboundarylist[which_bound];

	for(i=0;i<curboundary->nelems;i++)
	{
		/*add to the control point list*/
		add_to_shapeCtrPtsList(curboundary->pts[i]->x,curboundary->pts[i]->y, 
			get_cellID_givencoords(curboundary->pts[i]->x, curboundary->pts[i]->y));

	}
		
	resolution = get_Resolution_adp();
	//curboundary->i
	CalOpenHermiteCurve();

	/*calculate the cell strip that contains the curve*/
	//init_dis_verts();
	init_dis_cells();
	get_cellstrip_curve_quad(boundarycells, nboundarycells);

	/////////////////////
	//nboundverts=0;
	for(i=0; i<nboundarycells; i++)
	{
		face=quadmesh->quadcells[boundarycells[i]];

		for(j=0; j<face->nverts; j++)
			quadmesh->quad_verts[face->verts[j]]->visited=false;
	}

	for(i=0; i<nboundarycells; i++)
	{
		face=quadmesh->quadcells[boundarycells[i]];

		for(j=0; j<face->nverts; j++)
		{
			v=quadmesh->quad_verts[face->verts[j]];

			if(!v->visited)
			{
				boundvertlist[nboundverts]=v->index;
				nboundverts++;
			}
			v->visited=true;

			curve_pos = 0;

			if(!find_pos_curve(face->index, curve_pos, curve_pos))
			{
				curve_pos = pre_pos;
			}
			else
				pre_pos = curve_pos;

			if(curve_pos>0)
			{
				start[0]=out_pts[curve_pos-1].x;
				start[1]=out_pts[curve_pos-1].y;
				end[0]=out_pts[curve_pos].x;
				end[1]=out_pts[curve_pos].y;
				DistanceFromLine(v->x, v->y,
					start[0],start[1], end[0],end[1], distanceLine, disttoline);

				if(distanceLine<v->distance)
				{
					v->distance=distanceLine;
					//v->which_region=which_bound+1;   //11/17/2007
				}
			}

			if(curve_pos<num_curvepts_output-1)
			{
				start[0]=out_pts[curve_pos].x;
				start[1]=out_pts[curve_pos].y;
				end[0]=out_pts[curve_pos+1].x;
				end[1]=out_pts[curve_pos+1].y;
				DistanceFromLine(v->x, v->y,
					start[0],start[1], end[0],end[1], distanceLine, disttoline);
				
				if(distanceLine<v->distance)
				{
					v->distance=distanceLine;
					//v->which_region=which_bound+1;   //11/17/2007
				}
			}

		}
	}

	/*propagate to certain distance to get the vertices on the front*/
	//cal_ten_with_fast_marching_quad(disthred, which_bound);

}


void cal_ten_with_fast_marching_quad(double disthred, int which_bound)
{
	/*first, initialization*/
	/*since some vertices may contain the distance from other boundary*/
	/*11/17/2007*/
	//init_dis_verts();  
	//cal_init_bound_verts();

	reset_narrowband();

	init_narrow_band();

	///*second, loop as follows*/
	//	/*
	//	while (min-heap is not empty and no more far away points)
	//	{
	//		remove the root of the heap (with smallest distance value), 
	//		mark the corresponding vertex as "known";
	//		update the distance values of its neighbors if they are not "known",
	//		update the heap if they are "in the narrow band" with the new distance values;
	//		if any of them are not "in the narrow band" (i.e. "far away"), insert them to the heap and 
	//		mark them as "in the narrow band"
	//	}
	//	*/

	int nfinished = nboundverts;
	while(!narrowband->is_empty() && nfinished < quadmesh->nverts)
	{
		int cur_v = narrowband->FindSmallest();
		quadmesh->quad_verts[cur_v]->type = 2;

		/*we need to get the tensor on this vertex*/
		double t[4];

		if(!sharedvars.UseAllBoundsOn)
		{
			if(quadmesh->quad_verts[cur_v]->which_region>0)
			{
				/*we need to blend with previous tensor there*/
				get_tensor(quadmesh->quad_verts[cur_v]->x, quadmesh->quad_verts[cur_v]->y, t);
				quadmesh->quad_verts[cur_v]->Jacobian.entry[0][0]+=t[0];
				quadmesh->quad_verts[cur_v]->Jacobian.entry[0][1]+=t[1];
				quadmesh->quad_verts[cur_v]->Jacobian.entry[1][0]+=t[2];
				quadmesh->quad_verts[cur_v]->Jacobian.entry[1][1]+=t[3];
			}
			else
			{
				/*just compute the tensor there*/
				get_tensor(quadmesh->quad_verts[cur_v]->x, quadmesh->quad_verts[cur_v]->y, t);
				quadmesh->quad_verts[cur_v]->Jacobian.entry[0][0]=t[0];
				quadmesh->quad_verts[cur_v]->Jacobian.entry[0][1]=t[1];
				quadmesh->quad_verts[cur_v]->Jacobian.entry[1][0]=t[2];
				quadmesh->quad_verts[cur_v]->Jacobian.entry[1][1]=t[3];
			}
		}
		else
		{
				get_tensor(quadmesh->quad_verts[cur_v]->x, quadmesh->quad_verts[cur_v]->y, t);
				quadmesh->quad_verts[cur_v]->Jacobian.entry[0][0]=t[0];
				quadmesh->quad_verts[cur_v]->Jacobian.entry[0][1]=t[1];
				quadmesh->quad_verts[cur_v]->Jacobian.entry[1][0]=t[2];
				quadmesh->quad_verts[cur_v]->Jacobian.entry[1][1]=t[3];
		}

		//quadmesh->quad_verts[cur_v]->which_region=which_bound+1;  //11/17/2007

		boundvertlist[nboundverts]=cur_v;
		nboundverts++;
		
		/*update neighbor distances of cur_v*/
		if(quadmesh->quad_verts[cur_v]->distance<disthred)
			update_neighbor_Dis(cur_v);

		nfinished ++;
	}
}

void cal_ten_with_fast_marching_allbounds(double disthred)
{
	/*first, initialization*/
	/*since some vertices may contain the distance from other boundary*/
	/*11/17/2007*/
	reset_narrowband();

	init_narrow_band();

	///*second, loop as follows*/
	//	/*
	//	while (min-heap is not empty and no more far away points)
	//	{
	//		remove the root of the heap (with smallest distance value), 
	//		mark the corresponding vertex as "known";
	//		update the distance values of its neighbors if they are not "known",
	//		update the heap if they are "in the narrow band" with the new distance values;
	//		if any of them are not "in the narrow band" (i.e. "far away"), insert them to the heap and 
	//		mark them as "in the narrow band"
	//	}
	//	*/

	int nfinished = nboundverts;
	while(!narrowband->is_empty() && nfinished < quadmesh->nverts)
	{
		if(narrowband->elems==NULL)
		{
			int test=0;
			return;
		}

		int cur_v = narrowband->FindSmallest();
		quadmesh->quad_verts[cur_v]->type = 2;

		/*we need to get the tensor on this vertex*/
		double t[4];

		get_tensor(quadmesh->quad_verts[cur_v]->x, quadmesh->quad_verts[cur_v]->y, t);
		quadmesh->quad_verts[cur_v]->Jacobian.entry[0][0]=t[0];
		quadmesh->quad_verts[cur_v]->Jacobian.entry[0][1]=t[1];
		quadmesh->quad_verts[cur_v]->Jacobian.entry[1][0]=t[2];
		quadmesh->quad_verts[cur_v]->Jacobian.entry[1][1]=t[3];


		boundvertlist[nboundverts]=cur_v;
		nboundverts++;
		
		/*update neighbor distances of cur_v*/
		if(quadmesh->quad_verts[cur_v]->distance<disthred)
			update_neighbor_Dis(cur_v);

		nfinished ++;
	}

	/*  we set the brush affecting region as 0-indexed   */
	int i;
	for(i=0;i<nboundverts;i++)
	{
		quadmesh->quad_verts[boundvertlist[i]]->inbrushregion=true;
	}

	/*  update the region indexes of all the cells  */
	for(i=0;i<quadmesh->nfaces;i++)
		quadmesh->quadcells[i]->which_region=get_region_id_for_cell(i);
}



extern int *region_quadverts;                ////mesh vertices inside user selected region
extern int nregion_quadverts ;
extern int curMaxRegionQuadVerts;

/*
After setting the initial conditions of the multiboundary design,
we find out all the inner vertices
*/
void find_innerverts_multibounds()
{
	int i;

	nregion_quadverts=0;

	for(i=0; i<quadmesh->nverts; i++)
	{
		if(quadmesh->quad_verts[i]->which_region>0)
		{
			quadmesh->quad_verts[i]->OnBoundary=true;
			quadmesh->quad_verts[i]->InRegion=false;
			quadmesh->quad_verts[i]->RegionListID=-1;
		}
		else
		{
			quadmesh->quad_verts[i]->OnBoundary=false;
			quadmesh->quad_verts[i]->InRegion=true;
			quadmesh->quad_verts[i]->RegionListID=nregion_quadverts;

			region_quadverts[nregion_quadverts]=i;
			nregion_quadverts++;
		}
	}
}

void find_innerverts_multibounds_2()
{
	int i;

	nregion_quadverts=0;

	for(i=0; i<quadmesh->nverts; i++)
	{
		if(!quadmesh->quad_verts[i]->OnBoundary)
		{
			quadmesh->quad_verts[i]->InRegion=true;
			quadmesh->quad_verts[i]->RegionListID=nregion_quadverts;

			region_quadverts[nregion_quadverts]=i;
			nregion_quadverts++;
		}
	}
}


/*
Obtain the inner vertices under multiple boundary design
NOTE: we assume that the boundaries are stored in "mapboundarylist" variable
11/17/2007
*/
void obtain_smooth_region_multibounds(double widtheachbound)
{
	int i, j;

	/*initial all the vertices only once !!*/
	init_dis_verts();  
	init_quad_regionsmooth();
	//for(i=0;i<quadmesh->nverts;i++)
	//	quadmesh->quad_verts[i]->inbrushregion=false;

	/*set the regular elements using all the boundary information*/
	if(sharedvars.UseAllBoundsOn)
		get_imgbound_approDir();

	nboundverts = 0;
	nboundarycells=0;


	for(i=0;i<mapboundarylist->nmapboundaries;i++)
	{
		init_boundvertlist();

		cal_init_multibound_verts_onebound(widtheachbound, i);

		/*set the boundary orientation as the set of regular elements*/
		if(!sharedvars.UseAllBoundsOn)
			get_imgbound_approRegElems_onebound(i);

		/*  compute the tensor values on the vertices in current "narrow band"  */

		/*get the tensor values for each vertex inside the band with the user
		specified width*/
		//cal_ten_with_fast_marching_quad(widtheachbound, i);
	}

	int cur_v;
	double t[4]={0.};
	for(i=0;i<nboundverts;i++)
	{
		cur_v=boundvertlist[i];
		get_tensor(quadmesh->quad_verts[cur_v]->x, quadmesh->quad_verts[cur_v]->y, t);
		quadmesh->quad_verts[cur_v]->Jacobian.entry[0][0]=t[0];
		quadmesh->quad_verts[cur_v]->Jacobian.entry[0][1]=t[1];
		quadmesh->quad_verts[cur_v]->Jacobian.entry[1][0]=t[2];
		quadmesh->quad_verts[cur_v]->Jacobian.entry[1][1]=t[3];
	}

	/*   we currently disable the fast marching method 1/1/2008  */

	//cal_ten_with_fast_marching_allbounds(widtheachbound);

	/* use constrained optimization to set up the sparse linear system*
	   NOTE: we currently generate the tensor field everywhere inside the domain
	   1/1/2008
	*/

	for(i=0;i<quadmesh->nverts;i++)
	{
		if(quadmesh->quad_verts[i]->inland)
			quadmesh->quad_verts[i]->OnBoundary=false;
		else
			quadmesh->quad_verts[i]->OnBoundary=true;
	}

	/*  propagate one more neighborhood   */
	int k;
	for(k=0; k<3; k++)
	{
		int oldnboundverts=nboundverts;
		for(i=0;i<oldnboundverts;i++)
		{
			QuadVertex *v=quadmesh->quad_verts[boundvertlist[i]];
			QuadEdge *edge;
			for(j=0;j<v->Num_edge;j++)
			{
				edge=v->edges[j];
				cur_v=edge->verts[0];
				if(cur_v==boundvertlist[i])
					cur_v=edge->verts[1];
				if(!quadmesh->quad_verts[cur_v]->OnBoundary)
				{
					boundvertlist[nboundverts]=cur_v;
					quadmesh->quad_verts[cur_v]->OnBoundary=true;
					get_tensor(quadmesh->quad_verts[cur_v]->x, quadmesh->quad_verts[cur_v]->y, t);
					quadmesh->quad_verts[cur_v]->Jacobian.entry[0][0]=t[0];
					quadmesh->quad_verts[cur_v]->Jacobian.entry[0][1]=t[1];
					quadmesh->quad_verts[cur_v]->Jacobian.entry[1][0]=t[2];
					quadmesh->quad_verts[cur_v]->Jacobian.entry[1][1]=t[3];
					nboundverts++;
				}
			}
		}
	}

	for(i=0;i<nboundverts;i++)
	{
		quadmesh->quad_verts[boundvertlist[i]]->OnBoundary=true;
	}

	find_innerverts_multibounds_2();
	smooth_Jac_quadregion();

	cal_all_eigenvecs_quad();

	delete narrowband;
	narrowband=NULL;
}

/*
    We use traditional basis field approach to generate the tensor field from the 
	map boundaries  1/1/2008
*/

void obtain_field_basis()
{
	int i;
	double t[4]={0.};
	if(sharedvars.UseAllBoundsOn)
		get_imgbound_approDir();

	for(i=0;i<quadmesh->nverts;i++)
	{
		if(!quadmesh->quad_verts[i]->inland) 
			continue;

		get_tensor(quadmesh->quad_verts[i]->x, quadmesh->quad_verts[i]->y, t);
		quadmesh->quad_verts[i]->Jacobian.entry[0][0]=t[0];
		quadmesh->quad_verts[i]->Jacobian.entry[0][1]=t[1];
		quadmesh->quad_verts[i]->Jacobian.entry[1][0]=t[2];
		quadmesh->quad_verts[i]->Jacobian.entry[1][1]=t[3];
	}

	cal_all_eigenvecs_quad();
}


/*
The following routine finds the collection of seeds along the boundaries
NOTE: we make sure that the seeds are all "inland"
*/

SeedList *seedsalongbounds=NULL;

void cal_seeds_basedon_multibounds()
{
	if(seedsalongbounds != NULL)
		delete seedsalongbounds;

	seedsalongbounds=new SeedList(500);

	/*  compute the proper seeds based on the obtained boundary points  */
	int i,j;
	MapBoundary *curboundary;
	double start[2], end[2];
	icVector2 line_dir;

	for(i=0;i<mapboundarylist->nmapboundaries;i++)
	{
		curboundary=&mapboundarylist->mapboundarylist[i];
		for(j=0;j<curboundary->nelems-1;j++)
		{
			if(j%2!=0 && j!=curboundary->nelems-1) continue;

			//start[0]=curboundary->pts[j-1]->x;
			//start[1]=curboundary->pts[j-1]->y;
			start[0]=curboundary->pts[j]->x;
			start[1]=curboundary->pts[j]->y;
			end[0]=curboundary->pts[j+1]->x;
			end[1]=curboundary->pts[j+1]->y;

			line_dir.entry[0]=end[0]-start[0];
			line_dir.entry[1]=end[1]-start[1];

			//cal_one_seed_from_one_sample(start, end);

			/*  the moving distance should be a parameter adjusted by user!  */
			/*  NOTE: we don't always want the seeds move inward for the same distance!  */
			cal_one_boundSeed_inward(start, line_dir, DistAwayBoundaries*quadmesh->xinterval
				/*quadmesh->xinterval/0.9*/);
		}
	}

	/*   we need to randomly generate some seeds in land 12/31/2007  */
	for(i=0;i<mapboundarylist->nmapboundaries;i++)
	{
		curboundary=&mapboundarylist->mapboundarylist[i];
		for(j=0;j<curboundary->nelems-1;j++)
		{
			if(j%10!=0) continue;

			//start[0]=curboundary->pts[j-1]->x;
			//start[1]=curboundary->pts[j-1]->y;
			start[0]=curboundary->pts[j]->x;
			start[1]=curboundary->pts[j]->y;
			end[0]=curboundary->pts[j+1]->x;
			end[1]=curboundary->pts[j+1]->y;

			line_dir.entry[0]=end[0]-start[0];
			line_dir.entry[1]=end[1]-start[1];

			//cal_one_seed_from_one_sample(start, end);

			/*  the moving distance should be a parameter adjusted by user!  */
			cal_one_boundSeed_inward(start, line_dir, quadmesh->xinterval*10);
		}
	}

	for(i=0;i<mapboundarylist->nmapboundaries;i++)
	{
		curboundary=&mapboundarylist->mapboundarylist[i];
		for(j=0;j<curboundary->nelems-1;j++)
		{
			if(j%10!=0) continue;

			//start[0]=curboundary->pts[j-1]->x;
			//start[1]=curboundary->pts[j-1]->y;
			start[0]=curboundary->pts[j]->x;
			start[1]=curboundary->pts[j]->y;
			end[0]=curboundary->pts[j+1]->x;
			end[1]=curboundary->pts[j+1]->y;

			line_dir.entry[0]=end[0]-start[0];
			line_dir.entry[1]=end[1]-start[1];

			//cal_one_seed_from_one_sample(start, end);

			/*  the moving distance should be a parameter adjusted by user!  */
			cal_one_boundSeed_inward(start, line_dir, quadmesh->xinterval*20);
		}
	}

	for(i=0;i<mapboundarylist->nmapboundaries;i++)
	{
		curboundary=&mapboundarylist->mapboundarylist[i];
		for(j=0;j<curboundary->nelems-1;j++)
		{
			if(j%20!=0) continue;

			//start[0]=curboundary->pts[j-1]->x;
			//start[1]=curboundary->pts[j-1]->y;
			start[0]=curboundary->pts[j]->x;
			start[1]=curboundary->pts[j]->y;
			end[0]=curboundary->pts[j+1]->x;
			end[1]=curboundary->pts[j+1]->y;

			line_dir.entry[0]=end[0]-start[0];
			line_dir.entry[1]=end[1]-start[1];

			//cal_one_seed_from_one_sample(start, end);

			/*  the moving distance should be a parameter adjusted by user!  */
			cal_one_boundSeed_inward(start, line_dir, quadmesh->xinterval*40);
		}
	}

	/*   try to find the seed in the middle of the land  */
	for(i=0;i<mapboundarylist->nmapboundaries;i++)
	{
		curboundary=&mapboundarylist->mapboundarylist[i];
		for(j=0;j<curboundary->nelems-1;j++)
		{
			if(j%quadmesh->XDIM!=0) continue;

			//start[0]=curboundary->pts[j-1]->x;
			//start[1]=curboundary->pts[j-1]->y;
			start[0]=curboundary->pts[j]->x;
			start[1]=curboundary->pts[j]->y;
			end[0]=curboundary->pts[j+1]->x;
			end[1]=curboundary->pts[j+1]->y;

			line_dir.entry[0]=end[0]-start[0];
			line_dir.entry[1]=end[1]-start[1];

			//cal_one_seed_from_one_sample(start, end);

			/*  the moving distance should be a parameter adjusted by user!  */
			cal_one_boundSeed_inward(start, line_dir, quadmesh->xinterval*quadmesh->XDIM/2.);
		}
	}

	/* update the maximal weight and minmal weight  */
	seedsalongbounds->update_max_min_weights();
}

void cal_one_seed_from_one_sample(double start[2], double end[2])
{
	/*get the vector*/
	icVector2 line_dir, trace_dir;
	line_dir.entry[0]=end[0]-start[0];
	line_dir.entry[1]=end[1]-start[1];
	normalize(line_dir);

	trace_dir.entry[0]=-line_dir.entry[1];
	trace_dir.entry[1]=line_dir.entry[0];

	/*move start point along trace_dir with certain distance*/
	double newpt[2];
	newpt[0]=start[0]+quadmesh->xinterval*trace_dir.entry[0];
	newpt[1]=start[1]+quadmesh->xinterval*trace_dir.entry[1];
	int cellid=get_cellID_givencoords(newpt[0], newpt[1]);

	if(!is_not_inland(cellid))
	{
		/*add to the sample list*/
		Seed *newseed=(Seed*)malloc(sizeof(Seed));
		newseed->pos[0]=newpt[0];
		newseed->pos[1]=newpt[1];
		newseed->triangle=cellid;
		seedsalongbounds->append(newseed);
		return;
	}

	newpt[0]=start[0]-quadmesh->xinterval*trace_dir.entry[0];
	newpt[1]=start[1]-quadmesh->xinterval*trace_dir.entry[1];
	cellid=get_cellID_givencoords(newpt[0], newpt[1]);
	if(!is_not_inland(cellid))
	{
		/*add to the sample list*/
		Seed *newseed=(Seed*)malloc(sizeof(Seed));
		newseed->pos[0]=newpt[0];
		newseed->pos[1]=newpt[1];
		newseed->triangle=cellid;
		seedsalongbounds->append(newseed);
	}
}


/*
    We obtain one seed along the boundary, but move slightly away the water
	12/29/2007
	NOTE: this is not accurate to use the flag of a grid cell to judge
	whether a point is "inland" or not. But if we use pixel level, it will
	be slow as well.
*/

void cal_one_boundSeed_inward(double pt[2], icVector2 line_dir,
							  double move_dist)
{
	double xstart=quadmesh->xstart;
	double ystart=quadmesh->ystart;
	double xrang=quadmesh->xend-quadmesh->xstart;
	double yrang=quadmesh->yend-quadmesh->ystart;
	double dx=xrang/511;
	double dy=yrang/511;


	icVector2 norm;
	norm.entry[0]=-line_dir.entry[1];
	norm.entry[1]=line_dir.entry[0];

	normalize(norm);

	double newpt[2];
	newpt[0]=pt[0]+move_dist*norm.entry[0];
	newpt[1]=pt[1]+move_dist*norm.entry[1];
	int cellid=get_cellID_givencoords(newpt[0], newpt[1]);


	if(cellid>=0&&cellid<=quadmesh->nfaces&&/*is_inland_pixel(newpt[0],newpt[1], xstart, xrang, ystart, yrang, dx, dy,
			fittedmap1, 512)*/!is_not_inland(cellid))
	{
		/*add to the sample list*/
		Seed *newseed=(Seed*)malloc(sizeof(Seed));
		newseed->pos[0]=newpt[0];
		newseed->pos[1]=newpt[1];
		newseed->triangle=cellid;

		newseed->weight=cal_weight_for_a_seed(newseed->pos[0], newseed->pos[1],
			norm, move_dist);
		seedsalongbounds->sorted_add(newseed);
		//newseed->weight=0;
		//seedsalongbounds->append(newseed);
	}
	else
	{
		newpt[0]=pt[0]+0.2*move_dist*norm.entry[0];
		newpt[1]=pt[1]+0.2*move_dist*norm.entry[1];
		cellid=get_cellID_givencoords(newpt[0], newpt[1]);
		if(cellid>=0&&cellid<=quadmesh->nfaces&&/*is_inland_pixel(newpt[0],newpt[1], xstart, xrang, ystart, yrang, dx, dy,
			fittedmap1, 512)*/!is_not_inland(cellid))
		{
			/*we try a smaller amount of distance*/
			Seed *newseed=(Seed*)malloc(sizeof(Seed));
			newseed->pos[0]=newpt[0];
			newseed->pos[1]=newpt[1];
			newseed->triangle=cellid;

			newseed->weight=cal_weight_for_a_seed(newseed->pos[0], newseed->pos[1],
				norm, move_dist);
			seedsalongbounds->sorted_add(newseed);
		}

		else
		{
			/*  we consider the opposite diretion  */

			newpt[0]=pt[0]-move_dist*norm.entry[0];
			newpt[1]=pt[1]-move_dist*norm.entry[1];
			cellid=get_cellID_givencoords(newpt[0], newpt[1]);
			if(cellid>=0&&cellid<=quadmesh->nfaces&&/*is_inland_pixel(newpt[0],newpt[1], xstart, xrang, ystart, yrang, dx, dy,
			fittedmap1, 512)*/!is_not_inland(cellid))
			{
				/*add to the sample list*/
				Seed *newseed=(Seed*)malloc(sizeof(Seed));
				newseed->pos[0]=newpt[0];
				newseed->pos[1]=newpt[1];
				newseed->triangle=cellid;

				newseed->weight=cal_weight_for_a_seed(newseed->pos[0], newseed->pos[1],
					norm, -move_dist);
				seedsalongbounds->sorted_add(newseed);
			}

			else
			{
				/*  we try a smaller amount of the distance  */
				newpt[0]=pt[0]-0.2*move_dist*norm.entry[0];
				newpt[1]=pt[1]-0.2*move_dist*norm.entry[1];
				cellid=get_cellID_givencoords(newpt[0], newpt[1]);
				if(cellid>=0&&cellid<=quadmesh->nfaces&&!is_not_inland(cellid))
				{
					/*we try a smaller amount of distance*/
					Seed *newseed=(Seed*)malloc(sizeof(Seed));
					newseed->pos[0]=newpt[0];
					newseed->pos[1]=newpt[1];
					newseed->triangle=cellid;

					newseed->weight=cal_weight_for_a_seed(newseed->pos[0], newseed->pos[1],
						norm, move_dist);
					seedsalongbounds->sorted_add(newseed);
				}
			}
		}
	}
}

/*
    judge whether a point is "inland or not" in pixel level
*/

bool is_inland_pixel(double x, double y, double xstart, double xrang, double ystart, double yrang,
				 double dx, double dy,  unsigned char *map, int width)
{
	int c=(x-xstart)/dx;
	int r=(y-ystart)/dy;

	int id=(r*(width)+c);

	if(r>=width || c>=width || r<0 || c<0) 
		return false;

	if(map[3*id]<125)
		return true;
	else
		return false;
}

/*
    Define a weighting scheme based on the position of the seeds
	NOTE: for simplicity, we are using discrete weights only now
*/

double cal_weight_for_a_seed(double sx, double sy, icVector2 norm, double disc_radius)
{
	/**/

	double weight=0;
	double tp[2];
	double low, high, curvalue;
	//double mov_dist;

	low=0; high=2;
	double xstart=quadmesh->xstart;
	double ystart=quadmesh->ystart;
	double xrang=quadmesh->xend-quadmesh->xstart;
	double yrang=quadmesh->yend-quadmesh->ystart;
	double dx=xrang/511;
	double dy=yrang/511;

	tp[0]=sx+2.*disc_radius*norm.entry[0];
	tp[1]=sy+2.*disc_radius*norm.entry[1];
	
	if(is_inland_pixel(tp[0],tp[1], xstart, xrang, ystart, yrang, dx, dy,
			fittedmap1, 512))
			return (1+.25);

	int level=0;

	while(level<6)
	{
		curvalue=(low+high)/2.;
		//mov_dist=pow(2.f, curvalue);

		tp[0]=sx+curvalue*disc_radius*norm.entry[0];
		tp[1]=sy+curvalue*disc_radius*norm.entry[1];

		/*   use the loaded map to judge whether a point is "inland" or not  */
		if(is_inland_pixel(tp[0],tp[1], xstart, xrang, ystart, yrang, dx, dy,
			fittedmap1, 512))
		{
			low=curvalue;
		}
		else
		{
			high=curvalue;
		}
		
		level++;
	}

	return (1+2./pow(2.,curvalue));

}

void place_tensorlines_consider_bounds()
{
	if(seedsalongbounds!=NULL)
	{
		major->reset_placement_quad();
		major->init_major_minor_line_info(false);

		major->place_streamlines(0, brushinterfaceOn, seedsalongbounds);
		
		minor->reset_placement_quad();
		minor->init_major_minor_line_info(true);
		minor->place_streamlines(1, brushinterfaceOn, seedsalongbounds);
	}
	else
		place();
}