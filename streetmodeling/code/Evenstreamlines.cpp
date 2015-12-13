
/*
This file contains the routines for even streamline placement.

Created and modified by Guoning Chen
copyright @2007
*/

#include "stdafx.h"
#include "vfdatastructure.h"

#include "EvenlyStreamlines.h"

#include "tensoranalysis.h"
#include "tensordesign.h"
#include "tensoranalysis.h"

#include "computeroadvis.h"

#include "regionsmooth_quad.h"

#include "shareinterfacevars.h"
extern SharedInterfaceVars sharedvars;

double road_width[]={
	0.0026,
	0.0035,
	0.0042,
	0.005,
	0.005
};

extern QuadMesh *quadmesh;
extern DegeneratePt *degpts;
extern int globalface;
extern icVector2 tenline_dir_global, tenline_dir_global_p;
extern double euler_stepsize;
extern double predict_stepsize;/* = quadmesh->xinterval;*/

extern int *Extend_link(int *edge_link, int Num_edges);
extern bool IsRepeated(int *a, int b, int num);
extern void cal_eigen_vector_sym(icMatrix2x2 m, icVector2 ev[2]);

extern bool StoreToGlobalList(CurvePoints *temp, int num);
extern int cal_intersect(double PointA[2], double PointB[2], double PointC[2], 
						 double PointD[2], double t[2]);

extern void sample_along_tensorline_from_to(Trajectory *traj, int start_lineseg, double startp[2],
									 int end_lineseg, double endp[2], 
									 int edgeid, StreetNet *net);



EvenStreamlinePlace *major = NULL;
EvenStreamlinePlace *minor = NULL;

/*  for dead end computation 11/24/2007  */
TrajectoryList *t_major = NULL;
TrajectoryList *t_minor = NULL;

/*  Major roads (or Highways) */
int which_level=1;
bool majorroadsexisted=false;
EvenStreamlinePlace *major_level1 = NULL;
EvenStreamlinePlace *minor_level1 = NULL;

StreetNet *majRoadnet=NULL;
TensorLineIntersectionInfoList **majRoadintersectinfo_maj=NULL;
TensorLineIntersectionInfoList **majRoadintersectinfo_min=NULL;
int pre_nmajRoads_maj, pre_nmajRoads_min;

int prev_nmajors = 0;
int prev_nminors = 0;

bool brushinterfaceOn=false;


//double majorDensity = 1.5;
//double minorDensity = 1.5;
//double mintenline_length=1.5;

double majorDensity = 13.5;
double minorDensity = 11.5;
double mintenline_length=5.9;

extern unsigned char *fittedmap1;
extern bool is_inland_pixel(double x, double y, double xstart, double xrang, double ystart, double yrang,
				 double dx, double dy,  unsigned char *map, int width);

//double min_waterwidth;
//double dist_follow_bound;

StreetNet *streetnet = NULL;
TensorLineIntersectionInfoList **majorintersectinfo = NULL;
TensorLineIntersectionInfoList **minorintersectinfo = NULL;

extern SeedList *seedsalongbounds;
extern double hstep;

/*  for the tracing crossing water region  */
extern double MinDistCrossRiver, MinAngCrossRiver, MaxDistFollowBoundary;

extern Degenerate_Design *ten_designelems ;


/*     This variable is used to show whether highway network 
		has been specified
*/
#include "sketchdesign.h"

extern bool highwayexisted;
extern TrajectoryList *sketchlines;

extern bool flag_loadmap;

//extern int temp_count;


/*  implement the sorting of the intersections locating at the same line segment 1/12/2008 
    We only consider this for minor road placement
*/
int TensorLineIntersectionInfoList::cal_pos_on_same_lineseg(IntersectionInfo *newinfo, 
															int start, int end)
{
	if(major==NULL || minor==NULL || trajID<0 || streetnet==NULL) 
		return -1;  /* bug */

	Trajectory *curtraj;
	if(!majormin)
	{
		curtraj=major->evenstreamlines->trajs[trajID];
	}
	else
	{
		curtraj=minor->evenstreamlines->trajs[trajID];
	}

	Intersection *curintersect;
	int i;
	double newDist, curDist;
	icVector2 DistVec;

	/* compute the distance from the new intersection to the gstart point of the line segment */
	LineSeg *theline=&curtraj->linesegs[newinfo->lineseg_id];
	curintersect=streetnet->nodelist->intersects[newinfo->intersect_id];
	DistVec.entry[0]=curintersect->gpos[0]-theline->gstart[0];
	DistVec.entry[1]=curintersect->gpos[1]-theline->gstart[1];
	newDist=length(DistVec);

	for(i=start; i<=end; i++)
	{
		curintersect=streetnet->nodelist->intersects[infolist[i]->intersect_id];
		DistVec.entry[0]=curintersect->gpos[0]-theline->gstart[0];
		DistVec.entry[1]=curintersect->gpos[1]-theline->gstart[1];
		curDist=length(DistVec);

		if(newDist<curDist)
		{
			return i;
		}
	}

	return (end+1);
}


int *extend_link(int *edge_link, int Num_edges)
{
    int *temp = edge_link;
	edge_link = (int *) malloc(sizeof(int)*(Num_edges + 1));
	if(Num_edges > 0)
	{
		for(int i = 0; i < Num_edges; i++)
			edge_link[i] = temp[i];
		//delete temp;
		free (temp);
	}

	return edge_link;
}
void init_evenplace_ten()
{
	FILE *fp;

	if(major != NULL)
		delete major;

	major = new EvenStreamlinePlace(false);
	major->init();

	if(minor != NULL)
		delete minor;

	minor = new EvenStreamlinePlace(true);
	minor->init();

	if(which_level==1)
	{
		//if(major_level1 != NULL)
		//	delete major_level1;
		//major_level1=new EvenStreamlinePlace(false);
		//major_level1->init();

		//if(minor_level1 != NULL)
		//	delete minor_level1;
		//minor_level1=new EvenStreamlinePlace(true);
		//minor_level1->init();
		init_level1_placement();
	}

}

void place()
{
	FILE *fp;

	major->reset_placement_quad();
	major->init_major_minor_line_info(false);
	major->place_streamlines(0, brushinterfaceOn);
	
	minor->reset_placement_quad();
	minor->init_major_minor_line_info(true);
	minor->place_streamlines(1, brushinterfaceOn);
}




/*********************************************************************************/

/*                   For major road placement                     */
void init_level1_placement()
{
	if(major_level1 != NULL)
		delete major_level1;
	major_level1=new EvenStreamlinePlace(false);
	major_level1->init();
	//major_level1->set_default_parameters(false);
	//major_level1->reset_placement_quad();
	//major_level1->init_major_minor_line_info(false);

	if(minor_level1 != NULL)
		delete minor_level1;
	minor_level1=new EvenStreamlinePlace(true); /* could change many things!!! 1/21/2008 */
	minor_level1->init();
	//minor_level1->set_default_parameters(true);
	//minor_level1->reset_placement_quad();
	//minor_level1->init_major_minor_line_info(true);
}

void place_alternative_level1_makeup_seeds()
{
	SeedList *inputseeds=new SeedList(1);

	Seed *aseed=(Seed*)malloc(sizeof(Seed));

	QuadCell *face=quadmesh->quadcells[(int)(quadmesh->nfaces/2)];
	aseed->pos[0]=face->x_start_coord+quadmesh->xinterval/2;
	aseed->pos[1]=face->y_start_coord+quadmesh->yinterval/2;
	aseed->triangle=(int)(quadmesh->nfaces/2);
	inputseeds->append(aseed);

	if(!sharedvars.JobardMethodOn)
		place_alternative_level1(inputseeds);
	else
	{
		major_level1->place_streamlines(0, false, inputseeds);
		minor_level1->place_streamlines(1, false, inputseeds);
	}

	/*temporary solution here*/
	//major = major_level1;
	//minor = minor_level1;

	/*   we finally need to copy the obtained major roads to the regular tracing
	*/

	delete inputseeds;

	/*   we compute the major road network here   */
	init_majRoadnet();
	compute_majRoad_intersects();
	search_for_connection_majRoad();
}


/*
     For placing level 1 tensor lines under the loading map mode
*/
void place_alternative_level1_loadmap()
{
	SeedList *inputseeds=new SeedList();

	/*   copy seeds from the boundary    */
	inputseeds->copy_otherseedList(seedsalongbounds);

	//MinDistCrossRiver=(MinDistCrossRiver/map_xrang*(quadmesh->xend-quadmesh->xstart))/
	//	quadmesh->xinterval;
	//MinAngCrossRiver=MinAngCrossRiver/180*M_PI;
	//MaxDistFollowBoundary=(MaxDistFollowBoundary/map_xrang*(quadmesh->xend-quadmesh->xstart))/
	//	quadmesh->xinterval;

	if(!sharedvars.JobardMethodOn)
		place_alternative_level1(inputseeds);
	else
	{
		major_level1->place_streamlines(0, false, inputseeds);
		minor_level1->place_streamlines(1, false, inputseeds);
	}


	/*   we finally need to copy the obtained major roads to the regular tracing
	*/

	delete inputseeds;

	/*   we compute the major road network here   */
	init_majRoadnet();
	compute_majRoad_intersects();
	search_for_connection_majRoad();

	/*  we now extend the dead end if user desires  */

}

/*
In this routine, we perform alternative tracing between major and minor fields
*/
void place_alternative_level1(SeedList *inputseeds)
{
	/*initialization*/
	if(major_level1 != NULL)
		delete major_level1;
	major_level1=new EvenStreamlinePlace(false);
	major_level1->init();
	major_level1->set_default_parameters(false);
	major_level1->reset_placement_quad();
	major_level1->init_major_minor_line_info(false);

	if(minor_level1 != NULL)
		delete minor_level1;
	minor_level1=new EvenStreamlinePlace(true);
	minor_level1->init();
	minor_level1->set_default_parameters(true);
	minor_level1->reset_placement_quad();
	minor_level1->init_major_minor_line_info(true);

	SeedList *majorseeds=new SeedList(100);
	SeedList *minorseeds=new SeedList(100);

	majorseeds->copy_otherseedList(inputseeds);
	minorseeds->copy_otherseedList(inputseeds);

	/*The following is the alternative process to generate
	the level 1 tensor line placement*/
	int majormin_whoseturn = 0;  /*initially, major field*/
	int test_counter=0;
	while(1)
	{
		if(majorseeds->nseeds==0 && minorseeds->nseeds==0)
			break;

		SeedList *outputseeds=new SeedList(100);
		if(majormin_whoseturn==0)
		{
			if(place_onetensorline_level1(false, majorseeds, outputseeds))
			{
				/*copy the obtain new seeds*/
				minorseeds->copyandappend_otherseedList(outputseeds);
			}

			majormin_whoseturn = (majormin_whoseturn+1)%2;
		}
		else
		{
			if(place_onetensorline_level1(true, minorseeds, outputseeds))
			{
				/*copy the obtain new seeds*/
				majorseeds->copyandappend_otherseedList(outputseeds);
			}
			majormin_whoseturn = (majormin_whoseturn+1)%2;
		}

		delete outputseeds;
		test_counter++;

	}

	/*  NOTE: we also need to delete the sample point list associated with each major road  
	    1/2/2008
	*/

	/*release the sample point list*/
	int i;
	for(i = 0; i < major_level1->evenstreamlines->curMaxNumTrajs; i++)
	{
		delete major_level1->samplepts[i];
		major_level1->samplepts[i] = NULL;
	}
	delete [] major_level1->samplepts;
	major_level1->samplepts = NULL;

	for(i = 0; i < minor_level1->evenstreamlines->curMaxNumTrajs; i++)
	{
		delete minor_level1->samplepts[i];
		minor_level1->samplepts[i] = NULL;
	}
	delete [] minor_level1->samplepts;
	minor_level1->samplepts = NULL;

	minor_level1->reset_placement_quad();
}

/*
This routine will place one and only one tensor line according to the
specified type and the input seed list
*/
bool place_onetensorline_level1(bool type, SeedList *inputseeds, SeedList *outputseeds)
{
	EvenStreamlinePlace *curplace;
	double startpt[2];
	int startcell;
	int fieldtype;
	if(!type) /*major direction*/
	{
		curplace=major_level1;
		fieldtype=0;
	}
	else /*minor direction*/
	{
		curplace=minor_level1;
		fieldtype=1;
	}

	if(inputseeds == NULL)
	{
		/*obtain the seeds from previous computed tensor lines*/
		if(curplace->evenstreamlines->ntrajs==0)
		{
			/*need to randomly generate a seed to start tracing*/
			startcell = (int)((double)(rand())/RAND_MAX *(quadmesh->nfaces-1));
			QuadCell *face=quadmesh->quadcells[startcell];
			startpt[0] = face->x_start_coord+quadmesh->xinterval/2.;
			startpt[1] = face->y_start_coord+quadmesh->yinterval/2.;
		}
		else
		{
			if(curplace->seedpts==NULL)
			{
				if(curplace->cur_traj==curplace->evenstreamlines->ntrajs-1)
				{
					/*currently, we don't generate new seeds for this case*/
					return false;
				}
				else
				{
					/*compute new seeds*/
					curplace->seedpts=new SeedList();
					curplace->cal_seeds(curplace->cur_traj, curplace->dsep,
						curplace->every_nsample, fieldtype, brushinterfaceOn);
					/*choose a proper seed*/
					if(!get_a_proper_seed(curplace, startpt, startcell, fieldtype))
						return false;
				}
			}
			else
			{
				if(!get_a_proper_seed(curplace, startpt, startcell, fieldtype))
					return false;
			}
		}
	}

	else  /*trace use one of the input seeds*/
	{
		if(curplace->close_to_cur_streamline(inputseeds->seeds[0]->pos,
			inputseeds->seeds[0]->triangle, &curplace->evenstreamlines->ntrajs, 0, 
			curplace->minstartdist, 1., 0))
		{
			inputseeds->del_Node_byindex(0);
			return false;
		}
		else
		{
			startcell=inputseeds->seeds[0]->triangle;
			startpt[0]=inputseeds->seeds[0]->pos[0];
			startpt[1]=inputseeds->seeds[0]->pos[1];
		}
	}

	bool trace_flag;

	//if(!flag_loadmap)
	//{
	//	trace_flag = curplace->grow_a_tensorline(startpt, startcell, curplace->percentage_dsep*curplace->dsep, 
	//	curplace->discsize, curplace->sample_interval, curplace->loopdsep, curplace->dist2sing, 
	//	curplace->streamlinelength, fieldtype, brushinterfaceOn);
	//}

	//else
	//{
		trace_flag = curplace->grow_a_majRoad(startpt, startcell, curplace->percentage_dsep*curplace->dsep, 
				curplace->discsize, curplace->sample_interval, curplace->loopdsep, curplace->dist2sing, 
				curplace->streamlinelength, fieldtype, brushinterfaceOn);
	//}

	if(trace_flag)
	{
		/*generate a sequence of seeds along this tensor line as the return*/

		int i;
		SamplePtList *cursamplist=curplace->samplepts[curplace->evenstreamlines->ntrajs-1];

		outputseeds->nseeds=0;
		for(i=0;i<cursamplist->nsamples;i++)
		{
			if(i!=0 && i%curplace->every_nsample==0 || i==cursamplist->nsamples-1)
			{
				Seed *newseed=(Seed*)malloc(sizeof(Seed));
				newseed->pos[0]=cursamplist->samples[i]->gpt[0];
				newseed->pos[1]=cursamplist->samples[i]->gpt[1];
				newseed->triangle=cursamplist->samples[i]->triangle;
				outputseeds->append(newseed);
			}
		}

		/*we may need to extend the trajectory list*/
		if(curplace->evenstreamlines->isFull())
			if(!curplace->extend_trajList())
				exit(-1);
		/*remove the used seed*/
		inputseeds->del_Node_byindex(0);
		
		return true;
	}
	else
	{
		/*remove the used seed*/
		inputseeds->del_Node_byindex(0);

		return false;
	}
}


/* 
   After computing the major roads, we try to construct a major road network
*/
void init_majRoad_intersectionlists()
{
	int i;
	if(majRoadintersectinfo_maj != NULL)
	{
		for(i=0;i<pre_nmajRoads_maj; i++)
			if(majRoadintersectinfo_maj[i]!=NULL)
			{
				delete majRoadintersectinfo_maj[i];
				majRoadintersectinfo_maj[i]=NULL;
			}
		delete [] majRoadintersectinfo_maj;
		majRoadintersectinfo_maj = NULL;
	}

	if(majRoadintersectinfo_maj == NULL)
	{
		majRoadintersectinfo_maj = new TensorLineIntersectionInfoList *[major_level1->evenstreamlines->ntrajs];
		for(i=0; i<major_level1->evenstreamlines->ntrajs; i++)
			majRoadintersectinfo_maj[i] = new TensorLineIntersectionInfoList();
		pre_nmajRoads_maj=major_level1->evenstreamlines->ntrajs;
	}
	
	if(majRoadintersectinfo_min != NULL)
	{
		for(i=0;i<pre_nmajRoads_min; i++)
			if(majRoadintersectinfo_min[i]!=NULL)
			{
				delete majRoadintersectinfo_min[i];
				majRoadintersectinfo_min[i] = NULL;
			}
		delete [] majRoadintersectinfo_min;
		majRoadintersectinfo_min = NULL;
	}
	
	if(majRoadintersectinfo_min == NULL)
	{
		majRoadintersectinfo_min = new TensorLineIntersectionInfoList *[minor_level1->evenstreamlines->ntrajs];
		for(i=0; i<minor_level1->evenstreamlines->ntrajs; i++)
			majRoadintersectinfo_min[i] = new TensorLineIntersectionInfoList();
		pre_nmajRoads_min=minor_level1->evenstreamlines->ntrajs;
	}
}


void init_majRoadnet()
{
	if(majRoadnet == NULL)
		//delete streetnet;
		majRoadnet=new StreetNet(20, 80);
	majRoadnet->reset_streetnet();
	init_majRoad_intersectionlists();
}




bool compute_majRoad_intersect_between_twolines(int majtraj, int majstart, int majend,
										int mintraj, int minstart, int minend,
										double intersect[2], 
										int &majlinesegid, int &minlinesegid)
{
	int i, j;
	Trajectory *major1, *minor1;
	double A[2], B[2], C[2], D[2], t[2];
	major1=major_level1->evenstreamlines->trajs[majtraj];
	minor1=minor_level1->evenstreamlines->trajs[mintraj];
	
	//for(i=majstart; i<=majend; i++)
	for(i=max(0, majstart-1); i<=min(major1->nlinesegs-1, majend+1); i++)
	{
		A[0]=major1->linesegs[i].gstart[0];
		A[1]=major1->linesegs[i].gstart[1];
		B[0]=major1->linesegs[i].gend[0];
		B[1]=major1->linesegs[i].gend[1];

		//for(j=minstart; j<=minend; j++)
		for(j=max(0, minstart-1); j<=min(minor1->nlinesegs-1, minend+1); j++)
		{
			C[0]=minor1->linesegs[j].gstart[0];
			C[1]=minor1->linesegs[j].gstart[1];
			D[0]=minor1->linesegs[j].gend[0];
			D[1]=minor1->linesegs[j].gend[1];

			//if(GetIntersection2(A, B, C, D, t)==1)
			if(cal_intersect(A, B, C, D, t)==1)
			{
				intersect[0]=A[0]+t[0]*(B[0]-A[0]);
				intersect[1]=A[1]+t[0]*(B[1]-A[1]);
				majlinesegid = i;
				minlinesegid = j;
				return true;
			}
		}
	}
	return false;
}

void compute_majRoad_intersects_in_cell(int cellid)
{
	QuadCell *face = quadmesh->quadcells[cellid];

	/*method 1: we use brute force method 10/02/2007*/
	int i, j;
	int majlinesegid, minlinesegid;
	double intersect[2] = {0.};

	/*Test file*/
	FILE *fp;

		for(i=0; i<face->majorlines->nlines; i++)
		{
			for(j=0; j<face->minorlines->nlines; j++)
			{
				if(compute_majRoad_intersect_between_twolines(face->majorlines->lines[i]->whichtraj,
					face->majorlines->lines[i]->start, face->majorlines->lines[i]->end,
					face->minorlines->lines[j]->whichtraj,
					face->minorlines->lines[j]->start, face->minorlines->lines[j]->end,
					intersect, majlinesegid, minlinesegid))
				{
					/*create a new intersection data*/
					Intersection *newintersect=(Intersection*)malloc(sizeof(Intersection));
					newintersect->gpos[0]=intersect[0];
					newintersect->gpos[1]=intersect[1];
					cellid=get_cellID_givencoords(intersect[0], intersect[1]);
					newintersect->cellid=cellid;
					newintersect->majorline_id=face->majorlines->lines[i]->whichtraj;
					newintersect->minorline_id=face->minorlines->lines[j]->whichtraj;
					newintersect->majlineseg=majlinesegid;
					newintersect->minlineseg=minlinesegid;
					newintersect->nadjedges=0;
					newintersect->adj_edges=NULL;
					newintersect->endpt=false;
					newintersect->inside_region=false;
					newintersect->deleted=false;

					majRoadnet->nodelist->addNew(newintersect);

					/*   We need to save the intersection into the list of the intersections
					     in the triangle that contains this intersection.
					*/
					add_to_cell_intersectlist(cellid, newintersect->index, false);
					
					IntersectionInfo *majinfo=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
					majinfo->intersect_id=majRoadnet->nodelist->nelems-1;
					majinfo->lineseg_id = majlinesegid;
					majRoadintersectinfo_maj[newintersect->majorline_id]->sorted_add(majinfo);

					IntersectionInfo *mininfo=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
					mininfo->intersect_id=majRoadnet->nodelist->nelems-1;
					mininfo->lineseg_id = minlinesegid;
					majRoadintersectinfo_min[newintersect->minorline_id]->sorted_add(mininfo);
				}
			}
		}
}



void compute_majRoad_intersects()
{
	/*   Initialize the intersection list in each cell   */
	int i;
	QuadCell *face;
	for(i=0;i<quadmesh->nfaces;i++)
	{
		face=quadmesh->quadcells[i];

		/*   initialize the intersection lists  */
		if(face->intersectlist!=NULL)
		{
			free(face->intersectlist);
			face->intersectlist=NULL;
		}
		face->nintersects=0;

		/*   initialize the edge lists   */
		if(face->streetgraphedgelist != NULL)
		{
			free(face->streetgraphedgelist);
			face->streetgraphedgelist=NULL;
		}
		face->nstreetgraphedges=0;
	}

	/*we first locate the cells that may contain intersections*/
	for(i=0; i<quadmesh->nfaces; i++)
	{
		face = quadmesh->quadcells[i];

		if(face->hasmajor && face->hasminor)
			compute_majRoad_intersects_in_cell(i);
	}

	/*    second, 
	      record the end point of the tensor lines
	*/
	
	Trajectory *temp;
	icVector2 dist;
	Intersection *intersect;

	/*          major tensor lines          */

	for(i=0; i<major_level1->evenstreamlines->ntrajs; i++)
	{
		/*record the start point*/
		temp = major_level1->evenstreamlines->trajs[i];

		/*if it is exactly (close enough to) the first intersection*/
		if(majRoadintersectinfo_maj[i]->nelems>0)
		{
			intersect=majRoadnet->nodelist->intersects[
				majRoadintersectinfo_maj[i]->infolist[0]->intersect_id];
			dist.entry[0]=temp->linesegs[0].gstart[0]-intersect->gpos[0];
			dist.entry[1]=temp->linesegs[0].gstart[1]-intersect->gpos[1];
			if(length(dist)<1e-3) continue;
		}

		Intersection *newintersect=(Intersection*)malloc(sizeof(Intersection));
		newintersect->gpos[0]=temp->linesegs[0].gstart[0];
		newintersect->gpos[1]=temp->linesegs[0].gstart[1];
		newintersect->cellid=temp->linesegs[0].Triangle_ID;
		newintersect->majorline_id=temp->index;
		newintersect->minorline_id=-1;
		newintersect->majlineseg=0;
		newintersect->nadjedges=0;
		newintersect->adj_edges=NULL;
		newintersect->endpt=true;
		newintersect->inside_region=false;
		newintersect->deleted=false;
		newintersect->intersect_type=0;

		majRoadnet->nodelist->addNew(newintersect);
		add_to_cell_intersectlist(newintersect->cellid, newintersect->index, false);

		/*compute the edges incident to it*/
		IntersectionInfo *majinfo1=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
		majinfo1->intersect_id=majRoadnet->nodelist->nelems-1;
		majinfo1->lineseg_id = 0;
		majRoadintersectinfo_maj[i]->sorted_add(majinfo1);
		
		/*record the end point*/
				
		/*if it is exactly (close enough to) the first intersection*/
		if(majRoadintersectinfo_maj[i]->nelems>0)
		{
			intersect=majRoadnet->nodelist->intersects[
				majRoadintersectinfo_maj[i]->infolist[
					majRoadintersectinfo_maj[i]->nelems-1]->intersect_id];
			dist.entry[0]=temp->linesegs[temp->nlinesegs-1].gend[0]-intersect->gpos[0];
			dist.entry[1]=temp->linesegs[temp->nlinesegs-1].gend[1]-intersect->gpos[1];
			if(length(dist)<1e-3) continue;
		}


		temp = major_level1->evenstreamlines->trajs[i];
		Intersection *newintersect2=(Intersection*)malloc(sizeof(Intersection));
		newintersect2->gpos[0]=temp->linesegs[temp->nlinesegs-1].gend[0];
		newintersect2->gpos[1]=temp->linesegs[temp->nlinesegs-1].gend[1];
		newintersect2->cellid=temp->linesegs[temp->nlinesegs-1].Triangle_ID;
		newintersect2->majorline_id=temp->index;
		newintersect2->minorline_id=-1;
		newintersect2->majlineseg=temp->nlinesegs-1;
		newintersect2->nadjedges=0;
		newintersect2->adj_edges=NULL;
		newintersect2->endpt=true;
		newintersect2->inside_region=false;
		newintersect2->deleted=false;
		newintersect2->intersect_type=0;

		majRoadnet->nodelist->addNew(newintersect2);

		add_to_cell_intersectlist(newintersect2->cellid, newintersect2->index, false);

		/*compute the edges incident to it*/
		IntersectionInfo *majinfo2=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
		majinfo2->intersect_id=majRoadnet->nodelist->nelems-1;
		majinfo2->lineseg_id = temp->nlinesegs-1;
		majRoadintersectinfo_maj[i]->sorted_add(majinfo2);
	}

	/*          minor tensor lines          */
	for(i=0; i<minor_level1->evenstreamlines->ntrajs; i++)
	{
		/*record the start point*/
		temp = minor_level1->evenstreamlines->trajs[i];
		/*if it is exactly (close enough to) the first intersection*/

		if(majRoadintersectinfo_min[i]->nelems>0)
		{
			intersect=majRoadnet->nodelist->intersects[
				majRoadintersectinfo_min[i]->infolist[0]->intersect_id];
			dist.entry[0]=temp->linesegs[0].gstart[0]-intersect->gpos[0];
			dist.entry[1]=temp->linesegs[0].gstart[1]-intersect->gpos[1];
			if(length(dist)<1e-3) continue;
		}

		Intersection *newintersect=(Intersection*)malloc(sizeof(Intersection));
		newintersect->gpos[0]=temp->linesegs[0].gstart[0];
		newintersect->gpos[1]=temp->linesegs[0].gstart[1];
		newintersect->cellid=temp->linesegs[0].Triangle_ID;
		newintersect->minorline_id=temp->index;
		newintersect->majorline_id=-1;
		newintersect->minlineseg=0;
		newintersect->nadjedges=0;
		newintersect->adj_edges=NULL;
		newintersect->endpt=true;
		newintersect->inside_region=false;
		newintersect->deleted=false;


		majRoadnet->nodelist->addNew(newintersect);
		
		add_to_cell_intersectlist(newintersect->cellid, newintersect->index, false);

		/*compute the edges incident to it*/
		IntersectionInfo *mininfo1=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
		mininfo1->intersect_id=majRoadnet->nodelist->nelems-1;
		mininfo1->lineseg_id = 0;
		majRoadintersectinfo_min[i]->sorted_add(mininfo1);
		
		/*record the end point*/
		/*if it is exactly (close enough to) the first intersection*/

		if(majRoadintersectinfo_min[i]->nelems>0)
		{
			intersect=majRoadnet->nodelist->intersects[
				majRoadintersectinfo_min[i]->infolist[
					majRoadintersectinfo_min[i]->nelems-1]->intersect_id];
			dist.entry[0]=temp->linesegs[temp->nlinesegs-1].gend[0]-intersect->gpos[0];
			dist.entry[1]=temp->linesegs[temp->nlinesegs-1].gend[1]-intersect->gpos[1];
			if(length(dist)<1e-3) continue;
		}

		temp = minor_level1->evenstreamlines->trajs[i];
		Intersection *newintersect2=(Intersection*)malloc(sizeof(Intersection));
		newintersect2->gpos[0]=temp->linesegs[temp->nlinesegs-1].gend[0];
		newintersect2->gpos[1]=temp->linesegs[temp->nlinesegs-1].gend[1];
		newintersect2->cellid=temp->linesegs[temp->nlinesegs-1].Triangle_ID;
		newintersect2->minorline_id=temp->index;
		newintersect2->majorline_id=-1;
		newintersect2->minlineseg=temp->nlinesegs-1;
		newintersect2->nadjedges=0;
		newintersect2->adj_edges=NULL;
		newintersect2->endpt=true;
		newintersect2->inside_region=false;
		newintersect2->deleted=false;

		majRoadnet->nodelist->addNew(newintersect2);
		
		add_to_cell_intersectlist(newintersect2->cellid, newintersect2->index, false);

		/*compute the edges incident to it*/
		IntersectionInfo *mininfo2=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
		mininfo2->intersect_id=majRoadnet->nodelist->nelems-1;
		mininfo2->lineseg_id = temp->nlinesegs-1;
		majRoadintersectinfo_min[i]->sorted_add(mininfo2);
	}
}


/*
   Judge whether we already have an edge between these two intersections or not
   This routine will be used to avoid the small loop in the major road network
*/
bool has_edge_between_majRoad(int intersect1, int intersect2)
{
	int i;
	StreetGraphEdge *edge;
	for(i=0;i<majRoadnet->edgelist->nedges;i++)
	{
		edge=majRoadnet->edgelist->edges[i];

		if((edge->node_index1==intersect1 && edge->node_index2==intersect2)
			||(edge->node_index1==intersect2 && edge->node_index2==intersect1))
			return true;
	}
	return false;
}


bool has_edge_between_minRoad(int intersect1, int intersect2)
{
	int i;
	StreetGraphEdge *edge;
	for(i=0;i<streetnet->edgelist->nedges;i++)
	{
		edge=streetnet->edgelist->edges[i];

		if((edge->node_index1==intersect1 && edge->node_index2==intersect2)
			||(edge->node_index1==intersect2 && edge->node_index2==intersect1))
			return true;
	}
	return false;
}


void search_for_connection_majRoad()
{
	/*we can simply search the two groups of the intersection lists once
	and construct the edges*/
	int i, j;
	TensorLineIntersectionInfoList *infolist;
	IntersectionInfo *info, *infonext;

	int start_lineseg, end_lineseg;
	double startp[2], endp[2];
	Trajectory *traj;

	/*search major lines first*/
	for(i=0; i<major_level1->evenstreamlines->ntrajs; i++)
	{
		infolist=majRoadintersectinfo_maj[i];
		traj=major_level1->evenstreamlines->trajs[i];

		for(j=0; j<infolist->nelems-1; j++)
		{
			info = infolist->infolist[j];
			infonext = infolist->infolist[j+1];

			/*  we are NOT considering isolated edges now 12/27/2007  */
			if(majRoadnet->nodelist->intersects[info->intersect_id]->endpt
				&&majRoadnet->nodelist->intersects[infonext->intersect_id]->endpt)
				continue;

			/*  If this is a repeated edge, don't add it 1/10/2008
			    May contain bugs
			*/
			if(has_edge_between_majRoad(info->intersect_id, infonext->intersect_id))
				continue;

			/*we use new data structure for the street network 11/06/2007*/
			StreetGraphEdge *newedge = (StreetGraphEdge*)malloc(sizeof(StreetGraphEdge));
			newedge->node_index1 = info->intersect_id;
			newedge->node_index2 = infonext->intersect_id;
			newedge->cancel = newedge->visited = false;
			newedge->inter_pts = NULL; //we will add the intermediate points later
			newedge->ninter_pts = 0;
			newedge->roadtype=traj->roadtype;
			//newedge->deleted=false;
			majRoadnet->edgelist->append(newedge);
			
			/* resample the line and record the information into the 
			corresponding line  11/26/2007 */
			start_lineseg=info->lineseg_id;
			end_lineseg=infonext->lineseg_id;

			startp[0]=majRoadnet->nodelist->intersects[info->intersect_id]->gpos[0];
			startp[1]=majRoadnet->nodelist->intersects[info->intersect_id]->gpos[1];
			endp[0]=majRoadnet->nodelist->intersects[infonext->intersect_id]->gpos[0];
			endp[1]=majRoadnet->nodelist->intersects[infonext->intersect_id]->gpos[1];

			sample_along_tensorline_from_to(traj, start_lineseg, startp,
									 end_lineseg, endp, 
									 newedge->index, majRoadnet);

			/**/
			add_edge_to_majRoad_intersectnode(info->intersect_id, newedge->index);
			add_edge_to_majRoad_intersectnode(infonext->intersect_id, newedge->index);
		}

		/*  if the tensor line is closed, add one more edge */
		if(major_level1->evenstreamlines->trajs[i]->closed)
		{
			info = infolist->infolist[infolist->nelems-1];
			infonext = infolist->infolist[0];

			if(!has_edge_between_majRoad(info->intersect_id, infonext->intersect_id))
			{
				/*we use new data structure for the street network 11/06/2007*/
				StreetGraphEdge *onenewedge = (StreetGraphEdge*)malloc(sizeof(StreetGraphEdge));
				onenewedge->node_index1 = info->intersect_id;
				onenewedge->node_index2 = infonext->intersect_id;
				onenewedge->cancel = onenewedge->visited = false;
				onenewedge->inter_pts = NULL; //we will add the intermediate points later
				onenewedge->ninter_pts = 0;
				onenewedge->roadtype=traj->roadtype;
				majRoadnet->edgelist->append(onenewedge);
				
				/* resample the line and record the information into the 
				corresponding line  11/26/2007 */
				start_lineseg=info->lineseg_id;
				end_lineseg=infonext->lineseg_id;

				startp[0]=majRoadnet->nodelist->intersects[info->intersect_id]->gpos[0];
				startp[1]=majRoadnet->nodelist->intersects[info->intersect_id]->gpos[1];
				int start_cell=get_cellID_givencoords(startp[0], startp[1]);
				endp[0]=majRoadnet->nodelist->intersects[infonext->intersect_id]->gpos[0];
				endp[1]=majRoadnet->nodelist->intersects[infonext->intersect_id]->gpos[1];
				int end_cell=get_cellID_givencoords(endp[0], endp[1]);

				Trajectory *tempTraj=new Trajectory(-1);
				get_linesegs_anytwopts(startp,start_cell, endp,end_cell, tempTraj, 0, 10);
				onenewedge->inter_pts=(Point **)malloc(sizeof(Point *)* (tempTraj->nlinesegs+1));

				int m;
				int pre_cell=-1;
				for(m=0;m<tempTraj->nlinesegs;m++)
				{
					onenewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
					onenewedge->inter_pts[m]->x=tempTraj->linesegs[m].gstart[0];
					onenewedge->inter_pts[m]->y=tempTraj->linesegs[m].gstart[1];
					onenewedge->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge->inter_pts[m]->x,
						onenewedge->inter_pts[m]->y);
					if(onenewedge->inter_pts[m]->cellid!=pre_cell)
					{
						add_to_edgelist_one_cell(onenewedge->inter_pts[m]->cellid, i);
						pre_cell=onenewedge->inter_pts[m]->cellid;
					}
				}

				onenewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
				onenewedge->inter_pts[m]->x=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[0];
				onenewedge->inter_pts[m]->y=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[1];
				onenewedge->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge->inter_pts[m]->x,
					onenewedge->inter_pts[m]->y);
				onenewedge->ninter_pts=tempTraj->nlinesegs+1;
				delete tempTraj;

				/**/
				add_edge_to_majRoad_intersectnode(info->intersect_id, onenewedge->index);
				add_edge_to_majRoad_intersectnode(infonext->intersect_id, onenewedge->index);
			}
		}
	}

	/*search minor lines first*/
	for(i=0; i<minor_level1->evenstreamlines->ntrajs; i++)
	{
		infolist=majRoadintersectinfo_min[i];
		traj=minor_level1->evenstreamlines->trajs[i];

		for(j=0; j<infolist->nelems-1; j++)
		{
			info = infolist->infolist[j];
			infonext = infolist->infolist[j+1];

			
			/*  we are NOT considering isolated edges now 12/27/2007  */
			if(majRoadnet->nodelist->intersects[info->intersect_id]->endpt
				&&majRoadnet->nodelist->intersects[infonext->intersect_id]->endpt)
				continue;
			
			/*  If this is a repeated edge, don't add it 1/10/2008
			    May contain bugs
			*/
			if(has_edge_between_majRoad(info->intersect_id, infonext->intersect_id))
				continue;

			/*we use new data structure for the street network 11/06/2007*/
			StreetGraphEdge *newedge = (StreetGraphEdge*)malloc(sizeof(StreetGraphEdge));
			newedge->node_index1 = info->intersect_id;
			newedge->node_index2 = infonext->intersect_id;
			newedge->cancel = newedge->visited = false;
			newedge->inter_pts = NULL; //we will add the intermediate points later
			newedge->ninter_pts = 0;
			newedge->roadtype=traj->roadtype;
			majRoadnet->edgelist->append(newedge);

			/* resample the line and record the information into the 
			corresponding line  11/26/2007 */
			start_lineseg=info->lineseg_id;
			end_lineseg=infonext->lineseg_id;

			startp[0]=majRoadnet->nodelist->intersects[info->intersect_id]->gpos[0];
			startp[1]=majRoadnet->nodelist->intersects[info->intersect_id]->gpos[1];
			endp[0]=majRoadnet->nodelist->intersects[infonext->intersect_id]->gpos[0];
			endp[1]=majRoadnet->nodelist->intersects[infonext->intersect_id]->gpos[1];

			sample_along_tensorline_from_to(traj, start_lineseg, startp,
									 end_lineseg, endp, 
									 newedge->index, majRoadnet);

			/**/
			add_edge_to_majRoad_intersectnode(info->intersect_id, newedge->index);
			add_edge_to_majRoad_intersectnode(infonext->intersect_id, newedge->index);
		}
		/*  if the tensor line is closed, add one more edge */
		if(minor_level1->evenstreamlines->trajs[i]->closed)
		{
			info = infolist->infolist[infolist->nelems-1];
			infonext = infolist->infolist[0];

			if(!has_edge_between_majRoad(info->intersect_id, infonext->intersect_id))
			{
				/*we use new data structure for the street network 11/06/2007*/
				StreetGraphEdge *onenewedge = (StreetGraphEdge*)malloc(sizeof(StreetGraphEdge));
				onenewedge->node_index1 = info->intersect_id;
				onenewedge->node_index2 = infonext->intersect_id;
				onenewedge->cancel = onenewedge->visited = false;
				onenewedge->inter_pts = NULL; //we will add the intermediate points later
				onenewedge->ninter_pts = 0;
				onenewedge->roadtype=traj->roadtype;
				majRoadnet->edgelist->append(onenewedge);
				
				/* resample the line and record the information into the 
				corresponding line  11/26/2007 */
				start_lineseg=info->lineseg_id;
				end_lineseg=infonext->lineseg_id;

				startp[0]=majRoadnet->nodelist->intersects[info->intersect_id]->gpos[0];
				startp[1]=majRoadnet->nodelist->intersects[info->intersect_id]->gpos[1];
				int start_cell=get_cellID_givencoords(startp[0], startp[1]);
				endp[0]=majRoadnet->nodelist->intersects[infonext->intersect_id]->gpos[0];
				endp[1]=majRoadnet->nodelist->intersects[infonext->intersect_id]->gpos[1];
				int end_cell=get_cellID_givencoords(endp[0], endp[1]);

				Trajectory *tempTraj=new Trajectory(-1);
				get_linesegs_anytwopts(startp,start_cell, endp,end_cell, tempTraj, 0, 10);
				onenewedge->inter_pts=(Point **)malloc(sizeof(Point *)* (tempTraj->nlinesegs+1));

				int m;
				int pre_cell=-1;
				for(m=0;m<tempTraj->nlinesegs;m++)
				{
					onenewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
					onenewedge->inter_pts[m]->x=tempTraj->linesegs[m].gstart[0];
					onenewedge->inter_pts[m]->y=tempTraj->linesegs[m].gstart[1];
					onenewedge->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge->inter_pts[m]->x,
						onenewedge->inter_pts[m]->y);
					if(onenewedge->inter_pts[m]->cellid!=pre_cell)
					{
						add_to_edgelist_one_cell(onenewedge->inter_pts[m]->cellid, i);
						pre_cell=onenewedge->inter_pts[m]->cellid;
					}
				}

				onenewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
				onenewedge->inter_pts[m]->x=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[0];
				onenewedge->inter_pts[m]->y=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[1];
				onenewedge->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge->inter_pts[m]->x,
					onenewedge->inter_pts[m]->y);
				onenewedge->ninter_pts=tempTraj->nlinesegs+1;
				delete tempTraj;

				/**/
				add_edge_to_majRoad_intersectnode(info->intersect_id, onenewedge->index);
				add_edge_to_majRoad_intersectnode(infonext->intersect_id, onenewedge->index);
			}
		}
	}
}




void init_majRoadinfo_incells()
{
	int i, j;
	QuadCell *face;
	for(i=0; i<quadmesh->nfaces;i++)
	{
		/**/
		face=quadmesh->quadcells[i];

		if(face->majorlines!=NULL)
		{
			delete face->majorlines;
			face->majorlines=NULL;
		}
		
		if(face->minorlines!=NULL)
		{
			delete face->minorlines;
			face->minorlines=NULL;
		}
	}
}

void add_edge_to_majRoad_intersectnode(int node, int edgeindex)
{
	if(node < 0)
		return;

	majRoadnet->nodelist->intersects[node]->adj_edges = 
		extend_link(majRoadnet->nodelist->intersects[node]->adj_edges,
		majRoadnet->nodelist->intersects[node]->nadjedges);

	majRoadnet->nodelist->intersects[node]->adj_edges[
		majRoadnet->nodelist->intersects[node]->nadjedges] = edgeindex;
	majRoadnet->nodelist->intersects[node]->nadjedges++;
}


/*
    we connect one major road to existing road(s) using the straight line direction
	NOTE: this routine will be called during tracing!
	The basic idea is we search a small disc whose radius is determined by
	"search_dist" and center is "pt".
*/
void cal_one_euclidean_disc(int triangle, double p[3], double dsep, double discsize, 
					DynList_Int *trianglelist)
{
	////we suppose the global point always falls in the triangle
	QuadCell *face = quadmesh->quadcells[triangle];
	QuadVertex *vert;
	QuadEdge *cur_edge;
	icVector3 dis;

	////Find all the triangle
	trianglelist->nelems = 0;
	trianglelist->add_New(triangle);
	int cur_id = 0;
	int i, j;

	while(cur_id < trianglelist->nelems)
	{
		face = quadmesh->quadcells[trianglelist->elems[cur_id]];

		for(i = 0; i < face->nverts; i++)
		{
			vert = quadmesh->quad_verts[face->verts[i]];

			if(vert->visited)
				continue;

			dis.entry[0] = vert->x - p[0];
			dis.entry[1] = vert->y - p[1];

			vert->distance = length(dis); /*compute the Euclidean distance*/

			vert->visited = true; //set the visited flag
			if(vert->distance > discsize*dsep)
				continue;

			/*if the distance between the vertex and the input point is smaller than the threshold
			we need to add all its adjacent triangles into the triangle list
			*/
			for(j = 0; j < vert->ncells; j++)
			{
				QuadCell *c = vert->cells[j];

				if(c->index < 0)   //reach the boundary!
					continue;
				trianglelist->add_New(c->index);
			}
		}

		cur_id++;
	}

	for(i=0; i<trianglelist->nelems; i++)
	{
		face = quadmesh->quadcells[trianglelist->elems[i]];
		face->visited = false;
		for(j=0; j<face->nverts; j++)
		{
			vert = quadmesh->quad_verts[face->verts[j]];
			vert->visited = false;
		}
	}
}

bool is_at_validDir(int cell, double pt[2], icVector2 trajDir)
{
	int i;
	icVector2 otherDir;
	QuadCell *face=quadmesh->quadcells[cell];
	QuadVertex *v;
	//int ninvalidverts=0;
	for(i=0;i<face->nverts;i++)
	{
		v=quadmesh->quad_verts[face->verts[i]];
		otherDir.entry[0]=v->x-pt[0];
		otherDir.entry[1]=v->y-pt[1];

		if(dot(otherDir, trajDir)>=0)
			return true;
	}

	return false;
}

void cal_dist_to_lines(double pt[2], int curCell, icVector2 trajDir, double cosang,
					   double appro_inter[2], double &dist, bool majormin)
{
	int i, j;
	dist=1.e50;
	double distToLine, distToLineSeg;
	QuadCell *face=quadmesh->quadcells[curCell];
	LineInfo *lineinfo;
	Trajectory *traj;
	double a[2], b[2];
	icVector2 Dir;
	//normalize(trajDir);
	//double cosang=cos(ang);

	if(!majormin)
	{
		lineinfo=face->majorlines;
	}
	else
	{
		lineinfo=face->minorlines;
	}

	for(i=0;i<lineinfo->nlines;i++)
	{
		/*  obtain the corresponding trajectory  */
		if(!majormin)
			traj=major_level1->evenstreamlines->trajs[lineinfo->lines[i]->whichtraj];
		else
			traj=minor_level1->evenstreamlines->trajs[lineinfo->lines[i]->whichtraj];

		for(j=max(0,lineinfo->lines[i]->start-1); 
			j<=min(lineinfo->lines[i]->end, traj->nlinesegs-1); j++)
		{
			/*  compute the distance to the corresponding line segment */
			a[0]=traj->linesegs[j].gstart[0];
			a[1]=traj->linesegs[j].gstart[1];
			b[0]=traj->linesegs[j].gend[0];
			b[1]=traj->linesegs[j].gend[1];

			DistanceFromLine(pt[0],pt[1],  a[0],a[1],  b[0],b[1],
				distToLineSeg, distToLine);

			/*  judge whether it can be a potential intersection  */

			if(distToLineSeg<dist)
			{
				double temp[2]={(a[0]+b[0])/2., (a[1]+b[1])/2.};
				Dir.entry[0]=temp[0]-pt[0];
				Dir.entry[1]=temp[1]-pt[1];
				normalize(Dir);

				/*  after adding "|| distToLineSeg<quadmesh->xinterval"
				    this may contain bugs
				*/
				if(dot(Dir, trajDir)>=cosang /*|| distToLineSeg<quadmesh->xinterval*/)
				{
					dist=distToLineSeg;
					/*  modified at 1/10/2008  bug  */
					if(j==0)
					{
						appro_inter[0]=traj->linesegs[0].gstart[0];
						appro_inter[1]=traj->linesegs[0].gstart[1];
					}
					else if(j==traj->nlinesegs-1)
					{
						appro_inter[0]=traj->linesegs[j].gend[0];
						appro_inter[1]=traj->linesegs[j].gend[1];
					}
					else
					{
						appro_inter[0]=temp[0];
						appro_inter[1]=temp[1];
					}
				}

			}

		}
	}
}

bool find_other_nearby_traj(double pt[2], int cell, icVector2 trajDir, bool majormin,
							int &other_trajID, int &other_cellID, double otherpt[2],
							bool &other_majormin, double search_dist, double ang)
{
	DynList_Int *onedisc=new DynList_Int();
	cal_one_euclidean_disc(cell, pt, search_dist, 1.1, onedisc);

	/*  Search the cells in the disc and see whether there is another major roads */
	int i;
	QuadCell *cur_face;
	double smallest_dist=1.e50;
	double temp_dist=1.e50;
	double cosang=cos(ang);
	normalize(trajDir);
	double temp_otherpt[2];
	bool found=false;

	for(i=0;i<onedisc->nelems;i++)
	{
		if(!is_at_validDir(cell, pt, trajDir))
			continue;

		cur_face=quadmesh->quadcells[onedisc->elems[i]];

		//if(!cur_face->hasmajor && !cur_face->hasminor)
		//	continue;

		//if(majormin) /* we may according to the type of current trajectory to decide the order */

		/*  consider the major direction first  */
		if(cur_face->hasmajor)
		{
			/**/
			cal_dist_to_lines(pt, cur_face->index, trajDir, cosang, temp_otherpt, temp_dist,
				false);
			if(temp_dist<smallest_dist)
			{
				smallest_dist=temp_dist;
				otherpt[0]=temp_otherpt[0];
				otherpt[1]=temp_otherpt[1];
				other_cellID=onedisc->elems[i];
				other_majormin=false;
				found=true;
			}
		}
		if(cur_face->hasminor)
		{
			cal_dist_to_lines(pt, cur_face->index, trajDir, cosang, temp_otherpt, temp_dist,
				true);
			if(temp_dist<smallest_dist)
			{
				smallest_dist=temp_dist;
				otherpt[0]=temp_otherpt[0];
				otherpt[1]=temp_otherpt[1];
				other_cellID=onedisc->elems[i];
				other_majormin=true;
				found=true;
			}
		}

	}

	delete onedisc;
	return found;
}


/*
    find the nearby existing tensor line through tracing
*/

//bool trace_to_comp_nearby_traj(double pt[2], int cell, bool majormin,
//							int &other_trajID, int &other_cellID, double otherpt[2],
//							bool &other_majormin, double search_dist, int MaxCells,
//							double min_search,	Trajectory *trajSeg)
bool trace_to_comp_nearby_traj(double pt[2], int cell, bool majormin,
							double search_dist, int MaxCells, double min_waterwidth,
							double min_search,	double &dist_did, Trajectory *trajSeg
							/*int &flag*/)
{
	int i;
	double cur_trace_dist=0;
	int pre_face, cur_face;
	pre_face=cur_face=cell;
	int flag=-1;

	double globalp[2]={pt[0], pt[1]};
	int other_trajID;
	bool other_majormin;
	int other_cellID;
	double otherpt[2]={0.};


	/*  variables for pixel level judgement of "inland" or not  */
	double xstart=quadmesh->xstart;
	double ystart=quadmesh->ystart;
	double xrang=quadmesh->xend-quadmesh->xstart;
	double yrang=quadmesh->yend-quadmesh->ystart;
	double dx=xrang/511;
	double dy=yrang/511;

	//DynList_Int *onedisc=new DynList_Int();

	dist_did=0;

	euler_stepsize = quadmesh->xinterval/10;

	for(i=0;i<MaxCells;i++)
	{
		if(cur_face < 0 || cur_face >= quadmesh->nfaces)
		{
			flag=0;
			return true;
		}


		/*it reaches any boundary, we should stop as well*/
		if(globalp[0]>=quadmesh->xend-1.e-8||globalp[0]<=quadmesh->xstart+1.e-8
			||globalp[1]>=quadmesh->yend-1.e-8||globalp[1]<=quadmesh->ystart+1.e-8)
		{
			if(globalp[0]>=quadmesh->xend-1.e-8) globalp[0]=quadmesh->xend;
			else if(globalp[0]<=quadmesh->xstart+1.e-8) globalp[0]=quadmesh->xstart;
			if(globalp[1]>=quadmesh->yend-1.e-8) globalp[1]=quadmesh->yend;
			else if(globalp[1]<=quadmesh->ystart+1.e-8) globalp[1]=quadmesh->ystart;
			//break;
			return true;
		}

		/*  search the nearby triangles  */

		//if(majormin&&quadmesh->quadcells[cur_face]->hasmajor
		//{
		//}
		//else if(!majormin&&quadmesh
		
		icVector2 trajDir = tenline_dir_global;
		//globalp[0]+=1.e-7*trajDir.entry[0];
		//globalp[1]+=1.e-7*trajDir.entry[1];
		if(find_other_nearby_traj(globalp, cur_face, trajDir, majormin,
							other_trajID, other_cellID, otherpt,
							other_majormin, /*search_dist/3.*/min_search, M_PI/6.))
		{
			/*  connect to the neary by trajectory  */
		//otherpt[0]+=1.e-7*trajDir.entry[0];
		//otherpt[1]+=1.e-7*trajDir.entry[1];
			get_linesegs_anytwopts(globalp, cur_face,   otherpt, other_cellID,   trajSeg, 0, 200);
			dist_did=trajSeg->get_length();
			return true;
		}

		/**************************************************************/

		if(is_inveg_cell_weak(cur_face))
			break;

		/*for loading the map 10/24/2007*/
		//if(is_not_inland(cur_face))
		if(is_not_inland_cell_weak(cur_face))
		//if(!is_inland_pixel(globalp[0],globalp[1], xstart, xrang, ystart, yrang, dx, dy,
		//		fittedmap1, 512))
		{
			if(!sharedvars.AllowMajRoadCrossRiverOn /*&& !sharedvars.AllowMajRoadFollowBoundaryOn*/)
				break;

			/*   deal with water region   */
			double ce[2];
			icVector2 norm;
			//norm.entry[0]=-tenline_dir_global.entry[1];
			//norm.entry[1]=tenline_dir_global.entry[0];
			norm=tenline_dir_global;
			normalize(norm);
			ce[0]=globalp[0]+min_waterwidth*quadmesh->xinterval*norm.entry[0];
			ce[1]=globalp[1]+min_waterwidth*quadmesh->xinterval*norm.entry[1];
			int ce_cell=get_cellID_givencoords(ce[0], ce[1]);

			if(sharedvars.AllowMajRoadCrossRiverOn &&
				is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
					fittedmap1, 512)/*&& ang_traceDir_boundDir<cross_angle*/)
			{
				get_linesegs_anytwopts(globalp, cur_face,   ce, ce_cell,   trajSeg, 0, 200);
				globalp[0]=ce[0];
				globalp[1]=ce[1];
				cur_face=ce_cell;
			}

			/*if(!is_inland_pixel(globalp[0],globalp[1], xstart, xrang, ystart, yrang, dx, dy,
					fittedmap1, 512))*/
			else
				break;
		}
		
		pre_face = cur_face;

		if(!majormin)
			cur_face = trajSeg->trace_in_quad(cur_face, globalp, 0, flag);
		else
			cur_face = trajSeg->trace_in_quad(cur_face, globalp, 1, flag);

		//if this is the first line segment, there may be a case that the start point is on the edge !!!
		
		cur_trace_dist = trajSeg->get_length();

		//if(cur_trace_dist>search_dist)
		//	break;

		if(flag == 1 || flag == 4 || flag == 2  || pre_face == cur_face ) 
		{
			break;
		}
	}

	//delete onedisc;
	dist_did=trajSeg->get_length();
	return false;
}





/**************************************************************************************/
/*
   We try to implement a similar function as major road placement
   for placing minor roads with as fewer dead ends as possible
*/

void cal_dist_to_lines_minRoad(double pt[2], int curCell, icVector2 trajDir, double cosang,
					   double appro_inter[2], double &dist, bool majormin)
{
	int i, j;
	dist=1.e50;
	double distToLine, distToLineSeg;
	QuadCell *face=quadmesh->quadcells[curCell];
	LineInfo *lineinfo;
	Trajectory *traj;
	double a[2], b[2];
	icVector2 Dir;
	//normalize(trajDir);
	//double cosang=cos(ang);

	if(!majormin)
	{
		lineinfo=face->majorlines;
	}
	else
	{
		lineinfo=face->minorlines;
	}

	for(i=0;i<lineinfo->nlines;i++)
	{
		/*  obtain the corresponding trajectory  */
		if(!majormin)
			traj=major->evenstreamlines->trajs[lineinfo->lines[i]->whichtraj];
		else
			traj=minor->evenstreamlines->trajs[lineinfo->lines[i]->whichtraj];

		for(j=max(0,lineinfo->lines[i]->start-1); 
			j<=min(lineinfo->lines[i]->end, traj->nlinesegs-1); j++)
		{
			/*  compute the distance to the corresponding line segment */
			a[0]=traj->linesegs[j].gstart[0];
			a[1]=traj->linesegs[j].gstart[1];
			b[0]=traj->linesegs[j].gend[0];
			b[1]=traj->linesegs[j].gend[1];

			DistanceFromLine(pt[0],pt[1],  a[0],a[1],  b[0],b[1],
				distToLineSeg, distToLine);

			/*  judge whether it can be a potential intersection  */

			if(distToLineSeg<dist)
			{
				double temp[2]={(a[0]+b[0])/2., (a[1]+b[1])/2.};
				Dir.entry[0]=temp[0]-pt[0];
				Dir.entry[1]=temp[1]-pt[1];
				normalize(Dir);

				/*  after adding "|| distToLineSeg<quadmesh->xinterval"
				    this may contain bugs
				*/
				if(dot(Dir, trajDir)>=cosang 
					/*|| (distToLineSeg<quadmesh->xinterval&&dot(Dir, trajDir)>0)*/)
				{
					dist=distToLineSeg;
					/*  modified at 1/10/2008  bug  */
					if(j==0)
					{
						appro_inter[0]=traj->linesegs[0].gstart[0];
						appro_inter[1]=traj->linesegs[0].gstart[1];
					}
					else if(j==traj->nlinesegs-1)
					{
						appro_inter[0]=traj->linesegs[j].gend[0];
						appro_inter[1]=traj->linesegs[j].gend[1];
					}
					else
					{
						appro_inter[0]=temp[0];
						appro_inter[1]=temp[1];
					}
				}

			}

		}
	}
}

bool find_other_nearby_minRoad(double pt[2], int cell, icVector2 trajDir, bool majormin,
							int &other_trajID, int &other_cellID, double otherpt[2],
							bool &other_majormin, double search_dist, double ang)
{
	DynList_Int *onedisc=new DynList_Int();
	cal_one_euclidean_disc(cell, pt, search_dist, 1.1, onedisc);

	/*  Search the cells in the disc and see whether there is another major roads */
	int i;
	QuadCell *cur_face;
	double smallest_dist=1.e50;
	double temp_dist=1.e50;
	double cosang=cos(ang);
	normalize(trajDir);
	double temp_otherpt[2];
	bool found=false;

	for(i=0;i<onedisc->nelems;i++)
	{
		if(!is_at_validDir(cell, pt, trajDir))
			continue;

		cur_face=quadmesh->quadcells[onedisc->elems[i]];

		//if(!cur_face->hasmajor && !cur_face->hasminor)
		//	continue;

		//if(majormin) /* we may according to the type of current trajectory to decide the order */

		/*  consider the major direction first  */
		if(cur_face->hasmajor)
		{
			/**/
			cal_dist_to_lines_minRoad(pt, cur_face->index, trajDir, cosang, temp_otherpt, temp_dist,
				false);
			if(temp_dist<smallest_dist)
			{
				smallest_dist=temp_dist;
				otherpt[0]=temp_otherpt[0];
				otherpt[1]=temp_otherpt[1];
				other_cellID=onedisc->elems[i];
				other_majormin=false;
				found=true;
			}
		}
		if(cur_face->hasminor)
		{
			cal_dist_to_lines_minRoad(pt, cur_face->index, trajDir, cosang, temp_otherpt, temp_dist,
				true);
			if(temp_dist<smallest_dist)
			{
				smallest_dist=temp_dist;
				otherpt[0]=temp_otherpt[0];
				otherpt[1]=temp_otherpt[1];
				other_cellID=onedisc->elems[i];
				other_majormin=true;
				found=true;
			}
		}

	}

	delete onedisc;
	return found;
}



bool trace_to_comp_nearby_minRoad(double pt[2], int cell, bool majormin,
							int MaxCells, double max_traceDist,
							double min_search,	double &dist_did, Trajectory *trajSeg)
{
	int i;
	double cur_trace_dist=0;
	int pre_face, cur_face;
	pre_face=cur_face=cell;
	int flag=-1;

	double globalp[2]={pt[0], pt[1]};
	int other_trajID;
	bool other_majormin;
	int other_cellID;
	double otherpt[2]={0.};


	dist_did=0;

	euler_stepsize = quadmesh->xinterval/10;


	for(i=0;i<MaxCells;i++)
	{
		if(cur_face < 0 || cur_face >= quadmesh->nfaces)
		{
			dist_did=trajSeg->get_length();
			return true;
		}

		/*it reaches any boundary, we should stop as well*/
		if(globalp[0]>=quadmesh->xend-1.e-8||globalp[0]<=quadmesh->xstart+1.e-8
			||globalp[1]>=quadmesh->yend-1.e-8||globalp[1]<=quadmesh->ystart+1.e-8)
		{
			if(globalp[0]>=quadmesh->xend-1.e-8) globalp[0]=quadmesh->xend;
			else if(globalp[0]<=quadmesh->xstart+1.e-8) globalp[0]=quadmesh->xstart;
			if(globalp[1]>=quadmesh->yend-1.e-8) globalp[1]=quadmesh->yend;
			else if(globalp[1]<=quadmesh->ystart+1.e-8) globalp[1]=quadmesh->ystart;
			dist_did=trajSeg->get_length();
			return true;
		}

		/*  search the nearby triangles  */
		
		icVector2 trajDir = tenline_dir_global;
		if(find_other_nearby_minRoad(globalp, cur_face, trajDir, majormin,
							other_trajID, other_cellID, otherpt,
							other_majormin, min_search, M_PI/10.))
		{
			/*  connect to the neary by trajectory  */
			get_linesegs_anytwopts(globalp, cur_face,   otherpt, other_cellID,   trajSeg, 0, 20);
			dist_did=trajSeg->get_length();
			return true;
		}

		/**************************************************************/

		if(is_inveg_cell_weak(cur_face))
			break;

		/*for loading the map 10/24/2007*/
		//if(is_not_inland(cur_face))
		if(is_not_inland_cell_weak(cur_face))
		{
			break;
		}
		
		pre_face = cur_face;

		if(!majormin)
			cur_face = trajSeg->trace_in_quad(cur_face, globalp, 0, flag);
		else
			cur_face = trajSeg->trace_in_quad(cur_face, globalp, 1, flag);

		//if this is the first line segment, there may be a case that the start point is on the edge !!!
		
		cur_trace_dist = trajSeg->get_length();

		if(cur_trace_dist>max_traceDist)
			break;

		if(flag == 1 || flag == 4 || flag == 2  || pre_face == cur_face ) 
		{
			break;
		}
	}

	dist_did=trajSeg->get_length();
	return false;
}




bool connect_one_majRoad_to_exist(int trajID, bool majormin, icVector2 trajDir,
								  double pt[2], int cell, Trajectory *trajSeg, 
								  double search_dist, double ang)
{
	Trajectory *curtraj;
	if(!majormin) /*  major direction  */
	{
		curtraj=major_level1->evenstreamlines->trajs[trajID];
	}
	else  /*  minor direction  */
	{
		curtraj=minor_level1->evenstreamlines->trajs[trajID];
	}

	/*  start tracing from pt[2] along the direction trajDir  
	    The tracing process will stop when it reaches a cell containing other existing major road
	*/
	double otherpt[2];
	bool other_majormin;
	int other_cellID;
	int other_trajID;


	if(find_other_nearby_traj(pt, cell, trajDir, majormin, other_trajID, other_cellID, 
		otherpt, other_majormin, search_dist, ang))
	{
		/*  call the routine of getting the line segments between two points  */
		/*  we adjust the start and end points a little bit  */
		pt[0]+=1.e-7*trajDir.entry[0];
		pt[1]+=1.e-7*trajDir.entry[1];
		otherpt[0]+=1.e-7*trajDir.entry[0];
		otherpt[1]+=1.e-7*trajDir.entry[1];
		get_linesegs_anytwopts(pt, cell,   otherpt, other_cellID,   trajSeg, 0, 200);
		return true;
	}

	return false;
}

/*
   We need a post process to connect the dead ends of those major roads again.
   We can borrow some idea of the removal of the dead ends of the minor roads
   Contain bugs!
*/

void connect_one_majRoad(int id, bool majormin, bool startorend)
{
	int i;

	EvenStreamlinePlace *curplace;
	Trajectory *curtraj;
	TensorLineIntersectionInfoList *theinfolist;
	Intersection *theintersect;
	LineSeg *curline;
	double pt[2];
	int cell;
	icVector2 trajDir;
	double otherpt[2];
	double distToIntersect, dist_did;
	/*  variables for pixel level judgement of "inland" or not  */
	double xstart=quadmesh->xstart;
	double ystart=quadmesh->ystart;
	double xrang=quadmesh->xend-quadmesh->xstart;
	double yrang=quadmesh->yend-quadmesh->ystart;
	double dx=xrang/511;
	double dy=yrang/511;

	icVector2 loc_dist;

	double intersect_coords[2];
	int which_cell;
	int intersect_id;

	int del_till;

	//bool notuseIntersect=false;
	bool notuseIntersect=true;  //testing

	hstep = quadmesh->xinterval/2.;
	predict_stepsize = quadmesh->xinterval/2.;
	euler_stepsize = quadmesh->xinterval/5.;

	/*create a new line segment*/
	icVector2 line_dir;


	//int nlines_backneedtoremove=curtraj->nlinesegs-del_back_till;

	/*    intialization    */

	if(!majormin)  /*  major direction  */
	{
		curplace=major_level1;
		curtraj=major_level1->evenstreamlines->trajs[id];
		theinfolist=majRoadintersectinfo_maj[id];
	}
	else           /*  minor direction  */
	{
		curplace=minor_level1;
		curtraj=minor_level1->evenstreamlines->trajs[id];
		theinfolist=majRoadintersectinfo_min[id];
	}

	/*  judge whether it is a dead end or not  */
	if(!startorend) /* start point  */
	{
		if(!majRoadnet->nodelist->intersects[theinfolist->infolist[0]->intersect_id]
			->endpt)
			return;
	}
	else
	{
		if(!majRoadnet->nodelist->intersects[theinfolist->infolist[theinfolist->nelems-1]->intersect_id]
			->endpt)
			return;
	}

	if(!startorend)/*  start point  */
	{
		curline=&curtraj->linesegs[0];
		//trajDir.entry[0]=curline->gstart[0]-curline->gend[0];
		//trajDir.entry[1]=curline->gstart[1]-curline->gend[1];
		int later_line=min(curtraj->nlinesegs-1, 1);
		trajDir.entry[0]=curline->gstart[0]-curtraj->linesegs[later_line].gend[0];
		trajDir.entry[1]=curline->gstart[1]-curtraj->linesegs[later_line].gend[1];
		
		pt[0]=curtraj->linesegs[0].gstart[0];
		pt[1]=curtraj->linesegs[0].gstart[1];
		
		if(theinfolist->nelems==1) /*  a bug here need to be fixed  */
			notuseIntersect=true;
		else
			intersect_id=theinfolist->infolist[1]->intersect_id; //avoid itself
	}
	else           /*  end point  */
	{
		curline=&curtraj->linesegs[curtraj->nlinesegs-1];
		//trajDir.entry[0]=curline->gend[0]-curline->gstart[0];
		//trajDir.entry[1]=curline->gend[1]-curline->gstart[1];
		int earlier_line=max(0, curtraj->nlinesegs-2);
		trajDir.entry[0]=curline->gend[0]-curtraj->linesegs[earlier_line].gstart[0];
		trajDir.entry[1]=curline->gend[1]-curtraj->linesegs[earlier_line].gstart[1];

		pt[0]=curtraj->linesegs[curtraj->nlinesegs-1].gend[0];
		pt[1]=curtraj->linesegs[curtraj->nlinesegs-1].gend[1];

		if(theinfolist->nelems==1) /*  a bug here need to be fixed  */
			notuseIntersect=true;
		else
			intersect_id=theinfolist->infolist[theinfolist->nelems-2]->intersect_id; //avoid itself
	}

	if(!notuseIntersect)
	theintersect=majRoadnet->nodelist->intersects[intersect_id];

	cell=get_cellID_givencoords(pt[0], pt[1]);

	tenline_dir_global=trajDir;
	normalize(tenline_dir_global);

	distToIntersect=1.e50;
	int nlines_needtoremove;

	if(!majormin)  /*  major direction  */
	{
		del_till = majRoadnet->nodelist->intersects[intersect_id]->majlineseg;
	
	}
	else           /*  minor direction  */
	{
		del_till = majRoadnet->nodelist->intersects[intersect_id]->minlineseg;
	}

	LineSeg *one_new=NULL;
	/*  create new line segment for the back track to the preceding intersection  */
	if(!startorend && !notuseIntersect)  /*  start point */
	{
		nlines_needtoremove=del_till;

		one_new=(LineSeg*)malloc(sizeof(LineSeg));
		one_new->gstart[0]=theintersect->gpos[0];
		one_new->gstart[1]=theintersect->gpos[1];
		one_new->start[0]=theintersect->gpos[0];
		one_new->start[1]=theintersect->gpos[1];
		one_new->gend[0]=curtraj->linesegs[del_till].gend[0];
		one_new->gend[1]=curtraj->linesegs[del_till].gend[1];
		one_new->end[0]=curtraj->linesegs[del_till].gend[0];
		one_new->end[1]=curtraj->linesegs[del_till].gend[1];
		one_new->Triangle_ID=curtraj->linesegs[del_till].Triangle_ID;
		line_dir.entry[0]=one_new->gend[0]-one_new->gstart[0];
		line_dir.entry[1]=one_new->gend[1]-one_new->gstart[1];
		one_new->length=length(line_dir);
	}
	else if(startorend && !notuseIntersect)
	{
		nlines_needtoremove=curtraj->nlinesegs-del_till;

		one_new=(LineSeg*)malloc(sizeof(LineSeg));
		one_new->gstart[0]=curtraj->linesegs[del_till].gstart[0];
		one_new->gstart[1]=curtraj->linesegs[del_till].gstart[1];
		one_new->start[0]=curtraj->linesegs[del_till].gstart[0];
		one_new->start[1]=curtraj->linesegs[del_till].gstart[1];
		one_new->gend[0]=theintersect->gpos[0];
		one_new->gend[1]=theintersect->gpos[1];
		one_new->end[0]=theintersect->gpos[0];
		one_new->end[1]=theintersect->gpos[1];
		one_new->Triangle_ID=curtraj->linesegs[del_till].Triangle_ID;
		line_dir.entry[0]=one_new->gend[0]-one_new->gstart[0];
		line_dir.entry[1]=one_new->gend[1]-one_new->gstart[1];
		one_new->length=length(line_dir);
	}


	if(pt[0]<=quadmesh->xstart+1.e-8||pt[0]>=quadmesh->xend-1.e-8
		||pt[1]<=quadmesh->ystart+1.e-8||pt[1]>=quadmesh->yend-1.e-8)
	{
		/* do nothing */
	}
	else
	{
		if(!notuseIntersect&&!theintersect->endpt)
		{
			/*  calculate the distance to its preceding intersection  */

			/*obtain the coordinates of this intersection*/
			intersect_coords[0]=theintersect->gpos[0];
			intersect_coords[1]=theintersect->gpos[1];
			which_cell=majRoadnet->nodelist->intersects[intersect_id]->cellid;

			loc_dist.entry[0]=intersect_coords[0]-pt[0];
			loc_dist.entry[1]=intersect_coords[1]-pt[1];

			distToIntersect=length(loc_dist);

			/*  if the distance is really small, move the end point back to 
				the intersection
			*/
			if(distToIntersect<1.e-2)
			{
				/*  merge the end point with the intersection  
				*/
				if(!startorend)
				{
					curtraj->remove_front_nlines(nlines_needtoremove);
					/*and then add a new one*/
					curtraj->add_front_nlines(one_new, 1);
				}
				else
				{
					curtraj->remove_last_nlines(nlines_needtoremove);
					/*and then add a new one*/
					curtraj->add_last_nlines(one_new, 1);
				}
				return;
			}
		}

		/*  see whether it is besides the water region  */
		if(is_inland_pixel(pt[0],pt[1], xstart, xrang, ystart, yrang, dx, dy,
				fittedmap1, 512))
		{

			/*
				otherwise, we trace from this end point for certain distance
				until it hits other road, record the tracing distance and
				compare this distance with the distance to the intersection
				choose the smallest one
			*/

			/*  remember we set the tracing distance as trajDir  */
			Trajectory *temp_traj=new Trajectory(-1);

			if(trace_to_comp_nearby_traj(pt, cell, majormin, curplace->dsep, 
				100, MinDistCrossRiver, 2*quadmesh->xinterval, dist_did,  temp_traj))
			{
				/*  if the distance to the intersection is smaller than half of the 
					distance to the traveling
				*/
				if(distToIntersect<dist_did/2.)
				{
					/*  retreat to the intersection */
					if(!notuseIntersect&&!startorend)
					{
						curtraj->remove_front_nlines(nlines_needtoremove);
						/*and then add a new one*/
						curtraj->add_front_nlines(one_new, 1);
					}
					else if(!notuseIntersect && startorend)
					{
						curtraj->remove_last_nlines(nlines_needtoremove);
						/*and then add a new one*/
						curtraj->add_last_nlines(one_new, 1);
					}
				}
				else
				{
					if(!startorend) /*  start point  */
					{
						temp_traj->reverse_lines(); //reverse the obtain short lines
						curtraj->add_front_nlines(temp_traj->linesegs, temp_traj->nlinesegs);
					}
					else
						curtraj->add_last_nlines(temp_traj->linesegs, temp_traj->nlinesegs);

					if(one_new!=NULL)
						delete one_new;
				}
			}
			else if(!notuseIntersect)
			{
				/*  retreat to the intersection */
				if(!startorend)
				{
					curtraj->remove_front_nlines(nlines_needtoremove);
					/*and then add a new one*/
					curtraj->add_front_nlines(one_new, 1);
				}
				else
				{
					curtraj->remove_last_nlines(nlines_needtoremove);
					/*and then add a new one*/
					curtraj->add_last_nlines(one_new, 1);
				}
			}

			delete temp_traj;

		}

		else
		{
			if(sharedvars.AllowMajRoadCrossRiverOn)
			{
				/* make it cross river first and continue the tracing from there */
				double ce[2]={0.};
				icVector2 norm;
				norm=trajDir;
				normalize(norm);
				ce[0]=pt[0]+MinDistCrossRiver*quadmesh->xinterval*norm.entry[0];
				ce[1]=pt[1]+MinDistCrossRiver*quadmesh->xinterval*norm.entry[1];
				int ce_cell=get_cellID_givencoords(ce[0], ce[1]);

				if(is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
						fittedmap1, 512))
				{
					/*  we need to do a bindary search to find out the proper 
					    distance to cross the river
					*/
					double smaller_dist=MinDistCrossRiver/2.;
					ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
					ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
					ce_cell=get_cellID_givencoords(ce[0], ce[1]);
					if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
							fittedmap1, 512))
					{
						smaller_dist=.75*MinDistCrossRiver;
						ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
						ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
						ce_cell=get_cellID_givencoords(ce[0], ce[1]);
						if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
								fittedmap1, 512))
						{
							ce[0]=pt[0]+MinDistCrossRiver*quadmesh->xinterval*norm.entry[0];
							ce[1]=pt[1]+MinDistCrossRiver*quadmesh->xinterval*norm.entry[1];
							ce_cell=get_cellID_givencoords(ce[0], ce[1]);
						}
					}
					else
					{
						smaller_dist=MinDistCrossRiver/4.;
						ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
						ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
						ce_cell=get_cellID_givencoords(ce[0], ce[1]);
						if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
							fittedmap1, 512))
						{
							smaller_dist=MinDistCrossRiver*.375;
							ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
							ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
							ce_cell=get_cellID_givencoords(ce[0], ce[1]);
							if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
								fittedmap1, 512))
							{
								smaller_dist=MinDistCrossRiver/2.;
								ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
								ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
								ce_cell=get_cellID_givencoords(ce[0], ce[1]);
							}
						}
					}



					Trajectory *temp_traj=new Trajectory(-1);

					get_linesegs_anytwopts(pt, cell,   ce, ce_cell,   temp_traj, 0, 200);

					if(trace_to_comp_nearby_traj(ce, ce_cell, majormin, curplace->dsep, 
						100, MinDistCrossRiver, 2*quadmesh->xinterval, dist_did,  temp_traj))
					{
						/*  if the distance to the intersection is smaller than half of the 
							distance to the traveling
						*/
						if(!notuseIntersect&& distToIntersect<dist_did/2.)
						{
							/*  retreat to the intersection */
							if(!startorend)
							{
								curtraj->remove_front_nlines(nlines_needtoremove);
								/*and then add a new one*/
								curtraj->add_front_nlines(one_new, 1);
							}
							else
							{
								curtraj->remove_last_nlines(nlines_needtoremove);
								/*and then add a new one*/
								curtraj->add_last_nlines(one_new, 1);
							}
						}
						else
						{
							if(!startorend) /*  start point  */
							{
								temp_traj->reverse_lines(); //reverse the obtain short lines
								curtraj->add_front_nlines(temp_traj->linesegs, temp_traj->nlinesegs);
							}
							else
								curtraj->add_last_nlines(temp_traj->linesegs, temp_traj->nlinesegs);

							if(one_new!=NULL)
								delete one_new;
						}
					}
					else if(!notuseIntersect)
					{
						/*  retreat to the intersection */
						if(!startorend)
						{
							curtraj->remove_front_nlines(nlines_needtoremove);
							/*and then add a new one*/
							curtraj->add_front_nlines(one_new, 1);
						}
						else
						{
							curtraj->remove_last_nlines(nlines_needtoremove);
							/*and then add a new one*/
							curtraj->add_last_nlines(one_new, 1);
						}
					}

					delete temp_traj;
				}
			}
			else if(!notuseIntersect)
			{
				/*   retreat to the preceding intersection   */
			}
		}
	}
}


void connect_majRoads_postproc()
{
	if(major_level1==NULL || minor_level1 == NULL 
		|| majRoadintersectinfo_maj == NULL || majRoadintersectinfo_min == NULL)
		return;

	int i;
	/*   we first process the major roads along major direction   */
	/*consider major tensor lines*/
	for(i=0; i<major_level1->evenstreamlines->ntrajs; i++)
	{
		/*  we leave out the highway  12/24/2007  */
		connect_one_majRoad(i, false, false);
		connect_one_majRoad(i, false, true);
	}

	FILE *fp;
	//fp=fopen("remove_deadend_error.txt", "w");
	//fprintf(fp, "finish removing dead end from major roads\n");
	//fclose(fp);

	major_level1->rebuild_all_lineinfo();

	/*   we then rebuild the tensor line information and rebuild the graph   */
	init_majRoadnet();
	compute_majRoad_intersects();
	search_for_connection_majRoad();

	/*   second, we process the major roads along minor direction   */
	for(i=0; i<minor_level1->evenstreamlines->ntrajs; i++)
	{
		/*  we leave out the highway  12/24/2007  */
		connect_one_majRoad(i, true, false);
		connect_one_majRoad(i, true, true);
	}

	init_majRoadnet();
	compute_majRoad_intersects();
	search_for_connection_majRoad();
	

	/*   Re-consider the major direction again since minor direction has been changed  */
	for(i=0; i<major_level1->evenstreamlines->ntrajs; i++)
	{
		/*  we leave out the highway  12/24/2007  */
		connect_one_majRoad(i, false, false);
		connect_one_majRoad(i, false, true);
	}
	init_majRoadnet();
	compute_majRoad_intersects();
	search_for_connection_majRoad();
}


void connect_majRoads_graph(bool majormin, bool startorend)
{
	icVector2 dir;
	int i, j;
	double A[2], ang, dist_intersect, dist_edge;
	Trajectory *traj;
	int cell_id;
	int which_intersect, which_edge, which_samp, which_cell;
	double the_intersect[2];

	int cur_intersect_id;
	EvenStreamlinePlace *cur_place;
	TensorLineIntersectionInfoList **infolist;

	int except_edge;

	if(!majormin)
	{
		cur_place=major_level1;
		infolist=majRoadintersectinfo_maj;
	}
	else{
		cur_place=minor_level1;
		infolist=majRoadintersectinfo_min;
	}

	for(i=cur_place->evenstreamlines->ntrajs-cur_place->nnewlines_inReg;
		i<cur_place->evenstreamlines->ntrajs;i++)
	{
		traj=cur_place->evenstreamlines->trajs[i];

		dist_intersect=dist_edge=1.e50;
		bool closetointersect=false;
		bool closetoedge=false;
		double temp_dist=0;

		/*  get the start point  */

		if(!startorend)
		{
			cur_intersect_id = infolist[i]->infolist[0]->intersect_id;

			A[0]=traj->linesegs[0].gstart[0];
			A[1]=traj->linesegs[0].gstart[1];

			/*  if it is already on the domain boundary, do nothing */
			if(A[0]<=quadmesh->xstart+1.e-8 || A[0]>=quadmesh->xend-1.e-8
				|| A[1]<=quadmesh->ystart+1.e-8 || A[1]>=quadmesh->yend-1.e-8)
				return;

			//if(traj->nlinesegs>2)
			//{
			//	dir.entry[0]=traj->linesegs[0].gstart[0]-traj->linesegs[2].gend[0];
			//	dir.entry[1]=traj->linesegs[0].gstart[1]-traj->linesegs[2].gend[1];
			//}
			//else{
			//	dir.entry[0]=traj->linesegs[0].gstart[0]-traj->linesegs[1].gend[0];
			//	dir.entry[1]=traj->linesegs[0].gstart[1]-traj->linesegs[1].gend[1];
			//}
			dir.entry[0]=traj->linesegs[0].gstart[0]-traj->linesegs[0].gend[0];
			dir.entry[1]=traj->linesegs[0].gstart[1]-traj->linesegs[0].gend[1];
			normalize(dir);
			cell_id=traj->linesegs[0].Triangle_ID;
		}
		else
		{
			cur_intersect_id = infolist[i]->infolist[infolist[i]->nelems-1]->intersect_id;

			A[0]=traj->linesegs[traj->nlinesegs-1].gend[0];
			A[1]=traj->linesegs[traj->nlinesegs-1].gend[1];

			/*  if it is already on the domain boundary, do nothing */
			if(A[0]<=quadmesh->xstart+1.e-8 || A[0]>=quadmesh->xend-1.e-8
				|| A[1]<=quadmesh->ystart+1.e-8 || A[1]>=quadmesh->yend-1.e-8)
				return;

			//if(traj->nlinesegs>2)
			//{
			//	dir.entry[0]=traj->linesegs[traj->nlinesegs-1].gend[0]
			//		-traj->linesegs[traj->nlinesegs-3].gstart[0];
			//	dir.entry[1]=traj->linesegs[traj->nlinesegs-1].gend[1]
			//		-traj->linesegs[traj->nlinesegs-3].gstart[1];
			//}
			//else
			//{
			//	dir.entry[0]=traj->linesegs[traj->nlinesegs-1].gend[0]
			//		-traj->linesegs[traj->nlinesegs-2].gstart[0];
			//	dir.entry[1]=traj->linesegs[traj->nlinesegs-1].gend[1]
			//		-traj->linesegs[traj->nlinesegs-2].gstart[1];
			//}
			dir.entry[0]=traj->linesegs[traj->nlinesegs-1].gend[0]
				-traj->linesegs[traj->nlinesegs-1].gstart[0];
			dir.entry[1]=traj->linesegs[traj->nlinesegs-1].gend[1]
				-traj->linesegs[traj->nlinesegs-1].gstart[1];
			normalize(dir);
			cell_id=traj->linesegs[traj->nlinesegs-2].Triangle_ID;
		}

		/*  obtain the edge that is associated with the intersection  */
		except_edge=streetnet->nodelist->intersects[cur_intersect_id]->adj_edges[0];

		ang=atan2(dir.entry[1], dir.entry[0]);
		if(ang<0) ang+=(2.*M_PI);

		/*  compute the distance between A and other nearby intersections  
		    The following routines still use streetnet inside !
		*/
		//if(cal_smallest_dist_to_intersects(cell_id, A, ang, dir, dist_intersect, which_intersect,
		//	cur_intersect_id))
		{
			//dist_intersect=temp_dist;
			if(which_intersect!=cur_intersect_id)
			closetointersect=true;

			/*  we may want to record the distance to this intersection */
		}

		//int result=trace_for_intersection(A, dir, which_intersect, which_edge, which_samp, 
		//	the_intersect, which_cell, except_edge, cur_intersect_id, 5);
	}
}
/************************************************************/

bool get_a_proper_seed(EvenStreamlinePlace *curplace, double samp[2], int &cell, int fieldtype)
{
	int i;

	while(1)
	{
		for(i = curplace->cur_seed_pos; i < curplace->seedpts->nseeds; i++)
		{
			////Judge whether this is a valid seed or not
			/*meaning that it is not close to any exiting tensor lines*/
			if(!curplace->close_to_cur_streamline(curplace->seedpts->seeds[i]->pos,
				curplace->seedpts->seeds[i]->triangle, &curplace->evenstreamlines->ntrajs, 1, 
				curplace->minstartdist, curplace->discsize, 0))
			{
				////if we find a valid seed point, break
				curplace->cur_seed_pos = i+1;
				samp[0]=curplace->seedpts->seeds[i]->pos[0];
				samp[1]=curplace->seedpts->seeds[i]->pos[1];
				cell=curplace->seedpts->seeds[i]->triangle;

				return true;
			}
	
			curplace->seedpts->seeds[i]->state = 2;  //reject before starting
		}

		if(i >= curplace->seedpts->nseeds ) // it means no more seeds available for current streamline
		{
			curplace->cur_seed_pos = curplace->seedpts->nseeds;
			curplace->cur_traj ++;  //let next streamline as current streamline

			if(curplace->cur_traj >= curplace->evenstreamlines->ntrajs) //no more streamlines available
			{
				/*release the seed point list*/
				return false;
			}
			
			//Get seeds for the current streamline
			curplace->cal_seeds(curplace->cur_traj, curplace->seeddist, 
				curplace->every_nsample, fieldtype, brushinterfaceOn);
			continue;
		}
	}
}


/*
Judge whether an element has already been stored in an integer array
This can be extend to a template!
*/
bool is_repeated_elem(int *a, int b, int num)
{
	int i;
	for(i = 0; i < num; i++)
	{
		if(a[i] == b)
			return true;
	}
	return false;
}


/*
We use the following routine to get rid of all the "dead ends"
This routine will be called *after* the placement of both major and minor tensor lines
The basic idea is: we move the two ends of the each tensor line along the 
directions determined by the two end line segments, and try to find out the 
closest tensor line (not the same time as current tensor line).
And we compute the intersection and connect them (could move forward or backward)

NOTE: we also need to rebuild the majlineinfo and minlineinfo after the 
reconfiguration of the obtained tensor lines
11/18/2007

*/
extern LineSeg **trajectories;                 //trajectories' list
extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory
extern int MaxNumTrajectories;
extern int MaxNumLinesegsPerTraj;

void remove_all_deadends()
{
	int i;
	if(major==NULL || minor == NULL || majorintersectinfo == NULL || minorintersectinfo == NULL)
		return;

	/*                                                          */
	/*     we need to copy the major tensor lines 11/25/2007    */
	
	if(t_major != NULL)
		delete t_major;
		
	t_major=new TrajectoryList(major->evenstreamlines->ntrajs);

	for(i=0;i<t_major->curMaxNumTrajs;i++)
	{
		if(t_major->trajs[i]==NULL)
			t_major->trajs[i]=new Trajectory(i, major->evenstreamlines->trajs[i]->curMaxNumLinesegs);

		t_major->trajs[i]->nlinesegs=0;
	}

	for(i=0;i<t_major->ntrajs;i++)
	{
		t_major->trajs[i]->add_front_nlines(major->evenstreamlines->trajs[i]->linesegs,
			major->evenstreamlines->trajs[i]->nlinesegs);
	}

	/*consider major tensor lines*/
	for(i=0; i<major->evenstreamlines->ntrajs; i++)
	{
		/*  we leave out the highway  12/24/2007  */
		if(major->evenstreamlines->trajs[i]->roadtype==HIGHWAY
			||major->evenstreamlines->trajs[i]->roadtype==MAJOR)
		{
			if(!sharedvars.RemoveMajDeadEndsOn)
				continue;
		}

		if(!sharedvars.UseNewMedRemDeadEndOn)
			adj_deadend_onetensorline(0, i, major);
		else
		{
			remove_one_deadend_minRoad(i, false, false, major->dsep);
			remove_one_deadend_minRoad(i, false, true, major->dsep);
		}
	}

	FILE *fp;
	//fp=fopen("remove_deadend_error.txt", "w");
	//fprintf(fp, "finish removing dead end from major roads\n");
	//fclose(fp);

	major->rebuild_all_lineinfo();
	
	/*rebuild the network, since some information has changed*/
	reset_roadnetvis();
	compute_intersects();
	search_for_connection();

	/*                                                          */
	/*     we need to copy the major tensor lines 11/25/2007    */

	/*     since we are not dealing with minor and minor intersection
	       we don't need to duplicate the information for minor lines?
	*/
	
	//fp=fopen("remove_deadend_error.txt", "w");
	//fprintf(fp, "start removing dead end from minor roads\n");
	//fclose(fp);

	/*consider minor tensor lines*/
	for(i=0; i<minor->evenstreamlines->ntrajs; i++)
	{
		/*  we leave out the highway  12/24/2007  */
		if(minor->evenstreamlines->trajs[i]->roadtype==HIGHWAY
			||minor->evenstreamlines->trajs[i]->roadtype==MAJOR)
		{
			if(!sharedvars.RemoveMajDeadEndsOn)
				continue;
		}

		if(!sharedvars.UseNewMedRemDeadEndOn)
			adj_deadend_onetensorline(1, i, minor);
		else
		{
			remove_one_deadend_minRoad(i, true, false, minor->dsep);
			remove_one_deadend_minRoad(i, true, true, minor->dsep);
		}
	}
	
	//fp=fopen("remove_deadend_error.txt", "w");
	//fprintf(fp, "start reconstructing minor roads info\n");
	//fclose(fp);

	/*we also need to reconstruct the major and minor line information list
	to obtain the proper intersection calculation*/
	minor->rebuild_all_lineinfo();


	/*  release the memeory  */
	delete t_major;
	t_major=NULL;
	//delete t_minor;
	//t_minor=NULL;
}

/*
choose a proper end :)
11/18/2007
*/
void adj_deadend_onetensorline(int fieldtype, int id, EvenStreamlinePlace *curplace)
{
	/*we need not change the closed tensor line*/
	if(curplace->evenstreamlines->trajs[id]->closed)
		return;

	EvenStreamlinePlace *othertypeplace;
	if(fieldtype==0)
		othertypeplace=minor;
	else
		othertypeplace=major;

	/*                                                    */
	/*              process the start position            */

	icVector2 line_dir, line_dir2;
	Trajectory *curtraj=curplace->evenstreamlines->trajs[id];
	LineSeg *theline=&curplace->evenstreamlines->trajs[id]->linesegs[0];
	int nlines=curplace->evenstreamlines->trajs[id]->nlinesegs;
	line_dir.entry[0]=theline->gend[0]-theline->gstart[0];
	line_dir.entry[1]=theline->gend[1]-theline->gstart[1];


	line_dir=-line_dir;

	int cell_id=theline->Triangle_ID;
	double start[2]={theline->gstart[0], theline->gstart[1]};
	bool firstdirok=false;
	double intersect1[2], intersect2[2];
	double dist_front_forward=1.e50;
	LineSeg *temp_linesegs=NULL;
	int temp_lines;

	bool noavailable_intersect_start=false;
	bool noavailable_intersect_end=false;
	double dist_front_backward=1.e50;
	int intersect_id=0;
	double intersect_coords[2]={0.};
	int which_cell=0;

	//if(theline->gstart[0]>=quadmesh->xend-1.e-7 || theline->gstart[0]<=quadmesh->xstart+1.e-7
	//	||theline->gstart[1]>=quadmesh->yend-1.e-7 || theline->gstart[1]<=quadmesh->ystart+1.e-7)
	//	goto LEND;

	FILE *fp;

	if(find_cell_contain_bothmajmin(start, cell_id, curplace->dsep,
		line_dir, fieldtype))
	{
		/*compute the intersection*/
		compute_extended_intersect(start, cell_id, intersect1, 
			othertypeplace, fieldtype);
		line_dir2.entry[0]=intersect1[0]-theline->gstart[0];
		line_dir2.entry[1]=intersect1[1]-theline->gstart[1];
		dist_front_forward=length(line_dir2);
		if(dot(line_dir2, line_dir)<0)
			firstdirok=false;

		else{
			firstdirok=true;

			/*Then, we rearrange the line segment list
			NOTE: if we change the arrangement of these line segments,
			the following computation will be wrong
			We can simply record them and modify the trajectory later
			*/
			int i;
			temp_linesegs=(LineSeg*)malloc(sizeof(LineSeg)*(num_linesegs_curtraj[cur_traj_index]+2));
			temp_lines=num_linesegs_curtraj[cur_traj_index];
			for(i=0; i<temp_lines; i++)
			{
				temp_linesegs[i].gstart[0]=trajectories[cur_traj_index][temp_lines-1-i].gend[0];
				temp_linesegs[i].gstart[1]=trajectories[cur_traj_index][temp_lines-1-i].gend[1];
				temp_linesegs[i].gend[0]=trajectories[cur_traj_index][temp_lines-1-i].gstart[0];
				temp_linesegs[i].gend[1]=trajectories[cur_traj_index][temp_lines-1-i].gstart[1];
				
				temp_linesegs[i].start[0]=trajectories[cur_traj_index][temp_lines-1-i].end[0];
				temp_linesegs[i].start[1]=trajectories[cur_traj_index][temp_lines-1-i].end[1];
				temp_linesegs[i].end[0]=trajectories[cur_traj_index][temp_lines-1-i].start[0];
				temp_linesegs[i].end[1]=trajectories[cur_traj_index][temp_lines-1-i].start[1];

				temp_linesegs[i].length=trajectories[cur_traj_index][temp_lines-1-i].length;
				temp_linesegs[i].Triangle_ID=trajectories[cur_traj_index][temp_lines-1-i].Triangle_ID;
			}

			/*we need to add one more line segment here !*/
			temp_linesegs[temp_lines].gstart[0]=intersect1[0];
			temp_linesegs[temp_lines].gstart[1]=intersect1[1];
			temp_linesegs[temp_lines].gend[0]=start[0];
			temp_linesegs[temp_lines].gend[1]=start[1];
			temp_linesegs[temp_lines].Triangle_ID=cell_id;
			line_dir.entry[0]=intersect1[0]-start[0];
			line_dir.entry[1]=intersect1[1]-start[1];
			temp_linesegs[i].length=length(line_dir);
			temp_lines++;
		}
	}
	

	/*search the other direction along the tensor line to find the closest
	intersection to the start point
	We may need to use the majorintersectinfo/minorintersectinfo*/

	/**/
	if(fieldtype==0) /*it is major tensor line*/
	{
		intersect_id=majorintersectinfo[id]->infolist[0]->intersect_id;
		if(streetnet->nodelist->intersects[intersect_id]->endpt)
		{
			if(majorintersectinfo[id]->nelems>1)
				intersect_id=majorintersectinfo[id]->infolist[1]->intersect_id;
			else
				noavailable_intersect_start=true;
		}
		else
			return;  // we don't need to extend it! 1/10/2008
	}
	else  /* it is minor direction */
	{
	
		intersect_id=minorintersectinfo[id]->infolist[0]->intersect_id;

		if(streetnet->nodelist->intersects[intersect_id]->endpt)
		{
			if(minorintersectinfo[id]->nelems>1)
				intersect_id=minorintersectinfo[id]->infolist[1]->intersect_id;
			else
				noavailable_intersect_start=true;
		}
		else
			return;  // we don't need to extend it!  1/10/2008

	}
	

	/*   may contain bug here 1/8/2008 */
	if(fieldtype==0)
	{
		if(majorintersectinfo[id]->nelems<=1) /*  a bug here need to be fixed  */
			return;
	}
	else
	{
		if(minorintersectinfo[id]->nelems<=1) /*  a bug here need to be fixed  */
			return;
	}

	/*  it will contain some issues (bugs) after the new intersection computation  */

	/*obtain the coordinates of this intersection*/
	intersect_coords[0]=streetnet->nodelist->intersects[intersect_id]->gpos[0];
	intersect_coords[1]=streetnet->nodelist->intersects[intersect_id]->gpos[1];
	which_cell=streetnet->nodelist->intersects[intersect_id]->cellid;
	int del_front_till;

	if(fieldtype==0)
	{
		if(streetnet->nodelist->intersects[intersect_id]->intersect_type==0)
			del_front_till = streetnet->nodelist->intersects[intersect_id]->majlineseg;
		else
		{
			if(streetnet->nodelist->intersects[intersect_id]->majorline_id==id)
				del_front_till = streetnet->nodelist->intersects[intersect_id]->majlineseg;
			else if(streetnet->nodelist->intersects[intersect_id]->minorline_id==id)
				del_front_till = streetnet->nodelist->intersects[intersect_id]->minlineseg;
			else
			{
				int test=0;
			}
		}
	}
	else{
		if(streetnet->nodelist->intersects[intersect_id]->intersect_type==0)
			del_front_till = streetnet->nodelist->intersects[intersect_id]->minlineseg;
		else
		{
			if(streetnet->nodelist->intersects[intersect_id]->majorline_id==id)
				del_front_till = streetnet->nodelist->intersects[intersect_id]->majlineseg;
			else if(streetnet->nodelist->intersects[intersect_id]->minorline_id==id)
				del_front_till = streetnet->nodelist->intersects[intersect_id]->minlineseg;
			else
			{
				int test=0;
			}
		}
	}
	

	/*create a new line segment*/
	LineSeg *front_new=(LineSeg*)malloc(sizeof(LineSeg));

	//if(del_front_till 
	front_new->gstart[0]=intersect_coords[0];
	front_new->gstart[1]=intersect_coords[1];
	front_new->start[0]=intersect_coords[0];
	front_new->start[1]=intersect_coords[1];
	front_new->gend[0]=curtraj->linesegs[del_front_till].gend[0];
	front_new->gend[1]=curtraj->linesegs[del_front_till].gend[1];
	front_new->end[0]=curtraj->linesegs[del_front_till].gend[0];
	front_new->end[1]=curtraj->linesegs[del_front_till].gend[1];
	front_new->Triangle_ID=curtraj->linesegs[del_front_till].Triangle_ID;
	line_dir.entry[0]=front_new->gend[0]-front_new->gstart[0];
	line_dir.entry[1]=front_new->gend[1]-front_new->gstart[1];
	front_new->length=length(line_dir);

	/*we then choose the intersection that closest to the current end point*/
	line_dir.entry[0]=intersect_coords[0]-theline->gstart[0];
	line_dir.entry[1]=intersect_coords[1]-theline->gstart[1];

	dist_front_backward=length(line_dir);

	/////////////////////////////////////////////////////////////
	/*                                                    */
	/*                 process the end point              */
LEND:

	theline=&curplace->evenstreamlines->trajs[id]->linesegs[nlines-1];
	line_dir.entry[0]=theline->gend[0]-theline->gstart[0];
	line_dir.entry[1]=theline->gend[1]-theline->gstart[1];
	
	//if(theline->gend[0]>=quadmesh->xend-1.e-7 || theline->gend[0]<=quadmesh->xstart+1.e-7
	//	||theline->gend[1]>=quadmesh->yend-1.e-7 || theline->gend[1]<=quadmesh->ystart+1.e-7)
	//{
	//	if(temp_linesegs!=NULL)
	//	{
	//		free(temp_linesegs);
	//	}
	//	if(front_new!=NULL)
	//		free(front_new);
	//	//free(back_new);
	//	return;
	//}

	cell_id=theline->Triangle_ID;
	start[0]=theline->gend[0];
	start[1]=theline->gend[1];
	bool enddirok=false;
	double dist_back_forward=1.e50;

	if(find_cell_contain_bothmajmin(start, cell_id, curplace->dsep,
		line_dir, fieldtype))
	{
		/*compute the intersection*/
		compute_extended_intersect(start, cell_id, intersect1, 
			othertypeplace, fieldtype);
		line_dir2.entry[0]=intersect1[0]-theline->gstart[0];
		line_dir2.entry[1]=intersect1[1]-theline->gstart[1];
		dist_back_forward=length(line_dir2);
		if(dot(line_dir2, line_dir)<0) enddirok=false;

		else{
			enddirok=true;

			/*we need one more line segment here !*/
			trajectories[cur_traj_index][num_linesegs_curtraj[cur_traj_index]].gend[0]=intersect1[0];
			trajectories[cur_traj_index][num_linesegs_curtraj[cur_traj_index]].gend[1]=intersect1[1];
			trajectories[cur_traj_index][num_linesegs_curtraj[cur_traj_index]].gstart[0]=start[0];
			trajectories[cur_traj_index][num_linesegs_curtraj[cur_traj_index]].gstart[1]=start[1];
			trajectories[cur_traj_index][num_linesegs_curtraj[cur_traj_index]].Triangle_ID=cell_id;
			line_dir.entry[0]=intersect1[0]-start[0];
			line_dir.entry[1]=intersect1[1]-start[1];
			trajectories[cur_traj_index][num_linesegs_curtraj[cur_traj_index]].length=length(line_dir);
			num_linesegs_curtraj[cur_traj_index]++;
		}
	}
	

	/*search the other direction along the tensor line to find the closest
	intersection to the start point
	We may need to use the majorintersectinfo/minorintersectinfo*/

	/**/
	if(fieldtype==0) /*it is major tensor line*/
	{
		intersect_id=majorintersectinfo[id]->infolist[majorintersectinfo[id]->nelems-1]->intersect_id;
		if(streetnet->nodelist->intersects[intersect_id]->endpt)
		{
			if(majorintersectinfo[id]->nelems>1)
				intersect_id=majorintersectinfo[id]->infolist[majorintersectinfo[id]->nelems-2]->intersect_id;
			else
				noavailable_intersect_end=true;
		}
		else
			return;  // we don't need to extend it!  1/10/2008
	}
	else
	{
		intersect_id=minorintersectinfo[id]->infolist[minorintersectinfo[id]->nelems-1]->intersect_id;
		if(streetnet->nodelist->intersects[intersect_id]->endpt)
		{
			if(minorintersectinfo[id]->nelems>1)
				intersect_id=minorintersectinfo[id]->infolist[minorintersectinfo[id]->nelems-2]->intersect_id;
			else
				noavailable_intersect_end=true;
		}
		else 
			return;  // we don't need to extend it!  1/10/2008
	}

	/*obtain the coordinates of this intersection*/
	intersect_coords[0]=streetnet->nodelist->intersects[intersect_id]->gpos[0];
	intersect_coords[1]=streetnet->nodelist->intersects[intersect_id]->gpos[1];
	which_cell=streetnet->nodelist->intersects[intersect_id]->cellid;
	int del_back_till;

	if(fieldtype==0)
	{
		if(streetnet->nodelist->intersects[intersect_id]->intersect_type==0)
			del_back_till = streetnet->nodelist->intersects[intersect_id]->majlineseg;
		else
		{
			if(streetnet->nodelist->intersects[intersect_id]->majorline_id==id)
				del_back_till = streetnet->nodelist->intersects[intersect_id]->majlineseg;
			else if(streetnet->nodelist->intersects[intersect_id]->minorline_id==id)
				del_back_till = streetnet->nodelist->intersects[intersect_id]->minlineseg;
			else
			{
				int test=0;
			}
		}
	}
	else
	{
		if(streetnet->nodelist->intersects[intersect_id]->intersect_type==0)
			del_back_till = streetnet->nodelist->intersects[intersect_id]->minlineseg;
		else
		{
			if(streetnet->nodelist->intersects[intersect_id]->majorline_id==id)
				del_back_till = streetnet->nodelist->intersects[intersect_id]->majlineseg;
			else if(streetnet->nodelist->intersects[intersect_id]->minorline_id==id)
				del_back_till = streetnet->nodelist->intersects[intersect_id]->minlineseg;
			else
			{
				int test=0;
			}
		}
	}

	int nlines_backneedtoremove=curtraj->nlinesegs-del_back_till;
	
	/*create a new line segment*/
	LineSeg *back_new=(LineSeg*)malloc(sizeof(LineSeg));
	back_new->gstart[0]=curtraj->linesegs[del_back_till].gstart[0];
	back_new->gstart[1]=curtraj->linesegs[del_back_till].gstart[1];
	back_new->start[0]=curtraj->linesegs[del_back_till].gstart[0];
	back_new->start[1]=curtraj->linesegs[del_back_till].gstart[1];
	back_new->gend[0]=intersect_coords[0];
	back_new->gend[1]=intersect_coords[1];
	back_new->end[0]=intersect_coords[0];
	back_new->end[1]=intersect_coords[1];
	back_new->Triangle_ID=curtraj->linesegs[del_back_till].Triangle_ID;
	line_dir.entry[0]=back_new->gend[0]-back_new->gstart[0];
	line_dir.entry[1]=back_new->gend[1]-back_new->gstart[1];
	back_new->length=length(line_dir);

	/*we then choose the intersection that closest to the current end point*/
	line_dir.entry[0]=intersect_coords[0]-theline->gstart[0];
	line_dir.entry[1]=intersect_coords[1]-theline->gstart[1];

	double dist_back_backward=length(line_dir);
	

	/***************************************************************************/
	/*we rearrange the whole line segment list including the head and the tail*/

	/*re-arrange the first n line segments*/
	if((dist_front_forward<dist_front_backward && firstdirok) \
		|| noavailable_intersect_start && firstdirok /*&& fieldtype==0*/)
	{
		/*we add new line segment(s) in the front*/
		curtraj->add_front_nlines(temp_linesegs, temp_lines);
	}
	else if(!noavailable_intersect_start)
	{
		/*we remove front line n segment(s)*/
		curtraj->remove_front_nlines(del_front_till);

		/*and then add a new one*/
		curtraj->add_front_nlines(front_new, 1);
	}

	/*re-arrange the last n line segments*/
	if((dist_back_forward<dist_back_backward && enddirok)
		|| noavailable_intersect_end && enddirok /*&& fieldtype==0*/)
	{
		/*we add new line segment(s) in the front*/
		curtraj->add_last_nlines(trajectories[cur_traj_index],num_linesegs_curtraj[cur_traj_index]);
	}
	else if(!noavailable_intersect_end)
	{
		/*we remove back n line segment(s)*/
		curtraj->remove_last_nlines(nlines_backneedtoremove);

		/*and then add a new one*/
		curtraj->add_last_nlines(back_new, 1);
	}

	/*release pre-allocated variables*/
	if(temp_linesegs!=NULL)
	{
		free(temp_linesegs);
	}
	free(front_new);
	free(back_new);
}


/*
    Remove one dead end given the type of the trajectory and its id
	and whether it is start point or end point

	could contain bug
*/

void remove_one_deadend_minRoad(int id, bool majormin, bool startorend,
								double search_dist)
{
	int i;

	EvenStreamlinePlace *curplace;
	Trajectory *curtraj;
	TensorLineIntersectionInfoList *theinfolist;
	Intersection *theintersect;
	LineSeg *curline;
	double pt[2];
	int cell;
	icVector2 trajDir;
	double otherpt[2];
	double distToIntersect, dist_did;

	///*  variables for pixel level judgement of "inland" or not  */
	//double xstart=quadmesh->xstart;
	//double ystart=quadmesh->ystart;
	//double xrang=quadmesh->xend-quadmesh->xstart;
	//double yrang=quadmesh->yend-quadmesh->ystart;
	//double dx=xrang/511;
	//double dy=yrang/511;

	icVector2 loc_dist;

	double intersect_coords[2];
	int which_cell;
	int intersect_id;

	int del_till;

	bool notuseIntersect=false;

	hstep = quadmesh->xinterval/2.;
	predict_stepsize = quadmesh->xinterval/2.;
	euler_stepsize = quadmesh->xinterval/5.;

	/*create a new line segment*/
	icVector2 line_dir;


	/*    intialization    */

	if(!majormin)  /*  major direction  */
	{
		curplace=major;
		curtraj=major->evenstreamlines->trajs[id];
		theinfolist=majorintersectinfo[id];
	}
	else           /*  minor direction  */
	{
		curplace=minor;
		curtraj=minor->evenstreamlines->trajs[id];
		theinfolist=minorintersectinfo[id];
	}

	/*  judge whether it is a dead end or not  */
	if(!startorend) /* start point  */
	{
		if(!streetnet->nodelist->intersects[theinfolist->infolist[0]->intersect_id]
			->endpt)
			return;
	}
	else
	{
		if(!streetnet->nodelist->intersects[theinfolist->infolist[theinfolist->nelems-1]->intersect_id]
			->endpt)
			return;
	}

	if(!startorend)/*  start point  */
	{
		curline=&curtraj->linesegs[0];
		//trajDir.entry[0]=curline->gstart[0]-curline->gend[0];
		//trajDir.entry[1]=curline->gstart[1]-curline->gend[1];
		int later_line=min(curtraj->nlinesegs-1, 1);
		trajDir.entry[0]=curline->gstart[0]-curtraj->linesegs[later_line].gend[0];
		trajDir.entry[1]=curline->gstart[1]-curtraj->linesegs[later_line].gend[1];
		
		pt[0]=curtraj->linesegs[0].gstart[0];
		pt[1]=curtraj->linesegs[0].gstart[1];
		
		if(theinfolist->nelems==1) /*  a bug here need to be fixed  */
			notuseIntersect=true;
		else
			intersect_id=theinfolist->infolist[1]->intersect_id; //avoid itself
	}
	else           /*  end point  */
	{
		curline=&curtraj->linesegs[curtraj->nlinesegs-1];
		//trajDir.entry[0]=curline->gend[0]-curline->gstart[0];
		//trajDir.entry[1]=curline->gend[1]-curline->gstart[1];
		int earlier_line=max(0, curtraj->nlinesegs-2);
		trajDir.entry[0]=curline->gend[0]-curtraj->linesegs[earlier_line].gstart[0];
		trajDir.entry[1]=curline->gend[1]-curtraj->linesegs[earlier_line].gstart[1];

		pt[0]=curtraj->linesegs[curtraj->nlinesegs-1].gend[0];
		pt[1]=curtraj->linesegs[curtraj->nlinesegs-1].gend[1];

		if(theinfolist->nelems==1) /*  a bug here need to be fixed  */
			notuseIntersect=true;
		else
			intersect_id=theinfolist->infolist[theinfolist->nelems-2]->intersect_id; //avoid itself
	}

	if(!notuseIntersect)
		theintersect=streetnet->nodelist->intersects[intersect_id];

	cell=get_cellID_givencoords(pt[0], pt[1]);

	tenline_dir_global=trajDir;
	normalize(tenline_dir_global);

	distToIntersect=1.e50;
	int nlines_needtoremove;

	if(!majormin)  /*  major direction  */
	{
		if(streetnet->nodelist->intersects[intersect_id]->intersect_type==1)
		{
			if(streetnet->nodelist->intersects[intersect_id]->majorline_id==id)
				del_till=streetnet->nodelist->intersects[intersect_id]->majlineseg;
			else
				del_till=streetnet->nodelist->intersects[intersect_id]->minlineseg;
		}
		else
			del_till = streetnet->nodelist->intersects[intersect_id]->majlineseg;
	
	}
	else           /*  minor direction  */
	{
		if(streetnet->nodelist->intersects[intersect_id]->intersect_type==2)
		{
			if(streetnet->nodelist->intersects[intersect_id]->majorline_id==id)
				del_till=streetnet->nodelist->intersects[intersect_id]->majlineseg;
			else
				del_till=streetnet->nodelist->intersects[intersect_id]->minlineseg;
		}
		else
		    del_till = streetnet->nodelist->intersects[intersect_id]->minlineseg;
	}

	LineSeg *one_new=NULL;
		
	/*  create new line segment for the back track to the preceding intersection  */
	if(!startorend && !notuseIntersect)  /*  start point */
	{
		nlines_needtoremove=del_till;

		one_new=(LineSeg*)malloc(sizeof(LineSeg));
		one_new->gstart[0]=theintersect->gpos[0];
		one_new->gstart[1]=theintersect->gpos[1];
		one_new->start[0]=theintersect->gpos[0];
		one_new->start[1]=theintersect->gpos[1];
		one_new->gend[0]=curtraj->linesegs[del_till].gend[0];
		one_new->gend[1]=curtraj->linesegs[del_till].gend[1];
		one_new->end[0]=curtraj->linesegs[del_till].gend[0];
		one_new->end[1]=curtraj->linesegs[del_till].gend[1];
		one_new->Triangle_ID=curtraj->linesegs[del_till].Triangle_ID;
		line_dir.entry[0]=one_new->gend[0]-one_new->gstart[0];
		line_dir.entry[1]=one_new->gend[1]-one_new->gstart[1];
		one_new->length=length(line_dir);
	}
	else if(startorend && !notuseIntersect)
	{
		nlines_needtoremove=curtraj->nlinesegs-del_till;

		one_new=(LineSeg*)malloc(sizeof(LineSeg));
		one_new->gstart[0]=curtraj->linesegs[del_till].gstart[0];
		one_new->gstart[1]=curtraj->linesegs[del_till].gstart[1];
		one_new->start[0]=curtraj->linesegs[del_till].gstart[0];
		one_new->start[1]=curtraj->linesegs[del_till].gstart[1];
		one_new->gend[0]=theintersect->gpos[0];
		one_new->gend[1]=theintersect->gpos[1];
		one_new->end[0]=theintersect->gpos[0];
		one_new->end[1]=theintersect->gpos[1];
		one_new->Triangle_ID=curtraj->linesegs[del_till].Triangle_ID;
		line_dir.entry[0]=one_new->gend[0]-one_new->gstart[0];
		line_dir.entry[1]=one_new->gend[1]-one_new->gstart[1];
		one_new->length=length(line_dir);
	}

	if(one_new->Triangle_ID<0||one_new->Triangle_ID>=quadmesh->nfaces)
	{
		int test=0;
	}

	if(pt[0]<=quadmesh->xstart+1.e-8||pt[0]>=quadmesh->xend-1.e-8
		||pt[1]<=quadmesh->ystart+1.e-8||pt[1]>=quadmesh->yend-1.e-8)
	{
		/* do nothing */
	}

	else{
		if(!notuseIntersect&&!theintersect->endpt)
		{
			/*  calculate the distance to its preceding intersection  */

			/*obtain the coordinates of this intersection*/
			intersect_coords[0]=theintersect->gpos[0];
			intersect_coords[1]=theintersect->gpos[1];
			which_cell=streetnet->nodelist->intersects[intersect_id]->cellid;

			loc_dist.entry[0]=intersect_coords[0]-pt[0];
			loc_dist.entry[1]=intersect_coords[1]-pt[1];

			distToIntersect=length(loc_dist);

			/*  if the distance is really small, move the end point back to 
				the intersection
			*/
			if(distToIntersect<1.e-2)
			{
				/*  merge the end point with the intersection  
				*/
				if(!startorend)
				{
					curtraj->remove_front_nlines(nlines_needtoremove);
					/*and then add a new one*/
					curtraj->add_front_nlines(one_new, 1);
				}
				else
				{
					curtraj->remove_last_nlines(nlines_needtoremove);
					/*and then add a new one*/
					curtraj->add_last_nlines(one_new, 1);
				}
				return;
			}
			else
			{
				/*  Let's trace from the end point for certain distance and see what we can get */
				Trajectory *temp_traj=new Trajectory(-1);
				double dist_did=0;
				if(trace_to_comp_nearby_minRoad(pt, cell, majormin, 20,
					search_dist, search_dist/2., dist_did, temp_traj))
				{
					for(int k=0;k<temp_traj->nlinesegs;k++)
					{
						if(temp_traj->linesegs[k].Triangle_ID<0
							||temp_traj->linesegs[k].Triangle_ID>=quadmesh->nfaces)
						{
							int test=0;
						}
					}
					if(!startorend && dist_did < 2*distToIntersect)
					{
						temp_traj->reverse_lines();
						curtraj->add_front_nlines(temp_traj->linesegs, temp_traj->nlinesegs);
					}
					else if(startorend && dist_did < 2*distToIntersect)
						curtraj->add_last_nlines(temp_traj->linesegs, temp_traj->nlinesegs);
				}
				delete temp_traj;
			}
		}
		else
		{
			Trajectory *temp_traj=new Trajectory(-1);
			double dist_did=0;
			if(trace_to_comp_nearby_minRoad(pt, cell, majormin, 20,
				search_dist, search_dist/2., dist_did, temp_traj))
			{
				//testing
				for(int k=0;k<temp_traj->nlinesegs;k++)
				{
					if(temp_traj->linesegs[k].Triangle_ID<0
						||temp_traj->linesegs[k].Triangle_ID>=quadmesh->nfaces)
					{
						int test=0;
					}
				}

				if(!startorend && dist_did < 2*distToIntersect)
				{
					temp_traj->reverse_lines();
					curtraj->add_front_nlines(temp_traj->linesegs, temp_traj->nlinesegs);
				}
				else if(startorend && dist_did < 2*distToIntersect)
					curtraj->add_last_nlines(temp_traj->linesegs, temp_traj->nlinesegs);
			}
			delete temp_traj;
		}
	}

}

/*
*/
void compute_extended_intersect(double p[2], int cellid, double intersect[2], 
								EvenStreamlinePlace *curplace, int fieldtype)
{
	int i;
	double smallest_dist=1.e50;
	icVector2 dist;
	int chosen_line=0;
	LineInfo *lineinfo;
	Trajectory *traj;
	int startid, endid, middle;

	if(fieldtype == 0) /*major try to find the intersection with an existing minor*/
	{
		lineinfo=quadmesh->quadcells[cellid]->minorlines;
	}
	else
	{
		lineinfo=quadmesh->quadcells[cellid]->majorlines;
	}


	normalize(tenline_dir_global);

	double A[2],B[2],C[2],D[2], t[2];

	A[0]=p[0];
	A[1]=p[1];
	B[0]=p[0]+quadmesh->xinterval*tenline_dir_global.entry[0];
	B[1]=p[1]+quadmesh->xinterval*tenline_dir_global.entry[1];
	int j;

	for(i=0;i<lineinfo->nlines;i++)
	{
		traj=curplace->evenstreamlines->trajs[lineinfo->lines[i]->whichtraj];

		if(fieldtype==1)
			traj=t_major->trajs[lineinfo->lines[i]->whichtraj];

		startid=lineinfo->lines[i]->start;
		endid=lineinfo->lines[i]->end;

		for(j=startid; j<=endid; j++)
		{
			C[0]=traj->linesegs[j].gstart[0];
			C[1]=traj->linesegs[j].gstart[1];
			D[0]=traj->linesegs[j].gend[0];
			D[1]=traj->linesegs[j].gend[1];

			if(cal_intersect(A, B, C, D, t)==1)
			{
				intersect[0]=A[0]+t[0]*(B[0]-A[0]);
				intersect[1]=A[1]+t[0]*(B[1]-A[1]);
				return;
			}
		}
	}

	/*if we cannot find the more accurate intersection using above method, we
	can use the following approximation*/
	for(i=0;i<lineinfo->nlines;i++)
	{
		traj=curplace->evenstreamlines->trajs[lineinfo->lines[i]->whichtraj];
		
		if(fieldtype==1)
			traj=t_major->trajs[lineinfo->lines[i]->whichtraj];

		startid=lineinfo->lines[i]->start;
		endid=lineinfo->lines[i]->end;
		double temp_intersect[2];

		for(j=startid;j<=endid;j++)
		{
			temp_intersect[0]=(traj->linesegs[j].gstart[0]+traj->linesegs[j].gend[0])/2;
			temp_intersect[1]=(traj->linesegs[j].gstart[1]+traj->linesegs[j].gend[1])/2;

			dist.entry[0]=temp_intersect[0]-p[0];
			dist.entry[1]=temp_intersect[1]-p[1];
			double temp_dist=length(dist);
			if(temp_dist<smallest_dist)
			{
				smallest_dist=temp_dist;
				chosen_line=j;
				intersect[0]=temp_intersect[0];
				intersect[1]=temp_intersect[1];
			}
		}
	}
}

/*
In the following routine, we search only one direction
11/18/2007
*/

bool find_cell_contain_bothmajmin(double start[2], int &startcell, double dsep, 
								  icVector2 line_dir, int fieldtype)
{
	int i;
	int flag = 0;
	double globalp[2];
	int pre_face, cur_face;

	pre_face = cur_face = startcell;
	globalp[0] = start[0];   globalp[1] = start[1];
	
	euler_stepsize = quadmesh->xinterval/4.;
	predict_stepsize = quadmesh->xinterval/2.;

	////If current trajectories list is not enough to store all the trajectories
	////extend it!
	if(cur_traj_index >= MaxNumTrajectories)
	{
		MaxNumTrajectories += 100;

		num_linesegs_curtraj = (int*)realloc(num_linesegs_curtraj, sizeof(int)*MaxNumTrajectories);
		trajectories = (LineSeg **)realloc(trajectories, sizeof(LineSeg *) * MaxNumTrajectories);

		////extend the old trajectories
		for(i = 0; i < MaxNumTrajectories-100; i++)
			trajectories[i] = (LineSeg *)realloc(trajectories[i], sizeof(LineSeg ) * MaxNumLinesegsPerTraj);

		////allocate the new trajectories
		for( ; i < MaxNumTrajectories; i++)
		{
			trajectories[i] = (LineSeg *)malloc(sizeof(LineSeg) * MaxNumLinesegsPerTraj);
			num_linesegs_curtraj[i] = 0;
		}
	}

	num_linesegs_curtraj[cur_traj_index] = 0;

	/*we need to pick a direction for it first*/

	tenline_dir_global = line_dir;  /*obtain the major eigen vector*/

	int NUMTRACINGCELLS = (int)(dsep/quadmesh->xinterval+2);
	for(i = 0; i < NUMTRACINGCELLS; i++)
	{
		if(cur_face == -1 || cur_face>=quadmesh->nfaces) /*reach boundary*/
		{
			return false;
		}

		/*if it reaches any boundary, we should stop as well*/
		if(globalp[0]>quadmesh->xend-1.e-8||globalp[0]<quadmesh->xstart+1.e-8
			||globalp[1]>quadmesh->yend-1.e-8||globalp[1]<quadmesh->ystart+1.e-8)
			return false;
		
		/*major searching for minor*/
		if(fieldtype == 0 && quadmesh->quadcells[cur_face]->hasminor)
		{
			startcell=cur_face;
			start[0]=globalp[0];
			start[1]=globalp[1];
			return true;
		}

		else if(fieldtype ==1 && quadmesh->quadcells[cur_face]->hasmajor)
			/*minor searching for major*/
		{
			/**/
			start[0]=globalp[0];
			start[1]=globalp[1];
			startcell=cur_face;
			return true;
		}

		pre_face = cur_face;
		cur_face = trace_ten_in_one_tri_quad(cur_face, globalp, flag, fieldtype); 


		if(flag == 1 || flag == 2 || pre_face == cur_face) 
		{
			return false;
		}
	}

	return false;
}

/************************************************************************************/
/************************************************************************************/
/************************************************************************************/
/************************************************************************************/
/*The definitions of the member functions of class Trajectory 09/29/2007*/

Trajectory::Trajectory(int index, int curMaxNum)
{
	this->index = index;

	if(curMaxNum==0)
	{
		linesegs=NULL;
		curMaxNumLinesegs=0;
		return;
	}

	linesegs = (LineSeg *)malloc(sizeof(LineSeg)*curMaxNum);

	if(linesegs == NULL)
	{
		exit(-1);
	}

	int i;
	for(i=0;i<curMaxNum;i++)
	{
		linesegs[i].gend[0]=linesegs[i].end[0]=linesegs[i].gstart[0]=linesegs[i].start[0]=
			linesegs[i].gend[1]=linesegs[i].end[1]=linesegs[i].gstart[1]=linesegs[i].start[1]=0.;
		linesegs[i].length=0;
		linesegs[i].Triangle_ID=0;
	}
	
	curMaxNumLinesegs = curMaxNum;
	nlinesegs = 0;

	//eulerstep_scalar = 	0.004 * object->radius; 

}


bool Trajectory::extend_line_segments(int add_size)
{

	LineSeg *extendlist=(LineSeg*)malloc(sizeof(LineSeg)*(curMaxNumLinesegs+add_size));

	if(extendlist == NULL)
	//if(linesegs == NULL)
	{
		return false;
	}

	int i;
	for(i = 0; i < curMaxNumLinesegs; i++)
	{
		extendlist[i].end[0] = linesegs[i].end[0];
		extendlist[i].end[1] = linesegs[i].end[1];

		extendlist[i].start[0] = linesegs[i].start[0];
		extendlist[i].start[1] = linesegs[i].start[1];

		extendlist[i].gend[0] = linesegs[i].gend[0];
		extendlist[i].gend[1] = linesegs[i].gend[1];

		extendlist[i].gstart[0] = linesegs[i].gstart[0];
		extendlist[i].gstart[1] = linesegs[i].gstart[1];

		extendlist[i].length = linesegs[i].length;
		extendlist[i].Triangle_ID = linesegs[i].Triangle_ID;
		
	}

	FILE *fp;

	free(linesegs);
	
	linesegs = extendlist;

	for(i=curMaxNumLinesegs;i<curMaxNumLinesegs+add_size;i++)
	{
		linesegs[i].gend[0]=linesegs[i].end[0]=linesegs[i].gstart[0]=linesegs[i].start[0]=
			linesegs[i].gend[1]=linesegs[i].end[1]=linesegs[i].gstart[1]=linesegs[i].start[1]=0.;
		linesegs[i].length=0;
		linesegs[i].Triangle_ID=0;
	}

	curMaxNumLinesegs += add_size;
	return true;
}

/*get the flow length of the streamline*/
double Trajectory::get_length()
{
	int i;
	double len = 0;
	for(i = 0 ; i < nlinesegs; i++)
		len += linesegs[i].length;
	return len;
}


//remove the front n line segments
bool Trajectory::remove_front_nlines(int n)
{
	if(nlinesegs-n<0) return false;
	/*move the content forward*/
	int i;
	for(i=0;i<nlinesegs-n;i++)
	{
		linesegs[i].gstart[0]=linesegs[i+n].gstart[0];
		linesegs[i].gstart[1]=linesegs[i+n].gstart[1];
		linesegs[i].gend[0]=linesegs[i+n].gend[0];
		linesegs[i].gend[1]=linesegs[i+n].gend[1];
		
		linesegs[i].start[0]=linesegs[i+n].start[0];
		linesegs[i].start[1]=linesegs[i+n].start[1];
		linesegs[i].end[0]=linesegs[i+n].end[0];
		linesegs[i].end[1]=linesegs[i+n].end[1];

		linesegs[i].length=linesegs[i+n].length;
		linesegs[i].Triangle_ID=linesegs[i+n].Triangle_ID;
	}
	nlinesegs-=n;
	return true;
}

//add n new line segments in the front
bool Trajectory::add_front_nlines(LineSeg *otherlinesegs, int n)
{
	if(nlinesegs+n>=curMaxNumLinesegs)
	{
		if(!extend_line_segments(nlinesegs+n-curMaxNumLinesegs))
			exit(-1);
	}

	/*move backward n elements*/
	int i;
	if(nlinesegs>0)
	{
		for(i=nlinesegs-1;i>=0;i--)
		{
			linesegs[i+n].gstart[0]=linesegs[i].gstart[0];
			linesegs[i+n].gstart[1]=linesegs[i].gstart[1];
			linesegs[i+n].gend[0]=linesegs[i].gend[0];
			linesegs[i+n].gend[1]=linesegs[i].gend[1];
			
			linesegs[i+n].start[0]=linesegs[i].start[0];
			linesegs[i+n].start[1]=linesegs[i].start[1];
			linesegs[i+n].end[0]=linesegs[i].end[0];
			linesegs[i+n].end[1]=linesegs[i].end[1];

			linesegs[i+n].length=linesegs[i].length;
			linesegs[i+n].Triangle_ID=linesegs[i].Triangle_ID;
		}
	}

	/*copy the new n line segments to the front*/
	for(i=0;i<n;i++)
	{
		linesegs[i].gstart[0]=otherlinesegs[i].gstart[0];
		linesegs[i].gstart[1]=otherlinesegs[i].gstart[1];
		linesegs[i].gend[0]=otherlinesegs[i].gend[0];
		linesegs[i].gend[1]=otherlinesegs[i].gend[1];
		
		linesegs[i].start[0]=otherlinesegs[i].start[0];
		linesegs[i].start[1]=otherlinesegs[i].start[1];
		linesegs[i].end[0]=otherlinesegs[i].end[0];
		linesegs[i].end[1]=otherlinesegs[i].end[1];

		linesegs[i].length=otherlinesegs[i].length;
		linesegs[i].Triangle_ID=otherlinesegs[i].Triangle_ID;
	}
	nlinesegs+=n;
	return true;
}

//remove the last n line segments
bool Trajectory::remove_last_nlines(int n)
{
	if(nlinesegs-n<0) return false;
	nlinesegs-=n;
	return true;
}

//add n new line segments at the end
bool Trajectory::add_last_nlines(LineSeg *otherlinesegs, int n)
{
	if(nlinesegs+n>=curMaxNumLinesegs)
	{
		if(!extend_line_segments(nlinesegs+n-curMaxNumLinesegs))
			exit(-1);
	}

	/*copy the content of "linesegs" to the end of current list*/
	int i;
	for(i=nlinesegs;i<n+nlinesegs;i++)
	{
		linesegs[i].gstart[0]=otherlinesegs[i-nlinesegs].gstart[0];
		linesegs[i].gstart[1]=otherlinesegs[i-nlinesegs].gstart[1];
		linesegs[i].gend[0]=otherlinesegs[i-nlinesegs].gend[0];
		linesegs[i].gend[1]=otherlinesegs[i-nlinesegs].gend[1];
		
		linesegs[i].start[0]=otherlinesegs[i-nlinesegs].start[0];
		linesegs[i].start[1]=otherlinesegs[i-nlinesegs].start[1];
		linesegs[i].end[0]=otherlinesegs[i-nlinesegs].end[0];
		linesegs[i].end[1]=otherlinesegs[i-nlinesegs].end[1];

		linesegs[i].length=otherlinesegs[i-nlinesegs].length;
		linesegs[i].Triangle_ID=otherlinesegs[i-nlinesegs].Triangle_ID;
	}
	nlinesegs+=n;
	return true;
}


bool Trajectory::reverse_lines()
{
	int i;

	LineSeg *temp = (LineSeg *)malloc(sizeof(LineSeg)*(this->nlinesegs+1));
		
	if(temp == NULL)
	{
		return false;
	}

	int newnum_lines = 0;

	////store the line segment in reversed order

	for(i = this->nlinesegs-1; i >= 0; i--)
	{
		if(this->linesegs[i].Triangle_ID < 0
			|| this->linesegs[i].Triangle_ID >= quadmesh->nfaces  
			|| this->linesegs[i].length < 0)
		{
			continue;
		}

		temp[newnum_lines].gstart[0] = this->linesegs[i].gend[0];
		temp[newnum_lines].gstart[1] = this->linesegs[i].gend[1];
		
		temp[newnum_lines].gend[0] = this->linesegs[i].gstart[0];
		temp[newnum_lines].gend[1] = this->linesegs[i].gstart[1];
		
		temp[newnum_lines].start[0] = this->linesegs[i].end[0];
		temp[newnum_lines].start[1] = this->linesegs[i].end[1];

		temp[newnum_lines].end[0] = this->linesegs[i].start[0];
		temp[newnum_lines].end[1] = this->linesegs[i].start[1];

		temp[newnum_lines].length = this->linesegs[i].length;
		temp[newnum_lines].Triangle_ID = this->linesegs[i].Triangle_ID;


		newnum_lines++;
	}

	////Copy it back to the origin array
	for(i = 0; i < newnum_lines; i++)
	{
		this->linesegs[i].gstart[0] = temp[i].gstart[0];
		this->linesegs[i].gstart[1] = temp[i].gstart[1];
		
		this->linesegs[i].gend[0] = temp[i].gend[0];
		this->linesegs[i].gend[1] = temp[i].gend[1];

		this->linesegs[i].start[0] = temp[i].start[0];
		this->linesegs[i].start[1] = temp[i].start[1];
		
		this->linesegs[i].end[0] = temp[i].end[0];
		this->linesegs[i].end[1] = temp[i].end[1];

		this->linesegs[i].length = temp[i].length;
		this->linesegs[i].Triangle_ID = temp[i].Triangle_ID;
	}

	this->nlinesegs = newnum_lines;

	free(temp);
}


/*****************************************************************
store temp curve point array to global line segment array
*****************************************************************/

bool Trajectory::store_to_global_line_segs(CurvePoints *temp, int num)
{
	int i;
	int tempid = nlinesegs;
	icVector3 dis_vec;

	////if the number of the line segements over the maximum number of the line segments each trajectory can store
	////extend the space for each trajectory
	if(tempid + num - 1 >= curMaxNumLinesegs)
	{
		//if(curMaxNumLinesegs>1000) 
		//	return false; // possible bug here! 12/27/2007

		if(!extend_line_segments(200))
		{
			return false;
		}
	}

	/*save to the global list*/

	for( i = 0; i < num-1; i++)
	{
		////Build the line segment
		linesegs[tempid+i].gstart[0] = temp[i].gpx;
		linesegs[tempid+i].gstart[1] = temp[i].gpy;
		linesegs[tempid+i].start[0] = temp[i].lpx;
		linesegs[tempid+i].start[1] = temp[i].lpy;

		linesegs[tempid+i].gend[0] = temp[i+1].gpx;
		linesegs[tempid+i].gend[1] = temp[i+1].gpy;
		linesegs[tempid+i].end[0] = temp[i+1].lpx;
		linesegs[tempid+i].end[1] = temp[i+1].lpy;

		////Use local coordinates to calculate the length
		dis_vec.entry[0] = temp[i+1].gpx - temp[i].gpx;
		dis_vec.entry[1] = temp[i+1].gpy - temp[i].gpy;
		dis_vec.entry[2] = 0;

		linesegs[tempid+i].length = length(dis_vec);

		linesegs[tempid+i].Triangle_ID = temp[i].triangleid;
	}

	nlinesegs = tempid + num - 1;
	return true;
}



int Trajectory::trace_in_quad(int &face_id, double globalp[2], int type, int &flag)
{

	int i;
	double pre_point[2];

	/*  will this be a good solution? 1/9/2008 */
	if(!is_in_cell(face_id, globalp[0], globalp[1]))
	{
		face_id = get_cellID_givencoords(globalp[0], globalp[1]);
	}

	if(face_id < 0 || face_id>=quadmesh->nfaces)
		return -1;

	QuadCell *face = quadmesh->quadcells[face_id];

	QuadCell *pre_f = face;
	
	////Temporary curve point array

	CurvePoints *temp_point_list = (CurvePoints*) malloc(sizeof(CurvePoints) * 200);

	if(temp_point_list == NULL)
	{
		exit(-1);
	}

	int NumPoints = 0;
	
	/*the tracing will be performed under the global frame*/
	globalface = face_id;

	pre_point[0] = globalp[0];
	pre_point[1] = globalp[1];


	////////////////////////////////////////////////////
    for(i = 0; i < 200; i++)
	{

		////2. if current point is inside current triangle
		if(is_in_cell(face_id, globalp[0], globalp[1]))
		{
			////store the point into the temp curve points list

			temp_point_list[NumPoints].gpx = globalp[0];
			temp_point_list[NumPoints].gpy = globalp[1];
			temp_point_list[NumPoints].triangleid = face->index;  
			NumPoints++;

			pre_point[0] = globalp[0];
			pre_point[1] = globalp[1];
			
			/*change to use other integration scheme 07/09/07*/
			//if(compute_next_pt_tensor_quad_global(pre_point, globalp, face_id))
			if(get_nextpt_2ndeuler_ten_quad(pre_point, globalp, face_id, type))
			//if(get_nextpt_RK23_ten_quad(pre_point, globalp, face_id, type))
			//if(get_nextpt_RK45_ten_quad(pre_point, globalp, face_id, type))
			{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global_p = tenline_dir_global;

				tenline_dir_global.entry[0] = globalp[0] - pre_point[0];
				tenline_dir_global.entry[1] = globalp[1] - pre_point[1];

			}

			else{  ////the curve reach a singularity/degenerate point
				flag = 1;

				////Store the record into global line segment array
                
				if(store_to_global_line_segs(temp_point_list, NumPoints))
				{
					////Not enough memory
					flag = 4;
					free(temp_point_list);
					return face_id;
				}

				free(temp_point_list);

				return face_id;
			}
		}


		////3. if the point is out of current cell
		else{

			/*!!!!!!need to judge which cell it will enter!!!!!*/
			int PassVertornot = 0;

			get_next_cell_2(face_id, pre_point, globalp, PassVertornot, type);
			
			if(PassVertornot>0)  /*cross a vertex*/
			{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global_p = tenline_dir_global;

				tenline_dir_global.entry[0] = pre_point[0] - globalp[0];
				tenline_dir_global.entry[1] = pre_point[1] - globalp[1];

				/**/
				temp_point_list[NumPoints].gpx = globalp[0];
				temp_point_list[NumPoints].gpy = globalp[1];
				temp_point_list[NumPoints].triangleid = face->index;  ////cause problem 05/25/05
				NumPoints++;

				temp_point_list[NumPoints].gpx = pre_point[0];
				temp_point_list[NumPoints].gpy = pre_point[1];
				temp_point_list[NumPoints].triangleid = face_id;  ////cause problem 05/25/05
				NumPoints++;
			}
			else{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = globalp[0] - pre_point[0];
				tenline_dir_global.entry[1] = globalp[1] - pre_point[1];

				////Add the intersection point to the temporary points' list
				temp_point_list[NumPoints].gpx = globalp[0];
				temp_point_list[NumPoints].gpy = globalp[1];
				temp_point_list[NumPoints].triangleid = face->index;  ////cause problem 05/25/05
				NumPoints++;
			}

			//if(face_id<0||face_id>=quadmesh->nfaces)
			//	{
			//		int test=0;
			//	}

			if(NumPoints > 1){
 				////Store the record into global line segment array
				if(!store_to_global_line_segs(temp_point_list, NumPoints))
			   {   ////Not enough memory
				   flag = 4;
				   free(temp_point_list);
				   return face_id;
			   }
			}

			free(temp_point_list);
			return face_id;
		}
	
	}

	if(NumPoints > 0)
		if(!store_to_global_line_segs(temp_point_list, NumPoints))
			flag=4;

	free(temp_point_list);

	return face_id;
}



/************************************************************************************/
/************************************************************************************/
/************************************************************************************/
/************************************************************************************/
/************************************************************************************/

#include ".\glview.h"

/*test the seeds and samples*/
void EvenStreamlinePlace::display_sampts_seeds()  
{
	//SeedList *seedpts;
	//SamplePtList **samplepts;
	int i;

	/*display seeds as green dots*/
	glPointSize(3.);
	glColor3f(0, 1, 0);
	glBegin(GL_POINTS);
	for(i=0; i<seedpts->nseeds; i++)
	{
		glVertex2f(seedpts->seeds[i]->pos[0],seedpts->seeds[i]->pos[1]);
	}
	glEnd();

	/*display samples as blue dots*/
	int j;
	SamplePtList *samplist;
	glColor3f(0, 0, 1);
	glBegin(GL_POINTS);
	for(i=0; i<evenstreamlines->ntrajs; i++)
	{
		samplist=samplepts[i];
		for(j=0; j<samplist->nsamples; j++)
			glVertex2f(samplist->samples[j]->gpt[0],samplist->samples[j]->gpt[1]);
	}
	glEnd();
}

void EvenStreamlinePlace::reset_placement_quad()
{
	int i;
	QuadCell *face;

	for(i = 0; i < quadmesh->nfaces; i++)
	{
		face = quadmesh->quadcells[i];
		face->visited = false;

		//face->reset_sampleList();
	}
	init_samplelist_in_cell(majororminor);

	for(i = 0; i < quadmesh->nverts; i++)
	{
		quadmesh->quad_verts[i]->distance = 1e49;
		quadmesh->quad_verts[i]->visited = false;
	}
}

/*initialize the information of the tensor lines for the cells that
they passed.
This information will be used to compute the intersections of
major lines and minor lines efficiently*/
void EvenStreamlinePlace::init_major_minor_line_info(bool type)
{
	int i;
	QuadCell *face;
	for(i=0; i<quadmesh->nfaces; i++)
	{
		face = quadmesh->quadcells[i];
		if(type == false)
		{
			if(face->majorlines != NULL)
			{
				delete face->majorlines;
				face->majorlines = NULL;
			}
			face->hasmajor = false;
		}
		else{
			if(face->minorlines != NULL)
			{
				delete face->minorlines;
				face->minorlines = NULL;
			}
			face->hasminor = false;
		}
	}
}


extern ctr_point *out_pts; /*the Hermite output curve*/

void EvenStreamlinePlace::set_default_parameters(bool fieldtype)
{
	if(!fieldtype)
		dsep = majorDensity*quadmesh->xinterval;    //using the radius of the object instead of the edge
	else
		dsep = minorDensity*quadmesh->xinterval;    //using the radius of the object instead of the edge

	if(which_level==1)
		percentage_dsep = 0.5;
	else
		percentage_dsep = 0.9;

	discsize = 2.;
	if(which_level==1)
		sample_interval = min(0.05*dsep, quadmesh->xinterval/2.);
	else
		sample_interval = 0.25*dsep;

	every_nsample = 4;
	loopdsep = 0.4*quadmesh->xinterval; /*we now allow the loops to be closed*/
	dist2sing = 0.1*quadmesh->xinterval;

	streamlinelength = mintenline_length*quadmesh->xinterval;
	seeddist = 0.9*dsep;
	minstartdist = 0.9*dsep;

	euler_stepsize = quadmesh->xinterval/5.;
	predict_stepsize = quadmesh->xinterval;

	cur_traj=0;
	cur_seed_pos=0;
}


/*         Copy the user sketches to be part of the obtained roads          */

void EvenStreamlinePlace::copy_sketchlines()
{
	int i;
	int cur_line, movetonext;
	double cur_length;
	evenstreamlines->ntrajs=0;
	for(i=0;i<sketchlines->ntrajs;i++)
	{
		evenstreamlines->trajs[i]->nlinesegs=0;
		evenstreamlines->trajs[i]->add_front_nlines(sketchlines->trajs[i]->linesegs,
			sketchlines->trajs[i]->nlinesegs);

		/*do the sampling*/
		cur_line = 0;
		cur_length = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].length;
		samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[0] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
		samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[1] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
		samplepts[evenstreamlines->ntrajs]->samples[0]->triangle = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;
		samplepts[evenstreamlines->ntrajs]->nsamples = 1;

		movetonext = 0;
		cal_samplepts_when_tracing(evenstreamlines->ntrajs, sample_interval, cur_line, movetonext, cur_length,
				samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
			
		update_samples_in_cell(evenstreamlines->ntrajs, samplepts[evenstreamlines->ntrajs]->samples,
			samplepts[evenstreamlines->ntrajs]->nsamples);


		/*  compute the line information */
		rebuild_lineinfo_for_a_tensorline(i);

		evenstreamlines->trajs[i]->roadtype=HIGHWAY;

		evenstreamlines->ntrajs++;
	}
	
	/*     We need to copy the sample points to a global seed point list which will be
	       used by the minor line placement
	*/
	if(seedsalongbounds != NULL)
	{
		delete seedsalongbounds;
		seedsalongbounds=NULL;
	}
	seedsalongbounds=new SeedList(1000);

	for(i=0;i<sketchlines->ntrajs;i++)
	{
		int j;
		for(j=0;j<samplepts[i]->nsamples;j++)
		{
			if(j%every_nsample!=0)
				continue;
						
			/*   create a new seed   */
			Seed *s = (Seed *) malloc(sizeof(Seed));
			if(s == NULL)
			{
				return;
			}

			/*  NOTE: if it is boundary sketch, we need to move away from
			    water! 12/29/2007
			*/

			s->pos[0] = samplepts[i]->samples[j]->gpt[0];
			s->pos[1] = samplepts[i]->samples[j]->gpt[1];
			
			s->triangle = samplepts[i]->samples[j]->triangle;
			s->state = 0; //set it is active 

			if(samplepts[i]->samples[j]->triangle<0||
				samplepts[i]->samples[j]->triangle>=quadmesh->nfaces)
			{
				//int test=0;
				continue;
			}

			seedsalongbounds->append(s);
		}
	}

}


/*      Copy the obtained major roads as part of the roads      */

void EvenStreamlinePlace::copy_majorRoads(bool majormin)
{
	int i;
	int cur_line, movetonext;
	double cur_length;
	//evenstreamlines->ntrajs=0;
	EvenStreamlinePlace *curplace;

	if(!majormin)
	{
		curplace=major_level1;  /* level 1 tracing results:  major direction */
	}
	else
	{
		curplace=minor_level1;  /* level 1 tracing results:  minor direction */
	}

	for(i=0;i<curplace->evenstreamlines->ntrajs;i++)
	{
		evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs=0;

		if(curplace->evenstreamlines->trajs[i/*evenstreamlines->ntrajs*/]->nlinesegs<=0)
		{
			//int test = 0;
			continue;
		}

		//if(curplace->evenstreamlines->trajs[i]->is_mapboundary
		//	&&!sharedvars.UseBoundsAsRoadsOn)
		//	continue;

		evenstreamlines->trajs[evenstreamlines->ntrajs]->add_front_nlines(
			curplace->evenstreamlines->trajs[i/*evenstreamlines->ntrajs*/]->linesegs,
			curplace->evenstreamlines->trajs[i/*evenstreamlines->ntrajs*/]->nlinesegs);



		/*do the sampling*/
		cur_line = 0;
		cur_length = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].length;
		samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[0] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
		samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[1] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
		samplepts[evenstreamlines->ntrajs]->samples[0]->triangle = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;
		samplepts[evenstreamlines->ntrajs]->nsamples = 1;

		movetonext = 0;
		cal_samplepts_when_tracing(evenstreamlines->ntrajs, sample_interval, cur_line, movetonext, cur_length,
				samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
			
		update_samples_in_cell(evenstreamlines->ntrajs, samplepts[evenstreamlines->ntrajs]->samples,
			samplepts[evenstreamlines->ntrajs]->nsamples);


		/*  compute the line information */
		//if(curplace->evenstreamlines->trajs[i]->is_mapboundary
		//	&&!sharedvars.UseBoundsAsRoadsOn)
		//{
		//}
		//else
			rebuild_lineinfo_for_a_tensorline(evenstreamlines->ntrajs);

		evenstreamlines->trajs[evenstreamlines->ntrajs]->roadtype=MAJOR;
		evenstreamlines->trajs[evenstreamlines->ntrajs]->is_mapboundary
			=curplace->evenstreamlines->trajs[i]->is_mapboundary;

		evenstreamlines->trajs[evenstreamlines->ntrajs]->closed=
			curplace->evenstreamlines->trajs[i]->closed;

		evenstreamlines->ntrajs++;

	}
	
	/*     We need to copy the sample points to a global seed point list which will be
	       used by the minor line placement
	*/
	//if(seedsalongbounds != NULL)
	//{
	//	delete seedsalongbounds;
	//	seedsalongbounds=NULL;
	//}
	//seedsalongbounds=new SeedList(1000);

	//for(i=0;i<sketchlines->ntrajs;i++)
	//{
	//	int j;
	//	for(j=0;j<samplepts[i]->nsamples;j++)
	//	{
	//		if(j%every_nsample!=0)
	//			continue;
	//					
	//		/*   create a new seed   */
	//		Seed *s = (Seed *) malloc(sizeof(Seed));
	//		if(s == NULL)
	//		{
	//			return;
	//		}

	//		s->pos[0] = samplepts[i]->samples[j]->gpt[0];
	//		s->pos[1] = samplepts[i]->samples[j]->gpt[1];
	//		
	//		s->triangle = samplepts[i]->samples[j]->triangle;
	//		s->state = 0; //set it is active 

	//		if(samplepts[i]->samples[j]->triangle<0||
	//			samplepts[i]->samples[j]->triangle>=quadmesh->nfaces)
	//		{
	//			int test=0;
	//			continue;
	//		}

	//		seedsalongbounds->append(s);
	//	}
	//}

}


void EvenStreamlinePlace::sample_majorRoads_seeds(SeedList *seedlist, int samp_rate)
{
	/*  for the minor direction, we need to sample along the sketched major roads
	    to obtain the initial seeds
	*/

	int i, j;
	Trajectory *curtraj;
	for(i=0;i<major_level1->evenstreamlines->ntrajs;i++)
	{
		curtraj=major_level1->evenstreamlines->trajs[i];

		for(j=0;j<curtraj->nlinesegs;j++)
		{
			if(j%samp_rate==0)
			{
				Seed *s=(Seed*)malloc(sizeof(Seed));
				s->pos[0]=curtraj->linesegs[j].gend[0];
				s->pos[1]=curtraj->linesegs[j].gend[1];
				s->triangle=curtraj->linesegs[j].Triangle_ID;
				s->state=0;

				if(s->triangle<0||s->triangle>=quadmesh->nfaces)
				{
					free(s);
					continue;
				}

				seedlist->append(s);
			}
		}
	}
}



///*
//Use sample points to control the termination of streamlines
//*/
////void EvenStreamlinePlace::place_streamlines(int num_initial, double dsep, double percentage_dsep, 
////				  double discsize, double sample_interval, int every_nsample, 
////				  double loopdsep, double dist2sing, double streamlinelength, 
////				  double seeddist, double minstartdist, int flag)

void EvenStreamlinePlace::place_streamlines(int type, bool brushon)
{
	int i;
	int locate_preseed = 0;
	int begin_sample = 0;
    
	/*-----Set a set of default initial values for testing 5/2/06-----*/
	/*--- this setting should be provided by user ---*/

	//dsep = /*0.04 * quadmesh->radius*/1.5*quadmesh->xinterval;    //using the radius of the object instead of the edge
	
	//dsep = 1.5*quadmesh->xinterval;    //using the radius of the object instead of the edge

	if(type == 0)
		dsep = majorDensity*quadmesh->xinterval;    //using the radius of the object instead of the edge
	else
		dsep = minorDensity*quadmesh->xinterval;    //using the radius of the object instead of the edge

	if(which_level==1)
		percentage_dsep = 0.15;
	else
		percentage_dsep = 0.84;
	discsize = 2.;
	sample_interval = 0.25*dsep;
	every_nsample = 4;
	loopdsep = min(0.4*dsep, quadmesh->xinterval/2.);
	//dist2sing = 0.1*dsep;
	//loopdsep = 0.2*quadmesh->xinterval; /*we now allow the loops to be closed*/
	dist2sing = 0.1*quadmesh->xinterval;

	//streamlinelength = /*0.03**/5*quadmesh->xinterval;
	streamlinelength = mintenline_length*quadmesh->xinterval;
	seeddist = 0.9*dsep;
	minstartdist = 0.85*dsep;

	euler_stepsize = quadmesh->xinterval/5.;
	predict_stepsize = quadmesh->xinterval;

	double adj_dsep=dsep;
	double adj_seeddist=seeddist;
	double adj_startdist=minstartdist;


	/**************************************************************************/
	/*allocate space for seed points*/
	seedpts = new SeedList();

	cur_traj = 0; //set the first streamline as current streamline

	evenstreamlines->ntrajs = 0;


	/****************************************************************************/
	/*            Here, we consider to include the user specified highways
	*/

	if(highwayexisted && sharedvars.rdSketchMajorHighways)
	{
		/*    we copy the highway to the major *only*    */
		if(!majororminor)
		{
			copy_sketchlines();

			seedpts->nseeds = 0; //reset the seedpts list

			cal_seeds(cur_traj, seeddist, every_nsample, type, brushon);
		}
		else
		{
			/* for minor line placement, we make use of the seed points 
			   sampled from the highway
			*/
			if(seedsalongbounds!=NULL && seedsalongbounds->nseeds>0)
			{
				/*   we copy the seeds from the "seedsalongbounds"   */
				seedpts->copy_otherseedList(seedsalongbounds);
			}
			else{
				cal_init_streamlines(1, streamlinelength, percentage_dsep*dsep,
					discsize, sample_interval, loopdsep, dist2sing, type, brushon);

				seedpts->nseeds = 0; //reset the seedpts list

				cal_seeds(cur_traj, seeddist, every_nsample, type, brushon);
			}
		}
	}


	/****************************************************************************/
	/*            Here, we consider to include the level 1 tracing results
	*/

	if(majorroadsexisted)
	{
		/*    we copy the highway to the major *only*    */
		if(!majororminor)
		{
			copy_majorRoads(false);

			//seedpts->nseeds = 0; //reset the seedpts list

			/*  we need to sample along the major roads (use fewer seeds, bug) */
			cal_seeds(cur_traj, seeddist, 100*every_nsample, type, brushon);

			/*   release the majorlevel1   */
		}
		else
		{
			copy_majorRoads(true);

			//seedpts->nseeds = 0; //reset the seedpts list

			cal_seeds(cur_traj, seeddist, every_nsample, type, brushon);

			/*  bug 1/10/2008  */
			if(!sharedvars.rdSketchMajorHighways && !flag_loadmap 
				&& minor_level1->evenstreamlines->ntrajs==0
				/*sharedvars.EnableSketchBasedDesign*/)
			{
				/*  we need to sample along the major roads (use fewer seeds, bug) */
				sample_majorRoads_seeds(seedpts, 100*every_nsample);
			}

			/*    release the minorlevel1    */
		}
	}

	/*   finally, we perform regular tracing if no highway and major roads existed  */

	else{
		cal_init_streamlines(1, streamlinelength, percentage_dsep*dsep,
			discsize, sample_interval, loopdsep, dist2sing, type, brushon);

		seedpts->nseeds = 0; //reset the seedpts list

		cal_seeds(cur_traj, seeddist, every_nsample, type, brushon);
	}

	/************************************/
	/*       The regular tracing        */
	/************************************/

	locate_preseed = 0;
	while(1 ) //if there still are streamlines not being processed
	{
		////at each step, we grow streamlines for all the seed points associated with current streamline
		for(i = locate_preseed; i < seedpts->nseeds; i++)
		{
			/*   combined with density map  11/25/2007   */
			if(sharedvars.CombinePopDensityOn)
			{
				double approx_den=cal_approx_density_at(seedpts->seeds[i]->pos,
					seedpts->seeds[i]->triangle);
				adj_startdist=minstartdist/approx_den;
			}


			////Judge whether this is a valid seed or not
			//if(!close_to_cur_streamline(seedpts->seeds[i]->pos,
			//	seedpts->seeds[i]->triangle, &cur_traj, 1, minstartdist, discsize, 0))

			if(!close_to_cur_streamline(seedpts->seeds[i]->pos,
				seedpts->seeds[i]->triangle, &cur_traj, 1, adj_startdist, discsize, 0))
			{
				////if we find a valid seed point, break
				locate_preseed = i+1;
				break;
			}
	
			seedpts->seeds[i]->state = 2;  //reject before starting
		}

		if(i >= seedpts->nseeds ) // it means no more seeds available for current streamline
		{
			locate_preseed = seedpts->nseeds;
			cur_traj ++;  //let next streamline as current streamline

			if(cur_traj >= evenstreamlines->ntrajs) //no more streamlines available
			{
				goto LL;
			}
			
			/*        Get seeds for the current streamline      */
			//cal_seeds(cur_traj, dsep, every_nsample, type, brushon);
			cal_seeds(cur_traj, seeddist, every_nsample, type, brushon);
			continue;
		}

		////else, we find a valid seed point, then grow a new streamline from this seed

		if(grow_a_tensorline(seedpts->seeds[i]->pos, seedpts->seeds[i]->triangle,
			percentage_dsep*dsep, discsize, sample_interval, loopdsep, dist2sing, 
			streamlinelength, type, brushon))
		{
			seedpts->seeds[i]->state = 1;

			if(evenstreamlines->isFull())
			{
				int oldnum = evenstreamlines->curMaxNumTrajs;
				if(!evenstreamlines->extend())
				{
					/*release the seed point list*/
					return;
				}

				/*allocate memeory for the elements*/
				for(i = oldnum; i < evenstreamlines->curMaxNumTrajs; i++)
				{
					evenstreamlines->trajs[i] = new Trajectory(i);
					if(evenstreamlines->trajs[i] == NULL)
					{
						return;
					}
				}

				/*increase the sample point lists as well*/
				SamplePtList **temp = samplepts;
				samplepts = (SamplePtList **)malloc(sizeof(SamplePtList*)*evenstreamlines->curMaxNumTrajs);
				if(samplepts == NULL)
				{
					/*report error*/
					if(samplepts == NULL)
					{
						return;
					}

					/*release the seed point list*/
					delete  seedpts;
					seedpts=NULL;
					return;
				}

				/*copy the information from the old list*/
				for(i = 0; i < oldnum; i++)
					samplepts[i] = temp[i];
				delete [] temp;

				/*allocate the pointers for the real samples*/
				for(i = oldnum; i < evenstreamlines->curMaxNumTrajs; i++)
				{
					samplepts[i] = new SamplePtList(100);
					if(samplepts[i] == NULL)
					{
						return;
					}

					for(int j = 0; j < samplepts[i]->curMaxNumSamplePts; j++)
					{
						samplepts[i]->samples[j] = (SamplePt *)malloc(sizeof(SamplePt));

						if(samplepts[i]->samples[j] == NULL)
						{
							return;
						}
					}
				}

				/*extend the list of seedposition_ineachtensorline*/
				seedposition_ineachtensorline=(int*)realloc(seedposition_ineachtensorline,
					sizeof(int)*evenstreamlines->curMaxNumTrajs);
				if(seedposition_ineachtensorline==NULL)
					return;
			}
		}

		else
		{
			seedpts->seeds[i]->state = 3;   //the seed rejected due to the short streamline
		}

	}

	/*release the seed point list*/
LL:	
	delete  seedpts;
	seedpts=NULL;

	/*release the sample point list*/
	for(i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
	{
		delete samplepts[i];
		samplepts[i] = NULL;
	}
	delete [] samplepts;
	samplepts = NULL;
	reset_placement_quad();
}



/*
The placement function with user specified seeds as one of the input
*/

void EvenStreamlinePlace::place_streamlines(int type, bool brushon, SeedList *preseeds)
{
	int i;
	int locate_preseed = 0;
	int begin_sample = 0;
    
	/*-----Set a set of default initial values for testing 5/2/06-----*/
	/*--- this setting should be provided by user ---*/

	if(type == 0)
		dsep = majorDensity*quadmesh->xinterval;    //using the radius of the object instead of the edge
	else
		dsep = minorDensity*quadmesh->xinterval;    //using the radius of the object instead of the edge

	//percentage_dsep = 0.8;
	if(which_level==1)
	{
		if(!sharedvars.JobardMethodOn)
			percentage_dsep = 0.15;
		else
			percentage_dsep = 0.84;
	}
	else
		percentage_dsep = 0.84;
	discsize = 2.;
	sample_interval = 0.25*dsep;
	every_nsample = 4;
	//loopdsep = 0.4*quadmesh->xinterval; /*we now allow the loops to be closed*/
	loopdsep = min(0.4*dsep, quadmesh->xinterval/2.);
	dist2sing = 0.1*quadmesh->xinterval;

	//streamlinelength = 5*quadmesh->xinterval;
	streamlinelength = mintenline_length*quadmesh->xinterval;
	
	//seeddist = 1.*dsep;
	//minstartdist = 0.9*dsep;
	seeddist = 0.9*dsep;
	minstartdist = 0.85*dsep;

	euler_stepsize = quadmesh->xinterval/5.;
	predict_stepsize = quadmesh->xinterval;

	double adj_dsep=dsep;
	double adj_seeddist=seeddist;
	double adj_startdist=minstartdist;

	/*******************************************************************/
	/*allocate space for seed points*/
	seedpts = new SeedList();

	/*copy the seed points here*/
	seedpts->copy_otherseedList(preseeds);

	cur_traj = 0; //set the first streamline as current streamline

	/*Record execution information into file*/
	FILE *fp;


	evenstreamlines->ntrajs = 0;

	/****************************************************************************/
	/*            Here, we consider to include the user specified highways
	*/

	if(highwayexisted && sharedvars.rdSketchMajorHighways)
	{
		/*    we copy the highway to the major *only*    */
		if(!majororminor)
		{
			copy_sketchlines();

			seedpts->nseeds = 0; //reset the seedpts list

			cal_seeds(cur_traj, seeddist, every_nsample, type, brushon);
		}
		else
		{
			/* for minor line placement, we make use of the seed points 
			   sampled from the highway
			*/
			if(seedsalongbounds!=NULL && seedsalongbounds->nseeds>0)
			{
				/*   we copy the seeds from the "seedsalongbounds"   */
				//seedpts->copy_otherseedList(seedsalongbounds);
			}
			else{
				cal_init_streamlines(1, streamlinelength, percentage_dsep*dsep,
					discsize, sample_interval, loopdsep, dist2sing, type, brushon);

				seedpts->nseeds = 0; //reset the seedpts list

				cal_seeds(cur_traj, seeddist, every_nsample, type, brushon);
			}
		}
	}


	/****************************************************************************/
	/*            Here, we consider to include the level 1 tracing results
	*/

	if(majorroadsexisted)
	{
		/*    we copy the highway to the major *only*    */
		if(!majororminor)
		{
			copy_majorRoads(false);

			cal_seeds(cur_traj, seeddist, 20*every_nsample, type, brushon);

			/*   release the majorlevel1   */
		}
		else
		{
			copy_majorRoads(true);

			cal_seeds(cur_traj, seeddist, every_nsample, type, brushon);

			if(sharedvars.rdSketchMajorHighways==0&&sharedvars.EnableSketchBasedDesign)
			{
				/*  we need to sample along the major roads  */
				sample_majorRoads_seeds(seedpts, 20*every_nsample);
			}

			/*    release the minorlevel1    */
		}
	}

	/*  Jobard's method 1/21/2008  */
	if(sharedvars.JobardMethodOn && which_level ==1 && preseeds->nseeds<=1 )
	{
		cal_init_streamlines(1, streamlinelength, percentage_dsep*dsep,
			discsize, sample_interval, loopdsep, dist2sing, type, brushon);

		seedpts->nseeds = 0; //reset the seedpts list

		cal_seeds(cur_traj, seeddist, every_nsample, type, brushon);
	}

	

	locate_preseed = 0;
	while(1) //if there still are streamlines not being processed
	{
		////at each step, we grow streamlines for all the seed points associated with current streamline
		for(i = locate_preseed; i < seedpts->nseeds; i++)
		{
			/*   combined with density map  11/25/2007   */
			if(sharedvars.CombinePopDensityOn)
			{
				double approx_den=cal_approx_density_at(seedpts->seeds[i]->pos,
					seedpts->seeds[i]->triangle);
				adj_startdist=minstartdist/approx_den;
			}


			////Judge whether this is a valid seed or not
			//if(!close_to_cur_streamline(seedpts->seeds[i]->pos,
			//	seedpts->seeds[i]->triangle, &cur_traj, 1, minstartdist, discsize, 0))

			if(!close_to_cur_streamline(seedpts->seeds[i]->pos,
				seedpts->seeds[i]->triangle, &cur_traj, 1, adj_startdist, discsize, 0))
			{
				////if we find a valid seed point, break
				locate_preseed = i+1;
				break;
			}
	
			seedpts->seeds[i]->state = 2;  //reject before starting
		}

		if(i >= seedpts->nseeds ) // it means no more seeds available for current streamline
		{
			locate_preseed = seedpts->nseeds;
			cur_traj ++;  //let next streamline as current streamline

			if(cur_traj >= evenstreamlines->ntrajs) //no more streamlines available
			{
				/*release the seed point list*/
				//delete seedpts;
				//return;
				goto LL;
			}
			
			//Get seeds for the current streamline
			cal_seeds(cur_traj, seeddist, every_nsample, type, brushon);

			//fp=fopen("lineseg_mem_error.txt", "w");
			//fprintf(fp, "finish finding seeds\n");
			//fprintf(fp, "We have %d seeds now", seedpts->nseeds);
			//fclose(fp);

			continue;
		}

		////else, we find a valid seed point, then grow a new streamline from this seed

		//fp=fopen("lineseg_mem_error.txt", "w");
		//if(!majororminor)
		//	fprintf(fp, "start tracing traj %d under major direction\n", evenstreamlines->ntrajs);
		//else
		//	fprintf(fp, "start tracing traj %d under minor direction\n", evenstreamlines->ntrajs);
		//fclose(fp);

		if(grow_a_tensorline(seedpts->seeds[i]->pos, seedpts->seeds[i]->triangle,
			percentage_dsep*dsep, discsize, sample_interval, loopdsep, dist2sing, 
			streamlinelength, type, brushon))
		{
			seedpts->seeds[i]->state = 1;

			//fp=fopen("lineseg_mem_error.txt", "w");
			//if(!majororminor)
			//	fprintf(fp, "finish tracing traj %d under major direction\n", evenstreamlines->ntrajs-1);
			//else
			//	fprintf(fp, "finish tracing traj %d under minor direction\n", evenstreamlines->ntrajs-1);
			//fprintf(fp, "it contains %d line segments.\n", evenstreamlines->trajs[
			//	evenstreamlines->ntrajs-1]->nlinesegs);
			//fclose(fp);

			if(evenstreamlines->isFull())
			{
				int oldnum = evenstreamlines->curMaxNumTrajs;
				if(!evenstreamlines->extend())
				{
					/*release the seed point list*/
					return;
				}

				/*allocate memeory for the elements*/
				for(i = oldnum; i < evenstreamlines->curMaxNumTrajs; i++)
				{
					evenstreamlines->trajs[i] = new Trajectory(i);
					if(evenstreamlines->trajs[i] == NULL)
					{
						return;
					}
				}

				/*increase the sample point lists as well*/
				SamplePtList **temp = samplepts;
				samplepts = (SamplePtList **)malloc(sizeof(SamplePtList*)*evenstreamlines->curMaxNumTrajs);
				if(samplepts == NULL)
				{
					/*report error*/
					if(samplepts == NULL)
					{
						return;
					}

					/*release the seed point list*/
					delete  seedpts;
					seedpts=NULL;
					return;
				}

				/*copy the information from the old list*/
				for(i = 0; i < oldnum; i++)
					samplepts[i] = temp[i];
				delete [] temp;

				/*allocate the pointers for the real samples*/
				for(i = oldnum; i < evenstreamlines->curMaxNumTrajs; i++)
				{
					samplepts[i] = new SamplePtList(100);
					if(samplepts[i] == NULL)
					{
						return;
					}

					for(int j = 0; j < samplepts[i]->curMaxNumSamplePts; j++)
					{
						samplepts[i]->samples[j] = (SamplePt *)malloc(sizeof(SamplePt));

						if(samplepts[i]->samples[j] == NULL)
						{
							return;
						}
					}
				}
				
				/*extend the list of seedposition_ineachtensorline*/
				seedposition_ineachtensorline=(int*)realloc(seedposition_ineachtensorline,
					sizeof(int)*evenstreamlines->curMaxNumTrajs);
				if(seedposition_ineachtensorline==NULL)
					return;
			}
		}

		else
		{
			seedpts->seeds[i]->state = 3;   //the seed rejected due to the short streamline
		}

	}

	/*release the seed point list*/
LL:	delete  seedpts;
	seedpts=NULL;

	/*release the sample point list*/
	for(i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
	{
		delete samplepts[i];
		samplepts[i] = NULL;
	}
	delete [] samplepts;
	samplepts = NULL;
	reset_placement_quad();
}



/*   Place tensor lines inside a user specified closed *region*
     This is useful for local street network editing
*/

void EvenStreamlinePlace::place_tensorlines_inReg(int type, SeedList *ini_seeds)
{
	int i;
	int locate_preseed = 0;
	int begin_sample = 0;
    
	/*-----Set a set of default initial values for testing 5/2/06-----*/
	/*--- this setting should be provided by user ---*/

	if(type == 0)
		dsep = majorDensity*quadmesh->xinterval;    //using the radius of the object instead of the edge
	else
		dsep = minorDensity*quadmesh->xinterval;    //using the radius of the object instead of the edge

	//percentage_dsep = 0.8;
	if(which_level==1)
		percentage_dsep = 0.55;
	else
		percentage_dsep = 0.84;
	discsize = 2.;
	sample_interval = 0.25*dsep;
	every_nsample = 4;

	loopdsep = min(0.5*dsep, quadmesh->xinterval);/*we now allow the loops to be closed*/
	dist2sing = 0.1*quadmesh->xinterval;

	//streamlinelength = 5*quadmesh->xinterval;
	streamlinelength = mintenline_length*quadmesh->xinterval;
	
	seeddist = 0.9*dsep;
	minstartdist = 0.85*dsep;

	euler_stepsize = quadmesh->xinterval/5.;
	predict_stepsize = quadmesh->xinterval;

	double adj_dsep=dsep;
	double adj_seeddist=seeddist;
	double adj_startdist=minstartdist;

	bool brushon=false;
	nnewlines_inReg=0;


	/*******************************************************************/
	/*allocate space for seed points*/
	seedpts = new SeedList();

	/* we also need to reallocate the space for sampling points */
		samplepts = (SamplePtList **)malloc(sizeof(SamplePtList*)*evenstreamlines->curMaxNumTrajs);
		for(i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
		{
			samplepts[i] = new SamplePtList(50);

			if(samplepts[i] == NULL)
			{
				exit(-1);
			}

			/*allocate memory for the elements*/
			for(int j = 0 ; j < samplepts[i]->curMaxNumSamplePts; j++)
			{
				samplepts[i]->samples[j] = (SamplePt *)malloc(sizeof(SamplePt));
			}
		}

	/*copy the seed points here*/
	seedpts->copy_otherseedList(ini_seeds);

	cur_traj = evenstreamlines->ntrajs-1; //set the first streamline as current streamline

	locate_preseed = 0;
	while(1) //if there still are streamlines not being processed
	{
		////at each step, we grow streamlines for all the seed points associated with current streamline
		for(i = locate_preseed; i < seedpts->nseeds; i++)
		{
			/*   combined with density map  11/25/2007   */
			if(sharedvars.CombinePopDensityOn)
			{
				double approx_den=cal_approx_density_at(seedpts->seeds[i]->pos,
					seedpts->seeds[i]->triangle);
				adj_startdist=minstartdist/approx_den;
			}

			////Judge whether this is a valid seed or not
			if(!quadmesh->quadcells[seedpts->seeds[i]->triangle]->in_region)
			{
				locate_preseed = i+1;
				break;
			}
			//if(!close_to_cur_streamline(seedpts->seeds[i]->pos,
			//	seedpts->seeds[i]->triangle, &cur_traj, 1, minstartdist, discsize, 0))
			if(!close_to_cur_streamline(seedpts->seeds[i]->pos,
				seedpts->seeds[i]->triangle, &cur_traj, 1, adj_startdist, discsize, 0))
			{
				////if we find a valid seed point, break
				locate_preseed = i+1;
				break;
			}
	
			seedpts->seeds[i]->state = 2;  //reject before starting
		}

		if(i >= seedpts->nseeds ) // it means no more seeds available for current streamline
		{
			locate_preseed = seedpts->nseeds;
			cur_traj ++;  //let next streamline as current streamline

			if(cur_traj >= evenstreamlines->ntrajs) //no more streamlines available
			{
				/*release the seed point list*/
				goto LL;
			}
			
			//Get seeds for the current streamline
			cal_seeds(cur_traj, seeddist, every_nsample, type, brushon);
			continue;
		}

		////else, we find a valid seed point, then grow a new streamline from this seed

		if(grow_a_tensorline_inReg(seedpts->seeds[i]->pos, seedpts->seeds[i]->triangle,
			percentage_dsep*dsep, discsize, sample_interval, loopdsep, dist2sing, 
			streamlinelength, type, brushon))
		{
			seedpts->seeds[i]->state = 1;

			if(evenstreamlines->isFull())
			{
				int oldnum = evenstreamlines->curMaxNumTrajs;
				if(!evenstreamlines->extend())
				{
					/*release the seed point list*/
					return;
				}

				/*allocate memeory for the elements*/
				for(i = oldnum; i < evenstreamlines->curMaxNumTrajs; i++)
				{
					evenstreamlines->trajs[i] = new Trajectory(i);
					if(evenstreamlines->trajs[i] == NULL)
					{
						return;
					}
				}

				/*increase the sample point lists as well*/
				SamplePtList **temp = samplepts;
				samplepts = (SamplePtList **)malloc(sizeof(SamplePtList*)*evenstreamlines->curMaxNumTrajs);
				if(samplepts == NULL)
				{
					/*report error*/
					if(samplepts == NULL)
					{
						return;
					}

					/*release the seed point list*/
					delete  seedpts;
					seedpts=NULL;
					return;
				}

				/*copy the information from the old list*/
				for(i = 0; i < oldnum; i++)
					samplepts[i] = temp[i];
				delete [] temp;

				/*allocate the pointers for the real samples*/
				for(i = oldnum; i < evenstreamlines->curMaxNumTrajs; i++)
				{
					samplepts[i] = new SamplePtList(100);
					if(samplepts[i] == NULL)
					{
						return;
					}

					for(int j = 0; j < samplepts[i]->curMaxNumSamplePts; j++)
					{
						samplepts[i]->samples[j] = (SamplePt *)malloc(sizeof(SamplePt));

						if(samplepts[i]->samples[j] == NULL)
						{
							return;
						}
					}
				}
				
				/*extend the list of seedposition_ineachtensorline*/
				seedposition_ineachtensorline=(int*)realloc(seedposition_ineachtensorline,
					sizeof(int)*evenstreamlines->curMaxNumTrajs);
				if(seedposition_ineachtensorline==NULL)
					return;
			}
		}

		else
		{
			seedpts->seeds[i]->state = 3;   //the seed rejected due to the short streamline
		}

	}

	/*release the seed point list*/
LL:	delete  seedpts;
	seedpts=NULL;

	/*release the sample point list*/
	for(i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
	{
		delete samplepts[i];
		samplepts[i] = NULL;
	}
	delete [] samplepts;
	samplepts = NULL;
	reset_placement_quad();
}

int EvenStreamlinePlace::trace_in_quad_inReg(int &face_id, double globalp[2], int type, 
					double dtest, double loopsep, double dist2sing, 
					double sample_interval, double discsize, int &flag)
{
	int i;
	double pre_point[2];

	double origin_dtest=dtest;
	double origin_dist2sing=dist2sing;
	double origin_loopsep=loopsep;
	
	/*  will this be a good solution? 1/9/2008 */
	if(!is_in_cell(face_id, globalp[0], globalp[1]))
	{
		face_id = get_cellID_givencoords(globalp[0], globalp[1]);
	}

	if(face_id < 0 || face_id>=quadmesh->nfaces)
		return -1;

	QuadCell *face = quadmesh->quadcells[face_id];

	QuadCell *pre_f = face;
	
	////Temporary curve point array

	CurvePoints *temp_point_list = (CurvePoints*) malloc(sizeof(CurvePoints) * 200);

	if(temp_point_list == NULL)
	{
		exit(-1);
	}

	int NumPoints = 0;

	/*   NOTE:  if this is a boundary cell, we need to be very careful   */
	bool boundarytest=false;
	if(quadmesh->quadcells[face_id]->OnBoundary)
		boundarytest=true;
	
	/*the tracing will be performed under the global frame*/
	globalface = face_id;

	pre_point[0] = globalp[0];
	pre_point[1] = globalp[1];

	////////////////////////////////////////////////////
    for(i = 0; i < 200; i++)
	{

		////2. if current point is inside current triangle
		if(is_in_cell(face_id, globalp[0], globalp[1]))
		{
			////store the point into the temp curve points list

			temp_point_list[NumPoints].gpx = globalp[0];
			temp_point_list[NumPoints].gpy = globalp[1];
			temp_point_list[NumPoints].triangleid = face->index;  
			NumPoints++;

			pre_point[0] = globalp[0];
			pre_point[1] = globalp[1];
			
			/*change to use other integration scheme 07/09/07*/
			//if(compute_next_pt_tensor_quad_global(pre_point, globalp, face_id))
			//if(get_nextpt_2ndeuler_ten_quad(pre_point, globalp, face_id, type))
			if(get_nextpt_RK23_ten_quad(pre_point, globalp, face_id, type))
			//if(get_nextpt_RK45_ten_quad(pre_point, globalp, face_id, type))
			{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = globalp[0] - pre_point[0];
				tenline_dir_global.entry[1] = globalp[1] - pre_point[1];

				/*      we need to combine the density map to change 
						the separation distance automatically
						11/25/2007
				*/
				if(sharedvars.CombinePopDensityOn)
				{
					double approx_den=cal_approx_density_at(globalp, face_id);
					dtest = origin_dtest/approx_den;
					dist2sing=origin_dist2sing/approx_den;
					loopsep=origin_loopsep/approx_den;
				}

				if(boundarytest&&!is_inregion(globalp[0], globalp[1]))
				{
					/*  if it reaches the region boundary  */
					flag = 1;
					if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
						(temp_point_list, NumPoints))
					{
						////Not enough memory
						flag = 4;
						free(temp_point_list);
						return face_id;
					}

					free(temp_point_list);
				    return face_id;
				}


				////using distance to judge whether a point close to current degenerate points
				if(close_to_degPt_except(globalp, face_id, -1, dist2sing, discsize))
				{
					flag = 1;

					//if(!StoreToGlobalList(temp_point_list, NumPoints))
					if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
						(temp_point_list, NumPoints))
					{
						////Not enough memory
						flag = 4;
						free(temp_point_list);
						return face_id;
					}

					free(temp_point_list);
				    return face_id;
				}

				////Judge whether it is too close to other existing streamlines
				if(evenstreamlines->ntrajs > 0 &&
					close_to_cur_streamline(globalp, face_id, &evenstreamlines->ntrajs, 1, dtest, discsize, 0)) //scale the separate distance
				{
					//if(!StoreToGlobalList(temp_point_list, NumPoints))
					if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
						(temp_point_list, NumPoints))
					{
						////Not enough memory
						flag = 4;
						free(temp_point_list);
						return face_id;
					}

					free(temp_point_list);
					//flag = 3;
					flag = 2;

					return face_id;
				}

				////We may also need to compare the current point with the sampling point on itself!
				if(close_to_cur_samplePt(globalp, face_id, samplepts[evenstreamlines->ntrajs]->samples,
					samplepts[evenstreamlines->ntrajs]->nsamples, loopsep, discsize, sample_interval)) //scale the separate distance
				{
					//if(!StoreToGlobalList(temp_point_list, NumPoints))
					if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
						(temp_point_list, NumPoints))
					{
						////Not enough memory
						flag = 4;
						free(temp_point_list);
						return face_id;
					}

					flag = 3; //form a loop! 4/18/06

					free(temp_point_list);
					return face_id;
				}

			}

			else{  ////the curve reach a singularity/degenerate point
				flag = 1;

				////Store the record into global line segment array
                
				//if(!StoreToGlobalList(temp_point_list, NumPoints))
				if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
					(temp_point_list, NumPoints))
				{
					////Not enough memory
					flag = 4;
					free(temp_point_list);
					return face_id;
				}

				free(temp_point_list);

				return face_id;
			}
		}



		////3. if the point is out of current cell
		else{

			/*!!!!!!need to judge which cell it will enter!!!!!*/
			int PassVertornot = 0;
			get_next_cell_2(face_id, pre_point, globalp, PassVertornot, type);
			
			if(PassVertornot>0)  /*cross a vertex*/
			{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = pre_point[0] - globalp[0];
				tenline_dir_global.entry[1] = pre_point[1] - globalp[1];

				/**/
				temp_point_list[NumPoints].gpx = globalp[0];
				temp_point_list[NumPoints].gpy = globalp[1];
				temp_point_list[NumPoints].triangleid = face_id/*face->index*/;  ////cause problem 05/25/05
				NumPoints++;

				temp_point_list[NumPoints].gpx = pre_point[0];
				temp_point_list[NumPoints].gpy = pre_point[1];
				temp_point_list[NumPoints].triangleid = face_id;  ////cause problem 05/25/05
				NumPoints++;
			}
			else{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = globalp[0] - pre_point[0];
				tenline_dir_global.entry[1] = globalp[1] - pre_point[1];

				////Add the intersection point to the temporary points' list
				temp_point_list[NumPoints].gpx = globalp[0];
				temp_point_list[NumPoints].gpy = globalp[1];
				temp_point_list[NumPoints].triangleid = face->index;  ////cause problem 05/25/05
				NumPoints++;
			}

			//if(globalp[1] > 2)
			//{
			//	int test=0;
			//}

			/*obtain the global direction of current tensor line 09/20/2007*/

			if(NumPoints > 1){
 				////Store the record into global line segment array
               //if(!StoreToGlobalList(temp_point_list, NumPoints))
				if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
					(temp_point_list, NumPoints))
			   {   ////Not enough memory
				   flag = 4;
				   free(temp_point_list);
				   return face_id;
			   }
			}

			free(temp_point_list);
			return face_id;
		}
	
	}

	if(NumPoints > 0)
		if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
						(temp_point_list, NumPoints))
						flag=4;

	free(temp_point_list);

	return face_id;
}

/*
grow a tensor line from a certain seed
*/
bool EvenStreamlinePlace::grow_a_tensorline_inReg(double seed_p[2], int triangle, double dtest, 
											double discsize, double Sample_interval, 
											double loopdsep, double dist2sing, 
											double streamlinelength, int type, bool brushon)
{
	int i;
	int flag = -1;

	int pre_face, cur_face;
	double globalp[3] = {0.};
	int cur_line = 0;
	double cur_length = 0;
	int movetonext = 0;
		
	if(triangle < 0 || triangle >= quadmesh->nfaces)
		return false;

	pre_face = cur_face = triangle;
	globalp[0] = seed_p[0];
	globalp[1] = seed_p[1];
	
	FILE *fp;
	
	icMatrix2x2 ten;
	double t[4]={0.};

	if(!is_in_reg_cell(triangle, seed_p[0], seed_p[1]))
	{
		triangle = get_cellID_givencoords(seed_p[0], seed_p[1]);
	}

	compute_tensor_at_quad(triangle, seed_p[0], seed_p[1], ten);

	double evalues[2] = {0.};
	icVector2 ev[2], startdir;
	cal_eigen_vector_sym(ten, ev);

	if(type == 1) /*use minor eigen vector field*/
		ev[0] = ev[1];

	tenline_dir_global = startdir = -ev[0];  /*obtain the major eigen vector*/


	evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs = 0;
	evenstreamlines->trajs[evenstreamlines->ntrajs]->closed=false;


	samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[0] = seed_p[0];
	samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[1] = seed_p[1];
	samplepts[evenstreamlines->ntrajs]->samples[0]->triangle = triangle;
	samplepts[evenstreamlines->ntrajs]->nsamples=1;
	
extern double hstep;

	hstep = quadmesh->xinterval/2.;
	predict_stepsize = quadmesh->xinterval/2.;
	euler_stepsize = quadmesh->xinterval/8.;


	cur_line = 0;
	cur_length = 0;

		/*record the postion of the original seed point in the tensor line 11/18/2007*/
	seedposition_ineachtensorline[evenstreamlines->ntrajs]=0;
	bool firsttri=true;
	//int start_tri=triangle;

	//////////////////////////////////////////////////////////////////////////


	////Backward tracing
	int NUMTRACETRIS = (int)sqrt((double)quadmesh->nfaces);

	for(i = 0; i < 3*NUMTRACETRIS; i++)
	{
		if(globalp[0]>quadmesh->xend-1.e-7||globalp[0]<quadmesh->xstart+1.e-7
			||globalp[1]>quadmesh->yend-1.e-7||globalp[1]<quadmesh->ystart+1.e-7)
		{
			if(globalp[0]>quadmesh->xend-1.e-7) globalp[0]=quadmesh->xend;
			else if(globalp[0]<quadmesh->xstart+1.e-7) globalp[0]=quadmesh->xstart;
			if(globalp[1]>quadmesh->yend-1.e-7) globalp[1]=quadmesh->yend;
			else if(globalp[1]<quadmesh->ystart+1.e-7) globalp[1]=quadmesh->ystart;
			break;
		}

		////The cell does not exist. Something is wrong!
		if(cur_face < 0 || cur_face >= quadmesh->nfaces)
			break;

		if(cur_face != triangle)
			firsttri=false;

		/*for brush interface 10/10/2007*/
		if(!quadmesh->quadcells[cur_face]->in_region)
			break;

		if(!firsttri && quadmesh->quadcells[cur_face]->OnBoundary)
			break;

		/*for load the map 10/24/2007*/
		if(is_not_inland(cur_face))
			break;
		if(is_inveg_cell_weak(cur_face))
			break;
		
		/*it reaches any boundary, we should stop as well*/
		//if(globalp[0]>quadmesh->xend-1.e-8||globalp[0]<quadmesh->xstart+1.e-8
		//	||globalp[1]>quadmesh->yend-1.e-8||globalp[1]<quadmesh->ystart+1.e-8)
		//	break;
	
		pre_face = cur_face;
        cur_face = trace_in_quad_inReg(cur_face, globalp, type, dtest, loopdsep, dist2sing,
			Sample_interval, discsize, flag);

		//if this is the first line segment, there may be a case that the start point is on the edge !!!
		if(pre_face == triangle && evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
		{
			cur_length = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].length;
		}

		
		////We need to select the sampling points from current trajectory  4/22/06
		cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
			samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
		

		if(flag == 1 || flag == 4 || flag == 2 /*|| flag == 3  || pre_face == cur_face*/ ) 
		{
			break;
		}

		if(flag == 3 && is_valid_loop(evenstreamlines->ntrajs, loopdsep) 
			&& sharedvars.CloseLoopOn) /*the tensor line forms a loop here*/
		{
			/*connect the whole tensor line*/
			int cur_line_index=evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs;
			if(cur_line_index<=0) break;

			
			double startp[2], endp[2];
			Trajectory *temp=new Trajectory(-1);
			startp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[0];
			startp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[1];
			int cell1=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].Triangle_ID;
			endp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
			endp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
			int cell2=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;

			temp->nlinesegs=0;

			get_linesegs_anytwopts(startp,cell1,  endp,cell2, temp, 0, 10);
			
			evenstreamlines->trajs[evenstreamlines->ntrajs]->add_last_nlines(
				temp->linesegs, temp->nlinesegs);
			delete temp;

			evenstreamlines->trajs[evenstreamlines->ntrajs]->closed=true;
		   
			cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
				samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);

			goto L10;
		}

	}

 	////Reverse the order of the obtained line segments
	if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
	{
		reverse_streamline(evenstreamlines->ntrajs);

		//////Resample the streamline after the reversion
		cur_line = 0;
		cur_length = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].length;
		samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[0] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
		samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[1] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
		samplepts[evenstreamlines->ntrajs]->samples[0]->triangle = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;
		samplepts[evenstreamlines->ntrajs]->nsamples = 1;

		movetonext = 0;
	    cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
			samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);

		/*record the postion of the original seed point in the tensor line 11/18/2007*/
		seedposition_ineachtensorline[evenstreamlines->ntrajs]=evenstreamlines->trajs[
			evenstreamlines->ntrajs]->nlinesegs-1;
	}
	

	//////////Forward tracing
	pre_face = cur_face = triangle;
	globalp[0] = seed_p[0];
	globalp[1] = seed_p[1];
	flag = -1;
	tenline_dir_global = -startdir;
	firsttri=true;

	for(i = 0; i < 3*NUMTRACETRIS; i++)
	{

		if(globalp[0]>quadmesh->xend-1.e-7||globalp[0]<quadmesh->xstart+1.e-7
			||globalp[1]>quadmesh->yend-1.e-7||globalp[1]<quadmesh->ystart+1.e-7)
		{
			if(globalp[0]>quadmesh->xend-1.e-7) globalp[0]=quadmesh->xend;
			else if(globalp[0]<quadmesh->xstart+1.e-7) globalp[0]=quadmesh->xstart;
			if(globalp[1]>quadmesh->yend-1.e-7) globalp[1]=quadmesh->yend;
			else if(globalp[1]<quadmesh->ystart+1.e-7) globalp[1]=quadmesh->ystart;
			break;
		}

		if(cur_face < 0 || cur_face >= quadmesh->nfaces)
			break;

		if(cur_face != triangle)
			firsttri=false;
		
		/*for brush interface 10/10/2007*/
		if(!quadmesh->quadcells[cur_face]->in_region)
			break;
			
		if(!firsttri && quadmesh->quadcells[cur_face]->OnBoundary)
			break;
	
		/*for load the map 10/24/2007*/
		if(is_not_inland(cur_face))
			break;
		if(is_inveg_cell_weak(cur_face))
			break;
		
		/*it reaches any boundary, we should stop as well*/
		//if(globalp[0]>quadmesh->xend-1.e-8||globalp[0]<quadmesh->xstart+1.e-8
		//	||globalp[1]>quadmesh->yend-1.e-8||globalp[1]<quadmesh->ystart+1.e-8)
		//	break;
		

		pre_face = cur_face;
        cur_face = trace_in_quad_inReg(cur_face, globalp, type, dtest, loopdsep, dist2sing,
			Sample_interval, discsize, flag);
		

		////We need to select the sampling points from current trajectory  3/9/06
	    cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
			samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
		
		
	    if(flag == 1 || flag == 4 || flag == 2 /*|| flag == 3 || pre_face == cur_face*/ ) 
		{
			break;
		}

		if(flag == 3 && is_valid_loop(evenstreamlines->ntrajs, loopdsep) 
			&& sharedvars.CloseLoopOn) /*the tensor line forms a loop here*/
		{
			/*connect the whole tensor line*/
			/*we need to extend if not enough memory*/
			int cur_line_index=evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs;
			if(cur_line_index<=0) break;

			double startp[2], endp[2];
			Trajectory *temp=new Trajectory(-1);
			startp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[0];
			startp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[1];
			int cell1=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].Triangle_ID;
			endp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
			endp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
			int cell2=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;

			temp->nlinesegs=0;

			get_linesegs_anytwopts(startp,cell1,  endp,cell2, temp, 0, 10);
			
			evenstreamlines->trajs[evenstreamlines->ntrajs]->add_last_nlines(
				temp->linesegs, temp->nlinesegs);
			delete temp;
		   
			cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
				samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
			
			evenstreamlines->trajs[evenstreamlines->ntrajs]->closed=true;

			break;
		}
	}
	

    
L10:	if((/*evenstreamlines->ntrajs == 0 ||*/
		evenstreamlines->trajs[evenstreamlines->ntrajs]->get_length() > streamlinelength)
		&& evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
	{
		update_samples_in_cell(evenstreamlines->ntrajs, samplepts[evenstreamlines->ntrajs]->samples,
			samplepts[evenstreamlines->ntrajs]->nsamples);

		//ndisplay_trajs = evenstreamlines->ntrajs-1;

		/*we need to update the line information for each cell of the mesh 10/02/2007*/
		Trajectory *cur_traj = evenstreamlines->trajs[evenstreamlines->ntrajs];
		int pre_cell = cur_traj->linesegs[0].Triangle_ID;
		int start, end;
		start = end = 0;
		for(i=1; i<cur_traj->nlinesegs; i++)
		{
			if(cur_traj->linesegs[i].Triangle_ID == pre_cell)
			{
				end = i;
				continue;
			}
			else{
				LinesInOneCell *newline = (LinesInOneCell*)malloc(sizeof(LinesInOneCell));
				newline->whichtraj = evenstreamlines->ntrajs;
				newline->start = start;
				newline->end = end;
				QuadCell *curc = quadmesh->quadcells[pre_cell];

				//if it is major field tracing, add it to the major line info
				if(!majororminor)
				{
					if(curc->majorlines == NULL)
						curc->majorlines = new LineInfo(1);
					curc->majorlines->addNew(newline);
					curc->hasmajor=true;
				}

				//else, add to the minor line info
				else
				{
					if(curc->minorlines == NULL)
						curc->minorlines = new LineInfo(1);
					curc->minorlines->addNew(newline);
					curc->hasminor=true;
				}
				start = end = i;
				pre_cell = cur_traj->linesegs[i].Triangle_ID;
			}
		}

		/*we need to handle the last line segments*/
		start = end = cur_traj->nlinesegs-1;
		pre_cell = cur_traj->linesegs[start].Triangle_ID;
		if(pre_cell>=0&&pre_cell<quadmesh->nfaces)
		{
			for(i=cur_traj->nlinesegs-2; i>=0; i--)
			{
				if(cur_traj->linesegs[i].Triangle_ID == pre_cell)
				{
					start = i;
					continue;
				}
				else{
					LinesInOneCell *newline = (LinesInOneCell*)malloc(sizeof(LinesInOneCell));
					newline->whichtraj = evenstreamlines->ntrajs;
					newline->start = start;
					newline->end = end;
					QuadCell *curc = quadmesh->quadcells[pre_cell];

					//if it is major field tracing, add it to the major line info
					if(!majororminor)
					{
						if(curc->majorlines == NULL)
							curc->majorlines = new LineInfo(1);
						curc->majorlines->addNew(newline);
						curc->hasmajor=true;
					}

					//else, add to the minor line info
					else
					{
						if(curc->minorlines == NULL)
							curc->minorlines = new LineInfo(1);
						curc->minorlines->addNew(newline);
						curc->hasminor=true;
					}
					start = end = i;
					pre_cell = cur_traj->linesegs[i].Triangle_ID;
					break;
				}
			}
		}
		else
		{
			/*  we need to handle this!  */
			int test=0;
			cur_traj->nlinesegs--;
			while(cur_traj->linesegs[cur_traj->nlinesegs-1].Triangle_ID<0
				||cur_traj->linesegs[cur_traj->nlinesegs-1].Triangle_ID>=quadmesh->nfaces)
				cur_traj->nlinesegs--;
		}

		/*set the corresponding road type as default one 10/03/2007*/
		if(majororminor)
			evenstreamlines->trajs[evenstreamlines->ntrajs]->roadtype = MINOR;
		else
			evenstreamlines->trajs[evenstreamlines->ntrajs]->roadtype = MINOR;
		
		evenstreamlines->trajs[evenstreamlines->ntrajs]->traj_len=
			evenstreamlines->trajs[evenstreamlines->ntrajs]->get_length();

		/*reorder the line segments of the streamline 10/13/2007*/
		evenstreamlines->ntrajs ++;
		nnewlines_inReg++;

		return true;
	}

	return false;

}




bool EvenStreamlinePlace::extend_trajList()
{
	int i;
	int oldnum = evenstreamlines->curMaxNumTrajs;
	if(!evenstreamlines->extend())
	{
		/*release the seed point list*/
		return false;
	}

	/*allocate memeory for the elements*/
	for(i = oldnum; i < evenstreamlines->curMaxNumTrajs; i++)
	{
		evenstreamlines->trajs[i] = new Trajectory(i);
		if(evenstreamlines->trajs[i] == NULL)
		{
			return false;
		}
	}

	/*increase the sample point lists as well*/
	SamplePtList **temp = samplepts;
	samplepts = (SamplePtList **)malloc(sizeof(SamplePtList*)*evenstreamlines->curMaxNumTrajs);
	if(samplepts == NULL)
	{
		/*report error*/
		if(samplepts == NULL)
		{
			return false;
		}

		/*release the seed point list*/
		delete  seedpts;
		seedpts=NULL;
		return false;
	}

	/*copy the information from the old list*/
	for(i = 0; i < oldnum; i++)
		samplepts[i] = temp[i];
	delete [] temp;

	/*allocate the pointers for the real samples*/
	for(i = oldnum; i < evenstreamlines->curMaxNumTrajs; i++)
	{
		samplepts[i] = new SamplePtList(100);
		if(samplepts[i] == NULL)
		{
			return false;
		}

		for(int j = 0; j < samplepts[i]->curMaxNumSamplePts; j++)
		{
			samplepts[i]->samples[j] = (SamplePt *)malloc(sizeof(SamplePt));

			if(samplepts[i]->samples[j] == NULL)
			{
				return false;
			}
		}
	}
	
	/*extend the list of seedposition_ineachtensorline*/
	seedposition_ineachtensorline=(int*)realloc(seedposition_ineachtensorline,
		sizeof(int)*evenstreamlines->curMaxNumTrajs);
	if(seedposition_ineachtensorline==NULL)
		return false;

	return true;
}


/*
This routine implement the alternative placement of tensor lines
NOTE: each call, the routine will trace at most one tensor line
*/
void EvenStreamlinePlace::place_tensorlines_alternatively(
	int type, bool brushon, SeedList *preseeds)
{
	int i;
	int locate_preseed = 0;
	int begin_sample = 0;
    
	/*-----Set a set of default initial values for testing 5/2/06-----*/
	/*--- this setting should be provided by user ---*/

	if(type == 0)
		dsep = majorDensity*quadmesh->xinterval;    //using the radius of the object instead of the edge
	else
		dsep = minorDensity*quadmesh->xinterval;    //using the radius of the object instead of the edge

	//percentage_dsep = 0.8;
	if(which_level==1)
		percentage_dsep = 0.15;
	else
		percentage_dsep = 0.6;

	discsize = 2.;
	sample_interval = 0.25*dsep;
	every_nsample = 4;
	loopdsep = 0.4*quadmesh->xinterval; /*we now allow the loops to be closed*/
	dist2sing = 0.1*quadmesh->xinterval;

	streamlinelength = mintenline_length*quadmesh->xinterval;
	
	seeddist = 1.2*dsep;
	minstartdist = 0.9*dsep;

	euler_stepsize = quadmesh->xinterval/5.;
	predict_stepsize = 2*quadmesh->xinterval;


	/*******************************************************************/
	/*allocate space for seed points*/
	if(seedpts == NULL)
		seedpts = new SeedList();

	/*copy the seed points here*/
	seedpts->copy_otherseedList(preseeds);
}



/*
calculate a set of initial streamlines without considering the obtained separatrices and periodic orbits
*/
void EvenStreamlinePlace::cal_init_streamlines(int num, double streamlinelength, double dtest, double discsize,
				double Sample_interval, double loopdsep, double dist2sing, int type, bool brushon)
{
	double begin_p[2] = {0.};
	int pre_triangle = 0;
	int triangle;

	int cur_traj_index;

	//triangle = 0;
	//triangle = 1;  /*need to avoid the triangle containing fixed point*/
	//triangle = (rand())%quadmesh->nfaces;  //we should use farthest scheme to generate
	if(!brushon)
		triangle = (int)(quadmesh->nfaces/2.);
		//triangle = (int)(quadmesh->nfaces/3.);
	else
		triangle = out_pts[0].cellid;

	/////////////////////////

	while(1)
	{
		if(evenstreamlines->ntrajs > 0)
			triangle = (pre_triangle + rand())%quadmesh->nfaces;  //we should use farthest scheme to generate

		cur_traj_index = evenstreamlines->ntrajs;

		begin_p[0] = quadmesh->quadcells[triangle]->x_start_coord+quadmesh->xinterval/2.;
		begin_p[1] = quadmesh->quadcells[triangle]->y_start_coord+quadmesh->yinterval/2.;


		if(evenstreamlines->ntrajs > 0 && close_to_cur_streamline(begin_p, triangle, &cur_traj_index, 1, 10*dtest, discsize, 0))
		{
			pre_triangle = triangle;
			continue;
		}


		if(grow_a_tensorline(begin_p, triangle,
			dtest, discsize, Sample_interval, loopdsep, dist2sing, streamlinelength, type, brushon))
		{
			pre_triangle = triangle;
		}
		

		if(evenstreamlines->ntrajs >= num)
			return;
	}
}


/*
calculate a set of initial streamlines according to user specified
seeds
*/
void EvenStreamlinePlace::cal_init_streamlines(double streamlinelength, double dtest, double discsize,
				double Sample_interval, double loopdsep, double dist2sing, int type, bool brushon,
				SeedList *userspecifiedseeds)
{
	double begin_p[2] = {0.};
	int pre_triangle = 0;
	int triangle;

	int cur_traj_index;

	/////////////////////////
	int i;

	for(i=0;i<userspecifiedseeds->nseeds; i++)
	{
		cur_traj_index = evenstreamlines->ntrajs;

		begin_p[0] = userspecifiedseeds->seeds[i]->pos[0];
		begin_p[1] = userspecifiedseeds->seeds[i]->pos[1];

		triangle = userspecifiedseeds->seeds[i]->triangle;


		if(evenstreamlines->ntrajs > 0 && close_to_cur_streamline(begin_p, triangle, 
			&cur_traj_index, 1, 3*dtest, discsize, 0))
		{
			continue;
		}

		grow_a_tensorline(begin_p, triangle,
			dtest, discsize, Sample_interval, loopdsep, dist2sing, streamlinelength, type, brushon);
	}
}



bool testflag = false;

bool is_not_inland(int cellid)
{
	if(cellid<0 || cellid>=quadmesh->nfaces) return false;
	QuadCell *face = quadmesh->quadcells[cellid];
	QuadVertex *v;

	int i;

	for(i=0; i<face->nverts; i++)
	{
		v = quadmesh->quad_verts[face->verts[i]];

		if(!v->inland)
			return true;
	}
	return false;
}

/*
   We define that the cell is not an "inland" cell iff it contains more than
   one not inland vertices
*/
bool is_not_inland_cell_weak(int cellid)
{
	if(cellid<0 || cellid>=quadmesh->nfaces) return false;

	QuadCell *face = quadmesh->quadcells[cellid];
	QuadVertex *v;

	int i;
	int NnotinlandVerts=0;

	for(i=0; i<face->nverts; i++)
	{
		v = quadmesh->quad_verts[face->verts[i]];

		if(!v->inland)
		{
			NnotinlandVerts++;
			if(NnotinlandVerts>=2)
				return true;
		}
	}
	return false;
}



bool is_inveg(int cellid)
{
	QuadCell *face = quadmesh->quadcells[cellid];
	QuadVertex *v;

	int i;

	for(i=0; i<face->nverts; i++)
	{
		v = quadmesh->quad_verts[face->verts[i]];

		if(v->inveg)
			return true;
	}
	return false;
}

/*
   We define that the cell is not an "inland" cell iff it contains more than
   one not inland vertices
*/
bool is_inveg_cell_weak(int cellid)
{
	QuadCell *face = quadmesh->quadcells[cellid];
	QuadVertex *v;

	int i;
	int NinvegVerts=0;

	for(i=0; i<face->nverts; i++)
	{
		v = quadmesh->quad_verts[face->verts[i]];

		if(v->inveg)
		{
			NinvegVerts++;
			if(NinvegVerts>=2)
				return true;
		}
	}
	return false;
}


/*
This routine should be called after rearrange the previously obtained
tensor lines
11/19/2007
*/

void EvenStreamlinePlace::rebuild_all_lineinfo()
{
	/*we first need to release the previous information list*/
	int i, j;
	QuadCell *face;
	for(i=0; i<quadmesh->nfaces;i++)
	{
		/**/
		face=quadmesh->quadcells[i];

		if(!majororminor)
		{
			if(face->majorlines!=NULL)
			{
				delete face->majorlines;
				face->majorlines=NULL;
			}
			face->hasmajor=false;
		}
		else
		{
			if(face->minorlines!=NULL)
			{
				//for(j=0;j<face->minorlines->nlines;j++)
				//{
				//	if(face->minorlines[j]!=NULL)
				//	{
				//		delete face->minorlines[j]);
				//		face->minorlines[j]=NULL;
				//	}
				//}
				//free(face->minorlines);
				delete face->minorlines;
				face->minorlines=NULL;
			}
			face->hasminor=false;
		}
	}

	/*rebuild the information list for each cell according to the current tensor lines*/
	for(i=0;i<evenstreamlines->ntrajs;i++)
	{
		rebuild_lineinfo_for_a_tensorline(i);
	}
}

void EvenStreamlinePlace::rebuild_lineinfo_for_a_tensorline(int id)
{
	/*this is not necessary right now*/
	//update_samples_in_cell(evenstreamlines->ntrajs, samplepts[evenstreamlines->ntrajs]->samples,
	//	samplepts[evenstreamlines->ntrajs]->nsamples);

	/*we need to update the line information for each cell of the mesh 10/02/2007*/
	int i;
	Trajectory *cur_traj = evenstreamlines->trajs[id];
	int pre_cell = cur_traj->linesegs[0].Triangle_ID;
	int start, end;
	start = end = 0;
	for(i=1; i<cur_traj->nlinesegs; i++)
	{
		if(cur_traj->linesegs[i].Triangle_ID == pre_cell)
		{
			end = i;
			continue;
		}
		else{
			LinesInOneCell *newline = (LinesInOneCell*)malloc(sizeof(LinesInOneCell));
			newline->whichtraj = id;
			newline->start = start;
			newline->end = end;
			QuadCell *curc = quadmesh->quadcells[pre_cell];

			//if it is major field tracing, add it to the major line info
			if(!majororminor)
			{
				if(curc->majorlines == NULL)
					curc->majorlines = new LineInfo(1);
				curc->majorlines->addNew(newline);
				curc->hasmajor=true;
			}

			//else, add to the minor line info
			else
			{
				if(curc->minorlines == NULL)
					curc->minorlines = new LineInfo(1);
				curc->minorlines->addNew(newline);
				curc->hasminor=true;
			}
			start = end = i;
			pre_cell = cur_traj->linesegs[i].Triangle_ID;
		}
	}

	/*we need to handle the last line segments*/
	start = end = cur_traj->nlinesegs-1;
	pre_cell = cur_traj->linesegs[start].Triangle_ID;
	if(pre_cell>=0&&pre_cell<quadmesh->nfaces)
	{
		for(i=cur_traj->nlinesegs-2; i>=0; i--)
		{
			if(cur_traj->linesegs[i].Triangle_ID == pre_cell)
			{
				start = i;
				continue;
			}
			else{
				LinesInOneCell *newline = (LinesInOneCell*)malloc(sizeof(LinesInOneCell));
				newline->whichtraj = id;
				newline->start = start;
				newline->end = end;
				QuadCell *curc = quadmesh->quadcells[pre_cell];

				//if it is major field tracing, add it to the major line info
				if(!majororminor)
				{
					if(curc->majorlines == NULL)
						curc->majorlines = new LineInfo(1);
					if(!curc->majorlines->is_repeated(newline)){
						curc->majorlines->addNew(newline);
						curc->hasmajor=true;
					}
					else
						free(newline);

				}

				//else, add to the minor line info
				else
				{
					if(curc->minorlines == NULL)
						curc->minorlines = new LineInfo(1);
					if(!curc->minorlines->is_repeated(newline)){
						curc->minorlines->addNew(newline);
						curc->hasminor=true;
					}
					else
						free(newline);
				}
				start = end = i;
				pre_cell = cur_traj->linesegs[i].Triangle_ID;
				break;
			}
		}
	}
	else
	{
		/*  we need to handle this!  */
		int test=0;
		cur_traj->nlinesegs--;
		while(cur_traj->linesegs[cur_traj->nlinesegs-1].Triangle_ID<0
			||cur_traj->linesegs[cur_traj->nlinesegs-1].Triangle_ID>=quadmesh->nfaces)
			cur_traj->nlinesegs--;
	}

	/*set the corresponding road type as default one 10/03/2007*/
	if(majororminor)
		evenstreamlines->trajs[evenstreamlines->ntrajs]->roadtype = MINOR;
	else
		evenstreamlines->trajs[evenstreamlines->ntrajs]->roadtype = MINOR;
}


bool EvenStreamlinePlace::is_valid_loop(int trajid, double sep)
{
	Trajectory *traj=evenstreamlines->trajs[trajid];
	icVector2 startvec, endvec, dist;
	startvec.entry[0]=traj->linesegs[0].gstart[0]-traj->linesegs[0].gend[0];
	startvec.entry[1]=traj->linesegs[0].gstart[1]-traj->linesegs[0].gend[1];
	endvec.entry[0]=traj->linesegs[traj->nlinesegs-1].gend[0]
	-traj->linesegs[traj->nlinesegs-1].gstart[0];
	endvec.entry[1]=traj->linesegs[traj->nlinesegs-1].gend[1]
	-traj->linesegs[traj->nlinesegs-1].gstart[1];

	dist.entry[0]=traj->linesegs[traj->nlinesegs-1].gend[0]-traj->linesegs[0].gstart[0];
	dist.entry[1]=traj->linesegs[traj->nlinesegs-1].gend[1]-traj->linesegs[0].gstart[1];

	if(length(dist)<=sep+1.e-2)
		return true;
	return false;
}


/*
grow a tensor line from a certain seed
*/
bool EvenStreamlinePlace::grow_a_tensorline(double seed_p[2], int triangle, double dtest, 
											double discsize, double Sample_interval, 
											double loopdsep, double dist2sing, 
											double streamlinelength, int type, bool brushon)
{
	int i;
	int flag = -1;

	int pre_face, cur_face;
	double globalp[3] = {0.};
	int cur_line = 0;
	double cur_length = 0;
	int movetonext = 0;
		

	pre_face = cur_face = triangle;
	globalp[0] = seed_p[0];
	globalp[1] = seed_p[1];
	
	FILE *fp;
	
	icMatrix2x2 ten;
	double t[4]={0.};

	if(!/*is_in_reg_cell*/is_in_cell(triangle, seed_p[0], seed_p[1]))
	{
		triangle = get_cellID_givencoords(seed_p[0], seed_p[1]);
	}

	if(triangle < 0 || triangle >= quadmesh->nfaces)
		return false;
	
	compute_tensor_at_quad(triangle, seed_p[0], seed_p[1], ten);

	double evalues[2] = {0.};
	icVector2 ev[2], startdir;
	cal_eigen_vector_sym(ten, ev);

	if(type == 1) /*use minor eigen vector field*/
		ev[0] = ev[1];

	tenline_dir_global = startdir = -ev[0];  /*obtain the major eigen vector*/


	evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs = 0;
	evenstreamlines->trajs[evenstreamlines->ntrajs]->closed=false;


	samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[0] = seed_p[0];
	samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[1] = seed_p[1];
	samplepts[evenstreamlines->ntrajs]->samples[0]->triangle = triangle;
	samplepts[evenstreamlines->ntrajs]->nsamples=1;
	
extern double hstep;

	hstep = quadmesh->xinterval/2.;
	predict_stepsize = quadmesh->xinterval/2.;
	euler_stepsize = quadmesh->xinterval/2.;


	cur_line = 0;
	cur_length = 0;

		/*record the postion of the original seed point in the tensor line 11/18/2007*/
	seedposition_ineachtensorline[evenstreamlines->ntrajs]=0;

	//////////////////////////////////////////////////////////////////////////

	////Backward tracing
	int NUMTRACETRIS = (int)sqrt((double)quadmesh->nfaces);

	for(i = 0; i < 3*NUMTRACETRIS; i++)
	{

		////The cell does not exist. Something is wrong!
		if(cur_face < 0 || cur_face >= quadmesh->nfaces)
		{
			break;
		}

		/*it reaches any boundary, we should stop as well*/
		if(globalp[0]>quadmesh->xend-1.e-7||globalp[0]<quadmesh->xstart+1.e-7
			||globalp[1]>quadmesh->yend-1.e-7||globalp[1]<quadmesh->ystart+1.e-7)
		{
			if(globalp[0]>quadmesh->xend-1.e-7) globalp[0]=quadmesh->xend;
			else if(globalp[0]<quadmesh->xstart+1.e-7) globalp[0]=quadmesh->xstart;
			if(globalp[1]>quadmesh->yend-1.e-7) globalp[1]=quadmesh->yend;
			else if(globalp[1]<quadmesh->ystart+1.e-7) globalp[1]=quadmesh->ystart;
			break;
		}

		
		/*for brush interface 10/10/2007*/
		if(brushon && !quadmesh->quadcells[cur_face]->OnBoundary)
			break;

		/*for load the map 10/24/2007*/
		if(is_not_inland(cur_face))
			break;
		
		if(is_inveg_cell_weak(cur_face))
			break;

	
		pre_face = cur_face;
        cur_face = trace_in_quad_even(cur_face, globalp, type, dtest, loopdsep, dist2sing,
			Sample_interval, discsize, flag);

		//if this is the first line segment, there may be a case that the start point is on the edge !!!
		if(pre_face == triangle && evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
		{
			cur_length = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].length;
		}

		
		////We need to select the sampling points from current trajectory  4/22/06
		cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
			samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
		

		if(flag == 1 || flag == 4 || flag == 2 /*|| flag == 3  || pre_face == cur_face*/ ) 
		{
			if(sharedvars.RemDeadEndsTraceOn)
			{
				if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
				{
					globalp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs
						[evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs-1].gend[0];
					globalp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs
						[evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs-1].gend[1];
					cur_face=get_cellID_givencoords(globalp[0], globalp[1]);
				}
				Trajectory *temp_traj=new Trajectory(-1);
				double dist_did=0;
				if(trace_to_comp_nearby_minRoad(globalp, cur_face, majororminor, 20,
					dsep, dsep/2., dist_did, temp_traj))
				{
					evenstreamlines->trajs[evenstreamlines->ntrajs]->add_last_nlines(
						temp_traj->linesegs, temp_traj->nlinesegs);
					cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
						samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
				}
				delete temp_traj;
			}
			break;
		}

		if(flag == 3 && is_valid_loop(evenstreamlines->ntrajs, loopdsep)
			&& sharedvars.CloseLoopOn) /*the tensor line forms a loop here*/
		{
			/*connect the whole tensor line*/
			int cur_line_index=evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs;
			if(cur_line_index<=0) break;

			/*we need to extend if not enough memory*/
			//if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->curMaxNumLinesegs)
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->extend_line_segments(1);

			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gstart[0]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[0];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gstart[1]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[1];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].start[0]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].end[0];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].start[1]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].end[1];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].Triangle_ID=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].Triangle_ID;
			//
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gend[0]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gend[1]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].end[0]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].start[0];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].end[1]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].start[1];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].Triangle_ID=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;

			//evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs++;
			
			double startp[2], endp[2];
			Trajectory *temp=new Trajectory(-1);
			startp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[0];
			startp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[1];
			int cell1=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].Triangle_ID;
			endp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
			endp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
			int cell2=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;

			temp->nlinesegs=0;

			get_linesegs_anytwopts(startp,cell1,  endp,cell2, temp, 0, 5);
			
			evenstreamlines->trajs[evenstreamlines->ntrajs]->add_last_nlines(
				temp->linesegs, temp->nlinesegs);
			delete temp;

			evenstreamlines->trajs[evenstreamlines->ntrajs]->closed=true;
		   
			cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
				samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);

			goto L10;
		}

		//if(evenstreamlines->ntrajs==154 && !majororminor)
		//{
		//	fp = fopen("lineseg_mem_error.txt", "w");
		//	fprintf(fp, "# of line segments of traj %d is %d (backward).\n", evenstreamlines->ntrajs,
		//		evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs);
		//	fclose(fp);
		//}

	}

 	////Reverse the order of the obtained line segments
	if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
	{
		reverse_streamline(evenstreamlines->ntrajs);

		//////Resample the streamline after the reversion
		cur_line = 0;
		cur_length = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].length;
		samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[0] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
		samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[1] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
		samplepts[evenstreamlines->ntrajs]->samples[0]->triangle = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;
		samplepts[evenstreamlines->ntrajs]->nsamples = 1;

		movetonext = 0;
	    cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
			samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);

		/*record the postion of the original seed point in the tensor line 11/18/2007*/
		seedposition_ineachtensorline[evenstreamlines->ntrajs]=evenstreamlines->trajs[
			evenstreamlines->ntrajs]->nlinesegs-1;
	}
	

	//////////Forward tracing
	pre_face = cur_face = triangle;
	globalp[0] = seed_p[0];
	globalp[1] = seed_p[1];
	flag = -1;
	tenline_dir_global = -startdir;
		
	//if(evenstreamlines->ntrajs==51/*&&cur_face==7276*/)
	//{
	//	int test=0;
	//}

	for(i = 0; i < 3*NUMTRACETRIS; i++)
	{

		if(cur_face < 0 || cur_face >= quadmesh->nfaces)
			break;

		if(globalp[0]>quadmesh->xend-1.e-7||globalp[0]<quadmesh->xstart+1.e-7
			||globalp[1]>quadmesh->yend-1.e-7||globalp[1]<quadmesh->ystart+1.e-7)
		{
			if(globalp[0]>quadmesh->xend-1.e-7) globalp[0]=quadmesh->xend;
			else if(globalp[0]<quadmesh->xstart+1.e-7) globalp[0]=quadmesh->xstart;
			if(globalp[1]>quadmesh->yend-1.e-7) globalp[1]=quadmesh->yend;
			else if(globalp[1]<quadmesh->ystart+1.e-7) globalp[1]=quadmesh->ystart;
			break;
		}

		
		/*for brush interface 10/10/2007*/
		if(brushon && !quadmesh->quadcells[cur_face]->OnBoundary)
			break;
		
		/*for load the map 10/24/2007*/
		if(is_not_inland(cur_face))
			break;
		if(is_inveg_cell_weak(cur_face))
			break;
		
		///*it reaches any boundary, we should stop as well*/
		//if(globalp[0]>quadmesh->xend-1.e-8||globalp[0]<quadmesh->xstart+1.e-8
		//	||globalp[1]>quadmesh->yend-1.e-8||globalp[1]<quadmesh->ystart+1.e-8)
		//	break;
		

		pre_face = cur_face;
        cur_face = trace_in_quad_even(cur_face, globalp, type, dtest, loopdsep, dist2sing,
			Sample_interval, discsize, flag);
		

		////We need to select the sampling points from current trajectory  3/9/06
	    cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
			samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
		
		
	    if(flag == 1 || flag == 4 || flag == 2 /*|| flag == 3|| pre_face == cur_face*/ ) 
		{
			if(sharedvars.RemDeadEndsTraceOn)
			{
				if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
				{
					globalp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs
						[evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs-1].gend[0];
					globalp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs
						[evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs-1].gend[1];
					cur_face=get_cellID_givencoords(globalp[0], globalp[1]);
				}

				Trajectory *temp_traj=new Trajectory(-1);
				double dist_did=0;
				if(trace_to_comp_nearby_minRoad(globalp, cur_face, majororminor, 20,
					dsep, dsep/2., dist_did, temp_traj))
				{
					evenstreamlines->trajs[evenstreamlines->ntrajs]->add_last_nlines(
						temp_traj->linesegs, temp_traj->nlinesegs);
					cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
						samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
				}
				delete temp_traj;
			}

			break;
		}

		if(flag == 3 && is_valid_loop(evenstreamlines->ntrajs, loopdsep)
			&& sharedvars.CloseLoopOn) /*the tensor line forms a loop here*/
		{
			/*connect the whole tensor line*/
			/*we need to extend if not enough memory*/
			int cur_line_index=evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs;
			if(cur_line_index<=0) break;

			//if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->curMaxNumLinesegs)
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->extend_line_segments(5);

			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gstart[0]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[0];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gstart[1]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[1];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].start[0]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].end[0];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].start[1]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].end[1];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].Triangle_ID=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].Triangle_ID;
			//
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gend[0]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gend[1]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].end[0]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].start[0];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].end[1]=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].start[1];
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].Triangle_ID=
			//	evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;
			//
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs++;
			
			double startp[2], endp[2];
			Trajectory *temp=new Trajectory(-1);
			startp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[0];
			startp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[1];
			int cell1=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].Triangle_ID;
			endp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
			endp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
			int cell2=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;

			temp->nlinesegs=0;

			get_linesegs_anytwopts(startp,cell1,  endp,cell2, temp, 0, 5);
			
			evenstreamlines->trajs[evenstreamlines->ntrajs]->add_last_nlines(
				temp->linesegs, temp->nlinesegs);
			delete temp;
		   
			cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
				samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
			
			evenstreamlines->trajs[evenstreamlines->ntrajs]->closed=true;

			break;
		}
	
		//if(evenstreamlines->ntrajs==154 && !majororminor)
		//{
		//	fp = fopen("lineseg_mem_error.txt", "w");
		//	fprintf(fp, "# of line segments of traj %d is %d (forward).\n", evenstreamlines->ntrajs,
		//		evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs);
		//	fclose(fp);
		//}
	}
	

    
L10:	if((evenstreamlines->ntrajs == 0 ||
		evenstreamlines->trajs[evenstreamlines->ntrajs]->get_length() > streamlinelength)
		&& evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
	{
		update_samples_in_cell(evenstreamlines->ntrajs, samplepts[evenstreamlines->ntrajs]->samples,
			samplepts[evenstreamlines->ntrajs]->nsamples);

		//ndisplay_trajs = evenstreamlines->ntrajs-1;

		//if(evenstreamlines->ntrajs==154 && !majororminor)
		//{
		//	fp = fopen("lineseg_mem_error.txt", "w");
		//	fprintf(fp, "start putting information to cells :).");
		//	fclose(fp);
		//}

		/*we need to update the line information for each cell of the mesh 10/02/2007*/
		Trajectory *cur_traj = evenstreamlines->trajs[evenstreamlines->ntrajs];
		int pre_cell = cur_traj->linesegs[0].Triangle_ID;
		int start, end;
		start = end = 0;
		for(i=1; i<cur_traj->nlinesegs; i++)
		{
			if(cur_traj->linesegs[i].Triangle_ID == pre_cell)
			{
				end = i;
				continue;
			}
			else{
				LinesInOneCell *newline = (LinesInOneCell*)malloc(sizeof(LinesInOneCell));
				newline->whichtraj = evenstreamlines->ntrajs;
				newline->start = start;
				newline->end = end;
				QuadCell *curc = quadmesh->quadcells[pre_cell];

				//if it is major field tracing, add it to the major line info
				if(!majororminor)
				{
					if(curc->majorlines == NULL)
						curc->majorlines = new LineInfo(1);
					curc->majorlines->addNew(newline);
					curc->hasmajor=true;
				}

				//else, add to the minor line info
				else
				{
					if(curc->minorlines == NULL)
						curc->minorlines = new LineInfo(1);
					curc->minorlines->addNew(newline);
					curc->hasminor=true;
				}
				start = end = i;
				pre_cell = cur_traj->linesegs[i].Triangle_ID;
			}
		}

		/*we need to handle the last line segments*/
		start = end = cur_traj->nlinesegs-1;
		pre_cell = cur_traj->linesegs[start].Triangle_ID;
		if(pre_cell>=0&&pre_cell<quadmesh->nfaces)
		{
			for(i=cur_traj->nlinesegs-2; i>=0; i--)
			{
				if(cur_traj->linesegs[i].Triangle_ID == pre_cell)
				{
					start = i;
					continue;
				}
				else{
					LinesInOneCell *newline = (LinesInOneCell*)malloc(sizeof(LinesInOneCell));
					newline->whichtraj = evenstreamlines->ntrajs;
					newline->start = start;
					newline->end = end;
					QuadCell *curc = quadmesh->quadcells[pre_cell];

					//if it is major field tracing, add it to the major line info
					if(!majororminor)
					{
						if(curc->majorlines == NULL)
							curc->majorlines = new LineInfo(1);
						curc->majorlines->addNew(newline);
						curc->hasmajor=true;
					}

					//else, add to the minor line info
					else
					{
						if(curc->minorlines == NULL)
							curc->minorlines = new LineInfo(1);
						curc->minorlines->addNew(newline);
						curc->hasminor=true;
					}
					start = end = i;
					pre_cell = cur_traj->linesegs[i].Triangle_ID;
					break;
				}
			}
		}
		else
		{
			/*  we need to handle this!  */
			int test=0;
			cur_traj->nlinesegs--;
			while(cur_traj->linesegs[cur_traj->nlinesegs-1].Triangle_ID<0
				||cur_traj->linesegs[cur_traj->nlinesegs-1].Triangle_ID>=quadmesh->nfaces)
				cur_traj->nlinesegs--;
		}

		/*set the corresponding road type as default one 10/03/2007*/
		if(majororminor)
			evenstreamlines->trajs[evenstreamlines->ntrajs]->roadtype = MINOR;
		else
			//evenstreamlines->trajs[evenstreamlines->ntrajs]->roadtype = rand()%(FREEWAY+1);
			evenstreamlines->trajs[evenstreamlines->ntrajs]->roadtype = MINOR;

		/*reorder the line segments of the streamline 10/13/2007*/
		//reorder_streamline(evenstreamlines->ntrajs);

		evenstreamlines->trajs[evenstreamlines->ntrajs]->traj_len=
			evenstreamlines->trajs[evenstreamlines->ntrajs]->get_length();

		evenstreamlines->ntrajs ++;

		return true;
	}

	return false;

}


#include "ImgBoundaryExtract.h"
extern MapBoundaryList *mapboundarylist;

void follow_boundary(int boundID, int startpos, icVector2 orient, double followDist,
					 Trajectory *traj, double endp[2], int &endcell)
{
	/*  variables for pixel level judgement of "inland" or not  */
	double xstart=quadmesh->xstart;
	double ystart=quadmesh->ystart;
	double xrang=quadmesh->xend-quadmesh->xstart;
	double yrang=quadmesh->yend-quadmesh->ystart;
	double dx=xrang/511;
	double dy=yrang/511;

	MapBoundary *thebound=&mapboundarylist->mapboundarylist[boundID];
	icVector2 boundDir;
	if(startpos==thebound->nelems-1)
	{
		boundDir.entry[0]=thebound->pts[startpos]->x-thebound->pts[startpos-1]->x;
		boundDir.entry[1]=thebound->pts[startpos]->y-thebound->pts[startpos-1]->y;
	}
	else
	{
		boundDir.entry[0]=thebound->pts[startpos+1]->x-thebound->pts[startpos]->x;
		boundDir.entry[1]=thebound->pts[startpos+1]->y-thebound->pts[startpos]->y;
	}

	int i;
	icVector2 cur_lineDir;
	if(dot(boundDir,orient)>=0)
	{
		/*   we follow the boundary direction   */

		/*   we may still need to check whether it is too close to other existing 
		     major roads
		*/
		for(i=startpos;i<thebound->nelems-1;i++)
		{
			double p1[2],p2[2];
			p1[0]=thebound->pts[i]->x;
			p1[1]=thebound->pts[i]->y;
			int cell1=get_cellID_givencoords(p1[0], p1[1]);
			p2[0]=thebound->pts[i+1]->x;
			p2[1]=thebound->pts[i+1]->y;
			int cell2=get_cellID_givencoords(p2[0], p2[1]);
			if(cell1<0||cell1>=quadmesh->nfaces
				||cell2<0||cell2>=quadmesh->nfaces)
				return;

			/* judge whether it is too close to existing major road  
			   it it is, connect them! 1/1/2008
			*/

			cur_lineDir.entry[0]=p2[0]-p1[0];
			cur_lineDir.entry[1]=p2[1]-p1[1];

			if(dot(cur_lineDir, tenline_dir_global)<0)
			{
				/*  need to move away from the water!  */
				normalize(cur_lineDir);
				endp[0]=p1[0]+2*quadmesh->xinterval*cur_lineDir.entry[0];
				endp[1]=p1[1]+2*quadmesh->xinterval*cur_lineDir.entry[1];
				endcell=get_cellID_givencoords(endp[0], endp[1]);
				//if(is_not_inland(endcell))
				if(!is_inland_pixel(endp[0],endp[1], xstart, xrang, ystart, yrang, dx, dy,
						fittedmap1, 512))
				{
					endp[0]=p1[0]-2*quadmesh->xinterval*cur_lineDir.entry[0];
					endp[1]=p1[1]-2*quadmesh->xinterval*cur_lineDir.entry[1];
					endcell=get_cellID_givencoords(endp[0], endp[1]);

					//if(is_not_inland(endcell))
					if(!is_inland_pixel(endp[0],endp[1], xstart, xrang, ystart, yrang, dx, dy,
							fittedmap1, 512))
						return;

					get_linesegs_anytwopts(p1, cell1,  endp, endcell, traj, 0, 10);
				}
				else
				{
					get_linesegs_anytwopts(p1, cell1,  endp, endcell, traj, 0, 10);
				}
				return;
			}

			get_linesegs_anytwopts(p1, cell1,  p2, cell2, traj, 0, 10);
		}
	}

	else
	{
		/*   follow the inversed boundary direction  */
		for(i=startpos;i>=1;i--)
		{
			double p1[2],p2[2];
			p1[0]=thebound->pts[i]->x;
			p1[1]=thebound->pts[i]->y;
			int cell1=get_cellID_givencoords(p1[0], p1[1]);
			p2[0]=thebound->pts[i-1]->x;
			p2[1]=thebound->pts[i-1]->y;
			int cell2=get_cellID_givencoords(p2[0], p2[1]);
			if(cell1<0||cell1>=quadmesh->nfaces
				||cell2<0||cell2>=quadmesh->nfaces)
				return;
			
			/* judge whether it is too close to existing major road  
			   it it is, connect them! 1/1/2008
			*/

			cur_lineDir.entry[0]=p2[0]-p1[0];
			cur_lineDir.entry[1]=p2[1]-p1[1];

			if(dot(cur_lineDir, tenline_dir_global)<0)
			{
				/*  need to move away from the water!  */
				normalize(cur_lineDir);
				endp[0]=p1[0]+2*quadmesh->xinterval*cur_lineDir.entry[0];
				endp[1]=p1[1]+2*quadmesh->xinterval*cur_lineDir.entry[1];
				endcell=get_cellID_givencoords(endp[0], endp[1]);
				//if(is_not_inland(endcell))
				if(!is_inland_pixel(endp[0],endp[1], xstart, xrang, ystart, yrang, dx, dy,
						fittedmap1, 512))
				{
					endp[0]=p1[0]-2*quadmesh->xinterval*cur_lineDir.entry[0];
					endp[1]=p1[1]-2*quadmesh->xinterval*cur_lineDir.entry[1];
					endcell=get_cellID_givencoords(endp[0], endp[1]);
					//if(is_not_inland(endcell))
					//	return;
					get_linesegs_anytwopts(p1, cell1,  endp, endcell, traj, 0, 10);
				}
				else
				{
					get_linesegs_anytwopts(p1, cell1,  endp, endcell, traj, 0, 10);
				}
				return;
			}

			get_linesegs_anytwopts(p1, cell1,  p2, cell2, traj, 0, 10);
		}
	}
}

	//double xstart=quadmesh->xstart;
	//double ystart=quadmesh->ystart;
	//double xrang=quadmesh->xend-quadmesh->xstart;
	//double yrang=quadmesh->yend-quadmesh->ystart;
	//double dx=xrang/511;
	//double dy=yrang/511;

	//tp[0]=sx+2.*disc_radius*norm.entry[0];
	//tp[1]=sy+2.*disc_radius*norm.entry[1];
	//
	//if(is_inland_pixel(tp[0],tp[1], xstart, xrang, ystart, yrang, dx, dy,
	//		fittedmap1, 512))

/*
    Deal with the major roads closing to the water region
*/

void EvenStreamlinePlace::trace_water_region(double pt[2], int &cell, double Sample_interval,
						double loopdsep, double min_waterwidth, double cross_angle,
						double dist_follow_bound, int type)
{
	/*   first, obtain the boundary (index) crossing the cell (or nearby)   */
	double xstart=quadmesh->xstart;
	double ystart=quadmesh->ystart;
	double xrang=quadmesh->xend-quadmesh->xstart;
	double yrang=quadmesh->yend-quadmesh->ystart;
	double dx=xrang/511;
	double dy=yrang/511;

	int closest_bound;
	double intersect[2];
	icVector2 appro_bound_dir, norm;
	int startpos;

	/*  if we can not find a boundary across this cell, return now  */
	if(!find_closest_bound(pt, cell, closest_bound, intersect, startpos, appro_bound_dir))
		return;

	/*  if we find the boundary and the intersection, we need to evaluate the tangent
	    direction along the boundary and obtain its normal direction 
	*/
	normalize(appro_bound_dir);
	if(dot(appro_bound_dir, tenline_dir_global)<0)
		appro_bound_dir=-appro_bound_dir;

	norm.entry[0]=-appro_bound_dir.entry[1];
	norm.entry[1]=appro_bound_dir.entry[0];

	icVector2 traceDir=tenline_dir_global;

	normalize(traceDir);
	double ang_traceDir_boundDir=acos(dot(norm, traceDir));

	/*  go along the normal direction for min_waterwidth/2*/
	double ce[2];
	//ce[0]=pt[0]+min_waterwidth/2.*norm.entry[0];
	//ce[1]=pt[1]+min_waterwidth/2.*norm.entry[1];
	int ce_cell/*=get_cellID_givencoords(ce[0], ce[1])*/;

	///*  obtain a disc center at "ce"  */
	//trianglelist = new DynList_Int();

	//trianglelist->nelems = 0;
	//cal_euclidean_dist_2(ce_cell, ce, min_waterwidth/2., 1.2, trianglelist);

	//bool cross_river=false;
	//int i,j;
	//QuadCell *face;
	//QuadVertex *v;
	//icVector2 pos_dir;
	//for(i=0;i<trianglelist->nelems;i++)
	//{
	//	//face=quadmesh->quadcells[trianglelist->elems[i]];
	//	//for(j=0;j<face->nverts;j++)
	//	//{
	//	//	v=quadmesh->quad_verts[face->verts[j]];
	//	//	pos_dir.entry[0]=v->x-ce[0];
	//	//	pos_dir.entry[1]=v->y-ce[1];
	//	//	
	//	//	if(dot(pos_dir, norm)<0)
	//	//		continue;

	//	//	if(v->inland && ang_traceDir_boundDir<cross_angle)
	//	//	{
	//	//		cross_river=true;
	//	//		break;
	//	//	}
	//	//}
	//	
	//	if(!is_not_inland_cell_weak(trianglelist->elems[i])
	//&& ang_traceDir_boundDir<cross_angle)
	//	{
	//		cross_river=true;
	//		break;
	//	}
	//}

	Trajectory *tempTraj=new Trajectory(-1);
	ce[0]=pt[0]+1.0*min_waterwidth*quadmesh->xinterval*norm.entry[0];
	ce[1]=pt[1]+1.0*min_waterwidth*quadmesh->xinterval*norm.entry[1];
	ce_cell=get_cellID_givencoords(ce[0], ce[1]);
	if(sharedvars.AllowMajRoadCrossRiverOn &&
		is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
			fittedmap1, 512)&& ang_traceDir_boundDir<cross_angle)
	{
		/*  find a point in the other side of the water  */

		/*  we need to do a bindary search to find out the proper 
			distance to cross the river
		*/
		double smaller_dist=MinDistCrossRiver/2.;
		ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
		ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
		ce_cell=get_cellID_givencoords(ce[0], ce[1]);
		if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
				fittedmap1, 512))
		{
			smaller_dist=.75*MinDistCrossRiver;
			ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
			ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
			ce_cell=get_cellID_givencoords(ce[0], ce[1]);
			if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
					fittedmap1, 512))
			{
				ce[0]=pt[0]+MinDistCrossRiver*quadmesh->xinterval*norm.entry[0];
				ce[1]=pt[1]+MinDistCrossRiver*quadmesh->xinterval*norm.entry[1];
				ce_cell=get_cellID_givencoords(ce[0], ce[1]);
			}
		}
		else
		{
			smaller_dist=MinDistCrossRiver/4.;
			ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
			ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
			ce_cell=get_cellID_givencoords(ce[0], ce[1]);
			if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
				fittedmap1, 512))
			{
				smaller_dist=MinDistCrossRiver*.375;
				ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
				ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
				ce_cell=get_cellID_givencoords(ce[0], ce[1]);
				if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
					fittedmap1, 512))
				{
					smaller_dist=MinDistCrossRiver/2.;
					ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
					ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
					ce_cell=get_cellID_givencoords(ce[0], ce[1]);
				}
			}
		}

		get_linesegs_anytwopts(pt, cell,  ce, ce_cell, tempTraj, 0, 50);
		tenline_dir_global.entry[0]=ce[0]-pt[0];
		tenline_dir_global.entry[1]=ce[1]-pt[1];
		pt[0]=ce[0];
		pt[1]=ce[1];
		cell=ce_cell;
	}
	else
	{
		norm=-norm;
		ang_traceDir_boundDir=acos(dot(norm, traceDir));
		ce[0]=pt[0]+1.0*min_waterwidth*quadmesh->xinterval*norm.entry[0];
		ce[1]=pt[1]+1.0*min_waterwidth*quadmesh->xinterval*norm.entry[1];
		ce_cell=get_cellID_givencoords(ce[0], ce[1]);
		if(sharedvars.AllowMajRoadCrossRiverOn &&
			is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
				fittedmap1, 512)&& ang_traceDir_boundDir<cross_angle)
		{
			/*  find a point in the other side of the water  */
			/*  we need to do a bindary search to find out the proper 
				distance to cross the river
			*/
			double smaller_dist=MinDistCrossRiver/2.;
			ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
			ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
			ce_cell=get_cellID_givencoords(ce[0], ce[1]);
			if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
					fittedmap1, 512))
			{
				smaller_dist=.75*MinDistCrossRiver;
				ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
				ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
				ce_cell=get_cellID_givencoords(ce[0], ce[1]);
				if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
						fittedmap1, 512))
				{
					ce[0]=pt[0]+MinDistCrossRiver*quadmesh->xinterval*norm.entry[0];
					ce[1]=pt[1]+MinDistCrossRiver*quadmesh->xinterval*norm.entry[1];
					ce_cell=get_cellID_givencoords(ce[0], ce[1]);
				}
			}
			else
			{
				smaller_dist=MinDistCrossRiver/4.;
				ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
				ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
				ce_cell=get_cellID_givencoords(ce[0], ce[1]);
				if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
					fittedmap1, 512))
				{
					smaller_dist=MinDistCrossRiver*.375;
					ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
					ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
					ce_cell=get_cellID_givencoords(ce[0], ce[1]);
					if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
						fittedmap1, 512))
					{
						smaller_dist=MinDistCrossRiver/2.;
						ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
						ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
						ce_cell=get_cellID_givencoords(ce[0], ce[1]);
					}
				}
			}

			get_linesegs_anytwopts(pt, cell,  ce, ce_cell, tempTraj, 0, 50);
			tenline_dir_global.entry[0]=ce[0]-pt[0];
			tenline_dir_global.entry[1]=ce[1]-pt[1];
			pt[0]=ce[0];
			pt[1]=ce[1];
			cell=ce_cell;
		}

		else if(sharedvars.AllowMajRoadFollowBoundaryOn)
			/*  follow the "closest_bound" boundary for "dist_follow_bound" distance  */
			follow_boundary(closest_bound, startpos, appro_bound_dir, dist_follow_bound*quadmesh->xinterval, 
				tempTraj, pt, cell);

		else
		{
			/*  let's try to move away from the boundary a little bit */
			double cosang=dot(appro_bound_dir, traceDir);
			if(cosang>cos(M_PI/4.))
			{
				icVector2 newdir;
				double rot_ang=M_PI/6;
				//newdir.entry[0]=cos(rot_ang)*appro_bound_dir.entry[0]-sin(rot_ang)*appro_bound_dir.entry[1];
				//newdir.entry[1]=sin(rot_ang)*appro_bound_dir.entry[0]+cos(rot_ang)*appro_bound_dir.entry[1];
				newdir.entry[0]=cos(rot_ang)*traceDir.entry[0]-sin(rot_ang)*traceDir.entry[1];
				newdir.entry[1]=sin(rot_ang)*traceDir.entry[0]+cos(rot_ang)*traceDir.entry[1];
				ce[0]=pt[0]+/*smaller_dist**/quadmesh->xinterval*newdir.entry[0]*1.5;
				ce[1]=pt[1]+/*smaller_dist**/quadmesh->xinterval*newdir.entry[1]*1.5;
				ce_cell=get_cellID_givencoords(ce[0], ce[1]);
				//if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
				//	fittedmap1, 512))
				if(is_not_inland_cell_weak(ce_cell))
				{
					//newdir.entry[0]=cos(rot_ang)*appro_bound_dir.entry[0]+sin(rot_ang)*appro_bound_dir.entry[1];
					//newdir.entry[1]=-sin(rot_ang)*appro_bound_dir.entry[0]+cos(rot_ang)*appro_bound_dir.entry[1];
					newdir.entry[0]=cos(rot_ang)*traceDir.entry[0]+sin(rot_ang)*traceDir.entry[1];
					newdir.entry[1]=-sin(rot_ang)*traceDir.entry[0]+cos(rot_ang)*traceDir.entry[1];
					ce[0]=pt[0]+/*smaller_dist**/quadmesh->xinterval*newdir.entry[0]*1.5;
					ce[1]=pt[1]+/*smaller_dist**/quadmesh->xinterval*newdir.entry[1]*1.5;
					ce_cell=get_cellID_givencoords(ce[0], ce[1]);
					//if(is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
					//	fittedmap1, 512))
					if(!is_not_inland_cell_weak(ce_cell))
					{
						get_linesegs_anytwopts(pt, cell,  ce, ce_cell, tempTraj, 0, 50);
						tenline_dir_global.entry[0]=ce[0]-pt[0];
						tenline_dir_global.entry[1]=ce[1]-pt[1];
						pt[0]=ce[0];
						pt[1]=ce[1];
						cell=ce_cell;
					}
				}
				else
				{
					get_linesegs_anytwopts(pt, cell,  ce, ce_cell, tempTraj, 0, 50);
					tenline_dir_global.entry[0]=ce[0]-pt[0];
					tenline_dir_global.entry[1]=ce[1]-pt[1];
					pt[0]=ce[0];
					pt[1]=ce[1];
					cell=ce_cell;
			}
			}
		}
	}

	/*  copy the line segment in tempTraj to current trajectory data structure  */
	evenstreamlines->trajs[evenstreamlines->ntrajs]->add_last_nlines(
		tempTraj->linesegs, tempTraj->nlinesegs);

	delete tempTraj;
}


void EvenStreamlinePlace::comp_smallest_dist_to_one_mapBound(double pt[2], int cell, int boundID, 
															 double &dist, double intersect[2],
															 int &pos, icVector2 &appro_bound_dir)
{
	int i;
	MapBoundary *thebound=&mapboundarylist->mapboundarylist[boundID];
	icVector2 dist_vect;
	double temp_dist, Distline;
	pos=-1;
	for(i=0;i<thebound->nelems-1;i++)
	{
		if(thebound->pts[i]==NULL)
			continue;

		if(get_cellID_givencoords(thebound->pts[i]->x, thebound->pts[i]->y)!=cell)
			continue;

		DistanceFromLine(pt[0],pt[1],  thebound->pts[i]->x,thebound->pts[i]->y,
			thebound->pts[i+1]->x,thebound->pts[i+1]->y,   temp_dist, Distline);
		if(temp_dist<dist)
		{
			intersect[0]=(thebound->pts[i]->x+thebound->pts[i+1]->x)/2.;
			intersect[1]=(thebound->pts[i]->y+thebound->pts[i+1]->y)/2.;
			dist=temp_dist;
			if(i==0)
			{
				appro_bound_dir.entry[0]=thebound->pts[i+1]->x-thebound->pts[i]->x;
				appro_bound_dir.entry[1]=thebound->pts[i+1]->y-thebound->pts[i]->y;
			}
			else
			{
				appro_bound_dir.entry[0]=thebound->pts[i+1]->x-thebound->pts[i-1]->x;
				appro_bound_dir.entry[1]=thebound->pts[i+1]->y-thebound->pts[i-1]->y;
			}
			//appro_bound_dir.entry[0]=thebound->pts[i+1]->x-thebound->pts[i]->x;
			//appro_bound_dir.entry[1]=thebound->pts[i+1]->y-thebound->pts[i]->y;
			pos=i;
		}

		if(i==thebound->nelems-1)
		{
			/*  compute the Euclidean distance  */
			dist_vect.entry[0]=pt[0]-thebound->pts[i]->x;
			dist_vect.entry[1]=pt[1]-thebound->pts[i]->y;
			temp_dist=length(dist_vect);

			if(temp_dist<dist)
			{
				intersect[0]=thebound->pts[i]->x;
				intersect[1]=thebound->pts[i]->y;
				dist=temp_dist;
				appro_bound_dir.entry[0]=thebound->pts[i]->x-thebound->pts[i-1]->x;
				appro_bound_dir.entry[1]=thebound->pts[i]->y-thebound->pts[i-1]->y;
				pos=i;
			}
		}
	}
}

bool EvenStreamlinePlace::find_closest_bound(double pt[2], int cell, int &boundID, 
											 double intersect[2], int &pos, icVector2 &appro_bound_dir)
{
	int i;
	QuadCell *face=quadmesh->quadcells[cell];
	if(face->mapbounds==NULL || face->nmapbounds==0) return false;

	double smallest_dist=1.e50;
	double t_intersect[2];
	boundID=0;
	bool found=false;
	icVector2 t_appro_bound_dir;
	int t_pos=-1;
	for(i=0;i<face->nmapbounds;i++)
	{
		double temp_dist=1.e50;
		comp_smallest_dist_to_one_mapBound(pt, cell, face->mapbounds[i], temp_dist, 
			t_intersect, t_pos, t_appro_bound_dir);

		if(temp_dist<smallest_dist)
		{
			smallest_dist=temp_dist;
			intersect[0]=t_intersect[0];
			intersect[1]=t_intersect[1];
			boundID=face->mapbounds[i];
			appro_bound_dir=t_appro_bound_dir;
			pos=t_pos;
			found=true;
		}
	}
	return found;
}


/*
   Adaptively adjust the step size to cross the river 
*/
bool adp_crossRiver_test(double pt[2], int in_cell, double ce[2], int &ce_cell, icVector2 norm, 
						 double max_crossRiverDist)
{
	ce[0]=pt[0]+1.0*max_crossRiverDist*quadmesh->xinterval*norm.entry[0];
	ce[1]=pt[1]+1.0*max_crossRiverDist*quadmesh->xinterval*norm.entry[1];
	ce_cell=get_cellID_givencoords(ce[0], ce[1]);

	if(ce_cell<0 || ce_cell>=quadmesh->nfaces)
		return false;

	if(sharedvars.AllowMajRoadCrossRiverOn &&
		//is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
		//	fittedmap1, 512)&& ang_traceDir_boundDir<cross_angle)
		!is_not_inland_cell_weak(ce_cell))
	{
		/*  find a point in the other side of the water  */

		/*  we need to do a bindary search to find out the proper 
			distance to cross the river
		*/
		double smaller_dist=max_crossRiverDist/2.;
		ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
		ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
		ce_cell=get_cellID_givencoords(ce[0], ce[1]);
		//if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
		//		fittedmap1, 512))
		if(is_not_inland_cell_weak(ce_cell))
		{
			smaller_dist=.75*max_crossRiverDist;
			ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
			ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
			ce_cell=get_cellID_givencoords(ce[0], ce[1]);
			//if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
			//		fittedmap1, 512))
			if(is_not_inland_cell_weak(ce_cell))
			{
				ce[0]=pt[0]+max_crossRiverDist*quadmesh->xinterval*norm.entry[0];
				ce[1]=pt[1]+max_crossRiverDist*quadmesh->xinterval*norm.entry[1];
				ce_cell=get_cellID_givencoords(ce[0], ce[1]);
			}
		}
		else
		{
			smaller_dist=max_crossRiverDist/4.;
			ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
			ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
			ce_cell=get_cellID_givencoords(ce[0], ce[1]);
			//if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
			//	fittedmap1, 512))
			if(is_not_inland_cell_weak(ce_cell));
			{
				smaller_dist=max_crossRiverDist*.375;
				ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
				ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
				ce_cell=get_cellID_givencoords(ce[0], ce[1]);
				//if(!is_inland_pixel(ce[0],ce[1], xstart, xrang, ystart, yrang, dx, dy,
				//	fittedmap1, 512))
				if(is_not_inland_cell_weak(ce_cell));
				{
					smaller_dist=max_crossRiverDist/2.;
					ce[0]=pt[0]+smaller_dist*quadmesh->xinterval*norm.entry[0];
					ce[1]=pt[1]+smaller_dist*quadmesh->xinterval*norm.entry[1];
					ce_cell=get_cellID_givencoords(ce[0], ce[1]);
				}
			}
		}

	}
	else
		return false;

	return true;
}


/*
    compute a major road
*/
bool EvenStreamlinePlace::grow_a_majRoad(double seed_p[2], int triangle, double dtest, 
											double discsize, double Sample_interval, 
											double loopdsep, double dist2sing, 
											double streamlinelength, 
											int type, bool brushon)
{
	int i;
	int flag = -1;

	int pre_face, cur_face;
	double globalp[3] = {0.};
	int cur_line = 0;
	double cur_length = 0;
	int movetonext = 0;
		
	/*  variables for pixel level judgement of "inland" or not  */
	double xstart=quadmesh->xstart;
	double ystart=quadmesh->ystart;
	double xrang=quadmesh->xend-quadmesh->xstart;
	double yrang=quadmesh->yend-quadmesh->ystart;
	double dx=xrang/511;
	double dy=yrang/511;
	double dist_did;


	pre_face = cur_face = triangle;
	globalp[0] = seed_p[0];
	globalp[1] = seed_p[1];
	
	FILE *fp;
	
	icMatrix2x2 ten;
	//double t[4]={0.};

	if(!is_in_cell(triangle, seed_p[0], seed_p[1]))
	{
		triangle = get_cellID_givencoords(seed_p[0], seed_p[1]);
	}

	if(triangle < 0 || triangle >= quadmesh->nfaces)
		return false;
	
	compute_tensor_at_quad(triangle, seed_p[0], seed_p[1], ten);

	double evalues[2] = {0.};
	icVector2 ev[2], startdir;
	cal_eigen_vector_sym(ten, ev);

	if(type == 1) /*use minor eigen vector field*/
		ev[0] = ev[1];

	tenline_dir_global = startdir = -ev[0];  /*obtain the major eigen vector*/


	evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs = 0;
	evenstreamlines->trajs[evenstreamlines->ntrajs]->closed=false;


	samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[0] = seed_p[0];
	samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[1] = seed_p[1];
	samplepts[evenstreamlines->ntrajs]->samples[0]->triangle = triangle;
	samplepts[evenstreamlines->ntrajs]->nsamples=1;
	

	hstep = quadmesh->xinterval/2.;
	predict_stepsize = quadmesh->xinterval/2.;
	euler_stepsize = quadmesh->xinterval/10.;


	cur_line = 0;
	cur_length = 0;

		/*record the postion of the original seed point in the tensor line 11/18/2007*/
	seedposition_ineachtensorline[evenstreamlines->ntrajs]=0;

	//min_waterwidth=quadmesh->xinterval*2;

	//////////////////////////////////////////////////////////////////////////

	////Backward tracing
	int NUMTRACETRIS = (int)sqrt((double)quadmesh->nfaces);

	for(i = 0; i < 3*NUMTRACETRIS; i++)
	{

		////The cell does not exist. Something is wrong!
		if(cur_face < 0 || cur_face >= quadmesh->nfaces)
			break;


		/*it reaches any boundary, we should stop as well*/
		if(globalp[0]>=quadmesh->xend-1.e-8||globalp[0]<=quadmesh->xstart+1.e-8
			||globalp[1]>=quadmesh->yend-1.e-8||globalp[1]<=quadmesh->ystart+1.e-8)
		{
			if(globalp[0]>=quadmesh->xend-1.e-8) globalp[0]=quadmesh->xend;
			else if(globalp[0]<=quadmesh->xstart+1.e-8) globalp[0]=quadmesh->xstart;
			if(globalp[1]>=quadmesh->yend-1.e-8) globalp[1]=quadmesh->yend;
			else if(globalp[1]<=quadmesh->ystart+1.e-8) globalp[1]=quadmesh->ystart;
			break;
		}

		/*for brush interface 10/10/2007*/
		if(brushon && !quadmesh->quadcells[cur_face]->OnBoundary)
			break;

		if(is_inveg_cell_weak(cur_face))
			break;

		/*for loading the map 10/24/2007*/
		//if(is_not_inland(cur_face))
		if(is_not_inland_cell_weak(cur_face))
		//if(!is_inland_pixel(globalp[0],globalp[1], xstart, xrang, ystart, yrang, dx, dy,
		//		fittedmap1, 512))
		{
			if(!sharedvars.AllowMajRoadCrossRiverOn && !sharedvars.AllowMajRoadFollowBoundaryOn)
				break;

			/*   deal with water region   */
			//trace_water_region(globalp, cur_face, Sample_interval, loopdsep, min_waterwidth,
			//	M_PI/6, dist_follow_bound, type);
			trace_water_region(globalp, cur_face, Sample_interval, loopdsep, MinDistCrossRiver,
				MinAngCrossRiver, MaxDistFollowBoundary, type);

			if(cur_face<0||cur_face>=quadmesh->nfaces)
				break;

			//if(is_not_inland(cur_face))
			if(is_not_inland_cell_weak(cur_face))
			//if(!is_inland_pixel(globalp[0],globalp[1], xstart, xrang, ystart, yrang, dx, dy,
			//		fittedmap1, 512))
			{
				/*  How about we do one more test to make sure it will not cross river  */
				if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs==0)
					break;

				icVector2 trajDir;
				int curLineSegs=evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs;
				Trajectory *curTraj=evenstreamlines->trajs[evenstreamlines->ntrajs];
				trajDir.entry[0]=curTraj->linesegs[curLineSegs-1].gend[0]-
					curTraj->linesegs[curLineSegs-1].gstart[0];
				trajDir.entry[1]=curTraj->linesegs[curLineSegs-1].gend[1]-
					curTraj->linesegs[curLineSegs-1].gstart[1];

				if(length(trajDir)<1.e-8)
					break;

				normalize(trajDir);
				//double ce_test[2];
				//ce_test[0]=globalp[0]+MinDistCrossRiver*quadmesh->xinterval*trajDir.entry[0];
				//ce_test[1]=globalp[1]+MinDistCrossRiver*quadmesh->xinterval*trajDir.entry[1];
				//int cell_test=get_cellID_givencoords(ce_test[0], ce_test[1]);

				//if(is_not_inland_cell_weak(cell_test))
				//	break;
				//if(cell_test < 0 || cell_test >= quadmesh->nfaces)
				//	break;

				/*  use adpative test here */
				double ce_test[2];
				int cell_test;
				
				if(!adp_crossRiver_test(globalp, cur_face, ce_test, cell_test, 
					trajDir, MinDistCrossRiver))
					break;

				Trajectory *tempTraj=new Trajectory(-1);
				get_linesegs_anytwopts(globalp, cur_face,  ce_test, cell_test, tempTraj, 0, 50);
				curTraj->add_last_nlines(tempTraj->linesegs, tempTraj->nlinesegs);
				globalp[0]=ce_test[0];
				globalp[1]=ce_test[1];
				cur_face=cell_test;
				delete tempTraj;
			}
		}
		
		pre_face = cur_face;
        cur_face = trace_majRoad_in_quad(cur_face, globalp, type, dtest, loopdsep, dist2sing,
			Sample_interval, discsize, flag);

		//if this is the first line segment, there may be a case that the start point is on the edge !!!
		if(pre_face == triangle && evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
		{
			cur_length = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].length;
		}

		
		////We need to select the sampling points from current trajectory  4/22/06
		cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
			samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
	

		if(flag == 1 || flag == 4 || flag == 2  || (flag != 3 && pre_face == cur_face)/**/ ) 
		{
			/*   we still need to connect to the other kind of major roads   */
			Trajectory *temp_traj=new Trajectory(-1);
			Trajectory *curtraj=evenstreamlines->trajs[evenstreamlines->ntrajs];

			//icVector2 trajDir;
			//trajDir.entry[0]=curtraj->linesegs[curtraj->nlinesegs-1].gend[0]-
			//	curtraj->linesegs[curtraj->nlinesegs-1].gstart[0];
			//trajDir.entry[1]=curtraj->linesegs[curtraj->nlinesegs-1].gend[1]-
			//	curtraj->linesegs[curtraj->nlinesegs-1].gstart[1];

			//temp_traj->nlinesegs=0;
			//if(connect_one_majRoad_to_exist(evenstreamlines->ntrajs, majororminor,
			//	trajDir, globalp, cur_face, temp_traj, dsep, M_PI/2.))
			//{
			//	curtraj->add_last_nlines(temp_traj->linesegs, temp_traj->nlinesegs);
			//	cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
			//		samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
			//}

			if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
			{
				globalp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs
					[evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs-1].gend[0];
				globalp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs
					[evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs-1].gend[1];
				cur_face=get_cellID_givencoords(globalp[0], globalp[1]);
			}


			trace_to_comp_nearby_traj(globalp, cur_face, majororminor, dsep, 
				100, MinDistCrossRiver, 1.5*quadmesh->xinterval, dist_did, temp_traj);
			curtraj->add_last_nlines(temp_traj->linesegs, temp_traj->nlinesegs);
			cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
				samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
			
			delete temp_traj;

			break;
		}

		if(flag == 3 && is_valid_loop(evenstreamlines->ntrajs, loopdsep) 
			&& sharedvars.CloseLoopOn) /*the tensor line forms a loop here*/
		{
			/*connect the whole tensor line*/
			int cur_line_index=evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs;
			if(cur_line_index<=0) break;

			double startp[2], endp[2];
			Trajectory *temp=new Trajectory(-1);
			startp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[0];
			startp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[1];
			int cell1=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].Triangle_ID;
			endp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
			endp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
			int cell2=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;

			temp->nlinesegs=0;

			get_linesegs_anytwopts(startp,cell1,  endp,cell2, temp, 0, 10);
			
			evenstreamlines->trajs[evenstreamlines->ntrajs]->add_last_nlines(
				temp->linesegs, temp->nlinesegs);
			delete temp;

			evenstreamlines->trajs[evenstreamlines->ntrajs]->closed=true;
		   
			cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
				samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);

			goto L20;
		}

	}

 	////Reverse the order of the obtained line segments
	if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
	{
		reverse_streamline(evenstreamlines->ntrajs);

		//////Resample the streamline after the reversion
		cur_line = 0;
		cur_length = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].length;
		samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[0] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
		samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[1] = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
		samplepts[evenstreamlines->ntrajs]->samples[0]->triangle = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;
		samplepts[evenstreamlines->ntrajs]->nsamples = 1;

		movetonext = 0;
	    cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
			samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);

		/*record the postion of the original seed point in the tensor line 11/18/2007*/
		seedposition_ineachtensorline[evenstreamlines->ntrajs]=evenstreamlines->trajs[
			evenstreamlines->ntrajs]->nlinesegs-1;
	}
	

	//////////Forward tracing
	pre_face = cur_face = triangle;
	globalp[0] = seed_p[0];
	globalp[1] = seed_p[1];
	flag = -1;
	tenline_dir_global = -startdir;

	for(i = 0; i < 3*NUMTRACETRIS; i++)
	{

		if(cur_face < 0 || cur_face >= quadmesh->nfaces)
			break;

		if(globalp[0]>=quadmesh->xend-1.e-8||globalp[0]<=quadmesh->xstart+1.e-8
			||globalp[1]>=quadmesh->yend-1.e-8||globalp[1]<=quadmesh->ystart+1.e-8)
		{
			if(globalp[0]>=quadmesh->xend-1.e-8) globalp[0]=quadmesh->xend;
			else if(globalp[0]<=quadmesh->xstart+1.e-8) globalp[0]=quadmesh->xstart;
			if(globalp[1]>=quadmesh->yend-1.e-8) globalp[1]=quadmesh->yend;
			else if(globalp[1]<=quadmesh->ystart+1.e-8) globalp[1]=quadmesh->ystart;
			break;
		}

		
		if(is_inveg_cell_weak(cur_face))
			break;
		
		/*for brush interface 10/10/2007*/
		if(brushon && !quadmesh->quadcells[cur_face]->OnBoundary)
			break;
		
		/*for load the map 10/24/2007*/
		//if(is_not_inland(cur_face))
		if(is_not_inland_cell_weak(cur_face))
		//if(!is_inland_pixel(globalp[0],globalp[1], xstart, xrang, ystart, yrang, dx, dy,
		//		fittedmap1, 512))
		{
			if(!sharedvars.AllowMajRoadCrossRiverOn && !sharedvars.AllowMajRoadFollowBoundaryOn)
				break;

			/*   deal with water region   */
			//trace_water_region(globalp, cur_face, Sample_interval, loopdsep, min_waterwidth,
			//	M_PI/6, dist_follow_bound, type);
			trace_water_region(globalp, cur_face, Sample_interval, loopdsep, MinDistCrossRiver,
				MinAngCrossRiver, MaxDistFollowBoundary, type);

			if(cur_face<0||cur_face>=quadmesh->nfaces)
				break;

			//if(is_not_inland(cur_face))
			//if(is_not_inland_cell_weak(cur_face))
			if(!is_inland_pixel(globalp[0],globalp[1], xstart, xrang, ystart, yrang, dx, dy,
					fittedmap1, 512))
			{
				/*  How about we do one more test to make sure it will not cross river  */
				if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs==0)
					break;

				icVector2 trajDir;
				int curLineSegs=evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs;
				Trajectory *curTraj=evenstreamlines->trajs[evenstreamlines->ntrajs];
				trajDir.entry[0]=curTraj->linesegs[curLineSegs-1].gend[0]-
					curTraj->linesegs[curLineSegs-1].gstart[0];
				trajDir.entry[1]=curTraj->linesegs[curLineSegs-1].gend[1]-
					curTraj->linesegs[curLineSegs-1].gstart[1];

				if(length(trajDir)<1.e-8)
					break;

				normalize(trajDir);
				//double ce_test[2];
				//ce_test[0]=globalp[0]+MinDistCrossRiver*quadmesh->xinterval*trajDir.entry[0];
				//ce_test[1]=globalp[1]+MinDistCrossRiver*quadmesh->xinterval*trajDir.entry[1];
				//int cell_test=get_cellID_givencoords(ce_test[0], ce_test[1]);

				//if(is_not_inland_cell_weak(cell_test))
				//	break;

				//if(cell_test < 0 || cell_test >= quadmesh->nfaces)
				//	break;

				/*  use adpative test here */
				double ce_test[2];
				int cell_test;
				
				if(!adp_crossRiver_test(globalp, cur_face, ce_test, cell_test, 
					trajDir, MinDistCrossRiver))
					break;

				Trajectory *tempTraj=new Trajectory(-1);
				get_linesegs_anytwopts(globalp, cell_test,  ce_test, cell_test, tempTraj, 0, 50);
				curTraj->add_last_nlines(tempTraj->linesegs, tempTraj->nlinesegs);
				globalp[0]=ce_test[0];
				globalp[1]=ce_test[1];
				cur_face=cell_test;
				delete tempTraj;
			}
		}
		
		///*it reaches any boundary, we should stop as well*/
		//if(globalp[0]>quadmesh->xend-1.e-8||globalp[0]<quadmesh->xstart+1.e-8
		//	||globalp[1]>quadmesh->yend-1.e-8||globalp[1]<quadmesh->ystart+1.e-8)
		//	break;
		

		pre_face = cur_face;
        //cur_face = trace_in_quad_even(cur_face, globalp, type, dtest, loopdsep, dist2sing,
        cur_face = trace_majRoad_in_quad(cur_face, globalp, type, dtest, loopdsep, dist2sing,
			Sample_interval, discsize, flag);
		

		////We need to select the sampling points from current trajectory  3/9/06
	    cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
			samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
		
		
	    if(flag == 1 || flag == 4 || flag == 2 || (flag != 3&& pre_face == cur_face)/**/ ) 
		{
			/*   we still need to connect to the other kind of major roads   */
			Trajectory *temp_traj=new Trajectory(-1);
			Trajectory *curtraj=evenstreamlines->trajs[evenstreamlines->ntrajs];

			//icVector2 trajDir;
			//trajDir.entry[0]=curtraj->linesegs[curtraj->nlinesegs-1].gend[0]-
			//	curtraj->linesegs[curtraj->nlinesegs-1].gstart[0];
			//trajDir.entry[1]=curtraj->linesegs[curtraj->nlinesegs-1].gend[1]-
			//	curtraj->linesegs[curtraj->nlinesegs-1].gstart[1];

			//temp_traj->nlinesegs=0;
			//if(connect_one_majRoad_to_exist(evenstreamlines->ntrajs, majororminor,
			//	trajDir, globalp, cur_face, temp_traj, dsep, M_PI/2.))
			//{
			//	curtraj->add_last_nlines(temp_traj->linesegs, temp_traj->nlinesegs);
			//	cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
			//		samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
			//}

			if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
			{
				globalp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs
					[evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs-1].gend[0];
				globalp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs
					[evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs-1].gend[1];
				cur_face=get_cellID_givencoords(globalp[0], globalp[1]);
			}

			trace_to_comp_nearby_traj(globalp, cur_face, majororminor, dsep, 
				100, MinDistCrossRiver, 1.5*quadmesh->xinterval, dist_did, temp_traj);
			curtraj->add_last_nlines(temp_traj->linesegs, temp_traj->nlinesegs);
			cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
				samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
			
			delete temp_traj;

			break;
		}

		if(flag == 3 && is_valid_loop(evenstreamlines->ntrajs, loopdsep)
			&& sharedvars.CloseLoopOn) /*the tensor line forms a loop here*/
		{
			/*connect the whole tensor line*/
			/*we need to extend if not enough memory*/
			int cur_line_index=evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs;
			if(cur_line_index<=0) break;
			
			double startp[2], endp[2];
			Trajectory *temp=new Trajectory(-1);
			startp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[0];
			startp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[1];
			int cell1=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].Triangle_ID;
			endp[0]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
			endp[1]=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
			int cell2=evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;

			temp->nlinesegs=0;

			get_linesegs_anytwopts(startp,cell1,  endp,cell2, temp, 0, 10);
			
			evenstreamlines->trajs[evenstreamlines->ntrajs]->add_last_nlines(
				temp->linesegs, temp->nlinesegs);
			delete temp;
		   
			cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
				samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
			
			evenstreamlines->trajs[evenstreamlines->ntrajs]->closed=true;

			break;
		}
	}
	
    
L20:	if((evenstreamlines->ntrajs == 0 ||
		evenstreamlines->trajs[evenstreamlines->ntrajs]->get_length() > streamlinelength)
		&& evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
	{
		update_samples_in_cell(evenstreamlines->ntrajs, samplepts[evenstreamlines->ntrajs]->samples,
			samplepts[evenstreamlines->ntrajs]->nsamples);

		Trajectory *cur_traj = evenstreamlines->trajs[evenstreamlines->ntrajs];
		int pre_cell = cur_traj->linesegs[0].Triangle_ID;
		int start, end;
		start = end = 0;
		for(i=1; i<cur_traj->nlinesegs; i++)
		{
			/*  this contains a bug  */
			if(cur_traj->linesegs[i].Triangle_ID <0 
				||cur_traj->linesegs[i].Triangle_ID >=quadmesh->nfaces)
				continue;

			if(cur_traj->linesegs[i].Triangle_ID == pre_cell)
			{
				end = i;
				continue;
			}
			else{
				LinesInOneCell *newline = (LinesInOneCell*)malloc(sizeof(LinesInOneCell));
				newline->whichtraj = evenstreamlines->ntrajs;
				newline->start = start;
				newline->end = end;
				QuadCell *curc = quadmesh->quadcells[pre_cell];

				//if it is major field tracing, add it to the major line info
				if(!majororminor)
				{
					if(curc->majorlines == NULL)
						curc->majorlines = new LineInfo(1);
					curc->majorlines->addNew(newline);
					curc->hasmajor=true;
				}

				//else, add to the minor line info
				else
				{
					if(curc->minorlines == NULL)
						curc->minorlines = new LineInfo(1);
					curc->minorlines->addNew(newline);
					curc->hasminor=true;
				}
				start = end = i;
				pre_cell = cur_traj->linesegs[i].Triangle_ID;
			}
		}

		/*we need to handle the last line segments*/
		start = end = cur_traj->nlinesegs-1;
		pre_cell = cur_traj->linesegs[start].Triangle_ID;
		if(pre_cell>=0&&pre_cell<quadmesh->nfaces)
		{
			for(i=cur_traj->nlinesegs-2; i>=0; i--)
			{
				/*  this contains a bug  */
				if(cur_traj->linesegs[i].Triangle_ID <0 
					||cur_traj->linesegs[i].Triangle_ID >=quadmesh->nfaces)
					continue;

				if(cur_traj->linesegs[i].Triangle_ID == pre_cell)
				{
					start = i;
					continue;
				}
				else{
					LinesInOneCell *newline = (LinesInOneCell*)malloc(sizeof(LinesInOneCell));
					newline->whichtraj = evenstreamlines->ntrajs;
					newline->start = start;
					newline->end = end;
					QuadCell *curc = quadmesh->quadcells[pre_cell];

					//if it is major field tracing, add it to the major line info
					if(!majororminor)
					{
						if(curc->majorlines == NULL)
							curc->majorlines = new LineInfo(1);
						curc->majorlines->addNew(newline);
						curc->hasmajor=true;
					}

					//else, add to the minor line info
					else
					{
						if(curc->minorlines == NULL)
							curc->minorlines = new LineInfo(1);
						curc->minorlines->addNew(newline);
						curc->hasminor=true;
					}
					start = end = i;
					pre_cell = cur_traj->linesegs[i].Triangle_ID;
					break;
				}
			}
		}
		else
		{
			/*  we need to handle this!  */
			int test=0;
			cur_traj->nlinesegs--;
			while(cur_traj->linesegs[cur_traj->nlinesegs-1].Triangle_ID<0
				||cur_traj->linesegs[cur_traj->nlinesegs-1].Triangle_ID>=quadmesh->nfaces)
				cur_traj->nlinesegs--;
		}

		/*set the corresponding road type as default one 10/03/2007*/
		if(majororminor)
			evenstreamlines->trajs[evenstreamlines->ntrajs]->roadtype = MINOR;
		else
			evenstreamlines->trajs[evenstreamlines->ntrajs]->roadtype = MINOR;

		/*reorder the line segments of the streamline 10/13/2007*/

		evenstreamlines->trajs[evenstreamlines->ntrajs]->traj_len=
			evenstreamlines->trajs[evenstreamlines->ntrajs]->get_length();

		evenstreamlines->ntrajs ++;

		return true;
	}

	return false;

}

bool EvenStreamlinePlace::grow_a_tensorline_withoutinversing
		(double seed_p[2], int triangle, double dtest, double discsize, double Sample_interval, 
		double loopdsep, double dist2sing, double streamlinelength, int type, bool brushon)
{
	int i;
	int flag = -1;

	int pre_face, cur_face;
	double globalp[3] = {0.};
	int cur_line = 0;
	double cur_length = 0;
	int movetonext = 0;
		
	if(triangle < 0 || triangle > quadmesh->nfaces)
		return false;

	pre_face = cur_face = triangle;
	globalp[0] = seed_p[0];
	globalp[1] = seed_p[1];
	
	FILE *fp;
	
	icMatrix2x2 ten;
	double t[4]={0.};

	if(!is_in_reg_cell(triangle, seed_p[0], seed_p[1]))
	{
		triangle = get_cellID_givencoords(seed_p[0], seed_p[1]);
	}

	compute_tensor_at_quad(triangle, seed_p[0], seed_p[1], ten);

	double evalues[2] = {0.};
	icVector2 ev[2], startdir;
	cal_eigen_vector_sym(ten, ev);

	if(type == 1) /*use minor eigen vector field*/
		ev[0] = ev[1];

	tenline_dir_global = startdir = -ev[0];  /*obtain the major eigen vector*/


	evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs = 0;

	samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[0] = seed_p[0];
	samplepts[evenstreamlines->ntrajs]->samples[0]->gpt[1] = seed_p[1];
	samplepts[evenstreamlines->ntrajs]->samples[0]->triangle = triangle;
	samplepts[evenstreamlines->ntrajs]->nsamples=1;
	
extern double hstep;

	hstep = quadmesh->xinterval/2.;
	predict_stepsize = quadmesh->xinterval/2.;
	euler_stepsize = quadmesh->xinterval/4.;


	cur_line = 0;
	cur_length = 0;

		/*record the postion of the original seed point in the tensor line 11/18/2007*/
	seedposition_ineachtensorline[evenstreamlines->ntrajs]=0;

	//////////////////////////////////////////////////////////////////////////


	////Backward tracing
	int NUMTRACETRIS = (int)sqrt((double)quadmesh->nfaces);

	for(i = 0; i < 3*NUMTRACETRIS; i++)
	{
		////The cell does not exist. Something is wrong!
		if(cur_face < 0 || cur_face >= quadmesh->nfaces)
		{
			break;
		}

		/*for brush interface 10/10/2007*/
		if(brushon && !quadmesh->quadcells[cur_face]->OnBoundary)
			break;

		/*for load the map 10/24/2007*/
		if(is_not_inland(cur_face))
			break;
		
		/*it reaches any boundary, we should stop as well*/
		if(globalp[0]>quadmesh->xend-1.e-8||globalp[0]<quadmesh->xstart+1.e-8
			||globalp[1]>quadmesh->yend-1.e-8||globalp[1]<quadmesh->ystart+1.e-8)
			break;

	
		pre_face = cur_face;
        cur_face = trace_in_quad_even(cur_face, globalp, type, dtest, loopdsep, dist2sing,
			Sample_interval, discsize, flag);

		//if this is the first line segment, there may be a case that the start point is on the edge !!!
		if(pre_face == triangle && evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>0)
		{
			cur_length = evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].length;
		}

		
		////We need to select the sampling points from current trajectory  4/22/06
		   cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
			samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
		

		if(flag == 1 || flag == 4 || flag == 2 /*|| flag == 3 || pre_face == cur_face*/ ) 
		{
			break;
		}

		if(flag == 3) /*the tensor line forms a loop here*/
		{
			/*connect the whole tensor line*/
			/*we need to extend if not enough memory*/
			if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->curMaxNumLinesegs)
				evenstreamlines->trajs[evenstreamlines->ntrajs]->extend_line_segments(1);

			int cur_line_index=evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs;
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gstart[0]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[0];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gstart[1]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[1];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].start[0]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].end[0];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].start[1]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].end[1];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].Triangle_ID=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].Triangle_ID;
			
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gend[0]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gend[1]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].end[0]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].start[0];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].end[1]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].start[1];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].Triangle_ID=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;

			evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs++;
		   
			cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
				samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);

			goto L10;
		}

	}

	//////////Forward tracing
	pre_face = cur_face = triangle;
	globalp[0] = seed_p[0];
	globalp[1] = seed_p[1];
	flag = -1;
	tenline_dir_global = -startdir;

	for(i = 0; i < 3*NUMTRACETRIS; i++)
	{

		if(cur_face < 0 || cur_face >= quadmesh->nfaces)
		{
			break;
		}
		
		/*for brush interface 10/10/2007*/
		if(brushon && !quadmesh->quadcells[cur_face]->OnBoundary)
			break;
		
		/*for load the map 10/24/2007*/
		if(is_not_inland(cur_face))
			break;
		
		/*it reaches any boundary, we should stop as well*/
		if(globalp[0]>quadmesh->xend-1.e-8||globalp[0]<quadmesh->xstart+1.e-8
			||globalp[1]>quadmesh->yend-1.e-8||globalp[1]<quadmesh->ystart+1.e-8)
			break;
		

		pre_face = cur_face;
        cur_face = trace_in_quad_even(cur_face, globalp, type, dtest, loopdsep, dist2sing,
			Sample_interval, discsize, flag);
		

		////We need to select the sampling points from current trajectory  3/9/06
	    cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
			samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);
		
		
	    if(flag == 1 || flag == 4 || flag == 2 /*|| flag == 3|| pre_face == cur_face*/ ) 
		{
			break;
		}

		if(flag == 3) /*the tensor line forms a loop here*/
		{
			/*connect the whole tensor line*/
			/*we need to extend if not enough memory*/
			if(evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs>=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->curMaxNumLinesegs)
				evenstreamlines->trajs[evenstreamlines->ntrajs]->extend_line_segments(1);

			int cur_line_index=evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs;
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gstart[0]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[0];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gstart[1]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].gend[1];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].start[0]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].end[0];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].start[1]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].end[1];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].Triangle_ID=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index-1].Triangle_ID;
			
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gend[0]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[0];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].gend[1]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].gstart[1];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].end[0]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].start[0];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].end[1]=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].start[1];
			evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[cur_line_index].Triangle_ID=
				evenstreamlines->trajs[evenstreamlines->ntrajs]->linesegs[0].Triangle_ID;
			
			evenstreamlines->trajs[evenstreamlines->ntrajs]->nlinesegs++;
		   
			cal_samplepts_when_tracing(evenstreamlines->ntrajs, Sample_interval, cur_line, movetonext, cur_length,
				samplepts[evenstreamlines->ntrajs]->samples, samplepts[evenstreamlines->ntrajs]->nsamples);

			break;
		}
	}
	

    
L10:	if(evenstreamlines->ntrajs == 0 ||
		evenstreamlines->trajs[evenstreamlines->ntrajs]->get_length() > streamlinelength)
	{
		update_samples_in_cell(evenstreamlines->ntrajs, samplepts[evenstreamlines->ntrajs]->samples,
			samplepts[evenstreamlines->ntrajs]->nsamples);

		//ndisplay_trajs = evenstreamlines->ntrajs-1;

		/*we need to update the line information for each cell of the mesh 10/02/2007*/
		Trajectory *cur_traj = evenstreamlines->trajs[evenstreamlines->ntrajs];
		int pre_cell = cur_traj->linesegs[0].Triangle_ID;
		int start, end;
		start = end = 0;
		for(i=1; i<cur_traj->nlinesegs; i++)
		{
			if(cur_traj->linesegs[i].Triangle_ID == pre_cell)
			{
				end = i;
				continue;
			}
			else{
				LinesInOneCell *newline = (LinesInOneCell*)malloc(sizeof(LinesInOneCell));
				newline->whichtraj = evenstreamlines->ntrajs;
				newline->start = start;
				newline->end = end;
				QuadCell *curc = quadmesh->quadcells[pre_cell];

				//if it is major field tracing, add it to the major line info
				if(!majororminor)
				{
					if(curc->majorlines == NULL)
						curc->majorlines = new LineInfo(1);
					curc->majorlines->addNew(newline);
					curc->hasmajor=true;
				}

				//else, add to the minor line info
				else
				{
					if(curc->minorlines == NULL)
						curc->minorlines = new LineInfo(1);
					curc->minorlines->addNew(newline);
					curc->hasminor=true;
				}
				start = end = i;
				pre_cell = cur_traj->linesegs[i].Triangle_ID;
			}
		}

		/*we need to handle the last line segments*/
		start = end = cur_traj->nlinesegs-1;
		pre_cell = cur_traj->linesegs[start].Triangle_ID;
		for(i=cur_traj->nlinesegs-2; i>=0; i--)
		{
			if(cur_traj->linesegs[i].Triangle_ID == pre_cell)
			{
				start = i;
				continue;
			}
			else{
				LinesInOneCell *newline = (LinesInOneCell*)malloc(sizeof(LinesInOneCell));
				newline->whichtraj = evenstreamlines->ntrajs;
				newline->start = start;
				newline->end = end;
				QuadCell *curc = quadmesh->quadcells[pre_cell];

				//if it is major field tracing, add it to the major line info
				if(!majororminor)
				{
					if(curc->majorlines == NULL)
						curc->majorlines = new LineInfo(1);
					curc->majorlines->addNew(newline);
					curc->hasmajor=true;
				}

				//else, add to the minor line info
				else
				{
					if(curc->minorlines == NULL)
						curc->minorlines = new LineInfo(1);
					curc->minorlines->addNew(newline);
					curc->hasminor=true;
				}
				start = end = i;
				pre_cell = cur_traj->linesegs[i].Triangle_ID;
				break;
			}
		}

		/*set the corresponding road type as default one 10/03/2007*/
		if(majororminor)
			evenstreamlines->trajs[evenstreamlines->ntrajs]->roadtype = MINOR;
		else
			evenstreamlines->trajs[evenstreamlines->ntrajs]->roadtype = MINOR;

		evenstreamlines->ntrajs ++;

		return true;
	}

	return false;

}



void EvenStreamlinePlace::cal_seeds(int traj, double dsep, int every_nsample, int type, bool brushon)
{
	int i;
	double sample_p[2], newseed[2];
	int triangle, new_triangle;

	//Trace along the minor (major) field !!!
	icMatrix2x2 ten;
	double t[4]={0.};
	double evalues[2] = {0.};
	icVector2 ev[2], startdir;

	double origin_dsep=dsep;

	icVector2 tenline_dir;    /*  tensor line direction  */
	

	//Get the seed points
	for(i = 0; i < samplepts[traj]->nsamples; i++)
	{
		if(i % every_nsample != 0)
		{
			if(i != samplepts[traj]->nsamples-1) //we consider the last sample point
			    continue;
		}

		sample_p[0] = samplepts[traj]->samples[i]->gpt[0];
		sample_p[1] = samplepts[traj]->samples[i]->gpt[1];

		/*   obtain the direction of the tensor line segment  */
		if(i==0)
		{
			if(samplepts[traj]->nsamples>1)
			{
				tenline_dir.entry[0]=samplepts[traj]->samples[1]->gpt[0]-
					samplepts[traj]->samples[0]->gpt[0];
				tenline_dir.entry[1]=samplepts[traj]->samples[1]->gpt[1]-
					samplepts[traj]->samples[0]->gpt[1];
			}
			else
			{
				/*  use tensor line orientation  */
				tenline_dir.entry[0]=evenstreamlines->trajs[traj]->linesegs[evenstreamlines->trajs[traj]->nlinesegs-1].gstart[0]-
					evenstreamlines->trajs[traj]->linesegs[0].gend[0];
				tenline_dir.entry[1]=evenstreamlines->trajs[traj]->linesegs[evenstreamlines->trajs[traj]->nlinesegs-1].gstart[1]-
					evenstreamlines->trajs[traj]->linesegs[0].gend[1];
			}
		}
		else if(i==samplepts[traj]->nsamples-1)
		{
			if(samplepts[traj]->nsamples>1)
			{
				tenline_dir.entry[0]=samplepts[traj]->samples[i]->gpt[0]-
					samplepts[traj]->samples[samplepts[traj]->nsamples-1]->gpt[0];
				tenline_dir.entry[1]=samplepts[traj]->samples[i]->gpt[1]-
					samplepts[traj]->samples[samplepts[traj]->nsamples-1]->gpt[1];
			}
			else
			{
				/*  use tensor line orientation  */
				tenline_dir.entry[0]=evenstreamlines->trajs[traj]->linesegs[evenstreamlines->trajs[traj]->nlinesegs-1].gstart[0]-
					evenstreamlines->trajs[traj]->linesegs[0].gend[0];
				tenline_dir.entry[1]=evenstreamlines->trajs[traj]->linesegs[evenstreamlines->trajs[traj]->nlinesegs-1].gstart[1]-
					evenstreamlines->trajs[traj]->linesegs[0].gend[1];
			}
		}
		else
		{
			tenline_dir.entry[0]=samplepts[traj]->samples[i+1]->gpt[0]-
				samplepts[traj]->samples[i-1]->gpt[0];
			tenline_dir.entry[1]=samplepts[traj]->samples[i+1]->gpt[1]-
				samplepts[traj]->samples[i-1]->gpt[1];
		}

		triangle = samplepts[traj]->samples[i]->triangle;

		if(triangle < 0 || triangle >= quadmesh->nfaces)
			continue;

		compute_tensor_at_quad(triangle, sample_p[0], sample_p[1], ten);
		//get_tensor(sample_p[0], sample_p[1], t);
		//ten.entry[0][0]=t[0];
		//ten.entry[0][1]=t[1];
		//ten.entry[1][0]=t[2];
		//ten.entry[1][1]=t[3];

		cal_eigen_vector_sym(ten, ev);
		if(type == 0)  /*we place the tensor lines according to the major field*/
		{
			ev[0] = ev[1];  /*for choosing seeds, we use minor field*/
		}

		tenline_dir_global = startdir = ev[0];  /*obtain the major eigen vector*/

		/*  combined with the density map 11/25/2007  */
		if(sharedvars.CombinePopDensityOn)
		{
			double approx_den=cal_approx_density_at(sample_p, triangle);
			dsep = origin_dsep/approx_den;
		}

		////Get the forward seed, we scale the separate distance 4/18/06
		//if(get_a_seed(sample_p, triangle, traj, 1-type, newseed, new_triangle, dsep, brushon))
		if(get_a_seed(sample_p, tenline_dir, newseed, new_triangle, dsep, brushon))
		{
			if(new_triangle < 0 || new_triangle >= quadmesh->nfaces)
				continue;

			////Add to the seed point list

			/*use new method to create and link the element of the seeds*/
			Seed *s = (Seed *) malloc(sizeof(Seed));
			if(s == NULL)
			{
				return;
			}

			s->pos[0] = newseed[0];
			s->pos[1] = newseed[1];
			//s->pos[2] = 0;
			
			s->triangle = new_triangle;
			s->state = 0; //set it is active 
			s->weight=0;

			seedpts->append(s);
 		}
		
		////Get the backward seed, we scale the separate distance 
		tenline_dir_global = -startdir;

		tenline_dir=-tenline_dir;

		//if(get_a_seed(sample_p, triangle, traj, 1-type, newseed, new_triangle, dsep, brushon))
		if(get_a_seed(sample_p, tenline_dir, newseed, new_triangle, dsep, brushon))
		{
			if(new_triangle < 0 || new_triangle >= quadmesh->nfaces)
				continue;

			////Add to the seed point list

			/*use new method to create and link the element of the seeds*/
			Seed *s = (Seed *) malloc(sizeof(Seed));
			if(s == NULL)
			{
				return;
			}

			s->pos[0] = newseed[0];
			s->pos[1] = newseed[1];
			//s->pos[2] = newseed[2];
			
			s->triangle = new_triangle;
			s->state = 0; //set it is active 
			s->weight=0;

			seedpts->append(s);
		}

		////Extend the space for seeds if needed
		
		if(seedpts->isFull())
		{
			int oldnum = seedpts->curMaxNumSeeds;
			if(!seedpts->extend())
			{
				return;
			}

		}
	}
	
}


/*
   This may not work well in high curvature region
*/
bool EvenStreamlinePlace::get_a_seed(double sample_p[2], icVector2 line_dir,
			double end_p[2], int &end_triangle, double dsep, bool brushon)
{
	normalize(line_dir);

	icVector2 norm;
	norm.set(-line_dir.entry[1], line_dir.entry[0]);

	end_p[0]=sample_p[0]+dsep*norm.entry[0];
	end_p[1]=sample_p[1]+dsep*norm.entry[1];

	if(end_p[0]<quadmesh->xstart+1.e-8||end_p[0]>quadmesh->xend-1.e-8
		||end_p[1]<quadmesh->ystart+1.e-8||end_p[1]>quadmesh->yend-1.e-8)
		return false;

	end_triangle=get_cellID_givencoords(end_p[0], end_p[1]);

	if(end_triangle<0 || end_triangle>=quadmesh->nfaces)
		return false;

	return true;
}




bool EvenStreamlinePlace::get_a_seed(double sample_p[2], int begin_triangle, int cur_traj, int type,
			double end_p[2], int &end_triangle, double dsep, bool brushon)
{
	int i;
	int flag = -1;
	double globalp[2];
	int pre_face, cur_face;
	double pre_p[2] = {0.};
	double cur_p[2] = {0.};
	double cur_length = 0;
	double smallest_dist = dsep + 1.;
	int closest_traj = -1;

	pre_face = cur_face = begin_triangle;

	globalp[0] = sample_p[0];   
	globalp[1] = sample_p[1];

	tracing_points = (CurvePoints*) malloc(sizeof(CurvePoints) * 805);
	if(tracing_points == NULL)
	{
		return false;
	}

	num_tracingpoints = 0;

	int NUMTRACETRIS = (int)sqrt((double)quadmesh->nfaces);

	for(i = 0; i < NUMTRACETRIS; i++)
	{
		if(cur_face < 0 || cur_face >= quadmesh->nfaces)
		{
			free(tracing_points);
			//tracing_points=NULL;
			return false;
		}
			
		/*for brush interface 10/10/2007*/
		if(brushon && !quadmesh->quadcells[cur_face]->OnBoundary)
			break;
		
		///*for load the map 10/24/2007*/
		//if(is_not_inland(cur_face))
		//	break;
	
		/*it reaches any boundary, we should stop as well*/
		if(globalp[0]>quadmesh->xend-1.e-8||globalp[0]<quadmesh->xstart+1.e-8
			||globalp[1]>quadmesh->yend-1.e-8||globalp[1]<quadmesh->ystart+1.e-8)
			return false;

		pre_face = cur_face;
		cur_face = trace_in_triangle_seed_quad(cur_face, globalp, type, flag, 
			pre_p, cur_p, dsep, cur_length, cur_traj); ////0 means always forward tracing here


		if(flag == 2 || flag == 3 || cur_face < 0 ) //flag = 2--reach singularity;   flag = 3--reach maximum linesegments 
		{
			free(tracing_points);
			return false;
		}

		if(flag == 1) //reach the threshold, which means the length >= dsep
		{
			////Now we need to get the exact seed point
			double extra_length = cur_length - dsep;
			cal_exact_seed(cur_p, pre_p, cur_face, extra_length, end_p);

			end_triangle = pre_face;  //modified at 5/8/06

			if(end_triangle < 0)
			{
				free(tracing_points);
				return false;
			}

			free(tracing_points);
			return true;
		}
	}
}


double EvenStreamlinePlace::cal_approx_density_at(double p[2], int cell)
{
	/*we just need to call the bilinear interpolation routine, 
	since the density is scalar value*/
	QuadCell *qc = quadmesh->quadcells[cell];

	/*get the x coeff and y coeff*/
	double a = (p[0]-qc->x_start_coord)/quadmesh->xinterval;
	double b = (p[1]-qc->y_start_coord)/quadmesh->yinterval;

	if(fabs(a)<1e-6)
		a = 0;
	if(fabs(b)<1e-6)
		b = 0;

	/*obtain the vertices of this cell in the order of
	  v00 v01
	  v10 v11
	*/
	QuadVertex *v00 = quadmesh->quad_verts[qc->verts[0]];
	QuadVertex *v01 = quadmesh->quad_verts[qc->verts[3]];
	QuadVertex *v10 = quadmesh->quad_verts[qc->verts[1]];
	QuadVertex *v11 = quadmesh->quad_verts[qc->verts[2]];

	return bilinear_interpolate(a,b, v00->density,v01->density,v10->density,v11->density);
}



int EvenStreamlinePlace::trace_in_quad_even(int &face_id, double globalp[2], int type, 
					double dtest, double loopsep, double dist2sing, 
					double sample_interval, double discsize, int &flag)
{

	int i;
	double pre_point[2];

	double origin_dtest=dtest;
	double origin_dist2sing=dist2sing;
	double origin_loopsep=loopsep;

	/*  will this be a good solution? 1/9/2008 */
	if(!is_in_cell(face_id, globalp[0], globalp[1]))
	{
		face_id = get_cellID_givencoords(globalp[0], globalp[1]);
	}

	if(face_id < 0 || face_id>=quadmesh->nfaces)
		return -1;

	QuadCell *face = quadmesh->quadcells[face_id];

	QuadCell *pre_f = face;
	
	////Temporary curve point array

	CurvePoints *temp_point_list = (CurvePoints*) malloc(sizeof(CurvePoints) * 150);

	if(temp_point_list == NULL)
	{
		exit(-1);
	}

	int NumPoints = 0;
	
	/*the tracing will be performed under the global frame*/
	globalface = face_id;

	pre_point[0] = globalp[0];
	pre_point[1] = globalp[1];

	/*  can we let it go for the first point?  12/25/2007  */
			//temp_point_list[NumPoints].gpx = globalp[0];
			//temp_point_list[NumPoints].gpy = globalp[1];
			//temp_point_list[NumPoints].triangleid = face->index;  
			//NumPoints++;

			//get_nextpt_RK23_ten_quad(pre_point, globalp, face_id, type);

	////////////////////////////////////////////////////
    for(i = 0; i < 150; i++)
	{

		////2. if current point is inside current triangle
		if(is_in_cell(face_id, globalp[0], globalp[1]))
		{
			////store the point into the temp curve points list

			temp_point_list[NumPoints].gpx = globalp[0];
			temp_point_list[NumPoints].gpy = globalp[1];
			temp_point_list[NumPoints].triangleid = face->index;  
			NumPoints++;

			pre_point[0] = globalp[0];
			pre_point[1] = globalp[1];
			
			/*change to use other integration scheme 07/09/07*/
			//if(compute_next_pt_tensor_quad_global(pre_point, globalp, face_id))
			//if(get_nextpt_2ndeuler_ten_quad(pre_point, globalp, face_id, type))
			if(get_nextpt_RK23_ten_quad(pre_point, globalp, face_id, type))
			//if(get_nextpt_RK45_ten_quad(pre_point, globalp, face_id, type))
			{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = globalp[0] - pre_point[0];
				tenline_dir_global.entry[1] = globalp[1] - pre_point[1];

				/*      we need to combine the density map to change 
						the separation distance automatically
						11/25/2007
				*/
				if(sharedvars.CombinePopDensityOn)
				{
					double approx_den=cal_approx_density_at(globalp, face_id);
					dtest = origin_dtest/approx_den;
					dist2sing=origin_dist2sing/approx_den;
					loopsep=origin_loopsep/approx_den;
				}


				////using distance to judge whether a point close to current degenerate points
				if(close_to_degPt_except(globalp, face_id, -1, dist2sing, discsize))
				{
					flag = 1;

					//if(!StoreToGlobalList(temp_point_list, NumPoints))
					if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
						(temp_point_list, NumPoints))
					{
						////Not enough memory
						flag = 4;
						free(temp_point_list);
						return face_id;
					}

					free(temp_point_list);
				    return face_id;
				}

				////Judge whether it is too close to other existing streamlines
				//if(evenstreamlines->ntrajs > 0 &&
				//	close_to_cur_streamline(globalp, face_id, &evenstreamlines->ntrajs, 1, dtest, discsize, 0)) //scale the separate distance
				if(evenstreamlines->ntrajs > 0 &&
					close_to_cur_streamline(globalp, face_id, tenline_dir_global,
					&evenstreamlines->ntrajs, 1, dtest, discsize)) //scale the separate distance
				{
					//if(!StoreToGlobalList(temp_point_list, NumPoints))
					if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
						(temp_point_list, NumPoints))
					{
						////Not enough memory
						flag = 4;
						free(temp_point_list);
						return face_id;
					}

					free(temp_point_list);
					//flag = 3;
					flag = 2;

					return face_id;
				}

				////We may also need to compare the current point with the sampling point on itself!
				if(close_to_cur_samplePt(globalp, face_id, samplepts[evenstreamlines->ntrajs]->samples,
					samplepts[evenstreamlines->ntrajs]->nsamples, loopsep, discsize, sample_interval)) //scale the separate distance
				{
					//if(!StoreToGlobalList(temp_point_list, NumPoints))
					if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
						(temp_point_list, NumPoints))
					{
						////Not enough memory
						flag = 4;
						free(temp_point_list);
						return face_id;
					}

					flag = 3; //form a loop! 4/18/06

					free(temp_point_list);
					return face_id;
				}

			}

			else{  ////the curve reach a singularity/degenerate point
				flag = 1;

				////Store the record into global line segment array
                
				//if(!StoreToGlobalList(temp_point_list, NumPoints))
				if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
					(temp_point_list, NumPoints))
				{
					////Not enough memory
					flag = 4;
					free(temp_point_list);
					return face_id;
				}

				free(temp_point_list);

				return face_id;
			}
		}



		////3. if the point is out of current cell
		else{

			/*!!!!!!need to judge which cell it will enter!!!!!*/
			int PassVertornot = 0;
			//get_next_cell(face_id, pre_point, globalp, PassVertornot, type);
			get_next_cell_2(face_id, pre_point, globalp, PassVertornot, type);
			
			if(PassVertornot>0)  /*cross a vertex*/
			{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = pre_point[0] - globalp[0];
				tenline_dir_global.entry[1] = pre_point[1] - globalp[1];

				/**/
				temp_point_list[NumPoints].gpx = globalp[0];
				temp_point_list[NumPoints].gpy = globalp[1];
				temp_point_list[NumPoints].triangleid = face_id/*face->index*/;  ////cause problem 05/25/05
				NumPoints++;

				temp_point_list[NumPoints].gpx = pre_point[0];
				temp_point_list[NumPoints].gpy = pre_point[1];
				temp_point_list[NumPoints].triangleid = face_id;  ////cause problem 05/25/05
				NumPoints++;
			}
			else{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = globalp[0] - pre_point[0];
				tenline_dir_global.entry[1] = globalp[1] - pre_point[1];

				////Add the intersection point to the temporary points' list
				temp_point_list[NumPoints].gpx = globalp[0];
				temp_point_list[NumPoints].gpy = globalp[1];
				temp_point_list[NumPoints].triangleid = face->index;  ////cause problem 05/25/05
				NumPoints++;
			}

			//if(globalp[1] > 2)
			//{
			//	int test=0;
			//}

			/*obtain the global direction of current tensor line 09/20/2007*/
			//tenline_dir_global.entry[0] = globalp[0] - pre_point[0];
			//tenline_dir_global.entry[1] = globalp[1] - pre_point[1];

			//////Add the intersection point to the temporary points' list
			//temp_point_list[NumPoints].gpx = globalp[0];
			//temp_point_list[NumPoints].gpy = globalp[1];
			//temp_point_list[NumPoints].triangleid = face->index;  ////cause problem 05/25/05
			//NumPoints++;

			if(NumPoints > 1){
 				////Store the record into global line segment array
               //if(!StoreToGlobalList(temp_point_list, NumPoints))
				if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
					(temp_point_list, NumPoints))
			   {   ////Not enough memory
				   flag = 4;
				   free(temp_point_list);
				   return face_id;
			   }
			}

			free(temp_point_list);
			return face_id;
		}
	
	}

	if(NumPoints > 0)
		//StoreToGlobalList(temp_point_list, NumPoints);
		if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
						(temp_point_list, NumPoints))
						flag=4;

	free(temp_point_list);

	return face_id;
}


int EvenStreamlinePlace::trace_majRoad_in_quad(int &face_id, double globalp[2], int type, 
					double dtest, double loopsep, double dist2sing, 
					double sample_interval, double discsize, int &flag)
{

	int i;
	double pre_point[2];

	double origin_dtest=dtest;
	double origin_dist2sing=dist2sing;
	double origin_loopsep=loopsep;
	
	/*  will this be a good solution? 1/9/2008 */
	if(!is_in_cell(face_id, globalp[0], globalp[1]))
	{
		face_id = get_cellID_givencoords(globalp[0], globalp[1]);
	}

	if(face_id < 0 || face_id>=quadmesh->nfaces)
		return -1;

	QuadCell *face = quadmesh->quadcells[face_id];
	QuadCell *pre_f = face;

	
	////Temporary curve point array

	CurvePoints *temp_point_list = (CurvePoints*) malloc(sizeof(CurvePoints) * 200);

	if(temp_point_list == NULL)
	{
		exit(-1);
	}

	int NumPoints = 0;
	
	/*the tracing will be performed under the global frame*/
	globalface = face_id;

	pre_point[0] = globalp[0];
	pre_point[1] = globalp[1];

	////////////////////////////////////////////////////
    for(i = 0; i < 200; i++)
	{

		////2. if current point is inside current triangle
		if(is_in_cell(face_id, globalp[0], globalp[1]))
		{
			////store the point into the temp curve points list

			temp_point_list[NumPoints].gpx = globalp[0];
			temp_point_list[NumPoints].gpy = globalp[1];
			temp_point_list[NumPoints].triangleid = face->index;  
			NumPoints++;

			pre_point[0] = globalp[0];
			pre_point[1] = globalp[1];
			
			/*change to use other integration scheme 07/09/07*/
			if(get_nextpt_RK23_ten_quad(pre_point, globalp, face_id, type))
			//if(get_nextpt_2ndeuler_ten_quad(pre_point, globalp, face_id, type))
			{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = globalp[0] - pre_point[0];
				tenline_dir_global.entry[1] = globalp[1] - pre_point[1];

				/*      we need to combine the density map to change 
						the separation distance automatically
						11/25/2007
				*/
				if(sharedvars.CombinePopDensityOn)
				{
					double approx_den=cal_approx_density_at(globalp, face_id);
					dtest = origin_dtest/approx_den;
					dist2sing=origin_dist2sing/approx_den;
					loopsep=origin_loopsep/approx_den;
				}

				if(!sharedvars.AllowCrossSingularitiesOn && 
					close_to_degPt_except(globalp, face_id, -1, dist2sing, discsize))
				{
					flag = 1;

					//if(!StoreToGlobalList(temp_point_list, NumPoints))
					if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
						(temp_point_list, NumPoints))
					{
						////Not enough memory
						flag = 4;
						free(temp_point_list);
						return face_id;
					}

					free(temp_point_list);
				    return face_id;
				}

			////Judge whether it is too close to other existing streamlines
				if(evenstreamlines->ntrajs > 0 &&
					close_to_cur_streamline(globalp, face_id, &evenstreamlines->ntrajs, 1, dtest, discsize, 0)) //scale the separate distance
				{
					if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
						(temp_point_list, NumPoints))
					{
						////Not enough memory
						flag = 4;
						free(temp_point_list);
						return face_id;
					}

					flag = 2;
					free(temp_point_list);

			/* judge whether it is too close to existing major road  
			   it it is, connect them! 1/1/2008
			*/
					return face_id;
				}

				////We may also need to compare the current point with the sampling point on itself!
				if(close_to_cur_samplePt(globalp, face_id, samplepts[evenstreamlines->ntrajs]->samples,
					samplepts[evenstreamlines->ntrajs]->nsamples, loopsep, discsize, sample_interval)) //scale the separate distance
				{
					if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
						(temp_point_list, NumPoints))
					{
						////Not enough memory
						flag = 4;
						free(temp_point_list);
						return face_id;
					}

					flag = 3; //form a loop! 4/18/06

					free(temp_point_list);
					return face_id;
				}

			}

			else{  ////the curve reach a singularity/degenerate point
				flag = 1;

				////Store the record into global line segment array
                
				if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
					(temp_point_list, NumPoints))
				{
					////Not enough memory
					flag = 4;
					free(temp_point_list);
					return face_id;
				}

				free(temp_point_list);

				return face_id;
			}
		}



		////3. if the point is out of current cell
		else{

			/*!!!!!!need to judge which cell it will enter!!!!!*/
			int PassVertornot = 0;
			get_next_cell_2(face_id, pre_point, globalp, PassVertornot, type);

			//if(face_id==pre_f->index)
			//{
			//	int test=0;
			//}
			
			if(PassVertornot>0)  /*cross a vertex*/
			{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = pre_point[0] - globalp[0];
				tenline_dir_global.entry[1] = pre_point[1] - globalp[1];

				/**/
				temp_point_list[NumPoints].gpx = globalp[0];
				temp_point_list[NumPoints].gpy = globalp[1];
				temp_point_list[NumPoints].triangleid = face_id/*face->index*/;  ////cause problem 05/25/05
				NumPoints++;

				temp_point_list[NumPoints].gpx = pre_point[0];
				temp_point_list[NumPoints].gpy = pre_point[1];
				temp_point_list[NumPoints].triangleid = face_id;  ////cause problem 05/25/05
				NumPoints++;
			}
			else{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = globalp[0] - pre_point[0];
				tenline_dir_global.entry[1] = globalp[1] - pre_point[1];

				////Add the intersection point to the temporary points' list
				temp_point_list[NumPoints].gpx = globalp[0];
				temp_point_list[NumPoints].gpy = globalp[1];
				temp_point_list[NumPoints].triangleid = face->index;  ////cause problem 05/25/05
				NumPoints++;
			}


			if(NumPoints > 1){
 				////Store the record into global line segment array
				if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
					(temp_point_list, NumPoints))
			   {   ////Not enough memory
				   flag = 4;
				   free(temp_point_list);
				   return face_id;
			   }
			}

			free(temp_point_list);
			return face_id;
		}
	
	}

	if(NumPoints > 0)
		if(!evenstreamlines->trajs[evenstreamlines->ntrajs]->store_to_global_line_segs
						(temp_point_list, NumPoints))
						flag=4;

	free(temp_point_list);
	return face_id;
}





SamplePt **EvenStreamlinePlace::cal_samplepts_when_tracing(int traj, double interval, int &cur_line, int &movetonext, double &cur_length, 
							SamplePt **samples, int &num_samples)
{
	double cur_p[2] = {0.};
	QuadCell *face;

	while(cur_line < evenstreamlines->trajs[traj]->nlinesegs)
	{
		if(cal_a_sample_of_streamline(traj, cur_line, movetonext,
			cur_p, interval, cur_length))
		{

			FILE *fp;

			if(evenstreamlines->trajs[traj]->linesegs[cur_line].Triangle_ID<0
				||evenstreamlines->trajs[traj]->linesegs[cur_line].Triangle_ID>=quadmesh->nfaces)
				return samplepts[traj]->samples;
			
			face = quadmesh->quadcells[evenstreamlines->trajs[traj]->linesegs[cur_line].Triangle_ID];
			/*since we use global frame here, we need not convert to global frame here*/
			samplepts[traj]->samples[num_samples]->gpt[0] = cur_p[0];
			samplepts[traj]->samples[num_samples]->gpt[1] = cur_p[1];

			samplepts[traj]->samples[num_samples]->triangle = face->index;

			if(!is_in_reg_cell(face->index, cur_p[0], cur_p[1]))
			{
				int test = 0;

				samplepts[traj]->samples[num_samples]->triangle = 
					get_cellID_givencoords(cur_p[0], cur_p[1]);
			}
			
			if(samplepts[traj]->samples[num_samples]->triangle < 0 
				|| samplepts[traj]->samples[num_samples]->triangle >= quadmesh->nfaces)
			{
				continue;
			}

			num_samples++;

			if(num_samples >= samplepts[traj]->curMaxNumSamplePts)
			{

				int oldnum = samplepts[traj]->curMaxNumSamplePts;

					if(!samplepts[traj]->extend(200))
					{
						return NULL;
					}

				/*allocate memory for the elements*/
				for(int i = oldnum ; i < samplepts[traj]->curMaxNumSamplePts; i++)
				{
					samplepts[traj]->samples[i] = (SamplePt *)malloc(sizeof(SamplePt));
				}
				
			}
		}
	}

	cur_line--;

	if(cur_line < 0)
		cur_line = 0;

	return samplepts[traj]->samples;
}


void EvenStreamlinePlace::reverse_streamline(int streamlineid)
{
	////
	int i;
	int num_lines = evenstreamlines->trajs[streamlineid]->nlinesegs;
	Trajectory *traj = evenstreamlines->trajs[streamlineid];

	LineSeg *temp = (LineSeg *)malloc(sizeof(LineSeg)*(num_lines+1));
		
	if(temp == NULL)
	{
		return;
	}

	
	int newnum_lines = 0;

	////store the line segment in reversed order

	for(i = num_lines-1; i >= 0; i--)
	{
		if(traj->linesegs[i].Triangle_ID < 0
			|| traj->linesegs[i].Triangle_ID >= quadmesh->nfaces  
			|| traj->linesegs[i].length < 0)
		{
			continue;
		}

		temp[newnum_lines].gstart[0] = traj->linesegs[i].gend[0];
		temp[newnum_lines].gstart[1] = traj->linesegs[i].gend[1];
		
		temp[newnum_lines].gend[0] = traj->linesegs[i].gstart[0];
		temp[newnum_lines].gend[1] = traj->linesegs[i].gstart[1];
		
		temp[newnum_lines].start[0] = traj->linesegs[i].end[0];
		temp[newnum_lines].start[1] = traj->linesegs[i].end[1];

		temp[newnum_lines].end[0] = traj->linesegs[i].start[0];
		temp[newnum_lines].end[1] = traj->linesegs[i].start[1];

		temp[newnum_lines].length = traj->linesegs[i].length;
		temp[newnum_lines].Triangle_ID = traj->linesegs[i].Triangle_ID;


		newnum_lines++;
	}

	////Copy it back to the origin array
	for(i = 0; i < newnum_lines; i++)
	{
		traj->linesegs[i].gstart[0] = temp[i].gstart[0];
		traj->linesegs[i].gstart[1] = temp[i].gstart[1];
		
		traj->linesegs[i].gend[0] = temp[i].gend[0];
		traj->linesegs[i].gend[1] = temp[i].gend[1];

		traj->linesegs[i].start[0] = temp[i].start[0];
		traj->linesegs[i].start[1] = temp[i].start[1];
		
		traj->linesegs[i].end[0] = temp[i].end[0];
		traj->linesegs[i].end[1] = temp[i].end[1];

		traj->linesegs[i].length = temp[i].length;
		traj->linesegs[i].Triangle_ID = temp[i].Triangle_ID;
	}

	traj->nlinesegs = newnum_lines;

	free(temp);
}


void EvenStreamlinePlace::reorder_streamline(int streamlineid)
{
	////
	int i;
	int num_lines = evenstreamlines->trajs[streamlineid]->nlinesegs;
	Trajectory *traj = evenstreamlines->trajs[streamlineid];
	
	icVector2 globaldir, linedir;
	double re;

	LineSeg *temp = (LineSeg *)malloc(sizeof(LineSeg)*(num_lines+1));
		
	if(temp == NULL)
	{
		return;
	}

	
	int newnum_lines = 0;

	globaldir.entry[0]=traj->linesegs[0].gend[0]-traj->linesegs[0].gstart[0];
	globaldir.entry[1]=traj->linesegs[0].gend[1]-traj->linesegs[0].gstart[1];

	////store the line segment in reversed order

	for(i = 0; i < num_lines; i++)
	{
		if(traj->linesegs[i].Triangle_ID < 0
			|| traj->linesegs[i].Triangle_ID >= quadmesh->nfaces  
			|| traj->linesegs[i].length < 0)
		{
			continue;
		}

		linedir.entry[0]=traj->linesegs[i].gend[0]-traj->linesegs[0].gstart[0];
		linedir.entry[1]=traj->linesegs[i].gend[1]-traj->linesegs[0].gstart[1];

		if(dot(globaldir, linedir)>=0)
		{
			temp[newnum_lines].gstart[0] = traj->linesegs[i].gstart[0];
			temp[newnum_lines].gstart[1] = traj->linesegs[i].gstart[1];
			
			temp[newnum_lines].gend[0] = traj->linesegs[i].gend[0];
			temp[newnum_lines].gend[1] = traj->linesegs[i].gend[1];
			
			temp[newnum_lines].start[0] = traj->linesegs[i].start[0];
			temp[newnum_lines].start[1] = traj->linesegs[i].start[1];

			temp[newnum_lines].end[0] = traj->linesegs[i].end[0];
			temp[newnum_lines].end[1] = traj->linesegs[i].end[1];
			
			globaldir=linedir;
		}
		else
		{
			temp[newnum_lines].gstart[0] = traj->linesegs[i].gend[0];
			temp[newnum_lines].gstart[1] = traj->linesegs[i].gend[1];
			
			temp[newnum_lines].gend[0] = traj->linesegs[i].gstart[0];
			temp[newnum_lines].gend[1] = traj->linesegs[i].gstart[1];
			
			temp[newnum_lines].start[0] = traj->linesegs[i].end[0];
			temp[newnum_lines].start[1] = traj->linesegs[i].end[1];

			temp[newnum_lines].end[0] = traj->linesegs[i].start[0];
			temp[newnum_lines].end[1] = traj->linesegs[i].start[1];
			
			globaldir=-linedir;
		}

			
		temp[newnum_lines].length = traj->linesegs[i].length;
		temp[newnum_lines].Triangle_ID = traj->linesegs[i].Triangle_ID;


		newnum_lines++;
	}

	////Copy it back to the origin array
	for(i = 0; i < newnum_lines; i++)
	{
		traj->linesegs[i].gstart[0] = temp[i].gstart[0];
		traj->linesegs[i].gstart[1] = temp[i].gstart[1];
		
		traj->linesegs[i].gend[0] = temp[i].gend[0];
		traj->linesegs[i].gend[1] = temp[i].gend[1];

		traj->linesegs[i].start[0] = temp[i].start[0];
		traj->linesegs[i].start[1] = temp[i].start[1];
		
		traj->linesegs[i].end[0] = temp[i].end[0];
		traj->linesegs[i].end[1] = temp[i].end[1];

		traj->linesegs[i].length = temp[i].length;
		traj->linesegs[i].Triangle_ID = temp[i].Triangle_ID;
	}

	traj->nlinesegs = newnum_lines;

	free(temp);
}


void EvenStreamlinePlace::update_samples_in_cell(int traj, SamplePt **samples, int num_samples)
{
	int i;

	for(i = 0; i < num_samples; i++)
	{
		add_sample_to_cell(samples[i]->triangle, traj, i, majororminor);
	}
}


int EvenStreamlinePlace::trace_in_triangle_seed_quad(int &face_id, double globalp[3], int type, 
							int &flag, double pre_p[3], double cur_p[3], 
							double dsep, double &cur_length, int cur_traj)
{
	if(face_id < 0)
		return -1;

	int i;
	double pre_point[2];
	icVector2 dist;

	QuadCell *face = quadmesh->quadcells[face_id];

	QuadCell *pre_f = face;


	double smallest_dist = dsep+1.;
	int closest_traj = -1;
	
	////initialize
	globalface = face_id;

	pre_point[0] = globalp[0];
	pre_point[1] = globalp[1];


	////////////////////////////////////////////////////
    for(i = 0; i < 100; i++)
	{
		if(is_in_cell(face_id, globalp[0], globalp[1]))
		{
			////store the point into the temp curve points list

			tracing_points[num_tracingpoints].gpx = globalp[0];
			tracing_points[num_tracingpoints].gpy = globalp[1];
			tracing_points[num_tracingpoints].triangleid = face->index;  

			////Get the length for each line segment
			dist.entry[0] = globalp[0] - pre_point[0];
			dist.entry[1] = globalp[1] - pre_point[1];
			cur_length += tracing_points[num_tracingpoints].length = length(dist); //sum the lengths
            num_tracingpoints++;

			if(cur_length > dsep) //if the total length of the tracing curve is larger than the threshold
			{
				flag = 1;
				//globalv = pre_point[0] * face->LX + pre_point[1] * face->LY;
				pre_p[0] = pre_point[0];
				pre_p[1] = pre_point[1];
				//pre_p[2] = vert0[2] + globalv.entry[2];
				//
				//globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;
				cur_p[0] = globalp[0];
				cur_p[1] = globalp[1];
				//cur_p[2] = vert0[2] + globalv.entry[2];
				return face_id;
			}

			////Testing codes 4/17/06
			if(num_tracingpoints>=800) //can not get enough length
			{
				flag = 3;

				return face_id;
			}

			pre_point[0] = globalp[0];
			pre_point[1] = globalp[1];

			/*get next point*/
			//if(compute_next_pt_tensor_quad_global(pre_point, globalp, face_id))
			//if(get_nextpt_2ndeuler_ten_quad(pre_point, globalp, face_id, type))
			if(get_nextpt_RK23_ten_quad(pre_point, globalp, face_id, type))
			//if(get_nextpt_RK45_ten_quad(pre_point, globalp, face_id, type))
			{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = globalp[0] - pre_point[0];
				tenline_dir_global.entry[1] = globalp[1] - pre_point[1];
			}

			else{  ////the curve reach a singularity/degenerate point
				flag = 2;
				return face_id;
			}
		}

		////3. if the point is out of current triangle
		else{

			/*!!!!!!need to judge which cell it will enter!!!!!*/
			int PassVertornot = 0;
			double t_p[2]={pre_point[0], pre_point[1]};

			//get_next_cell(face_id, pre_point, globalp, PassVertornot, type);
			get_next_cell_2(face_id, pre_point, globalp, PassVertornot, type);

			if(PassVertornot>0)  /*cross a vertex*/
			{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = pre_point[0] - globalp[0];
				tenline_dir_global.entry[1] = pre_point[1] - globalp[1];

				/**/
				tracing_points[num_tracingpoints].gpx = globalp[0];
				tracing_points[num_tracingpoints].gpy = globalp[1];
				tracing_points[num_tracingpoints].lpx = globalp[0];
				tracing_points[num_tracingpoints].lpy = globalp[1];
				tracing_points[num_tracingpoints].triangleid = face->index;  ////cause problem 05/25/05
				num_tracingpoints++;

				tracing_points[num_tracingpoints].gpx = pre_point[0];
				tracing_points[num_tracingpoints].gpy = pre_point[1];
				tracing_points[num_tracingpoints].lpx = pre_point[0];
				tracing_points[num_tracingpoints].lpy = pre_point[1];
				tracing_points[num_tracingpoints].triangleid = face_id;  ////cause problem 05/25/05
				num_tracingpoints++;
				
				//////Get the length for each line segment
				dist.entry[0] = t_p[0] - pre_point[0];
				dist.entry[1] = t_p[1] - pre_point[1];
				cur_length += tracing_points[num_tracingpoints].length = length(dist);
				num_tracingpoints++;
			}
			else{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = globalp[0] - pre_point[0];
				tenline_dir_global.entry[1] = globalp[1] - pre_point[1];

				////Add the intersection point to the temporary points' list
				tracing_points[num_tracingpoints].gpx = globalp[0];
				tracing_points[num_tracingpoints].gpy = globalp[1];
				tracing_points[num_tracingpoints].lpx = globalp[0];
				tracing_points[num_tracingpoints].lpy = globalp[1];
				tracing_points[num_tracingpoints].triangleid = face->index;  ////cause problem 05/25/05
				num_tracingpoints++;

				//////Get the length for each line segment
				dist.entry[0] = globalp[0] - pre_point[0];
				dist.entry[1] = globalp[1] - pre_point[1];
				cur_length += tracing_points[num_tracingpoints].length = length(dist);
				num_tracingpoints++;
			}

			/*obtain the global direction of current tensor line 09/20/2007*/
			//tenline_dir_global.entry[0] = globalp[0] - pre_point[0];
			//tenline_dir_global.entry[1] = globalp[1] - pre_point[1];
   //     
			//////Add the intersection point to the temporary points' list
			//tracing_points[num_tracingpoints].gpx = globalp[0];
			//tracing_points[num_tracingpoints].gpy = globalp[1];
			//tracing_points[num_tracingpoints].lpx = globalp[0];
			//tracing_points[num_tracingpoints].lpy = globalp[1];
			//tracing_points[num_tracingpoints].triangleid = pre_f->index;  
			//

			////if the length of the tracing curve is larger than the threshold
			if(cur_length > dsep)
			{
				flag = 1;

				if(PassVertornot>0)
				{
					pre_p[0] = t_p[0];
					pre_p[1] = t_p[1];
					
					cur_p[0] = pre_point[0];
					cur_p[1] = pre_point[1];
				}
				else
				{
					pre_p[0] = pre_point[0];
					pre_p[1] = pre_point[1];

					cur_p[0] = globalp[0];
					cur_p[1] = globalp[1];
				}
			}
			
			////Testing codes 4/17/06
			if(num_tracingpoints>=800) //can not get enough length
			{
				flag = 3;

				return face_id;
			}

			return face_id;
		}
	}

	flag = 3;   //can not go outside of the triangle

	return face_id;
}


void EvenStreamlinePlace::cal_exact_seed(double cur_p[2], double pre_p[2], int triangle, 
						double extra_length, double exact_seed[2])
{
	icVector3 line;
	line.entry[0] = cur_p[0] - pre_p[0];
	line.entry[1] = cur_p[1] - pre_p[1];

	double t = 1-(extra_length/length(line)); //it seems that this can get more even separate distance

	exact_seed[0] = (1-t) * pre_p[0] + t * cur_p[0];
	exact_seed[1] = (1-t) * pre_p[1] + t * cur_p[1];

}


bool EvenStreamlinePlace::close_to_degPt_except(double p[3], int triangle, 
												  int singid, double threshold, double discsize)
{
	int i;
	int *NearbyTriangles = NULL;
	int num_triangles = 0;
	QuadCell *face;
	//QuadVertex *v;
	icVector3 VP/*, v0*/;
	double dist/*, a, b, alpha[3]*/;

	////For unfolding
	int edgeid;
	double rotmat[16] = {0.};
	double pt[2];
	//set_Indentity16(rotmat);

	////Get the triangles that fall into the circle region with p as the center and 3*separate_dist as radius
	//NearbyTriangles = cal_euclidean_dist(triangle, p, threshold, discsize, NearbyTriangles, num_triangles);
	
		trianglelist = new DynList_Int();
	trianglelist->nelems = 0;
	cal_euclidean_dist_2(triangle, p, threshold, discsize, trianglelist);
	NearbyTriangles = trianglelist->elems;
	num_triangles = trianglelist->nelems;


	for(i = 0; i < num_triangles; i++)
	{
		face = quadmesh->quadcells[NearbyTriangles[i]];

		if(face->degpt_index < 0) //this triangle does not contain fixed point
			continue;

		if(face->degpt_index == singid)
			continue;

		/*   we allow the tensor lines go through centers and nodes */
		if(ten_designelems[face->degpt_index].type == 2 /*node*/
			|| ten_designelems[face->degpt_index].type == 3 /*center*/)
			continue;
		
		pt[0]= degpts[face->degpt_index].gcx;
		pt[1]= degpts[face->degpt_index].gcy;
			
		VP.entry[0] = pt[0] - p[0];
		VP.entry[1] = pt[1] - p[1];

		dist = length(VP); /*compute the Euclidean distance*/

		if(dist <= threshold)
		{
			//reset_dist(NearbyTriangles,num_triangles);
			//free(NearbyTriangles);

			reset_dist(trianglelist);
			delete trianglelist;
			return true;
		}
	}
			
	//reset_dist(NearbyTriangles,num_triangles);
	//free(NearbyTriangles);

	reset_dist(trianglelist);
	delete trianglelist;
	return false;
}


bool EvenStreamlinePlace::close_to_cur_streamline(double p[3], int triangle, 
						int *Except_trajs, int num_trajs, 
						double separate_dist, double discsize, int Init)
						//int &which_triangle, double samp[2])
{
	int *NearbyTriangles = NULL;
	int num_triangles = 0;

	int i, j, k;
	int traj;
	QuadCell *face;
	//QuadVertex *v1, *v2, *v3;
	double dis/*, alpha[3]*/;
	int sampleid;

	icVector2 VP;

	double pt[2] = {0.};

	trianglelist = new DynList_Int();
	trianglelist->nelems = 0;
	cal_euclidean_dist_2(triangle, p, separate_dist, discsize, trianglelist);
	NearbyTriangles = trianglelist->elems;
	num_triangles = trianglelist->nelems;

	int nsamplepts;
	SampleListInTriangle *samplelist;

	//////then we test all the curve points in the nearby triangles
	for(i = 0; i < num_triangles; i++)
	{
		face = quadmesh->quadcells[NearbyTriangles[i]];

		if(!majororminor) /*trace major*/
		{
			nsamplepts=face->maj_nsamplepts;
			samplelist=face->maj_samplepts;
		}
		else /*trace minor*/
		{
			nsamplepts=face->min_nsamplepts;
			samplelist=face->min_samplepts;
		}

		//if(face->num_samplepts <= 0)

		if(nsamplepts <=0)
			continue;

		/////////

		//for(j = 0; j < face->num_samplepts; j++)
		for(j = 0; j < nsamplepts; j++)
		{
			//traj = face->samplepts[j].which_traj;
			traj = samplelist[j].which_traj;

			if(is_repeated_elem(Except_trajs, traj, num_trajs))
				continue;

			/*   we allow it close to the major or highway 12/26/2007 */
			if(sharedvars.AllowMinCloseToMajOn&&
				(evenstreamlines->trajs[traj]->roadtype==MAJOR
				||evenstreamlines->trajs[traj]->roadtype==HIGHWAY))
				continue;
			
			//sampleid = face->samplepts[j].which_sample;
			sampleid = samplelist[j].which_sample;
		
			////For the center 4 triangles we unfold them and use the Euclidean distance directly

			pt[0] = samplepts[traj]->samples[sampleid]->gpt[0];
			pt[1] = samplepts[traj]->samples[sampleid]->gpt[1];
			
			VP.entry[0] = pt[0] - p[0];
			VP.entry[1] = pt[1] - p[1];

			dis = length(VP);  /*calculate the Euclidean distance only*/

			if(dis <= separate_dist) //scale the distance a little bit
			{
				reset_dist(trianglelist);
				delete trianglelist;

				/*we may need to find out a point that the connection
				is perpandicular to both tensor lines 11/17/2007*/
				int cur_cellid=face->index;
				//if(cal_approx_perpendicular_pt(traj, cur_cellid, p, samp, majororminor))
				//	which_triangle=cur_cellid;
				//else{
				//	which_triangle=face->index;
				//	samp[0]=pt[0];
				//	samp[1]=pt[1];
				//}

				return true;
			}
		}
	}

	if(NearbyTriangles != NULL)
	{
		reset_dist(trianglelist);
		delete trianglelist;
	}
	return false;
}



/*
    In this version of the overloaded routine, 
	we allow current tensor line to go as close to the existing major roads (or highways) as possible
	but we need to judge whether it is parallel to the nearby  major roads (or highways),
	if it is, we stop the tracing; otherwise, keep going :)
*/
bool EvenStreamlinePlace::close_to_cur_streamline(double p[2], int triangle, icVector2 go_dir,
						int *Except_trajs, int num_trajs, 
						double separate_dist, double discsize/*, int Init*/)
{
	int *NearbyTriangles = NULL;
	int num_triangles = 0;

	int i, j, k;
	int traj;
	QuadCell *face;
	double dis;
	int sampleid;

	icVector2 VP;

	double pt[2] = {0.};

	trianglelist = new DynList_Int();
	trianglelist->nelems = 0;
	cal_euclidean_dist_2(triangle, p, separate_dist, discsize, trianglelist);
	NearbyTriangles = trianglelist->elems;
	num_triangles = trianglelist->nelems;

	int nsamplepts;
	SampleListInTriangle *samplelist;
	normalize(go_dir);

	//////then we test all the curve points in the nearby triangles
	for(i = 0; i < num_triangles; i++)
	{
		face = quadmesh->quadcells[NearbyTriangles[i]];

		if(!majororminor) /* trace major direction */
		{
			nsamplepts=face->maj_nsamplepts;
			samplelist=face->maj_samplepts;
		}
		else              /* trace minor direction */
		{
			nsamplepts=face->min_nsamplepts;
			samplelist=face->min_samplepts;
		}

		if(nsamplepts <=0)
			continue;

		/////////

		for(j = 0; j < nsamplepts; j++)
		{
			traj = samplelist[j].which_traj;

			if(is_repeated_elem(Except_trajs, traj, num_trajs))
				continue;

			sampleid = samplelist[j].which_sample;
			pt[0] = samplepts[traj]->samples[sampleid]->gpt[0];
			pt[1] = samplepts[traj]->samples[sampleid]->gpt[1];

			/*   we allow it close to the major or highway 12/26/2007 */
			if(evenstreamlines->trajs[traj]->roadtype==MAJOR
				||evenstreamlines->trajs[traj]->roadtype==HIGHWAY)
			{
				/*judge whether the tensor line goes parallely to the major road (or highway)*/

				icVector2 maj_dir;

				/*  first, obtain the local major road direction maj_dir  */
				if(samplepts[traj]->nsamples==1)
					continue;

				if(sampleid==0)
				{
					maj_dir.entry[0]=samplepts[traj]->samples[sampleid+1]->gpt[0]-pt[0];
					maj_dir.entry[1]=samplepts[traj]->samples[sampleid+1]->gpt[1]-pt[1];
				}
				else
				{
					maj_dir.entry[0]=pt[0]-samplepts[traj]->samples[sampleid-1]->gpt[0];
					maj_dir.entry[1]=pt[1]-samplepts[traj]->samples[sampleid-1]->gpt[1];
				}
				normalize(maj_dir);

				if(fabs(dot(maj_dir, go_dir))<0.2)
				continue;
			}
			
		
			////For the center 4 triangles we unfold them and use the Euclidean distance directly

			
			VP.entry[0] = pt[0] - p[0];
			VP.entry[1] = pt[1] - p[1];

			dis = length(VP);  /*calculate the Euclidean distance only :)*/

			if(dis <= separate_dist) //scale the distance a little bit
			{
				reset_dist(trianglelist);
				delete trianglelist;

				/*we may need to find out a point that the connection
				is perpandicular to both tensor lines 11/17/2007*/
				int cur_cellid=face->index;
				return true;
			}
		}

		/*  for minor direction tracing with major road existing, we also need to judge whether 
		    it is too close and parallel to the nearby major road (or highway)! 12/26/2007
		*/
		if(majororminor)/*  trace minor direction  */
		{
			nsamplepts=face->maj_nsamplepts;
			samplelist=face->maj_samplepts;
		}
		else            /*  trace major direction  */
		{
			nsamplepts=face->min_nsamplepts;
			samplelist=face->min_samplepts;
		}

		for(j = 0; j < nsamplepts; j++)
		{
			traj = samplelist[j].which_traj;

			//if(is_repeated_elem(Except_trajs, traj, num_trajs))
			//	continue;

			if(evenstreamlines->trajs[traj]->roadtype!=MAJOR
				&&evenstreamlines->trajs[traj]->roadtype!=HIGHWAY)
				continue;

			sampleid = samplelist[j].which_sample;
			pt[0] = samplepts[traj]->samples[sampleid]->gpt[0];
			pt[1] = samplepts[traj]->samples[sampleid]->gpt[1];
			
			VP.entry[0] = pt[0] - p[0];
			VP.entry[1] = pt[1] - p[1];
			dis = length(VP);  /*calculate the Euclidean distance only :)*/

			if(dis > separate_dist/*/2.*/) //scale the distance a little bit
				continue;


			/*   we allow it close to the major or highway 12/26/2007 */
			/*judge whether the tensor line goes parallely to the major road (or highway)*/
			icVector2 maj_dir;

			/*  first, obtain the local major road direction maj_dir  */
			if(samplepts[traj]->nsamples==1)
				continue;

			if(sampleid==0)
			{
				maj_dir.entry[0]=samplepts[traj]->samples[sampleid+1]->gpt[0]-pt[0];
				maj_dir.entry[1]=samplepts[traj]->samples[sampleid+1]->gpt[1]-pt[1];
			}
			else
			{
				maj_dir.entry[0]=pt[0]-samplepts[traj]->samples[sampleid-1]->gpt[0];
				maj_dir.entry[1]=pt[1]-samplepts[traj]->samples[sampleid-1]->gpt[1];
			}
			normalize(maj_dir);

			if(fabs(dot(maj_dir, go_dir))<0.2)
				continue;
			return true;
		}
	}

	if(NearbyTriangles != NULL)
	{
		reset_dist(trianglelist);
		delete trianglelist;
	}
	return false;
}


/*
   For major road tracing.
   Judge whether current major road is too close to existing roads
*/

bool EvenStreamlinePlace::close_to_cur_streamline_maj(double p[3], int triangle, 
						int *Except_trajs, int num_trajs, 
						double separate_dist, double discsize, int &which_traj, int &which_samp)
{
	int *NearbyTriangles = NULL;
	int num_triangles = 0;

	int i, j, k;
	int traj;
	QuadCell *face;
	double dis;
	int sampleid;

	icVector2 VP;

	double pt[2] = {0.};

	trianglelist = new DynList_Int();
	trianglelist->nelems = 0;
	cal_euclidean_dist_2(triangle, p, separate_dist, discsize, trianglelist);
	NearbyTriangles = trianglelist->elems;
	num_triangles = trianglelist->nelems;

	int nsamplepts;
	SampleListInTriangle *samplelist;

	double smallest_dist=1.e50;
	bool found=false;

	//////then we test all the curve points in the nearby triangles
	for(i = 0; i < num_triangles; i++)
	{
		face = quadmesh->quadcells[NearbyTriangles[i]];

		if(!majororminor) /*trace major*/
		{
			nsamplepts=face->maj_nsamplepts;
			samplelist=face->maj_samplepts;
		}
		else /*trace minor*/
		{
			nsamplepts=face->min_nsamplepts;
			samplelist=face->min_samplepts;
		}

		if(nsamplepts <=0)
			continue;

		/////////

		for(j = 0; j < nsamplepts; j++)
		{
			traj = samplelist[j].which_traj;

			if(is_repeated_elem(Except_trajs, traj, num_trajs))
				continue;

			/*   we allow it close to the major or highway 12/26/2007 */
			//if(evenstreamlines->trajs[traj]->roadtype==MAJOR
			//	||evenstreamlines->trajs[traj]->roadtype==HIGHWAY)
			//	continue;
			
			sampleid = samplelist[j].which_sample;
		
			////For the center 4 triangles we unfold them and use the Euclidean distance directly

			pt[0] = samplepts[traj]->samples[sampleid]->gpt[0];
			pt[1] = samplepts[traj]->samples[sampleid]->gpt[1];
			
			VP.entry[0] = pt[0] - p[0];
			VP.entry[1] = pt[1] - p[1];

			dis = length(VP);  /*calculate the Euclidean distance only*/

			if(dis <= separate_dist) 
			{
				//reset_dist(trianglelist);
				//delete trianglelist;
				found=true;
				if(dis<smallest_dist)
				{
					which_traj=traj;
					which_samp=sampleid;
				}

				/*we may need to find out a point that the connection
				is perpandicular to both tensor lines 11/17/2007*/
				int cur_cellid=face->index;
			}
		}
	}

	if(NearbyTriangles != NULL)
	{
		reset_dist(trianglelist);
		delete trianglelist;
	}

	return found;
}




/*NOTE that we already have the information of the line segments of the existing
tensor lines for each triangle
*/
bool EvenStreamlinePlace::cal_approx_perpendicular_pt(int traj, int &cellid, double p[2],
													  double pt[2], bool type)
{
	QuadCell *face=quadmesh->quadcells[cellid];
	int i;
	//face->majorlines->lines[0]->
	LineInfo *lineinfo;
	int startlineseg, endlineseg;
	Trajectory *tt=evenstreamlines->trajs[traj];

	if(!type) lineinfo=face->majorlines;
	else lineinfo=face->minorlines;
	if(lineinfo == NULL) return false;
	for(i=0; i<lineinfo->nlines; i++)
	{
		if(lineinfo->lines[i]->whichtraj==traj)
		{
			startlineseg=lineinfo->lines[i]->start;
			endlineseg=lineinfo->lines[i]->end;
			break;
		}

	}

	if(i>=lineinfo->nlines) return false;

	double distolineseg, distoline;
	DistanceFromLine(p[0],p[1],tt->linesegs[startlineseg].gstart[0], tt->linesegs[startlineseg].gstart[1],
		tt->linesegs[startlineseg].gend[0],tt->linesegs[startlineseg].gend[1],distolineseg, distoline);

	double smallestdis=distolineseg;
	int smallestlineseg=startlineseg;

	//void DistanceFromLine(double cx, double cy, double ax, double ay ,
	//				  double bx, double by, double &distanceSegment,
	//				  double &distanceLine);

	/*find out a point on the trajectory who will 
	provide us an approximate perpendicular connection*/
	for(i=max(0, startlineseg-20); i<=min(endlineseg+20, evenstreamlines->trajs[traj]->nlinesegs-1);
		i++)
	{
		DistanceFromLine(p[0],p[1],tt->linesegs[i].gstart[0],tt->linesegs[i].gstart[1],
			tt->linesegs[i].gend[0],tt->linesegs[i].gend[1],distolineseg, distoline);

		if(distolineseg<smallestdis)
		{
			smallestdis=distolineseg;
			smallestlineseg=i;
		}
	}

	/*use the start point of the line segment "smallestlineseg" as the return point*/
	pt[0]=(tt->linesegs[smallestlineseg].gend[0]+tt->linesegs[smallestlineseg].gstart[0])/2;
	pt[1]=(tt->linesegs[smallestlineseg].gend[1]+tt->linesegs[smallestlineseg].gstart[1])/2;
	cellid=tt->linesegs[smallestlineseg].Triangle_ID;

	return true;
}



bool EvenStreamlinePlace::close_to_cur_samplePt(double p[2], int triangle, SamplePt **samples, int num_samples,
						double separate_dist, double discsize, double sample_interval)
{
	int i;
	double dist/*, alpha[3], a, b*/;
	int stop_locate = floor(separate_dist/sample_interval)+3;
	QuadCell *face;
	QuadVertex *v;
	icVector2 VP;

	int *NearbyTriangles = NULL;
	int num_triangles = 0;
	

	////Get the triangles that fall into the circle region with p as the center and 3*separate_dist as radius
	trianglelist = new DynList_Int();

	trianglelist->nelems = 0;
	cal_euclidean_dist_2(triangle, p, separate_dist, discsize, trianglelist);
	NearbyTriangles = trianglelist->elems;
	num_triangles = trianglelist->nelems;


	for(i = 0 ; i < num_samples-stop_locate; i++)
	{
		if(!is_repeated_elem(NearbyTriangles, samples[i]->triangle, num_triangles))
			continue;

		face = quadmesh->quadcells[samples[i]->triangle];
		
		VP.entry[0] = p[0]-samples[i]->gpt[0];
		VP.entry[1] = p[1]-samples[i]->gpt[1];

		dist = length(VP);  /*calculate the Euclidean distance*/

		if(dist < separate_dist)
		{

			reset_dist(trianglelist);
			delete trianglelist;

			return true;
		}
	}

	reset_dist(trianglelist);
	delete trianglelist;
	return false;
}


void EvenStreamlinePlace::reset_dist(int *NearbyTriangles, int num)
{
	int i, j;
	QuadCell *face;
	QuadVertex *v;
	for(i = 0; i < num; i++)
	{
		face = quadmesh->quadcells[NearbyTriangles[i]];
		face->visited = false;
		for(j = 0; j < 3; j++)
		{
			v = quadmesh->quad_verts[face->verts[j]];
			v->distance = 1e50;
		}
	}
}



void EvenStreamlinePlace::reset_dist(DynList_Int *trianglelist)
{
	int i, j;
	QuadCell *face;
	QuadVertex *v;
	for(i = 0; i < trianglelist->nelems; i++)
	{
		face = quadmesh->quadcells[trianglelist->elems[i]];
		face->visited = false;
		for(j = 0; j < 3; j++)
		{
			v = quadmesh->quad_verts[face->verts[j]];
			v->distance = 1e50;
		}
	}
}


bool EvenStreamlinePlace::cal_a_sample_of_streamline(int traj, int &cur_lineindex, int &movetonext,
					double curpt[2], double interval, double &cur_length)
{
	int i;
	icVector2 len_vec;
	double alpha;

	int num_lines = evenstreamlines->trajs[traj]->nlinesegs;
	
	curpt[0] = curpt[1] = -1;

	if(cur_length >= interval)
	{
		alpha = (cur_length-interval)/evenstreamlines->trajs[traj]->linesegs[cur_lineindex].length;
		curpt[0] = alpha*evenstreamlines->trajs[traj]->linesegs[cur_lineindex].gstart[0] 
			+ (1-alpha)*evenstreamlines->trajs[traj]->linesegs[cur_lineindex].gend[0];
		curpt[1] = alpha*evenstreamlines->trajs[traj]->linesegs[cur_lineindex].gstart[1] 
			+ (1-alpha)*evenstreamlines->trajs[traj]->linesegs[cur_lineindex].gend[1];

			//if(alpha>1 || alpha<0)
			//{
			//	int test =0;
			//}

		cur_length -= interval;
		return true;
	}

	else{
		cur_lineindex++;
		if(cur_lineindex >= num_lines)
		{
			return false;
		}
		cur_length += evenstreamlines->trajs[traj]->linesegs[cur_lineindex].length;
		return false;
	}
}


SampleListInTriangle * EvenStreamlinePlace::extend_cell_samplelist(SampleListInTriangle *samplepts, int nsamples)
{
	if(samplepts == NULL || nsamples == 0)
	{
		samplepts = (SampleListInTriangle*)malloc(sizeof(SampleListInTriangle));
		return samplepts;
	}

	SampleListInTriangle *temp = samplepts;
	samplepts = (SampleListInTriangle*)malloc(sizeof(SampleListInTriangle)*(nsamples+1));
	for(int i=0; i<nsamples; i++)
		samplepts[i] = temp[i];
	return samplepts;
}


void EvenStreamlinePlace::add_sample_to_cell(int triangle, int which_traj, 
											 int which_sample, bool fieldtype)
{
	if(triangle < 0 || triangle >= quadmesh->nfaces)
		return;

	QuadCell *face = quadmesh->quadcells[triangle];
	
	//face->samplepts = extend_cell_samplelist(face->samplepts, face->num_samplepts);
	//
	//face->samplepts[face->num_samplepts].which_traj = which_traj;
	//face->samplepts[face->num_samplepts].which_sample = which_sample;

	//face->num_samplepts++;

	if(!fieldtype) /*major*/
	{
		face->maj_samplepts = extend_cell_samplelist(face->maj_samplepts, 
			face->maj_nsamplepts);
		
		face->maj_samplepts[face->maj_nsamplepts].which_traj = which_traj;
		face->maj_samplepts[face->maj_nsamplepts].which_sample = which_sample;

		face->maj_nsamplepts++;
	}
	else /*minor*/
	{
		face->min_samplepts = extend_cell_samplelist(face->min_samplepts, 
			face->min_nsamplepts);
		
		face->min_samplepts[face->min_nsamplepts].which_traj = which_traj;
		face->min_samplepts[face->min_nsamplepts].which_sample = which_sample;

		face->min_nsamplepts++;
	}
}


void EvenStreamlinePlace::init_samplelist_in_cell(bool fieldtype)
{
	int i;
	QuadCell *face;
	//for(i=0; i<quadmesh->nfaces; i++)
	//{
	//	face = quadmesh->quadcells[i];
	//	if(face->samplepts != NULL)
	//		free(face->samplepts);
	//	face->samplepts = NULL;
	//	face->num_samplepts = 0;
	//	//face->sampleindex=-1;
	//}

	if(!fieldtype) /*major*/
	{
		for(i=0; i<quadmesh->nfaces; i++)
		{
			face = quadmesh->quadcells[i];
			if(face->maj_samplepts != NULL)
				free(face->maj_samplepts);
			face->maj_samplepts = NULL;
			face->maj_nsamplepts = 0;
			face->MAJMaxSampNum=0;
		}
	}
	else  /*minor*/
	{
		for(i=0; i<quadmesh->nfaces; i++)
		{
			face = quadmesh->quadcells[i];
			if(face->min_samplepts != NULL)
				free(face->min_samplepts);
			face->min_samplepts = NULL;
			face->min_nsamplepts = 0;
			face->MINMaxSampNum=0;
		}
	}
}

int *EvenStreamlinePlace::cal_euclidean_dist(int triangle, double p[2], double dsep, double discsize, 
					int *NearbyTriangles, int &num_triangles)
{
	////we suppose the global point always falls in the triangle
	QuadCell *face = quadmesh->quadcells[triangle];
	QuadVertex *vert;
	QuadEdge *cur_edge;
	icVector3 dis;

	////Find all the triangle
	num_triangles = 0;
	NearbyTriangles = Extend_link(NearbyTriangles, num_triangles);
	NearbyTriangles[0] = triangle;
	num_triangles++;
	int cur_id = 0;
	int i, j;

	while(cur_id < num_triangles)
	{
		face = quadmesh->quadcells[NearbyTriangles[cur_id]];

		for(i = 0; i < face->nverts; i++)
		{
			vert = quadmesh->quad_verts[face->verts[i]];

			if(vert->visited)
				continue;

			dis.entry[0] = vert->x - p[0];
			dis.entry[1] = vert->y - p[1];

			vert->distance = length(dis); /*compute the Euclidean distance*/

			vert->visited = true; //set the visited flag
			if(vert->distance > discsize*dsep)
				continue;

			////if the distance between the vertex and the input point is smaller than the threshold
			////we need to add all its adjacent triangles into the triangle list
			for(j = 0; j < vert->ncells; j++)
			{
				QuadCell *c = vert->cells[j];

				if(c->index < 0) //reach the boundary!
					continue;

				if(IsRepeated(NearbyTriangles, c->index, num_triangles))
					continue;

				NearbyTriangles = Extend_link(NearbyTriangles, num_triangles);
				NearbyTriangles[num_triangles] = c->index;
				num_triangles++;
			}

		}

		cur_id++;
	}
	return NearbyTriangles;
}



/*
Use new data structure to save the geodesic disc
*/
void EvenStreamlinePlace::cal_euclidean_dist_2(int triangle, double p[3], double dsep, double discsize, 
					DynList_Int *trianglelist)
{
	////we suppose the global point always falls in the triangle
	QuadCell *face = quadmesh->quadcells[triangle];
	QuadVertex *vert;
	QuadEdge *cur_edge;
	icVector3 dis;

	////Find all the triangle
	trianglelist->nelems = 0;
	trianglelist->add_New(triangle);
	int cur_id = 0;
	int i, j;

	while(cur_id < trianglelist->nelems)
	{
		face = quadmesh->quadcells[trianglelist->elems[cur_id]];

		for(i = 0; i < face->nverts; i++)
		{
			vert = quadmesh->quad_verts[face->verts[i]];

			if(vert->visited)
				continue;

			dis.entry[0] = vert->x - p[0];
			dis.entry[1] = vert->y - p[1];

			vert->distance = length(dis); /*compute the Euclidean distance*/

			vert->visited = true; //set the visited flag
			if(vert->distance > discsize*dsep)
				continue;

			////if the distance between the vertex and the input point is smaller than the threshold
			////we need to add all its adjacent triangles into the triangle list
			for(j = 0; j < vert->ncells; j++)
			{
				QuadCell *c = vert->cells[j];

				if(c->index < 0) //reach the boundary!
					continue;
				trianglelist->add_New(c->index);
			}
		}

		cur_id++;
	}

	for(i=0; i<trianglelist->nelems; i++)
	{
		face = quadmesh->quadcells[trianglelist->elems[i]];
		face->visited = false;
		for(j=0; j<face->nverts; j++)
		{
			vert = quadmesh->quad_verts[face->verts[j]];
			vert->visited = false;
		}
	}
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*****************************************************************************/
/*we compute the intersections according to the information stored in
each cell*/


void init_streetnet()
{
	if(streetnet == NULL)
		streetnet=new StreetNet(200, 800);
	streetnet->reset_streetnet();
	init_tensorline_intersectionlists();
}


/*this will be called after the tensor line placement*/
void init_tensorline_intersectionlists()
{
	int i;
	if(majorintersectinfo != NULL)
	{
		for(i=0;i<prev_nmajors; i++)
			if(majorintersectinfo[i]!=NULL)
			{
				delete majorintersectinfo[i];
				majorintersectinfo[i]=NULL;
			}
		delete [] majorintersectinfo;
		majorintersectinfo = NULL;
	}

	if(majorintersectinfo == NULL)
	{
		majorintersectinfo = new TensorLineIntersectionInfoList *[major->evenstreamlines->ntrajs];
		for(i=0; i<major->evenstreamlines->ntrajs; i++)
			majorintersectinfo[i] = new TensorLineIntersectionInfoList(false,i);
		prev_nmajors=major->evenstreamlines->ntrajs;
	}
	
	if(minorintersectinfo != NULL)
	{
		for(i=0;i<prev_nminors; i++)
			if(minorintersectinfo[i]!=NULL)
			{
				delete minorintersectinfo[i];
				minorintersectinfo[i] = NULL;
			}
		delete [] minorintersectinfo;
		minorintersectinfo = NULL;
	}
	
	if(minorintersectinfo == NULL)
	{
		minorintersectinfo = new TensorLineIntersectionInfoList *[minor->evenstreamlines->ntrajs];
		for(i=0; i<minor->evenstreamlines->ntrajs; i++)
			minorintersectinfo[i] = new TensorLineIntersectionInfoList(true,i);
		prev_nminors=minor->evenstreamlines->ntrajs;
	}
}


void init_streetnetwork_info_in_cells()
{
	int i;
	QuadCell *face;
	for(i=0;i<quadmesh->nfaces;i++)
	{
		face=quadmesh->quadcells[i];

		/*   initialize the intersection lists  */
		if(face->intersectlist!=NULL)
		{
			free(face->intersectlist);
			face->intersectlist=NULL;
		}
		face->nintersects=0;

		/*   initialize the edge lists   */
		if(face->streetgraphedgelist != NULL)
		{
			free(face->streetgraphedgelist);
			face->streetgraphedgelist=NULL;
		}
		face->nstreetgraphedges=0;
	}
}


void compute_intersects()
{
	/*   Initialize the intersection list in each cell   */
	int i;
	QuadCell *face;
	//for(i=0;i<quadmesh->nfaces;i++)
	//{
	//	face=quadmesh->quadcells[i];

	//	/*   initialize the intersection lists  */
	//	if(face->intersectlist!=NULL)
	//	{
	//		free(face->intersectlist);
	//		face->intersectlist=NULL;
	//	}
	//	face->nintersects=0;

	//	/*   initialize the edge lists   */
	//	if(face->streetgraphedgelist != NULL)
	//	{
	//		free(face->streetgraphedgelist);
	//		face->streetgraphedgelist=NULL;
	//	}
	//	face->nstreetgraphedges=0;
	//}
	init_streetnetwork_info_in_cells();

	/*we first locate the cells that may contain intersections*/
	for(i=0; i<quadmesh->nfaces; i++)
	{
		face = quadmesh->quadcells[i];

		if(face->hasmajor && face->hasminor)
			compute_intersects_in_cell(i);

		if(face->hasmajor && face->majorlines!=NULL && face->majorlines->nlines>=2)
			compute_intersects_in_cell_sametype(i, false);
		
		if(face->hasminor && face->minorlines!=NULL && face->minorlines->nlines>=2)
			compute_intersects_in_cell_sametype(i, true);
	}

	/*    second, 
	      record the end point of the tensor lines
	*/
	
	Trajectory *temp;
	icVector2 dist;
	Intersection *intersect;

	/*          major tensor lines          */

	for(i=0; i<major->evenstreamlines->ntrajs; i++)
	{
		/*record the start point*/
		temp = major->evenstreamlines->trajs[i];

		if(temp->closed)
			continue;

		/*if it is exactly (close enough to) the first intersection*/
		if(majorintersectinfo[i]->nelems>0)
		{
			intersect=streetnet->nodelist->intersects[majorintersectinfo[i]->infolist[0]->intersect_id];
			dist.entry[0]=temp->linesegs[0].gstart[0]-intersect->gpos[0];
			dist.entry[1]=temp->linesegs[0].gstart[1]-intersect->gpos[1];
			if(length(dist)<1e-3) continue;
		}

		Intersection *newintersect=(Intersection*)malloc(sizeof(Intersection));
		newintersect->gpos[0]=temp->linesegs[0].gstart[0];
		newintersect->gpos[1]=temp->linesegs[0].gstart[1];
		newintersect->cellid=temp->linesegs[0].Triangle_ID;
		newintersect->majorline_id=temp->index;
		newintersect->minorline_id=-1;
		newintersect->majlineseg=0;
		newintersect->nadjedges=0;
		newintersect->adj_edges=NULL;
		newintersect->endpt=true;
		newintersect->inside_region=false;
		newintersect->deleted=false;

		//streetnet->danglepts->addNew(newintersect);
		//newintersect->index += streetnet->nodelist->nelems; /*change the index of it*/
		streetnet->nodelist->addNew(newintersect);

		//if(newintersect->cellid<0)
		//{
		//	int test=0;
		//}
		add_to_cell_intersectlist(newintersect->cellid, newintersect->index, false);

		/*compute the edges incident to it*/
		IntersectionInfo *majinfo1=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
		majinfo1->intersect_id=streetnet->nodelist->nelems-1;
		majinfo1->lineseg_id = 0;
		majorintersectinfo[i]->sorted_add(majinfo1);
		
		/*record the end point*/
				
		/*if it is exactly (close enough to) the first intersection*/
		if(majorintersectinfo[i]->nelems>0)
		{
			intersect=streetnet->nodelist->intersects[
				majorintersectinfo[i]->infolist[majorintersectinfo[i]->nelems-1]->intersect_id];
			dist.entry[0]=temp->linesegs[temp->nlinesegs-1].gend[0]-intersect->gpos[0];
			dist.entry[1]=temp->linesegs[temp->nlinesegs-1].gend[1]-intersect->gpos[1];
			if(length(dist)<1e-3) continue;
		}


		temp = major->evenstreamlines->trajs[i];
		Intersection *newintersect2=(Intersection*)malloc(sizeof(Intersection));
		newintersect2->gpos[0]=temp->linesegs[temp->nlinesegs-1].gend[0];
		newintersect2->gpos[1]=temp->linesegs[temp->nlinesegs-1].gend[1];
		newintersect2->cellid=temp->linesegs[temp->nlinesegs-1].Triangle_ID;
		newintersect2->majorline_id=temp->index;
		newintersect2->minorline_id=-1;
		newintersect2->majlineseg=temp->nlinesegs-1;
		newintersect2->nadjedges=0;
		newintersect2->adj_edges=NULL;
		newintersect2->endpt=true;
		newintersect2->inside_region=false;
		newintersect2->deleted=false;

		//streetnet->danglepts->addNew(newintersect2);
		//newintersect2->index += streetnet->nodelist->nelems; /*change the index of it*/
		streetnet->nodelist->addNew(newintersect2);

		//if(newintersect2->cellid<0)
		//{
		//	int test=0;
		//}
		add_to_cell_intersectlist(newintersect2->cellid, newintersect2->index, false);

		/*compute the edges incident to it*/
		IntersectionInfo *majinfo2=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
		majinfo2->intersect_id=streetnet->nodelist->nelems-1;
		majinfo2->lineseg_id = temp->nlinesegs-1;
		majorintersectinfo[i]->sorted_add(majinfo2);
	}

	/*minor tensor lines*/
	for(i=0; i<minor->evenstreamlines->ntrajs; i++)
	{
		/*record the start point*/
		temp = minor->evenstreamlines->trajs[i];
		if(temp->closed)
			continue;

		/*if it is exactly (close enough to) the first intersection*/
		if(minorintersectinfo[i]->nelems>0)
		{
			intersect=streetnet->nodelist->intersects[minorintersectinfo[i]->infolist[0]->intersect_id];
			dist.entry[0]=temp->linesegs[0].gstart[0]-intersect->gpos[0];
			dist.entry[1]=temp->linesegs[0].gstart[1]-intersect->gpos[1];
			if(length(dist)<1e-3) continue;
		}

		Intersection *newintersect=(Intersection*)malloc(sizeof(Intersection));
		newintersect->gpos[0]=temp->linesegs[0].gstart[0];
		newintersect->gpos[1]=temp->linesegs[0].gstart[1];
		newintersect->cellid=temp->linesegs[0].Triangle_ID;
		newintersect->minorline_id=temp->index;
		newintersect->majorline_id=-1;
		newintersect->minlineseg=0;
		newintersect->nadjedges=0;
		newintersect->adj_edges=NULL;
		newintersect->endpt=true;
		newintersect->inside_region=false;
		newintersect->deleted=false;
		newintersect->intersect_type=0;


		//streetnet->danglepts->addNew(newintersect);
		//newintersect2->index += streetnet->nodelist->nelems; /*change the index of it*/
		streetnet->nodelist->addNew(newintersect);
		
		//if(newintersect->cellid<0)
		//{
		//	int test=0;
		//}
		add_to_cell_intersectlist(newintersect->cellid, newintersect->index, false);

		/*compute the edges incident to it*/
		IntersectionInfo *mininfo1=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
		mininfo1->intersect_id=streetnet->nodelist->nelems-1;
		mininfo1->lineseg_id = 0;
		minorintersectinfo[i]->sorted_add(mininfo1);
		
		/*record the end point*/
		/*if it is exactly (close enough to) the first intersection*/

		if(minorintersectinfo[i]->nelems>0)
		{
			intersect=streetnet->nodelist->intersects[
				minorintersectinfo[i]->infolist[minorintersectinfo[i]->nelems-1]->intersect_id];
			dist.entry[0]=temp->linesegs[temp->nlinesegs-1].gend[0]-intersect->gpos[0];
			dist.entry[1]=temp->linesegs[temp->nlinesegs-1].gend[1]-intersect->gpos[1];
			if(length(dist)<1e-3) continue;
		}

		temp = minor->evenstreamlines->trajs[i];
		Intersection *newintersect2=(Intersection*)malloc(sizeof(Intersection));
		newintersect2->gpos[0]=temp->linesegs[temp->nlinesegs-1].gend[0];
		newintersect2->gpos[1]=temp->linesegs[temp->nlinesegs-1].gend[1];
		newintersect2->cellid=temp->linesegs[temp->nlinesegs-1].Triangle_ID;
		newintersect2->minorline_id=temp->index;
		newintersect2->majorline_id=-1;
		newintersect2->minlineseg=temp->nlinesegs-1;
		newintersect2->nadjedges=0;
		newintersect2->adj_edges=NULL;
		newintersect2->endpt=true;
		newintersect2->inside_region=false;
		newintersect2->deleted=false;
		newintersect2->intersect_type=0;

		//streetnet->danglepts->addNew(newintersect2);
		//newintersect2->index += streetnet->nodelist->nelems; /*change the index of it*/
		streetnet->nodelist->addNew(newintersect2);
		
		//if(newintersect2->cellid<0)
		//{
		//	int test=0;
		//}

		add_to_cell_intersectlist(newintersect2->cellid, newintersect2->index, false);

		/*compute the edges incident to it*/
		IntersectionInfo *mininfo2=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
		mininfo2->intersect_id=streetnet->nodelist->nelems-1;
		mininfo2->lineseg_id = temp->nlinesegs-1;
		minorintersectinfo[i]->sorted_add(mininfo2);
	}
}

extern int GetIntersection2(double PointA[2], double PointB[2], 
					 double PointC[2], double PointD[2], double t[2]);


/*we now use brute force method to find out the intersections*/
bool compute_intersect_between_twolines(int majtraj, int majstart, int majend,
										int mintraj, int minstart, int minend,
										double intersect[2], 
										int &majlinesegid, int &minlinesegid)
{
	int i, j;
	Trajectory *major1, *minor1;
	double A[2], B[2], C[2], D[2], t[2];
	major1=major->evenstreamlines->trajs[majtraj];
	minor1=minor->evenstreamlines->trajs[mintraj];
	for(i=majstart; i<=majend; i++)
	{
		A[0]=major1->linesegs[i].gstart[0];
		A[1]=major1->linesegs[i].gstart[1];
		B[0]=major1->linesegs[i].gend[0];
		B[1]=major1->linesegs[i].gend[1];

		for(j=minstart; j<=minend; j++)
		{
			C[0]=minor1->linesegs[j].gstart[0];
			C[1]=minor1->linesegs[j].gstart[1];
			D[0]=minor1->linesegs[j].gend[0];
			D[1]=minor1->linesegs[j].gend[1];

			//if(GetIntersection2(A, B, C, D, t)==1)
			if(cal_intersect(A, B, C, D, t)==1)
			{
				intersect[0]=A[0]+t[0]*(B[0]-A[0]);
				intersect[1]=A[1]+t[0]*(B[1]-A[1]);
				majlinesegid = i;
				minlinesegid = j;
				return true;
			}
		}
	}
	return false;
}

/*we compute the intersections in a particular cell*/
void compute_intersects_in_cell(int cellid)
{
	QuadCell *face = quadmesh->quadcells[cellid];

	/*method 1: we use brute force method 10/02/2007*/
	int i, j;
	int majlinesegid, minlinesegid;
	double intersect[2] = {0.};

	/*Test file*/
	FILE *fp;

		for(i=0; i<face->majorlines->nlines; i++)
		{
			for(j=0; j<face->minorlines->nlines; j++)
			{
				if(compute_intersect_between_twolines(face->majorlines->lines[i]->whichtraj,
					face->majorlines->lines[i]->start, face->majorlines->lines[i]->end,
					face->minorlines->lines[j]->whichtraj,
					face->minorlines->lines[j]->start, face->minorlines->lines[j]->end,
					intersect, majlinesegid, minlinesegid))
				{
					/*create a new intersection data*/
					Intersection *newintersect=(Intersection*)malloc(sizeof(Intersection));
					newintersect->gpos[0]=intersect[0];
					newintersect->gpos[1]=intersect[1];
					cellid=get_cellID_givencoords(intersect[0], intersect[1]);
					newintersect->cellid=cellid;
					newintersect->majorline_id=face->majorlines->lines[i]->whichtraj;
					newintersect->minorline_id=face->minorlines->lines[j]->whichtraj;
					newintersect->majlineseg=majlinesegid;
					newintersect->minlineseg=minlinesegid;
					newintersect->nadjedges=0;
					newintersect->adj_edges=NULL;
					newintersect->endpt=false;
					newintersect->inside_region=false;
					newintersect->deleted=false;
					newintersect->intersect_type=0;

					streetnet->nodelist->addNew(newintersect);

					/*   We need to save the intersection into the list of the intersections
					     in the triangle that contains this intersection.
					*/
					add_to_cell_intersectlist(cellid, newintersect->index, false);
					
					IntersectionInfo *majinfo=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
					majinfo->intersect_id=streetnet->nodelist->nelems-1;
					majinfo->lineseg_id = majlinesegid;
					majorintersectinfo[newintersect->majorline_id]->sorted_add(majinfo);

					IntersectionInfo *mininfo=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
					mininfo->intersect_id=streetnet->nodelist->nelems-1;
					mininfo->lineseg_id = minlinesegid;
					minorintersectinfo[newintersect->minorline_id]->sorted_add(mininfo);
				}
			}
		}
}

/**
     We now consider the cases that two major (minor) roads intersect 
**/

bool compute_intersect_between_twolines_sametype(int majtraj, int majstart, int majend,
										int mintraj, int minstart, int minend,
										double intersect[2], 
										int &majlinesegid, int &minlinesegid, bool majormin)
{
	int i, j;
	Trajectory *major1, *minor1;
	double A[2], B[2], C[2], D[2], t[2];
	if(!majormin)  /*  two major roads intersect  */
	{
		major1=major->evenstreamlines->trajs[majtraj];
		minor1=major->evenstreamlines->trajs[mintraj];
	}
	else           /*  two minor roads intersect  */
	{
		major1=minor->evenstreamlines->trajs[majtraj];
		minor1=minor->evenstreamlines->trajs[mintraj];
	}

	//for(i=majstart; i<=majend; i++)
	for(i=max(0, majstart-1); i<=min(major1->nlinesegs-1, majend+1); i++)
	{
		A[0]=major1->linesegs[i].gstart[0];
		A[1]=major1->linesegs[i].gstart[1];
		B[0]=major1->linesegs[i].gend[0];
		B[1]=major1->linesegs[i].gend[1];

		//for(j=minstart; j<=minend; j++)
		for(j=max(0, minstart-1); j<=min(minor1->nlinesegs-1, minend+1); j++)
		{
			C[0]=minor1->linesegs[j].gstart[0];
			C[1]=minor1->linesegs[j].gstart[1];
			D[0]=minor1->linesegs[j].gend[0];
			D[1]=minor1->linesegs[j].gend[1];

			if(cal_intersect(A, B, C, D, t)==1)
			{
				intersect[0]=A[0]+t[0]*(B[0]-A[0]);
				intersect[1]=A[1]+t[0]*(B[1]-A[1]);
				majlinesegid = i;
				minlinesegid = j;
				return true;
			}
		}
	}
	return false;
}


/*
   compute the intersections between two major (or minor) roads
*/
void compute_intersects_in_cell_sametype(int cellid, bool majormin)
{
	QuadCell *face = quadmesh->quadcells[cellid];

	/*method 1: we use brute force method 10/02/2007*/
	int i, j;
	int majlinesegid, minlinesegid;
	double intersect[2] = {0.};

	/*Test file*/
	FILE *fp;

	LineInfo *thelineinfo;

	if(!majormin) /*  two major roads intersect  */
	{
		thelineinfo=face->majorlines;
	}
	else          /*  two minor roads intersect  */
	{
		thelineinfo=face->minorlines;
	}

	for(i=0; i<thelineinfo->nlines; i++)
	{
		for(j=i+1; j<thelineinfo->nlines; j++)
		{
			if(compute_intersect_between_twolines_sametype(thelineinfo->lines[i]->whichtraj,
				thelineinfo->lines[i]->start, thelineinfo->lines[i]->end,
				thelineinfo->lines[j]->whichtraj,
				thelineinfo->lines[j]->start, thelineinfo->lines[j]->end,
				intersect, majlinesegid, minlinesegid, majormin))
			{
				/*create a new intersection data*/
				Intersection *newintersect=(Intersection*)malloc(sizeof(Intersection));
				newintersect->gpos[0]=intersect[0];
				newintersect->gpos[1]=intersect[1];
				cellid=get_cellID_givencoords(intersect[0], intersect[1]);
				newintersect->cellid=cellid;
				newintersect->majorline_id=thelineinfo->lines[i]->whichtraj;
				newintersect->minorline_id=thelineinfo->lines[j]->whichtraj;
				newintersect->majlineseg=majlinesegid;
				newintersect->minlineseg=minlinesegid;
				newintersect->nadjedges=0;
				newintersect->adj_edges=NULL;
				newintersect->endpt=false;
				newintersect->inside_region=false;
				newintersect->deleted=false;

				if(!majormin) /* two major roads intersect */
					newintersect->intersect_type=1;
				else          /* two minor roads intersect */
					newintersect->intersect_type=2;

				streetnet->nodelist->addNew(newintersect);

				/*   We need to save the intersection into the list of the intersections
					    in the triangle that contains this intersection.
				*/
				add_to_cell_intersectlist(cellid, newintersect->index, false);
				
				if(!majormin)
				{
					IntersectionInfo *majinfo=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
					majinfo->intersect_id=streetnet->nodelist->nelems-1;
					majinfo->lineseg_id = majlinesegid;
					majorintersectinfo[newintersect->majorline_id]->sorted_add(majinfo);

					if(major->evenstreamlines->trajs[newintersect->majorline_id]->nlinesegs-1
						<majinfo->lineseg_id)
					{
						int test=0;
					}

					IntersectionInfo *majinfo2=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
					majinfo2->intersect_id=streetnet->nodelist->nelems-1;
					majinfo2->lineseg_id = minlinesegid;
					majorintersectinfo[newintersect->minorline_id]->sorted_add(majinfo2);

					if(major->evenstreamlines->trajs[newintersect->minorline_id]->nlinesegs-1
						<majinfo2->lineseg_id)
					{
						int test=0;
					}
				}
				else
				{
					IntersectionInfo *mininfo=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
					mininfo->intersect_id=streetnet->nodelist->nelems-1;
					mininfo->lineseg_id = minlinesegid;
					minorintersectinfo[newintersect->minorline_id]->sorted_add(mininfo);
					
					if(minor->evenstreamlines->trajs[newintersect->minorline_id]->nlinesegs-1
						<mininfo->lineseg_id)
					{
						int test=0;
					}
					
					IntersectionInfo *mininfo2=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
					mininfo2->intersect_id=streetnet->nodelist->nelems-1;
					mininfo2->lineseg_id = majlinesegid;
					minorintersectinfo[newintersect->majorline_id]->sorted_add(mininfo2);
					
					if(minor->evenstreamlines->trajs[newintersect->majorline_id]->nlinesegs-1
						<mininfo2->lineseg_id)
					{
						int test=0;
					}
				}
			}
		}
	}
}



/*
     !!!!NOTE: we consider only those new computed tensor lines
	 We need to avoid the previous tensor lines
*/
//void compute_intersect_in_boundarycell(int cell_id)
//{
//	QuadCell *face = quadmesh->quadcells[cell_id];
//
//	/*method 1: we use brute force method 10/02/2007*/
//	int i, j;
//	int majlinesegid, minlinesegid;
//	double intersect[2] = {0.};
//
//	/*Test file*/
//	FILE *fp;
//
//		for(i=0; i<face->majorlines->nlines; i++)
//		{
//			for(j=0; j<face->minorlines->nlines; j++)
//			{
//				if(compute_intersect_between_twolines(face->majorlines->lines[i]->whichtraj,
//					face->majorlines->lines[i]->start, face->majorlines->lines[i]->end,
//					face->minorlines->lines[j]->whichtraj,
//					face->minorlines->lines[j]->start, face->minorlines->lines[j]->end,
//					intersect, majlinesegid, minlinesegid))
//				{
//					/*create a new intersection data*/
//					Intersection *newintersect=(Intersection*)malloc(sizeof(Intersection));
//					newintersect->gpos[0]=intersect[0];
//					newintersect->gpos[1]=intersect[1];
//					cellid=get_cellID_givencoords(intersect[0], intersect[1]);
//					newintersect->cellid=cellid;
//					newintersect->majorline_id=face->majorlines->lines[i]->whichtraj;
//					newintersect->minorline_id=face->minorlines->lines[j]->whichtraj;
//					newintersect->majlineseg=majlinesegid;
//					newintersect->minlineseg=minlinesegid;
//					newintersect->nadjedges=0;
//					newintersect->adj_edges=NULL;
//					newintersect->endpt=false;
//					newintersect->inside_region=false;
//					newintersect->deleted=false;
//
//					streetnet->nodelist->addNew(newintersect);
//
//					/*   We need to save the intersection into the list of the intersections
//					     in the triangle that contains this intersection.
//					*/
//					add_to_cell_intersectlist(cellid, newintersect->index, false);
//					
//					IntersectionInfo *majinfo=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
//					majinfo->intersect_id=streetnet->nodelist->nelems-1;
//					majinfo->lineseg_id = majlinesegid;
//					majorintersectinfo[newintersect->majorline_id]->sorted_add(majinfo);
//
//					IntersectionInfo *mininfo=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
//					mininfo->intersect_id=streetnet->nodelist->nelems-1;
//					mininfo->lineseg_id = minlinesegid;
//					minorintersectinfo[newintersect->minorline_id]->sorted_add(mininfo);
//				}
//			}
//		}
//}




/*add the index of an intersection into the corresponding cell*/
void add_to_cell_intersectlist(int cellid, int intersectid, bool endpt)
{
	QuadCell *face = quadmesh->quadcells[cellid];

	if(face->intersectlist==NULL||face->nintersects==0)
		face->intersectlist=(int*)malloc(sizeof(int));
	else
		face->intersectlist = extend_link(face->intersectlist, face->nintersects);
	face->intersectlist[face->nintersects]=intersectid;
	face->nintersects++;
}


void add_edge_to_intersectnode(int node, int edgeindex)
{
	if(node < 0)
		return;

	streetnet->nodelist->intersects[node]->adj_edges = 
		extend_link(streetnet->nodelist->intersects[node]->adj_edges,
		streetnet->nodelist->intersects[node]->nadjedges);

	streetnet->nodelist->intersects[node]->adj_edges[
		streetnet->nodelist->intersects[node]->nadjedges] = edgeindex;
	streetnet->nodelist->intersects[node]->nadjedges++;
}



/*search for the connectivity information using the intersection lists for 
the major and minor tensor lines, respectively*/
void search_for_connection()
{
	/*we can simply search the two groups of the intersection lists once
	and construct the edges*/
	int i, j;
	TensorLineIntersectionInfoList *infolist;
	IntersectionInfo *info, *infonext;

	int start_lineseg, end_lineseg;
	double startp[2], endp[2];
	Trajectory *traj;

	/*search major lines first*/
	for(i=0; i<major->evenstreamlines->ntrajs; i++)
	{
		infolist=majorintersectinfo[i];
		traj=major->evenstreamlines->trajs[i];

		for(j=0; j<infolist->nelems-1; j++)
		{
			info = infolist->infolist[j];
			infonext = infolist->infolist[j+1];

			/*build connection between it and its succeeding*/
			//Graph_Edge *newedge = (Graph_Edge*)malloc(sizeof(Graph_Edge));
			//newedge->node_index1 = info->intersect_id;
			//newedge->node_index2 = infonext->intersect_id;
			//newedge->cancel = newedge->visited = false;
			//streetnet->edgelist->append(newedge);

			/*  we are NOT considering isolated edges now 12/27/2007  */
			if(streetnet->nodelist->intersects[info->intersect_id]->endpt
				&&streetnet->nodelist->intersects[infonext->intersect_id]->endpt)
				continue;
			if(has_edge_between_minRoad(info->intersect_id, infonext->intersect_id))
				continue;

			/*we use new data structure for the street network 11/06/2007*/
			StreetGraphEdge *newedge = (StreetGraphEdge*)malloc(sizeof(StreetGraphEdge));
			newedge->node_index1 = info->intersect_id;
			newedge->node_index2 = infonext->intersect_id;
			newedge->cancel = newedge->visited = false;
			newedge->inter_pts = NULL; //we will add the intermediate points later
			newedge->ninter_pts = 0;
			newedge->roadtype=traj->roadtype;
			streetnet->edgelist->append(newedge);
			
			/* resample the line and record the information into the 
			corresponding line  11/26/2007 */
			start_lineseg=info->lineseg_id;
			end_lineseg=infonext->lineseg_id;

			startp[0]=streetnet->nodelist->intersects[info->intersect_id]->gpos[0];
			startp[1]=streetnet->nodelist->intersects[info->intersect_id]->gpos[1];
			endp[0]=streetnet->nodelist->intersects[infonext->intersect_id]->gpos[0];
			endp[1]=streetnet->nodelist->intersects[infonext->intersect_id]->gpos[1];

			sample_along_tensorline_from_to(traj, start_lineseg, startp,
									 end_lineseg, endp, 
									 newedge->index, streetnet);

			/**/
			add_edge_to_intersectnode(info->intersect_id, newedge->index);
			add_edge_to_intersectnode(infonext->intersect_id, newedge->index);
		}

		/*  if the tensor line is closed, add one more edge */
		if(major->evenstreamlines->trajs[i]->closed)
		{
			info = infolist->infolist[infolist->nelems-1];
			infonext = infolist->infolist[0];

			//if(streetnet->nodelist->intersects[info->intersect_id]->endpt
			//	&&streetnet->nodelist->intersects[infonext->intersect_id]->endpt)
			//	continue;
			if(!has_edge_between_minRoad(info->intersect_id, infonext->intersect_id))
			{
				/*we use new data structure for the street network 11/06/2007*/
				StreetGraphEdge *onenewedge = (StreetGraphEdge*)malloc(sizeof(StreetGraphEdge));
				onenewedge->node_index1 = info->intersect_id;
				onenewedge->node_index2 = infonext->intersect_id;
				onenewedge->cancel = onenewedge->visited = false;
				onenewedge->inter_pts = NULL; //we will add the intermediate points later
				onenewedge->ninter_pts = 0;
				onenewedge->roadtype=traj->roadtype;
				streetnet->edgelist->append(onenewedge);
				
				/* resample the line and record the information into the 
				corresponding line  11/26/2007 */
				start_lineseg=info->lineseg_id;
				end_lineseg=infonext->lineseg_id;

				startp[0]=streetnet->nodelist->intersects[info->intersect_id]->gpos[0];
				startp[1]=streetnet->nodelist->intersects[info->intersect_id]->gpos[1];
				int start_cell=get_cellID_givencoords(startp[0], startp[1]);
				endp[0]=streetnet->nodelist->intersects[infonext->intersect_id]->gpos[0];
				endp[1]=streetnet->nodelist->intersects[infonext->intersect_id]->gpos[1];
				int end_cell=get_cellID_givencoords(endp[0], endp[1]);

				Trajectory *tempTraj=new Trajectory(-1);
				get_linesegs_anytwopts(startp,start_cell, endp,end_cell, tempTraj, 0, 10);
				onenewedge->inter_pts=(Point **)malloc(sizeof(Point *)* (tempTraj->nlinesegs+1));

				int m;
				int pre_cell=-1;
				for(m=0;m<tempTraj->nlinesegs;m++)
				{
					onenewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
					onenewedge->inter_pts[m]->x=tempTraj->linesegs[m].gstart[0];
					onenewedge->inter_pts[m]->y=tempTraj->linesegs[m].gstart[1];
					onenewedge->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge->inter_pts[m]->x,
						onenewedge->inter_pts[m]->y);
					if(onenewedge->inter_pts[m]->cellid!=pre_cell)
					{
						add_to_edgelist_one_cell(onenewedge->inter_pts[m]->cellid, i);
						pre_cell=onenewedge->inter_pts[m]->cellid;
					}
				}

				onenewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
				onenewedge->inter_pts[m]->x=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[0];
				onenewedge->inter_pts[m]->y=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[1];
				onenewedge->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge->inter_pts[m]->x,
					onenewedge->inter_pts[m]->y);
				onenewedge->ninter_pts=tempTraj->nlinesegs+1;
				delete tempTraj;

				/**/
				add_edge_to_intersectnode(info->intersect_id, onenewedge->index);
				add_edge_to_intersectnode(infonext->intersect_id, onenewedge->index);
			}
		}
	}

	/*search minor lines first*/
	for(i=0; i<minor->evenstreamlines->ntrajs; i++)
	{
		infolist=minorintersectinfo[i];
		traj=minor->evenstreamlines->trajs[i];

		for(j=0; j<infolist->nelems-1; j++)
		{
			info = infolist->infolist[j];
			infonext = infolist->infolist[j+1];

			if(has_edge_between_minRoad(info->intersect_id, infonext->intersect_id))
				continue;
		
			/*  we are NOT considering isolated edges now 12/27/2007  */
			if(streetnet->nodelist->intersects[info->intersect_id]->endpt
				&&streetnet->nodelist->intersects[infonext->intersect_id]->endpt)
				continue;

			/*we use new data structure for the street network 11/06/2007*/
			StreetGraphEdge *newedge2 = (StreetGraphEdge*)malloc(sizeof(StreetGraphEdge));
			newedge2->node_index1 = info->intersect_id;
			newedge2->node_index2 = infonext->intersect_id;
			newedge2->cancel = newedge2->visited = false;
			newedge2->inter_pts = NULL; //we will add the intermediate points later
			newedge2->ninter_pts = 0;
			newedge2->roadtype=traj->roadtype;
			streetnet->edgelist->append(newedge2);

			/* resample the line and record the information into the 
			corresponding line  11/26/2007 */
			start_lineseg=info->lineseg_id;
			end_lineseg=infonext->lineseg_id;

			startp[0]=streetnet->nodelist->intersects[info->intersect_id]->gpos[0];
			startp[1]=streetnet->nodelist->intersects[info->intersect_id]->gpos[1];
			endp[0]=streetnet->nodelist->intersects[infonext->intersect_id]->gpos[0];
			endp[1]=streetnet->nodelist->intersects[infonext->intersect_id]->gpos[1];

			sample_along_tensorline_from_to(traj, start_lineseg, startp,
									 end_lineseg, endp, 
									 newedge2->index, streetnet);

			/**/
			add_edge_to_intersectnode(info->intersect_id, newedge2->index);
			add_edge_to_intersectnode(infonext->intersect_id, newedge2->index);
		}

		/*  if the tensor line is closed, add one more edge */
		if(minor->evenstreamlines->trajs[i]->closed)
		{
			info = infolist->infolist[infolist->nelems-1];
			infonext = infolist->infolist[0];

			//if(streetnet->nodelist->intersects[info->intersect_id]->endpt
			//	&&streetnet->nodelist->intersects[infonext->intersect_id]->endpt)
			//	continue;
			if(!has_edge_between_minRoad(info->intersect_id, infonext->intersect_id))
			{
				/*we use new data structure for the street network 11/06/2007*/
				StreetGraphEdge *onenewedge = (StreetGraphEdge*)malloc(sizeof(StreetGraphEdge));
				onenewedge->node_index1 = info->intersect_id;
				onenewedge->node_index2 = infonext->intersect_id;
				onenewedge->cancel = onenewedge->visited = false;
				onenewedge->inter_pts = NULL; //we will add the intermediate points later
				onenewedge->ninter_pts = 0;
				onenewedge->roadtype=traj->roadtype;
				streetnet->edgelist->append(onenewedge);
				
				/* resample the line and record the information into the 
				corresponding line  11/26/2007 */
				start_lineseg=info->lineseg_id;
				end_lineseg=infonext->lineseg_id;

				startp[0]=streetnet->nodelist->intersects[info->intersect_id]->gpos[0];
				startp[1]=streetnet->nodelist->intersects[info->intersect_id]->gpos[1];
				int start_cell=get_cellID_givencoords(startp[0], startp[1]);
				endp[0]=streetnet->nodelist->intersects[infonext->intersect_id]->gpos[0];
				endp[1]=streetnet->nodelist->intersects[infonext->intersect_id]->gpos[1];
				int end_cell=get_cellID_givencoords(endp[0], endp[1]);

				Trajectory *tempTraj=new Trajectory(-1);
				get_linesegs_anytwopts(startp,start_cell, endp,end_cell, tempTraj, 0, 10);
				onenewedge->inter_pts=(Point **)malloc(sizeof(Point *)* (tempTraj->nlinesegs+1));

				int m;
				int pre_cell=-1;
				for(m=0;m<tempTraj->nlinesegs;m++)
				{
					onenewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
					onenewedge->inter_pts[m]->x=tempTraj->linesegs[m].gstart[0];
					onenewedge->inter_pts[m]->y=tempTraj->linesegs[m].gstart[1];
					onenewedge->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge->inter_pts[m]->x,
						onenewedge->inter_pts[m]->y);
					if(onenewedge->inter_pts[m]->cellid!=pre_cell)
					{
						add_to_edgelist_one_cell(onenewedge->inter_pts[m]->cellid, i);
						pre_cell=onenewedge->inter_pts[m]->cellid;
					}
				}

				onenewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
				onenewedge->inter_pts[m]->x=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[0];
				onenewedge->inter_pts[m]->y=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[1];
				onenewedge->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge->inter_pts[m]->x,
					onenewedge->inter_pts[m]->y);
				onenewedge->ninter_pts=tempTraj->nlinesegs+1;
				delete tempTraj;

				/**/
				add_edge_to_intersectnode(info->intersect_id, onenewedge->index);
				add_edge_to_intersectnode(infonext->intersect_id, onenewedge->index);
			}
		}
	}


	/*  finally mark those nodes with only one edge as dead end  */
	for(i=0;i<streetnet->nodelist->nelems;i++)
	{
		if(streetnet->nodelist->intersects[i]->nadjedges<=1)
			streetnet->nodelist->intersects[i]->endpt=true;
	}
}


/*                                                                          */
/*   routines for saving/loading the obtained street network into a file    */

extern void update_street_network();
void save_cur_street_network(char *filename)
{
	if(streetnet==NULL || streetnet->nodelist->nelems==0
		||streetnet->edgelist->nedges==0)
		return;

	/*  remove isolated nodes first 
	   possible bug: 1/21/2008
	*/
	int i;
	//for(i=0;i<streetnet->nodelist->nelems;i++)
	//{
	//	if(streetnet->nodelist->intersects[i]->nadjedges==0)
	//		streetnet->nodelist->intersects[i]->deleted=true;
	//}
	
	update_street_network();

	FILE *fp=fopen(filename, "w");

	if(fp==NULL)
		return;

	int counter=0;
	Intersection *cur_intersect;

	/*    Save the street intersection information     */
	for(i=0;i<streetnet->nodelist->nelems;i++)
	{
		cur_intersect=streetnet->nodelist->intersects[i];
		if(cur_intersect->deleted)
			continue;
		counter++;
	}

	fprintf(fp, "#VERTICES: %d\n", counter/*streetnet->nodelist->nelems*/);
	for(i=0;i<streetnet->nodelist->nelems;i++)
	{
		cur_intersect=streetnet->nodelist->intersects[i];
		if(cur_intersect->deleted)
			continue;

		fprintf(fp, "%d %f %f\n", cur_intersect->nadjedges, 
			cur_intersect->gpos[0], cur_intersect->gpos[1]);
	}
	
	/*    Save the street line segment information     */
	counter=0;
	StreetGraphEdge *cur_edge;
	for(i=0;i<streetnet->edgelist->nedges;i++)
	{
		cur_edge=streetnet->edgelist->edges[i];
		if(cur_edge->cancel)
			continue;
		counter++;
	}
	fprintf(fp, "#EDGES: %d\n", counter /*streetnet->edgelist->nedges*/);
	for(i=0;i<streetnet->edgelist->nedges;i++)
	{
		cur_edge=streetnet->edgelist->edges[i];
		if(cur_edge->cancel)
			continue;

		if(cur_edge->node_index1>=streetnet->nodelist->nelems
			||cur_edge->node_index2>=streetnet->nodelist->nelems)
		{
			int test=0;
		}

		int streetwidth=2;
		if(cur_edge->roadtype==MINOR)
			streetwidth=2;
		else if(cur_edge->roadtype==MAJOR)
			streetwidth=4;
		else if(cur_edge->roadtype==HIGHWAY)
			streetwidth=6;
		fprintf(fp, "%d %d %d %d %c\n", cur_edge->node_index1, cur_edge->node_index2,
			0, streetwidth, 'S'); /*  we use default setting right now  11/30/2007*/
	}
	fclose(fp);
}


/*  we need this functionality to fullfil the upstream editing */
bool load_a_street_network(char *filename)
{
	FILE *fp=fopen(filename, "r");

	if(fp==NULL)
		return false;

	int i;
	int counter=0;
	int nadjedges;
	float posx, posy;
	int streetwidth;
	char streetType;
	int node1, node2;
	int notsure;
	Intersection *theintersect, *theintersect2;

	FILE *fpt;

	/*  initialize the street network variable  */
	init_streetnet();
	init_streetnetwork_info_in_cells();

	fscanf(fp, "#VERTICES: %d\n", &counter);


	//fpt=fopen("graphload_test.txt", "w");
	//fprintf(fpt, "start loading the vertices");
	//fclose(fpt);

	for(i=0;i<counter;i++)
	{
		fscanf(fp, "%d %f %f\n", &nadjedges, &posx, &posy);

		/*  create a new intersection */
		Intersection *newintersect=(Intersection *)malloc(sizeof(Intersection));
		newintersect->gpos[0]=posx;
		newintersect->gpos[1]=posy;
		newintersect->nadjedges=0; /* a bug */
		newintersect->adj_edges=(int*)malloc(sizeof(int)*(nadjedges));
		newintersect->index=i;
		newintersect->cellid=get_cellID_givencoords(posx, posy);
		newintersect->deleted=false;

		streetnet->nodelist->addNew(newintersect);
	}
	
	/*    Save the street line segment information     */
	counter=0;
	fscanf(fp, "#EDGES: %d\n", &counter);
	for(i=0;i<counter;i++)
	{
		fscanf(fp, "%d %d %d %d %c\n", &node1, &node2, &notsure,
			&streetwidth, &streetType); /*  we use default setting right now  11/30/2007*/

		/*  create a new edge */
		StreetGraphEdge *newedge=(StreetGraphEdge*)malloc(sizeof(StreetGraphEdge));
		newedge->index=i;
		newedge->node_index1=node1;
		newedge->node_index2=node2;
		newedge->cancel=false;
		newedge->ninter_pts=0;
		newedge->inter_pts=NULL;

		if(streetwidth==2)
			newedge->roadtype=MINOR;
		else if(streetwidth ==4)
			newedge->roadtype=MAJOR;
		else
			newedge->roadtype=HIGHWAY;

		newedge->nregionblocks=0;
		newedge->visited=false;

		streetnet->edgelist->append(newedge);
		
		/*  add the edge to the corresponding edge lists of the two intersections  */
		theintersect=streetnet->nodelist->intersects[node1];
		theintersect->adj_edges[theintersect->nadjedges]=i;
		theintersect->nadjedges++;
		
		theintersect2=streetnet->nodelist->intersects[node2];
		theintersect2->adj_edges[theintersect2->nadjedges]=i;
		theintersect2->nadjedges++;

		
		/*  perform sub-sampling along the edge */
		Trajectory *temptraj=new Trajectory(-1);
		int cell1=get_cellID_givencoords(theintersect->gpos[0], theintersect->gpos[1]);
		int cell2=get_cellID_givencoords(theintersect2->gpos[0], theintersect2->gpos[1]);
		get_linesegs_anytwopts(theintersect->gpos, cell1, theintersect2->gpos, cell2,
			temptraj, 0, 50);
		
		newedge->inter_pts=(Point **)malloc(sizeof(Point *)* (temptraj->nlinesegs+1));
		int j;
		for(j=0;j<temptraj->nlinesegs;j++)
		{
			newedge->inter_pts[j]=(Point*)malloc(sizeof(Point));
			newedge->inter_pts[j]->x=temptraj->linesegs[j].gstart[0];
			newedge->inter_pts[j]->y=temptraj->linesegs[j].gstart[1];
			newedge->inter_pts[j]->cellid=get_cellID_givencoords(newedge->inter_pts[j]->x,
				newedge->inter_pts[j]->y);
			
			add_to_edgelist_one_cell(newedge->inter_pts[j]->cellid, i);

		}

		newedge->inter_pts[j]=(Point*)malloc(sizeof(Point));
		newedge->inter_pts[j]->x=temptraj->linesegs[temptraj->nlinesegs-1].gend[0];
		newedge->inter_pts[j]->y=temptraj->linesegs[temptraj->nlinesegs-1].gend[1];
		newedge->inter_pts[j]->cellid=get_cellID_givencoords(newedge->inter_pts[j]->x,
			newedge->inter_pts[j]->y);
		newedge->ninter_pts=temptraj->nlinesegs+1;

		//fpt=fopen("graphload_test.txt", "w");
		//fprintf(fpt, "finish loading edge %d", i);
		//fclose(fpt);
		delete temptraj;

	}
	fclose(fp);
	
	fpt=fopen("graphload_test.txt", "w");
	fprintf(fpt, "finish loading the edges");
	fclose(fpt);

	/*  set the 'end point' flag for each intersection of the graph */
	for(i=0;i<streetnet->nodelist->nelems;i++)
	{
		theintersect=streetnet->nodelist->intersects[i];
		theintersect->endpt=false;
		if(theintersect->nadjedges<=1)
			theintersect->endpt=true;
		
		add_to_cell_intersectlist(theintersect->cellid, theintersect->index, theintersect->endpt);
	}

	return true;
}

/***********************************************************************************/
/*
   saving functionality
*/

/*
    save the major roads as tensor lines
*/

bool save_obtained_majRoads_tenLines(char *filename)
{
	if(major_level1==NULL || minor_level1==NULL)
		return false;

	FILE *fp=fopen(filename, "w");
	if(fp==NULL)
		return false;

	int i, j;
	Trajectory *curtraj;

	/*  major direction  */
	fprintf(fp, "#major lines: %d\n", major_level1->evenstreamlines->ntrajs);
	for(i=0;i<major_level1->evenstreamlines->ntrajs;i++)
	{
		curtraj=major_level1->evenstreamlines->trajs[i];
		fprintf(fp, "#line %d: %d\n", i, curtraj->nlinesegs);
		if(curtraj->closed)
			fprintf(fp, "#closed: %c\n", 'y');
		else
			fprintf(fp, "#closed: %c\n", 'n');
		for(j=0;j<curtraj->nlinesegs;j++)
		{
			fprintf(fp, "%f, %f\n", curtraj->linesegs[j].gstart[0], 
				curtraj->linesegs[j].gstart[1]);
		}
		fprintf(fp, "%f, %f\n", curtraj->linesegs[curtraj->nlinesegs-1].gend[0],
			curtraj->linesegs[curtraj->nlinesegs-1].gend[1]);
	}
	
	/*  minor direction  */
	fprintf(fp, "#minor lines: %d\n", minor_level1->evenstreamlines->ntrajs);
	for(i=0;i<minor_level1->evenstreamlines->ntrajs;i++)
	{
		curtraj=minor_level1->evenstreamlines->trajs[i];
		fprintf(fp, "#line %d: %d\n", i, curtraj->nlinesegs);
		if(curtraj->closed)
			fprintf(fp, "#closed: %c\n", 'y');
		else
			fprintf(fp, "#closed: %c\n", 'n');
		for(j=0;j<curtraj->nlinesegs;j++)
		{
			fprintf(fp, "%f, %f\n", curtraj->linesegs[j].gstart[0], 
				curtraj->linesegs[j].gstart[1]);
		}
		fprintf(fp, "%f, %f\n", curtraj->linesegs[curtraj->nlinesegs-1].gend[0],
			curtraj->linesegs[curtraj->nlinesegs-1].gend[1]);
	}
	fclose(fp);
	return true;
}

bool load_majRoads_tenLines(char *filename)
{
	FILE *fp=fopen(filename, "r");
	if(fp==NULL)
		return false;

	int nlines, nlinesegs, lineindex;
	float x, y;
	int i, j;
	Trajectory *curtraj;
	icVector2 dist;
		
	char ch;

	/*  load the major direction  */
	fscanf(fp, "#major lines: %d\n", &nlines);

	if(major_level1 != NULL)
		delete major_level1;
	major_level1=new EvenStreamlinePlace(false, nlines);

	if(major_level1 == NULL)
		return false;

	for(i=0; i<nlines;i++)
	{
		fscanf(fp, "#line %d: %d\n", &lineindex, &nlinesegs);
		major_level1->evenstreamlines->trajs[i]=new Trajectory(i, nlinesegs);
		if(major_level1->evenstreamlines->trajs[i]==NULL)
			return false;

		fscanf(fp, "#closed: %c\n", &ch);
		if(ch == 'y')
			major_level1->evenstreamlines->trajs[i]->closed=true;
		else
			major_level1->evenstreamlines->trajs[i]->closed=false;

		curtraj=major_level1->evenstreamlines->trajs[i];
		fscanf(fp, "%f, %f\n", &x, &y);
		curtraj->linesegs[0].gstart[0]=x;
		curtraj->linesegs[0].gstart[1]=y;
		curtraj->linesegs[0].Triangle_ID=get_cellID_givencoords(x, y);
		for(j=1;j<nlinesegs;j++)
		{
			fscanf(fp, "%f, %f\n", &x, &y);

			curtraj->linesegs[j].gstart[0]=
			curtraj->linesegs[j-1].gend[0]=x;

			curtraj->linesegs[j].gstart[1]=
			curtraj->linesegs[j-1].gend[1]=y;

			dist.entry[0]=curtraj->linesegs[j-1].gend[0]-curtraj->linesegs[j-1].gstart[0];
			dist.entry[1]=curtraj->linesegs[j-1].gend[1]-curtraj->linesegs[j-1].gstart[1];
			curtraj->linesegs[j-1].length=length(dist);

			curtraj->linesegs[j].Triangle_ID=get_cellID_givencoords(x, y);

			//if(curtraj->linesegs[j].Triangle_ID!=curtraj->linesegs[j-1].Triangle_ID)
		}
		fscanf(fp, "%f, %f\n", &x, &y);
		curtraj->linesegs[nlinesegs-1].gend[0]=x;
		curtraj->linesegs[nlinesegs-1].gend[1]=y;
		dist.entry[0]=curtraj->linesegs[nlinesegs-1].gend[0]-curtraj->linesegs[nlinesegs-1].gstart[0];
		dist.entry[1]=curtraj->linesegs[nlinesegs-1].gend[1]-curtraj->linesegs[nlinesegs-1].gstart[1];
		curtraj->linesegs[nlinesegs-1].length=length(dist);
		curtraj->linesegs[nlinesegs-1].Triangle_ID=get_cellID_givencoords(x, y);
		curtraj->nlinesegs=nlinesegs;
	}
	major_level1->evenstreamlines->ntrajs=nlines;

	/*  load the minor direction  */
	fscanf(fp, "#minor lines: %d\n", &nlines);

	if(minor_level1 != NULL)
		delete minor_level1;
	minor_level1=new EvenStreamlinePlace(true, nlines);

	if(minor_level1 == NULL)
		return false;

	for(i=0; i<nlines;i++)
	{
		fscanf(fp, "#line %d: %d\n", &lineindex, &nlinesegs);
		minor_level1->evenstreamlines->trajs[i]=new Trajectory(i, nlinesegs);
		if(minor_level1->evenstreamlines->trajs[i]==NULL)
			return false;

		fscanf(fp, "#closed: %c\n", &ch);
		if(ch == 'y')
			minor_level1->evenstreamlines->trajs[i]->closed=true;
		else
			minor_level1->evenstreamlines->trajs[i]->closed=false;

		curtraj=minor_level1->evenstreamlines->trajs[i];
		fscanf(fp, "%f, %f\n", &x, &y);
		curtraj->linesegs[0].gstart[0]=x;
		curtraj->linesegs[0].gstart[1]=y;
		curtraj->linesegs[0].Triangle_ID=get_cellID_givencoords(x, y);
		for(j=1;j<nlinesegs;j++)
		{
			fscanf(fp, "%f, %f\n", &x, &y);

			curtraj->linesegs[j].gstart[0]=
			curtraj->linesegs[j-1].gend[0]=x;

			curtraj->linesegs[j].gstart[1]=
			curtraj->linesegs[j-1].gend[1]=y;

			curtraj->linesegs[j].Triangle_ID=get_cellID_givencoords(x, y);

			//if(curtraj->linesegs[j].Triangle_ID!=curtraj->linesegs[j-1].Triangle_ID)
		}
		fscanf(fp, "%f, %f\n", &x, &y);
		curtraj->linesegs[nlinesegs-1].gend[0]=x;
		curtraj->linesegs[nlinesegs-1].gend[1]=y;
		curtraj->linesegs[nlinesegs-1].Triangle_ID=get_cellID_givencoords(x, y);
		curtraj->nlinesegs=nlinesegs;
	}
	minor_level1->evenstreamlines->ntrajs=nlines;

	return true;
}

/*
    save the major roads as a street network
*/

bool save_obtained_majRoads_network(char *filename)
{
	if(majRoadnet==NULL || majRoadnet->nodelist->nelems==0
		||majRoadnet->edgelist->nedges==0)
		return false;
	
	FILE *fp=fopen(filename, "w");

	if(fp==NULL)
		return false;

	int i;
	int counter=0;
	Intersection *cur_intersect;

	/*    Save the street intersection information     */
	for(i=0;i<majRoadnet->nodelist->nelems;i++)
	{
		cur_intersect=majRoadnet->nodelist->intersects[i];
		if(cur_intersect->deleted)
			continue;
		counter++;
	}

	fprintf(fp, "#VERTICES: %d\n", counter/*streetnet->nodelist->nelems*/);
	for(i=0;i<majRoadnet->nodelist->nelems;i++)
	{
		cur_intersect=majRoadnet->nodelist->intersects[i];
		if(cur_intersect->deleted)
			continue;

		fprintf(fp, "%d %f %f\n", cur_intersect->nadjedges, 
			cur_intersect->gpos[0], cur_intersect->gpos[1]);
	}
	
	/*    Save the street line segment information     */
	counter=0;
	StreetGraphEdge *cur_edge;
	for(i=0;i<majRoadnet->edgelist->nedges;i++)
	{
		cur_edge=majRoadnet->edgelist->edges[i];
		if(cur_edge->cancel)
			continue;
		counter++;
	}
	fprintf(fp, "#EDGES: %d\n", counter /*streetnet->edgelist->nedges*/);
	for(i=0;i<majRoadnet->edgelist->nedges;i++)
	{
		cur_edge=majRoadnet->edgelist->edges[i];
		if(cur_edge->cancel)
			continue;

		int streetwidth=4;
		fprintf(fp, "%d %d %d %d %c\n", cur_edge->node_index1, cur_edge->node_index2,
			0, streetwidth, 'S'); /*  we use default setting right now  11/30/2007*/
	}
	fclose(fp);
	return true;
}

bool load_majRoads_network(char *filename)
{
	FILE *fp=fopen(filename, "r");

	if(fp==NULL)
		return false;

	int i;
	int counter=0;
	int nadjedges;
	float posx, posy;
	int streetwidth;
	char streetType;
	int node1, node2;
	int notsure;
	Intersection *theintersect, *theintersect2;

	FILE *fpt;

	/*  initialize the street network variable  */
	init_majRoadnet();
	//init_streetnetwork_info_in_cells();
	init_majRoadinfo_incells();
	init_majRoad_intersectionlists();

	fscanf(fp, "#VERTICES: %d\n", &counter);


	for(i=0;i<counter;i++)
	{
		fscanf(fp, "%d %f %f\n", &nadjedges, &posx, &posy);

		/*  create a new intersection */
		Intersection *newintersect=(Intersection *)malloc(sizeof(Intersection));
		newintersect->gpos[0]=posx;
		newintersect->gpos[1]=posy;
		newintersect->nadjedges=0; /* a bug */
		newintersect->adj_edges=(int*)malloc(sizeof(int)*(nadjedges));
		newintersect->index=i;
		newintersect->cellid=get_cellID_givencoords(posx, posy);
		newintersect->deleted=false;

		majRoadnet->nodelist->addNew(newintersect);
	}
	
	/*    Save the street line segment information     */
	counter=0;
	fscanf(fp, "#EDGES: %d\n", &counter);
	for(i=0;i<counter;i++)
	{
		fscanf(fp, "%d %d %d %d %c\n", &node1, &node2, &notsure,
			&streetwidth, &streetType); /*  we use default setting right now  11/30/2007*/

		/*  create a new edge */
		StreetGraphEdge *newedge=(StreetGraphEdge*)malloc(sizeof(StreetGraphEdge));
		newedge->index=i;
		newedge->node_index1=node1;
		newedge->node_index2=node2;
		newedge->cancel=false;
		newedge->ninter_pts=0;
		newedge->inter_pts=NULL;

		if(streetwidth==2)
			newedge->roadtype=MINOR;
		else if(streetwidth ==4)
			newedge->roadtype=MAJOR;
		else
			newedge->roadtype=HIGHWAY;

		newedge->nregionblocks=0;
		newedge->visited=false;

		majRoadnet->edgelist->append(newedge);
		
		/*  add the edge to the corresponding edge lists of the two intersections  */
		theintersect=majRoadnet->nodelist->intersects[node1];
		theintersect->adj_edges[theintersect->nadjedges]=i;
		theintersect->nadjedges++;
		
		theintersect2=majRoadnet->nodelist->intersects[node2];
		theintersect2->adj_edges[theintersect2->nadjedges]=i;
		theintersect2->nadjedges++;

		
		/*  perform sub-sampling along the edge */
		Trajectory *temptraj=new Trajectory(-1);
		int cell1=get_cellID_givencoords(theintersect->gpos[0], theintersect->gpos[1]);
		int cell2=get_cellID_givencoords(theintersect2->gpos[0], theintersect2->gpos[1]);
		get_linesegs_anytwopts(theintersect->gpos, cell1, theintersect2->gpos, cell2,
			temptraj, 0, 50);
		
		newedge->inter_pts=(Point **)malloc(sizeof(Point *)* (temptraj->nlinesegs+1));
		int j;
		for(j=0;j<temptraj->nlinesegs;j++)
		{
			newedge->inter_pts[j]=(Point*)malloc(sizeof(Point));
			newedge->inter_pts[j]->x=temptraj->linesegs[j].gstart[0];
			newedge->inter_pts[j]->y=temptraj->linesegs[j].gstart[1];
			newedge->inter_pts[j]->cellid=get_cellID_givencoords(newedge->inter_pts[j]->x,
				newedge->inter_pts[j]->y);
			
			add_to_edgelist_one_cell(newedge->inter_pts[j]->cellid, i);

		}

		newedge->inter_pts[j]=(Point*)malloc(sizeof(Point));
		newedge->inter_pts[j]->x=temptraj->linesegs[temptraj->nlinesegs-1].gend[0];
		newedge->inter_pts[j]->y=temptraj->linesegs[temptraj->nlinesegs-1].gend[1];
		newedge->inter_pts[j]->cellid=get_cellID_givencoords(newedge->inter_pts[j]->x,
			newedge->inter_pts[j]->y);
		newedge->ninter_pts=temptraj->nlinesegs+1;

		delete temptraj;

	}
	fclose(fp);

	/*  set the 'end point' flag for each intersection of the graph */
	for(i=0;i<majRoadnet->nodelist->nelems;i++)
	{
		theintersect=majRoadnet->nodelist->intersects[i];
		theintersect->endpt=false;
		if(theintersect->nadjedges<=1)
			theintersect->endpt=true;
		
		add_to_cell_intersectlist(theintersect->cellid, theintersect->index, theintersect->endpt);
	}

	return true;
}


/*
Get the resampled points for each edge in the street network
11/06/2007
*/
double street_sample_interval;


void resample_edge_streetnetwork(int edgeindex)
{
	StreetGraphEdge *edge = streetnet->edgelist->edges[edgeindex];

	int node1 = edge->node_index1;
	int node2 = edge->node_index2;

	/*find the corresponding intersections*/
	Intersection *intersect1=streetnet->nodelist->intersects[node1];
	Intersection *intersect2=streetnet->nodelist->intersects[node2];

	/*obtain the information of the tensor lines that the two intersections belong to*/
	
	/*first consider major line*/
	int line_id=intersect1->majorline_id;
	int line_start=intersect1->majlineseg;
	int line_end=intersect2->majlineseg;
	if(intersect2->majlineseg<line_start)
	{
		line_start=intersect2->majlineseg;
		line_end=intersect1->majlineseg;
		Intersection *temp=intersect1;
		intersect1=intersect2;
		intersect2=temp;
	}

	if(line_id != intersect2->majorline_id)
	{
		/*this pair of intersections is on a minor line*/
		line_id=intersect1->minorline_id;

		//set a debug here
		if(line_id != intersect2->minorline_id)
		{
			int test = 0; //something is wrong here
		}

		line_start=intersect1->minlineseg;
		line_end=intersect2->minlineseg;
		if(intersect2->minlineseg<line_start)
		{
			line_start=intersect2->minlineseg;
			line_end=intersect1->minlineseg;
			Intersection *temp=intersect1;
			intersect1=intersect2;
			intersect2=temp;
		}
	}


	/*re-sample the line segments between the two intersections*/

	int cur_line = line_start;
	icVector2 line_dir;
	LineSeg *cur_lineseg=&major->evenstreamlines->trajs[line_id]->linesegs[cur_line];
	line_dir.entry[0]=cur_lineseg->gend[0]-intersect1->gpos[0];
	line_dir.entry[1]=cur_lineseg->gend[1]-intersect1->gpos[1];
	double cur_length=length(line_dir);

	int count = 0;

	while(count<50) /*avoid infinite loop*/
	{
		if(cur_length>=street_sample_interval)
		{
			/*generate a new sample point*/
			cur_lineseg=&major->evenstreamlines->trajs[line_id]->linesegs[cur_line];

			double alpha = (cur_length-street_sample_interval)/cur_lineseg->length;
			double x = alpha*cur_lineseg->gstart[0] + (1-alpha)*cur_lineseg->gend[0];
			double y = alpha*cur_lineseg->gstart[1] + (1-alpha)*cur_lineseg->gend[1];

			cur_length -= street_sample_interval;
			Point *newpt=(Point*)malloc(sizeof(Point));
			newpt->x=x;
			newpt->y=y;
			newpt->cellid=cur_lineseg->Triangle_ID;

			/*need to add to the edge sample point list*/

			if(edge->ninter_pts==0)
			{
				edge->inter_pts = (Point**)malloc(sizeof(Point*));
				edge->inter_pts[0]=newpt;
			}
			else
			{
				//StreetGraphEdge **temp_e=edge->inter_pts;
				edge->inter_pts=(Point**)realloc(edge->inter_pts, sizeof(Point*)*(edge->ninter_pts+1));
				edge->inter_pts[edge->ninter_pts]=newpt;
			}

			edge->ninter_pts++;
		}
		else
		{
			//move to next line segment
			if(cur_line == line_end)
				break;

			cur_line += 1;
			cur_lineseg=&major->evenstreamlines->trajs[line_id]->linesegs[cur_line];

			if(cur_line == line_end)
			{
				line_dir.entry[0]=cur_lineseg->gstart[0]-intersect2->gpos[0];
				line_dir.entry[1]=cur_lineseg->gstart[1]-intersect2->gpos[1];
				cur_length+=length(line_dir);
			}
			else
				cur_length+=cur_lineseg->length;

			//if(cur_line == line_end && cur_length < street_sample_interval)
			//	break;
		}

		count++;
	}

}

void resample_alledges_streetnetwork()
{
	int i;
	for(i=0; i<streetnet->edgelist->nedges; i++)
		resample_edge_streetnetwork(i);
}


void StreetNet::add_edge_to_node(int node, int edgeindex)
{
	if(node < 0)
		return;

	nodelist->intersects[node]->adj_edges = extend_link(nodelist->intersects[node]->adj_edges,
		nodelist->intersects[node]->nadjedges);

	nodelist->intersects[node]->adj_edges[nodelist->intersects[node]->nadjedges] = edgeindex;
	nodelist->intersects[node]->nadjedges++;
}













/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/
/*Fast Marching for quad mesh 09/30/2007*/
MinHeap::MinHeap(int initsize)
{
	/*allocate the Heap according to user specification*/
	if(initsize==0)
	{
		elems=NULL;
		curMaxNum=initsize=0;
		return;
	}

	curMaxNum = initsize;
	elems=(HeapElem**)malloc(sizeof(HeapElem*)*curMaxNum);
	nelems = 0;

	for(int i=0; i<curMaxNum; i++)
		elems[i]=(HeapElem*)malloc(sizeof(HeapElem));
}


MinHeap::~MinHeap()
{
	if(elems!=NULL)
	{
		int i;
		for(i=0;i<curMaxNum;i++)
		{
			if(elems[i]!=NULL)
			{
				free(elems[i]);
				elems[i]=NULL;
			}
		}
		free(elems);
	}
}


void MinHeap::reset()
{
	nelems = 0;
}

bool MinHeap::is_empty()
{
	if(nelems == 0)
		return true;
	return false;
}


int MinHeap::get_pos(int vert)
{
	int i;
	for(i=0; i<nelems; i++)
	{
		if(elems[i]->vertid == vert)
			return i;
	}
	return -1;
}

//Operations
/*return the root of the heap, and restore the heap*/
int MinHeap::FindSmallest()
{
	int smallest_vert = elems[0]->vertid;

	if(smallest_vert<0)
	{
		int test = 0;
	}

	/*remove the first element and restore the heap property*/
	DownHeap(0);

	nelems --;

	return smallest_vert; //return the corresponding vertex index of the original first element of the heap
}


/*move forward one step from the given position,
it is equivalent to remove one element at the postion "pos" */
void MinHeap::move_forward(int pos)
{
	int i;
	for(i=pos; i<nelems-1; i++)
	{
		elems[i]->vertid=elems[i+1]->vertid;
		elems[i]->T=elems[i+1]->T;
	}
}

/*restore the min-heap property after removing one element
The order is downward along the heap
*/
void MinHeap::DownHeap(int ielem)
{
	/**/
	if(ielem == nelems-1) /*if it is the last element*/
	{
		return;
	}
	/*if it is one of the leafs*/
	if(ielem>=nelems/2)
	{
		/*move one step forward for the rest of the elements*/
		move_forward(ielem);

		/*restore the heap property for all the moved elements: UpHeap*/
		int i;
		for(i=ielem; i<nelems-1; i++)
			UpHeap(i);
	}
	else
	{
		/*find the child that has smallest value, move it upward,
		repeat until reaching the leaf level*/
		int child = get_leftchild(ielem);
		while(ielem < nelems/2)
		{
			if(child+1>=nelems||elems[child+1]->vertid<0)
				break;

			if(elems[child]->T<elems[child+1]->T)
			{
				/*move the left child upward*/
				elems[ielem]->vertid=elems[child]->vertid;
				elems[ielem]->T = elems[child]->T;
				ielem=child;
			}
			else
			{
				/*move the right child upward*/
				elems[ielem]->vertid=elems[child+1]->vertid;
				elems[ielem]->T = elems[child+1]->T;
				ielem=child+1;
			}
			child = get_leftchild(ielem);

			if(child<0) break;
		}

		/*for the leaf level:*/
		/*move one step forward for the rest of the elements*/
		move_forward(ielem);


		/*restore the heap property for all the moved elements: UpHeap*/
		int i;
		for(i=ielem; i<nelems-1; i++)
			UpHeap(i);

	}
}


/*start from the given position and restore the min-heap property upward
The order is upward along the heap
*/
void MinHeap::UpHeap(int ielem)
{
	/*start from the position ielem*/
	if(ielem==0) return;
	int parent = get_parent(ielem);
	if(elems[parent]->T<elems[ielem]->T) return;

	while(elems[parent]->T>elems[ielem]->T)
	{
		/*we swap them*/
		int tempi = elems[ielem]->vertid;
		double tempv = elems[ielem]->T;
		elems[ielem]->vertid = elems[parent]->vertid;
		elems[ielem]->T = elems[parent]->T;
		elems[parent]->vertid = tempi;
		elems[parent]->T = tempv;

		/*keep upward*/
		ielem = parent;
		parent = get_parent(ielem);
		if(parent<0) return;  /*should not reach here!*/
	}
}


/*Insert a new element at the end of the array,
and use UpHeap to restore the heap property*/
bool MinHeap::Insert(int vert, double T)
{
	/*extend the array if necessary*/
	if(nelems>=curMaxNum)
		if(!extend(100))
			return false;

	//if(vert<0)
	//{
	//	int test = 0;
	//}

	/*add the new element*/
	elems[nelems]->vertid = vert;
	elems[nelems]->T = T;
	nelems++;

	/*UpHeap to maintain the heap property*/
	UpHeap(nelems-1);
	return true;
}


/*extend the space of the array*/
bool MinHeap::extend(int step)
{
	HeapElem **temp = elems;
	elems=(HeapElem**)malloc(sizeof(HeapElem*)*(curMaxNum+step));

	if(elems == NULL)
		return false;

	int i;
	for(i=0; i<curMaxNum; i++)
		elems[i]=temp[i];
	for(i=curMaxNum; i<curMaxNum+step; i++)
	{
		elems[i]=(HeapElem*)malloc(sizeof(HeapElem));
		if(elems[i]==NULL)
		{
			return false;
		}
	}

	curMaxNum += step;
	free(temp);
	return true;
}


/*get the parent index of the given element*/
int MinHeap::get_parent(int child)
{
	if(child <=0 ) return -1;
	return (int)((child-1)/2);
}

/*get the left child index of the given element*/
int MinHeap::get_leftchild(int parent)
{
	if((parent*2+1)>nelems) return -1;

	return parent*2+1;
}

/*get the right child index of the given element*/
int MinHeap::get_rightchild(int parent)
{
	if((parent*2+2)>nelems) return -1;
	return parent*2+2;
}




/*The following we implement the fast marching method on quad mesh*/

/*we need a data structure to store the curve (zero level set),
but we allow open curve now*/
typedef struct BoundCurve{
	LineSeg **lines; /*we use short line segments to represent the boundary curve*/
	int nlines;
	int curMaxNum;
}BoundCurve;


/*NOTE: we could have several disconnected boundary curves*/
//BoundCurve *curves;
class BoundCurveList
{
public:
	BoundCurve **curves;
	int ncurves;
	int curMaxNum;

	BoundCurveList(int initsize = 1);
};


/* we implement the brush interface here 10/10/2007*/
#include "LimitCycleCreator.h "
#include ".\glview.h"
#include "regionsmooth_quad.h"
#include "tensorvis.h"

extern int *boundarycells;
extern int nboundarycells;
extern DesignTriangleCycle myCycle;
extern ctr_point *control_pts;        // allocate our control point array
extern int resolution;    // how many points in our output array
extern int num_shapecontrol_pts;
extern int num_curvepts_output;
extern int MaxNumShapeControlPts;
extern int HermiteStep;
//extern void display_quadmesh(GLenum mode);
extern void display_quadmesh_select(GLenum mode);


BoundCurveList *boundlist = NULL;
int *boundvertlist = NULL;
int nboundverts;
MinHeap *narrowband = NULL;
double brush_width;


/*This routine initialize the distances on 
all vertices to be as large as possible*/
void init_dis_verts()
{
	int i;
	for(i=0; i<quadmesh->nverts; i++)
	{
		quadmesh->quad_verts[i]->distance = 1.e50;
		quadmesh->quad_verts[i]->type = 0;
		////quadmesh->quad_verts[i]->which_region = 0;  //11/17/2007
		quadmesh->quad_verts[i]->Jacobian.set(0.);
	}

	//if(boundvertlist!=NULL)
	//	free(boundvertlist);
	//boundvertlist=(int*)malloc(sizeof(int)*quadmesh->nverts);
	//nboundverts=0;

	init_boundvertlist();

	//if(narrowband == NULL)
	//	narrowband = new MinHeap(100);
	//narrowband->reset();
	reset_narrowband();
}

void init_verts_for_brush()
{
	int i;
	for(i=0; i<quadmesh->nverts; i++)
	{
		quadmesh->quad_verts[i]->distance = 1.e50;
		quadmesh->quad_verts[i]->type = 0;
		////quadmesh->quad_verts[i]->which_region = 0;  //11/17/2007
		//quadmesh->quad_verts[i]->Jacobian.set(0.);
	}

	//if(boundvertlist!=NULL)
	//	free(boundvertlist);
	//boundvertlist=(int*)malloc(sizeof(int)*quadmesh->nverts);
	//nboundverts=0;

	init_boundvertlist();

	//if(narrowband == NULL)
	//	narrowband = new MinHeap(100);
	//narrowband->reset();
	reset_narrowband();
}


/*
In this routine, we will reset the distance related variable for each vertex only
*/

void reset_vert_dis()
{
	int i;
	for(i=0; i<quadmesh->nverts; i++)
	{
		quadmesh->quad_verts[i]->distance = 1.e50;
		quadmesh->quad_verts[i]->type = 0;
	}

	init_boundvertlist();

	reset_narrowband();
}

void init_verts_all()
{
	int i;
	for(i=0; i<quadmesh->nverts; i++)
	{
		quadmesh->quad_verts[i]->which_region = 0;  //11/17/2007
		quadmesh->quad_verts[i]->inbrushregion=false;  //11/17/2007
		quadmesh->quad_verts[i]->phi=0;
	}

	init_dis_verts();
}


void init_boundvertlist()
{
	if(boundvertlist!=NULL)
		free(boundvertlist);
	boundvertlist=(int*)malloc(sizeof(int)*quadmesh->nverts);
	nboundverts=0;
}


void reset_narrowband()
{
	if(narrowband == NULL)
		//narrowband = new MinHeap(quadmesh->nverts);
		narrowband = new MinHeap(100);
	narrowband->reset();
}

void init_dis_cells()
{
	int i;
	for(i=0; i<quadmesh->nfaces; i++)
		quadmesh->quadcells[i]->OnBoundary=false;
}

/*This routine calculates the initial distances of the boundary vertices,
and save them into an array (this can be simply an integer array storing the indices
of the boundary vertices)
NOTE: we define boundary vertices as the vertices of the cells that the
curve passes through.*/
void cal_init_bound_verts()
{
	int i, j, k;
	int curve_pos=0;
	QuadCell *face;
	QuadVertex *v;
	double start[2], end[2];
	double distanceLine, disttoline;
	int pre_pos = 0;

	nboundverts=0;
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
		
			//for(k=curve_pos; k<num_curvepts_output; k++)
			//{
				//if(find_pos_curve(face->index, curve_pos, curve_pos))
				//{
				//}

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
							v->distance=distanceLine;
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
							v->distance=distanceLine;
					}

			//}

			//v->type=2;
		}
	}
}


/*
The problem is that the two consecutive points could be in two neighboring cells!!!
Therefore, when move to surfaces, problems may appear
10/10/2007
*/
bool find_pos_curve(int cellid, int startid, int &pos)
{
	int i;
	for(i=startid; i<num_curvepts_output; i++)
	{
		if(out_pts[i].cellid==cellid)
		{
			pos = i;
			return true;
		}
	}
	return false;
}

/*10/08/2007 initialize the intial boundary and the boundary vertices
using a straight line segment*/

void cal_init_bound_verts_test()
{
	int i, j;

	//use cell 12x12 and 30x30
	double start[2], end[2];
	int start_cell, end_cell;

	QuadCell *face = quadmesh->quadcells[20*20];
	start[0]=(quadmesh->quad_verts[face->verts[0]]->x
		+quadmesh->quad_verts[face->verts[1]]->x
		+quadmesh->quad_verts[face->verts[2]]->x
		+quadmesh->quad_verts[face->verts[3]]->x)/4;
	start[1]=(quadmesh->quad_verts[face->verts[0]]->y
		+quadmesh->quad_verts[face->verts[1]]->y
		+quadmesh->quad_verts[face->verts[2]]->y
		+quadmesh->quad_verts[face->verts[3]]->y)/4;
	start_cell=20*20;

	face = quadmesh->quadcells[40*40];
	end[0]=(quadmesh->quad_verts[face->verts[0]]->x
		+quadmesh->quad_verts[face->verts[1]]->x
		+quadmesh->quad_verts[face->verts[2]]->x
		+quadmesh->quad_verts[face->verts[3]]->x)/4;
	end[1]=(quadmesh->quad_verts[face->verts[0]]->y
		+quadmesh->quad_verts[face->verts[1]]->y
		+quadmesh->quad_verts[face->verts[2]]->y
		+quadmesh->quad_verts[face->verts[3]]->y)/4;
	end_cell=40*40;

	nboundarycells=0;
	find_boundarycells_oneline(start, start_cell, end,
								end_cell, boundarycells, nboundarycells);

	/*calculate the distance from the vertices of the boundary cells to
	the line segment*/
	//QuadCell *face;
	QuadVertex *v;
	double distanceLine;
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
			if(v->visited)
				continue;
			v->visited=true;

			DistanceFromLine(v->x, v->y,
				start[0],start[1], end[0],end[1], v->distance, distanceLine);

			boundvertlist[nboundverts]=v->index;
			nboundverts++;
		}
	}
}



/*the following code is borrowed from 
http://www.codeguru.com/forum/printthread.php?t=194400
*/
void DistanceFromLine(double cx, double cy, double ax, double ay ,
					  double bx, double by, double &distanceSegment,
					  double &distanceLine)
{

	//
	// find the distance from the point (cx,cy) to the line
	// determined by the points (ax,ay) and (bx,by)
	//
	// distanceSegment = distance from the point to the line segment
	// distanceLine = distance from the point to the line (assuming
	//					infinite extent in both directions
	//

/*

Subject 1.02: How do I find the distance from a point to a line?


    Let the point be C (Cx,Cy) and the line be AB (Ax,Ay) to (Bx,By).
    Let P be the point of perpendicular projection of C on AB.  The parameter
    r, which indicates P's position along AB, is computed by the dot product 
    of AC and AB divided by the square of the length of AB:
    
    (1)     AC dot AB
        r = ---------  
            ||AB||^2
    
    r has the following meaning:
    
        r=0      P = A
        r=1      P = B
        r<0      P is on the backward extension of AB
        r>1      P is on the forward extension of AB
        0<r<1    P is interior to AB
    
    The length of a line segment in d dimensions, AB is computed by:
    
        L = sqrt( (Bx-Ax)^2 + (By-Ay)^2 + ... + (Bd-Ad)^2)

    so in 2D:   
    
        L = sqrt( (Bx-Ax)^2 + (By-Ay)^2 )
    
    and the dot product of two vectors in d dimensions, U dot V is computed:
    
        D = (Ux * Vx) + (Uy * Vy) + ... + (Ud * Vd)
    
    so in 2D:   
    
        D = (Ux * Vx) + (Uy * Vy) 
    
    So (1) expands to:
    
            (Cx-Ax)(Bx-Ax) + (Cy-Ay)(By-Ay)
        r = -------------------------------
                          L^2

    The point P can then be found:

        Px = Ax + r(Bx-Ax)
        Py = Ay + r(By-Ay)

    And the distance from A to P = r*L.

    Use another parameter s to indicate the location along PC, with the 
    following meaning:
           s<0      C is left of AB
           s>0      C is right of AB
           s=0      C is on AB

    Compute s as follows:

            (Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay)
        s = -----------------------------
                        L^2


    Then the distance from C to P = |s|*L.

*/


	double r_numerator = (cx-ax)*(bx-ax) + (cy-ay)*(by-ay);
	double r_denomenator = (bx-ax)*(bx-ax) + (by-ay)*(by-ay);
	double r = r_numerator / r_denomenator;
//
    double px = ax + r*(bx-ax);
    double py = ay + r*(by-ay);
//     
    double s =  ((ay-cy)*(bx-ax)-(ax-cx)*(by-ay) ) / r_denomenator;

	distanceLine = fabs(s)*sqrt(r_denomenator);

//
// (xx,yy) is the point on the lineSegment closest to (cx,cy)
//
	double xx = px;
	double yy = py;

	if ( (r >= 0) && (r <= 1) )
	{
		distanceSegment = distanceLine;
	}
	else
	{

		double dist1 = (cx-ax)*(cx-ax) + (cy-ay)*(cy-ay);
		double dist2 = (cx-bx)*(cx-bx) + (cy-by)*(cy-by);
		if (dist1 < dist2)
		{
			xx = ax;
			yy = ay;
			distanceSegment = sqrt(dist1);
		}
		else
		{
			xx = bx;
			yy = by;
			distanceSegment = sqrt(dist2);
		}


	}

	return;
}



//void quadratic_solver(double a, double b, double c)
//{
//}
extern int solve_ten_quadratic(double a, double b, double c, double solutions[2]);

bool quadratic_2d(int vert, double &result)
{
	QuadVertex *v = quadmesh->quad_verts[vert];
	int right, left, up, bottom;
	bool useRow, useCol;
	double rowDis, colDis, Dis;
	double a, b, c;

	if(v->x>=quadmesh->xend-1.e-8) /*right boundary*/
	{
		left=vert-1;
		if(quadmesh->quad_verts[left]->type == 2)
		{
			useRow = true;
			rowDis = quadmesh->quad_verts[left]->distance;
		}
		else{
			useRow = false;
			rowDis = 0;
		}
	}
	else if(v->x<=quadmesh->xstart+1.e-8) /*right boundary*/
	{
		right=vert+1;
		if(quadmesh->quad_verts[right]->type == 2)
		{
			useRow = true;
			rowDis = quadmesh->quad_verts[right]->distance;
		}
		else{
			useRow = false;
			rowDis = 0;
		}
	}
	else
	{
		right = vert+1;
		left = vert-1;
		rowDis = min(quadmesh->quad_verts[left]->distance, quadmesh->quad_verts[right]->distance);

		if(quadmesh->quad_verts[left]->type==2 || quadmesh->quad_verts[right]->type==2)
			useRow = true;
		else{
			useRow = false;
			rowDis = 0;
		}
	}

	if(v->y>=quadmesh->yend-1.e-8)  /*upper boundar*/
	{
		bottom = vert-quadmesh->XDIM;
		if(quadmesh->quad_verts[bottom]->type == 2)
		{
			useCol = true;
			colDis = quadmesh->quad_verts[bottom]->distance;
		}
		else{
			useCol = false;
			colDis = 0;
		}
	}
	else if(v->y<=quadmesh->ystart+1.e-8) /*bottom boundary*/
	{
		up = vert+quadmesh->XDIM;
		if(quadmesh->quad_verts[up]->type == 2)
		{
			useCol = true;
			colDis = quadmesh->quad_verts[up]->distance;
		}
		else{
			useCol = false;
			colDis = 0;
		}
	}
	else
	{
		up = vert+quadmesh->XDIM;
		bottom = vert-quadmesh->YDIM;
		colDis = min(quadmesh->quad_verts[up]->distance, quadmesh->quad_verts[bottom]->distance);


		if(quadmesh->quad_verts[up]->type == 2 || quadmesh->quad_verts[bottom]->type == 2)
			useCol = true;
		else{
			useCol = false;
			colDis = 0;
		}
	}

	if((useRow && !useCol) || (useCol && !useRow))
	{
		Dis = max(rowDis, colDis);
		//Dis = min(rowDis, colDis);
		a = 1;
		b = -2*Dis;
		c = Dis*Dis-1;
	}
	else
	{
		a = 2;
		b = -(2*colDis)-(2*rowDis);
		c = rowDis*rowDis + colDis*colDis-1;
	}

	double solutions[2] = {0.};
	/*solve the quadratic function*/
	if(solve_ten_quadratic(a, b, c, solutions)==0)
	{
		return false;
	}
	else
	{
		result = solutions[0];
		return true;
	}
}




void update_oneVer(int cur_v, int next_v, double &result, bool &sflag)
{
	sflag = false;
	if(quadratic_2d(next_v, result))
	{
		if(result<quadmesh->quad_verts[next_v]->distance)
		{
			quadmesh->quad_verts[next_v]->distance = result;
			sflag = true;
		}
	}
	else  //use Dijkstra distance instead
	{
		double temp_dis = quadmesh->quad_verts[cur_v]->distance+
			quadmesh->xinterval;
		if(temp_dis<quadmesh->quad_verts[next_v]->distance)
		{
			quadmesh->quad_verts[next_v]->distance=temp_dis;
			sflag = true;
		}
	}

	if(quadmesh->quad_verts[next_v]->type==0) /*avoid repeated*/
	{
		if(!narrowband->Insert(next_v, quadmesh->quad_verts[next_v]->distance))
		{
			/*exit(-1)*/return;
		}
		quadmesh->quad_verts[next_v]->type = 1;
	}

	else if(quadmesh->quad_verts[next_v]->type==1 && sflag)
	{
		/*we need to update the corresponding Min Heap*/
		int pos = narrowband->get_pos(next_v);

		if(pos < 0) return;

		narrowband->elems[pos]->T = quadmesh->quad_verts[next_v]->distance;
		narrowband->UpHeap(pos); /*the updating always keeps the distance minimum*/
	}
}


void init_narrow_band()
{
	int i, cur_v, next_v;
	double result;
	bool sflag = false;

	/*mark the boundary vertices as "known"*/
	for(i=0; i<nboundverts; i++)
	{
		quadmesh->quad_verts[boundvertlist[i]]->type = 2;
	}
	
	QuadVertex *v;
	for(i=0; i<nboundverts; i++)
	{
		/*propagate the distance to the one ring neighborhood of each of these "known" points*/
		cur_v = boundvertlist[i];

		v=quadmesh->quad_verts[cur_v];
		//if(cur_v%quadmesh->XDIM==0) /*left boundary*/
		if(v->x<=quadmesh->xstart+1.e-8) /*left boundary*/
		{
			/*check right only*/
			next_v = cur_v+1;

			if(quadmesh->quad_verts[next_v]->type!=2)
				update_oneVer(cur_v, next_v, result, sflag);

		}
		//else if((cur_v+1)%quadmesh->XDIM == 0) /*right boundary*/
		else if(v->x>=quadmesh->xend-1.e-8) /*right boundary*/
		{
			/*check left only*/
			next_v = cur_v-1;
			if(quadmesh->quad_verts[next_v]->type!=2)
				update_oneVer(cur_v, next_v, result, sflag);
		}
		else /*check both left and right*/
		{
			/*check right*/
			next_v = cur_v+1;
			if(quadmesh->quad_verts[next_v]->type!=2)
				update_oneVer(cur_v, next_v, result, sflag);
			
			/*check left*/
			next_v = cur_v-1;
			if(quadmesh->quad_verts[next_v]->type!=2)
				update_oneVer(cur_v, next_v, result, sflag);
		}

		//if(cur_v/quadmesh->YDIM == 0)  //bottom
		if(v->y <=quadmesh->ystart+1.e-8)  //bottom
		{
			/*check above only*/
			next_v = cur_v+quadmesh->XDIM;
			if(quadmesh->quad_verts[next_v]->type!=2)
				update_oneVer(cur_v, next_v, result, sflag);
		}
		//else if(cur_v/(quadmesh->YDIM-1) == 0)  //upper
		else if(v->y>=quadmesh->yend-1.e-8)  //upper
		{
			/*check below only*/
			next_v = cur_v-quadmesh->XDIM;
			if(quadmesh->quad_verts[next_v]->type!=2)
				update_oneVer(cur_v, next_v, result, sflag);
		}
		else /*check both above and below*/
		{
			/*check above*/
			next_v = cur_v+quadmesh->XDIM;
			if(quadmesh->quad_verts[next_v]->type!=2)
				update_oneVer(cur_v, next_v, result, sflag);
			
			/*check below*/
			next_v = cur_v-quadmesh->XDIM;
			if(quadmesh->quad_verts[next_v]->type!=2)
				update_oneVer(cur_v, next_v, result, sflag);
		}

		/*mark these neighbors as "in the narrow band", insert them into the min-heap accordingly*/

		/*mark the remaining points of the mesh as "far away"*/
	}

}


/*Update the distances of the neighbors of the given vertex*/
void update_neighbor_Dis(int cur_v)
{
	int next_v;
	double result;
	QuadVertex *v = quadmesh->quad_verts[cur_v];
	bool sflag = false;

	//if(cur_v%quadmesh->XDIM==0) /*left boundary*/
	if(v->x<=quadmesh->xstart+1.e-8) /*left boundary*/
	{
		/*check right only*/
		next_v = cur_v+1;
		if(quadmesh->quad_verts[next_v]->type!=2)
			update_oneVer(cur_v, next_v, result, sflag);

	}
	//else if((cur_v+1)%quadmesh->XDIM == 0) /*right boundary*/
    else if(v->x>=quadmesh->xend-1.e-8)	
	{
		/*check left only*/
		next_v = cur_v-1;
		if(quadmesh->quad_verts[next_v]->type!=2)
			update_oneVer(cur_v, next_v, result, sflag);
	}
	else /*check both left and right*/
	{
		/*check right*/
		next_v = cur_v+1;
		if(quadmesh->quad_verts[next_v]->type!=2)
			update_oneVer(cur_v, next_v, result, sflag);
		
		/*check left*/
		next_v = cur_v-1;
		if(quadmesh->quad_verts[next_v]->type!=2)
			update_oneVer(cur_v, next_v, result, sflag);
	}

	//if(cur_v/quadmesh->YDIM == 0)  //bottom
	if(v->y<=quadmesh->ystart+1.e-8)  //bottom
	{
		/*check above only*/
		next_v = cur_v+quadmesh->XDIM;
		if(quadmesh->quad_verts[next_v]->type!=2)
			update_oneVer(cur_v, next_v, result, sflag);
	}
	//else if(cur_v/(quadmesh->YDIM-1) == 0)  //upper
	else if(v->y>=quadmesh->yend-1.e-8)  //upper
	{
		/*check below only*/
		next_v = cur_v-quadmesh->XDIM;
		if(quadmesh->quad_verts[next_v]->type!=2)
			update_oneVer(cur_v, next_v, result, sflag);
	}
	else /*check both above and below*/
	{
		/*check above*/
		next_v = cur_v+quadmesh->XDIM;
		if(quadmesh->quad_verts[next_v]->type!=2)
			update_oneVer(cur_v, next_v, result, sflag);
		
		/*check below*/
		next_v = cur_v-quadmesh->XDIM;
		if(quadmesh->quad_verts[next_v]->type!=2)
			update_oneVer(cur_v, next_v, result, sflag);
	}
}

/*After given the boundary vertices with known values on them,
we start the process of fast marching
NOTE: we assume that the boundary vertices have been added to the array and
the distances on them have also been computed
*/
void fast_marching_quad()
{
	/*first, initialization*/
	init_dis_verts();

	//cal_init_bound_verts_test(); /*testing code 10/08/2007*/
	cal_init_bound_verts();

	init_narrow_band();


	int nfinished = nboundverts;
	while(!narrowband->is_empty() && nfinished < quadmesh->nverts)
	{
		int cur_v = narrowband->FindSmallest();

		quadmesh->quad_verts[cur_v]->type = 2;

			boundvertlist[nboundverts]=cur_v;
			nboundverts++;

		update_neighbor_Dis(cur_v);
		nfinished ++;
	}
}


void fast_marching_quad_withDis(double disthred)
{
	/*first, initialization*/
	//init_dis_verts();
	reset_vert_dis();

	//cal_init_bound_verts_test(); /*testing code 10/08/2007*/
	cal_init_bound_verts();

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

		boundvertlist[nboundverts]=cur_v;
		nboundverts++;
		
		/*update neighbor distances of cur_v*/
		if(quadmesh->quad_verts[cur_v]->distance<disthred)
			update_neighbor_Dis(cur_v);

		nfinished ++;
	}
}



extern void  HsvRgb( float hsv[3], float rgb[3] );

/*use color coding to visualize the distance calculation result*/
void vis_distance()
{
	double fartest, closest;
	int i,j;
	QuadVertex *v;
	fartest=closest=quadmesh->quad_verts[0]->distance;
	for(i=1; i<quadmesh->nverts; i++)
	{
		v=quadmesh->quad_verts[i];
		if(v->distance>fartest) fartest=v->distance;
		if(v->distance<closest) closest=v->distance;
	}

	QuadCell *face;
	glShadeModel(GL_SMOOTH);

	float hsv[3]={0.,1.,1.};
	float rgb[3]={0.};
	for(i=0; i<quadmesh->nfaces; i++)
	{
		face=quadmesh->quadcells[i];

		glBegin(GL_POLYGON);
		for(j=0; j<face->nverts; j++)
		{
			v=quadmesh->quad_verts[face->verts[j]];
			hsv[0]=240.*(v->distance-closest)/(fartest-closest);
			HsvRgb(hsv, rgb);
			glColor3fv((GLfloat*)rgb);
			glVertex2f(v->x, v->y);
		}
		glEnd();
	}

}

extern int *region_quadverts;
extern int nregion_quadverts;

void display_zero_level()
{
	int i;
	QuadVertex *v;
	glColor3f(1, 0, 1);
	glPointSize(3.);
	glBegin(GL_POINTS);
	for(i=0; i<nboundverts; i++)
	{
		v=quadmesh->quad_verts[boundvertlist[i]];
		glVertex2f(v->x, v->y);
	}
	glEnd();

	glColor3f(0, 0, 0);
	glBegin(GL_POINTS);
	for(i=0; i<narrowband->nelems; i++)
	{
		v=quadmesh->quad_verts[narrowband->elems[i]->vertid];
		glVertex2f(v->x, v->y);
	}
	glEnd();
	
	glColor3f(1, 1, 0);
	//glBegin(GL_POINTS);
	//for(i=0; i<quadmesh->nverts; i++)
	//{
	//	v=quadmesh->quad_verts[i];
	//	if(!v->InRegion) continue;
	//	glVertex2f(v->x, v->y);
	//}
	//glEnd();

	/*visualize the inner vertices*/
	glBegin(GL_POINTS);
	for(i=0; i<nregion_quadverts; i++)
	{
		v=quadmesh->quad_verts[region_quadverts[i]];
		if(!v->InRegion) continue;
		glVertex2f(v->x, v->y);
	}
	glEnd();
}



/**************************************************************/
/**************************************************************/
/**************************************************************/
/**************************************************************/


int get_cellID_picking(double x, double y)
{
	GLuint  selectBuffer[SELECTBUFFERSIZE];
	GLint	vp[4] = {0, 0 ,1, 1};
	int hits;
	int selectedTriangle = -1;


	////Build the selection buffer here
	glSelectBuffer(SELECTBUFFERSIZE, selectBuffer);

   	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);

	gluPickMatrix(x, y, 1e-9, 1e-9, vp );
		
	glOrtho(0,1,  0,1,  0,50);

	////If one of the element has been selected for being edited

	display_quadmesh_select(GL_SELECT);

	hits = glRenderMode(GL_RENDER);

	if(hits>0)
	{
		////Because this is 2D plane, we need not worry about the overlap of objects
		////It will have at most one object being selected at one click (it may not be true)
		int objname = selectBuffer[3];

		selectedTriangle = objname - NAMEOFTRIANGLE;

	}

	else{
		selectedTriangle = -1;
	}

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode( GL_MODELVIEW );

	return selectedTriangle;
}


/*
   The following routine may contain bug
*/

int get_cellID_givencoords(double x, double y)
{
	int i=(x-quadmesh->xstart)/quadmesh->xinterval;
	int j=(y-quadmesh->ystart)/quadmesh->yinterval;

	//if(i<0) i=0;
	//if(j<0) j=0;

	if(i==quadmesh->XDIM-1) i=quadmesh->XDIM-2;
	if(j==quadmesh->YDIM-1) j=quadmesh->YDIM-2;

	return (j*(quadmesh->XDIM-1)+i);
}


/*
The function deal with the curve representation rather than
a single line segment. The curve is stored in the variable

ctr_point *out_pts;

*/
//void find_boundarycells_curve_quad(double start[2], int start_cell, double end[2],
//								int end_cell, int *celllist, int &ncells)
void get_cellstrip_curve_quad(int *celllist, int &ncells)
{
	double pre_p[2]={out_pts[0].x, out_pts[0].y};
	double cur_p[2]={out_pts[1].x, out_pts[1].y};
	int i;
	int cur_cell /*= start_cell*/;
	QuadCell *face;
	icVector2 linedir;
	linedir.set((cur_p[0]-pre_p[0]),(cur_p[1]-pre_p[1]));
	normalize(linedir);
	tenline_dir_global = linedir;

	/*add the first cell*/
	//celllist[ncells]=out_pts[0].cellid;
	/*or*/
	//cur_cell=celllist[ncells]=get_cellID_picking(pre_p[0], pre_p[1]);get_cellID_givencoords
	cur_cell=celllist[ncells]=get_cellID_givencoords(pre_p[0], pre_p[1]);
	ncells=1;
	quadmesh->quadcells[cur_cell]->OnBoundary=true;
	out_pts[0].cellid=cur_cell;

	icVector2 t_major[4];
	int count = 1;

	while(count <= resolution)
	{
		if(!is_in_cell(cur_cell, cur_p[0], cur_p[1]))
		{
			/*find the next cell the curve will enter*/

			double temp_p[2]={cur_p[0], cur_p[1]};

			/*we first store the major vectors of the vertices
			of current cell*/
			face = quadmesh->quadcells[cur_cell];
			for(i=0; i<4; i++)
				t_major[i]=quadmesh->quad_verts[face->verts[i]]->major;

			/*replace them with the line segment direction*/
			for(i=0; i<4; i++)
				quadmesh->quad_verts[face->verts[i]]->major=linedir;

			int passvertornot = 0;
			//get_next_cell(cur_cell, pre_p, cur_p, passvertornot, 0);
			get_next_cell_2(cur_cell, pre_p, cur_p, passvertornot, 0);

			if(cur_cell<0 || cur_cell>quadmesh->nfaces)
				return;

			/*add to the cell list*/
			if(!quadmesh->quadcells[cur_cell]->OnBoundary)
			{
				celllist[ncells]=cur_cell;
				ncells++;

				quadmesh->quadcells[cur_cell]->OnBoundary=true;
			}

			if(passvertornot == 0)
			{
				pre_p[0] = cur_p[0];
				pre_p[1] = cur_p[1];
			}
			cur_p[0] = temp_p[0];
			cur_p[1] = temp_p[1];

			linedir.set((cur_p[0]-pre_p[0]),(cur_p[1]-pre_p[1]));
			normalize(linedir);

			/*store back the original vectors*/
			for(i=0; i<4; i++)
				quadmesh->quad_verts[face->verts[i]]->major=t_major[i];
		}
		else
		{
			cur_p[0] = out_pts[count].x;
			cur_p[1] = out_pts[count].y;
			out_pts[count-1].cellid=cur_cell;
			count++;
		}
	}
	
	/*as a remedy, we force the finding of the cellid*/
	//for(i=1; i<resolution; i++)
	//	out_pts[i].cellid=get_cellID_picking(out_pts[i].x, out_pts[i].y);
}


int get_Resolution_adp()
{
	double edgelen;
	int resolute = (int) 2* GetWholeLengthofDesignCurve()/quadmesh->xinterval;

	return resolute;
}


////Add the user selected point into the points list
void add_to_shapeCtrPtsList(double x, double y, int cellid)
{
	//int vert;

	if(num_shapecontrol_pts >= MaxNumShapeControlPts)
	{
		MessageBox(NULL, "Can not add more control points!", "Error", MB_OK);
		return;
	}
	control_pts[num_shapecontrol_pts].x = x;
	control_pts[num_shapecontrol_pts].y = y;
	control_pts[num_shapecontrol_pts].cellid = cellid;


	num_shapecontrol_pts ++;
}


////
/*
In this routine, we use the control points (not the output shape points) to approximate 
the direction of the boundary curve
*/
void get_bound_approDir()
{
	int i;
	//for(i=0; i<num_shapecontrol_pts-1; i++)
	//{
	//	set_ten_regBasis(control_pts[i].x, control_pts[i].y, 0);
	//	set_ten_regDir(control_pts[i+1].x, control_pts[i+1].y);
	//}
	
	for(i=0; i<resolution-1; i++)
	{
		set_ten_regBasis(out_pts[i].x, out_pts[i].y, 0);
		set_ten_regDir(out_pts[i+1].x, out_pts[i+1].y);
	}
}


/*
Get the tensor field in the interior of the brush
*/
void get_brushinterior_ten()
{
	int i;
	double t[4]={0.};
	QuadVertex *v;
	for(i=0; i<quadmesh->nverts; i++)
	{
		v=quadmesh->quad_verts[i];
		//v->Jacobian.set(0.);
		//v->major.set(0,0);
		//v->minor.set(0,0);

		if(v->type!=2 || v->InRegion || !v->inland) /*consider the geographics map*/
			continue;

		get_tensor(v->x, v->y, t);

		/*we also need to normalize the tensor here*/
		v->Jacobian.entry[0][0]=t[0]/50.;
		v->Jacobian.entry[0][1]=t[1]/50.;
		v->Jacobian.entry[1][0]=t[2]/50.;
		v->Jacobian.entry[1][1]=t[3]/50.;


		/*get the major and minor*/
		cal_eigenvecs_onevert_quad(i);
	}

	
	/*normalize the major and minor field*/
	normalized_tensorfield_quad();
			
	init_degpts();
	render_alpha_map_quad(false);
	render_alpha_map_quad(true);


	QuadCell *face;
	int num_nonzeroverts = 0;
	for(i=0; i<quadmesh->nfaces; i++)
	{
		num_nonzeroverts = 0;
		face=quadmesh->quadcells[i];

		for(int j=0; j<face->nverts; j++)
		{
			if(quadmesh->quad_verts[face->verts[j]]->type==2)
				num_nonzeroverts++;
		}

		if(num_nonzeroverts > 2)
			face->OnBoundary = true;
	}

	for(i=0; i<quadmesh->nverts; i++)
		quadmesh->quad_verts[i]->type=0;
}

bool closedbrush=false;

/*
save the brush curve into a file
*/
void save_brush_curve(char *filename)
{
	int i;
	FILE *fp = fopen(filename, "w");
	
	////write control points
	fprintf(fp, "Control points:%d\n", num_shapecontrol_pts);
	fprintf(fp, "type:%d\n", 0);
	if(closedbrush)
		fprintf(fp, "closed:%d\n", 1);
	else
		fprintf(fp, "closed:%d\n", 0);
	for(i = 0; i < num_shapecontrol_pts; i++)
	{
		fprintf(fp, "%f,%f\n", control_pts[i].x, control_pts[i].y);
	}
	fclose(fp);
}

/*
load the brush curve
*/
void load_brush_curve(char *filename)
{
	int i;
	int type;
	int brushclosed;
	FILE *fp = fopen(filename, "r");

	if(fp == NULL)
	{
		MessageBox(NULL, "Can not open file!", "Error", MB_OK);
		return;
	}

	//int num_shapecontrol_pts;
	float x, y;

	////write control points
	fscanf(fp, "Control points:%d\n", &num_shapecontrol_pts);
	fscanf(fp, "type:%d\n", &type);
	fscanf(fp, "closed:%d\n", &brushclosed);
	for(i = 0; i < num_shapecontrol_pts; i++)
	{
		fscanf(fp, "%f,%f\n", &x, &y);
		control_pts[i].x = x;
		control_pts[i].y = y;
	}
	fclose(fp);
		
	resolution = get_Resolution_adp();

	if(brushclosed == 0)
	{
		CalOpenHermiteCurve();
	}
	else
	{
		CalHermiteCurve();
	}
}







/**********************************************************************************/
/* The following we are trying to implement a functionality that allow 
the user to merge the brush stroke to the background tensor if there is one 10/20/2007*/
/*
obtain a buffer region around the brush stroke 10/20/2007
*/
void get_brushbuffer_quad_withDis(double disthred)
{
	int i;
	nboundverts = 0;
	for(i=0; i<quadmesh->nverts; i++)
	{
		quadmesh->quad_verts[i]->InRegion = false;
		quadmesh->quad_verts[i]->OnBoundary = true;
		quadmesh->quad_verts[i]->RegionListID = -1;
		if(quadmesh->quad_verts[i]->type != 2)
			continue;

		boundvertlist[nboundverts]=i;
		nboundverts++;
	}

	init_narrow_band();

	int nfinished = nboundverts;

	nregion_quadverts = 0;

	/*we make use of the idea of fast marching but setting a maximal propagated distance
	when the distance exceeds this threshold, we stop propagation*/

	/*
	To obtain the smoothing region, we also mark the obtained vertices accordingly!
	*/

	while(!narrowband->is_empty() && nfinished < quadmesh->nverts)
	{
		int cur_v = narrowband->FindSmallest();
		quadmesh->quad_verts[cur_v]->type = 2;

		/*save it to the smoothing region*/
		if(quadmesh->quad_verts[cur_v]->inland) /*consider the geographics map 10/24/2007*/
		{
			region_quadverts[nregion_quadverts] = cur_v;
			quadmesh->quad_verts[region_quadverts[nregion_quadverts]]->RegionListID 
				= nregion_quadverts;
			quadmesh->quad_verts[cur_v]->InRegion = true;
			quadmesh->quad_verts[cur_v]->OnBoundary = false;
			nregion_quadverts++;
		}
		/*********************/

		boundvertlist[nboundverts]=cur_v;
		nboundverts++;
		
		/*update neighbor distances of cur_v*/
		if(quadmesh->quad_verts[cur_v]->distance<disthred)
			update_neighbor_Dis(cur_v);

		nfinished ++;
	}
}



/*
We need to save the previous brush information  10/20/2007
*/
DynList_Int *vertsinbrush = NULL;  //the vertices falling in previous specified brushes!

void save_cur_brushverts()
{
	int i;
	if(vertsinbrush == NULL)
		return;

	for(i=0; i<quadmesh->nverts; i++)
	{
		if(quadmesh->quad_verts[i]->type!=2)
			continue;

		vertsinbrush->add_New(i);
	}
}






/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
/*              Routines for finding all the blocks             */

/*
This routine initializes all the edges in the street graph
*/

RegionBlockList *regionblocklist=NULL;

void init_regionblocklist()
{
	if(regionblocklist!=NULL)
		delete regionblocklist;
	regionblocklist=new RegionBlockList();
}

void init_all_streetgraph_edges()
{
	int i;
	StreetGraphEdge *cur_sedge;
	if(streetnet==NULL) return;
	for(i=0; i<streetnet->edgelist->nedges;i++)
	{
		cur_sedge=streetnet->edgelist->edges[i];
		cur_sedge->nregionblocks=0;
	}
}

bool select_a_valid_edge(int nodeid, icVector2 normal, int &nextedge, int exceptedge)
{
	Intersection *pnode=streetnet->nodelist->intersects[nodeid];
	StreetGraphEdge *cur_e;
	double sim=-50.;
	icVector2 line_dir;
	int i;
	nextedge=-1;

	if(pnode->nadjedges==2)
	{
		for(i=0;i<pnode->nadjedges;i++)
		{
			cur_e=streetnet->edgelist->edges[pnode->adj_edges[i]];
			if(pnode->adj_edges[i]==exceptedge)
				continue;
			if(streetnet->nodelist->intersects[cur_e->node_index1]->endpt
				||streetnet->nodelist->intersects[cur_e->node_index2]->endpt)
				continue;
			if(cur_e->nregionblocks==2)
				continue;

			nextedge=pnode->adj_edges[i];
			return true;
		}

		return false;
	}

	for(i=0;i<pnode->nadjedges;i++)
	{
		if(pnode->adj_edges[i]==exceptedge)
			continue;

		cur_e=streetnet->edgelist->edges[pnode->adj_edges[i]];

		if(cur_e->nregionblocks==2)
			continue;

		/*we don't consider danggling edge here as well*/
		if(streetnet->nodelist->intersects[cur_e->node_index1]->endpt
			||streetnet->nodelist->intersects[cur_e->node_index2]->endpt)
			continue;

		if(cur_e->node_index1==nodeid)
		{
			line_dir.entry[0]=streetnet->nodelist->intersects[cur_e->node_index2]->gpos[0]
				-streetnet->nodelist->intersects[cur_e->node_index1]->gpos[0];
			line_dir.entry[1]=streetnet->nodelist->intersects[cur_e->node_index2]->gpos[1]
				-streetnet->nodelist->intersects[cur_e->node_index1]->gpos[1];
		}
		else
		{
			line_dir.entry[0]=streetnet->nodelist->intersects[cur_e->node_index1]->gpos[0]
				-streetnet->nodelist->intersects[cur_e->node_index2]->gpos[0];
			line_dir.entry[1]=streetnet->nodelist->intersects[cur_e->node_index1]->gpos[1]
				-streetnet->nodelist->intersects[cur_e->node_index2]->gpos[1];
		}

		normalize(line_dir);

		double temp=dot(line_dir, normal);
		if(temp>sim)
		{
			sim=temp;
			nextedge=pnode->adj_edges[i];
		}
	}

	if(nextedge==-1)
		return false;
	else
		return true;
}


bool select_a_valid_edge_anglebased(int nodeid, icVector2 curdir, int &nextedge, int exceptedge)
{
	Intersection *pnode=streetnet->nodelist->intersects[nodeid];
	StreetGraphEdge *cur_e;
	double smallest_ag=50.;
	icVector2 line_dir;
	int i;
	nextedge=-1;

	curdir=-curdir;
	double ag1=atan2(curdir.entry[1], curdir.entry[0]);
	if(ag1<0) ag1+=2*M_PI;
	double ag2, ag_diff;

	if(pnode->nadjedges==2)
	{
		for(i=0;i<pnode->nadjedges;i++)
		{
			cur_e=streetnet->edgelist->edges[pnode->adj_edges[i]];
			if(pnode->adj_edges[i]==exceptedge)
				continue;
			if(streetnet->nodelist->intersects[cur_e->node_index1]->endpt
				||streetnet->nodelist->intersects[cur_e->node_index2]->endpt)
				continue;
			if(cur_e->nregionblocks==2)
				continue;

			nextedge=pnode->adj_edges[i];
			return true;
		}

		return false;
	}

	for(i=0;i<pnode->nadjedges;i++)
	{
		if(pnode->adj_edges[i]==exceptedge)
			continue;

		cur_e=streetnet->edgelist->edges[pnode->adj_edges[i]];

		if(cur_e->nregionblocks==2)
			continue;

		/*we don't consider danggling edge here as well*/
		if(streetnet->nodelist->intersects[cur_e->node_index1]->endpt
			||streetnet->nodelist->intersects[cur_e->node_index2]->endpt)
			continue;

		if(cur_e->node_index1==nodeid)
		{
			line_dir.entry[0]=streetnet->nodelist->intersects[cur_e->node_index2]->gpos[0]
				-streetnet->nodelist->intersects[cur_e->node_index1]->gpos[0];
			line_dir.entry[1]=streetnet->nodelist->intersects[cur_e->node_index2]->gpos[1]
				-streetnet->nodelist->intersects[cur_e->node_index1]->gpos[1];
		}
		else
		{
			line_dir.entry[0]=streetnet->nodelist->intersects[cur_e->node_index1]->gpos[0]
				-streetnet->nodelist->intersects[cur_e->node_index2]->gpos[0];
			line_dir.entry[1]=streetnet->nodelist->intersects[cur_e->node_index1]->gpos[1]
				-streetnet->nodelist->intersects[cur_e->node_index2]->gpos[1];
		}

		//normalize(line_dir);
		ag2=atan2(line_dir.entry[1], line_dir.entry[0]);
		if(ag2<0) ag2+=2*M_PI;

		ag_diff=ag1-ag2;
		if( ag_diff < -M_PI)
			ag_diff += 2 * M_PI;
		
		//if( ag_diff > M_PI)
		//	ag_diff -= 2 * M_PI;

		if(ag_diff>0&&ag_diff<smallest_ag)
		{
			smallest_ag=ag_diff;
			nextedge=pnode->adj_edges[i];
		}
	}

	if(nextedge==-1)
		return false;
	else
		return true;
}


/*
After obtaining the lists of edges and nodes that form a region block,
we call the following routine to store them into the corresponding data structure
*/
void construct_a_regionblock(int *nodes, int nnodes, int *edges, int nedges)
{
	int i;

	/*create a new block*/
	RegionBlock *block=(RegionBlock*)malloc(sizeof(RegionBlock));
	block->edgelist=(int *)malloc(sizeof(int)*nedges);
	block->nodelist=(int *)malloc(sizeof(int)*nnodes);

	for(i=0; i<nnodes; i++)
	{
		block->nodelist[i]=nodes[i];
	}
	block->num_nodes=nnodes;

	for(i=0; i<nedges; i++)
	{
		block->edgelist[i]=edges[i];
		streetnet->edgelist->edges[edges[i]]->regionblocks[
			streetnet->edgelist->edges[edges[i]]->nregionblocks]=regionblocklist->nelems;
		streetnet->edgelist->edges[edges[i]]->nregionblocks++;
	}
	block->num_edges=nedges;

	block->closed_curved=false;
	block->trajid=-1;

	/*add to the block list*/
	regionblocklist->append(block);
}



/*
search a block from current input edge based on the specified orientation
orient: false--follow the edge direction;  true--follow the inverse edge direction
*/
bool form_a_block(int edgeid, bool orient)
{
	StreetGraphEdge *start_sedge=streetnet->edgelist->edges[edgeid];
	StreetGraphEdge *cur_sedge=start_sedge;
	int startnode, endnode, cur_node;
	Intersection *pendnode, *pcur_node;
	icVector2 line_dir;
	icVector2 cur_normal;
	int next_edge;

	int curMaxSelectedEdge=10;
	int *sel_edges=(int*)malloc(sizeof(int)*curMaxSelectedEdge);
	int *sel_nodes=(int*)malloc(sizeof(int)*curMaxSelectedEdge);
	sel_edges[0]=edgeid;
	int nsel_edges=1;

	
	if(!orient) /*it is following the edge direction*/
	{
		endnode=start_sedge->node_index1;
		startnode=cur_node=start_sedge->node_index2;
		sel_nodes[0]=cur_node;
	}

	else /*following the inverse edge direction*/
	{
		endnode=start_sedge->node_index2;
		startnode=cur_node=start_sedge->node_index1;
		sel_nodes[0]=cur_node;
	}

	pendnode=streetnet->nodelist->intersects[endnode];
	pcur_node=streetnet->nodelist->intersects[cur_node];
	line_dir.entry[0]=pcur_node->gpos[0]-pendnode->gpos[0];
	line_dir.entry[1]=pcur_node->gpos[1]-pendnode->gpos[1];
	/*obtain the inner normal*/
	cur_normal.entry[0]=-line_dir.entry[1];
	cur_normal.entry[1]=line_dir.entry[0];
	normalize(cur_normal);

	int cur_edgeid=edgeid;

	while(1)
	{
		//if(!select_a_valid_edge(cur_node, cur_normal, next_edge, cur_edgeid))
		if(!select_a_valid_edge_anglebased(cur_node, line_dir, next_edge, cur_edgeid))
		/*if we run out of all the choices*/
		{
			free(sel_edges);
			free(sel_nodes);
			return false;
		}

		cur_edgeid=next_edge;
		
		/*add the next_edge to the sel_edges list*/
		if(nsel_edges>=curMaxSelectedEdge)
		{
			sel_edges=(int*)realloc(sel_edges, sizeof(int)*(curMaxSelectedEdge+10));
			sel_nodes=(int*)realloc(sel_nodes, sizeof(int)*(curMaxSelectedEdge+10));
			curMaxSelectedEdge+=10;
		}

		if(is_repeated_elem(sel_edges, next_edge, nsel_edges))
		{
			free(sel_edges);
			free(sel_nodes);
			return false;
		}

		sel_edges[nsel_edges]=next_edge;
		nsel_edges++;

		/*update the normal*/
		/*find the next node to continue on*/
		cur_sedge=streetnet->edgelist->edges[next_edge];
		pcur_node=streetnet->nodelist->intersects[cur_node];
		if(cur_node==cur_sedge->node_index1) /*we should follow the edge direction*/
		{
			pendnode=streetnet->nodelist->intersects[cur_sedge->node_index2];
			cur_sedge->inversed_orient=false;
			cur_node=cur_sedge->node_index2;
		}
		else if(cur_node==cur_sedge->node_index2)
		{
			pendnode=streetnet->nodelist->intersects[cur_sedge->node_index1];
			cur_sedge->inversed_orient=true;
			cur_node=cur_sedge->node_index1;
		}
		else
			return false;

		sel_nodes[nsel_edges-1]=cur_node;

		if(cur_node==endnode) /*we find a block*/
		{
			/*construct a new block*/
			construct_a_regionblock(sel_nodes, nsel_edges, sel_edges, nsel_edges);

			free(sel_edges);
			free(sel_nodes);
			return true;
		}

		line_dir.entry[0]=pendnode->gpos[0]-pcur_node->gpos[0];
		line_dir.entry[1]=pendnode->gpos[1]-pcur_node->gpos[1];

		cur_normal.entry[0]=-line_dir.entry[1];
		cur_normal.entry[1]=line_dir.entry[0];

		normalize(cur_normal);
	}
}

				  
/*
construct block(s) for an input edge
Return: 1 block if the edge has become one edge of an existing block
        2 blocks if the edge has never be included by existing block
		0 block if this is an daggling edge
*/

int construct_blocks_for_an_edge(int edgeid)
{
	int obtained_blocks=0;

	StreetGraphEdge *cur_sedge=streetnet->edgelist->edges[edgeid];

	if(cur_sedge->nregionblocks==0) /*never be included by existing block*/
	{
		cur_sedge=streetnet->edgelist->edges[edgeid];
		cur_sedge->inversed_orient=false;
		if(form_a_block(edgeid, cur_sedge->inversed_orient))
			obtained_blocks++;

		cur_sedge->inversed_orient=true;
		if(form_a_block(edgeid, cur_sedge->inversed_orient))
			obtained_blocks++;
	}

	else if(cur_sedge->nregionblocks==1)/*has become one edge of an existing block*/
	{
		/*we need to figure out which orientation it has been used in that block*/
		cur_sedge->inversed_orient=!cur_sedge->inversed_orient;
		if(form_a_block(edgeid, cur_sedge->inversed_orient))
			obtained_blocks++;
		else
		{
			if(form_a_block(edgeid, !cur_sedge->inversed_orient))
				obtained_blocks++;
		}
	}

	return obtained_blocks;
}

/*
construct the regions by traversing the edge list
*/
void construct_regionblocks_edgewise()
{
	if(streetnet==NULL) return;

	init_regionblocklist();
	init_all_streetgraph_edges();
	int i;
	StreetGraphEdge *cur_sedge;

	for(i=0;i<streetnet->edgelist->nedges;i++)
	{
		cur_sedge=streetnet->edgelist->edges[i];
		if(cur_sedge->nregionblocks==2)
			continue;

		if(streetnet->nodelist->intersects[cur_sedge->node_index1]->endpt
			||streetnet->nodelist->intersects[cur_sedge->node_index2]->endpt)
			continue;

		/*cur_sedge->nregionblocks+=*/construct_blocks_for_an_edge(i);
	}
}

extern void getcolor_scc_frac(int num, int frac, float rgb[3]);

void vis_regionblocks()
{
	if(regionblocklist==NULL) return;

	int i, j;
	float rgb[3]={0.};
	Intersection *cur_n;
	for(i=0;i<regionblocklist->nelems;i++)
	{
		//if(regionblocklist->blocks[i]->num_nodes>5)
		//	continue;

		getcolor_scc_frac(i, 9, rgb);
		glColor3fv(rgb);
		glBegin(GL_POLYGON);
		for(j=0;j<regionblocklist->blocks[i]->num_nodes; j++)
		{
			cur_n=streetnet->nodelist->intersects[regionblocklist->blocks[i]->nodelist[j]];
			glVertex2f(cur_n->gpos[0], cur_n->gpos[1]);
		}
		glEnd();
	}

	/*visualize those edges having only one block sharing them*/
	StreetGraphEdge *cur_sedge;
	glLineWidth(3.);
	glColor3f(1,1,1);
	for(i=0;i<streetnet->edgelist->nedges;i++)
	{
		cur_sedge=streetnet->edgelist->edges[i];

		if(cur_sedge->nregionblocks<2)
		{
			Intersection *n1, *n2;
			n1=streetnet->nodelist->intersects[cur_sedge->node_index1];
			n2=streetnet->nodelist->intersects[cur_sedge->node_index2];
			glBegin(GL_LINES);
				glVertex2f(n1->gpos[0], n1->gpos[1]);
				glVertex2f(n2->gpos[0], n2->gpos[1]);
			glEnd();
		}

	}
	glLineWidth(1.);
}


