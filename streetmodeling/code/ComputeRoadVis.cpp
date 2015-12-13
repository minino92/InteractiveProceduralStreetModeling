/*
This file implements the computation of the road network for visualization purpose 
NOTE: the data structure has combined with current data structure
if the data structure is changed, according changes are necessary.
10/03/2007
*/
#include "stdafx.h"

#include "VFDataStructure.h"

#include "tensorvis.h"

#include "tensoranalysis.h"

#include "tensordesign.h"

#include "evenlystreamlines.h"

#include "computeroadvis.h"

extern EvenStreamlinePlace *major, *minor;
extern StreetNet *streetnet ;
extern TensorLineIntersectionInfoList **majorintersectinfo ;
extern TensorLineIntersectionInfoList **minorintersectinfo ;
extern double road_width[];

StreetVis *roadnetvis = NULL;

extern int prev_nmajors ;
extern int prev_nminors ;


/*intialize the road network visualization*/
void init_roadnetvis()
{
	if(roadnetvis != NULL)
	{
		delete roadnetvis;
		roadnetvis = NULL;
	}
	roadnetvis=new StreetVis(major->evenstreamlines->ntrajs, minor->evenstreamlines->ntrajs);
	roadnetvis->init_lineseglist();
}

void reset_roadnetvis()
{
	if(roadnetvis != NULL)
	{
		roadnetvis->nmajDirRoads=0;
		roadnetvis->nminDirRoads=0;
	}

	int i;
	if(majorintersectinfo!=NULL)
	{
		for(i=0; i<prev_nmajors; i++)
			majorintersectinfo[i]->nelems=0;
	}
	if(minorintersectinfo!=NULL)
	{
		for(i=0; i<prev_nminors; i++)
			minorintersectinfo[i]->nelems=0;
	}

	if(streetnet!=NULL)
	{
		streetnet->edgelist->nedges = 0;
		streetnet->nodelist->nelems=0;
	}
}

/*
compute the four visual points of each intersection for the purpose of 
street network visualization 
*/
void compute_visual_points_intersects()
{
	int i;
	for(i=0; i<streetnet->nodelist->nelems; i++)
	{
		if(streetnet->nodelist->intersects[i]->endpt)
			break;
		compute_visualpts_oneintersect(i);
	}
}


/*
Very crucial function:
In this function, we compute the four neighboring points of each intersection.
This is crucial for finding the parcel of each block and visualizing the 
street properly.
*/
void compute_visualpts_oneintersect(int intersectid)
{
	Intersection *intersect=streetnet->nodelist->intersects[intersectid];
	Trajectory *majtraj=major->evenstreamlines->trajs[intersect->majorline_id];
	Trajectory *mintraj=minor->evenstreamlines->trajs[intersect->minorline_id];
	int maj_lineseg=intersect->majlineseg;
	int min_lineseg=intersect->minlineseg;
	icVector2 majdir, mindir;
	majdir.entry[0]=majtraj->linesegs[maj_lineseg].gend[0]-majtraj->linesegs[maj_lineseg].gstart[0];
	majdir.entry[1]=majtraj->linesegs[maj_lineseg].gend[1]-majtraj->linesegs[maj_lineseg].gstart[1];
	normalize(majdir);
	

	mindir.entry[0]=mintraj->linesegs[min_lineseg].gend[0]-mintraj->linesegs[min_lineseg].gstart[0];
	mindir.entry[1]=mintraj->linesegs[min_lineseg].gend[1]-mintraj->linesegs[min_lineseg].gstart[1];
	normalize(mindir);

	icVector2 temp;
	temp.entry[0]=-majdir.entry[1];
	temp.entry[1]=majdir.entry[0];

	if(dot(temp, mindir)<0) /*flip it*/
	{
		mindir=-mindir;
	}

	double im[2];
	im[0]=intersect->gpos[0]-(road_width[majtraj->roadtype])*mindir.entry[0];
	im[1]=intersect->gpos[1]-(road_width[majtraj->roadtype])*mindir.entry[1];
	intersect->majpt1[0]=im[0]-road_width[mintraj->roadtype]*majdir.entry[0];
	intersect->majpt1[1]=im[1]-road_width[mintraj->roadtype]*majdir.entry[1];
	intersect->majpt2[0]=im[0]+road_width[mintraj->roadtype]*majdir.entry[0];
	intersect->majpt2[1]=im[1]+road_width[mintraj->roadtype]*majdir.entry[1];
	//im[0]=intersect->gpos[0]-(road_width[majtraj->roadtype]+0.001)*mindir.entry[0];
	//im[1]=intersect->gpos[1]-(road_width[majtraj->roadtype]+0.001)*mindir.entry[1];
	//intersect->majpt1[0]=im[0]-(road_width[mintraj->roadtype]+0.002)*majdir.entry[0];
	//intersect->majpt1[1]=im[1]-(road_width[mintraj->roadtype]+0.002)*majdir.entry[1];
	//intersect->majpt2[0]=im[0]+(road_width[mintraj->roadtype]+0.002)*majdir.entry[0];
	//intersect->majpt2[1]=im[1]+(road_width[mintraj->roadtype]+0.002)*majdir.entry[1];
	
	im[0]=intersect->gpos[0]+(road_width[majtraj->roadtype])*mindir.entry[0];
	im[1]=intersect->gpos[1]+(road_width[majtraj->roadtype])*mindir.entry[1];
	intersect->majpt3[0]=im[0]-road_width[mintraj->roadtype]*majdir.entry[0];
	intersect->majpt3[1]=im[1]-road_width[mintraj->roadtype]*majdir.entry[1];
	intersect->majpt4[0]=im[0]+road_width[mintraj->roadtype]*majdir.entry[0];
	intersect->majpt4[1]=im[1]+road_width[mintraj->roadtype]*majdir.entry[1];
	//intersect->majpt3[0]=im[0]-(road_width[mintraj->roadtype]+0.002)*majdir.entry[0];
	//intersect->majpt3[1]=im[1]-(road_width[mintraj->roadtype]+0.002)*majdir.entry[1];
	//intersect->majpt4[0]=im[0]+(road_width[mintraj->roadtype]+0.002)*majdir.entry[0];
	//intersect->majpt4[1]=im[1]+(road_width[mintraj->roadtype]+0.002)*majdir.entry[1];

	/*compute the corresponding order for the minor road, should be 2,4 and 1,3*/
}


/*given the directions of two lines, compute the average of the normal.
NOTE: we assume the two lines have the same orientation*/
void StreetVis::ave_normal(double l1[2], double l2[2], double ave_n[2])
{
	double n1[2]={-l1[1], l1[0]};
	double n2[2]={-l2[1], l2[0]};

	ave_n[0]=(n1[0]+n2[0])/2.;
	ave_n[1]=(n1[1]+n2[1])/2.;
}


/*check whether there are any intersections at the specified line segment 
of the specified tensor line*/
bool contain_intersect(int linesegid, int &startpos, int &nextlineseg, 
					   TensorLineIntersectionInfoList *lineinfo, int &intersect_id)
{
	int i;
	for(i=startpos; i<lineinfo->nelems; i++)
	{
		if(lineinfo->infolist[i]->lineseg_id == linesegid)
		{
			intersect_id = lineinfo->infolist[i]->intersect_id;
			startpos = i;
			if(i<lineinfo->nelems-1)
				nextlineseg=lineinfo->infolist[i+1]->lineseg_id;
			else
				nextlineseg=-2;
			return true;
		}
	}
	//startpos = lineinfo->nelems;
	nextlineseg = -2;
	return false;
}


void StreetVis::construct_roads_vis(TrajectoryList *major, TrajectoryList *minor)
{
	int i, j;
	Trajectory *curtraj;
	TensorLineIntersectionInfoList *infolist;
	int intersect_linesegid;

	/*record the wider points for road network*/
	double pcur_s[2], ncur_s[2], ppre_s[2], npre_s[2];
	double l1[2], l2[2];
	icVector2 ln;
	int infolistpos = 0;
	int intersect_id;
	FILE *fp;
	int nextlineseg = -1;

	/*major direction*/
	for(i=0; i<major->ntrajs; i++)
	{
		curtraj=major->trajs[i];
		infolist = majorintersectinfo[i];
		infolistpos = 0;

		/*obtain the starting points*/
		l1[0]=curtraj->linesegs[0].gend[0]-curtraj->linesegs[0].gstart[0];
		l1[1]=curtraj->linesegs[0].gend[1]-curtraj->linesegs[0].gstart[1];
		ln.set(-l1[1], l1[0]);
		normalize(ln);
		/*positive width*/
		ppre_s[0]=curtraj->linesegs[0].gstart[0]+road_width[curtraj->roadtype]*ln.entry[0];
		ppre_s[1]=curtraj->linesegs[0].gstart[1]+road_width[curtraj->roadtype]*ln.entry[1];

		/*negative width*/
		npre_s[0]=curtraj->linesegs[0].gstart[0]-road_width[curtraj->roadtype]*ln.entry[0];
		npre_s[1]=curtraj->linesegs[0].gstart[1]-road_width[curtraj->roadtype]*ln.entry[1];

		nextlineseg = -1;

		for(j=1; j<curtraj->nlinesegs; j++)
		{
			/*get the direction of the current line segment*/
			l2[0]=curtraj->linesegs[j].gend[0]-curtraj->linesegs[j].gstart[0];
			l2[1]=curtraj->linesegs[j].gend[1]-curtraj->linesegs[j].gstart[1];

			/*obtain the average normal*/
			ave_normal(l1, l2, ln.entry);
			normalize(ln);

			/*positive: build a line segment for road_line_1*/
			pcur_s[0]=curtraj->linesegs[j].gstart[0]+road_width[curtraj->roadtype]*ln.entry[0];
			pcur_s[1]=curtraj->linesegs[j].gstart[1]+road_width[curtraj->roadtype]*ln.entry[1];
			RoadLineSeg *newlineseg1=(RoadLineSeg*)malloc(sizeof(RoadLineSeg));
			newlineseg1->start[0]=ppre_s[0];
			newlineseg1->start[1]=ppre_s[1];
			newlineseg1->end[0]=pcur_s[0];
			newlineseg1->end[1]=pcur_s[1];
			newlineseg1->cell_id=curtraj->linesegs[j].Triangle_ID;
			majDirRoads[i]->addNew(newlineseg1, false);


			/*negative: build a line segment for road_line_2*/
			ncur_s[0]=curtraj->linesegs[j].gstart[0]-road_width[curtraj->roadtype]*ln.entry[0];
			ncur_s[1]=curtraj->linesegs[j].gstart[1]-road_width[curtraj->roadtype]*ln.entry[1];
			RoadLineSeg *newlineseg2=(RoadLineSeg*)malloc(sizeof(RoadLineSeg));
			newlineseg2->start[0]=npre_s[0];
			newlineseg2->start[1]=npre_s[1];
			newlineseg2->end[0]=ncur_s[0];
			newlineseg2->end[1]=ncur_s[1];
			newlineseg2->cell_id=curtraj->linesegs[j].Triangle_ID;
			majDirRoads[i]->addNew(newlineseg2, true);
			
			/*save the previously obtained points*/
			ppre_s[0]=pcur_s[0];
			ppre_s[1]=pcur_s[1];
			npre_s[0]=ncur_s[0];
			npre_s[1]=ncur_s[1];

			l1[0]=l2[0];
			l1[1]=l2[1];

			/*check whether this line segment contain intersection or not*/
			if(infolistpos>=infolist->nelems-1 ) 
				continue;

			nextlineseg = -1;
L1:			if(contain_intersect(j, infolistpos, nextlineseg, infolist, intersect_id))
			{
				/*if it is end point, continue*/
				if(streetnet->nodelist->intersects[intersect_id]->endpt)
					continue;

				/*buile a new line segment and break the continuous line at 
				the position of the intersection*/
				//pcur_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt1[0];
				//pcur_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt1[1];
				pcur_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt3[0];
				pcur_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt3[1];
				RoadLineSeg *newlineseg3=(RoadLineSeg*)malloc(sizeof(RoadLineSeg));
				newlineseg3->start[0]=ppre_s[0];
				newlineseg3->start[1]=ppre_s[1];
				newlineseg3->end[0]=pcur_s[0];
				newlineseg3->end[1]=pcur_s[1];
				newlineseg3->cell_id=curtraj->linesegs[j].Triangle_ID;
				majDirRoads[i]->addNew(newlineseg3, false);

				//ncur_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt3[0];
				//ncur_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt3[1];
				ncur_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt1[0];
				ncur_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt1[1];
				RoadLineSeg *newlineseg4=(RoadLineSeg*)malloc(sizeof(RoadLineSeg));
				newlineseg4->start[0]=npre_s[0];
				newlineseg4->start[1]=npre_s[1];
				newlineseg4->end[0]=ncur_s[0];
				newlineseg4->end[1]=ncur_s[1];
				newlineseg4->cell_id=curtraj->linesegs[j].Triangle_ID;
				majDirRoads[i]->addNew(newlineseg4, true);

				/*save the previously obtained points*/
				//ppre_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt2[0];
				//ppre_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt2[1];
				//npre_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt4[0];
				//npre_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt4[1];

				ppre_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt4[0];
				ppre_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt4[1];
				npre_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt2[0];
				npre_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt2[1];

				//if(nextlineseg == j)/*have another intersection on this linesegment*/
				//	goto L1;
			}
		}
	}

	/*minor direction*/
	for(i=0; i<minor->ntrajs; i++)
	{
		curtraj=minor->trajs[i];
		infolist = minorintersectinfo[i];
		infolistpos = 0;

		/*obtain the starting points*/
		l1[0]=curtraj->linesegs[0].gend[0]-curtraj->linesegs[0].gstart[0];
		l1[1]=curtraj->linesegs[0].gend[1]-curtraj->linesegs[0].gstart[1];
		ln.set(-l1[1], l1[0]);
		normalize(ln);
		/*positive width*/
		ppre_s[0]=curtraj->linesegs[0].gstart[0]+road_width[curtraj->roadtype]*ln.entry[0];
		ppre_s[1]=curtraj->linesegs[0].gstart[1]+road_width[curtraj->roadtype]*ln.entry[1];

		/*negative width*/
		npre_s[0]=curtraj->linesegs[0].gstart[0]-road_width[curtraj->roadtype]*ln.entry[0];
		npre_s[1]=curtraj->linesegs[0].gstart[1]-road_width[curtraj->roadtype]*ln.entry[1];

		nextlineseg = -1;

		for(j=1; j<curtraj->nlinesegs; j++)
		{
			/*get the direction of the current line segment*/
			l2[0]=curtraj->linesegs[j].gend[0]-curtraj->linesegs[j].gstart[0];
			l2[1]=curtraj->linesegs[j].gend[1]-curtraj->linesegs[j].gstart[1];

			/*obtain the average normal*/
			ave_normal(l1, l2, ln.entry);
			normalize(ln);

			/*positive: build a line segment for road_line_1*/
			pcur_s[0]=curtraj->linesegs[j].gstart[0]+road_width[curtraj->roadtype]*ln.entry[0];
			pcur_s[1]=curtraj->linesegs[j].gstart[1]+road_width[curtraj->roadtype]*ln.entry[1];
			RoadLineSeg *newlineseg1=(RoadLineSeg*)malloc(sizeof(RoadLineSeg));
			newlineseg1->start[0]=ppre_s[0];
			newlineseg1->start[1]=ppre_s[1];
			newlineseg1->end[0]=pcur_s[0];
			newlineseg1->end[1]=pcur_s[1];
			newlineseg1->cell_id=curtraj->linesegs[j].Triangle_ID;
			minDirRoads[i]->addNew(newlineseg1, false);


			/*negative: build a line segment for road_line_2*/
			ncur_s[0]=curtraj->linesegs[j].gstart[0]-road_width[curtraj->roadtype]*ln.entry[0];
			ncur_s[1]=curtraj->linesegs[j].gstart[1]-road_width[curtraj->roadtype]*ln.entry[1];
			RoadLineSeg *newlineseg2=(RoadLineSeg*)malloc(sizeof(RoadLineSeg));
			newlineseg2->start[0]=npre_s[0];
			newlineseg2->start[1]=npre_s[1];
			newlineseg2->end[0]=ncur_s[0];
			newlineseg2->end[1]=ncur_s[1];
			newlineseg2->cell_id=curtraj->linesegs[j].Triangle_ID;
			minDirRoads[i]->addNew(newlineseg2, true);
			
			/*save the previously obtained points*/
			ppre_s[0]=pcur_s[0];
			ppre_s[1]=pcur_s[1];
			npre_s[0]=ncur_s[0];
			npre_s[1]=ncur_s[1];
			
			l1[0]=l2[0];
			l1[1]=l2[1];

			/*check whether this line segment contain intersection or not*/
			nextlineseg = -1;
			
			if(infolistpos>=infolist->nelems-1) 
				continue;
			
L2:			if(contain_intersect(j, infolistpos, nextlineseg, infolist, intersect_id))
			{
				/*if it is end point, continue*/
				if(streetnet->nodelist->intersects[intersect_id]->endpt)
					continue;

				Intersection *intersect=streetnet->nodelist->intersects[intersect_id];
				Trajectory *majtraj=major->trajs[intersect->majorline_id];
				Trajectory *mintraj=minor->trajs[intersect->minorline_id];
				int maj_lineseg=intersect->majlineseg;
				int min_lineseg=intersect->minlineseg;
				icVector2 temp_maj, temp_min;
				temp_maj.entry[0]=majtraj->linesegs[maj_lineseg].gend[0]
					-majtraj->linesegs[maj_lineseg].gstart[0];
				temp_maj.entry[1]=majtraj->linesegs[maj_lineseg].gend[1]
					-majtraj->linesegs[maj_lineseg].gstart[1];
				normalize(temp_maj);
				

				temp_min.entry[0]=mintraj->linesegs[min_lineseg].gend[0]
					-mintraj->linesegs[min_lineseg].gstart[0];
				temp_min.entry[1]=mintraj->linesegs[min_lineseg].gend[1]
					-mintraj->linesegs[min_lineseg].gstart[1];
				normalize(temp_min);

				icVector2 temp;
				temp.entry[0]=-temp_maj.entry[1];
				temp.entry[1]=temp_maj.entry[0];

				double re=dot(temp, temp_min);

				/*buile a new line segment and break the continuous line at 
				the position of the intersection*/
				if(re>=0)
				{
					pcur_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt1[0];
					pcur_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt1[1];
					RoadLineSeg *newlineseg3=(RoadLineSeg*)malloc(sizeof(RoadLineSeg));
					newlineseg3->start[0]=ppre_s[0];
					newlineseg3->start[1]=ppre_s[1];
					newlineseg3->end[0]=pcur_s[0];
					newlineseg3->end[1]=pcur_s[1];
					newlineseg3->cell_id=curtraj->linesegs[j].Triangle_ID;
					minDirRoads[i]->addNew(newlineseg3, false);

					ncur_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt2[0];
					ncur_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt2[1];
					RoadLineSeg *newlineseg4=(RoadLineSeg*)malloc(sizeof(RoadLineSeg));
					newlineseg4->start[0]=npre_s[0];
					newlineseg4->start[1]=npre_s[1];
					newlineseg4->end[0]=ncur_s[0];
					newlineseg4->end[1]=ncur_s[1];
					newlineseg4->cell_id=curtraj->linesegs[j].Triangle_ID;
					minDirRoads[i]->addNew(newlineseg4, true);

					/*save the previously obtained points*/
					ppre_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt3[0];
					ppre_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt3[1];
					npre_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt4[0];
					npre_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt4[1];
				}
				else
				{
					pcur_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt4[0];
					pcur_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt4[1];
					RoadLineSeg *newlineseg3=(RoadLineSeg*)malloc(sizeof(RoadLineSeg));
					newlineseg3->start[0]=ppre_s[0];
					newlineseg3->start[1]=ppre_s[1];
					newlineseg3->end[0]=pcur_s[0];
					newlineseg3->end[1]=pcur_s[1];
					newlineseg3->cell_id=curtraj->linesegs[j].Triangle_ID;
					minDirRoads[i]->addNew(newlineseg3, false);

					ncur_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt3[0];
					ncur_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt3[1];
					RoadLineSeg *newlineseg4=(RoadLineSeg*)malloc(sizeof(RoadLineSeg));
					newlineseg4->start[0]=npre_s[0];
					newlineseg4->start[1]=npre_s[1];
					newlineseg4->end[0]=ncur_s[0];
					newlineseg4->end[1]=ncur_s[1];
					newlineseg4->cell_id=curtraj->linesegs[j].Triangle_ID;
					minDirRoads[i]->addNew(newlineseg4, true);

					/*save the previously obtained points*/
					ppre_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt2[0];
					ppre_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt2[1];
					npre_s[0]=streetnet->nodelist->intersects[intersect_id]->majpt1[0];
					npre_s[1]=streetnet->nodelist->intersects[intersect_id]->majpt1[1];
				}

				//if(nextlineseg == j)/*have another intersection on this linesegment*/
				//	goto L2;
			}
		}
	}
}


void StreetVis::init_lineseglist()
{
	int i;
	for(i=0; i<nmajDirRoads; i++)
	{
		majDirRoads[i]=new OneRoadVis(major->evenstreamlines->trajs[i]->nlinesegs);
	}
	for(i=0; i<nminDirRoads; i++)
	{
		minDirRoads[i]=new OneRoadVis(minor->evenstreamlines->trajs[i]->nlinesegs);
	}
}
