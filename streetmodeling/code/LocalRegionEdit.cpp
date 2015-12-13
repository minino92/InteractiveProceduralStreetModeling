
/*  LocalRegionEdit.cpp */


#include "stdafx.h"
#include "RegionSmoothing.h"
#include "VFDataStructure.h"
#include "LocalTracing.h"
#include "tensoranalysis.h"
#include "evenlystreamlines.h"

#include "LocalRegionEdit.h"
#include "regionsmooth_quad.h"
#include "SketchDesign.h"

#include "Numerical.h"

extern int stack[200000];              /////stack for flood search inside a closed region
extern int stack_top;
extern QuadMesh *quadmesh;

extern icVector2 tenline_dir_global;

extern int *region_quadverts;                ////mesh vertices inside user selected region
extern int nregion_quadverts;
extern int curMaxRegionQuadVerts;

extern Point *point;                       ////we may initial it as 50 points, over 50, we can extend it
extern int Num_SmoothRegionpoints;                     ////Number of points that user selected
extern int firstVertID_quad;


extern int *boundarycells;
extern int nboundarycells;

extern StreetNet *streetnet;

int *innercells=NULL;
int ninnercells=0;

int *innerintersections=NULL;
int ninnerintersections=0;

/*  record the intersections on the boundary  */
int *boundaryintersections=NULL;
int nboundaryintersections=0;

Trajectory *regionboundary=NULL;
extern unsigned char cur_max_reg_index;  /*   current maximum region index   */

extern EvenStreamlinePlace *major, *minor;
extern TensorLineIntersectionInfoList **majorintersectinfo;
extern TensorLineIntersectionInfoList **minorintersectinfo;

//int *boundaryintersects=NULL;
//int nboundaryintersects=0;
extern int prev_nmajors ;
extern int prev_nminors ;

int *contour_cells=NULL;
int ncontour_cells=0;


extern bool is_inregion(double x, double y);
extern int cal_intersect(double PointA[2], double PointB[2], 
				  double PointC[2], double PointD[2], double t[2]);
extern void sample_along_tensorline_from_to(Trajectory *traj, int start_lineseg, double startp[2],
									 int end_lineseg, double endp[2], 
									 int edgeid, StreetNet *net);
extern int *extend_link(int *edge_link, int Num_edges);

bool is_on_local_editing=false;

/*   After obtaining the inner vertices, we need to find out the inner cells */
void find_inner_cells()
{
	int i, j;
	QuadVertex *v;

	if(innercells==NULL)
		innercells=(int*)malloc(sizeof(int)*quadmesh->nfaces);
	ninnercells=0;

	for(i=0;i<quadmesh->nfaces;i++)
		quadmesh->quadcells[i]->visited=false;

	for(i=0;i<nregion_quadverts;i++)
	{
		v=quadmesh->quad_verts[region_quadverts[i]];
		for(j=0;j<v->ncells;j++)
		{
			if(v->cells[j]->visited||v->cells[j]->OnBoundary)
				continue;

			/*   add to the inner cell list   */
			innercells[ninnercells]=v->cells[j]->index;
			ninnercells++;

			/*   set it as visited   */
			v->cells[j]->visited=true;
		}
	}
}



/*  After finding the inner cells, we can locate the intersections that fall in the 
    user specified region.
	We now need to locate the *boundary* intersections, they are the intersections
	that have at least one edge connecting with the intersection outside the region
*/

void find_inner_intersections()
{
	int i, j;
	QuadCell *face;

	/*  initialization   */
	if(innerintersections!=NULL)
	{
		free(innerintersections);
		innerintersections=NULL;
	}
	innerintersections=(int*)malloc(sizeof(int)*streetnet->nodelist->nelems);
	ninnerintersections=0;

	for(i=0;i<streetnet->nodelist->nelems;i++)
	{
		streetnet->nodelist->intersects[i]->inside_region=false;
	}

	/*we first deal with the inner cells, since all the intersections of 
	the inner cells should be inner intersections
	*/
	for(i=0;i<ninnercells;i++)
	{
		face=quadmesh->quadcells[innercells[i]];

		for(j=0;j<face->nintersects;j++)
		{
			innerintersections[ninnerintersections]=face->intersectlist[j];
			ninnerintersections++;

			streetnet->nodelist->intersects[face->intersectlist[j]]->inside_region=true;
		}
	}

	/*second, we deal with the boundary cells, because the intersections in these
	cell may not be inner intersections
	*/
	for(i=0;i<nboundarycells;i++)
	{
		face=quadmesh->quadcells[boundarycells[i]];

		for(j=0;j<face->nintersects;j++)
		{
			/*  need to judge whether the intersection is inside the region  */
			Intersection *intersect=streetnet->nodelist->intersects[face->intersectlist[j]];
			if(is_inregion(intersect->gpos[0],intersect->gpos[1]))
			{
				innerintersections[ninnerintersections]=face->intersectlist[j];
				ninnerintersections++;
				intersect->inside_region=true;
			}
		}
	}
}


/*
    Obtain the length of the sample list associated with the specified edge 
	of the given network (graph)
*/

double get_samplength_given_edge(int edge_id, StreetNet *net)
{
	int i;
	StreetGraphEdge *edge=net->edgelist->edges[edge_id];

	double len=0;

	if(edge->ninter_pts<=1) return len;

	icVector2 dist;
	for(i=0;i<edge->ninter_pts-1;i++)
	{
		dist.entry[0]=edge->inter_pts[i+1]->x-edge->inter_pts[i]->x;
		dist.entry[1]=edge->inter_pts[i+1]->y-edge->inter_pts[i]->y;

		len+=length(dist);
	}

	return len;
}



/*
    After determining the inner intersections, we now can compute the new intersections
	on the region boundary and rebuild the graph 
*/

void recal_intersections_on_RegBoundary()
{
	int i, j;
	Intersection *intersect, *other_intersect;
	StreetGraphEdge *cur_edge;
	double intersect_p[2];
	int which_cell;

	/* initialize edge list */
	for(i=0;i<streetnet->edgelist->nedges;i++)
		streetnet->edgelist->edges[i]->visited=false;

	if(boundaryintersections!=NULL)
		free(boundaryintersections);
	boundaryintersections=(int*)malloc(sizeof(int)*streetnet->nodelist->nelems);
	nboundaryintersections=0;

	for(i=0;i<ninnerintersections;i++)
	{
		intersect=streetnet->nodelist->intersects[innerintersections[i]];

		for(j=0;j<intersect->nadjedges;j++)
		{
			cur_edge=streetnet->edgelist->edges[intersect->adj_edges[j]];
			cur_edge->visited=true;

			if(cur_edge->node_index1!=innerintersections[i])
				other_intersect=streetnet->nodelist->intersects[cur_edge->node_index1];
			else
				other_intersect=streetnet->nodelist->intersects[cur_edge->node_index2];

			if(other_intersect->inside_region)
			{
				cur_edge->cancel=true;  // mark this edge as "cancelled" !
				continue;
			}

			/*  we now need to search for the new intersection along current edge   */
			search_for_intersection_along_edge(cur_edge->index, intersect_p, which_cell);

			/*  we make use of the line segment list associated with this edge to 
			    compute the new intersection, we then update this edge with new intersection
			*/

			Intersection *newintersect=(Intersection*)malloc(sizeof(Intersection));
			newintersect->gpos[0]=intersect_p[0];
			newintersect->gpos[1]=intersect_p[1];
			newintersect->cellid=which_cell;
			newintersect->majorline_id=-1;   /*   we need to figure out how to update these   */
			newintersect->minorline_id=-1;
			newintersect->majlineseg=0;
			newintersect->nadjedges=0;
			newintersect->adj_edges=NULL;
			newintersect->endpt=true;
			newintersect->inside_region=false;
			newintersect->deleted=false;

			streetnet->nodelist->addNew(newintersect);

			add_to_cell_intersectlist(newintersect->cellid, newintersect->index, false);

			/*   update its edge list   */

			if(cur_edge->node_index1==innerintersections[i])
			{
				cur_edge->node_index1 = newintersect->index;
				
				/*  Update the sample list associated with this edge  */
				int k, pos=-1;
				for(k=cur_edge->ninter_pts-1; k>=0; k--)
				{
					if(cur_edge->inter_pts[k]->cellid==newintersect->cellid)
					{
						cur_edge->inter_pts[k]->x=intersect_p[0];
						cur_edge->inter_pts[k]->y=intersect_p[1];
						pos=k;

						break;
					}
				}
				/* release the others */
				if(pos>=0)
				{
					for(k=0;k<pos;k++)
					{
						free(cur_edge->inter_pts[k]);
						cur_edge->inter_pts[k]=NULL;
					}

					for(k=pos;k<cur_edge->ninter_pts;k++)
					{
						cur_edge->inter_pts[k-pos]=cur_edge->inter_pts[k];
					}
					cur_edge->ninter_pts-=pos;
				}
			}
			else
			{
				cur_edge->node_index2 = newintersect->index;

				int k, pos=-1;

				/*  Update the sample list associated with this edge  */
				for(k=0; k<cur_edge->ninter_pts; k++)
				{
					if(cur_edge->inter_pts[k]->cellid==newintersect->cellid)
					{
						cur_edge->inter_pts[k]->x=intersect_p[0];
						cur_edge->inter_pts[k]->y=intersect_p[1];
						pos=k;

						break;
					}
				}

				/* release the others */
				if(pos>=0 && pos<cur_edge->ninter_pts)
				{
					for(k=pos+1;k<cur_edge->ninter_pts;k++)
					{
						free(cur_edge->inter_pts[k]);
						cur_edge->inter_pts[k]=NULL;
					}
					cur_edge->ninter_pts=pos+1;
				}
			}

			/*  obtain the length  of the new edge  */
			if(get_samplength_given_edge(cur_edge->index, streetnet)<quadmesh->xinterval/10.)
			{
				cur_edge->cancel=true;
				newintersect->deleted=true;
			}
			else{
				newintersect->adj_edges=(int*)malloc(sizeof(int));
				newintersect->adj_edges[0]=cur_edge->index;
				newintersect->nadjedges=1;

				/*  record the new intersection  */
				boundaryintersections[nboundaryintersections]=newintersect->index;
				nboundaryintersections++;
			}
		}


		/*  After finishing dealing with this intersection, we need to remove this intersection
		    and update the street network graph as well (hard)
		*/

		intersect->deleted=true;
	}
}


/*
   Compute the distance to the boundary of the region.
   If it is really small, then we say it is "inside" the region 
*/
bool is_inregion_appro_dis(double x, double y, double threshold)
{
	int i;
	double distSeg, distLine;

	////Calculate the sum of the angle
	for( i = 0; i < Num_SmoothRegionpoints; i++)
	{

		DistanceFromLine(x,y, point[i].x, point[i].y, 
			point[(i+1)%Num_SmoothRegionpoints].x, point[(i+1)%Num_SmoothRegionpoints].y,
			distSeg, distLine);

		if(distSeg<threshold)
			return true;

	}

	return false;
}


/*
     A new approach to compute the intersections between the graph edges and the 
	 boundary of the user specified region
*/
void recal_intersections_on_RegBoundary_2()
{
	int i, j;
	Intersection *intersect, *other_intersect;
	StreetGraphEdge *cur_edge;
	double intersect_p[2];
	int which_cell;

	/* initialize edge list */
	for(i=0;i<streetnet->edgelist->nedges;i++)
		streetnet->edgelist->edges[i]->visited=false;

	if(boundaryintersections!=NULL)
		free(boundaryintersections);
	boundaryintersections=(int*)malloc(sizeof(int)*streetnet->nodelist->nelems);
	nboundaryintersections=0;

	/*   mark all the inner edges of the street graph that are completely  
	     inside the user specified region 
		 NOTE: we are not considering the curved edges yet by using this 
		 naive method
	*/
	for(i=0;i<ninnerintersections;i++)
	{
		intersect=streetnet->nodelist->intersects[innerintersections[i]];

		for(j=0;j<intersect->nadjedges;j++)
		{
			cur_edge=streetnet->edgelist->edges[intersect->adj_edges[j]];

			if(cur_edge->visited) continue;

			cur_edge->visited=true;

			if(cur_edge->node_index1!=innerintersections[i])
				other_intersect=streetnet->nodelist->intersects[cur_edge->node_index1];
			else
				other_intersect=streetnet->nodelist->intersects[cur_edge->node_index2];

			if(other_intersect->inside_region)
			{
				/*  we use the middle point to handle the second order(quadratic) curved edge  */
				if(cur_edge->ninter_pts>2)
				{
					/*  choose the middle point of the sample points associated with the edge  */
					double middle[2];
					middle[0]=cur_edge->inter_pts[(int)(cur_edge->ninter_pts/2)]->x;
					middle[1]=cur_edge->inter_pts[(int)(cur_edge->ninter_pts/2)]->y;

					if(is_inregion(middle[0], middle[1]))
					{
						cur_edge->cancel=true;
					}
				}
				else
					cur_edge->cancel=true;  // mark this edge as "cancelled" !
			}
			else
			{
				/*  if other_intersect is also really close (need a threshold) to the boundary  */
				if(is_inregion_appro(other_intersect->gpos[0], other_intersect->gpos[1],
						0.05)
					||is_inregion_appro_dis(other_intersect->gpos[0], other_intersect->gpos[1],
						5.e-4))
				{
					cur_edge->cancel=true;
					boundaryintersections[nboundaryintersections]=other_intersect->index;
					nboundaryintersections++;
				}
			}
		}
		intersect->deleted=true;  /*  delete this inner intersection  */
	}

	/*
	    We now search along the converted boundary (i.e. a trajectory)
	*/
	int center_cell;
	QuadCell *face;
	int which_samp;
	double A[2], B[2];

	int dual;

	for(i=0; i<regionboundary->nlinesegs; i++)
	{
		center_cell=regionboundary->linesegs[i].Triangle_ID;
		A[0]=regionboundary->linesegs[i].gstart[0];
		A[1]=regionboundary->linesegs[i].gstart[1];
		B[0]=regionboundary->linesegs[i].gend[0];
		B[1]=regionboundary->linesegs[i].gend[1];

		for(dual=0;dual<2;dual++)
		{
			if(dual==0)
		center_cell=get_cellID_givencoords(A[0]+1.e-7, A[1]+1.e-7);
			else
		center_cell=get_cellID_givencoords(A[0]-1.e-7, A[1]-1.e-7);

		/*  method 1:  consider only this cell  */
		face=quadmesh->quadcells[center_cell];
		for(j=0;j<face->nstreetgraphedges;j++)
		{
			cur_edge=streetnet->edgelist->edges[face->streetgraphedgelist[j]];
			/*  search intersection along the edge */
			if(cal_intersect_at_graph_edge(face->streetgraphedgelist[j],
				intersect_p, which_samp, A, B))
			{
				/*   create a new intersection here   */
				Intersection *newintersect=(Intersection*)malloc(sizeof(Intersection));
				newintersect->gpos[0]=intersect_p[0];
				newintersect->gpos[1]=intersect_p[1];
				newintersect->cellid=get_cellID_givencoords(intersect_p[0], intersect_p[1]);
				newintersect->majorline_id=-1;   /*   we need to figure out how to update these   */
				newintersect->minorline_id=-1;
				newintersect->majlineseg=0;
				newintersect->nadjedges=0;
				newintersect->adj_edges=NULL;
				newintersect->endpt=true;
				newintersect->inside_region=false;
				newintersect->deleted=false;

				streetnet->nodelist->addNew(newintersect);

				add_to_cell_intersectlist(newintersect->cellid, newintersect->index, false);

				/*   update its edge list   */

				if(streetnet->nodelist->intersects[cur_edge->node_index1]->inside_region)
				{
					cur_edge->node_index1 = newintersect->index;
					
					/*  Update the sample list associated with this edge  */
					int k, pos=-1;
					for(k=cur_edge->ninter_pts-1; k>=0; k--)
					{
						if(cur_edge->inter_pts[k]->cellid==newintersect->cellid)
						{
							cur_edge->inter_pts[k]->x=intersect_p[0];
							cur_edge->inter_pts[k]->y=intersect_p[1];
							pos=k;

							break;
						}
					}
					/* release the others */
					if(pos>=0)
					{
						if(pos==0)
						{
							cur_edge->cancel=true;
							newintersect->deleted=true;
						}
						else
						{
							for(k=0;k<pos;k++)
							{
								free(cur_edge->inter_pts[k]);
								cur_edge->inter_pts[k]=NULL;
							}

							for(k=pos;k<cur_edge->ninter_pts;k++)
							{
								cur_edge->inter_pts[k-pos]=cur_edge->inter_pts[k];
							}
							cur_edge->ninter_pts-=pos;
						}
					}

					/*  new method to update the sample list of this edge
					    could contain bug  1/17/2008
					*/
					//double start[2], end[2];
					//int start_cell, end_cell;
					//start[0]=newintersect->gpos[0];
					//start[1]=newintersect->gpos[1];
					//start_cell=get_cellID_givencoords(start[0], start[1]);
					//end[0]=cur_edge->inter_pts[cur_edge->ninter_pts-1]->x;
					//end[1]=cur_edge->inter_pts[cur_edge->ninter_pts-1]->y;
					//end_cell=get_cellID_givencoords(end[0], end[1]);
					//Trajectory *tempTraj=new Trajectory(-1);
					//get_linesegs_anytwopts(start,start_cell, end,end_cell, tempTraj, 0, 5);
					//free(cur_edge->inter_pts);
					//cur_edge->inter_pts=NULL;
					//cur_edge->inter_pts=(Point **)malloc(sizeof(Point *)* (tempTraj->nlinesegs+1));

					//int m;
					//for(m=0;m<tempTraj->nlinesegs;m++)
					//{
					//	cur_edge->inter_pts[m]=(Point*)malloc(sizeof(Point));
					//	cur_edge->inter_pts[m]->x=tempTraj->linesegs[m].gstart[0];
					//	cur_edge->inter_pts[m]->y=tempTraj->linesegs[m].gstart[1];
					//	cur_edge->inter_pts[m]->cellid=get_cellID_givencoords(cur_edge->inter_pts[m]->x,
					//		cur_edge->inter_pts[m]->y);
					//	
					//	//add_to_edgelist_one_cell(cur_edge->inter_pts[j]->cellid, i);
					//}

					//cur_edge->inter_pts[m]=(Point*)malloc(sizeof(Point));
					//cur_edge->inter_pts[m]->x=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[0];
					//cur_edge->inter_pts[m]->y=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[1];
					//cur_edge->inter_pts[m]->cellid=get_cellID_givencoords(cur_edge->inter_pts[m]->x,
					//	cur_edge->inter_pts[m]->y);
					//cur_edge->ninter_pts=tempTraj->nlinesegs+1;

					//delete tempTraj;
				}
				else
				{
					cur_edge->node_index2 = newintersect->index;

					int k, pos=-1;

					/*  Update the sample list associated with this edge  */
					for(k=0; k<cur_edge->ninter_pts; k++)
					{
						if(cur_edge->inter_pts[k]->cellid==newintersect->cellid)
						{
							cur_edge->inter_pts[k]->x=intersect_p[0];
							cur_edge->inter_pts[k]->y=intersect_p[1];
							pos=k;

							break;
						}
					}

					/* release the others */
					if(pos>=0 && pos<cur_edge->ninter_pts)
					{
						if(pos==0)
						{
							cur_edge->cancel=true;
							newintersect->deleted=true;
						}
						else
						{
							for(k=pos+1;k<cur_edge->ninter_pts;k++)
							{
								free(cur_edge->inter_pts[k]);
								cur_edge->inter_pts[k]=NULL;
							}
							cur_edge->ninter_pts=pos+1;
						}
					}

					/*  new method to update the sample list of this edge
					    could contain bug  1/17/2008
					*/
					//double start[2], end[2];
					//int start_cell, end_cell;
					//end[0]=newintersect->gpos[0];
					//end[1]=newintersect->gpos[1];
					//end_cell=get_cellID_givencoords(end[0], end[1]);
					//start[0]=cur_edge->inter_pts[cur_edge->ninter_pts-1]->x;
					//start[1]=cur_edge->inter_pts[cur_edge->ninter_pts-1]->y;
					//start_cell=get_cellID_givencoords(start[0], start[1]);
					//Trajectory *tempTraj=new Trajectory(-1);
					//get_linesegs_anytwopts(start,start_cell, end,end_cell, tempTraj, 0, 5);
					//free(cur_edge->inter_pts);
					//cur_edge->inter_pts=NULL;
					//cur_edge->inter_pts=(Point **)malloc(sizeof(Point *)* (tempTraj->nlinesegs+1));

					//int m;
					//for(m=0;m<tempTraj->nlinesegs;m++)
					//{
					//	cur_edge->inter_pts[m]=(Point*)malloc(sizeof(Point));
					//	cur_edge->inter_pts[m]->x=tempTraj->linesegs[m].gstart[0];
					//	cur_edge->inter_pts[m]->y=tempTraj->linesegs[m].gstart[1];
					//	cur_edge->inter_pts[m]->cellid=get_cellID_givencoords(cur_edge->inter_pts[m]->x,
					//		cur_edge->inter_pts[m]->y);
					//	
					//	//add_to_edgelist_one_cell(cur_edge->inter_pts[j]->cellid, i);
					//}

					//cur_edge->inter_pts[m]=(Point*)malloc(sizeof(Point));
					//cur_edge->inter_pts[m]->x=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[0];
					//cur_edge->inter_pts[m]->y=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[1];
					//cur_edge->inter_pts[m]->cellid=get_cellID_givencoords(cur_edge->inter_pts[m]->x,
					//	cur_edge->inter_pts[m]->y);
					//cur_edge->ninter_pts=tempTraj->nlinesegs+1;

					//delete tempTraj;
				}

				/*  obtain the length of the new edge  */
				if(get_samplength_given_edge(cur_edge->index, streetnet)<quadmesh->xinterval/2.)
				{
					cur_edge->cancel=true;
					newintersect->deleted=true;
				}
				else{
					newintersect->adj_edges=(int*)malloc(sizeof(int));
					newintersect->adj_edges[0]=cur_edge->index;
					newintersect->nadjedges=1;

					/*  record the new intersection  */
					boundaryintersections[nboundaryintersections]=newintersect->index;
					nboundaryintersections++;
				}
			}

			else
			{
				/*  Can we use some approximation here?  */

				/*  How about we get rid of this edge completely  */
				//cur_edge->cancel=true;
			}
		}

		}
		/*  method 2:  consider a set of cells with this cell as the center one  */
		/* if one of the end point of the line segment is (close to) a vertex of the mesh 
		   we need to use this method!
		*/
	}

	/*  search the edge list and mark those artifact edges  */
	for(i=0;i<streetnet->edgelist->nedges;i++)
	{
		cur_edge=streetnet->edgelist->edges[i];

		if(cur_edge->cancel)
		{
			if(streetnet->nodelist->intersects[cur_edge->node_index1]->endpt)
				streetnet->nodelist->intersects[cur_edge->node_index1]->deleted=true;
			if(streetnet->nodelist->intersects[cur_edge->node_index2]->endpt)
				streetnet->nodelist->intersects[cur_edge->node_index2]->deleted=true;
			continue;
		}

		if(streetnet->nodelist->intersects[cur_edge->node_index1]->deleted
			|| streetnet->nodelist->intersects[cur_edge->node_index2]->deleted)
		{
			cur_edge->cancel=true;
			if(streetnet->nodelist->intersects[cur_edge->node_index1]->endpt)
				streetnet->nodelist->intersects[cur_edge->node_index1]->deleted=true;
			if(streetnet->nodelist->intersects[cur_edge->node_index2]->endpt)
				streetnet->nodelist->intersects[cur_edge->node_index2]->deleted=true;
		}

		/*  remove dangling edges  */
		if(streetnet->nodelist->intersects[cur_edge->node_index1]->endpt
			&&streetnet->nodelist->intersects[cur_edge->node_index1]->endpt)
			cur_edge->cancel=true;
	}

	/*  remove artifact nodes  */
	for(i=0;i<streetnet->nodelist->nelems;i++)
	{
		intersect=streetnet->nodelist->intersects[i];

		if(intersect->deleted)
			continue;

		//if(intersect->nadjedges==1)
		//{
		//	if(streetnet->edgelist->edges[intersect->adj_edges[0]]->cancel)
		//		intersect->deleted=true;
		//}

		int nondeleted=0;
		for(j=0;j<intersect->nadjedges;j++)
		{
			if(streetnet->edgelist->edges[intersect->adj_edges[j]]->cancel)
				continue;
			nondeleted++;
		}

		if(nondeleted==0)
		{
			intersect->deleted=true;
		}
	}


	/*   update the whole street network here   */
	update_street_network();
}





void recal_intersections_on_RegBoundary_3()
{
	int i, j;
	Intersection *intersect, *other_intersect;
	StreetGraphEdge *cur_edge;
	double intersect_p[2];
	int which_cell;

	/* initialize edge list */
	for(i=0;i<streetnet->edgelist->nedges;i++)
		streetnet->edgelist->edges[i]->visited=false;

	if(boundaryintersections!=NULL)
		free(boundaryintersections);
	boundaryintersections=(int*)malloc(sizeof(int)*streetnet->nodelist->nelems);
	nboundaryintersections=0;

	/*   mark all the inner edges of the street graph that are completely  
	     inside the user specified region 
		 NOTE: we are not considering the curved edges yet by using this 
		 naive method
	*/
	for(i=0;i<ninnerintersections;i++)
	{
		intersect=streetnet->nodelist->intersects[innerintersections[i]];

		for(j=0;j<intersect->nadjedges;j++)
		{
			cur_edge=streetnet->edgelist->edges[intersect->adj_edges[j]];
			cur_edge->cancel=true;
			if(cur_edge->node_index1!=innerintersections[i])
				other_intersect=streetnet->nodelist->intersects[cur_edge->node_index1];
			else
				other_intersect=streetnet->nodelist->intersects[cur_edge->node_index2];

			/*  if other_intersect is also really close (need a threshold) to the boundary  */
			if(is_inregion_appro(other_intersect->gpos[0], other_intersect->gpos[1],
					0.05)
				||is_inregion_appro_dis(other_intersect->gpos[0], other_intersect->gpos[1],
					5.e-4))
			{
				boundaryintersections[nboundaryintersections]=other_intersect->index;
				nboundaryintersections++;
			}
		}
		intersect->deleted=true;  /*  delete this inner intersection  */
	}

	/*
	    We now search along the converted boundary (i.e. a trajectory)
	*/
	int center_cell;
	QuadCell *face;
	int which_samp;
	double A[2], B[2];

	int dual;

	for(i=0; i<regionboundary->nlinesegs; i++)
	{
		center_cell=regionboundary->linesegs[i].Triangle_ID;
		A[0]=regionboundary->linesegs[i].gstart[0];
		A[1]=regionboundary->linesegs[i].gstart[1];
		B[0]=regionboundary->linesegs[i].gend[0];
		B[1]=regionboundary->linesegs[i].gend[1];

		for(dual=0;dual<2;dual++)
		{
			if(dual==0)
		center_cell=get_cellID_givencoords(A[0]+1.e-7, A[1]+1.e-7);
			else
		center_cell=get_cellID_givencoords(A[0]-1.e-7, A[1]-1.e-7);

		/*  method 1:  consider only this cell  */
		face=quadmesh->quadcells[center_cell];
		for(j=0;j<face->nstreetgraphedges;j++)
		{
			cur_edge=streetnet->edgelist->edges[face->streetgraphedgelist[j]];
			/*  search intersection along the edge */
			if(cal_intersect_at_graph_edge(face->streetgraphedgelist[j],
				intersect_p, which_samp, A, B))
			{
				/*   create a new intersection here   */
				Intersection *newintersect=(Intersection*)malloc(sizeof(Intersection));
				newintersect->gpos[0]=intersect_p[0];
				newintersect->gpos[1]=intersect_p[1];
				newintersect->cellid=get_cellID_givencoords(intersect_p[0], intersect_p[1]);
				newintersect->majorline_id=-1;   /*   we need to figure out how to update these   */
				newintersect->minorline_id=-1;
				newintersect->majlineseg=0;
				newintersect->nadjedges=0;
				newintersect->adj_edges=NULL;
				newintersect->endpt=true;
				newintersect->inside_region=false;
				newintersect->deleted=false;

				streetnet->nodelist->addNew(newintersect);

				add_to_cell_intersectlist(newintersect->cellid, newintersect->index, false);

				/*   update its edge list   */

				/*  new method: possible bug 1/19/2008
				    we create a new edge in the global edge list
					and cancel the old one!
				*/

				cur_edge->cancel=true;
				StreetGraphEdge *onenewedge=(StreetGraphEdge *)
					malloc(sizeof(StreetGraphEdge));
				onenewedge->cancel=onenewedge->visited=false;
				onenewedge->node_index1=newintersect->index;
				if(streetnet->nodelist->intersects[cur_edge->node_index1]->inside_region)
					onenewedge->node_index2=cur_edge->node_index2;
				else
					onenewedge->node_index2=cur_edge->node_index1;
				onenewedge->roadtype = 
					streetnet->edgelist->edges[face->streetgraphedgelist[j]]->roadtype;

				streetnet->edgelist->append(onenewedge);

				double start[2], end[2];
				int start_cell, end_cell;
				start[0]=newintersect->gpos[0];
				start[1]=newintersect->gpos[1];
				start_cell=get_cellID_givencoords(start[0], start[1]);
				end[0]=streetnet->nodelist->intersects[onenewedge->node_index2]->gpos[0];
				end[1]=streetnet->nodelist->intersects[onenewedge->node_index2]->gpos[1];
				end_cell=get_cellID_givencoords(end[0], end[1]);
				Trajectory *tempTraj=new Trajectory(-1);
				get_linesegs_anytwopts(start,start_cell, end,end_cell, tempTraj, 0, 10);
				onenewedge->inter_pts=(Point **)malloc(sizeof(Point *)* (tempTraj->nlinesegs+1));

				int m;
				for(m=0;m<tempTraj->nlinesegs;m++)
				{
					onenewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
					onenewedge->inter_pts[m]->x=tempTraj->linesegs[m].gstart[0];
					onenewedge->inter_pts[m]->y=tempTraj->linesegs[m].gstart[1];
					onenewedge->inter_pts[m]->cellid=get_cellID_givencoords(cur_edge->inter_pts[m]->x,
						cur_edge->inter_pts[m]->y);
					
					//add_to_edgelist_one_cell(cur_edge->inter_pts[j]->cellid, i);
				}

				onenewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
				onenewedge->inter_pts[m]->x=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[0];
				onenewedge->inter_pts[m]->y=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[1];
				onenewedge->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge->inter_pts[m]->x,
					onenewedge->inter_pts[m]->y);
				onenewedge->ninter_pts=tempTraj->nlinesegs+1;
				delete tempTraj;
				if(get_samplength_given_edge(onenewedge->index, streetnet)<quadmesh->xinterval/2.)
				{
					onenewedge->cancel=true;
					newintersect->deleted=true;
				}
				else{
					newintersect->adj_edges=(int*)malloc(sizeof(int));
					newintersect->adj_edges[0]=onenewedge->index;
					newintersect->nadjedges=1;

					/*  add to the other intersection  */
					add_edge_to_intersectnode(onenewedge->node_index2, onenewedge->index);

					/*  record the new intersection  */
					boundaryintersections[nboundaryintersections]=newintersect->index;
					nboundaryintersections++;
				}

			}
		}

		}
		/*  method 2:  consider a set of cells with this cell as the center one  */
		/* if one of the end point of the line segment is (close to) a vertex of the mesh 
		   we need to use this method!
		*/
	}

	/*  search the edge list and mark those artifact edges  */
	//for(i=0;i<streetnet->edgelist->nedges;i++)
	//{
	//	cur_edge=streetnet->edgelist->edges[i];

	//	if(cur_edge->cancel)
	//	{
	//		if(streetnet->nodelist->intersects[cur_edge->node_index1]->endpt)
	//			streetnet->nodelist->intersects[cur_edge->node_index1]->deleted=true;
	//		if(streetnet->nodelist->intersects[cur_edge->node_index2]->endpt)
	//			streetnet->nodelist->intersects[cur_edge->node_index2]->deleted=true;
	//		continue;
	//	}

	//	if(streetnet->nodelist->intersects[cur_edge->node_index1]->deleted
	//		|| streetnet->nodelist->intersects[cur_edge->node_index2]->deleted)
	//	{
	//		cur_edge->cancel=true;
	//		if(streetnet->nodelist->intersects[cur_edge->node_index1]->endpt)
	//			streetnet->nodelist->intersects[cur_edge->node_index1]->deleted=true;
	//		if(streetnet->nodelist->intersects[cur_edge->node_index2]->endpt)
	//			streetnet->nodelist->intersects[cur_edge->node_index2]->deleted=true;
	//	}

	//	/*  remove dangling edges  */
	//	if(streetnet->nodelist->intersects[cur_edge->node_index1]->endpt
	//		&&streetnet->nodelist->intersects[cur_edge->node_index2]->endpt)
	//		cur_edge->cancel=true;
	//}

	/*  remove artifact nodes  */
	for(i=0;i<streetnet->nodelist->nelems;i++)
	{
		intersect=streetnet->nodelist->intersects[i];

		if(intersect->deleted)
			continue;

		int nondeleted=0;
		for(j=0;j<intersect->nadjedges;j++)
		{
			if(streetnet->edgelist->edges[intersect->adj_edges[j]]->cancel)
				continue;
			nondeleted++;
		}

		if(nondeleted==0)
		{
			intersect->deleted=true;
		}
	}


	/*   update the whole street network here   */

	/*  a test 1/19/2008 */
	//for(i=0;i<streetnet->edgelist->nedges;i++)
	//{
	//	if(streetnet->edgelist->edges[i]->cancel)
	//		continue;

	//	if(streetnet->nodelist->intersects[streetnet->edgelist->edges[i]->node_index1]->deleted
	//		||streetnet->nodelist->intersects[streetnet->edgelist->edges[i]->node_index2]->deleted)
	//	{
	//		int test=0;
	//	}
	//}
	update_street_network();
}





/*
*/
bool cal_intersect_at_graph_edge(int edgeid, double intersect[2], int &which_samp,
								 double A[2], double B[2])
{
	double C[2], D[2], t[2];
	int i;
	StreetGraphEdge *edge=streetnet->edgelist->edges[edgeid];
	for(i=0;i<edge->ninter_pts-1;i++)
	{
		C[0]=edge->inter_pts[i]->x;
		C[1]=edge->inter_pts[i]->y;
		D[0]=edge->inter_pts[i+1]->x;
		D[1]=edge->inter_pts[i+1]->y;
		if(cal_intersect_2(A, B, C, D, t)==1)
		{
			intersect[0]=A[0]+t[0]*(B[0]-A[0]);
			intersect[1]=A[1]+t[0]*(B[1]-A[1]);
			which_samp=i;
			return true;
		}
	}
	return false;
}



/*  We need to convert the boundary curve into a tensor line, such that
    we can compute the intersection more conveniently
*/
extern void get_linesegs_anytwopts(double p1[2], int cell1, 
						double p2[2], int cell2,
						Trajectory *traj, int type, int MaxIters);

void convert_to_tensorline()
{
	if(regionboundary!=NULL)
	{
		delete regionboundary;
		regionboundary=NULL;
	}
	regionboundary=new Trajectory(-1);
	regionboundary->nlinesegs=0;

	/*  we now consider the input region using the "Region smoothing" interface  */
	int i;
	double start[2], end[2];
	int cell1, cell2;
	for(i=0; i<Num_SmoothRegionpoints; i++)
	{
		start[0] = point[i].x;
		start[1] = point[i].y;
		cell1=point[i].cellid=get_cellID_givencoords(start[0], start[1]);
		end[0] = point[(i+1)%Num_SmoothRegionpoints].x;
		end[1] = point[(i+1)%Num_SmoothRegionpoints].y;
		cell2=point[(i+1)%Num_SmoothRegionpoints].cellid=get_cellID_givencoords(end[0], end[1]);

		get_linesegs_anytwopts(start, cell1, end, cell2, regionboundary, 0, quadmesh->nfaces);
	}
}

/*
    compute the intersection between an edge of the street network graph
	and the boundary of the user specified region.
	We make use of the obtained boundary cell list to locate the intersection
*/
void search_for_intersection_along_edge(int edgeid, double intersect[2], int &which_cell)
{
	int i, pos=-1;
	StreetGraphEdge *edge=streetnet->edgelist->edges[edgeid];

	for(i=0;i<edge->ninter_pts;i++)
	{
		if(is_repeated_elem(boundarycells, edge->inter_pts[i]->cellid, nboundarycells))
		{
			pos=i;
			break;
		}
	}

	/*   we compute the approximate intersection   */
	if(pos>=0)
	{
		int cell=edge->inter_pts[pos]->cellid;

		int start_line_edge, end_line_edge;
		int start_line_bound, end_line_bound;
		bool found=false;
		double A[2], B[2], C[2], D[2], t[2];

		start_line_edge=end_line_edge=
			start_line_bound=end_line_bound=-1;

		/*  search the sample list of the edge  */
		for(i=0;i<edge->ninter_pts;i++)
		{
			if(cell==edge->inter_pts[i]->cellid && !found)
			{
				end_line_edge=start_line_edge=i;
				found=true;
			}
			else if(cell==edge->inter_pts[i]->cellid && found)
			{
				end_line_edge=i;
			}
		}

		/*  search the region boundary tensor line  */
		found=false;
		for(i=0;i<regionboundary->nlinesegs;i++)
		{
			if(cell==regionboundary->linesegs[i].Triangle_ID && !found)
			{
				end_line_bound=start_line_bound=i;
				found=true;
			}
			else if(cell==regionboundary->linesegs[i].Triangle_ID && found)
			{
				end_line_bound=i;
			}
		}

		if(start_line_edge<0 || start_line_bound<0)
			goto LA;

		/*  compute the intersection here  */
		int j;
		//for(i=max(0, start_line_edge-1); i<=min(end_line_edge+1, edge->ninter_pts-1)-1; i++)
		for(i=0; i<edge->ninter_pts-1; i++)
		{
			A[0]=edge->inter_pts[i]->x;
			A[1]=edge->inter_pts[i]->y;
			B[0]=edge->inter_pts[i+1]->x;
			B[1]=edge->inter_pts[i+1]->y;

			for(j=max(0, start_line_bound-2); 
				j<=min(end_line_bound+2, regionboundary->nlinesegs-1); j++)
			{
				C[0]=regionboundary->linesegs[j].gstart[0];
				C[1]=regionboundary->linesegs[j].gstart[1];
				D[0]=regionboundary->linesegs[j].gend[0];
				D[1]=regionboundary->linesegs[j].gend[1];

				if(cal_intersect(A, B, C, D, t)==1)
				{
					intersect[0]=A[0]+t[0]*(B[0]-A[0]);
					intersect[1]=A[1]+t[0]*(B[1]-A[1]);
					which_cell=cell;
					return;
				}
			}
		}

		/*  if we can't find the proper intersection, just use the approximate one  */

LA:		if(!is_inregion(edge->inter_pts[start_line_edge]->x, edge->inter_pts[start_line_edge]->y))
		{
			intersect[0]=edge->inter_pts[start_line_edge]->x;
			intersect[1]=edge->inter_pts[start_line_edge]->y;
		}
		else
		{
			intersect[0]=edge->inter_pts[end_line_edge]->x;
			intersect[1]=edge->inter_pts[end_line_edge]->y;
		}

		which_cell=get_cellID_givencoords(intersect[0],intersect[1]);
	}
	else
	{
		//currently, we leave it there, and return the end point as the intersection
		if(!streetnet->nodelist->intersects[edge->node_index1]->inside_region)
		{
			intersect[0]=streetnet->nodelist->intersects[edge->node_index1]->gpos[0];
			intersect[1]=streetnet->nodelist->intersects[edge->node_index1]->gpos[1];
		}
		else
		{
			intersect[0]=streetnet->nodelist->intersects[edge->node_index2]->gpos[0];
			intersect[1]=streetnet->nodelist->intersects[edge->node_index2]->gpos[1];
		}

		which_cell=get_cellID_givencoords(intersect[0],intersect[1]);

		/*  can we try brute force method?  */
	}
}


/*
	We need to reset the tensor line related information of the inner
	and boundary cells
*/
void reset_inner_and_boundary_cells()
{
	int i;
	QuadCell *face;
	for(i=0;i<ninnercells;i++)
	{
		face=quadmesh->quadcells[i];
		if(face->majorlines!=NULL)
		{
			//free(face->majorlines);
			face->majorlines=NULL;
		}

		face->hasmajor=false;
		face->hasminor=false;
	}
}

/*
    Remove the sub street network inside the specified region
*/

void remove_street_in_reg()
{
	init_quad_regionsmooth();
	find_boundarycells();
	mark_boundVerts();
	find_innerVerts();
	find_inner_cells();

	convert_to_tensorline();
	find_inner_intersections();


	//recal_intersections_on_RegBoundary();
	//recal_intersections_on_RegBoundary_2();
	recal_intersections_on_RegBoundary_3();

	is_on_local_editing=true;
}





/*
    Mark the inner vertices using "cur_max_reg_index+1" index number
*/

void mark_inner_verts_with_new_RegIndex()
{
	int i, j;
	QuadVertex *v;


	for(i=0;i<nregion_quadverts;i++)
	{
		v=quadmesh->quad_verts[region_quadverts[i]];
		v->which_region = cur_max_reg_index+1;
	}

	cur_max_reg_index++;
}


void update_RegIndex_inner_and_boundary_cells()
{
	int i;
	QuadCell *face;

	for(i=0;i<ninnercells;i++)
	{
		face=quadmesh->quadcells[innercells[i]];
		face->which_region=get_region_id_for_cell(innercells[i]);
	}

	for(i=0;i<nboundarycells;i++)
	{
		face=quadmesh->quadcells[boundarycells[i]];
		face->which_region=get_region_id_for_cell(boundarycells[i]);
	}
}


/*
    Retrace the tensor lines inside the user specified region 
*/

void replace_inReg()
{
	/*  construct a seed list based on the intersections on the boundary  */

	SeedList *ini_seeds=new SeedList(nboundaryintersections);
	int i;

	QuadCell *face;
	for(i=0;i<quadmesh->nfaces;i++)
		quadmesh->quadcells[i]->in_region=false;

	ini_seeds->nseeds=0;
	for(i=0;i<nboundaryintersections;i++)
	{
		if(boundaryintersections[i]<0
			||boundaryintersections[i]>=streetnet->nodelist->nelems)
			continue;

		if(!is_inregion(streetnet->nodelist->intersects[boundaryintersections[i]]->gpos[0],
			streetnet->nodelist->intersects[boundaryintersections[i]]->gpos[1]))
			continue;

		Seed *s=(Seed *)malloc(sizeof(Seed));
		s->pos[0]=streetnet->nodelist->intersects[boundaryintersections[i]]->gpos[0];
		s->pos[1]=streetnet->nodelist->intersects[boundaryintersections[i]]->gpos[1];
		s->triangle=streetnet->nodelist->intersects[boundaryintersections[i]]->cellid;
		s->state=0;

		ini_seeds->append(s);

		quadmesh->quadcells[s->triangle]->in_region=true;
	}

	/*  NOTE: we need at least one extra seed in the middle of the region 
	    1/19/2008
	*/
	face=quadmesh->quadcells[innercells[(int)(ninnercells/2)]];
	Seed *s2 = (Seed *)malloc(sizeof(Seed));
	s2->pos[0]=face->x_start_coord+quadmesh->xinterval/2.;
	s2->pos[1]=face->y_start_coord+quadmesh->yinterval/2.;
	s2->triangle=innercells[(int)(ninnercells/2)];
	s2->state=0;
	ini_seeds->append(s2);

	face=quadmesh->quadcells[innercells[(ninnercells-1)]];
	Seed *s3 = (Seed *)malloc(sizeof(Seed));
	s3->pos[0]=face->x_start_coord+quadmesh->xinterval/2.;
	s3->pos[1]=face->y_start_coord+quadmesh->yinterval/2.;
	s3->triangle=innercells[ninnercells-1];
	s3->state=0;
	ini_seeds->append(s3);

	/*  initialize the region  */
	for(i=0;i<quadmesh->nverts;i++)
	{
		quadmesh->quad_verts[i]->InRegion=false;
	}
	for(i=0;i<nregion_quadverts;i++)
		quadmesh->quad_verts[region_quadverts[i]]->InRegion=true;

	/*    reset the information of the inner cells    */
	for(i=0;i<ninnercells;i++)
	{
		face=quadmesh->quadcells[innercells[i]];
		face->in_region=true;
		face->OnBoundary=false;

		/*  intersection list  */
		if(face->intersectlist!=NULL)
		{
			free(face->intersectlist);
			face->intersectlist=NULL;
		}
		face->nintersects=0;

		/*  graph edge list  */
		if(face->streetgraphedgelist!=NULL
			||face->nstreetgraphedges>0)
		{
			free(face->streetgraphedgelist);
			face->streetgraphedgelist=NULL;
		}
		face->nstreetgraphedges=0;

		/*   major tensor line sample point list  */
		if(face->maj_samplepts!=NULL)
		{
			free(face->maj_samplepts);
			face->maj_samplepts=NULL;
		}
		face->maj_nsamplepts=0;
		face->hasmajor=false;
		
		/*   minor tensor line sample point list  */
		if(face->min_samplepts!=NULL)
		{
			free(face->min_samplepts);
			face->min_samplepts=NULL;
		}
		face->min_nsamplepts=0;
		face->hasminor=false;

		if(face->majorlines != NULL)
		{
			//free(face->majorlines);
			delete face->majorlines;
			face->majorlines=NULL;
		}

		if(face->minorlines != NULL)
		{
			//free(face->minorlines);
			delete face->minorlines;
			face->minorlines=NULL;
		}

	}


	/*    reset the information of the boundary cells    */
	//for(i=0;i<nboundarycells;i++)
	//	quadmesh->quadcells[boundarycells[i]]->in_region=true;
	for(i=0;i<nboundarycells;i++)
	{
		face=quadmesh->quadcells[boundarycells[i]];
		face->in_region=true;
		face->OnBoundary=true;

		if(face->maj_samplepts!=NULL)
		{
			free(face->maj_samplepts);
			face->maj_samplepts=NULL;
		}
		face->maj_nsamplepts=0;
		face->hasmajor=false;
		
		if(face->min_samplepts!=NULL)
		{
			free(face->min_samplepts);
			face->min_samplepts=NULL;
		}
		face->min_nsamplepts=0;
		face->hasminor=false;

		if(face->majorlines != NULL)
		{
			//free(face->majorlines);
			delete face->majorlines;
			face->majorlines=NULL;
		}

		if(face->minorlines != NULL)
		{
			//free(face->minorlines);
			delete face->minorlines;
			face->minorlines=NULL;
		}

	}

	/*  start tracing with the ini_seeds as the intial seeds  */

	//if(upstreaming_edit)
	//{
	/*
	possible bug:
	currently, I reset the major and minor tensor line data structure  */
		major->evenstreamlines->ntrajs=minor->evenstreamlines->ntrajs=0;
	//}

	major->place_tensorlines_inReg(0, ini_seeds);

	minor->place_tensorlines_inReg(1, ini_seeds);

}


/*
    We try to compute the sub-graph inside the region
*/

void construct_sub_graph_inReg()
{
	/*   we first find out all the new intersections inside the region   */
	compute_intersections_inReg();

	/*   based on the new obtained tensor lines and their corresponding  
	     line information lists, we search the connections of the sub-graph
	*/
	search_connection_inReg();

	/*  we need to connect the sub graph with existing street network */
}




bool is_repeated_intersect_inCell(int cell, double inter[2], int &which_intersect)
{
	int i;
	QuadCell *face=quadmesh->quadcells[cell];
	icVector2 dist;
	for(i=0;i<face->nintersects;i++)
	{
		dist.entry[0]=inter[0]-streetnet->nodelist->intersects[face->intersectlist[i]]->gpos[0];
		dist.entry[1]=inter[1]-streetnet->nodelist->intersects[face->intersectlist[i]]->gpos[1];

		if(length(dist)==0.0)
		{
			which_intersect=face->intersectlist[i];
			return true;
		}
	}
	return false;
}


/*
   Compute the new intersections inside the region
   Note: we make use of the set of the inner cell list
*/
void compute_intersections_inReg()
{
	int i;
	QuadCell *face;

	/*   Initialize the intersection list in the inner cells   */
	for(i=0;i<ninnercells;i++)
	{
		face=quadmesh->quadcells[innercells[i]];

		if(face->intersectlist!=NULL)
		{
			free(face->intersectlist);
			face->intersectlist=NULL;
		}
		face->nintersects=0;
	}

	/*   Extend the majorintersectinfo list   */
	TensorLineIntersectionInfoList **maj_temp=majorintersectinfo;
	majorintersectinfo = new TensorLineIntersectionInfoList *[major->evenstreamlines->ntrajs];
	for(i=0;i<major->evenstreamlines->ntrajs-major->nnewlines_inReg;i++)
		majorintersectinfo[i]=maj_temp[i];
	delete [] maj_temp;

	for(i=major->evenstreamlines->ntrajs-major->nnewlines_inReg; 
		i<major->evenstreamlines->ntrajs; i++)
		majorintersectinfo[i] = new TensorLineIntersectionInfoList();
	prev_nmajors=major->evenstreamlines->ntrajs;
	
	/*   Extend the minorintersectinfo list   */
	TensorLineIntersectionInfoList **min_temp=minorintersectinfo;
	minorintersectinfo = new TensorLineIntersectionInfoList *[minor->evenstreamlines->ntrajs];
	for(i=0;i<minor->evenstreamlines->ntrajs-minor->nnewlines_inReg;i++)
		minorintersectinfo[i]=min_temp[i];
	delete [] min_temp;

	for(i=minor->evenstreamlines->ntrajs-minor->nnewlines_inReg; 
		i<minor->evenstreamlines->ntrajs; i++)
		minorintersectinfo[i] = new TensorLineIntersectionInfoList();
	prev_nminors=minor->evenstreamlines->ntrajs;

	/*   first, compute the intersections   */
	for(i=0;i<ninnercells;i++)
	{
		face = quadmesh->quadcells[innercells[i]];

		if(face->hasmajor && face->hasminor)
			compute_intersects_in_cell(innercells[i]);
	}

	/*  we also need to consider the intersection inside the boundary cells  */
	for(i=0;i<nboundarycells;i++)
	{
		face=quadmesh->quadcells[boundarycells[i]];

		//if((face->hasmajor && face->majorlines==NULL)
		//	||(face->hasminor && face->minorlines==NULL))
		//{
		//	int test=0;
		//}

		if(face->hasmajor&&face->hasminor)
			//compute_intersect_in_boundarycell(boundarycells[i]);
			compute_intersects_in_cell(boundarycells[i]);
	}


	/*    second, 
	      record the end point of the tensor lines
	*/
	
	Trajectory *temp;
	icVector2 dist;
	Intersection *intersect;

	/*          major tensor lines          */

	for(i=major->evenstreamlines->ntrajs-major->nnewlines_inReg; 
		i<major->evenstreamlines->ntrajs; i++)
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
			int which_intersect;
			if(is_repeated_intersect_inCell(temp->linesegs[0].Triangle_ID, 
				temp->linesegs[0].gstart, which_intersect))
				continue;
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

		add_to_cell_intersectlist(newintersect->cellid, newintersect->index, true);

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

			int which_intersect;
			if(is_repeated_intersect_inCell(temp->linesegs[temp->nlinesegs-1].Triangle_ID, 
				temp->linesegs[temp->nlinesegs-1].gstart, which_intersect))
				continue;
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

		add_to_cell_intersectlist(newintersect2->cellid, newintersect2->index, false);

		/*compute the edges incident to it*/
		IntersectionInfo *majinfo2=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
		majinfo2->intersect_id=streetnet->nodelist->nelems-1;
		majinfo2->lineseg_id = temp->nlinesegs-1;
		majorintersectinfo[i]->sorted_add(majinfo2);
	}

	/*minor tensor lines*/
	for(i=minor->evenstreamlines->ntrajs-minor->nnewlines_inReg; 
		i<minor->evenstreamlines->ntrajs; i++)
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
			int which_intersect;
			if(is_repeated_intersect_inCell(temp->linesegs[0].Triangle_ID, 
				temp->linesegs[0].gstart, which_intersect))
				continue;
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


		//streetnet->danglepts->addNew(newintersect);
		//newintersect2->index += streetnet->nodelist->nelems; /*change the index of it*/
		streetnet->nodelist->addNew(newintersect);
		
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
			dist.entry[0]=temp->linesegs[temp->nlinesegs-1].gstart[0]-intersect->gpos[0];
			dist.entry[1]=temp->linesegs[temp->nlinesegs-1].gstart[1]-intersect->gpos[1];
			if(length(dist)<1e-3) continue;
			int which_intersect;
			if(is_repeated_intersect_inCell(temp->linesegs[temp->nlinesegs-1].Triangle_ID, 
				temp->linesegs[temp->nlinesegs-1].gstart, which_intersect))
				continue;
		}

		temp = minor->evenstreamlines->trajs[i];
		Intersection *newintersect2=(Intersection*)malloc(sizeof(Intersection));
		newintersect2->gpos[0]=temp->linesegs[temp->nlinesegs-1].gstart[0];
		newintersect2->gpos[1]=temp->linesegs[temp->nlinesegs-1].gstart[1];
		newintersect2->cellid=temp->linesegs[temp->nlinesegs-1].Triangle_ID;
		newintersect2->minorline_id=temp->index;
		newintersect2->majorline_id=-1;
		newintersect2->minlineseg=temp->nlinesegs-1;
		newintersect2->nadjedges=0;
		newintersect2->adj_edges=NULL;
		newintersect2->endpt=true;
		newintersect2->inside_region=false;
		newintersect2->deleted=false;

		//streetnet->danglepts->addNew(newintersect2);
		//newintersect2->index += streetnet->nodelist->nelems; /*change the index of it*/
		streetnet->nodelist->addNew(newintersect2);
		
		add_to_cell_intersectlist(newintersect2->cellid, newintersect2->index, false);

		/*compute the edges incident to it*/
		IntersectionInfo *mininfo2=(IntersectionInfo*)malloc(sizeof(IntersectionInfo));
		mininfo2->intersect_id=streetnet->nodelist->nelems-1;
		mininfo2->lineseg_id = temp->nlinesegs-1;
		minorintersectinfo[i]->sorted_add(mininfo2);
	}
}


/*
    Search for the connections based on the obtained new intersections and tensor lines
*/

void search_connection_inReg()
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
	for(i=major->evenstreamlines->ntrajs-major->nnewlines_inReg; 
		i<major->evenstreamlines->ntrajs; i++)
	{
		infolist=majorintersectinfo[i];
		traj=major->evenstreamlines->trajs[i];

		for(j=0; j<infolist->nelems-1; j++)
		{
			info = infolist->infolist[j];
			infonext = infolist->infolist[j+1];

			/*build connection between it and its succeeding*/

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
	for(i=minor->evenstreamlines->ntrajs-minor->nnewlines_inReg; 
		i<minor->evenstreamlines->ntrajs; i++)
	{
		infolist=minorintersectinfo[i];
		traj=minor->evenstreamlines->trajs[i];

		for(j=0; j<infolist->nelems-1; j++)
		{
			info = infolist->infolist[j];
			infonext = infolist->infolist[j+1];

			/*build connection between it and its succeeding*/
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
}



/*
    Compute the level set contour based on the level set propagation
    NOTE: here, we assume always *one* closed contour will be obtained 
*/

void compute_level_set_contour(double dist)
{
	/*  first, we find out the triangles that the contour will pass throught  */
	mark_contour_cells(dist);

	/*  second, we sort the set of the triangles (borrow the idea from 
	    the boundary extraction of the image)
	*/
	get_sorted_contour_cells();

	/*  third, we compute the intersections and connect them to form a curve  */
	compute_the_contour(dist);
}


bool is_contour_cell(int id, double dist)
{
	QuadCell *face=quadmesh->quadcells[id];
	int i;
	bool smaller=false;
	bool larger=false;

	for(i=0;i<face->nverts;i++)
	{
		if(quadmesh->quad_verts[face->verts[i]]->distance>dist)
			larger=true;
		if(quadmesh->quad_verts[face->verts[i]]->distance<=dist)
			smaller=true;

	}

	if(smaller&&larger)
		return true;
	else
		return false;
}


void mark_contour_cells(double dist)
{
	int i;

	for(i=0;i<quadmesh->nfaces;i++)
	{
		quadmesh->quadcells[i]->is_contour_cell=false;
		if(!is_contour_cell(i, dist))
			continue;
		quadmesh->quadcells[i]->is_contour_cell=true;
		//contour_cells[ncontour_cells]=i;
		//ncontour_cells++;
	}
}


/*
   We order the contour cells
*/

void get_sorted_contour_cells()
{
	if(contour_cells!=NULL)
		free(contour_cells);
	contour_cells=(int*)malloc(sizeof(int)*(int)(quadmesh->nfaces/5));
	ncontour_cells=0;

	int i;
	int cur;

	/*   intialization   */
	for(i=0;i<quadmesh->nfaces;i++)
	{
		quadmesh->quadcells[i]->visited=false;
	}

	/*   obtain the first cell   */

	for(i=0;i<quadmesh->nfaces;i++)
	{
		if(quadmesh->quadcells[i]->is_contour_cell)
		{
			contour_cells[ncontour_cells]=i;
			ncontour_cells++;
			break;
		}
	}

	if(ncontour_cells==0) return;

	/*   start searching and sorting the remaining cells   */

	cur=contour_cells[0];
	/*   the first searching order is:
	     right->up->upperright->lower->lowerright->left->lowerleft->upperleft
	*/
	while(1)
	{
		/*rigth*/
		if(quadmesh->quadcells[cur]->x_start_coord+quadmesh->xinterval<quadmesh->xend-1.e-6
			&& quadmesh->quadcells[cur+1]->is_contour_cell
			&& !quadmesh->quadcells[cur+1]->visited)
		{
			cur = cur+1;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		/*up*/
		else if(cur<quadmesh->nfaces-quadmesh->XDIM+1
			&& quadmesh->quadcells[cur+quadmesh->XDIM-1]->is_contour_cell
			&& !quadmesh->quadcells[cur+quadmesh->XDIM-1]->visited)
		{
			cur = cur+quadmesh->XDIM-1;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		/*upper right*/
		else if(quadmesh->quadcells[cur]->x_start_coord+quadmesh->xinterval<quadmesh->xend-1.e-6
			&& quadmesh->quadcells[cur]->y_start_coord+quadmesh->yinterval<quadmesh->yend-1.e-6
			&& quadmesh->quadcells[cur+quadmesh->XDIM]->is_contour_cell
			&& !quadmesh->quadcells[cur+quadmesh->XDIM]->visited)
		{
			cur=cur+quadmesh->XDIM;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		/*lower*/
		else if(cur>=quadmesh->XDIM-1
			&& quadmesh->quadcells[cur-quadmesh->XDIM+1]->is_contour_cell
			&& !quadmesh->quadcells[cur-quadmesh->XDIM+1]->visited)
		{
			cur=cur-quadmesh->XDIM+1;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		/*lowerright*/
		else if(quadmesh->quadcells[cur]->x_start_coord+quadmesh->xinterval<quadmesh->xend-1.e-6
			&& quadmesh->quadcells[cur]->y_start_coord>quadmesh->ystart+1.e-6
			&& quadmesh->quadcells[cur-quadmesh->XDIM+2]->is_contour_cell
			&& !quadmesh->quadcells[cur-quadmesh->XDIM+2]->visited)
		{
			cur=cur-quadmesh->XDIM+2;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		/*left*/
		else if(quadmesh->quadcells[cur]->x_start_coord>quadmesh->xstart+1.e-6
			&& quadmesh->quadcells[cur-1]->is_contour_cell
			&& !quadmesh->quadcells[cur-1]->visited)
		{
			cur=cur-1;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		/*lowerleft*/
		else if(quadmesh->quadcells[cur]->x_start_coord>quadmesh->xstart+1.e-6
			&& quadmesh->quadcells[cur]->y_start_coord>quadmesh->ystart+1.e-6
			&& quadmesh->quadcells[cur-quadmesh->XDIM]->is_contour_cell
			&& !quadmesh->quadcells[cur-quadmesh->XDIM]->visited)
		{
			cur=cur-quadmesh->XDIM;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		/*upperleft*/
		else if(quadmesh->quadcells[cur]->x_start_coord>quadmesh->xstart+1.e-6
			&& quadmesh->quadcells[cur]->y_start_coord+quadmesh->yinterval<quadmesh->yend-1.e-6
			&& quadmesh->quadcells[cur+quadmesh->XDIM-2]->is_contour_cell
			&& !quadmesh->quadcells[cur+quadmesh->XDIM-2]->visited)
		{
			cur=cur+quadmesh->XDIM-2;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		else
			break;
	}

	/*   reverse the obtained cells   */
	int *temp=(int*)malloc(sizeof(int)*ncontour_cells);
	for(i=ncontour_cells-1;i>=0;i--)
	{
		temp[ncontour_cells-1-i]=contour_cells[i];
	}
	for(i=0;i<ncontour_cells;i++)
	{
		contour_cells[i]=temp[i];
	}
	free(temp);

	/*   start from the last cell in the list and follow the pseudo-inversed searching order:
	     left->lower->lowerleft->up->upperleft->right->upperright->lowerright
	*/

	cur=contour_cells[ncontour_cells-1];
	while(1)
	{
		/*left*/
		if(quadmesh->quadcells[cur]->x_start_coord>quadmesh->xstart+1.e-6
			&&quadmesh->quadcells[cur-1]->is_contour_cell
			&& !quadmesh->quadcells[cur-1]->visited)
		{
			cur=cur-1;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		/*lower*/
		else if(cur-quadmesh->XDIM+1>0
			&& quadmesh->quadcells[cur-quadmesh->XDIM+1]->is_contour_cell
			&& !quadmesh->quadcells[cur-quadmesh->XDIM+1]->visited)
		{
			cur=cur-quadmesh->XDIM+1;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		/*lowerleft*/
		else if(quadmesh->quadcells[cur]->x_start_coord>quadmesh->xstart+1.e-6
			&&quadmesh->quadcells[cur]->y_start_coord>quadmesh->ystart+1.e-6
			&& quadmesh->quadcells[cur-quadmesh->XDIM]->is_contour_cell
			&& !quadmesh->quadcells[cur-quadmesh->XDIM]->visited)
		{
			cur=cur-quadmesh->XDIM;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		/*up*/
		else if(cur<quadmesh->nfaces-quadmesh->XDIM+1
			&& quadmesh->quadcells[cur+quadmesh->XDIM-1]->is_contour_cell
			&& !quadmesh->quadcells[cur+quadmesh->XDIM-1]->visited)
		{
			cur = cur+quadmesh->XDIM-1;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		/*upperleft*/
		else if(quadmesh->quadcells[cur]->x_start_coord>quadmesh->xstart+1.e-6
			&& quadmesh->quadcells[cur]->y_start_coord+quadmesh->yinterval<quadmesh->yend-1.e-6
			&& quadmesh->quadcells[cur+quadmesh->XDIM-2]->is_contour_cell
			&& !quadmesh->quadcells[cur+quadmesh->XDIM-2]->visited)
		{
			cur=cur+quadmesh->XDIM-2;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		/*rigth*/
		else if(quadmesh->quadcells[cur]->x_start_coord+quadmesh->xinterval<quadmesh->xend-1.e-6
			&& quadmesh->quadcells[cur+1]->is_contour_cell
			&& !quadmesh->quadcells[cur+1]->visited)
		{
			cur = cur+1;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		/*upper right*/
		else if(quadmesh->quadcells[cur]->x_start_coord+quadmesh->xinterval<quadmesh->xend-1.e-6
			&&quadmesh->quadcells[cur]->y_start_coord+quadmesh->yinterval<quadmesh->yend-1.e-6
			&&quadmesh->quadcells[cur+quadmesh->XDIM]->is_contour_cell
			&& !quadmesh->quadcells[cur+quadmesh->XDIM]->visited)
		{
			cur=cur+quadmesh->XDIM;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		/*lowerright*/
		else if(quadmesh->quadcells[cur]->x_start_coord+quadmesh->xinterval<quadmesh->xend-1.e-6
			&& quadmesh->quadcells[cur]->y_start_coord>quadmesh->ystart+1.e-6
			&&quadmesh->quadcells[cur-quadmesh->XDIM+2]->is_contour_cell
			&& !quadmesh->quadcells[cur-quadmesh->XDIM+2]->visited)
		{
			cur=cur-quadmesh->XDIM+2;
			quadmesh->quadcells[cur]->visited=true;
			contour_cells[ncontour_cells]=cur;
			ncontour_cells++;
		}
		else
			break;
	}
}


/*
    After sorting the cells that contain the level set contour,
	we compute the real contour approximately
	NOTE: we assume that the "contour" cells have been stored at "contour_cells" 
*/

void compute_the_contour(double dist)
{
	int i;
	int from_edge, to_edge;
	double from[2], to[2];

	/*   Initialization   */
	if(regionboundary!=NULL)
	{
		delete regionboundary;
		regionboundary=NULL;
	}
	regionboundary=new Trajectory(-1);
	regionboundary->nlinesegs=0;

	for(i=0;i<ncontour_cells;i++)
	{
		int num_intersects=compute_contour_in_one_cell(contour_cells[i], dist,
			from, from_edge, to, to_edge);
		if(num_intersects==2)
		{
			/*  connect with previous line segment  */
			if(regionboundary->nlinesegs>=regionboundary->curMaxNumLinesegs)
			{
				if(!regionboundary->extend_line_segments(100))
					exit(-1);
			}

			if(regionboundary->nlinesegs>0)
			{
				if((quadmesh->edgelist[from_edge]->tris[0]==
					regionboundary->linesegs[regionboundary->nlinesegs-1].Triangle_ID)
					||(quadmesh->edgelist[from_edge]->tris[1]==
					regionboundary->linesegs[regionboundary->nlinesegs-1].Triangle_ID))
				{
					regionboundary->linesegs[regionboundary->nlinesegs].gstart[0]=from[0];
					regionboundary->linesegs[regionboundary->nlinesegs].gstart[1]=from[1];
					regionboundary->linesegs[regionboundary->nlinesegs].gend[0]=to[0];
					regionboundary->linesegs[regionboundary->nlinesegs].gend[1]=to[1];
					regionboundary->linesegs[regionboundary->nlinesegs].Triangle_ID=contour_cells[i];
				}
				else if((quadmesh->edgelist[to_edge]->tris[0]==
					regionboundary->linesegs[regionboundary->nlinesegs-1].Triangle_ID)
					||(quadmesh->edgelist[to_edge]->tris[1]==
					regionboundary->linesegs[regionboundary->nlinesegs-1].Triangle_ID))
				{
					regionboundary->linesegs[regionboundary->nlinesegs].gstart[0]=to[0];
					regionboundary->linesegs[regionboundary->nlinesegs].gstart[1]=to[1];
					regionboundary->linesegs[regionboundary->nlinesegs].gend[0]=from[0];
					regionboundary->linesegs[regionboundary->nlinesegs].gend[1]=from[1];
					regionboundary->linesegs[regionboundary->nlinesegs].Triangle_ID=contour_cells[i];
				}
				else
				{
					regionboundary->linesegs[regionboundary->nlinesegs].gstart[0]=from[0];
					regionboundary->linesegs[regionboundary->nlinesegs].gstart[1]=from[1];
					regionboundary->linesegs[regionboundary->nlinesegs].gend[0]=to[0];
					regionboundary->linesegs[regionboundary->nlinesegs].gend[1]=to[1];
					regionboundary->linesegs[regionboundary->nlinesegs].Triangle_ID=contour_cells[i];
				}
			}
			else
			{
					regionboundary->linesegs[0].gstart[0]=from[0];
					regionboundary->linesegs[0].gstart[1]=from[1];
					regionboundary->linesegs[0].gend[0]=to[0];
					regionboundary->linesegs[0].gend[1]=to[1];
					regionboundary->linesegs[0].Triangle_ID=contour_cells[i];
			}

			regionboundary->nlinesegs++;

		}
		else
		{
			/*   we can't deal with that right now*/
			int test=0;
		}
	}
}


int compute_contour_in_one_cell(int cellid, double dist, double from[2], int &from_edge,
								 double to[2], int &to_edge)
{
	int i;
	QuadCell *face=quadmesh->quadcells[cellid];
	

	QuadEdge *edge;
	double alpha;
	from_edge=-1;
	to_edge=-1;
	int num_intersects=0;

	for(i=0;i<4;i++)
	{
		edge=face->edges[i];

		if(quadmesh->quad_verts[edge->verts[0]]->distance<dist&&
			quadmesh->quad_verts[edge->verts[1]]->distance>dist)
		{
			alpha=(dist-quadmesh->quad_verts[edge->verts[0]]->distance)/
				(quadmesh->quad_verts[edge->verts[1]]->distance-
				quadmesh->quad_verts[edge->verts[0]]->distance);
			if(from_edge<0)
			{
				from[0]=quadmesh->quad_verts[edge->verts[0]]->x
					+alpha*(quadmesh->quad_verts[edge->verts[1]]->x
					-quadmesh->quad_verts[edge->verts[0]]->x);
				from[1]=quadmesh->quad_verts[edge->verts[0]]->y
					+alpha*(quadmesh->quad_verts[edge->verts[1]]->y
					-quadmesh->quad_verts[edge->verts[0]]->y);
				from_edge=edge->index;
			}
			else
			{
				to[0]=quadmesh->quad_verts[edge->verts[0]]->x
					+alpha*(quadmesh->quad_verts[edge->verts[1]]->x
					-quadmesh->quad_verts[edge->verts[0]]->x);
				to[1]=quadmesh->quad_verts[edge->verts[0]]->y
					+alpha*(quadmesh->quad_verts[edge->verts[1]]->y
					-quadmesh->quad_verts[edge->verts[0]]->y);
				to_edge=edge->index;
			}

			num_intersects++;
		}

		else if(quadmesh->quad_verts[edge->verts[0]]->distance>dist&&
			quadmesh->quad_verts[edge->verts[1]]->distance<dist)
		{
			alpha=(dist-quadmesh->quad_verts[edge->verts[1]]->distance)/
				(quadmesh->quad_verts[edge->verts[0]]->distance-
				quadmesh->quad_verts[edge->verts[1]]->distance);
			if(from_edge<0)
			{
				from[0]=quadmesh->quad_verts[edge->verts[1]]->x
					+alpha*(quadmesh->quad_verts[edge->verts[0]]->x
					-quadmesh->quad_verts[edge->verts[1]]->x);
				from[1]=quadmesh->quad_verts[edge->verts[1]]->y
					+alpha*(quadmesh->quad_verts[edge->verts[0]]->y-
					quadmesh->quad_verts[edge->verts[1]]->y);
				from_edge=edge->index;
			}
			else
			{
				to[0]=quadmesh->quad_verts[edge->verts[1]]->x
					+alpha*(quadmesh->quad_verts[edge->verts[0]]->x
					-quadmesh->quad_verts[edge->verts[1]]->x);
				to[1]=quadmesh->quad_verts[edge->verts[0]]->y
					+alpha*(quadmesh->quad_verts[edge->verts[1]]->y
					-quadmesh->quad_verts[edge->verts[1]]->y);
				to_edge=edge->index;
			}
			num_intersects++;
		}
	}

	return 	num_intersects;
}


/*  We update the node list and edge list of the street network according to the
    "deleted" or "cancel" flag of each node or edge
*/
int *result_node_indices=NULL;
int *result_edge_indices=NULL;

void update_street_network()
{
	int i, j;
	int counter=0;
	StreetGraphEdge *edge;
	Intersection *intersect;

	/*  first, assign the new index for each survival node and edge  */
	if(result_node_indices!=NULL)
		free(result_node_indices);
	result_node_indices=(int*)malloc(sizeof(int)*streetnet->nodelist->nelems);
	for(i=0;i<streetnet->nodelist->nelems;i++)
		result_node_indices[i]=i;
	
	if(result_edge_indices!=NULL)
		free(result_edge_indices);
	result_edge_indices=(int*)malloc(sizeof(int)*streetnet->edgelist->nedges);
	for(i=0;i<streetnet->edgelist->nedges;i++)
		result_edge_indices[i]=i;
	
	/*  update node index  */
	counter=0;
	for(i=0;i<streetnet->nodelist->nelems;i++)
	{
		if(streetnet->nodelist->intersects[i]->deleted)
			continue;
		streetnet->nodelist->intersects[i]->index=counter;
		result_node_indices[i]=counter;
		counter++;
	}

	/*  update edge index  */
	counter=0;
	for(i=0;i<streetnet->edgelist->nedges;i++)
	{
		if(streetnet->edgelist->edges[i]->cancel)
			continue;
		streetnet->edgelist->edges[i]->index=counter;
		result_edge_indices[i]=counter;
		counter++;
	}

	/*  copy nodes according to their new indices update the edgelist in each node  */

	counter=0;
	for(i=0;i<streetnet->nodelist->nelems;i++)
	{
		/**/
		if(streetnet->nodelist->intersects[i]->deleted)
		{
			if(streetnet->nodelist->intersects[i]->adj_edges!=NULL)
			{
				free(streetnet->nodelist->intersects[i]->adj_edges);
				streetnet->nodelist->intersects[i]->adj_edges=NULL;
			}
			free(streetnet->nodelist->intersects[i]);
			streetnet->nodelist->intersects[i]=NULL;
			continue;
		}

		counter++;

		streetnet->nodelist->intersects[result_node_indices[i]]=
			streetnet->nodelist->intersects[i];

		/* update the edge list */
		intersect=streetnet->nodelist->intersects[i];

		//if(/*intersect->nadjedges==0 || */intersect->nadjedges>4)
		//{
		//	int test=0;
		//}

		int *temp_edge_indices=(int*)malloc(sizeof(int)*intersect->nadjedges);
		for(j=0;j<intersect->nadjedges;j++)
			temp_edge_indices[j]=-1;

		int nedge_counter=0;
		for(j=0;j<intersect->nadjedges;j++)
		{
			edge=streetnet->edgelist->edges[intersect->adj_edges[j]];
			if(edge->cancel) 
				continue;

			temp_edge_indices[nedge_counter]=intersect->adj_edges[j];
			nedge_counter++;
		}


		//int *temp_edgelist=intersect->adj_edges;
		free(intersect->adj_edges);
		intersect->adj_edges=(int*)malloc(sizeof(int)*nedge_counter);

		for(j=0;j<nedge_counter;j++)
		{
			/*  obtain the new index of the survival edge   */
			intersect->adj_edges[j]=
				result_edge_indices[temp_edge_indices[j]];

		}
		intersect->nadjedges=nedge_counter;

	}

	/*  reset the remaining space of the intersection list  */
	for(i=counter;i<streetnet->nodelist->nelems;i++)
	{
		streetnet->nodelist->intersects[i]=NULL;
	}
	streetnet->nodelist->nelems=counter;

	/*   !! Update the boundary intersection list  */
	for(i=0;i<nboundaryintersections;i++)
	{
		if(boundaryintersections[i]==result_node_indices[boundaryintersections[i]]
		&&boundaryintersections[i]>=streetnet->nodelist->nelems)
			boundaryintersections[i]=-1;
		else
			boundaryintersections[i]=result_node_indices[boundaryintersections[i]];
	}

	/*  !!Update the inner intersection list  */
	for(i=0;i<ninnerintersections;i++)
	{
		if(innerintersections[i]==result_node_indices[innerintersections[i]]
		&&innerintersections[i]>=streetnet->nodelist->nelems)
			innerintersections[i]=-1;
		else
		innerintersections[i]=result_node_indices[innerintersections[i]];
	}


	/*   update the edge list   */
	counter=0;
	for(i=0;i<streetnet->edgelist->nedges;i++)
	{
		if(streetnet->edgelist->edges[i]->cancel)
		{
			if(streetnet->edgelist->edges[i]->inter_pts!=NULL)
			{
				free(streetnet->edgelist->edges[i]->inter_pts);
				streetnet->edgelist->edges[i]->inter_pts=NULL;
			}
			free(streetnet->edgelist->edges[i]);
			streetnet->edgelist->edges[i]=NULL;
			continue;
		}

		streetnet->edgelist->edges[result_edge_indices[i]]=
			streetnet->edgelist->edges[i];

		/*  update the two end points of the edge  */
		edge=streetnet->edgelist->edges[i];

		edge->node_index1=result_node_indices[edge->node_index1];
		edge->node_index2=result_node_indices[edge->node_index2];

		//if(edge->node_index1>=streetnet->nodelist->nelems
		//	||edge->node_index2>=streetnet->nodelist->nelems)
		//{
		//	int test=0;
		//}
		counter++;
	}
	/*  reset the remaining space of the edge list  */
	for(i=counter;i<streetnet->edgelist->nedges;i++)
	{
		streetnet->edgelist->edges[i]=NULL;
	}
	streetnet->edgelist->nedges=counter;

	/*  Let's make a test here  */
	//for(i=0;i<streetnet->nodelist->nelems;i++)
	//{
	//	Intersection *intersect=streetnet->nodelist->intersects[i];
	//	for(j=0;j<intersect->nadjedges;j++)
	//	{
	//		if(intersect->adj_edges[j]>=streetnet->edgelist->nedges)
	//		{
	//			int test=0;
	//			intersect->adj_edges[j]=result_edge_indices[intersect->adj_edges[j]];
	//		}
	//	}
	//}

	/*-------------------------------------------------------------*/

	/*  update the intersection lists of the cells in the mesh  */

	QuadCell *face;
	for(i=0;i<quadmesh->nfaces;i++)
	{
		face=quadmesh->quadcells[i];

		if(face->nintersects==0 || face->intersectlist==NULL) 
			continue;

		free(face->intersectlist);
		face->intersectlist=NULL;
		face->nintersects=0;

	}

	for(i=0;i<streetnet->nodelist->nelems;i++)
	{
		intersect=streetnet->nodelist->intersects[i];
		add_to_cell_intersectlist(intersect->cellid, i, intersect->endpt);
	}

	/*  update majorlineinfo list  */
	for(i=0;i<major->evenstreamlines->ntrajs;i++)
	{
		for(j=0;j<majorintersectinfo[i]->nelems;j++)
		{
			majorintersectinfo[i]->infolist[j]->intersect_id=
				result_node_indices[majorintersectinfo[i]->infolist[j]->intersect_id];
		}
	}

	/*  update minorlineinfo list  */
	for(i=0;i<minor->evenstreamlines->ntrajs;i++)
	{
		for(j=0;j<minorintersectinfo[i]->nelems;j++)
		{
			minorintersectinfo[i]->infolist[j]->intersect_id=
				result_node_indices[minorintersectinfo[i]->infolist[j]->intersect_id];
		}
	}

	/*-------------------------------------------------------------*/

	/*  update the edge lists of the cells in the mesh   */

	for(i=0;i<quadmesh->nfaces;i++)
	{
		face=quadmesh->quadcells[i];

		if(face->nstreetgraphedges==0 || face->streetgraphedgelist==NULL ) continue;

		free(face->streetgraphedgelist);
		face->streetgraphedgelist=NULL;
		face->nstreetgraphedges=0;

	}

	int pre_cell=-1;
	for(i=0;i<streetnet->edgelist->nedges;i++)
	{
		edge=streetnet->edgelist->edges[i];
		pre_cell=-1;

		for(j=0;j<edge->ninter_pts;j++)
		{
			int test_cell=edge->inter_pts[j]->cellid;
			if(test_cell<0 ||test_cell>=quadmesh->nfaces)
			{
				test_cell=edge->inter_pts[j]->cellid=
					get_cellID_givencoords(edge->inter_pts[j]->x,
					edge->inter_pts[j]->y);
				if(test_cell<0 ||test_cell>=quadmesh->nfaces)
					continue;
			}

			if(pre_cell==edge->inter_pts[j]->cellid)
				continue;

			/*   add to the cell edge list   */
			pre_cell=edge->inter_pts[j]->cellid;
			if(pre_cell<0 ||pre_cell>=quadmesh->nfaces)
			{
				pre_cell=get_cellID_givencoords(edge->inter_pts[j]->x,
					edge->inter_pts[j]->y);
				if(pre_cell<0 ||pre_cell>=quadmesh->nfaces)
					continue;
			}

			QuadCell *face=quadmesh->quadcells[pre_cell];

			if(face->streetgraphedgelist==NULL
				||face->nstreetgraphedges==0)
				face->streetgraphedgelist=
				(int*)malloc(sizeof(int));
			else
				face->streetgraphedgelist=
					extend_link(face->streetgraphedgelist,
					face->nstreetgraphedges);

			face->streetgraphedgelist[
				face->nstreetgraphedges]=i;
				face->nstreetgraphedges++;
		}

	}


	/*  for undo operation, we need to save the address information of those being deleted
	    intersections and edges  12/15/2007
	*/
}

/*
   Merge the intersections that are on (close to) the region boundary
void merge_close_intersections_on_RegBoundary(double threshold)
*/

void merge_subgraph_to_streetnet(bool majormin, bool startorend)
{
	icVector2 dir;
	int i;
	double A[2], ang, dist_intersect, dist_edge;
	Trajectory *traj;
	int cell_id;
	int which_intersect, which_edge, which_samp;

	int cur_intersect_id;
	EvenStreamlinePlace *cur_place;
	TensorLineIntersectionInfoList **infolist;

	int except_edge;

	if(!majormin)
	{
		cur_place=major;
		infolist=majorintersectinfo;
	}
	else{
		cur_place=minor;
		infolist=minorintersectinfo;
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

			dir.entry[0]=traj->linesegs[traj->nlinesegs-1].gend[0]
				-traj->linesegs[traj->nlinesegs-1].gstart[0];
			dir.entry[1]=traj->linesegs[traj->nlinesegs-1].gend[1]
				-traj->linesegs[traj->nlinesegs-1].gstart[1];
			normalize(dir);
			cell_id=traj->linesegs[traj->nlinesegs-1].Triangle_ID;
		}
		if(!quadmesh->quadcells[cell_id]->OnBoundary)
			return;

		/*  obtain the edge that is associated with the intersection  */
		except_edge=streetnet->nodelist->intersects[cur_intersect_id]->adj_edges[0];

		//if(find_closest_intersection(A, dist_intersect, which_intersect))
		{
			/*  extend to the intersection  */
			/*  merge to the intersection "which_intersect"  */
			/*  update "cur_intersect_id"  */
			streetnet->nodelist->intersects[cur_intersect_id]->deleted=true;
			StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];

			if(edge->node_index1==cur_intersect_id)
			{
				edge->node_index1=which_intersect;
				/*update the sampling list*/
				edge->inter_pts[0]->x=streetnet->nodelist->intersects[which_intersect]->gpos[0];
				edge->inter_pts[0]->y=streetnet->nodelist->intersects[which_intersect]->gpos[1];
				edge->inter_pts[0]->cellid=streetnet->nodelist->intersects[which_intersect]->cellid;
			}
			else
			{
				edge->node_index2=which_intersect;
				/*update the sampling list*/
				edge->inter_pts[edge->ninter_pts-1]->x
					=streetnet->nodelist->intersects[which_intersect]->gpos[0];
				edge->inter_pts[edge->ninter_pts-1]->y
					=streetnet->nodelist->intersects[which_intersect]->gpos[1];
				edge->inter_pts[edge->ninter_pts-1]->cellid
					=streetnet->nodelist->intersects[which_intersect]->cellid;
			}

			streetnet->nodelist->intersects[which_intersect]->endpt=false;
			continue;
		}

		ang=atan2(dir.entry[1], dir.entry[0]);
		if(ang<0) ang+=(2.*M_PI);


		/*  compute the distance between A and other nearby intersections  */
		if(cal_smallest_dist_to_intersects(cell_id, A, ang, dir, dist_intersect, which_intersect,
			cur_intersect_id))
		{
			//dist_intersect=temp_dist;
			if(which_intersect!=cur_intersect_id)
			closetointersect=true;
		}

		/*  compute the distance between A and other nearby edges of the street network  */
		if(cal_smallest_dist_to_edges(cell_id, A, ang, dir, dist_edge, which_edge, which_samp, except_edge))
		{
			//dist_edge=temp_dist;
			if(which_edge != except_edge)
			closetoedge=true;
		}

		//if(closetoedge && streetnet->edgelist->edges[which_edge]->ninter_pts<which_samp)
		//{
		//	int test=0;
		//}

		if(closetointersect && closetoedge)
		{
			/**/
			if(dist_intersect<=dist_edge)
			{
				/*  extend to the intersection  */
				/*  merge to the intersection "which_intersect"  */
				/*  update "cur_intersect_id"  */
				streetnet->nodelist->intersects[cur_intersect_id]->deleted=true;
				StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];

				if(edge->node_index1==cur_intersect_id)
				{
					edge->node_index1=which_intersect;
					/*update the sampling list*/
					edge->inter_pts[0]->x=streetnet->nodelist->intersects[which_intersect]->gpos[0];
					edge->inter_pts[0]->y=streetnet->nodelist->intersects[which_intersect]->gpos[1];
					edge->inter_pts[0]->cellid=streetnet->nodelist->intersects[which_intersect]->cellid;
				}
				else
				{
					edge->node_index2=which_intersect;
					/*update the sampling list*/
					edge->inter_pts[edge->ninter_pts-1]->x
						=streetnet->nodelist->intersects[which_intersect]->gpos[0];
					edge->inter_pts[edge->ninter_pts-1]->y
						=streetnet->nodelist->intersects[which_intersect]->gpos[1];
					edge->inter_pts[edge->ninter_pts-1]->cellid
						=streetnet->nodelist->intersects[which_intersect]->cellid;
				}

				streetnet->nodelist->intersects[which_intersect]->endpt=false;
			}
			else
			{
				/*  extend to the edge  */
				
				streetnet->nodelist->intersects[cur_intersect_id]->endpt=false;
				//streetnet->nodelist->intersects[cur_intersect_id]->gpos[0]=
				//	(streetnet->edgelist->edges[which_edge]->inter_pts[which_samp]->x+
				//	streetnet->edgelist->edges[which_edge]->inter_pts[which_samp+1]->x)/2.;
				//streetnet->nodelist->intersects[cur_intersect_id]->gpos[1]=
				//	(streetnet->edgelist->edges[which_edge]->inter_pts[which_samp]->y+
				//	streetnet->edgelist->edges[which_edge]->inter_pts[which_samp+1]->y)/2.;

				/*   use the first point of the line segment   */
				streetnet->nodelist->intersects[cur_intersect_id]->gpos[0]=
					streetnet->edgelist->edges[which_edge]->inter_pts[which_samp]->x;
				streetnet->nodelist->intersects[cur_intersect_id]->gpos[1]=
					streetnet->edgelist->edges[which_edge]->inter_pts[which_samp]->y;
				StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];

				if(edge->node_index1==cur_intersect_id)
				{
					edge->inter_pts[0]->x=
						streetnet->nodelist->intersects[cur_intersect_id]->gpos[0];
					edge->inter_pts[0]->y=
						streetnet->nodelist->intersects[cur_intersect_id]->gpos[1];
				}
				else
				{
					edge->inter_pts[edge->ninter_pts-1]->x=
						streetnet->nodelist->intersects[cur_intersect_id]->gpos[0];
					edge->inter_pts[edge->ninter_pts-1]->y=
						streetnet->nodelist->intersects[cur_intersect_id]->gpos[1];
				}

				/*  update the original edge and   
					generate one new edges in the street network (graph)  
				*/
				edge=streetnet->edgelist->edges[which_edge];
				StreetGraphEdge *newedge=(StreetGraphEdge*)malloc(sizeof(StreetGraphEdge));
				newedge->node_index1=cur_intersect_id;
				newedge->node_index2=edge->node_index2;
				newedge->cancel = newedge->visited = false;
				newedge->inter_pts=(Point**)malloc(sizeof(Point*)*(edge->ninter_pts-which_samp));
				int j;
				for(j=which_samp;j<edge->ninter_pts;j++)
					newedge->inter_pts[j-which_samp]=edge->inter_pts[j];
				newedge->ninter_pts=edge->ninter_pts-which_samp;
				newedge->roadtype=edge->roadtype;
				
				streetnet->edgelist->append(newedge);

				/*  update the old edge  */
				edge->node_index2=cur_intersect_id;

				Point **temp=edge->inter_pts;
				edge->inter_pts=(Point**)malloc(sizeof(Point*)*(which_samp+1));
				for(j=0;j<=which_samp;j++)
					edge->inter_pts[j]=temp[j];
				edge->ninter_pts=which_samp+1;
				free(temp);

				/*  update the edge list of the intersection "cur_intersect_id"  */
				add_edge_to_intersectnode(cur_intersect_id, edge->index);
				add_edge_to_intersectnode(cur_intersect_id, newedge->index);

			}
		}

		else if(closetointersect)
		{
			/*  extend to the intersection  */
			/*  merge to the intersection "which_intersect"  */
			/*  update "cur_intersect_id"  */
			streetnet->nodelist->intersects[cur_intersect_id]->deleted=true;
			StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];
			if(edge->node_index1==cur_intersect_id)
			{
				edge->node_index1=which_intersect;
				/*update the sampling list*/
				edge->inter_pts[0]->x=streetnet->nodelist->intersects[which_intersect]->gpos[0];
				edge->inter_pts[0]->y=streetnet->nodelist->intersects[which_intersect]->gpos[1];
				edge->inter_pts[0]->cellid=streetnet->nodelist->intersects[which_intersect]->cellid;
			}
			else
			{
				edge->node_index2=which_intersect;
				/*update the sampling list*/
				edge->inter_pts[edge->ninter_pts-1]->x
					=streetnet->nodelist->intersects[which_intersect]->gpos[0];
				edge->inter_pts[edge->ninter_pts-1]->y
					=streetnet->nodelist->intersects[which_intersect]->gpos[1];
				edge->inter_pts[edge->ninter_pts-1]->cellid
					=streetnet->nodelist->intersects[which_intersect]->cellid;
			}
				streetnet->nodelist->intersects[which_intersect]->endpt=false;

		}

		else if(closetoedge)
		{
			/*  extend to the edge   */

				streetnet->nodelist->intersects[cur_intersect_id]->endpt=false;
				//streetnet->nodelist->intersects[cur_intersect_id]->gpos[0]=
				//	(streetnet->edgelist->edges[which_edge]->inter_pts[which_samp]->x+
				//	streetnet->edgelist->edges[which_edge]->inter_pts[which_samp+1]->x)/2.;
				//streetnet->nodelist->intersects[cur_intersect_id]->gpos[1]=
				//	(streetnet->edgelist->edges[which_edge]->inter_pts[which_samp]->y+
				//	streetnet->edgelist->edges[which_edge]->inter_pts[which_samp+1]->y)/2.;

				/*   use the first point of the line segment   */
				streetnet->nodelist->intersects[cur_intersect_id]->gpos[0]=
					streetnet->edgelist->edges[which_edge]->inter_pts[which_samp]->x;
				streetnet->nodelist->intersects[cur_intersect_id]->gpos[1]=
					streetnet->edgelist->edges[which_edge]->inter_pts[which_samp]->y;
				StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];

				if(edge->node_index1==cur_intersect_id)
				{
					edge->inter_pts[0]->x=
						streetnet->nodelist->intersects[cur_intersect_id]->gpos[0];
					edge->inter_pts[0]->y=
						streetnet->nodelist->intersects[cur_intersect_id]->gpos[1];
				}
				else
				{
					edge->inter_pts[edge->ninter_pts-1]->x=
						streetnet->nodelist->intersects[cur_intersect_id]->gpos[0];
					edge->inter_pts[edge->ninter_pts-1]->y=
						streetnet->nodelist->intersects[cur_intersect_id]->gpos[1];
				}

				/*  update the original edge and   
					generate one new edges in the street network (graph)  
				*/
				edge=streetnet->edgelist->edges[which_edge];
				StreetGraphEdge *newedge=(StreetGraphEdge*)malloc(sizeof(StreetGraphEdge));
				newedge->node_index1=cur_intersect_id;
				newedge->node_index2=edge->node_index2;
				newedge->cancel = newedge->visited = false;
				newedge->inter_pts=(Point**)malloc(sizeof(Point*)*(edge->ninter_pts-which_samp));
				int j;
				for(j=which_samp;j<edge->ninter_pts;j++)
					newedge->inter_pts[j-which_samp]=edge->inter_pts[j];
				newedge->ninter_pts=edge->ninter_pts-which_samp;
				newedge->roadtype=edge->roadtype;
				
				streetnet->edgelist->append(newedge);

				/*  update the old edge  */
				edge->node_index2=cur_intersect_id;

				Point **temp=edge->inter_pts;
				edge->inter_pts=(Point**)malloc(sizeof(Point*)*(which_samp+1));
				for(j=0;j<=which_samp;j++)
					edge->inter_pts[j]=temp[j];
				edge->ninter_pts=which_samp+1;
				free(temp);

				/*  update the edge list of the intersection "cur_intersect_id"  */
				add_edge_to_intersectnode(cur_intersect_id, edge->index);
				add_edge_to_intersectnode(cur_intersect_id, newedge->index);
		}

		else
		{
			/*  just leave it there  */
		}

		/*   deal with end point  */
	}
}


/*
	calculate the distance between A and other nearby intersections
	using the intersection list stored in each nearby cells
*/
bool cal_smallest_dist_to_intersects(int cellid, double A[2], double ang, icVector2 dir,
									 double &dist, int &which_intersect, int except_intersect)
{
	dist=1.e50;
	which_intersect=-1;
	bool flag=false;
	if(ang>=0 && ang<=M_PI/2.)
	{
		/*  we consider right, upper and upper right cells  */
		if(quadmesh->quadcells[cellid]->x_start_coord+quadmesh->xinterval<quadmesh->xend-1.e-8)
		{
			/*right*/
			double temp_dist=1.e50;
			int temp_intersect;
			if(cal_smallest_dist_from_one_cell(cellid+1, A, dir, temp_dist, temp_intersect,except_intersect))
			{
				if(dist>temp_dist)
				{
					dist=temp_dist;
					which_intersect=temp_intersect;
				}
				//if(dist>quadmesh->xinterval*2)
				//{
				//	int test=0;
				//}
				flag=true;
			}
		}
		if(cellid<quadmesh->nfaces-(quadmesh->XDIM-1))
		{
			/*up*/
			double temp_dist;
			int temp_intersect;
			if(cal_smallest_dist_from_one_cell(cellid+quadmesh->XDIM-1, 
				A, dir, temp_dist, temp_intersect,except_intersect))
			{
				if(dist>temp_dist)
				{
					dist=temp_dist;
					which_intersect=temp_intersect;
				}

				//if(dist>quadmesh->xinterval*2)
				//{
				//	int test=0;
				//}

				flag=true;
			}
		}
		if((quadmesh->quadcells[cellid]->x_start_coord+quadmesh->xinterval<quadmesh->xend-1.e-6)
			&&quadmesh->quadcells[cellid]->y_start_coord+quadmesh->yinterval<quadmesh->yend-1.e-6)
		{
			/*upper right*/
			double temp_dist;
			int temp_intersect;
			if(cal_smallest_dist_from_one_cell(cellid+quadmesh->XDIM, 
				A, dir, temp_dist, temp_intersect,except_intersect))
			{
				if(dist>temp_dist)
				{
					dist=temp_dist;
					which_intersect=temp_intersect;
				}
				//if(dist>quadmesh->xinterval*2)
				//{
				//	int test=0;
				//}
				flag=true;
			}
		}
	}


	else if(ang>M_PI/2 && ang<=M_PI)
	{
		/*  we consider up, left and upper left cells  */
		if(cellid<quadmesh->nfaces-(quadmesh->XDIM-1))
		{
			/*up*/
			double temp_dist;
			int temp_intersect;
			if(cal_smallest_dist_from_one_cell(cellid+quadmesh->XDIM-1, 
				A, dir, temp_dist, temp_intersect,except_intersect))
			{
				if(dist>temp_dist)
				{
					dist=temp_dist;
					which_intersect=temp_intersect;
				}
				//if(dist>quadmesh->xinterval*2)
				//{
				//	int test=0;
				//}
				flag=true;
			}
		}
		if(quadmesh->quadcells[cellid]->x_start_coord>quadmesh->xstart+1.e-6)
		{
			/*left*/
			double temp_dist;
			int temp_intersect;
			if(cal_smallest_dist_from_one_cell(cellid-1, 
				A, dir, temp_dist, temp_intersect,except_intersect))
			{
				if(dist>temp_dist)
				{
					dist=temp_dist;
					which_intersect=temp_intersect;
				}
				//if(dist>quadmesh->xinterval*2)
				//{
				//	int test=0;
				//}
				flag=true;
			}
		}
		if((quadmesh->quadcells[cellid]->x_start_coord>quadmesh->xstart+1.e-6)
			&& cellid<quadmesh->nfaces-(quadmesh->XDIM-1))
		{
			/* upper left */
			double temp_dist;
			int temp_intersect;
			if(cal_smallest_dist_from_one_cell(cellid+quadmesh->XDIM-2, 
				A, dir, temp_dist, temp_intersect,except_intersect))
			{
				if(dist>temp_dist)
				{
					dist=temp_dist;
					which_intersect=temp_intersect;
				}
				//if(dist>quadmesh->xinterval*2)
				//{
				//	int test=0;
				//}
				flag=true;
			}
		}

	}


	else if(ang>M_PI && ang<M_PI/2*3)
	{
		/*  we consider left, lower and lower left cells  */ 
		if(quadmesh->quadcells[cellid]->x_start_coord>quadmesh->xstart+1.e-6)
		{
			/*left*/
			double temp_dist;
			int temp_intersect;
			if(cal_smallest_dist_from_one_cell(cellid-1, 
				A, dir, temp_dist, temp_intersect,except_intersect))
			{
				if(dist>temp_dist)
				{
					dist=temp_dist;
					which_intersect=temp_intersect;
				}
				//if(dist>quadmesh->xinterval*2)
				//{
				//	int test=0;
				//}
				flag=true;
			}
		}

		if(cellid>=quadmesh->XDIM-1)
		{
			/*  lower  */
			double temp_dist;
			int temp_intersect;
			if(cal_smallest_dist_from_one_cell(cellid-quadmesh->XDIM+1, 
				A, dir, temp_dist, temp_intersect,except_intersect))
			{
				if(dist>temp_dist)
				{
					dist=temp_dist;
					which_intersect=temp_intersect;
				}
				//if(dist>quadmesh->xinterval*2)
				//{
				//	int test=0;
				//}
				flag=true;
			}
		}
		if((quadmesh->quadcells[cellid]->x_start_coord>quadmesh->xstart+1.e-6)
			&&cellid>=quadmesh->XDIM-1)
		{
			/*  lower left  */
			double temp_dist;
			int temp_intersect;
			if(cal_smallest_dist_from_one_cell(cellid-quadmesh->XDIM, 
				A, dir, temp_dist, temp_intersect,except_intersect))
			{
				if(dist>temp_dist)
				{
					dist=temp_dist;
					which_intersect=temp_intersect;
				}
				//if(dist>quadmesh->xinterval*2)
				//{
				//	int test=0;
				//}
				flag=true;
			}
		}
	}


	else
	{
		/*  we consider lower, right and lower right cells  */
		if(cellid>=quadmesh->XDIM-1)
		{
			/*  lower  */
			double temp_dist;
			int temp_intersect;
			if(cal_smallest_dist_from_one_cell(cellid-quadmesh->XDIM+1, 
				A, dir, temp_dist, temp_intersect,except_intersect))
			{
				if(dist>temp_dist)
				{
					dist=temp_dist;
					which_intersect=temp_intersect;
				}
				//if(dist>quadmesh->xinterval*2)
				//{
				//	int test=0;
				//}
				flag=true;
			}
		}
		if(quadmesh->quadcells[cellid]->x_start_coord+quadmesh->xinterval<quadmesh->xend-1.e-6)
		{
			/*right*/
			double temp_dist;
			int temp_intersect;
			if(cal_smallest_dist_from_one_cell(cellid+1, A, dir, temp_dist, temp_intersect,except_intersect))
			{
				if(dist>temp_dist)
				{
					dist=temp_dist;
					which_intersect=temp_intersect;
				}
				//if(dist>quadmesh->xinterval*2)
				//{
				//	int test=0;
				//}
				flag=true;
			}
		}
		if((quadmesh->quadcells[cellid]->x_start_coord+quadmesh->xinterval<quadmesh->xend-1.e-6)
			&&cellid>=quadmesh->XDIM-1)
		{
			/*  lower right */
			double temp_dist;
			int temp_intersect;
			if(cal_smallest_dist_from_one_cell(cellid-quadmesh->XDIM+2, 
				A, dir, temp_dist, temp_intersect,except_intersect))
			{
				if(dist>temp_dist)
				{
					dist=temp_dist;
					which_intersect=temp_intersect;
				}
				//if(dist>quadmesh->xinterval*2)
				//{
				//	int test=0;
				//}
				flag=true;
			}
		}
	}

	return flag;
}


/*
   Calculate the smallest distance to the intersections of *one* nearby cell  
*/
bool cal_smallest_dist_from_one_cell(int cellid, double A[2], icVector2 dir,
									 double &smallestdist, int &which_intersect, int except_intersect)
{
	int i;
	QuadCell *face=quadmesh->quadcells[cellid];
	if(face->nintersects==0) return false;

	icVector2 dist;
	smallestdist=1.e50;
	which_intersect=-1;
	for(i=0;i<face->nintersects;i++)
	{
		if(face->intersectlist[i]==except_intersect) 
			continue;

		if(streetnet->nodelist->intersects[face->intersectlist[i]]->deleted)
			continue;

		dist.entry[0]=streetnet->nodelist->intersects[face->intersectlist[i]]->gpos[0]-A[0];
		dist.entry[1]=streetnet->nodelist->intersects[face->intersectlist[i]]->gpos[1]-A[1];
		double len=length(dist);

		//if(len==0.0)
		//{
		//	int test=0;
		//}

		if(smallestdist>len)
		{
			normalize(dist);
			if(len>1.e-2 && dot(dist, dir)<0.96) 
				continue;

			smallestdist=len;
			which_intersect=face->intersectlist[i];
		}
	}

	if(which_intersect>=0) return true;
	return false;
}


//
/*
   Calculate the smallest distance to one nearby edge of the street network
void DistanceFromLine(double cx, double cy, double ax, double ay ,
//					  double bx, double by, double &distanceSegment,
//					  double &distanceLine)
*/
double cal_smallest_dist_to_one_edge(int edgeid, double p[2], int &which_samp)
{
	int i;
	StreetGraphEdge *edge=streetnet->edgelist->edges[edgeid];
	double smallestdist=1.e50;
	which_samp=-1;
	double A[2], B[2];
	double distSeg, distLine;
	for(i=0;i<edge->ninter_pts-1;i++)
	{
		A[0]=edge->inter_pts[i]->x;
		A[1]=edge->inter_pts[i]->y;
		B[0]=edge->inter_pts[i+1]->x;
		B[1]=edge->inter_pts[i+1]->y;
		DistanceFromLine(p[0],p[1],  A[0],A[1],  B[0],B[1],  distSeg, distLine);
		if(distSeg<smallestdist)
		{
			//if(distSeg<1.e-6)
			//{
			//	int test=0;
			//	continue;
			//}
			smallestdist=distSeg;
			which_samp=i;
		}
	}

	return smallestdist;
}


bool cal_smallest_dist_to_one_edge_in_a_cell(int cellid, double A[2], icVector2 dir, double &dist,
											   int &which_edge, int &which_samp, int except_edge)
{
	
	int i;
	QuadCell *face;
	double temp_dist;
	int temp_samp;
	icVector2 check_dir_m, check_dir_s, check_dir_e;
	StreetGraphEdge *edge;

	face=quadmesh->quadcells[cellid];

	if(face->nstreetgraphedges==0||face->streetgraphedgelist==NULL) return false;

	bool updateflag=false;

	for(i=0;i<face->nstreetgraphedges;i++)
	{
		if(face->streetgraphedgelist[i]==except_edge) continue;

		edge=streetnet->edgelist->edges[face->streetgraphedgelist[i]];
		if(edge->cancel)
			continue;

		check_dir_m.entry[0]=(edge->inter_pts[0]->x+edge->inter_pts[edge->ninter_pts-1]->x)/2.-
			A[0];
		check_dir_m.entry[1]=(edge->inter_pts[0]->y+edge->inter_pts[edge->ninter_pts-1]->y)/2.-
			A[1];

		check_dir_s.entry[0]=edge->inter_pts[0]->x-A[0];
		check_dir_s.entry[1]=edge->inter_pts[0]->y-A[1];
		
		check_dir_e.entry[0]=edge->inter_pts[edge->ninter_pts-1]->x-A[0];
		check_dir_e.entry[1]=edge->inter_pts[edge->ninter_pts-1]->y-A[1];

		normalize(check_dir_m);
		normalize(check_dir_s);
		normalize(check_dir_e);

		double dot1=dot(check_dir_m, dir);
		double dot2=dot(check_dir_s, dir);
		double dot3=dot(check_dir_e, dir);

		if(dot1<dot2) dot1=dot2;
		if(dot1<dot3) dot1=dot3;

		if(dot1<0.96) 
			continue;

		temp_dist=cal_smallest_dist_to_one_edge(face->streetgraphedgelist[i], 
			A, temp_samp);
		if(dist>temp_dist)
		{
			dist=temp_dist;
			which_edge=face->streetgraphedgelist[i];
			which_samp=temp_samp;

			if(which_samp==0)
				which_samp=1;

			updateflag=true;
		}
	}

	return updateflag;
}

/*
	calculate the distance between A and the nearby edges of the street network
	in the nearby cells
*/
bool cal_smallest_dist_to_edges(int &cellid, double A[2], double ang, icVector2 dir,
									 double &dist, int &which_edge, int &which_samp, int except_edge)
{
	dist=1.e50;
	which_edge=-1;
	which_samp=-1;
	bool flag=false;
	int which_cell;

	if(ang>=0 && ang<=M_PI/2.)
	{
		/*  we consider right, upper and upper right cells  */
		if(quadmesh->quadcells[cellid]->x_start_coord+quadmesh->xinterval<quadmesh->xend-1.e-6)
		{
			/*right*/
			which_cell=cellid+1;
			if(cal_smallest_dist_to_one_edge_in_a_cell(which_cell, A, dir, dist, which_edge, which_samp,
				except_edge))
				flag=true;
		}
		if(cellid<quadmesh->nfaces-(quadmesh->XDIM-1))
		{
			/*up*/
			which_cell=cellid+quadmesh->XDIM-1;
			if(cal_smallest_dist_to_one_edge_in_a_cell(which_cell, A, dir, dist, which_edge, which_samp,
				except_edge))
				flag=true;
		}
		if((quadmesh->quadcells[cellid]->x_start_coord+quadmesh->xinterval<quadmesh->xend-1.e-6)
			&&quadmesh->quadcells[cellid]->y_start_coord+quadmesh->yinterval<quadmesh->yend-1.e-6)
		{
			/*upper right*/
			which_cell=cellid+quadmesh->XDIM;
			if(cal_smallest_dist_to_one_edge_in_a_cell(which_cell, A, dir, dist, which_edge, which_samp,
				except_edge))
				flag=true;
		}
	}


	else if(ang>M_PI/2 && ang<=M_PI)
	{
		/*  we consider up, left and upper left cells  */
		if(cellid<quadmesh->nfaces-(quadmesh->XDIM-1))
		{
			/*up*/
			which_cell=cellid+quadmesh->XDIM-1;
			if(cal_smallest_dist_to_one_edge_in_a_cell(which_cell, A, dir, dist, which_edge, which_samp,
				except_edge))
				flag=true;
		}
		if(quadmesh->quadcells[cellid]->x_start_coord>quadmesh->xstart+1.e-6)
		{
			/*left*/
			which_cell=cellid-1;
			if(cal_smallest_dist_to_one_edge_in_a_cell(which_cell, A, dir, dist, which_edge, which_samp,
				except_edge))
				flag=true;
		}
		if((quadmesh->quadcells[cellid]->x_start_coord>quadmesh->xstart+1.e-6)
			&& cellid<quadmesh->nfaces-(quadmesh->XDIM-1))
		{
			/* upper left */
			which_cell=cellid+quadmesh->XDIM-2;
			if(cal_smallest_dist_to_one_edge_in_a_cell(which_cell, A, dir, dist, which_edge, which_samp,
				except_edge))
				flag=true;
		}

	}


	else if(ang>M_PI && ang<M_PI/2*3)
	{
		/*  we consider left, lower and lower left cells  */ 
		if(quadmesh->quadcells[cellid]->x_start_coord>quadmesh->xstart+1.e-6)
		{
			/*left*/
			which_cell=cellid-1;
			if(cal_smallest_dist_to_one_edge_in_a_cell(which_cell, A, dir, dist, which_edge, which_samp,
				except_edge))
				flag=true;
		}

		if(cellid>=quadmesh->XDIM-1)
		{
			/*  lower  */
			which_cell=cellid-quadmesh->XDIM+1;
			if(cal_smallest_dist_to_one_edge_in_a_cell(which_cell, A, dir, dist, which_edge, which_samp,
				except_edge))
				flag=true;
		}
		if((quadmesh->quadcells[cellid]->x_start_coord>quadmesh->xstart+1.e-6)
			&&cellid>=quadmesh->XDIM-1)
		{
			/*  lower left  */
			which_cell=cellid-quadmesh->XDIM;
			if(cal_smallest_dist_to_one_edge_in_a_cell(which_cell, A, dir, dist, which_edge, which_samp,
				except_edge))
				flag=true;
		}
	}


	else
	{
		/*  we consider lower, right and lower right cells  */
		if(cellid>=quadmesh->XDIM-1)
		{
			/*  lower  */
			which_cell=cellid-quadmesh->XDIM+1;
			if(cal_smallest_dist_to_one_edge_in_a_cell(which_cell, A, dir, dist, which_edge, which_samp,
				except_edge))
				flag=true;
		}
		if(quadmesh->quadcells[cellid]->x_start_coord+quadmesh->xinterval<quadmesh->xend-1.e-6)
		{
			/*right*/
			which_cell=cellid+1;
			if(cal_smallest_dist_to_one_edge_in_a_cell(which_cell, A, dir, dist, which_edge, which_samp,
				except_edge))
				flag=true;
		}
		if((quadmesh->quadcells[cellid]->x_start_coord+quadmesh->xinterval<quadmesh->xend-1.e-6)
			&&cellid>=quadmesh->XDIM-1)
		{
			/*  lower right */
			which_cell=cellid-quadmesh->XDIM+2;
			if(cal_smallest_dist_to_one_edge_in_a_cell(which_cell, A, dir, dist, which_edge, which_samp,
				except_edge))
				flag=true;
		}
	}

	return flag;
}


/*  
    Search the closest intersection in the boundary intersection list
*/

bool find_closest_intersection(double p[2], double &dist, int &which_intersect, 
							   icVector2 goDir, double cosang)
{
	int i;
	icVector2 dis;
	Intersection *intersect;
	

	for(i=0;i<nboundaryintersections;i++)
	{
		if(boundaryintersections[i]==which_intersect) 
			continue;

		intersect=streetnet->nodelist->intersects[boundaryintersections[i]];
		if(intersect->deleted)
			continue;

		dis.entry[0]=p[0]-intersect->gpos[0];
		dis.entry[1]=p[1]-intersect->gpos[1];

		double tempDist=length(dis);

		normalize(dis);

		if(dot(dis, goDir)<cosang)
			continue;

		if(length(dis)<1.e-2)
		{
			dist=length(dis);
			which_intersect=boundaryintersections[i];
			return true;
		}
	}

	return false;
}


/*
*/
bool find_closest_intersection(double p[2], double &dist, 
							   int &which_intersect, int &which_cell)
{
	int i;
	icVector2 dis;
	Intersection *intersect;

	for(i=0;i<nboundaryintersections;i++)
	{
		if(boundaryintersections[i]==which_intersect)
			continue;

		if(boundaryintersections[i]<0||
			boundaryintersections[i]>=streetnet->nodelist->nelems)
			continue;

		intersect=streetnet->nodelist->intersects[boundaryintersections[i]];
		dis.entry[0]=p[0]-intersect->gpos[0];
		dis.entry[1]=p[1]-intersect->gpos[1];

		//if(length(dis)<1.e-2)
		if(length(dis)<quadmesh->xinterval/4.)
		{
			dist=length(dis);
			which_intersect=boundaryintersections[i];
			which_cell=get_cellID_givencoords(intersect->gpos[0], intersect->gpos[1]);
			return true;
		}
	}

	return false;
}


/*
    compute the intersection that is closest to the point "p" on the edge "edgeid"
	inside a specified cell
*/
bool cal_smallest_intersect_in_one_cell(int cellid, double p1[2], double p2[2], int except_edge, 
										int &which_edge, int &which_samp,
										  double intersect[2])
{
	int i;
	QuadCell *face=quadmesh->quadcells[cellid];
	for(i=0;i<face->nstreetgraphedges;i++)
	{
		if(streetnet->edgelist->edges[face->streetgraphedgelist[i]]->cancel)
			continue;

		if(except_edge==face->streetgraphedgelist[i])
			continue;

		if(cal_intersect_at_graph_edge(face->streetgraphedgelist[i],
			intersect, which_samp, p1, p2))
		{
			which_edge=face->streetgraphedgelist[i];
			return true;
		}
	}

	return false;
}


/*
    This routine resorts to the tracing to find out the possible intersection
*/

unsigned char trace_for_intersection(double start[2], icVector2 linedir, int &which_intersection,
									int &which_edge, int &which_samp, 
									double intersect[2], int &cell_intersect, 
									int except_edge, int except_intersect, int MaxIter)
{
	double distance_intersect, distance_edge;

	/*  we first consider the boundary intersections in the list  */
	normalize(linedir);

	if(find_closest_intersection(start, distance_intersect, which_intersection, cell_intersect))
	{
		return 1;
	}

	/*  second, we trace along the "dir" direction and see which intersection or edge it will 
	touch*/
	double pre_p[2]={start[0], start[1]};
	double cur_p[2]={start[0]+linedir.entry[0], start[1]+linedir.entry[1]};
	int i, j;
	int cur_cell = get_cellID_givencoords(start[0], start[1]);
	QuadCell *face;
	tenline_dir_global = linedir;

	icVector2 t_major[4];


	for(i=0;i<MaxIter;i++)
	{
		if(cur_cell <0 || cur_cell>=quadmesh->nfaces)
			return 0;

		/*   check the intersection with an edge in current cell if there are any  */
		if(cal_smallest_intersect_in_one_cell(cur_cell, pre_p, cur_p, except_edge, 
			which_edge, which_samp, intersect))
		{
			return 2;  /* intersect with an edge */
		}
		
		if(!is_in_cell(cur_cell, cur_p[0], cur_p[1]))
		{
			/*  search for next cell  */
			face = quadmesh->quadcells[cur_cell];

			for(j=0; j<4; j++)
				t_major[j]=quadmesh->quad_verts[face->verts[j]]->major;

			/*replace them with the line segment direction*/
			for(j=0; j<4; j++)
				quadmesh->quad_verts[face->verts[j]]->major=tenline_dir_global;

			int passvertornot = 0;
			double tp[2]={pre_p[0], pre_p[1]};
			get_next_cell_2(cur_cell, pre_p, cur_p, passvertornot, 0);

			if(passvertornot==0)
			{
				//cur_p[0]=pre_p[0];
				//cur_p[1]=pre_p[1];
				//pre_p[0]=tp[0];
				//pre_p[1]=tp[1];
				pre_p[0]=cur_p[0];
				pre_p[1]=cur_p[1];
			}

			//pre_p[0]=cur_p[0];
			//pre_p[1]=cur_p[1];
			cur_p[0]=start[0]+linedir.entry[0];
			cur_p[1]=start[1]+linedir.entry[1];

			/*store back the original vectors*/
			for(j=0; j<4; j++)
				quadmesh->quad_verts[face->verts[j]]->major=t_major[j];
		}
	}

	return 0;
}

/*
   We search only the boundary intersections 
   and the edges
*/
void merge_subgraph_to_streetnet_2(bool majormin, bool startorend)
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
		cur_place=major;
		infolist=majorintersectinfo;
	}
	else{
		cur_place=minor;
		infolist=minorintersectinfo;
	}

	for(i=cur_place->evenstreamlines->ntrajs-cur_place->nnewlines_inReg;
		i<cur_place->evenstreamlines->ntrajs;i++)
	{
		traj=cur_place->evenstreamlines->trajs[i];

		if(traj->closed)
			continue;

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
			cell_id=traj->linesegs[traj->nlinesegs-1].Triangle_ID;
		}

		//if(!quadmesh->quadcells[cell_id]->OnBoundary)
		//	return;

		/*  obtain the edge that is associated with the intersection  */
		except_edge=streetnet->nodelist->intersects[cur_intersect_id]->adj_edges[0];

		//for(j=0;j<2;j++)
		//{
		//	if(j==0){
		//	A[0]+=1.e-7;
		//	A[1]+=1.e-7;
		//	cell_id=get_cellID_givencoords(A[0],A[1]);
		//	}
		//	else
		//	{
		//	A[0]-=2.e-7;
		//	A[1]-=2.e-7;
		//	cell_id=get_cellID_givencoords(A[0],A[1]);
		//	}
		
		ang=atan2(dir.entry[1], dir.entry[0]);
		if(ang<0) ang+=(2.*M_PI);


		/*  compute the distance between A and other nearby intersections  */
		if(cal_smallest_dist_to_intersects(cell_id, A, ang, dir, dist_intersect, which_intersect,
			cur_intersect_id))
		{
			//dist_intersect=temp_dist;
			if(which_intersect!=cur_intersect_id)
			closetointersect=true;

			//if(which_intersect>=streetnet->nodelist->nelems)
			//{
			//	int test=0;
			//}
		}

			int result=trace_for_intersection(A, dir, which_intersect, which_edge, which_samp, 
				the_intersect, which_cell, except_edge, cur_intersect_id, 5);

			//if(which_intersect>=streetnet->nodelist->nelems)
			//{
			//	int test=0;
			//}

			if(result==1)
			{
				/*   merge with existing intersection  */
					/*  extend to the intersection  */
					/*  merge to the intersection "which_intersect"  */
					/*  update "cur_intersect_id"  */

					streetnet->nodelist->intersects[cur_intersect_id]->deleted=true;
					StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];

					Trajectory *temp_traj = new Trajectory(-1);
					get_linesegs_anytwopts(the_intersect, which_cell,
						streetnet->nodelist->intersects[cur_intersect_id]->gpos, 
						streetnet->nodelist->intersects[cur_intersect_id]->cellid,
						temp_traj, 0, 10);

					int edge_cell;

					if(edge->node_index1==cur_intersect_id)
					{
						edge->node_index1=which_intersect;
						/*update the sampling list*/
						edge->inter_pts[0]->x=streetnet->nodelist->intersects[which_intersect]->gpos[0];
						edge->inter_pts[0]->y=streetnet->nodelist->intersects[which_intersect]->gpos[1];
						edge_cell=edge->inter_pts[0]->cellid=streetnet->nodelist->intersects[which_intersect]->cellid;
					}
					else
					{
						edge->node_index2=which_intersect;
						/*update the sampling list*/
						edge->inter_pts[edge->ninter_pts-1]->x
							=streetnet->nodelist->intersects[which_intersect]->gpos[0];
						edge->inter_pts[edge->ninter_pts-1]->y
							=streetnet->nodelist->intersects[which_intersect]->gpos[1];
						edge_cell=edge->inter_pts[edge->ninter_pts-1]->cellid
							=streetnet->nodelist->intersects[which_intersect]->cellid;
					}

					if(!is_repeated_elem(quadmesh->quadcells[edge_cell]->streetgraphedgelist,
						edge->index, quadmesh->quadcells[edge_cell]->nstreetgraphedges))
					{
						add_to_edgelist_one_cell(edge_cell, edge->index);
					}

					int pre_cell=-1;
					int k;
					for(k=0;k<temp_traj->nlinesegs;k++)
					{
						if(pre_cell==temp_traj->linesegs[k].Triangle_ID) 
							continue;
						pre_cell = temp_traj->linesegs[k].Triangle_ID;
						QuadCell *face=quadmesh->quadcells[pre_cell];
						if(!is_repeated_elem(face->streetgraphedgelist, pre_cell,
							face->nstreetgraphedges))
						{
							add_to_edgelist_one_cell(pre_cell, edge->index);
						}
					}

					streetnet->nodelist->intersects[which_intersect]->endpt=false;

					/*  possible bug:  remember to add this edge to the "which_intersect"  
							1/19/2008
					*/
					add_edge_to_intersectnode(which_intersect, except_edge);

					delete temp_traj;
			}


			else if(result==2)
			{
				/*  compute the distance to the obtained intersection  */
				dir.entry[0]=streetnet->nodelist->intersects[cur_intersect_id]->gpos[0]-
					the_intersect[0];
				dir.entry[1]=streetnet->nodelist->intersects[cur_intersect_id]->gpos[1]-
					the_intersect[1];
				dist_edge=length(dir);

				if(closetointersect&&dist_edge>dist_intersect)
				{
					streetnet->nodelist->intersects[cur_intersect_id]->deleted=true;

					StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];

					Trajectory *temp_traj = new Trajectory(-1);
					get_linesegs_anytwopts(streetnet->nodelist->intersects[which_intersect]->gpos, 
						streetnet->nodelist->intersects[which_intersect]->cellid,
						streetnet->nodelist->intersects[cur_intersect_id]->gpos, 
						streetnet->nodelist->intersects[cur_intersect_id]->cellid,
						temp_traj, 0, 5);

					int edge_cell;

					if(edge->node_index1==cur_intersect_id)
					{
						edge->node_index1=which_intersect;
						/*update the sampling list*/
						edge->inter_pts[0]->x=streetnet->nodelist->intersects[which_intersect]->gpos[0];
						edge->inter_pts[0]->y=streetnet->nodelist->intersects[which_intersect]->gpos[1];
						edge_cell=edge->inter_pts[0]->cellid=streetnet->nodelist->intersects[which_intersect]->cellid;
					}
					else
					{
						edge->node_index2=which_intersect;
						/*update the sampling list*/
						edge->inter_pts[edge->ninter_pts-1]->x
							=streetnet->nodelist->intersects[which_intersect]->gpos[0];
						edge->inter_pts[edge->ninter_pts-1]->y
							=streetnet->nodelist->intersects[which_intersect]->gpos[1];
						edge_cell=edge->inter_pts[edge->ninter_pts-1]->cellid
							=streetnet->nodelist->intersects[which_intersect]->cellid;
					}

					if(!is_repeated_elem(quadmesh->quadcells[edge_cell]->streetgraphedgelist,
						edge->index, quadmesh->quadcells[edge_cell]->nstreetgraphedges))
					{
						add_to_edgelist_one_cell(edge_cell, edge->index);
					}

					int pre_cell=-1;
					int k;
					for(k=0;k<temp_traj->nlinesegs;k++)
					{
						if(pre_cell==temp_traj->linesegs[k].Triangle_ID) continue;
						pre_cell = temp_traj->linesegs[k].Triangle_ID;
						QuadCell *face=quadmesh->quadcells[pre_cell];
						if(!is_repeated_elem(face->streetgraphedgelist, pre_cell,
							face->nstreetgraphedges))
						{
							add_to_edgelist_one_cell(pre_cell, edge->index);
						}
					}

					streetnet->nodelist->intersects[which_intersect]->endpt=false;
					delete temp_traj;
				}
				else
				{
					/*  extend to the edge  */
					
					streetnet->nodelist->intersects[cur_intersect_id]->endpt=false;

					/*   use the first point of the line segment   */
					//streetnet->nodelist->intersects[cur_intersect_id]->gpos[0]=
					//	streetnet->edgelist->edges[which_edge]->inter_pts[which_samp]->x;
					//streetnet->nodelist->intersects[cur_intersect_id]->gpos[1]=
					//	streetnet->edgelist->edges[which_edge]->inter_pts[which_samp]->y;
					//StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];

					StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];
					Trajectory *temp_traj = new Trajectory(-1);
					get_linesegs_anytwopts(the_intersect, get_cellID_givencoords(the_intersect[0], the_intersect[1]),
						streetnet->nodelist->intersects[cur_intersect_id]->gpos, 
						streetnet->nodelist->intersects[cur_intersect_id]->cellid,
						temp_traj, 0, 5);
				
					streetnet->nodelist->intersects[cur_intersect_id]->gpos[0]=
						the_intersect[0];
					streetnet->nodelist->intersects[cur_intersect_id]->gpos[1]=
						the_intersect[1];


					if(edge->node_index1==cur_intersect_id)
					{
						edge->inter_pts[0]->x=
							streetnet->nodelist->intersects[cur_intersect_id]->gpos[0];
						edge->inter_pts[0]->y=
							streetnet->nodelist->intersects[cur_intersect_id]->gpos[1];

						edge->inter_pts[0]->cellid=get_cellID_givencoords(
							edge->inter_pts[0]->x, edge->inter_pts[0]->y);
					}
					else
					{
						edge->inter_pts[edge->ninter_pts-1]->x=
							streetnet->nodelist->intersects[cur_intersect_id]->gpos[0];
						edge->inter_pts[edge->ninter_pts-1]->y=
							streetnet->nodelist->intersects[cur_intersect_id]->gpos[1];
						edge->inter_pts[edge->ninter_pts-1]->cellid=get_cellID_givencoords(
							edge->inter_pts[edge->ninter_pts-1]->x, 
							edge->inter_pts[edge->ninter_pts-1]->y);
					}

					/*  update the cells that this edge pass!!  */
					int k;
					int pre_cell=-1;
					for(k=0;k<edge->ninter_pts;k++)
					{
						if(pre_cell == edge->inter_pts[k]->cellid) continue;
						pre_cell=edge->inter_pts[k]->cellid;
						if(pre_cell<0||pre_cell>=quadmesh->nfaces)
						{
							pre_cell=edge->inter_pts[k]->cellid=
								get_cellID_givencoords(edge->inter_pts[k]->x,
								edge->inter_pts[k]->y);
							if(pre_cell<0||pre_cell>=quadmesh->nfaces)
								continue;
						}
						QuadCell *face=quadmesh->quadcells[edge->inter_pts[k]->cellid];
						if(!is_repeated_elem(face->streetgraphedgelist, pre_cell,
							face->nstreetgraphedges))
						{
							add_to_edgelist_one_cell(pre_cell, edge->index);
						}
					}
					pre_cell=-1;
					for(k=0;k<temp_traj->nlinesegs;k++)
					{
						if(pre_cell==temp_traj->linesegs[k].Triangle_ID) continue;
						pre_cell = temp_traj->linesegs[k].Triangle_ID;
						QuadCell *face=quadmesh->quadcells[pre_cell];
						if(!is_repeated_elem(face->streetgraphedgelist, pre_cell,
							face->nstreetgraphedges))
						{
							add_to_edgelist_one_cell(pre_cell, edge->index);
						}
					}

					/*  update the original edge and   
						generate one new edges in the street network (graph)  
					*/
					edge=streetnet->edgelist->edges[which_edge];

					StreetGraphEdge *newedge=(StreetGraphEdge*)malloc(sizeof(StreetGraphEdge));
					newedge->node_index1=cur_intersect_id;
					newedge->node_index2=edge->node_index2;
					newedge->cancel = newedge->visited = false;
					newedge->inter_pts=(Point**)malloc(sizeof(Point*)*(edge->ninter_pts-which_samp));
					int j;
					for(j=which_samp;j<edge->ninter_pts;j++)
					{
						if(j==which_samp)
						{
							newedge->inter_pts[j-which_samp]=(Point*)malloc(sizeof(Point));
							newedge->inter_pts[j-which_samp]->x=the_intersect[0];
							newedge->inter_pts[j-which_samp]->y=the_intersect[1];
							newedge->inter_pts[j-which_samp]->cellid=get_cellID_givencoords(
								the_intersect[0], the_intersect[1]);
							continue;
						}
						newedge->inter_pts[j-which_samp]=(Point*)malloc(sizeof(Point));
						newedge->inter_pts[j-which_samp]->x=edge->inter_pts[j]->x;
						newedge->inter_pts[j-which_samp]->y=edge->inter_pts[j]->y;
						newedge->inter_pts[j-which_samp]->cellid=edge->inter_pts[j]->cellid;

						if(j!=which_samp)
							free(edge->inter_pts[j]);
					}
					newedge->ninter_pts=edge->ninter_pts-which_samp;
					newedge->roadtype=edge->roadtype;
					
					streetnet->edgelist->append(newedge);

					/*  update the old edge  */
					edge->node_index2=cur_intersect_id;

					Point **temp=edge->inter_pts;
					edge->inter_pts=(Point**)malloc(sizeof(Point*)*(which_samp+2));
					for(j=0;j<=which_samp;j++)
						edge->inter_pts[j]=temp[j];
					edge->inter_pts[j]=(Point*)malloc(sizeof(Point));
					edge->inter_pts[j]->x=the_intersect[0];
					edge->inter_pts[j]->y=the_intersect[1];

					edge->ninter_pts=which_samp+2;
					free(temp);

					/*   should we remove the edge from those cells after "the_intersect"  */

					/*  update the edge list of the intersection "cur_intersect_id"  */
					add_edge_to_intersectnode(cur_intersect_id, edge->index);
					add_edge_to_intersectnode(cur_intersect_id, newedge->index);

					delete temp_traj;
				}
			}
		//}
	}
}


/*  1/20/2008   */
void merge_subgraph_to_streetnet_3(bool majormin, bool startorend)
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
		cur_place=major;
		infolist=majorintersectinfo;
	}
	else{
		cur_place=minor;
		infolist=minorintersectinfo;
	}

	for(i=cur_place->evenstreamlines->ntrajs-cur_place->nnewlines_inReg;
		i<cur_place->evenstreamlines->ntrajs;i++)
	{
		traj=cur_place->evenstreamlines->trajs[i];

		if(traj->closed)
			continue;

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

				dir.entry[0]=traj->linesegs[traj->nlinesegs-1].gend[0]
					-traj->linesegs[traj->nlinesegs-1].gstart[0];
				dir.entry[1]=traj->linesegs[traj->nlinesegs-1].gend[1]
					-traj->linesegs[traj->nlinesegs-1].gstart[1];
			normalize(dir);
			cell_id=traj->linesegs[traj->nlinesegs-1].Triangle_ID;
		}

		/*  obtain the edge that is associated with the intersection  */
		except_edge=streetnet->nodelist->intersects[cur_intersect_id]->adj_edges[0];

		ang=atan2(dir.entry[1], dir.entry[0]);
		if(ang<0) ang+=(2.*M_PI);


		/*  compute the distance between A and other nearby intersections  */
		if(cal_smallest_dist_to_intersects(cell_id, A, ang, dir, dist_intersect, which_intersect,
			cur_intersect_id))
		{
			//dist_intersect=temp_dist;
			if(which_intersect!=cur_intersect_id)
				closetointersect=true;

			/*   merge with existing intersection  */
				/*  extend to the intersection  */
				/*  merge to the intersection "which_intersect"  */
				/*  update "cur_intersect_id"  */


				StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];

				the_intersect[0]=streetnet->nodelist->intersects[which_intersect]->gpos[0];
				the_intersect[1]=streetnet->nodelist->intersects[which_intersect]->gpos[1];
				which_cell=streetnet->nodelist->intersects[which_intersect]->cellid;

				Trajectory *temp_traj = new Trajectory(-1);
				get_linesegs_anytwopts(
					streetnet->nodelist->intersects[cur_intersect_id]->gpos, 
					streetnet->nodelist->intersects[cur_intersect_id]->cellid,
					the_intersect, which_cell,
					temp_traj, 0, 10);

				/* we add a new edge  */
				if(has_edge_between_minRoad(cur_intersect_id, which_intersect))
				{
					delete temp_traj;
					continue;
				}

				StreetGraphEdge *anewedge=(StreetGraphEdge *)malloc(sizeof(StreetGraphEdge));
				anewedge->node_index1=cur_intersect_id;
				anewedge->node_index2=which_intersect;
				anewedge->roadtype=edge->roadtype;
				anewedge->cancel=anewedge->visited=false;
				anewedge->inter_pts=(Point**)malloc(sizeof(Point*)*(temp_traj->nlinesegs+1));
				streetnet->edgelist->append(anewedge);


				/*  possible bug:  remember to add this edge to the "which_intersect"  
						1/19/2008
				*/
				int m;
				int pre_cell=-1;
				for(m=0;m<temp_traj->nlinesegs;m++)
				{
					anewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
					anewedge->inter_pts[m]->x=temp_traj->linesegs[m].gstart[0];
					anewedge->inter_pts[m]->y=temp_traj->linesegs[m].gstart[1];
					anewedge->inter_pts[m]->cellid=get_cellID_givencoords(anewedge->inter_pts[m]->x,
						anewedge->inter_pts[m]->y);
					
					if(pre_cell!=anewedge->inter_pts[m]->cellid)
					{
						pre_cell=anewedge->inter_pts[m]->cellid;
						add_to_edgelist_one_cell(anewedge->inter_pts[m]->cellid, 
							anewedge->index);
					}
				}

				anewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
				anewedge->inter_pts[m]->x=temp_traj->linesegs[temp_traj->nlinesegs-1].gend[0];
				anewedge->inter_pts[m]->y=temp_traj->linesegs[temp_traj->nlinesegs-1].gend[1];
				anewedge->inter_pts[m]->cellid=get_cellID_givencoords(anewedge->inter_pts[m]->x,
					anewedge->inter_pts[m]->y);
				anewedge->ninter_pts=temp_traj->nlinesegs+1;
				if(pre_cell!=anewedge->inter_pts[m]->cellid)
				{
					pre_cell=anewedge->inter_pts[m]->cellid;
					add_to_edgelist_one_cell(anewedge->inter_pts[m]->cellid, 
						anewedge->index);
				}

				streetnet->nodelist->intersects[which_intersect]->endpt=false;
				streetnet->nodelist->intersects[cur_intersect_id]->endpt=false;

				/*   add to intersection list  */
				add_edge_to_intersectnode(which_intersect, anewedge->index);
				add_edge_to_intersectnode(cur_intersect_id, anewedge->index);

				delete temp_traj;
		}
		else
		{

			int result=trace_for_intersection(A, dir, which_intersect, which_edge, which_samp, 
				the_intersect, which_cell, except_edge, cur_intersect_id, 5);

			if(result==1)  /*  merge with existing intersections  */
			{
				/*   merge with existing intersection  */
					/*  extend to the intersection  */
					/*  merge to the intersection "which_intersect"  */
					/*  update "cur_intersect_id"  */

				StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];

				the_intersect[0]=streetnet->nodelist->intersects[which_intersect]->gpos[0];
				the_intersect[1]=streetnet->nodelist->intersects[which_intersect]->gpos[1];
				which_cell=streetnet->nodelist->intersects[which_intersect]->cellid;

				Trajectory *temp_traj = new Trajectory(-1);
				get_linesegs_anytwopts(
					streetnet->nodelist->intersects[cur_intersect_id]->gpos, 
					streetnet->nodelist->intersects[cur_intersect_id]->cellid,
					the_intersect, which_cell,
					temp_traj, 0, 10);

				/* we add a new edge  */
				StreetGraphEdge *anewedge=(StreetGraphEdge *)malloc(sizeof(StreetGraphEdge));
				anewedge->node_index1=cur_intersect_id;
				anewedge->node_index2=which_intersect;
				anewedge->roadtype=edge->roadtype;
				anewedge->cancel=anewedge->visited=false;
				anewedge->inter_pts=(Point**)malloc(sizeof(Point*)*(temp_traj->nlinesegs+1));
				streetnet->edgelist->append(anewedge);


				/*  possible bug:  remember to add this edge to the "which_intersect"  
						1/19/2008
				*/
				int m;
				int pre_cell=-1;
				for(m=0;m<temp_traj->nlinesegs;m++)
				{
					anewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
					anewedge->inter_pts[m]->x=temp_traj->linesegs[m].gstart[0];
					anewedge->inter_pts[m]->y=temp_traj->linesegs[m].gstart[1];
					anewedge->inter_pts[m]->cellid=get_cellID_givencoords(anewedge->inter_pts[m]->x,
						anewedge->inter_pts[m]->y);
					
					if(pre_cell!=anewedge->inter_pts[m]->cellid)
					{
						pre_cell=anewedge->inter_pts[m]->cellid;
						add_to_edgelist_one_cell(anewedge->inter_pts[m]->cellid, 
							anewedge->index);
					}
				}

				anewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
				anewedge->inter_pts[m]->x=temp_traj->linesegs[temp_traj->nlinesegs-1].gend[0];
				anewedge->inter_pts[m]->y=temp_traj->linesegs[temp_traj->nlinesegs-1].gend[1];
				anewedge->inter_pts[m]->cellid=get_cellID_givencoords(anewedge->inter_pts[m]->x,
					anewedge->inter_pts[m]->y);
				anewedge->ninter_pts=temp_traj->nlinesegs+1;
				if(pre_cell!=anewedge->inter_pts[m]->cellid)
				{
					pre_cell=anewedge->inter_pts[m]->cellid;
					add_to_edgelist_one_cell(anewedge->inter_pts[m]->cellid, 
						anewedge->index);
				}

				streetnet->nodelist->intersects[which_intersect]->endpt=false;
				streetnet->nodelist->intersects[cur_intersect_id]->endpt=false;

				/*   add to intersection list  */
				add_edge_to_intersectnode(which_intersect, anewedge->index);
				add_edge_to_intersectnode(cur_intersect_id, anewedge->index);

				delete temp_traj;

			}


			//else if(result==2)
			//{
			//	/*  compute the distance to the obtained point on the edge  */
			//	//dir.entry[0]=streetnet->nodelist->intersects[cur_intersect_id]->gpos[0]-
			//	//	the_intersect[0];
			//	//dir.entry[1]=streetnet->nodelist->intersects[cur_intersect_id]->gpos[1]-
			//	//	the_intersect[1];
			//	//dist_edge=length(dir);

			//	//if(closetointersect&&dist_edge>dist_intersect)
			//	//{
			//	//	streetnet->nodelist->intersects[cur_intersect_id]->deleted=true;

			//	//	StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];

			//	//	Trajectory *temp_traj = new Trajectory(-1);
			//	//	get_linesegs_anytwopts(streetnet->nodelist->intersects[which_intersect]->gpos, 
			//	//		streetnet->nodelist->intersects[which_intersect]->cellid,
			//	//		streetnet->nodelist->intersects[cur_intersect_id]->gpos, 
			//	//		streetnet->nodelist->intersects[cur_intersect_id]->cellid,
			//	//		temp_traj, 0, 5);

			//	//	int edge_cell;

			//	//	if(edge->node_index1==cur_intersect_id)
			//	//	{
			//	//		edge->node_index1=which_intersect;
			//	//		/*update the sampling list*/
			//	//		edge->inter_pts[0]->x=streetnet->nodelist->intersects[which_intersect]->gpos[0];
			//	//		edge->inter_pts[0]->y=streetnet->nodelist->intersects[which_intersect]->gpos[1];
			//	//		edge_cell=edge->inter_pts[0]->cellid=streetnet->nodelist->intersects[which_intersect]->cellid;
			//	//	}
			//	//	else
			//	//	{
			//	//		edge->node_index2=which_intersect;
			//	//		/*update the sampling list*/
			//	//		edge->inter_pts[edge->ninter_pts-1]->x
			//	//			=streetnet->nodelist->intersects[which_intersect]->gpos[0];
			//	//		edge->inter_pts[edge->ninter_pts-1]->y
			//	//			=streetnet->nodelist->intersects[which_intersect]->gpos[1];
			//	//		edge_cell=edge->inter_pts[edge->ninter_pts-1]->cellid
			//	//			=streetnet->nodelist->intersects[which_intersect]->cellid;
			//	//	}

			//	//	if(!is_repeated_elem(quadmesh->quadcells[edge_cell]->streetgraphedgelist,
			//	//		edge->index, quadmesh->quadcells[edge_cell]->nstreetgraphedges))
			//	//	{
			//	//		add_to_edgelist_one_cell(edge_cell, edge->index);
			//	//	}

			//	//	int pre_cell=-1;
			//	//	int k;
			//	//	for(k=0;k<temp_traj->nlinesegs;k++)
			//	//	{
			//	//		if(pre_cell==temp_traj->linesegs[k].Triangle_ID) continue;
			//	//		pre_cell = temp_traj->linesegs[k].Triangle_ID;
			//	//		QuadCell *face=quadmesh->quadcells[pre_cell];
			//	//		if(!is_repeated_elem(face->streetgraphedgelist, pre_cell,
			//	//			face->nstreetgraphedges))
			//	//		{
			//	//			add_to_edgelist_one_cell(pre_cell, edge->index);
			//	//		}
			//	//	}

			//	//	streetnet->nodelist->intersects[which_intersect]->endpt=false;
			//		/*  possible bug:  remember to add this edge to the "which_intersect"  
			//				1/19/2008
			//		*/
			//		//add_edge_to_intersectnode(which_intersect, except_edge);
			//	//	delete temp_traj;
			//	//}
			//	//else
			//	//{
			//		/*  extend to the edge  */
			//		
			//		streetnet->nodelist->intersects[cur_intersect_id]->endpt=false;
			//		streetnet->nodelist->intersects[cur_intersect_id]->deleted=false;

			//		/*   use the first point of the line segment   */

			//		StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];
			//		Trajectory *temp_traj = new Trajectory(-1);
			//		get_linesegs_anytwopts(the_intersect, get_cellID_givencoords(the_intersect[0], the_intersect[1]),
			//			streetnet->nodelist->intersects[cur_intersect_id]->gpos, 
			//			streetnet->nodelist->intersects[cur_intersect_id]->cellid,
			//			temp_traj, 0, 5);

			//		if(which_samp==0)
			//		{
			//			the_intersect[0]=streetnet->edgelist->edges[which_edge]->inter_pts[0]->x;
			//			the_intersect[1]=streetnet->edgelist->edges[which_edge]->inter_pts[0]->y;
			//		}
			//		else if(which_samp==streetnet->edgelist->edges[which_edge]->ninter_pts-1)
			//		{
			//			int tempID=streetnet->edgelist->edges[which_edge]->ninter_pts-1;
			//			the_intersect[0]=streetnet->edgelist->edges[which_edge]->inter_pts[tempID]->x;
			//			the_intersect[1]=streetnet->edgelist->edges[which_edge]->inter_pts[tempID]->y;
			//		}
			//	
			//		streetnet->nodelist->intersects[cur_intersect_id]->gpos[0]=
			//			the_intersect[0];
			//		streetnet->nodelist->intersects[cur_intersect_id]->gpos[1]=
			//			the_intersect[1];


			//		if(edge->node_index1==cur_intersect_id)
			//		{
			//			edge->inter_pts[0]->x=
			//				streetnet->nodelist->intersects[cur_intersect_id]->gpos[0];
			//			edge->inter_pts[0]->y=
			//				streetnet->nodelist->intersects[cur_intersect_id]->gpos[1];

			//			edge->inter_pts[0]->cellid=get_cellID_givencoords(
			//				edge->inter_pts[0]->x, edge->inter_pts[0]->y);
			//		}
			//		else
			//		{
			//			edge->inter_pts[edge->ninter_pts-1]->x=
			//				streetnet->nodelist->intersects[cur_intersect_id]->gpos[0];
			//			edge->inter_pts[edge->ninter_pts-1]->y=
			//				streetnet->nodelist->intersects[cur_intersect_id]->gpos[1];
			//			edge->inter_pts[edge->ninter_pts-1]->cellid=get_cellID_givencoords(
			//				edge->inter_pts[edge->ninter_pts-1]->x, 
			//				edge->inter_pts[edge->ninter_pts-1]->y);
			//		}

			//		/*  update the cells that this edge passes!!  */
			//		int k;
			//		int pre_cell=-1;
			//		for(k=0;k<edge->ninter_pts;k++)
			//		{
			//			if(pre_cell == edge->inter_pts[k]->cellid) continue;
			//			pre_cell=edge->inter_pts[k]->cellid;
			//			if(pre_cell<0||pre_cell>=quadmesh->nfaces)
			//			{
			//				pre_cell=edge->inter_pts[k]->cellid=
			//					get_cellID_givencoords(edge->inter_pts[k]->x,
			//					edge->inter_pts[k]->y);
			//				if(pre_cell<0||pre_cell>=quadmesh->nfaces)
			//					continue;
			//			}
			//			QuadCell *face=quadmesh->quadcells[edge->inter_pts[k]->cellid];
			//			if(!is_repeated_elem(face->streetgraphedgelist, pre_cell,
			//				face->nstreetgraphedges))
			//			{
			//				add_to_edgelist_one_cell(pre_cell, edge->index);
			//			}
			//		}
			//		pre_cell=-1;
			//		for(k=0;k<temp_traj->nlinesegs;k++)
			//		{
			//			if(pre_cell==temp_traj->linesegs[k].Triangle_ID) continue;
			//			pre_cell = temp_traj->linesegs[k].Triangle_ID;
			//			QuadCell *face=quadmesh->quadcells[pre_cell];
			//			if(!is_repeated_elem(face->streetgraphedgelist, pre_cell,
			//				face->nstreetgraphedges))
			//			{
			//				add_to_edgelist_one_cell(pre_cell, edge->index);
			//			}
			//		}

			//		/*  update the original edge and   
			//			generate one new edges in the street network (graph)  
			//		*/
			//		edge=streetnet->edgelist->edges[which_edge];
			//		/*  delete the old edge  */
			//		edge->cancel=true;

			//		/*  create two new edges  */

			//		/* 1: */
			//		StreetGraphEdge *onenewedge=(StreetGraphEdge *)
			//			malloc(sizeof(StreetGraphEdge));
			//		onenewedge->cancel=
			//			onenewedge->visited=false;
			//		onenewedge->node_index1=cur_intersect_id;
			//		onenewedge->node_index2=edge->node_index1;
			//		onenewedge->roadtype = edge->roadtype;
			//		onenewedge->cancel=onenewedge->visited=false;

			//		streetnet->edgelist->append(onenewedge);

			//		double start[2], end[2];
			//		int start_cell, end_cell;
			//		start[0]=streetnet->nodelist->intersects[onenewedge->node_index1]->gpos[0];
			//		start[1]=streetnet->nodelist->intersects[onenewedge->node_index1]->gpos[1];
			//		start_cell=get_cellID_givencoords(start[0], start[1]);
			//		end[0]=streetnet->nodelist->intersects[onenewedge->node_index2]->gpos[0];
			//		end[1]=streetnet->nodelist->intersects[onenewedge->node_index2]->gpos[1];
			//		end_cell=get_cellID_givencoords(end[0], end[1]);
			//		Trajectory *tempTraj=new Trajectory(-1);
			//		get_linesegs_anytwopts(start,start_cell, end,end_cell, tempTraj, 0, 10);
			//		onenewedge->inter_pts=(Point **)malloc(sizeof(Point *)* (tempTraj->nlinesegs+1));

			//		int m;
			//		pre_cell=-1;
			//		for(m=0;m<tempTraj->nlinesegs;m++)
			//		{
			//			onenewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
			//			onenewedge->inter_pts[m]->x=tempTraj->linesegs[m].gstart[0];
			//			onenewedge->inter_pts[m]->y=tempTraj->linesegs[m].gstart[1];
			//			onenewedge->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge->inter_pts[m]->x,
			//				onenewedge->inter_pts[m]->y);
			//			
			//			if(pre_cell!=onenewedge->inter_pts[m]->cellid)
			//			{
			//				pre_cell=onenewedge->inter_pts[m]->cellid;
			//				add_to_edgelist_one_cell(onenewedge->inter_pts[m]->cellid, 
			//					onenewedge->index);
			//			}
			//		}

			//		onenewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
			//		onenewedge->inter_pts[m]->x=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[0];
			//		onenewedge->inter_pts[m]->y=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[1];
			//		onenewedge->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge->inter_pts[m]->x,
			//			onenewedge->inter_pts[m]->y);
			//		onenewedge->ninter_pts=tempTraj->nlinesegs+1;
			//		if(pre_cell!=onenewedge->inter_pts[m]->cellid)
			//		{
			//			pre_cell=onenewedge->inter_pts[m]->cellid;
			//			add_to_edgelist_one_cell(onenewedge->inter_pts[m]->cellid, 
			//				onenewedge->index);
			//		}
			//		delete tempTraj;

			//		/*  add to the adjacent edge lists of the two end intersections  */
			//		add_edge_to_intersectnode(onenewedge->node_index1, onenewedge->index);
			//		add_edge_to_intersectnode(onenewedge->node_index2, onenewedge->index);


			//		/* 2: */
			//		StreetGraphEdge *onenewedge2=(StreetGraphEdge *)
			//			malloc(sizeof(StreetGraphEdge));
			//		onenewedge2->cancel=
			//			onenewedge2->visited=false;
			//		onenewedge2->node_index1=cur_intersect_id;
			//		onenewedge2->node_index2=edge->node_index2;
			//		onenewedge2->roadtype = edge->roadtype;
			//		onenewedge2->cancel=onenewedge2->visited=false;

			//		streetnet->edgelist->append(onenewedge2);

			//		start[0]=streetnet->nodelist->intersects[onenewedge2->node_index1]->gpos[0];
			//		start[1]=streetnet->nodelist->intersects[onenewedge2->node_index1]->gpos[1];
			//		start_cell=get_cellID_givencoords(start[0], start[1]);
			//		end[0]=streetnet->nodelist->intersects[onenewedge2->node_index2]->gpos[0];
			//		end[1]=streetnet->nodelist->intersects[onenewedge2->node_index2]->gpos[1];
			//		end_cell=get_cellID_givencoords(end[0], end[1]);
			//		tempTraj=new Trajectory(-1);
			//		get_linesegs_anytwopts(start,start_cell, end,end_cell, tempTraj, 0, 10);
			//		onenewedge2->inter_pts=(Point **)malloc(sizeof(Point *)* (tempTraj->nlinesegs+1));

			//		pre_cell=-1;
			//		for(m=0;m<tempTraj->nlinesegs;m++)
			//		{
			//			onenewedge2->inter_pts[m]=(Point*)malloc(sizeof(Point));
			//			onenewedge2->inter_pts[m]->x=tempTraj->linesegs[m].gstart[0];
			//			onenewedge2->inter_pts[m]->y=tempTraj->linesegs[m].gstart[1];
			//			onenewedge2->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge2->inter_pts[m]->x,
			//				onenewedge2->inter_pts[m]->y);
			//			
			//			if(pre_cell!=onenewedge2->inter_pts[m]->cellid)
			//			{
			//				pre_cell=onenewedge2->inter_pts[m]->cellid;
			//				add_to_edgelist_one_cell(onenewedge2->inter_pts[m]->cellid, 
			//					onenewedge2->index);
			//			}
			//		}

			//		onenewedge2->inter_pts[m]=(Point*)malloc(sizeof(Point));
			//		onenewedge2->inter_pts[m]->x=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[0];
			//		onenewedge2->inter_pts[m]->y=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[1];
			//		onenewedge2->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge2->inter_pts[m]->x,
			//			onenewedge2->inter_pts[m]->y);
			//		onenewedge2->ninter_pts=tempTraj->nlinesegs+1;
			//		if(pre_cell!=onenewedge2->inter_pts[m]->cellid)
			//		{
			//			pre_cell=onenewedge2->inter_pts[m]->cellid;
			//			add_to_edgelist_one_cell(onenewedge2->inter_pts[m]->cellid, 
			//				onenewedge2->index);
			//		}
			//		delete tempTraj;

			//		/*  add to the adjacent edge lists of the two end intersections  */
			//		add_edge_to_intersectnode(onenewedge2->node_index1, onenewedge2->index);
			//		add_edge_to_intersectnode(onenewedge2->node_index2, onenewedge2->index);
			//		
			//	//}
			//}
			else
			{
				streetnet->nodelist->intersects[cur_intersect_id]->deleted=true;
				streetnet->edgelist->edges[except_edge]->cancel=true;
			}
		}
	}
}




void merge_subgraph_to_streetnet_4(bool majormin, bool startorend)
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
		cur_place=major;
		infolist=majorintersectinfo;
	}
	else{
		cur_place=minor;
		infolist=minorintersectinfo;
	}

	for(i=cur_place->evenstreamlines->ntrajs-cur_place->nnewlines_inReg;
		i<cur_place->evenstreamlines->ntrajs;i++)
	{
		traj=cur_place->evenstreamlines->trajs[i];

		if(traj->closed)
			continue;

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

				dir.entry[0]=traj->linesegs[traj->nlinesegs-1].gend[0]
					-traj->linesegs[traj->nlinesegs-1].gstart[0];
				dir.entry[1]=traj->linesegs[traj->nlinesegs-1].gend[1]
					-traj->linesegs[traj->nlinesegs-1].gstart[1];
			normalize(dir);
			cell_id=traj->linesegs[traj->nlinesegs-1].Triangle_ID;
		}

		/*  obtain the edge that is associated with the intersection  */
		except_edge=streetnet->nodelist->intersects[cur_intersect_id]->adj_edges[0];


		int result=trace_comp_nearby_intersect(A, cell_id, dir, 
								cur_intersect_id, except_edge, 
								which_intersect, the_intersect,
								which_edge, which_samp, 5);

		/*  compute the distance between A and other nearby intersections  */
		if(result==1)
		{
			/*   merge with existing intersection  */
				/*  extend to the intersection  */
				/*  merge to the intersection "which_intersect"  */
				/*  update "cur_intersect_id"  */


				StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];

				the_intersect[0]=streetnet->nodelist->intersects[which_intersect]->gpos[0];
				the_intersect[1]=streetnet->nodelist->intersects[which_intersect]->gpos[1];
				which_cell=streetnet->nodelist->intersects[which_intersect]->cellid;

				Trajectory *temp_traj = new Trajectory(-1);
				get_linesegs_anytwopts(
					streetnet->nodelist->intersects[cur_intersect_id]->gpos, 
					streetnet->nodelist->intersects[cur_intersect_id]->cellid,
					the_intersect, which_cell,
					temp_traj, 0, 10);

				/* we add a new edge  */
				StreetGraphEdge *anewedge=(StreetGraphEdge *)malloc(sizeof(StreetGraphEdge));
				anewedge->node_index1=cur_intersect_id;
				anewedge->node_index2=which_intersect;
				anewedge->roadtype=edge->roadtype;
				anewedge->cancel=anewedge->visited=false;
				anewedge->inter_pts=(Point**)malloc(sizeof(Point*)*(temp_traj->nlinesegs+1));
				streetnet->edgelist->append(anewedge);


				/*  possible bug:  remember to add this edge to the "which_intersect"  
						1/19/2008
				*/
				int m;
				int pre_cell=-1;
				for(m=0;m<temp_traj->nlinesegs;m++)
				{
					anewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
					anewedge->inter_pts[m]->x=temp_traj->linesegs[m].gstart[0];
					anewedge->inter_pts[m]->y=temp_traj->linesegs[m].gstart[1];
					anewedge->inter_pts[m]->cellid=get_cellID_givencoords(anewedge->inter_pts[m]->x,
						anewedge->inter_pts[m]->y);
					
					if(pre_cell!=anewedge->inter_pts[m]->cellid)
					{
						pre_cell=anewedge->inter_pts[m]->cellid;
						add_to_edgelist_one_cell(anewedge->inter_pts[m]->cellid, 
							anewedge->index);
					}
				}

				anewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
				anewedge->inter_pts[m]->x=temp_traj->linesegs[temp_traj->nlinesegs-1].gend[0];
				anewedge->inter_pts[m]->y=temp_traj->linesegs[temp_traj->nlinesegs-1].gend[1];
				anewedge->inter_pts[m]->cellid=get_cellID_givencoords(anewedge->inter_pts[m]->x,
					anewedge->inter_pts[m]->y);
				anewedge->ninter_pts=temp_traj->nlinesegs+1;
				if(pre_cell!=anewedge->inter_pts[m]->cellid)
				{
					pre_cell=anewedge->inter_pts[m]->cellid;
					add_to_edgelist_one_cell(anewedge->inter_pts[m]->cellid, 
						anewedge->index);
				}

				streetnet->nodelist->intersects[which_intersect]->endpt=false;
				streetnet->nodelist->intersects[cur_intersect_id]->endpt=false;

				/*   add to intersection list  */
				add_edge_to_intersectnode(which_intersect, anewedge->index);
				add_edge_to_intersectnode(cur_intersect_id, anewedge->index);

				delete temp_traj;
		}

		else if(result == 2)
		{
			/*  compute the distance to the obtained point on the edge  */
			/*  extend to the edge  */
			
			streetnet->nodelist->intersects[cur_intersect_id]->endpt=false;
			streetnet->nodelist->intersects[cur_intersect_id]->deleted=false;

			/*   use the first point of the line segment   */

			StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];

			if(which_samp==0)
			{
				the_intersect[0]=edge->inter_pts[0]->x;
				the_intersect[1]=edge->inter_pts[0]->y;
			}
			else if(which_samp==edge->ninter_pts-1)
			{
				the_intersect[0]=edge->inter_pts[edge->ninter_pts-1]->x;
				the_intersect[1]=edge->inter_pts[edge->ninter_pts-1]->y;
			}

			Trajectory *temp_traj = new Trajectory(-1);
			get_linesegs_anytwopts(the_intersect, get_cellID_givencoords(the_intersect[0], the_intersect[1]),
				streetnet->nodelist->intersects[cur_intersect_id]->gpos, 
				streetnet->nodelist->intersects[cur_intersect_id]->cellid,
				temp_traj, 0, 5);
		
			streetnet->nodelist->intersects[cur_intersect_id]->gpos[0]=
				the_intersect[0];
			streetnet->nodelist->intersects[cur_intersect_id]->gpos[1]=
				the_intersect[1];


			if(edge->node_index1==cur_intersect_id)
			{
				edge->inter_pts[0]->x=
					streetnet->nodelist->intersects[cur_intersect_id]->gpos[0];
				edge->inter_pts[0]->y=
					streetnet->nodelist->intersects[cur_intersect_id]->gpos[1];

				edge->inter_pts[0]->cellid=get_cellID_givencoords(
					edge->inter_pts[0]->x, edge->inter_pts[0]->y);
			}
			else
			{
				edge->inter_pts[edge->ninter_pts-1]->x=
					streetnet->nodelist->intersects[cur_intersect_id]->gpos[0];
				edge->inter_pts[edge->ninter_pts-1]->y=
					streetnet->nodelist->intersects[cur_intersect_id]->gpos[1];
				edge->inter_pts[edge->ninter_pts-1]->cellid=get_cellID_givencoords(
					edge->inter_pts[edge->ninter_pts-1]->x, 
					edge->inter_pts[edge->ninter_pts-1]->y);
			}

			/*  update the cells that this edge pass!!  */
			int k;
			int pre_cell=-1;
			for(k=0;k<edge->ninter_pts;k++)
			{
				if(pre_cell == edge->inter_pts[k]->cellid) continue;
				pre_cell=edge->inter_pts[k]->cellid;
				if(pre_cell<0||pre_cell>=quadmesh->nfaces)
				{
					pre_cell=edge->inter_pts[k]->cellid=
						get_cellID_givencoords(edge->inter_pts[k]->x,
						edge->inter_pts[k]->y);
					if(pre_cell<0||pre_cell>=quadmesh->nfaces)
						continue;
				}
				QuadCell *face=quadmesh->quadcells[edge->inter_pts[k]->cellid];
				if(!is_repeated_elem(face->streetgraphedgelist, pre_cell,
					face->nstreetgraphedges))
				{
					add_to_edgelist_one_cell(pre_cell, edge->index);
				}
			}
			pre_cell=-1;
			for(k=0;k<temp_traj->nlinesegs;k++)
			{
				if(pre_cell==temp_traj->linesegs[k].Triangle_ID) continue;
				pre_cell = temp_traj->linesegs[k].Triangle_ID;
				QuadCell *face=quadmesh->quadcells[pre_cell];
				if(!is_repeated_elem(face->streetgraphedgelist, pre_cell,
					face->nstreetgraphedges))
				{
					add_to_edgelist_one_cell(pre_cell, edge->index);
				}
			}

			/*  update the original edge and   
				generate one new edges in the street network (graph)  
			*/
			edge=streetnet->edgelist->edges[which_edge];
			/*  delete the old edge  */
			edge->cancel=true;

			/*  create two new edges  */

			/* 1: */
			StreetGraphEdge *onenewedge=(StreetGraphEdge *)
				malloc(sizeof(StreetGraphEdge));
			onenewedge->cancel=
				onenewedge->visited=false;
			onenewedge->node_index1=cur_intersect_id;
			onenewedge->node_index2=edge->node_index1;
			onenewedge->roadtype = edge->roadtype;
			onenewedge->cancel=onenewedge->visited=false;

			streetnet->edgelist->append(onenewedge);

			double start[2], end[2];
			int start_cell, end_cell;
			start[0]=streetnet->nodelist->intersects[onenewedge->node_index1]->gpos[0];
			start[1]=streetnet->nodelist->intersects[onenewedge->node_index1]->gpos[1];
			start_cell=get_cellID_givencoords(start[0], start[1]);
			end[0]=streetnet->nodelist->intersects[onenewedge->node_index2]->gpos[0];
			end[1]=streetnet->nodelist->intersects[onenewedge->node_index2]->gpos[1];
			end_cell=get_cellID_givencoords(end[0], end[1]);
			Trajectory *tempTraj=new Trajectory(-1);
			get_linesegs_anytwopts(start,start_cell, end,end_cell, tempTraj, 0, 10);
			onenewedge->inter_pts=(Point **)malloc(sizeof(Point *)* (tempTraj->nlinesegs+1));

			int m;
			pre_cell=-1;
			for(m=0;m<tempTraj->nlinesegs;m++)
			{
				onenewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
				onenewedge->inter_pts[m]->x=tempTraj->linesegs[m].gstart[0];
				onenewedge->inter_pts[m]->y=tempTraj->linesegs[m].gstart[1];
				onenewedge->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge->inter_pts[m]->x,
					onenewedge->inter_pts[m]->y);
				
				if(pre_cell!=onenewedge->inter_pts[m]->cellid)
				{
					pre_cell=onenewedge->inter_pts[m]->cellid;
					add_to_edgelist_one_cell(onenewedge->inter_pts[m]->cellid, 
						onenewedge->index);
				}
			}

			onenewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
			onenewedge->inter_pts[m]->x=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[0];
			onenewedge->inter_pts[m]->y=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[1];
			onenewedge->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge->inter_pts[m]->x,
				onenewedge->inter_pts[m]->y);
			onenewedge->ninter_pts=tempTraj->nlinesegs+1;
			if(pre_cell!=onenewedge->inter_pts[m]->cellid)
			{
				pre_cell=onenewedge->inter_pts[m]->cellid;
				add_to_edgelist_one_cell(onenewedge->inter_pts[m]->cellid, 
					onenewedge->index);
			}
			delete tempTraj;

			/*  add to the adjacent edge lists of the two end intersections  */
			add_edge_to_intersectnode(onenewedge->node_index1, onenewedge->index);
			add_edge_to_intersectnode(onenewedge->node_index2, onenewedge->index);


			/* 2: */
			StreetGraphEdge *onenewedge2=(StreetGraphEdge *)
				malloc(sizeof(StreetGraphEdge));
			onenewedge2->cancel=
				onenewedge2->visited=false;
			onenewedge2->node_index1=cur_intersect_id;
			onenewedge2->node_index2=edge->node_index2;
			onenewedge2->roadtype = edge->roadtype;
			onenewedge2->cancel=onenewedge->visited=false;

			streetnet->edgelist->append(onenewedge2);

			start[0]=streetnet->nodelist->intersects[onenewedge2->node_index1]->gpos[0];
			start[1]=streetnet->nodelist->intersects[onenewedge2->node_index1]->gpos[1];
			start_cell=get_cellID_givencoords(start[0], start[1]);
			end[0]=streetnet->nodelist->intersects[onenewedge2->node_index2]->gpos[0];
			end[1]=streetnet->nodelist->intersects[onenewedge2->node_index2]->gpos[1];
			end_cell=get_cellID_givencoords(end[0], end[1]);
			tempTraj=new Trajectory(-1);
			get_linesegs_anytwopts(start,start_cell, end,end_cell, tempTraj, 0, 10);
			onenewedge2->inter_pts=(Point **)malloc(sizeof(Point *)* (tempTraj->nlinesegs+1));

			pre_cell=-1;
			for(m=0;m<tempTraj->nlinesegs;m++)
			{
				onenewedge2->inter_pts[m]=(Point*)malloc(sizeof(Point));
				onenewedge2->inter_pts[m]->x=tempTraj->linesegs[m].gstart[0];
				onenewedge2->inter_pts[m]->y=tempTraj->linesegs[m].gstart[1];
				onenewedge2->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge2->inter_pts[m]->x,
					onenewedge2->inter_pts[m]->y);
				
				if(pre_cell!=onenewedge2->inter_pts[m]->cellid)
				{
					pre_cell=onenewedge2->inter_pts[m]->cellid;
					add_to_edgelist_one_cell(onenewedge2->inter_pts[m]->cellid, 
						onenewedge2->index);
				}
			}

			onenewedge2->inter_pts[m]=(Point*)malloc(sizeof(Point));
			onenewedge2->inter_pts[m]->x=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[0];
			onenewedge2->inter_pts[m]->y=tempTraj->linesegs[tempTraj->nlinesegs-1].gend[1];
			onenewedge2->inter_pts[m]->cellid=get_cellID_givencoords(onenewedge2->inter_pts[m]->x,
				onenewedge2->inter_pts[m]->y);
			onenewedge2->ninter_pts=tempTraj->nlinesegs+1;
			if(pre_cell!=onenewedge2->inter_pts[m]->cellid)
			{
				pre_cell=onenewedge2->inter_pts[m]->cellid;
				add_to_edgelist_one_cell(onenewedge2->inter_pts[m]->cellid, 
					onenewedge2->index);
			}
			delete tempTraj;

			/*  add to the adjacent edge lists of the two end intersections  */
			add_edge_to_intersectnode(onenewedge2->node_index1, onenewedge2->index);
			add_edge_to_intersectnode(onenewedge2->node_index2, onenewedge2->index);
					
		}

		else
		{
			streetnet->nodelist->intersects[cur_intersect_id]->deleted=true;
			streetnet->edgelist->edges[except_edge]->cancel=true;
		}
	}
}

/*  we still need a routine to connect the dead ends of the original network
along the boundary 
*/

void connect_outer_deadends_Reg()
{
	/*  search all the boundary intersections  */
	int i;
	double A[2];
	icVector2 dir;
	double ang;
	Intersection *cur_intersect, *other_intersect;
	int except_edge;
	StreetGraphEdge *cur_edge;
	int cell_id;
	double dist_intersect;
	int which_intersect=-1;
	double the_intersect[2];
	int which_cell;

	for(i=0;i<nboundaryintersections;i++)
	{
		if(boundaryintersections[i]<0||
			boundaryintersections[i]>=streetnet->nodelist->nelems)
			continue;
		
		if(!is_inregion(streetnet->nodelist->intersects[boundaryintersections[i]]->gpos[0],
			streetnet->nodelist->intersects[boundaryintersections[i]]->gpos[1]))
			continue;

		cur_intersect=streetnet->nodelist->intersects[boundaryintersections[i]];

		if(!cur_intersect->endpt)
			continue;

		except_edge=cur_intersect->adj_edges[0];
		cur_edge=streetnet->edgelist->edges[except_edge];
		if(cur_intersect->nadjedges>1)
		{
			int test=0;
		}

		cell_id=cur_intersect->cellid;
		A[0]=cur_intersect->gpos[0];
		A[1]=cur_intersect->gpos[1];

		if(cur_edge->node_index1==cur_intersect->index)
		{
			other_intersect=streetnet->nodelist->intersects[cur_edge->node_index2];
			dir.entry[0]=A[0]-cur_edge->inter_pts[1]->x;
			dir.entry[1]=A[1]-cur_edge->inter_pts[1]->y;
		}
		else
		{
			other_intersect=streetnet->nodelist->intersects[cur_edge->node_index1];
			dir.entry[0]=A[0]-cur_edge->inter_pts[cur_edge->ninter_pts-1]->x;
			dir.entry[1]=A[1]-cur_edge->inter_pts[cur_edge->ninter_pts-1]->y;
		}

		if(!cur_intersect->endpt)
			continue;

		/**/
		
		ang=atan2(dir.entry[1], dir.entry[0]);
		if(ang<0) ang+=(2.*M_PI);
		
		if(cal_smallest_dist_to_intersects(cell_id, A, ang, dir, dist_intersect, which_intersect,
			cur_intersect->index))
		{
				StreetGraphEdge *edge=streetnet->edgelist->edges[except_edge];

				the_intersect[0]=streetnet->nodelist->intersects[which_intersect]->gpos[0];
				the_intersect[1]=streetnet->nodelist->intersects[which_intersect]->gpos[1];
				which_cell=streetnet->nodelist->intersects[which_intersect]->cellid;

				Trajectory *temp_traj = new Trajectory(-1);
				get_linesegs_anytwopts(
					cur_intersect->gpos, 
					cur_intersect->cellid,
					the_intersect, which_cell,
					temp_traj, 0, 10);

				/* we add a new edge  */
				StreetGraphEdge *anewedge=(StreetGraphEdge *)malloc(sizeof(StreetGraphEdge));
				anewedge->node_index1=cur_intersect->index;
				anewedge->node_index2=which_intersect;
				anewedge->roadtype=edge->roadtype;
				anewedge->cancel=anewedge->visited=false;
				anewedge->inter_pts=(Point**)malloc(sizeof(Point*)*(temp_traj->nlinesegs+1));
				streetnet->edgelist->append(anewedge);


				/*  possible bug:  remember to add this edge to the "which_intersect"  
						1/19/2008
				*/
				int m;
				int pre_cell=-1;
				for(m=0;m<temp_traj->nlinesegs;m++)
				{
					anewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
					anewedge->inter_pts[m]->x=temp_traj->linesegs[m].gstart[0];
					anewedge->inter_pts[m]->y=temp_traj->linesegs[m].gstart[1];
					anewedge->inter_pts[m]->cellid=get_cellID_givencoords(anewedge->inter_pts[m]->x,
						anewedge->inter_pts[m]->y);
					
					if(pre_cell!=anewedge->inter_pts[m]->cellid)
					{
						pre_cell=anewedge->inter_pts[m]->cellid;
						add_to_edgelist_one_cell(anewedge->inter_pts[m]->cellid, 
							anewedge->index);
					}
				}

				anewedge->inter_pts[m]=(Point*)malloc(sizeof(Point));
				anewedge->inter_pts[m]->x=temp_traj->linesegs[temp_traj->nlinesegs-1].gend[0];
				anewedge->inter_pts[m]->y=temp_traj->linesegs[temp_traj->nlinesegs-1].gend[1];
				anewedge->inter_pts[m]->cellid=get_cellID_givencoords(anewedge->inter_pts[m]->x,
					anewedge->inter_pts[m]->y);
				anewedge->ninter_pts=temp_traj->nlinesegs+1;
				if(pre_cell!=anewedge->inter_pts[m]->cellid)
				{
					pre_cell=anewedge->inter_pts[m]->cellid;
					add_to_edgelist_one_cell(anewedge->inter_pts[m]->cellid, 
						anewedge->index);
				}

				streetnet->nodelist->intersects[which_intersect]->endpt=false;
				cur_intersect->endpt=false;

				/*   add to intersection list  */
				add_edge_to_intersectnode(which_intersect, anewedge->index);
				add_edge_to_intersectnode(cur_intersect->index, anewedge->index);

				delete temp_traj;
		}

		else
		{
			/*   we just simply back track the edges having dangling dead ends   */
			cur_intersect->deleted=true;
			cur_edge->cancel=true;
		}

	}
}


extern bool is_at_validDir(int cell, double pt[2], icVector2 trajDir);

void obtain_neighboringCells(int cell, int *neighbors, int &Nneighbors)
{
	int i, j;

	/**/
	for(i=-1;i<=1;i++)
	{
		/*  different row  */
		if(cell+i*(quadmesh->XDIM-1)<0
			||cell+i*(quadmesh->XDIM-1)>=quadmesh->nfaces)
			continue;

		for(j=-1;j<=1;j++)
		{
			if((cell+j)/(quadmesh->XDIM-1)!=cell/(quadmesh->XDIM-1))
				continue;

			neighbors[Nneighbors]=cell+i*(quadmesh->XDIM-1)+j;
			Nneighbors++;
		}
	}

}


bool search_nearby_intersection_around_oneCell(int cell, double threshold, double pt[2], int except_intersect,
											   icVector2 goDir, double cosang, int &which_intersect)
{
	/*  search the neighboring */
	int neighbors[9];
	int Nneighbors=0;

	int i, j;
	double dist=1.e50;
	obtain_neighboringCells(cell, neighbors, Nneighbors);
	QuadCell *face;
	icVector2 dir_to_intersect;
	double tempDist=0;
	Intersection *intersect=NULL;
	bool found=false;

	for(i=0;i<Nneighbors;i++)
	{
		if(!is_at_validDir(neighbors[i], pt, goDir))
			continue;

		face=quadmesh->quadcells[neighbors[i]];

		for(j=0;j<face->nintersects;j++)
		{
			if(face->intersectlist[j]==except_intersect)
				continue;

			intersect=streetnet->nodelist->intersects[face->intersectlist[j]];
			dir_to_intersect.entry[0]=intersect->gpos[0]-pt[0];
			dir_to_intersect.entry[1]=intersect->gpos[1]-pt[1];

			tempDist=length(dir_to_intersect);

			if(tempDist==0.0)  /*  could they be the same intersection?  */
			{
				int test=0;
			}

			normalize(dir_to_intersect);

			if(dot(goDir, dir_to_intersect)<cosang || tempDist>threshold) /* not in the valid direction */
				continue;

			if(tempDist<dist)
			{
				dist=tempDist;
				which_intersect=face->intersectlist[j];
				found=true;
			}
		}
	}

	return found;
}


/*  compute the intersection with the nearby edge  */
bool search_nearby_intersection_onEdge_oneCell(int cell,/*double threshold,*/double p1[2], double p2[2],
											   icVector2 goDir, double cosang, int except_edge, 
											   int &which_edge, int &which_samp, double intersect[2])
{
	QuadCell *face=quadmesh->quadcells[cell];
	int i;
	icVector2 dir_to_intersect;
	double tempDist;
	double dist=1.e50;
	int temp_samp;
	for(i=0;i<face->nstreetgraphedges;i++)
	{
		if(streetnet->edgelist->edges[face->streetgraphedgelist[i]]->cancel)
			continue;

		if(except_edge==face->streetgraphedgelist[i])
			continue;

		if(cal_intersect_at_graph_edge(face->streetgraphedgelist[i],
			intersect, temp_samp, p1, p2))
		{
			dir_to_intersect.entry[0]=intersect[0]-p1[0];
			dir_to_intersect.entry[1]=intersect[1]-p1[1];
			tempDist=length(dir_to_intersect);

			normalize(dir_to_intersect);


			if(dot(dir_to_intersect, goDir)<cosang)
				continue;

			if(tempDist<dist)
			{
				dist=tempDist;
				which_samp=temp_samp;
				which_edge=face->streetgraphedgelist[i];
				return true;
			}
		}
	}

	return false;
}

/*
   start from specify end points trace along the input direction 
*/

int trace_comp_nearby_intersect(double startp[2], int start_cell, icVector2 goDir, 
								int except_intersect, int except_edge, 
								int &which_intersect, double intersect_p[2],
								int &which_edge, int &which_samp, int MaxSearchCells)
{
	int i,j;

	int cur_cell=start_cell;
	double pre_p[2],cur_p[2];

	normalize(goDir);
	tenline_dir_global = goDir;

	pre_p[0]=startp[0];
	pre_p[1]=startp[1];

	cur_p[0]=pre_p[0]+goDir.entry[0];
	cur_p[1]=pre_p[1]+goDir.entry[1];

	icVector2 t_major[4];
	QuadCell *face;

	for(i=0;i<MaxSearchCells;i++)
	{
		/*  if we find a near by valid intersection, return 1  */
		if(search_nearby_intersection_around_oneCell(cur_cell, quadmesh->xinterval, startp/*pre_p*/, 
			except_intersect, goDir, 0, which_intersect))
		{
			intersect_p[0]=pre_p[0];
			intersect_p[1]=pre_p[1];
			return 1;  /*  we find a nearby intersection  */
		}

		/*  if we find an intersection with a nearby edge, return 2  */
		if(search_nearby_intersection_onEdge_oneCell(cur_cell, startp, cur_p, goDir, 0.,
			except_edge, which_edge, which_samp, intersect_p))
		{
			return 2;
		}

		if(cur_cell<0 || cur_cell>=quadmesh->nfaces)
			return 0;

		if(!is_in_cell(cur_cell, cur_p[0], cur_p[1]))
		{
			face = quadmesh->quadcells[cur_cell];

			for(j=0; j<4; j++)
				t_major[j]=quadmesh->quad_verts[face->verts[j]]->major;

			/*replace them with the line segment direction*/
			for(j=0; j<4; j++)
				quadmesh->quad_verts[face->verts[j]]->major=tenline_dir_global;

			int passvertornot = 0;
			double tp[2]={pre_p[0], pre_p[1]};
			get_next_cell_2(cur_cell, pre_p, cur_p, passvertornot, 0);
			
			if(passvertornot==0)
			{
				pre_p[0]=cur_p[0];
				pre_p[1]=cur_p[1];
			}

			cur_p[0]=startp[0]+goDir.entry[0];
			cur_p[1]=startp[1]+goDir.entry[1];

			/*store back the original vectors*/
			for(j=0; j<4; j++)
				quadmesh->quad_verts[face->verts[j]]->major=t_major[j];
		}
	}

	return 0;  //doesn't find suitable intersection
}



/*
   for brush based region selection 
*/
