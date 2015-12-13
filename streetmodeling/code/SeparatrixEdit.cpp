////SeparatrixEdit.cpp

#include "stdafx.h"

#include "SeparatrixEdit.h"

#include "topologyedit.h"

#include "VFDataStructure.h"

#include "HermiteCurve.h"

#include "VFSynthesis.h"

#include "LimitCycleCreator.h"

#include "RegionSmoothing.h"

////variables for singularities pair cancellation and movement
extern TriangularRegion repellerRegion;       ////region containing a repeller
extern TriangularRegion attractorRegion;      ////region containing an attractor
extern RegionBoundary repellerBoundary;
extern RegionBoundary attractorBoundary;

extern TriangularRegion intersectRegion;     ////The intersect region
extern RegionBoundary intersectBoundary;     ////The intersect boundary 
extern InnerVertices intersectInnerverts;    ////The inner vertices inside the intersect region

extern Separatrices *separatrices;             //array for group of separatrices
extern Singularities *singularities;           //being captured singularites' list

extern LineSeg **trajectories;                 //trajectories' list
extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory

extern int resolution;                // how many points in our output array of the controlling curve
extern ctr_point *out_pts;

extern ctr_point *control_pts;        // allocate our control point array
extern int num_shapecontrol_pts;
extern int num_curvepts_output;

extern void Hermitecurve_open(int n, ctr_point *interpts, icVector2 *T, ctr_point *output, int step, int &num_output);


//////////////////////////
extern Edge **Cycle_edge;
extern int num_cycleedges;
extern int MaxEdgeInCycle;

extern int *DesignCurveCellCycle;
extern int num_triangles_designcurve;
extern int MaxNumTrianglesDesignCurve;


////Variable to store the points and the corresponding triangles that contain the points
extern LineSeg *designcurve;
extern int num_lineseg_designcurve;

//////////////////////////
extern Polygon3D Object;

extern double dmax;

////Extern routines
extern void AddToInnerVerts(Vertex *, int);
extern int GetOppositeTriangle(Edge *, int);
extern bool AttractorExitEdgePending(Edge *);
extern bool RepellerExitEdgePending(Edge *);

extern bool TriangleSearch(int *acycle, int num_cycletriangles, int oneTriangle, int &position);

extern void UnionRegion(TriangularRegion &source1, TriangularRegion &source2, TriangularRegion &dest);

///////////
TriangularRegion newTriangleStrip;

////we may need the third variable to tell the routine what kind of separatrix it is, an incoming or outgoing one
void SeparatrixModify(int saddleID, int trajID)
{
	////Grow region
	GetInitSeparatrixStrip(saddleID, trajID);
	
	int inorout;
	////Get the incoming or outgoing feature of the separatrix
	if(trajID == separatrices[singularities[saddleID].separtices].sep1
		|| trajID == separatrices[singularities[saddleID].separtices].sep3)
		inorout = 0;  ////it is outgoing separatrix
	else
		inorout = 1;  ////it is incoming separatrix
	
	////Grow the region for the separatrix modification
	GetTheRegionForSeparatrixModification(inorout);


	////Get the new triangle strip that cover the new separatrix
	GetTheNewTriangleStrip();

	////Set the vectors on the boundary of the new triangle strip
	SetVectorsOnNewBoundary(inorout);

	////Get the inner vertices of the region 
	GetInnerVertices();

	////Build the final region and smooth it
	Cancel_RegionSmooth();

}


////This routine should be called after finding the intersection region 12/28/05
void GetInnerVertices()
{
	int i, j;
	Face *cur_f;
	Vertex *vert;
	
	//Initialization part, reset all the relative flags
	for(i = 0; i < Object.nverts; i++)
	{
		vert = Object.vlist[i];
		if(vert->OnBoundary == 1)
		{
			vert->InRegion = 0;
			vert->RegionListID = -1;
			continue;
		}

		vert->OnBoundary = 0;
		vert->InRegion = 0;
		vert->RegionListID = -1;
	}

	UpdateBoundary(2);

	////Set the vertices on the boundary as 'OnBoundary'
	for(i = 0; i < intersectBoundary.num; i++)
	{
		Object.vlist[intersectBoundary.edgelist[i]->verts[0]]->OnBoundary = 1;
		Object.vlist[intersectBoundary.edgelist[i]->verts[0]]->InRegion = 0;
		Object.vlist[intersectBoundary.edgelist[i]->verts[1]]->OnBoundary = 1;
		Object.vlist[intersectBoundary.edgelist[i]->verts[1]]->InRegion = 0;
	}

	////Get the intersection of inner vertices
	intersectInnerverts.num = 0;
	for(i = 0; i < intersectRegion.num; i++)
	{
		cur_f = Object.flist[intersectRegion.trianglelist[i]];

		for(j = 0; j < cur_f->nverts; j++)
		{
			vert = Object.vlist[cur_f->verts[j]];
			if(vert->OnBoundary == 0 && vert->InRegion == 0) ////The vertex is inside the region
			{
				intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
				vert->InRegion = 1;
				vert->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
				intersectInnerverts.num ++;
			}
		}
	}

}


////Initialize the triangle strip according to the original separatrix
void GetInitSeparatrixStrip(int saddleID, int trajindex)
{
	int i;
	int position;
	int cur_t;

	int numlinesegs = num_linesegs_curtraj[trajindex];

	//Initialize part
	repellerRegion.num = 0;
	attractorRegion.num = 0;

	//Search for the initial triangle strip for both region
	for(i = 0; i < numlinesegs; i++)
	{
		cur_t = trajectories[trajindex][i].Triangle_ID;

		if(cur_t < 0) //probably reaches boundary
			break;

		////Do not cover the saddle and the other singularity at the other end of the separatrix
		if(Object.flist[cur_t]->contain_singularity == 1)
			continue;

		if(!TriangleSearch(repellerRegion.trianglelist, repellerRegion.num, cur_t, position))
		{
			repellerRegion.trianglelist[repellerRegion.num] = cur_t;
			repellerRegion.num ++;

			attractorRegion.trianglelist[attractorRegion.num] = cur_t;
			attractorRegion.num ++;
		}
	}

	////Get boundaries for the two regions
	UpdateBoundary(1);
	GetRegionNormals(1); ////Get the outward normal of the region for each boundary edge

	UpdateBoundary(0);
	GetRegionNormals(0); ////Get the outward normal of the region for each boundary edge
}

////Get the initial smoothing region that is a little bit larger
void GetTheRegionForSeparatrixModification(int inorout)
{
	//1. Grow the regions for repeller and attractor respectively
	//Separatrix_Growing(0);
	//Separatrix_Growing(1);

	////2. Union the above two regions to get the final region
	//UnionRegion(repellerRegion, attractorRegion, intersectRegion);

	Separatrix_Growing(1-inorout);

	if(inorout == 1)
	{
		CopyRegion(repellerRegion.trianglelist, intersectRegion.trianglelist, repellerRegion.num);
		intersectRegion.num = repellerRegion.num;
	}
	else
	{
		CopyRegion(attractorRegion.trianglelist, intersectRegion.trianglelist, attractorRegion.num);
		intersectRegion.num = attractorRegion.num;
	}
}


void GetTheNewTriangleStrip()
{
	GetDesignCellCycle_new();
}


/////////////////////////////////////////////////////////////////////////
///There are so many repeated codes here as GetBoundary() routine

///Which edge list should store them
void GetBoundaryofTriangleStrip()
{
	////First we need to perform manifold testing and correct the original region
	int i;
	int EndVertID;
	Vertex *cur_vert;
	Edge *cur_edge;

	Edge **temp_edgelist = (Edge **) malloc(sizeof(Edge *) * MaxEdgeInCycle);

	int temp_num_edges, other_num_edges;
	temp_num_edges = 0;

	////rebuild boundary again
	for(i = 0; i < num_cycleedges; i++)
	{
		Cycle_edge[i]->OnBoundary = 0;
	}

	BuildBoundaryEdgeList(DesignCurveCellCycle, num_triangles_designcurve);

	////Initial the flag of all the edges on current boundaries
	for(i = 0; i < num_cycleedges; i++)
	{
		Cycle_edge[i]->BoundaryVisited = 0;
	}

	temp_num_edges = 0;
	temp_edgelist[0] = Cycle_edge[0];
	temp_num_edges++;
	EndVertID = Cycle_edge[0]->verts[0];
	temp_edgelist[0]->BoundaryVisited = 1;

	cur_vert = Object.vlist[Cycle_edge[0]->verts[1]];

	while(cur_vert->VertID != EndVertID)   ////Not form a closed edges list
	{
		for(i = 0; i < cur_vert->Num_edge; i++)
		{
			cur_edge = cur_vert->edges[i];
			if(cur_edge->OnBoundary == 1 && cur_edge->BoundaryVisited == 0)
			{
				temp_edgelist[temp_num_edges] = cur_edge;
				temp_num_edges ++;
				cur_edge->BoundaryVisited = 1;
	 		}
	    }
		
		////Get next testing vertex
		if(temp_edgelist[temp_num_edges-1]->verts[0] != cur_vert->VertID)
			cur_vert = Object.vlist[temp_edgelist[temp_num_edges-1]->verts[0]];
		else
			cur_vert = Object.vlist[temp_edgelist[temp_num_edges-1]->verts[1]];
    }

	free(temp_edgelist);

}

////Set the vectors on the boundary of the triangle strip that contains the new separatrix
void SetVectorsOnNewBoundary(int inorout)
{

	////Do we need to extend the triangle strip to generate a smoother result ??? 12/28/05

	//1. Get the boundary
    GetBoundaryofTriangleStrip();

	//2. Get the inward normals of the edges on the boundary
    GetNormalsForBoundaryEdges(1);

	//3. Assign the vectors according to these normals
	SetVectorOnVertices(inorout);

	//4. Normalize the vectors on the boundary
	int i;
	Vertex *cur_v;
	double r;

	for(i = 0; i < Object.nverts; i++)
	{
		cur_v = Object.vlist[i];
		if(cur_v->OnBoundary == 1)
		{
			r = length(cur_v->vec);
			r *= r;
						
			if (r < DistanceThreshold) 
			{
				r = DistanceThreshold;
				cur_v->vec *= dmax/r; 
			}

			r = length(cur_v->vec);
			r *= r;

			if (r > dmax*dmax) { 
				r  = sqrt(r); 
				cur_v->vec *= dmax/r; 
			}
		}
	}
}



////Find which triangle in the design triangle strip containing the edge
int WhichCellTriangleContainTheEdge(Edge* edge)
{
	int i, j;
	Face *cur_f;
	Edge *cur_e;

	for(i = 0; i < num_triangles_designcurve; i++)
	{
		cur_f = Object.flist[DesignCurveCellCycle[i]];
		for(j = 0; j < 3; j++)
		{
			cur_e = cur_f->edges[j];
			if(edge == cur_e)
				return cur_f->index;
		}
	}
}

void SetVectorOnVertices(int inorout)
{
	int i, j, k;
	int cur_triangle;
	Vertex *cur_v;
	Edge *cur_e, *next_e;
	icVector2 vert_normal, curvedirect;
	double theta = (45./90.) * (M_PI/2.);
	
	for(i = 0; i < num_cycleedges - 1; i++)
	{
		cur_e = Cycle_edge[i];
		next_e = Cycle_edge[i+1];

		//1. Find the common vertex of the two edges
		for(j = 0; j < 2; j++)
		{
			for(k = 0; k < 2; k++)
			{
				if(cur_e->verts[j] == next_e->verts[k])
				{
					cur_v = Object.vlist[cur_e->verts[j]];
					break;
				}
			}
		}

		//2. we need to get the normal for the shared vertex
		vert_normal = cur_e->normal + next_e->normal;
		normalize(vert_normal);

		//3. Get the current direction of the design curve
		cur_triangle = WhichCellTriangleContainTheEdge(cur_e);
		curvedirect = GetCurveDirection(cur_triangle);

		//4. According to the incoming or outgoing feature of the separatrix to set the vector on the vertex
		cur_v->vec = GetAVector(curvedirect, vert_normal, theta, inorout);
		cur_v->OnBoundary = 1;
	}
}


////Get the design curve direction inside the specific triangle  12/28/05
icVector2 GetCurveDirection(int triangleID)
{
	int i;
	icVector2 result;

	for(i = 0; i < num_lineseg_designcurve; i++)
	{
		if(designcurve[i].Triangle_ID == triangleID)
			break;
	}

	if(i == num_lineseg_designcurve)  ////the triangle does not contain any controlling points
	{
		return result;   ////return zero vector;
	}

	////Calculate the vector
	result.entry[0] = designcurve[i].gend[0] - designcurve[i].gstart[0];
	result.entry[1] = designcurve[i].gend[1] - designcurve[i].gstart[1];

	return result;
}


////according to the curve direction and the feature of the separatrix, we design the new vector
////by rotating the input normal
icVector2 GetAVector(icVector2 curveorient, icVector2 normal, double ang, int inorout)
{
	icVector2 ccw_vec, cw_vec;

	////calculate the vector by rotating along clockwise orientation
	cw_vec.entry[0] = cos(ang)*normal.entry[0]+sin(ang)*normal.entry[1];
	cw_vec.entry[1] = -sin(ang)*normal.entry[0]+cos(ang)*normal.entry[1];

	////calculate the vector by rotating along counter clockwise orientation
	ccw_vec.entry[0] = cos(ang)*normal.entry[0]-sin(ang)*normal.entry[1];
	ccw_vec.entry[1] = sin(ang)*normal.entry[0]+cos(ang)*normal.entry[1];

	////normalize the vector
	normalize(cw_vec);
	normalize(ccw_vec);
	normalize(curveorient);

	////
	if(inorout == 0)  //it is an incoming separatrix
	{
		if(dot(cw_vec, curveorient) < dot(ccw_vec, curveorient))
		{
			return cw_vec;  //cw_vec more opposite to the curve direction
		}

		else
			return ccw_vec;
	}

	else
	{
		if(dot(cw_vec, curveorient) > dot(ccw_vec, curveorient))
			return cw_vec;  //cw_vec more close to the curve direction
		else
			return ccw_vec;
	}
}


/************************************************************
The main routine for region growing for separatrix editing
************************************************************/
void Separatrix_Growing(int type)
{
	int i;
	bool exitornot = false;
	int num_edges;
	Edge **edgelist;
	Edge *cur_edge;
    int oppositeTriangle;
	int num_newaddedtriangle;

	if(type == 0)
	{
		edgelist = repellerBoundary.edgelist;
		num_edges = repellerBoundary.num;
	}
	else{
		edgelist = attractorBoundary.edgelist;
		num_edges = attractorBoundary.num;
	}

	while(1)
	{
		num_newaddedtriangle = 0;
		for(i = 0; i < num_edges; i++)
		{
			cur_edge = edgelist[i];

			if(type == 0)  ////repeller
			{
				exitornot = RepellerExitEdgePending(cur_edge);
			}
			else{          ////attractor
				exitornot = AttractorExitEdgePending(cur_edge);
			}

			oppositeTriangle = GetOppositeTriangle(cur_edge, type);

			if(oppositeTriangle < 0) ////not suitable for 3D surfaces
				continue;
			
			////if the opposite triangle containing singularity or separatrices
			if(Object.flist[oppositeTriangle]->contain_singularity == 1
				|| Object.flist[oppositeTriangle]->contain_separatrix == 1) 
			{
				continue;
			}

			if(exitornot && oppositeTriangle >= 0)  ////it is an exiting edge 
			{
				////if the opposite triangle is not inside the region!!!
				////add the adjacent triangle into the region

				AddToRegionTriangles(oppositeTriangle, type);
				num_newaddedtriangle++;

				//////add the two ending points of current edge to the inner vertices list
				AddToInnerVerts(Object.vlist[cur_edge->verts[0]], type);
				AddToInnerVerts(Object.vlist[cur_edge->verts[1]], type);

			}

		}

		if(num_newaddedtriangle == 0) ////No more exiting edges can be found
			return;

		////Update the boundary list under new region
		UpdateBoundary(type);
		GetRegionNormals(type);

		if(type == 0)
			num_edges = repellerBoundary.num;
		else
			num_edges = attractorBoundary.num;
	}
}


//////////////////////////////////////////////////////////////
void DiscretizeSeparatrix(int trajID)
{
	int i;
	////Using default 10 controlling points here now
	int num_curvepoints = num_linesegs_curtraj[trajID];

	////Note that: do not close the controlling curve
	if(num_curvepoints - 2 <= 10)
	{
		for(i = 0; i < 10; i++)
		{
			control_pts[i].x = trajectories[trajID][i+1].gstart[0];
			control_pts[i].y = trajectories[trajID][i+1].gstart[1];
		}
	}

	else
	{
		int interval = (int)(num_curvepoints/10.);
		for(i = 0; i < 10; i++)
		{
			control_pts[i].x = trajectories[trajID][i*interval+1].gstart[0];
			control_pts[i].y = trajectories[trajID][i*interval+1].gstart[1];
		}
	}

	num_shapecontrol_pts = 10;

}



void AllocSeparatrixEdit()
{
}


void InitSeparatrixEdit()
{
	num_cycleedges = 0;
}