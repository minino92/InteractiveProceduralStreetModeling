////topologyedit.cpp


#include "stdafx.h"
#include "topologyedit.h"
#include "VFDataStructure.h"
#include "LocalTracing.h"
#include "RegionSmoothing.h"
#include "VFSynthesis.h"
#include "LimitCycleEdit.h"
#include "LimitCycleCreator.h"
#include "Numerical.h"

#include "gl/glut.h"

#include "ConleyRelationGraph.h"
#include "PairCancellationModule.h"

////variables for singularities pair cancellation and movement
TriangularRegion repellerRegion;       ////region containing a repeller
TriangularRegion attractorRegion;      ////region containing an attractor
RegionBoundary repellerBoundary;
RegionBoundary attractorBoundary;
InnerVertices  repellerInnerverts;
InnerVertices  attractorInnerverts;


int MaxNumTriangle ;
int MaxNumBoundaryEdges;
int MaxNumInnerVerts;

////How many variables we need to store all the possible region growing? 11/06/05
TriangularRegion Source_re1, Source_re2, Source_re;
TriangularRegion Sink_re1, Sink_re2, Sink_re;

TriangularRegion intersectRegion;     ////The intersect region
RegionBoundary intersectBoundary;     ////The intersect boundary 
InnerVertices intersectInnerverts;    ////The inner vertices inside the intersect region

icVector2 *VerticalField;             ////Tempory field for singularity movement which is perpendicular to original field
icMatrix3x3 pretransform_matrix;
icMatrix3x3 compensate_matrix;        ////matrix for compensated transformation of pair cancellation
/*-------------------------------------------------*/


extern int SelectTriangleID;
extern Polygon3D Object;
extern int TotalEdgesNum;
extern double dmax;
extern Singularities *singularities;           //being captured singularites' list
extern int cur_singularity_index;

extern LineSeg **trajectories;                 //trajectories' list
extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory
extern Separatrices *separatrices;             //array for group of separatrices
extern int cur_separatrices_index;

extern int InitialStateOn;
/*-------------------------------------------------*/

////variables for numerical calculation routines
extern Vec_INT *ija_p;
extern Vec_DP *sa_p;

/*-------------------------------------------------*/
extern GraphNode *graphnodes ;
extern GraphEdge *graphedges ;
extern int cur_node_index;
extern int cur_graphedge_index;
extern int *MediaNodes;
extern int Num_MediaNodes;
extern int MaxMediaNodes;

extern int *repeller_nodes;     //the indices of the selected repellers in the conley graph
extern int NumRepellers;        //the number of the being selected repellers
extern int *attractor_nodes;    //the indices of the selected attractors in the conley graph
extern int NumAttractors;       //the number of the being selected attractors

/*-------------------------------------------------*/
extern Edge **Cycle_edge;
extern int num_cycleedges;
extern int *DesignCurveCellCycle;
extern int num_triangles_designcurve;
extern int MaxNumTrianglesDesignCurve;


/*-------------------------------------------------*/
extern void GetLocalVector();
extern void CaptureSing();
extern void ReflectTheWholeField(double);
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

/*****************************************************************************/
////routins for singularities pair cancellation and movement
void UnionRegion(TriangularRegion &, TriangularRegion &, TriangularRegion &);
void IntersectRegion(TriangularRegion &source1, TriangularRegion &source2, TriangularRegion &dest);
void SetBoundaryFlag_Ver(Edge **boundaryedgelist, int num_edges);


bool IsRepeatedEdge(Edge **edgelist, int num_edges, Edge *anedge)
{
	int i;

	for(i = 0; i < num_edges; i++)
	{
		if(edgelist[i] == anedge)
			return true;
	}
	return false;
}

////Add the input edge into the boundary list
void AddToBoundaryEdges(Edge *cur_edge, int type)
{
	////if the space of the boundary is not enough, extend it
	int cur_boundarynum;
	if(type == 0)
		cur_boundarynum = repellerBoundary.num;
	else if(type == 1)
		cur_boundarynum = attractorBoundary.num;
	else
		cur_boundarynum = intersectBoundary.num;

	if(cur_boundarynum >= MaxNumBoundaryEdges - 1)
	{
		MaxNumBoundaryEdges += 50;
		repellerBoundary.edgelist = (Edge**)realloc(repellerBoundary.edgelist, sizeof(Edge*)*MaxNumBoundaryEdges);
		attractorBoundary.edgelist = (Edge**)realloc(attractorBoundary.edgelist, sizeof(Edge*)*MaxNumBoundaryEdges);
		intersectBoundary.edgelist = (Edge**)realloc(intersectBoundary.edgelist, sizeof(Edge*)*MaxNumBoundaryEdges);
	}

	if(type == 0)
	{
		if(!IsRepeatedEdge(repellerBoundary.edgelist, repellerBoundary.num, cur_edge))
		{
			repellerBoundary.edgelist[repellerBoundary.num] = cur_edge;
			repellerBoundary.num ++;
		}
	}
	else if(type == 1){
		if(!IsRepeatedEdge(attractorBoundary.edgelist, attractorBoundary.num, cur_edge))
		{
			attractorBoundary.edgelist[attractorBoundary.num] = cur_edge;
			attractorBoundary.num ++;
		}
	}
	else{  ////add to the intersected boundary
		if(!IsRepeatedEdge(intersectBoundary.edgelist, intersectBoundary.num, cur_edge))
		{
			intersectBoundary.edgelist[intersectBoundary.num] = cur_edge;
			intersectBoundary.num ++;
		}
	}
}


////Test whether the input triangle is inside current list or not
////Here we use the flag to quickly judge whether the triangle is in the list or not
bool InsideTrianglesListOrnot(int triangleID, int type)
{
	if(type == 0)
	{
		if(Object.flist[triangleID]->repell_inregion == 1)
			return true;
	}

	else{
		if(Object.flist[triangleID]->attract_inregion == 1)
			return true;
	}

	return false;
}

////Add the triangle to the region list
void AddToRegionTriangles(int triangleID, int type)
{
	if(triangleID < 0)  ////if the triangle does not exist
		return;

	int cur_trianglenum;
	if(type == 0)
		cur_trianglenum = repellerRegion.num;
	else
		cur_trianglenum = attractorRegion.num;

	////if the space for the triangles list is not enough, extend it
	if(cur_trianglenum >= MaxNumTriangle - 1)
	{
		MaxNumTriangle += 100;
		repellerRegion.trianglelist = (int*)realloc(repellerRegion.trianglelist, sizeof(int) * MaxNumTriangle);
		attractorRegion.trianglelist = (int*)realloc(attractorRegion.trianglelist, sizeof(int) * MaxNumTriangle);
		intersectRegion.trianglelist = (int*)realloc(intersectRegion.trianglelist, sizeof(int) * MaxNumTriangle);
	}

	////Make sure that the triangle not inside current triangles' list, very important!!!!

	if(InsideTrianglesListOrnot(triangleID, type))
		return;

	////set the flag of the triangle to mark it as a triangle inside the region
	if(type == 0)
	{
		Object.flist[triangleID]->repell_inregion = 1;
		repellerRegion.trianglelist[repellerRegion.num] = triangleID;
		repellerRegion.num++;
	}
	else{
		Object.flist[triangleID]->attract_inregion = 1;
		attractorRegion.trianglelist[attractorRegion.num] = triangleID;
		attractorRegion.num++;
	}
}

////Add the vertex to the inner vertices' list
void AddToInnerVerts(Vertex *vert, int type)
{
	int cur_innervertsnum;
	if(type == 0)
		cur_innervertsnum = repellerInnerverts.num;
	else
		cur_innervertsnum = attractorInnerverts.num;

	if(cur_innervertsnum >= MaxNumInnerVerts - 1)
	{
		MaxNumInnerVerts += 100;
		repellerInnerverts.vertslist = (Vertex**)realloc(repellerInnerverts.vertslist, sizeof(Vertex *) * MaxNumInnerVerts);
		attractorInnerverts.vertslist = (Vertex**)realloc(attractorInnerverts.vertslist, sizeof(Vertex *) * MaxNumInnerVerts);
		intersectInnerverts.vertslist = (Vertex**)realloc(intersectInnerverts.vertslist, sizeof(Vertex *) * MaxNumInnerVerts);
	}
	
	////add the input vertex to the inner vertices list
	if(type == 0)
	{
		if(vert->repell_flag == 2)
			return;

		repellerInnerverts.vertslist[repellerInnerverts.num] = vert;
		repellerInnerverts.num ++;

		vert->repell_flag = 2;
	}
	else{
		if(vert->attract_flag == 2)
			return;

		attractorInnerverts.vertslist[attractorInnerverts.num] = vert;
		attractorInnerverts.num ++;

		vert->attract_flag = 2;
	}

	////Update the flags of the being added vertices
	//if(type == 0)
	//	vert->repell_flag = 2;   ////set it as an inner vertex
	//else
	//    vert->attract_flag = 2;

}

////Get the opposite triangle for current triangle that sharing the specific edge
int GetOppositeTriangle(Edge *cur_edge, int type)
{
	////we just return the triangle that is not inside the region
	if(cur_edge->tris[0] < 0) //meet the boundary
		return  -1;

	if(type == 0)
	{
		if(Object.flist[cur_edge->tris[0]]->repell_inregion == 0)
			return cur_edge->tris[0];
		else
			return cur_edge->tris[1];
	}
	else
	{
		if(Object.flist[cur_edge->tris[0]]->attract_inregion == 0)
			return cur_edge->tris[0];
		else
			return cur_edge->tris[1];
	}
}


////There are only 6 cases for each triangle
int EdgeVerticesCase(Edge *cur_edge, Face *face)
{
	if(cur_edge->verts[0] == face->verts[0] && cur_edge->verts[1] == face->verts[1])
		return 0;           ////case 0->1, the same as triangle's orientation
	else if(cur_edge->verts[0] == face->verts[1] && cur_edge->verts[1] == face->verts[0])
		return 1;           ////case 1->0, the inverse to triangle's orientation
	else if(cur_edge->verts[0] == face->verts[1] && cur_edge->verts[1] == face->verts[2])
		return 2;           ////case 1->2, the same as triangle's orientation
	else if(cur_edge->verts[0] == face->verts[2] && cur_edge->verts[1] == face->verts[1])
		return 3;           ////case 2->1, the inverse to triangle's orientation
	else if(cur_edge->verts[0] == face->verts[2] && cur_edge->verts[1] == face->verts[0])
		return 4;           ////case 2->0, the same as triangle's orientation
	else if(cur_edge->verts[0] == face->verts[0] && cur_edge->verts[1] == face->verts[2])
		return 5;           ////case 0->2, the inverse to triangle's orientation
	else 
		return -1;          ////wrong
}


////Calculate the region normal at a specific edge
void CalNormalAtEdge(Edge *cur_edge, Face *face, int type)
{
	icVector2 normal;
	int edgecase;
	int orientation;

	////Judge the orientation of the triangle, we need it to get the correct normal
	if(face->xy[2][1] > 0) ////It is counter clockwise orientation
	{
		orientation = 0;
	}
	else{                  ////It is clockwise orientation
		orientation = 1;
	}

	/*if the face is boundary face, we don't call following routine*/

	edgecase = EdgeVerticesCase(cur_edge, face);
	////calculate the outward normal for the edge
	if(orientation == 0) ////counter clockwise orientation
	{
		switch(edgecase)
		{
		case 0:   ////the same as triangle's orientation
		case 2:
		case 4:   normal.entry[0] = Object.vlist[cur_edge->verts[1]]->y - Object.vlist[cur_edge->verts[0]]->y;
			normal.entry[1] = -(Object.vlist[cur_edge->verts[1]]->x - Object.vlist[cur_edge->verts[0]]->x);
			normalize(normal);
			break;

		case 1:   ////the inverse orientation
		case 3:
		case 5:   normal.entry[0] = -(Object.vlist[cur_edge->verts[1]]->y - Object.vlist[cur_edge->verts[0]]->y);
			normal.entry[1] = Object.vlist[cur_edge->verts[1]]->x - Object.vlist[cur_edge->verts[0]]->x;
			normalize(normal);
			break;

		default:
			MessageBox(NULL, "wrong edge!", "Error", MB_OK); 
		}
	}

	else{
		switch(edgecase)
		{
		case 0:   ////the same as triangle's orientation
		case 2:
		case 4:  normal.entry[0] = -(Object.vlist[cur_edge->verts[1]]->y - Object.vlist[cur_edge->verts[0]]->y);
			normal.entry[1] = Object.vlist[cur_edge->verts[1]]->x - Object.vlist[cur_edge->verts[0]]->x;
			normalize(normal);
			break; 

		case 1:   ////the inverse orientation
		case 3:
		case 5:   normal.entry[0] = Object.vlist[cur_edge->verts[1]]->y - Object.vlist[cur_edge->verts[0]]->y;
			normal.entry[1] = -(Object.vlist[cur_edge->verts[1]]->x - Object.vlist[cur_edge->verts[0]]->x);
			normalize(normal);
			break;

		default:
			MessageBox(NULL, "wrong edge!", "Error", MB_OK); 
		}
	}

	////We need to store the normal for separatrix editing and limit cycle design
	if(type == 0) //repeller, calculate the outward normal
	{
		cur_edge->normal = normal;
	}
	else  //attractor,
	{
		cur_edge->normal = -normal;  //separatrix editing will need this
	}

	////Store to the edge structure
	if(type == 0)
		cur_edge->repell_normal = normal;
	else if(type == 1)
		cur_edge->attract_normal = normal;
	else{
		cur_edge->attract_normal = cur_edge->repell_normal = normal;
	}
}

////Get the normal for current edges on the boundary
////This routine should be called after we build the new boundary list
void GetRegionNormals(int type)
{
	int i;
	Edge *cur_edge;
	Face *face;
	int EdgeNum;
	Edge **edgelist;

	if(type == 0)
	{
		EdgeNum = repellerBoundary.num;
		edgelist = repellerBoundary.edgelist;
	}
	else if(type == 1){
		EdgeNum = attractorBoundary.num;
		edgelist = attractorBoundary.edgelist;
	}
	else{
		EdgeNum = intersectBoundary.num;
		edgelist = intersectBoundary.edgelist;
	}


	for(i = 0; i < EdgeNum; i ++)
	{
		cur_edge = edgelist[i];

		////Get the face that inside the region
		face = Object.flist[cur_edge->tris[0]];

		if(type == 0){
    		if(face->repell_inregion == 0)
				face = Object.flist[cur_edge->tris[1]];
		}
		else{
			if(face->attract_inregion == 0)
				face = Object.flist[cur_edge->tris[1]];
		}

		CalNormalAtEdge(cur_edge, face, type);
	}
}


////build the edges list for the region boundary
void UpdateBoundary(int type)
{
	int i, j;
	Face *face, *face2;
	Edge *cur_edge;

	int TriangleNum;
	int *trianglelist;

	if(type == 0){
		TriangleNum = repellerRegion.num;
		trianglelist = repellerRegion.trianglelist;
		repellerBoundary.num = 0; ////Reset the boundary list
	}
	else if(type == 1){
		TriangleNum = attractorRegion.num;
		trianglelist = attractorRegion.trianglelist;
		attractorBoundary.num = 0;
	}

	else{   ////calculate the intersected region
		TriangleNum = intersectRegion.num;
		trianglelist = intersectRegion.trianglelist;
		intersectBoundary.num = 0;
	}


	////Initialization part

	////1. Reset the edges' flags
	for(i = 0; i < TriangleNum; i++)
	{
		face = Object.flist[trianglelist[i]];

		for(j = 0; j < 3; j++)
		{
			cur_edge = face->edges[j];
			if(cur_edge == NULL)
				continue;
			cur_edge->attract_visited = cur_edge->repell_visited = 0;
		}
	}

	////2. Reset region flags
	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];
		if(type == 0)
			face->repell_inregion = 0;
		else
		    face->attract_inregion = 0;
	}

	for(i = 0; i < TriangleNum; i++)
	{
		face = Object.flist[trianglelist[i]];
		if(type == 0)
			face->repell_inregion = 1;
		else
			face->attract_inregion = 1;
	}

	//////////////////////////////////////////////
	for(i = 0; i < TriangleNum; i++)
	{
		face = Object.flist[trianglelist[i]];
	

		for(j = 0; j < 3; j++)
		{
			cur_edge = face->edges[j];

			if(cur_edge == NULL || cur_edge->index < 0 || cur_edge->index > TotalEdgesNum)
				continue;

			////testing the boundary of the whole field  

			//We need to let the boundary edge be included in the region boundary
			if( cur_edge->tris[0] < 0 || cur_edge->tris[1] < 0)
			{
				AddToBoundaryEdges(cur_edge, type);
				continue;
			}

			if(type == 0)
			{
				if(cur_edge->repell_visited == 1)
					continue;
			}
			else{
				if(cur_edge->attract_visited == 1)
					continue;
			}

			////Test the two faces that share current edge
			if(cur_edge->tris[0] != trianglelist[i])
				face2 = Object.flist[cur_edge->tris[0]];
			else
				face2 = Object.flist[cur_edge->tris[1]];


			if(type == 0)
			{
				cur_edge->repell_visited = 1;
				if(face2->repell_inregion == 1) ////The other triangle is inside the region, this is not on boundary
					continue;
				else{
					////Add to the boundary list
					AddToBoundaryEdges(cur_edge, type);
				}
			}

			else{
				cur_edge->attract_visited = 1;
				if(face2->attract_inregion == 1)////The other triangle is inside the region, this is not on boundary
					continue;
				else{
					////Add to the boundary list
					AddToBoundaryEdges(cur_edge, type);
				}
			}
		}//for
	}//for
}

bool RepellerExitEdgePending(Edge *cur_edge)
{
	double dot1, dot2, dotresult = 0;
	Vertex *v1, *v2;
	v1 = Object.vlist[cur_edge->verts[0]];
	v2 = Object.vlist[cur_edge->verts[1]];

	dot1 = dot(cur_edge->repell_normal, v1->vec);
	dot2 = dot(cur_edge->repell_normal, v2->vec);

	dotresult = max(dot1, dot2);   //growing with mixed edges
	//dotresult = min(dot1, dot2); //pure exit edge growing

	//if(dotresult >= 0) return true;
	if(dotresult > 0) return true;

	return false;
}

bool AttractorExitEdgePending(Edge *cur_edge)
{
	double dot1, dot2, dotresult = 0;
	Vertex *v1, *v2;
	v1 = Object.vlist[cur_edge->verts[0]];
	v2 = Object.vlist[cur_edge->verts[1]];

	dot1 = dot(cur_edge->attract_normal, v1->vec);
	dot2 = dot(cur_edge->attract_normal, v2->vec);

	dotresult = min(dot1, dot2);   //growing with mixed edges
	//dotresult = max(dot1, dot2); //pure entrance edge growing

	//if(dotresult <= 0) return true;
	if(dotresult < 0) return true;

	return false;
}

/******************************************************************************/
/*************************************************************/
////Old version for the region growing 

/********************************************************
The main routine for region growing 
********************************************************/
void Cancel_Growing(int type, int target_triangle)
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

			if(oppositeTriangle < 0)
				continue;

			if(Object.flist[oppositeTriangle]->contain_singularity == 1) ////if the opposite triangle containing singularity
			{
					////If the singularity is the intermediate saddle point, add the triangle to the region
					if(IstheMediaNode(Object.flist[oppositeTriangle]->singularity_index)
						||oppositeTriangle == target_triangle)
					{
						////add the triangle to the region
						AddToRegionTriangles(oppositeTriangle, type);
						num_newaddedtriangle++;
					}

					else
						continue;
			}

			if(/*Object.flist[oppositeTriangle]->contain_separatrix == 1 &&*/
				Object.flist[oppositeTriangle]->fence_flag == 1) ////set fence 2/16/06
				continue;

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


/*
Grow region for periodic orbit detection 4/15/06
*/

bool PeriodicDetect_Growing(int type)
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

			if(oppositeTriangle < 0) //reach the mesh boundary
				return false;  

			if(Object.flist[oppositeTriangle]->contain_singularity == 1) ////if the opposite triangle containing singularity
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
			return true;

		////Update the boundary list under new region
		UpdateBoundary(type);
		GetRegionNormals(type);


		if(type == 0)
			num_edges = repellerBoundary.num;
		else
			num_edges = attractorBoundary.num;
	}

	return true;
}


/*
New region growing for limit cycle detection.
The modifications are that we return the reason of the region growing stop
*/
bool PeriodicDetect_Growing(int type, int centersing, int &flag)
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

			if(oppositeTriangle < 0) //reach the mesh boundary
			{
				flag = 1;            //set the flag
				continue;
			}

			if(Object.flist[oppositeTriangle]->contain_singularity == 1) ////if the opposite triangle containing singularity
			{
				////if this is a intermediate saddle, add it! changed at 06/18/06!
				if(singularities[Object.flist[oppositeTriangle]->singularity_index].type != SADDLE)
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
		{
			flag = 0;
			return true;
		}

		////Update the boundary list under new region
		UpdateBoundary(type);
		GetRegionNormals(type);

		if(type == 0)
			num_edges = repellerBoundary.num;
		else
			num_edges = attractorBoundary.num;
	}

	flag = 0;
	return true;
}


/*************************************************************
Grow region for a repeller
*************************************************************/

//void Cancel_GrowRepellerRegion(int which_triangle, int target_triangle, int singularID, int initsaddlelength)
//{
//	////initialization part
//
//	////Intialize the region related flags
//	int i;
//	for(i = 0; i < Object.nfaces; i++)
//	{
//		Object.flist[i]->repell_inregion = 0;
//	}
//
//	for(i = 0; i < Object.nverts; i++)
//	{
//		Object.vlist[i]->repell_flag = 0;
//	}
//
//	////initialize repeller region with first input triangle
//	repellerRegion.trianglelist[0] = which_triangle;
//	Object.flist[which_triangle]->repell_inregion = 1;
//	repellerRegion.num = 1;
//
//	////Initialize a larger region for saddle
//	if(singularities[singularID].type == SADDLE)
//	{
//		Cancel_InitSaddleRegion(singularID, 1, initsaddlelength);
//		UpdateBoundary(0);
//	}
//	else{
//		////initialize repeller boundary with the 3 edges of first triangle
//		repellerBoundary.edgelist[0] = Object.flist[which_triangle]->edges[0];
//		repellerBoundary.edgelist[1] = Object.flist[which_triangle]->edges[1];
//		repellerBoundary.edgelist[2] = Object.flist[which_triangle]->edges[2];
//		repellerBoundary.num = 3;
//	}
//
//	GetRegionNormals(0); ////Get the outward normal of the region for each boundary edge
//
//	////initialize the inner vertices list
//	repellerInnerverts.num = 0;    ////no inner vertex at present
//
//	////region grow processing
//	Cancel_Growing(0, target_triangle);
//}
//
//
/*************************************************************
Grow region for attractor
*************************************************************/
//void Cancel_GrowAttractorRegion(int which_triangle, int target_triangle, int singularID, int initsaddlelength)
//{
//	////initialization part
//
//	////Intialize the region related flags
//	int i;
//	Face *face;
//
//	for(i = 0; i < Object.nfaces; i++)
//	{
//		face = Object.flist[i];
//		face->attract_inregion = 0;
//	}
//
//
//	////initialize attractor region with first input triangle
//	attractorRegion.trianglelist[0] = which_triangle;
//	Object.flist[which_triangle]->attract_inregion = 1;
//	attractorRegion.num = 1;
//
//	////Initialize a larger region for saddle
//	if(singularities[singularID].type == SADDLE)
//	{
//		Cancel_InitSaddleRegion(singularID, 0, initsaddlelength);
//		UpdateBoundary(1);
//	}
//	else{
//		////initialize attractor boundary with the 3 edges of first triangle
//		attractorBoundary.edgelist[0] = Object.flist[which_triangle]->edges[0];
//		attractorBoundary.edgelist[1] = Object.flist[which_triangle]->edges[1];
//		attractorBoundary.edgelist[2] = Object.flist[which_triangle]->edges[2];
//		attractorBoundary.num = 3;
//	}
//
//	GetRegionNormals(1); ////Get the outward normal of the region for each boundary edge
//
//	////initialize the inner vertices list
//	attractorInnerverts.num = 0;    ////no inner vertex at present
//
//	////region grow processing
//	Cancel_Growing(1, target_triangle);
//}



/*************************************************************
Grow region for attractor
*************************************************************/
void Cancel_GrowAttractorRegion(int which_triangle, int target_triangle, int singularID, int initsaddlelength,
								int MultiCancel, double length_percentage)
{
	////initialization part

	////Intialize the region related flags
	int i;
	Face *face;

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];
		face->attract_inregion = 0;
	}


	////initialize attractor region with first input triangle
	attractorRegion.trianglelist[0] = which_triangle;
	Object.flist[which_triangle]->attract_inregion = 1;
	attractorRegion.num = 1;

	////Initialize a larger region for saddle
	if(singularities[singularID].type == SADDLE)
	{
		if(MultiCancel == 0)
			Cancel_InitSaddleRegion(singularID, 0, length_percentage);
		else
		{
			/* -- Only suitable for one indirectly connected pair -- */
			int attract = singularities[singularID].node_index;

			InitSaddleGrowforMultRegion(singularID, 1, length_percentage, &attract, 1);
		}

		UpdateBoundary(1);
	}
	else{
		////initialize attractor boundary with the 3 edges of first triangle
		attractorBoundary.edgelist[0] = Object.flist[which_triangle]->edges[0];
		attractorBoundary.edgelist[1] = Object.flist[which_triangle]->edges[1];
		attractorBoundary.edgelist[2] = Object.flist[which_triangle]->edges[2];
		attractorBoundary.num = 3;

		////
		attractorBoundary.edgelist[0]->OnAttractBoundary = 1;
		attractorBoundary.edgelist[1]->OnAttractBoundary = 1;
		attractorBoundary.edgelist[2]->OnAttractBoundary = 1;
	}
		
	GetRegionNormals(1); ////Get the outward normal of the region for each boundary edge

	////initialize the inner vertices list
	attractorInnerverts.num = 0;    ////no inner vertex at present

	////region grow processing
	Cancel_Growing(1, target_triangle);
}


/*************************************************************
Grow region for a repeller
*************************************************************/

void Cancel_GrowRepellerRegion(int which_triangle, int target_triangle, int singularID, int initsaddlelength,
								int MultiCancel, double length_percentage)
{
	////initialization part

	////Intialize the region related flags
	int i;
	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->repell_inregion = 0;
	}

	for(i = 0; i < Object.nverts; i++)
	{
		Object.vlist[i]->repell_flag = 0;
	}

	////initialize repeller region with first input triangle
	repellerRegion.trianglelist[0] = which_triangle;
	Object.flist[which_triangle]->repell_inregion = 1;
	repellerRegion.num = 1;

	////Initialize a larger region for saddle
	if(singularities[singularID].type == SADDLE)
	{
		if(MultiCancel == 0)
			Cancel_InitSaddleRegion(singularID, 1, length_percentage);
		else
		{
			/* -- Only suitable for one indirectly connected pair -- */
			int repell = singularities[singularID].node_index;

			InitSaddleGrowforMultRegion(singularID, 0, length_percentage, &repell, 1);
		}

		UpdateBoundary(0);
	}
	else{
		////initialize repeller boundary with the 3 edges of first triangle
		repellerBoundary.edgelist[0] = Object.flist[which_triangle]->edges[0];
		repellerBoundary.edgelist[1] = Object.flist[which_triangle]->edges[1];
		repellerBoundary.edgelist[2] = Object.flist[which_triangle]->edges[2];
		repellerBoundary.num = 3;

		repellerBoundary.edgelist[0]->OnRepellBoundary = 1;
		repellerBoundary.edgelist[1]->OnRepellBoundary = 1;
		repellerBoundary.edgelist[2]->OnRepellBoundary = 1;
	}
	
	GetRegionNormals(0); ////Get the outward normal of the region for each boundary edge

	////initialize the inner vertices list
	repellerInnerverts.num = 0;    ////no inner vertex at present

	////region grow processing
	Cancel_Growing(0, target_triangle);
}


/*********************************************************************
 Get the intersected region of repeller region and attractor region
 Note that here we need not worry about the space issue!:)
*********************************************************************/
void Cancel_GetIntersectedRegion(int repellID, int attractID)
{
	int i, j;
	Face *face;
	Edge *cur_edge;
	Vertex *vert;

	////Get the intersect triangular region
	intersectRegion.num = 0;
	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];
		if(face->repell_inregion == 1 && face->attract_inregion == 1)
		{
			intersectRegion.trianglelist[intersectRegion.num] = face->index;
			intersectRegion.num ++;
		}
		else{
			face->repell_inregion = 0;
			face->attract_inregion = 0;
		}
	}

	////Add the triangles that contain the two singularities into the intersectRegion 2/23/06
	if(!IsRepeated(intersectRegion.trianglelist, singularities[repellID].Triangle_ID, intersectRegion.num))
	{
		intersectRegion.trianglelist[intersectRegion.num] = singularities[repellID].Triangle_ID;
		intersectRegion.num ++;
	}

	if(!IsRepeated(intersectRegion.trianglelist, singularities[attractID].Triangle_ID, intersectRegion.num))
	{
		intersectRegion.trianglelist[intersectRegion.num] = singularities[attractID].Triangle_ID;
		intersectRegion.num ++;
	}

	/*----------------------------------------------------------*/
	////Add the intermediary points into the region
	////11/19/05
	//for(i = 0; i < Num_MediaNodes; i++)
	//{
	//	int repeated = 0;
	//	for(j = 0; j < intersectRegion.num; j++)
	//	{
	//		if(singularities[MediaNodes[i]].Triangle_ID == intersectRegion.trianglelist[j])
	//		{
	//			repeated = 1;
	//			break;
	//		}
	//	}

	//	if(repeated == 0)
	//	{
	//		intersectRegion.trianglelist[intersectRegion.num] = singularities[MediaNodes[i]].Triangle_ID;
	//		intersectRegion.num ++;
	//	}
	//}
	/*-----------------------------------------------------------*/

	////Reuse the UpdateBoundary routine to build the boundary edge list
	////First, we need to reset the flag of edge visting
	for(i = 0; i < intersectRegion.num; i++)
	{
		face = Object.flist[intersectRegion.trianglelist[i]];
		for(j = 0; j < 3; j++)
		{
			cur_edge =  face->edges[j];
			cur_edge->repell_visited = 0;
			cur_edge->attract_visited = 0;
		}
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
	for(i = 0; i < Object.nverts; i++)
	{
		vert = Object.vlist[i];
		if(vert->attract_flag == 2 && vert->repell_flag == 2 && vert->OnBoundary == 0) ////The vertex is inside both region
		{
			intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
			vert->InRegion = 1;
			vert->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
			intersectInnerverts.num ++;
		}
	}

	////Add the vertices of the triangles that contain the two singularities to the innerverts list
	////Testing codes here
	face = Object.flist[singularities[repellID].Triangle_ID];
	for(j = 0; j < face->nverts; j++)
	{
		vert = Object.vlist[face->verts[j]];

		if(vert->InRegion != 1)
		{
			intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
			vert->InRegion = 1;
			vert->OnBoundary = 0;
			vert->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
			intersectInnerverts.num ++;
		}
	}
	
	face = Object.flist[singularities[attractID].Triangle_ID];
	for(j = 0; j < face->nverts; j++)
	{
		vert = Object.vlist[face->verts[j]];

		if(vert->InRegion != 1)
		{
			intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
			vert->InRegion = 1;
			vert->OnBoundary = 0;
			vert->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
			intersectInnerverts.num ++;
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////
}


////What version it is?
////we may get the final intersected region here
//void Cancel_GetRegion(int repell_triangle, int attract_triangle, int repellID, int attractID)
//{
//	////first, grow the repeller region
//	Cancel_GrowRepellerRegion(repell_triangle, attract_triangle, repellID, InitSaddleRegionLength);
//
//	////second, grow the attractor region
//	Cancel_GrowAttractorRegion(attract_triangle, repell_triangle, attractID, InitSaddleRegionLength);
//
//	////third, calculate the intersection of the two regions obtained above
//	Cancel_GetIntersectedRegion(repellID, attractID);
//}


////Adaptive pair cancellation for sink and saddle or source and saddle
bool AdaptivePairCancel(int repell_triangle, int attract_triangle, int repellID, int attractID)
{
	int i;
	double cur_length, small_length, large_length;
	int success = 0;

	int traid[2], min_triangles;

	small_length = 0;
	large_length = 2.;

	if(singularities[attractID].type == SADDLE) ////saddle acts as an attractor
	{
		traid[0] = separatrices[singularities[attractID].separtices].sep1;
		traid[1] = separatrices[singularities[attractID].separtices].sep3;
	}
	else                                        ////saddle acts as a repeller
	{
		traid[0] = separatrices[singularities[repellID].separtices].sep2;
		traid[1] = separatrices[singularities[repellID].separtices].sep4;
	}

	min_triangles = min(num_linesegs_curtraj[traid[0]]-1, num_linesegs_curtraj[traid[1]]-1);

    //// Loop until the difference between current saddle initial length and previous one is equal to 1
	int count = 0;
	do{
		InitCancellationAndMovement();

		cur_length = /*(int)*/(small_length + large_length)/2;
			
		//Grow region according to current length
		if(singularities[attractID].type == SADDLE) ////saddle acts as an attractor
		{
			//Cancel_GrowRepellerRegion(repell_triangle, attract_triangle, repellID, cur_length, 0, 1.);
			//Cancel_GrowAttractorRegion(attract_triangle, repell_triangle, attractID, cur_length, 0, 1.);
			Cancel_GrowRepellerRegion(repell_triangle, attract_triangle, repellID, cur_length, 0, cur_length);
			Cancel_GrowAttractorRegion(attract_triangle, repell_triangle, attractID, cur_length, 0, cur_length);
		}

		else{
			//Cancel_GrowAttractorRegion(attract_triangle, repell_triangle, attractID, cur_length, 0, 1.);
			//Cancel_GrowRepellerRegion(repell_triangle, attract_triangle, repellID, cur_length, 0, 1.);
			Cancel_GrowAttractorRegion(attract_triangle, repell_triangle, attractID, cur_length, 0, cur_length);
			Cancel_GrowRepellerRegion(repell_triangle, attract_triangle, repellID, cur_length, 0, cur_length);
		}

		////3c) Get intersected region of the two region
		Cancel_GetIntersectedRegion(repellID, attractID);

		////Perform previous conditions testing before smoothing
		////The euler number should be decided by the two components we want to delete!!!
		////which we call them the conley boundary condition (Homology) 06/09/06
		if(CalEulerValue(intersectRegion.trianglelist, intersectRegion.num) == 1)
		{
			success = 1;
			if(count == 0)
			{
				break;
			}
			else
				small_length = cur_length;  //get a bigger region
		}

		else
		{
			large_length = cur_length;      //get a smaller region
		}

		count++;

	} while(large_length - small_length > 1e-5 && (double)min_triangles*(large_length - small_length) > 0.5
		&& (int)min_triangles*cur_length >= 1);

	if(success == 1)
	{
		////perform smoothing on the intersected region
		Cancel_RegionSmooth();

		if(!IntersectedRegionSingCapture())
			return true;
		else
		{
			Undo();
			return false;
		}
	}

	else
		return false;

}


/********************************************************************************/
/********************************************************************************/
////New version routines for region smoothing based on Conley index theory

/*---------------------------------------------------------------------
Type: 0 -- repeller;  1 -- attractor
target_triangle
---------------------------------------------------------------------*/
void Cancel_Growing_adv(int type, int target_triangle)
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

			if(Object.flist[oppositeTriangle]->contain_singularity == 1) ////if the opposite triangle containing singularity
			{
				if(IstheMediaNode(Object.flist[oppositeTriangle]->singularity_index)
					/*||oppositeTriangle == target_triangle*/)
				{
					////add the triangle to the region
					AddToRegionTriangles(oppositeTriangle, type);
					num_newaddedtriangle++;
				}

				else //// it is not the triangle containing the counterpart singularity
					continue;
			}

							
			if(Object.flist[oppositeTriangle]->contain_separatrix == 1&&
				Object.flist[oppositeTriangle]->fence_flag == 1) ////set fence 2/16/06
				continue;


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

/*************************************************************
Grow region for a repeller
*************************************************************/

void Cancel_GrowRepellerRegion_adv(int which_triangle, int target_triangle, int singularID,
								   double length_percentage)
{
	////initialization part

	////Intialize the region related flags
	int i;

	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->repell_inregion = 0;
	}

	for(i = 0; i < Object.nverts; i++)
	{
		Object.vlist[i]->repell_flag = 0;
	}

	////initialize repeller region with first input triangle
	repellerRegion.trianglelist[0] = which_triangle;
	Object.flist[which_triangle]->repell_inregion = 1;
	repellerRegion.num = 1;

	////Initialize a larger region for saddle
	if(singularities[singularID].type == SADDLE)
	{
		InitSaddleGrowforMultRegion_new(singularID, 0, 600, length_percentage,
			attractor_nodes, NumAttractors);
		UpdateBoundary(0);
	}
	else{
		////initialize repeller boundary with the 3 edges of first triangle
		repellerBoundary.edgelist[0] = Object.flist[which_triangle]->edges[0];
		repellerBoundary.edgelist[1] = Object.flist[which_triangle]->edges[1];
		repellerBoundary.edgelist[2] = Object.flist[which_triangle]->edges[2];
		repellerBoundary.num = 3;
	}

	GetRegionNormals(0); ////Get the outward normal of the region for each boundary edge

	////initialize the inner vertices list
	repellerInnerverts.num = 0;    ////no inner vertex at present

	////region grow processing
	Cancel_Growing_adv(0, target_triangle);
}


/*************************************************************
Grow region for attractor
*************************************************************/
void Cancel_GrowAttractorRegion_adv(int which_triangle, int target_triangle, int singularID, 
									double length_percentage)
{
	////initialization part

	////Intialize the region related flags
	int i;
	Face *face;

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];
		face->attract_inregion = 0;
	}


	////initialize attractor region with first input triangle
	attractorRegion.trianglelist[0] = which_triangle;
	Object.flist[which_triangle]->attract_inregion = 1;
	attractorRegion.num = 1;

	////Initialize a larger region for saddle
	if(singularities[singularID].type == SADDLE)
	{
		InitSaddleGrowforMultRegion_new(singularID, 1,600 , length_percentage,
			attractor_nodes, NumAttractors);
		UpdateBoundary(1);
	}
	else{
		////initialize attractor boundary with the 3 edges of first triangle
		attractorBoundary.edgelist[0] = Object.flist[which_triangle]->edges[0];
		attractorBoundary.edgelist[1] = Object.flist[which_triangle]->edges[1];
		attractorBoundary.edgelist[2] = Object.flist[which_triangle]->edges[2];
		attractorBoundary.num = 3;
	}

	GetRegionNormals(1); ////Get the outward normal of the region for each boundary edge

	////initialize the inner vertices list
	attractorInnerverts.num = 0;    ////no inner vertex at present

	////region grow processing
	Cancel_Growing_adv(1, target_triangle);
}


void Cancel_GetIntersectedRegion_adv(int repellID, int attractID)
{
	int i, j;
	Face *face;
	Edge *cur_edge;
	Vertex *vert;

	////Get the intersect triangular region
	intersectRegion.num = 0;
	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];
		if(face->repell_inregion == 1 && face->attract_inregion == 1)
		{
			intersectRegion.trianglelist[intersectRegion.num] = face->index;
			intersectRegion.num ++;
		}
		else{
			face->repell_inregion = 0;
			face->attract_inregion = 0;
		}
	}

	/*----------------------------------------------------------*/
	////Add the intermediary points into the region
	////11/19/05
	for(i = 0; i < Num_MediaNodes; i++)
	{
		int repeated = 0;
		for(j = 0; j < intersectRegion.num; j++)
		{
			if(singularities[MediaNodes[i]].Triangle_ID == intersectRegion.trianglelist[j])
			{
				repeated = 1;
				break;
			}
		}

		if(repeated == 0)
		{
			intersectRegion.trianglelist[intersectRegion.num] = singularities[MediaNodes[i]].Triangle_ID;
			intersectRegion.num ++;
		}
	}
	/*-----------------------------------------------------------*/

	////Reuse the UpdateBoundary routine to build the boundary edges' list
	////First, we need to reset the flag of edge visting
	for(i = 0; i < intersectRegion.num; i++)
	{
		face = Object.flist[intersectRegion.trianglelist[i]];
		for(j = 0; j < 3; j++)
		{
			cur_edge =  face->edges[j];
			cur_edge->repell_visited = 0;
			cur_edge->attract_visited = 0;
		}
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
	
	/*----------------------------11/19/05-------------------------------------*/
	//// add the two median saddles into the region if there are not in the region now
	for(i = 0; i < Num_MediaNodes; i++)
	{
		face = Object.flist[singularities[MediaNodes[i]].Triangle_ID];
		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];
			vert->OnBoundary = 0;
			vert->attract_flag = 2;
			vert->repell_flag = 2;
		}
	}
	/*------------------------------------------------------------------------*/

	////Get the intersection of inner vertices
	intersectInnerverts.num = 0;
	for(i = 0; i < Object.nverts; i++)
	{
		vert = Object.vlist[i];
		if(vert->attract_flag == 2 && vert->repell_flag == 2 /*&& vert->OnBoundary == 0*/) ////The vertex is inside both region
		{
			intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
			vert->InRegion = 1;
			vert->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
			intersectInnerverts.num ++;
		}
	}

	////Add the vertices of the triangles that contain the two singularities to the innerverts list
	////Testing codes here
	face = Object.flist[singularities[repellID].Triangle_ID];
	for(j = 0; j < face->nverts; j++)
	{
		vert = Object.vlist[face->verts[j]];

		if(vert->InRegion != 1)
		{
			intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
			vert->InRegion = 1;
			vert->OnBoundary = 0;
			vert->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
			intersectInnerverts.num ++;
		}
	}
	
	face = Object.flist[singularities[attractID].Triangle_ID];
	for(j = 0; j < face->nverts; j++)
	{
		vert = Object.vlist[face->verts[j]];

		if(vert->InRegion != 1)
		{
			intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
			vert->InRegion = 1;
			vert->OnBoundary = 0;
			vert->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
			intersectInnerverts.num ++;
		}
	}

}

/*----------------------------------------------------------------------------
Adaptive pair cancellation for two source and sink cancellation
----------------------------------------------------------------------------*/
//void AdaptivePairCancel_adv(int repeller, int attractor, int saddle1, int saddle2)
//{
//	int i, j;
//	////Do not use adaptive region adjustment here now!!! 11/06/05
//
//	////1. get the intersected region (repeller ^ saddle1) and (repeller ^ saddle2)
//	
//	////1a) treat repeller + saddle1 as the pair that we want to cancel
//	////Now, saddle1 acts as an attractor
//	Cancel_GrowRepellerRegion(singularities[repeller].Triangle_ID, singularities[saddle1].Triangle_ID,\
//		repeller, 200);
//	Cancel_GrowAttractorRegion(singularities[saddle1].Triangle_ID, singularities[repeller].Triangle_ID,\
//		saddle1, 200);
//	Cancel_GetIntersectedRegion(repeller, saddle1);
//
//	////Copy the region from the intersection region to Source_re1
//	CopyRegion(intersectRegion.trianglelist, Source_re1.trianglelist, intersectRegion.num);
//	Source_re1.num = intersectRegion.num;
//	////Reset all the flag for next pair region building
//    InitCancellationAndMovement();
//
//	/**********************************************************************/
//	
//	////1b) treat repeller + saddle2 as the pair that we want to cancel
//	Cancel_GrowRepellerRegion(singularities[repeller].Triangle_ID, singularities[saddle2].Triangle_ID,\
//		repeller, 200);
//	Cancel_GrowAttractorRegion(singularities[saddle2].Triangle_ID, singularities[repeller].Triangle_ID,\
//		saddle2, 200);
//	Cancel_GetIntersectedRegion(repeller, saddle2);
//
//	////Copy the region from the intersection region to Source_re2
//	CopyRegion(intersectRegion.trianglelist, Source_re2.trianglelist, intersectRegion.num);
//	Source_re2.num = intersectRegion.num;
//
//	////Reset all the flag for next pair region building
//    InitCancellationAndMovement();
//
//	/*-----------------------------------------------------------*/
//
//	////2. Get the union of these two regions = source_re
//    UnionRegion(Source_re1, Source_re2, Source_re);
//
//	/*-----------------------------------------------------------*/
//	////3. Get the intersected region (attractor ^ saddle1) and (attractor ^ saddle2)
//
//	////3a) treat attractor + saddle1 as the pair that we want to cancel
//	////Now, saddle1 acts as a repeller
//	Cancel_GrowAttractorRegion(singularities[attractor].Triangle_ID, singularities[saddle1].Triangle_ID,\
//		attractor, 200);
//	Cancel_GrowRepellerRegion(singularities[saddle1].Triangle_ID, singularities[attractor].Triangle_ID,\
//		saddle1, 200);
//	Cancel_GetIntersectedRegion(saddle1, attractor);
//
//	////Copy the region from the intersection region to Sink_re1
//	CopyRegion(intersectRegion.trianglelist, Sink_re1.trianglelist, intersectRegion.num);
//	Sink_re1.num = intersectRegion.num;
//
//	////Reset all the flag for next pair region building
//    InitCancellationAndMovement();
//
//	/**********************************************************************/
//	
//	////3b) treat attractor + saddle2 as the pair that we want to cancel
//	Cancel_GrowAttractorRegion(singularities[attractor].Triangle_ID, singularities[saddle2].Triangle_ID,\
//		attractor, 200);
//	Cancel_GrowRepellerRegion(singularities[saddle2].Triangle_ID, singularities[attractor].Triangle_ID,\
//		saddle2, 200);
//	Cancel_GetIntersectedRegion(saddle2, attractor);
//
//	////Copy the region from the intersection region to Sink_re2
//	CopyRegion(intersectRegion.trianglelist, Sink_re2.trianglelist, intersectRegion.num);
//	Sink_re2.num = intersectRegion.num;
//
//	////Reset all the flag for next pair region building
//    InitCancellationAndMovement();
//
//	/*-----------------------------------------------------------*/
//	////4. Get the union of these two regions = sink_re
//    UnionRegion(Sink_re1, Sink_re2, Sink_re);
//
//	/*-----------------------------------------------------------*/
//	////5. Get the final smoothing region source_re U sink_re
//	UnionRegion(Source_re, Sink_re, intersectRegion);
//	
//	////It seems the region is too large
//	Face *cur_f;
//	Vertex * cur_v;
//	for(i = 0; i < intersectRegion.num; i++)
//	{
//		cur_f = Object.flist[intersectRegion.trianglelist[i]];
//
//		for(j = 0; j < cur_f->nverts; j++)
//		{
//			cur_v = Object.vlist[cur_f->verts[j]];
//			cur_v->repell_flag = 2;
//			cur_v->attract_flag = 2;
//		}
//	}
//
//	UpdateBoundary(2);
//
//
//	////Set the vertices on the boundary as 'OnBoundary'
//	for(int i = 0; i < intersectBoundary.num; i++)
//	{
//		Object.vlist[intersectBoundary.edgelist[i]->verts[0]]->OnBoundary = 1;
//		Object.vlist[intersectBoundary.edgelist[i]->verts[0]]->InRegion = 0;
//		Object.vlist[intersectBoundary.edgelist[i]->verts[1]]->OnBoundary = 1;
//		Object.vlist[intersectBoundary.edgelist[i]->verts[1]]->InRegion = 0;
//	}
//
//	////Get the intersection of inner vertices
//	intersectInnerverts.num = 0;
//	for(i = 0; i < Object.nverts; i++)
//	{
//		cur_v = Object.vlist[i];
//		if(cur_v->attract_flag == 2 && cur_v->repell_flag == 2 && cur_v->OnBoundary == 0) ////The vertex is inside both region
//		{
//			intersectInnerverts.vertslist[intersectInnerverts.num] = cur_v;
//			cur_v->InRegion = 1;
//			cur_v->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
//			intersectInnerverts.num ++;
//		}
//	}
//
//	Cancel_RegionSmooth();
//}
//
//
/*----------------------------------------------------------------------------
Adaptive pair cancellation for two source and sink cancellation
Add On 11/19/05
----------------------------------------------------------------------------*/
//void AdaptivePairCancel_adv2(int repeller, int attractor/*, int saddle1, int saddle2*/)
//{
//	//1. first, grow regions from repeller (for example, a source)
//	Cancel_GrowRepellerRegion(singularities[repeller].Triangle_ID,
//		singularities[attractor].Triangle_ID, repeller, 200);
//	Cancel_GrowRepellerRegion(singularities[repeller].Triangle_ID,
//		-1, repeller, 200);
//	Cancel_GrowRepellerRegion_adv(singularities[repeller].Triangle_ID,
//		-1, repeller, 200);
//
//	//2. second, grow region from attractor (for example, a sink)
//	Cancel_GrowAttractorRegion(singularities[attractor].Triangle_ID,
//		singularities[repeller].Triangle_ID, attractor, 200);
//	Cancel_GrowAttractorRegion(singularities[attractor].Triangle_ID,
//		-1, attractor, 200);
//	Cancel_GrowAttractorRegion_adv(singularities[attractor].Triangle_ID,
//		-1, attractor, 200);
//
//	//3. third, calculate the intersection of the two regions
//	Cancel_GetIntersectedRegion(repeller, attractor);
//	Cancel_GetIntersectedRegion_adv(repeller, attractor);
//
//	int i, j;
//	Face *face;
//	Vertex *vert;
//	////4. fourth, add the two median saddles into the region if there are not in the region now
//	for(i = 0 i < Num_MediaNodes; i++)
//	{
//		face = Object.flist[singularities[MediaNodes[i]].Triangle_ID];
//		for(j = 0; j < face->nverts; j++)
//		{
//			vert = Object.vlist[face->verts[j]];
//		}
//	}
//
//	//5. fifth, we may need to test the validation of the region before smoothing
//
//	//6. sixth, smooth the region
//	Cancel_RegionSmooth();
//
//	if(IntersectedRegionSingCapture()) ////contain singularity
//	{
//		MessageBox(NULL, "Can not cancel this pair of singularities", "Exception", MB_OK);
//		Undo();
//	}
//
//}
//









/*
Add at 11/22/05 Can be reused at many places
*/
bool IsRepeated(int *a, int b, int num)
{
	for(int i = 0; i < num; i++)
	{
		if(a[i] == b)
			return true;
	}

	return false;
}


/*
Set fence for the scenario of one saddle+one source or one saddle+one sink cancellation
2/16/06
*/

void SetFenceForASep(int index)
{
	int i;

	for(i = 0; i < num_linesegs_curtraj[index]; i++)
	{
		Object.flist[trajectories[index][i].Triangle_ID]->contain_separatrix = 1;
		Object.flist[trajectories[index][i].Triangle_ID]->fence_flag = 1;
	}
}

void SetFenceForOneSaddleandOtherSingCancel(int saddleNode, int singNode)
{
	////we want to set fence for the separatrix that does not connect to sing
	////and is not the region growing intial separatrix

	int i, j, sep, traj;
	sep = singularities[graphnodes[saddleNode].singularityID].separtices;
	int singID = graphnodes[singNode].singularityID;

	if(graphnodes[singNode].type == 0) //if the singularity is a repeller
	{
		////the region growing from the saddle should begin from sep1, sep3
		////so we just need to check sep2, and sep4
		traj = separatrices[sep].sep2;

		if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID == 
			singularities[graphnodes[singNode].singularityID].Triangle_ID)
		{
			//then, we need to set fence around the separatrix sep4
			traj = separatrices[sep].sep4;

		}
		else
		{
			//we need to set fence around the separatrix sep2
			traj = separatrices[sep].sep2;
		}
	}

	else{  //it is an attractor
		////the region growing from the saddle should begin from sep2, sep4
		////so we just need to check sep1, and sep3
		traj = separatrices[sep].sep1;

		if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID == 
			singularities[graphnodes[singNode].singularityID].Triangle_ID)
		{
			//then, we need to set fence around the separatrix sep4
			traj = separatrices[sep].sep3;

		}
		else
		{
			//we need to set fence around the separatrix sep2
			traj = separatrices[sep].sep1;
		}
	}
		
	////For two connection orbit cases
	if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID == 
		singularities[graphnodes[singNode].singularityID].Triangle_ID)
		traj = -1;

	////1. set fences for the separatrix of the input saddles
	if(traj > 0)
		SetFenceForASep(traj);

	////2. set fences for other separatrices

	for(i = 0; i < cur_separatrices_index; i++)
	{
		if(i == sep) continue;

		for(j = 0; j < 4; j++)
		{
			switch(j){
				case 0:
					traj = separatrices[i].sep1;
					break;
				case 1:
					traj = separatrices[i].sep2;
					break;
				case 2:
					traj = separatrices[i].sep3;
					break;
				case 3:
					traj = separatrices[i].sep4;
					break;
			}

			////Test code here 2/23/06
			////we do not set the fences on those seps that connect to the singNode
			if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID ==
				singularities[singID].Triangle_ID)
				continue;
			 
			////set fences for the separatrix of the input saddles
            SetFenceForASep(traj);
		}
	}
}


/*
Set fence for the generic scenario multiple repellers + attractors cancellation
2/16/06
//idea: each of the saddles in the "intervals" list connects to both a repeller
and an attractors in the lists at the same time, so the 4 separatrices associated
with the saddle should be set as fence, maybe the region will be little bit larger ???
*/

bool SepoftheSaddleIntheList(int *intervals, int num_intervals, int sepindex)
{
	int i;
	for(i = 0; i < num_intervals; i++)
	{
		if(sepindex == singularities[graphnodes[intervals[i]].singularityID].separtices)
			return true;
	}
	return false;
}



////We find all the separatrices that involve in the connections with the deleted set 2/16/06

void SetConnectedSeps(int *repellers, int num_repellers,
			   int *attractors, int num_attractors,
			   int *intervals, int num_intervals)
{
	int i, j, k;
	int sep, traj;

	////Reset the flags
	for(i = 0; i < cur_separatrices_index; i++)
	{
		separatrices[i].connect1 = 0;
		separatrices[i].connect2 = 0;
		separatrices[i].connect3 = 0;
		separatrices[i].connect4 = 0;
	}

	////For each intermediary saddle, we test its 4 seps to find whether they connect with
	////one of the element in the deleted set
	for(i = 0; i < num_intervals; i++)
	{
		sep = singularities[graphnodes[intervals[i]].singularityID].separtices;

		for(j = 0; j < 4; j++)
		{
			switch(j){
				case 0:
					traj = separatrices[sep].sep1;
					break;
				case 1:
					traj = separatrices[sep].sep2;
					break;
				case 2:
					traj = separatrices[sep].sep3;
					break;
				case 3:
					traj = separatrices[sep].sep4;
					break;
			}
			
			////Test the set of repellers 
			for(k = 0; k < num_repellers; k++)
			{
				if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID
					== singularities[graphnodes[repellers[k]].singularityID].Triangle_ID)
				{
					switch(j){
						case 0:
							separatrices[sep].connect1 = 1;
							goto LL;
						case 1:
							separatrices[sep].connect2 = 1;
							goto LL;
						case 2:
							separatrices[sep].connect3 = 1;
							goto LL;
						case 3:
							separatrices[sep].connect4 = 1;
							goto LL;
					}
				}
			}

			////Test the set of attractors
			for(k = 0; k < num_attractors; k++)
			{
				if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID
					== singularities[graphnodes[attractors[k]].singularityID].Triangle_ID)
				{
					switch(j){
						case 0:
							separatrices[sep].connect1 = 1;
							goto LL;
						case 1:
							separatrices[sep].connect2 = 1;
							goto LL;
						case 2:
							separatrices[sep].connect3 = 1;
							goto LL;
						case 3:
							separatrices[sep].connect4 = 1;
							goto LL;
					}
				}
			}
LL:         ;
		}

	}
}

////we should set fences at all undirectly connected separatrices 2/16/06
void SetFenceForMultiPairCancel(int *repellers, int num_repellers,
			   int *attractors, int num_attractors,
			   int *intervals, int num_intervals)
{
	//
	int i, j;
	int traj;

	for(i = 0; i < cur_separatrices_index; i++)
	{
		////if it is one the group of the separatrices of the saddle in the interval list
		////do not set fence
		if(SepoftheSaddleIntheList(intervals, num_intervals, i))
			continue;

		////set fence for other separatrices
		for(j = 0; j < 4; j++)
		{
			switch(j){
				case 0:
					traj = separatrices[i].sep1;
					break;
				case 1:
					traj = separatrices[i].sep2;
					break;
				case 2:
					traj = separatrices[i].sep3;
					break;
				case 3:
					traj = separatrices[i].sep4;
					break;
			}

			////set fences for the separatrix of the input saddles
            SetFenceForASep(traj);
		}
	}
}


/*------------------------------------------------------------------
Set the fences of the region  2/16/06
Input: repeller, attractor, intervals-the intermediary saddles, num_intervals
------------------------------------------------------------------*/

void SetFences(int *repellers, int num_repellers,
			   int *attractors, int num_attractors,
			   int *intervals, int num_intervals)
{
	//if there are only two nodes input and one of the repeller and attractor is saddle
	int i;
	int cur_r, cur_a;

	////reset the fence flag
	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->contain_separatrix = 0; ////this may affect the separatrix editing laterly 2/23/06
		Object.flist[i]->fence_flag = 0;
	}

	if(num_repellers == 1 && num_attractors == 1) //saddle + one other singularity scenario
	{
		cur_r = graphnodes[repellers[0]].singularityID;
		cur_a = graphnodes[attractors[0]].singularityID;

		if(singularities[cur_r].type == SADDLE)
			SetFenceForOneSaddleandOtherSingCancel(repellers[0], attractors[0]);

		if(singularities[cur_a].type == SADDLE)
			SetFenceForOneSaddleandOtherSingCancel(attractors[0], repellers[0]);

		////we may also need to set fences for other separatrices !!

		////we need to deal with source-sink scenario, it is also one repeller and one attractor case
	}

	else  //more general scenarios
	{
		SetFenceForMultiPairCancel(repellers, num_repellers,
			attractors, num_attractors,
			intervals, num_intervals);
        
	}
}





/*-------------------------------------------------------------------*/
////11/20/05
////Get region for multiple repellers and attractors cancellation
////This routine only considers the singularity cancellation now!!!
////Modified at 11/30/05 (Adpative length of separatrix)
//Need to be modified according to the surface cases
/*-------------------------------------------------------------------*/
void InitSaddleGrowforMultRegion(int singularID, int type, 
									double length_percentage,
								     int *sametypesings, int numsings)
{
	int i;
    int endsingID;
	int traid_1 = separatrices[singularities[singularID].separtices].sep1;

	//int dividen = pow(2, curtimes);  ////12.20

	if(type == 1)    ////Saddle as an attractor, forward tracing along outgoing direction
	{
		//going = singularities[singularID].outgoing;
		
		////tracing along both directions
		//1. positive outgoing direction

		////Get the singularity at the end of the separatrix
		endsingID = Object.flist[trajectories[traid_1][num_linesegs_curtraj[traid_1]-1].Triangle_ID]->singularity_index;
		if(endsingID >= 0 && IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
		{
			for(i = 0; i < num_linesegs_curtraj[traid_1]; i++)
			{
				////avoid to add the triangle containing singularity
				if(Object.flist[trajectories[traid_1][i].Triangle_ID]->contain_singularity != 1
					|| Object.flist[trajectories[traid_1][i].Triangle_ID]->singularity_index == singularID)
				AddToRegionTriangles(trajectories[traid_1][i].Triangle_ID, 1);
			}
		}

		if((endsingID >= 0 && !IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
			|| endsingID == -1) //endsingID == -1 means the separtrix ends at the boundary or closed orbit
		{
			for(i = 0; i < (int)((num_linesegs_curtraj[traid_1]-1)*length_percentage); i++)
			{
				////avoid to add the triangle containing singularity
				if(Object.flist[trajectories[traid_1][i].Triangle_ID]->contain_singularity != 1
					|| Object.flist[trajectories[traid_1][i].Triangle_ID]->singularity_index == singularID)
				AddToRegionTriangles(trajectories[traid_1][i].Triangle_ID, 1);
			}
		}

		//2. negative outgoing direction
		endsingID = Object.flist[trajectories[traid_1+2][num_linesegs_curtraj[traid_1+2]-1].Triangle_ID]->singularity_index;

		if(endsingID >= 0 && IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
		{
			for(i = 0; i < num_linesegs_curtraj[traid_1+2]; i++)
			{
				////avoid to add the triangle containing singularity
				if(Object.flist[trajectories[traid_1+2][i].Triangle_ID]->contain_singularity != 1
					|| Object.flist[trajectories[traid_1+2][i].Triangle_ID]->singularity_index == singularID)
				AddToRegionTriangles(trajectories[traid_1+2][i].Triangle_ID, 1);
			}
		}

		if((endsingID >= 0 && !IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
			|| endsingID == -1) //endsingID == -1 means the separtrix ends at the boundary or closed orbit
		{
			for(i = 0; i < (int)((num_linesegs_curtraj[traid_1+2]-1)*length_percentage); i++)
			{
				////avoid to add the triangle containing singularity
				if(Object.flist[trajectories[traid_1+2][i].Triangle_ID]->contain_singularity != 1
					|| Object.flist[trajectories[traid_1+2][i].Triangle_ID]->singularity_index == singularID)
				AddToRegionTriangles(trajectories[traid_1+2][i].Triangle_ID, 1);
			}
		}
	}

	/*---------------------------------------------------------------*/
	else{           
		////Saddle as a repeller, backward tracing along incoming direction
		
		////tracing along both directions
		//1. positive incoming direction
		endsingID = Object.flist[trajectories[traid_1+1][num_linesegs_curtraj[traid_1+1]-1].Triangle_ID]->singularity_index;
		if(endsingID >= 0 && IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
		{
			for(i = 0; i < num_linesegs_curtraj[traid_1+1]; i++)
			{
				////avoid to add the triangle containing singularity
				if(Object.flist[trajectories[traid_1+1][i].Triangle_ID]->contain_singularity != 1
					|| Object.flist[trajectories[traid_1+1][i].Triangle_ID]->singularity_index == singularID)
				AddToRegionTriangles(trajectories[traid_1+1][i].Triangle_ID, 0);
			}
		}

		if((endsingID >= 0 && !IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
			|| endsingID == -1) //endsingID == -1 means the separatrix ends at the boundary or closed orbit
		{
			for(i = 0; i < (int)((num_linesegs_curtraj[traid_1+1]-1)*length_percentage); i++)
			{
				////avoid to add the triangle containing singularity
				if(Object.flist[trajectories[traid_1+1][i].Triangle_ID]->contain_singularity != 1
					|| Object.flist[trajectories[traid_1+1][i].Triangle_ID]->singularity_index == singularID)
				AddToRegionTriangles(trajectories[traid_1+1][i].Triangle_ID, 0);
			}
		}


		//2. negative incoming direction
		endsingID = Object.flist[trajectories[traid_1+3][num_linesegs_curtraj[traid_1+3]-1].Triangle_ID]->singularity_index;

		if(endsingID >= 0 && IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
		{
			for(i = 0; i < num_linesegs_curtraj[traid_1+3]; i++)
			{
				////avoid to add the triangle containing singularity
				if(Object.flist[trajectories[traid_1+3][i].Triangle_ID]->contain_singularity != 1
					|| Object.flist[trajectories[traid_1+3][i].Triangle_ID]->singularity_index == singularID)
				AddToRegionTriangles(trajectories[traid_1+3][i].Triangle_ID, 0);
			}
		}

		if((endsingID >= 0 && !IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
			|| endsingID == -1) //endsingID == -1 means the separatrix ends at the boundary or closed orbit
		{
			for(i = 0; i < (int)((num_linesegs_curtraj[traid_1+3]-1)*length_percentage); i++)
			{
				////avoid to add the triangle containing singularity
				if(Object.flist[trajectories[traid_1+3][i].Triangle_ID]->contain_singularity != 1
					|| Object.flist[trajectories[traid_1+3][i].Triangle_ID]->singularity_index == singularID)
				AddToRegionTriangles(trajectories[traid_1+3][i].Triangle_ID, 0);
			}
		}
	}
}


/*-------------------------------------------------------------------*/
////11/20/05
////Get region for multiple repellers and attractors cancellation
////This routine only considers the singularity cancellation now!!!
////Modified at 11/30/05 (Adpative length of separatrix)
////Modified at 2/16/06 (Use percentage instead of the real length
/*-------------------------------------------------------------------*/
void InitSaddleGrowforMultRegion_new(int singularID, int type, 
									 int initsaddlelength, double length_percentage,
								     int *sametypesings, int numsings)
{
	int i;
    int endsingID;
	int num_addedtriangles = 0;
	int traid_1 = separatrices[singularities[singularID].separtices].sep1;

	//int dividen = pow(2, curtimes);  ////12.20

	if(type == 1)    ////Saddle as an attractor, forward tracing along outgoing direction
	{
		//going = singularities[singularID].outgoing;
		
		////tracing along both directions
		//1. positive outgoing direction

		////Get the singularity at the end of the separatrix
		endsingID = Object.flist[trajectories[traid_1][num_linesegs_curtraj[traid_1]-1].Triangle_ID]->singularity_index;
		if(endsingID >= 0 && IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
		{
			for(i = 0; i < num_linesegs_curtraj[traid_1]; i++)
			{
				AddToRegionforMultRegion(trajectories[traid_1][i].Triangle_ID, 1, num_addedtriangles);

				if(num_addedtriangles == initsaddlelength)
					break;
			}
		}

		if((endsingID >= 0 && !IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
			|| endsingID == -1) //endsingID == -1 means the separtrix ends at the boundary or closed orbit
		{
			//for(i = 0; i < (int)(num_linesegs_curtraj[traid_1]-1)/dividen; i++)
			for(i = 0; i < (int)((num_linesegs_curtraj[traid_1]-1)*length_percentage); i++)
			{
				AddToRegionforMultRegion(trajectories[traid_1][i].Triangle_ID, 1, num_addedtriangles);

				//if(num_addedtriangles == min(5,initsaddlelength))
				//if(num_addedtriangles == cursaddlelength)
				//	break;
			}
		}

		//2. negative outgoing direction
		endsingID = Object.flist[trajectories[traid_1+2][num_linesegs_curtraj[traid_1+2]-1].Triangle_ID]->singularity_index;
        num_addedtriangles = 0;

		if(endsingID >= 0 && IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
		{
			for(i = 0; i < num_linesegs_curtraj[traid_1+2]; i++)
			{
				AddToRegionforMultRegion(trajectories[traid_1+2][i].Triangle_ID, 1, num_addedtriangles);

				if(num_addedtriangles == initsaddlelength)
					break;
			}
		}

		if((endsingID >= 0 && !IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
			|| endsingID == -1) //endsingID == -1 means the separtrix ends at the boundary or closed orbit
		{
			for(i = 0; i < (int)((num_linesegs_curtraj[traid_1+2]-1)*length_percentage); i++)
			{
				AddToRegionforMultRegion(trajectories[traid_1+2][i].Triangle_ID, 1, num_addedtriangles);

			}
		}
	}

	/*---------------------------------------------------------------*/
	else{           
		////Saddle as a repeller, backward tracing along incoming direction
		
		////tracing along both directions
		//1. positive incoming direction
		endsingID = Object.flist[trajectories[traid_1+1][num_linesegs_curtraj[traid_1+1]-1].Triangle_ID]->singularity_index;
		if(endsingID >= 0 && IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
		{
			for(i = 0; i < num_linesegs_curtraj[traid_1+1]; i++)
			{
				AddToRegionforMultRegion(trajectories[traid_1+1][i].Triangle_ID, 0, num_addedtriangles);

				if(num_addedtriangles == initsaddlelength)
					break;
			}
		}

		if((endsingID >= 0 && !IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
			|| endsingID == -1) //endsingID == -1 means the separtrix ends at the boundary or closed orbit
		{
			for(i = 0; i < (int)((num_linesegs_curtraj[traid_1+1]-1)*length_percentage); i++)
			{
				AddToRegionforMultRegion(trajectories[traid_1+1][i].Triangle_ID, 0, num_addedtriangles);

				//if(num_addedtriangles == min(5,initsaddlelength))
				//if(num_addedtriangles == cursaddlelength)
				//	break;
			}
		}


		//2. negative incoming direction
		endsingID = Object.flist[trajectories[traid_1+3][num_linesegs_curtraj[traid_1+3]-1].Triangle_ID]->singularity_index;
		num_addedtriangles = 0;

		if(endsingID >= 0 && IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
		{
			for(i = 0; i < num_linesegs_curtraj[traid_1+3]; i++)
			{
				AddToRegionforMultRegion(trajectories[traid_1+3][i].Triangle_ID, 0, num_addedtriangles);

				if(num_addedtriangles == initsaddlelength)
					break;
			}
		}

		if((endsingID >= 0 && !IsRepeated(sametypesings, singularities[endsingID].node_index, numsings))
			|| endsingID == -1) //endsingID == -1 means the separtrix ends at the boundary or closed orbit
		{
			for(i = 0; i < (int)((num_linesegs_curtraj[traid_1+3]-1)*length_percentage); i++)
			{
				AddToRegionforMultRegion(trajectories[traid_1+3][i].Triangle_ID, 0, num_addedtriangles);

				//if(num_addedtriangles == min(5,initsaddlelength))
				//if(num_addedtriangles == cursaddlelength)
				//	break;
			}
		}
	}
}


/*
Count the number of triangles covering the specific separatrix
*/
int CountNumTriangleonSeparatrix(int index)
{
	int Num_triangles = 0;
	int *temp_triangles = new int[(int)sqrt((double)Object.nfaces)];

	int i;
	for(i = 0; i < num_linesegs_curtraj[index]; i++)
	{
		if(IsRepeated(temp_triangles, trajectories[index][i].Triangle_ID, Num_triangles))
			continue;

		//add to the temparary list of the triangles
		temp_triangles[Num_triangles] = trajectories[index][i].Triangle_ID;
		Num_triangles ++;
	}

	delete temp_triangles;
	return Num_triangles;
}


void GetMultRegion(int *repellers, int num_repellers, int *attractors, int num_attractors)
{
	int i;

	int singID;

	////Initialization part
	Source_re.num = 0;
	Source_re1.num = 0;
	Sink_re.num = 0;
	Sink_re1.num = 0;

	////1. Grow and calculate the union region of all repellers
	for(i = 0; i < num_repellers; i++)
	{
		////Reset
		repellerRegion.num = 0;
		Source_re1.num = 0;

		singID = graphnodes[repellers[i]].singularityID;

		////Grow region for each repeller
		Cancel_GrowRepellerRegion_adv(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength);

		////Union with a common variable storing the previous union result
		CopyRegion(Source_re.trianglelist, Source_re1.trianglelist, Source_re.num);
		Source_re1.num = Source_re.num;

		UnionRegion(repellerRegion, Source_re1, Source_re);
	}


	////2. Grow and calcualte the union region for all attractors
	for(i = 0; i < num_attractors; i++)
	{
		////Reset
		attractorRegion.num = 0;
		Sink_re1.num = 0;

		singID = graphnodes[attractors[i]].singularityID;

		////Grow region for each repeller
		Cancel_GrowAttractorRegion_adv(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength);

		////Union with a common variable storing the previous union result
		CopyRegion(Sink_re.trianglelist, Sink_re1.trianglelist, Sink_re.num);
		Sink_re1.num = Sink_re.num;

		UnionRegion(attractorRegion, Sink_re1, Sink_re);
	}

	////3. Get the intersected region of the previous two regions
	IntersectRegion(Source_re, Sink_re, intersectRegion);

	GetMultIntersection(repellers, num_repellers, attractors, num_attractors);
}


/*
Add at 11/22/05
*/
void AddToRegionforMultRegion(int triangleID, int type, int &num)
{
	if(triangleID < 0)  ////if the triangle does not exist
		return;

	int cur_trianglenum;
	if(type == 0)
		cur_trianglenum = repellerRegion.num;
	else
		cur_trianglenum = attractorRegion.num;

	////if the space for the triangles list is not enough, extend it
	if(cur_trianglenum >= MaxNumTriangle - 1)
	{
		MaxNumTriangle += 100;
		repellerRegion.trianglelist = (int*)realloc(repellerRegion.trianglelist, sizeof(int) * MaxNumTriangle);
		attractorRegion.trianglelist = (int*)realloc(attractorRegion.trianglelist, sizeof(int) * MaxNumTriangle);
		intersectRegion.trianglelist = (int*)realloc(intersectRegion.trianglelist, sizeof(int) * MaxNumTriangle);
	}

	////Make sure that the triangle not inside current triangles' list, very important!!!!

	if(InsideTrianglesListOrnot(triangleID, type))
		return;

	////set the flag of the triangle to mark it as a triangle inside the region
	if(type == 0)
	{
		Object.flist[triangleID]->repell_inregion = 1;
		repellerRegion.trianglelist[repellerRegion.num] = triangleID;
		repellerRegion.num++;
		num++;
	}
	else{
		Object.flist[triangleID]->attract_inregion = 1;
		attractorRegion.trianglelist[attractorRegion.num] = triangleID;
		attractorRegion.num++;
		num++;
	}
}

/*-------------------------------------------------------------------
////11/22/05
Get region for multiple repellers and attractors cancellation
This routine only considers the singularity cancellation now!!!
We tend to use new method to grow the two regions from both sources, sinks and saddles
-------------------------------------------------------------------*/

void GetMultRegion2(int *repellers, int num_repellers, int *attractors, int num_attractors)
{
	int i;

	int singID;

	////Initialization part
	Source_re.num = 0;
	Source_re1.num = 0;
	Sink_re.num = 0;
	Sink_re1.num = 0;

	////1. Grow and calculate the union region of all repellers
	for(i = 0; i < num_repellers; i++)
	{
		////Reset
		repellerRegion.num = 0;
		Source_re1.num = 0;

		singID = graphnodes[repellers[i]].singularityID;

		////Grow region for each repeller
		Cancel_GrowRepellerRegion_adv(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength);

		////Union with a common variable storing the previous union result
		CopyRegion(Source_re.trianglelist, Source_re1.trianglelist, Source_re.num);
		Source_re1.num = Source_re.num;

		UnionRegion(repellerRegion, Source_re1, Source_re);
	}

	////Grow from those intermediary saddles, here saddles act as repellers
	for(i = 0; i < Num_MediaNodes; i++)
	{
		repellerRegion.num = 0;
		Source_re1.num = 0;
		
		singID = graphnodes[MediaNodes[i]].singularityID;

		////Grow region for each repeller
		Cancel_GrowRepellerRegion_adv(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength);

		////Union with a common variable storing the previous union result
		CopyRegion(Source_re.trianglelist, Source_re1.trianglelist, Source_re.num);
		Source_re1.num = Source_re.num;

		UnionRegion(repellerRegion, Source_re1, Source_re);
	}

	////2. Grow and calcualte the union region for all attractors
	for(i = 0; i < num_attractors; i++)
	{
		////Reset
		attractorRegion.num = 0;
		Sink_re1.num = 0;

		singID = graphnodes[attractors[i]].singularityID;

		////Grow region for each repeller
		Cancel_GrowAttractorRegion_adv(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength);

		////Union with a common variable storing the previous union result
		CopyRegion(Sink_re.trianglelist, Sink_re1.trianglelist, Sink_re.num);
		Sink_re1.num = Sink_re.num;

		UnionRegion(attractorRegion, Sink_re1, Sink_re);
	}
	
	////Grow from those intermediary saddles, here saddles act as repellers
	for(i = 0; i < Num_MediaNodes; i++)
	{
		attractorRegion.num = 0;
		Sink_re1.num = 0;
		
		singID = graphnodes[MediaNodes[i]].singularityID;

		////Grow region for each repeller
		Cancel_GrowAttractorRegion_adv(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength);

		////Union with a common variable storing the previous union result
		CopyRegion(Sink_re.trianglelist, Sink_re1.trianglelist, Sink_re.num);
		Sink_re1.num = Sink_re.num;

		UnionRegion(attractorRegion, Sink_re1, Sink_re);
	}

	////3. Get the intersected region of the previous two regions
	IntersectRegion(Source_re, Sink_re, intersectRegion);

	GetMultIntersection(repellers, num_repellers, attractors, num_attractors);
}

/*-------------------------------------------------------------------
////11/22/05
Get region for multiple repellers and attractors cancellation
This routine only considers the singularity cancellation now!!!
We tend to use new method to grow the two regions from both sources, sinks and saddles
-------------------------------------------------------------------*/

//void AdaptiveGetMultRegion(int *repellers, int num_repellers, int *attractors, int num_attractors)
//{
//	int i;
//
//	int singID;
//	
//	//variables for adpative region growing
//	int small_length = 0;
//	int large_length = 2*InitSaddleRegionLength;
//	int cur_length;
//
//	////Initialization part
//
//	//variables for adaptive saddle regions
//    TriangularRegion temp_rsaddle_re, temp_asaddle_re;
//	temp_rsaddle_re.trianglelist = (int *)malloc(sizeof(int) * Object.nfaces);
//	temp_asaddle_re.trianglelist = (int *)malloc(sizeof(int) * Object.nfaces);
//
//	icVector2 *temp_vec = (icVector2 *)malloc(sizeof(icVector2) * Object.nverts);
//
//	temp_asaddle_re.num = 0;
//	temp_rsaddle_re.num = 0;
//
//	////Note that here, Source_re1 and Sink_re1 can be local variables! 11/30/05
//	Source_re.num = 0;
//	Source_re1.num = 0;
//	Sink_re.num = 0;
//	Sink_re1.num = 0;
//
//	////1. Save current whole field into a temporary array
//	for(i = 0; i < Object.nverts; i++)
//	{
//		VerticalField[i] = Object.vlist[i]->vec;
//	}
//
//	////Here Source_re and Sink_re will not change during the adpative steps
//	////2. Grow and calculate the union region of all repellers (sources only)
//	for(i = 0; i < num_repellers; i++)
//	{
//		////Reset
//		repellerRegion.num = 0;
//		Source_re1.num = 0;
//
//		singID = graphnodes[repellers[i]].singularityID;
//
//		////Grow region for each repeller
//		if(singularities[singID].type == SADDLE) //To avoid saddle growing break down 1/22/06
//			Cancel_GrowRepellerRegion_adv(singularities[singID].Triangle_ID, -1, singID, 1);
//		else
//		    Cancel_GrowRepellerRegion_adv(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength);
//
//		////Union with a common variable storing the previous union result
//		CopyRegion(Source_re.trianglelist, Source_re1.trianglelist, Source_re.num);
//		Source_re1.num = Source_re.num;
//
//		UnionRegion(repellerRegion, Source_re1, Source_re);
//	}
//
//	////Grow and calcualte the union region for all attractors(sinks only)
//	for(i = 0; i < num_attractors; i++)
//	{
//		////Reset
//		attractorRegion.num = 0;
//		Sink_re1.num = 0;
//
//		singID = graphnodes[attractors[i]].singularityID;
//
//		////Grow region for each attractor
//		if(singularities[singID].type == SADDLE) //To avoid saddle growing break down 1/22/06
//			Cancel_GrowAttractorRegion_adv(singularities[singID].Triangle_ID, -1, singID, 1);
//		else
//			Cancel_GrowAttractorRegion_adv(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength);
//
//		////Union with a common variable storing the previous union result
//		CopyRegion(Sink_re.trianglelist, Sink_re1.trianglelist, Sink_re.num);
//		Sink_re1.num = Sink_re.num;
//
//		UnionRegion(attractorRegion, Sink_re1, Sink_re);
//	}
//
//	/*--------------------------------------------------------------------------*/
//	//3. the followings are adaptive steps
//	int count = 0;
//    cur_length = InitSaddleRegionLength;
//	do{
//		//store the previous smoothing field for latter restore
//		for(i = 0; i < Object.nverts; i++)
//			temp_vec[i] = Object.vlist[i]->vec;
//
//		////1)restore the original field
//		for(i = 0; i < Object.nverts; i++)
//			Object.vlist[i]->vec = VerticalField[i];
//
//		////2)Grow from those intermediary saddles, here saddles act as repellers
//		temp_rsaddle_re.num = 0;
//
//		for(i = 0; i < Num_MediaNodes; i++)
//		{
//			repellerRegion.num = 0;
//			Source_re1.num = 0;
//			
//			singID = graphnodes[MediaNodes[i]].singularityID;
//
//			////Grow region for each repeller
//			//Cancel_GrowRepellerRegion_adv(singularities[singID].Triangle_ID, -1, singID, cur_length);
//			
//			Cancel_GrowRepellerRegion_adv(singularities[singID].Triangle_ID, -1, singID, count);
//
//			////Union with a common variable storing the previous union result
//			CopyRegion(temp_rsaddle_re.trianglelist, Source_re1.trianglelist, temp_rsaddle_re.num);
//			Source_re1.num = temp_rsaddle_re.num;
//
//			UnionRegion(repellerRegion, Source_re1, temp_rsaddle_re);
//		}
//			
//		//union with the previous sources' region
//		CopyRegion(temp_rsaddle_re.trianglelist, Source_re1.trianglelist, temp_rsaddle_re.num);
//		Source_re1.num = temp_rsaddle_re.num;
//
//		UnionRegion(Source_re, Source_re1, temp_rsaddle_re);
//
//		
//		////3)Grow from those intermediary saddles, here saddles act as repellers
//		temp_asaddle_re.num = 0;
//		for(i = 0; i < Num_MediaNodes; i++)
//		{
//			attractorRegion.num = 0;
//			Sink_re1.num = 0;
//			
//			singID = graphnodes[MediaNodes[i]].singularityID;
//
//			////Grow region for each repeller
//			//Cancel_GrowAttractorRegion_adv(singularities[singID].Triangle_ID, -1, singID, cur_length);
//			
//			Cancel_GrowAttractorRegion_adv(singularities[singID].Triangle_ID, -1, singID, count);
//
//			////Union with a common variable storing the previous union result
//			CopyRegion(temp_asaddle_re.trianglelist, Sink_re1.trianglelist, temp_asaddle_re.num);
//			Sink_re1.num = temp_asaddle_re.num;
//
//			UnionRegion(attractorRegion, Sink_re1, temp_asaddle_re);
//		}
//
//		//union with previous sinks' region
//		CopyRegion(temp_asaddle_re.trianglelist, Sink_re1.trianglelist, temp_asaddle_re.num);
//		Sink_re1.num = temp_asaddle_re.num;
//
//		UnionRegion(Sink_re, Sink_re1, temp_asaddle_re);
//
//		////4)Get the intersected region of the previous two regions
//		IntersectRegion(temp_rsaddle_re, temp_asaddle_re, intersectRegion);
//
//		GetMultIntersection(repellers, num_repellers, attractors, num_attractors);
//
//		////5) validate the euler character and manifold here (undone) 11/30/05
//
//		////6) perform one smoothing
//		Cancel_RegionSmooth();
//
//		////7) Get rid of those singularities?
//		if(IntersectedRegionSingCount() != abs(num_repellers + num_attractors - Num_MediaNodes))
//		{
//			//Undo();
//			break;
//		}
//        
//		large_length = cur_length;
//		cur_length = (int)(large_length + small_length)/2;
//
//		count++;
//
//	}while(large_length - small_length > 1);
//
//	////4. perform one more smoothing here if necessary
//	////if there is singularities after smoothing of the first region, stop! 
//	if(count == 0){  
//		MessageBox(NULL, "Can not find proper region!", "Error", MB_OK);
//	    Undo();
//		return;
//	}
//
//	//store the previous smoothing field for previous smoothing
//	for(i = 0; i < Object.nverts; i++)
//		Object.vlist[i]->vec = temp_vec[i];
//
//	free(temp_vec);
//	free(temp_rsaddle_re.trianglelist);
//	free(temp_asaddle_re.trianglelist);
//}
//
//
//
//
////Use new idea to perform adaptive region growing
int GetMinSepForSinkSourceCancel(int source, int sink, int *Media, int num_media)
{
	int min = Object.nfaces;

	int i, j;
	int sep;

	int traj;

	for(i = 0; i < num_media; i++)
	{
		//
		sep = singularities[Media[i]].separtices;

		for(j = 0; j < 4; j++)
		{
			switch(j){
				case 0:
					traj = separatrices[sep].sep1;
					break;
				case 1:
					traj = separatrices[sep].sep2;
					break;
				case 2:
					traj = separatrices[sep].sep3;
					break;
				case 3:
					traj = separatrices[sep].sep4;
					break;
			}

			if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID == singularities[source].Triangle_ID
				|| trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID == singularities[sink].Triangle_ID)
				continue;

			if(num_linesegs_curtraj[traj] < min) min = num_linesegs_curtraj[traj];
		}
	}

	return min;
}

////New idea: if the region grow fail in larger separatrix length, make it shorten
////if it succeed, make it longer until it fail, then shorten it, until 
////the length finally converges 2/16/06

////How about we cut the following routine into several smaller pieces?

/*
Input: repellers and attractors are all node list in the Conley partial order
*/

bool AdaptiveGetMultRegion_new(int *repellers, int num_repellers, int *attractors, int num_attractors)
{
	int i;
	int singID;
	int min_triangles = 0;
	int count = 0;
	int success = 0;
	
	//variables for adpative region growing, modified at 2/16/06
	double small_length = 0;
	double large_length = 1;
	double cur_length;

	////Initialization part

	//variables for adaptive saddle regions
    TriangularRegion temp_rsaddle_re, temp_asaddle_re;
	temp_rsaddle_re.trianglelist = (int *)malloc(sizeof(int) * Object.nfaces);
	temp_asaddle_re.trianglelist = (int *)malloc(sizeof(int) * Object.nfaces);

	temp_asaddle_re.num = 0;
	temp_rsaddle_re.num = 0;

	////Note that here, Source_re1 and Sink_re1 can be local variables! 11/30/05
	Source_re.num = 0;
	Source_re1.num = 0;
	Sink_re.num = 0;
	Sink_re1.num = 0;

	//we need to get the shortest separatrix
	//We assume only sink and source pair being considered
	min_triangles = GetMinSepForSinkSourceCancel(graphnodes[repellers[0]].singularityID,
		graphnodes[attractors[0]].singularityID, MediaNodes, Num_MediaNodes);

    
	cur_length = 1.;

	/*   Adaptive step   */

	do{
		InitCancellationAndMovement();
		Source_re.num = 0;
		Source_re1.num = 0;
		Sink_re.num = 0;
		Sink_re1.num = 0;

		//// Grow and calculate the union region of all repellers and attractors respectively
		for(i = 0; i < num_repellers; i++)
		{
			////Reset
			repellerRegion.num = 0;
			Source_re1.num = 0;

			singID = graphnodes[repellers[i]].singularityID;

			//we still assume one sink and one source pair now 07/02/06
			int target_t = singularities[graphnodes[attractors[0]].singularityID].Triangle_ID;

			Cancel_GrowRepellerRegion(singularities[singID].Triangle_ID, target_t, singID, InitSaddleRegionLength,
				0, 1.);

			////Union with a common variable storing the previous union result
			CopyRegion(Source_re.trianglelist, Source_re1.trianglelist, Source_re.num);
			Source_re1.num = Source_re.num;

			UnionRegion(repellerRegion, Source_re1, Source_re);
		}

		for(i = 0; i < num_attractors; i++)
		{
			////Reset
			attractorRegion.num = 0;
			Sink_re1.num = 0;

			singID = graphnodes[attractors[i]].singularityID;

			//we still assume one sink and one source pair now 07/02/06
			int target_t = singularities[graphnodes[attractors[0]].singularityID].Triangle_ID;

			Cancel_GrowAttractorRegion(singularities[singID].Triangle_ID, target_t, singID, InitSaddleRegionLength,
				0, 1.);

			////Union with a common variable storing the previous union result
			CopyRegion(Sink_re.trianglelist, Sink_re1.trianglelist, Sink_re.num);
			Sink_re1.num = Sink_re.num;

			UnionRegion(attractorRegion, Sink_re1, Sink_re);
		}


		/*   Adpative step   */
		//// Grow from those intermediary saddles, here saddles act as repellers
		temp_rsaddle_re.num = 0;

		for(i = 0; i < Num_MediaNodes; i++)
		{
			repellerRegion.num = 0;
			Source_re1.num = 0;
			
			singID = graphnodes[MediaNodes[i]].singularityID;

			Cancel_GrowRepellerRegion(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength,
				0, cur_length);

			////Union with a common variable storing the previous union result
			CopyRegion(temp_rsaddle_re.trianglelist, Source_re1.trianglelist, temp_rsaddle_re.num);
			Source_re1.num = temp_rsaddle_re.num;

			UnionRegion(repellerRegion, Source_re1, temp_rsaddle_re);
		}
			
		//union with the previous sources' region
		CopyRegion(temp_rsaddle_re.trianglelist, Source_re1.trianglelist, temp_rsaddle_re.num);
		Source_re1.num = temp_rsaddle_re.num;

		UnionRegion(Source_re, Source_re1, temp_rsaddle_re);

		
		////3)Grow from those intermediary saddles, here saddles act as repellers
		temp_asaddle_re.num = 0;
		for(i = 0; i < Num_MediaNodes; i++)
		{
			attractorRegion.num = 0;
			Sink_re1.num = 0;
			
			singID = graphnodes[MediaNodes[i]].singularityID;

			////Grow region for each repeller
		
			Cancel_GrowAttractorRegion(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength,
				0, cur_length);

			////Union with a common variable storing the previous union result
			CopyRegion(temp_asaddle_re.trianglelist, Sink_re1.trianglelist, temp_asaddle_re.num);
			Sink_re1.num = temp_asaddle_re.num;

			UnionRegion(attractorRegion, Sink_re1, temp_asaddle_re);
		}

		//union with previous sinks' region
		CopyRegion(temp_asaddle_re.trianglelist, Sink_re1.trianglelist, temp_asaddle_re.num);
		Sink_re1.num = temp_asaddle_re.num;

		UnionRegion(Sink_re, Sink_re1, temp_asaddle_re);

		//// Get the intersected region of the previous two regions
		IntersectRegion(temp_rsaddle_re, temp_asaddle_re, intersectRegion);

		GetMultIntersection(repellers, num_repellers, attractors, num_attractors);

		//// validate the euler character and manifold here (Conley boundary condition validate here)
		if(CalEulerValue(intersectRegion.trianglelist, intersectRegion.num) == 1)
		{
			success = 1;
			if(count == 0)
				break;                  //it is the possible largest region, stop here

			small_length = cur_length;  //try a larger regoin
		}

		else
			large_length = cur_length;  //try a smaller region

        
		cur_length = (large_length + small_length)/2;

		count++;

	}while(large_length - small_length > 1e-5 && (double)min_triangles*(large_length - small_length) > 0.5
		&& (int)min_triangles*cur_length >= 1 && count < 3);

	if(success == 1)
	{
		Cancel_RegionSmooth();

		if(IntersectedRegionSingCapture())
		{
			Undo();
			free(temp_rsaddle_re.trianglelist);
			free(temp_asaddle_re.trianglelist);
			return false;
		}

		free(temp_rsaddle_re.trianglelist);
		free(temp_asaddle_re.trianglelist);
		return true;
	}

	else
	{
		free(temp_rsaddle_re.trianglelist);
		free(temp_asaddle_re.trianglelist);
		return false;
	}

}




////Get the intersetion of multiple singularities cancellation
void GetMultIntersection(int *repellers, int num_repellers, int *attractors, int num_attractors)
{
	int i, j;
	Face *face;
	Edge *cur_edge;
	Vertex *vert;

	////Initial part
	for(i = 0; i < Object.nverts; i++)
	{
		vert = Object.vlist[i];

		vert->InRegion = 0;
		vert->OnBoundary = 0;
		vert->attract_flag = 0;
		vert->repell_flag = 0;
		vert->RegionListID = -1;
	}

	
	////Do we need to add the vertices of the triangles that contain the singularities we want to cancel
	////Yes, we need to do this
	////we do not test the repeatation here!!! we still consider singularity only
	for(i = 0; i < num_repellers; i++)
	{
		if(!IsRepeated(intersectRegion.trianglelist, 
			singularities[graphnodes[repellers[i]].singularityID].Triangle_ID,intersectRegion.num))
		{
			intersectRegion.trianglelist[intersectRegion.num]
				= singularities[graphnodes[repellers[i]].singularityID].Triangle_ID;
			intersectRegion.num++;
		}
	}
	
	for(i = 0; i < num_attractors; i++)
	{
		if(!IsRepeated(intersectRegion.trianglelist, 
			singularities[graphnodes[attractors[i]].singularityID].Triangle_ID,intersectRegion.num))
		{
			intersectRegion.trianglelist[intersectRegion.num]
				= singularities[graphnodes[attractors[i]].singularityID].Triangle_ID;
			intersectRegion.num++;
		}
	}

	////Set the flag for region vertices (not only inner vertices)
	for(i = 0; i < intersectRegion.num; i++)
	{
		face = Object.flist[intersectRegion.trianglelist[i]];

		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];

			vert->attract_flag = 2;
			vert->repell_flag = 2;
		}
	}

	////Reuse the UpdateBoundary routine to build the boundary edges' list
	////First, we need to reset the flag of edge visting
	for(i = 0; i < intersectRegion.num; i++)
	{
		face = Object.flist[intersectRegion.trianglelist[i]];
		for(j = 0; j < face->nverts; j++)
		{
			cur_edge =  face->edges[j];
			cur_edge->repell_visited = 0;
			cur_edge->attract_visited = 0;
		}
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
	for(i = 0; i < Object.nverts; i++)
	{
		vert = Object.vlist[i];
		if(vert->attract_flag == 2 && vert->repell_flag == 2 && vert->OnBoundary == 0) ////The vertex is inside both region
		{
			intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
			vert->InRegion = 1;
			vert->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
			intersectInnerverts.num ++;
		}
	}

	////Set the vertices of the triangle containing the being wanted to be cancelled singularities
	////as inner vertices 2/16/06

	for(i = 0; i < num_repellers; i++)
	{
		face = Object.flist[singularities[graphnodes[repellers[i]].singularityID].Triangle_ID];

		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];

			if(vert->InRegion == 1 || vert->RegionListID >= 0)
			{
				vert->OnBoundary = 0;
				continue;
			}
			
			intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
			vert->OnBoundary = 0;
			vert->InRegion = 1;
			vert->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
			intersectInnerverts.num ++;
		}
	}

	for(i = 0; i < num_attractors; i++)
	{
		face = Object.flist[singularities[graphnodes[attractors[i]].singularityID].Triangle_ID];

		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];

			if(vert->InRegion == 1 || vert->RegionListID >= 0)
			{
				vert->OnBoundary = 0;
				continue;
				}
			
			intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
			vert->OnBoundary = 0;
			vert->InRegion = 1;
			vert->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
			intersectInnerverts.num ++;
		}
	}

}


////Copy one region to a new array
//dest <- source
void CopyRegion(int *source, int *dest, int num)
{
	int i;

	for(i = 0; i < num; i++)
	{
	    dest[i] = source[i];	
	}
}

bool IstheMediaNode(int index)
{
	int i;

	for(i = 0; i < Num_MediaNodes; i++)
	{
		if(graphnodes[MediaNodes[i]].singularityID >= 0 
			&& graphnodes[MediaNodes[i]].singularityID == index)
			return true;
		else if(graphnodes[MediaNodes[i]].LimitCycleID >= 0
			&& graphnodes[MediaNodes[i]].LimitCycleID == index)
			return true;
	}

	return false;
}


/*------------------------------------------------------------------
Get the union of two regions
dest = source1 U source2
------------------------------------------------------------------*/
void UnionRegion(TriangularRegion &source1, TriangularRegion &source2, TriangularRegion &dest)
{
	int i, j;
	int cur_index;
	int cur_triangle;
	int repeat_flag = 0;

	for(i = 0; i < source1.num; i++)
	{
		dest.trianglelist[i] = source1.trianglelist[i];
	}

	cur_index = source1.num;

	for(i = 0; i < source2.num; i++)
	{
		cur_triangle = source2.trianglelist[i];
		repeat_flag = 0;

		for(j = 0; j < source1.num; j++)
		{
			if(cur_triangle == source1.trianglelist[j])
			{
				repeat_flag = 1;
				break;
			}
		}

		if(repeat_flag != 1)
		{
			dest.trianglelist[cur_index] = cur_triangle;
			cur_index++;
		}
	}

	dest.num = cur_index;
}


/*------------------------------------------------------------------
Get the intersections of two regions, quite slow when the mesh is huge!!!
dest = source1 ^ source2   11/20/05
------------------------------------------------------------------*/
void IntersectRegion(TriangularRegion &source1, TriangularRegion &source2, TriangularRegion &dest)
{
	int i, j;
	int cur_index = 0;
	int cur_triangle;


	for(i = 0; i < source1.num; i++)
	{
		cur_triangle = source1.trianglelist[i];

		for(j = 0; j < source2.num; j++)
		{
			if(cur_triangle == source2.trianglelist[j])
			{
				dest.trianglelist[cur_index] = cur_triangle;
				cur_index++;
			}
		}

	}

	dest.num = cur_index;
}





/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

/*****************************************************************************************/

////Routines for singularity movement

////Grow the virtual saddle region using higher threshold exit edges pending condition
bool Move_AttractorExitEdgePending(Edge *cur_edge)
{
	double dot1, dot2, dotresult = 0;
	Vertex *v1, *v2;
	v1 = Object.vlist[cur_edge->verts[0]];
	v2 = Object.vlist[cur_edge->verts[1]];

	dot1 = dot(cur_edge->attract_normal, v1->vec);
	dot2 = dot(cur_edge->attract_normal, v2->vec);

	dotresult = min(dot1, dot2);

	if(dotresult < 0/*2e-3*/) return true;   ////using higher threshold to get a proper region

	return false;

}

bool Move_RepellerExitEdgePending(Edge *cur_edge)
{
	double dot1, dot2, dotresult = 0;
	Vertex *v1, *v2;
	v1 = Object.vlist[cur_edge->verts[0]];
	v2 = Object.vlist[cur_edge->verts[1]];

	dot1 = dot(cur_edge->repell_normal, v1->vec);
	dot2 = dot(cur_edge->repell_normal, v2->vec);

	dotresult = max(dot1, dot2);

	if(dotresult > 0/*1e-4*/) return true;   ////using higher threshold to get a proper region

	return false;

}

////growing the source(repeller) region for movement, testing 07/14/05
void Move_Growing(int type, int target_triangle)
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
				exitornot = Move_RepellerExitEdgePending(cur_edge);
			}
			else{          ////attractor
				exitornot = Move_AttractorExitEdgePending(cur_edge);
			}

			oppositeTriangle = GetOppositeTriangle(cur_edge, type);

			if(Object.flist[oppositeTriangle]->contain_singularity == 1) ////if the opposite triangle containing singularity
			{
				if(oppositeTriangle != target_triangle) ////not the triangle containing the counterpart singularity
					continue;
			}

			if(exitornot && oppositeTriangle > 0)  ////it is an exiting edge
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

////Remove specific vertex from the Innerverts list
void Move_RemoveVertex(Vertex *del_vert, InnerVertices &innerverts)
{
	int i;
	////Find the vertex in the list
	for(i = 0; i < innerverts.num; i++)
	{
		if(innerverts.vertslist[i] == del_vert)
			break;
	}

	if(i >= innerverts.num)
		return;

	for( ; i < innerverts.num-1; i++)
	{
		innerverts.vertslist[i] = innerverts.vertslist[i+1];
		innerverts.vertslist[i]->RegionListID = i;                  ////Update the index of the vertex in the list
	}

	innerverts.num--;
}


//void Move_RemoveTriangle(int del_triangle, TriangularRegion &region)
//{
//}


void Move_GetIntersectedRegion(int sourceID, int target_triangle)
{
	int i, j;
	Face *face;
	Edge *cur_edge;
	Vertex *vert;

	////Get the intersect triangular region
	intersectRegion.num = 0;
	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];
		if(face->repell_inregion == 1 && face->attract_inregion == 1)
		{
			intersectRegion.trianglelist[intersectRegion.num] = face->index;
			intersectRegion.num ++;
		}
		else{
			face->repell_inregion = 0;
			face->attract_inregion = 0;
		}
	}

	////Reuse the UpdateBoundary routine to build the boundary edges' list
	////First, we need to reset the flag of edge visting
	for(i = 0; i < intersectRegion.num; i++)
	{
		face = Object.flist[intersectRegion.trianglelist[i]];
		for(j = 0; j < 3; j++)
		{
			cur_edge =  face->edges[j];
			cur_edge->repell_visited = 0;
			cur_edge->attract_visited = 0;
		}
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

	//// Set the second boundary of the region and  the vectors of the target triangle
	face = Object.flist[target_triangle];
	for(i = 0; i < 3; i++)
	{
		Object.vlist[face->verts[i]]->OnBoundary = 1;
		Object.vlist[face->verts[i]]->InRegion = 0;
	}

	////Get the intersection of inner vertices
	intersectInnerverts.num = 0;
	for(i = 0; i < Object.nverts; i++)
	{
		vert = Object.vlist[i];
		if(vert->attract_flag == 2 && vert->repell_flag == 2 && vert->OnBoundary == 0) ////The vertex is inside both region
		{
			intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
			vert->InRegion = 1;
			vert->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
			intersectInnerverts.num ++;
		}
	}

	////Add the vertices of the triangles that contain the source to the innerverts list
	////Testing codes here
	face = Object.flist[singularities[sourceID].Triangle_ID];
	for(j = 0; j < face->nverts; j++)
	{
		vert = Object.vlist[face->verts[j]];

		if(vert->InRegion != 1)
		{
			intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
			vert->InRegion = 1;
			vert->OnBoundary = 0;
			vert->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
			intersectInnerverts.num ++;
		}
	}
	
	///////////////////////////////////////////////////////////////////////////////////////////
}


////Get the smooth region for singularity movement
void Move_GetRegion(int triangleID, int singID, int target_triangle, double newx, double newy)
{
	////1. Grow the repeller region covering the original source
	Move_GrowRepellerRegion(triangleID, target_triangle);

	////2. Grow the virtual saddle region covering the new position
	Move_GrowAttractorRegion(triangleID, target_triangle, singID, newx, newy);

	////3. Get the final smoothing region
	//Move_BuildSmoothRegion();
	Move_GetIntersectedRegion(singID, target_triangle);


	////Set the three vectors on the vertices of the second boundary
	Move_SetVectorsAtNewPos(triangleID, target_triangle); 

}


void Move_SetVectorsAtNewPos(int source_triangle, int target_triangle)
{
	int i;
	icVector2 globalv, localv;
	Face *face1, *face2;
	face1 = Object.flist[source_triangle];
	face2 = Object.flist[target_triangle];


    ////New method to set the three vectors
	face2->direct_vec[0].entry[0]=
		face2->xy[0][0] - (face2->xy[1][0]+face2->xy[2][0])/2;
	face2->direct_vec[0].entry[1]=
		face2->xy[0][1] - (face2->xy[1][1]+face2->xy[2][1])/2;
	normalize(face2->direct_vec[0]);
	
	face2->direct_vec[1].entry[0]=
		face2->xy[1][0] - (face2->xy[0][0]+face2->xy[2][0])/2;
	face2->direct_vec[1].entry[1]=
		face2->xy[1][1] - (face2->xy[0][1]+face2->xy[2][1])/2;
	normalize(face2->direct_vec[1]);

	face2->direct_vec[2].entry[0]=
		face2->xy[2][0] - (face2->xy[1][0]+face2->xy[0][0])/2;
	face2->direct_vec[2].entry[1]=
		face2->xy[2][1] - (face2->xy[1][1]+face2->xy[0][1])/2;
	normalize(face2->direct_vec[2]);

	for(i = 0; i < 3; i++)
	{
		globalv = face2->direct_vec[i].entry[0]*face2->LX + face2->direct_vec[i].entry[1]*face2->LY;
		Object.vlist[face2->verts[i]]->vec = 0.12 * globalv;
	}
}


////Get the incoming and outgoing eigen vectors for the virtual saddle
void Move_GetVirtualSaddle(int old_singID, double x, double y, icVector2 &in, icVector2 &out)
{
	in.entry[0] = x - singularities[old_singID].gcx;
	in.entry[1] = y - singularities[old_singID].gcy;

	normalize(in);

	out.entry[0] = -in.entry[1];
	out.entry[1] = in.entry[0];
}


///////////////////////////////////////
void Move_InitSaddleRegion(int triangleID, double newx, double newy)
{
	int i;
	int flag = 0;
	double globalp[2];
	int pre_face, cur_face;

	////first, begin from the snew position to perform forward tracing to get s'
	////Here we do not use optimized method to find s', just trace from snew through at most 4 triangle
	pre_face = cur_face = triangleID;
	globalp[0] = newx;   
	globalp[1] = newy;
	for(i = 0; i < 10; i++)   ////How long we should go from snew to find s'??????!!!!!!
	{
		if(cur_face == -1)
			break;

		pre_face = cur_face;
		cur_face = TraceInATriangle2(cur_face, globalp, 0, flag); ////0 means always forward tracing here

		if(flag == 3 || flag == 4 || pre_face == cur_face) 
			break;
	}
	////Now, globalp storing the position of s'
	
	
	////second, rotate the whole field with 90 degree
	for(i = 0; i < Object.nverts; i++)
	{
		VerticalField[i].entry[0] = Object.vlist[i]->vec.entry[0];
		VerticalField[i].entry[1] = Object.vlist[i]->vec.entry[1];

		Object.vlist[i]->vec.entry[0] = -VerticalField[i].entry[1];
		Object.vlist[i]->vec.entry[1] =  VerticalField[i].entry[0];
	}
	GetLocalVector();

	////Perform initisaddle region growing as usual here under the vertical field
	double saved_globalp[2];
	saved_globalp[0] = globalp[0];
	saved_globalp[1] = globalp[1];
	int saved_face = cur_face;           ////save the triangle that contains s'
	for(i = 0; i < 30/*InitSaddleRegionLength*/; i++)
	{
		if(cur_face == -1)
			break;
		
		AddToRegionTriangles(cur_face, 1);

		pre_face = cur_face;
		cur_face = TraceInATriangle2(cur_face, globalp, 0, flag); ////0 means always forward tracing here

		if(flag == 3 || flag == 4 || pre_face == cur_face) 
			break;
	}
	
	globalp[0] = saved_globalp[0];
	globalp[1] = saved_globalp[1];
	cur_face = saved_face;
	for(i = 0; i < 30/*InitSaddleRegionLength*/; i++)
	{
		if(cur_face == -1)
			break;

		pre_face = cur_face;
		cur_face = TraceInATriangle2(cur_face, globalp, 1, flag); ////1 means always backward tracing here

		if(flag == 3 || flag == 4 || pre_face == cur_face) 
			break;

		AddToRegionTriangles(cur_face, 1);
	}

	////finally, recovery the original field for following processing
	for(i = 0; i < Object.nverts; i++)
	{
		Object.vlist[i]->vec.entry[0] = VerticalField[i].entry[0];
		Object.vlist[i]->vec.entry[1] = VerticalField[i].entry[1];
	}
	GetLocalVector();
}



////Get the initial triangular strip for the virtual saddle for region growing
////Grow along the outgoing direction (straight line), judge which triangles it passes
void Move_GetInitRegionofVirtualSaddle(int Init_triangle, icVector2 out, double newx, double newy)
{
	int i, c, pre_triangle, cur_triangle;
	double testx, testy;

	////Track along the positive outgoing direction
	cur_triangle = Init_triangle;
	testx = newx + SEPARATRIXSTEP*out.entry[0];
	testy = newy + SEPARATRIXSTEP*out.entry[1];
	for(i = 0; i < InitSaddleRegionLength; i++)
	{
		pre_triangle = cur_triangle;
		////Track inside a triangle until it go outside of the triangle or reaches the maximum steps
		for(c = 0; c < 400; c++)
		{
			testx += 0.008 * out.entry[0];
			testy += 0.008 * out.entry[1];
			cur_triangle = TriangleDetect(testx, testy);

			if(cur_triangle != pre_triangle)
				break;
		}

		if(cur_triangle < 0 || cur_triangle > Object.nfaces)
			break;

		////else, add to the attractor region
		AddToRegionTriangles(cur_triangle, 1);
	}
	
	////Track along the negative outgoing direction
	cur_triangle = Init_triangle;
	testx = newx - SEPARATRIXSTEP*out.entry[0];
	testy = newy - SEPARATRIXSTEP*out.entry[1];
	for(i = 0; i < InitSaddleRegionLength; i++)
	{
		pre_triangle = cur_triangle;
		////Track inside a triangle until it go outside of the triangle or reaches the maximum steps
		for(c = 0; c < 400; c++)
		{
			testx -= 0.008 * out.entry[0];
			testy -= 0.008 * out.entry[1];
			cur_triangle = TriangleDetect(testx, testy);

			if(cur_triangle != pre_triangle)
				break;
		}

		if(cur_triangle < 0 || cur_triangle > Object.nfaces)
			break;

		////else, add to the attractor region
		AddToRegionTriangles(cur_triangle, 1);
	}
}



////
void Move_GrowRepellerRegion(int triangleID, int target_triangle)
{
	////initialize repeller region with first input triangle
	repellerRegion.trianglelist[0] = triangleID;
	Object.flist[triangleID]->repell_inregion = 1;
	repellerRegion.num = 1;

	////initialize repeller boundary with the 3 edges of first triangle
	repellerBoundary.edgelist[0] = Object.flist[triangleID]->edges[0];
	repellerBoundary.edgelist[1] = Object.flist[triangleID]->edges[1];
	repellerBoundary.edgelist[2] = Object.flist[triangleID]->edges[2];
	repellerBoundary.num = 3;

	GetRegionNormals(0); ////Get the outward normal of the region for each boundary edge

	////initialize the inner vertices list
	repellerInnerverts.num = 0;    ////no inner vertex at present

	////
	Move_Growing(0, target_triangle);

}


////
void Move_GrowAttractorRegion(int triangleID, int target_triangle, int singularID, double newx, double newy)
{
	icVector2 in, out;
	
	////Calculate the incoming and outgoing eigen vector of the virtual saddle
	//Move_GetVirtualSaddle(singularID, newx, newy, in, out);

	//double s_primex = newx + 0.05*in.entry[0];
	//double s_primey = newy + 0.05*in.entry[1];

	//triangleID = TriangleDetect(s_primex, s_primey);
	//////Get the initial region for attractor growing from the s'
	//Move_GetInitRegionofVirtualSaddle(target_triangle, out, s_primex, s_primey);
    Move_InitSaddleRegion(target_triangle, newx, newy);

	UpdateBoundary(1);
	GetRegionNormals(1);

	////
	Move_Growing(1, target_triangle);
}
/*****************************************************************************************/

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

//////Initialize the beginning region of a saddle as a triangular strip
//////The first triangle that contain the saddle has been added to the list before calling this routine
//void Cancel_InitSaddleRegion(int singularID, int type, int initsaddlelength)
//{
//	icVector2 going;
//
//	int i;
//	int flag = 0;
//	double globalp[2];
//	int pre_face, cur_face;
//	int loop_flag = 0;  ////for detect a closed loop from the beginning triangle
//	int *triangles = (int*)malloc(sizeof(int) * NUMTRACINGTRIANGLE);
//	int num_triangles = 0;
//
//
//	if(type == 0)    ////Saddle as an attractor, forward tracing along outgoing direction
//	{
//		going = singularities[singularID].outgoing;
//		
//		////tracing along both directions
//		//1. positive outgoing direction
//		pre_face = cur_face = singularities[singularID].Triangle_ID;
//		globalp[0] = singularities[singularID].gcx + SEPARATRIXSTEP * going.entry[0];   
//		globalp[1] = singularities[singularID].gcy + SEPARATRIXSTEP * going.entry[1];
//		for(i = 0; i < initsaddlelength; i++)
//		{
//			if(cur_face < 0)
//			{
//				break;
//			}
//
//			pre_face = cur_face;
//			cur_face = TraceInATriangle2(cur_face, globalp, type, flag); ////0 means always forward tracing here
//
//			if(flag == 3 || flag == 4 || pre_face == cur_face || loop_flag == 1) 
//			{
//				break;
//			}
//
//			AddToRegionTriangles(cur_face, 1);
//			
//			if(LoopDetect(triangles, num_triangles, cur_face))
//				loop_flag = 1;
//			else{
//				triangles[num_triangles] = cur_face;
//				num_triangles++;
//			}
//		}
//
//		//2. negative outgoing direction
//		flag = 0;
//		loop_flag = 0;
//		pre_face = cur_face = singularities[singularID].Triangle_ID;
//		globalp[0] = singularities[singularID].gcx - SEPARATRIXSTEP * going.entry[0];   
//		globalp[1] = singularities[singularID].gcy - SEPARATRIXSTEP * going.entry[1];
//		for(i = 0; i < initsaddlelength; i++)
//		{
//			if(cur_face < 0)
//			{
//				break;
//			}
//
//			pre_face = cur_face;
//			cur_face = TraceInATriangle2(cur_face, globalp, type, flag); ////0 means always forward tracing here
//
//			if(flag == 3 || flag == 4 || pre_face == cur_face || loop_flag == 1) 
//			{
//				break;
//			}
//			AddToRegionTriangles(cur_face, 1);
//			//attractorRegion.trianglelist[attractorRegion.num] = cur_face;
//			//attractorRegion.num ++;
//			if(LoopDetect(triangles, num_triangles, cur_face))
//				loop_flag = 1;
//			else{
//				triangles[num_triangles] = cur_face;
//				num_triangles++;
//			}
//		}
//	}
//	else{            ////Saddle as a repeller, backward tracing along incoming direction
//		going = singularities[singularID].incoming;
//		
//		////tracing along both directions
//		//1. positive outgoing direction
//		pre_face = cur_face = singularities[singularID].Triangle_ID;
//		globalp[0] = singularities[singularID].gcx + SEPARATRIXSTEP * going.entry[0];   
//		globalp[1] = singularities[singularID].gcy + SEPARATRIXSTEP * going.entry[1];
//		for(i = 0; i < InitSaddleRegionLength; i++)
//		{
//			if(cur_face < 0)
//			{
//				break;
//			}
//
//			pre_face = cur_face;
//			cur_face = TraceInATriangle2(cur_face, globalp, 1, flag); ////1 means always backward tracing here
//
//			if(flag == 3 || flag == 4 || pre_face == cur_face || loop_flag == 1) 
//			{
//				break;
//			}
//
//
//			AddToRegionTriangles(cur_face, 0);
//			
//			if(LoopDetect(triangles, num_triangles, cur_face))
//				loop_flag = 1;
//			else{
//				triangles[num_triangles] = cur_face;
//				num_triangles++;
//			}
//		}
//
//		//2. negative outgoing direction
//		flag = 0;
//		loop_flag = 0;
//		pre_face = cur_face = singularities[singularID].Triangle_ID;
//		globalp[0] = singularities[singularID].gcx - SEPARATRIXSTEP * going.entry[0];   
//		globalp[1] = singularities[singularID].gcy - SEPARATRIXSTEP * going.entry[1];
//		for(i = 0; i < InitSaddleRegionLength; i++)
//		{
//			if(cur_face < 0)
//			{
//				break;
//			}
//
//			pre_face = cur_face;
//			cur_face = TraceInATriangle2(cur_face, globalp, 1, flag); ////0 means always forward tracing here
//
//			if(flag == 3 || flag == 4 || pre_face == cur_face || loop_flag == 1) 
//			{
//				break;
//			}
//
//			AddToRegionTriangles(cur_face, 0);
//			
//			if(LoopDetect(triangles, num_triangles, cur_face))
//				loop_flag = 1;
//			else{
//				triangles[num_triangles] = cur_face;
//				num_triangles++;
//			}
//		}
//	}
//}


/*------------------------------------------------------------------------------*/
////Initialize the beginning region of a saddle as a triangular strip
////The first triangle that contain the saddle has been added to the list before calling this routine
////Reuse the calculation result of separatrices
////sep1 and sep3 are outgoing direction; sep2 and sep4 are incoming direction

////Now we can not adjust the length !!  11/22/05
////We need to find some way to count the number of triangles that has been added to the region

void Cancel_InitSaddleRegion(int saddle, int type, double percentage)
{

	int i;

	int sep = singularities[saddle].separtices;

	int traj;

	int pre_triangle = -1;

	if(type == 0)    ////Saddle as an attractor, forward tracing along outgoing direction
	{

		////tracing along both directions
		//1. positive outgoing direction
		traj = separatrices[sep].sep1;
		for(i = 0; i < (int)num_linesegs_curtraj[traj]*percentage; i++)
		{
			if(trajectories[traj][i].Triangle_ID != singularities[saddle].Triangle_ID
				&& Object.flist[trajectories[traj][i].Triangle_ID]->contain_singularity == 1)
				break;

			if(trajectories[traj][i].Triangle_ID == pre_triangle)
				continue;

			AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 1);
			pre_triangle = trajectories[traj][i].Triangle_ID;
		}

		//2. negative outgoing direction
		traj = separatrices[sep].sep3;
		pre_triangle = -1;
		for(i = 0; i < (int)num_linesegs_curtraj[traj]*percentage; i++)
		{
			if(trajectories[traj][i].Triangle_ID != singularities[saddle].Triangle_ID
				&& Object.flist[trajectories[traj][i].Triangle_ID]->contain_singularity == 1)
				break;

			if(trajectories[traj][i].Triangle_ID == pre_triangle)
				continue;

			AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 1);
			pre_triangle = trajectories[traj][i].Triangle_ID;
		}
	}

	else{           
		////Saddle as a repeller, backward tracing along incoming direction
		
		////tracing along both directions
		//1. positive incoming direction
		traj = separatrices[sep].sep2;
		for(i = 0; i <  (int)num_linesegs_curtraj[traj]*percentage; i++)
		{
			if(trajectories[traj][i].Triangle_ID != singularities[saddle].Triangle_ID
				&& Object.flist[trajectories[traj][i].Triangle_ID]->contain_singularity == 1)
				break;

			if(trajectories[traj][i].Triangle_ID == pre_triangle)
				continue;

			AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 0);
			pre_triangle = trajectories[traj][i].Triangle_ID;
		}

		//2. negative incoming direction
		for(i = 0; i < (int)num_linesegs_curtraj[traj]*percentage; i++)
		{
			if(trajectories[traj][i].Triangle_ID != singularities[saddle].Triangle_ID
				&& Object.flist[trajectories[traj][i].Triangle_ID]->contain_singularity == 1)
				break;

			if(trajectories[traj][i].Triangle_ID == pre_triangle)
				continue;

			AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 0);
			pre_triangle = trajectories[traj][i].Triangle_ID;
		}
	}
}





/*-------------------------------------------------------------------------------------*/
////Routines for two connection orbit cases 2/23/06

//judge whether there are two connections orbits between the input singularities
//one of this singularity must be saddle
bool TwoConnectOrbits(int saddle, int singID)
{
	if(singularities[saddle].type != SADDLE)
		return false;

	int sep = singularities[saddle].separtices;
	int traj;
	int count_connection = 0;
	int i;

	for(i = 0; i < 4; i++)
	{
		switch(i)
		{
		case 0:
			traj = separatrices[sep].sep1;
			break;
		case 1:
			traj = separatrices[sep].sep2;
			break;
		case 2:
			traj = separatrices[sep].sep3;
			break;
		case 3:
			traj = separatrices[sep].sep4;
			break;
		}

		if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID == singularities[singID].Triangle_ID)
			count_connection ++;
	}

	if(count_connection == 2)
		return true;
	return false;
}


/*--------------------------------------------------------------*/
/*Get the direction of the triangle containing the saddle
Input: saddle index
       and the separatrix that user want to reverse
	   we always assume that user input the correct trajectory index
*/
icVector2 GetSaddleTriangleDirection(int saddle, int traj)
{
	int i, sep;

	sep = singularities[saddle].separtices;
	for(i = 0; i < 4; i++)
	{
		switch(i){
			case 0:
				if(traj == separatrices[sep].sep1)
				{
					return singularities[saddle].outgoing;
				}
			case 1:
				if(traj == separatrices[sep].sep2)
				{
					return singularities[saddle].incoming;
				}
			case 2:
				if(traj == separatrices[sep].sep3)
				{
					return (-singularities[saddle].outgoing);
				}
			case 3:
				if(traj == separatrices[sep].sep4)
				{
					return (-singularities[saddle].incoming);
				}
		}
	}
}


////Set the extra constrains for the two connection orbit cancelllation
void SetExtraBoundary(int saddle, int singID, icVector2 saddle_dir, int traj)
{
	Face *face;
	Vertex *vert;
	int i;

	////set the triangle containing the saddle
	face = Object.flist[singularities[saddle].Triangle_ID];;
	for(i = 0; i < face->nverts; i++)
	{
		vert = Object.vlist[face->verts[i]];
		vert->vec = -0.01*saddle_dir;
						
		if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
			vert->vec *= -1;

		vert->OnBoundary = 1;
		vert->InRegion = 0;
		vert->repell_flag = 0;
		vert->attract_flag = 0;
	}


	////set the triangle containing the other singularity
	face = Object.flist[singularities[singID].Triangle_ID];
	for(i = 0; i < face->nverts; i++)
	{
		vert = Object.vlist[face->verts[i]];
		vert->vec = 0.01*saddle_dir;
		
		if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
			vert->vec *= -1;
		
		vert->OnBoundary = 1;
		vert->InRegion = 0;
		vert->repell_flag = 0;
		vert->attract_flag = 0;
	}

	////do we need to reset the vectors on the triangles containing the selected separatrix
	////2/27/06
	int numsegs = num_linesegs_curtraj[traj];
	icVector2 line_dir, n_vec1, n_vec2;
	icVector2 eva_dis1, eva_dis2;
	double theta = (10./90.)*(M_PI / 2);
	for(i = 0; i < (int)numsegs/1.6; i++)
	{
		if(trajectories[traj][i].Triangle_ID == singularities[saddle].Triangle_ID
			|| trajectories[traj][i].Triangle_ID == singularities[singID].Triangle_ID)
			continue;

		face = Object.flist[trajectories[traj][i].Triangle_ID];

		if(face->attract_inregion == 1 || face->repell_inregion == 1)
		{
			face->attract_inregion = face->repell_inregion = 0;

			////Get the approximate direction of the separatrix in the triangle
			line_dir.entry[0] = trajectories[traj][i].gend[0] - trajectories[traj][i].gstart[0];
			line_dir.entry[1] = trajectories[traj][i].gend[1] - trajectories[traj][i].gstart[1];
			
			if(singularities[singID].type == SINK || singularities[singID].type == AFOCUS)
				line_dir *= -1;
			
			normalize(line_dir);
			
			n_vec1.entry[0] = line_dir.entry[0]*cos(theta) - line_dir.entry[1]*sin(theta);
			n_vec1.entry[1] = line_dir.entry[0]*sin(theta) + line_dir.entry[1]*cos(theta);
			
			n_vec2.entry[0] = line_dir.entry[0]*cos(theta) + line_dir.entry[1]*sin(theta);
			n_vec2.entry[1] = -line_dir.entry[0]*sin(theta) + line_dir.entry[1]*cos(theta);

			for(int j = 0; j < face->nverts; j++)
			{
				vert = Object.vlist[face->verts[j]];

				if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
				{
					////we may create a repelling limit cycle, it means that all the vectors 
					////on the boundary should separate away from the separatrix
					eva_dis1.entry[0] = vert->x+0.005*n_vec1.entry[0] - trajectories[traj][i].gend[0];
					eva_dis1.entry[1] = vert->y+0.005*n_vec1.entry[1] - trajectories[traj][i].gend[1];

					eva_dis2.entry[0] = vert->x+0.005*n_vec2.entry[0] - trajectories[traj][i].gend[0];
					eva_dis2.entry[1] = vert->y+0.005*n_vec2.entry[1] - trajectories[traj][i].gend[1];

					if(length(eva_dis1) < length(eva_dis2))
					{
						vert->vec = 0.03*n_vec2;
					}
					else
						vert->vec = 0.03*n_vec1;
				}
				else if(singularities[singID].type == SINK || singularities[singID].type == AFOCUS)
				{
					////we may create an attractor limit cycle, it means that all the vectors
					////on the boundary should converge to the separatrix
					
					eva_dis1.entry[0] = vert->x+0.005*n_vec1.entry[0] - trajectories[traj][i].gend[0];
					eva_dis1.entry[1] = vert->y+0.005*n_vec1.entry[1] - trajectories[traj][i].gend[1];

					eva_dis2.entry[0] = vert->x+0.005*n_vec2.entry[0] - trajectories[traj][i].gend[0];
					eva_dis2.entry[1] = vert->y+0.005*n_vec2.entry[1] - trajectories[traj][i].gend[1];

					if(length(eva_dis1) < length(eva_dis2))
					{
						vert->vec = 0.03*n_vec1;
					}
					else
						vert->vec = 0.03*n_vec2;
				}

				vert->OnBoundary = 1;
				vert->InRegion = 0;
				vert->repell_flag = 0;
				vert->attract_flag = 0;
			}
		}
	}
}



////New method to reset the second boundary 3/9/06
////We use adaptive method to set the rotation angle
void SetExtraBoundary2(int saddle, int singID, icVector2 saddle_dir, int traj)
{
	Face *face;
	Vertex *vert;
	int i;

	/*-----------------------------------------------------------------------*/
	////set the triangle containing the saddle
	face = Object.flist[singularities[saddle].Triangle_ID];;
	for(i = 0; i < face->nverts; i++)
	{
		vert = Object.vlist[face->verts[i]];
		vert->vec = -0.01*saddle_dir;
						
		if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
			vert->vec *= -1;

		vert->OnBoundary = 1;
		vert->InRegion = 0;
		vert->repell_flag = 0;
		vert->attract_flag = 0;
	}


	////set the triangle containing the other singularity
	face = Object.flist[singularities[singID].Triangle_ID];
	for(i = 0; i < face->nverts; i++)
	{
		vert = Object.vlist[face->verts[i]];
		vert->vec = 0.01*saddle_dir;
		
		if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
			vert->vec *= -1;
		
		vert->OnBoundary = 1;
		vert->InRegion = 0;
		vert->repell_flag = 0;
		vert->attract_flag = 0;
	}
	/*-----------------------------------------------------------------------*/

	////reset the vectors on the triangles containing the selected separatrix
	////3/9/06
	////1. We first find all the triangles
	icVector2 line_dir, n_vec1, n_vec2;
	icVector2 eva_dis1, eva_dis2;
	int *triangles = (int*)malloc(sizeof(int)*Object.nfaces/20);
	int num_triangles = 0;
	int numsegs = num_linesegs_curtraj[traj];
	int pre_triangle = -1;
	int j;
	double largest_dis = 0;

	////Build the triangle list
	for(i = 0; i < numsegs; i++)
	{
		if(pre_triangle != trajectories[traj][i].Triangle_ID)
		{
			triangles[num_triangles] = trajectories[traj][i].Triangle_ID;
			pre_triangle = trajectories[traj][i].Triangle_ID;
			num_triangles++;
		}
	}

	////Reset the distance for each vertex
	for(i = 0; i < num_triangles; i++)
	{
		face = Object.flist[triangles[i]];

		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];
			vert->distance = 1e40;
		    vert->visited = 0;
		}
	}

	////calculate the distances for all the vertices on the boundaries
	for(i = 0; i < num_triangles; i++)
	{
		face = Object.flist[triangles[i]];

		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];

			if(vert->distance > 1e20 && vert->visited == 0)
			{
				FindTheDisForOneVert(triangles, num_triangles, vert->VertID);
				vert->visited = 1;
			}
		}
	}

	////Find the maximum distance
	for(i = 0; i < num_triangles; i++)
	{
		face = Object.flist[triangles[i]];

		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];
			if(vert->distance > largest_dis)
				largest_dis = vert->distance;
			vert->visited = 0;
		}
	}

	////Calculate the vector values according to the distance
	double theta, magnitude = 0;
	for(i = 0; i < (int)numsegs/1.5; i++)
	{
		if(trajectories[traj][i].Triangle_ID == singularities[saddle].Triangle_ID
			|| trajectories[traj][i].Triangle_ID == singularities[singID].Triangle_ID)
			continue;

		face = Object.flist[trajectories[traj][i].Triangle_ID];

		if(face->attract_inregion == 1 || face->repell_inregion == 1)
		{
			face->attract_inregion = face->repell_inregion = 0;

			////Get the approximate direction of the separatrix in the triangle
			line_dir.entry[0] = trajectories[traj][i].gend[0] - trajectories[traj][i].gstart[0];
			line_dir.entry[1] = trajectories[traj][i].gend[1] - trajectories[traj][i].gstart[1];
			
			if(singularities[singID].type == SINK || singularities[singID].type == AFOCUS)
			//if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
				line_dir *= -1;
			
			normalize(line_dir);

			for(j = 0; j < face->nverts; j++)
			{
				vert = Object.vlist[face->verts[j]];

				theta = 30 + (vert->distance/largest_dis)*30;
				magnitude = 0.01 + (vert->distance/largest_dis)*0.02;

				n_vec1.entry[0] = line_dir.entry[0]*cos(theta) - line_dir.entry[1]*sin(theta);
				n_vec1.entry[1] = line_dir.entry[0]*sin(theta) + line_dir.entry[1]*cos(theta);
				
				n_vec2.entry[0] = line_dir.entry[0]*cos(theta) + line_dir.entry[1]*sin(theta);
				n_vec2.entry[1] = -line_dir.entry[0]*sin(theta) + line_dir.entry[1]*cos(theta);

				eva_dis1.entry[0] = vert->x+0.005*n_vec1.entry[0] - trajectories[traj][i].gend[0];
				eva_dis1.entry[1] = vert->y+0.005*n_vec1.entry[1] - trajectories[traj][i].gend[1];

				eva_dis2.entry[0] = vert->x+0.005*n_vec2.entry[0] - trajectories[traj][i].gend[0];
				eva_dis2.entry[1] = vert->y+0.005*n_vec2.entry[1] - trajectories[traj][i].gend[1];

				if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
				{
					////we may create a repelling limit cycle, it means that all the vectors 
					////on the boundary should separate away from the separatrix

					if(length(eva_dis1) < length(eva_dis2))
						vert->vec = magnitude*n_vec2;
						//vert->vec = 0.01*n_vec2;
					else
						vert->vec = magnitude*n_vec1;
						//vert->vec = 0.01*n_vec1;
				}

				else if(singularities[singID].type == SINK || singularities[singID].type == AFOCUS)
				{
					////we may create an attractor limit cycle, it means that all the vectors
					////on the boundary should converge to the separatrix

					if(length(eva_dis1) < length(eva_dis2))
						vert->vec = magnitude*n_vec1;
						//vert->vec = 0.01*n_vec1;
					else
						vert->vec = magnitude*n_vec2;
						//vert->vec = 0.01*n_vec2;
				}

				vert->OnBoundary = 1;
				vert->InRegion = 0;
				vert->repell_flag = 0;
				vert->attract_flag = 0;
			}
		}
	}
}






void GrowAttractRegionForTwoConnections(int singularID, int target_triangle)
{
	////Intialize the region related flags
	int i;
	Face *face;

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];
		face->attract_inregion = 0;
	}

	for(i = 0; i < Object.nverts; i++)
	{
		Object.vlist[i]->attract_flag = 0;
	}

	////initialize attractor region with first input triangle

	if(singularities[singularID].type != SADDLE)
	{
		attractorRegion.trianglelist[0] = singularities[singularID].Triangle_ID;
		Object.flist[singularities[singularID].Triangle_ID]->attract_inregion = 1;
		attractorRegion.num = 1;

		////initialize attractor boundary with the 3 edges of first triangle
		attractorBoundary.edgelist[0] = Object.flist[singularities[singularID].Triangle_ID]->edges[0];
		attractorBoundary.edgelist[1] = Object.flist[singularities[singularID].Triangle_ID]->edges[1];
		attractorBoundary.edgelist[2] = Object.flist[singularities[singularID].Triangle_ID]->edges[2];
		attractorBoundary.num = 3;
	}
	else
	{
		InitSaddleRegionForTwoConnections(singularID, 0, 
			Object.flist[target_triangle]->singularity_index, 1);
	}

	GetRegionNormals(1); ////Get the outward normal of the region for each boundary edge

	////initialize the inner vertices list
	attractorInnerverts.num = 0;    ////no inner vertex at present

	////region grow processing
	Cancel_Growing(1, target_triangle);
	//Cancel_Growing_adv(1, target_triangle);
}



void GrowRepellRegionForTwoConnections(int singularID, int target_triangle)
{
	////initialization part

	////Intialize the region related flags
	int i;
	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->repell_inregion = 0;
	}

	for(i = 0; i < Object.nverts; i++)
	{
		Object.vlist[i]->repell_flag = 0;
	}

	////initialize repeller region with first input triangle
	if(singularities[singularID].type != SADDLE)
	{
		repellerRegion.trianglelist[0] = singularities[singularID].Triangle_ID;
		Object.flist[singularities[singularID].Triangle_ID]->repell_inregion = 1;
		repellerRegion.num = 1;

		////initialize repeller boundary with the 3 edges of first triangle
		repellerBoundary.edgelist[0] = Object.flist[singularities[singularID].Triangle_ID]->edges[0];
		repellerBoundary.edgelist[1] = Object.flist[singularities[singularID].Triangle_ID]->edges[1];
		repellerBoundary.edgelist[2] = Object.flist[singularities[singularID].Triangle_ID]->edges[2];
		repellerBoundary.num = 3;
	}
	
    else
	{
		InitSaddleRegionForTwoConnections(singularID, 1, 
			Object.flist[target_triangle]->singularity_index, 1);
	}

	GetRegionNormals(0); ////Get the outward normal of the region for each boundary edge

	////initialize the inner vertices list
	repellerInnerverts.num = 0;    ////no inner vertex at present

	////region grow processing
	Cancel_Growing(0, target_triangle);
	//Cancel_Growing_adv(0, target_triangle);
}



void InitSaddleRegionForTwoConnections(int saddle, int type, int connect_sing, int initsaddlelength)
{
	int i, j, sep, traj;
	//int dualsep = 0;
	////it seems that we just need to add one triangle that contains the saddle
	if(type == 0) //saddle acts as attractor
	{
		attractorRegion.trianglelist[0] = singularities[saddle].Triangle_ID;
		attractorRegion.num = 1;
		Object.flist[singularities[saddle].Triangle_ID]->attract_inregion = 1;

		//attractorBoundary.edgelist[0] = Object.flist[singularities[saddle].Triangle_ID]->edges[0];
		//attractorBoundary.edgelist[1] = Object.flist[singularities[saddle].Triangle_ID]->edges[1];
		//attractorBoundary.edgelist[2] = Object.flist[singularities[saddle].Triangle_ID]->edges[2];
		//attractorBoundary.num = 3;

		sep = singularities[saddle].separtices;

		for(i = 0; i < 4; i++)
		{
			switch(i)
			{
			case 0:
				traj = separatrices[sep].sep1;
				break;
			case 1:
				traj = separatrices[sep].sep2;
				break;
			case 2:
				traj = separatrices[sep].sep3;
				break;
			case 3:
				traj = separatrices[sep].sep4;
				break;
			}

			////we just add one more at each direction and see what happen
			if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID!=singularities[connect_sing].Triangle_ID)
			{
				for(j = 0; j < num_linesegs_curtraj[traj]; j++)
				{
					if(trajectories[traj][j].Triangle_ID != singularities[saddle].Triangle_ID)
					{
						//attractorRegion.trianglelist[attractorRegion.num ] = 
						//	trajectories[traj][j].Triangle_ID;
						//attractorRegion.num ++;

						//if(attractorRegion.num > 4)
						//	break;
						if(Object.flist[trajectories[traj][j].Triangle_ID]->contain_singularity != 1)
							AddToRegionTriangles(trajectories[traj][j].Triangle_ID, 1);
					}
				}
			}
		}

		UpdateBoundary(1);
	}

	else{ //saddle acts as repeller
		repellerRegion.trianglelist[0] = singularities[saddle].Triangle_ID;
		Object.flist[singularities[saddle].Triangle_ID]->repell_inregion = 1;
		repellerRegion.num = 1;

		//repellerBoundary.edgelist[0] = Object.flist[singularities[saddle].Triangle_ID]->edges[0];
		//repellerBoundary.edgelist[1] = Object.flist[singularities[saddle].Triangle_ID]->edges[1];
		//repellerBoundary.edgelist[2] = Object.flist[singularities[saddle].Triangle_ID]->edges[2];
		//repellerBoundary.num = 3;
		
		sep = singularities[saddle].separtices;

		for(i = 0; i < 4; i++)
		{
			switch(i)
			{
			case 0:
				traj = separatrices[sep].sep1;
				break;
			case 1:
				traj = separatrices[sep].sep2;
				break;
			case 2:
				traj = separatrices[sep].sep3;
				break;
			case 3:
				traj = separatrices[sep].sep4;
				break;
			}

			////we just add one more at each direction and see what happen
			if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID!=singularities[connect_sing].Triangle_ID)
			{
				for(j = 0; j < num_linesegs_curtraj[traj]; j++)
				{
					if(trajectories[traj][j].Triangle_ID != singularities[saddle].Triangle_ID)
					{
						//repellerRegion.trianglelist[repellerRegion.num ] = 
						//	trajectories[traj][j].Triangle_ID;
						//repellerRegion.num ++;
						if(Object.flist[trajectories[traj][j].Triangle_ID]->contain_singularity != 1)
							AddToRegionTriangles(trajectories[traj][j].Triangle_ID, 0);
					}
				}
				//dualsep = 1;
			}
		}

		UpdateBoundary(0);
	}
}


/*************************************************************************************/

////Routines for building and solving a sparse linear system

////build the sparse linear matrix in compact form
void BuildTheSparseLinearSystem2(Vec_DP &tsa, Vec_INT &tija, Mat_DP &bound_v)
{
	int i, j; 
	Vertex *adj_v, *cur_v = NULL;
	Edge *adj_e;
	int *RegionIndex;
	int num_nonzeroarow = 0;
	int num_elements = 0;
	
	////fill the diagnal elements first
	for(i = 0; i < intersectInnerverts.num; i++)
		//tsa[i] = 1 + SMOOTHSTEP;
		tsa[i] = 1.;    ////11/06/05

	num_elements = intersectInnerverts.num + 1;

	////difficult part to store other non-zero, non-diagnal elements
	for(i = 0; i < intersectInnerverts.num; i++)
	{
		cur_v = intersectInnerverts.vertslist[i];
		RegionIndex = new int[cur_v->Num_edge];

		num_nonzeroarow = 0;

		for(j = 0; j < cur_v->Num_edge; j++)
		{
			adj_e = cur_v->edges[j];

			////get the adjacent vertex on the other side of current edge
			if(Object.vlist[adj_e->verts[0]] != cur_v)
				adj_v = Object.vlist[adj_e->verts[0]];
			else
				adj_v = Object.vlist[adj_e->verts[1]];

			////if the adjacent vertex is on boundary, put it to the right side
			if(adj_v->OnBoundary == 1 || adj_v->RegionListID < 0)
			{
				RegionIndex[j] = -1;

				//bound_v[i][0] += (SMOOTHSTEP/cur_v->Num_edge)*adj_v->vec.entry[0];
				//bound_v[i][1] += (SMOOTHSTEP/cur_v->Num_edge)*adj_v->vec.entry[1];
				bound_v[i][0] += (1./cur_v->Num_edge)*adj_v->vec.entry[0];  
				bound_v[i][1] += (1./cur_v->Num_edge)*adj_v->vec.entry[1];
			}
			else
			{
				RegionIndex[j] = adj_v->RegionListID;
				num_nonzeroarow++;   ////add one more non-zero, non-diagnal element in current row
			}
		}

		if(num_nonzeroarow == 0)
		{
			tija[i] = num_elements;
			
		}

		else{
			////sorting the regionindex array
			BubbleSorting(RegionIndex, cur_v->Num_edge);

			////move all non -1 element in the 'RegionIndex' to the front
			if(num_nonzeroarow < cur_v->Num_edge)
			{
				int firstnonfuyielement = cur_v->Num_edge - num_nonzeroarow;

				for(j = 0; j < num_nonzeroarow; j++)
				{
					RegionIndex[j] = RegionIndex[firstnonfuyielement + j];
				}
			}

			////Add elements to the corresponding positions of sa and ija array
			tija[i] = num_elements;

			for(j = 0; j < num_nonzeroarow; j++)
			{
				//tsa[num_elements + j] = -SMOOTHSTEP/cur_v->Num_edge;
				tsa[num_elements + j] = -1./cur_v->Num_edge;   ////11/06/05
				tija[num_elements + j] = RegionIndex[j];
			}

			num_elements += num_nonzeroarow;
		}

		delete [] RegionIndex;
	}

	tsa[intersectInnerverts.num] = 0.;
	tija[intersectInnerverts.num] = num_elements-1;
}

/**************************************************************************************************
////Calling the routine from numerical recipe to solve the sparse line system for region smoothing
**************************************************************************************************/
void Cancel_RegionSmooth()
{
	if(intersectInnerverts.num < 2)
	{
		MessageBox(NULL, "Can not find enough inner vertices", "error", MB_OK);
		return;
	}

	int i;
	Vertex *cur_v;
    int NMAX=10*intersectInnerverts.num;

    Vec_INT tempija(NMAX);
    Vec_DP tempsa(NMAX);

    DP err;
    Vec_DP b(intersectInnerverts.num),bcmp(intersectInnerverts.num),x(intersectInnerverts.num);

	Mat_DP bound_vert(0.0, intersectInnerverts.num, 2);  ////store the vectors that on the boundary vertices

    const int ITOL=1,ITMAX=100;
    const DP TOL=1.0e-10;
    int iter;

	//icVector2 tempv;

	BuildTheSparseLinearSystem2(tempsa, tempija, bound_vert);

    Vec_INT ija(tempija[intersectInnerverts.num]);
    Vec_DP sa(tempija[intersectInnerverts.num]);

	for(i = 0; i < tempija[intersectInnerverts.num]; i++)
	{
		ija[i] = tempija[i];
		sa[i] = tempsa[i];
	}

    ija_p = &ija;
	sa_p = &sa;
     

	////Smoothing the x component of the vectors
	for (i = 0; i < intersectInnerverts.num; i++) 
	{
		cur_v = intersectInnerverts.vertslist[i];
        x[i]=0.0;

		b[i]=bound_vert[i][0];
    }
    linbcg(b,x,ITOL,TOL,ITMAX,iter,err);

	////Store the results back to the vertices
	for (i = 0; i < intersectInnerverts.num; i++) 
	{
		cur_v = intersectInnerverts.vertslist[i];
		cur_v->vec.entry[0] = x[i];
    }

	////Smoothing the y component of the vectors
	for (i = 0; i < intersectInnerverts.num; i++) 
	{
		cur_v = intersectInnerverts.vertslist[i];
        x[i]=0.0;

		b[i]=bound_vert[i][1];
    }
    linbcg(b,x,ITOL,TOL,ITMAX,iter,err);

	////Store the results back to the vertices
	for (i = 0; i < intersectInnerverts.num; i++) 
	{
		cur_v = intersectInnerverts.vertslist[i];
		cur_v->vec.entry[1] = x[i];
    }

	////For correct output, You need to normalize the result vector
    double r;

	for(i = 0; i < intersectInnerverts.num; i++)
	{
		cur_v = intersectInnerverts.vertslist[i];
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
	
	////Store the field after performing smoothing operation
	////May be not enough memory !!!
	//for(i = 0; i < intersectInnerverts.num; i++)
	//{
	//	aftercancel[i] = intersectInnerverts.vertslist[i]->vec;
	//}

}

//void UndoPairCancel()
//{
//	int i;
//
//	for(i = 0; i < intersectInnerverts.num; i++)
//	{
//		intersectInnerverts.vertslist[i]->vec = beforecancel[i];
//	}
//
//}
//
//
//void RedoPairCancel()
//{
//	int i;
//
//	for(i = 0; i < intersectInnerverts.num; i++)
//	{
//		intersectInnerverts.vertslist[i]->vec = aftercancel[i];
//	}
//
//}
//
//

////Testing routine 07/10/05
////Test from only one triangle constituted neighborhood
////The routine should be called after singularities extraction

//void TestingRoutine(int sing1, int sing2)
//{
//	////Here we assume only 2 singularities existing in the field
//	////One is repeller (a source for example), the other is an attractor ( a saddle)
//
//	////From the two triangles that contain the two singularities repectively
//	////Get the intersection region to perform smoothing
//
//	////1.Get the beginning triangles
//	int repell_triangle, attract_triangle;
//	int repell_index, attract_index;
//	int compensate_flag = 0;
//
//	if(singularities[sing1].type == SOURCE && singularities[sing2].type == SADDLE)
//	{
//		repell_triangle = singularities[sing1].Triangle_ID;
//		repell_index = sing1;
//		attract_triangle = singularities[sing2].Triangle_ID;
//		attract_index = sing2;
//	}
//
//	else if(singularities[sing1].type == SADDLE && singularities[sing2].type == SOURCE)
//	{
//		repell_triangle = singularities[sing2].Triangle_ID;
//		repell_index = sing2;
//		attract_triangle = singularities[sing1].Triangle_ID;
//		attract_index = sing1;
//	}
//
//	else if(singularities[sing1].type == SINK && singularities[sing2].type == SADDLE)
//	{
//		repell_triangle = singularities[sing2].Triangle_ID;
//		repell_index = sing2;
//		attract_triangle = singularities[sing1].Triangle_ID;
//		attract_index = sing1;
//	}
//
//	else if(singularities[sing1].type == SADDLE && singularities[sing2].type == SINK)
//	{
//		repell_triangle = singularities[sing1].Triangle_ID;
//		repell_index = sing1;
//		attract_triangle = singularities[sing2].Triangle_ID;
//		attract_index = sing2;
//	}
//
//	else if(singularities[sing1].type == SADDLE && singularities[sing2].type != SADDLE)
//	{
//		repell_triangle = singularities[sing2].Triangle_ID;
//		repell_index = sing2;
//		attract_triangle = singularities[sing1].Triangle_ID;
//		attract_index = sing1;
//
//		RotateFieldBasedOnSing(sing2);
//		compensate_flag = 1;
//	}
//
//	else if(singularities[sing1].type != SADDLE && singularities[sing2].type == SADDLE)
//	{
//		repell_triangle = singularities[sing1].Triangle_ID;
//		repell_index = sing1;
//		attract_triangle = singularities[sing2].Triangle_ID;
//		attract_index = sing2;
//
//		RotateFieldBasedOnSing(sing1);
//		compensate_flag = 1;
//	}
//
//	else
//	{
//		MessageBox(NULL, "not suitable singularities' pair!", "Error", MB_OK);
//		return;
//	}
//
//
//	Cancel_GetRegion(repell_triangle, attract_triangle, repell_index, attract_index);
//
//	UpdateBoundary(0);
//	UpdateBoundary(1);
//
//	////3. Smooth the region
//	Cancel_RegionSmooth();
//
//	if(compensate_flag == 1)
//	{
//		CompensateRotate(compensate_matrix);
//	}
//
//}




/*******************************************************************************/

////Routines for region validation


////Singularity extraction on the intersected region
bool IntersectedRegionSingCapture()
{
	int i, j;
	Face *face;
	int *verts;
	icVector2 vec[3];      //vectors at three vertices
    double  theta[3];      //for storing the angles between two vector for Gauss circle judgement

	double  vec_ang[3];  //the angle for each vector under the polar frame of current triangle
	double  ang_sum;

	////Calculate the Gaussian angle
	for (i = 0; i < intersectRegion.num; i++) {
		face = Object.flist[intersectRegion.trianglelist[i]];
		verts = face->verts;

		ang_sum = 0;

		//For each triangle, calculate the vector for all vertices
		for (j=0; j < face->nverts; j++) {
			////using the vectors stored in the vertices
			vec[j] =  Object.vlist[verts[j]]->vec;
			normalize(vec[j]);

			vec_ang[j] = atan2(vec[j].entry[1], vec[j].entry[0]);
		}

		for(j = 0; j < face->nverts; j++)
		{

			if(j == 0)
				theta[j] = vec_ang[face->nverts-1] - vec_ang[0];
			else
				theta[j] = vec_ang[j-1] - vec_ang[j];

			if( theta[j] < -M_PI)
				theta[j] += 2 * M_PI;
			
			if( theta[j] > M_PI)
				theta[j] -= 2 * M_PI;

			ang_sum += theta[j];
		}

		if(fabs(ang_sum) >= (2 * M_PI - 0.01))
			return true; ////find a singularity
	}

    return false;  ////No singularity inside the region
}


////Count the number of the singularities inside the region
int IntersectedRegionSingCount(int &totalindex)
{
	int i, j;
	int num_sings = 0;
	Face *face;
	int *verts;
	icVector2 vec[3];      //vectors at three vertices
    double  theta[3];      //for storing the angles between two vector for Gauss circle judgement

	double  vec_ang[3];  //the angle for each vector under the polar frame of current triangle
	double  ang_sum;

	////Calculate the Gaussian angle
	for (i = 0; i < intersectRegion.num; i++) {
		face = Object.flist[intersectRegion.trianglelist[i]];
		face->contain_singularity = 0;
		verts = face->verts;

		ang_sum = 0;

		//For each triangle, calculate the vector for all vertices
		for (j=0; j < face->nverts; j++) {
			////using the vectors stored in the vertices
			vec[j] =  Object.vlist[verts[j]]->vec;
			normalize(vec[j]);

			vec_ang[j] = atan2(vec[j].entry[1], vec[j].entry[0]);
			
			////
			if(vec_ang[j] < 0) vec_ang[j] += 2 * M_PI;
		}

		for(j = 0; j < face->nverts; j++)
		{
			
			if(j == 2)
			{
				theta[j] = vec_ang[0] - vec_ang[2];
			}
			else{
				theta[j] = vec_ang[j+1] - vec_ang[j];
			}

			if( theta[j] < -M_PI)
				theta[j] += 2 * M_PI;
			
			if( theta[j] > M_PI)
				theta[j] -= 2 * M_PI;

			ang_sum += theta[j];
		}

		if(fabs(ang_sum) >= (2 * M_PI - 0.1))
		{
			if(ang_sum > 0)
				totalindex++;
			else
				totalindex--;

			num_sings ++;
			//face->contain_singularity = 1;
		}
	}

    return num_sings;  //// return the number of the singularity inside the region
}


////Entry routine for pair cancellation
void PairCancellation(int sing1, int sing2)
{
	////From the two triangles that contain the two singularities repectively
	////Get the intersection region to perform smoothing

	int repell_triangle, attract_triangle;
	int repell_index, attract_index;
	int compensate_flag = 0;
	int repellers, attractors; 

	////initialization part

	Num_MediaNodes = 0;

	////we need to set fences here
	SetFence_LimitCyclesexcept(-1);  //set fences for all limit cycles 08/06/06

	////1.Get the beginning triangles
	if(singularities[sing1].type == SOURCE && singularities[sing2].type == SADDLE)
	{
		repell_triangle = singularities[sing1].Triangle_ID;
		repell_index = sing1;
		attract_triangle = singularities[sing2].Triangle_ID;
		attract_index = sing2;

		////
		repellers  = singularities[sing1].node_index;
		attractors = singularities[sing2].node_index;

		if(TwoConnectOrbits(sing2, sing1))
		{
			PairCancelWithTwoConnections(sing2, sing1);
			return;
		}
	}

	else if(singularities[sing1].type == SADDLE && singularities[sing2].type == SOURCE)
	{
		repell_triangle = singularities[sing2].Triangle_ID;
		repell_index = sing2;
		attract_triangle = singularities[sing1].Triangle_ID;
		attract_index = sing1;

		////
		repellers  = singularities[sing2].node_index;
		attractors = singularities[sing1].node_index;
		
		if(TwoConnectOrbits(sing1, sing2))
		{
			PairCancelWithTwoConnections(sing1, sing2);
			return;
		}
	}

	else if(singularities[sing1].type == SINK && singularities[sing2].type == SADDLE)
	{
		repell_triangle = singularities[sing2].Triangle_ID;
		repell_index = sing2;
		attract_triangle = singularities[sing1].Triangle_ID;
		attract_index = sing1;

		////
		repellers  = singularities[sing2].node_index;
		attractors = singularities[sing1].node_index;
		
		if(TwoConnectOrbits(sing2, sing1))
		{
			PairCancelWithTwoConnections(sing2, sing1);
			return;
		}
	}

	else if(singularities[sing1].type == SADDLE && singularities[sing2].type == SINK)
	{
		repell_triangle = singularities[sing1].Triangle_ID;
		repell_index = sing1;
		attract_triangle = singularities[sing2].Triangle_ID;
		attract_index = sing2;

		////
		repellers  = singularities[sing1].node_index;
		attractors = singularities[sing2].node_index;
		
		if(TwoConnectOrbits(sing1, sing2))
		{
			PairCancelWithTwoConnections(sing1, sing2);
			return;
		}
	}

	else if(singularities[sing1].type == SADDLE && singularities[sing2].type != SADDLE)
	{
		repell_triangle = singularities[sing2].Triangle_ID;
		repell_index = sing2;
		attract_triangle = singularities[sing1].Triangle_ID;
		attract_index = sing1;

		RotateFieldBasedOnSing(sing2);

		////We somehow need to update the information of the being cancelled saddle
		CaptureSing();  //we assume it will not change the order of detected singularites
		CalSeparatrices();
		compensate_flag = 1;
		////
		repellers  = singularities[sing2].node_index;
		attractors = singularities[sing1].node_index;
	}

	else if(singularities[sing1].type != SADDLE && singularities[sing2].type == SADDLE)
	{
		repell_triangle = singularities[sing1].Triangle_ID;
		repell_index = sing1;
		attract_triangle = singularities[sing2].Triangle_ID;
		attract_index = sing2;

		RotateFieldBasedOnSing(sing1);
		
		////We somehow need to update the information of the being cancelled saddle
		//CaptureSing();  //we assume it will not change the order of detected singularites
		CalSeparatrices();

		compensate_flag = 1;
		////
		repellers  = singularities[sing1].node_index;
		attractors = singularities[sing2].node_index;
	}

	////we may also can cancel two centers after reflection!!! 11/06/05
	
	/*-------------------------------------------------------------------------------*/

	else
	{
		MessageBox(NULL, "not suitable singularities' pair, try to rotate the field first!", "Error", MB_OK);
		return;
	}


    if(!AdaptivePairCancel(repell_triangle, attract_triangle, repell_index, attract_index))
	{
		MessageBox(NULL, "can not cancel this pair!", "", MB_OK);
		if(compensate_flag == 1)
			goto LL;
		return;
	}

	//we need to mark the singularity pair as cancelled
	MarkCancel(&repellers, 1,
		&attractors, 1,
		NULL, 0);

	if(compensate_flag == 1)
	{
LL:		CompensateRotate(compensate_matrix);
	}
}


/*-----------------------------------------------------------------------------------------*/



int GetOneConnection(int saddle, int singID)
{
	int i, traj;
	int sep = singularities[saddle].separtices;

	for(i = 0; i < 4; i++)
	{
		switch(i)
		{
		case 0:
			traj = separatrices[sep].sep1;
			break;

		case 1:
			traj = separatrices[sep].sep2;
			break;
			
		case 2:
			traj = separatrices[sep].sep3;
			break;

		case 3:
			traj = separatrices[sep].sep4;
			break;
		}

		if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID == singularities[singID].Triangle_ID)
			return traj;
	}
}


void GetIntersectRegionForTwoConnections(int repellID, int attractID)
{
	int i, j;
	Face *face;
	Edge *cur_edge;
	Vertex *vert;

	////Remove the triangles that contain the two singularities from the intersected region 2/27/06
	//Object.flist[singularities[repellID].Triangle_ID]->attract_inregion = 0;
	//Object.flist[singularities[repellID].Triangle_ID]->repell_inregion = 0;
	//Object.flist[singularities[attractID].Triangle_ID]->attract_inregion = 0;
	//Object.flist[singularities[attractID].Triangle_ID]->repell_inregion = 0;
	
	////Test, we only remove the triangle containing the saddle
	if(singularities[repellID].type == SADDLE)
	{
		Object.flist[singularities[repellID].Triangle_ID]->attract_inregion = 0;
		Object.flist[singularities[repellID].Triangle_ID]->repell_inregion = 0;
		Object.flist[singularities[attractID].Triangle_ID]->attract_inregion = 1;
		Object.flist[singularities[attractID].Triangle_ID]->repell_inregion = 1;
	}
	else
	{
		Object.flist[singularities[attractID].Triangle_ID]->attract_inregion = 0;
		Object.flist[singularities[attractID].Triangle_ID]->repell_inregion = 0;
		Object.flist[singularities[repellID].Triangle_ID]->attract_inregion = 1;
		Object.flist[singularities[repellID].Triangle_ID]->repell_inregion = 1;
	}

	////Get the intersect triangular region
	intersectRegion.num = 0;
	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];
		if(face->repell_inregion == 1 && face->attract_inregion == 1)
		{
			intersectRegion.trianglelist[intersectRegion.num] = face->index;
			intersectRegion.num ++;
		}
		else{
			face->repell_inregion = 0;
			face->attract_inregion = 0;
		}
	}

	////Add the triangles that contain the two singularities into the intersectRegion 2/23/06
	//if(!IsRepeated(intersectRegion.trianglelist, singularities[repellID].Triangle_ID, intersectRegion.num))
	//{
	//	intersectRegion.trianglelist[intersectRegion.num] = singularities[repellID].Triangle_ID;
	//	intersectRegion.num ++;
	//}

	//if(!IsRepeated(intersectRegion.trianglelist, singularities[attractID].Triangle_ID, intersectRegion.num))
	//{
	//	intersectRegion.trianglelist[intersectRegion.num] = singularities[attractID].Triangle_ID;
	//	intersectRegion.num ++;
	//}

	

	/*----------------------------------------------------------*/
	////Add the intermediary points into the region
	////11/19/05
	//for(i = 0; i < Num_MediaNodes; i++)
	//{
	//	int repeated = 0;
	//	for(j = 0; j < intersectRegion.num; j++)
	//	{
	//		if(singularities[MediaNodes[i]].Triangle_ID == intersectRegion.trianglelist[j])
	//		{
	//			repeated = 1;
	//			break;
	//		}
	//	}

	//	if(repeated == 0)
	//	{
	//		intersectRegion.trianglelist[intersectRegion.num] = singularities[MediaNodes[i]].Triangle_ID;
	//		intersectRegion.num ++;
	//	}
	//}
	/*-----------------------------------------------------------*/

	////Reuse the UpdateBoundary routine to build the boundary edge list
	////First, we need to reset the flag of edge visting
	for(i = 0; i < intersectRegion.num; i++)
	{
		face = Object.flist[intersectRegion.trianglelist[i]];
		for(j = 0; j < 3; j++)
		{
			cur_edge =  face->edges[j];
			cur_edge->repell_visited = 0;
			cur_edge->attract_visited = 0;
		}
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
	
	////Set vertices of the triangles that contain the two singularities as other boundaries

	face = Object.flist[singularities[repellID].Triangle_ID];
	for(j = 0; j < face->nverts; j++)
	{
		vert = Object.vlist[face->verts[j]];
		vert->InRegion = 0;
		vert->OnBoundary = 1;
		vert->attract_flag = vert->repell_flag = 0;
	}
	
	face = Object.flist[singularities[attractID].Triangle_ID];
	for(j = 0; j < face->nverts; j++)
	{
		vert = Object.vlist[face->verts[j]];
		vert->InRegion = 0;
		vert->OnBoundary = 1;
		vert->attract_flag = vert->repell_flag = 0;
	}

	////Get the intersection of inner vertices
	intersectInnerverts.num = 0;
	for(i = 0; i < Object.nverts; i++)
	{
		vert = Object.vlist[i];
		if(vert->attract_flag == 2 && vert->repell_flag == 2 && vert->OnBoundary == 0) ////The vertex is inside both region
		{
			intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
			vert->InRegion = 1;
			vert->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
			intersectInnerverts.num ++;
		}
	}

}


/*------------------------------------------------------------------------------*/
/*
Set the boundary triangles
*/
void SetBoundaryTriangles()
{
	int i;
	Face *face;

	////Update the boundary for the intersect region
	intersectBoundary.num = 0;
	UpdateBoundary(2);

	////Initial the whole face
	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->discard = 0;
		Object.flist[i]->inDesignCellCycle = 0;
	}

	////
	for(i = 0; i < intersectBoundary.num; i++)
	{
		face = Object.flist[intersectBoundary.edgelist[i]->tris[0]];

		if(face->attract_inregion == 1 && face->repell_inregion == 1)
		{
			face->discard = 1;
			continue;
		}

		face = Object.flist[intersectBoundary.edgelist[i]->tris[0]];
		if(face->attract_inregion == 1 && face->repell_inregion == 1)
		{
			face->discard = 1;
		}
	}
}

/*
The routine tries to extend the triangle strip containing the "being reversed" separatrix
traj: the index of the separatrix being reversed
length: the percentage of the length of the separatrix being reversed
*/
void ExtendIniStrip(int traj, double length)
{
	int i, j, num_lines;
	int pre_triangle = -1;
	Edge *cur_e;
	Vertex *cur_v;
	Corner *cur_c;

	////Copy the triangles that containing the "being reversed" separatrix
	////Here we reuse the DesignCellCycle variable
	num_lines = (int)num_linesegs_curtraj[traj]*length;
	num_triangles_designcurve = 0;
	for(i = 0; i < num_lines; i++)
	{
		if(pre_triangle == trajectories[traj][i].Triangle_ID)
			continue;

		DesignCurveCellCycle[num_triangles_designcurve] = trajectories[traj][i].Triangle_ID;
		num_triangles_designcurve++;
		Object.flist[trajectories[traj][i].Triangle_ID]->inDesignCellCycle = 1;
		pre_triangle = trajectories[traj][i].Triangle_ID;
	}

	////Using the similar method as we did in the limit cycle creation, 
	////but we try to avoid the boundary triangles
	BuildBoundaryEdgeList(DesignCurveCellCycle, num_triangles_designcurve); ////get current boundary edges list

	////Initialize. Added at 3/9/06
	for(i = 0; i < num_cycleedges; i++)
	{
		cur_e = Cycle_edge[i];

		cur_v = Object.vlist[cur_e->verts[0]];
		cur_v->repell_flag = 0;

		cur_v = Object.vlist[cur_e->verts[1]];
		cur_v->repell_flag = 0;
	}

	////Get the vertices on the current boundaries
	repellerInnerverts.num = 0;
	for(i = 0; i < num_cycleedges; i++)
	{
		cur_e = Cycle_edge[i];

		cur_v = Object.vlist[cur_e->verts[0]];
		AddToInnerVerts(cur_v, 0);
		
		cur_v = Object.vlist[cur_e->verts[1]];
		AddToInnerVerts(cur_v, 0);
	}

	////Extend the cell cycle by adding all the adjacent triangles of the vertices into the cycle list
	////Suppose current vertices on the boundaries are stored in the variable "repellerInnerverts"
	for(i = 0; i < repellerInnerverts.num; i++)
	{
		cur_v = repellerInnerverts.vertslist[i];

		for(j = 0; j < cur_v->Num_corners; j++)
		{
			cur_c = Object.clist[cur_v->Corners[j]];

			////if it is a boundary triangle, skip it!
			if(Object.flist[cur_c->t]->discard == 1)
				continue;

			if(!IsRepeated(DesignCurveCellCycle, cur_c->t, num_triangles_designcurve))
			{
				if(num_triangles_designcurve >= MaxNumTrianglesDesignCurve - 1)
				{
					MaxNumTrianglesDesignCurve += 50;
					DesignCurveCellCycle = (int*)realloc(DesignCurveCellCycle, sizeof(int)*MaxNumTrianglesDesignCurve);
				}

				DesignCurveCellCycle[num_triangles_designcurve] = cur_c->t;
				Object.flist[cur_c->t]->inDesignCellCycle = 1;
				num_triangles_designcurve ++;
			}
		}
	}

	////reset
	for(i = 0; i < num_cycleedges; i++)
	{
		Cycle_edge[i]->OnBoundary = 0;
		Cycle_edge[i]->BoundaryVisited = 0;
	}
	num_cycleedges = 0;

	for(i = 0; i < repellerInnerverts.num; i++)
	{
		////reset the flag
		repellerInnerverts.vertslist[i]->repell_flag = 0;
	}
	repellerInnerverts.num = 0;
}


void SetExtraBoundary3(int saddle, int singID, icVector2 saddle_dir, int traj, double len)
{
	Face *face;
	Vertex *vert, *other_v;
	Edge *cur_e, *temp_e;
	int i, j;
	
	////set the triangle containing the saddle
	face = Object.flist[singularities[saddle].Triangle_ID];;
	for(i = 0; i < face->nverts; i++)
	{
		vert = Object.vlist[face->verts[i]];
		vert->vec = -0.01*saddle_dir;
						
		if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
			vert->vec *= -1;

		vert->OnBoundary = 1;
		vert->InRegion = 0;
		vert->repell_flag = 0;
		vert->attract_flag = 0;
	}


	////set the triangle containing the other singularity
	face = Object.flist[singularities[singID].Triangle_ID];
	for(i = 0; i < face->nverts; i++)
	{
		vert = Object.vlist[face->verts[i]];
		vert->vec = 0.01*saddle_dir;
		
		if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
			vert->vec *= -1;
		
		vert->OnBoundary = 1;
		vert->InRegion = 0;
		vert->repell_flag = 0;
		vert->attract_flag = 0;
	}

	////Initial
	for(i = 0; i < Object.nverts; i++)
	{
		Object.vlist[i]->visited = 0;
	}

	////we need to reset the vectors on the triangles containing the selected separatrix
	////2/27/06
	int numsegs = (int)num_linesegs_curtraj[traj]*len;
	icVector2 line_dir, n_vec1, n_vec2;
	icVector2 eva_dis1, eva_dis2;
	double theta = (30./180.)*M_PI;

	//attractorInnerverts.num = 0; //to store the vertices on the boundary of original triangle strip

	for(i = 0; i < numsegs; i++)
	{
		if(trajectories[traj][i].Triangle_ID == singularities[saddle].Triangle_ID
			|| trajectories[traj][i].Triangle_ID == singularities[singID].Triangle_ID)
			continue;

		face = Object.flist[trajectories[traj][i].Triangle_ID];

		if(face->attract_inregion == 1 || face->repell_inregion == 1)
		{
			face->attract_inregion = face->repell_inregion = 0;

			////Get the approximate direction of the separatrix in the triangle
			line_dir.entry[0] = trajectories[traj][i].gend[0] - trajectories[traj][i].gstart[0];
			line_dir.entry[1] = trajectories[traj][i].gend[1] - trajectories[traj][i].gstart[1];
			
			if(singularities[singID].type == SINK || singularities[singID].type == AFOCUS)
				line_dir *= -1;
			
			normalize(line_dir);
			
			n_vec1.entry[0] = line_dir.entry[0]*cos(theta) - line_dir.entry[1]*sin(theta);
			n_vec1.entry[1] = line_dir.entry[0]*sin(theta) + line_dir.entry[1]*cos(theta);
			
			n_vec2.entry[0] = line_dir.entry[0]*cos(theta) + line_dir.entry[1]*sin(theta);
			n_vec2.entry[1] = -line_dir.entry[0]*sin(theta) + line_dir.entry[1]*cos(theta);

			for(int j = 0; j < face->nverts; j++)
			{
				vert = Object.vlist[face->verts[j]];

				if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
				{
					////we may create a repelling limit cycle, it means that all the vectors 
					////on the boundary should separate away from the separatrix
					eva_dis1.entry[0] = vert->x+0.005*n_vec1.entry[0] - trajectories[traj][i].gend[0];
					eva_dis1.entry[1] = vert->y+0.005*n_vec1.entry[1] - trajectories[traj][i].gend[1];

					eva_dis2.entry[0] = vert->x+0.005*n_vec2.entry[0] - trajectories[traj][i].gend[0];
					eva_dis2.entry[1] = vert->y+0.005*n_vec2.entry[1] - trajectories[traj][i].gend[1];

					if(length(eva_dis1) < length(eva_dis2))
					{
						vert->vec = 0.03*n_vec2;
					}
					else
						vert->vec = 0.03*n_vec1;
				}
				else if(singularities[singID].type == SINK || singularities[singID].type == AFOCUS)
				{
					////we may create an attractor limit cycle, it means that all the vectors
					////on the boundary should converge to the separatrix
					
					eva_dis1.entry[0] = vert->x+0.005*n_vec1.entry[0] - trajectories[traj][i].gend[0];
					eva_dis1.entry[1] = vert->y+0.005*n_vec1.entry[1] - trajectories[traj][i].gend[1];

					eva_dis2.entry[0] = vert->x+0.005*n_vec2.entry[0] - trajectories[traj][i].gend[0];
					eva_dis2.entry[1] = vert->y+0.005*n_vec2.entry[1] - trajectories[traj][i].gend[1];

					if(length(eva_dis1) < length(eva_dis2))
					{
						vert->vec = 0.03*n_vec1;
					}
					else
						vert->vec = 0.03*n_vec2;
				}

				vert->visited = 1;
			}
		}
	}

	//// Get the vector values by average the values on the vertices that are on the boundary of 
	////original triangle strip
	int vec_count = 0;
	icVector2 temp_v;

	BuildBoundaryEdgeList(DesignCurveCellCycle, num_triangles_designcurve); ////get current boundary edges list
	for(i = 0; i < num_cycleedges; i++)
	{
		cur_e = Cycle_edge[i];

		vert = Object.vlist[cur_e->verts[0]];

		vec_count = 0;
		temp_v.set(0,0);

		////We still need to worry about that if the vertex itself is on 
		////the boundary of the original triangle strip
		if(vert->OnBoundary == 0)
		{
			if(vert->visited == 0)
			{
				for(j = 0; j < vert->Num_edge; j++)
				{
					temp_e = vert->edges[j];

					if(temp_e->verts[0] != vert->VertID)
						other_v = Object.vlist[temp_e->verts[0]];
					else
						other_v = Object.vlist[temp_e->verts[1]];

					if(other_v->visited == 1)
					{
						temp_v += other_v->vec;
						vec_count++;
					}
				}

				vert->vec = (1./vec_count)*temp_v; //get the average vector values
			}

			vert->OnBoundary = 1;
			vert->InRegion = 0;
		}

		vert = Object.vlist[cur_e->verts[1]];

		vec_count = 0;
		temp_v.set(0,0);

		////We still need to worry about that if the vertex itself is on 
		////the boundary of the original triangle strip
		if(vert->OnBoundary == 0)
		{
			if(vert->visited == 0)
			{
				for(j = 0; j < vert->Num_edge; j++)
				{
					temp_e = vert->edges[j];

					if(temp_e->verts[0] != vert->VertID)
						other_v = Object.vlist[temp_e->verts[0]];
					else
						other_v = Object.vlist[temp_e->verts[1]];

					if(other_v->visited == 1)
					{
						temp_v += other_v->vec;
						vec_count++;
					}
				}

				vert->vec = (1./vec_count)*temp_v; //get the average vector values
			}

			vert->OnBoundary = 1;
			vert->InRegion = 0;
		}
	}

}

/*------------------------------------------------------------------------------*/


////Pair cancellation for the case that there are two connections between the saddle and 
////the other singularity 2/23/06
void PairCancelWithTwoConnections(int saddle, int singID)
{
	int type;
	int repellID, attractID;
	int repell_triangle, attract_triangle;

	if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
		type = 0; //repeller, thus, saddle should act as attractor

	else if(singularities[singID].type == SINK || singularities[singID].type == AFOCUS)
		type = 1; //attractor, thus, saddle should act as repeller
	////we do not consider other case now

	if(type == 0)
	{
		repellID = singID;
		repell_triangle = singularities[singID].Triangle_ID;

		attractID = saddle;
		attract_triangle = singularities[saddle].Triangle_ID;
	}

	else{
		attractID = singID;
		attract_triangle = singularities[singID].Triangle_ID;

		repellID = saddle;
		repell_triangle = singularities[saddle].Triangle_ID;
	}
	
	////For setting the fences
	int repellers[1] = {singularities[repellID].node_index};
	int num_repellers = 1;
	int attractors[1] = {singularities[attractID].node_index};
	int num_attractors = 1;

 	SetFences(repellers, num_repellers,
		attractors, num_attractors,
		MediaNodes, Num_MediaNodes);


	GrowAttractRegionForTwoConnections(attractID, repell_triangle);
	GrowRepellRegionForTwoConnections(repellID ,attract_triangle);

	int traj = GetOneConnection(saddle, singID);
	icVector2 saddle_dir = GetSaddleTriangleDirection(saddle, traj);

	SetExtraBoundary(saddle, singID, saddle_dir, traj);
	//SetExtraBoundary2(saddle, singID, saddle_dir, traj);

	/*----------------------------------------------------------------------*/
	////Try wider extra boundary
	//SetBoundaryTriangles();
	//ExtendIniStrip(traj, 1);
	//SetExtraBoundary3(saddle, singID, saddle_dir, traj, 1);
	/*----------------------------------------------------------------------*/

	GetIntersectRegionForTwoConnections(repellID, attractID);
	UpdateBoundary(2);

	////At this moment, we just randomly pick a connection, later we can add
	////the interface for user to specify which connection he wants


	Cancel_RegionSmooth();

	NormalizeField();

	MarkCancel(repellers,1, attractors,1, MediaNodes, Num_MediaNodes);
}



/////////////////////////////////////////////

////Multiple cancellation
void MultCancellation(int *repellers, int num_repellers, int *attractors, int num_attractors)
{
	////1. Initialization Part
	if(MediaNodes != NULL)
	{
		free(MediaNodes);
	}
	MediaNodes = (int *)malloc(sizeof(int)*MaxMediaNodes);
	Num_MediaNodes = 0;

	////2. Find the connected intermediary components
	SearchConnectComponents_adv(repellers, num_repellers,
		attractors, num_attractors,
		MediaNodes, Num_MediaNodes);

	
	////New added case for sink+source with double connections between them 07/02/06
	if(Num_MediaNodes == 1)
	{
		//
		IndirectlyDoublyConnectedSingularityPair(graphnodes[repellers[0]].singularityID,
			graphnodes[attractors[0]].singularityID, MediaNodes, Num_MediaNodes);
		goto LL;
	}


	//set fence for all separatrices
	//SetFences(repellers, num_repellers,
	//	attractors, num_attractors,
	//	MediaNodes, Num_MediaNodes);

	//set fences according to the intermediate saddles
	SetFenceForSeps(-1);
	SetFence_LimitCyclesexcept(-1);
	RemoveFenceForIndirectlyDoublyPair(graphnodes[repellers[0]].singularityID, 
		graphnodes[attractors[0]].singularityID);

	for(int i = 0; i < Num_MediaNodes; i++)
	{
		RemoveFenceForSaddle(graphnodes[MediaNodes[i]].singularityID);
	}

	//we need to remove the fences of those separatrices that connect to the two singularities

	//AdaptiveGetMultRegion(repellers, num_repellers, attractors, num_attractors);	
	if(AdaptiveGetMultRegion_new(repellers, num_repellers, attractors, num_attractors))

		////if cancellation is successful, mark those nodes 'cancelled'
		MarkCancel(repellers, num_repellers,
				attractors, num_attractors,
				MediaNodes, Num_MediaNodes);

LL:	free(MediaNodes);
	MediaNodes = NULL;
	Num_MediaNodes = 0;
}


/////////////////////////////////////////////////////////////////
////Routines for region validation

////
bool InEdgesList(Edge **edges, int num, Edge *anedge)
{
	int i;
	for(i = 0; i < num; i++)
	{
		if(edges[i] == anedge)
			return true;
	}
	return false;
}

bool InVertsList(int *verts, int num, int vert_index)
{
	int i;
	for(i = 0; i < num; i++)
	{
		if(verts[i] == vert_index)
			return true;
	}
	return false;
}

////calculate the euler charateristices value of a specific region
int CalEulerValue(int *trianglelist, int num)
{
	int i, j;
	Face *face;
	Edge *cur_edge, **temp_edges;
	//Vertex *vert;
	int *temp_verts;
	int num_edges, num_verts;
	int MaxEdges = 200;
	int MaxVerts = 1000;

	////Initialize part
	num_edges = num_verts = 0;

	////allocate memory for the edges' list and vertices' list
	temp_edges = (Edge **) malloc(sizeof(Edge *) * MaxEdges);
	temp_verts = (int *)malloc(sizeof(int) * MaxVerts);

	for(i = 0; i < num; i++)
	{
		face = Object.flist[trianglelist[i]];
		////count the number of the edges and vertices
		for(j = 0; j < 3; j ++)
		{
			cur_edge = face->edges[j];
			////judge and count the edge
			if(num_edges >= MaxEdges - 1)
			{
				MaxEdges += 100;
				temp_edges = (Edge **)realloc(temp_edges, sizeof(Edge*) * MaxEdges);
			}
			if(!InEdgesList(temp_edges, num_edges, cur_edge))
			{
				////add to the list and count the number of edges
				temp_edges[num_edges] = cur_edge;
				num_edges ++;
			}


			//vert = Object.vlist[face->verts[j]];
			////judge and count the vertices;
			if(num_verts >= MaxVerts - 1)
			{
				MaxVerts += 200;
				temp_verts = (int *)realloc(temp_verts, sizeof(int ) *MaxVerts);
			}
			if(!InVertsList(temp_verts, num_verts, face->verts[j]))
			{
				temp_verts[num_verts] = face->verts[j];
				num_verts ++;
			}
		}
	}

	free(temp_edges);
	free(temp_verts);

	return (num_verts + num - num_edges);
}

/*************************************************************************************/

////Routines for field rotation and reflection

////Routines for automatic rotation for pair cancellation and singularity movement
void RotateFieldBasedOnSing(int sing_index)
{
	int i;
	double s[2];

	////our goal is to rotate the singularity to become a source
	double tanvalue = singularities[sing_index].Jacobian.entry[0][1]/singularities[sing_index].Jacobian.entry[1][1];
	double rot_angle = atan(tanvalue);


	if((singularities[sing_index].Jacobian.entry[0][0]*cos(rot_angle)-singularities[sing_index].Jacobian.entry[1][0]*sin(rot_angle)) < 0
		|| (singularities[sing_index].Jacobian.entry[1][1]*cos(rot_angle)+singularities[sing_index].Jacobian.entry[0][1]*sin(rot_angle)) < 0)
		rot_angle += M_PI;

	pretransform_matrix.setIdentity();
	pretransform_matrix.set(cos(rot_angle), -sin(rot_angle), 0,\
		sin(rot_angle), cos(rot_angle), 0,\
		0, 0, 1);


	compensate_matrix.setIdentity();
	compensate_matrix.set(cos(rot_angle), sin(rot_angle), 0,\
		-sin(rot_angle), cos(rot_angle), 0,\
		0, 0, 1);

	////rotate the whole field
	for(i = 0; i < Object.nverts; i++)
	{
		s[0] = Object.vlist[i]->vec.entry[0];
		s[1] = Object.vlist[i]->vec.entry[1];


		Object.vlist[i]->vec.entry[0] = pretransform_matrix.entry[0][0] * s[0] \
			+ pretransform_matrix.entry[0][1] * s[1];
		Object.vlist[i]->vec.entry[1] = pretransform_matrix.entry[1][0] * s[0] \
			+ pretransform_matrix.entry[1][1] * s[1];
	}

    GetLocalVector();

	////can we just update the singularities we want to cancel??? Yes, we can, think about it!
	CaptureSing();
}

void CompensateRotate(icMatrix3x3 compense_matrix)
{
	////rotate the whole field
	int i;
	double s[2];

	for(i = 0; i < Object.nverts; i++)
	{
		s[0] = Object.vlist[i]->vec.entry[0];
		s[1] = Object.vlist[i]->vec.entry[1];


		Object.vlist[i]->vec.entry[0] = compensate_matrix.entry[0][0] * s[0] \
			+ compensate_matrix.entry[0][1] * s[1];
		Object.vlist[i]->vec.entry[1] = compensate_matrix.entry[1][0] * s[0] \
			+ compensate_matrix.entry[1][1] * s[1];
	}

	compensate_matrix.setIdentity();
	pretransform_matrix.setIdentity();

}

/************************************************************************************************/



////////////////////////////////////////////////////////////////////////////////////
////Testing routine for singularity movement
void Move_TestingRoutine(int source_triangle, int target_triangle, int singularityID, double newx, double newy)
{
	int compensate_flag, reflect_flag;
	compensate_flag = reflect_flag = 0;
	int totalindex = 0;

	////According to the type of the singularity to perform corresponding pre-process
	if(singularities[singularityID].type != SOURCE)
	{
		compensate_flag = 1;

		if(singularities[singularityID].type == SADDLE){
			reflect_flag = 1;
			////Need to reflect the whole field first
			ReflectTheWholeField(0);
			CaptureSing();

			////Then, rotate the field to change it to become a source
			RotateFieldBasedOnSing(singularityID);
		}

		else{
			////Rotate the field
			RotateFieldBasedOnSing(singularityID);
		}
	}


	////1. Get the smoothing region

	////we may need to perform region validation
	////if after movement, more singularities appear, do not allow the smoothing

	Move_GetRegion(source_triangle, singularityID, target_triangle, newx, newy);

	////2. Perform smoothing in the region
	Cancel_RegionSmooth();

	////3. we need to perform region validation here
	////if after smoothing, the region contains more than one singularity, do not allow the movement!
	if(IntersectedRegionSingCount(totalindex)>1) ////contain singularity
	{
		MessageBox(NULL,"Can not move this singularity!", "Error", MB_OK);
		Undo();
		return;
	}


	////4. Normalize the vectors on the second boundary
	Vertex *cur_v;
	double r;
	for(int i = 0; i < 3; i++)
	{
		cur_v = Object.vlist[Object.flist[target_triangle]->verts[i]];
	    r = length(cur_v->vec);
		r *= r;
					
		if (r < 0.0001) 
		{
			r = 0.0001;
			cur_v->vec *= dmax/r; 
		}

	    r = length(cur_v->vec);
		r *= r;

		if (r > dmax*dmax) { 
			r  = sqrt(r); 
			cur_v->vec *= dmax/r; 
		}
	}

	/////Compensate the tranformation of the field
	if(compensate_flag == 1)
	{
		if(reflect_flag == 1)
		{
			CompensateRotate(compensate_matrix);
			ReflectTheWholeField(0);
		}

		else{
			CompensateRotate(compensate_matrix);
		}
	}
}


extern bool IsExitEdge(Edge *cur_edge);
extern bool IsEntranceEdge(Edge *cur_edge);
extern bool IsMixedEdge(Edge *cur_edge);


void TestDisplayRegion(GLenum mode)
{
		Face *face;
		Vertex *vert;
		int i, j;

		/*display all the triangle mesh first*/
		glColor3f(.8,.6, 0);
		for(i=0;i<Object.nfaces;i++)
		{
			face = Object.flist[i];
			glBegin(GL_LINE_LOOP);
			for(j = 0; j < face->nverts; j++)
			{
				vert = Object.vlist[face->verts[j]];
				glVertex2f(vert->x, vert->y);
			}
			glEnd();
		}
		
		////Draw repeller region mesh
		//glColor3f(0, 1, 0);
		//for(i = 0; i < repellerRegion.num; i++)
		//{
		//	face = Object.flist[repellerRegion.trianglelist[i]];
		//	glBegin(GL_LINE_LOOP);
		//	for(j = 0; j < face->nverts; j++)
		//	{
		//		vert = Object.vlist[face->verts[j]];
		//		glVertex2f(vert->x, vert->y);
		//	}
		//	glEnd();
		//}

		//for(i = 0; i < Source_re.num; i++)
		//{
		//	face = Object.flist[Source_re.trianglelist[i]];
		//	glBegin(GL_LINE_LOOP);
		//	for(j = 0; j < face->nverts; j++)
		//	{
		//		vert = Object.vlist[face->verts[j]];
		//		glVertex2f(vert->x, vert->y);
		//	}
		//	glEnd();
		//}
		
		////Draw attractor region mesh
		//glColor3f(0, 0, 1);
		//for(i = 0; i < attractorRegion.num; i++)
		//{
		//	face = Object.flist[attractorRegion.trianglelist[i]];
		//	glBegin(GL_LINE_LOOP);
		//	for(j = 0; j < face->nverts; j++)
		//	{
		//		vert = Object.vlist[face->verts[j]];
		//		glVertex2f(vert->x, vert->y);
		//	}
		//	glEnd();
		//}
		
		//for(i = 0; i < Sink_re.num; i++)
		//{
		//	face = Object.flist[Sink_re.trianglelist[i]];
		//	glBegin(GL_LINE_LOOP);
		//	for(j = 0; j < face->nverts; j++)
		//	{
		//		vert = Object.vlist[face->verts[j]];
		//		glVertex2f(vert->x, vert->y);
		//	}
		//	glEnd();
		//}

		////Draw the testing arrows on vertices of the saddle region
		//glColor3f(1, 1, 0);
		//for(i = 0; i < attractorRegion.num; i++)
		//{
		//	face = Object.flist[attractorRegion.trianglelist[i]];

		//	for(j = 0; j < face->nverts; j++)
		//	{
		//		vert = Object.vlist[face->verts[j]];
		//		glPushMatrix();
		//		glTranslatef(vert->x, vert->y, 0);
		//		glRotatef(atan2(vert->vec.entry[1],vert->vec.entry[0])*360/(2*M_PI), 0, 0, 1);
		//		glScalef(ARROWSCALE/5, ARROWSCALE/5, 1);
		//			glBegin(GL_LINES);
		//			glVertex2f(0, 0);
		//			glVertex2f(1, 0);
		//			glEnd();

		//			////Draw the wings of the arrow
		//			glBegin(GL_LINES);
		//			glVertex2f(1, 0);
		//			glVertex2f(0.8, 0.16);

		//			glVertex2f(1, 0);
		//			glVertex2f(0.8, -0.16);
		//			glEnd();
		//		glPopMatrix();
		//	}
		//}

		////Draw the normal of the edges on boundary
		//UpdateBoundary(/*attractorBoundary, attractorRegion, */1);
		//GetRegionNormals(/*attractorBoundary,*/ 1);
		Edge *cur_edge;
		//glColor3f(0.7, 1, 1);
		//for(i = 0; i < attractorBoundary.num; i++)
		//{
		//	cur_edge = attractorBoundary.edgelist[i];

		//	vert = Object.vlist[cur_edge->verts[0]];
		//	 
		//	glPushMatrix();
		//	glTranslatef(vert->x, vert->y, 0);
		//	glRotatef(atan2(cur_edge->attract_normal.entry[1],\
		//		cur_edge->attract_normal.entry[0])*360/(2*M_PI), 0, 0, 1);
		//	glScalef(ARROWSCALE/5, ARROWSCALE/5, 1);
		//		glBegin(GL_LINES);
		//		glVertex2f(0, 0);
		//		glVertex2f(1, 0);
		//		glEnd();

		//		////Draw the wings of the arrow
		//		glBegin(GL_LINES);
		//		glVertex2f(1, 0);
		//		glVertex2f(0.8, 0.16);

		//		glVertex2f(1, 0);
		//		glVertex2f(0.8, -0.16);
		//		glEnd();
		//	glPopMatrix();
		//}

		////Draw the normal of the edges on boundary
		//Edge *cur_edge;
		//glColor3f(0.7, 1, 1);
		//for(i = 0; i < repellerBoundary.num; i++)
		//{
		//	cur_edge = repellerBoundary.edgelist[i];

		//	vert = Object.vlist[cur_edge->verts[0]];
		//	 
		//	glPushMatrix();
		//	glTranslatef(vert->x, vert->y, 0);
		//	glRotatef(atan2(cur_edge->repell_normal.entry[1],\
		//		cur_edge->repell_normal.entry[0])*360/(2*M_PI), 0, 0, 1);
		//	glScalef(ARROWSCALE/5, ARROWSCALE/5, 1);
		//		glBegin(GL_LINES);
		//		glVertex2f(0, 0);
		//		glVertex2f(1, 0);
		//		glEnd();

		//		////Draw the wings of the arrow
		//		glBegin(GL_LINES);
		//		glVertex2f(1, 0);
		//		glVertex2f(0.8, 0.16);

		//		glVertex2f(1, 0);
		//		glVertex2f(0.8, -0.16);
		//		glEnd();
		//	glPopMatrix();
		//}

       ////Draw the intersected region
		glColor3f(1, 1, 0);
		for(i = 0; i < intersectRegion.num; i++)
		{
			face = Object.flist[intersectRegion.trianglelist[i]];
			glBegin(GL_LINE_LOOP);
			for(j = 0; j < face->nverts; j++)
			{
				vert = Object.vlist[face->verts[j]];
				glVertex2f(vert->x, vert->y);
			}
			glEnd();
		}

		////Drawing boundary
		glColor3f(1, 1, 1);
		GetRegionNormals(2);
		//Edge *cur_edge;
		for(i = 0; i < intersectBoundary.num; i++)
		{
			cur_edge = intersectBoundary.edgelist[i];
			if(RepellerExitEdgePending(cur_edge))
				glColor3f(0, 1, 0);
			else if(AttractorExitEdgePending(cur_edge))
				glColor3f(0, 0, 1);
			else
				glColor3f(1, 0, 0);
			
			glBegin(GL_LINES);
			vert = Object.vlist[intersectBoundary.edgelist[i]->verts[0]];
			glVertex2f(vert->x, vert->y);
			vert = Object.vlist[intersectBoundary.edgelist[i]->verts[1]];
			glVertex2f(vert->x, vert->y);
			glEnd();
		}


		//Draw the boundary of repeller region
		glColor3f(0, 1, 0);
		for(i = 0; i < repellerBoundary.num; i++)
		{
			glBegin(GL_LINES);
			vert = Object.vlist[repellerBoundary.edgelist[i]->verts[0]];
			glVertex2f(vert->x, vert->y);
			vert = Object.vlist[repellerBoundary.edgelist[i]->verts[1]];
			glVertex2f(vert->x, vert->y);
			glEnd();
		}
		
		//Draw the boundary of attractor region
		glColor3f(0, 0, 1);
		for(i = 0; i < attractorBoundary.num; i++)
		{
			glBegin(GL_LINES);
			vert = Object.vlist[attractorBoundary.edgelist[i]->verts[0]];
			glVertex2f(vert->x, vert->y);
			vert = Object.vlist[attractorBoundary.edgelist[i]->verts[1]];
			glVertex2f(vert->x, vert->y);
			glEnd();
		}

		////Mark the inner vertices
		glColor3f(1, 0, 0);
		glPointSize(2.);
		glBegin(GL_POINTS);
		for(i = 0; i < Object.nverts; i++)
		{
			vert = Object.vlist[i];
			if(vert->InRegion == 1)
				glVertex2f(vert->x, vert->y);
		}
		glEnd();

		////Draw arrows on the boundary
		//glColor3f(0.5, 1, 1);
		//for(i = 0; i < intersectBoundary.num; i++)
		//{
		//	cur_edge = intersectBoundary.edgelist[i];

		//	for(j = 0; j < 2; j++)
		//	{
		//		vert = Object.vlist[cur_edge->verts[j]];
		//		 
  //   			glPushMatrix();
		//		glTranslatef(vert->x, vert->y, 0);
		//		glRotatef(atan2(vert->vec.entry[1],\
		//			vert->vec.entry[0])*360/(2*M_PI), 0, 0, 1);
		//		glScalef(ARROWSCALE/5, ARROWSCALE/5, 1);
		//			glBegin(GL_LINES);
		//			glVertex2f(0, 0);
		//			glVertex2f(1, 0);
		//			glEnd();

		//			////Draw the wings of the arrow
		//			glBegin(GL_LINES);
		//			glVertex2f(1, 0);
		//			glVertex2f(0.8, 0.16);

		//			glVertex2f(1, 0);
		//			glVertex2f(0.8, -0.16);
		//			glEnd();
		//		glPopMatrix();
		//	}
		//}

	    ////Draw the triangle containing singularity  2/16/06
		//glColor3f(1, 0, 1);
		//for(i = 0; i < intersectRegion.num; i++)
		//{
		//	face = Object.flist[intersectRegion.trianglelist[i]];

		//	if(face->contain_singularity == 0)
		//		continue;

		//	glBegin(GL_LINE_LOOP);
		//	for(j = 0; j < face->nverts; j++)
		//	{
		//		vert = Object.vlist[face->verts[j]];
		//		glVertex2f(vert->x, vert->y);
		//	}
		//	glEnd();
		//}

		////Draw the boundary of attractor region
		//glColor3f(1, 1, 1);
		//for(i = 0; i < attractorBoundary.num; i++)
		//{
		//	glBegin(GL_LINES);
		//	vert = Object.vlist[attractorBoundary.edgelist[i]->verts[0]];
		//	glVertex2f(vert->x, vert->y);
		//	vert = Object.vlist[attractorBoundary.edgelist[i]->verts[1]];
		//	glVertex2f(vert->x, vert->y);
		//	glEnd();
		//}


		////Draw the boundary of repeller region
		//glColor3f(1, 1, 1);
		//for(i = 0; i < repellerBoundary.num; i++)
		//{
		//	glBegin(GL_LINES);
		//	vert = Object.vlist[repellerBoundary.edgelist[i]->verts[0]];
		//	glVertex2f(vert->x, vert->y);
		//	vert = Object.vlist[repellerBoundary.edgelist[i]->verts[1]];
		//	glVertex2f(vert->x, vert->y);
		//	glEnd();
		//}
}


/******************************************************************************/

////Routines for initialization of the variables for region smoothing

////Initialize the singularities pair cancellation and movement
void InitCancellationAndMovement()
{
	repellerRegion.num = 0;
	repellerBoundary.num = 0;
	repellerInnerverts.num = 0;
	
	attractorRegion.num = 0;
	attractorBoundary.num = 0;
	attractorInnerverts.num = 0;

	intersectRegion.num = 0;
	intersectBoundary.num = 0;
	intersectInnerverts.num = 0;

	////Reset the flags
	int i, j;
	Face *face;
	Edge *cur_edge;

	for(i = 0; i < Object.nverts; i++)
	{
		Object.vlist[i]->attract_flag = 0;
		Object.vlist[i]->repell_flag = 0;
		Object.vlist[i]->RegionListID = -1;
	}

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];
		face->attract_inregion = 0;
		face->repell_inregion = 0;

		//face->fence_flag = 0; //reset all the fence flags

		for(j = 0; j < 3; j++)
		{
			cur_edge = face->edges[j];
			cur_edge->attract_visited = 0;
			cur_edge->repell_visited = 0;
		}
	}
		
	ResetSmoothVars();
}


/*
Remove the fences everywhere
*/
void ClearAllFences()
{
	int i;
	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->fence_flag = 0;
	}
}


void InitMultiRegionVars()
{
	Source_re1.num = 0;
	Source_re2.num = 0;
	Source_re.num = 0;

	Sink_re1.num = 0;
	Sink_re2.num = 0;
	Sink_re.num = 0;
}


/*
Double pointer reallocation may have problem
realloc() routine may not return a correct double pointer
*/
void AllocateVarforTopologyEdit()
{
	////Allocate space for singularities pair cancellation and movement
	MaxNumTriangle = Object.nfaces;   ////modified on 11/06/05
    MaxNumBoundaryEdges = 1000;
    MaxNumInnerVerts = 4000;
	
	repellerRegion.trianglelist = (int *)malloc(sizeof(int) * MaxNumTriangle);
	repellerBoundary.edgelist = (Edge **)malloc(sizeof(Edge *) * MaxNumBoundaryEdges);
	repellerInnerverts.vertslist = (Vertex **)malloc(sizeof(Vertex *) * MaxNumInnerVerts);
	
	attractorRegion.trianglelist = (int *)malloc(sizeof(int) * MaxNumTriangle);
	attractorBoundary.edgelist = (Edge **)malloc(sizeof(Edge *) * MaxNumBoundaryEdges);
	attractorInnerverts.vertslist = (Vertex **)malloc(sizeof(Vertex *) * MaxNumInnerVerts);

	intersectRegion.trianglelist = (int *)malloc(sizeof(int) * MaxNumTriangle);
	intersectBoundary.edgelist = (Edge **)malloc(sizeof(Edge *) * MaxNumBoundaryEdges);
	intersectInnerverts.vertslist = (Vertex **)malloc(sizeof(Vertex *) * MaxNumInnerVerts);


	VerticalField = (icVector2 *)malloc(sizeof(icVector2) * Object.nverts);

	////
}


void AllocateVarforMultiRegion()
{
	Source_re1.trianglelist = (int *)malloc(sizeof(int) * Object.nfaces);
	//Source_re2.trianglelist = (int *)malloc(sizeof(int) * Object.nfaces);
	Source_re.trianglelist = (int *)malloc(sizeof(int) * Object.nfaces);

	Sink_re1.trianglelist = (int *)malloc(sizeof(int) * Object.nfaces);
	//Sink_re2.trianglelist = (int *)malloc(sizeof(int) * Object.nfaces);
	Sink_re.trianglelist = (int *)malloc(sizeof(int) * Object.nfaces);
}


void FinalizeTopologyEdit()
{
	free(repellerRegion.trianglelist);
	free(repellerBoundary.edgelist);
	free(repellerInnerverts.vertslist);

	free(attractorRegion.trianglelist);
	free(attractorBoundary.edgelist);
	free(attractorInnerverts.vertslist);

	free(intersectRegion.trianglelist);
	free(intersectBoundary.edgelist);
	free(intersectInnerverts.vertslist);
}

void FianlizeMultiRegion()
{
	free(Source_re1.trianglelist);
	free(Source_re.trianglelist);

	free(Sink_re1.trianglelist);
	free(Sink_re.trianglelist);
}

/***************************************************************
Extend the number of the edge incident to a specific vertex
The function extend one space each time
***************************************************************/

int *Extend_link(int *edge_link, int Num_edges)
{
    int *temp = edge_link;
	int *temp_edge_link = (int *) malloc(sizeof(int)*(Num_edges + 1));
	if(Num_edges > 0)
	{
		for(int i = 0; i < Num_edges; i++)
			temp_edge_link[i] = temp[i];
		//delete temp;
		free (temp);
	}

	return temp_edge_link;
}


/*
Set boundary flag for the vertices on the input boundary
*/

void SetBoundaryFlag_Ver(Edge **boundaryedgelist, int num_edges)
{
	int i;
	Edge *cur_e;

	for(i = 0; i < num_edges; i++)
	{
		cur_e = boundaryedgelist[i];

		Object.vlist[cur_e->verts[0]]->OnBoundary = 1;
		Object.vlist[cur_e->verts[1]]->OnBoundary = 1;
	}
}


