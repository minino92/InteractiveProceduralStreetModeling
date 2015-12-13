////LimitCycleDetect.cpp

////Defect of this file: too much global variables!!!!!!!!!!Not good

////This module implements the detection of limit cycle using the closed streamline detection method

#include "stdafx.h"

#include "LimitCycleDetect.h"

#include "VFDataStructure.h"

#include "LocalTracing.h"

#include "numerical.h"

#include "topologyedit.h"

extern Polygon3D Object;
extern Singularities *singularities;           //being captured singularites' list
extern int cur_singularity_index;
extern LineSeg **trajectories;                 //trajectories' list
extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory

////Actually, we need to find another efficient data structure to store the whole information of a limit cycle
////for future operations, such as reshape, moving
extern LimitCycle *limitcycles;
extern int cur_limitcycle_index;
extern int MaxNumLimitCycles;


extern double avg_edge_length; //record the average edge length 1/8/06


extern void getVector(double, double, icVector2&, double&);

//////extern routines for outer limit cycle detection
extern void GetPointsOnStreamline(int limitcycle_index, int num_pts);
extern bool InRegion(int vert_id);
extern bool InRegion(double, double);


//////////////////////////////////////////////////////////////////////////////////
////extern variables from topologyedit 07/27/05
extern TriangularRegion repellerRegion;       ////region containing a repeller
extern TriangularRegion attractorRegion;      ////region containing an attractor
extern TriangularRegion intersectRegion;      ////region containing an attractor
extern RegionBoundary repellerBoundary;
extern RegionBoundary attractorBoundary;
extern RegionBoundary intersectBoundary;
extern InnerVertices  repellerInnerverts;
extern InnerVertices  attractorInnerverts;
extern InnerVertices  intersectInnerverts;

extern int MaxNumTriangle ;
extern int MaxNumBoundaryEdges;
extern int MaxNumInnerVerts;


extern void AddToBoundaryEdges(Edge *cur_edge, int type);
extern void AddToInnerVerts(Vertex *vert, int type);
extern bool RepellerExitEdgePending(Edge *cur_edge);
extern bool AttractorExitEdgePending(Edge *cur_edge);
extern int GetOppositeTriangle(Edge *cur_edge, int type);


/////////////////////////////////////////////////////////////////////////////////
extern double problemx, problemy;

/////////for the flow length
extern double sum_flow_length;


////variables for limit cycle extraction
int *TriangleList;
int *cellcycle;
int num_trianglelist;                ////number of triangles in the triangle list
int num_celltriangle;                ////number of triangles in the found cell cycle
Edge **Cycle_edge;                   ////mesh edges of the boundary of the cell cycle
Vertex **PotentialVerts;             ////Potential exit points (vertices) on the boundary of cell cycle
double realexit[2];                  ////we may need the information of vertex!!!!05/30/05
Vertex *realexitvert;                ////currently founded real exit
int num_cycleedges;                  ////number of edges consist of the cell cycle 
int num_potentialverts;              ////number of vertices consist of the cell cycle
Edge *OneSharedEdge;                 ////one shared edge in the cell cycle for finding fixed point
double theFixPoint[2];               ////the fix point on the closed streamline
Edge *SecondSharedEdge;              ////another shared edge for finding next level beginning point

int MaxTriangleInList;
int MaxTriangleInCellCycle;
int MaxEdgeInCycle;
int MaxVertOnBoundary;

////To store multiple limit cycles detection results
int **CellCycleList;
int MaxNumCellCycle;                       ////Maximum number of limit cycles that can be stored
int Cur_CellCycleIndex;                    ////current cell cycle index
int *NumTriangleInEachCellCycle;           ////store the number of triangles in each cell cycle

int pre_cur_traj_index;

//////Variable for limit cycle tangent curve detection
double tang_pre[3], tang_cur[3];  ////They are 3D global coordinates
icVector2 intersect_edge;     // record the direction of the current intersected edge


//double gbeginx, gbeginy;              ////global variables for storing the beginning point

/////////////////
////Testing global variables
double MarkNextBx[2], MarkNextBy[2];
int embedlevel;


int test_numcelledges;

////Global thread variable
HANDLE	thread;

/*-------------------------------------------------------------------------------------------------------------*/
////Some global functions for limit cycle detection that are not suitable to define in limitcycledetect.h file


////may have problem here 05/30/05
bool GetTangentPoint(Edge *e, double &Ex, double &Ey)
{
	double alpha;
	double Tx, Ty;  ////coordinates of the tangent point
	////if we can not find the point, return false
	icVector2 v0, v1, v0andv1;

	v0 = Object.vlist[e->verts[0]]->vec;

	v1 = Object.vlist[e->verts[1]]->vec;

	normalize(v0);
	normalize(v1);

	v0andv1 = v0 - v1;

    normalize(v0andv1);

	double dot1 = dot(v1, v0andv1);

	if(dot1 < 0 )
	{
		v0andv1 *= -1;
	}


	if(fabs(v0.entry[0] - v1.entry[0]) >= 1e-8)
	{
		alpha = (v0andv1.entry[0] - v0.entry[0])/(v0.entry[0] - v1.entry[0]);
	}
	else
		alpha = (v0andv1.entry[1] - v0.entry[1])/(v0.entry[1] - v1.entry[1]);

	if(alpha < 0 || alpha > 1)
	    return false;

	////get the coordinates of tangent point on the edge
	Tx = Object.vlist[e->verts[0]]->x + alpha*(Object.vlist[e->verts[1]]->x - Object.vlist[e->verts[0]]->x);
	Ty = Object.vlist[e->verts[0]]->y + alpha*(Object.vlist[e->verts[1]]->y - Object.vlist[e->verts[0]]->y);

	////Actually, we should let 'E' not the same as the tangent point 'T'
	Ex = Tx;
	Ey = Ty;

	return true;
}


////Tangential potential exits testing
bool TangentExitTest(Edge *e, int repellorattract)
{
	double Ex, Ey;
	Ex = Ey = 0;
	////first, we need to extract the potential exit point 'E' on the edge 'e'
	////if we can not find such a point, return false
	if(!GetTangentPoint(e, Ex, Ey))
	    return false;

	return backwardExitTest(Ex, Ey, repellorattract);
}


bool IsEdgeAlreadyInList(Edge *cur_e, Edge **Cycleedge, int NumofCellEdges)
{
	int i;
	for(i = 0; i < NumofCellEdges; i++)
	{
		if(Cycleedge[i] == cur_e)
			return true;
	}

	return false;
}

////Add the input edge to the cell boundary edge list
////This routine assume that the edge is always added to Cycle_edge list
void AddToBoundaryEdgeList(Edge *cur_e)
{
	if(IsEdgeAlreadyInList(cur_e, Cycle_edge, num_cycleedges))
		return;
	
	if(num_cycleedges >= MaxEdgeInCycle)
	{
		MaxEdgeInCycle += 100;
		Cycle_edge = (Edge **)realloc(Cycle_edge, sizeof(Edge *) * MaxEdgeInCycle);

		if(Cycle_edge == NULL)
		{
			MessageBox(NULL,"fail to reallocate memory for Cycle_edge!", "Error", MB_OK);
			exit(-1);
		}
	}

	Cycle_edge[num_cycleedges] = cur_e;
	num_cycleedges++;
}



bool IsVertAlreadyInList(Vertex *cur_v, Vertex **Potentialverts, int NumofPotentialVerts)
{
	int i;
	for(i = 0; i < NumofPotentialVerts; i++)
	{
		if(Potentialverts[i] == cur_v)
			return true;
	}

	return false;
}

/********************************************
********************************************/
void derivs(const DP t, Vec_I_DP &y, Vec_O_DP &dydx)
{
	icVector2 vec ;
	double mag = 0;

	////first, we need to get the change of next points
	getVector(y[0], y[1], vec, mag);

	dydx[0] = vec.entry[0];
	dydx[1] = vec.entry[1];

}

void inverse_derivs(const DP t, Vec_I_DP &y, Vec_O_DP &dydx)
{
	icVector2 vec ;
	double mag = 0;

	////first, we need to get the change of next points
	getVector(y[0], y[1], vec, mag);

	dydx[0] = -vec.entry[0];
	dydx[1] = -vec.entry[1];

}
/*-------------------------------------------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
////Routines and variables for limit cycle detection

bool TriangleSearch(int *acycle, int num_cycletriangles, int oneTriangle, int &position)
{
	int i;

	for(i = 0; i < num_cycletriangles; i++)
	{
		if(acycle[i] == oneTriangle)
		{
			position = i;
			return true;
		}
	}

	return false;
}


////This routine tends to find a closed cell cycle
bool FindCellCycle(int *trilist, int originNum, int *acycle, int &CellNum, int cur_triangleID)
{
	int position, j;

	position = -1;

	if(!TriangleSearch(trilist, originNum, cur_triangleID, position))
		return false;


	////we find a cell cycle
	////copy the triangles inside the cell cycle to build a new list
	CellNum = 0;
    for( j = 0; j < originNum-position; j++)
	{
		acycle[j] = trilist[position + j];
		CellNum ++;
	}

	////we may need to test whether the two intersections are too close
	////If they are, we may meet a tangent curve here
    //if this point is a fixed point, there may be problem here
	//especailly passing through a vertex 1/5/06
	if(CellNum <= 2)
	{
		icVector2 dis;
		dis.entry[0] = tang_cur[0] - tang_pre[0]; 
		dis.entry[1] = tang_cur[1] - tang_pre[1];

		////testing codes here, write to a temporary file
		//FILE *fp = fopen("tangentpoint.txt", "a");
		//
		//fprintf(fp, "previous point : %f, %f \n", tang_pre[0], tang_pre[1]);
		//fprintf(fp, "previous point : %f, %f \n", tang_cur[0], tang_cur[1]);
		//fprintf(fp, "the distance: %f\n", length(dis));

		//fclose(fp);

		icVector2 t_dis = dis;
		normalize(t_dis);

		///we can treat the two points as the same point similarly
		if(length(dis) < 1e-8 || fabs(dot(t_dis, intersect_edge)) > 1-0.01)  
		{
			CellNum = 0;
			return false;
		}
	}

	return true;
}


/**************************************************************************
Routine for beginning forward/backward tracing of repeller/attractor using local tracing
Entry: x, y are the coordinates of a point very close to the center repeller
       type ---tell routine what kind of tracing scheme it should use
Output: flag records the reason that we exit this routine
***************************************************************************/
	
void local_CellCycleDetect(double &x, double &y, int type, int &flag)
{
	int i;
	int pre_face, cur_face;
	double globalp[2];
	//Face *cur_face;

	////initialize
	////Capture the first triangle
LL:	pre_face = cur_face = TriangleDetect(x, y);

	if(pre_face < 0)
	{
		x += 1e-8;
		goto LL;
	}

	//cur_face = Object.flist[pre_triangle];
	globalp[0] = x;
	globalp[1] = y;

	for(i = 0; i < NUMTRACINGTRIANGLE; i++)
	{
		if(cur_face == -1)
		{
			flag = 2;
			return;
		}

		pre_face = cur_face;
		cur_face = TraceInATriangle(cur_face, globalp, type, flag); ////0 means always forward tracing here

		if(flag == 3 || flag == 4 || pre_face == cur_face) 
		{
			flag = 2;            ////reach a singularity or the boundary
			return;
		}

		if(FindCellCycle(TriangleList, num_trianglelist, cellcycle, num_celltriangle, cur_face))
		{
			flag = 1;

			x = globalp[0];  ////return current global point for possible next cell cycle testing
			y = globalp[1];

			return;
		}

		else{
			TriangleList[num_trianglelist] = cur_face;
			num_trianglelist++;
		}
	}
}



////Using backward tracing to test whether current cell cycle contain a closed streamline

bool closedStreamlineTest(int repellorattract)
{
	int i;
	Vertex *vert;
	Edge *cur_edge;

	////first, testing all the boundary vertices
	for(i = 0; i < num_potentialverts; i++)
	{
		vert = PotentialVerts[i];
		if(backwardExitTest(vert->x, vert->y, repellorattract))
		{
			realexitvert = vert;
			return false;
		}
	}

	////second, testing all possible tagent point on the boundary
	for( i = 0; i < num_cycleedges; i++)
	{
		cur_edge = Cycle_edge[i];
		if(TangentExitTest(cur_edge, repellorattract))
			return false;
	}

	return true;
}


bool backwardExitTest(double x, double y, int repellorattract)
{
	if(local_TravelCellCycleTest(x, y, repellorattract))
	{
		realexit[0] = x;
		realexit[1] = y;

		return true;
	}

	else
		return false;
}

	
////using local tracing to test whether this potential exit is a real exit
bool local_TravelCellCycleTest(double x, double y, int repellorattract)
{
	int i;
	int flag = -1;
	double cur_p[2];
	int pre_face, cur_face;
	//Face *cur_face;

	////variables for local cell cycle testing
	int *localCycle = new int[MaxTriangleInCellCycle];
	int num_localcycletriangles = 0;
	int position;

	////initialize
	cur_p[0] = x;
	cur_p[1] = y;

	////integrate one step to avoid judge vertex

LL:	OneBackwardRungeKutta(cur_p[0], cur_p[1], cur_p[0], cur_p[1], repellorattract); ////backward tracing
	pre_face = cur_face = TriangleDetect(cur_p[0], cur_p[1]); ////get the first triangle

	if(pre_face < 0)
		goto LL;


	/////
	for( i = 0 ; i < NUMTRACINGTRIANGLE ; i ++)
	{
		if(cur_face == -1)
			return false;

		pre_face = cur_face;
		cur_face = TraceInATriangle(cur_face, cur_p, repellorattract, flag);

		if(flag == 3 || flag == 4 || pre_face == cur_face) 
		{   ////reach a singularity or the boundary
			return true;
		}
		
		if(TriangleSearch(cellcycle, num_celltriangle, cur_face, position))
		{
			////store the current triangle into local triangle list
			if(num_localcycletriangles < MaxTriangleInCellCycle)
			{
				if(TriangleSearch(localCycle, num_localcycletriangles, cur_face, position))
				{
					////testing whether local triangle strip closed or not
					if(num_localcycletriangles <= num_celltriangle) ////has same number of triangles as current cycle
					{
						delete [] localCycle;
					    return true;
					}
				}
				////add to local list if not get a close cell cycle
				localCycle[num_localcycletriangles] = cur_face;
				num_localcycletriangles ++;
			
     		}
		}

		else{
			delete [] localCycle;
			return false;
		}
	}
	
	delete [] localCycle;  ////this should not reach!!!!!
	return false;
}



////after we locate a cell cycle, using binary search to find the 
////accurate closed streamline
void GettheClosedStreamline(int repellorattract)
{
	////first, we need to choose an inner edge in the cell cycle (one of the sharing edge of the triangle strip)
	////We have already stored it in the OneSharedEdge
	double v1[2], v2[2];
	Vertex *vert1, *vert0;

	////initialize
	vert0 = Object.vlist[OneSharedEdge->verts[0]];
	vert1 = Object.vlist[OneSharedEdge->verts[1]];

	v1[0] = vert0->x;
	v1[1] = vert0->y;

	v2[0] = vert1->x;
	v2[1] = vert1->y;

	////begin binary search
	BinaryFindtheFixedPoint(v1, v2, repellorattract);
}


	
void BinaryFindtheFixedPoint(double v1[2], double v2[2], int type)
{
	int i = 0;  
	double alpha;
	double midx, midy, interx, intery;
	icVector2 dis;
	int flag = 0;

	interx = intery = 0;

	int triangleID = OneSharedEdge->tris[0];  //modified at 1/5/06

	while(i < 200)    ////avoid dead loop
	{
		flag = 0;
		midx = (v1[0]+v2[0])/2.;
	    midy = (v1[1]+v2[1])/2.;

 		local_GetNextIntersection(midx, midy, triangleID, interx, intery, type, flag); //modified at 1/5/06

		if(flag == 1)
		{
		    MessageBox(NULL, "Can not find next intersection!", "error", MB_OK);
			return;
		}

		if(flag == 2) ////simply use current point as the fixed point
		{
			theFixPoint[0] = interx;
			theFixPoint[1] = intery;
			return;
		}

		dis.entry[0] = interx - midx;
		dis.entry[1] = intery - midy;

		if(length(dis) < INTERSECTIONERROR) ////we find the fixed point
		{
			theFixPoint[0] = interx;
			theFixPoint[1] = intery;
			return;
		}

		////calculate the alpha value
		////Modified at 1/5/06
		//if((v2[0] - v1[0]) != 0)
		//	alpha = (interx - v1[0])/(v2[0] - v1[0]);
		//else
		//	alpha = (intery - v1[1])/(v2[1] - v1[1]);
		icVector2 v1_m, v1_v2;
		v1_m.entry[0] = interx - v1[0];
		v1_m.entry[1] = intery - v1[1];
		v1_v2.entry[0] = v2[0] - v1[0];
		v1_v2.entry[1] = v2[1] - v1[1];

		if(length(v1_v2) < 1e-9)
		{
			theFixPoint[0] = midx;
			theFixPoint[1] = midy;
			return;
		}

		else{
			alpha = length(v1_m)/length(v1_v2);
		}

		////get next searching half line segment
		if(alpha > 0.5)
		{
			v1[0] = midx;
			v1[1] = midy;
		}

		else if(alpha < -1e-8)
		{
			//MessageBox(NULL, "Wrong test: not a closed streamline!", "error", MB_OK);
			theFixPoint[0] = midx;
			theFixPoint[1] = midy;
			return;
		}

		else{
			v2[0] = midx;
			v2[1] = midy;
		}

		i++;
	}

}

	
////To reduce the calculation of intersections, we may need to store previous intersection 
////and compare it with current intersection, here is the temporary global variable for this purpose
////1/6/06
double first_intersection[2] = {0.};
int first_intersect = 0;




////Using local tracing to detect the intersection of trajectory and the chosen edge
////Modified at 1/5/06
void  local_GetNextIntersection(double beginx, double beginy, int thetriangle,
								double &interx, double &intery, 
								int type, int &flag)
{
	int i, localflag = 0;
	double /*pre_p[2],*/ cur_p[2];
	int cur_face, pre_face;
	//double t[2];
	icVector2 v_length;
	icVector2 v0, v1; ////two end points of the edges

	cur_p[0] = cur_p[1] = 0;

	pre_face = cur_face = thetriangle;

	////previous tracing part
	for(i = 0; i < NUMTRACINGTRIANGLE ; i++)
	{
		if(cur_face == -1)
		{
			flag = 1;
			return;
		}

    	pre_face = cur_face;
		cur_face = TraceInATriangle2(cur_face, cur_p, type, localflag);

		if(localflag == 3 || localflag == 4 /*|| pre_face == cur_face*/)
		{
			flag = 1;
			return;
		}

		if(cur_face == thetriangle) //modified at 1/5/06
			break;

	}
	
	////calculate intersection part

	////here we use global tracing to get the intersection 
	v0.entry[0] = Object.vlist[OneSharedEdge->verts[0]]->x;
	v0.entry[1] = Object.vlist[OneSharedEdge->verts[0]]->y;
	v1.entry[0] = Object.vlist[OneSharedEdge->verts[1]]->x;
	v1.entry[1] = Object.vlist[OneSharedEdge->verts[1]]->y;

	double A, B, C;
	A = v0.entry[1] - v1.entry[1];
	B = v1.entry[0] - v0.entry[0];
	C = (v0.entry[0]*v1.entry[1] - v1.entry[0]*v0.entry[1]);

	double pending = A*cur_p[0] + B*cur_p[1] + C;

	/*----Modified at 1/5/06------*/
	////if the intersection falls in the shared edge
	////this is the intersection we want
	if(fabs(pending) < 1e-9) ////passing the vertex
	{
		interx = cur_p[0];
		intery = cur_p[1];
		goto LL;
	}
	////else, we perform one more tracing in current triangle till we go out of it
	////now the intersection is what we want
	cur_face = TraceInATriangle2(cur_face, cur_p, type, localflag);
	interx = cur_p[0];
	intery = cur_p[1];

	////New added codes 1/6/06
LL:	if(first_intersect == 0)
	{
		first_intersection[0] = interx;
		first_intersection[1] = intery;
		first_intersect++;
	}

	if(first_intersect == 1)
	{
		if(first_intersection[0]== interx && first_intersection[1] == intery)
		{
			flag = 3;  //we probably reach the fixed point
		}
		else
		{
			first_intersection[0] = interx;
			first_intersection[1] = intery;
		}
	}


}


////if previous cell cycle does not contain closed streamline
////Get next testing beginning point for finding new cell cycle
void GetNextTestBeginPoint(double &bx, double &by)
{
	bx = realexit[0];
	by = realexit[1];
}


////Judge whether it is a good shared edge, a temporary method 1/2/06
////because not all limit cycles will have a specific center !!!!
bool GoodSharedEdge(Edge *the_edge, double center[2])
{
	icVector2 vec0, vec1;
	//vec0 = Object.vlist[the_edge->verts[0]]->vec;
	//vec1 = Object.vlist[the_edge->verts[1]]->vec;
	vec0.entry[0] = Object.vlist[the_edge->verts[0]]->x - center[0];
	vec0.entry[1] = Object.vlist[the_edge->verts[0]]->y - center[1];
	
	vec1.entry[0] = Object.vlist[the_edge->verts[1]]->x - center[0];
	vec1.entry[1] = Object.vlist[the_edge->verts[1]]->y - center[1];

	normalize(vec0);
	normalize(vec1);
	
	double theta = acos(dot(vec0, vec1));

	if(theta < M_PI/18) return true;
	else
		return false;
}

////Building the edges list consist of the cell cycle
void boundaryBuilding()
{
	////search the direct former and later triangles of current triangle
	////find the sharing edges of these two pairs of triangle, mark the two edges of current triangle
	////the remaining edge must be on the boundary

	int i, j;
	int former, later;
	int adj_f;
	Face *curf;
	Edge *cur_e;

	for(i = 0; i < num_celltriangle; i++)
	{
		curf = Object.flist[cellcycle[i]];
		former = cellcycle[((i-1)+num_celltriangle)%num_celltriangle];
        later = cellcycle[(i+1)%num_celltriangle];

		for(j = 0; j < curf->nverts; j++)
		{
			cur_e = curf->edges[j];

			////Get the adjacent face to current face cellcycle[i]
			if(cur_e->tris[0] != cellcycle[i])
				adj_f = cur_e->tris[0];
			else
				adj_f = cur_e->tris[1];

			if(adj_f == former || adj_f == later)  ////this edge is an inner edge, not add to the boundary list
			{
				//if( i == 0)
				//if(GoodSharedEdge(cur_e) || i == 0)
				//    OneSharedEdge = cur_e;             ////recalled the first shared edge inside cell cycle

				//if( i == 4)
				//	SecondSharedEdge = cur_e;

				continue;
			}

			else{
				AddToBoundaryEdgeList(cur_e);
			}
		}
	}
}


void GetASharedEdge(double center[2])
{
	int i, j;
	int former, later;
	int adj_f;
	Face *curf;
	Edge *cur_e;

	for(i = 0; i < num_celltriangle; i++)
	{
		curf = Object.flist[cellcycle[i]];
		former = cellcycle[((i-1)+num_celltriangle)%num_celltriangle];
        later = cellcycle[(i+1)%num_celltriangle];

		for(j = 0; j < curf->nverts; j++)
		{
			cur_e = curf->edges[j];

			////Get the adjacent face to current face cellcycle[i]
			if(cur_e->tris[0] != cellcycle[i])
				adj_f = cur_e->tris[0];
			else
				adj_f = cur_e->tris[1];

			if(adj_f == former || adj_f == later)  ////this edge is an inner edge, not add to the boundary list
			{
				if(GoodSharedEdge(cur_e, center) || i == 0)
				    OneSharedEdge = cur_e;             ////recalled the first shared edge inside cell cycle

				//if( i == 4)
				//	SecondSharedEdge = cur_e;

				continue;
			}
		}
	}

}

	
////Building the vertices list consist of the boundary
////Note: this routine must be called after building the boundary edges list
void PotentialVertsBuilding()
{
	//if(num_celltriangle < 3)
	//	return;
   
	int i;

	for(i = 0; i < num_cycleedges; i++)
	{
		AddToPotentialVertList(Cycle_edge[i]->verts[0]);
		
		AddToPotentialVertList(Cycle_edge[i]->verts[1]);
	}
}




void AddToPotentialVertList(int cur_v_id)
{
	Vertex *cur_v = Object.vlist[cur_v_id];

	if(IsVertAlreadyInList(cur_v, PotentialVerts, num_potentialverts))
		return;

	if(num_potentialverts >= MaxVertOnBoundary - 1)
	{
		MaxVertOnBoundary += 100;
		realloc((Vertex **)PotentialVerts, sizeof(Vertex*) * MaxVertOnBoundary);
	}

	PotentialVerts[num_potentialverts] = cur_v;
	num_potentialverts++;
}



////routines for one step Runge Kutta using global driver!!!
void OneBackwardRungeKutta(double x, double y, double &nextx, double &nexty, int repellorattract)
{
	////calling Runge Kutta to get the next point
    const int N=2;
    int i,j;
    DP eps,hdid,hnext,htry, t = 1.;

    Vec_DP by(N),dydx(N),dysav(N),ysav(N),yscal(N);

	ysav[0] = x;
	ysav[1] = y;

	if(repellorattract == 1)
        derivs(t, ysav, dysav);
	else
		inverse_derivs(t, ysav, dysav);

    for (i=0;i<N;i++) yscal[i]=1.0;
    htry=0.7;


	/*--------------------------------------------------------------------*/
	////calling adaptive stepsize runge-kutta method here to get next step
	for (i=0;i<5;i++) {
		eps=exp(-DP(i+1));
		t = 1.0;
		for (j=0;j<N;j++) {
			by[j]=ysav[j];
			dydx[j]=dysav[j];
		}
					
		if(repellorattract == 1)
		    rkqs(by,dydx,t,htry,eps,yscal,hdid,hnext, derivs);
		else
		    rkqs(by,dydx,t,htry,eps,yscal,hdid,hnext, inverse_derivs);
			 
	}
	/////store some important information
	nextx = by[0];
	nexty = by[1];

	htry = hnext;
	if(htry > 2.)
		htry = 0.7;
	/*--------------------------------------------------------------------*/
}

void  StoreCurrentCellCycle(int *acellcycle, int num)
{
	//CellCycleList[Cur_CellCycleIndex] = new int[num_celltriangle];

	limitcycles[cur_limitcycle_index].cellcycle = (int*)malloc(sizeof(int) * num);

	for(int i = 0; i < num; i++)
	{
		//CellCycleList[Cur_CellCycleIndex][i] = cellcycle[i];

		limitcycles[cur_limitcycle_index].cellcycle[i] = acellcycle[i];
	}

	//NumTriangleInEachCellCycle[Cur_CellCycleIndex] = num_celltriangle;

	limitcycles[cur_limitcycle_index].num_triangles = num;

}
	
void StoreCurStreamline(LineSeg *streamline, int num)
{
	limitcycles[cur_limitcycle_index].closed_streamline = (LineSeg*)malloc(sizeof(LineSeg)*num);

	if(limitcycles[cur_limitcycle_index].closed_streamline == NULL)
	{
		MessageBox(NULL, "Not enough memory!", "Error", MB_OK);
		return;
	}

	for(int i = 0; i < num; i++)
	{
		limitcycles[cur_limitcycle_index].closed_streamline[i].gstart[0] = streamline[i].gstart[0];
		limitcycles[cur_limitcycle_index].closed_streamline[i].gstart[1] = streamline[i].gstart[1];
		limitcycles[cur_limitcycle_index].closed_streamline[i].gstart[2] = streamline[i].gstart[2];

		limitcycles[cur_limitcycle_index].closed_streamline[i].gend[0] = streamline[i].gend[0];
		limitcycles[cur_limitcycle_index].closed_streamline[i].gend[1] = streamline[i].gend[1];
		limitcycles[cur_limitcycle_index].closed_streamline[i].gend[2] = streamline[i].gend[2];

		limitcycles[cur_limitcycle_index].closed_streamline[i].start[0] = streamline[i].start[0];
		limitcycles[cur_limitcycle_index].closed_streamline[i].start[1] = streamline[i].start[1];
		limitcycles[cur_limitcycle_index].closed_streamline[i].start[2] = streamline[i].start[2];

		limitcycles[cur_limitcycle_index].closed_streamline[i].end[0] = streamline[i].end[0];
		limitcycles[cur_limitcycle_index].closed_streamline[i].end[1] = streamline[i].end[1];
		limitcycles[cur_limitcycle_index].closed_streamline[i].end[2] = streamline[i].end[2];

		limitcycles[cur_limitcycle_index].closed_streamline[i].length = streamline[i].length; //store the length

		limitcycles[cur_limitcycle_index].closed_streamline[i].Triangle_ID = streamline[i].Triangle_ID;
	}

	limitcycles[cur_limitcycle_index].num_linesegs = num;
}


////Judge the further vertex on the chosen shared edge
////return the vertex that much further to the specific singularity
Vertex * GetFurtherVer(Edge *sharededge, int singID)
{
	icVector2 dis1, dis2;
	dis1.entry[0] = Object.vlist[sharededge->verts[0]]->x - singularities[singID].gcx;
	dis1.entry[1] = Object.vlist[sharededge->verts[0]]->y - singularities[singID].gcy;

	dis2.entry[0] = Object.vlist[sharededge->verts[1]]->x - singularities[singID].gcx;
	dis2.entry[1] = Object.vlist[sharededge->verts[1]]->y - singularities[singID].gcy;

	if(length(dis1) > length(dis2)) return(Object.vlist[sharededge->verts[0]]);
	else return(Object.vlist[sharededge->verts[1]]);
}


////Main routine for limit cycle detection
void LimitCycleDetect()
{
	/*------------ variables for cell cycle detection ----------*/
	int i;
	double beginx, beginy;
	int type = 0;
	int flag = -1;
	Vertex *ver1;

	int inner = 0;

	int deadloopdetect = 0;

	int index_separatrices = 0;

	int Saddleornot = 0;

	double center[2] = {0.};  //variable for a good shared edge selection to calculate the fixed point 1/2/06

	/*---------- variables for closed streamline display ---------*/
	double streamlinex, streamliney;
	int triangleid;
	//Face *cur_face;

	////Initialize
	beginx = beginy = 0;
    cur_limitcycle_index = 0;
	Cur_CellCycleIndex = 0;
    
	cur_traj_index = pre_cur_traj_index;


	/////Testing codes here 07/25/05
	MarkNextBx[0] = MarkNextBy[0] = 0;
	MarkNextBx[1] = MarkNextBy[1] = 0;


	for(i = 0; i < cur_singularity_index; i++)
	{
        flag = -1;
		inner = 0;
		deadloopdetect = 0;

		/*----------------------------------------*/
		////11/10/05
		index_separatrices = 0;
		Saddleornot = 0;

		/*----------------------------------------*/
		////a limit cycle can be reached from saddle !!!!
		if( singularities[i].type == CWCENTER || singularities[i].type == CCWCENTER)
		{
			continue;
		}
		else{
			if(singularities[i].type == SOURCE || singularities[i].type == RFOCUS)
				type = 0;
			else if(singularities[i].type == SINK || singularities[i].type == AFOCUS)
				type = 1;
			else { //it is a saddle
				Saddleornot = 1;
			}
		}

		////The problem of this selection of beginning point!!!
		////we need to consider the eigenvector direction of saddles!!!! 11/09/05

L2:		if(Saddleornot == 1)
		{
			switch(index_separatrices)
			{
			case 0: beginx = singularities[i].gcx + singularities[i].outgoing.entry[0]*SEPARATRIXSTEP;
				    beginy = singularities[i].gcy + singularities[i].outgoing.entry[1]*SEPARATRIXSTEP; 
					type = 0;
					break;

			case 1: beginx = singularities[i].gcx - singularities[i].outgoing.entry[0]*SEPARATRIXSTEP;
				    beginy = singularities[i].gcy - singularities[i].outgoing.entry[1]*SEPARATRIXSTEP; 
					type = 0;
					break;
			
			case 2: beginx = singularities[i].gcx + singularities[i].incoming.entry[0]*SEPARATRIXSTEP;
				    beginy = singularities[i].gcy + singularities[i].incoming.entry[1]*SEPARATRIXSTEP; 
					type = 1;
					break;

			case 3: beginx = singularities[i].gcx - singularities[i].incoming.entry[0]*SEPARATRIXSTEP;
				    beginy = singularities[i].gcy - singularities[i].incoming.entry[1]*SEPARATRIXSTEP; 
					type = 1;
					break;
			}

			index_separatrices ++;
		}

		else{
			beginx = singularities[i].gcx + SEPARATRIXSTEP*1;  //you need to change this method
			beginy = singularities[i].gcy + SEPARATRIXSTEP*1;
		}

LL:		InitLimitCycleDetection();  ////intialize the cell cycle detection data structures
		flag = -1;

		triangleid = TriangleDetect(beginx, beginy);
		Double_CellCycleDetect(beginx, beginy, triangleid, type, flag);

		if(flag == 2) 
		{
			if(Saddleornot == 0)
			    continue;
			else if(index_separatrices < 3)
				goto L2;
				
		}

		if(flag == 1)  ////we find a closed cellcycle
		{
			////First, we need to check whether it is the limit cycle being found before!!
            ////We can simple compare the cell cycle with those of previous limit cycles
			int pre_limit_index = -1;

			if(IsPreviousLimitCycle(pre_limit_index))
			{ 
				//we need to get another beginning point or just simple move to next singularity?

				//if current testing singularity is a saddle, we may need to store the index of the 
				//previously founded limit cycle to the saddle and add the saddle to the limit cycle
				//variable 10/13/05
				if(singularities[i].type == SADDLE)
				{
					////1)Add the index of the limit cycle 'pre_limit_index' to the list of the saddle
                    UpdateListInSaddle(i, pre_limit_index);
					////2)Add the index of the saddle to the list of limit cycle
					UpdateListInLimitCycle(pre_limit_index, i);
				}

				if(inner > 0)
				{
					////Add the pre_limit_index to the connected list of limit cycle found at previous step
					UpdateCycleListInLimitCycle(i, pre_limit_index); ////not sure!!!!!! 10/14/05
				}

				if(Saddleornot == 1)
				{
					if(index_separatrices < 3)
						goto L2;
				}
				
				continue;

			}
			



			////build the boundary of the cell cycle, find out all the potential exits
			boundaryBuilding();


			////Testing code, save the current number of edges in the cell cycle 08/25/05
			test_numcelledges = num_cycleedges;

			////Deal with the streamline tangent to the shared edge of two triangles 08/28/05
			////Always falls into dead loop because of the choosing of next "beginx, beginy"
			//if(num_celltriangle < 3)  
			//{
			//	////if we have only two triangles in the cell cycle, it may be the tangential case!
			//	if(deadloopdetect >= 6)
			//		continue;

			//	PotentialVertsBuilding();
			//	if(!closedStreamlineTest(type))
			//	{
			//		////Get next beginning point
			//		//beginx = Object.vlist[OneSharedEdge->verts[0]]->x;
			//		//beginy = Object.vlist[OneSharedEdge->verts[0]]->y;

			//		//beginx = gbeginx;
			//		//beginy = gbeginy;

			//		if(type == 0)
			//		    OneBackwardRungeKutta(beginx, beginy, beginx, beginy, 1);
			//		else
			//		    OneBackwardRungeKutta(beginx, beginy, beginx, beginy, 0);

			//		deadloopdetect++;
			//		goto LL;
			//	}
			//}


			////we continue to test all the potential exits
			//if(closedStreamlineTest(type))
			//{
				////we need to locate the accurate position of the closed streamline and display it
			    ////At the same time, we need to store the related information of the found limit cycle

				center[0] = singularities[i].gcx;
				center[1] = singularities[i].gcy;
				GetASharedEdge(center);
			    
				////store the cell cycle to corresponding limit cycle data structure
				StoreCurrentCellCycle(cellcycle, num_celltriangle);  
				Cur_CellCycleIndex++;


				GettheClosedStreamline(type);

				/*------------------------------------------------*/

				////Using local tracing to calculate and store the streamline
				streamlinex = theFixPoint[0];
				streamliney = theFixPoint[1];
				//if(type == 0)
				//{
				//	OneBackwardRungeKutta(streamlinex, streamliney, streamlinex, streamliney, 1);
				//}
				//
				//else{
				//	OneBackwardRungeKutta(streamlinex, streamliney, streamlinex, streamliney, 0);
				//}


				if(streamlinex == 0&& streamliney == 0)
				{
					triangleid = TriangleDetect(streamlinex, streamliney);
				}
				else
					triangleid = OneSharedEdge->tris[0];

				CalLocalTracing(triangleid, streamlinex, streamliney, type);

				////store the streamline to corresponding limit cycle data structure
				StoreCurStreamline(trajectories[cur_traj_index], num_linesegs_curtraj[cur_traj_index]);     

				limitcycles[cur_limitcycle_index].singularID = i;  ////store the center singularity of current limit cycle
				
				limitcycles[cur_limitcycle_index].node_index = -1; ////10/16/05

				/*---------------------------------------------------------------*/
				////Initialize the connected list for the graph
				limitcycles[cur_limitcycle_index].connected_limitcycle = NULL;
				limitcycles[cur_limitcycle_index].num_connectedcycles = 0;
				limitcycles[cur_limitcycle_index].connected_saddle = NULL;
				limitcycles[cur_limitcycle_index].num_connectedsaddles = 0;
				/*---------------------------------------------------------------*/

				////store the type of the limit cycle
				if(type == 0)
					limitcycles[cur_limitcycle_index].type = 1;
				else
					limitcycles[cur_limitcycle_index].type = 0;

				/*---------------------------------------------------------------*/
				////build/add to the connected list
				if(singularities[i].type == SADDLE)
				{
					////1)Add the index of the limit cycle 'pre_limit_index' to the list of the saddle
                    UpdateListInSaddle(i, cur_limitcycle_index);
					////2)Add the index of the saddle to the list of limit cycle
					UpdateListInLimitCycle(cur_limitcycle_index, i);
				}

				if(inner > 0)
				{
					////Add the pre_limit_index to the connected list of limit cycle found at previous step
					UpdateCycleListInLimitCycle(i, pre_limit_index); ////not sure!!!!!! 10/14/05
				}
				/*---------------------------------------------------------------*/


				/*---------------------------------------------------------------*/
				////create the handle for the limit cycle for later interative operation
			    ////Using the next level beginning point as the legend center 08/25/05
				ver1 = Object.vlist[OneSharedEdge->verts[0]];
				limitcycles[cur_limitcycle_index].legend_center[0] = ver1->x;
				limitcycles[cur_limitcycle_index].legend_center[1] = ver1->y;

				////Get the base of the legend
				if(ver1->VertID != OneSharedEdge->verts[0])
					ver1 = Object.vlist[OneSharedEdge->verts[0]];
				else
					ver1 = Object.vlist[OneSharedEdge->verts[1]];
				limitcycles[cur_limitcycle_index].legend_base[0] = ver1->x;
				limitcycles[cur_limitcycle_index].legend_base[1] = ver1->y;
				
				cur_traj_index ++;
				cur_limitcycle_index++;



				/*----------------------------------------------------------*/
				////if we have judged all the separatrices of the saddle keep going
                if(index_separatrices < 3 && singularities[i].type == SADDLE)
					goto L2;

				/*------------------------------------------------*/
                ////Another method to judge next beginning point
               
				ver1 = GetFurtherVer(OneSharedEdge, i);
				beginx = ver1->x;
				beginy = ver1->y;
				

				icVector2 temp_v;
				temp_v.entry[0] = beginx - singularities[i].gcx;
				temp_v.entry[1] = beginy - singularities[i].gcy;
				normalize(temp_v);
				beginx += 0.01 * temp_v.entry[0];  //we should use adaptive stepsize 1/2/06
				beginy += 0.01 * temp_v.entry[1];

                OneBackwardRungeKutta(beginx, beginy, beginx, beginy, type);



			/*------------------------------------------------*/


 				////Inverse the type of the limit cycle
				if(type == 0)
					type = 1;
				else
					type = 0;

				////Build the Conley relation graph 10/13/05
				if(inner == 0)
				{
					if(singularities[i].type == SADDLE) 
					{////if this is a single limit cycle can be reached from the saddle
						
						////1)Add the index of the limit cycle to the list of the saddle
						UpdateListInSaddle(i, cur_limitcycle_index);

					    ////2)Add the index of the saddle to the list of limit cycle
						UpdateListInLimitCycle(cur_limitcycle_index, i);
				    }
				}
				else{
					////if this is a limit cycle that contains previous limit cycle
					////we need to build the connection between them (update the informaiton of both limit cycle)
					UpdateCycleListInLimitCycle(cur_limitcycle_index-1, cur_limitcycle_index-2);
					UpdateCycleListInLimitCycle(cur_limitcycle_index-2, cur_limitcycle_index-1);
				}

				////Perform cell cycle detection in the outer space
				inner ++;
				goto LL;
			//}
		}
	}

	if(cur_limitcycle_index < 1){
		num_linesegs_curtraj[cur_traj_index] = 0;
		MessageBox(NULL, "No limit cycle has been found!", "error", MB_OK);
	}
}


bool IsPreviousLimitCycle(int &limit_index)
{
	if(cur_limitcycle_index == 0)
		return false;

	int i, j;
	int position;
	int cur_triangle;
	int count;

	for(i = 0; i < cur_limitcycle_index; i++)
	{
		////compare current cell cycle with previous cell cycles of the limit cycles being found before
		count = 0;

		for(j = 0; j < num_celltriangle; j++)
		{
			cur_triangle = cellcycle[j];

			if(!TriangleSearch(limitcycles[i].cellcycle, limitcycles[i].num_triangles, cur_triangle, position))
				goto LL;   ////continue to test next limit cycle

			////if it can be found
			count ++;
		}

		//if(count == limitcycles[i].num_triangles)
		if(count == limitcycles[i].num_triangles - 1) ////using a loosen pending condition
		{
			limit_index = i;
			return true;
		}


		else
LL:			continue;
	}

	return false;
}


////using region growing to locate limit cycle
void LimitCycleDetectNew()
{
	/*------------ variables for cell cycle detection ----------*/
	int i;
	double beginx, beginy;
	int type = 0;
	int flag = -1;
	Edge *chosen_edge;
	Vertex *ver1, *ver2;

	/*---------- variables for closed streamline display ---------*/
	double streamlinex, streamliney;
	int triangleid;
	//Face *cur_face;

	////Initialize
	beginx = beginy = 0;
    cur_limitcycle_index = 0;
	Cur_CellCycleIndex = 0;
    
	cur_traj_index = pre_cur_traj_index;


	/////Testing codes here 07/25/05
	MarkNextBx[0] = MarkNextBy[0] = 0;
	MarkNextBx[1] = MarkNextBy[1] = 0;


	for(i = 0; i < cur_singularity_index; i++)
	{
        flag = -1;

		if(singularities[i].type == SADDLE || singularities[i].type == CWCENTER || singularities[i].type == CCWCENTER)
		{
			continue;
		}
		else{
			if(singularities[i].type == SOURCE || singularities[i].type == RFOCUS)
				type = 0;
			else
				type = 1;
		}

		beginx = singularities[i].gcx + SEPARATRIXSTEP*10;
		beginy = singularities[i].gcy + SEPARATRIXSTEP*10;
        
		InitLimitCycleDetection();  ////intialize the cell cycle detection data structures
		flag = -1;

		////perform region growing here

		GetLimitCycleRegion(i, type, 0);
		////we still need some mechanism to judge whether it is a center!!!!


		////choose the middle point of one edge on the boundary of the region to perform cell cycle detection
		UpdateBoundary(type);

		if(type == 0)
		{
			chosen_edge = repellerBoundary.edgelist[0];
		}
		else{
			chosen_edge = attractorBoundary.edgelist[0];
		}
		ver1 = Object.vlist[chosen_edge->verts[0]];
		ver2 = Object.vlist[chosen_edge->verts[1]];
		
		beginx = (ver1->x + ver2->x) / 2;
		beginy = (ver1->y + ver2->y) / 2;

		////Move away from the edge a little bit
		//if(type == 0)
		//	OneBackwardRungeKutta(beginx, beginy, beginx, beginy, 0);
		//else
		//	OneBackwardRungeKutta(beginx, beginy, beginx, beginy, 1);

		////////////////////////////////////
		////New method to get the beginning point
		icVector2 temp_v;
		temp_v.entry[0] = beginx - singularities[i].gcx;
		temp_v.entry[1] = beginy - singularities[i].gcy;
		beginx += 0.01 * temp_v.entry[0];
		beginy += 0.01 * temp_v.entry[1];

		problemx = beginx;
		problemy = beginy;

		//Double_CellCycleDetect(beginx, beginy, type, flag);
		//if(type == 0)
		//    Double_CellCycleDetect(beginx, beginy, 1, flag);
		//else
		//    Double_CellCycleDetect(beginx, beginy, 0, flag);

		if(flag == 2) 
		{
			continue;  ////reaches singularity or boundary
		}

		if(flag == 1)  ////we find a closed cellcycle
		{
			//////build the boundary of the cell cycle, find out all the potential exits
			boundaryBuilding();  ////we need this step to get a shared edge

			////we need to locate the accurate position of the closed streamline and display it
			////At the same time, we need to store the related information of the found limit cycle

			////store the cell cycle to corresponding limit cycle data structure
			StoreCurrentCellCycle(cellcycle, num_celltriangle);  
			Cur_CellCycleIndex++;

			GettheClosedStreamline(type);

			/*------------------------------------------------*/

			////Using local tracing to calculate and store the streamline
			if(type == 0)
			{
				streamlinex = theFixPoint[0];
				streamliney = theFixPoint[1];

				OneBackwardRungeKutta(streamlinex, streamliney, streamlinex, streamliney, 1);
				OneBackwardRungeKutta(streamlinex, streamliney, streamlinex, streamliney, 1);
				triangleid = TriangleDetect(streamlinex, streamliney);

				CalLocalTracing(triangleid, streamlinex, streamliney, type);
			}
			
			else{
				streamlinex = theFixPoint[0];
				streamliney = theFixPoint[1];

				OneBackwardRungeKutta(streamlinex, streamliney, streamlinex, streamliney, 0);
				OneBackwardRungeKutta(streamlinex, streamliney, streamlinex, streamliney, 0);
				triangleid = TriangleDetect(streamlinex, streamliney);

				CalLocalTracing(triangleid, streamlinex, streamliney, type);
			}

			////store the streamline to corresponding limit cycle data structure
			StoreCurStreamline(trajectories[cur_traj_index], num_linesegs_curtraj[cur_traj_index]);     

			limitcycles[cur_limitcycle_index].singularID = i;  ////store the center singularity of current limit cycle
			////store the type of the limit cycle
			if(type == 0)
				limitcycles[cur_limitcycle_index].type = 1;
			else
				limitcycles[cur_limitcycle_index].type = 0;
			
			cur_traj_index ++;
			cur_limitcycle_index++;
		}
	}

	if(cur_limitcycle_index < 1){
		num_linesegs_curtraj[cur_traj_index] = 0;
		MessageBox(NULL, "No limit cycle has been found!", "error", MB_OK);
	}
}


	
////Multiple rounds testing for locate a closed cell cycle accurately	
//Modified at 1/4/06
void  Double_CellCycleDetect(double &x, double &y,  int &begin_triangle, int type, int &flag)
{
	int *local_cellcycle = new int[MaxTriangleInCellCycle];
	int num_localcell = 0;

	//int thirdtest = 0;

	int loopcontrol = 0;

	double gx = x;
	double gy = y;

	//double pre_gx, pre_gy;
	icVector2 dis;

	////first cell cycle test
LL:	First_CellCycleDetect(gx, gy, begin_triangle, local_cellcycle, num_localcell, type, flag);


	if(flag == 2)
	{
		return;
	}
 
	if(flag == 1)  ////found a cell cycle, retest it
	{
    	if(Second_CellCycleTest(gx, gy, type, local_cellcycle, num_localcell, flag))
		{
    		flag = 1;

			////store it to global cellcycle
			for(int k = 0; k < num_localcell; k++)
				cellcycle[k] = local_cellcycle[k];
			num_celltriangle = num_localcell;

			x = gx;
			y = gy;
			return;
		}

		else{ ////if it can not pass second test, it means that the streamline will leave previous cell cycle

			if(flag == 2) ////if reach a singularity or boundary of the domain, need not to trace any more
			{
				x = gx;
				y = gy;
				return;
			}

			if(loopcontrol > 5)
			{
				x = gx;
				y = gy;
				return;
			}

			InitLimitCycleDetection();
			loopcontrol++;
			goto LL;
		}
	}
}
	

void  First_CellCycleDetect(double &x, double &y, int &begin_triangle, 
							int *acycle, int &num_localcell, int type, int &flag)
{
	int i;
	int pre_face, cur_face;
	double globalp[2];
	int position;
	//Face *cur_face;

	////initialize

	////Capture the first triangle
LL:	pre_face = cur_face = TriangleDetect(x, y);

	if(pre_face == -1)
	{
		x += 1e-7;
		goto LL;
	}
	globalp[0] = x;
	globalp[1] = y;

	//pre_face = cur_face = begin_triangle;

	////We need to add the first triangle
	TriangleList[num_trianglelist] = cur_face;
	num_trianglelist++;

	for(i = 0; i < NUMTRACINGTRIANGLE; i++)
	{
		if(cur_face == -1)
			return;

		pre_face = cur_face ;
		cur_face = TraceInATriangle2(cur_face, globalp, type, flag);
 
		if(flag == 3 || flag == 4 || pre_face == cur_face) 
		{
			flag = 2;            ////reach a singularity or the boundary
			return;
		}

		if(FindCellCycle(TriangleList, num_trianglelist, acycle, num_localcell, cur_face))
		{
			flag = 1;
			x = globalp[0];
			y = globalp[1];
			return;
		}

		else{
			if(TriangleSearch(TriangleList, num_trianglelist, cur_face, position))
		      continue;

			if(num_trianglelist >= MaxTriangleInList - 1)
			{
				MaxTriangleInList += 100;
				TriangleList = (int*)realloc(TriangleList, sizeof(int)*MaxTriangleInList);
                cellcycle = (int*)realloc(cellcycle, sizeof(int)*MaxTriangleInList); 
			}

			TriangleList[num_trianglelist] = cur_face;
			num_trianglelist++;
		}
	}
}


bool  Second_CellCycleTest(double &x, double &y, int type, int *localcycle, int num_localcell, int &flag)
{
	int i;
	double globalp[2];
	int pre_face, cur_face;
	int num_passtriangle = 0;

	//// initialize
	pre_face = cur_face = localcycle[0];  ////get the first triangle
	num_passtriangle++;

	globalp[0] = x;
	globalp[1] = y;

	for(i = 0; i < NUMTRACINGTRIANGLE; i++)
	{
		////Match triangle one by one
		if(cur_face != localcycle[i])
		{
			x = globalp[0];  //Modified at 1/5/06
			y = globalp[1];
			return false;
		}

		if(cur_face == -1)
			return false;

		pre_face = cur_face ;
		cur_face = TraceInATriangle2(cur_face, globalp, type, flag);
 
		if(flag == 3 || flag == 4 || pre_face == cur_face) 
		{
			flag = 2;            ////reach a singularity or the boundary
			return false;
		}

		num_passtriangle++;

		if(num_passtriangle == num_localcell)
			return true;
 
	}
}


/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/
////The new cell cycle detection basic routines
////The difference between the new method from the original one is
////the detection of triangle, we hope that we can remove the TriangleDetect() routine
////These routines will be associated with the region growing method
void  First_CellCycleDetect_2(double &x, double &y, int &begin_triangle, 
							int *acycle, int &num_localcell, int type, int &flag)
{
	int i;
	int pre_face, cur_face;
	double globalp[2];
	int position;

	globalp[0] = x;
	globalp[1] = y;
	pre_face = cur_face = begin_triangle;

	////We need to add the first triangle
	TriangleList[num_trianglelist] = cur_face;
	num_trianglelist++;

	for(i = 0; i < NUMTRACINGTRIANGLE; i++)
	{
		if(cur_face == -1)
			return;

		pre_face = cur_face ;
		cur_face = TraceInATriangle2(cur_face, globalp, type, flag);
 
		if(flag == 3 || flag == 4 || pre_face == cur_face) 
		{
			flag = 2;            ////reach a singularity or the boundary
			return;
		}

		if(FindCellCycle(TriangleList, num_trianglelist, acycle, num_localcell, cur_face))
		{
			flag = 1;
			x = globalp[0];
			y = globalp[1];
			return;
		}

		else{
			if(TriangleSearch(TriangleList, num_trianglelist, cur_face, position))
		      continue;

			if(num_trianglelist >= MaxTriangleInList - 1)
			{
				MaxTriangleInList += 100;
				TriangleList = (int*)realloc(TriangleList, sizeof(int)*MaxTriangleInList);
                cellcycle = (int*)realloc(cellcycle, sizeof(int)*MaxTriangleInList); 
			}

			TriangleList[num_trianglelist] = cur_face;
			num_trianglelist++;
		}
	}
}


bool  Second_CellCycleTest_2(double &x, double &y, int type, int *localcycle, 
						   int num_localcell, int &cur_Triangle, int &flag)
{
	int i;
	double globalp[2];
	int pre_face, cur_face;
	int num_passtriangle = 0;

	//// initialize
	pre_face = cur_face = localcycle[0];  ////get the first triangle
	num_passtriangle++;

	globalp[0] = x;
	globalp[1] = y;

	for(i = 0; i < NUMTRACINGTRIANGLE; i++)
	{
		////Match triangle one by one
		if(cur_face != localcycle[i])
		{
			x = globalp[0];  //Modified at 1/5/06
			y = globalp[1];
			cur_Triangle = cur_face;   //Modified at 1/17/06
			return false;
		}

		if(cur_face == -1)
			return false;

		pre_face = cur_face ;
		cur_face = TraceInATriangle2(cur_face, globalp, type, flag);
 
		if(flag == 3 || flag == 4 || pre_face == cur_face) 
		{
			flag = 2;            ////reach a singularity or the boundary
			return false;
		}

		num_passtriangle++;

		if(num_passtriangle == num_localcell)
			return true;
 
	}
}


////New routine to test whether the two arrays have the same elements in the same order 1/19/06
////The two arrays should have the same number of elements
bool HavetheSameElems(int *a1, int *a2, int num)
{
	int i;
	int a2_beginpos = 0;

	for(i = 0; i < num; i++)
	{
		if(a1[0] == a2[i])
		{
			a2_beginpos = i;
			break;
		}
	}

	if(a2_beginpos >= num) ////we can not find such element in the array a2
		return false;

	////ok, now let's compare each element in the two arrays
	for(i = 0; i < num; i++)
	{
		if(a1[i] != a2[(i+a2_beginpos)%num])
			return false;
	}

	return true;
}

void  Double_CellCycleDetect_2(double &x, double &y,  int &begin_triangle, int type, int &flag)
{
	int *local_cellcycle = new int[MaxTriangleInCellCycle];
	int num_localcell = 0;

	int loopcontrol = 0;

	double gx = x;
	double gy = y;

	//double pre_gx, pre_gy;
	icVector2 dis;

	////New added variables for dealing with the numerical issue 1/19/06
	int *pre_cellcycle = new int[MaxTriangleInCellCycle];
	int pre_num_localcell = 0;

	////first cell cycle test
LL:	
	/*----------------------------------------------------*/
	////save the previous cell cycle firstly 1/19/06
	for(int i = 0; i < num_localcell; i++)
	{
		pre_cellcycle[i] = local_cellcycle[i];
	}
	pre_num_localcell = num_localcell;
	/*----------------------------------------------------*/
	
	num_localcell = 0;
	First_CellCycleDetect_2(gx, gy, begin_triangle, local_cellcycle, num_localcell, type, flag);


	if(flag == 2)
	{
		return;
	}
 
	if(flag == 1)  ////found a cell cycle, retest it
	{
    	if(Second_CellCycleTest_2(gx, gy, type, local_cellcycle, num_localcell, begin_triangle, flag))
		{
L2:    		flag = 1;

			////store it to global cellcycle
			for(int k = 0; k < num_localcell; k++)
				cellcycle[k] = local_cellcycle[k];
			num_celltriangle = num_localcell;

			x = gx;
			y = gy;
			return;
		}

		else{ ////if it can not pass second test, it means that the streamline will leave previous cell cycle

			////if it is just some numerical issue, perform another test here 1/19/06
			if(num_localcell == pre_num_localcell)
			{
				////if it is the same cell cycle as before, return the previously detected cell cycle
				if(HavetheSameElems(pre_cellcycle, local_cellcycle, num_localcell))
				{
					goto L2;
				}
			}

			if(flag == 2) ////if reach a singularity or boundary of the domain, need not to trace any more
			{
				x = gx;
				y = gy;
				return;
			}

			if(loopcontrol > 5)
			{
				x = gx;
				y = gy;
				return;
			}

			InitLimitCycleDetection();
			loopcontrol++;
			goto LL;
		}
	}
}

/*-------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/

////Get the next detection beginning point
void GetNextLevelBeginPoint(double &next_x, double &next_y, int limitcycle_index, int &flag)
{
	Vertex *vert, *other_vert;
	double test_x, test_y;
	icVector2 edgevec;
	double scaler = 0.01;
	
	test_x = test_y = 0;

	////1. build the the testing region of current found limit cycle
	int num_pts = 4;
LL:	GetPointsOnStreamline(limitcycle_index, num_pts);

	//if()
	//{
	//	vert = Object.vlist[OneSharedEdge->verts[0]];
	//	other_vert = Object.vlist[OneSharedEdge->verts[1]];

	//	edgevec.entry[0] = other_vert->x - vert->x;
	//	edgevec.entry[1] = other_vert->y - vert->y;

	//	next_x = vert->x + scaler * edgevec.entry[0]; ////pick the point between two ending point
	//	next_y = vert->y + scaler * edgevec.entry[1]; ////it must be not inside the region

	//	return;
	//}
	if(SecondSharedEdge == NULL && OneSharedEdge == NULL)
		return;

	if(SecondSharedEdge == NULL)
		SecondSharedEdge = OneSharedEdge;

	////2. Test the two ending points of the stored shared edge to find the vertex not inside the region
	if((InRegion(SecondSharedEdge->verts[0]) && InRegion(SecondSharedEdge->verts[1]))
		||(!InRegion(SecondSharedEdge->verts[0]) && !InRegion(SecondSharedEdge->verts[1])))
	{
		if(num_pts >= limitcycles[limitcycle_index].num_linesegs)
		{
			flag = 1;
			MessageBox(NULL, "Wrong edge! Can not find the next beginning point!", "Error", MB_OK);
			return;
		}
		num_pts ++;
		goto LL;
	}

	else if(InRegion(SecondSharedEdge->verts[0])&& !InRegion(SecondSharedEdge->verts[1]))
	{
		vert = Object.vlist[SecondSharedEdge->verts[1]];
		other_vert = Object.vlist[SecondSharedEdge->verts[0]];
	}

	else if(InRegion(SecondSharedEdge->verts[1])&& !InRegion(SecondSharedEdge->verts[0]))
	{
		vert = Object.vlist[SecondSharedEdge->verts[0]];
		other_vert = Object.vlist[SecondSharedEdge->verts[1]];
	}

	////3. Pick a point on the edge and close to the vertex found above as the next beginning point
	edgevec.entry[0] = other_vert->x - vert->x;
	edgevec.entry[1] = other_vert->y - vert->y;

	test_x = vert->x - scaler * edgevec.entry[0];
	test_y = vert->y - scaler * edgevec.entry[1];

	//if(InRegion(test_x, test_y))
	//{
	//	scaler *= 0.5;
	//	goto LP;
	//}

	next_x = test_x;
	next_y = test_y;
}


////intialize the limit cycle detection
////this routine will not allocate new memory for the cell cycle detector
////just clear all the previous information
void  InitLimitCycleDetection()
{
	num_trianglelist = 0;
	num_celltriangle = 0;
	num_cycleedges = 0;
    num_potentialverts = 0;
	Cur_CellCycleIndex = 0;

	OneSharedEdge = NULL;
	theFixPoint[0] = theFixPoint[1] = 0.;

	repellerRegion.num = 0;
	attractorRegion.num = 0;
	intersectRegion.num = 0;
	repellerBoundary.num = 0;
	attractorBoundary.num = 0;
	intersectBoundary.num = 0;
	repellerInnerverts.num = 0;
	attractorInnerverts.num = 0;
	intersectInnerverts.num = 0;

	///Reset the region growing flags
	for(int i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->attract_inregion = 0;
		Object.flist[i]->repell_inregion = 0;
		Object.flist[i]->inDesignCellCycle = 0;
	}

}


void ResetLimitCycles()
{
	num_trianglelist = 0;
	num_celltriangle = 0;
	num_cycleedges = 0;
    num_potentialverts = 0;
	Cur_CellCycleIndex = 0;

	OneSharedEdge = NULL;
	theFixPoint[0] = theFixPoint[1] = 0.;

	repellerRegion.num = 0;
	attractorRegion.num = 0;
	intersectRegion.num = 0;
	repellerBoundary.num = 0;
	attractorBoundary.num = 0;
	intersectBoundary.num = 0;
	repellerInnerverts.num = 0;
	attractorInnerverts.num = 0;
	intersectInnerverts.num = 0;

	///Reset the region growing flags
	for(int i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->attract_inregion = 0;
		Object.flist[i]->repell_inregion = 0;
		Object.flist[i]->inDesignCellCycle = 0;
	}

	////Initialize the variables in the limit cycle data structure
	////There are some memeory issues here if I set the pointers as NULL
	//for(int i = 0; i < MaxNumLimitCycles; i++)
	for(int i = 0; i < cur_limitcycle_index; i++)
	{
		if(limitcycles[i].cellcycle != NULL)
		{
			free(limitcycles[i].cellcycle);
			limitcycles[i].cellcycle = NULL;
		}
		limitcycles[i].num_triangles = 0;

		if(limitcycles[i].closed_streamline != NULL)
		{
			free(limitcycles[i].closed_streamline);
			limitcycles[i].closed_streamline = NULL;
		}
		limitcycles[i].num_linesegs = 0;
		limitcycles[i].singularID = -1;
		limitcycles[i].type = -1;

		if(limitcycles[i].connected_limitcycle != NULL)
		{
			free(limitcycles[i].connected_limitcycle);
			limitcycles[i].connected_limitcycle = NULL;
		}
		limitcycles[i].num_connectedcycles = 0;
		
		if(limitcycles[i].connected_saddle != NULL)
		{
			free(limitcycles[i].connected_saddle);
			limitcycles[i].connected_saddle = NULL;
		}
		limitcycles[i].num_connectedsaddles = 0;
	}
}

void InitLimitCycleStructure()
{
	cur_limitcycle_index = 0;

	//Initialize the variables in the limit cycle data structure
	for(int i = 0; i < MaxNumLimitCycles; i++)
	{
		limitcycles[i].cellcycle = NULL;
		limitcycles[i].num_triangles = 0;
		limitcycles[i].closed_streamline = NULL;
		limitcycles[i].num_linesegs = 0;
		limitcycles[i].singularID = -1;
		limitcycles[i].type = -1;
	}
}



void AllocVarforLimitCycleDetect()
{
	////put to limit cycle initial codes
	MaxTriangleInList = 300;
	MaxTriangleInCellCycle = 300;
	MaxEdgeInCycle = 500;
	MaxVertOnBoundary = 1000;
    MaxNumCellCycle = 50;

	TriangleList = (int*) malloc(sizeof(int) * MaxTriangleInList);
	cellcycle = (int*) malloc(sizeof(int) * MaxTriangleInList);
	num_trianglelist = 0;
	num_celltriangle = 0;
	Cycle_edge = (Edge **) malloc(sizeof(Edge*)*MaxEdgeInCycle);                  
	//PotentialVerts = (Vertex**) malloc(sizeof(Vertex*)*MaxVertOnBoundary);           
	num_cycleedges = 0;
    num_potentialverts = 0;
	cur_limitcycle_index = 0;

	pre_cur_traj_index = 0;

	CellCycleList = (int **) malloc(sizeof(int *) * MaxNumCellCycle);
    NumTriangleInEachCellCycle = (int*) malloc(sizeof(int) * MaxNumCellCycle);


	Cur_CellCycleIndex = 0;
	OneSharedEdge = NULL;
	theFixPoint[0] = theFixPoint[1] = 0;
	
	//TypeofLimitCycles = new int[MaxNumCellCycle];   ////the type of corresponding limit cycles

}

////release memory allocate to variables of limit cycle detection
void finalizeLimitCycleDetect()
{
	free(TriangleList);
	free(cellcycle);

	free(Cycle_edge);
	//free(PotentialVerts);

	free(CellCycleList);
	free(NumTriangleInEachCellCycle);
}






//////////////////////////////////////////////////////////////////
////New method to detect limit cycle using Conley region growing
void GrowLimitCycleRegion(int type)
{
	int i;
	bool exitornot = false;
	int num_edges;
	int oppositeTriangle;
	Edge **edgelist;
	Edge *cur_edge;
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
				return;
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

void GetLimitCycleRegion(int singID, int type, int inner)
{

	////we can not always grow from the triangle that contains the center singularity of the limit cycle
	if(inner == 0){
		////Initialize the region
		repellerRegion.num = attractorRegion.num = 0;

		if(type == 0) ////the limit cycle is an attractor, so the center of it is a repeller
		{
			////set the first triangle of the region which contains the center of the limit cycle
			repellerRegion.trianglelist[0] = singularities[singID].Triangle_ID;
			repellerRegion.num ++;
		}
		else          ////the limit cycle is an repeller, so the center of it is a attractor
		{
			////set the first triangle of the region which contains the center of the limit cycle
			attractorRegion.trianglelist[0] = singularities[singID].Triangle_ID;
			attractorRegion.num ++;
		}
	}
	////Begin region growing here
	UpdateBoundary(type);
	GrowLimitCycleRegion(type);
}



//////////////////////////////////////////////////////////
////Not just for saddle, other kind of singularities can also connect to limit cycle!!!
void UpdateListInSaddle(int saddleID, int limitcycle)
{
	////1.Extend current list
	////2.Add to the end of current list
	////3.Update the counter

	int i;
	//int *temp_list = limitcycles[limitcycle].connected_saddle;
	int *temp_list = singularities[saddleID].connected_limitcycles;

	if(temp_list == NULL)
	{
		singularities[saddleID].connected_limitcycles = (int*)malloc(sizeof(int));
		singularities[saddleID].connected_limitcycles[0] = limitcycle;
		singularities[saddleID].num_connected_limitcycles = 1;
	}

	else{
		singularities[saddleID].connected_limitcycles =
			(int*)malloc(sizeof(int)*(singularities[saddleID].num_connected_limitcycles+1));

		for(i = 0; i < singularities[saddleID].num_connected_limitcycles; i++)
		{
			singularities[saddleID].connected_limitcycles[i] = temp_list[i];
		}

		singularities[saddleID].connected_limitcycles[i] = limitcycle;

		singularities[saddleID].num_connected_limitcycles += 1;

		free(temp_list);
	}
}


////It seems that following codes do not check the repeating limit cycles in the list
////This is what we want for double connections (3/11/06)
void UpdateListInLimitCycle(int limitcycle, int saddleID)
{
	////Extend current list	and add to the end of current list

	int i;
	int *temp_list = limitcycles[limitcycle].connected_saddle;
	double *temp_length = limitcycles[limitcycle].flow_length_connectedsing;

	if(temp_list == NULL)
	{
		limitcycles[limitcycle].connected_saddle = (int*)malloc(sizeof(int));
		limitcycles[limitcycle].connected_saddle[0] = saddleID;
		limitcycles[limitcycle].num_connectedsaddles = 1;

		////Allocate space to save the flow length of the connection
		limitcycles[limitcycle].flow_length_connectedsing = 
			(double*)malloc(sizeof(double)*	limitcycles[limitcycle].num_connectedsaddles); //08/10/06
		limitcycles[limitcycle].flow_length_connectedsing[0] = sum_flow_length;  //08/10/06
	}

	else{
		if(IsRepeated(limitcycles[limitcycle].connected_saddle,
			saddleID, limitcycles[limitcycle].num_connectedsaddles))
			return;

		limitcycles[limitcycle].connected_saddle =
			(int*)malloc(sizeof(int)*(limitcycles[limitcycle].num_connectedsaddles+1));
		
		limitcycles[limitcycle].flow_length_connectedsing =
			(double*)malloc(sizeof(double)*(limitcycles[limitcycle].num_connectedsaddles+1)); //08/10/06

		for(i = 0; i < limitcycles[limitcycle].num_connectedsaddles; i++)
		{
			limitcycles[limitcycle].connected_saddle[i] = temp_list[i];
			limitcycles[limitcycle].flow_length_connectedsing[i] = temp_length[i]; //08/10/06
		}

		limitcycles[limitcycle].connected_saddle[i] = saddleID;
		limitcycles[limitcycle].flow_length_connectedsing[i] = sum_flow_length;  //08/10/06

		limitcycles[limitcycle].num_connectedsaddles += 1;

		free(temp_list);
	}

}


void UpdateCycleListInLimitCycle(int limitcycle1, int limitcycle2)
{
	////1. Extend current list
	////2. Add to the end of current list
	////3. Update the counter

	
	int i;
	int *temp_list = limitcycles[limitcycle1].connected_limitcycle;
	double *temp_length = limitcycles[limitcycle1].flow_length_connectedcycle;

	if(temp_list == NULL)
	{
		limitcycles[limitcycle1].connected_limitcycle = (int*)malloc(sizeof(int));
		limitcycles[limitcycle1].connected_limitcycle[0] = limitcycle2;
		limitcycles[limitcycle1].num_connectedcycles = 1;

		////08/10/06 store the flow length of this connection
		limitcycles[limitcycle1].flow_length_connectedcycle = (double*)malloc(sizeof(double));
		limitcycles[limitcycle1].flow_length_connectedcycle[0] = sum_flow_length;
	}

	else{
		//if the connection has already been built, return
		if(IsRepeated(limitcycles[limitcycle1].connected_limitcycle,
			limitcycle2, limitcycles[limitcycle1].num_connectedcycles))
			return;

		limitcycles[limitcycle1].connected_limitcycle =
			(int*)malloc(sizeof(int)*(limitcycles[limitcycle1].num_connectedcycles+1));
		
		////08/10/06 store the flow length of this connection
		limitcycles[limitcycle1].flow_length_connectedcycle = 
			(double*)malloc(sizeof(double)*(limitcycles[limitcycle1].num_connectedcycles+1));

		for(i = 0; i < limitcycles[limitcycle1].num_connectedcycles; i++)
		{
			limitcycles[limitcycle1].connected_limitcycle[i] = temp_list[i];

			limitcycles[limitcycle1].flow_length_connectedcycle[i] = temp_length[i]; //08/10/06
		}

		limitcycles[limitcycle1].connected_limitcycle[i] = limitcycle2;

		limitcycles[limitcycle1].flow_length_connectedcycle[i] = sum_flow_length;  ////08/10/06

		limitcycles[limitcycle1].num_connectedcycles += 1;

		free(temp_list);
	}
}
