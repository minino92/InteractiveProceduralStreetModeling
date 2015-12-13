//NewLimitCycleDetect.cpp  1/24/06


#include "stdafx.h"

#include "NewLimitCycleDetect.h"

#include "LimitCycleDetect.h"

#include "VFDataStructure.h"

#include "LocalTracing.h"

#include "numerical.h"

#include "topologyedit.h"

#include "Test_LimitcycleDetect.h"   ////when we change to new method, we may copy some useful routines to here 1/24/06

#include "LimitCycleCreator.h"

/*For testing*/
#include "GL/glut.h"

extern Polygon3D Object;
extern Singularities *singularities;           //being captured singularites' list
extern int cur_singularity_index;
extern LineSeg **trajectories;                 //trajectories' list
extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory
extern Separatrices *separatrices;             //array for group of separatrices


extern int MaxNumTrajectories;

extern int MaxNumLinesegsPerTraj;


////Actually, we need to find another efficient data structure to store the whole information of a limit cycle
////for future operations, such as reshape, moving
extern LimitCycle *limitcycles;
extern int cur_limitcycle_index;
extern int MaxNumLimitCycles;

extern int MaxTriangleInCellCycle;
extern int *TriangleList;
extern int *cellcycle;
extern int num_trianglelist;                ////number of triangles in the triangle list
extern int num_celltriangle;                ////number of triangles in the found cell cycle
extern Edge **Cycle_edge;                   ////mesh edges of the boundary of the cell cycle
extern Vertex *realexitvert;                ////currently founded real exit
extern int num_cycleedges;                  ////number of edges consist of the cell cycle 
extern int num_potentialverts;              ////number of vertices consist of the cell cycle
extern Edge *OneSharedEdge;                 ////one shared edge in the cell cycle for finding fixed point
extern double theFixPoint[2];               ////the fix point on the closed streamline
extern Edge *SecondSharedEdge;              ////another shared edge for finding next level beginning point

extern int Cur_CellCycleIndex;                    ////current cell cycle index


extern int *DesignCurveCellCycle;
extern int num_triangles_designcurve;


//// Use region growing to detect the limit cycle 1/11/2006
extern TriangularRegion repellerRegion;       ////region containing a repeller
extern TriangularRegion attractorRegion;      ////region containing an attractor
extern RegionBoundary repellerBoundary;
extern RegionBoundary attractorBoundary;
extern InnerVertices repellerInnerverts;
extern InnerVertices attractorInnerverts;

extern TriangularRegion intersectRegion;     ////The intersect region
extern RegionBoundary intersectBoundary;     ////The intersect boundary 
extern InnerVertices intersectInnerverts;    ////The inner vertices inside the intersect region

extern TriangularRegion Source_re;
extern TriangularRegion Sink_re;

//extern void IntersectRegion(TriangularRegion &source1, TriangularRegion &source2, TriangularRegion &dest);

extern Vertex * GetFurtherVer(Edge *sharededge, int singID);
extern void StoreCurStreamline(LineSeg *streamline, int num);
extern void IntersectRegion(TriangularRegion &source1, TriangularRegion &source2, TriangularRegion &dest);




/////Testing variables 01/5/06
extern double MarkNextBx[2];

////To reduce the calculation of intersections, we may need to store previous intersection 
////and compare it with current intersection, here is the temporary global variable for this purpose
////1/6/06
extern double first_intersection[2] ;
extern int first_intersect ;


/////////////////////////////////////////////////////////////
//Testing variables here
extern int test_aboundtriangle;
extern int second_begintriangle;
extern Edge *test_sel_edge ;


/////////Boundary list 4/11/06
BoundaryList aboundarylist;

/*------------------------------------------------------------------------*/

///////////////////////////////////////////////////////////////////////////
extern bool IsGoodSharedEdge(Edge *cur_e);

//////////////////////////////////////////////////////////////////////////



void CellCycleDetect_new(double &x, double &y, int &begin_triangle, int type, int &flag, int &stop_face)
{
	int i;
	int pre_face, cur_face;
	double globalp[2];
	//int position;

	////Initialization part
	globalp[0] = x;
	globalp[1] = y;
	pre_face = cur_face = begin_triangle;

	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->access_count = 0;
		Object.flist[i]->discard = 0;
	}

	////Perform local tracing and mark the triangles the curve passes through 
	Object.flist[cur_face]->access_count++;

	////Modified at 5/25/06
	for(i = 0; i < /*2*NUMTRACINGTRIANGLE*/3*Object.nfaces; i++)
	{
		if(cur_face == -1)
			return;

		//if(i == NUMTRACINGTRIANGLE-1)
		//{
		//	i = i;
		//}

		pre_face = cur_face ;
		cur_face = TraceInATriangle2(cur_face, globalp, type, flag);
 
		if(flag == 3 || flag == 4 /*|| pre_face == cur_face*/) 
		{
			flag = 2;            ////reach a singularity or the boundary
			return;
		}

		else if(Object.flist[cur_face]->access_count >= 5)
		{
			Object.flist[cur_face]->access_count-= 4;
			stop_face = cur_face;
			x = globalp[0];
			y = globalp[1];
			flag = 1;
			return;
		}

		else{
			Object.flist[cur_face]->access_count++;
		}
	}
}

int FindTheMostAccessTriangle()
{
	int i;
	int theone = -1;
	int cur_accessnum = 0;

	for(i = 0; i < Object.nfaces; i++)
	{
		if(Object.flist[i]->access_count > cur_accessnum && Object.flist[i]->discard != 1)
		{
			theone = i;
			cur_accessnum = Object.flist[i]->access_count;
		}
	}
    return theone;
}

////Test whether the triangle has an edge that tangent to current field
////Not finished 1/24/06
bool TangentialTest(int triangle)
{
	return false;
}


////Judge whether it is the edge we want to calculate the intersection
bool IsSameEdgeandTriangle(int sel_triangle, Edge *sel_edge, int pre_face)
{
	if((sel_edge->tris[0] == sel_triangle && sel_edge->tris[1] == pre_face)
		|| (sel_edge->tris[0] == pre_face && sel_edge->tris[1] == sel_triangle))
		return true;
	else
		return false;
}


////Judge whether the intersection is convergent or not
bool ConvergeIntersect(icVector2 IntersectList[], int size, int lhead, int lend)
{
	int i;
	//icVector2 diff;
	double len1, len2;

	if(lhead == lend || (lend-lhead+size)%size <= 2) 
		return true;

	len1 = length(IntersectList[(lhead+1)%size] - IntersectList[lhead]);

	for(i = 1; i < (lend-lhead+size)%size; i++)
	{
		len2 = length(IntersectList[(lhead+i+1)%size] - IntersectList[(lhead+i)%size]);
		if( len2 > len1)
		{
			return false;  ////it is not convergent
		}
		len1 = len2;
	}
	return true;
}


////Judge whether the intersection tends to converge to the same point
bool ReachThreshold(icVector2 IntersectList[], int size, int lhead, int lend)
{
	double len = length(IntersectList[(lend)%size] - IntersectList[(lend-1+size)%size]);
	if(len < 1e-8)
		return true;
	return false;
}



////Get the cell cycle after locating the closed streamline of the limit cycle
void GetCellCycleandStreamline(int triangle, double fix_pt[2], int type, \
							   int *cellcycle, int &num_cells, int &MaxTriangleInCellCycle)
{
	int i;
	int flag = 0;
	double globalp[2];
	int pre_face, cur_face;

	////Initialization part

	pre_face = cur_face = triangle;

	globalp[0] = fix_pt[0];   
	globalp[1] = fix_pt[1];

	num_cells = 0;

	////If current trajectories list is not enough to store all the trajectories
	////extend it!
	if(cur_traj_index >= MaxNumTrajectories - 1)
	{
		MaxNumTrajectories += 100;

		num_linesegs_curtraj = (int*)realloc(num_linesegs_curtraj, sizeof(int)*MaxNumTrajectories);
		trajectories = (LineSeg **)realloc(trajectories, sizeof(LineSeg *) * MaxNumTrajectories);

		////extend the old trajectories
		for(i = 0; i < MaxNumTrajectories-100; i++)
			trajectories[i] = (LineSeg *)realloc(trajectories[i], sizeof(LineSeg ) * MaxNumLinesegsPerTraj);

		////allocate the new trajectories
		for( ; i < MaxNumTrajectories; i++)
			trajectories[i] = (LineSeg *)malloc(sizeof(LineSeg) * MaxNumLinesegsPerTraj);

	}

	num_linesegs_curtraj[cur_traj_index] = 0;

	cellcycle[0] = cur_face;
	num_cells = 1;


	for(i = 0; i < NUMTRACINGTRIANGLE; i++)
	{

		////Set the contain_separatrix flag 1/3/06, may be useful in limit cycle destortion
		//Object.flist[cur_face]->contain_separatrix = 1;
		if(cur_face < 0)
			return;

		pre_face = cur_face;
		cur_face = TraceInATriangle(cur_face, globalp, type, flag); ////0 means always forward tracing here

		if(flag == 3 || flag == 4 || pre_face == cur_face ) 
		{
			return;
		}

		if(!IsRepeated(cellcycle, cur_face, num_cells))
		{
			////we may need to extend the cell cycle here firstly
            if(num_cells >= MaxTriangleInCellCycle - 1)
			{
				MaxTriangleInCellCycle += 50;
				cellcycle = (int*)realloc((int *)cellcycle, sizeof(int)*MaxTriangleInCellCycle);
			}

			cellcycle[num_cells] = cur_face;
			num_cells++;
		}

	}

	////saved for limit cycle detection
	//pre_cur_traj_index = cur_traj_index;


}




const int IntersectListsize = 3;

/*
Get the index of the sharing edge of two given triangle with respect to the first triangle

*/
int GetEdgeIndexofTriangle(int triangle1, int triangle2)
{
	Edge *cur_e;

	for(int i = 0; i < 3; i++)
	{
		cur_e = Object.flist[triangle1]->edges[i];

		if((cur_e->tris[0] == triangle1 && cur_e->tris[1] == triangle2)
			||(cur_e->tris[0] == triangle2 && cur_e->tris[1] == triangle1))
			return i;
	}
	return -1;
}


////
////we may record the cell cycle after locating the closed streamline
bool IsALimitCycle(int sel_triangle, Edge *sel_edge, double fix_point[2], int type)
{
	int i;
	int flag = -1;
	icVector2 IntersectList[IntersectListsize];
	int lhead, lend;
	double globalp[2];

	int pre_face, cur_face;
	
	////Initialization part
	fix_point[0] = globalp[0] = (Object.vlist[sel_edge->verts[0]]->x + Object.vlist[sel_edge->verts[1]]->x)/2;
	fix_point[1] = globalp[1] = (Object.vlist[sel_edge->verts[0]]->y + Object.vlist[sel_edge->verts[1]]->y)/2;

	////store the first intersection
	IntersectList[0].entry[0] = globalp[0];
	IntersectList[0].entry[1] = globalp[1];
	lhead = 0;
	lend = 1;

	pre_face = cur_face = sel_triangle;

	for(i = 0; i < Object.nfaces/*2*NUMTRACINGTRIANGLE*/; i++)
	{
		pre_face = cur_face;
		cur_face = TraceInATriangle2(cur_face, globalp, type, flag); ////0 means always forward tracing here

		if(flag == 3 || flag == 4 || pre_face == cur_face ) ////if it reach the boundary or singularity
			////or stop at somewhere, it can not be a limit cycle
		{
			return false;
		}

		else
		{
			if(cur_face == sel_triangle) ////we may need to calculate and compare the intersections here
			{
				if(IsSameEdgeandTriangle(sel_triangle, sel_edge, pre_face))
				{
					////first, we need to store the intersection
					if(lend == lhead)
						lhead = (lhead+1)%IntersectListsize;

					IntersectList[lend].entry[0] = globalp[0];
					IntersectList[lend].entry[1] = globalp[1];
					lend = (lend+1)%IntersectListsize;

					////If the intersection is not converge, return false
					if(!ConvergeIntersect(IntersectList, IntersectListsize, lhead, lend))
						return false;
					
					if(ReachThreshold(IntersectList, IntersectListsize, lhead, lend))
					{
						fix_point[0] = globalp[0];
						fix_point[1] = globalp[1];

						GetCellCycleandStreamline(sel_triangle, fix_point, type, \
							cellcycle, num_celltriangle, MaxTriangleInCellCycle);
						return true;
					}
				}

				else //update with the new edge and begin searching from the middle of the new edge
				{
					sel_edge = Object.flist[sel_triangle]->edges[GetEdgeIndexofTriangle(sel_triangle, pre_face)];

					//reset the fixed point as the middle point of the new edge
					fix_point[0] = globalp[0] = (Object.vlist[sel_edge->verts[0]]->x + Object.vlist[sel_edge->verts[1]]->x)/2;
					fix_point[1] = globalp[1] = (Object.vlist[sel_edge->verts[0]]->y + Object.vlist[sel_edge->verts[1]]->y)/2;
				}
			}
		}

	}

	GetCellCycleandStreamline(sel_triangle, fix_point, type, \
		cellcycle, num_celltriangle, MaxTriangleInCellCycle);
	return true;

}


bool OntheEdge(double x1, double y1, double x2, double y2, double pt[2])
{
	double A, B, C;
	A = y1 - y2;
	B = x2 - x1;
	C = (x1*y2 - x2*y1);

	double pending = A*pt[0] + B*pt[1] + C;

	if(fabs(pending) < 1e-8) ////passing the vertex
	{
		return true;
	}
}


///////////////////////////////////////////////////////
////Use the point and the triangle that contains the point to get the edge that the point falls on
Edge* GetTheEdge(int triangle, double globalp[2], Edge *sel_edge, int type, int &passverflag)
{
	Vertex *v1, *v2;
	Edge *edge;
    int i;

	int cur_face, pre_face;
	int cur_access_count;
	int flag;
	double sec_globalp[2] = {globalp[0], globalp[1]};

	////perform passing vertex testing
	icVector2 temp_d;
	for(i = 0; i < 3; i++)
	{
		v1 = Object.vlist[Object.flist[triangle]->verts[i]];
		temp_d.entry[0] = v1->x - globalp[0];
		temp_d.entry[1] = v1->y - globalp[1];

		if(length(temp_d) < 1e-8)  ////it means that the given point happens to be a vertex
		{
			pre_face = cur_face = triangle;
			cur_access_count = Object.flist[triangle]->access_count;
			Object.flist[triangle]->access_count++;

			////we need to perform one more round local tracing
			for(int j = 0; j < NUMTRACINGTRIANGLE; j++)
			{
				pre_face = cur_face ;
				cur_face = TraceInATriangle2(cur_face, sec_globalp, type, flag);
		 
				if(flag == 3 || flag == 4 || pre_face == cur_face) 
				{
					////reach a singularity or the boundary
					break;
				}

				else if(Object.flist[cur_face]->access_count >= cur_access_count+2)
				{
					Object.flist[cur_face]->access_count-= 1;
					
					temp_d.entry[0] = sec_globalp[0] - globalp[0];
					temp_d.entry[1] = sec_globalp[1] - globalp[1];
					if(length(temp_d) < 1e-8)  ////it always passes this vertex, then it can be the fixed point
					{
						passverflag = 1;
						return NULL;
					}
					
					break;  ////just perform one round tracing here, not any more!!
				}

				else{
					Object.flist[cur_face]->access_count++;
				}
			}
		}
	}

	for(i = 0; i < 3; i++)
	{
		edge = Object.flist[triangle]->edges[i];

		v1 = Object.vlist[edge->verts[0]];
		v2 = Object.vlist[edge->verts[1]];

		if(OntheEdge(v1->x, v1->y, v2->x, v2->y, globalp))
		{
			sel_edge = edge;
			return edge;
		}
	}
}



/*----------------------------------------------------------------------------*/

bool FindCellCycle_2(int *trilist, int originNum, int *acycle, int &CellNum, int cur_triangleID, 
					 int &position, int &flag)
{
	int j;

	if(!TriangleSearch(trilist, originNum, cur_triangleID, position))
	{
		return false;
	}

	else{  ////at this moment, we do not consider the number of cell cycle  <= 2
		if(position == originNum - 1)
		{
			flag = 1;
			return false;
		}
	}

	////we find a cell cycle
	////copy the triangles inside the cell cycle to build a new list
	CellNum = 0;
    for( j = 0; j < originNum-position; j++)
	{
		acycle[j] = trilist[position + j];
		CellNum ++;
	}

	return true;
}

void GetPartialArray(int *a, int &originNum, int beginning_pos, int ending_pos)
{
	int i;

	for(i = 0; i < ending_pos - beginning_pos; i++)
	{
		a[i] = a[i+beginning_pos];
	}

	originNum = ending_pos - beginning_pos;
}


bool IsPreviousLimitCycleForSaddle(int &limit_index)
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
				continue;   ////continue to test next limit cycle

			limit_index = i;
			return true;
		}

	}

	return false;
}



////Detect the four separatrices of the saddle to find the cell cycles, and find the connections
////between the saddle and the limit cycles
////This routine will be really really important one! 1/31/06

void SearchFromSaddle(int saddle)
{
	int i, j, sep;

	int MaxNumTriangles = 600;

	int *trianglelist = (int*)malloc(sizeof(int)*MaxNumTriangles);
	int originNum = 0;

	int flag = -1;
	int position = -1;

	////for each separatrix, we perform following operations
	for(i = 0; i < 4; i++)
	{
		originNum = 0;
		
		if(singularities[saddle].separtices < 0)
			continue;

		switch(i)
		{
		case 0:
			sep = separatrices[singularities[saddle].separtices].sep1;
			break;
		case 1:
			sep = separatrices[singularities[saddle].separtices].sep2;
			break;
		case 2:
			sep = separatrices[singularities[saddle].separtices].sep3;
			break;
		case 3:
			sep = separatrices[singularities[saddle].separtices].sep4;
			break;
		}

		if(sep < 0)
			continue;

		////search for a close cycle
		for(j = 0; j < num_linesegs_curtraj[sep]; j++)
		{
			flag = -1;
			if(FindCellCycle_2(trianglelist, originNum, cellcycle, num_celltriangle, trajectories[sep][j].Triangle_ID,
				position, flag))
			{
				////We find a cell cycle, test whether it is the limit cycle we obtained before
				int pre_limit_index = -1;
				if(IsPreviousLimitCycleForSaddle(pre_limit_index))
				{
					if(pre_limit_index == -1)
						continue;
					
					UpdateListInSingularity(saddle, pre_limit_index);
					UpdateListInLimitCycle(pre_limit_index, saddle);
					break;
				}

				GetPartialArray(trianglelist, originNum, position+1, originNum-1);
			}

			else
			{
				if(originNum >= MaxNumTriangles - 1)
				{
					MaxNumTriangles += 50;
					trianglelist = (int*)realloc(trianglelist, sizeof(int)*MaxNumTriangles);
				}

				if(flag != 1)
				{
					trianglelist[originNum] = trajectories[sep][j].Triangle_ID;
					originNum ++;
				}
			}
		}
	}

	free(trianglelist);
}




//////////////////////////////////////////////////////

void DetectALimitCycle_new(double begin_p[2], int begin_triangle, int type, int inner, int singID, int &flag)
{
	//int flag = -1;
	//int fix_triangleID;

	int chosen_triangle;
	Edge *chosen_edge = NULL;
	double fix_pt[2];

	////if the input singularity is a saddle, we need to search all its separatrices, so we have to use 
	////the old method to find their connections

	CellCycleDetect_new(begin_p[0], begin_p[1], begin_triangle, type, flag, chosen_triangle); 
	
	if(flag == 2) //reach a singularity or the boundary of the mesh
	{
		flag = 2;
		return;
	}

	if(flag == 1) //we may find a cell cycle, but need to perform further testing
	{
		////We may use the last intersected edge in previous step as the chosen_edge if it satisfies the
		////condition
		int throughver = 0;
		chosen_edge = GetTheEdge(chosen_triangle, begin_p, chosen_edge, type, throughver);
        OneSharedEdge = chosen_edge;

		if(throughver == 1)
		{
			GetCellCycleandStreamline(chosen_triangle, begin_p, type,\
				cellcycle, num_celltriangle, MaxTriangleInCellCycle);
			goto L2;
		}

		/*-----------------------------------------------*/
		//////Testing code here only
		second_begintriangle = chosen_triangle;
		test_sel_edge = chosen_edge;


		if(singularities[singID].type == SINK)
		{
			singID = singID;

		}
		/*-----------------------------------------------*/

		////perform local tracing from the middle point of the edge and calculate the intersection
		if(IsALimitCycle(chosen_triangle, chosen_edge, fix_pt, type))
		{

		    int pre_limit_index = -1;

			//if it is a previously detected limit cycle
L2:			if(IsPreviousLimitCycleForSaddle(pre_limit_index))
			//if(IsPreviousLimitCycle(pre_limit_index)) //it relies on the cell cycle
			{
				if(inner == 0) //the center singularity can only connect with the most inner limit cycle
				{
					UpdateListInSingularity(singID, pre_limit_index);
					UpdateListInLimitCycle(pre_limit_index, singID);
				}
				else //if it is an outer limit cycle, think about that here !!
				{
					UpdateCycleListInLimitCycle(cur_limitcycle_index-1, pre_limit_index);
					UpdateCycleListInLimitCycle(pre_limit_index, cur_limitcycle_index-1);
				}

				return;
			}

			/*---------------------------------------------------------------*/
			
			/*---------------------------------------------------------------*/
			////Store the information of the new limit cycle
					
			StoreCurrentCellCycle(cellcycle, num_celltriangle);  

			limitcycles[cur_limitcycle_index].singularID = singID;  ////store the center singularity of current limit cycle
			
			limitcycles[cur_limitcycle_index].node_index = -1; ////10/16/05

			StoreCurStreamline(trajectories[cur_traj_index], num_linesegs_curtraj[cur_traj_index]);     

			////store the type of the limit cycle
			if(type == 0)
				limitcycles[cur_limitcycle_index].type = 1;
			else
				limitcycles[cur_limitcycle_index].type = 0;

			////Initialize the connected list for the graph
			limitcycles[cur_limitcycle_index].connected_limitcycle = NULL;
			limitcycles[cur_limitcycle_index].num_connectedcycles = 0;
			limitcycles[cur_limitcycle_index].connected_saddle = NULL;
			limitcycles[cur_limitcycle_index].num_connectedsaddles = 0;

			/*---------------------------------------------------------------*/

			//Build the conley connection for the limit cycle
			
			if(inner == 0) //if it is the most inner limit cycle
			{
				UpdateListInSingularity(singID, cur_limitcycle_index);
				UpdateListInLimitCycle(cur_limitcycle_index, singID);
			}

			else //if it is an outer limit cycle
			{
				UpdateCycleListInLimitCycle(cur_limitcycle_index-1, cur_limitcycle_index);
				UpdateCycleListInLimitCycle(cur_limitcycle_index, cur_limitcycle_index-1);
			}
			/*---------------------------------------------------------------*/

			////Build the graphical handler for the limit cycle for user interaction
			
			//BuildHandlerforLimitCycle(cur_limitcycle_index); //it depends on the cell cycle also
			BuildHandlerforLimitCycle(cur_limitcycle_index, fix_pt, chosen_edge->verts[0]);
			
			cur_traj_index ++;
			cur_limitcycle_index++;
		}
	}
}

////Find the connected edge that shares the given vertex 2/4/06
Edge *FindConnectedEdge(Edge **edgelist, int num_edges, int vert, int &zeroorone)
{
	Edge *theedge = NULL;

	for(int i = 0; i < num_edges; i++)
	{
		if(edgelist[i]->visited == 1)
			continue;

		if( edgelist[i]->verts[0] == vert || edgelist[i]->verts[1] == vert)
		{
			theedge = edgelist[i];
			return theedge;
		}
	}
 
	return theedge;
}


////sort the boundary edge list 2/4/06
void SortEdgelist(int type, int &flag)
{
	Edge **temp_edgelist, **secondtemp_edgelist;
	Edge *findtheedge;
	int temp_edgenums;
	int cur_edgenums;
	int zeroorone = 0;

	int end_vert, cur_vert;
	int i;

	////Initialization part

	if(type == 0)
	{
		temp_edgelist = repellerBoundary.edgelist;
		temp_edgenums = repellerBoundary.num;
	}
	else if(type == 1)
	{
		temp_edgelist = attractorBoundary.edgelist;
		temp_edgenums = attractorBoundary.num;
	}
	else
	{
		temp_edgelist = intersectBoundary.edgelist;
		temp_edgenums = intersectBoundary.num;
	}

	////set the visited flag here
	for(i = 0; i < temp_edgenums; i++)
	{
		temp_edgelist[i]->visited = 0;
	}

	secondtemp_edgelist = (Edge **)malloc(sizeof(Edge*)*temp_edgenums);

	secondtemp_edgelist[0] = temp_edgelist[0];
	end_vert = temp_edgelist[0]->verts[0];
	cur_vert = temp_edgelist[0]->verts[1];
	cur_edgenums = 1;
	temp_edgelist[0]->visited = 1;

	flag = 0;
	////connect the edges in order
	while(cur_edgenums != temp_edgenums )
	{
		////try to deal with multiple boundaries here, may be wrong
		if(cur_vert == end_vert && cur_edgenums != temp_edgenums)
		{
			for(i = 0; i < temp_edgenums; i++)
			{
				if(temp_edgelist[i]->visited == 0)  ////find any unvisited edge
				{
					secondtemp_edgelist[cur_edgenums] = temp_edgelist[i];
					cur_edgenums ++;
					temp_edgelist[i]->visited = 1;
					end_vert = temp_edgelist[i]->verts[0];
					cur_vert = temp_edgelist[i]->verts[1];
					flag ++;
					break;
				}
			}
		}


		findtheedge = FindConnectedEdge(temp_edgelist, temp_edgenums, cur_vert, zeroorone);

		if(findtheedge == NULL)
		{
			////something is wrong in the boundary edge list, they should all connected
			//MessageBox(NULL, "Wrong edge list!", "Error", MB_OK);

			///This always means that the region is separated, which may not be valid any more!!!! 2/4/06
			return;
		}

		secondtemp_edgelist[cur_edgenums] = findtheedge;
		cur_edgenums ++;
		findtheedge ->visited = 1;

		if(findtheedge->verts[0] != cur_vert)
			cur_vert = findtheedge->verts[0];
		else
			cur_vert = findtheedge->verts[1];
	}

	////set back to the orginal edge list
	if(type == 0)
	{
		for(i = 0; i < cur_edgenums; i++)
		{
			repellerBoundary.edgelist[i] = secondtemp_edgelist[i];
		}
		repellerBoundary.num = cur_edgenums;
	}

	else if(type == 1)
	{
		for(i = 0; i < cur_edgenums; i++)
		{
			attractorBoundary.edgelist[i] = secondtemp_edgelist[i];
		}
		attractorBoundary.num = cur_edgenums;
	}

	else
	{
		for(i = 0; i < cur_edgenums; i++)
		{
			intersectBoundary.edgelist[i] = secondtemp_edgelist[i];
		}
		intersectBoundary.num = cur_edgenums;
	}

	free(secondtemp_edgelist);
}


////The routine for returning an outer boundary triangle 2/4/06
int OuterBoundaryTriangle(int type, int &anintriangle)
{
	int thetriangle = -1;
	int i;
	int EndVertID;
	Vertex *cur_vert;
	Edge *cur_edge, *other_e;
	Edge **temp_boundaryedgelist;
	int temp_boundaryedgenums;

	Edge **temp_edgelist;

	int zeroorone;

	int temp_num_edges, other_num_edges;
	
	////First build the boundary edge list for the region
	//UpdateBoundary(type);

	if(type == 0)
	{
		temp_boundaryedgelist = repellerBoundary.edgelist;
		temp_boundaryedgenums = repellerBoundary.num;
	}

	else
	{
		temp_boundaryedgelist = attractorBoundary.edgelist;
		temp_boundaryedgenums = attractorBoundary.num;
	}

	temp_edgelist = (Edge **) malloc(sizeof(Edge *) * temp_boundaryedgenums);


	////Initial the flag of all the edges on current boundaries
	for(i = 0; i < temp_boundaryedgenums; i++)
	{
		temp_boundaryedgelist[i]->visited = 0;
	}

	temp_num_edges = 0;
	temp_edgelist[0] = temp_boundaryedgelist[0];
	temp_num_edges++;
	EndVertID = temp_boundaryedgelist[0]->verts[0];
	temp_edgelist[0]->visited = 1;

	cur_vert = Object.vlist[temp_boundaryedgelist[0]->verts[1]];

	while(cur_vert->VertID != EndVertID)   ////Not form a closed edges list
	{
		cur_edge = FindConnectedEdge(temp_boundaryedgelist, temp_boundaryedgenums, cur_vert->VertID, zeroorone);

		temp_edgelist[temp_num_edges] = cur_edge;
		temp_num_edges ++;
		cur_edge->visited = 1;
		
		////Get next testing vertex
		if(cur_edge->verts[0] != cur_vert->VertID)
			cur_vert = Object.vlist[cur_edge->verts[0]];
		else
			cur_vert = Object.vlist[cur_edge->verts[1]];
	}

    other_num_edges = temp_boundaryedgenums - temp_num_edges;

	if(other_num_edges > temp_num_edges)
	{
		////get an edge from the outer boundary
		for(i = 0; i < temp_boundaryedgenums; i++)
		{
			cur_edge = temp_boundaryedgelist[i];

			if(cur_edge->visited == 0)
				break;
		}

		////Get an edge on the inner boundary
		other_e = temp_edgelist[0];
	}

	else
	{
		////Get an edge on outer boundary
		cur_edge = temp_edgelist[0];
		
		////get an edge on the inner boundary
		for(i = 0; i < temp_boundaryedgenums; i++)
		{
			other_e = temp_boundaryedgelist[i];

			if(other_e->visited == 0)
				break;
		}
	}

	////find the triangle that inside the region 
	if(type == 0)
	{
		if(Object.flist[cur_edge->tris[0]]->repell_inregion == 1)
		{
			thetriangle = cur_edge->tris[0];
		}
		else
		{
			thetriangle = cur_edge->tris[1];
		}

		////Get a triangle at the inner boundary 3/22/06
		if(Object.flist[other_e->tris[0]]->repell_inregion == 1)
			anintriangle = other_e->tris[0];
		else
			anintriangle = other_e->tris[1];
	}

	else
	{
		if(Object.flist[cur_edge->tris[0]]->attract_inregion == 1)
		{
			thetriangle = cur_edge->tris[0];
		}
		else
		{
			thetriangle = cur_edge->tris[1];
		}
		
		////Get a triangle at the inner boundary 3/22/06
		if(Object.flist[other_e->tris[0]]->attract_inregion == 1)
			anintriangle = other_e->tris[0];
		else
			anintriangle = other_e->tris[1];
	}

	free(temp_edgelist);
	return thetriangle;
 
}



/*-----------------------------------------------------------------------------*/
////3/22/06
/*
Get an inner edge after getting the cell cycle
Suppose we have already get a cell cycle stored in "cellcycle"
*/
Edge* GetAnEdgeForFixedPoint(int *cellcycle, int num_celltriangle)
{
	int i, j;
	Edge *cur_e = NULL;
	Face *face;
	icVector2 edge_vec;
	icVector2 vec1, vec2;
	Vertex *v1, *v2;

	for(i = 0; i < num_celltriangle; i++)
	{
		face = Object.flist[cellcycle[i]];

		for(j = 0; j < 3; j++)
		{
			cur_e = face->edges[j];
			v1 = Object.vlist[cur_e->verts[0]];
			v2 = Object.vlist[cur_e->verts[1]];
			vec1 = v1->vec;
			vec2 = v2->vec;

			edge_vec.entry[0] = v1->x - v2->x;
			edge_vec.entry[1] = v1->y - v2->y;

			normalize(edge_vec);
			normalize(vec1);
			normalize(vec2);
	
			if(dot(edge_vec, vec1) <= 0.5 && dot(edge_vec, vec1) <= 0.5)
				return cur_e;
		}
	}
    
	return cur_e; //just pick the last edge
}

/*
Search the fixed point using the old method
*/
void GettheFixedPoint(Edge *theedge, int type, double fixedpt[2])
{
	////first, we need to choose an inner edge in the cell cycle (one of the sharing edge of the triangle strip)
	////We have already stored it in the OneSharedEdge
	double v1[2], v2[2];
	Vertex *vert1, *vert0;

	////initialize
	vert0 = Object.vlist[theedge->verts[0]];
	vert1 = Object.vlist[theedge->verts[1]];

	v1[0] = vert0->x;
	v1[1] = vert0->y;

	v2[0] = vert1->x;
	v2[1] = vert1->y;

	////begin binary search
	BinaryFindtheFixedPoint(v1, v2, type);

	fixedpt[0] = theFixPoint[0];
	fixedpt[1] = theFixPoint[1];
}

/*
Get a limit cycle
*/
void GetALimitCycle(int out_t, double out_gp[2],
					int in_t, double in_gp[2],
					int *cellcycle, int &num_cells, 
					int singID, int type, int inner)
{
	double fix_pt[2] = {0.};

	ParallelCellCycleLocate(out_t, out_gp,
		in_t, in_gp, cellcycle, num_cells, type);

	OneSharedEdge = GetAnEdgeForFixedPoint(cellcycle,num_cells);
	GettheFixedPoint(OneSharedEdge, type, fix_pt);

	if(OneSharedEdge->tris[0] < 0 || OneSharedEdge->tris[1] < 0)
		return;

	GetCellCycleandStreamline(OneSharedEdge->tris[0], fix_pt, type, \
							  cellcycle, num_cells, MaxTriangleInCellCycle);

	int pre_limit_index = -1;

	//if it is a previously detected limit cycle
	if(IsPreviousLimitCycleForSaddle(pre_limit_index))
	//if(IsPreviousLimitCycle(pre_limit_index)) //it relies on the cell cycle
	{
		if(inner == 0) //the center singularity can only connect with the most inner limit cycle
		{
			UpdateListInSingularity(singID, pre_limit_index);
			UpdateListInLimitCycle(pre_limit_index, singID);
		}
		else //if it is an outer limit cycle, think about that here !!
		{
			UpdateCycleListInLimitCycle(cur_limitcycle_index-1, pre_limit_index);
			UpdateCycleListInLimitCycle(pre_limit_index, cur_limitcycle_index-1);
		}

		return;
	}
	
	/*---------------------------------------------------------------*/
	////Store the information of the new limit cycle
			
	StoreCurrentCellCycle(cellcycle, num_celltriangle);  

	limitcycles[cur_limitcycle_index].singularID = singID;  ////store the center singularity of current limit cycle
	
	limitcycles[cur_limitcycle_index].node_index = -1; ////10/16/05

	StoreCurStreamline(trajectories[cur_traj_index], num_linesegs_curtraj[cur_traj_index]);     

	////store the type of the limit cycle
	if(type == 0)
		limitcycles[cur_limitcycle_index].type = 1;
	else
		limitcycles[cur_limitcycle_index].type = 0;

	////Initialize the connected list for the graph
	limitcycles[cur_limitcycle_index].connected_limitcycle = NULL;
	limitcycles[cur_limitcycle_index].num_connectedcycles = 0;
	limitcycles[cur_limitcycle_index].connected_saddle = NULL;
	limitcycles[cur_limitcycle_index].num_connectedsaddles = 0;

	/*---------------------------------------------------------------*/

	//Build the conley connection for the limit cycle
	
	if(inner == 0) //if it is the most inner limit cycle
	{
		UpdateListInSingularity(singID, cur_limitcycle_index);
		UpdateListInLimitCycle(cur_limitcycle_index, singID);
	}

	else //if it is an outer limit cycle
	{
		UpdateCycleListInLimitCycle(cur_limitcycle_index-1, cur_limitcycle_index);
		UpdateCycleListInLimitCycle(cur_limitcycle_index, cur_limitcycle_index-1);
	}
	/*---------------------------------------------------------------*/

	////Build the graphical handler for the limit cycle for user interaction
	
	BuildHandlerforLimitCycle(cur_limitcycle_index); //it depends on the cell cycle also

	cur_limitcycle_index++;

	/*---------------------------------------------------------------*/
}


/* ------------------4/6/06---------------------*/
void GetInitTriangleStrip(int type)
{
	int i;
	Edge *cur_e;

	if(type == 0)
	{
		UpdateBoundary(0);
		repellerRegion.num = 0;
		for(i = 0; i < repellerBoundary.num; i++)
		{
			cur_e = repellerBoundary.edgelist[i];

			if(IsRepeated(Source_re.trianglelist, cur_e->tris[0], Source_re.num))
			{
				repellerRegion.trianglelist[repellerRegion.num] = cur_e->tris[0];
				repellerRegion.num++;
			}
			else
			{
				repellerRegion.trianglelist[repellerRegion.num] = cur_e->tris[1];
				repellerRegion.num++;
			}
		}
	}

	else
	{
		UpdateBoundary(1);
		attractorRegion.num = 0;
		for(i = 0; i < attractorBoundary.num; i++)
		{
			cur_e = attractorBoundary.edgelist[i];

			if(IsRepeated(Sink_re.trianglelist, cur_e->tris[0], Sink_re.num))
			{
				attractorRegion.trianglelist[attractorRegion.num] = cur_e->tris[0];
				attractorRegion.num++;
			}
			else
			{
				attractorRegion.trianglelist[attractorRegion.num] = cur_e->tris[1];
				attractorRegion.num++;
			}
		}
	}
}


/*-----------------------------------------------------------------------------*/


void GetACellCycle_new2(int singID, int type)
{
	//int i;
	Edge **temp_edgelist;
	int num_boundedges = 0;
    Edge *oneboundedge;
	int aboundtriangle = singularities[singID].Triangle_ID;
	int multiboundaryornot = 0;
	
	//a triangle at the inner boundary if exists 3/22/06
	//this will fail if the limit cycle cross the center triangle
	int aninerboundtriangle = Object.clist[3*singularities[singID].Triangle_ID]->ot; 
	
	double begin_p_out[2]/*, begin_p_in[2]*/;
	int inner = 0;

	int test_count = 0;

	int flag = -1;

	repellerRegion.num = 0;
	attractorRegion.num = 0;

	////first, grow a region from the triangle containing the center singularity
LL:	if(type == 0)  
	{
		//the type of the singularity at the center is a repeller, hence the limit cycle 
		//should be an attractor
		if(inner == 0)
		{
			repellerRegion.trianglelist[repellerRegion.num] = aboundtriangle;
			repellerRegion.num ++;
		}

		UpdateBoundary(0);
		GetRegionNormals(0);
		Cancel_Growing(0, -1);

		temp_edgelist = repellerBoundary.edgelist;
		num_boundedges = repellerBoundary.num;
		//oneboundedge = repellerBoundary.edgelist[0];
		oneboundedge = repellerBoundary.edgelist[num_boundedges-1];
		aboundtriangle = repellerRegion.trianglelist[repellerRegion.num-1]; //this is a boundary triangle

		////Copy the region
		CopyRegion(repellerRegion.trianglelist, Source_re.trianglelist, repellerRegion.num);
		Source_re.num = repellerRegion.num;
	}

	else
	{
		if(inner == 0)
		{
			attractorRegion.trianglelist[attractorRegion.num] = aboundtriangle;
			attractorRegion.num ++;
		}

		UpdateBoundary(1);
		GetRegionNormals(1);
		Cancel_Growing(1, -1);

		temp_edgelist = attractorBoundary.edgelist;
		num_boundedges = attractorBoundary.num;
		//oneboundedge = attractorBoundary.edgelist[0];
		oneboundedge = attractorBoundary.edgelist[num_boundedges-1];
		aboundtriangle = attractorRegion.trianglelist[attractorRegion.num-1]; //this is a boundary triangle
		
		////Copy the region
		CopyRegion(attractorRegion.trianglelist, Sink_re.trianglelist, attractorRegion.num);
		Sink_re.num = attractorRegion.num;
	}

	/*-------------------------------------------------------------------------*/
	////grow the backward region 3/22/06

	GetInitTriangleStrip(type);

	if(type == 0)
	{
		//repellerRegion.num = 0;
		//repellerRegion.trianglelist[0] = aboundtriangle;
		//repellerRegion.num = 1;
		
		UpdateBoundary(0);
		GetRegionNormals(0);
		Cancel_Growing(0, -1);

		////calculate the intersection region
		IntersectRegion(repellerRegion, Source_re, intersectRegion);
	}

	else
	{
		//attractorRegion.num = 0;
		//attractorRegion.trianglelist[0] = aboundtriangle;
		//attractorRegion.num = 1;

		UpdateBoundary(1);
		GetRegionNormals(1);
		Cancel_Growing(1, -1);
		
		////calculate the intersection region
		IntersectRegion(attractorRegion, Sink_re, intersectRegion);
	}
	////After calculate the intersection of two regions, we get a ring-like region
	////That is a little bit larger than the cell cycle containing the limit cycle
	/*---------------------------------------------------------------------------*/

	////pick one triangle at the boundary of the region obtained above

	////An important problem! for embedded limit cycle detection, we can not use random triangle 
	////on the boundary, here we need to use the triangle on the outer boundary!!!

	////Note that there is no inner or outer boundary on surfaces, especially on torus 3/22/06
	
	//UpdateBoundary(type);
	//SortEdgelist(type, multiboundaryornot);

	//if(multiboundaryornot > 0)
	//	aboundtriangle = OuterBoundaryTriangle(type, aninerboundtriangle);

	//else{
	//	if(IntheRegion(oneboundedge->tris[0], type))
	//		aboundtriangle = oneboundedge->tris[0];
	//	else
	//		aboundtriangle = oneboundedge->tris[1];
	//}

	////Here just a testing code
	//test_aboundtriangle = aboundtriangle;
	//second_begintriangle = aninerboundtriangle;

	if(intersectRegion.num == 0)
		return;

	aboundtriangle = intersectRegion.trianglelist[rand()%intersectRegion.num];

	////Locate the real cell cycle and closed orbit
	GetCenterofATriangle(aboundtriangle, begin_p_out);
	//GetCenterofATriangle(aninerboundtriangle, begin_p_in);

	DetectALimitCycle_new(begin_p_out, aboundtriangle, type, inner, singID, flag);

	////Testing on 5/24/06
		//if(FindALimitCycle(intersectRegion.trianglelist, intersectRegion.num,
		//	begin_p_out, aboundtriangle, inner, singID, type))
		//{
		//	//save the limit cycle information and build the connection relation
		//}

		//else //not detect a limit cycle, abort
		//{
		//	return;
		//}

	///////////////////////////////////////////////////////////////////////////////////////


	//////For embedded limit cycles, we should begin from the previous selected triangle and perform 
	//////backward region growing!!(which means repeller growing will become attractor growing)

	//////first, we need to reverse the type
	//////if the selected boundary triangle is not at the boundary of the whole field
	if(IsAtFieldBoundary(aboundtriangle) || test_count >= 4 || flag == 2)
		return;

	////Interlace copy here !!! 1/11/06
	if(type == 0)
	{
		CopyRegion(repellerRegion.trianglelist, attractorRegion.trianglelist, repellerRegion.num);
		attractorRegion.num = repellerRegion.num;
	}

	else
	{
		CopyRegion(attractorRegion.trianglelist, repellerRegion.trianglelist, attractorRegion.num);
		repellerRegion.num = attractorRegion.num;
	}

	type = 1 - type;

	test_count ++;

	inner ++;

	goto LL;
}


//#define TESTLIMITDETECT


/*
Main routine for limit cycle detection
*/
void DetectLimitCycle_new2()
{
	//declare local variables
    int i, Saddleornot, type;

	int *saddlelist = NULL;
	int Num_saddles = 0;

	//Initialization part
    cur_limitcycle_index = 0;

	/* -- test -- 06/18/06 -- */
//#ifdef TESTLIMITDETECT
//	int onesing = 0;
//#endif


	for(i = 0; i < cur_singularity_index; i++)
	{
		Saddleornot = 0;

		InitLimitCycleDetection(); //re-initialized the limit cycle detection each time starts from a new singularity

		if(singularities[i].type == CWCENTER || singularities[i].type == CCWCENTER)
			continue;

		if(singularities[i].type == SOURCE || singularities[i].type == RFOCUS)//it is a repeller
		{
			type = 0;
		}

		else if(singularities[i].type == SINK || singularities[i].type == AFOCUS) //it is an attractor
		{
			type = 1;
		}

		else   //it is a saddle
		{
			type = 0;
			Saddleornot = 1;
			saddlelist = Extend_link(saddlelist, Num_saddles);
			saddlelist[Num_saddles] = i;
			Num_saddles++;
			continue;
		}

		//GetACellCycle_new2(i, type);


		////GetACellCycle_3(i, type);  //method 2 4/11/06

		GetCellCycles(i, type);	 //method 3 5/23/06 use intersection not so good!
		
		/* -- test -- 06/18/06 -- */
//#ifdef TESTLIMITDETECT
//		onesing = 1; //we just test a detection from one single singularity
//#endif
//		/* -- test -- 06/18/06 -- */
//#ifdef TESTLIMITDETECT
//		if (onesing == 1)//we just test a detection from one single singularity
//			return;
//#endif
	}

	for(i = 0; i < Num_saddles; i++)
	{
		if(singularities[saddlelist[i]].type == SADDLE)
		SearchFromSaddle(saddlelist[i]);
	}

	if(saddlelist != NULL)
		free(saddlelist);
}









/******************************************************************************/
////Using two local tracing from the outer and inner boundary to 
////local a cell cycle
bool ParallelCellCycleLocate(int out_t, double out_gp[2],
							 int in_t, double in_gp[2],
							 int *cellcycle, int &num_cells,
							 int type)
{
	int *local_cell_out = (int*)malloc(sizeof(int)*1000);
	int num_localcellout = 0;
	int *local_cell_in = (int*)malloc(sizeof(int)*1000);
	int num_localcellin = 0;

	int flag = -1;

	int origin_out_t = out_t;
	int origin_in_t = in_t;

	int test_rount = 0;
	//int i;

	First_CellCycleDetect_2(out_gp[0], out_gp[1], out_t,
			cellcycle, num_cells, type, flag);


	free(local_cell_out);
	free(local_cell_in);
	return false;
}

bool IsSameIntArray(int *a1, int num_a1, int *a2, int num_a2)
{
	int i, j;
	int position = 0;
	
	//for(i = 0; i < num_a1; i++)
	//{
	//	////for each element in a1, we find the corresponding element in a2
	//	yes = 0;
	//	for(j = 0; j < num_a2; j++)
	//	{
	//		if(a1[i] == a2[j])
	//		{
	//			yes = 1;
	//			break;
	//		}
	//	}

	//	if(yes == 0)
	//		return false;
	//}

	//for(i = 0; i < num_a2; i++)
	//{
	//	////for each element in a1, we find the corresponding element in a2
	//	yes = 0;
	//	for(j = 0; j < num_a1; j++)
	//	{
	//		if(a2[i] == a1[j])
	//		{
	//			yes = 1;
	//			break;
	//		}
	//	}

	//	if(yes == 0)
	//		return false;
	//}

	if(num_a1 != num_a2) return false;

	for(i = 0; i < num_a2; i++)
	{
		if(a1[0] == a2[i])
		{
			position = i;
			break;
		}
	}
  
	if(i == num_a2) return false; //element a1[0] does not appear in array a2;

	for(i = 1; i < num_a1; i++)
	{
		for(j = 1; j < num_a2; j++)
		{
			if(a1[i] == a2[(position+j)%num_a2])
				continue;

			else
				return false;
		}
	}

	return true;
}



/*
Judge whether two integer arrays are similar
*/
bool IsSimilarIntArray(int *a1, int num_a1, int *a2, int num_a2, double percent)
{
	int i, j;
	int position = 0;
	int sametriangles = 0;
	//int a2_position;

	for(i = 0; i < num_a1; i++)
	{
		for(j = 0; j < num_a2; j++)
		{
			if(a1[i] == a2[j])
			{
				position = j;
				sametriangles++;
				break;
			}
		}
	}
  
	//a2_position = 1;
	j = 1;
	for(i++; i < num_a1; i++)
	{
		for(; j < num_a2; j++)
		{
			if(a1[i] == a2[(position+j)%num_a2])
			{
				sametriangles++;
				//a2_position = j;
				break;
			}
		}
	}

	////similar
	if((double)sametriangles/(double)min(num_a1, num_a2) >= percent)
		return true;

	else
		return false;
}







/******************************************************************************************/


/////New boundary extraction 4/11/06 these following routines can be put to a library
void AllocBoundaryList()
{
	////set the maximum number of current boundaries in the list
	aboundarylist.MaxNumBoundaries = 10;

	////Build the boundary list
	aboundarylist.boundarylist = (RegionBoundary *)malloc(sizeof(RegionBoundary)*
		aboundarylist.MaxNumBoundaries);

	////
	int i;
	for(i = 0; i < aboundarylist.MaxNumBoundaries; i++)
	{
		aboundarylist.boundarylist[i].MaxNumEdgesOnBoundary = 1000;
		aboundarylist.boundarylist[i].edgelist = 
			(Edge**)malloc(sizeof(Edge*)*aboundarylist.boundarylist[i].MaxNumEdgesOnBoundary);
		aboundarylist.boundarylist[i].num = 0;
	}
}



/*
Remove one element from the given integer array if the element is in the array
Return true if success; otherwise, return false
*/
bool RemoveOneElem(int *a, int b, int &num)
{
	int i, pos = -1;

	for(i = 0; i < num; i++)
	{
		if(a[i] == b)
		{
			pos = i;
			break;
		}
	}

	if(pos == -1) return false;  //can not find the specific element

	//remove the element and reorganize the array
	for(i = pos; i < num-1; i++)
	{
		a[i] = a[i+1];
	}

	return true;
}

/*
Extract all the boundaries of the input region 
Input: region -- a triangulation mesh
       num -- number of triangles in this region
Output: flag -- if the region contain non-manifold vertex (singularity), set it as 1; otherwise 0
*/
void GetAllBoundaries(int *region, int &num)
{
	int i/*, j*/;
	//Face *face;
	Edge *cur_e;
	int EndVertID;
    Vertex *cur_vert;

	int num_visitededges = 0;
	int num_trianinregion = num;

	//flag = 0;  //set the initial flag

	////Reset the boundary flag 
LL:	aboundarylist.cur_boundary_num = 0;
	for(i = 0; i < aboundarylist.MaxNumBoundaries; i++)
	{
		aboundarylist.boundarylist[i].num = 0;
	}

	////Rebuild the boundary edge list
	BuildBoundaryEdgeList(region, num_trianinregion, 0);

	////Search all the boundaries of the region
	aboundarylist.boundarylist[aboundarylist.cur_boundary_num].edgelist[0]=
		Cycle_edge[0];
	aboundarylist.boundarylist[aboundarylist.cur_boundary_num].num = 1;

	EndVertID = Cycle_edge[0]->verts[0];
	Cycle_edge[0]->BoundaryVisited = 1;
	num_visitededges ++;

	cur_vert = Object.vlist[Cycle_edge[0]->verts[1]];

	while(num_visitededges <= num_cycleedges)   ////Not visit all the boundary edges
	{
		int num_edgesonboundary = 0;

		for(i = 0; i < cur_vert->Num_edge; i++)
		{
			cur_e = cur_vert->edges[i];
			if(cur_e->OnBoundary == 1 && cur_e->BoundaryVisited == 0)
			{
				if(num_edgesonboundary == 0)
				{
					aboundarylist.boundarylist[aboundarylist.cur_boundary_num].edgelist
						[aboundarylist.boundarylist[aboundarylist.cur_boundary_num].num]= cur_e;
					aboundarylist.boundarylist[aboundarylist.cur_boundary_num].num ++;
					cur_e->BoundaryVisited = 1; //set the visited flag

					num_visitededges++;
				}
				num_edgesonboundary++;
	 		}
	    }

		if(num_edgesonboundary % 2 == 0) ////there are odd number of edges on the boundary, wrong!!
		{
			//MessageBox(NULL, "singular vertex in boundary extraction", "Error", MB_OK);

			//flag = 1; //set flag

			//remove one of those triangles adjacent to currrent vertex and in side the region
			int oneinside = 0;
			Corner *c;

			for(i = 0; i < cur_vert->Num_corners; i++)
			{
				c = Object.clist[cur_vert->Corners[i]];

				if(!IsRepeated(region, c->t, num))
					continue;

				if(oneinside == 0)
				{
					oneinside ++;
					continue;
				}

				////remove all the other adjacent triangles that are inside the region
				if(!RemoveOneElem(region, c->t, num))
				{
					MessageBox(NULL, "Wrong removement!", "Error", MB_OK);
				}

				//Object.flist[c->t]->i
			}

			//flag = 1;  //we meet the case of non-manifold vertex, and we perform removing
			goto LL;
			//return;
		}

		////Get next testing vertex
		if(aboundarylist.boundarylist[aboundarylist.cur_boundary_num].edgelist
			[aboundarylist.boundarylist[aboundarylist.cur_boundary_num].num-1]->verts[0] != cur_vert->VertID)
		{
			cur_vert = Object.vlist[aboundarylist.boundarylist[aboundarylist.cur_boundary_num].edgelist
			                          [aboundarylist.boundarylist[aboundarylist.cur_boundary_num].num-1]->verts[0]];
		}
		else
		{
			cur_vert = Object.vlist[aboundarylist.boundarylist[aboundarylist.cur_boundary_num].edgelist
			                          [aboundarylist.boundarylist[aboundarylist.cur_boundary_num].num-1]->verts[1]];
		}
 
		if(cur_vert->VertID == EndVertID) ////form a loop
		{
			////we found a boundary
			aboundarylist.cur_boundary_num++;

			////Get an unvisited edge and continue the process

			for(i = 1; i < num_cycleedges; i++)
			{
				cur_e = Cycle_edge[i];

				if(cur_e->BoundaryVisited == 1)
					continue;

				aboundarylist.boundarylist[aboundarylist.cur_boundary_num].edgelist[0]=
					cur_e;
				aboundarylist.boundarylist[aboundarylist.cur_boundary_num].num = 1;

				EndVertID = cur_e->verts[0];
				cur_e->BoundaryVisited = 1;

				cur_vert = Object.vlist[cur_e->verts[1]];

				break;
			}

			if(i == num_cycleedges) //no more edges not being visited
				return;
		}
	}

}


////Suppose we have already know the boundary edges (not form boundaries) 4/15/06
////Boundary edges are stored in Cycle_edge data structure
bool ReachMeshBoundary()
{
	int i;
	Edge *cur_e;

	for(i = 0; i < num_cycleedges; i++)
	{
		cur_e = Cycle_edge[i];
		if(cur_e->tris[0] < 0 || cur_e->tris[1] < 0)
			return true;
	}
	return false;
}


////Suppose we have already extracted all the boundaries
////The following routine return a triangle at the outer boundary
////This only works for 2D planar case 4/15/06
int GetAnOuterBoundaryTriangle(int *regiontriangles, int nums)
{
	int i;
	int maxnum = aboundarylist.boundarylist[0].num;
	int outerboundaryid = 0;

	////we assume the outer boundary always has the maximum number of boundary edges
	for(i = 1; i < aboundarylist.cur_boundary_num; i++)
	{
		if(aboundarylist.boundarylist[i].num > maxnum)
		{
			maxnum = aboundarylist.boundarylist[i].num;
			outerboundaryid = i;
		}
	}

	////Get a triangle at the outer boudary

	Edge *cur_e = aboundarylist.boundarylist[outerboundaryid].edgelist
		[rand()%aboundarylist.boundarylist[outerboundaryid].num];

	if(IsRepeated(regiontriangles, cur_e->tris[0], nums))
		return cur_e->tris[0];

	else
		return cur_e->tris[1];	
}




/* ------------------4/6/06---------------------*/
void GetInitTriangleStrip(RegionBoundary aboundary, int type)
{
	int i;
	Edge *cur_e;

	if(type == 0)
	{
		UpdateBoundary(0);
		repellerRegion.num = 0;
		for(i = 0; i < aboundary.num; i++)
		{
			cur_e = aboundary.edgelist[i];

			if(IsRepeated(Source_re.trianglelist, cur_e->tris[0], Source_re.num))
			{
				repellerRegion.trianglelist[repellerRegion.num] = cur_e->tris[0];
				repellerRegion.num++;
			}
			else
			{
				repellerRegion.trianglelist[repellerRegion.num] = cur_e->tris[1];
				repellerRegion.num++;
			}
		}
	}

	else
	{
		UpdateBoundary(1);
		attractorRegion.num = 0;
		for(i = 0; i < aboundary.num; i++)
		{
			cur_e = aboundary.edgelist[i];

			if(IsRepeated(Sink_re.trianglelist, cur_e->tris[0], Sink_re.num))
			{
				attractorRegion.trianglelist[attractorRegion.num] = cur_e->tris[0];
				attractorRegion.num++;
			}
			else
			{
				attractorRegion.trianglelist[attractorRegion.num] = cur_e->tris[1];
				attractorRegion.num++;
			}
		}
	}
}


//Using new method to determine the embedded periodic orbits  4/11/06
void GetACellCycle_3(int singID, int type)
{
	int i;
	Edge **temp_edgelist;
	int num_boundedges = 0;
    Edge *oneboundedge;
	int aboundtriangle = singularities[singID].Triangle_ID;
	int multiboundaryornot = 0;
	int anouterboundtriangle = -1;

	////temporary region variables for boundary growing 4/11/06
	TriangularRegion temp_r1, temp_r2;
	temp_r1.MaxNumTrianglesInRegion = Object.nfaces;
	temp_r2.MaxNumTrianglesInRegion = Object.nfaces;
	temp_r1.trianglelist = (int*)malloc(sizeof(int)*temp_r1.MaxNumTrianglesInRegion);
	temp_r2.trianglelist = (int*)malloc(sizeof(int)*temp_r2.MaxNumTrianglesInRegion);
	////

	//a triangle at the inner boundary if exists 3/22/06
	//this will fail if the limit cycle cross the center triangle
	int aninerboundtriangle = Object.clist[3*singularities[singID].Triangle_ID]->ot; 
	
	double begin_p_out[2]/*, begin_p_in[2]*/;
	int inner = 0;
	int test_count = 0;
	int flag = -1;

	////Initialization part

	repellerRegion.num = 0;
	attractorRegion.num = 0;

	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->attract_inregion = 0;
		Object.flist[i]->repell_inregion = 0;
	}

	////first, grow a region from the triangle containing the center singularity
LL:	if(type == 0)  
	{
		//the type of the singularity at the center is a repeller, hence the limit cycle 
		//should be an attractor
		if(inner == 0)
		{
			repellerRegion.trianglelist[repellerRegion.num] = aboundtriangle;
			repellerRegion.num ++;
		}

		UpdateBoundary(0);
		GetRegionNormals(0);
		//Cancel_Growing(0, -1);
		if(!PeriodicDetect_Growing(type))
			goto L2;

		temp_edgelist = repellerBoundary.edgelist;
		num_boundedges = repellerBoundary.num;
		oneboundedge = repellerBoundary.edgelist[num_boundedges-1];
		aboundtriangle = repellerRegion.trianglelist[repellerRegion.num-1]; //this is a boundary triangle

		////Copy the region
		CopyRegion(repellerRegion.trianglelist, Source_re.trianglelist, repellerRegion.num);
		Source_re.num = repellerRegion.num;
	}

	else
	{
		if(inner == 0)
		{
			attractorRegion.trianglelist[attractorRegion.num] = aboundtriangle;
			attractorRegion.num ++;
		}

		UpdateBoundary(1);
		GetRegionNormals(1);
		//Cancel_Growing(1, -1);
		if(!PeriodicDetect_Growing(type))
			goto L2;

		temp_edgelist = attractorBoundary.edgelist;
		num_boundedges = attractorBoundary.num;
		oneboundedge = attractorBoundary.edgelist[num_boundedges-1];
		aboundtriangle = attractorRegion.trianglelist[attractorRegion.num-1]; //this is a boundary triangle
		
		////Copy the region
		CopyRegion(attractorRegion.trianglelist, Sink_re.trianglelist, attractorRegion.num);
		Sink_re.num = attractorRegion.num;
	}

	/***************************************************************************/
	////Temporary judgement 4/15/06
	//if(type == 0)
	//	BuildBoundaryEdgeList(repellerRegion.trianglelist, repellerRegion.num, type);
	//else
	//	BuildBoundaryEdgeList(attractorRegion.trianglelist, attractorRegion.num, type);

	//if(ReachMeshBoundary())
	//	return;
	/***************************************************************************/

	if(type == 0)
	{
		//extract all the boundaries of the obtained region 4/11/06
		GetAllBoundaries(repellerRegion.trianglelist, repellerRegion.num); 
	}
	else
	{
		//extract all the boundaries of the obtained region 4/11/06
		GetAllBoundaries(attractorRegion.trianglelist, attractorRegion.num); 
	}

	if(type == 0)
		anouterboundtriangle = 
			GetAnOuterBoundaryTriangle(repellerRegion.trianglelist, repellerRegion.num);
	else
		anouterboundtriangle = 
			GetAnOuterBoundaryTriangle(attractorRegion.trianglelist, attractorRegion.num);

	if(aboundarylist.cur_boundary_num == 0) //we reach the whole field/mesh boundary
		goto L2;

	/*--------------------------------------------------------------------------------*/
	////grow the backward regions from all extracted boundaries respectively 4/11/06

	////We need to grow the region from the first boundary before performing intersection operation
	GetInitTriangleStrip(aboundarylist.boundarylist[0], type);
	if(type == 0)
	{
		UpdateBoundary(0);
		GetRegionNormals(0);
		//Cancel_Growing(0, -1);
		if(!PeriodicDetect_Growing(type))
			goto L2;

		CopyRegion(repellerRegion.trianglelist, temp_r1.trianglelist, repellerRegion.num);
		temp_r1.num = repellerRegion.num;
	}
	else 
	{
		UpdateBoundary(1);
		GetRegionNormals(1);
		//Cancel_Growing(1, -1);
		if(!PeriodicDetect_Growing(type))
			goto L2;
		CopyRegion(attractorRegion.trianglelist, temp_r1.trianglelist, attractorRegion.num);
		temp_r1.num = attractorRegion.num;
	}

	////Then we grow regions from other boundaries and intersect them with previously obtained region 
	for(i = 1; i < aboundarylist.cur_boundary_num; i++)
	{
		GetInitTriangleStrip(aboundarylist.boundarylist[i], type);

		if(type == 0)
		{
			UpdateBoundary(0);
			GetRegionNormals(0);
			//Cancel_Growing(0, -1);
			if(!PeriodicDetect_Growing(type))
				goto L2;

			IntersectRegion(temp_r1, repellerRegion, temp_r2);
		}
		else
		{
			UpdateBoundary(1);
			GetRegionNormals(1);
			//Cancel_Growing(1, -1);
			if(!PeriodicDetect_Growing(type))
				goto L2;

			IntersectRegion(temp_r1, attractorRegion, temp_r2);
		}

		CopyRegion(temp_r2.trianglelist, temp_r1.trianglelist, temp_r2.num);
		temp_r1.num = temp_r2.num;
	}
	/*--------------------------------------------------------------------------------*/


	/*--------------------------------------------------------------------------------*/
	////Intersect with previous main region to get the ring-shaped region 4/11/06

	if(type == 0)
	{
		////calculate the intersection region
		IntersectRegion(temp_r1, Source_re, intersectRegion);
	}

	else
	{
		////calculate the intersection region
		IntersectRegion(temp_r1, Sink_re, intersectRegion);
	}


	/*--------------------------------------------------------------------------------*/
	////After calculate the intersection of two regions, we get a ring-shaped region
	////That is a little bit larger than the cell cycle containing the limit cycle

	if(intersectRegion.num == 0)
		goto L2;

	//select a random triangle inside the ring-shaped region to find the fixed point
	aboundtriangle = intersectRegion.trianglelist[rand()%intersectRegion.num];

	////Locate the real cell cycle and closed orbit

	GetCenterofATriangle(aboundtriangle, begin_p_out);
	//GetCenterofATriangle(aninerboundtriangle, begin_p_in);

	DetectALimitCycle_new(begin_p_out, aboundtriangle, type, inner, singID, flag);
	
	//////For embedded limit cycles, we should begin from the previous selected triangle and perform 
	//////backward region growing!!(which means repeller growing will become attractor growing)

	//////first, we need to reverse the type
	//////if the selected boundary triangle is not at the boundary of the whole field
	if(IsAtFieldBoundary(anouterboundtriangle) || test_count == 4 || flag == 2)
	{
		goto L2;
	}

	////Interlace copy here !!! 1/11/06
	if(type == 0)
	{
		CopyRegion(Source_re.trianglelist, attractorRegion.trianglelist, Source_re.num);
		attractorRegion.num = Source_re.num;
	}

	else
	{
		CopyRegion(Sink_re.trianglelist, repellerRegion.trianglelist, Sink_re.num);
		repellerRegion.num = Sink_re.num;
	}

	type = 1 - type;

	test_count ++;

	inner ++;

	goto LL;

L2:	/* release the allocated memory */	
	free(temp_r1.trianglelist);
	free(temp_r2.trianglelist);
}



/*Testing routine to show all the boundaries*/
void ShowAllBoundaries()
{
	int i, j;
	Edge *cur_e;
	Vertex *v;

	glLineWidth(2.8);
	for(i = 0; i < aboundarylist.cur_boundary_num; i++)
	{
		if(i % 3 == 0)
			glColor3f(1, 0, 0);
		else if(i % 3 == 1)
			glColor3f(0, 1, 0);
		else
			glColor3f(0, 0, 1);

		for(j = 0; j < aboundarylist.boundarylist[i].num; j++)
		{
			cur_e = aboundarylist.boundarylist[i].edgelist[j];

			glBegin(GL_LINES);
			v = Object.vlist[cur_e->verts[0]];
			glVertex2f(v->x, v->y);
			v = Object.vlist[cur_e->verts[1]];
			glVertex2f(v->x, v->y);
			glEnd();
		}
	}
	glLineWidth(1.);
}





/******************************************************************************/


/*
Grow regions from all boundaries, return the intersect region of all these region
*/
bool GetOneRingShapedRegion(int *ring, int &num, int type)
{
	int i;
	int flag = 0;

	/*Temporary region variable*/
	TriangularRegion temp_r1, temp_r2;
	temp_r1.MaxNumTrianglesInRegion = Object.nfaces;
	temp_r2.MaxNumTrianglesInRegion = Object.nfaces;
	temp_r1.trianglelist = (int*)malloc(sizeof(int)*temp_r1.MaxNumTrianglesInRegion);
	temp_r2.trianglelist = (int*)malloc(sizeof(int)*temp_r2.MaxNumTrianglesInRegion);

	temp_r1.num = temp_r2.num = 0;

	////We need to grow the region from the first boundary before performing intersection operation
	GetInitTriangleStrip(aboundarylist.boundarylist[0], type);
	if(type == 0)
	{
		UpdateBoundary(0);
		GetRegionNormals(0);
		PeriodicDetect_Growing(0, -1, flag);
		if(flag == 1)
		{
			free(temp_r1.trianglelist);
			free(temp_r2.trianglelist);
			return false;
		}

		CopyRegion(repellerRegion.trianglelist, temp_r1.trianglelist, repellerRegion.num);
		temp_r1.num = repellerRegion.num;
	}
	else 
	{
		UpdateBoundary(1);
		GetRegionNormals(1);
		PeriodicDetect_Growing(1, -1, flag);
		if(flag == 1)
		{
			free(temp_r1.trianglelist);
			free(temp_r2.trianglelist);
			return false;
		}

		CopyRegion(attractorRegion.trianglelist, temp_r1.trianglelist, attractorRegion.num);
		temp_r1.num = attractorRegion.num;
	}

	////Then we grow regions from other boundaries and intersect them with previously obtained region 
	for(i = 1; i < aboundarylist.cur_boundary_num; i++)
	{
		flag = 0;

		GetInitTriangleStrip(aboundarylist.boundarylist[i], type);

		if(type == 0)
		{
			UpdateBoundary(0);
			GetRegionNormals(0);
			PeriodicDetect_Growing(0, -1, flag);
			if(flag == 1)
			{
				free(temp_r1.trianglelist);
				free(temp_r2.trianglelist);
				return false;
			}

			IntersectRegion(temp_r1, repellerRegion, temp_r2);
		}
		else
		{
			UpdateBoundary(1);
			GetRegionNormals(1);
			PeriodicDetect_Growing(1, -1, flag);
			if(flag == 1)
			{
				free(temp_r1.trianglelist);
				free(temp_r2.trianglelist);
				return false;
			}

			IntersectRegion(temp_r1, attractorRegion, temp_r2);
		}

		CopyRegion(temp_r2.trianglelist, temp_r1.trianglelist, temp_r2.num);
		temp_r1.num = temp_r2.num;
	}
	
	CopyRegion(temp_r2.trianglelist, ring, temp_r2.num);
	num = temp_r2.num;

	free(temp_r1.trianglelist);
	free(temp_r2.trianglelist);
	return true;
}



/*
*/
void GetBackwardRegion(int *region, int num, int type, int *backwardregion, int num_backward)
{
	//Get all the possible boundaries of the input region
	GetAllBoundaries(region, num);

	///*The following can be changed according to the method we will use*/
	//Grow regions from all the boundaries using the same region growing (exit or entrance) 
	//For multiple region intersection method,
	//We intersect all the obtained regions, and return the obtain region
	//
	//For single boundary method,
	//We pick the region that has ring-shaped, and return this region.
	GetOneRingShapedRegion(backwardregion, num_backward, type);
}

void GetCellCycles(int singID, int type)
{
	//int i;
	int aboundtriangle = singularities[singID].Triangle_ID;
	
	double begin_p_out[2]/*, begin_p_in[2]*/;

	int inner = 0;

	int flag = 0;

	int reachboundary = 0;  //judge whether the region growing reach boundary or not

	repellerRegion.num = 0;
	attractorRegion.num = 0;

	//set the initial region as the triangle containing the input singularity
	if(type == 0)
	{
		repellerRegion.trianglelist[0] = singularities[singID].Triangle_ID;
		repellerRegion.num++;
	}
	else
	{
		attractorRegion.trianglelist[0] = singularities[singID].Triangle_ID;
		attractorRegion.num++;
	}

	//get the regions and the limit cycles
	while(inner <= 6)
	{

		////first, grow a region from the initial region according to current type of region growing
		if(type == 0)  
		{
			UpdateBoundary(0);
			GetRegionNormals(0);
			PeriodicDetect_Growing(0, singID, reachboundary);

			////Copy the region
			CopyRegion(repellerRegion.trianglelist, Source_re.trianglelist, repellerRegion.num);
			Source_re.num = repellerRegion.num;
		}

		else
		{
			UpdateBoundary(1);
			GetRegionNormals(1);
			PeriodicDetect_Growing(1, singID, reachboundary);

			////Copy the region
			CopyRegion(attractorRegion.trianglelist, Sink_re.trianglelist, attractorRegion.num);
			Sink_re.num = attractorRegion.num;
		}




		//Call GetBackwardRegion() to get a backward growing region
		if(type == 0)
		{
			GetBackwardRegion(Source_re.trianglelist, Source_re.num, 0,	
				repellerRegion.trianglelist, repellerRegion.num);
			
			//calculate the intersection region
			IntersectRegion(repellerRegion, Source_re, intersectRegion);
		}
		else
		{
			GetBackwardRegion(Sink_re.trianglelist, Sink_re.num, 1,
				attractorRegion.trianglelist, attractorRegion.num);
			
			//calculate the intersection region
			IntersectRegion(attractorRegion, Sink_re, intersectRegion);
		}

		//Get one triangle from the intersection region, get the center point of it
		if(intersectRegion.num == 0) //region growing failed
			return;

		aboundtriangle = intersectRegion.trianglelist[rand()%intersectRegion.num];

		GetCenterofATriangle(aboundtriangle, begin_p_out);

		//Call FindALimitCycle() to get the limit cycle

		//if(FindALimitCycle(intersectRegion.trianglelist, intersectRegion.num,
		//	begin_p_out, aboundtriangle, inner, singID, type))
		//{
		//	//save the limit cycle information and build the connection relation
		//}

		//else //not detect a limit cycle, abort
		//{
		//	return;
		//}

		////Now we try to call the old routine to get the limit cycle 5/28/06
		flag = 0;
		DetectALimitCycle_new(begin_p_out, aboundtriangle, type, inner, singID, flag);

		if(flag == 2)
			return;

		if(reachboundary == 1) /* 06/18/06 */
			return;  //already reach boundary ! 06/18/06

		//////////////////////////////05/28/06////////

		//Set previously obtained region as the current initial region, switch the type of region growing
		//This part may have problem
		//should we copy only the ring-shaped region or the whole obtained region ? 06/18/06
		if(type == 0)
		{
			//CopyRegion(repellerRegion.trianglelist, attractorRegion.trianglelist, repellerRegion.num);
			//attractorRegion.num = repellerRegion.num;
			CopyRegion(Source_re.trianglelist, attractorRegion.trianglelist, Source_re.num);
			attractorRegion.num = Source_re.num;
		}

		else
		{
			//CopyRegion(attractorRegion.trianglelist, repellerRegion.trianglelist, attractorRegion.num);
			//repellerRegion.num = attractorRegion.num;
			CopyRegion(Sink_re.trianglelist, repellerRegion.trianglelist, Sink_re.num);
			repellerRegion.num = Sink_re.num;
		}

		type = 1 - type; //switch the type

		inner ++; //not an inner limit cycle any more

	}
}





/*
Grow regions from all possible boundaries, return the region that has ring-shaped
*/
bool GetOneRingShapedRegion2(int *ring, int num, int type)
{
	return false;
}

/*
One tracing step testing for finding a good triangle associated with the edge
*/
bool OneStepTest(int v1, int v2, int triangle, int type)
{
	double middle[2];
	icVector2 vec;
	double alpha[3], a, b;
	Face *face = Object.flist[triangle];

	//Get the middle point of the edge v1v2
	middle[0] = (Object.vlist[v1]->x + Object.vlist[v2]->x)/2.;
	middle[1] = (Object.vlist[v1]->y + Object.vlist[v2]->y)/2.;

	vec = 0.5*Object.vlist[v1]->vec + 0.5*Object.vlist[v2]->vec;

	//trace from the middle one step and get a new point
	if(type == 0)
	{
		middle[0] += vec.entry[0]/8.;
		middle[1] += vec.entry[1]/8.;
	}
	else
	{
		middle[0] -= vec.entry[0]/8.;
		middle[1] -= vec.entry[1]/8.;
	}

	//test whether the point still falls in the triangle
	vec.entry[0] = middle[0] - Object.vlist[face->verts[0]]->x;
	vec.entry[1] = middle[1] - Object.vlist[face->verts[0]]->y;

	a = dot(vec, face->LX);
	b = dot(vec, face->LY);

	Get2DBarycentricFacters(triangle, a, b, alpha);

	////2. if current point is inside current triangle
	if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1 
		&& (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1 
		&& (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
	{
		return true;
	}

	return false;
}

Edge *GetAGoodEdge(int *cellcycle, int num_cells, int type, int &triangle)
{
	int i,j;
	Face *face;
	Edge *cur_e;

	for(i = 0; i < num_cells; i++)
	{
		face = Object.flist[cellcycle[i]];

		for(j = 0; j < 3; j++)
		{
			cur_e = face->edges[j];

			if(!IsRepeated(cellcycle, cur_e->tris[0], num_cells)
				|| !IsRepeated(cellcycle, cur_e->tris[1], num_cells))
				continue; //not the inner edge

			else
			{
				//there are two possible triangles can be associated with the edge
				//we want the tracing curve starts from the middle point of the edge
				//will enter the triangle rather than leave the triangle 
				//this will provide convenience for the fixed point finding
				triangle = cur_e->tris[0];
				if(!OneStepTest(cur_e->verts[0], cur_e->verts[1], cur_e->tris[0], type))
					triangle = cur_e->tris[1];

				return cur_e;
			}
		}
	}
}


/*
*/
bool OneCellCycleDetect(int *trianglestrip, int &num_triangles, double bgp[2], int &begin_triangle, int type)
{
	int i;
	int pre_face, cur_face;
	double globalp[2];
	//int position;
	int flag;

	////
	int *cellcycle = (int*)malloc(sizeof(int)*num_triangles);
	int num_cells = 0;

	////Initialization part
	globalp[0] = bgp[0];
	globalp[1] = bgp[1];
	pre_face = cur_face = begin_triangle;

	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->access_count = 0;
		Object.flist[i]->discard = 0;
	}

	////Perform local tracing and mark the triangles the curve passes through 
	Object.flist[cur_face]->access_count++;

	for(i = 0; i < /*2*NUMTRACINGTRIANGLE*/3*Object.nfaces; i++)
	{
		if(cur_face == -1)
		{
			free(cellcycle);
			return false;
		}

		pre_face = cur_face ;
		cur_face = TraceInATriangle2(cur_face, globalp, type, flag);
 
		if(flag == 3 || flag == 4 || !IsRepeated(trianglestrip, cur_face, num_triangles)) 
		{
			free(cellcycle);
			return false;	////reach a singularity or the boundary
		}

		else if(Object.flist[cur_face]->access_count >= 4)
		{
			Object.flist[cur_face]->access_count-= 1;
			bgp[0] = globalp[0];
			bgp[1] = globalp[1];
			begin_triangle = cur_face;

			//we need to store the cell cycle here (not do this now!!5/23/06)

			//call the old method to perform one more round tracing
			First_CellCycleDetect(bgp[0], bgp[1], cur_face, cellcycle, num_cells, type, flag);

			//copy the obtained cell cycle to the triangle strip
			CopyRegion(cellcycle, trianglestrip, num_cells);
			num_triangles = num_cells;

			free(cellcycle);
			return true;
		}

		else{
			Object.flist[cur_face]->access_count++;
		}
	}

	free(cellcycle);
	return false;
}



////Build the graphical handler for the limit cycle for user interaction
void BuildHandlerforLimitCycle(int cur_limitcycle_index, double fixedpoint[2], int ver)
{
	Vertex *ver1;
	ver1 = Object.vlist[ver];

	////Using the input fixed point as the legend center 05/23/06
	limitcycles[cur_limitcycle_index].legend_center[0] = fixedpoint[0];
	limitcycles[cur_limitcycle_index].legend_center[1] = fixedpoint[1];

	////Get the base of the legend
	limitcycles[cur_limitcycle_index].legend_base[0] = ver1->x;
	limitcycles[cur_limitcycle_index].legend_base[1] = ver1->y;
}


////Build the graphical handler for the limit cycle based on the fixed point
//void BuildHandlerforLimitCycle(int cur_limitcycle_index)
//{
//	Vertex *ver1;
//	////Using the next level beginning point as the legend center 08/25/05
//	ver1 = Object.vlist[OneSharedEdge->verts[0]];
//	limitcycles[cur_limitcycle_index].legend_center[0] = ver1->x;
//	limitcycles[cur_limitcycle_index].legend_center[1] = ver1->y;
//
//	////Get the base of the legend
//	if(ver1->VertID != OneSharedEdge->verts[0])
//		ver1 = Object.vlist[OneSharedEdge->verts[0]];
//	else
//		ver1 = Object.vlist[OneSharedEdge->verts[1]];
//
//	limitcycles[cur_limitcycle_index].legend_base[0] = ver1->x;
//	limitcycles[cur_limitcycle_index].legend_base[1] = ver1->y;
//}

/*
*/
bool FindALimitCycle(int *trianglestrip, int num_triangles, double bgp[2], int begin_triangle, 
					 int inner, int singID, int type)
{
	int pre_limitcycle_id = -1;
	double fixedpoint[2] = {0.};
	int fix_triangle = -1;
	//Edge *theEdge;
	int OneEndVert;
	
	//Call OneCellCycleDetect() to locate a one-triangle width triangle strip
	if(!OneCellCycleDetect(trianglestrip, num_triangles, bgp, begin_triangle, type))
		return false; //not detect any limit cycle

	//Compare with previous limit cycle
	if(IsPreviousDetectedLimitCycle(trianglestrip, num_triangles, pre_limitcycle_id))
	{
		//need to build connection with other elements
		if(inner == 0) //the center singularity can only connect with the most inner limit cycle
		{
			UpdateListInSingularity(singID, pre_limitcycle_id);
			UpdateListInLimitCycle(pre_limitcycle_id, singID);
		}
		else //if it is an outer limit cycle, think about that here !!
		{
			UpdateCycleListInLimitCycle(cur_limitcycle_index-1, pre_limitcycle_id);
			UpdateCycleListInLimitCycle(pre_limitcycle_id, cur_limitcycle_index-1);
		}

		return true;
	}

	else
	{
		if(FindFixedPoint(trianglestrip, num_triangles, type, fixedpoint, fix_triangle, OneEndVert))
		{
			//Get the closed orbit from the fixed point and associated triangle
			CalLocalTracing(fix_triangle, fixedpoint[0], fixedpoint[1], type);

			//Get other information for the limit cycle
			
			/*---------------------------------------------------------------*/
			////Store the information of the new limit cycle
					
			StoreCurrentCellCycle(trianglestrip, num_triangles);  

			limitcycles[cur_limitcycle_index].singularID = singID;  ////store the center singularity of current limit cycle
			
			limitcycles[cur_limitcycle_index].node_index = -1; ////10/16/05

			StoreCurStreamline(trajectories[cur_traj_index], num_linesegs_curtraj[cur_traj_index]);     

			////store the type of the limit cycle
			if(type == 0)
				limitcycles[cur_limitcycle_index].type = 1;
			else
				limitcycles[cur_limitcycle_index].type = 0;

			////Initialize the connected list for the graph
			limitcycles[cur_limitcycle_index].connected_limitcycle = NULL;
			limitcycles[cur_limitcycle_index].num_connectedcycles = 0;
			limitcycles[cur_limitcycle_index].connected_saddle = NULL;
			limitcycles[cur_limitcycle_index].num_connectedsaddles = 0;

			/*---------------------------------------------------------------*/

			//Build the conley connection for the limit cycle
			
			if(inner == 0) //if it is the most inner limit cycle
			{
				UpdateListInSingularity(singID, cur_limitcycle_index);
				UpdateListInLimitCycle(cur_limitcycle_index, singID);
			}

			else //if it is an outer limit cycle
			{
				UpdateCycleListInLimitCycle(cur_limitcycle_index-1, cur_limitcycle_index);
				UpdateCycleListInLimitCycle(cur_limitcycle_index, cur_limitcycle_index-1);
			}
			/*---------------------------------------------------------------*/

			////Build the graphical handler for the limit cycle for user interaction
			
			BuildHandlerforLimitCycle(cur_limitcycle_index, fixedpoint, OneEndVert); //it depends on the cell cycle also

			//cur_traj_index ++;
			cur_limitcycle_index++;
		}

		else //we can not find the fixed point for this edge, probably we can choose another edge?
		{
		}
	}
}


/*
Test whether the new detected limit cycle is one of the previously detected one
Here, we still use cell cycle to decide. In future, we need to consider to 
use the sample points on the closed orbits
*/
bool IsPreviousDetectedLimitCycle(int *cycle, int num, int &pre_limitcycle_id)
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

		for(j = 0; j < num; j++)
		{
			cur_triangle = cycle[j];

			if(!TriangleSearch(limitcycles[i].cellcycle, limitcycles[i].num_triangles, cur_triangle, position))
				goto LL;   ////continue to test next limit cycle

			////if it can be found
			count ++;
		}

		if(count == limitcycles[i].num_triangles - 1) ////using a loosen pending condition
		{
			pre_limitcycle_id = i;
			return true;
		}


		else
LL:			continue;
	}

	return false;
}

/*
*/
bool FindFixedPoint(int *cellcycle, int num_cells, int type, double fixedpoint[2], 
					int &fix_triangle, int &OneEndVert)
{
	icVector2 dist;
	int cur_triangle;
	double pre_intersection[2] = {0.};
	fix_triangle = -1;
	int flag;

	//Find the edge inside the cell cycle
	Edge *theEdge = GetAGoodEdge(cellcycle, num_cells, type, fix_triangle);

	pre_intersection[0] = fixedpoint[0] = (Object.vlist[theEdge->verts[0]]->x + Object.vlist[theEdge->verts[1]]->x)/2.;
	pre_intersection[1] = fixedpoint[1] = (Object.vlist[theEdge->verts[0]]->y + Object.vlist[theEdge->verts[1]]->y)/2.;

	cur_triangle = fix_triangle;
	//begin from the middle point of the edge
	do{
		//trace one round until it returns 'fix_triangle' again

		do{
			cur_triangle = TraceInATriangle(cur_triangle, fixedpoint, type, flag);
			if(flag == 3 || flag == 4)
				return false;
		}while(cur_triangle != fix_triangle);

		dist.entry[0] = fixedpoint[0] - pre_intersection[0];
		dist.entry[1] = fixedpoint[1] - pre_intersection[1];

		pre_intersection[0] = fixedpoint[0];
		pre_intersection[1] = fixedpoint[1];
		
	}while(length(dist)>1e-6);
	
	OneEndVert = theEdge->verts[0]; //return one vertex of the selected edge for making the legend of the limit cycle
	return true;
}

