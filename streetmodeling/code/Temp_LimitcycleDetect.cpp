

////A temporary file for limit cycle detection 1/4/06
////Temp_LimitcycleDetect.cpp 
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

extern int *TriangleList;
extern int *cellcycle;
extern int num_trianglelist;                ////number of triangles in the triangle list
extern int num_celltriangle;                ////number of triangles in the found cell cycle
extern Edge **Cycle_edge;                   ////mesh edges of the boundary of the cell cycle
//extern Vertex **PotentialVerts;             ////Potential exit points (vertices) on the boundary of cell cycle
//extern double realexit[2];                  ////we may need the information of vertex!!!!05/30/05
extern Vertex *realexitvert;                ////currently founded real exit
extern int num_cycleedges;                  ////number of edges consist of the cell cycle 
extern int num_potentialverts;              ////number of vertices consist of the cell cycle
extern Edge *OneSharedEdge;                 ////one shared edge in the cell cycle for finding fixed point
extern double theFixPoint[2];               ////the fix point on the closed streamline
extern Edge *SecondSharedEdge;              ////another shared edge for finding next level beginning point

extern int Cur_CellCycleIndex;                    ////current cell cycle index

extern Vertex * GetFurtherVer(Edge *sharededge, int singID);
extern void StoreCurStreamline(LineSeg *streamline, int num);




/////Testing variables 01/5/06
extern double MarkNextBx[2];

////To reduce the calculation of intersections, we may need to store previous intersection 
////and compare it with current intersection, here is the temporary global variable for this purpose
////1/6/06
extern double first_intersection[2] ;
extern int first_intersect ;

/*------------------------------------------------------------------------------*/
////The following two routines are for the calculation of the closed orbit
void BinaryFindtheFixedPoint_new(double v1[2], double v2[2], int type, double FixPoint[2], int &triangleID)
{
	int i = 0;  
	double alpha;
	double midx, midy, interx, intery;
	icVector2 dis;
	int flag = 0;

	interx = intery = 0;

	triangleID = OneSharedEdge->tris[0];  //we pick any one of the two faces

	first_intersect = 0;   ////new added 1/6/06

	while(i < 200)    ////avoid dead loop
	{
		flag = 0;
		midx = (v1[0]+v2[0])/2.;
	    midy = (v1[1]+v2[1])/2.;

 		local_GetNextIntersection(midx, midy, triangleID, interx, intery, type, flag);

		if(flag == 1)
		{
		    MessageBox(NULL, "Can not find next intersection!", "error", MB_OK);
			return;
		}

		if(flag == 2) ////simply use current point as the fixed point
		{
			FixPoint[0] = interx;
			FixPoint[1] = intery;
			return;
		}

		if(flag == 3) ////we have reached the fixed point during calculating the intersection
		{
			FixPoint[0] = interx;
			FixPoint[1] = intery;
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

		else if(alpha < -1e-9)
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



////after we locate a cell cycle, using binary search to find the 
////accurate closed streamline
void GettheClosedStreamline_new(int repellorattract, double FixPoint[2], int &triangleID)
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
	BinaryFindtheFixedPoint_new(v1, v2, repellorattract, FixPoint, triangleID);
}

///////////////
//Testing variables 1/5/06
int firstlimit = 0;

bool CloseToSingularity(double point[2], int singID)
{
	if(fabs(point[0] - singularities[singID].gcx) < 1e-4
		&& fabs(point[1] - singularities[singID].gcy) < 1e-4)
		return true;
	else
		return false;
}


bool IsASingularity(double alpha[3], int triangle)
{
    icVector2 v = GetVectorAtPoints(triangle, alpha);
	
	if(length(v) < 1e-9) return true;

	return false;
}

void ChooseBeginningPointForRegularSing(int singID, double begin_p[2], int &begin_triangle)
{
	int i;
	int triangle = singularities[singID].Triangle_ID;
	Face *face = Object.flist[triangle];
	Edge *cur_e;
	Vertex *v1, *v2;

	////we pick the middle point of one the edges of the triangle
	////make sure it is not too close to the singularity
	//for(i = 0; i < 3; i++)
	//{
	//	cur_e = face->edges[i];
	//	v1 = Object.vlist[cur_e->verts[0]];
	//	v2 = Object.vlist[cur_e->verts[1]];

	//	begin_p[0] = (v1->x + v2->x)/2.;
	//	begin_p[1] = (v1->y + v2->y)/2.;

	//	if(!CloseToSingularity(begin_p, singID))
	//		return;
	//}


	////$$second way$$, we pick the center of any one of the 3 neighboring triangles of current triangle
	Face *other_f;
	Vertex *v3;
	double alpha[3] = {0.3, 0.3, 0.4};
	for(i = 0; i < 3; i++)
	{
		cur_e = face->edges[i];

		if(cur_e->tris[0] != face->index)
			other_f = Object.flist[cur_e->tris[0]];
		else
			other_f = Object.flist[cur_e->tris[1]];

		v1 = Object.vlist[other_f->verts[0]];
		v2 = Object.vlist[other_f->verts[1]];
		v3 = Object.vlist[other_f->verts[2]];

		begin_p[0] = 0.3*v1->x + 0.3*v2->x + 0.4*v3->x;
		begin_p[1] = 0.3*v1->y + 0.3*v2->y + 0.4*v3->y;

		if(!IsASingularity(alpha, other_f->index))
		{
			begin_triangle = other_f->index;
			return;
		}
	}
}

////Get the beginning point in 2D plane
void GetBeginPoint(int singID, int type, double begin_p[2], int &begin_triangle, int Saddleornot, icVector2 sep_vec)
{
	double sing_center[2] = {singularities[singID].gcx, singularities[singID].gcy};
	//icVector

	if(Saddleornot == 0)
	{
		
		//begin_p[0] = sing_center[0] + SEPARATRIXSTEP/2;
		//begin_p[1] = sing_center[1] + SEPARATRIXSTEP/2;
		
		//sep_vec.entry[0] = SEPARATRIXSTEP;
		//sep_vec.entry[1] = SEPARATRIXSTEP;
		//normalize(sep_vec);
		//FixSepBeginningPos(singularities[singID].Triangle_ID, sep_vec, sing_center, begin_p);
		ChooseBeginningPointForRegularSing(singID, begin_p, begin_triangle);
	}

	else
	{
		FixSepBeginningPos(singularities[singID].Triangle_ID, sep_vec, sing_center, begin_p);
	}
	//FixSepBeginningPos(singularities[singID].Triangle_ID, sep_vec, sing_center, begin_p);

	////Testing codes 1/5/06
	if(firstlimit == 0)
	{
		MarkNextBx[0] = begin_p[0];
		MarkNextBx[1] = begin_p[1];
		firstlimit++;
	}
}


////This method may not be suitable for surfaces !!!!!!!!
void GetNextLevelBeginPoint(int singID, double begin_p[2], int type)
{
	Vertex *ver1;

	ver1 = GetFurtherVer(OneSharedEdge, singID);
	begin_p[0] = ver1->x;
	begin_p[1] = ver1->y;

	icVector2 temp_v;
	temp_v.entry[0] = begin_p[0] - singularities[singID].gcx;
	temp_v.entry[1] = begin_p[1] - singularities[singID].gcy;
	normalize(temp_v);
	
	begin_p[0] += 0.01 * temp_v.entry[0];  //we should use adaptive stepsize 1/2/06
	begin_p[1] += 0.01 * temp_v.entry[1];

    OneBackwardRungeKutta(begin_p[0], begin_p[1], begin_p[0], begin_p[1], type);
}

////Build the conley connection after detection of a limit cycle
void BuildConleyConnectionforLimitCycle(int index1, int index2, int saddlelimitornot)
{
}


//////////////////////////////////////////////////////////
////Not just for saddle, other kind of singularities can also connect to limit cycle!!!
////It seems that following codes do not check the repeating limit cycles in the list
////This is what we want for double connections (3/11/06)
void UpdateListInSingularity(int singID, int limitcycle)
{
	////1.Extend current list
	////2.Add to the end of current list
	////3.Update the counter

	int i;
	//int *temp_list = limitcycles[limitcycle].connected_saddle;
	int *temp_list = singularities[singID].connected_limitcycles;

	if(temp_list == NULL)
	{
		singularities[singID].connected_limitcycles = (int*)malloc(sizeof(int));
		singularities[singID].connected_limitcycles[0] = limitcycle;
		singularities[singID].num_connected_limitcycles = 1;
	}

	else{
		if(IsRepeated(singularities[singID].connected_limitcycles,
			limitcycle, singularities[singID].num_connected_limitcycles))
			return;

		singularities[singID].connected_limitcycles =
			(int*)malloc(sizeof(int)*(singularities[singID].num_connected_limitcycles+1));

		for(i = 0; i < singularities[singID].num_connected_limitcycles; i++)
		{
			singularities[singID].connected_limitcycles[i] = temp_list[i];
		}

		singularities[singID].connected_limitcycles[i] = limitcycle;

		singularities[singID].num_connected_limitcycles += 1;

		free(temp_list);
	}
}



////Build the graphical handler for the limit cycle for user interaction
void BuildHandlerforLimitCycle(int cur_limitcycle_index)
{
	Vertex *ver1;
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
}


//Detect from a specified singularity
void DetectFromASingularity(int singID, int type, int Saddleornot)
{
	double begin_p[2];
	icVector2 sep_vec;
	int index_separatrix = 0;
	int flag;
	int inner = 0;

	int fix_triangleID, begin_triangle;

	double center[2] = {singularities[singID].gcx, singularities[singID].gcy};
	begin_triangle = singularities[singID].Triangle_ID;
	
	//Get the beginning point according to the type and the center of the singularity
	//Note that we now tend to make sure that the begining point falls in the triangle containing
	//the singularity
	if(Saddleornot == 0)
	{
		GetBeginPoint(singID, type, begin_p, begin_triangle, Saddleornot, sep_vec);
	}

	else   //if it is a saddle, they are somewhat different 
	{
L1:		switch(index_separatrix)
		{
			case 0: sep_vec = singularities[singID].outgoing; 
					type = 0;
					break;

			case 1: sep_vec = -singularities[singID].outgoing; 
					type = 0;
					break;
			
			case 2: sep_vec = singularities[singID].incoming; 
					type = 1;
					break;

			case 3: sep_vec = -singularities[singID].incoming; 
					type = 1;
					break;

			default:
				   return;
		}
		GetBeginPoint(singID, type, begin_p, begin_triangle, Saddleornot, sep_vec);

		index_separatrix ++;
	}

	
	// Initialize for the detection
L2:	InitLimitCycleDetection();

	flag = -1;

	//Find cell cycle, it may be an unstable part! 
	//Double_CellCycleDetect(begin_p[0], begin_p[1], begin_trinagle, type, flag); 
	Double_CellCycleDetect(begin_p[0], begin_p[1], begin_triangle, type, flag); 

	if(flag == 2) //reach a singularity or the boundary of the mesh
	{
		if(Saddleornot == 0)
			return;
		else      //we need to test the other separatrix if there is still one not being visited
		{
			goto L1; //?????
		}
	}

	if(flag == 1) //we found a cell cycle
	{
		int pre_limit_index = -1;

		if(num_celltriangle == 0)
			return;

		//if it is a previously detected limit cycle
		if(IsPreviousLimitCycle(pre_limit_index))
		{
			UpdateListInSingularity(singID, pre_limit_index);
			UpdateListInLimitCycle(pre_limit_index, singID);

			if(Saddleornot == 1 && index_separatrix < 3)
				goto L1;

			return;
		}

		/*---------------------------------------------------------------*/
		//if it is a new extracted limit cycle
		boundaryBuilding();
       
		GetASharedEdge(center);

		StoreCurrentCellCycle(cellcycle, num_celltriangle);  
		Cur_CellCycleIndex++;   //do we have multiple cell cycle during one limit cycle detection process

		//Get the closed orbit, this is the unstable part of the whole detection
		GettheClosedStreamline_new(type, theFixPoint, fix_triangleID);

		////since some times, the program can not find the fixed point, you need to deal with 
		////the beginning point is origin 1/5/06
		if((theFixPoint[0] == 0 && theFixPoint[1] == 0)||fix_triangleID < 0)
		{
			fix_triangleID = TriangleDetect(theFixPoint[0], theFixPoint[0]);
		}
			

		//calculate the closed streamline
		CalLocalTracing(fix_triangleID, theFixPoint[0], theFixPoint[1], type);

		/*---------------------------------------------------------------*/
		////Store the information of the new limit cycle

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
		BuildHandlerforLimitCycle(cur_limitcycle_index);

		cur_traj_index ++;
		cur_limitcycle_index++;

		//calculate the next level beginning point
		GetNextLevelBeginPoint(singID, begin_p, type);
		//the next extracted limit cycle may be an outer limit cycle, but it is not always correct
		//we need to judge wether the new limit cycle contain the center singularity!!!!!
		inner ++;  

		//reverse the type of the limit cycle detection
		if(type == 0) type = 1;
		else type = 0;
		goto L2;

		if(Saddleornot == 1 && index_separatrix < 3)
			goto L1;

	}
}




//Main routine

void LimitCycleDetection()
{
	//declare local variables
    int i, Saddleornot, type;

	//Initialization part
    cur_limitcycle_index = 0;

	for(i = 0; i < cur_singularity_index; i++)
	{
		Saddleornot = 0;

		if(singularities[i].type == CWCENTER || singularities[i].type == CCWCENTER)
			continue;

		if(singularities[i].type == SOURCE || singularities[i].type == RFOCUS)//it is a repeller
			type = 0;

		else if(singularities[i].type == SINK || singularities[i].type == AFOCUS) //it is an attractor
			type = 1;

		else   //it is a saddle
		{
			type = 0;
			Saddleornot = 1;
		}

		DetectFromASingularity(i, type, Saddleornot);
	}

	if(cur_limitcycle_index == 0)
	{
		num_linesegs_curtraj[cur_traj_index] = 0;
		MessageBox(NULL, "No limit cycle has been found!", "error", MB_OK);
	}
}






////////////////////////////////////////////////////////////////////////////////////
/**********************************************************************************/

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

extern void IntersectRegion(TriangularRegion &source1, TriangularRegion &source2, TriangularRegion &dest);

bool IntheRegion(int triangle, int type)
{
	int num_regtriangles;
	int *trianglelist;

	if(type == 0)
	{
		trianglelist = repellerRegion.trianglelist;
		num_regtriangles = repellerRegion.num;
	}
	else{
		trianglelist = attractorRegion.trianglelist;
		num_regtriangles = attractorRegion.num;
	}

	if(IsRepeated(trianglelist, triangle, num_regtriangles))
		return true;
	else 
		return false;
}


void GetTheIntersectRegion()
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
}



bool IsGoodSharedEdge(Edge *cur_e)
{
	Vertex *v1, *v2;
	icVector2 vec1, vec2, edge_vec;
	v1 = Object.vlist[cur_e->verts[0]];
	v2 = Object.vlist[cur_e->verts[1]];

	vec1 = v1->vec;
	vec2 = v2->vec;
	edge_vec.entry[0] = v2->x - v1->x;
	edge_vec.entry[1] = v2->y - v1->y;

	normalize(vec1);
	normalize(vec2);
	normalize(edge_vec);

	if(dot(vec1, vec2) < 0.5) return false;

	if(fabs(dot(vec1, edge_vec)) > 0.5 || fabs(dot(vec2, edge_vec)) > 0.5)
		return false;

	return true;
}





void GetAShareEdge_new()
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
				if(IsGoodSharedEdge(cur_e) || i == 0)
				    OneSharedEdge = cur_e;             ////recalled the first shared edge inside cell cycle

				continue;
			}
		}
	}

}


extern int *cellcycle;
extern int num_celltriangle;                ////number of triangles in the triangle list


////Here begin_triangle is the triangle on the boundary
////begin_p[] is the center of the triangle
void DetectALimitCycle(double begin_p[2], int begin_triangle, int type, int inner)
{
	int flag = -1;
	int fix_triangleID;

	//double center[2] = {}

	Double_CellCycleDetect_2(begin_p[0], begin_p[1], begin_triangle, type, flag); 
	//First_CellCycleDetect_2(begin_p[0], begin_p[1], begin_triangle, cellcycle, num_celltriangle, type, flag);
	
	if(flag == 2) //reach a singularity or the boundary of the mesh
	{
		return;
	}

	if(flag == 1) //we found a cell cycle
	{
		int pre_limit_index = -1;

		if(num_celltriangle == 0)
			return;

		//if it is a previously detected limit cycle
		if(IsPreviousLimitCycle(pre_limit_index))
		{
			//UpdateListInSingularity(singID, pre_limit_index);
			//UpdateListInLimitCycle(pre_limit_index, singID);

			//if(Saddleornot == 1 && index_separatrix < 3)
			//	goto L1;

			return;
		}

		/*---------------------------------------------------------------*/
		//if it is a new extracted limit cycle
		boundaryBuilding();
       
		//GetASharedEdge(center);
        GetAShareEdge_new();

		StoreCurrentCellCycle(cellcycle, num_celltriangle);  
		Cur_CellCycleIndex++;   //do we have multiple cell cycle during one limit cycle detection process

		//Get the closed orbit, this is the unstable part of the whole detection
		GettheClosedStreamline_new(type, theFixPoint, fix_triangleID);

		////since some times, the program can not find the fixed point, you need to deal with 
		////the beginning point is origin 1/5/06
		if((theFixPoint[0] == 0 && theFixPoint[1] == 0)||fix_triangleID < 0)
		{
			fix_triangleID = TriangleDetect(theFixPoint[0], theFixPoint[0]);
		}
			

		//calculate the closed streamline
		CalLocalTracing(fix_triangleID, theFixPoint[0], theFixPoint[1], type);

		/*---------------------------------------------------------------*/
		////Store the information of the new limit cycle

		//limitcycles[cur_limitcycle_index].singularID = singID;  ////store the center singularity of current limit cycle
		
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
		
		//if(inner == 0) //if it is the most inner limit cycle
		//{
		//	UpdateListInSingularity(singID, cur_limitcycle_index);
		//	UpdateListInLimitCycle(cur_limitcycle_index, singID);
		//}

		//else //if it is an outer limit cycle
		//{
		//	UpdateCycleListInLimitCycle(cur_limitcycle_index-1, cur_limitcycle_index);
		//	UpdateCycleListInLimitCycle(cur_limitcycle_index, cur_limitcycle_index-1);
		//}
		/*---------------------------------------------------------------*/

		////Build the graphical handler for the limit cycle for user interaction
		BuildHandlerforLimitCycle(cur_limitcycle_index);

		cur_traj_index ++;
		cur_limitcycle_index++;
	}
}





void GetCenterofATriangle(int triangle, double center[2])
{
	Face *other_f = Object.flist[triangle];
	Vertex *v1 = Object.vlist[other_f->verts[0]];
	Vertex *v2 = Object.vlist[other_f->verts[1]];
	Vertex *v3 = Object.vlist[other_f->verts[2]];

	center[0] = 0.3*v1->x + 0.3*v2->x + 0.4*v3->x;
	center[1] = 0.3*v1->y + 0.3*v2->y + 0.4*v3->y;

}


bool IsAtFieldBoundary(int triangle)
{
	Corner *c;
	int i;

	for(i = 0; i < 3; i++)
	{
		c = Object.clist[3*triangle + i];

		if(c->t == -1)
			return true;
	}

	return false;
}



int test_aboundtriangle;
int second_begintriangle;

void GetACellCycle_new(int singID, int type)
{
	int i;
	Edge **temp_edgelist;
	int num_boundedges = 0;
    Edge *oneboundedge;
	int aboundtriangle = singularities[singID].Triangle_ID;
	
	double begin_p[2];
	int inner = 0;

	int test_count = 0;

	repellerRegion.num = 0;
	attractorRegion.num = 0;

	////first, grow a region from the triangle containing the center singularity
LL:	if(type == 0)  
	{
		//the type of the singularity at the center is a repeller, hence the limit cycle 
		//should be an attractor
		repellerRegion.trianglelist[repellerRegion.num] = aboundtriangle;
		repellerRegion.num ++;

		UpdateBoundary(0);
		GetRegionNormals(0);
		Cancel_Growing(0, -1);

		temp_edgelist = repellerBoundary.edgelist;
		num_boundedges = repellerBoundary.num;
		oneboundedge = repellerBoundary.edgelist[0];
		//oneboundedge = repellerBoundary.edgelist[num_boundedges-1];

		////Copy the region
		CopyRegion(repellerRegion.trianglelist, Source_re.trianglelist, repellerRegion.num);
		Source_re.num = repellerRegion.num;
	}

	else
	{
		attractorRegion.trianglelist[attractorRegion.num] = aboundtriangle;
		attractorRegion.num ++;

		UpdateBoundary(1);
		GetRegionNormals(1);
		Cancel_Growing(1, -1);

		temp_edgelist = attractorBoundary.edgelist;
		num_boundedges = attractorBoundary.num;
		oneboundedge = attractorBoundary.edgelist[0];
		//oneboundedge = attractorBoundary.edgelist[num_boundedges-1];
		
		////Copy the region
		CopyRegion(attractorRegion.trianglelist, Sink_re.trianglelist, repellerRegion.num);
		Sink_re.num = repellerRegion.num;
	}

	////pick one triangle at the boundary of the region obtained above

	////A important problem! for embedded limit cycle detection, we can not use random triangle 
	////on the boundary, here we need to use the triangle on the outer boundary!!!
	if(IntheRegion(oneboundedge->tris[0], type))
		aboundtriangle = oneboundedge->tris[0];
	else
		aboundtriangle = oneboundedge->tris[1];
	
	////grow the backward region

	//if(type == 0)
	//{
	//	repellerRegion.num = 0;
	//	repellerRegion.trianglelist[0] = aboundtriangle;
	//	repellerRegion.num = 1;

	//	UpdateBoundary(0);
	//	GetRegionNormals(0);
	//	Cancel_Growing(0, -1);

	//	////calculate the intersection region
	//	IntersectRegion(repellerRegion, Source_re, intersectRegion);
	//}

	//else
	//{
	//	attractorRegion.num = 0;
	//	attractorRegion.trianglelist[0] = aboundtriangle;
	//	attractorRegion.num = 1;

	//	UpdateBoundary(1);
	//	GetRegionNormals(1);
	//	Cancel_Growing(1, -1);
	//	
	//	////calculate the intersection region
	//	IntersectRegion(attractorRegion, Sink_re, intersectRegion);
	//}

	////then we can perform local tracing from the outer boundary of the region to locate the cell cycle!
	////Now we can judge the tangent case easily, since we have already know the number of the
	////cell triangles roughly!!! so if the number of detected cell triangles is not equal to the rough 
	////number of cell triangles, we think that maybe a numerical error 1/11/06

	/*-------------------------------------------------------*/
	////Here just a testing code
	test_aboundtriangle = aboundtriangle;

	if(test_count == 0)
		second_begintriangle = aboundtriangle;
	/*-------------------------------------------------------*/

	////Locate the real cell cycle and closed orbit
	GetCenterofATriangle(aboundtriangle, begin_p);
	DetectALimitCycle(begin_p, aboundtriangle, type, inner);


	////For embedded limit cycles, we should begin from the previous selected triangle and perform 
	////backward region growing!!(which means repeller growing will become attractor growing)

	////first, we need to reverse the type
	////if the selected boundary triangle is not at the boundary of the whole field
	if(IsAtFieldBoundary(aboundtriangle) || test_count == 1)
		return;

	////Interlace copy here !!! 1/11/06
	if(type == 0)
	{
		//CopyRegion(Source_re.trianglelist, attractorRegion.trianglelist, Source_re.num);
		//attractorRegion.num = Source_re.num;
		CopyRegion(repellerRegion.trianglelist, attractorRegion.trianglelist, repellerRegion.num);
		attractorRegion.num = repellerRegion.num;
	}

	else
	{
		//CopyRegion(Sink_re.trianglelist, repellerRegion.trianglelist, Sink_re.num);
		//repellerRegion.num = Sink_re.num;
		CopyRegion(attractorRegion.trianglelist, repellerRegion.trianglelist, attractorRegion.num);
		repellerRegion.num = attractorRegion.num;
	}

	type = 1 - type;

	test_count ++;

	goto LL;
}


void DetectLimitCycle_new()
{
	//declare local variables
    int i, Saddleornot, type;

	//Initialization part
    cur_limitcycle_index = 0;

	for(i = 0; i < cur_singularity_index; i++)
	{
		Saddleornot = 0;

		if(singularities[i].type == CWCENTER || singularities[i].type == CCWCENTER)
			continue;

		if(singularities[i].type == SOURCE || singularities[i].type == RFOCUS)//it is a repeller
			type = 0;

		else if(singularities[i].type == SINK || singularities[i].type == AFOCUS) //it is an attractor
			type = 1;

		else   //it is a saddle
		{
			type = 0;
			Saddleornot = 1;
		}

		GetACellCycle_new(i, type);
	}

	//if(cur_limitcycle_index == 0)
	//{
	//	num_linesegs_curtraj[cur_traj_index] = 0;
	//	MessageBox(NULL, "No limit cycle has been found!", "error", MB_OK);
	//}
}