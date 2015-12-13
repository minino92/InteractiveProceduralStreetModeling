////LimitCycleEdit.cpp

#include "stdafx.h"

#include "LimitCycleEdit.h"

#include "VFDataStructure.h"
#include "VFSynthesis.h"

#include "LimitCycleDetect.h"

#include "LimitCycleCreator.h"

#include "Regionsmoothing.h"

#include "topologyedit.h"

#include "LocalTracing.h"

#include "ConleyRelationGraph.h"

#include "PairCancellationModule.h"

extern Polygon3D Object;

extern Singularities *singularities;           //being captured singularites' list
extern int cur_singularity_index;
extern LineSeg **trajectories;                 //trajectories' list
extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory

extern Separatrices *separatrices;             //array for group of separatrices

extern LimitCycle *limitcycles;                 //limit cycles data structure
extern int cur_limitcycle_index;
extern int MaxNumLimitCycles;

////variable from region smoothing
extern Vertex **regionverts;                  ////vertices inside user selected region
extern int Num_verts;                         ////number of inner vertices
extern int MaxNumVerts;

extern Point *point;                       ////we may initial it as 50 points, over 50, we can extend it
extern int Num_SmoothRegionpoints;                     ////Number of points that user selected
extern int MaxNumPoints;                   ////Maximum element in the point array, we can extend it if needed


////Variables for finding the Intermediary saddles
extern int *MediaNodes;
extern int Num_MediaNodes;
extern int MaxMediaNodes;

extern GraphNode *graphnodes;


////we may use the region and boundary variables of the topology edit modula
extern TriangularRegion repellerRegion;       ////region containing a repeller
extern TriangularRegion attractorRegion;      ////region containing an attractor
extern RegionBoundary repellerBoundary;
extern RegionBoundary attractorBoundary;
extern InnerVertices repellerInnerverts;
extern InnerVertices attractorInnerverts;

extern TriangularRegion intersectRegion;     ////The intersect region
extern RegionBoundary intersectBoundary;     ////The intersect boundary 
extern InnerVertices intersectInnerverts;    ////The inner vertices inside the intersect region

extern int MaxNumInnerVerts;

////Extern routines
extern void AddToInnerVerts(Vertex *, int);
extern int GetOppositeTriangle(Edge *, int);
extern bool AttractorExitEdgePending(Edge *);
extern bool RepellerExitEdgePending(Edge *);
extern void UnionRegion(TriangularRegion &source1, TriangularRegion &source2, TriangularRegion &dest);
extern void IntersectRegion(TriangularRegion &source1, TriangularRegion &source2, TriangularRegion &dest);


////Testing variables 07/24/05
extern int **CellCycleList;
extern int Cur_CellCycleIndex;                    ////current cell cycle index
extern int *NumTriangleInEachCellCycle;           ////store the number of triangles in each cell cycle


////The external variables are introduced for limit cycle relocation 3/8/06
extern int resolution;    // how many points in our output array
extern ctr_point *out_pts;
extern ctr_point *control_pts;        // allocate our control point array
extern int num_shapecontrol_pts;

extern int *DesignCurveCellCycle;
extern int num_triangles_designcurve;
extern int MaxNumTrianglesDesignCurve;

extern LineSeg *designcurve;
extern int num_lineseg_designcurve;


extern Edge **Cycle_edge;
extern int num_cycleedges;
extern int num_curvepts_output;
extern int num_innertriangles;

extern int *Boundaryverts;
extern int num_boundverts;
extern DesignTriangleCycle myCycle;



extern void SetBoundaryFlag_Ver(Edge **boundaryedgelist, int num_edges);


/***************************************************************/
////Routines for limit cycle cancellation

bool InVertsList(Vertex **verts, int num, Vertex *vert_index)
{
	int i;
	for(i = 0; i < num; i++)
	{
		if(verts[i] == vert_index)
			return true;
	}
	return false;
}

void GetPointsOnStreamline(int limitcycle_index, int num_pts)
{
	////try to avoid selecting a vertex !!!!
	int interval = (int)limitcycles[limitcycle_index].num_linesegs/num_pts;

	int i;
	for(i = 0; i < num_pts; i++)
	{
	//if(not a vertex)
		point[i].x = limitcycles[limitcycle_index].closed_streamline[i*interval].gstart[0];
		point[i].y = limitcycles[limitcycle_index].closed_streamline[i*interval].gstart[1];
	}

	Num_SmoothRegionpoints = num_pts;
}

////Get the vertex that most close to the center of the limit cycle
////But how about the limit cycle without any singularity inside it??????
int GetBeginningVert(int limitcycle_index)
{
	//int vert_index;

	int i, j;
	int num_pts = 4;
	Face *face;
	Vertex *vert;

	////pick at least four points on the closed streamline of the limit cycle
	////construct a region according to the order of the points, perform the similar first vertex finding
	////as we perform in region smoothing
LL:	GetPointsOnStreamline(limitcycle_index, num_pts);   ////pick 4 points intially

	////Need to optimize the following algorithm, we may access a vertex many times!!!!! 
	for(i = 0; i < limitcycles[limitcycle_index].num_triangles; i++)
	{
		face = Object.flist[limitcycles[limitcycle_index].cellcycle[i]];

		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];

			if(InRegion(vert->VertID))
				return (vert->VertID);
		}
	}

	////if we find a vertex, return the vertex
	////otherwise, we need to change the region(increase the number of points) to perform searching again
	if(num_pts >= limitcycles[limitcycle_index].num_triangles)
	{
		MessageBox(NULL, "Can not find a valid inner vertex", "Error", MB_OK);
		return -1;
	}

	num_pts++;
	goto LL;


}

bool GetInnerLimitCycleRegion(int limitcycle_index, int *trianglelist, int num)
{
	int begin_vert = -1;
	int i, j;
	Face *face;
	Vertex *vert;

	////Add all the vertices in the cell cycle into the inner vertices list
	Num_verts = 0;
	for(i = 0; i < num; i++)
	{
		face = Object.flist[trianglelist[i]];

		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];

			////If the vertex is not in the list, add it to the list
			if(!InVertsList(regionverts, Num_verts, vert))
			{
				if(Num_verts >= MaxNumVerts-1)
				{
					MaxNumVerts += 100;
					regionverts = (Vertex**)realloc(regionverts, sizeof(Vertex*) * MaxNumVerts);
					
				}
				regionverts[Num_verts] = vert;
				vert->RegionListID = Num_verts;
				vert->InRegion = 1;
				vert->OnBoundary = 0;
				Num_verts ++;
			}

		}
	}

	////Find the true inner vertices!! (Diffcult part)
	begin_vert = GetBeginningVert(limitcycle_index);  ////Testing the first limit cycle

	if(begin_vert < 0)
		return false;

	////Call the same inner vertices searching routine of region smoothing to find all the other inner vertices
	InnerLimitCycleSearch(begin_vert);

	return true;
}


////The following will be somewhat different from the original one in region smoothing
void InnerLimitCycleSearch(int seed)
{
	int i, j;

	//regionverts[Num_verts] = Object.vlist[seed];
	//Num_verts++;

	Vertex *vert; 
	Vertex *adj_v;

	Edge *temp_e;

	i = Num_verts;

	////Perform searching on seed vertex here
	vert = Object.vlist[seed];

	for(j = 0; j < vert->Num_edge; j++)
	{
		temp_e = vert->edges[j];

		if(Object.vlist[temp_e->verts[0]] != vert)
			adj_v = Object.vlist[temp_e->verts[0]];
		else
			adj_v = Object.vlist[temp_e->verts[1]];

		////if the vertex is not on boundary and has not been marked, add it to the list
		if(adj_v->InRegion != 1 && adj_v->OnBoundary != 1)
		{
			if(Num_verts >= MaxNumVerts - 1)
			{
				MaxNumVerts += 100;
				regionverts = (Vertex**)realloc(regionverts, sizeof(Vertex*) * MaxNumVerts);
				
			}

			regionverts[Num_verts] = adj_v;
			regionverts[Num_verts]->RegionListID = Num_verts;
			adj_v->InRegion = 1;
			adj_v->OnBoundary = 0;
			Num_verts ++;
		}
	}



	////search all the other inner vertices here
	while(i < Num_verts)
	{
		vert = regionverts[i];  ////get current source vertex
		 
		////search all its adjacent vertices and judge whether it should be added to the inner vertices list or not
		for(j = 0; j < vert->Num_edge; j++)
		{
			temp_e = vert->edges[j];

			if(Object.vlist[temp_e->verts[0]] != vert)
				adj_v = Object.vlist[temp_e->verts[0]];
			else
				adj_v = Object.vlist[temp_e->verts[1]];

			////if the vertex is not on boundary and has not been marked, add it to the list
			if(adj_v->InRegion != 1 && adj_v->OnBoundary != 1)
			{
				if(Num_verts >= MaxNumVerts - 1)
				{
					MaxNumVerts += 100;
					regionverts = (Vertex**)realloc(regionverts, sizeof(Vertex*) * MaxNumVerts);
					
				}

				regionverts[Num_verts] = adj_v;
				regionverts[Num_verts]->RegionListID = Num_verts;
				adj_v->InRegion = 1;
				adj_v->OnBoundary = 0;
				Num_verts ++;
			}
		}

		i++;
	}
}


/////testing routine for limit cycle removement/cancellation (a limit cycle + a singularity)
void CancelLimitCycle(int index)
{
	GetInnerLimitCycleRegion(index, limitcycles[index].cellcycle, limitcycles[index].num_triangles);

	SavePreviousField();
	RegionSmooth();
	SavePostField();
}


////Note that  we can use limit cycle region growing (Conley region) to get a better region for smoothing
void CancelOneLimitCycle(int index)
{
	////Perform region growing
	GetLimitCycleRegion(index, limitcycles[index].type, 0);

	///Get inner vertices
    SavePreviousField();
	RegionSmooth();
	SavePostField();
}


/*---------------------------------------------------------------------*/
////Saddle and singularity connection judgement
bool SaddleSingConnectionJudge(int singID, int saddleID)
{
	int i;
	double globalp[2];
	int pre_face, cur_face;
	int type = -1;
	int flag = 0;

	int loop_flag = 0;  ////for detect a closed loop from the beginning triangle
	int *triangles = (int*)malloc(sizeof(int) * NUMTRACINGTRIANGLE);
	int num_triangles = 0;

	////
	icVector2 dist;

	////We need to calculate the first point using the similar way as calculating separatrices

	////Positive direction
	if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
	{
		globalp[0] = singularities[saddleID].gcx + SEPARATRIXSTEP * singularities[saddleID].incoming.entry[0];   
		globalp[1] = singularities[saddleID].gcy + SEPARATRIXSTEP * singularities[saddleID].incoming.entry[1];
		type = 1;   ////perform backward tracing
	}
	else if(singularities[singID].type == SINK || singularities[singID].type == AFOCUS)
	{
		globalp[0] = singularities[saddleID].gcx + SEPARATRIXSTEP * singularities[saddleID].outgoing.entry[0];   
		globalp[1] = singularities[saddleID].gcy + SEPARATRIXSTEP * singularities[saddleID].outgoing.entry[1];
		type = 0;   ////perform forward tracing
	}
	pre_face = cur_face = TriangleDetect(globalp[0], globalp[1]);


	for(i = 0; i < NUMTRACINGTRIANGLE; i++)
	{

		if(cur_face <= -1)
		{
			break;
		}

		pre_face = cur_face;
		cur_face = TraceInATriangle(cur_face, globalp, type, flag); //

		if(flag == 3)
		{
			//if(cur_face == singularities[singID].Triangle_ID)
			int lineseg_index = num_linesegs_curtraj[cur_traj_index]-1;
			dist.entry[0] = trajectories[cur_traj_index][lineseg_index].gend[0] - singularities[singID].gcx;
			dist.entry[1] = trajectories[cur_traj_index][lineseg_index].gend[1] - singularities[singID].gcy;

			if(length(dist) < 1e-3)
			{
				free(triangles);
				return true;
			}
			else
				break;
		}

		if(flag == 4 /*|| pre_face == cur_face || loop_flag == 1*/) 
		{
			break;
		}

		////Not accurate to use triangle loop detection only
		////after we find a triangle loop, we need to test whether the curve constitutues a "closed" orbit!!!
		////we can use a small neighbor hood to measure the distance !

		//if(LoopDetect(triangles, num_triangles, cur_face))
		//	loop_flag = 1;
		//else{
		//	triangles[num_triangles] = cur_face;
		//	num_triangles++;
		//}
	}


	////Negative direction
	flag = 0;
	num_triangles = 0;
	if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
	{
		globalp[0] = singularities[saddleID].gcx - SEPARATRIXSTEP * singularities[saddleID].incoming.entry[0];   
		globalp[1] = singularities[saddleID].gcy - SEPARATRIXSTEP * singularities[saddleID].incoming.entry[1];
		type = 1;  ////perform backward tracing
	}
	else if(singularities[singID].type == SINK || singularities[singID].type == AFOCUS)
	{
		globalp[0] = singularities[saddleID].gcx - SEPARATRIXSTEP * singularities[saddleID].outgoing.entry[0];   
		globalp[1] = singularities[saddleID].gcy - SEPARATRIXSTEP * singularities[saddleID].outgoing.entry[1];
		type = 0;  ////perform forward tracing
	}
	pre_face = cur_face = TriangleDetect(globalp[0], globalp[1]);


	for(i = 0; i < NUMTRACINGTRIANGLE; i++)
	{

		if(cur_face <= -1)
		{
			free(triangles);
			return false;
		}

		pre_face = cur_face;
		cur_face = TraceInATriangle(cur_face, globalp, type, flag); //

		if(flag == 3)
		{
			//if(cur_face == singularities[singID].Triangle_ID)
			int lineseg_index = num_linesegs_curtraj[cur_traj_index]-1;
			dist.entry[0] = trajectories[cur_traj_index][lineseg_index].gend[0] - singularities[singID].gcx;
			dist.entry[1] = trajectories[cur_traj_index][lineseg_index].gend[1] - singularities[singID].gcy;

			if(length(dist) < 1e-3)
			{
				free(triangles);
				return true;
			}
			else
			{
				free(triangles);
				return false;
			}
		}

		if(flag == 4 /*|| pre_face == cur_face || loop_flag == 1*/) 
		{
			free(triangles);
			return false;
		}

		////Not accurate to use triangle loop detection only
		////after we find a triangle loop, we need to test whether the curve constitutues a "closed" orbit!!!
		////we can use a small neighbor hood to measure the distance !

		//if(LoopDetect(triangles, num_triangles, cur_face))
		//	loop_flag = 1;
		//else{
		//	triangles[num_triangles] = cur_face;
		//	num_triangles++;
		//}
	}

	free(triangles);
	return false;
}

////Find the saddles that have connecting trajectory with specific singularity
////we suppose only two saddles adjacent to the current singularity
////Modified on 09/06/05
bool FindTheConnectedSaddle(int singID, int &saddle1, int &saddle2)
{
	int i;

	int count = 0;

	////Suppose we don't know which saddles are inside the limit cycle
	for(i = 0; i < cur_singularity_index; i++)
	{
		if(singularities[i].type != SADDLE)
			continue;
		////if this is one of the saddle we want
		if(SaddleSingConnectionJudge(singID, i))  ////Here we store the two adjacent saddles at the same time
		{
			if(count == 0)
			{
				saddle1 = i;
				count ++;
			}
			else
				saddle2 = i;
		}
	}

	if(count == 0)
		return false;
	else
		return true;
}


////Grow initial triangular strip for the limit cycle
////Modified on 09/06/05

//void GrowInitLimitRegion(int SaddleID, int otherSaddleID, int type)
//{
//	////These two saddles should be included into the smoothing region according to Conley index theory
//	if(type == 0)
//	{
//		if(otherSaddleID >= 0)
//		{
//			repellerRegion.trianglelist[0] = singularities[SaddleID].Triangle_ID;
//			repellerRegion.trianglelist[1] = singularities[otherSaddleID].Triangle_ID;
//			repellerRegion.num = 2;
//			Object.flist[repellerRegion.trianglelist[0]]->repell_inregion = 1;
//			Object.flist[repellerRegion.trianglelist[0]]->attract_inregion = 1;
//			Object.flist[repellerRegion.trianglelist[1]]->repell_inregion = 1;
//			Object.flist[repellerRegion.trianglelist[1]]->attract_inregion = 1;
//		}
//		else{
//			repellerRegion.trianglelist[0] = singularities[SaddleID].Triangle_ID;
//			repellerRegion.num = 1;
//			Object.flist[repellerRegion.trianglelist[0]]->repell_inregion = 1;
//			Object.flist[repellerRegion.trianglelist[0]]->attract_inregion = 1;
//		}
//		Cancel_InitSaddleRegion(SaddleID, 1, InitSaddleRegionLength);
//	}
//	else{
//		if(otherSaddleID >= 0)
//		{
//			attractorRegion.trianglelist[0] = singularities[SaddleID].Triangle_ID;
//			attractorRegion.trianglelist[1] = singularities[otherSaddleID].Triangle_ID;
//			attractorRegion.num = 2;
//			Object.flist[attractorRegion.trianglelist[0]]->attract_inregion = 1;
//			Object.flist[attractorRegion.trianglelist[0]]->repell_inregion = 1;
//			Object.flist[attractorRegion.trianglelist[1]]->attract_inregion = 1;
//			Object.flist[attractorRegion.trianglelist[1]]->repell_inregion = 1;
//		}
//		else{
//			attractorRegion.trianglelist[0] = singularities[SaddleID].Triangle_ID;
//			attractorRegion.num = 1;
//			Object.flist[attractorRegion.trianglelist[0]]->attract_inregion = 1;
//			Object.flist[attractorRegion.trianglelist[0]]->repell_inregion = 1;
//		}
//		Cancel_InitSaddleRegion(SaddleID, 0, InitSaddleRegionLength);
//	}
//
//	UpdateBoundary(type);
//	GetRegionNormals(type);
//}
//
//
//bool GetLimitRegionFromSaddle(int singID, int type)
//{
//	int theSaddle1, theSaddle2;
//	theSaddle1 = theSaddle2 = -1;
//
//	////1. Find the saddle that has connecting trajectory between it and the other singularity
//	if(!FindTheConnectedSaddle(singID, theSaddle1, theSaddle2))
//		return false;
//
//
//	////2. Begin from the saddle captured above, perform initial triangular strip growing along other eigen vector
//	////direction, until it reaches other singularity, whole field boundary or form a cell cycle
//	GrowInitLimitRegion(theSaddle1, theSaddle2, type);   
//
//	////3. perform Conley region growing as usual for the limit cycle
//	Region_Growing(type);
//
//	return true;
//}


////Important notes (10/15/05)
////The region of the limit cycle should grow from the cell cycle of the limit cycle not from the saddle!!!!!
////The reigon of the singularity should grow from the singularity itself
////The region and the singularities that may be covered should follow the Conley index theory
////(We can get the information from the Conley relation/connection graph)

////Get the region for performing smoothing to cancel a limit cycle and a singularity inside it
////The variable "limitID" may be used to judge whether the singularities are inside the limit cycle
//void GetLimitSingSmoothingRegion(int singID, int limitID)
//{
//	int type;
//	int saddleID;
//
//	
//	////we may first test whether the singularity is inside the limit cycle or not
//	////we can use the inner boundary edges to build the judging region, but we need to sort them! 
//
//	////if current selected saddle is proved to be inside the limit cycle
//	////then, we can make sure that there are more than one singularities inside currrent limit cycle
//	
//	////A saddle can be used to cancel with a limit cycle!
//	if(singularities[singID].type == SADDLE)
//	{
//		////call another routine
//		GetLimitSaddlePairRegion(singID, limitID);
//		return;
//	}
//
//	if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
//		type = 0;
//	else
//		type = 1;
//
//	////We do not consider the same type cancellation now 12/29/05
//	if(type == limitcycles[limitID].type)
//	{
//		return;  
//	}
//
//
//	////1. grow region for the singularity
//    if(type == 0)
//	{
//		Cancel_GrowRepellerRegion(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength);
//	}
//	else{
//		Cancel_GrowAttractorRegion(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength);
//	}
//
//	////2. grow region for limit cycle
//	//GetInnerLimitCycleRegion(limitID, limitcycles[limitID].cellcycle, limitcycles[limitID].num_triangles);
//	GrowRegionforALimitCycle(limitID);
//
//
//	////3. Grow regions from all the intermediary singularities between them  12/29/05
//	//This should be based on the Conley connection/topology graph!
//	////Note we need to add the triangles containing the two adjacent saddles into the region!!
//	if(!GetLimitRegionFromSaddle(singID, 1-type))
//	{
//		MessageBox(NULL, "Can not find the connected saddle!", "Error", MB_OK);
//        return;
//	}
//
//	////4. Calculate the intersected region
//	////We need a new get intersect routine here
//	GetLimitSingIntersectRegion(singID/*, saddleID*/);
//}
//
//
////2/07/06
//void GetLimitSingSmoothingRegion2(int singID, int limitID, int *MediaNodes, int Num_medianodes)
//{
//	if(Num_medianodes == 0) //we suppose this case is the limit cycle and center singularity cancellation
//	{
//		if(limitcycles[limitID].type == 0) // limit cycle is a repeller, so the center is a attractor
//		{
//			Cancel_GrowAttractorRegion(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength);
//			CopyRegion(attractorRegion.trianglelist, intersectRegion.trianglelist, attractorRegion.num);
//			intersectRegion.num = attractorRegion.num;
//		}
//
//		else // limit cycle is an attractor, so the center is a repeller
//		{
//		    Cancel_GrowRepellerRegion(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength);
//			CopyRegion(repellerRegion.trianglelist, intersectRegion.trianglelist, repellerRegion.num);
//			intersectRegion.num = repellerRegion.num;
//		}
//		return;
//	}
//
//	////otherwise, there are some intervals (saddles) between them
//
//	////similarly as the multiple singularities cancellation, we perform region growing from saddles
//	////as repeller and attractor respectively, and combine them with the region of limit cycle!
//	int i;
//	int cur_saddle;
//	TriangularRegion Temp_re;
//	Temp_re.trianglelist = (int*)malloc(sizeof(int)*Object.nfaces);
//	Temp_re.num = 0;
//
//
//	if(limitcycles[limitID].type == 0) //limit cycle is a repeller
//	{
//		//1.1) Grow region for limit cycle(Note that it seems that we need not grow region from limit cycle!)
//
//		////Grow region for the singularity
//		Cancel_GrowAttractorRegion(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength);
//		UnionRegion(attractorRegion, intersectRegion, Temp_re);
//		CopyRegion(Temp_re.trianglelist, intersectRegion.trianglelist, Temp_re.num);
//		intersectRegion.num = Temp_re.num;
//
//		//Grow region from each saddle as attractors
//
//		for(i = 0; i < Num_medianodes; i++)
//		{
//			//after growing a region for one saddle, we need to union it with previous saddle regions
//			attractorRegion.num = 0;
//
//			cur_saddle = graphnodes[MediaNodes[i]].singularityID;
//			GetLimitSaddlePairRegion2(cur_saddle, limitID, singID);
//			AddToRegionTriangles(singularities[singID].Triangle_ID, 1);
//
//			UnionRegion(attractorRegion, intersectRegion, Temp_re);
//			CopyRegion(Temp_re.trianglelist, intersectRegion.trianglelist, Temp_re.num);
//			intersectRegion.num = Temp_re.num;
//
//		}
//
//		//The final region is the region for smoothing
//		//Make sure that the final region cover the singularity
//	}
//
//	else  //limit cycle is an attractor
//	{
//		////Grow region for the singularity
//		Cancel_GrowRepellerRegion(singularities[singID].Triangle_ID, -1, singID, InitSaddleRegionLength);
//		UnionRegion(repellerRegion, intersectRegion, Temp_re);
//		CopyRegion(Temp_re.trianglelist, intersectRegion.trianglelist, Temp_re.num);
//		intersectRegion.num = Temp_re.num;
//
//		//Grow region from each saddle as repellers
//
//		for(i = 0; i < Num_medianodes; i++)
//		{
//			//after growing a region for one saddle, we need to union it with previous saddle regions
//			repellerRegion.num = 0;
//
//			cur_saddle = graphnodes[MediaNodes[i]].singularityID;
//			GetLimitSaddlePairRegion2(cur_saddle, limitID, singID);
//			AddToRegionTriangles(singularities[singID].Triangle_ID, 0);
//
//			UnionRegion(repellerRegion, intersectRegion, Temp_re);
//			CopyRegion(Temp_re.trianglelist, intersectRegion.trianglelist, Temp_re.num);
//			intersectRegion.num = Temp_re.num;
//		}
//
//		//The final region is the region for smoothing
//		//Make sure that the final region cover the singularity
//	}
//
//	free(Temp_re.trianglelist);
//}


////Get the region for the limit cycle and a saddle cancellation
void GetLimitSaddlePairRegion(int saddleID, int limitID)
{
	int i, j;
	int type = limitcycles[limitID].type;
	Edge **edgelist;
	int edgesnum;
	Vertex **vertslist;
	int vertsnum;
	int *trianglelist;
	int trianglesnum;
	Vertex *cur_v;
	Face *face;

	if(type == 0) ////the saddle acts as an attractor
	{
		attractorRegion.trianglelist[0] = singularities[saddleID].Triangle_ID;
		Object.flist[singularities[saddleID].Triangle_ID]->attract_inregion = 1;
		attractorRegion.num = 1;
	}
	else{
		repellerRegion.trianglelist[0] = singularities[saddleID].Triangle_ID;
		Object.flist[singularities[saddleID].Triangle_ID]->repell_inregion = 1;
		repellerRegion.num = 1;
	}

	//Cancel_InitLimitSaddleRegion(saddleID, type, 200);
    Cancel_InitLimitSaddleRegion2(saddleID, type, 200);
    UpdateBoundary(1-type);
	GetRegionNormals(1-type);
	Region_Growing(1-type);

	UpdateBoundary(1-type);

	if(type == 0)
	{
		edgelist = attractorBoundary.edgelist;
		vertslist = attractorInnerverts.vertslist;
		edgesnum = attractorBoundary.num;
		vertsnum = attractorInnerverts.num;
		
		trianglelist = attractorRegion.trianglelist;
		trianglesnum = attractorRegion.num;
	}
	else{
		edgelist =  repellerBoundary.edgelist;
		vertslist = repellerInnerverts.vertslist;
		edgesnum = repellerBoundary.num;
		vertsnum = repellerInnerverts.num;
		
		trianglelist = repellerRegion.trianglelist;
		trianglesnum = repellerRegion.num;
	}

	////Get the all the inner vertices
	////Set the vertices on the boundary as 'OnBoundary'
	for(i = 0; i < edgesnum; i++)
	{
		Object.vlist[edgelist[i]->verts[0]]->OnBoundary = 1;
		Object.vlist[edgelist[i]->verts[0]]->InRegion = 0;
		Object.vlist[edgelist[i]->verts[1]]->OnBoundary = 1;
		Object.vlist[edgelist[i]->verts[1]]->InRegion = 0;
	}

	for(i = 0; i < trianglesnum; i++)
	{
		face = Object.flist[trianglelist[i]];
		if(face->contain_singularity == 1)
		{
			if(face->singularity_index == saddleID) ////containing the saddle need to be cancelled
				continue;

			for(j = 0; j < face->nverts; j++)
			{
				cur_v = Object.vlist[face->verts[j]];
				cur_v->OnBoundary = 1;
			}
		}
	}

	////Get the intersection of inner vertices
	intersectInnerverts.num = 0;
	for(i = 0; i < vertsnum; i++)
	{
		cur_v = vertslist[i];
		
		if(cur_v->OnBoundary != 1)
		{
			intersectInnerverts.vertslist[intersectInnerverts.num] = cur_v;
			cur_v->InRegion = 1;
			cur_v->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
			intersectInnerverts.num++;
		}
	}
}

////New method to grow region from saddle 2/6/06
void GetLimitSaddlePairRegion2(int saddleID, int limitID, int singID)
{
	int type = limitcycles[limitID].type;

	if(type == 0) ////the saddle acts as an attractor
	{
		attractorRegion.trianglelist[0] = singularities[saddleID].Triangle_ID;
		Object.flist[singularities[saddleID].Triangle_ID]->attract_inregion = 1;
		attractorRegion.num = 1;
	}
	else{
		repellerRegion.trianglelist[0] = singularities[saddleID].Triangle_ID;
		Object.flist[singularities[saddleID].Triangle_ID]->repell_inregion = 1;
		repellerRegion.num = 1;
	}

    Cancel_InitLimitSaddleRegion2(saddleID, type, singID); //we need the 'singID' instead of initlength now
    UpdateBoundary(1-type);
	GetRegionNormals(1-type);
	Region_Growing(1-type);

	UpdateBoundary(1-type);

}


////This routine will be called after building the intersection region 2/7/06
void GetInnerVerts(int singID, int *Medianodes, int Num_medianodes)
{
	Face *face;
	Edge *edge;
	Vertex *v;
	int i, j;

	UpdateBoundary(2);

	for(i = 0; i < Object.nverts; i++)
	{
		Object.vlist[i]->OnBoundary = 0;
		Object.vlist[i]->RegionListID = -1;
		Object.vlist[i]->InRegion = 0;
	}

	////Set the boundary vertices
	for(i = 0; i < intersectBoundary.num; i++)
	{
		edge = intersectBoundary.edgelist[i];
		v = Object.vlist[edge->verts[0]];
		v->OnBoundary = 1;
		
		v = Object.vlist[edge->verts[1]];
		v->OnBoundary = 1;
	}

	////Get the intersection of inner vertices
	intersectInnerverts.num = 0;

	for(i = 0; i < intersectRegion.num; i++)
	{
		face = Object.flist[intersectRegion.trianglelist[i]];

		for(j = 0; j < face->nverts; j++)
		{
			v = Object.vlist[face->verts[j]];

			if(v->RegionListID >= 0 || v->OnBoundary == 1)
				continue;

			if(intersectInnerverts.num >= MaxNumInnerVerts -1)
			{
				MaxNumInnerVerts += 200;
				intersectInnerverts.vertslist = (Vertex**)realloc(intersectInnerverts.vertslist,
					sizeof(Vertex*)*MaxNumInnerVerts);
			}

			intersectInnerverts.vertslist[intersectInnerverts.num] = v;
			v->RegionListID = intersectInnerverts.num;
			v->InRegion = 1;
			intersectInnerverts.num ++;
		}
	}

	////Add the vertices of the triangles that contain the singularity to the innerverts list
	if(singID >= 0)
	{
		face = Object.flist[singularities[singID].Triangle_ID];
		for(j = 0; j < face->nverts; j++)
		{
			v = Object.vlist[face->verts[j]];

			if(v->InRegion != 1)
			{
				intersectInnerverts.vertslist[intersectInnerverts.num] = v;
				v->InRegion = 1;
				v->OnBoundary = 0;
				v->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
				intersectInnerverts.num ++;
			}
		}
	}

	////Add the intervals (saddles)
	for(i = 0; i < Num_medianodes; i++)
	{
		face = Object.flist[singularities[graphnodes[Medianodes[i]].singularityID].Triangle_ID];
		for(j = 0; j < face->nverts; j++)
		{
			v = Object.vlist[face->verts[j]];

			if(v->InRegion != 1)
			{
				intersectInnerverts.vertslist[intersectInnerverts.num] = v;
				v->InRegion = 1;
				v->OnBoundary = 0;
				v->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
				intersectInnerverts.num ++;
			}
		}
	}
}


/////We should change to use more efficient method
////need not perform local tracing again! We just track the triangles along the separatrices

void Cancel_InitLimitSaddleRegion2(int singularID, int type, int singID)
{
	int i;
	int sep, traj;

	sep = singularities[singularID].separtices;


	if(type == 0)    ////Saddle acts as an attractor, forward tracing along outgoing directions
	{
		////1. Positive direction
		traj = separatrices[sep].sep1;
		if(singID == -1 
			|| trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID == singularities[singID].Triangle_ID)
		{
			for(i = 0; i < num_linesegs_curtraj[traj]; i++)
			{
				if(trajectories[traj][i].Triangle_ID != singularities[singularID].Triangle_ID
					&& Object.flist[trajectories[traj][i].Triangle_ID]->contain_singularity == 1)
					break;

				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 1);
			}
		}
		else{ //we pick half of the separatrix
			for(i = 0; i < (int)num_linesegs_curtraj[traj]/2; i++)
			{
				if(trajectories[traj][i].Triangle_ID != singularities[singularID].Triangle_ID
					&& Object.flist[trajectories[traj][i].Triangle_ID]->contain_singularity == 1)
					break;

				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 1);
			}
		}

		////2. Negative direction
		traj = separatrices[sep].sep3;
		if(singID == -1 
			|| trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID == singularities[singID].Triangle_ID)
		{
			for(i = 0; i < num_linesegs_curtraj[traj]; i++)
			{
				if(trajectories[traj][i].Triangle_ID != singularities[singularID].Triangle_ID
					&& Object.flist[trajectories[traj][i].Triangle_ID]->contain_singularity == 1)
					break;
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 1);
			}
		}
		else{ //we pick half of the separatrix
			for(i = 0; i < (int)num_linesegs_curtraj[traj]/2; i++)
			{
				if(trajectories[traj][i].Triangle_ID != singularities[singularID].Triangle_ID
					&& Object.flist[trajectories[traj][i].Triangle_ID]->contain_singularity == 1)
					break;

				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 1);
			}
		}
	}

	else  ////Saddle acts as a repeller, follow the incoming directions
	{
		////1. Positive direction
		traj = separatrices[sep].sep2;
		if(singID == -1 
			|| trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID == singularities[singID].Triangle_ID)
		{
			for(i = 0; i < num_linesegs_curtraj[traj]; i++)
			{
				if(trajectories[traj][i].Triangle_ID != singularities[singularID].Triangle_ID
					&& Object.flist[trajectories[traj][i].Triangle_ID]->contain_singularity == 1)
					break;
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 0);
			}
		}
		else{ //we pick half of the separatrix
			for(i = 0; i < (int)num_linesegs_curtraj[traj]/2; i++)
			{
				if(trajectories[traj][i].Triangle_ID != singularities[singularID].Triangle_ID
					&& Object.flist[trajectories[traj][i].Triangle_ID]->contain_singularity == 1)
					break;

				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 1);
			}
		}

		////2. Negative direction
		traj = separatrices[sep].sep4;
		if(singID == -1 
			|| trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID == singularities[singID].Triangle_ID)
		{
			for(i = 0; i < num_linesegs_curtraj[traj]; i++)
			{
				if(trajectories[traj][i].Triangle_ID != singularities[singularID].Triangle_ID
					&& Object.flist[trajectories[traj][i].Triangle_ID]->contain_singularity == 1)
					break;
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 0);
			}
		}
		else{ //we pick half of the separatrix
			for(i = 0; i < (int)num_linesegs_curtraj[traj]/2; i++)
			{
				if(trajectories[traj][i].Triangle_ID != singularities[singularID].Triangle_ID
					&& Object.flist[trajectories[traj][i].Triangle_ID]->contain_singularity == 1)
					break;

				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 1);
			}
		}
	}
}


////Grow intial triangle strip for limit cycle pair cancellation 2/8/06
void InitSaddleRegionforLimitPairCancel(int saddle, int type, int initlength)
{
	int i;
	int sep, traj;

	sep = singularities[saddle].separtices;

	if(type == 0)  //saddle acts as a repeller, we need to grow along incoming directions
	{
		////Positive incoming direction
		traj = separatrices[sep].sep2;
		for(i = 0; i < num_linesegs_curtraj[traj]/10; i++)
		{
			AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 0);
		}

		////Negative incoming direction
		traj = separatrices[sep].sep4;
		for(i = 0; i < num_linesegs_curtraj[traj]/10; i++)
		{
			AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 0);
		}
	}

	else           //saddle acts as an attractor, we need to grow along outgoing directions
	{
		////Positive outgoing direction
		traj = separatrices[sep].sep1;
		for(i = 0; i < num_linesegs_curtraj[traj]/10; i++)
		{
			AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 1);
		}

		////Negative outgoing direction
		traj = separatrices[sep].sep3;
		for(i = 0; i < num_linesegs_curtraj[traj]/10; i++)
		{
			AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 1);
		}
	}
}


void Cancel_InitLimitSaddleRegion(int singularID, int type, int initsaddlelength)
{
	icVector2 going;

	int i;
	int flag = 0;
	double globalp[2];
	int pre_face, cur_face;
	int loop_flag = 0;  ////for detect a closed loop from the beginning triangle
	int *triangles = (int*)malloc(sizeof(int) * NUMTRACINGTRIANGLE);
	int num_triangles = 0;


	if(type == 0)    ////Saddle as an attractor, forward tracing along outgoing direction
	{
		going = singularities[singularID].outgoing;
		
		////tracing along both directions
		//1. positive outgoing direction
		pre_face = cur_face = singularities[singularID].Triangle_ID;
		globalp[0] = singularities[singularID].gcx + SEPARATRIXSTEP * going.entry[0];   
		globalp[1] = singularities[singularID].gcy + SEPARATRIXSTEP * going.entry[1];
		for(i = 0; i < initsaddlelength; i++)
		{
			if(cur_face < 0)
			{
				break;
			}

			pre_face = cur_face;
			cur_face = TraceInATriangle2(cur_face, globalp, type, flag); ////0 means always forward tracing here

			if(flag == 3 || flag == 4 || pre_face == cur_face || loop_flag == 1) 
			{
				if(flag == 3)  ////reached a singularity, we should not include this separatrix!
				{
					if(type == 1)
						repellerRegion.num -= 2;
					else
						attractorRegion.num -= 2;
				}
				break;

			}

			AddToRegionTriangles(cur_face, 1);
			
			if(LoopDetect(triangles, num_triangles, cur_face))
				loop_flag = 1;
			else{
				triangles[num_triangles] = cur_face;
				num_triangles++;
			}
		}

		//2. negative outgoing direction
		flag = 0;
		loop_flag = 0;
		pre_face = cur_face = singularities[singularID].Triangle_ID;
		globalp[0] = singularities[singularID].gcx - SEPARATRIXSTEP * going.entry[0];   
		globalp[1] = singularities[singularID].gcy - SEPARATRIXSTEP * going.entry[1];
		for(i = 0; i < initsaddlelength; i++)
		{
			if(cur_face < 0)
			{
				break;
			}

			pre_face = cur_face;
			cur_face = TraceInATriangle2(cur_face, globalp, type, flag); ////0 means always forward tracing here

			if(flag == 3 || flag == 4 || pre_face == cur_face || loop_flag == 1) 
			{
				if(flag == 3)  ////reached a singularity, we should not include this separatrix!
				{
					if(type == 1)
						repellerRegion.num -= 2;
					else
						attractorRegion.num -= 2;
				}
				break;
			}
			AddToRegionTriangles(cur_face, 1);
			//attractorRegion.trianglelist[attractorRegion.num] = cur_face;
			//attractorRegion.num ++;
			if(LoopDetect(triangles, num_triangles, cur_face))
				loop_flag = 1;
			else{
				triangles[num_triangles] = cur_face;
				num_triangles++;
			}
		}
	}
	else{            ////Saddle as a repeller, backward tracing along incoming direction
		going = singularities[singularID].incoming;
		
		////tracing along both directions
		//1. positive outgoing direction
		pre_face = cur_face = singularities[singularID].Triangle_ID;
		globalp[0] = singularities[singularID].gcx + SEPARATRIXSTEP * going.entry[0];   
		globalp[1] = singularities[singularID].gcy + SEPARATRIXSTEP * going.entry[1];
		for(i = 0; i < InitSaddleRegionLength; i++)
		{
			if(cur_face < 0)
			{
				break;
			}

			pre_face = cur_face;
			cur_face = TraceInATriangle2(cur_face, globalp, 1, flag); ////1 means always backward tracing here

			if(flag == 3 || flag == 4 || pre_face == cur_face || loop_flag == 1) 
			{
				if(flag == 3)  ////reached a singularity, we should not include this separatrix!
				{
					if(type == 1)
						repellerRegion.num -= 2;
					else
						attractorRegion.num -= 2;
				}
				break;
			}

			//repellerRegion.trianglelist[repellerRegion.num] = cur_face;
			//repellerRegion.num ++;
			AddToRegionTriangles(cur_face, 0);
			
			if(LoopDetect(triangles, num_triangles, cur_face))
				loop_flag = 1;
			else{
				triangles[num_triangles] = cur_face;
				num_triangles++;
			}
		}

		//2. negative outgoing direction
		flag = 0;
		loop_flag = 0;
		pre_face = cur_face = singularities[singularID].Triangle_ID;
		globalp[0] = singularities[singularID].gcx - SEPARATRIXSTEP * going.entry[0];   
		globalp[1] = singularities[singularID].gcy - SEPARATRIXSTEP * going.entry[1];
		for(i = 0; i < InitSaddleRegionLength; i++)
		{
			if(cur_face < 0)
			{
				break;
			}

			pre_face = cur_face;
			cur_face = TraceInATriangle2(cur_face, globalp, 1, flag); ////0 means always forward tracing here

			if(flag == 3 || flag == 4 || pre_face == cur_face || loop_flag == 1) 
			{
				if(flag == 3)  ////reached a singularity, we should not include this separatrix!
				{
					if(type == 1)
						repellerRegion.num -= 2;
					else
						attractorRegion.num -= 2;
				}
				break;
			}
			//repellerRegion.trianglelist[repellerRegion.num] = cur_face;
			//repellerRegion.num ++;
			AddToRegionTriangles(cur_face, 0);
			
			if(LoopDetect(triangles, num_triangles, cur_face))
				loop_flag = 1;
			else{
				triangles[num_triangles] = cur_face;
				num_triangles++;
			}
		}
	}

	free(triangles);
}


////Get the intersect region for smoothing
////Modified on 09/06/05 
////Suppose we have two adjacent saddles to current singularity

void GetLimitSingIntersectRegion(int singID/*, int saddleID*/)
{
	int i, j;
	Face *face, *face1, *face2;
	Edge *cur_edge;
	Vertex *vert;

	////Get the intersect triangular region
	intersectRegion.trianglelist[0] = singularities[singID].Triangle_ID;
	intersectRegion.num = 1;

	////The following adding may have problem! 09/04/05
	//intersectRegion.trianglelist[1] = singularities[saddleID].Triangle_ID;
	//intersectRegion.num = 2;

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
	////First, we need to reset the flag of edge visiting
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

	////Get the intersection of inner vertices
	intersectInnerverts.num = 0;

    ////we need to add the vertices on the boundary into the vertices list? why?
	////Modified on 09/04/05

	for(i = 0; i < intersectBoundary.num; i++)
	{
		Object.vlist[intersectBoundary.edgelist[i]->verts[0]]->OnBoundary = 0;
		Object.vlist[intersectBoundary.edgelist[i]->verts[0]]->attract_flag = 2;
		Object.vlist[intersectBoundary.edgelist[i]->verts[0]]->repell_flag = 2;

		Object.vlist[intersectBoundary.edgelist[i]->verts[1]]->OnBoundary = 0;
		Object.vlist[intersectBoundary.edgelist[i]->verts[1]]->attract_flag = 2;
		Object.vlist[intersectBoundary.edgelist[i]->verts[1]]->repell_flag = 2;
	}

	for(i = 0; i < intersectBoundary.num; i++)
	{
		cur_edge = intersectBoundary.edgelist[i];

		face1 = Object.flist[cur_edge->tris[0]];
		face2 = Object.flist[cur_edge->tris[1]];

		if(face1->contain_singularity == 1 || face2->contain_singularity == 1)
		{
			if(face1->contain_singularity == 1 && face1->singularity_index != singID 
				&& singularities[face1->singularity_index].type != SADDLE)  ////not the cancelled sing and its neighboring saddles
			{
				for(j = 0; j < face1->nverts; j++)
				{
					vert = Object.vlist[face1->verts[j]];
					vert->OnBoundary = 1;
				}
			}
			if(face2->contain_singularity == 1 && face2->singularity_index != singID
				&& singularities[face2->singularity_index].type != SADDLE)  ////not the cancelled sing and its neighboring saddles
			{
				for(j = 0; j < face2->nverts; j++)
				{
					vert = Object.vlist[face2->verts[j]];
					vert->OnBoundary = 1;
				}
			}
		}

	}

	////Set all the vertices of the triangles that do not contain the adjacent saddles as boundary!!! 09/04/05
	//for(i = 0; i < intersectRegion.num; i++)
	//{
	//	face = Object.flist[intersectRegion.trianglelist[i]];
	//	if(face->contain_singularity == 1)
	//	{
	//		if(singularities[singID].Triangle_ID == face->index) ////containing the singularity need to be cancelled
	//			continue;
	//		if(singularities[face->singularity_index].type == SADDLE) ////containing the adjacent saddles
	//			continue;
	//		for(j = 0; j < face->nverts; j++)
	//		{
	//			vert = Object.vlist[face->verts[j]];
	//			vert->OnBoundary = 1;
	//		}
	//	}
	//}


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
	face = Object.flist[singularities[singID].Triangle_ID];
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



void GetLimitSingIntersectRegion2(int singID, int limitID)
{
}


////Remove the fence for limit cyclc involved cancellation 06/19/06

//void RemoveLimitCycleFence(int limitID)
//{
//	int i;
//	int triangle;
//
//	for(i = 0; i < limitcycles[limitID].num_triangles; i++)
//	{
//		triangle = limitcycles[limitID].cellcycle[i];
//		Object.flist[triangle]->fence_flag = 0;
//	}
//}


bool LimitSaddleCancel(int saddle, int limitID)
{
	int temp_repeller[1], temp_attractor[1];
	double cur_length, large_length, small_length, success_length;
	int saddletype;
	int i;
	
	int num_t1, num_t2;

	//initialize
	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->attract_inregion = 0;
		Object.flist[i]->repell_inregion = 0;
	}

	repellerRegion.num = 0;
	attractorRegion.num = 0;
	intersectRegion.num = 0;

	//1. Grow region for the limit cycle (this will be the same during optimization)
	//GrowRegionforALimitCycle(limitID);
	
	////set the node for updating the C-graph
	if(limitcycles[limitID].type == 0) //saddle acts as an attractor
	{
		temp_repeller[0] = limitcycles[limitID].node_index;
		temp_attractor[0] = singularities[saddle].node_index;
		saddletype = 1;
	}
	else //saddle acts as a repeller
	{
		temp_repeller[0] = singularities[saddle].node_index;
		temp_attractor[0] = limitcycles[limitID].node_index;
		saddletype = 0;
	}

	//2. Optimization process
	small_length = 0;
	large_length = 1.;
	cur_length = large_length;
	int count = 0;
	int success = 0;

	if(saddletype == 0) //saddle acts as a repeller
	{
		//get the two incoming directions
		num_t1 = num_linesegs_curtraj[separatrices[singularities[saddle].separtices].sep2];
		num_t2 = num_linesegs_curtraj[separatrices[singularities[saddle].separtices].sep4];
	}
	else
	{
		//get the two outgoing directions
		num_t1 = num_linesegs_curtraj[separatrices[singularities[saddle].separtices].sep1];
		num_t2 = num_linesegs_curtraj[separatrices[singularities[saddle].separtices].sep3];
	}

	int minnum = min(num_t1, num_t2);

	///we need to set fences for all the other saddles
	SetFenceForSeps(saddle);
	SetFence_LimitCyclesexcept(limitID);
	RemoveLimitCycleFence(limitID);
	RemoveFenceForSaddle(saddle);     //remove the fences for the separatrices of the saddle(temporary)

	do{
	
		//1. Grow region for the limit cycle (this will be the same during optimization)
		InitCancellationAndMovement();

		GrowRegionforALimitCycle(limitID);

		//2.1 ) Get the region for the saddle

		////initialize the previous saddle growing region
		for(i = 0; i < intersectRegion.num; i++)
		{
			if(saddletype == 0)
				Object.flist[intersectRegion.trianglelist[i]]->repell_inregion = 0;
			else
				Object.flist[intersectRegion.trianglelist[i]]->attract_inregion = 0;
		}

		if(saddletype == 0)
			repellerRegion.num = 0;
		else
			attractorRegion.num = 0;

		intersectRegion.num = 0;

		InitSaddleGrowforMultRegion(saddle, saddletype, cur_length, NULL, 0);

		UpdateBoundary(saddletype);
		GetRegionNormals(saddletype);
		Cancel_Growing(saddletype, -1);
		
		//2.2 ) Get the intersection region of the saddle and the limit cycle
		IntersectRegion(repellerRegion, attractorRegion, intersectRegion);

		////Add the triangle containing the saddle
		if(!IsRepeated(intersectRegion.trianglelist, singularities[saddle].Triangle_ID, intersectRegion.num))
		{
			intersectRegion.trianglelist[intersectRegion.num] = 
				singularities[saddle].Triangle_ID;
			intersectRegion.num ++;
		}

		//2.3 ) Judge whether the intersection region satisfies the Conley boundary condition or not
		if(CalEulerValue(intersectRegion.trianglelist, intersectRegion.num) <= 0) //a ring-shaped
		{
			success = 1; //we do find a valid region to perform smoothing
			success_length = cur_length;
			if(count == 0)
			{
				break;
			}

			small_length = cur_length;  //get a bigger region
		}

		else
		{
			large_length = cur_length;  //get a smaller region
		}

		cur_length = (small_length + large_length)/2.;  //binary search here
		
		count ++;

	}while(large_length - small_length > 1e-5 && (double)minnum*(large_length - small_length) > 0.5
		&& (int)minnum*cur_length >= 1/* && count <= 0*/);

	if(success == 1)
	{
		Face *face;
		Vertex *vert;
		Edge *cur_e;

		if(count != 0)
		{
			//need to use the successful length to grow the region again!
			for(i = 0; i < intersectRegion.num; i++)
			{
				if(saddletype == 0)
					Object.flist[intersectRegion.trianglelist[i]]->repell_inregion = 0;
				else
					Object.flist[intersectRegion.trianglelist[i]]->attract_inregion = 0;
			}

			if(saddletype == 0)
				repellerRegion.num = 0;
			else
				attractorRegion.num = 0;

			intersectRegion.num = 0;

			InitSaddleGrowforMultRegion(saddle, saddletype, success_length, NULL, 0);
			UpdateBoundary(saddletype);
			GetRegionNormals(saddletype);
			Cancel_Growing(saddletype, -1);
			
			IntersectRegion(repellerRegion, attractorRegion, intersectRegion);
		}

		////Add the triangle containing the saddle
		if(!IsRepeated(intersectRegion.trianglelist, singularities[saddle].Triangle_ID, intersectRegion.num))
		{
			intersectRegion.trianglelist[intersectRegion.num] = 
				singularities[saddle].Triangle_ID;
			intersectRegion.num ++;
		}

		//get the inner vertices of the region
		for(i = 0; i < Object.nfaces; i++)
		{
			face = Object.flist[i];

			for(int j = 0; j < 3; j++)
			{
				cur_e = face->edges[j];
				cur_e->OnBoundary = 0;

				Object.vlist[face->verts[j]]->RegionListID = -1; //reset the region index of the vertex
				Object.vlist[face->verts[j]]->InRegion = 0;
			}
		}
		UpdateBoundary(2);  //update the boundary of the intersect region

		for(i = 0; i < intersectBoundary.num; i++)
		{
			cur_e = intersectBoundary.edgelist[i];
			Object.vlist[cur_e->verts[0]]->OnBoundary = 1;
			Object.vlist[cur_e->verts[1]]->OnBoundary = 1;
		}

		intersectInnerverts.num = 0; //reset the number of inner vertices
		//get the inner vertices of the valid region
		for(i = 0; i < intersectRegion.num; i++)
		{
			face = Object.flist[intersectRegion.trianglelist[i]];

			for(int j = 0; j < face->nverts; j++)
			{
				vert = Object.vlist[face->verts[j]];

				//add to the inner vertex list!
				if(vert->OnBoundary == 0 && vert->RegionListID < 0)
				{
					intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
					vert->InRegion = 1;
					vert->RegionListID = intersectInnerverts.num;
					intersectInnerverts.num++;
				}
			}
		}

		////Add the 3 vertices of the triangle containing the saddle
		face = Object.flist[singularities[saddle].Triangle_ID];
		for(i = 0; i < face->nverts; i++)
		{
			vert = Object.vlist[face->verts[i]];

			if(vert->InRegion != 1 && vert->RegionListID < 0)
			{
				intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
				vert->InRegion = 1;
				vert->RegionListID = intersectInnerverts.num;
				intersectInnerverts.num++;
			}
		}

		//Perform one smoothing

		//save the previous field
		Cancel_RegionSmooth();

		//Calculate the Poincare index of the region
		int totalindex = 0;
		IntersectedRegionSingCount(totalindex);

		//we also need to perform limit cycle detection inside this region 06/22/06

		if(totalindex == -1)
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





////Modified at 2/06/06 and 2/07/06

void CancelLimitSingPair(int singID, int limitID)
{
	////1. Initialize the cancellation
	//InitCancellationAndMovement();

	////2. Find all the intervals between the being cancelled components
	////Here we suppose that we perform only single limit cycle and single singularity cancellation
	int temp_repeller[1];
	int temp_attractor[1];
	
	MediaNodes = (int *)malloc(sizeof(int)*MaxMediaNodes);
	Num_MediaNodes = 0;

	////
	TriangularRegion Temp_re;
	Temp_re.trianglelist = (int*)malloc(sizeof(int)*Object.nfaces);
	Temp_re.num = 0;

	////Limit cycle and saddle with two connections
	if(singularities[singID].type == SADDLE && IsTwoConnections(limitID, singID))
	{

		//set the being cancelled nodes for updating the C-graph
		if(limitcycles[limitID].type == 0)
		{
			temp_repeller[0] = limitcycles[limitID].node_index;
			temp_attractor[0] = singularities[singID].node_index;
		}
		else
		{
			temp_attractor[0] = limitcycles[limitID].node_index;
			temp_repeller[0] = singularities[singID].node_index;
		}

		LimitSaddleCancelWithTwoOrbits(limitID, singID);
		goto LL;
	}
	

	////Remove the fence of the being cancelled limit cycle 06/19/06
	RemoveLimitCycleFence(limitID);

	if(limitcycles[limitID].type == 0) //the limit cycle is a repeller
	{
		if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
		{
			MessageBox(NULL, "Not a correct cancelled pair!", "Error", MB_OK);
			goto L2;
		}


		if(singularities[singID].type == SADDLE)
		{
			////for saddle, we need not find the intermediary components

			temp_repeller[0] = limitcycles[limitID].node_index;
			temp_attractor[0] = singularities[singID].node_index;

			if(LimitSaddleCancel(singID, limitID))
			{
				goto LL;
			}

			else
			{
				MessageBox(NULL, "fail to cancel the pair!", "Fail", MB_OK);
				goto L2;
			}

		}
		else{
			temp_repeller[0] = limitcycles[limitID].node_index;
			temp_attractor[0] = singularities[singID].node_index;

			//search the intermediate components (saddles)
			SearchConnectComponents_adv(temp_repeller, 1, temp_attractor, 1, MediaNodes, Num_MediaNodes);
			if(Num_MediaNodes == 0)
			{
				if(DirectlyConnectedSingandCyclePair(singID, limitID))
					goto LL;
			}

			if(IndirectlyConnectedSingandCyclePair(singID, limitID)) //indirectly linked limit cycle and singularity cancellation
				goto LL;
		}
	}

	else //the limit cycle is an attractor
	{
		if(singularities[singID].type == SINK || singularities[singID].type == AFOCUS)
		{
			MessageBox(NULL, "Not a correct cancelled pair!", "Error", MB_OK);
			goto L2;
		}

		if(singularities[singID].type == SADDLE)
		{
			////for saddle, we need not find the intermediary components
			
			temp_repeller[0] = singularities[singID].node_index;
			temp_attractor[0] = limitcycles[limitID].node_index;
			
			if(LimitSaddleCancel(singID, limitID))
			{
				goto LL;
			}

			else
			{
				MessageBox(NULL, "fail to cancel the pair!", "Fail", MB_OK);
				goto L2;
			}
		}

		else
		{
			temp_repeller[0] = singularities[singID].node_index;
			temp_attractor[0] = limitcycles[limitID].node_index;
			
			//search the intermediate components (saddles)
			SearchConnectComponents_adv(temp_repeller, 1, temp_attractor, 1, MediaNodes, Num_MediaNodes);
			if(Num_MediaNodes == 0)
			{
				if(DirectlyConnectedSingandCyclePair(singID, limitID))
					goto LL;
			}

			if(IndirectlyConnectedSingandCyclePair(singID, limitID)) //indirectly linked limit cycle and singularity cancellation
			{
				goto LL;
			}
		}
	}


	////if successful(Satisfies Conley index/Poincare index), mark those nodes 'cancelled'
LL:	MarkCancel(temp_repeller, 1, temp_attractor, 1, MediaNodes, Num_MediaNodes);

L2:	free(MediaNodes);
	Num_MediaNodes = 0;
	free(Temp_re.trianglelist);
}


////The main routine for region growing 
void Region_Growing(int type)
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

			if(Object.flist[oppositeTriangle]->contain_singularity == 1 || oppositeTriangle < 0) ////if the opposite triangle containing singularity
			{
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


/*-----------------------------------------------------------------------*/
////Routines for limit cycle pair cancellation
void GetLimitPairIntersectRegion()
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

	////Get the intersection of inner vertices
	intersectInnerverts.num = 0;

    ////we need to add the vertices on the boundary into the vertices list(why?)

	for(i = 0; i < intersectBoundary.num; i++)
	{
		Object.vlist[intersectBoundary.edgelist[i]->verts[0]]->OnBoundary = 1;
		//Object.vlist[intersectBoundary.edgelist[i]->verts[0]]->attract_flag = 2;
		//Object.vlist[intersectBoundary.edgelist[i]->verts[0]]->repell_flag = 2;

		Object.vlist[intersectBoundary.edgelist[i]->verts[1]]->OnBoundary = 1;
		//Object.vlist[intersectBoundary.edgelist[i]->verts[1]]->attract_flag = 2;
		//Object.vlist[intersectBoundary.edgelist[i]->verts[1]]->repell_flag = 2;
	}

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

void CancelEmbededLimitCycles(int inner, int outer)
{
	if(limitcycles[inner].type == limitcycles[outer].type)
	{
		MessageBox(NULL, "Not a valid pair of limit cycles!", "Error", MB_OK);
		return;
	}

	////if two limit cycles are not embeded, they can not be cancelled!

	int i;
	int temp_repeller[1], temp_attractor[1];

	int cur_saddle;
	TriangularRegion Source_re, Sink_re;
	Source_re.trianglelist = (int*)malloc(sizeof(int)*Object.nfaces);
	Source_re.num = 0;
	Sink_re.trianglelist = (int*)malloc(sizeof(int)*Object.nfaces);
	Sink_re.num = 0;

	if(MediaNodes != NULL)
	{
		free(MediaNodes);
	}
	MediaNodes = (int*)malloc(sizeof(int)*MaxMediaNodes);
	Num_MediaNodes = 0;

	InitCancellationAndMovement();

	SetFence_LimitCycles(inner, outer, MediaNodes, Num_MediaNodes);


	////get smoothing region
	////Copy the cell cycle to the corresponding variables
	if(limitcycles[inner].type == 0)  ////the limit cycle is a repeller
	{
		////Get the interval components through the Conley connection graph
		temp_repeller[0] = limitcycles[inner].node_index;
		temp_attractor[0] = limitcycles[outer].node_index;
		SearchConnectComponents_adv(temp_repeller, 1, temp_attractor, 1, MediaNodes, Num_MediaNodes);

		////Grow inner limit cycle region
		repellerRegion.num = 0;
		for(i = 0; i < limitcycles[inner].num_triangles; i++)
			repellerRegion.trianglelist[i] = limitcycles[inner].cellcycle[i];
		repellerRegion.num = limitcycles[inner].num_triangles;
		UpdateBoundary(0);
		GetRegionNormals(0);
		Cancel_Growing(0, -1);

		CopyRegion(repellerRegion.trianglelist, Source_re.trianglelist, repellerRegion.num);
		Source_re.num = repellerRegion.num;

		////Grow region from the saddle as repeller

		for(i = 0; i < Num_MediaNodes; i++)
		{
			repellerRegion.num = 0;

			cur_saddle = graphnodes[MediaNodes[i]].singularityID;
			InitSaddleRegionforLimitPairCancel(cur_saddle, 0, 20);
			Cancel_Growing(0, -1);

			CopyRegion(Source_re.trianglelist, intersectRegion.trianglelist, Source_re.num);
			intersectRegion.num = Source_re.num;

			UnionRegion(repellerRegion, intersectRegion, Source_re);
		}


		////Grow outer limit cycle region
		attractorRegion.num = 0;
		for(i = 0; i < limitcycles[outer].num_triangles; i++)
			attractorRegion.trianglelist[i] = limitcycles[outer].cellcycle[i];
		attractorRegion.num = limitcycles[outer].num_triangles;
		UpdateBoundary(1);
		GetRegionNormals(1);
		Cancel_Growing(1, -1);

		CopyRegion(attractorRegion.trianglelist, Sink_re.trianglelist, attractorRegion.num);
		Sink_re.num = attractorRegion.num;

		////Grow region from the saddle as attractor
		for(i = 0; i < Num_MediaNodes; i++)
		{
			attractorRegion.num = 0;

			cur_saddle = graphnodes[MediaNodes[i]].singularityID;
			InitSaddleRegionforLimitPairCancel(cur_saddle, 1, 20);
			Cancel_Growing(1, -1);

			CopyRegion(Sink_re.trianglelist, intersectRegion.trianglelist, Sink_re.num);
			intersectRegion.num = Sink_re.num;

			UnionRegion(attractorRegion, intersectRegion, Sink_re);
		}
	}

	else{
		////Get the interval components through the Conley connection graph
		temp_repeller[0] = limitcycles[outer].node_index;
		temp_attractor[0] = limitcycles[inner].node_index;
		SearchConnectComponents_adv(temp_repeller, 1, temp_attractor, 1, MediaNodes, Num_MediaNodes);

		////Grow inner limit cycle region
		attractorRegion.num = 0;
		for(i = 0; i < limitcycles[inner].num_triangles; i++)
			attractorRegion.trianglelist[i] = limitcycles[inner].cellcycle[i];
		attractorRegion.num = limitcycles[inner].num_triangles;
		UpdateBoundary(1);
		GetRegionNormals(1);
		Cancel_Growing(1, -1);
		
		CopyRegion(attractorRegion.trianglelist, Sink_re.trianglelist, attractorRegion.num);
		Sink_re.num = attractorRegion.num;

		////Grow region from the saddle as attractor
		for(i = 0; i < Num_MediaNodes; i++)
		{
			attractorRegion.num = 0;

			cur_saddle = graphnodes[MediaNodes[i]].singularityID;
			InitSaddleRegionforLimitPairCancel(cur_saddle, 1, 20);
			Cancel_Growing(1, -1);

			CopyRegion(Sink_re.trianglelist, intersectRegion.trianglelist, Sink_re.num);
			intersectRegion.num = Sink_re.num;

			UnionRegion(attractorRegion, intersectRegion, Sink_re);
		}

		////Grow outer limit cycle region
		repellerRegion.num = 0;
		for(i = 0; i < limitcycles[outer].num_triangles; i++)
			repellerRegion.trianglelist[i] = limitcycles[outer].cellcycle[i];
		repellerRegion.num = limitcycles[outer].num_triangles;
		UpdateBoundary(0);
		GetRegionNormals(0);
		Cancel_Growing(0, -1);

		CopyRegion(repellerRegion.trianglelist, Source_re.trianglelist, repellerRegion.num);
		Source_re.num = repellerRegion.num;

		////Grow region from the saddle as repeller

		for(i = 0; i < Num_MediaNodes; i++)
		{
			repellerRegion.num = 0;

			cur_saddle = graphnodes[MediaNodes[i]].singularityID;
			InitSaddleRegionforLimitPairCancel(cur_saddle, 0, 20);
			Cancel_Growing(0, -1);

			CopyRegion(Source_re.trianglelist, intersectRegion.trianglelist, Source_re.num);
			intersectRegion.num = Source_re.num;

			UnionRegion(repellerRegion, intersectRegion, Source_re);
		}
	}

	////Get the intersect region of these two limit cycles
	//GetLimitPairIntersectRegion();
	IntersectRegion(Source_re, Sink_re, intersectRegion);

	////Add the triangles that contain the saddles into the region
	int position = -1;
	int cur_t;
	for(int i = 0; i < Num_MediaNodes; i++)
	{
		cur_t = singularities[graphnodes[MediaNodes[i]].singularityID].Triangle_ID;

		if(!TriangleSearch(intersectRegion.trianglelist, intersectRegion.num, cur_t, position))
		{
			intersectRegion.trianglelist[intersectRegion.num] = cur_t;
			intersectRegion.num++;
		}
	}


	////Get the inner vertices of the intersect region
	GetInnerVerts(-1, MediaNodes, Num_MediaNodes);

	////remove the center triangle of the inner limit cycle
	if(limitcycles[inner].singularID >= 0)
	{
		Face *face = Object.flist[singularities[limitcycles[inner].singularID].Triangle_ID];
		for(int i = 0; i < face->nverts; i++)
		{
			////Remove from the inner vertices list and update the number
			Object.vlist[face->verts[i]]->OnBoundary = 1;
			Object.vlist[face->verts[i]]->InRegion = 0;
		}
	}

	////Get inner vertices
    SavePreviousField();
	Cancel_RegionSmooth();
	SavePostField();

	////if success
	MarkCancel(temp_repeller, 1, temp_attractor, 1, MediaNodes, Num_MediaNodes);

	free(MediaNodes);
	Num_MediaNodes = 0;
	free(Source_re.trianglelist);
	free(Sink_re.trianglelist);
}




////Cancel two embeded limit cycles
void CancelPairLimitCycles(int inner, int outer)
{
	InitCancellationAndMovement();

	////Get the smoothing region from the cell cycle of the inner limit cycle
	GetInnerLimitCycleRegion(outer, limitcycles[outer].cellcycle, limitcycles[outer].num_triangles);

	////Set the center triangle as second boundary 
	Face *face = Object.flist[singularities[limitcycles[inner].singularID].Triangle_ID];
	for(int i = 0; i < face->nverts; i++)
	{
		////Remove from the inner vertices list and update the number
		Object.vlist[face->verts[i]]->OnBoundary = 1;
		Object.vlist[face->verts[i]]->InRegion = 0;
	}
    
	SavePreviousField();
	RegionSmooth();
	SavePostField();
}




/*--------------------------------------------------------------------------------*/
////12/29/05

////Grow the region for any limit cycle based on Conley index theory
void GrowRegionforALimitCycle(int limitID)
{
	if(limitcycles[limitID].type == 0) //this is a repell limit cycle
		GrowRepellLimitCycle_Region(limitID);
	else
		GrowAttractLimitCycle_Region(limitID);
}


////Grow region for limit cycle with attractor characteristics
void GrowAttractLimitCycle_Region(int limitID)
{
	int i;

	attractorRegion.num = 0;
	for(i = 0; i < limitcycles[limitID].num_triangles; i++)
		attractorRegion.trianglelist[i] = limitcycles[limitID].cellcycle[i];
	attractorRegion.num = limitcycles[limitID].num_triangles;

	UpdateBoundary(1);
	GetRegionNormals(1);
	Cancel_Growing(1, -1);

}


////Grow region for limit cycle with repeller characteristics
void GrowRepellLimitCycle_Region(int limitID)
{
	int i;

	repellerRegion.num = 0;
	for(i = 0; i < limitcycles[limitID].num_triangles; i++)
		repellerRegion.trianglelist[i] = limitcycles[limitID].cellcycle[i];
	repellerRegion.num = limitcycles[limitID].num_triangles;

	UpdateBoundary(0);
	GetRegionNormals(0);
	Cancel_Growing(0, -1);
}
















/****-------------------------------------------------------------------------****/
/****-------------------------------------------------------------------------****/

/****-------------------------------------------------------------------------****/
////Routines for limit cycle relocation 3/8/06

void GetALimitCycleRegion(int limitID)
{
	int type = limitcycles[limitID].type;
	attractorRegion.num = 0;
	repellerRegion.num = 0;
	intersectRegion.num = 0;
	intersectBoundary.num = 0;
	attractorBoundary.num = 0;
	repellerBoundary.num = 0;

	SetFence_LimitCyclesexcept(limitID);

	if(type == 0) //it is a repelling limit cycle
	{
		GrowRepellLimitCycle_Region(limitID);
	}

	else{ //it is an attracting limit cycle
		GrowAttractLimitCycle_Region(limitID);
	}
}

/*
Calculate the length of the whole limit cycle, avoid repeating line segments
*/
double GetTheLengthofWholeCycle(int limitID)
{
	int num_triangles = 0;
	int pre_triangle = -1;
	int i;

	double result_len = 0;

	int num_linesegs = limitcycles[limitID].num_linesegs;
	int max_triangles = limitcycles[limitID].num_triangles;

	for(i = 0; i < num_linesegs; i++)
	{
		if(num_triangles == max_triangles)
			return result_len;

		result_len += limitcycles[limitID].closed_streamline[i].length;

		if(pre_triangle != limitcycles[limitID].closed_streamline[i].Triangle_ID)
		{
			num_triangles ++;
			pre_triangle = limitcycles[limitID].closed_streamline[i].Triangle_ID;
		}
	}

	return result_len;
}


/*
Get the interval length between those sampling controlling points
*/
double GetIntervalLength(double whole_len, int num_controlpts)
{
	return whole_len/num_controlpts;
}


/*
Get one sampling control point on the cycle
*/
void GetOneSamplePoint(int limitID, int cur_lineindex, 
					   ctr_point prectrlpt, ctr_point &curctrlpt, 
					   int &new_lineindex, double interval)
{
	int i;
	double leng = 0;

	double temp_len1, temp_len2;

	icVector2 len_vec;

	LineSeg *closedorbit = limitcycles[limitID].closed_streamline;
	int num_lines = limitcycles[limitID].num_linesegs;


	////Here we try global coordinates 3/08/06
	len_vec.entry[0] = prectrlpt.x - closedorbit[cur_lineindex].gend[0];
	len_vec.entry[1] = prectrlpt.y - closedorbit[cur_lineindex].gend[1];

	leng = length(len_vec);

	for(i = cur_lineindex+1; i < num_lines; i++)
	{
		leng += closedorbit[i].length;

		if(leng == interval)
		{
			curctrlpt.x = closedorbit[i].gend[0];
			curctrlpt.y = closedorbit[i].gend[1];
			new_lineindex = i+1;
			return;
		}

		if(leng > interval)
		{
			temp_len1 = closedorbit[i].length; //get the length of the line segment
			temp_len2 = leng - interval;       //get the extra length

			////calculate the linear interpolated coefficient
			double alpha = 1 - (temp_len2/temp_len1);
			//get the global coordinates of the control point
			curctrlpt.x = (1-alpha)*closedorbit[i].gstart[0]+alpha*closedorbit[i].gend[0];
			curctrlpt.y = (1-alpha)*closedorbit[i].gstart[1]+alpha*closedorbit[i].gend[1];
			new_lineindex = i;
			return;
		}
	}
}


/*
Sample the new control points for the limit cycle
*/
void SampleControPts(int num_controlpts, int limitID, ctr_point *control_pts)
{
	double whole_len = GetTheLengthofWholeCycle(limitID);

	double interval = GetIntervalLength(whole_len, num_controlpts);


	int i, cur_lineindex = 0;

	////set the first control point as the first point on the curve
    num_shapecontrol_pts = 0;
	control_pts[num_shapecontrol_pts].x = limitcycles[limitID].closed_streamline[0].gstart[0];
	control_pts[num_shapecontrol_pts].y = limitcycles[limitID].closed_streamline[0].gstart[1];
	num_shapecontrol_pts++;

	for(i = 1; i < num_controlpts; i++)
	{
		GetOneSamplePoint(limitID, cur_lineindex,
			control_pts[num_shapecontrol_pts-1], control_pts[num_shapecontrol_pts],
			cur_lineindex, interval);
		num_shapecontrol_pts++;
	}
}


/*
Set the boundary flags of the vertices on the boundaries of the ring-like region
*/
void SetRingRegionBoundary(int type)
{
	////Suppose we have already got the boundaries of the region
	Edge *cur_e;
	Vertex *cur_v;
	Edge **edgelist;
	int num_edges;
	int i;

	if(type == 0)
	{
		edgelist = repellerBoundary.edgelist;
		num_edges = repellerBoundary.num;
	}
	else
	{
		edgelist = attractorBoundary.edgelist;
		num_edges = attractorBoundary.num;
	}

	for(i = 0 ; i < num_edges; i++)
	{
		cur_e = edgelist[i];
		//Set the first vertices
		cur_v = Object.vlist[cur_e->verts[0]];
		cur_v->OnBoundary = 1;
		cur_v->InRegion = 0;

		cur_v = Object.vlist[cur_e->verts[1]];
		cur_v->OnBoundary = 1;
		cur_v->InRegion = 0;
	}
}


/*
Build the ring-like region for smoothing
*/
void BuildSmoothRingRegion(int type)
{
	int *triangles;
	int num_triangles;
	int i, j;
	Face *face;
	Vertex *cur_v;

	if(type == 0) //use repeller region
	{
		triangles = repellerRegion.trianglelist;
		num_triangles = repellerRegion.num;
	}
	else
	{
		triangles = attractorRegion.trianglelist;
		num_triangles = attractorRegion.num;
	}
	
	////Reset the InRegion flags for all vertices
	for(i = 0; i < Object.nverts; i++)
	{
		Object.vlist[i]->InRegion = 0;
	}

	////Get the inner vertices
	Num_verts = 0;

	for(i = 0 ; i < num_triangles; i++)
	{
		face = Object.flist[triangles[i]];

		for(j = 0; j < face->nverts; j++)
		{
			cur_v = Object.vlist[face->verts[j]];
			if(cur_v->OnBoundary == 0 && cur_v->InRegion != 1)
			{
				regionverts[Num_verts] = cur_v;
				cur_v->RegionListID = Num_verts;
				cur_v->InRegion = 1;
				Num_verts ++;
			}
		}
	}
}


/*
After move the sampling points to the new positions, we regenerate the limit cycle in the new position
using the similar method we use to generate a limit cycle before!
*/
void GeneNewLimitCycle(int type)
{

	////We need some initial process
	InitUnderneathMesh();

	////Get the triangle strip containing the limit cycle

	//GetDesignCellCycle_new3();
	//if(CalEulerValue(DesignCurveCellCycle, num_triangles_designcurve) != 0)
	//{
	//	MessageBox(NULL, "not form a ring-shaped region, try it again!", "Error", MB_OK);
	//	num_shapecontrol_pts = 0;
	//    num_curvepts_output = 0;

	//	num_cycleedges = 0;
	//	num_triangles_designcurve = 0;
	//	num_innertriangles = 0;
	//	return;
	//}

	//////Get the vertices on current boundary and extend the cell cycle from these vertices
	//ExtendDesignCellCycle();
	//if(CalEulerValue(DesignCurveCellCycle, num_triangles_designcurve) != 0)
	//{
	//	MessageBox(NULL, "not form a ring-shaped region, try it again!", "Error", MB_OK);
	//	num_shapecontrol_pts = 0;
	//    num_curvepts_output = 0;

	//	num_cycleedges = 0;
	//	num_triangles_designcurve = 0;
	//	num_innertriangles = 0;
	//	return;
	//}

	//////Set a four boundary region
 //   GetBoundary();  //Get the boundary of the new triangle strip 
	//GetNormalsForBoundaryEdges(type);


	//////set the vector values on the two inner boundaries
	//GetVectorsOnBoundaries(type);
	
	GetDesignCellCycle_new4();

	if(CalEulerValue(myCycle.DesignCellCycle, myCycle.num_triangles_designcurve) != 0)
	{
		MessageBox(NULL, "not form a ring-shaped region, try it again!", "Error", MB_OK);
		num_shapecontrol_pts = 0;
	    num_curvepts_output = 0;

		num_cycleedges = 0;
		myCycle.num_triangles_designcurve = 0;
		num_innertriangles = 0;
		return;
	}


    GetBoundaryVerts(myCycle.DesignCellCycle, myCycle.num_triangles_designcurve);

	GetVectorsOnBoundaries(Boundaryverts, num_boundverts, type);
	
	SetRingRegionBoundary(type);

	////Perform Laplacian smoothing
	BuildSmoothRingRegion(type);
	RegionSmooth();
    NormalizeVectorsOnBoundary();
	//NormalizeField();
}



////New routines for assigning vector values on the boundary of the new limit cycle

/*
Calculate the distance between a point and a straight line
*/
double GetDistance(double x1, double y1,
				   double x2, double y2,
				   double outx, double outy)
{
	double dis;
	double root;

	root = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));

	dis = fabs(((x2-x1)*(y1-outy)-(x1-outx)*(y2-y1))/root);
	return dis;
}

/*
Find all the distance from each vertex to the curve inside the triangle strip
We may suppose that all the boundary vertices are stored in the
*/
void FindAllDistance()
{
	int i;
	Edge *cur_e;
	Vertex *cur_v;

	////Reset the distance for each vertex
	for(i = 0; i < num_cycleedges; i++)
	{
		cur_e = Cycle_edge[i];

		cur_v = Object.vlist[cur_e->verts[0]];
		cur_v->distance = 1e40;
		cur_v->visited = 0;

		cur_v = Object.vlist[cur_e->verts[1]];
		cur_v->distance = 1e40;
		cur_v->visited = 0;
	}

	////Calculate all the distances
	for(i = 0; i < num_cycleedges; i++)
	{
		cur_e = Cycle_edge[i];

		cur_v = Object.vlist[cur_e->verts[0]];

		if(cur_v->distance > 1e20 && cur_v->visited == 0)
		{
			FindTheDisForOneVert(DesignCurveCellCycle, num_triangles_designcurve, cur_v->VertID);
			cur_v->visited = 1;
		}
		
		cur_v = Object.vlist[cur_e->verts[1]];
		if(cur_v->distance > 1e20 && cur_v->visited == 0)
		{
			FindTheDisForOneVert(DesignCurveCellCycle, num_triangles_designcurve, cur_v->VertID);
			cur_v->visited = 1;
		}

		/*--------------Testing code here-------------*/
		if(cur_v->distance > 1e20)
		{
			cur_v->distance = 0;
		}
	}
}

/*
Calculate the distance from one vertex to the curve 
Here let's suppose we store the line segment in the variable designcurve
*/
void FindTheDisForOneVert(int *triangles, int num_triangles, int vertid)
{
	int i, j;
	Vertex *v = Object.vlist[vertid];
	Corner *c;
	Face *face;
	double temp_dis;

	for(i = 0; i < v->Num_corners; i++)
	{
		////First, find all the triangles in the triangle list that share the vertex "vertid"
		c = Object.clist[v->Corners[i]];

		if(!IsRepeated(triangles, c->t, num_triangles))
			continue;

		////Calculate all the distances between the vertex and the line segments in the triangles
		////only store the smallest distance;

		for(j = 0; j < num_lineseg_designcurve; j++)
		{
			if(designcurve[j].Triangle_ID != c->t)
				continue;

			temp_dis = GetDistance(designcurve[j].gstart[0], designcurve[j].gstart[1],
				designcurve[j].gend[0], designcurve[j].gend[1],
				v->x, v->y);

			if(temp_dis < v->distance)
			{
				v->distance = temp_dis;
				v->which_line = j;
			}
		}
	}
}

/*
Set the vector values on the boundary
Suppose the boundary edges are stored at "Cycle_edge"
*/
void ResetVectorOnBoundary(int type)
{
	int i, j;
	Vertex *cur_v;
	Edge *cur_e;
	double largest_dis = 0;
	icVector2 line_dir, n_vec1, n_vec2;
	icVector2 eva_dis1, eva_dis2;        // evaluate the distance
	Face *face;
	int pre_triangle = -1;
	double theta;           //the rotation degree, which should be between 10-70 degree
	double magnitude = 0.1;
	
	////First, we need to find the largest distance of the vertices along the boundary
	for(i = 0; i < num_cycleedges; i++)
	{
		cur_e = Cycle_edge[i];

		cur_v = Object.vlist[cur_e->verts[0]];
		if(cur_v->distance > largest_dis)
			largest_dis = cur_v->distance;
		cur_v->visited = 0;
		
		cur_v = Object.vlist[cur_e->verts[1]];
		if(cur_v->distance > largest_dis)
			largest_dis = cur_v->distance;
		cur_v->visited = 0;
	}

	////Set the vector values according to the distance of each vertex

	for(i = 0 ; i < num_lineseg_designcurve; i++)
	{
		if(pre_triangle == designcurve[i].Triangle_ID)
			continue;


		face = Object.flist[designcurve[i].Triangle_ID];
		for(j = 0; j < face->nverts; j++)
		{
			cur_v = Object.vlist[face->verts[j]];

			if(cur_v->visited == 1) //we don't want to reset the vector values again and again
				continue;

			////Use the line segment that is most close to the vertex
			////I.e. has the smallest distance to the vertex

			if(cur_v->distance == 0)
			{
				line_dir.entry[0] = designcurve[i].gend[0]-designcurve[i].gstart[0];
				line_dir.entry[1] = designcurve[i].gend[1]-designcurve[i].gstart[1];
			}
			else
			{
				line_dir.entry[0] = designcurve[cur_v->which_line].gend[0]-designcurve[cur_v->which_line].gstart[0];
				line_dir.entry[1] = designcurve[cur_v->which_line].gend[1]-designcurve[cur_v->which_line].gstart[1];
			}
			normalize(line_dir);

			cur_v->OnBoundary = 1; //remember to set the vertex on the boundary
			cur_v->InRegion = 0;

			theta = 20 + cur_v->distance/largest_dis * 40;
			theta = (theta/180.)*(M_PI);
			magnitude = 0.005 + (cur_v->distance/largest_dis)*0.005;

			n_vec1.entry[0] = line_dir.entry[0]*cos(theta) - line_dir.entry[1]*sin(theta);
			n_vec1.entry[1] = line_dir.entry[0]*sin(theta) + line_dir.entry[1]*cos(theta);
			
			n_vec2.entry[0] = line_dir.entry[0]*cos(theta) + line_dir.entry[1]*sin(theta);
			n_vec2.entry[1] = -line_dir.entry[0]*sin(theta) + line_dir.entry[1]*cos(theta);

			if(cur_v->distance == 0)
			{
				eva_dis1.entry[0] = cur_v->x+0.005*n_vec1.entry[0] - designcurve[i].gend[0];
				eva_dis1.entry[1] = cur_v->y+0.005*n_vec1.entry[1] - designcurve[i].gend[1];

				eva_dis2.entry[0] = cur_v->x+0.005*n_vec2.entry[0] - designcurve[i].gend[0];
				eva_dis2.entry[1] = cur_v->y+0.005*n_vec2.entry[1] - designcurve[i].gend[1];
			}

			else
			{
				eva_dis1.entry[0] = cur_v->x+0.005*n_vec1.entry[0] - designcurve[cur_v->which_line].gend[0];
				eva_dis1.entry[1] = cur_v->y+0.005*n_vec1.entry[1] - designcurve[cur_v->which_line].gend[1];

				eva_dis2.entry[0] = cur_v->x+0.005*n_vec2.entry[0] - designcurve[cur_v->which_line].gend[0];
				eva_dis2.entry[1] = cur_v->y+0.005*n_vec2.entry[1] - designcurve[cur_v->which_line].gend[1];
			}

			if(type == 0)
			{
				if(length(eva_dis1) < length(eva_dis2))
				{
					cur_v->vec = magnitude*n_vec2;
				}
				else
					cur_v->vec = magnitude*n_vec1;
			}
			else
			{
				if(length(eva_dis1) > length(eva_dis2))
				{
					cur_v->vec = magnitude*n_vec2;
				}
				else
					cur_v->vec = magnitude*n_vec1;
			}

			cur_v->visited = 1;
		}
	}
}




/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
////Routines for limit cycle and saddle pair cancellation with two connecting orbits
////3/11/06

/*
judge whether there are two connecting orbits between them
We search the connection list in the saddle to find the limit cycle
*/
bool IsTwoConnections(int limitID, int saddle)
{
	int count = 0;
	int i;

	for(i = 0; i < singularities[saddle].num_connected_limitcycles; i++)
	{
		if(singularities[saddle].connected_limitcycles[i] == limitID)
			count++;
	}

	if(count == 2)
		return true;

	else
		return false;
}



/*
Initial the region for saddle
*/
void InitSaddleRegionForTwoOrbits(int limitID, int saddle, double initsaddlelength)
{
	int i;
	int other_limit1, other_limit2;
	int sep = singularities[saddle].separtices;
	int traj;

	int num_lines;
	int pre_triangle = -1;

	other_limit1 = other_limit2 = -1;

	////Let's first find the other two limit cycles the saddle connecting with
	////Here we suppose the saddle connects to other two inner limit cycles 3/11/06
	////Can it connects to two singularities inside the limit cycle ??????
	for(i = 0; i < singularities[saddle].num_connected_limitcycles; i++)
	{
		if(singularities[saddle].connected_limitcycles[i] != limitID)
		{
			if(other_limit1 < 0)
				other_limit1 = singularities[saddle].connected_limitcycles[i];
			else
				other_limit2 = singularities[saddle].connected_limitcycles[i];
		}
	}

	int *triangles1 = limitcycles[other_limit1].cellcycle;
	int num_cells1 = limitcycles[other_limit1].num_triangles;
	
	int *triangles2 = limitcycles[other_limit2].cellcycle;
	int num_cells2 = limitcycles[other_limit2].num_triangles;

	////Grow the initial triangle strip along the separatrice connecting with the other
	////two limit cycles
	if(limitcycles[limitID].type == 0) //saddle will act as an attractor
	{
		//we need to grow the initial triangle strip along the other two directions
		traj = separatrices[sep].sep1;
		num_lines = num_linesegs_curtraj[traj];
		for(i = 0; i < num_lines; i++)
		{
			if(pre_triangle == trajectories[traj][i].Triangle_ID)
				continue;

			if(IsRepeated(triangles1, trajectories[traj][i].Triangle_ID, num_cells1)
				|| IsRepeated(triangles2, trajectories[traj][i].Triangle_ID, num_cells2))
				break;

			AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 1);
			pre_triangle = trajectories[traj][i].Triangle_ID;
		}

		traj = separatrices[sep].sep3;
		num_lines = num_linesegs_curtraj[traj];
		pre_triangle = trajectories[traj][0].Triangle_ID;

		for(i = 1; i < num_lines; i++)
		{
			if(pre_triangle == trajectories[traj][i].Triangle_ID)
				continue;
			
			if(IsRepeated(triangles1, trajectories[traj][i].Triangle_ID, num_cells1)
				|| IsRepeated(triangles2, trajectories[traj][i].Triangle_ID, num_cells2))
				break;

			AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 1);
			pre_triangle = trajectories[traj][i].Triangle_ID;
		}
	}
	else //saddle will act as a repeller
	{
		//we need to grow the initial triangle strip along the other two directions
		traj = separatrices[sep].sep2;
		num_lines = num_linesegs_curtraj[traj];
		for(i = 0; i < num_lines; i++)
		{
			if(pre_triangle == trajectories[traj][i].Triangle_ID)
				continue;

			if(IsRepeated(triangles1, trajectories[traj][i].Triangle_ID, num_cells1)
				|| IsRepeated(triangles2, trajectories[traj][i].Triangle_ID, num_cells2))
				break;

			AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 0);  //modified at 06/19/06
			pre_triangle = trajectories[traj][i].Triangle_ID;
		}

		traj = separatrices[sep].sep4;
		num_lines = num_linesegs_curtraj[traj];
		pre_triangle = trajectories[traj][0].Triangle_ID;

		for(i = 1; i < num_lines; i++)
		{
			if(pre_triangle == trajectories[traj][i].Triangle_ID)
				continue;
			
			if(IsRepeated(triangles1, trajectories[traj][i].Triangle_ID, num_cells1)
				|| IsRepeated(triangles2, trajectories[traj][i].Triangle_ID, num_cells2))
				break;

			AddToRegionTriangles(trajectories[traj][i].Triangle_ID, 0);  //modified at 06/19/06
			pre_triangle = trajectories[traj][i].Triangle_ID;
		}
	}
}


void SetFenceForALimitCycle(int limitID)
{
	int i;
	for(i = 0; i < limitcycles[limitID].num_triangles; i++)
	{
		Object.flist[limitcycles[limitID].cellcycle[i]]->contain_separatrix = 1;
		Object.flist[limitcycles[limitID].cellcycle[i]]->fence_flag = 1;
	}
}

/*
Set fences for those limit cycles not involving in the limit cycle pair cancellation
*/
void SetFence_LimitCycles(int inner, int outer, int *MediaNodes, int Num_MediaNodes)
{
	int i, j;
	int numlinesegs;
	Face *face;

	for(i = 0; i < cur_limitcycle_index; i++)
	{
		if(i == inner || i == outer)
			continue;

		numlinesegs = limitcycles[i].num_linesegs;

		for(j = 0; j < numlinesegs; j++)
		{
			face = Object.flist[limitcycles[i].closed_streamline[j].Triangle_ID];
			face->contain_separatrix = 1;
			face->fence_flag = 1;
		}
	}
}

/*
Set fences for all limit cycles except the specific one
*/
void SetFence_LimitCyclesexcept(int limitID)
{
	int i, j;
	int numlinesegs;
	Face *face;

	for(i = 0; i < cur_limitcycle_index; i++)
	{
		if(i == limitID)
			continue;

		numlinesegs = limitcycles[i].num_linesegs;

		for(j = 0; j < numlinesegs; j++)
		{
			face = Object.flist[limitcycles[i].closed_streamline[j].Triangle_ID];
			face->contain_separatrix = 1;
			face->fence_flag = 1;
		}
	}
}

/*
Grow the region for the saddle according to the type of the limit cycle
*/
void GrowSaddleRegionForTwoOrbits(int limitID, int saddle)
{
	InitSaddleRegionForTwoOrbits(limitID, saddle, 1);

	UpdateBoundary(1-limitcycles[limitID].type);
	
	GetRegionNormals(1-limitcycles[limitID].type); ////Get the outward normal of the region for each boundary edge
	
	////We should avoid those two inner limit cycles
	int sep = singularities[saddle].separtices;

	if(limitcycles[limitID].type == 0)
	{
		SetFenceForASep(separatrices[sep].sep1);
		SetFenceForASep(separatrices[sep].sep3);
	}
	
	else
	{
		SetFenceForASep(separatrices[sep].sep2);
		SetFenceForASep(separatrices[sep].sep4);
	}

	////region grow processing
	Cancel_Growing(1-limitcycles[limitID].type, -1);

}

/*
Intersect the limit cycle region and saddle region, set the boundary flags
*/
void GetIntersectRegionForLimitSaddleWithTwoOrbits(int limitID, int saddle)
{
	int i, j;
	Face *face;
	Edge *cur_edge;
	Vertex *vert;

	////Add the triangle containing the saddle into the region
	Object.flist[singularities[saddle].Triangle_ID]->attract_inregion = 1;
	Object.flist[singularities[saddle].Triangle_ID]->repell_inregion = 1;

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

	for(i = 0; i < intersectRegion.num; i++)
	{
		face = Object.flist[intersectRegion.trianglelist[i]];
		for(j = 0; j < 3; j++)
		{
			cur_edge = face->edges[j];
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
	for(i = 0; i < intersectRegion.num; i++)
	{
		face = Object.flist[intersectRegion.trianglelist[i]];

		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];
			if(vert->OnBoundary == 0 && vert->InRegion == 0) ////The vertex is inside both region
			{
				intersectInnerverts.vertslist[intersectInnerverts.num] = vert;
				vert->InRegion = 1;
				vert->RegionListID = intersectInnerverts.num;     ////set the id of the vertex for future region smoothing
				intersectInnerverts.num ++;
			}
		}
	}
}

/*
Temporary routine to get one of the two connection orbits for testing 
Suppose we can guarantee that the saddle has connections with the limit cycle
Later, we should allow user to specify the separatrix he wants to reverse!!
*/
int GetOneConnection(int limitID, int saddle, int type)
{
	int sep = singularities[saddle].separtices;

	if(limitcycles[limitID].type == 0) //one of the incoming separatrices of the saddle will be chosen
	{
		return separatrices[sep].sep2;
	}

	else
	{
		return separatrices[sep].sep1;
	}
}

///**
//Set fence for all separatrices except for the specific one
//06/21/06
//**/
//void SetFenceForSeps(int except_sa)
//{
//	int i, j, traj;
//	int sep;
//
//	for(i = 0; i < cur_singularity_index; i++)
//	{
//		if(singularities[i].type == SADDLE && i != except_sa)
//		{
//			sep = singularities[i].separtices;
//			////set fence for other separatrices
//			for(j = 0; j < 4; j++)
//			{
//				switch(j){
//					case 0:
//						traj = separatrices[sep].sep1;
//						break;
//					case 1:
//						traj = separatrices[sep].sep2;
//						break;
//					case 2:
//						traj = separatrices[sep].sep3;
//						break;
//					case 3:
//						traj = separatrices[sep].sep4;
//						break;
//				}
//
//				////set fences for the separatrix of the input saddles
//				if(traj < 0)
//				{
//					int test = 0;
//					continue;
//				}
//				SetFenceForASep(traj);
//			}
//		}
//	}
//}

/*
The main routine to perform limit cycle and saddle cancellation with two connecting orbits
*/
void LimitSaddleCancelWithTwoOrbits(int limitID, int saddle)
{
	//We need to set fence before growing region 06/25/06
	SetFenceForSeps(saddle);
	SetFence_LimitCyclesexcept(limitID);

	//first, we grow the region for limit cycle
	GetALimitCycleRegion(limitID);

	//second, we need to grow the region for saddle
	GrowSaddleRegionForTwoOrbits(limitID, saddle);

	//third,  set the extra boundary if necessary
    SetExtraBoundary(limitID, saddle, GetOneConnection(limitID, saddle, 0));

	//fourth, we get the intersection of these two region
	GetIntersectRegionForLimitSaddleWithTwoOrbits(limitID, saddle);

	//finally, perform Laplacian smoothing
    //UpdateBoundary(2);
	SavePreviousField();
	Cancel_RegionSmooth();
	SavePostField();
	NormalizeField();

}

/*
Test to set the extra boundary condition, make cause jaggy flow after smoothing !!!!
*/
void SetExtraBoundary(int limitID, int saddle, int traj)
{
	Face *face;
	Vertex *vert;
	int i;

	////
	int *triangles = limitcycles[limitID].cellcycle;
	int num_cells = limitcycles[limitID].num_triangles;

	////do we need to reset the vectors on the triangles containing the selected separatrix
	////2/27/06
	int numsegs = num_linesegs_curtraj[traj];
	icVector2 line_dir, n_vec1, n_vec2;
	icVector2 eva_dis1, eva_dis2;
	double theta = (10./90.)*(M_PI / 2);
	for(i = 0; i < (int)numsegs/1.6; i++)
	{
		face = Object.flist[trajectories[traj][i].Triangle_ID];

		////We do not set the limit cycle itself as part of the extra boundary!!!
		if(IsRepeated(triangles, trajectories[traj][i].Triangle_ID, num_cells))
			break;

		if(face->attract_inregion == 1 || face->repell_inregion == 1)
		{
			face->attract_inregion = face->repell_inregion = 0;

			////Get the approximate direction of the separatrix in the triangle
			line_dir.entry[0] = trajectories[traj][i].gend[0] - trajectories[traj][i].gstart[0];
			line_dir.entry[1] = trajectories[traj][i].gend[1] - trajectories[traj][i].gstart[1];
			
			if(limitcycles[limitID].type == 1)
				line_dir *= -1;
			
			normalize(line_dir);
			
			n_vec1.entry[0] = line_dir.entry[0]*cos(theta) - line_dir.entry[1]*sin(theta);
			n_vec1.entry[1] = line_dir.entry[0]*sin(theta) + line_dir.entry[1]*cos(theta);
			
			n_vec2.entry[0] = line_dir.entry[0]*cos(theta) + line_dir.entry[1]*sin(theta);
			n_vec2.entry[1] = -line_dir.entry[0]*sin(theta) + line_dir.entry[1]*cos(theta);

			for(int j = 0; j < face->nverts; j++)
			{
				vert = Object.vlist[face->verts[j]];

				eva_dis1.entry[0] = vert->x+0.005*n_vec1.entry[0] - trajectories[traj][i].gend[0];
				eva_dis1.entry[1] = vert->y+0.005*n_vec1.entry[1] - trajectories[traj][i].gend[1];

				eva_dis2.entry[0] = vert->x+0.005*n_vec2.entry[0] - trajectories[traj][i].gend[0];
				eva_dis2.entry[1] = vert->y+0.005*n_vec2.entry[1] - trajectories[traj][i].gend[1];
				if(limitcycles[limitID].type == 0 )
				{
					////we may create a repelling limit cycle, it means that all the vectors 
					////on the boundary should separate away from the separatrix

					if(length(eva_dis1) < length(eva_dis2))
						vert->vec = 0.03*n_vec2;
					else
						vert->vec = 0.03*n_vec1;
				}
				else if(limitcycles[limitID].type == 1)
				{
					////we may create an attractor limit cycle, it means that all the vectors
					////on the boundary should converge to the separatrix

					if(length(eva_dis1) < length(eva_dis2))
						vert->vec = 0.03*n_vec1;
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
