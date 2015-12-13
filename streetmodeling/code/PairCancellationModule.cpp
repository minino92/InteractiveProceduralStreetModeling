////PairCancellationModule.cpp
#include "stdafx.h"

#include "PairCancellationModule.h"

#include "vfdatastructure.h"
#include "topologyedit.h"

#include "limitcycleedit.h"
#include "ConleyRelationGraph.h"
#include "VFSynthesis.h"
#include "limitcycleedit.h"


extern Singularities *singularities;           //being captured singularites' list
extern int cur_singularity_index;
extern LineSeg **trajectories;                 //trajectories' list
extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory
extern Separatrices *separatrices;             //array for group of separatrices
extern int cur_separatrices_index;
extern LimitCycle *limitcycles;                 //limit cycles data structure
extern int cur_limitcycle_index;


////variables for singularities pair cancellation and movement
extern TriangularRegion repellerRegion;       ////region containing a repeller
extern TriangularRegion attractorRegion;      ////region containing an attractor
extern TriangularRegion intersectRegion;      ////region for intersection
extern RegionBoundary repellerBoundary;
extern RegionBoundary attractorBoundary;
extern InnerVertices  repellerInnerverts;
extern InnerVertices  attractorInnerverts;
extern TriangularRegion intersectRegion;     ////The intersect region
extern RegionBoundary intersectBoundary;     ////The intersect boundary 
extern InnerVertices intersectInnerverts;    ////The inner vertices inside the intersect region

extern TriangularRegion Source_re1, Source_re2, Source_re;
extern TriangularRegion Sink_re1, Sink_re2, Sink_re;

extern GraphNode *graphnodes ;
extern GraphEdge *graphedges ;
extern int cur_node_index;
extern int cur_graphedge_index;
extern int *MediaNodes;
extern int Num_MediaNodes;
extern int MaxMediaNodes;

extern Polygon3D Object;

/**************************/
double percentage = 1.;

extern void Undo();
extern void ReCalDiretionalVec();
extern void UnionRegion(TriangularRegion &source1, TriangularRegion &source2, TriangularRegion &dest);
extern void IntersectRegion(TriangularRegion &source1, TriangularRegion &source2, TriangularRegion &dest);
extern void SetBoundaryFlag_Ver(Edge **boundaryedgelist, int num_edges);



/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/
/*
A general structure for all pair cancellations
Input: 
      repeller -- the index of a singularity or a limit cycle acting as a repeller
	  attractor -- the index of a singularity or a limit cycle acting as an attractor
	            two nodes in C-graph corresponding to the two being cancelled elements
	  type      the type of cancellation
	            {0 -- singularity pair cancellation;
				 1 -- limit cycle + singularity pair cancellation;
				 2 -- limit cycle pair cancellation
				 }
*/

void Gen_PairCancellation(int repeller, int attractor, int type)
{
	if(type == 0)          //it is singularity pair cancellation
	{
		SingularityPairCancel(repeller, attractor);
	}

	else if(type == 1)     //it is singularity and limit cycle pair cancellation
	{
		// for this case, we assume that user pass singularity through first parameter /augment
		// and pass the index of limit cycle through seconde augment
		LimitCycleSingularityCancel(repeller, attractor);
	}

	else
	{
	}
}


/*
routine for singularity pair cancellation
include the simply directly connected cases, doubley connected cases and indirectly connected cases
*/
void SingularityPairCancel(int repeller, int attractor)
{
	int sing1, sing2;

	sing1 = graphnodes[repeller].singularityID;
	sing2 = graphnodes[attractor].singularityID;

	//directly connected
	if(singularities[sing1].type == SADDLE || singularities[sing2].type == SADDLE)
	{
		DirectlyConnectedSingularityPair(sing1, sing2);
	}

	else //indireclty connected
	{
		IndirectlyConnectedSingularityPair(sing1, sing2);
	}
}

/*
Routine for directly connected singularity pair cancellation
*/
void DirectlyConnectedSingularityPair(int sing1, int sing2)
{
	//
	if(singularities[sing1].type == SADDLE)
	{
		if(TwoConnectOrbits(sing1, sing2))
		{
			// doubley connected with sing1 is a saddle and sing2 is other kind of singuarity
			DoublyConnectedSingularityPair(sing1, sing2);
		}
	}

	else
	{
		if(TwoConnectOrbits(sing2, sing1))
		{
			// doubley connected with sing2 is a saddle and sing1 is other kind of singuarity
			DoublyConnectedSingularityPair(sing2, sing1);
		}
	}
}


/*
Doubly connected singularity pair cancellation
Input: saddle --  the index of the saddle in the singularity list
	   othersing -- the index of the other singularity in the list
*/

void DoublyConnectedSingularityPair(int saddle, int othersing)
{
	int type;
	int repellID, attractID;
	int repell_triangle, attract_triangle;

	if(singularities[othersing].type == SOURCE || singularities[othersing].type == RFOCUS)
		type = 0; //repeller, thus, saddle should act as attractor

	else if(singularities[othersing].type == SINK || singularities[othersing].type == AFOCUS)
		type = 1; //attractor, thus, saddle should act as repeller
	////we do not consider other case now

	if(type == 0)
	{
		repellID = othersing;
		repell_triangle = singularities[othersing].Triangle_ID;

		attractID = saddle;
		attract_triangle = singularities[saddle].Triangle_ID;
	}

	else{
		attractID = othersing;
		attract_triangle = singularities[othersing].Triangle_ID;

		repellID = saddle;
		repell_triangle = singularities[saddle].Triangle_ID;
	}
	
	////For setting the fences and updating the C-graph
	int repellers[1] = {singularities[repellID].node_index};
	int num_repellers = 1;
	int attractors[1] = {singularities[attractID].node_index};
	int num_attractors = 1;

 	SetFences(repellers, num_repellers,
		attractors, num_attractors,
		MediaNodes, Num_MediaNodes);

//we need to control the initial length of the saddle region here
	GrowAttractRegionForTwoConnections_test(attractID, repell_triangle, percentage);
	GrowRepellRegionForTwoConnections_test(repellID ,attract_triangle, percentage);

	int traj = GetOneConnection(saddle, othersing);
	icVector2 saddle_dir = GetSaddleTriangleDirection(saddle, traj);

	//SetExtraBoundary(saddle, singID, saddle_dir, traj);
	SetExtraBoundary2(saddle, othersing, saddle_dir, traj);

	GetIntersectRegionForTwoConnections(repellID, attractID);
	UpdateBoundary(2);

	////At this moment, we just randomly pick a connection, later we can add
	////the interface for user to specify which connection he wants

	//Cancel_RegionSmooth();

	//ProjectToTangentPlane();

	//MarkCancel(repellers,1, attractors,1, MediaNodes, Num_MediaNodes);
}

/* we need to modify the following routine to pass the percentage !!
Therefore, we also need to modify the two routines
	GrowAttractRegionForTwoConnections(attractID, repell_triangle, +percentage);
	GrowRepellRegionForTwoConnections(repellID ,attract_triangle, +percentage);
*/
void InitSaddleRegionForTwoConnections_test(int saddle, int type, int connect_sing, 
											int initsaddlelength, double percentage)
{
	int i, j, sep, traj;

	////it seems that we just need to add one triangle that contains the saddle
	int counter=0;
	if(type == 0) //saddle acts as an attractor
	{
		attractorRegion.trianglelist[0] = singularities[saddle].Triangle_ID;
		attractorRegion.num = 1;
		Object.flist[singularities[saddle].Triangle_ID]->attract_inregion = 1;

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

			////We just add the triangles containing the incoming separatrices 06/22/06
			
			if(i == 0 || i == 2)
			{
				//add only one more triangle along these directions
				////we just add one more at each direction and see what happen
				if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID!=singularities[connect_sing].Triangle_ID)
				{
					////use the 'percentage' to control the initial length of the separatrices
					for(j = 0; j < (int)num_linesegs_curtraj[traj]*percentage; j++)
					{
						if(trajectories[traj][j].Triangle_ID != singularities[saddle].Triangle_ID)
						{
							if(Object.flist[trajectories[traj][j].Triangle_ID]->contain_singularity != 1)
							{
								AddToRegionTriangles(trajectories[traj][j].Triangle_ID, 1);
								//break;
							}
						}
					}
				}

				continue;
			}

			////we just add one more at each direction and see what happen
			if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID==singularities[connect_sing].Triangle_ID)
			{
				for(j = 0; j < num_linesegs_curtraj[traj]; j++)
				{
					if(trajectories[traj][j].Triangle_ID != singularities[saddle].Triangle_ID)
					{
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

			////We just add the triangles containing the outgoing separatrices 06/22/06
			if(i == 1 || i == 3)
			{
				////use the 'percentage' to control the initial length of the separatrices
				for(j = 0; j < (int)num_linesegs_curtraj[traj]*percentage; j++)
				{
					if(trajectories[traj][j].Triangle_ID != singularities[saddle].Triangle_ID)
					{
						if(Object.flist[trajectories[traj][j].Triangle_ID]->contain_singularity != 1)
						{
							AddToRegionTriangles(trajectories[traj][j].Triangle_ID, 0);
							//break;
						}
					}
				}

				continue;
			}

			////we just add one more at each direction and see what happen
			if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID==singularities[connect_sing].Triangle_ID)
			{
				for(j = 0; j < num_linesegs_curtraj[traj]; j++)
				{
					if(trajectories[traj][j].Triangle_ID != singularities[saddle].Triangle_ID)
					{
						if(Object.flist[trajectories[traj][j].Triangle_ID]->contain_singularity != 1)
							AddToRegionTriangles(trajectories[traj][j].Triangle_ID, 0);
					}
				}
			}
		}

		UpdateBoundary(0);
	}
}


void GrowAttractRegionForTwoConnections_test(int singularID, int target_triangle, double percentage)
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
		
		////
		attractorBoundary.edgelist[0]->OnAttractBoundary = 1;
		attractorBoundary.edgelist[1]->OnAttractBoundary = 1;
		attractorBoundary.edgelist[2]->OnAttractBoundary = 1;
	}
	else
	{
		InitSaddleRegionForTwoConnections_test(singularID, 0, 
			Object.flist[target_triangle]->singularity_index, 1, percentage);
	}


	GetRegionNormals(1); ////Get the outward normal of the region for each boundary edge

	////initialize the inner vertices list
	attractorInnerverts.num = 0;    ////no inner vertex at present

	////region grow processing
	Cancel_Growing(1, target_triangle);
}



void GrowRepellRegionForTwoConnections_test(int singularID, int target_triangle, double percentage)
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
		
		////
		repellerBoundary.edgelist[0]->OnRepellBoundary = 1;
		repellerBoundary.edgelist[1]->OnRepellBoundary = 1;
		repellerBoundary.edgelist[2]->OnRepellBoundary = 1;
	}
	
    else
	{
		InitSaddleRegionForTwoConnections_test(singularID, 1, 
			Object.flist[target_triangle]->singularity_index, 1, percentage);
	}


	GetRegionNormals(0); ////Get the outward normal of the region for each boundary edge

	////initialize the inner vertices list
	repellerInnerverts.num = 0;    ////no inner vertex at present

	////region grow processing
	Cancel_Growing(0, target_triangle);
}




/*
Deal with the indirectly connected singularity pair cancellation
*/
void IndirectlyConnectedSingularityPair(int sing1, int sing2)
{
	//we need to search the intermediate components (saddles)
	MediaNodes = (int *)malloc(sizeof(int)*MaxMediaNodes);
	Num_MediaNodes = 0;

	//suppose we have already initialize the cancellation
	int temp_repeller[1], temp_attractor[1];

	if(singularities[sing1].type == SOURCE || singularities[sing1].type == RFOCUS)  //sing1 -- repeller, sing2 -- attractor
	{
		temp_repeller[0] = singularities[sing1].node_index;
		temp_attractor[0] = singularities[sing2].node_index;
	}
	else   //sing1 -- attractor, sing2 -- repeller
	{
		temp_attractor[0] = singularities[sing1].node_index;
		temp_repeller[0] = singularities[sing2].node_index;
	}

	//search the intermediate components (saddles)
	SearchConnectComponents_adv(temp_repeller, 1, temp_attractor, 1, MediaNodes, Num_MediaNodes);

	if(Num_MediaNodes == 1)  //the sink-source with double connections between their paths
	{
		IndirectlyDoublyConnectedSingularityPair(sing1, sing2, MediaNodes, Num_MediaNodes);
	}
	else{
		// Normal indirectly connected singularity pair

		//set fences according to the intermediate saddles
		SetFenceForSeps(-1);
		SetFence_LimitCyclesexcept(-1);
	}
}

/*
Deal with the indirectly connected cases but there exists double connections between the intermediate saddle
and one of the being cancelled singularities
Here we assume there is only one intermediate saddle
Input: sing1, sing2 -- two different non-saddle singularities
*/
void 
IndirectlyDoublyConnectedSingularityPair(int sing1, int sing2, int *Media, int num_media)
{
	//suppose that we have already judged the existing of double connections between 
	//one saddle and one singularity
	//suppose we have already initialize the cancellation
	int temp_repeller[1], temp_attractor[1];

	if(singularities[sing1].type == SOURCE || singularities[sing1].type == RFOCUS)  
	{
		temp_repeller[0] = singularities[sing1].node_index;
		temp_attractor[0] = singularities[sing2].node_index;
	}
	else if(singularities[sing1].type == SINK || singularities[sing1].type == AFOCUS)
	{
		temp_attractor[0] = singularities[sing1].node_index;
		temp_repeller[0] = singularities[sing2].node_index;
	}

	//set fences according to the intermediate saddles
	SetFence_LimitCyclesexcept(-1);
	SetFenceForSeps(Media[0]);
	RemoveFenceForIndirectlyDoublyPair(sing1, sing2);

	if(TwoConnectOrbits(Media[0], sing1))
	{
		if(!Ada_Getregion_IndirectlyDoublyConnectedSingularityPair(sing1, sing2, Media[0]))
		{
			return;
		}
	}

	else if(TwoConnectOrbits(Media[0], sing2))
	{
		if(!Ada_Getregion_IndirectlyDoublyConnectedSingularityPair(sing2, sing1, Media[0]))
		{
			return;
		}
	}

	////Get the inner vertices of the intersect region
	GetInnerVerts(sing1, Media, num_media);


	//perform one smoothing
	Cancel_RegionSmooth();

	//Judge whether the result field satisfies the Conley/Poincare index
	int totalindex = 0;
	if(IntersectedRegionSingCount(totalindex) != abs(2 - num_media))
	{
		Undo();
		return;
	}

	MarkCancel(temp_repeller, 1, temp_attractor, 1, Media, num_media);
}


/*
Assume that we always let 'sing1' have double connections with the 'inter_saddle'
*/

bool 
Ada_Getregion_IndirectlyDoublyConnectedSingularityPair(int sing1, int sing2, 
													   int inter_saddle)
{
	double small_length, large_length, cur_length;
	int min_trianglesnum;
	int success, count = 0;

	////Note that here, Source_re1 and Sink_re1 can be local variables! 11/30/05
	Source_re.num = 0;
	Source_re1.num = 0;
	Sink_re.num = 0;
	Sink_re1.num = 0;

	/*  Adaptive step  */
	small_length = 0;
	large_length = 2.;  //what is the best value ?
	min_trianglesnum = GetMinTriangleNumofSepsForDoublyConnectedSingularityPair(inter_saddle, sing1, sing2);

	success = 0;

	do{
		InitCancellationAndMovement(); //initialize all the flags except fences
		Source_re.num = 0;
		Source_re1.num = 0;
		Sink_re.num = 0;
		Sink_re1.num = 0;

		if(singularities[sing1].type == SOURCE || singularities[sing1].type == RFOCUS)
		{
			//Grow repelling region from sing1
			Cancel_GrowRepellerRegion(singularities[sing1].Triangle_ID, -1, sing1, InitSaddleRegionLength,
				0, 1.);

			//Grow attracting region from sing2
			Cancel_GrowAttractorRegion(singularities[sing2].Triangle_ID, -1, sing2, InitSaddleRegionLength,
				0, 1.);
		}

		else{
			//Grow attracting region from sing1
			Cancel_GrowAttractorRegion(singularities[sing1].Triangle_ID, -1, sing1, InitSaddleRegionLength,
				0, 1.);

			//Grow repelling region from sing2
			Cancel_GrowRepellerRegion(singularities[sing2].Triangle_ID, -1, sing2, InitSaddleRegionLength,
				0, 1.);
		}

		cur_length = (small_length + large_length)/2.;

		//Grow region from each saddle as attractors
		UnionRegion(attractorRegion, Sink_re1, Sink_re);

			//after growing a region for one saddle, we need to union it with previous saddle regions
			attractorRegion.num = 0;

			//Initialize the saddle region
			//InitSaddle_IndirectlyDoublyConnectedSingularityPair(inter_saddle,
			//	sing1, sing2, 1, percentage); //here 'percentage' is the global one!
			InitSaddle_IndirectlyDoublyConnectedSingularityPair(inter_saddle,
				sing1, sing2, 1, cur_length); //here 'cur_length' is the current percentage!

			//grow region from the initial region here
			UpdateBoundary(1);
			
			GetRegionNormals(1); ////Get the outward normal of the region for each boundary edge
			Cancel_Growing(1, -1);

			CopyRegion(Sink_re.trianglelist, Sink_re1.trianglelist, Sink_re.num);
			Sink_re1.num = Sink_re.num;
			UnionRegion(attractorRegion, Sink_re1, Sink_re);

		//Grow region from each saddle as repellers
		UnionRegion(repellerRegion, Source_re1, Source_re);

			//after growing a region for one saddle, we need to union it with previous saddle regions
			repellerRegion.num = 0;

			//Initialize the saddle region 
			//InitSaddle_IndirectlyDoublyConnectedSingularityPair(inter_saddle,
			//	sing1, sing2, 0, percentage); //here 'percentage' is the global one!
			InitSaddle_IndirectlyDoublyConnectedSingularityPair(inter_saddle,
				sing1, sing2, 0, cur_length); //here 'cur_length' is the current percentage!
			
			//grow region from the initial region here
			UpdateBoundary(0);
			
			GetRegionNormals(0); ////Get the outward normal of the region for each boundary edge
			Cancel_Growing(0, -1);

			CopyRegion(Source_re.trianglelist, Source_re1.trianglelist, Source_re.num);
			Source_re1.num = Source_re.num;
			UnionRegion(repellerRegion, Source_re1, Source_re);

		//get the intersection of the two region
		IntersectRegion(Source_re, Sink_re, intersectRegion);
		
		//Judge whether the intersection region satisfies the Conley boundary condition or not
		if(CalEulerValue(intersectRegion.trianglelist, intersectRegion.num) == 1) //a disc-shaped
		{
			success = 1; //we do find a valid region to perform smoothing
			//success_length = cur_length;
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

		count++;

	}while(large_length - small_length > 1e-5 && (double)min_trianglesnum*(large_length - small_length) > 0.5
		&& (int)min_trianglesnum*cur_length >= 1);

	if(success == 1)
		return true;
    
	return false;

}


/*
Assume sing1 has double connections with saddle
*/
int GetMinTriangleNumofSepsForDoublyConnectedSingularityPair(int saddle, int sing1, int sing2)
{
	int sep = singularities[saddle].separtices;
	int traj, result;

	if(singularities[sing1].type == SOURCE || singularities[sing1].type == RFOCUS)
	{
		//get the number of triangles of the outgoing separatrix not connecting to the 'sing2'

		traj = separatrices[sep].sep1;
		if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID != singularities[sing2].Triangle_ID)
		{
			traj = separatrices[sep].sep3;
		}
	}

	else if(singularities[sing1].type == SINK || singularities[sing1].type == AFOCUS)
	{
		//get the number of triangles of the incoming separatrix not connecting to the 'sing2'

		traj = separatrices[sep].sep2;
		if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID != singularities[sing2].Triangle_ID)
		{
			traj = separatrices[sep].sep4;
		}
	}

	result = num_linesegs_curtraj[traj];
	return result;
}


/*
Remove the fences from the separatrices that reach 'sing1' or 'sing2'
*/
void RemoveFenceForIndirectlyDoublyPair(int sing1, int sing2)
{
	int i, j, k;
	int sep, traj;
	int pre_triangle = -1;

	for(i = 0; i < cur_singularity_index; i++)
	{
		if(singularities[i].type != SADDLE)	
			continue;

		sep = singularities[i].separtices;

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

			if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID != singularities[sing1].Triangle_ID
				&& trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID != singularities[sing2].Triangle_ID)
				continue;

			pre_triangle = -1;
			for(k = 0; k < num_linesegs_curtraj[traj]; k++)
			{
				if(trajectories[traj][k].Triangle_ID != pre_triangle)
				{
					Object.flist[trajectories[traj][k].Triangle_ID]->fence_flag = 0;
					pre_triangle = trajectories[traj][k].Triangle_ID;
				}
			}
		}
	}
}


/*
Assume that we always let 'sing1' have double connections with the 'saddle'
*/
void
InitSaddle_IndirectlyDoublyConnectedSingularityPair(int saddle, int sing1, int sing2,
													int type, double percentage)
{
	int i;
	int sep = singularities[saddle].separtices;
	int traj, other_traj;
	int pre_triangle = -1;

	//Note that since saddle and sing1 have double connections

	if(singularities[sing1].type == SOURCE || singularities[sing1].type == RFOCUS)
	{
		//add the triangles containing the two incoming separatrices
		traj = separatrices[sep].sep2;
		for(i = 0; i < num_linesegs_curtraj[traj]; i++)
		{
			if(trajectories[traj][i].Triangle_ID != pre_triangle)
			{
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, type);
				pre_triangle = trajectories[traj][i].Triangle_ID;
			}
		}

		pre_triangle = -1;
		traj = separatrices[sep].sep4;
		for(i = 0; i < num_linesegs_curtraj[traj]; i++)
		{
			if(trajectories[traj][i].Triangle_ID != pre_triangle)
			{
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, type);
				pre_triangle = trajectories[traj][i].Triangle_ID;
			}
		}

		//add the triangles containing the outgoing separatrix connecting to the attracting 'sing2'
		traj = separatrices[sep].sep1;
		other_traj = separatrices[sep].sep3;
		if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID != singularities[sing2].Triangle_ID)
		{
			other_traj = traj;
			traj = separatrices[sep].sep3;
		}
		
		pre_triangle = -1;
		for(i = 0; i < num_linesegs_curtraj[traj]; i++)
		{
			if(trajectories[traj][i].Triangle_ID != pre_triangle)
			{
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, type);
				pre_triangle = trajectories[traj][i].Triangle_ID;
			}
		}

		//finally, add the left separatrix according to the optimization length
		pre_triangle = -1;
		for(i = 0; i < (int)num_linesegs_curtraj[other_traj]*percentage; i++)
		{
			//remember: do not contain the singularity other than 'sing' and the intermediate saddles
			if(Object.flist[trajectories[other_traj][i].Triangle_ID]->contain_singularity == 1
				&& Object.flist[trajectories[other_traj][i].Triangle_ID]->singularity_index != sing1
				&& Object.flist[trajectories[other_traj][i].Triangle_ID]->singularity_index != sing2 )
				continue ;

			if(trajectories[other_traj][i].Triangle_ID != pre_triangle)
			{
				AddToRegionTriangles(trajectories[other_traj][i].Triangle_ID, type);
				pre_triangle = trajectories[other_traj][i].Triangle_ID;
			}
		}
	}

	else
	{
		//add the triangles containing the two outgoing separatrices
		traj = separatrices[sep].sep1;
		for(i = 0; i < num_linesegs_curtraj[traj]; i++)
		{
			if(trajectories[traj][i].Triangle_ID != pre_triangle)
			{
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, type);
				pre_triangle = trajectories[traj][i].Triangle_ID;
			}
		}

		pre_triangle = -1;
		traj = separatrices[sep].sep3;
		for(i = 0; i < num_linesegs_curtraj[traj]; i++)
		{
			if(trajectories[traj][i].Triangle_ID != pre_triangle)
			{
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, type);
				pre_triangle = trajectories[traj][i].Triangle_ID;
			}
		}

		//add the triangles containing the incoming separatrix connecting to the attracting 'sing2'
		traj = separatrices[sep].sep2;
		other_traj = separatrices[sep].sep4;
		if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID != singularities[sing2].Triangle_ID)
		{
			other_traj = traj;
			traj = separatrices[sep].sep4;
		}
		
		pre_triangle = -1;
		for(i = 0; i < num_linesegs_curtraj[traj]; i++)
		{
			if(trajectories[traj][i].Triangle_ID != pre_triangle)
			{
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, type);
				pre_triangle = trajectories[traj][i].Triangle_ID;
			}
		}

		//finally, add the left separatrix according to the optimization length
		pre_triangle = -1;
		for(i = 0; i < (int)num_linesegs_curtraj[other_traj]*percentage; i++)
		{
			//remember: do not contain the singularity other than 'sing' and the intermediate saddles
			if(Object.flist[trajectories[other_traj][i].Triangle_ID]->contain_singularity == 1
				&& Object.flist[trajectories[other_traj][i].Triangle_ID]->singularity_index != sing1
				&& Object.flist[trajectories[other_traj][i].Triangle_ID]->singularity_index != sing2 )
				continue ;

			if(trajectories[other_traj][i].Triangle_ID != pre_triangle)
			{
				AddToRegionTriangles(trajectories[other_traj][i].Triangle_ID, type);
				pre_triangle = trajectories[other_traj][i].Triangle_ID;
			}
		}
	}
}


/*
*/
void RemoveFenceFromSepToCycle(int cycle)
{
	//Remove the fence of the separatrice that reach the limit cycle
	int i, j;
	int sep, traj;
	for(i = 0; i < cur_singularity_index; i++)
	{
		if(singularities[i].type == SADDLE && !IsRepeated(MediaNodes, singularities[i].node_index, Num_MediaNodes))
		{
			sep = singularities[i].separtices;
			if(limitcycles[cycle].type == 0)
			{
				traj = separatrices[sep].sep2;
				if(!IsRepeated(limitcycles[cycle].cellcycle,
					trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID, limitcycles[cycle].num_triangles))
				{
					traj = separatrices[sep].sep4;
					if(!IsRepeated(limitcycles[cycle].cellcycle,
						trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID, limitcycles[cycle].num_triangles))
						continue; //this saddle has no separatrices connecting to the limit cycle)
				}
			}

			else{
				traj = separatrices[sep].sep1;
				if(!IsRepeated(limitcycles[cycle].cellcycle,
					trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID, limitcycles[cycle].num_triangles))
				{
					traj = separatrices[sep].sep3;
					if(!IsRepeated(limitcycles[cycle].cellcycle,
						trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID, limitcycles[cycle].num_triangles))
						continue; //this saddle has no separatrices connecting to the limit cycle)
				}
			}
			
			for(j = 0; j < num_linesegs_curtraj[traj]; j++)
			{
				Object.flist[trajectories[traj][j].Triangle_ID]->fence_flag = 0;
			}
		}
	}
}


/*
Get the minimum number of triangles of the separatrices starting from the intermediate saddles but
not connecting to the specific singularity and limit cycle
*/
int GetMinTriangleNumofSeps(int *Media, int num, int sing, int cycle)
{
	int i, j;
	int sep, traj;
	int min;
	min = Object.nfaces;

	for(i = 0; i < num; i++)
	{
		sep = singularities[graphnodes[Media[i]].node_index].separtices;

		for(j = 0; j < 4; j++)
		{
			switch(j)
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
			
			//if the separatrix connects to the singularity or the limit cycle, skip it!
			if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID == singularities[sing].Triangle_ID)
				continue;

			if(IsRepeated(limitcycles[cycle].cellcycle, 
				trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID, limitcycles[cycle].num_triangles))
				continue;

			//get the maximum and minimum
			if(num_linesegs_curtraj[traj] < min) min = num_linesegs_curtraj[traj];
		}
	}

	return min;
}


/**
Set fence for all separatrices except for the specific one
06/21/06
**/
void SetFenceForSeps(int except_sa)
{
	int i, j, traj;
	int sep;

	for(i = 0; i < cur_singularity_index; i++)
	{
		if(singularities[i].type == SADDLE && i != except_sa)
		{
			sep = singularities[i].separtices;
			////set fence for other separatrices
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

				////set fences for the separatrix of the input saddles
				if(traj < 0) //?
				{
					int test = 0;
					continue;
				}
				SetFenceForASep(traj);
			}
		}
	}
}


////Remove the fence for limit cyclc involved cancellation 06/19/06

void RemoveLimitCycleFence(int limitID)
{
	int i;
	int triangle;

	for(i = 0; i < limitcycles[limitID].num_triangles; i++)
	{
		triangle = limitcycles[limitID].cellcycle[i];
		Object.flist[triangle]->fence_flag = 0;
	}
}


void RemoveFenceForSaddle(int saddle)
{
	int traj;
	int j, i;

	for(j = 0; j < 4; j ++)
	{
		switch(j){
			case 0:
				traj = separatrices[singularities[saddle].separtices].sep1;
				break;
			case 1:
				traj = separatrices[singularities[saddle].separtices].sep2;
				break;
			case 2:
				traj = separatrices[singularities[saddle].separtices].sep3;
				break;
			case 3:
				traj = separatrices[singularities[saddle].separtices].sep4;
				break;
		}

		////set fences for the separatrix of the input saddles
		if(traj < 0)
		{
			int test = 0;
			continue;
		}

		for(i = 0; i < num_linesegs_curtraj[traj]; i++)
		{
			Object.flist[trajectories[traj][i].Triangle_ID]->fence_flag = 0;
		}
	}
}



/* functions for different kinds of singularity and limit cycle pair cancellation */

bool 
IndirectlyConnectedSingandCyclePair(int sing, int cycle)
{
	//MediaNodes = (int *)malloc(sizeof(int)*MaxMediaNodes);
	//Num_MediaNodes = 0;

	//suppose we have already initialize the cancellation
	int temp_repeller[1], temp_attractor[1];

	if(limitcycles[cycle].type == 0)  //limit cycle -- repeller, singularity -- attractor
	{
		temp_repeller[0] = limitcycles[cycle].node_index;
		temp_attractor[0] = singularities[sing].node_index;
	}
	else   //limit cycle -- attractor, singularity -- repeller
	{
		temp_attractor[0] = limitcycles[cycle].node_index;
		temp_repeller[0] = singularities[sing].node_index;
	}

	//search the intermediate components (saddles)
	SearchConnectComponents_adv(temp_repeller, 1, temp_attractor, 1, MediaNodes, Num_MediaNodes);

	//set fences according to the intermediate saddles
	SetFenceForSeps(-1);
	SetFence_LimitCyclesexcept(cycle);
	RemoveLimitCycleFence(cycle);
	RemoveFenceFromSepToCycle(cycle);

	for(int i = 0; i < Num_MediaNodes; i++)
	{
		RemoveFenceForSaddle(graphnodes[MediaNodes[i]].singularityID);
	}


	if(Num_MediaNodes == 1)
	{
		//doubly connected pair!? 06/27/06 yes, it can happen
	}

	if(!Ada_Getregion_IndirectlyConnectedSingandCyclePair(sing, cycle, MediaNodes, Num_MediaNodes))
	{
		//free(MediaNodes);
		//Num_MediaNodes = 0;
		return false;
	}

	/* complementary process as follows */
	////Add the triangles that contain the saddles into the region
	int cur_t;
	for(int i = 0; i < Num_MediaNodes; i++)
	{
		cur_t = singularities[graphnodes[MediaNodes[i]].singularityID].Triangle_ID;

		if(!IsRepeated(intersectRegion.trianglelist, cur_t, intersectRegion.num))
		{
			intersectRegion.trianglelist[intersectRegion.num] = cur_t;
			intersectRegion.num++;
		}
	}

	////Add the being wanted to cancel singularity into the region
	cur_t = singularities[sing].Triangle_ID;

	if(!IsRepeated(intersectRegion.trianglelist, cur_t, intersectRegion.num))
	{
		intersectRegion.trianglelist[intersectRegion.num] = cur_t;
		intersectRegion.num++;
	}

	////Get the inner vertices of the intersect region
	GetInnerVerts(sing, MediaNodes, Num_MediaNodes);


	//perform one smoothing
	Cancel_RegionSmooth();

	//Judge whether the result field satisfies the Conley/Poincare index
	int totalindex = 0;
	if(!IntersectedRegionSingCount(totalindex) == abs(1 - Num_MediaNodes))
	{
		Undo();
		//free(MediaNodes);
		//Num_MediaNodes = 0;
		return false;
	}

	//MarkCancel(temp_repeller, 1, temp_attractor, 1, MediaNodes, Num_MediaNodes);
	//
	//free(MediaNodes);
	//Num_MediaNodes = 0;
	return true;
}


/*
We assume that the singularity and the limit cycle pair satisfy the basic requirement of cancellation
The shape of the smoothing region should be a ring
Modified at 06/29/06, Initial all the cancellation related flags in each adaptive step
*/
bool 
Ada_Getregion_IndirectlyConnectedSingandCyclePair
(int sing, int cycle, int *Media, int num_media)
{
	//first grow regions for the singularity and limit cycle respectively. These two regions are constant
	//during the optimization
	int i;
	int cur_saddle;

	double small_length, large_length, cur_length;
	int min_trianglesnum;
	int success, count = 0;
	
	success = 0;

	////Note that here, Source_re1 and Sink_re1 can be local variables! 11/30/05
	Source_re.num = 0;
	Source_re1.num = 0;
	Sink_re.num = 0;
	Sink_re1.num = 0;

	/*  Adaptive step  */
	small_length = 0;
	//large_length = 0.6;
	large_length = 2.;  //what is the best value ?
	min_trianglesnum = GetMinTriangleNumofSeps(Media, num_media, sing, cycle);


	do{
		InitCancellationAndMovement();
		Source_re.num = 0;
		Source_re1.num = 0;
		Sink_re.num = 0;
		Sink_re1.num = 0;


		if(limitcycles[cycle].type == 0) //limit cycle acts as a repeller
		{
			// Grow region for limit cycle (repeller)
			GrowRegionforALimitCycle(cycle);

			////Grow region for the singularity (attractor)
			Cancel_GrowAttractorRegion(singularities[sing].Triangle_ID, -1, sing, InitSaddleRegionLength,
				0, 1.);
		}

		else  //limit cycle acts as an attractor
		{
			// Grow region for limit cycle (attractor)
			GrowRegionforALimitCycle(cycle);

			// Grow region for the singularity (repeller)
			Cancel_GrowRepellerRegion(singularities[sing].Triangle_ID, -1, sing, InitSaddleRegionLength,
				0, 1.);
		}

		cur_length = (small_length + large_length)/2.;

		//Grow region from each saddle as attractors
		UnionRegion(attractorRegion, Sink_re1, Sink_re);

		for(i = 0; i < num_media; i++)
		{

			//after growing a region for one saddle, we need to union it with previous saddle regions
			attractorRegion.num = 0;
			cur_saddle = graphnodes[MediaNodes[i]].singularityID;

			//Initialize the saddle region
			//InitSaddle_IndirectlyConnectedSingandCycle(cur_saddle,
			//	1, sing, cycle, percentage); //here 'percentage' is the global one!
			InitSaddle_IndirectlyConnectedSingandCycle(cur_saddle,
				1, sing, cycle, cur_length); //here 'cur_length' is the current percentage!

			//grow region from the initial region here
			UpdateBoundary(1);
			
			GetRegionNormals(1); ////Get the outward normal of the region for each boundary edge
			Cancel_Growing(1, -1);

			CopyRegion(Sink_re.trianglelist, Sink_re1.trianglelist, Sink_re.num);
			Sink_re1.num = Sink_re.num;
			UnionRegion(attractorRegion, Sink_re1, Sink_re);
		}

		//Grow region from each saddle as repellers
		UnionRegion(repellerRegion, Source_re1, Source_re);
		for(i = 0; i < num_media; i++)
		{

			//after growing a region for one saddle, we need to union it with previous saddle regions
			repellerRegion.num = 0;
			cur_saddle = graphnodes[MediaNodes[i]].singularityID;

			//Initialize the saddle region 
			//InitSaddle_IndirectlyConnectedSingandCycle(cur_saddle,
			//	0, sing, cycle, percentage); //here 'percentage' is the global one!
			InitSaddle_IndirectlyConnectedSingandCycle(cur_saddle,
				0, sing, cycle, cur_length); //here 'cur_length' is the current percentage!
			
			//grow region from the initial region here
			UpdateBoundary(0);
			
			GetRegionNormals(0); ////Get the outward normal of the region for each boundary edge
			Cancel_Growing(0, -1);

			CopyRegion(Source_re.trianglelist, Source_re1.trianglelist, Source_re.num);
			Source_re1.num = Source_re.num;
			UnionRegion(repellerRegion, Source_re1, Source_re);
		}

		//get the intersection of the two region
		IntersectRegion(Source_re, Sink_re, intersectRegion);

		if(!IsRepeated(intersectRegion.trianglelist, singularities[sing].Triangle_ID,
			intersectRegion.num))
		{
			intersectRegion.trianglelist[intersectRegion.num] = singularities[sing].Triangle_ID;
			intersectRegion.num++;
		}
		
		//Judge whether the intersection region satisfies the Conley boundary condition or not
		if(CalEulerValue(intersectRegion.trianglelist, intersectRegion.num) == 0) //a ring-shaped
		{
			success = 1; //we do find a valid region to perform smoothing
			//success_length = cur_length;
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

		count++;

	}while(large_length - small_length > 1e-5 && (double)min_trianglesnum*(large_length - small_length) > 0.5
		&& (int)min_trianglesnum*cur_length >= 1/* && count < 1*/);

	if(success == 1)
		return true;

	return false;
}

/*
Initialize the length of the separatrices
*/
void 
InitSaddle_IndirectlyConnectedSingandCycle
(int saddle, int type, int sing, int cycle, double percentage)
{
	//Since both of the two cases will need to include the separatrices connecting 
	//to the attractor and repeller respectively
	//We can add them here!
	int i;
	int sep = singularities[saddle].separtices;
	int traj, other_traj;
	int pre_triangle = -1;

	if(limitcycles[cycle].type == 0) //the singularity will be an attractor
	{
		//follow the outgoing separatrix that leads to the singularity, add those triangles
		traj = separatrices[sep].sep1;
		other_traj = separatrices[sep].sep3;
		if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID != singularities[sing].Triangle_ID)
		{
			other_traj = traj;
			traj = separatrices[sep].sep3;	
			if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID != singularities[sing].Triangle_ID)
			{
				//some thing is wrong here
				MessageBox(NULL, "wrong connection!", "Error", MB_OK);
				return;
			}
		}
		
		//add those triangles
		for(i = 0; i < num_linesegs_curtraj[traj]; i++)
		{
			//add to attractor region
			if(trajectories[traj][i].Triangle_ID != pre_triangle)
			{
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, type);
				pre_triangle = trajectories[traj][i].Triangle_ID;
			}
		}

		//for the other sep, we probably need to add 30% of its triangles
		//pre_triangle = -1;
		//for(i = 0; i < (int)num_linesegs_curtraj[other_traj]*percentage/*0.3*/; i++)
		//{
		//	//remember: do not contain the singularity other than 'sing' and the intermediate saddles
		//	if(Object.flist[trajectories[other_traj][i].Triangle_ID]->contain_singularity == 1
		//		&& Object.flist[trajectories[other_traj][i].Triangle_ID]->singularity_index != sing)
		//		continue ;

		//	//add to attractor region
		//	if(trajectories[other_traj][i].Triangle_ID != pre_triangle)
		//	{
		//		AddToRegionTriangles(trajectories[other_traj][i].Triangle_ID, type);
		//		pre_triangle = trajectories[other_traj][i].Triangle_ID;
		//	}
		//}

		//follow the incoming separatrix that leads to the limit cycle, add those triangles
		traj = separatrices[sep].sep2;
		other_traj = separatrices[sep].sep4;

		if(!IsRepeated(limitcycles[cycle].cellcycle, trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID,
			limitcycles[cycle].num_triangles))
		{
			//if the last triangle of the sep falls in the cell cycle of the limit cycle, 
			//we say that the sep connects to the limit cycle
			//otherwise, choose the other one
			other_traj = traj;
			traj = separatrices[sep].sep4;
			if(!IsRepeated(limitcycles[cycle].cellcycle, trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID,
				limitcycles[cycle].num_triangles))
			{
				//some thing is wrong here
				MessageBox(NULL, "wrong connection!", "Error", MB_OK);
				return;
			}

		}
		
		pre_triangle = -1;
		for(i = 0; i < num_linesegs_curtraj[traj]; i++)
		{
			//add to repeller region

			if(trajectories[traj][i].Triangle_ID != pre_triangle)
			{
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, type);
				pre_triangle = trajectories[traj][i].Triangle_ID;
			}
		}

		//for the other sep, we probably need to add 30% of its triangles
		//pre_triangle = -1;
		//for(i = 0; i < (int)num_linesegs_curtraj[other_traj]*percentage/*0.3*/; i++)
		//{
		//	//add to attractor region
		//	//remember: do not contain the singularity other than 'sing' and the intermediate saddles
		//	if(Object.flist[trajectories[other_traj][i].Triangle_ID]->contain_singularity == 1
		//		&& Object.flist[trajectories[other_traj][i].Triangle_ID]->singularity_index != sing)
		//		continue ;

		//	if(trajectories[other_traj][i].Triangle_ID != pre_triangle)
		//	{
		//		AddToRegionTriangles(trajectories[other_traj][i].Triangle_ID, type);
		//		pre_triangle = trajectories[other_traj][i].Triangle_ID;
		//	}
		//}

	}

	else
	{
		//follow the outgoing separatrix that leads to the limit cycle, add those triangles
		traj = separatrices[sep].sep1;
		other_traj = separatrices[sep].sep3;

		if(!IsRepeated(limitcycles[cycle].cellcycle, trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID,
			limitcycles[cycle].num_triangles))
		{
			//if the last triangle of the sep falls in the cell cycle of the limit cycle, 
			//we say that the sep connects to the limit cycle
			//otherwise, choose the other one
			other_traj = traj;
			traj = separatrices[sep].sep3;
			if(!IsRepeated(limitcycles[cycle].cellcycle, trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID,
				limitcycles[cycle].num_triangles))
			{
				//some thing is wrong here
				MessageBox(NULL, "wrong connection!", "Error", MB_OK);
				return;
			}

		}
		
		pre_triangle = -1;
		for(i = 0; i < num_linesegs_curtraj[traj]; i++)
		{
			//add to repeller region

			if(trajectories[traj][i].Triangle_ID != pre_triangle)
			{
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, type);
				pre_triangle = trajectories[traj][i].Triangle_ID;
			}
		}
		
		//for the other sep, we probably need to add 30% of its triangles
		//pre_triangle = -1;
		//for(i = 0; i < (int)num_linesegs_curtraj[other_traj]*percentage/*0.3*/; i++)
		//{
		//	//remember: do not contain the singularity other than 'sing' and the intermediate saddles
		//	if(Object.flist[trajectories[other_traj][i].Triangle_ID]->contain_singularity == 1
		//		&& Object.flist[trajectories[other_traj][i].Triangle_ID]->singularity_index != sing)
		//		continue ;

		//	//add to attractor region
		//	if(trajectories[other_traj][i].Triangle_ID != pre_triangle)
		//	{
		//		AddToRegionTriangles(trajectories[other_traj][i].Triangle_ID, type);
		//		pre_triangle = trajectories[other_traj][i].Triangle_ID;
		//	}
		//}

		//follow the incoming separatrix that leads to the singularity, add those triangles
		traj = separatrices[sep].sep2;
		other_traj = separatrices[sep].sep4;

		if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID != singularities[sing].Triangle_ID)
		{
			other_traj = traj;
			traj = separatrices[sep].sep4;	
			if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID != singularities[sing].Triangle_ID)
			{
				//some thing is wrong here
				MessageBox(NULL, "wrong connection!", "Error", MB_OK);
				return;
			}
		}
		
		//add those triangles
		pre_triangle = -1;
		for(i = 0; i < num_linesegs_curtraj[traj]; i++)
		{
			//add to attractor region
			if(trajectories[traj][i].Triangle_ID != pre_triangle)
			{
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, type);
				pre_triangle = trajectories[traj][i].Triangle_ID;
			}
		}
		
		//for the other sep, we probably need to add 30% of its triangles
		//pre_triangle = -1;
		//for(i = 0; i < (int)num_linesegs_curtraj[other_traj]*percentage/*0.3*/; i++)
		//{
		//	//remember: do not contain the singularity other than 'sing' and the intermediate saddles
		//	if(Object.flist[trajectories[other_traj][i].Triangle_ID]->contain_singularity == 1
		//		&& Object.flist[trajectories[other_traj][i].Triangle_ID]->singularity_index != sing)
		//		continue ;

		//	//add to attractor region
		//	if(trajectories[other_traj][i].Triangle_ID != pre_triangle)
		//	{
		//		AddToRegionTriangles(trajectories[other_traj][i].Triangle_ID, type);
		//		pre_triangle = trajectories[other_traj][i].Triangle_ID;
		//	}
		//}
	}



	////The following codes add the triangles according to the adaptive length

	if(type == 0)  //saddle acts as repeller
	{
		pre_triangle = -1;

		//initialize the region including:
		//1) the triangles containing the sep(outgoing) connecting to the 'attractor'
		//2) the triangles containing the sep(incoming) connecting to the 'repeller'
		//3) the percentage * the triangles containing the other incoming sep
			
		traj = separatrices[sep].sep2;
		
		if(limitcycles[cycle].type == 0)
		{
			//optimize the incoming sep not connecting to the limit cycle
			
			if(IsRepeated(limitcycles[cycle].cellcycle, trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID,
				limitcycles[cycle].num_triangles))
			{
				//if the last triangle of the sep falls in the cell cycle of the limit cycle, 
				//we say that the sep connects to the limit cycle
				traj = separatrices[sep].sep4;
			}

		}

		else
		{
			//optimize the incoming sep not connecting to the singularity

			if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID == singularities[sing].Triangle_ID)
			{
				traj = separatrices[sep].sep4;	
			}
			
		}
		
		for(i = 0; i < (int)num_linesegs_curtraj[traj]*percentage; i++)
		{
			//remember: do not contain the singularity other than 'sing' and the intermediate saddles
			if(Object.flist[trajectories[traj][i].Triangle_ID]->contain_singularity == 1
				&& Object.flist[trajectories[traj][i].Triangle_ID]->singularity_index != sing)
				continue ;

			if(trajectories[traj][i].Triangle_ID != pre_triangle)
			{
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, type);
				pre_triangle = trajectories[traj][i].Triangle_ID;
			}
		}
	}

	else //saddle acts as attractor
	{
		//initialize the region including:
		//1) the triangles containing the sep(outgoing) connecting to the 'attractor'
		//2) the triangles containing the sep(incoming) connecting to the 'repeller'
		//3) the percentage * the triangles containing the other outgoing sep
		
		traj = separatrices[sep].sep1;

		if(limitcycles[cycle].type == 0)
		{
			//optimize the outgoing sep not connecting to the singularity

			if(trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID == singularities[sing].Triangle_ID)
			{
				traj = separatrices[sep].sep3;	
			}
		}

		else
		{
			if(IsRepeated(limitcycles[cycle].cellcycle, trajectories[traj][num_linesegs_curtraj[traj]-1].Triangle_ID,
				limitcycles[cycle].num_triangles))
			{
				//if the last triangle of the sep falls in the cell cycle of the limit cycle, 
				//we say that the sep connects to the limit cycle
				traj = separatrices[sep].sep3;
			}
		}
		
		for(i = 0; i < (int)num_linesegs_curtraj[traj]*percentage; i++)
		{
			//remember: do not contain the singularity other than 'sing' and the intermediate saddles
			if(Object.flist[trajectories[traj][i].Triangle_ID]->contain_singularity == 1
				&& Object.flist[trajectories[traj][i].Triangle_ID]->singularity_index != sing)
				continue ;

			if(trajectories[traj][i].Triangle_ID != pre_triangle)
			{
				AddToRegionTriangles(trajectories[traj][i].Triangle_ID, type);
				pre_triangle = trajectories[traj][i].Triangle_ID;
			}
		}
	}
}

/*
*/
bool IsDirectlyConnected(int node1, int node2)
{
	int i;
	int sec_node;

	for(i = 0; i < graphnodes[node1].nedges; i++)
	{
		sec_node = graphedges[graphnodes[node1].edges[i]].node_index1;
		if(sec_node == node1)
			sec_node = graphedges[graphnodes[node1].edges[i]].node_index2;

		if(sec_node == node2)
			return true;
	}

	return false;
}

/*
Deal with limit cycle and center singularity pair cancellation 07/05/06
*/
bool 
DirectlyConnectedSingandCyclePair(int sing, int cycle)
{
	SetFence_LimitCyclesexcept(cycle);
	SetFenceForSeps(-1);
	RemoveLimitCycleFence(cycle);

	////Grow the region for the limit cycle
	GrowRegionforALimitCycle(cycle);

	//
	if(limitcycles[cycle].type == 0) //the singularity acts as an attractor
	{
		Cancel_GrowAttractorRegion(singularities[sing].Triangle_ID, -1, sing, InitSaddleRegionLength, 0, 1.);
	}

	else
	{
		Cancel_GrowRepellerRegion(singularities[sing].Triangle_ID, -1, sing, InitSaddleRegionLength, 0, 1.);
	}

	//Get the intersection of the region
	IntersectRegion(repellerRegion, attractorRegion, intersectRegion);

	//Add the center triangle
	if(!IsRepeated(intersectRegion.trianglelist, singularities[sing].Triangle_ID, intersectRegion.num))
	{
		intersectRegion.trianglelist[intersectRegion.num] = singularities[sing].Triangle_ID;
		intersectRegion.num++;
	}

	if(CalEulerValue(intersectRegion.trianglelist, intersectRegion.num) != 1)
		return false;

	//Get the inner vertices of the intersect region
	GetInnerVerts(sing, NULL, 0);

	UpdateBoundary(2);
	SetBoundaryFlag_Ver(intersectBoundary.edgelist, intersectBoundary.num);
	Cancel_RegionSmooth(); //perform cancellation

	int totalindex = 0;
	if(IntersectedRegionSingCount(totalindex) == 1)
		return true;

	else
	{
		Undo();
		return false;
	}

}



/* Main subroutine for limit cycle and singularity pair cancellation */
void LimitCycleSingularityCancel(int sing, int cycle)
{
	int node1, node2;
	//if the singularity and the limit cycle are directly connected in C-graph (use a linear searching)
	if(singularities[sing].type == SADDLE)
	{
		//Judge whether they are doubly connected
	}

	else
	{
		//Judge whether they are doubly connected!!!!! (new finding 06/27/06)

		node1 = singularities[sing].node_index;
		node2 = limitcycles[cycle].node_index;
		if(IsDirectlyConnected(node1, node2))
		{
			// call the directly connected limit cycle and center singularity cancellation routine
		}

		else
		{
			IndirectlyConnectedSingandCyclePair(sing, cycle);
		}
	}
}


