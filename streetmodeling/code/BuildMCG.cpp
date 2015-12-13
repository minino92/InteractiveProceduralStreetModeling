/*
BuildMCG.cpp

This file implements the construction of the MCG (Morse set level Conley/Connection Graph).

The basic ideas are similar to the algorithm we used to build ECG.

Assume that we assign a flag for each boundary of each valid SCC, respectively.
(Here, a valid SCC is either a single triangle containing fixed point;
or a ring-liked triangle strip, or a disk triangle set containing fixed points or periodic orbits)

Definition of repelling elements: the flow on all boundaries of the elements is leaving the element
Definition of attracting elements: the flow of all boundaries of the elements is coming in to the element
Definition of saddle-like elements(need more thought): 
    the flow of the boundaries of the elements is either coming or leaving the element
	(Note that, saddle region may contain more than one triangle)
Definition of other elements (could they exist ?)

1. First, we get the connections bying tracing the traditional separatrices

2. Second, we start from the unmarked sides of repelling elements(but not including the single triangle element)
   We trace along the flow direction

3. Third, we start from the unmarked sides of attracting elements(but not including the single triangle element)
   we trace along the inverse of the flow

4. If there are still some unmarked SCC (usually they are single triangle elements), choose a point in the
residual of each of them to do tracing to find out remaining connections

For all of these tracing
First one, we find the starting points according to the saddle positions and their eigen vectors

2~4, we need to choose a point on the boundary of each element (we can simply choose the middle point
of a boundary edge)

The criterion of stop condition for tracing (difficult part?)

The tracing will stop when the tracing curve hits the boundary of one of the other elements.
This means that we need to perform triangle testing each time when we enter a new triangle, slow?
(It seems that we need another flag for each triangle to record the index of SCC it belongs to)


*************New way to find out the connections between two Morse sets***********************
02/21/07

We have reuse the directed graph we obtained before.
Start from the boundary nodes of each SCC, search along its edges until it hits the other SCC.
Do this for every boundary of each SCC
(Problem: how many nodes we should use for each single boundary?
can we prove that searching from boundary is sufficient? -- it seems so. If boundary searching and 
interior searching return different results, it will tell us that we have 
a intersect during the searching, which violates the rule of streamline in a vector field.
Otherwise, we may get a incorrect SCC, which we can guarantee will not happen using the Finding
SCC algorithm (reference).
)

To do that, we need to first judge the type of the Morse set.
Only three types can be exist in practice: 
One is repelling, which means all its boundary edges are exit or mixed edges;
Second is attracting, which means all its boundary edges are entrance edges or mixed edges 
All the other SCC's not satisfied previous two definitions belong to the third type
--interval (can we prove that in this type of SCC's, there are at least one saddle inside them?).

created by Guoning Chen 02/05/07
*/

#include "stdafx.h"

#include "FindSCC.h"

#include "LocalTracing.h"

#include "lib/icVector.h"

#include "VFDataStructure.h"


/*we may also use routines from Taustep.cpp */

extern SCCList scclist;
extern Polygon3D Object;
extern GraphEdge *sccedges;
extern int num_sccedges ;
extern int curMaxNumDirGraphEdges;
extern GraphNode2 *sccnodes;
extern int *cur_nodes_order;
extern int num_sccnodes;
extern int num_sccomps;               //record the current number of SCC components
extern Singularities *singularities;

extern int r_counter, a_counter, s_counter; //for labeling 3/14/06

void connect_From_Saddle_Elem(int SCC_id);

void connect_From_Repelling_Elem(int SCC_id, int which_boundary);



/*
First try to build the MCG
*/
#include "topologyedit.h"
extern TriangularRegion repellerRegion;       ////region containing a repeller
extern TriangularRegion attractorRegion;      ////region containing an attractor
extern RegionBoundary repellerBoundary;
extern RegionBoundary attractorBoundary;

extern bool RepellerExitEdgePending(Edge *cur_edge);
extern bool AttractorExitEdgePending(Edge *cur_edge);
extern int GetOppositeTriangle(Edge *, int);

void grow_from_region(int *triangles, int n, int scc, int type);
void grow_for_connect(int type, int scc_index);
bool is_saddle_region(int scc_index);
bool is_repelling_region(int scc_index);
bool is_attracting_region(int scc_index);
void grow_from_saddle_region();
void assign_mcgnodes();
void add_to_mcg_nodelist(int node);
void add_to_mcg_edgelist(int node1, int node2);
void add_edge_to_node(int node, int edgeindex);
void layout_mcg();
void init_growing();

void build_mcg_edges_1();

void build_mcg();

int get_region_type(int scc_index);
int get_region_type_2(int scc_index);
int make_correction_saddle(int scc_index);

void set_region_flags(int *triangle, int num, int type);
void remove_region_flags(int *triangle, int num);

void reset_sccnodes_flags();
void set_sccnodes_flags_for_SCC(int);

/*****************************************************************************/
/*the following routines implement the same idea but using
the obtained graph after \tau tracing instead
*/

/*
Here, we judge the type of each region using the obtained graph after \tau tracing
*/
int get_region_type_2(int scc_index);
int get_triangle_type(int tri, int scc_index);

void remove_redundant_edges();

void grow_repeller_region_graph_2(int, int);
void grow_all_mcgnodes();


void cal_mcgedge_regions();


/*---------------------------------------------------------------------------*/

extern GraphNode *graphnodes;
extern GraphEdge *graphedges;
extern int num_sccnodes;
//extern int cur_node_index;
extern int cur_graphedge_index;
//int curMaxNumMCGNodes;

/*dual graph*/
extern GraphEdge *edges_dual;
extern int num_edges_dual;
//int curMaxNumDirGraphEdges;
extern GraphNode2 *nodes_dual;
extern int num_nodes_dual;


MCGNode *mcgnodes = NULL;
MCGEdge *mcgedges = NULL;
int cur_mcgnode_index;
int curMaxNumMCGNodes;
int cur_mcgedge_index;
int curMaxNumMCGEdges;

//
extern int *searchstack;
extern int nelem_in_stack;
extern int MaxStackElems;

void reset_sccnodes_flags()
{
	int i;

	for(i = 0; i < num_sccnodes; i++)
	{
		sccnodes[i].visited = 0;
	}
}


void add_to_mcg_nodelist(int scc_index)
{
	if(cur_mcgnode_index >= curMaxNumMCGNodes)
	{
		/*extend the memory*/
		//graphnodes = (GraphNode *)realloc(graphnodes, sizeof(GraphNode)*(curMaxNumMCGNodes+50));
		mcgnodes = (MCGNode *)realloc(mcgnodes, sizeof(MCGNode)*(curMaxNumMCGNodes+50));
		if(mcgnodes == NULL)
		{
			MessageBox(NULL, "failed to allocate memory for mcgnodes!", "Error", MB_OK);
			/*later, we should right to a log file*/
			exit(-1);
		}

		curMaxNumMCGNodes+=50;
	}

	mcgnodes[cur_mcgnode_index].index = cur_mcgnode_index;
	mcgnodes[cur_mcgnode_index].scc_index = scc_index;
	mcgnodes[cur_mcgnode_index].edges = NULL;
	mcgnodes[cur_mcgnode_index].nedges = 0;
	mcgnodes[cur_mcgnode_index].visited = 0;
	mcgnodes[cur_mcgnode_index].cancelled = 0;

	/*we also need to get the type of the node!!!*/
	//mcgnodes[cur_mcgnode_index].type = get_region_type(scc_index);
	mcgnodes[cur_mcgnode_index].type = get_region_type_2(scc_index);

	/*we make the correction for the saddle region with triangles fewer than 5
	03/08/07*/
	if(mcgnodes[cur_mcgnode_index].type == 2
		&& scclist.scccomponents[scc_index].num_nodes < 5)
		mcgnodes[cur_mcgnode_index].type = make_correction_saddle(scc_index);

	/*------------------------------------------------------------------*/
	/*the inaccurate judgement could happen for larger Morse sets!!!!*/

	/*------------------------------------------------------------------*/

	if( mcgnodes[cur_mcgnode_index].type == 0)
	{
		mcgnodes[cur_mcgnode_index].labelindex = r_counter;
		r_counter++ ; 
	}
	else if(mcgnodes[cur_mcgnode_index].type == 1)
	{
		mcgnodes[cur_mcgnode_index].labelindex = a_counter;
		a_counter++;
	}
	else
	{
		mcgnodes[cur_mcgnode_index].labelindex = s_counter;
		s_counter++;
		//FILE *fp = fopen("found_saddle_region.txt", "a");
		//fprintf(fp, "current number of saddles: %d\n", s_counter);
		//fclose(fp);
	}

	/*
	File test
	*/
	//FILE *fp = fopen("error_build_mcg.txt", "w");
	//fprintf(fp, "the current number of nodes in mcg is: %d\n", cur_mcgnode_index);
	//fprintf(fp, "repellers: %d, attractors: %d, saddles: %d\n", r_counter, a_counter, s_counter);
	//fclose(fp);

	cur_mcgnode_index++;

}

void init_edgelist()
{
	curMaxNumMCGEdges = 200;
	mcgedges = (MCGEdge *)malloc(sizeof(MCGEdge)*curMaxNumMCGEdges);

	if(mcgedges == NULL)
	{
		MessageBox(NULL, "fail to initialize mcgedges!", "Error", MB_OK);
		exit(-1);
	}
}


/*
Add a new edge into the edge list
*/
void add_to_mcg_edgelist(int node1, int node2)
{
	if(cur_mcgedge_index >= curMaxNumMCGEdges)
	{
		mcgedges = (MCGEdge *)realloc(mcgedges, sizeof(MCGEdge)*(curMaxNumMCGEdges+50));

		if(mcgedges == NULL)
		{
			MessageBox(NULL, "fail to reallocate mcgedges!", "Error", MB_OK);
			exit(-1);
		}
		curMaxNumMCGEdges+=50;
	}

	mcgedges[cur_mcgedge_index].edge_index = cur_mcgedge_index;
	mcgedges[cur_mcgedge_index].node_index1 = node1;
	mcgedges[cur_mcgedge_index].node_index2 = node2;
	mcgedges[cur_mcgedge_index].cancel = false;
	mcgedges[cur_mcgedge_index].visited = 0;

	cur_mcgedge_index++;
}


void add_edge_to_node(int node, int edgeindex)
{
	if(node < 0)
		return;

	int i;
	int *temp_list = mcgnodes[node].edges;

	if(temp_list == NULL)
	{
		mcgnodes[node].edges = (int*)malloc(sizeof(int));
		mcgnodes[node].edges[0] = edgeindex;
		mcgnodes[node].nedges = 1;
	}

	else{
		mcgnodes[node].edges = (int*)malloc(sizeof(int)*(mcgnodes[node].nedges+1));
		if(mcgnodes[node].edges == NULL)
		{
			MessageBox(NULL, "fail to allocate edges for mcgnodes", "Error", MB_OK);
			exit(-1);
		}
		for(i = 0; i < mcgnodes[node].nedges; i++)
		{
			mcgnodes[node].edges[i] = temp_list[i];
		}
		mcgnodes[node].edges[i] = edgeindex;

		mcgnodes[node].nedges += 1;

		free(temp_list);
	}
}


/*routine to judge whether there is at least one interval(saddle)
between two nodes (usually a repeller and an attractor)
we always suppose the inputs are valid 03/08/07
*/
bool has_interval_mcg(int node1, int node2)
{
	int i, j;

	////we need to make sure node1 is the repeller
	if(mcgnodes[node1].type == 1)
	{
		int temp = node1;
		node1 = node2;
		node2 = temp;
	}


	////Begin to judge whether there is at least one interval connecting them
	int *temp_connect_saddlelist = new int[mcgnodes[node1].nedges];
	int num_connect_saddlelist = 0;
	int cur_edge, n2;

	////Begin searching from this repeller to search all the connected saddles
	for(j = 0; j < mcgnodes[node1].nedges; j++)
	{
		cur_edge = mcgnodes[node1].edges[j];
		n2 = mcgedges[cur_edge].node_index2;

		if(mcgnodes[n2].type == 2) //if this is a saddle
		{
			////we need to perform repeating test if we allow multiple connections between two nodes!!
			temp_connect_saddlelist[num_connect_saddlelist] = n2;
			num_connect_saddlelist++;
		}
	}

	for(j = 0; j < num_connect_saddlelist; j++)
	{
		for(i = 0; i < mcgnodes[temp_connect_saddlelist[j]].nedges; i++)
		{
			cur_edge = mcgnodes[temp_connect_saddlelist[j]].edges[i];
			n2 = mcgedges[cur_edge].node_index2;

			if(n2 == node2) //we find an interval
			{
				free(temp_connect_saddlelist);
				return true;
			}
		}
	}

	free(temp_connect_saddlelist);
	return false;
}


/*Judge whether there has already one edge between them
Note: not consider double connection now 03/08/07
*/
bool is_connected_mcg(int node1, int node2)
{
	int i, cur_edge;

	for(i = 0; i < mcgnodes[node1].nedges; i++)
	{
		cur_edge = mcgnodes[node1].edges[i];
		
		if(mcgedges[cur_edge].node_index1 == node1)
		{
			if(mcgedges[cur_edge].node_index2 == node2)
				return true;
		}

		else{
			if(mcgedges[cur_edge].node_index1 == node2)
				return true;
		}
	}

	return false;
}


/*lay out the nodes*/
extern double ConleyGraphWin_x ;
extern double ConleyGraphWin_y ;
extern double types_interval_y ;

void layout_mcg()
{
	int i;
	int num_repellers, num_attractors, num_saddles;
	int cur_repeller_index, cur_attractor_index, cur_saddle_index;
	double repeller_interval, attractor_interval, saddle_interval;
    double repeller_y, attractor_y, saddle_y;


	num_repellers = num_attractors = num_saddles = 0;
	cur_repeller_index = cur_attractor_index = cur_saddle_index = 0;

	////Actually, they may be counted during the building of the graph
	////1. First, count the number of the nodes of 3 types respectively
	for(i = 0; i < cur_mcgnode_index; i++)
	{
		if(mcgnodes[i].type == 0)
			num_repellers ++;
		else if(mcgnodes[i].type == 1)
			num_attractors ++;
		else
			num_saddles ++;
	}

	repeller_interval = ConleyGraphWin_x / (num_repellers+1);
	attractor_interval = ConleyGraphWin_x / (num_attractors+1);
	saddle_interval = ConleyGraphWin_x / (num_saddles+1);

	attractor_y = types_interval_y;
	saddle_y = 2 * types_interval_y;
	repeller_y = 3 * types_interval_y;
	
	FILE *fp = fopen("error_build_mcg.txt", "a");
	fprintf(fp, "repellers: %d;  attractors: %d;  saddles: %d\n", num_repellers, num_attractors, num_saddles);
	fprintf(fp, "start calculating the positions...\n", num_repellers, num_attractors, num_saddles);
	fclose(fp);

	////2. Then, evenly space the positions of these nodes
	for(i = 0; i < cur_mcgnode_index; i++)
	{
		if(mcgnodes[i].type == 0)
		{
			mcgnodes[i].pos_x = (cur_repeller_index+1)*repeller_interval;
			mcgnodes[i].pos_y = repeller_y;
			cur_repeller_index ++;
		}
		else if(mcgnodes[i].type == 1)
		{
			mcgnodes[i].pos_x = (cur_attractor_index+1)*attractor_interval;
			mcgnodes[i].pos_y = attractor_y;
			cur_attractor_index ++;
		}
		else{
			mcgnodes[i].pos_x = (cur_saddle_index+1)*saddle_interval;
			mcgnodes[i].pos_y = saddle_y;
			cur_saddle_index ++;
		}
	}

	/*file test*/
	//FILE *fp = fopen("error_build_mcg.txt", "a");
	//fprintf(fp, "repellers: %d;  attractors: %d;  saddles: %d\n", num_repellers, num_attractors, num_saddles);
	//fclose(fp);
	fp = fopen("error_build_mcg.txt", "a");
	fprintf(fp, "finish calculating the positions.\n", num_repellers, num_attractors, num_saddles);
	fclose(fp);
}


/*
Assign nodes for the Morse sets in MCG
*/
void assign_mcgnodes()
{
	////Release previous space
	//if(mcgnodes != NULL)
	//	free(graphnodes);

	////Allocate new space for nodes 
	if(mcgnodes == NULL)
	{
		curMaxNumMCGNodes = 200;
		mcgnodes = (MCGNode *)malloc(sizeof(MCGNode) * curMaxNumMCGNodes);

		if(mcgnodes == NULL)
		{
			MessageBox(NULL, "fail to allocate mcgnodes!", "Error", MB_OK);
			exit(-1);
		}
	}
	cur_mcgnode_index = 0;

	r_counter = a_counter = s_counter = 0;

	int i;

	for(i = 0; i < scclist.num_sccs; i++)
	{
		if(scclist.scccomponents[i].num_nodes < 2 && scclist.scccomponents[i].num_singularities < 1)
			continue;

		if(scclist.scccomponents[i].num_nodes >= 2 && scclist.scccomponents[i].valid == 0)
			continue;

		add_to_mcg_nodelist(i);
		scclist.scccomponents[i].node_index = cur_mcgnode_index-1;
	}

}

bool is_saddle_region(int scc_index)
{
	int i;

	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		if(Object.flist[scclist.scccomponents[scc_index].nodes[i]]->contain_singularity == 1)
		{
			if(singularities[Object.flist[scclist.scccomponents
				[scc_index].nodes[i]]->singularity_index].type == SADDLE)
				return true;
		}
	}
	return false;
}

void grow_from_saddle_region()
{
	int i;

	for(i = 0; i < scclist.num_sccs; i++)
	{
		if(scclist.scccomponents[i].num_nodes < 2 && scclist.scccomponents[i].num_singularities < 1)
			continue;

		if(scclist.scccomponents[i].num_nodes > 2 && scclist.scccomponents[i].valid == 0)
			continue;

		if(is_saddle_region(i))
		{
			/*first, find the connections to attracting regions */
			grow_from_region(scclist.scccomponents[i].nodes,
				scclist.scccomponents[i].num_nodes, i, 0);

			/*second, find the connections to repelling regions*/
			grow_from_region(scclist.scccomponents[i].nodes,
				scclist.scccomponents[i].num_nodes, i, 1);
		}
	}
}


void grow_from_region(int *triangles, int n, int scc_index, int type)
{
	/*According to the type, grow the region forward (type = 0) or backward (type = 1)*/

	if(type == 0) /*it is forward growing */
	{
		/*copy triangles to the repellerRegion*/
		CopyRegion(triangles, repellerRegion.trianglelist, n);
		repellerRegion.num = n;
	}

	else
	{
		/*copy triangles to the attractorRegion*/
		CopyRegion(triangles, attractorRegion.trianglelist, n);
		attractorRegion.num = n;
	}

	grow_for_connect(type, scc_index);
}


void grow_for_connect(int type, int scc_index)
{
	int i;
	bool exitornot = false;
	int num_edges;
	Edge **edgelist;
	Edge *cur_edge;
    int oppositeTriangle;
	int num_newaddedtriangle;
	int stop_flag = 0;
	int oneother_scc = scc_index;

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
				//if the triangle containing this fixed point has the same SCC index of the region
				if(sccnodes[oppositeTriangle].sscomp_index == scc_index)
				{
					////add the triangle to the region
					AddToRegionTriangles(oppositeTriangle, type);
					num_newaddedtriangle++;
				}

				else /*this means that we reach another Morse set*/
				{
					if(stop_flag >= 1)
						break; /*but, this can guarantee to find only one other Morse set */
					else
					{
						/*avoid reaching the same SCC twice*/
						oneother_scc = sccnodes[oppositeTriangle].sscomp_index;
						stop_flag ++;

						/*probably add a new edge to MCG*/
						if(type == 0)
							add_to_mcg_edgelist(scclist.scccomponents[scc_index].node_index,
							scclist.scccomponents[oneother_scc].node_index);
						else
							add_to_mcg_edgelist(scclist.scccomponents[oneother_scc].node_index,
							scclist.scccomponents[scc_index].node_index);
					}
				}
			}

			//if(Object.flist[oppositeTriangle]->fence_flag == 1) ////set fence 2/16/06
			//	continue;

			/*if it reach other Morse decomposed region*/
			if(sccnodes[oppositeTriangle].sscomp_index != scc_index)
			{
				if((scclist.scccomponents[sccnodes[oppositeTriangle].sscomp_index].num_singularities > 0
					||scclist.scccomponents[sccnodes[oppositeTriangle].sscomp_index].valid == 1)
					&& oneother_scc != sccnodes[oppositeTriangle].sscomp_index)
				{
					if(stop_flag >= 1)
						break; /*but, this can guarantee to find only one other Morse set */
					else{
						stop_flag ++;
						oneother_scc = sccnodes[oppositeTriangle].sscomp_index;
					}
				}
			}

			if(exitornot && oppositeTriangle >= 0)  ////it is an exiting edge 
			{
				////if the opposite triangle is not inside the region!!!
				////add the adjacent triangle into the region
				AddToRegionTriangles(oppositeTriangle, type);
				num_newaddedtriangle++;

			}

		}

		if(num_newaddedtriangle == 0) ////No more exiting/entrance edges can be found
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


void init_growing()
{
	int i;
	Face *face;

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];
		face->attract_inregion = 0;
		face->repell_inregion = 0;
	}

}
/*get the type of the specified region according to the flow behavior
of the boundaries of the region
*/
extern void build_Boundary_Edgelist(int *triangle, int ntris);
extern bool AttractorExitEdgePending(Edge *cur_edge);
extern bool RepellerExitEdgePending(Edge *cur_edge);
void CalNormalAtEdge(Edge *cur_edge, Face *face, int type);


bool is_exit_edge(Edge *cur_edge)
{
	double dot1, dot2, dotresult = 0;
	Vertex *v1, *v2;
	v1 = Object.vlist[cur_edge->verts[0]];
	v2 = Object.vlist[cur_edge->verts[1]];

	dot1 = dot(cur_edge->repell_normal, v1->vec);
	dot2 = dot(cur_edge->repell_normal, v2->vec);

	//dotresult = max(dot1, dot2);   //growing with mixed edges
	dotresult = min(dot1, dot2); //pure exit edge growing

	if(dotresult >= 0) return true;
	//if(dotresult > 0) return true;

	return false;
}

bool is_entrance_edge(Edge *cur_edge)
{
	double dot1, dot2, dotresult = 0;
	Vertex *v1, *v2;
	v1 = Object.vlist[cur_edge->verts[0]];
	v2 = Object.vlist[cur_edge->verts[1]];

	dot1 = dot(cur_edge->repell_normal, v1->vec);
	dot2 = dot(cur_edge->repell_normal, v2->vec);
	//dot1 = dot(cur_edge->attract_normal, v1->vec);
	//dot2 = dot(cur_edge->attract_normal, v2->vec);

	//dotresult = min(dot1, dot2);   //growing with mixed edges
	dotresult = max(dot1, dot2); //pure entrance edge growing

	if(dotresult <= 0) return true;
	//if(dotresult < 0) return true;

	return false;
}


extern Edge **Cycle_edge;
extern int num_cycleedges;

void get_boundary_normals(Edge **edges, int nedges, int type)
{
	int i;
	Edge *cur_edge;
	Face *face;
	int EdgeNum;
	Edge **edgelist;

	EdgeNum = nedges;
	edgelist = edges;

	for(i = 0; i < EdgeNum; i ++)
	{
		cur_edge = edgelist[i];

		////Get the face that inside the region
		if(cur_edge->tris[0] !=0)
			face = Object.flist[cur_edge->tris[0]];
		else
			face = Object.flist[cur_edge->tris[1]];

		if(type == 0){
			if(face->repell_inregion == 0 && cur_edge->tris[1] >=0)
				face = Object.flist[cur_edge->tris[1]];
		}
		else{
			if(face->attract_inregion == 0&& cur_edge->tris[1] >=0)
				face = Object.flist[cur_edge->tris[1]];
		}

		CalNormalAtEdge(cur_edge, face, type);
	}
}

int get_region_type(int scc_index)
{
	//if(is_saddle_region(scc_index)) /*I am cheating here!*/
	//{
	//	FILE *fp = fopen("found_saddle_region.txt", "a");
	//	fprintf(fp, "find a saddle region: %d\n", scc_index);
	//	fclose(fp);
	//	return 2;
	//}

	/*we need to evaluate the boundary edges of the region one by one*/
	set_region_flags(scclist.scccomponents[scc_index].nodes,
		scclist.scccomponents[scc_index].num_nodes, 0);

	/*first find out all the boundary edges*/
	build_Boundary_Edgelist(scclist.scccomponents[scc_index].nodes, 
		scclist.scccomponents[scc_index].num_nodes);

	/*calculate the outward normal*/
	get_boundary_normals(Cycle_edge, num_cycleedges, 0);
	//get_boundary_normals(Cycle_edge, num_cycleedges, 1);
	
	int i, r_flag, a_flag;
	r_flag = 0;
	a_flag = 0;
	Edge *e;
	/*judge the exit*/
	for(i = 0; i < num_cycleedges; i++)
	{
		e = Cycle_edge[i];
		if(is_exit_edge(e)) /*it is exit edge*/
		{
			r_flag = 1;

			if(a_flag == 1)
			{
				remove_region_flags(scclist.scccomponents[scc_index].nodes,
					scclist.scccomponents[scc_index].num_nodes);
				return 2;       //this is a saddle region
			}
		}
		else if(is_entrance_edge(e)) /*it is entrance edge*/
		{
			a_flag = 1;

			if(r_flag == 1)
			{
				remove_region_flags(scclist.scccomponents[scc_index].nodes,
					scclist.scccomponents[scc_index].num_nodes);
				return 2;       //this is a saddle region
			}
		}

		else  /*it is mixed edge*/
		{
			//a_flag = 1;

			//if(r_flag == 1)
			//{
				remove_region_flags(scclist.scccomponents[scc_index].nodes,
					scclist.scccomponents[scc_index].num_nodes);
				return 2;       //this is a saddle region
			//}
		}
	}
	
	/*judge the entrance*/
	//for(i = 0; i < num_cycleedges; i++)
	//{
	//	e = Cycle_edge[i];
	//	if(is_entrance_edge(e))
	//	{
	//		a_flag = 1;

	//		if(r_flag == 1)
	//		{
	//			remove_region_flags(scclist.scccomponents[scc_index].nodes,
	//				scclist.scccomponents[scc_index].num_nodes);
	//			return 2;       //this is a saddle region
	//		}
	//	}
	//	else
	//	{
	//		r_flag = 1;

	//		if(a_flag == 1)
	//		{
	//			remove_region_flags(scclist.scccomponents[scc_index].nodes,
	//				scclist.scccomponents[scc_index].num_nodes);
	//			return 2;       //this is a saddle region
	//		}
	//	}
	//}

	if(r_flag == 1 && a_flag == 0)
	{
		remove_region_flags(scclist.scccomponents[scc_index].nodes,
			scclist.scccomponents[scc_index].num_nodes);
		return 0;               //this is a repelling region
	}

	if(r_flag == 0 && a_flag == 1)
	{
		remove_region_flags(scclist.scccomponents[scc_index].nodes,
			scclist.scccomponents[scc_index].num_nodes);
		return 1;               //this is an attracting region
	}

}



void set_region_flags(int *triangle, int num, int type)
{
	int i;

	for(i = 0; i < num; i++)
	{
		if(type == 0)
			Object.flist[triangle[i]]->repell_inregion = 1;
		else
			Object.flist[triangle[i]]->attract_inregion = 1;
	}
}


void remove_region_flags(int *triangle, int num)
{
	int i;

	for(i = 0; i < num; i++)
	{
		Object.flist[triangle[i]]->repell_inregion = 1;
		Object.flist[triangle[i]]->attract_inregion = 1;
	}
}


void build_mcg_edges_1()
{
	int i;
	int scc_index;

	if(mcgedges == NULL)
	{
		curMaxNumMCGEdges = 50;
		mcgedges = (MCGEdge *)malloc(sizeof(MCGEdge)*curMaxNumMCGEdges);
		
		if(mcgedges == NULL)
		{
			MessageBox(NULL, "fail to allocate mcgedges!", "Error", MB_OK);
			exit(-1);
		}
	}

	cur_mcgedge_index=0;

	/*first, grow from saddle region*/

	for(i = 0; i < cur_mcgnode_index; i++)
	{
		if(mcgnodes[i].type == 2)
		{
			scc_index = mcgnodes[i].scc_index;
			init_growing();
			/*first, find the connections to attracting regions */
			set_region_flags(scclist.scccomponents[scc_index].nodes,
				scclist.scccomponents[scc_index].num_nodes, 0);

			grow_from_region(scclist.scccomponents[scc_index].nodes,
				scclist.scccomponents[scc_index].num_nodes, scc_index, 0);

			/*second, find the connections to repelling regions*/
			init_growing();

			set_region_flags(scclist.scccomponents[scc_index].nodes,
				scclist.scccomponents[scc_index].num_nodes, 1);

			grow_from_region(scclist.scccomponents[scc_index].nodes,
				scclist.scccomponents[scc_index].num_nodes, scc_index, 1);
		}
	}

	/* second, we grow from the repelling regions, but which ones???????*/
}


/*----------03/07/07--------------*/
/* The second group of routines constructing MCG using obtained directed graph
after \tau tracing*/

int *boundary_nodes = NULL;
int nboundarytris = 0;
/*
Get the boundary nodes of the specified SCC
*/
int get_boundary_nodes(int scc_index)
{
	int MaxTris;
	if(scclist.scccomponents[scc_index].num_nodes < 5)
	{
		if(boundary_nodes != NULL)
			free(boundary_nodes);
		MaxTris = (int)scclist.scccomponents[scc_index].num_nodes;
		boundary_nodes = (int*)malloc(sizeof(int)*MaxTris);
		nboundarytris = 0;
		int i;
		for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
		{
			boundary_nodes[nboundarytris]=scclist.scccomponents[scc_index].nodes[i];
			nboundarytris++;
		}
		return 1;
	}


	MaxTris = (int)scclist.scccomponents[scc_index].num_nodes;
	boundary_nodes = (int*)malloc(sizeof(int)*MaxTris);
	nboundarytris = 0;

	int i, j;
	Face *face;
	Edge *e;

	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		face = Object.flist[scclist.scccomponents[scc_index].nodes[i]];
		for(j = 0; j < 3; j++)
		{
			e = face->edges[j];

			if(IsRepeated(scclist.scccomponents[scc_index].nodes, e->tris[0],
				scclist.scccomponents[scc_index].num_nodes)
			 && IsRepeated(scclist.scccomponents[scc_index].nodes, e->tris[1],
				scclist.scccomponents[scc_index].num_nodes))
				continue;  //this is an inner edge

			/*this is a boundary triangle*/
			boundary_nodes[nboundarytris] = face->index;
			nboundarytris++;
			break;
		}
	}
	return nboundarytris;
}

/*
Get the type of the region according to the incoming and outgoing edges of the nodes in 
the SCC component
*/
int get_region_type_2(int scc_index)
{
	int i;
	int r_flag, a_flag;

	r_flag = a_flag = 0;


	if(get_boundary_nodes(scc_index)==1)
	{
		for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
		{
			if(get_triangle_type(scclist.scccomponents[scc_index].nodes[i], scc_index) == 0)
			{
				r_flag = 1;

				if(a_flag == 1)
					return 2;
			}

			else if(get_triangle_type(scclist.scccomponents[scc_index].nodes[i], scc_index) == 1)
			{
				a_flag = 1;

				if(r_flag == 1)
					return 2;
			}
			else
				return 2;
		}
	}

	else{
		for(i = 0; i < nboundarytris; i++)
		{
			if(get_triangle_type(boundary_nodes[i], scc_index) == 0)
			{
				r_flag = 1;

				if(a_flag == 1)
					return 2;
			}

			else if(get_triangle_type(boundary_nodes[i], scc_index) == 1)
			{
				a_flag = 1;

				if(r_flag == 1)
					return 2;
			}
			else
				return 2;
		}
	}

	if(r_flag == 1 && a_flag == 0)
		return 0;
	if(r_flag == 0 && a_flag == 1)
		return 1;
}

/*the routine return the type of the specified triangle "tri"
0-- repelling, means that all its incoming edges are inside the 
same SCC
1-- attracting,
2-- saddle

The problem is the triangle having tangential edge may be treated
as saddle triangle
*/
int get_triangle_type(int tri, int scc_index)
{
	//GraphNode2 *cur_node = &sccnodes[tri];
	int other_node;
	int i;
	int r_flag, a_flag;
	r_flag = a_flag = 0;

	for(i = 0; i < sccnodes[tri].nedges; i++)
	{
		/*we do not consider the edge pointing to itself*/
		if(sccedges[sccnodes[tri].edges[i]].node_index1 == tri
			&& sccedges[sccnodes[tri].edges[i]].node_index2 == tri)
			continue;

		other_node = sccedges[sccnodes[tri].edges[i]].node_index2;

		/*this is an outward edge*/
		if(other_node != tri)
		{
			if(sccnodes[other_node].sscomp_index != scc_index)
			{
				//r_flag = 1;  //it has at least one outgoing edge pointing to *other* SCC
				//if(a_flag == 1)
				//	return 2;
				r_flag ++;
			}
			continue;
		}

		/*this is an incoming edge*/
		other_node = sccedges[sccnodes[tri].edges[i]].node_index1;

		/*this edge is coming from the other node in the same SCC*/
		if(sccnodes[other_node].sscomp_index == scc_index)
			continue;
		else{
			//a_flag = 1;
			//if(r_flag == 1)
			//	return 2;   //this is a saddle triangle
			a_flag++;
		}
	}

	if(r_flag > 0 && a_flag == 0)
		return 0;
	else if(r_flag == 0 && a_flag > 0)
		return 1;
	else
		return 2;
}


/*
The number of the triangles in the saddle-like regions we want to make correction to
is 4 now. 03/08/07
Probably this should happen during the judgement
*/
int make_correction_saddle(int scc_index)
{
	int j;
	int source, sink, saddle;

	source=sink=saddle = 0;

	for(j = 0; j < scclist.scccomponents[scc_index].num_nodes; j++)
	{
		int t = scclist.scccomponents[scc_index].nodes[j];

		if(singularities[Object.flist[t]->singularity_index].type == SOURCE)
			source++;
		if(singularities[Object.flist[t]->singularity_index].type == SINK)
			sink++;
		if(singularities[Object.flist[t]->singularity_index].type == SADDLE)
			saddle++;
	}

	if(saddle == 0)
	{
		if(source > 0)
			return 0;
		else
			return 1;
	}
	else
		return 2;
}

/*implement region growing in the directed graph*/

void grow_saddle_region_graph();
void grow_repeller_region_graph(int scc_index, int);
void grow_attractor_region_graph(int scc_index);
void build_mcg_edges_graph();

int *repell_region = NULL;
int ntris_repell_region;
int curMaxNumTrisRepellRegion;

/*
set the nodes in the specified SCC as "visited"
*/
void set_sccnodes_flags_for_SCC(int scc_index)
{
	int i;

	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		sccnodes[scclist.scccomponents[scc_index].nodes[i]].visited = 1;
	}
}


/*
judge whether a new added link is proper or not.
The criterion for "proper" now is
source->saddle/sink
saddle->saddle/sink
The inputs are the node indices of the MCG
*/

bool is_valid_link(int from, int to)
{
	if(mcgnodes[from].type == 1 || mcgnodes[to].type == 0)
		return false;

	return true;
}


/*
grow the repeller region using the obtained directed graph
not sure the graph or the algorithm is correct here!!!
03/12/07
*/
void grow_repeller_region_graph(int scc_index, int inverse)
{
	/*probably we should use all the nodes in the subgraph here*/
	//get_boundary_nodes(scc_index);

	curMaxNumTrisRepellRegion = Object.nfaces;
	if(repell_region != NULL)
		free(repell_region);

	repell_region = (int*)malloc(sizeof(int)*curMaxNumTrisRepellRegion);
	ntris_repell_region = 0;

	/*copy the boundary nodes*/
	//CopyRegion(boundary_nodes,repell_region,nboundarytris);
	//ntris_repell_region = nboundarytris;

	
	/*
	use all the triangles in the region
	*/
	CopyRegion(scclist.scccomponents[scc_index].nodes, repell_region,
		scclist.scccomponents[scc_index].num_nodes);
	ntris_repell_region = scclist.scccomponents[scc_index].num_nodes;

	//int cur_index = 0;
	//int num_newnodes = 0;

	//do{
	//}while(num_newnodes > 0)

	//reset_sccnodes_flags();
	int i, j, cur_node, other_node;

	for(i = 0; i < ntris_repell_region; i++)
	{
		/*expand current node along its outgoing edges*/
		cur_node = repell_region[i];

		for(j = 0; j < sccnodes[cur_node].nedges; j++)
		{
			if(sccedges[sccnodes[cur_node].edges[j]].node_index2 == cur_node)
				continue;

			other_node = sccedges[sccnodes[cur_node].edges[j]].node_index2;

			if(sccnodes[other_node].sscomp_index != scc_index
				//&& sccnodes[other_node].visited == 0
				&& !IsRepeated(repell_region, other_node, ntris_repell_region))
			{
				/*judge whether it corresponds to other valid Morse set*/
				if((scclist.scccomponents[sccnodes[other_node].sscomp_index].num_nodes <= 2
					&& scclist.scccomponents[sccnodes[other_node].sscomp_index].num_singularities > 0)
					||(scclist.scccomponents[sccnodes[other_node].sscomp_index].num_nodes > 2
					&& scclist.scccomponents[sccnodes[other_node].sscomp_index].valid == 1))
				{
					/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
					/*add a new edge in mcg, will ignore double connection here !!!!!!!*/
					//if(!is_connected_mcg(cur_node, other_node)
					//	&& !has_interval_mcg(cur_node, other_node))
					//{
					//	add_to_mcg_edgelist(cur_node, other_node);

					//	/*add the edge to the corresponding nodes*/
					//	add_edge_to_node(cur_node, cur_mcgedge_index-1);
					//	add_edge_to_node(other_node, cur_mcgedge_index-1);

					//	FILE *fp = fopen("currentedges.txt", "w");
					//	fprintf(fp, "The MCG has %d edges now!\n", cur_mcgedge_index);
					//	fprintf(fp, "It is growing from SCC %d!\n", scc_index);
					//	fprintf(fp, "It is node %d in MCG!\n", scclist.scccomponents[scc_index].node_index);
					//	fclose(fp);
					//}
					
					if(!is_connected_mcg(scclist.scccomponents[scc_index].node_index, 
						scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index)
						//&& !has_interval_mcg(scclist.scccomponents[scc_index].node_index, 
						//scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index)
						)
					{
						if(inverse == 0)
						{
							if(is_valid_link(scclist.scccomponents[scc_index].node_index, 
								scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index))
							{
								add_to_mcg_edgelist(scclist.scccomponents[scc_index].node_index, 
									scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index);
								/*add the edge to the corresponding nodes*/
								add_edge_to_node(scclist.scccomponents[scc_index].node_index, cur_mcgedge_index-1);
								add_edge_to_node(scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index, 
									cur_mcgedge_index-1);
							}
						}
						else{
							if(is_valid_link(scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index, 
								scclist.scccomponents[scc_index].node_index))
							{
								add_to_mcg_edgelist(scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index, 
									scclist.scccomponents[scc_index].node_index);
								/*add the edge to the corresponding nodes*/
								add_edge_to_node(scclist.scccomponents[scc_index].node_index, cur_mcgedge_index-1);
								add_edge_to_node(scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index, 
									cur_mcgedge_index-1);
							}
						}

						/*add the edge to the corresponding nodes*/
						//add_edge_to_node(scclist.scccomponents[scc_index].node_index, cur_mcgedge_index-1);
						//add_edge_to_node(scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index, 
						//	cur_mcgedge_index-1);

						//FILE *fp = fopen("currentedges.txt", "w");
						//fprintf(fp, "The MCG has %d edges now!\n", cur_mcgedge_index);
						//fprintf(fp, "It is growing from SCC %d!\n", scc_index);
						//fprintf(fp, "It is node %d in MCG!\n", scclist.scccomponents[scc_index].node_index);
						//fclose(fp);

						//set_sccnodes_flags_for_SCC(sccnodes[other_node].sscomp_index);
						//continue;
					}
				}

				/*add the node into the array*/
				repell_region[ntris_repell_region] = other_node;
				ntris_repell_region++;
				
				//sccnodes[other_node].visited = 1;
			}
		}
	}
}




/*
grow the repelling region using the obtained directed graph
New method but still based BFS
03/27/07
*/
void grow_repeller_region_graph_2(int scc_index, int inverse)
{

	curMaxNumTrisRepellRegion = Object.nfaces;
	if(repell_region != NULL)
		free(repell_region);

	repell_region = (int*)malloc(sizeof(int)*curMaxNumTrisRepellRegion);
	ntris_repell_region = 0;

	/*
	use all the triangles in the region
	*/
	CopyRegion(scclist.scccomponents[scc_index].nodes, repell_region,
		scclist.scccomponents[scc_index].num_nodes);
	ntris_repell_region = scclist.scccomponents[scc_index].num_nodes;


	reset_sccnodes_flags();
	int i, j, k, cur_node, other_node;

	for(i = 0; i < ntris_repell_region; i++)
	{
		/*expand current node along its outgoing edges*/
		cur_node = repell_region[i];

		if(sccnodes[cur_node].visited == 1)
			continue;

		sccnodes[cur_node].visited = 1;

		for(j = 0; j < sccnodes[cur_node].nedges; j++)
		{
			if(sccedges[sccnodes[cur_node].edges[j]].node_index2 == cur_node)
				continue;

			other_node = sccedges[sccnodes[cur_node].edges[j]].node_index2;

			if(sccnodes[other_node].sscomp_index != scc_index
				&& sccnodes[other_node].visited == 0
				&& !IsRepeated(repell_region, other_node, ntris_repell_region))
			{
				/*judge whether it corresponds to other valid Morse set*/
				if((scclist.scccomponents[sccnodes[other_node].sscomp_index].num_nodes <= 2
					&& scclist.scccomponents[sccnodes[other_node].sscomp_index].num_singularities > 0)
					||(scclist.scccomponents[sccnodes[other_node].sscomp_index].num_nodes > 2
					&& scclist.scccomponents[sccnodes[other_node].sscomp_index].valid == 1))
				{
					if(!is_connected_mcg(scclist.scccomponents[scc_index].node_index, 
						scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index)
						//&& !has_interval_mcg(scclist.scccomponents[scc_index].node_index, 
						//scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index)
						)
					{
						if(inverse == 0)
						{
							if(is_valid_link(scclist.scccomponents[scc_index].node_index, 
								scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index))
							{
								add_to_mcg_edgelist(scclist.scccomponents[scc_index].node_index, 
									scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index);
								/*add the edge to the corresponding nodes*/
								add_edge_to_node(scclist.scccomponents[scc_index].node_index, cur_mcgedge_index-1);
								add_edge_to_node(scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index, 
									cur_mcgedge_index-1);
							}
						}
						else{
							if(is_valid_link(scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index, 
								scclist.scccomponents[scc_index].node_index))
							{
								add_to_mcg_edgelist(scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index, 
									scclist.scccomponents[scc_index].node_index);
								/*add the edge to the corresponding nodes*/
								add_edge_to_node(scclist.scccomponents[scc_index].node_index, cur_mcgedge_index-1);
								add_edge_to_node(scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index, 
									cur_mcgedge_index-1);
							}
						}

						/*add all the nodes in this Morse sets into the array*/
						for(k = 0; k < scclist.scccomponents[sccnodes[other_node].sscomp_index].num_nodes; k++)
						{
							if(sccnodes[scclist.scccomponents[sccnodes[other_node].sscomp_index].nodes[k]].visited
								== 1)
								continue;
							repell_region[ntris_repell_region] = 
								scclist.scccomponents[sccnodes[other_node].sscomp_index].nodes[k];
							ntris_repell_region++;
						}

					}
				}
				else{
					/*add the node into the array*/
					repell_region[ntris_repell_region] = other_node;
					ntris_repell_region++;
				}
			}
		}
	}
}




void grow_repeller_region_dual_graph(int scc_index, int inverse)
{
	get_boundary_nodes(scc_index);

	curMaxNumTrisRepellRegion = Object.nfaces;
	if(repell_region != NULL)
		free(repell_region);

	repell_region = (int*)malloc(sizeof(int)*curMaxNumTrisRepellRegion);

	/*copy the boundary nodes*/
	CopyRegion(boundary_nodes,repell_region,nboundarytris);
	ntris_repell_region = nboundarytris;

	//int cur_index = 0;
	//int num_newnodes = 0;

	//do{
	//}while(num_newnodes > 0)

	//reset_sccnodes_flags();

	int i, j, cur_node, other_node;
	for(i = 0; i < ntris_repell_region; i++)
	{
		/*expand current node along its outgoing edges*/
		cur_node = repell_region[i];

		for(j = 0; j < nodes_dual[cur_node].nedges; j++)
		{
			if(edges_dual[nodes_dual[cur_node].edges[j]].node_index2 == cur_node)
				continue;

			other_node = edges_dual[nodes_dual[cur_node].edges[j]].node_index2;

			if(nodes_dual[other_node].sscomp_index != scc_index
				//&& sccnodes[other_node].visited == 0
				&& !IsRepeated(repell_region, other_node, ntris_repell_region))
			{
				/*judge whether it corresponds to other valid Morse set*/
				if((scclist.scccomponents[nodes_dual[other_node].sscomp_index].num_nodes <= 2
					&& scclist.scccomponents[nodes_dual[other_node].sscomp_index].num_singularities > 0)
					||(scclist.scccomponents[nodes_dual[other_node].sscomp_index].num_nodes > 2
					&& scclist.scccomponents[nodes_dual[other_node].sscomp_index].valid == 1))
				{
					/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
					/*add a new edge in mcg, will ignore double connection here !!!!!!!*/
					//if(!is_connected_mcg(cur_node, other_node)
					//	&& !has_interval_mcg(cur_node, other_node))
					//{
					//	add_to_mcg_edgelist(cur_node, other_node);

					//	/*add the edge to the corresponding nodes*/
					//	add_edge_to_node(cur_node, cur_mcgedge_index-1);
					//	add_edge_to_node(other_node, cur_mcgedge_index-1);

					//	FILE *fp = fopen("currentedges.txt", "w");
					//	fprintf(fp, "The MCG has %d edges now!\n", cur_mcgedge_index);
					//	fprintf(fp, "It is growing from SCC %d!\n", scc_index);
					//	fprintf(fp, "It is node %d in MCG!\n", scclist.scccomponents[scc_index].node_index);
					//	fclose(fp);
					//}
					
					if(!is_connected_mcg(scclist.scccomponents[scc_index].node_index, 
						scclist.scccomponents[nodes_dual[other_node].sscomp_index].node_index)
						//&& !has_interval_mcg(scclist.scccomponents[scc_index].node_index, 
						//scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index)
						)
					{
						if(inverse == 0)
						{
							add_to_mcg_edgelist(scclist.scccomponents[scc_index].node_index, 
								scclist.scccomponents[nodes_dual[other_node].sscomp_index].node_index);
						}
						else
							add_to_mcg_edgelist(scclist.scccomponents[nodes_dual[other_node].sscomp_index].node_index, 
								scclist.scccomponents[scc_index].node_index);

						/*add the edge to the corresponding nodes*/
						add_edge_to_node(scclist.scccomponents[scc_index].node_index, cur_mcgedge_index-1);
						add_edge_to_node(scclist.scccomponents[nodes_dual[other_node].sscomp_index].node_index, 
							cur_mcgedge_index-1);

						//FILE *fp = fopen("currentedges.txt", "w");
						//fprintf(fp, "The MCG has %d edges now!\n", cur_mcgedge_index);
						//fprintf(fp, "It is growing from SCC %d!\n", scc_index);
						//fprintf(fp, "It is node %d in MCG!\n", scclist.scccomponents[scc_index].node_index);
						//fclose(fp);

						//set_sccnodes_flags_for_SCC(sccnodes[other_node].sscomp_index);
						//continue;
					}
				}

				/*add the node into the array*/
				repell_region[ntris_repell_region] = other_node;
				ntris_repell_region++;
				
				//sccnodes[other_node].visited = 1;
			}
		}
	}
}



/*grow the attracting region*/
void grow_attractor_region_graph(int scc_index)
{
	/*first reverse all the edges in the graph*/
	ReverseEdges();

	/*grow as repelling region*/
	grow_repeller_region_graph(scc_index, 1);
		//grow_repeller_region_dual_graph(scc_index, 1);

	/*reverse the edges back to original state*/
	ReverseEdges();
}


void grow_saddle_region_graph()
{
	int i;

	for(i = 0; i < cur_mcgnode_index; i++)
	{
		if(mcgnodes[i].type != 2)
			continue;

		/*first grow the repelling*/
		grow_repeller_region_graph(mcgnodes[i].scc_index, 0);
		//grow_repeller_region_dual_graph(mcgnodes[i].scc_index, 0);

		/*second grow the attracting*/
		//grow_attractor_region_graph(mcgnodes[i].scc_index);

		/*remove redundant*/
		//remove_redundant_edges();
	}
}






/*
Grow all the nodes forward without considering its type
03/27/07
*/
void grow_all_mcgnodes()
{
	int i;
	for(i = 0; i < cur_mcgnode_index; i++)
	{
		grow_repeller_region_graph_2(mcgnodes[i].scc_index, 0);
	}

	/*reverse the edges*/
 //   ReverseEdges();

	//
	///*grow based on reversed graph*/
	//for(i = 0; i < cur_mcgnode_index; i++)
	//{
	//	grow_repeller_region_graph_2(mcgnodes[i].scc_index, 1);
	//}

	//ReverseEdges();
}

extern bool del_one_edge_from(int *, int &, int);

/*
Remove the redundant edges (indirected).
Note: double connections may be removed by this routine as well
*/
void remove_redundant_edges()
{
	int i, j, k, l;

	int *expandnodes = new int[cur_mcgnode_index];
	int num_expanded_nodes;
	//int level = 0;
	int cur_node, other_node;

	for(i = 0; i < cur_mcgnode_index; i++)
	{
		/*search all the possible path between any two nodes in mcg*/
		if(mcgnodes[i].type == 1)
			continue;

		/*expand the node and remove redundant edges at the same time
		until it can not be expanded any more*/
		//level = 0;
		expandnodes[0] = i;
		num_expanded_nodes = 1;

		for(j = 0; j < num_expanded_nodes; j++)
		{
			cur_node = expandnodes[j];

			/*add all directly reachable nodes into the array*/
			for(k = 0; k < mcgnodes[cur_node].nedges; k++)
			{
				other_node = mcgedges[mcgnodes[cur_node].edges[k]].node_index2;

				/*we consider outgoing edges only*/
				if(other_node == cur_node)
					continue;

				if(IsRepeated(expandnodes, other_node, num_expanded_nodes))
				{
					/*if the node has been expanded before, it is a repeated node*/
					/*remove the edge between 'i' and this 'other_node'*/
					for(l = 0; l < mcgnodes[i].nedges; l++)
					{
						if(mcgedges[mcgnodes[i].edges[l]].node_index2 == other_node)
						{
							/*Remove this edge
							*/
							del_one_edge_from(mcgnodes[i].edges,
								mcgnodes[i].nedges, mcgnodes[i].edges[l]);
							
							del_one_edge_from(mcgnodes[other_node].edges,
								mcgnodes[other_node].nedges, mcgnodes[i].edges[l]);
							
							mcgedges[mcgnodes[i].edges[l]].cancel = true;
							
						    break;
						}
					}
				}

				else
				{
					expandnodes[num_expanded_nodes] = other_node;
					num_expanded_nodes++;
				}
			}

			//level ++;
		}
	}

	delete [] expandnodes;
}


int has_longer_path_in_mcg(int from, int to, int edgeindex)
{
	//the maximum searching number is the number of the edges in the graph
	int count = 0;
	int cur_node = from;
	int cur_node_id = 0;

	int i;
	//reset the 'visited' flags of all the edges
	for(i = 0; i < cur_mcgedge_index; i++)
	{
		mcgedges[i].visited = 0;
	}

	/*set the direct edge from 'from' to 'to' to be visited
	*/
	mcgedges[edgeindex].visited = 1;
	
	/* allocate the searching stack for DFS searching */
	int *searcharray = (int *)malloc(sizeof(int) * (cur_mcgnode_index*2));
	int MaxArrayElems = cur_mcgnode_index*2;
	int nelem_in_array = 0;

	searcharray[0] = from;
	nelem_in_array = 1;

	mcgnodes[from].parent = from;
	int num_triangles_in_path = 0;
	int pathlength = 0;

	/*Here we use BFS*/
	for(int j = 0; j < nelem_in_array; j++)
	{
		cur_node = searcharray[j];

		if(cur_node == to)
			goto LL;

		//pick one unvisited edges
		for(i = 0; i < mcgnodes[cur_node].nedges; i++)
		{
			if(mcgedges[mcgnodes[cur_node].edges[i]].visited == 1)
				continue;

			//set the node_index2 of the edge as the next cur_node
			//cur_node = sccedges[sccnodes[cur_node_id].edges[i]].node_index2;
			if(nelem_in_array >= MaxArrayElems) /*extend it*/
			{
				searcharray = (int *)realloc(searcharray, sizeof(int)*(MaxArrayElems+50));

				if(searcharray == NULL)
					exit(-1);

				MaxArrayElems += 50;
			}

			if(IsRepeated(searcharray,mcgedges[mcgnodes[cur_node].edges[i]].node_index2,nelem_in_array))
				continue;

			searcharray[nelem_in_array] = mcgedges[mcgnodes[cur_node].edges[i]].node_index2;
			mcgnodes[mcgedges[mcgnodes[cur_node].edges[i]].node_index2].parent = cur_node;
			nelem_in_array++;

			mcgedges[mcgnodes[cur_node].edges[i]].visited = 1;
		}
	}
			free(searcharray);
	return pathlength;

LL:		
		/*backtrack from the end triangle to the starting triangle*/
		cur_node = to;
		while(cur_node != from)
		{
			cur_node = mcgnodes[cur_node].parent;

			if(cur_node < 0)
			{
				free(searcharray);
				return -1;
			}
			pathlength ++;
		}


				free(searcharray);
		return pathlength;
}


bool has_longer_path_in_mcg_2(int n1, int n2, int edgeindex)
{
	//the maximum searching number is the number of the edges in the graph

	int i;
	//reset the 'visited' flags of all the edges
	for(i = 0; i < cur_mcgedge_index; i++)
	{
		mcgedges[i].visited = 0;
	}

	/*set the direct edge from 'from' to 'to' to be visited
	*/
	mcgedges[edgeindex].visited = 1;
	
	//the maximum searching number is the number of the edges in the graph
	int count = 0;
	int cur_node = n1;
	int cur_node_id = 0;

	if(searchstack != NULL)
		free(searchstack);

	searchstack = (int *)malloc(sizeof(int) * (cur_mcgnode_index*2));
	MaxStackElems = cur_mcgnode_index*2;

	searchstack[0] = n1;
	nelem_in_stack = 1;

	while(cur_node != n2 && count < cur_mcgedge_index)
	{
		/* pick the last element of the stack */
		cur_node = searchstack[nelem_in_stack-1];

		nelem_in_stack --;

		if(cur_node == n2)
			return true;

		/***********************************************************/

		//pick one unvisited edges
		for(i = 0; i < mcgnodes[cur_node].nedges; i++)
		{
			if(mcgedges[mcgnodes[cur_node].edges[i]].visited == 1
				||mcgedges[mcgnodes[cur_node].edges[i]].cancel)
				continue;

			/*consider outgoing edges only*/
			if(mcgedges[mcgnodes[cur_node].edges[i]].node_index2 == cur_node)
				continue;

			//set the node_index2 of the edge as the next cur_node
			//cur_node = sccedges[sccnodes[cur_node_id].edges[i]].node_index2;
			if(nelem_in_stack >= MaxStackElems) /*extend it*/
			{
				searchstack = (int *)realloc(searchstack, sizeof(int)*(MaxStackElems+50));

				if(searchstack == NULL)
					exit(-1);

				MaxStackElems += 50;
			}

			searchstack[nelem_in_stack] = mcgedges[mcgnodes[cur_node].edges[i]].node_index2;
			nelem_in_stack++;

			mcgedges[mcgnodes[cur_node].edges[i]].visited = 1;
			count++;
			//break;
		}

		if(nelem_in_stack == 0)
			return false;

	}

	if(count < cur_mcgedge_index)
	{
		return true;
	}
	return false;
}

/*
New routine of remove redundant edges
*/
void remove_redundant_edges_2()
{
	int i, j;
	int *expandnodes = new int[cur_mcgnode_index];
	int num_expanded_nodes;
	//int level = 0;
	int cur_node, other_node;

	for(i = 0; i < cur_mcgnode_index; i++)
	{
		for(j = 0; j < mcgnodes[i].nedges; j++)
		{
			other_node = mcgedges[mcgnodes[i].edges[j]].node_index2;
			if(other_node == i)
				continue;  /*consider outgoing edges only*/

			/*currently, we have remove double connections 03/14/07*/
			//if(has_longer_path_in_mcg(i, other_node, mcgnodes[i].edges[j]) > 1)
			if(has_longer_path_in_mcg_2(i, other_node, mcgnodes[i].edges[j]))
			{
					/*Remove this edge
					*/
					//del_one_edge_from(mcgnodes[i].edges,
					//	mcgnodes[i].nedges, mcgnodes[i].edges[j]);
					//
					//del_one_edge_from(mcgnodes[other_node].edges,
					//	mcgnodes[other_node].nedges, mcgnodes[i].edges[j]);
					
					mcgedges[mcgnodes[i].edges[j]].cancel = true;
			}
		}
	}

	/*update the edge list in the nodes*/
	for(i = 0; i < cur_mcgedge_index; i++)
	{
		if(mcgedges[i].cancel)
		{
			del_one_edge_from(mcgnodes[mcgedges[i].node_index1].edges,
				mcgnodes[mcgedges[i].node_index1].nedges, i);
			
			del_one_edge_from(mcgnodes[mcgedges[i].node_index2].edges,
				mcgnodes[mcgedges[i].node_index2].nedges, i);
		}
	}
	delete [] expandnodes;
}


/*
Grow the remaining regions
*/
void grow_remaining_graph()
{
	int i;

	for(i = 0; i < cur_mcgnode_index; i++)
	{
		if(mcgnodes[i].type == 2)
			continue;

		if(mcgnodes[i].type == 0)
		{
			grow_repeller_region_graph(mcgnodes[i].scc_index, 0);
		//grow_repeller_region_dual_graph(mcgnodes[i].scc_index, 0);
		}

		else
		{
			grow_attractor_region_graph(mcgnodes[i].scc_index);
		}
		
		/*remove redundant*/
		remove_redundant_edges_2();
	}
}

/*
build the edges in mcg according to the obtained directed graph
*/

extern void build_dual_DirGraph_no_Tau();

void build_mcg_edges_graph()
{
	//build_dual_DirGraph_no_Tau(); /*build the dual graph*/

	cur_mcgedge_index = 0;

	/*first get the edges from saddle regions*/
	//grow_saddle_region_graph();

	//remove_redundant_edges();

	/*for each repelling or attracing region, if it is a ring-shaped region
	and only have one edge incident to its corresponding node in MCG*/
	//grow_remaining_graph();

	/*grow all nodes without knowing their types*/
	grow_all_mcgnodes();

	/*remove any redundant edges*/
	remove_redundant_edges_2();
}

/*build the mcg according to the obtained Morse sets*/
void build_mcg()
{
	assign_mcgnodes();
	layout_mcg();
	//build_mcg_edges_1();
	build_mcg_edges_graph();


	/*calculate the regions containing the connectors of the MCG 08/29/2007*/
	cal_mcgedge_regions();
}


/************************************************************************************/

/*we find the regions that can visualize the connection relationship between Morse sets
08/29/2007*/

/*we will use the region intersection operation*/
//extern void 
void cal_forward_region_MFIG(int scc_index, int inverse, 
							 int *repell_region, int &ntris_repell_region, int &curMaxNumTrisRepellRegion);
void get_saddle_two_regions(int scc_index, 
	int *forward_region, int &nforward_tris, int &curMaxNumForward,
	int *backward_region, int &nbackward_tris, 	int &curMaxNumBackward);

extern void ReverseEdges();

void intersect_twoRegions(int *source1, int num1, 
					  int *source2, int num2,
					  int *dest, int &num_dest)
{
	int i, j;
	int cur_index = 0;
	int cur_triangle;


	for(i = 0; i < num1; i++)
	{
		cur_triangle = source1[i];

		for(j = 0; j < num2; j++)
		{
			if(cur_triangle == source2[j])
			{
				dest[cur_index] = cur_triangle;
				cur_index++;
			}
		}

	}

	num_dest = cur_index;
}


void remove_redundant(int *source1, int &num1, int *source2, int num2)
{
	int i, j;
	int *temp;
	int cur_elems = 0;
	if(num1 > num2) temp=(int*)malloc(sizeof(int)*num1);
	else temp = (int*)malloc(sizeof(int)*num2);
	for(i=0; i<num1; i++)
	{
		for(j=0; j<num2; j++)
		{
			if(source1[i]==source2[j])
				break;
		}

		if(j>=num2)
		{
			temp[cur_elems] = source1[i];
			cur_elems ++;
		}
	}

	for(i=0; i<cur_elems; i++)
		source1[i] = temp[i];
	num1 = cur_elems;
}


void cal_mcgedge_regions()
{
	/* The following provides the description of the algorithm:

	for each saddle like Morse set in MCG
		if it has edges associated with it
		compute the regions of it using forward and backward graph of MFIG
	      consider all its edges
			  case 1: if it is saddle to sink
			          grow the sink region using backward graph of MFIG
					  intersect the obtained region with the forward region of the saddle
					  save the intersect regions to the member variable *triangles of MCGEdge
			  case 2: if it is source to saddle
			          grow the source region using forward graph of MFIG
					  intersect the obtained region with the backward region of the saddle
					  save the intersect regions to the member variable *triangles of MCGEdge
			  case 3: if it is saddle to saddle
			          according to the direction of connection
					  use the forward region of the starting saddle to intersect
					  with the backward region of the ending saddle
					  save the intersect regions to the member variable *triangles of MCGEdge					  
	*/

	/*declare variables*/

	int *forward_region = NULL;
	int curMaxNumForward;
	int nforward_tris = 0;
	int *backward_region = NULL;
	int nbackward_tris = 0;
	int curMaxNumBackward;
	int *other_region = NULL;
	int curMaxNumOther;
	int nother_tris = 0;
	int *intersect_region = NULL;
	int curMaxNumIntersect;
	int nintersect_tris = 0;

	/*Allocate the initial space for those variables*/
	curMaxNumForward = Object.nfaces;
	forward_region = (int*)malloc(sizeof(int)*curMaxNumForward);
	
	curMaxNumBackward = Object.nfaces;
	backward_region = (int*)malloc(sizeof(int)*curMaxNumBackward);
	
	curMaxNumOther = Object.nfaces;
	other_region = (int*)malloc(sizeof(int)*curMaxNumOther);

	curMaxNumIntersect = Object.nfaces;
	intersect_region = (int*)malloc(sizeof(int)*curMaxNumIntersect);

	/*initialize the mcg edges*/
	int i, j;

	for(i=0; i<cur_mcgedge_index; i++)
	{
		mcgedges[i].visited = 0;
		mcgedges[i].ntris = 0;
	}

	/**/

	bool forward_backward = false;

	int othernode;

	for(i=0; i<cur_mcgnode_index; i++)
	{
		if(mcgnodes[i].type != 2)
			continue;

		/**/
		if(mcgnodes[i].nedges == 0)
			continue;

		nforward_tris = nbackward_tris = 0;

		/*grow two regions for the saddle*/
		get_saddle_two_regions(mcgnodes[i].scc_index, forward_region, nforward_tris, curMaxNumForward, 
			backward_region, nbackward_tris, curMaxNumBackward);

		for(j=0; j<mcgnodes[i].nedges; j++)
		{
			if(mcgedges[mcgnodes[i].edges[j]].visited == 1)
				continue;

			mcgedges[mcgnodes[i].edges[j]].visited = 1;

			othernode = mcgedges[mcgnodes[i].edges[j]].node_index2;
			if(othernode == i) /*if it is current saddle node*/
				othernode = mcgedges[mcgnodes[i].edges[j]].node_index1;

			nintersect_tris = nother_tris = 0;

			/*case 1: the other node is a sink*/
			if(mcgnodes[othernode].type == 1) 
			{
				/*grow backward*/

				cal_forward_region_MFIG(mcgnodes[othernode].scc_index, 1,
					other_region, nother_tris, curMaxNumOther);

				forward_backward = true;

			}

			/*case 2: the other node is a source*/
			else if(mcgnodes[othernode].type == 0)
			{
				/*grow forward*/
				cal_forward_region_MFIG(mcgnodes[othernode].scc_index, 0,
					other_region, nother_tris, curMaxNumOther);

				forward_backward = false;
			}

			/*case 3: the other node is a saddle*/
			else
			{
				/*use the direction to decide*/
				if(othernode == mcgedges[mcgnodes[i].edges[j]].node_index1)
				{
					/*othernode is a source*/
					/*grow forward from it*/
					cal_forward_region_MFIG(mcgnodes[othernode].scc_index, 0,
						other_region, nother_tris, curMaxNumOther);

					/*compute the intersection*/

					forward_backward = false;
				}
				else
				{
					/*othernode is a sink*/
					/*grow backward from it*/
					cal_forward_region_MFIG(mcgnodes[othernode].scc_index, 1,
						other_region, nother_tris, curMaxNumOther);
					
					/*compute the intersection*/
					forward_backward = true;
				}
			}

			/*compute the intersection*/

			if(!forward_backward) /*forward, should use the backward region of the saddle*/
					intersect_twoRegions(other_region, nother_tris,
						backward_region, nbackward_tris,
						intersect_region, nintersect_tris);
					//intersect_twoRegions(other_region, nother_tris,
					//	forward_region, nforward_tris,
					//	intersect_region, nintersect_tris);
			else
					intersect_twoRegions(other_region, nother_tris,
						forward_region, nforward_tris,
						intersect_region, nintersect_tris);
					//intersect_twoRegions(other_region, nother_tris,
					//	backward_region, nbackward_tris,
					//	intersect_region, nintersect_tris);


			/*remove the redundant triangles*/
			CopyRegion(scclist.scccomponents[mcgnodes[othernode].scc_index].nodes, other_region,
				scclist.scccomponents[mcgnodes[othernode].scc_index].num_nodes);
			nother_tris = scclist.scccomponents[mcgnodes[othernode].scc_index].num_nodes;

			remove_redundant(intersect_region, nintersect_tris,
				other_region, nother_tris);

			CopyRegion(scclist.scccomponents[mcgnodes[i].scc_index].nodes, other_region,
				scclist.scccomponents[mcgnodes[i].scc_index].num_nodes);
			nother_tris = scclist.scccomponents[mcgnodes[i].scc_index].num_nodes;

			remove_redundant(intersect_region, nintersect_tris,
				other_region, nother_tris);

			/*copy the intersection region to the edge*/

			mcgedges[mcgnodes[i].edges[j]].triangles = 
				(int *)malloc(sizeof(int)*nintersect_tris);
			CopyRegion(intersect_region, mcgedges[mcgnodes[i].edges[j]].triangles, nintersect_tris);
			mcgedges[mcgnodes[i].edges[j]].ntris = nintersect_tris;
		}
	}

	free(forward_region);
	free(backward_region);
	free(intersect_region);
	free(other_region);
}


/*
grow the repelling region using the obtained directed graph
New method but still based BFS
03/27/07
*/
void cal_forward_region_MFIG(int scc_index, int inverse, 
							 int *repell_region, int &ntris_repell_region, int &curMaxNumTrisRepellRegion)
{

	//curMaxNumTrisRepellRegion = Object.nfaces;
	//if(repell_region != NULL)
	//	free(repell_region);

	//repell_region = (int*)malloc(sizeof(int)*curMaxNumTrisRepellRegion);

	ntris_repell_region = 0;

	if(inverse == 1)
		ReverseEdges();


	/*
	use all the triangles in the region
	*/
	CopyRegion(scclist.scccomponents[scc_index].nodes, repell_region,
		scclist.scccomponents[scc_index].num_nodes);
	ntris_repell_region = scclist.scccomponents[scc_index].num_nodes;



	reset_sccnodes_flags();
	int i, j, k, cur_node, other_node;

	for(i = 0; i < ntris_repell_region; i++)
	{
		/*expand current node along its outgoing edges*/
		cur_node = repell_region[i];

		if(sccnodes[cur_node].visited == 1)
			continue;

		sccnodes[cur_node].visited = 1;

		for(j = 0; j < sccnodes[cur_node].nedges; j++)
		{
			if(sccedges[sccnodes[cur_node].edges[j]].node_index2 == cur_node)
				continue;

			other_node = sccedges[sccnodes[cur_node].edges[j]].node_index2;

			if(sccnodes[other_node].sscomp_index != scc_index
				&& sccnodes[other_node].visited == 0
				&& !IsRepeated(repell_region, other_node, ntris_repell_region))
			{
				///*judge whether it corresponds to other valid Morse set*/
				//if((scclist.scccomponents[sccnodes[other_node].sscomp_index].num_nodes <= 2
				//	&& scclist.scccomponents[sccnodes[other_node].sscomp_index].num_singularities > 0)
				//	||(scclist.scccomponents[sccnodes[other_node].sscomp_index].num_nodes > 2
				//	&& scclist.scccomponents[sccnodes[other_node].sscomp_index].valid == 1))
				//{
				//	if(!is_connected_mcg(scclist.scccomponents[scc_index].node_index, 
				//		scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index)
				//		//&& !has_interval_mcg(scclist.scccomponents[scc_index].node_index, 
				//		//scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index)
				//		)
				//	{
				//		//if(inverse == 0)
				//		//{
				//		//	if(is_valid_link(scclist.scccomponents[scc_index].node_index, 
				//		//		scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index))
				//		//	{
				//		//	}
				//		//}
				//		//else{
				//		//	if(is_valid_link(scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index, 
				//		//		scclist.scccomponents[scc_index].node_index))
				//		//	{
				//		//		//add_to_mcg_edgelist(scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index, 
				//		//		//	scclist.scccomponents[scc_index].node_index);
				//		//		///*add the edge to the corresponding nodes*/
				//		//		//add_edge_to_node(scclist.scccomponents[scc_index].node_index, cur_mcgedge_index-1);
				//		//		//add_edge_to_node(scclist.scccomponents[sccnodes[other_node].sscomp_index].node_index, 
				//		//		//	cur_mcgedge_index-1);
				//		//	}
				//		//}

				//		/*add all the nodes in this Morse sets into the array*/
						//for(k = 0; k < scclist.scccomponents[sccnodes[other_node].sscomp_index].num_nodes; k++)
						//{
						//	if(sccnodes[scclist.scccomponents[sccnodes[other_node].sscomp_index].nodes[k]].visited
						//		== 1 || IsRepeated(repell_region, 
						//		scclist.scccomponents[sccnodes[other_node].sscomp_index].nodes[k], ntris_repell_region))
						//		continue;
						//	repell_region[ntris_repell_region] = 
						//		scclist.scccomponents[sccnodes[other_node].sscomp_index].nodes[k];
						//	ntris_repell_region++;
						//}

				//	}
				//}
				//else{
					/*add the node into the array*/
					repell_region[ntris_repell_region] = other_node;
					ntris_repell_region++;
				//}
			}
		}
	}
	
	if(inverse == 1)
		ReverseEdges();
}



/**/
void get_saddle_two_regions(int scc_index, 
	int *forward_region, int &nforward_tris, int &curMaxNumForward, 
	int *backward_region, int &nbackward_tris, 	int &curMaxNumBackward)
{
	/*grow the forward region*/
	cal_forward_region_MFIG(scc_index, 0, forward_region, nforward_tris, curMaxNumForward);

	/*grow the backward region*/
	/*reverse the edges of MFIG*/

	cal_forward_region_MFIG(scc_index, 1, backward_region, nbackward_tris, curMaxNumForward);
	
	/*reverse the edges of MFIG*/
}