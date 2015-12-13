
/*
FindSCC.cpp
*/
#include "stdafx.h"

#include "FindSCC.h"

#include "topologyedit.h"

#include "VFDataStructure.h"

#include "FindSepAttPoints.h"

GraphEdge *sccedges;
int num_sccedges = 0;
int curMaxNumDirGraphEdges;
GraphNode2 *sccnodes;
int *cur_nodes_order;
int num_sccnodes = 0;
int num_sccomps = 0;               //record the current number of SCC components
int cur_sccnode_index = 0;
int cur_end_edgelist;                 /*the last index of the edge in the list*/


/*build a dual graph 03/12/07*/
GraphEdge *edges_dual = NULL;
int num_edges_dual = 0;
//int curMaxNumDirGraphEdges;
GraphNode2 *nodes_dual = NULL;
int num_nodes_dual = 0;

int num_morsesets = 0;

SCCList scclist;

double MixedEdgeRatio = 0.1;

extern TriangularRegion repellerRegion;       ////region containing a repeller

extern Polygon3D Object;



int *path = NULL;
int num_tris_inpath = 0;
int *searchstack = NULL;
int nelem_in_stack;
int MaxStackElems;
int glob_t1, glob_t2;


/**/
/*it is better to record the time*/
time_t g_rawtime;

extern void CalNormalAtEdge(Edge *cur_edge, Face *face, int type);

/* Test the \tau recurrence 01/10/07 */
extern void trace_Mesh_Tau(double);

extern void construct_tau_map(double tau);


bool IsExitEdge(Edge *cur_edge)
{
	double dot1, dot2, dotresult = 0;
	Vertex *v1, *v2;
	v1 = Object.vlist[cur_edge->verts[0]];
	v2 = Object.vlist[cur_edge->verts[1]];

	dot1 = dot(cur_edge->repell_normal, v1->vec);
	dot2 = dot(cur_edge->repell_normal, v2->vec);

	dotresult = min(dot1, dot2); //pure exit edge growing

	if(dotresult >= 0) return true;

	return false;
}

bool IsEntranceEdge(Edge *cur_edge)
{
	double dot1, dot2, dotresult = 0;
	Vertex *v1, *v2;
	v1 = Object.vlist[cur_edge->verts[0]];
	v2 = Object.vlist[cur_edge->verts[1]];

	dot1 = dot(cur_edge->repell_normal, v1->vec);
	dot2 = dot(cur_edge->repell_normal, v2->vec);

	dotresult = max(dot1, dot2); //pure entrance edge growing

	if(dotresult < 0) return true;

	return false;
}


bool IsMixedEdge(Edge *cur_edge)
{
	double dot1, dot2, dotresult = 0;
	Vertex *v1, *v2;
	v1 = Object.vlist[cur_edge->verts[0]];
	v2 = Object.vlist[cur_edge->verts[1]];

	dot1 = dot(cur_edge->repell_normal, v1->vec);
	dot2 = dot(cur_edge->repell_normal, v2->vec);
	
	double ratio;
	ratio = (cur_edge->repell_normal.entry[0]*v1->vec.entry[0]+v1->vec.entry[1]*cur_edge->repell_normal.entry[1])/
		((v1->vec.entry[0]-v2->vec.entry[0])*cur_edge->repell_normal.entry[0]+(v1->vec.entry[1]-v2->vec.entry[1])*cur_edge->repell_normal.entry[1]);

	//if(ratio < 0.2 || ratio > 0.8)
	if(ratio < MixedEdgeRatio || ratio > 1-MixedEdgeRatio)
		return false;

	if(dot1*dot2 < 0)  //do we need to consider tangent case here "=" 07/10/06
		return true;

	return false;
}

void build_DirGraph_no_Tau()
{
	int i;

	for(i = 0; i < Object.nfaces; i++)
	{
		build_DirGraph_Tri_no_Tau(i,0);
	}
}


void build_DirGraph_Tri_no_Tau(int tri, int backward)
{
	Face *face = Object.flist[tri];
	Edge *cur_e;
	Vertex *cur_v;
	
	int node1, node2;

	//Add the face to the node list
	//sccnodes[num_sccnodes].node_index = tri;
	//cur_nodes_order[num_sccnodes] = tri;
	//num_sccnodes ++;

	for(int j = 0; j < 3; j++)
	{
		cur_e = face->edges[j];
		//cur_e->mixed = 0;
		
		////Calculate the outward normal of the edge
		CalNormalAtEdge(cur_e, face, 0);

		//build the directed connection according to the type of the edge
		if(backward == 0)
		{
			if(IsExitEdge(cur_e)||IsMixedEdge(cur_e))
			{
				node1 = tri;   //go from currrent triangle
				node2 = cur_e->tris[0];
				if(node2 == node1)
					node2 = cur_e->tris[1];
				
				if(node2 < 0)  // we reach boundary
					continue;

				//add to the edge list
				SCC_AddToEdge(node1, node2, num_sccedges);
				//add the edge to the nodes
				SCC_AddEdgeToNode(node1, num_sccedges-1);
				SCC_AddEdgeToNode(node2, num_sccedges-1);

				if(IsMixedEdge(cur_e))
					cur_e->mixed = 1;  //mixed edge

				else
					cur_e->mixed = 2;  //exit edge
			}
		}
		else if(backward == 1)
		{
			if(IsEntranceEdge(cur_e)||IsMixedEdge(cur_e))
			{
				node1 = tri;   //go from currrent triangle
				node2 = cur_e->tris[0];
				if(node2 == node1)
					node2 = cur_e->tris[1];
				
				if(node2 < 0)  // we reach boundary
					continue;

				//add to the edge list
				SCC_AddToEdge(node1, node2, num_sccedges);
				//add the edge to the nodes
				SCC_AddEdgeToNode(node1, num_sccedges-1);
				SCC_AddEdgeToNode(node2, num_sccedges-1);

				if(IsMixedEdge(cur_e))
					cur_e->mixed = 1;  //mixed edge

				else
					cur_e->mixed = 2;  //exit edge
			}
		}

	}
	
	face->repell_inregion = 0;

}


/*********************************************************************/
/* build the dual graph for debugging here*/

/*-----------------------------------------------------------
* Add the edge "edgeindex" to the edge list of "node"
-----------------------------------------------------------*/
void add_EdgeToNode_dual(int node, int edgeindex)
{
	//we need to first extend the space of the edge list
	if(nodes_dual[node].nedges == 0)
	{
		if(nodes_dual[node].edges != NULL)
			free(nodes_dual[node].edges);
		nodes_dual[node].edges = (int*)malloc(sizeof(int)*1);
	}

	else
	{
		int *temp = nodes_dual[node].edges;

		if((nodes_dual[node].edges = 
			(int*)realloc((int*) nodes_dual[node].edges, sizeof(int) * (nodes_dual[node].nedges+1))) == NULL)
		{
			MessageBox(NULL, "fail to reallocate for edge list of node!", "Error", MB_OK);
			exit(-1);
			//sccnodes[node].edges = temp;
		}
	}

	nodes_dual[node].edges[nodes_dual[node].nedges] = edgeindex;
	nodes_dual[node].nedges ++;
}


/*-----------------------------------------------------------
* add node1 and node2 to the edge list
We always assume that the direction is node1 to node2
-----------------------------------------------------------*/
void add_ToEdge_dual(int node1, int node2, int &cur_index)
{
	if(cur_index >= curMaxNumDirGraphEdges)
	{
		GraphEdge *temp = edges_dual;
		edges_dual = (GraphEdge *)realloc(edges_dual, sizeof(GraphEdge)*(curMaxNumDirGraphEdges+Object.nedges/2));
		if(edges_dual == NULL)
		{
			MessageBox(NULL, "failed to reallocate memory!", "Error", MB_OK);
			exit(-1);
		}
		curMaxNumDirGraphEdges += Object.nedges/2;
	}

	edges_dual[cur_index].node_index1 = node1;
	edges_dual[cur_index].node_index2 = node2;
	edges_dual[cur_index].edge_index = cur_index;
	cur_index++;
}



void build_dual_DirGraph_Tri_no_Tau(int tri, int backward)
{
	Face *face = Object.flist[tri];
	Edge *cur_e;
	Vertex *cur_v;
	
	int node1, node2;

	//Add the face to the node list
	//sccnodes[num_sccnodes].node_index = tri;
	//cur_nodes_order[num_sccnodes] = tri;
	//num_sccnodes ++;

	for(int j = 0; j < 3; j++)
	{
		cur_e = face->edges[j];
		//cur_e->mixed = 0;
		
		////Calculate the outward normal of the edge
		CalNormalAtEdge(cur_e, face, 0);

		//build the directed connection according to the type of the edge
		if(backward == 0)
		{
			if(IsExitEdge(cur_e)||IsMixedEdge(cur_e))
			{
				node1 = tri;   //go from currrent triangle
				node2 = cur_e->tris[0];
				if(node2 == node1)
					node2 = cur_e->tris[1];
				
				if(node2 < 0)  // we reach boundary
					continue;

				//add to the edge list
				add_ToEdge_dual(node1, node2, num_edges_dual);
				//add the edge to the nodes
				add_EdgeToNode_dual(node1, num_edges_dual-1);
				add_EdgeToNode_dual(node2, num_edges_dual-1);

				if(IsMixedEdge(cur_e))
					cur_e->mixed = 1;  //mixed edge

				else
					cur_e->mixed = 2;  //exit edge
			}
		}
		else if(backward == 1)
		{
			if(IsEntranceEdge(cur_e)||IsMixedEdge(cur_e))
			{
				node1 = tri;   //go from currrent triangle
				node2 = cur_e->tris[0];
				if(node2 == node1)
					node2 = cur_e->tris[1];
				
				if(node2 < 0)  // we reach boundary
					continue;

				//add to the edge list
				add_ToEdge_dual(node1, node2, num_edges_dual);
				//add the edge to the nodes
				add_EdgeToNode_dual(node1, num_edges_dual-1);
				add_EdgeToNode_dual(node2, num_edges_dual-1);

				if(IsMixedEdge(cur_e))
					cur_e->mixed = 1;  //mixed edge

				else
					cur_e->mixed = 2;  //exit edge
			}
		}

	}
	
	face->repell_inregion = 0;

}


void build_dual_DirGraph_no_Tau()
{
	int i;

	if(nodes_dual != NULL)
	{
		for(i = 0; i < num_nodes_dual; i++)
			if(nodes_dual[i].edges!=NULL)
				free(nodes_dual[i].edges);
		free(nodes_dual);
	}
	nodes_dual = (GraphNode2*)malloc(sizeof(GraphNode2)*Object.nfaces);
	num_nodes_dual = 0;

	if(edges_dual != NULL)
		free(edges_dual);
	edges_dual = (GraphEdge *)malloc(sizeof(GraphEdge)*Object.nedges*2);
	num_edges_dual = 0;
		
	/*copy the SCC information for the nodes*/
	for(i = 0; i < Object.nfaces; i++)
	{
		nodes_dual[num_nodes_dual].node_index = i;
		nodes_dual[num_nodes_dual].sscomp_index = sccnodes[i].sscomp_index;
		nodes_dual[num_nodes_dual].edges = NULL;
		nodes_dual[num_nodes_dual].nedges = 0;
		num_nodes_dual ++;
	}


	for(i = 0; i < Object.nfaces; i++)
	{
		build_dual_DirGraph_Tri_no_Tau(i,0);
	}

}

/*************************************************************************/
extern void cal_Jacobian_All_Vers();

/*
Build the connectiong graph
*/
void BuildConnectedGraph(double tau)
{
	int i, j;

	Face *face;
	Edge *cur_e;
	Vertex *cur_v;

	int node1, node2;

	num_sccedges = 0;
	num_sccnodes = 0;

	/* calculate the Jacobian everywhere in the mesh */
	CalJacobianForWholeMesh();
	cal_Jacobian_All_Vers();
	
	/*The following means that we find all the separation and attachement features
	in the flow first
	This is proven not working well!!*/
	//CalAllSpecialPoints_alltraingles(); /*calculate all the special points first*/

	//initialization
	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		for(j = 0; j < 3; j++)
		{
			cur_e = face->edges[j];
			cur_e->OnBoundary = 0;
			cur_e->mixed = 0;
		}
	}

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		//Add the face to the node list
		sccnodes[num_sccnodes].node_index = i;
		cur_nodes_order[num_sccnodes] = i;
		num_sccnodes ++;
		cur_sccnode_index++;
	}

	if(tau < 1e-4)
	{
		build_DirGraph_no_Tau();
	}

	else{
		/* Test \tau recurrence mapping 01/10/07 */
		
		//trace_Mesh_Tau(tau); 

		construct_tau_map(tau);

		/* Here, we need to add more nodes corresponding to the boundary edges 02/09/07 */

		/* we need to extend the node list according to the number of boundary edges */

		/* But the problem is after we get the SCC, how can we deal with the nodes
		that correspond to boundary edges ?
		Can we simply ignore them ?
		*/

		/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
		/* The following codes cause memeory leakage issue 02/12/07 !*/

		//Edge *e;
		//int num_boundaryedges = 0;
		//for(i = 0; i < Object.nfaces; i++)
		//{
		//	face = Object.flist[i];

		//	for(j = 0; j < 3; j++)
		//	{
		//		e = face->edges[j];
		//		
		//		if(e->tris[0] < 0 || e->tris[1] < 0) /* this is a boundary edge*/
		//		{
		//			num_boundaryedges++;
		//		}
		//	}
		//}

		///* extend the node list */

		//sccnodes = (GraphNode2 *)realloc(sccnodes, sizeof(GraphNode2)*(Object.nfaces + num_boundaryedges));
		//if(sccnodes == NULL)
		//{
		//	MessageBox(NULL, "fail to realloc sccnodes!", "Error", MB_OK);
		//	exit(-1);
		//}

		///* allocate a node for each edge */
		//for(i = 0; i < Object.nfaces; i++)
		//{
		//	face = Object.flist[i];

		//	for(j = 0; j < 3; j++)
		//	{
		//		e = face->edges[j];
		//		
		//		if(e->tris[0] < 0 || e->tris[1] < 0) /* this is a boundary edge*/
		//		{
		//			sccnodes[num_sccnodes].node_index =Object.nfaces + e->index;//set the special index of the node
		//			cur_nodes_order[num_sccnodes] = num_sccnodes; 
		//			num_sccnodes ++;
		//		}
		//	}
		//}
		
		//trace_Mesh_Tau(30);   //we use \tau/2 now, after two tracing, it should equivalent to \tau


		//fp_n = fopen("num_edges.txt", "a");
		//fprintf(fp_n, "current number of edges in the directed graph is :%d. \n", num_sccedges);
		//fclose(fp_n);
	}
}


/*-----------------------------------------------------------
* Add the edge "edgeindex" to the edge list of "node"
-----------------------------------------------------------*/
void SCC_AddEdgeToNode(int node, int edgeindex)
{
	//we need to first extend the space of the edge list
	if(sccnodes[node].nedges == 0)
	{
		if(sccnodes[node].edges != NULL)
			free(sccnodes[node].edges);
		sccnodes[node].edges = (int*)malloc(sizeof(int)*1);
	}

	else
	{
		int *temp = sccnodes[node].edges;

		if((sccnodes[node].edges = 
			(int*)realloc((int*) sccnodes[node].edges, sizeof(int) * (sccnodes[node].nedges+1))) == NULL)
		{
			MessageBox(NULL, "fail to reallocate for edge list of node!", "Error", MB_OK);
			exit(-1);
			//sccnodes[node].edges = temp;
		}
	}

	sccnodes[node].edges[sccnodes[node].nedges] = edgeindex;
	sccnodes[node].nedges ++;
}


/*-----------------------------------------------------------
* add node1 and node2 to the edge list
We always assume that the direction is node1 to node2
if this routine is called, it means there is no
unused element in the edge list now 03/20/07
-----------------------------------------------------------*/
void SCC_AddToEdge(int node1, int node2, int &cur_index)
{
	if(cur_index >= curMaxNumDirGraphEdges)
	{
		GraphEdge *temp = sccedges;
		sccedges = (GraphEdge *)realloc(sccedges, sizeof(GraphEdge)*(curMaxNumDirGraphEdges+Object.nedges/2));
		if(sccedges == NULL)
		{
			MessageBox(NULL, "failed to reallocate memory!", "Error", MB_OK);
			exit(-1);
		}
		curMaxNumDirGraphEdges += Object.nedges/2;
	}

	sccedges[cur_index].node_index1 = node1;
	sccedges[cur_index].node_index2 = node2;
	sccedges[cur_index].edge_index = cur_index;
	sccedges[cur_index].cancelled = false;  /*03/15/07*/
	cur_index++;

	cur_end_edgelist = cur_index;   /*record the last index of the edge 03/20/07 */
}



void ReverseEdges()
{
	int i, temp;

	for(i = 0; i < num_sccedges; i++)
	{
		temp = sccedges[i].node_index1;
		sccedges[i].node_index1 = sccedges[i].node_index2;
		sccedges[i].node_index2 = temp;
	}
}


/* Sort the index of the nodes according to their finished times */
void SortNodes(GraphNode2 *nodes, int *nodeindices, int num_nodes)
{
    int i, j, temp, max;

	for (i = 0; i < num_nodes-1; i++)
	{
		max = i;
		for (j = i+1; j < num_nodes; j++)
		{
			//if (nodes[j].levels > nodes[max].levels)
			if (nodes[nodeindices[j]].levels > nodes[nodeindices[max]].levels)
				max = j;
		}
		temp = nodeindices[i];
		nodeindices[i] = nodeindices[max];
		nodeindices[max] = temp;
	}
}



/* Here we assume that u is always valid */
void DFS_visit(GraphNode2 *nodes, GraphEdge *edges, int u, int &time, int inverse)
{
	int i;
	GraphEdge cur_e;
	int adj_node;

	nodes[u].visited = 1;
	time ++;

	//nodes[u].levels = time;

	for(i = 0; i < nodes[u].nedges; i++)
	{
		cur_e = edges[nodes[u].edges[i]];

		//adj_node = cur_e.node_index1;               //search the transpose graph
		adj_node = cur_e.node_index2;               //modified at 8/12/06, since we have reversed the flow

		if(adj_node == u)
			continue;

		if(nodes[adj_node].visited == 0)
		{
			nodes[adj_node].sscomp_index = num_sccomps;  //mark 'adj_node' as the same sub-tree (the same scc)
			                                             //as u
			nodes[adj_node].parent = u;                  //set the parent
			DFS_visit(nodes, edges, adj_node, time, inverse); //recursive to find out all the nodes in the subtree
		}
	}

	nodes[u].visited = 2;
	time ++;
	nodes[u].levels = time;
	nodes[u].sscomp_index = num_sccomps;
	
	////Add the node into the list instead of sorting it later
	if(cur_sccnode_index >= 0)
	{
		cur_nodes_order[cur_sccnode_index] = u;
		cur_sccnode_index --;
	}
}


/*
The non-recursive DFS
*/
/*
Here we need a stack to keep track of the searching
*/
int *dfs_stack = NULL;
int curMaxDFSStack;
int top_dfsstack;



void add_to_dfsstack(int elem)
{
	if(top_dfsstack == curMaxDFSStack)
	{
		/*extend it*/
		dfs_stack = (int*)realloc(dfs_stack, sizeof(int)*(curMaxDFSStack+500));

		if(dfs_stack == NULL)
		{
			time ( &g_rawtime );
			FILE *fp = fopen("mem_error.txt", "a");
			fprintf(fp, "fail to re_alloc memory for dfs_stack in routine DFS_visit\n");
			fprintf(fp, "current date and time are :  %s. \n", ctime (&g_rawtime) );
			fclose(fp);
			exit(-1);
		}

		curMaxDFSStack+=800;
	}

	top_dfsstack++;
	dfs_stack[top_dfsstack] = elem;
}


void DFS_visit(GraphNode2 *nodes, int num_nodes, GraphEdge *edges, int num_edges, 
			   int u, int scc_index)
{
	int i;
	GraphEdge cur_e;
	int cur_node, adj_node;
	int nchildren = 0;

	//nodes[u].visited = 1;

	/*initialize the stack*/
	if(dfs_stack != NULL)
		free(dfs_stack);

	curMaxDFSStack = num_nodes*2;
	dfs_stack = (int*)malloc(sizeof(int)*curMaxDFSStack);
	if(dfs_stack == NULL)
	{
		time ( &g_rawtime );
		FILE *fp = fopen("mem_error.txt", "a");
		fprintf(fp, "fail to allocate memory for dfs_stack in routine DFS_visit\n");
		fclose(fp);
		fprintf(fp, "current date and time are :  %s. \n", ctime (&g_rawtime) );
		exit(-1);
	}
	dfs_stack[0] = u;
	top_dfsstack = 0;

	//nodes[u].levels = time;
	while(top_dfsstack >= 0)/*keep searching until the stack is empty*/
	{
		/*pop up the last node*/
		cur_node = dfs_stack[top_dfsstack];

		if(nodes[cur_node].visited == 2 ) 
		{
			/*remove it from the stack*/
			top_dfsstack --;
			continue;
		}

		nodes[cur_node].visited = 1; /*mark it as 'visited' but not 'finished' yet*/

		nchildren = 0;
		for(i = nodes[cur_node].nedges-1; i >= 0; i--)
		//for(i = 0; i < nodes[cur_node].nedges; i++)
		{
			cur_e = edges[nodes[cur_node].edges[i]];

			adj_node = cur_e.node_index2;

			if(adj_node == cur_node)/*consider outgoing edges only*/
				continue;

			if(nodes[adj_node].visited == 0)
			{
				//nodes[adj_node].sscomp_index = scc_index;  //mark 'adj_node' as the same sub-tree (the same scc)
				nodes[adj_node].parent = cur_node;         //set the parent

				/*add to the stack*/
				add_to_dfsstack(adj_node);
				nchildren ++;
			}
		}

		if(nchildren == 0)/*this is a leaf or no more 'unvisited' children*/
		{
			nodes[cur_node].visited = 2;  /*the searching is 'finished' for this node*/
			nodes[cur_node].sscomp_index = scc_index;

			/*put it to the cur_nodes_order array*/
			if(cur_sccnode_index >= 0)
			{
				cur_nodes_order[cur_sccnode_index] = cur_node;
				cur_sccnode_index --;
			}


			/*remove it from the stack*/
			top_dfsstack --;
		}
	}

	free(dfs_stack);
	dfs_stack = NULL;
}


void DFS(GraphNode2 *nodes, int num_nodes, GraphEdge *edges, int inverse)
{
	int i;
	int time = 0;
	num_sccomps = 0;


	//initialization
	for(i = 0; i < num_nodes; i++)
	{
		nodes[i].visited = 0;
		//nodes[i].levels = 0;
		nodes[i].sscomp_index = -1;
		nodes[i].parent = i;
	}

	//
	for(i = 0; i < num_nodes; i++)
	{
		if(inverse == 0)
		{
			if(nodes[i].visited == 0)
			{
				DFS_visit(nodes, sccedges, i, time, inverse);
				//DFS_visit(nodes, num_nodes, sccedges, num_sccedges, i, num_sccomps);
			}
		}

		else
		{
			if(nodes[cur_nodes_order[i]].visited == 0)
			{
				DFS_visit(nodes, sccedges, cur_nodes_order[i], time, inverse);
				//DFS_visit(nodes, num_nodes, sccedges, num_sccedges, cur_nodes_order[i], num_sccomps);
				num_sccomps ++;           //for each subtree, increase the index
			}
		}
	}
}


void count_MorseSets()
{
	num_morsesets = 0;
	int i;

	for(i = 0; i < scclist.num_sccs; i++)
	{

		if(scclist.scccomponents[i].num_nodes <= 2 && scclist.scccomponents[i].num_singularities <= 0)
			continue;

		else if(scclist.scccomponents[i].num_nodes > 2 && scclist.scccomponents[i].valid == 0)
			continue;

		num_morsesets++;
	}

}

void BuildSCCElemList()
{
	int i;

	//allocate space for the SCC components
	if(scclist.scccomponents != NULL)
	{
		//we need to first free the node list for each scc
		for(i = 0; i < scclist.num_sccs; i++)
		{
			if(scclist.scccomponents[i].nodes != NULL)
			{
				free(scclist.scccomponents[i].nodes);
				scclist.scccomponents[i].nodes = NULL;
				scclist.scccomponents[i].num_nodes = 0;
			}

			//Release the list for singularities inside the SCC 07/23/06
			if(scclist.scccomponents[i].limitcycles != NULL)
				free(scclist.scccomponents[i].limitcycles);

			//Release the list for limit cycles inside the SCC  07/23/06
			if(scclist.scccomponents[i].singular_tri != NULL)
				free(scclist.scccomponents[i].singular_tri);
		}
		free(scclist.scccomponents);
	}

	scclist.scccomponents = (SCComponent*)malloc(sizeof(SCComponent) * num_sccomps);
	scclist.num_sccs = num_sccomps;

	int *numtriangles_eachscc = (int*)malloc(sizeof(int)*num_sccomps);

	for(i = 0; i < num_sccomps; i++)
		numtriangles_eachscc[i] = 0;

	//first, we search the whole mesh once, and get the number of triangles in each scc components
	for(i = 0; i < Object.nfaces; i++)
	{
		numtriangles_eachscc[sccnodes[i].sscomp_index] += 1;
	}

	//allocate space for each scc
	for(i = 0; i < num_sccomps; i++)
	{

		scclist.scccomponents[i].nodes = (int*) malloc(sizeof(int)*numtriangles_eachscc[i]);
		scclist.scccomponents[i].num_nodes = 0;

		////initialize
		scclist.scccomponents[i].singular_tri = NULL;
		scclist.scccomponents[i].num_singularities = 0;
		scclist.scccomponents[i].limitcycles = NULL;
		scclist.scccomponents[i].num_limitcycles = 0;
		scclist.scccomponents[i].num_attpts = 0;
		scclist.scccomponents[i].num_seppts = 0;
	}

	//record the scc components

	for(i = 0; i < Object.nfaces; i++)
	{
		scclist.scccomponents[sccnodes[i].sscomp_index].nodes[
			scclist.scccomponents[sccnodes[i].sscomp_index].num_nodes] = i;
		scclist.scccomponents[sccnodes[i].sscomp_index].num_nodes++;

		/* count the number of fixed points in the SCC 02/21/07 */
		if(Object.flist[i]->contain_singularity > 0
			||Object.flist[i]->singularity_index >= 0)
		{
			scclist.scccomponents[sccnodes[i].sscomp_index].singular_tri
				=Extend_link(scclist.scccomponents[sccnodes[i].sscomp_index].singular_tri, 
				scclist.scccomponents[sccnodes[i].sscomp_index].num_singularities);
			scclist.scccomponents[sccnodes[i].sscomp_index].singular_tri[
				scclist.scccomponents[sccnodes[i].sscomp_index].num_singularities] = i;
			scclist.scccomponents[sccnodes[i].sscomp_index].num_singularities++;
		}
	}

	//Mark the boundary edges of each SCC component 07/18/06
	Edge *cur_e;
	Face *face;
	int j, k;
	int other_t;
	for(i = 0; i < scclist.num_sccs; i++)
	{
		for(j = 0; j < scclist.scccomponents[i].num_nodes; j++)
		{
			face = Object.flist[scclist.scccomponents[i].nodes[j]];

			for(k = 0; k < 3; k++)
			{
				cur_e = face->edges[k];

				other_t = cur_e->tris[0];

				if(other_t == face->index)
					other_t = cur_e->tris[1];

				//Test whether it is a boundary edge
				if(IsRepeated(scclist.scccomponents[i].nodes, other_t, scclist.scccomponents[i].num_nodes))
				{
					cur_e->OnBoundary = 0;
					continue;  //not a boundary edge
				}

				cur_e->OnBoundary = 1; 
			}
		}
	}

	free(numtriangles_eachscc);
}


extern void write_down_SCC(); //test code 02/04/07


/*
The main routine for finding the SCC
*/

void FindSCC()
{
	///
	cur_sccnode_index = num_sccnodes - 1;
	FILE *fp;

	DFS(sccnodes, num_sccnodes, sccedges, 0);
	//SortNodes(sccnodes, cur_nodes_order, num_sccnodes);
	ReverseEdges();
		
	//fp = fopen("detectprocess.txt", "a");
	//	fprintf(fp, "start second DFS\n");
	//	fclose(fp);	
	
	DFS(sccnodes, num_sccnodes, sccedges, 1);
	
	//fp = fopen("detectprocess.txt", "a");
	//	fprintf(fp, "start building SCC list\n");
	//	fclose(fp);

	BuildSCCElemList();
		
	ReverseEdges();  /*reverse the edges back to its original states */

		//fp = fopen("detectprocess.txt", "a");
		//fprintf(fp, "%d SCC's have been found. \n", num_sccomps);
		////fprintf(fp, "current date and time are :  %s. \n", ctime (&rawtime) );
		////fprintf(fp, "finding SCC spent %f s. \n", timespent);
		////fprintf(fp, "%d nodes and %d edges in the graph.\n", num_sccnodes, num_sccedges);
		//fclose(fp);	

		//write_down_SCC(); //test code 02/04/07
}

void AllocForSCC()
{
	sccedges = (GraphEdge*) malloc(sizeof(GraphEdge)* Object.nedges*2); //it requires the construction of the edge list
	num_sccedges = 0;
	curMaxNumDirGraphEdges = Object.nedges*2;

	sccnodes = (GraphNode2*) malloc(sizeof(GraphNode2) * Object.nfaces);
	cur_nodes_order = (int *) malloc(sizeof(int) * Object.nfaces);
	num_sccnodes = 0;

	int i;

	for(i = 0; i < Object.nfaces; i++)
	{
		sccnodes[i].edges = NULL;
		sccnodes[i].nedges = 0;
	}
}

void FinalizeSCC()
{
	int i;

	for(i = 0; i < Object.nfaces; i++)
	{
		if(sccnodes[i].edges!=NULL)
		{
			free(sccnodes[i].edges);
			sccnodes[i].edges = NULL;
		}
	}

	free(sccnodes);
	free(sccedges);
	free(cur_nodes_order);
}

void finalize_scclist()
{
	int i;

	for(i = 0; i < num_sccomps; i++)
	{
		if(scclist.scccomponents[i].nodes != NULL)
		{
			free(scclist.scccomponents[i].nodes);
			scclist.scccomponents[i].nodes = NULL;
		}
		scclist.scccomponents[i].num_nodes = 0;
	}

	free(scclist.scccomponents);
}

void InitForSCC()
{
	int i;

	//initialize the nodes
	for(i = 0; i < Object.nfaces; i++)
	{
		if(	sccnodes[i].edges != NULL)
		{
			free(sccnodes[i].edges);
			sccnodes[i].edges = NULL;
		}

		sccnodes[i].nedges = 0;
		sccnodes[i].visited = 0;
		sccnodes[i].parent = i;
		sccnodes[i].levels = 0;

		//Object.flist[i]->contain_singularity = 0;
	}
	
	//initialize the edges
	for(i = 0; i < curMaxNumDirGraphEdges/*2*Object.nedges*/; i++)
	{
		sccedges[i].visited = 0;
		sccedges[i].node_index1 = sccedges[i].node_index2 = -1;
	}

	//
	scclist.scccomponents = NULL;
	scclist.num_sccs = 0;

	num_sccnodes = 0;
	num_sccedges = 0;

}

/*
Mark all the valid SCC
*/
extern bool IsValidSCC(int);
void mark_All_Valid_SCCS()
{
	for(int i = 0; i < scclist.num_sccs; i++)
	{
		if(!IsValidSCC(i))
		{
			scclist.scccomponents[i].valid = 0;
			continue;
		}

		scclist.scccomponents[i].valid = 1;   //this is a valid SCC
	}
}



/**************************************************************************/


/* To visualize the path, we have to save the nodes along the path! 02/25/07
*/


bool path_From_T1_To_T2(int n1, int n2)
{
	//the maximum searching number is the number of the edges in the graph
	int count = 0;
	int cur_node = n1;
	int cur_node_id = 0;

	int i;
	//reset the 'visited' flags of all the edges
	for(i = 0; i < num_sccedges; i++)
	{
		sccedges[i].visited = 0;
	}
	
	/* allocate the searching stack for DFS searching */
	if(searchstack != NULL)
		free(searchstack);
	searchstack = (int *)malloc(sizeof(int) * (scclist.scccomponents[sccnodes[n1].sscomp_index].num_nodes));
	MaxStackElems = scclist.scccomponents[sccnodes[n1].sscomp_index].num_nodes+100;
	nelem_in_stack = 0;

	searchstack[0] = n1;
	nelem_in_stack = 1;

	if(path != NULL)
		free(path);
	num_tris_inpath = 0;
	sccnodes[n1].parent = n1;

	/*Here we use BFS*/
	for(int j = 0; j < nelem_in_stack; j++)
	{
		cur_node = searchstack[j];

		if(cur_node == n2)
			goto LL;

		//pick one unvisited edges
		for(i = 0; i < sccnodes[cur_node].nedges; i++)
		{
			if(sccedges[sccnodes[cur_node].edges[i]].visited == 1)
				continue;

			/*we need to make sure that they are in the same SCC!*/
			if(sccnodes[sccedges[sccnodes[cur_node].edges[i]].node_index2].sscomp_index!=
				sccnodes[cur_node].sscomp_index)
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

			if(IsRepeated(searchstack,sccedges[sccnodes[cur_node].edges[i]].node_index2,nelem_in_stack))
				continue;

			searchstack[nelem_in_stack] = sccedges[sccnodes[cur_node].edges[i]].node_index2;
			sccnodes[sccedges[sccnodes[cur_node].edges[i]].node_index2].parent = cur_node;
			nelem_in_stack++;

			sccedges[sccnodes[cur_node].edges[i]].visited = 1;
		}
	}
			free(searchstack);
			searchstack = NULL;
	return false;

	/*The following routine use DFS*/
	//while(cur_node != n2 && count < num_sccedges)
	//{
	//	/* pick the last element of the stack */
	//	cur_node = searchstack[nelem_in_stack-1];

	//	nelem_in_stack --;

	//	if(cur_node == n2)
	//		goto LL;

	//	/***********************************************************/

	//	//pick one unvisited edges
	//	for(i = 0; i < sccnodes[cur_node].nedges; i++)
	//	{
	//		if(sccedges[sccnodes[cur_node].edges[i]].visited == 1)
	//			continue;

	//		/*we need to make sure that they are in the same SCC!*/
	//		if(sccnodes[sccedges[sccnodes[cur_node].edges[i]].node_index2].sscomp_index!=
	//			sccnodes[cur_node].sscomp_index)
	//			continue;

	//		//set the node_index2 of the edge as the next cur_node
	//		//cur_node = sccedges[sccnodes[cur_node_id].edges[i]].node_index2;
	//		if(nelem_in_stack >= MaxStackElems) /*extend it*/
	//		{
	//			searchstack = (int *)realloc(searchstack, sizeof(int)*(MaxStackElems+50));

	//			if(searchstack == NULL)
	//				exit(-1);

	//			MaxStackElems += 50;
	//		}

	//		if(IsRepeated(searchstack,sccedges[sccnodes[cur_node].edges[i]].node_index2,nelem_in_stack))
	//			continue;

	//		searchstack[nelem_in_stack] = sccedges[sccnodes[cur_node].edges[i]].node_index2;
	//		sccnodes[sccedges[sccnodes[cur_node].edges[i]].node_index2].parent = cur_node;
	//		nelem_in_stack++;

	//		sccedges[sccnodes[cur_node].edges[i]].visited = 1;
	//		count++;
	//		//break;
	//	}

	//	if(nelem_in_stack == 0)
	//	{
	//		free(searchstack);
	//		searchstack = NULL;
	//		return false;
	//	}
	//}

	if(count < num_sccedges)
	{
LL:		FILE *fp = fopen("pathtest.txt","w");
		fprintf(fp, "Yes, find a path between %d and %d\n", n1, n2);
		fclose(fp);

		/*backtrack from the end triangle to the starting triangle*/
		int num_tris_in_path = 1;
		cur_node = n2;
		while(cur_node != n1)
		{
			cur_node = sccnodes[cur_node].parent;

			if(cur_node < 0)
			{
				free(searchstack);
			    searchstack = NULL;
				return false;
			}
			num_tris_in_path ++;
		}

		fp = fopen("pathtest.txt","a");
		fprintf(fp, "there are %d triangles in the path\n", num_tris_in_path);
		fclose(fp);

		/*allocate memeory for the path*/
		path = (int*)malloc(sizeof(int)*num_tris_in_path);
		cur_node = n2;
		path[0] = cur_node;
		num_tris_inpath = 1;
		fp = fopen("pathtest.txt","a");
		fprintf(fp, "they are: \n");
		while(cur_node != n1)
		{
			cur_node = sccnodes[cur_node].parent;
			path[num_tris_inpath] = cur_node;
			num_tris_inpath ++;
			fprintf(fp, "%d \n", cur_node);
		}
		fclose(fp);

			free(searchstack);
			searchstack = NULL;
		return true;
	}
}



