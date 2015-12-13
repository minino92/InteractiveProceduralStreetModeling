
/*
SCCVerify.cpp
*/
#include "stdafx.h"

#include "FindSCC.h"
#include "VFDataStructure.h"

extern GraphEdge *sccedges;
extern int num_sccedges;
extern int curMaxNumDirGraphEdges;
extern GraphNode2 *sccnodes;
extern int *cur_nodes_order;
extern int num_sccnodes;


extern int *searchstack;
extern int nelem_in_stack;
extern int MaxStackElems;



bool is_Node_in_Graph(int node_index)
{
	int i;
	for(i = 0; i < num_sccnodes; i++)
	{
		if(sccnodes[i].node_index == node_index)
			return true;
	}

	return false;
}


bool remove_Edge_from_Node(int node_index, int edge_index)
{
	int i;
	int node_id;

	for(i = 0; i < num_sccnodes; i++)
	{
		if(sccnodes[i].node_index == node_index)
		{
			node_id = i;
			break;
		}
	}

	//
	if(node_id >= num_sccnodes)
		return false;

	for(i = 0; i < sccnodes[node_id].nedges; i++)
	{
		if(sccnodes[node_id].edges[i] == edge_index)
		{
			//remove it
			if(i == sccnodes[node_id].nedges-1)
			{
				sccnodes[node_id].nedges--;
				return true;
			}

			for(int j = i; j < sccnodes[node_id].nedges-1; j++)
			{
				sccnodes[node_id].edges[j] = sccnodes[node_id].edges[j+1];
			}
			sccnodes[node_id].nedges--;
			return true;
		}
	}

	return false;
}


void load_Graph(char *filename)
{
	//load from file
	FILE *fp = fopen(filename, "r");

	if(fp == NULL)
	{
		printf(" Can't open file %s \n", filename);
		exit(-1);
	}

	int i, j;
	int nnodes, nedges;

	int node_id, countnodeedges = 0;
	
	//get the number of nodes and edges

	fscanf(fp, "#nodes:%d\n", &nnodes);
	fscanf(fp, "#edges:%d\n", &nedges);

	num_sccnodes = nnodes;
	num_sccedges = nedges;

	//allocate the memeory for node and edge lists, respectively
	//here, we consider only steady graph!
	sccnodes = (GraphNode2 *)malloc(sizeof(GraphNode2) * (nnodes+1));
	
	/* allocate the searching stack for DFS searching */
	searchstack = (int *)malloc(sizeof(int) * (nnodes+1));
	MaxStackElems = nnodes;

	//initialize the elements in the sccnodes list
	for(i = 0; i < nnodes; i++)
	{
		sccnodes[i].nedges = 0;
		sccnodes[i].edges = NULL;
		sccnodes[i].node_index = -1;
	}

	sccedges = (GraphEdge *)malloc(sizeof(GraphEdge )*(nedges +1));
	//initialize the elements in the sccedges list
	for(i = 0; i < nedges; i++)
	{
		sccedges[i].edge_index = -1;
	}


	//read nodes and edges
	int cur_sccedge_index = 0;
	for(i = 0; i < nnodes; i++)
	{
		fscanf(fp, "the node: %d has %d edges\n", &node_id, &countnodeedges);

		sccnodes[i].node_index = node_id;

		//allocate memeory for ith element of the sccnodes list
		sccnodes[i].edges = (int *)malloc(sizeof(int)*countnodeedges);
		sccnodes[i].nedges = countnodeedges;

		/* Note that some edges link to other nodes are not in the graph !! */
		for(j = 0; j < countnodeedges; j++)
		{
			int node2;

			fscanf(fp, "%d \n", &node2);

			SCC_AddToEdge(node_id, node2, cur_sccedge_index); //add a new edge
			sccnodes[i].edges[j] = cur_sccedge_index-1;  //add to node edge list
		}
	}

	num_sccedges = cur_sccedge_index;
	//remove those edges who point to the nodes that are not in the graph
	for(i = 0; i < num_sccedges; i++)
	{
		if(is_Node_in_Graph(sccedges[i].node_index2))
			continue;

		//remove the edge from the node "node_index1"
		remove_Edge_from_Node(sccedges[i].node_index1, i);
		
	}
}


/*
-------------------Recursive version---------
dfs(v)
    process(v)
    mark v as visited
    for all vertices i adjacent to v not visited
        dfs(i)

-------------------Another version-----------------------

dfs(graph G)
{
  list L = empty
  tree T = empty
  choose a starting vertex x
  search(x)
  while(L is not empty)
    remove edge (v, w) from end of L
    if w not yet visited
    {
      add (v, w) to T
      search(w)
    }
}

   
search(vertex v)
{
  visit v
  for each edge (v, w)
    add edge (v, w) to end of L
}


*/

/*
BFS algorithm from wikipedia 
*/

//struct Vertex {
//        ...
//        std::vector<int> out;
//        ...
//};
//
//std::vector<Vertex> graph(vertices);
//
//bool BFS(std::vector<Vertex>& graph, int start, int end) {
//  std::queue<int> next;
//  std::vector<int> parent(graph.size(), 0); // 0 means filled with zeros.
//  parent[start] = -1;
//  next.push(start);
//  while (!next.empty()) {
//    int u = next.front();
//    next.pop();
//    // Here is the point where you can examine the u th vertex of graph
//    // For example:
//    if (u == end) return true;
//    for (size_t j = 0; j < graph[u].out.size(); ++j) {
//      // Look through neighbors.
//      int v = graph[u].out[j];
//      if (parent[v] == 0) {
//        // If v is unvisited.
//        parent[v] = u;
//        next.push(v);
//      }
//    }
//  }
//  return false;
//}
//

/* To visualize the path, we have to save the nodes along the path! 02/25/07
*/

bool path_From_n1_To_n2(int n1, int n2)
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

	searchstack[0] = n1;
	nelem_in_stack = 1;

	while(cur_node != n2 && count < num_sccedges)
	{
		/* pick the last element of the stack */
		cur_node = searchstack[nelem_in_stack-1];

		nelem_in_stack --;

		if(cur_node == n2)
			return true;

		/***********************************************************/

		//pick one unvisited edges
		for(i = 0; i < sccnodes[cur_node].nedges; i++)
		{
			if(sccedges[sccnodes[cur_node].edges[i]].visited == 1)
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

			searchstack[nelem_in_stack] = sccedges[sccnodes[cur_node].edges[i]].node_index2;
			nelem_in_stack++;

			sccedges[sccnodes[cur_node].edges[i]].visited = 1;
			count++;
			//break;
		}

		if(nelem_in_stack == 0)
			return false;

	}

	if(count < num_sccedges)
	{
		return true;
	}
	return false;
}

bool twopaths_between_n1_n2(int n1, int n2)
{
	//search path from n1 to n2
	if(!path_From_n1_To_n2(n1, n2))
		return false;
	if(!path_From_n1_To_n2(n2, n1))
		return false;
	return true;
}


bool verify_SCC()
{
	//start detecting the path between two any nodes in the graph

	//if any given two nodes n1 and n2 in the graph, we can find a path from n1 to n2
	//and from n2 to n1, then we say the whole graph is one SCC
	int i, j;
	for(i = 0; i < num_sccnodes-1; i++)
	{
		for(j = i+1; j < num_sccnodes; j++)
		{
			if(!twopaths_between_n1_n2(sccnodes[i].node_index, sccnodes[j].node_index))
				return false;
		}
	}

	return true;
}


//void main()
//{
//	load_Graph("thescc.txt");
//
//	if(verify_SCC())
//		printf("great, the whole graph is one SCC !");
//	else
//		printf("sorry, the whole graph is not a SCC !");
//
//	getchar();
//}