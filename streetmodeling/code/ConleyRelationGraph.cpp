////ConleyRelationGraph.cpp

////This module implements the detection and creation of the Conley relation graph
////(i.e. a graph showing the connections of repeller, attractors and saddles)

#include "stdafx.h"

#include "ConleyRelationGraph.h"

#include "VFDataStructure.h"

#include "VFAnalysis.h"

#include "LocalTracing.h"

#include "LimitCycleDetect.h"


/*----------------------------------------------------*/
extern GraphNode *graphnodes;
extern GraphEdge *graphedges;
extern int cur_node_index;
extern int cur_graphedge_index;

extern LimitCycle *limitcycles;
extern int cur_limitcycle_index;
extern Singularities *singularities;
extern int cur_singularity_index;
extern Separatrices *separatrices;             //array for group of separatrices
extern int cur_separatrices_index;

extern LineSeg **trajectories;                 //trajectories' list
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory

extern Polygon3D Object;

/*----------------------------------------------------*/
int MaxNumGraphEdges = 50;

////Set the default window size for the Conley graph displaying (0.0f,2.0f,0.0f,4.0f)
double ConleyGraphWin_x = 4;
double ConleyGraphWin_y = 2;
double types_interval_y = ConleyGraphWin_y / 4;

/*----------------------------------------------------*/
extern int maxlength_searcharray ;
extern int *NodeSearchArray ;
extern int cur_endposition;
extern int cur_search_node_pos;
extern int *MediaNodes;
extern int Num_MediaNodes;
extern int MaxMediaNodes;

/*----------------------------------------------------*/
extern int *repeller_nodes;     //the indices of the selected repellers in the conley graph
extern int NumRepellers;        //the number of the being selected repellers
extern int *attractor_nodes;    //the indices of the selected attractors in the conley graph
extern int NumAttractors;       //the number of the being selected attractors

/*----------------------------------------------------*/
extern bool IsRepeated(int *a, int b, int num);


////a new variable to store the previous detected singularities 3/19/06
PreSingularityList *presingularityList = NULL;  //use a quick way now, need to extend to dynamic structure later
int num_presingularities = 0;
int r_counter, a_counter, s_counter; //for labeling 3/14/06

PreLimitCycle *prelimitcyclelist = NULL;
int num_prelimitcycles = 0;

////Urgent variables for video capture 3/29/06
int *Pre_CancelledNode = (int *)malloc(sizeof(int) * 20);
int num_cancellednode = 0;

/*----------------------------------------------------*/
////routine to judge whether there is at least one interval(saddle)
////between two nodes (usually a repeller and an attractor)
////we always suppose the inputs are valid 2/15/06
bool HasInterval(int node1, int node2)
{
	int i, j;
	////we need to make sure node1 is the repeller
	if(graphnodes[node1].type == 1)
	{
		int temp = node1;
		node1 = node2;
		node2 = temp;
	}


	////Begin to judge whether there is at least one interval connecting them
	int *temp_connect_saddlelist = new int[graphnodes[node1].nedges];
	int num_connect_saddlelist = 0;
	int cur_edge, n2;

	////Begin searching from this repeller to search all the connected saddles
	for(j = 0; j < graphnodes[node1].nedges; j++)
	{
		cur_edge = graphnodes[node1].edges[j];
		n2 = graphedges[cur_edge].node_index2;

		if(graphnodes[n2].type == 2) //if this is a saddle
		{
			////we need to perform repeating test if we allow multiple connections between two nodes!!
			temp_connect_saddlelist[num_connect_saddlelist] = n2;
			num_connect_saddlelist++;
		}
	}

	for(j = 0; j < num_connect_saddlelist; j++)
	{
		for(i = 0; i < graphnodes[temp_connect_saddlelist[j]].nedges; i++)
		{
			cur_edge = graphnodes[temp_connect_saddlelist[j]].edges[i];
			n2 = graphedges[cur_edge].node_index2;

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


/* An easy way for this is to set extra member for separatrix structure
to store the length of the separatrix
*/

/*
This routine return the length of the specific separatrix
sep_id: 1, 2, 3 or 4
*/
double GetLength_separatrix(int saddle, int sep_id)
{
	int sep = singularities[saddle].separtices;

	switch(sep_id)
	{
	case 1:	 return (separatrices[sep].length1);

	case 2:  return (separatrices[sep].length2);

	case 3:  return (separatrices[sep].length3);

	case 4:  return (separatrices[sep].length4);
	}
}

//double GetLength_saddletocycle(int saddle, int cycleID, int sep_id)
//{
//}

/* 
For the following non-separatrix connection, we can calculate their length during finding the 
connection between them
*/

double GetLength_cycletocycle(int cycle1, int cycle2)
{
	return 0.;
}

double GetLength_cycletosing(int cycleID, int singID)
{
	return 0.;
}


/*----------------------------------------------------*/
////Routine for building the Conley partial graph
void BuildGraph()
{
	int i, j;
	int trajectory_id;
	int lineseg_id;
	int last_triangle;

	int node1, node2;
	
	/*Important initialization part*/
	////It seems that it is not a good way to put the following codes here!!!! 10/14/05
	////perform singularities detection
	if(cur_singularity_index == 0)
		CaptureSing();

	////separatrices calculation
	if(cur_separatrices_index == 0)
		CalSeparatrices();

	////Limit cycle detection, not consider limit cycle now 10/15/05
	//if(cur_limitcycle_index == 0)
	//LimitCycleDetect();

	/*Build the graph*/

	//1. Build the nodes
	////Release previous space
	if(graphnodes != NULL)
		free(graphnodes);

	if(graphedges != NULL)
		free(graphedges);

	////Allocate new space for nodes according to the number of singularities and number of limit cycles 
	graphnodes = (GraphNode *)malloc(sizeof(GraphNode) * (cur_singularity_index + cur_limitcycle_index + 10));
	cur_node_index = 0;

	////Allocate a graph edge array with default length equal to MaxNumGraphEdges
	graphedges = (GraphEdge *)malloc(sizeof(GraphEdge) * MaxNumGraphEdges);
	cur_graphedge_index = 0;

	////Assign one node for each singularity and limit cycle
	r_counter = a_counter = s_counter = 0;

	for(i = 0; i < cur_singularity_index; i++)
	{
		/* we remove following condition at 02/10/2007 */
		//if(singularities[i].type != CWCENTER && singularities[i].type != CCWCENTER)
		//{
			graphnodes[cur_node_index].node_index = cur_node_index;
			graphnodes[cur_node_index].singularityID = i;
			graphnodes[cur_node_index].LimitCycleID = -1;
			graphnodes[cur_node_index].cancelled = 0;
			singularities[i].node_index = cur_node_index;
			
			if(singularities[i].type == SOURCE || singularities[i].type == RFOCUS)
			{
				graphnodes[cur_node_index].type = 0;   ////it is a repeller
				graphnodes[cur_node_index].labelindex = r_counter;
				r_counter++;
			}
			else if(singularities[i].type == SINK || singularities[i].type == AFOCUS)
			{
				graphnodes[cur_node_index].type = 1;   ////it is a attractor
				graphnodes[cur_node_index].labelindex = a_counter;
				a_counter++;
			}
			else
			{
				graphnodes[cur_node_index].type = 2;   ////it is a saddle
				graphnodes[cur_node_index].labelindex = s_counter;
				s_counter++;
			}

			graphnodes[cur_node_index].nedges = 0;
			graphnodes[cur_node_index].edges = NULL;

			cur_node_index ++;
		//}
	}

	////limit cycle
    for(i = 0; i < cur_limitcycle_index; i++)
	{
		graphnodes[cur_node_index].node_index = cur_node_index;
		graphnodes[cur_node_index].LimitCycleID = i;
		graphnodes[cur_node_index].singularityID = -1;
		graphnodes[cur_node_index].cancelled = 0;
		limitcycles[i].node_index = cur_node_index;
		if(limitcycles[i].type == 0)  ////it is a repeller
		{
			graphnodes[cur_node_index].type = 0;
			graphnodes[cur_node_index].labelindex = r_counter;
			r_counter++;
		}
		else                          ////it is an attractor
		{
			graphnodes[cur_node_index].type = 1;
			graphnodes[cur_node_index].labelindex = a_counter;
			a_counter++;
		}

		graphnodes[cur_node_index].nedges = 0;
		graphnodes[cur_node_index].edges = NULL;

		cur_node_index ++;
	}

	//2. Build the edges

	////The graph is a directed graph!!!! 11/19/05, the direction of an edge is node1->node2
	////I.E. node1 should be a repeller and node2 should be an attractor
	////Note that, during the building of each edge, we need to check whether it is a previous edge!

	////(a) Find the edges between saddles and other singularities or limit cycles
	for(i = 0; i < cur_singularity_index; i++)
	{

		////Begin from each saddle, check its separatrices
		////if the last triangle of the separatrix contains singularity, build an edge between the singularity and current saddle
		if(singularities[i].type == SADDLE)
		{
		////if the last triangle of the separatrix does not contain any singularity
		////it maybe forms a closed streamline (we do not consider the 'center' here now) or reach the boundary
		////Build an edge between the limit cycle and current saddle

			////Then, follow the separatrices of the saddle 
			for(j = 0; j < 4; j++)
			{
				switch(j)
				{
				case 0:
					trajectory_id = separatrices[singularities[i].separtices].sep1;
					break;
				case 1:
					trajectory_id = separatrices[singularities[i].separtices].sep2;
					break;
				case 2:
					trajectory_id = separatrices[singularities[i].separtices].sep3;
					break;
				case 3:
					trajectory_id = separatrices[singularities[i].separtices].sep4;
					break;
				}

				lineseg_id = num_linesegs_curtraj[trajectory_id];
				if(lineseg_id <= 0) continue;
				last_triangle = trajectories[trajectory_id][lineseg_id-1].Triangle_ID;
				if(last_triangle >= 0 && Object.flist[last_triangle]->contain_singularity > 0) ////contains singularity
				{
					////Get the node index for the other singularity
					node1 = singularities[i].node_index;
					node2 = singularities[Object.flist[last_triangle]->singularity_index].node_index;

					if(node2 < 0)
						continue;

					if(singularities[Object.flist[last_triangle]->singularity_index].type == SOURCE
						||singularities[Object.flist[last_triangle]->singularity_index].type == RFOCUS)
					{
						////here, saddle is the attractor
						node2 = node1;
						node1 = singularities[Object.flist[last_triangle]->singularity_index].node_index;
					}

					////Add to graph edge array
					AddToEdgeArray(node1, node2, cur_graphedge_index);

					////Save the edge length (the length of the separatrix) to the new added edge
					////08/10/06, here i happens to be the index of the saddle, j is the index of the sep of the saddle
					graphedges[cur_graphedge_index-1].flow_length = GetLength_separatrix(i, j+1);
				   
					////Need to add to the edge list of node1 and node2
					AddEdgeToNode(node1, cur_graphedge_index-1);
					AddEdgeToNode(node2, cur_graphedge_index-1);
				}
			}

			////(b) build the connections between saddles and limit cycles 2/15/06
			for(j = 0; j < singularities[i].num_connected_limitcycles; j++)
			{
				node1 = singularities[i].node_index;
				node2 = limitcycles[singularities[i].connected_limitcycles[j]].node_index;

				if(limitcycles[singularities[i].connected_limitcycles[j]].type == 0)
				{
					node2 = node1;
					node1 = limitcycles[singularities[i].connected_limitcycles[j]].node_index;
				}

				if(node2 < 0 || node1 < 0)
					continue;
				

				////If there is no connection between them, connect them
				if(!IsConnected(node1, node2))
				{
					////Add to graph edge array
					AddToEdgeArray(node1, node2, cur_graphedge_index);

					////We need to first find out the location of the singularity in the list 08/10/06
					int k;
					for(k = 0; k < limitcycles[singularities[i].connected_limitcycles[j]].num_connectedsaddles;
						k++)
					{
						if(limitcycles[singularities[i].connected_limitcycles[j]].connected_saddle[k]==i)
							break;
					}
					
					graphedges[cur_graphedge_index-1].flow_length = 
						limitcycles[singularities[i].connected_limitcycles[j]].flow_length_connectedsing[k];
					
					////Need to add to the edge list of node1 and node2
					AddEdgeToNode(node1, cur_graphedge_index-1);
					AddEdgeToNode(node2, cur_graphedge_index-1);
				}
			}
		}
	}

	////(c) Find the edges between limit cycles and other singularities
	////Here we build the saddle-limit cycle connections first
	////which can be fulfilled in previous step 2/15/06

	////Here we begin the searching from the connected list in limit cycle data structure
	////Note that we also need to consider the directions of the connections
	for(i = 0; i < cur_limitcycle_index; i++)
	{
		node1 = limitcycles[i].node_index;

		for(j = 0; j < limitcycles[i].num_connectedsaddles; j++)
		{
			node2 = singularities[limitcycles[i].connected_saddle[j]].node_index;
			
			if(limitcycles[i].type == 1)
			{
				////this is an attractive limit cycle, hence saddle becomes a repeller here
				node1 = node2;
				node2 = limitcycles[i].node_index;
			}

			if(node1<0 || node2<0)
				continue;

			if(HasInterval(node1, node2))  ////there is saddle connecting them, very important here 07/26/06
				continue;

			if(IsConnected(node1, node2)) ////avoid repeating  08/10/06
				continue;


			////Add to graph edge array
			if(AddToEdgeArray(node1, node2, cur_graphedge_index))
			{
				////Need to add to the edge list of node1 and node2
				AddEdgeToNode(node1, cur_graphedge_index-1);
				AddEdgeToNode(node2, cur_graphedge_index-1);

				////we need to save the edge length here 08/10/06
				////We need to first find out the location of the singularity in the list 08/10/06
				graphedges[cur_graphedge_index-1].flow_length = 
					limitcycles[i].flow_length_connectedsing[j];
			}
		}
	}

	////(d) Find the edges between embeded limit cycles
	for(i = 0; i < cur_limitcycle_index; i++)
	{
		node1 = limitcycles[i].node_index;

		for(j = 0; j < limitcycles[i].num_connectedcycles; j++)
		{
			////if there has already one edge between them, continue
			node2 = limitcycles[limitcycles[i].connected_limitcycle[j]].node_index;

			if(limitcycles[limitcycles[i].connected_limitcycle[j]].type == 
				limitcycles[i].type)
				continue;

			if(limitcycles[i].type == 1)
			{
				////this is an attractive limit cycle, hence the other becomes a repeller here
				node1 = node2;
				node2 = limitcycles[i].node_index;
		    }

			if(node1 < 0 || node2 < 0)
				continue;

			if(HasInterval(node1, node2))  ////there is saddel connecting them
				continue;

			if(IsConnected(node1, node2)) ////avoid repeating  08/10/06
				continue;

			if(!AddToEdgeArray(node1, node2, cur_graphedge_index))
				continue;

			else
			{
				AddEdgeToNode(node1, cur_graphedge_index-1);
				AddEdgeToNode(node2, cur_graphedge_index-1);

				////store the edge length 08/10/06
				graphedges[cur_graphedge_index-1].flow_length = 
					limitcycles[i].flow_length_connectedcycle[j];
			}
		}
	}


	////Make a test here 08/10/06, write the edge length into a file
	FILE *fp = fopen("edgelength.txt", "w");
	for(i = 0; i < cur_graphedge_index; i++)
	{
		fprintf(fp, "edge %d: %f\n", i, graphedges[i].flow_length);
	}
	fclose(fp);

}


/*-----------------------------------------------------------
* Create a new edge, and add node1 and node2 as its end points
-----------------------------------------------------------*/
bool AddToEdgeArray(int node1, int node2, int &cur_index)
{
	////Check whether there has been an edge between them
	if(!IsRepeatedEdge(node1, node2, cur_index))
	{
		////if there is no edge between them, add an edge
		if(cur_index >= MaxNumGraphEdges-1)
		{
			MaxNumGraphEdges += 50;
			graphedges = (GraphEdge*)realloc(graphedges, sizeof(GraphEdge)*MaxNumGraphEdges);
		}

		graphedges[cur_index].edge_index = cur_index;
		graphedges[cur_index].node_index1 = node1;
		graphedges[cur_index].node_index2 = node2;

		graphedges[cur_index].visited = 0;

		cur_index++;

		return true;
	}

	else
		return false;
}


/*-----------------------------------------------------------
* Add the edge "edgeindex" to the edge list of "node"
-----------------------------------------------------------*/
void AddEdgeToNode(int node, int edgeindex)
{
	if(node < 0)
		return;

	int i;
	int *temp_list = graphnodes[node].edges;

	if(temp_list == NULL)
	{
		graphnodes[node].edges = (int*)malloc(sizeof(int));
		graphnodes[node].edges[0] = edgeindex;
		graphnodes[node].nedges = 1;
	}

	else{
		graphnodes[node].edges = (int*)malloc(sizeof(int)*(graphnodes[node].nedges+1));
		for(i = 0; i < graphnodes[node].nedges; i++)
		{
			graphnodes[node].edges[i] = temp_list[i];
		}
		graphnodes[node].edges[i] = edgeindex;

		graphnodes[node].nedges += 1;

		free(temp_list);
	}
}



bool IsRepeatedEdge(int node1, int node2, int cur_index)
{
	int i, j;
	for(i = 0; i < cur_index; i++)
	{
	}
	return false;
}


////Judge whether there has already one edge between them
bool IsConnected(int node1, int node2)
{
	int i, cur_edge;

	for(i = 0; i < graphnodes[node1].nedges; i++)
	{
		cur_edge = graphnodes[node1].edges[i];
		
		if(graphedges[cur_edge].node_index1 == node1)
		{
			if(graphedges[cur_edge].node_index2 == node2)
				return true;
		}

		else{
			if(graphedges[cur_edge].node_index1 == node2)
				return true;
		}
	}

	return false;
}


////Calculate the positions of the graph nodes
////It will be called when we update the graph
////Modified at 1/31/06

void LayOutNodes()
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
	for(i = 0; i < cur_node_index; i++)
	{
		if(graphnodes[i].type == 0)
			num_repellers ++;
		else if(graphnodes[i].type == 1)
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

	////2. Then, evenly space the positions of these nodes
	for(i = 0; i < cur_node_index; i++)
	{
		if(graphnodes[i].type == 0)
		{
			graphnodes[i].pos_x = (cur_repeller_index+1)*repeller_interval;
			graphnodes[i].pos_y = repeller_y;
			cur_repeller_index ++;
		}
		else if(graphnodes[i].type == 1)
		{
			graphnodes[i].pos_x = (cur_attractor_index+1)*attractor_interval;
			graphnodes[i].pos_y = attractor_y;
			cur_attractor_index ++;
		}
		else{
			graphnodes[i].pos_x = (cur_saddle_index+1)*saddle_interval;
			graphnodes[i].pos_y = saddle_y;
			cur_saddle_index ++;
		}
	}
}


/*-----------------------------------------------------------
* Finding the connected nodes on the paths that connect
* the input two nodes in the graph (only for singularities)
-----------------------------------------------------------*/
bool SearchConnectComponents(int sourcenode, int destnode, int *media_nodes, int &num_medianodes)
{
	int i, j;
	GraphNode *cur_node;
	int nodewanttoadd = -1;
	int nodewanttocompare = -1;

	////Initial part
	maxlength_searcharray = 50;
	if(NodeSearchArray != NULL)
		free(NodeSearchArray);
    NodeSearchArray = (int*)malloc(sizeof(int) * maxlength_searcharray); 
	cur_endposition = 0;
    cur_search_node_pos = 0;

	////Reset the visited flag
	for(i = 0; i < cur_graphedge_index; i++)
		graphedges[i].visited = 0;

	for(i = 0; i < cur_node_index; i++)
		graphnodes[i].visitied = 0;

	////
	if((graphnodes[sourcenode].type == 0 && graphnodes[destnode].type == 1)
		||(graphnodes[sourcenode].type == 1 && graphnodes[destnode].type == 0))
	{
		////if the two nodes are not at the same level

        if(SearchTwoLevels(sourcenode, destnode, media_nodes, num_medianodes))
			return true;
		else return false;
	}

	else if((graphnodes[sourcenode].type == 0 && graphnodes[destnode].type == 0)
		||(graphnodes[sourcenode].type == 1 && graphnodes[destnode].type == 1))
	{
		////if the two nodes are at the same level
		////We allow the search go back to previous level, but can we cross saddle to another level?
		////At this moment, we do not allow it cross the saddle level!! 11/06/05
	}

	else{
		////Two saddle may be connected at time-varying vector field
	}

	return false;
}


/*------------------------------------------------------------------------------------
* Search only two levels of the graph to find the two saddles on the two paths
* that connect the repeller and the attractor
------------------------------------------------------------------------------------*/
bool SearchTwoLevels(int sourcenode, int destnode, int *media_nodes, int &num_medianodes)
{
	int i, j;
	GraphNode *cur_node;
	int nodewanttoadd = -1;
	int nodewanttocompare = -1;

	NodeSearchArray[0] = sourcenode;
	cur_search_node_pos = 0;
	cur_endposition = 1;

	////1. find the possible saddles that connected with the source node
	cur_node = &graphnodes[sourcenode];
	for(i = 0; i < cur_node->nedges; i++)
	{
		////Through each edge, we find another not visited node
		if(graphedges[cur_node->edges[i]].node_index1 != cur_node->node_index)
		{
			nodewanttoadd = graphedges[cur_node->edges[i]].node_index1;
		}
		else
			nodewanttoadd = graphedges[cur_node->edges[i]].node_index2;

		////Add to search array
		if(cur_endposition >= maxlength_searcharray-1)
		{
			maxlength_searcharray += 20;
			NodeSearchArray = (int*)realloc(NodeSearchArray, sizeof(int)*maxlength_searcharray);
		}

		////Because we just search two layers here, we do not perform repeated testing
		////But if we consider more than two layer searching, we need to test repeatation first!!!!

		NodeSearchArray[cur_endposition] = nodewanttoadd;
		cur_endposition++;

		graphedges[cur_node->edges[i]].visited = 1;
	}

	////2. find the two paths that can reach destnode with length of 2
	cur_search_node_pos ++;

	for(i = cur_search_node_pos; i < cur_endposition; i++)
	{
		cur_node = &graphnodes[NodeSearchArray[cur_search_node_pos]];
		
		for(j = 0; j < cur_node->nedges; j++)
		{
			if(graphedges[cur_node->edges[j]].visited == 1)
				continue;

			if(graphedges[cur_node->edges[j]].node_index1 != cur_node->node_index)
			{
				nodewanttocompare = graphedges[cur_node->edges[j]].node_index1;
			}
			else
				nodewanttocompare = graphedges[cur_node->edges[j]].node_index2;
			
			if(nodewanttocompare == destnode)
			{
				if(num_medianodes >= MaxMediaNodes-1)
				{
					MaxMediaNodes += 20;
					MediaNodes = (int*)realloc(MediaNodes, sizeof(int)*MaxMediaNodes);
				}

				media_nodes[num_medianodes] = cur_node->node_index;
				num_medianodes++;
				break;
			}

		}
        cur_search_node_pos ++;
	}
	if(num_medianodes >= 2) return true;
	else return false;

}


/*-----------------------------------------------------------
* Finding the connected nodes on the paths that connect
* the input two nodes in the graph (only for singularities)
-----------------------------------------------------------*/
void SearchConnectComponents_adv(
	int *repellers, int num_repellers, 
	int *attractors, int num_attractors, 
	int *media_nodes, int &num_medianodes)
{
	int i;
	int cur_node;

	////Initial part
	maxlength_searcharray = 50;
	if(NodeSearchArray != NULL)
		free(NodeSearchArray);
    NodeSearchArray = (int*)malloc(sizeof(int) * maxlength_searcharray); 
	cur_endposition = 0;
    cur_search_node_pos = 0;

	////Reset the visited flag
	for(i = 0; i < cur_graphedge_index; i++)
		graphedges[i].visited = 0;

	for(i = 0; i < cur_node_index; i++)
		graphnodes[i].visitied = 0;

	////Begin searching from each repeller in the "repellers" list 
	for(i = 0; i < num_repellers; i++)
	{
		cur_node = repellers[i];

		TwoLevelDBFS(cur_node, attractors, num_attractors, media_nodes, num_medianodes);
	}

}


/*------------------------------------------------------------------------------------
* Search only two levels of the graph to find the saddles on the paths
* that connect the input repeller and the attractors
------------------------------------------------------------------------------------*/
void TwoLevelDBFS(int source, int *attractors, int num_attractors,
				  int *media_nodes, int &num_medianodes)
{
	int i, j, k;
	GraphNode *cur_node;
	int nodewanttoadd = -1;
	int nodewanttocompare = -1;

	NodeSearchArray[0] = source;
	cur_search_node_pos = 0;
	cur_endposition = 1;

	////1. find all the saddles that connected with the source node
	////we may repeated adding the saddles
	cur_node = &graphnodes[source];

	for(i = 0; i < cur_node->nedges; i++)
	{
		////Through each edge, we find another not visited node
		nodewanttoadd = graphedges[cur_node->edges[i]].node_index2;
		
		if(graphnodes[nodewanttoadd].visitied == 1) //the saddle has been visited by other repellers
			continue;

		if(graphnodes[nodewanttoadd].type == 1) //we do not add attractor here
			continue;

		////Add to search array
		if(cur_endposition >= maxlength_searcharray-1)
		{
			maxlength_searcharray += 20;
			NodeSearchArray = (int*)realloc(NodeSearchArray, sizeof(int)*maxlength_searcharray);
		}

		if(!IsRepeated(NodeSearchArray, nodewanttoadd, cur_endposition))
		{
			NodeSearchArray[cur_endposition] = nodewanttoadd;
			cur_endposition++;
		}

		graphedges[cur_node->edges[i]].visited = 1;

		//graphnodes[nodewanttoadd].visitied = 1;
	}

	////2. find all the paths that can reach any nodes in the "attractors"

	////we may need to consider those "repeller-attractor" directly connected cases
	////Such as a simple limit cycle with a singularity in its center 11/20/05

	cur_search_node_pos ++;

	for(i = cur_search_node_pos; i < cur_endposition; i++)
	{
		cur_node = &graphnodes[NodeSearchArray[cur_search_node_pos]];

		if(cur_node->visitied == 1) //the saddle has been visited by other repellers
			continue;

		if(cur_node->type == 1)     //we do not consider an attractor as a middle point here! 11/20/05
			continue;
		
		for(j = 0; j < cur_node->nedges; j++)
		{
			if(graphedges[cur_node->edges[j]].visited == 1) //avoid to go back to repellers
				continue;

			//Always remember this is a directed graph
			nodewanttoadd = graphedges[cur_node->edges[j]].node_index2;

			if(nodewanttoadd == cur_node->node_index) 
				//the edge is coming from other repellers not current repeller "source"
				//we avoid it to go back to repeller layer
				continue;

			//compare with the nodes in the "attractors"
			for(k = 0; k < num_attractors; k++)
			{
				nodewanttocompare = attractors[k];

				if(nodewanttoadd == nodewanttocompare)
				{
					if(num_medianodes >= MaxMediaNodes-1)
					{
						MaxMediaNodes += 20;
						MediaNodes = (int*)realloc(MediaNodes, sizeof(int)*MaxMediaNodes);
					}

					if(!IsRepeated(MediaNodes, cur_node->node_index, num_medianodes))
					{
						media_nodes[num_medianodes] = cur_node->node_index;
						num_medianodes++;
						cur_node->visitied = 1;   //Tell the program that the saddle has been added
					}

					break;
				}
			}

		}
        cur_search_node_pos ++;
	}

}

/*------------------------------------------------------------------------------------
* Find the connected nodes on the paths that connect the input two nodes
* in the graph (one is a limit cycle, the other is a singularity)
------------------------------------------------------------------------------------*/
bool SearchbyDFS(int sourcenode, int destnode, int *media_nodes, int &num_medianodes)
{
	int i, j;
	GraphNode *cur_node;
	int nodewanttoadd = -1;
	int nodewanttocompare = -1;

	NodeSearchArray[0] = sourcenode;
	cur_search_node_pos = 0;
	cur_endposition = 1;

	////1. find the possible saddles that connected with the source node
	cur_node = &graphnodes[sourcenode];
	
	for(i = 0; i < cur_node->nedges; i++)
	{
		////Through each edge, we find another not visited node
		if(graphedges[cur_node->edges[i]].node_index1 != cur_node->node_index)
			nodewanttoadd = graphedges[cur_node->edges[i]].node_index1;
		else
			nodewanttoadd = graphedges[cur_node->edges[i]].node_index2;

		////Add to search array
		if(cur_endposition >= maxlength_searcharray-1)
		{
			maxlength_searcharray += 20;
			NodeSearchArray = (int*)realloc(NodeSearchArray, sizeof(int)*maxlength_searcharray);
		}

		////Because we just search two layers here, we do not perform repeated testing
		////But if we consider more than two layer searching, we need to test repeatation first!!!!

		NodeSearchArray[cur_endposition] = nodewanttoadd;
		cur_endposition++;

		graphedges[cur_node->edges[i]].visited = 1;
	}


	////2. find the two paths that can reach destnode with length of 2
	cur_search_node_pos ++;

	for(i = cur_search_node_pos; i < cur_endposition; i++)
	{
		cur_node = &graphnodes[NodeSearchArray[cur_search_node_pos]];
		
		for(j = 0; j < cur_node->nedges; j++)
		{
			if(graphedges[cur_node->edges[j]].visited == 1)
				continue;

			if(graphedges[cur_node->edges[j]].node_index1 != cur_node->node_index)
				nodewanttocompare = graphedges[cur_node->edges[j]].node_index1;
			else
				nodewanttocompare = graphedges[cur_node->edges[j]].node_index2;
			
			////Here we forbit the searching across the saddle level
			if(graphnodes[nodewanttoadd].type != graphnodes[sourcenode].type 
				&& graphnodes[nodewanttoadd].type != 2)
				continue;

			////If we meet the destinate node, it means that we find a path, stop the searching 11/06/05

			if(nodewanttocompare == destnode)
			{
				if(num_medianodes >= MaxMediaNodes-1)
				{
					MaxMediaNodes += 20;
					MediaNodes = (int*)realloc(MediaNodes, sizeof(int)*MaxMediaNodes);
				}

				media_nodes[num_medianodes] = cur_node->node_index;
				num_medianodes++;
				break;
			}

		}
        cur_search_node_pos ++;
	}

	return false;
}


/*------------------------------------------------------------------------------------
* Finding the connected nodes on the paths that connect the input two nodes
* in the graph (one is a limit cycle, the other is a singularity)
------------------------------------------------------------------------------------*/
bool SearchLimitCycleConnectComponents(int limitnode, int singnode, int *media_nodes, int &num_medianodes)
{

	////if the two nodes have direct connection and they have opposite feature
	////i.e., one is a repeller, the other is an attractor
	////just return true, and set num_mediannodes = 0;

	////else, we need to perform the similar search as above
	return false;
}




/*
Mark those being cancelled nodes as 'cancelled'
*/
void MarkCancel(int *repellers, int num_repellers, 
	int *attractors, int num_attractors, 
	int *media_nodes, int num_medianodes)
{
	int i;

//int *Pre_CancelledNode = (int *)malloc(sizeof(int) * 20);
//int num_cancellednode = 0;

	num_cancellednode = 0; //3/29/06 for video capture

	for(i = 0; i < num_repellers; i++)
	{
		graphnodes[repellers[i]].cancelled = 1;
		Pre_CancelledNode[num_cancellednode] = repellers[i];  //3/29/06 for video capture
		num_cancellednode ++;
	}
	
	for(i = 0; i < num_attractors; i++)
	{
		graphnodes[attractors[i]].cancelled = 1;
		Pre_CancelledNode[num_cancellednode] = attractors[i];  //3/29/06 for video capture
		num_cancellednode ++;
	}
	
	for(i = 0; i < num_medianodes; i++)
	{
		graphnodes[media_nodes[i]].cancelled = 1;
		Pre_CancelledNode[num_cancellednode] = media_nodes[i];  //3/29/06 for video capture
		num_cancellednode ++;
	}
}


/*
Save previous detected singularities
*/
void SavePreSingularities()
{
	int i;

	////Allocate the space for previous singularity list
	if(presingularityList != NULL)
		free(presingularityList);

	presingularityList = (PreSingularityList *)malloc(sizeof(PreSingularityList)*cur_singularity_index);
    num_presingularities = 0;

	////Store the information for singularities
	for(i = 0; i < cur_singularity_index; i++)
	{
		presingularityList[i].triangleID = singularities[i].Triangle_ID;
		presingularityList[i].type = singularities[i].type;
		presingularityList[i].nodeindex = singularities[i].node_index;
	}

	num_presingularities = cur_singularity_index;

}


/**/
void SavePreLimitCycle()
{
	int i, j;

	////Allocate the space
	if(prelimitcyclelist != NULL)
	{
		for(i = 0; i < num_prelimitcycles; i++)
		{
			if(prelimitcyclelist[i].cellcycle != NULL)
				free(prelimitcyclelist[i].cellcycle);
		}
	}

	////Hard code here, we need to extend to dynamic structure later
	prelimitcyclelist = (PreLimitCycle *)malloc(sizeof(PreLimitCycle)*cur_limitcycle_index);
	for(i = 0; i < cur_limitcycle_index; i++)
	{
		prelimitcyclelist[i].cellcycle = NULL;
	}

	////save the cell cycles for all the limit cycles
	for(i = 0; i < cur_limitcycle_index; i++)
	{
		prelimitcyclelist[i].cellcycle = (int*)malloc(sizeof(int)*limitcycles[i].num_triangles);

		for(j = 0; j < limitcycles[i].num_triangles; j++)
		{
			prelimitcyclelist[i].cellcycle[j] = limitcycles[i].cellcycle[j];
		}

		prelimitcyclelist[i].num_cells = limitcycles[i].num_triangles;
		prelimitcyclelist[i].type = limitcycles[i].type;
		prelimitcyclelist[i].nodeindex = limitcycles[i].node_index;
	}

	num_prelimitcycles = cur_limitcycle_index;
}

/*
Judge whether the captured singularity is an old one or not
*/
bool IsPreSingularity(int triangle, int type, int &org_nodeid)
{
	int i;

	for(i = 0; i < num_presingularities; i++)
	{
		if(triangle == presingularityList[i].triangleID && type == presingularityList[i].type)
		{
			org_nodeid = presingularityList[i].nodeindex;
			return true;
		}
	}

	return false;
}


/*
Judge whether the captured limit cycle is an old one or not
*/
bool IsPreLimitCycle(int limitID, int &org_nodeid)
{
	int i, j;

	int cur_t;
	int count = 0;

	for(i = 0; i < num_prelimitcycles; i++)
	{
		count++;
		for(j = 0; j < limitcycles[limitID].num_triangles; j++)
		{
			cur_t = limitcycles[limitID].cellcycle[j];
			
			if(IsRepeated(prelimitcyclelist[i].cellcycle, cur_t, prelimitcyclelist[i].num_cells))
			{
				count++;
			}
		}

		if(count >= (int)(0.95 * limitcycles[limitID].num_triangles))
		{
			org_nodeid = prelimitcyclelist[i].nodeindex;
			return true;
		}

	}

	return false;
}


/*
Update the C-graph, not rebuild it. We just add new nodes and edges for the new appearing singularities
Note that we should store the previous singularities before performing any cancellation
*/
void UpdateGraph()
{
	int i, j;
	int trajectory_id;
	int lineseg_id;
	int last_triangle;

	int node1, node2;

	int org_nodeid = -1;
	
	//int r_counter, a_counter, s_counter; //for labeling 3/14/06


	for(i = 0; i < cur_singularity_index; i++)
	{
		if(singularities[i].type != CWCENTER && singularities[i].type != CCWCENTER)
		{
			////if it is not previous singularity
			if(!IsPreSingularity(singularities[i].Triangle_ID, singularities[i].type, org_nodeid))
			{
				graphnodes[cur_node_index].node_index = cur_node_index;
				graphnodes[cur_node_index].singularityID = i;
				graphnodes[cur_node_index].LimitCycleID = -1;
				graphnodes[cur_node_index].cancelled = 0;
				singularities[i].node_index = cur_node_index;
				
				if(singularities[i].type == SOURCE || singularities[i].type == RFOCUS)
				{
					graphnodes[cur_node_index].type = 0;   ////it is a repeller
					graphnodes[cur_node_index].labelindex = r_counter;
					r_counter++;
				}
				else if(singularities[i].type == SINK || singularities[i].type == AFOCUS)
				{
					graphnodes[cur_node_index].type = 1;   ////it is a attractor
					graphnodes[cur_node_index].labelindex = a_counter;
					a_counter++;
				}
				else
				{
					graphnodes[cur_node_index].type = 2;   ////it is a saddle
					graphnodes[cur_node_index].labelindex = s_counter;
					s_counter++;
				}

				graphnodes[cur_node_index].nedges = 0;
				graphnodes[cur_node_index].edges = NULL;

				cur_node_index ++;
			}

			else
			{
				////if it is the previous singularity, update it with its original node index
				graphnodes[org_nodeid].singularityID = i; //is that possible that previuos node 
				                                         //corresponding to a limit cycle?
				singularities[i].node_index = org_nodeid;
			}
		}
	}

	//////Update limit cycle
    for(i = 0; i < cur_limitcycle_index; i++)
	{
		////If it is not a previous detected limit cycle, add the node
		if(!IsPreLimitCycle(i, org_nodeid))
		{
			graphnodes[cur_node_index].node_index = cur_node_index;
			graphnodes[cur_node_index].LimitCycleID = i;
			graphnodes[cur_node_index].singularityID = -1;
			graphnodes[cur_node_index].cancelled = 0;
			limitcycles[i].node_index = cur_node_index;
			if(limitcycles[i].type == 0)  ////it is a repeller
			{
				graphnodes[cur_node_index].type = 0;
				graphnodes[cur_node_index].labelindex = r_counter;
				r_counter++;
			}
			else                          ////it is an attractor
			{
				graphnodes[cur_node_index].type = 1;
				graphnodes[cur_node_index].labelindex = a_counter;
				a_counter++;
			}

			graphnodes[cur_node_index].nedges = 0;
			graphnodes[cur_node_index].edges = NULL;

			cur_node_index ++;
		}

		else
		{
			graphnodes[org_nodeid].LimitCycleID = i;
			limitcycles[i].node_index = org_nodeid;
		}
	}

	//2. Build the edges
	////The graph is a directed graph!!!! 11/19/05, the direction of an edge is node1->node2
	////I.E. node1 should be a repeller and node2 should be an attractor
	////Note that, during the building of each edge, we need to check whether it is a previous edge!

	////(a) Find the edges between saddles and other singularities or limit cycles
	for(i = 0; i < cur_singularity_index; i++)
	{

		////Begin from each saddle, check its separatrices
		////if the last triangle of the separatrix contains singularity, build an edge between the singularity and current saddle
		if(singularities[i].type == SADDLE)
		{
		////if the last triangle of the separatrix does not contain any singularity
		////it maybe forms a closed streamline (we do not consider the 'center' here now) or reach the boundary
		////Build an edge between the limit cycle and current saddle

			////Then, follow the separatrices of the saddle 
			for(j = 0; j < 4; j++)
			{
				switch(j)
				{
				case 0:
					trajectory_id = separatrices[singularities[i].separtices].sep1;
					break;
				case 1:
					trajectory_id = separatrices[singularities[i].separtices].sep2;
					break;
				case 2:
					trajectory_id = separatrices[singularities[i].separtices].sep3;
					break;
				case 3:
					trajectory_id = separatrices[singularities[i].separtices].sep4;
					break;
				}

				lineseg_id = num_linesegs_curtraj[trajectory_id];
				if(lineseg_id <= 0) continue;
				last_triangle = trajectories[trajectory_id][lineseg_id-1].Triangle_ID;
				if(last_triangle >= 0 && Object.flist[last_triangle]->contain_singularity > 0) ////contains singularity
				{
					org_nodeid = -1;
					int singid = Object.flist[last_triangle]->singularity_index;
					
					if(Object.flist[last_triangle]->singularity_index < 0)
						continue;

					if(IsPreSingularity(singularities[i].Triangle_ID, SADDLE, org_nodeid)
						&& IsPreSingularity(last_triangle, singularities[singid].type, org_nodeid))
						continue;

					////Get the node index for the other singularity
					node1 = singularities[i].node_index;
					node2 = singularities[Object.flist[last_triangle]->singularity_index].node_index;

					if(singularities[Object.flist[last_triangle]->singularity_index].type == SOURCE
						||singularities[Object.flist[last_triangle]->singularity_index].type == RFOCUS)
					{
						////here, saddle is the attractor
						node2 = node1;
						node1 = singularities[Object.flist[last_triangle]->singularity_index].node_index;
					}

					////Add to graph edge array
					AddToEdgeArray(node1, node2, cur_graphedge_index);
				   
					////Need to add to the edge list of node1 and node2
					AddEdgeToNode(node1, cur_graphedge_index-1);
					AddEdgeToNode(node2, cur_graphedge_index-1);
				}
			}

			////(b) build the connections between saddles and limit cycles 2/15/06
			for(j = 0; j < singularities[i].num_connected_limitcycles; j++)
			{
				org_nodeid = -1;

				if(IsPreSingularity(singularities[i].Triangle_ID, SADDLE, org_nodeid)
					&& IsPreLimitCycle(singularities[i].connected_limitcycles[j], org_nodeid))
					continue;

				node1 = singularities[i].node_index;
				node2 = limitcycles[singularities[i].connected_limitcycles[j]].node_index;

				if(node2 < 0)
					continue;

				if(limitcycles[singularities[i].connected_limitcycles[j]].type == 0)
				{
					node2 = node1;
					node1 = limitcycles[singularities[i].connected_limitcycles[j]].node_index;
				}
				

				////If there is no connection between them, connect them
				if(!IsConnected(node1, node2))
				{
					////Add to graph edge array
					AddToEdgeArray(node1, node2, cur_graphedge_index);
					
					////Need to add to the edge list of node1 and node2
					AddEdgeToNode(node1, cur_graphedge_index-1);
					AddEdgeToNode(node2, cur_graphedge_index-1);
				}
			}
		}
	}

	////(c) Find the edges between limit cycles and other singularities
	////Here we build the saddle-limit cycle connections first
	////which can be fulfilled in previous step 2/15/06

	////Here we begin the searching from the connected list in limit cycle data structure
	////Note that we also need to consider the directions of the connections
	for(i = 0; i < cur_limitcycle_index; i++)
	{
		node1 = limitcycles[i].node_index;

		for(j = 0; j < limitcycles[i].num_connectedsaddles; j++)
		{
			node2 = singularities[limitcycles[i].connected_saddle[j]].node_index;
			
			if(limitcycles[i].type == 1)
			{
				////this is an attractive limit cycle, hence saddle becomes a repeller here
				node1 = node2;
				node2 = limitcycles[i].node_index;
			}

			if(node1 < 0 || node2 < 0)
				continue;

			if(HasInterval(node1, node2))  ////there is saddel connecting them
				continue;

			////Add to graph edge array
			if(AddToEdgeArray(node1, node2, cur_graphedge_index))
			{
				////Need to add to the edge list of node1 and node2
				AddEdgeToNode(node1, cur_graphedge_index-1);
				AddEdgeToNode(node2, cur_graphedge_index-1);
			}
		}
	}

//int *Pre_CancelledNode = (int *)malloc(sizeof(int) * 20);
//int num_cancellednode = 0;
	////Set the cancelled node for previous cancellation!
	for(i = 0; i < cur_node_index; i++)
	{
		if(IsRepeated(Pre_CancelledNode, graphnodes[i].node_index, num_cancellednode))
			graphnodes[i].cancelled = 1;
	}


	////(c) Find the edges between embeded limit cycles
	//for(i = 0; i < cur_limitcycle_index; i++)
	//{
	//	node1 = limitcycles[i].node_index;

	//	for(j = 0; j < limitcycles[i].num_connectedcycles; j++)
	//	{
	//		////if there has already one edge between them, continue
	//		node2 = limitcycles[limitcycles[i].connected_limitcycle[j]].node_index;

	//		if(limitcycles[limitcycles[i].connected_limitcycle[j]].type == 
	//			limitcycles[i].type)
	//			continue;

	//		if(limitcycles[i].type == 1)
	//		{
	//			////this is an attractive limit cycle, hence the other becomes a repeller here
	//			node1 = node2;
	//			node2 = limitcycles[i].node_index;
	//	    }

	//		if(HasInterval(node1, node2))  ////there is saddel connecting them
	//			continue;

	//		if(!AddToEdgeArray(node1, node2, cur_graphedge_index))
	//			continue;

	//		else
	//		{
	//			AddEdgeToNode(node1, cur_graphedge_index-1);
	//			AddEdgeToNode(node2, cur_graphedge_index-1);
	//		}
	//	}
	//}

}