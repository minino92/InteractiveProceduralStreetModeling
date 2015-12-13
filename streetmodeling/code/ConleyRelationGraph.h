////ConleyRelationGraph.h

/* The Data structure for previous captured singularities*/
typedef struct PreSingularityList{
	int triangleID;
	int type;
	int nodeindex;
}PreSingularityList;

/* The Data structure for previous captured limit cycles*/
typedef struct PreLimitCycle{
	int *cellcycle;
	int num_cells;
	int type;
	int nodeindex;
}PreLimitCycle;

void BuildGraph();
bool AddToEdgeArray(int node1, int node2, int &cur_index);
bool IsRepeatedEdge(int node1, int node2, int cur_index);
bool IsConnected(int node1, int node2);

void AddEdgeToNode(int node, int edgeindex);

void LayOutNodes();

////Search the connected components
bool SearchConnectComponents(int sourcenode, int destnode, int *media_nodes, int &num_medianodes);
bool SearchLimitCycleConnectComponents(int limitnode, int singnode, int *media_nodes, int &num_medianodes);

bool SearchTwoLevels(int sourcenode, int destnode, int *media_nodes, int &num_medianodes);
bool SearchbyDFS(int sourcenode, int destnode, int *media_nodes, int &num_medianodes);


////More generic case
void TwoLevelDBFS(int source, int *attractors, int num_attractors,
				  int *media_nodes, int &num_medianodes);
void SearchConnectComponents_adv(
	int *repellers, int num_repellers, 
	int *attractors, int num_attractors, 
	int *media_nodes, int &num_medianodes);

void MarkCancel(int *repellers, int num_repellers, 
	int *attractors, int num_attractors, 
	int *media_nodes, int num_medianodes);

void UpdateGraph();
bool IsPreSingularity(int triangle, int type, int &org_nodeid);
void SavePreSingularities();
void SavePreLimitCycle();
bool IsPreLimitCycle(int limitID, int &org_nodeid);

