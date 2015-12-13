
/*
  FindSCC.h
*/

//#include "Polygon3D.h"

////Data structure for finding the strongly connected component 07/10/06
typedef struct GraphNode2{
	int node_index;                   //Node index (unique)
	int *edges;                       //An array of the edges incident to the node
	int nedges;
	int levels;                       //the level during DFS
	int parent;                       //store the parent of current nodes for backward tracking to build the SCC
	int sscomp_index;
	int visited;
} GraphNode2;


typedef struct SCComponent{
	int *nodes;
	int num_nodes;
	int valid;
	int num_seppts;
	int num_attpts;
	int num_boundaries;
	/*  Two lists for the limit cycles and singularities inside the component */
	int *singular_tri;           //the list of triangles containing singularity 
	int num_singularities;
	int *limitcycles;
	int num_limitcycles;
		
	bool bad_flag;              /*Mark whether the corresponding region is bad*/

	int node_index;             /*the index of the corresponding node in the MCG graph 03/03/07*/
}SCComponent;   //the data structure for only one strongly connected component

typedef struct SCCList{
	SCComponent *scccomponents;
	int num_sccs;
}SCCList;       //the list of all strongly connected components in a directed graph


/* The following is an positive integer stack */
class int_stack{
	int *stack;
	int size;
	int top;
public:
	int_stack(int size = 0) { this->size = size; top = 0;}
	int pop(){ top--; if(top < 0) return -1; return stack[top];}
	bool push(int elem)
	{
		if(top >= size)
		{
			int *temp = stack;
			if((stack = (int*)realloc(stack, sizeof(int)*(size+100))) == NULL)
			{
				stack = temp;
				return false;
			}
			size += 200;
		}
		stack[top] = elem;
		top++;
		return true;
	}
};


void BuildConnectedGraph(double);
void SCC_AddEdgeToNode(int node, int edgeindex);
void SCC_AddToEdge(int node1, int node2, int &cur_index);
void ReverseEdges();
void FindSCC();
void mark_All_Valid_SCCS();

//void DFS(GraphNode2, int, int);

void BuildSCCElemList();
void CalJacobianForWholeMesh();
void build_DirGraph_Tri_no_Tau(int tri, int backward);
void build_DirGraph_no_Tau();

//void remove_outedges_from(int node);
//bool del_one_edge_from(int *edges, int &n, int del_e);


void AllocForSCC();
void InitForSCC();
void finalize_scclist();


bool path_From_T1_To_T2(int n1, int n2);
void count_MorseSets();

