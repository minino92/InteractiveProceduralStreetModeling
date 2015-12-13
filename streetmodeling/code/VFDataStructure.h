//VFDataStructure.h

/* This file define some data structure that will be used by this tool for vector field design
The data structures defined here include singular element, regular element, captured singularities, limit cycles
underneath mesh structure, trajectories and separatrix
*/
#pragma once

#include "lib/icMatrix.h"
#include "lib/icVector.h"

#define M_PI 3.14159265358979323846264338327950288       ////Move to the datastructure.h file 06/23/05

#define SELECTBUFFERSIZE 128

////////////////////////////////////////////////////////
////Load name for object selection
#define NAMEOFSINGELEM      1           ////name of singular element for mouse selection
#define NAMEOFREGELEM       2001        ////name of regular element for mouse selection
#define NAMEOFSINGCONTROL   3000        ////name of singular element control points for mouse selection
#define NAMEOFREGCONTROL    4000        ////name of regular element control points for mouse selection
#define NAMEOFSINGULARITY   5000        ////name of singularities for mouse selection of topology editing
#define NAMEOFLIMITCYCLES   7000        ////name of limit cycles for limit cycle editing
#define NAMEOFSHAPECONTROL  8000        ////name of shape control point for limit cycle shape controlling
#define NAMEOFBRUSHES       10001       ////name of the control points on brushes (or sketches)
#define NAMEOFINTERSECTS    11001       ////name of the intersections of the streets 
#define NAMEOFSEEDS         20001       ////name of the initial seeds
#define NAMEOFMAJROADS      26001       ////name of the major roads
#define NAMEOFTRIANGLE      30001       ////name of triangles for mouse selection

////some constants
const int NUMTRACINGTRIANGLE = 700;
const double DistanceThreshold = 1e-5;
const double RegularStrength = 0.5;  //the scaler for the strength of regular element
const double ARROWSCALE = 0.07;
const double EDITBOXSIZE = 0.04;
const double INOROUTJUDGE = 0.015;     //for single separatrices
const double SEPARATRIXSTEP = 0.01;     //for separatrices, maybe a better setting for separatrices drawing
const double INTERSECTIONERROR = 1e-5; //the error threshold between beginning point and next intersection point
                                       //for closed streamline tracing
const double SMOOTHSTEP = 10000.;        //speed control for region smoothing
const double SCALESTEP = 15.;   //The scaler step control for scale
const double ROTATESTEP = 10.;   //The scaler step control for rotation

const int InitSaddleRegionLength = 200;  ////tracing n's triangles for each direction
const double ScaleforNewStrength = 100;

/*----------------------------------------------*/
//Edit box structure

typedef struct EditBox{
	icVector2 p1, p2, p3, p4, Up;  ////stands for the 4 points of the edit box
}EditBox;


/*----------------------------------------------*/
//singular element

typedef struct SingularElem{
	int ID;
	int Triangle_ID;                //the triangle that contains the singular element
    double centerx, centery;        //center under global frame
    int type;                       //singular type (0 –Null, 1—source, 2—sink, 3—saddle, 4—cwcenter, 
                                       //5—ccwcenter
    icMatrix3x3  transform_matrix;  //transformation matrix for user editing
	icMatrix3x3  Jacobian;
    //visual control icons ( edit box, control points)
	EditBox editbox;
	EditBox cur_editbox;
    //double rate of decreasing (optional) ---not consider at this moment

	////variables may be useful in editing
	double rotang;
	double sx, sy;
	double s;
	bool deleted;
} SingularElement;


/*----------------------------------------------*/
//regular element

typedef struct RegularElem{
	int ID;
    double base[2];                //base under global frame
    icVector2 Direct;                  //direction in global frame
    int type;                       //singular type (0 –regular, 1—convergent, 2—divergent)
    icMatrix3x3  transform_matrix;  //transformation matrix for user editing
    icMatrix3x3 transposeRot;   
	//visual control icons (control points)
    //double rate of decreasing (optional) ---not consider at this moment

	int basis_triangle;            //the triangle that contains the basis
	
	////variables may be useful in editing
	double rotang;
	double s;
} RegularElement;


/*----------------------------------------------*/
//singularities

////Upgrade: You may need to store the local coordinates of the singularity

typedef struct Singularities{
    double gcx, gcy;                 //global center
	float alpha[3];                  //the barycentric coordinates of the fixed point inside the triangle
    int Triangle_ID;                 //which triangle it belongs to
    double local_coordinates[3];     //we can choose to use 2D local coordinates or barycentric coordinates 
    int type;                        //singular type (0-–Null, 1—source, 2—sink, 3—saddle, 4—cwcenter, 
                                     //5—ccwcenter, 6—repel center, 7—attractive center )
    icMatrix3x3 Jacobian;            //local Jacobian matrix for this singularity
    icVector2 incoming, outgoing;    //for separatrices calculation beginning from a saddle
    //visual icons ( we need not store it, the display routine will use different colors to mark different kinds of singularities)
	
	/*---------Variables for Conley relation graph 10/14/05---------------*/
	int node_index;                  //the index of the node associates with this singularity
	int separtices;                  //the index of the separatrix group belongs to the saddle 10/13/05
	int selectedtraj;                //the index of selected separatrix for modification 12/23/05
	int *connected_limitcycles;      //An array of limit cycles that the saddle can reach
	int num_connected_limitcycles;   //number of the limit cycles inside the array 'connected_limitcycles'

	int connected;                   //1--connected with other element(s) 08/06/06
} Singularities;


/*----------------------------------------------*/
//line segment structure for trajectory

typedef struct LineSeg{
	int Triangle_ID;                    //which triangle this line segment locates in
	double start[2], end[2];            //local coordinates for start and end points
	double gstart[2], gend[2];          //global coordinates for start and end points
	double length;                      //we may need to store the length of the line segment
} LineSeg;

////For temporary trajetory 
typedef struct CurvePoints{
	double lpx, lpy;
	double gpx, gpy;
	double length;
	int triangleid;
}CurvePoints;


//single trajectory data structure
//typedef struct Trajectory{
//	LineSeg *line_segs;
//	int num_linesegs;
//	int cur_MaxNumLinesegs;
//} Trajectory;

/*----------------------------------------------*/
//structure for a group of separatrices (only 4 trajectories for each group)
typedef struct Separatrices{
	int sep1, sep2, sep3, sep4; //1,3 -- outgoing || 2,4 -- incoming
	int connect1, connect2, connect3, connect4; //whether the corresponding sep is a connection orbit during
	                                            //cancellation
	double length1, length2, length3, length4;  //save the flow lengthes of the separatrices 08/10/06

} Separatrices;

/*----------------------------------------------*/
//limit cycle
typedef struct LimitCycle{
	int *cellcycle;                   //triangles list of closed cell cycle
	int num_triangles;                //number of triangles in the cell cycle
	LineSeg *closed_streamline;       //closed streamline of the limit cycle
	int num_linesegs;                 //number of line segments in the streamline
	int singularID;                   //the index of the center singularity of the limit cycle if exists
	int type;                         //type of the limit cycle 0--repeller, 1--attractor

	////Operation legend for mouse selection and editing
	double legend_center[2];
	double legend_base[2];

	////member for Conley relation graph 10/13/05
	int node_index;
	int *connected_saddle;            //An array of all saddles that can reach the limit cycle
	int num_connectedsaddles;
	double *flow_length_connectedsing;

	int *connected_limitcycle;        //An array of other limit cycles that connect with current limit cycle
	int num_connectedcycles;	
	double *flow_length_connectedcycle;

	int connected_l, connected_r;    //1---connected with other element(s) 08/06/06

} LimitCycle;




/**************************************************************************/
////Data structure for the underneath meshes
///////////////////////////////////////////////////////////////////
extern struct Edge;

///////////////////////////////////////////////////////////////////
//Vertex structure

typedef struct Vertex {
	double x,y,z;             //the coordinates of the vertex under global frame
	double bx, by;        //jitter the coordinates a little bit
	double nx,ny,nz;
	double prob_on_path;
	void *other_props;        //other properties 
	icVector3 normal;         //normal at the vertex

	int *edges_id;            //edges incident to current vertex
	int Num_edge;             //number of edges
	Edge **edges;             //store the addresses of the edges incident to the vertex
	
	/*---------Variables for corner table---------------*/
	int *Corners;             //store the indexes of those corners associate with the vertexs
	int Num_corners;          //Number of the corners

	/*---------Variables for region smoothing---------------*/
	int OnBoundary; 
	int InRegion;
	int VertID;               //*Index in the whole object vertices list
	int RegionListID;         //Index in the inner vertices list
	
	/*---------Variables for pair cancellation---------------*/
	int repell_flag;          //for repeller neighborhood | 0--unknown, 1--on boundary, 2-- in region
	int attract_flag;         //for attractor neighborhood | 0--unknown, 1--on boundary, 2-- in region

	////For geodesic based search, now it is for limit cycle relocation
	int which_line;
	double distance;
	//int visited;
	bool visited;
	unsigned char tau_visited;

	/*----------Variable for finding the separation and attachment points --------------*/
	icVector2 vec;           //variable to store the vector on that vertex under global frame
	icVector2 vec_J;
	icVector2 evec[2];
	double length[3];

	icMatrix2x2 Jacobian;    /*we will calculate it according to the Jacobian around its one ring neighbor*/
	/*eigen vectors*/
	icVector2 major, minor;  /*for visualization only*/
	double major_ang, minor_ang;
	double tensor_major_ang;
	bool major_cos, major_sin, minor_cos, minor_sin;
	icVector2 tran_vec;      /*transfer eigen vector to a real vector field*/

	double mag_speed;

	/*------------------------------------------------------*/
	int *connected_limitcycle;
	int num_connectedlimitcycle;
	
	/*------------For local tracing----------------*/
	double *Anglist;   //store the angle allocated for each triangle on the tangent plane

	/*------------For adaptive \tau-----------*/ //03/15/06
	float tau[2];  /*previous two \tau for adaptive adjust \tau*/
	icVector2 endp[2];    /*the global coordinates of the end point*/
	int end_tri[2];       /*the triangle contains the end point*/
	bool end1_ornot;      /*tell program which tau it should use, 0--tau[0], 1--tau[1]*/
	bool done;           /*mark if it has been traced*/

	/*---- For asymmetric method----*/
	float gama_d, gama_r, gama_s;
	float hue, saturation;

	/*--- color mapping ---*/
	float vf_r, vf_g, vf_b;

	float f_color[3], b_color[3];

	float d_color[3];  /*for density visualization*/
	int s_count;       /*the number of the samples associated with this vertex*/


} Vertex;

///////////////////////////////////////////////////////
//Face structure
////A structure for streamline placement copy at 07/22/06
typedef struct SampleListInTriangle{
	int which_traj;
	int limitcycle;
	int which_sample;
}SampleListInTriangle;

typedef struct Face {
	unsigned char nverts;        // number of vertices in the list 
	int *verts;                  // indices of the vertices associated with current face list 
	void *other_props;           //other properties 
	double area;                 //area of the face
	double length[3], angle[3];  //???
	int nvisible;
	icVector3 normal;           //normal of the face
	//icVector3 center;           //

	Edge *edges[3];             //store the address of the associated edges

	//////The following members are used to define the local frame for each triangle
	int index;                  //index of current face

	icVector2 LX;               //local frame , x axis
	icVector2 LY;               //local frame,  y axis

	double xy[3][2];            //store the coordinates of its three vertices under local frame

	icVector2 direct_vec[3];    //store the directional vectors in local frame, we have 3 vectors for each triangle

	icMatrix2x2 Jacobian;       //every triangle mesh will have a local linear vector field (Jacobian Matrix)
	
	/*---------Variables for pair cancellation---------------*/
	int repell_inregion;        //whether the triangle is inside the repeller neighborhood
	int attract_inregion;       //whether the triangle is inside the attractor neighborhood

	int contain_singularity;    //flag to mark whether this triangle containing a singularity or not
	int singularity_index;      //id of the singularity contained by current triangle

	////Variables for limit cycle shape design
	int inDesignCellCycle;

	int contain_separatrix;     //for separatrix editing (grow region)1/3/06
	int fence_flag;             //to set fence to prevent the region growing crossing this triangle

	////Variables for new limit cycle detection 1/24/06
	int access_count;           //record the times that the local tracing accessing
	int discard;                //if it is not a good condition triangle, reuse it as boundary flag 3/16/06

	/*------------For streamline placement-----------*/ //2/18/06
	int *trajlist;              //a list of the trajectory indices that pass through the triangle
	int num_passing_trajs;      //number of the trajectories passing through the triangle

	/*----------------------------------------------*/
	int which_SSC, pre_SCC;              //which strongly component it belongs to

	/*----------For the Jacobian of the triangle--------- 07/17/06 */
	icVector2 eigen[2];
	double evalues[2];
	double center[2];           //store the center of the triangle

	/*----------For degenerate point detection (temp) ---------- 09/18/2007 */
	unsigned char degenerate_type;

	/*------------For streamline placement, using samplin point-----------*/ //copy at 07/22/06
	SampleListInTriangle *samplepts;
	int *sampleindex;
	int num_samplepts;
	int MaxSampNum;

} Face;


////////new structure to store the information of edge
////I use link to store the edge, because we donot know the number of the edge first
typedef struct Edge{
	int index;                //Id for specific edge
	int verts[2];             //The two points for specific edge
	int OppVerts[2];          //The two vertices that opposite to the edge
	int tris[2];              //Two neighbor faces that share this edge
	int visited;              //for my subdivision of the triangle mesh
	int MidPointID;           //The ID of the vertex of the middle point on the edge
	double length;            //store the length of the edge
	
	/*---------Variables for pair cancellation---------------*/
	int OnRepellBoundary;     //If the edge is at the boundary of repell region mark it as 1, otherwise 0 1/23/06
	int OnAttractBoundary;    //If the edge is at the boundary of repell region mark it as 1, otherwise 0 1/23/06
	int repell_visited;       //whether the edge has been visited during the boundary building
	int attract_visited;
	icVector2 repell_normal;  //outward normal of the repeller region at the edge 
	icVector2 attract_normal; //outward normal of the attractor region at the edge 

	icVector2 normal;

	////Variables for boundary edges list extraction
	int OnBoundary;
	int BoundaryVisited;

	/*----- To store the special points on the edge (two at most) 07/17/06 ----*/
	icVector2 attp, sep;
	int find_attp, find_sep;
	int att_visit, sep_visit;    //0 -- alive; 1--visited/dead;  2--too close/pending  07/23/06

	/***----For testing the SCC----***/
	int mixed;
	
	/***---- For finding the separation and attachment points ----***/
	icVector2 evec[2];
	icMatrix2x2 Jacobian;     //for calculate the decomposition of the Jacobian
	int valid;

	/**----- For recording the intersections -----**/
	icVector2 intersections[2];  //we store only the recent two intersections
	int num_intersections;
	//double pre_length;

	Edge *next;
}Edge;

///The following define the corner structure
///Because each triangle has 3 corner, so the number of the corner can be 
///decided by the number of the faces and times number of the vertices on each face
typedef struct Corner{
	int index;
	double angle;
	Edge *edge[2];    //two edges associated with this corner

	/////The corner operation 
	int v;            //the ID of the vertex of the corner
	int n;
	int p;
	int t;            //the triangle the corner belongs to
	Edge *e;          //the opposite edge of the corner
	int o;            //the index of its opposite corner
	int ot;           //the index of its opposite triangle for traversal 2/9/05

	int Edge_count;   //special variable for edges search 1/21/05
	
	/*----------------------------------------------------------*/
	////variables for singularities detection and local tracing
	double BeginAng, EndAng;
	double r;
	int orient;
}Corner;


//////////////////////////////////////////////////
////Object data structure for the triangular mesh
typedef struct Polygon3D {
	/*unsigned */int nverts, nfaces, nedges;   //number of vertices, faces, edges respectively
	Vertex **vlist;                        //list of vertices
	Face **flist;                          //list of faces
	//PlyOtherProp *vert_other,*face_other;  
	double area;
	double radius;                         //radius of the bounding sphere of the object
	double ave_edge_length;
	icVector3 center;                      //the center of the object

	Edge *elist;              //This is a different link from vertices and faces to store the unknown edges
	Edge **edgelist;

	////New added Corner information 08/17/05
	Corner **clist;
	int ncorners;

} Polygon3D;


/**************************************************************************/
////Data structures for singularities pair cancellation and movement, modified at 4/11/06
typedef struct TriangularRegion{
	int *trianglelist;
	int num;
	int MaxNumTrianglesInRegion;
} TriangularRegion;

typedef struct RegionBoundary{
	Edge **edgelist;
	int num;
	int MaxNumEdgesOnBoundary;
} RegionBoundary;

typedef struct InnerVertices{
	Vertex **vertslist;
	int num;
	int MaxNumVersInRegion;
} InnerVertices;

////Boundary list for boundary extraction of multiple boundary region 4/11/06
typedef struct BoundaryList{
	RegionBoundary *boundarylist;
	int cur_boundary_num;
	int MaxNumBoundaries;
}BoundaryList;

/**************************************************************************/
////Data structures for Conley relation graph
typedef struct GraphEdge{
	int edge_index;                   //Edge index (unique)
	int node_index1, node_index2;     //Indices of the two nodes this edge connects
	int visited;                      //For graph searching
	int trajID;                       //the trajetory index of the corresponding separatrix 2/23/06
	double flow_length;               //store the flow length of the corresponding separatrix 08/10/06
	//double geo_length;				  //store the geometry length of the two nodes 08/10/06
	bool cancelled;
} GraphEdge;

typedef struct GraphNode{
	int node_index;                   //Node index (unique)
	int singularityID;                //The singularity ID that the node associates with
	int LimitCycleID;                 //The limit cycle ID that the node associates with
	int type;                         //0~repeller; 1~attractor; 2~saddle
	int *edges;                       //An array of the edges incident to the node
	int nedges;
	int visitied;                     //For graph searching
	int labelindex;
	int cancelled;

	double pos_x, pos_y;              //the positions of the node in the graph
} GraphNode;

/*---------------------------------------------------------------------------*/
/*define two new data structures for MCG */
typedef struct MCGNode{
	int index;
	int scc_index;
	int *edges;
	int nedges;
	int parent;
	unsigned char type;               //0-repeller, 1-attractor, 3-saddle
	
	int visited;                     //For graph searching
	int labelindex;
	int cancelled;

	double pos_x, pos_y;              //the positions of the node in the graph
}MCGNode;

typedef struct MCGEdge{
	int edge_index;                   //Edge index (unique)
	int node_index1, node_index2;     //Indices of the two nodes this edge connects
	int visited;                      //For graph searching
	bool cancel;

	/*the variables to represent the region corresponds to this edge08/29/2007*/
	int *triangles;
	int ntris;

	/* two optional members for automatic simplification */
	//int trajID;                       //the trajetory index of the corresponding separatrix 2/23/06
	double flow_length;               //store the flow length of the corresponding separatrix 08/10/06
}MCGEdge;

/**************************************************************************/
///////////////////////////////////////////////////
//typedef double Point[2];

typedef struct Point
{
	double x, y;
	int cellid;
}Point;

typedef struct ctr_point {
  double x;
  double y;
  //double z;
  int cellid;
}ctr_point;

typedef struct Boundaryvert{
	int vertID;       //the index of the vertex
	icVector2 vec;    //the vector on the vertex
}Boundaryvert;

typedef struct Boundaryverts_List{
	Boundaryvert *boundaryverts;  //the list of the boundary vertices with the vectors on them
	int num_boundarverts;         //the number of the boundary vertices
}Boundaryverts_List;


///////////////////////////////////////////////////
enum Elementtype{
	NOTHING,
	SOURCE,
	SINK,
	SADDLE,
	CWCENTER,
	CCWCENTER,
	AFOCUS,
	RFOCUS,
	REGULAR
};

/*current color scheme  08/25/05
source: pure green (0, 1, 0)
repeller: light green (0, 1, 0.7)
sink:   pure red (1, 0, 0)
attractor:  orange (1, 0.5, 0)
saddle: pure blue (0, 0, 1)
center: light red (1, 0, 1)
*/

///////////////////////////////////////////////////
enum which_point{
	NON,
	LOWLEFT,
	UPPERLEFT,
	UPPERRIGHT,
	LOWRIGHT,
	LEFT,
	UPPER,
	RIGHT,
	BUTTOM,
	UPPERROTATE,
	ARROWBASE,
	ARROWHEAD
};

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*****************************************************************************/
/*for tensor field analysis 09/29/2007*/
typedef struct DegeneratePt{
	int degpt_index;
    double gcx, gcy;                 //global center
	float alpha[3];                  //the barycentric coordinates of the fixed point inside the triangle
    int Triangle_ID;                 //which triangle it belongs to
    double local_coordinates[3];     //we can choose to use 2D local coordinates or barycentric coordinates 
    int type;                        //degenerate type (0—wedge, 1—trisector, 2—node, 3—center, 
                                     //4—saddle, 5—other higher order)
    icMatrix2x2 ten;                 //local tensor for this degenerate point
	double s1_ang, s2_ang, s3_ang;   //the starting angles of the separatrices
	icVector2 s[3];
	unsigned char nseps;             //the number of separatrices
    //icVector2 s1, s2, s3;    //for separatrices calculation beginning from a saddle
    //visual icons ( we need not store it, the display routine will use different colors to mark different kinds of singularities)

	bool deleted;
} DegeneratePt;

/*the data structure for singular design element of a tensor field*/
typedef struct Degenerate_Design{
	int ID;
	int Triangle_ID;                //the triangle that contains the singular element
    double centerx, centery;        //center under global frame
    int type;                       //0 – wedge, 1—trisector, 2—node, 3—center, 4—saddle, 
                                       //5—ccwcenter
    icMatrix3x3  transform_matrix;  //transformation matrix for user editing
	//icMatrix3x3  Jacobian;
    /*visual control icons ( edit box, control points)*/
	EditBox editbox;
	EditBox cur_editbox;
    //double rate of decreasing (optional) ---not consider at this moment

	////variables may be useful in editing
	double rotang;
	double sx, sy;
	double s;
	bool deleted;

	////
	unsigned char which_region;    //for two level tensor field design 11/17/2007
}Degenerate_Design;

/*the data structure for regular design element of a tensor field*/
typedef struct TenRegularElem{
	int ID;
    double base[2], end[2];                //base under global frame
    icVector2 Direct;                  //direction in global frame
    unsigned char type;                       //singular type (0 –regular, 1—convergent, 2—divergent)
    icMatrix3x3  transform_matrix;  //transformation matrix for user editing
	icMatrix3x3  transpose_matrix;

    //icMatrix3x3 transposeRot;   
	//visual control icons (control points)
    //double rate of decreasing (optional) ---not consider at this moment

	int basis_triangle;            //the triangle that contains the basis
	
	////variables may be useful in editing
	//double originalang;
	double rotang;
	double s;
	
	bool deleted;
	
	////
	unsigned char which_region;    //for two level tensor field design 11/17/2007
}TenRegularElem;

/*************************************************************/
/*extra info for computing the intersections 10/02/2007*/
typedef struct LinesInOneCell{
	int whichtraj;    /*the index of the tensor line*/
	int start, end;  /*from the "start" line segement to the "end" one*/
}LinesInOneCell;

class LineInfo{
public:
	LinesInOneCell **lines;
	int nlines;
	int curMaxNum;

	LineInfo(int init_size = 0)
	{
		if (init_size == 0)
		{
			lines = NULL;
			nlines = curMaxNum = 0;
			return;
		}

		lines = (LinesInOneCell**)malloc(sizeof(LinesInOneCell*)*init_size);
		nlines = 0;
		if(lines == NULL)
			exit(-1);
		curMaxNum = init_size;
		for(int i=0; i<curMaxNum; i++)
			lines[i]=NULL;
	}

	~LineInfo()
	{
		int i;
		if(lines != NULL)
		{
			for(i=0; i<curMaxNum; i++)
			{
				if(lines[i] != NULL)
					free(lines[i]);
			}

			free(lines);
		}
	}

	void addNew(LinesInOneCell *l)
	{
		if(nlines>=curMaxNum)
		{
			if(!extend(1))
				exit(-1);  /*probably not enough memory*/
		}
		lines[nlines] = l;
		nlines++;
	}

	bool isFull()
	{
		if(nlines>=curMaxNum)
			return true;
		return false;
	}

	bool extend(int step)
	{
		LinesInOneCell **temp = lines;
		lines = (LinesInOneCell**)malloc(sizeof(LinesInOneCell*)*(curMaxNum+step));
		if(lines == NULL)
			return false;
		int i;
		for(i=0; i<curMaxNum; i++)
			lines[i]=temp[i];
		for(i=curMaxNum; i<curMaxNum+step; i++)
			lines[i]=NULL;
		curMaxNum += step;
		free(temp);
		return true;
	}

	bool is_repeated(LinesInOneCell *l)
	{
		int i;
		for(i=0;i<nlines;i++)
		{
			if(l->whichtraj==lines[i]->whichtraj)
				return true;
		}
		return false;
	}
};
/*************************************************************/

extern struct QuadEdge;
extern struct QuadCell;

typedef struct QuadVertex {
	double x,y/*,z*/;             //the coordinates of the vertex under global frame
	double vx,vy/*,vz*/;
	double prob_on_path;
	void *other_props;        //other properties 
	//icVector3 normal;         //normal at the vertex

	int *edges_id;            //edges incident to current vertex
	unsigned char Num_edge;             //number of edges
	QuadEdge **edges;             //store the addresses of the edges incident to the vertex

    unsigned char ncells;     //number of cells that share this vertex
	QuadCell **cells;         //the cells that share this vertex
	
	/*---------Variables for region smoothing---------------*/
	bool OnBoundary; 
	bool InRegion;
	int index;               //*Index in the whole object vertices list
	int RegionListID;         //Index in the inner vertices list
	
	/*---------Variables for pair cancellation---------------*/
	//unsigned char repell_flag;          //for repeller neighborhood | 0--unknown, 1--on boundary, 2-- in region
	//unsigned char attract_flag;         //for attractor neighborhood | 0--unknown, 1--on boundary, 2-- in region

	////For geodesic based search, now it is for limit cycle relocation
	int which_line;
	double distance;  /*for distance level set as well*/
	bool visited;
	unsigned char type;      /*0--far away; 1--in narrow band; 2--known 09/30/2007*/

	icMatrix2x2 Jacobian;    /*we will calculate it according to the Jacobian around its one ring neighbor*/
	icMatrix2x2 origin_ten;  /*the old Jacobian*/
	/*eigen vectors*/
	icVector2 major, minor;  /*for visualization only*/
	double major_ang, minor_ang;
	double tensor_major_ang;
	bool major_cos, major_sin, minor_cos, minor_sin;
	icVector2 tran_vec;      /*transfer eigen vector to a real vector field*/

	bool inland;             /*true -- in land; false -- in the water region*/
	bool inveg;              /*true -- in veg region;  false -- in regular land or water */

	unsigned char which_region;  /*for two level design 11/17/2007*/
	bool inbrushregion;

	/*   for multi-density tensor line placement   */
	double density;

	/*   for the asymmetric tensor field design 12/29/2007  */
	double phi;

} QuadVertex;


typedef struct QuadCell {
	unsigned char nverts;        // number of vertices in the list 
	int *verts;                  // indices of the vertices associated with current face list 
	int index;                  //index of current face
	void *other_props;           //other properties 
	//double area;                 //area of the face
	//icVector3 normal;           //normal of the face
	//icVector3 center;           //

	QuadEdge *edges[4];             //store the address of the associated edges
	int xstart, xend, ystart, yend;
	double x_start_coord, y_start_coord;

	/*---------Variables for pair cancellation---------------*/
	//bool repell_inregion;        //whether the triangle is inside the repeller neighborhood
	//bool attract_inregion;       //whether the triangle is inside the attractor neighborhood

	bool contain_degpt;    //flag to mark whether this triangle containing a singularity or not
	int degpt_index;      //id of the singularity contained by current triangle

	bool OnBoundary;

	/*----------For degenerate point detection (temp) ---------- 09/18/2007 */
	unsigned char degenerate_type;

	/*------------For streamline placement, using samplin point-----------*/ //copy at 07/22/06
	////SampleListInTriangle *samplepts;
	//////int *sampleindex;
	////int num_samplepts;
	////int MaxSampNum;
	
	SampleListInTriangle *maj_samplepts, *min_samplepts;
	int maj_nsamplepts, min_nsamplepts;
	int MAJMaxSampNum, MINMaxSampNum;

	bool visited;

	/*----- For calculating the intersections and constructing the street net*/
	LineInfo *majorlines, *minorlines;
	bool hasmajor, hasminor;

	int *intersectlist;   /*a list of intersections in this cell*/
	unsigned char nintersects;      /*the number of the intersections in the list*/

	int *streetgraphedgelist;
	unsigned char nstreetgraphedges;

	/* for region division 11/22/2007*/
	LineInfo *sketchlines;
	unsigned char which_region;

	bool in_region;
	//bool in_veg;                  /*true -- in veg region;  false -- in regular land or water */
	bool is_contour_cell;

	/*  record the water boundaries that cross this cell  */
	int *mapbounds;
	int nmapbounds;

} QuadCell;

typedef struct QuadEdge{
	int index;                //Id for specific edge
	int verts[2];             //The two points for specific edge
	int tris[2];              //Two neighbor faces that share this edge
	bool visited;              //for my subdivision of the triangle mesh
	//double length;            //store the length of the edge
	
	/*---------Variables for pair cancellation---------------*/
	bool OnRepellBoundary;     //If the edge is at the boundary of repell region mark it as 1, otherwise 0 1/23/06
	bool OnAttractBoundary;    //If the edge is at the boundary of repell region mark it as 1, otherwise 0 1/23/06
	bool repell_visited;       //whether the edge has been visited during the boundary building
	bool attract_visited;
	icVector2 repell_normal;  //outward normal of the repeller region at the edge 
	icVector2 attract_normal; //outward normal of the attractor region at the edge 

	//icVector2 normal;

	////Variables for boundary edges list extraction
	bool OnBoundary;
	bool BoundaryVisited;

	/***---- For finding the separation and attachment points ----***/
	icMatrix2x2 Jacobian;     //for calculate the decomposition of the Jacobian

	/**----- For recording the intersections -----**/
	icVector2 intersections[2];  //we store only the recent two intersections
	int num_intersections;
	//double pre_length;

	QuadEdge *next;
}QuadEdge;


////Object data structure for the triangular mesh
class QuadMesh {
public:	
	int nverts, nfaces, nedges;   //number of vertices, faces, edges respectively
	QuadCell **quadcells;                      //list of faces
	QuadVertex **quad_verts;				   //list of vertices
	int XDIM, YDIM;
	double xinterval, yinterval;
	double xstart, xend, ystart, yend;

	double area;
	double radius;                         //radius of the bounding sphere of the object
	double ave_edge_length;
	icVector3 center;                      //the center of the object

	QuadEdge *elist;              //This is a different link from vertices and faces to store the unknown edges
	QuadEdge **edgelist;

	QuadMesh();
	QuadMesh(int xdim, int ydim, double xstart, double xend, 
				   double ystart, double yend);
	bool gen_regquad_mesh(int xdim, int ydim, double xstart, double xend, 
					double ystart, double yend);
	bool gen_regquad_vertices(int xdim, int ydim, double xstart, double xend, 
					double ystart, double yend);
	void init_vertices();
	void gen_regquad_faces(int xdim, int ydim)  ;
	void finalize_quad_verts();
	void finalize_quad_cells();
	void construct_edges();
	QuadEdge  **Extend_Elist(QuadEdge **edge_link, int Num_edges);
	QuadCell  **extend_celllist_ver(QuadCell **cells, int ncells);
	void orient_edges_cells();
};


/** The min-heap structure/class for fast marching in quad mesh 09/30/2007**/

typedef struct HeapElem{
	int vertid;
	double T;
}HeapElem;


class MinHeap{
public:
	//Attributes
	HeapElem **elems;
	int nelems;
	int curMaxNum;

	MinHeap(int initsize = 0);
	~MinHeap();

	//Operations
	int FindSmallest();
	void DownHeap(int );
	void UpHeap(int );
	bool Insert(int, double);
	void move_forward(int );
	int get_parent(int);
	int get_leftchild(int);
	int get_rightchild(int);
	bool extend(int step=100);
	void reset();
	bool is_empty();
	int get_pos(int vert);
};


class Trajectory{
public:
	int index;
	int  nlinesegs;
	int  curMaxNumLinesegs;
	LineSeg *linesegs;
	unsigned char roadtype;  /*record the road type for street modeling*/
	bool closed;
	bool is_mapboundary;     /*  record whether it is the boundaries of a loaded map  */

	int saddleID;           /*which saddle this trajectory belongs to*/

	double eulerstep_scalar;

	double traj_len;

	/*Construct the trajectory*/
	Trajectory(int index, int curMaxNum = 200);

	~Trajectory()
	{
		if(curMaxNumLinesegs > 0)
		{
			free(linesegs);
			curMaxNumLinesegs=0;
		}
	}

	bool store_to_global_line_segs(CurvePoints *temp, int num);

	/*extend the line segment list if there is not enough space left*/
	bool extend_line_segments(int add_size);

	//get the length of the trajectory
	double get_length();

	//get any line segment according to the input index
	LineSeg *get_line_seg(int index);

	//remove the front n line segments
	bool remove_front_nlines(int n);

	//add n new line segments in the front
	bool add_front_nlines(LineSeg *, int);

	//remove the last n line segments
	bool remove_last_nlines(int n);

	//add n new line segments at the end
	bool add_last_nlines(LineSeg *, int);

	//copy line segment
	void copy_linesegs(LineSeg *, int);

	//reverse the trajectory
	bool reverse_lines();

	int trace_in_quad(int &face_id, double globalp[2], int type, int &flag);
}; //end of Trajectory class






class TrajectoryList{
public:
	Trajectory **trajs;         //the trajectory list
	int ntrajs;                  //current number of existing trajectories
	int curMaxNumTrajs;          //maximum number of trajectories can be stored
	double length;                //the flow length of the trajectory


	// The similar list operations
	TrajectoryList(int initsize = 1000) //construction
	{
		trajs = (Trajectory **)malloc(sizeof(Trajectory *)*initsize);
		curMaxNumTrajs = initsize;
		ntrajs = 0;

		if(trajs == NULL)
		{
			exit(-1);
		}
		
		for(int i = 0; i < initsize; i++)
			trajs[i] = NULL;
		curMaxNumTrajs = initsize;
	} 

	~TrajectoryList()
	{
		int i, j;

		for(i = 0; i < curMaxNumTrajs; i++)
		{
			if(trajs[i] != NULL)
			{
				free(trajs[i]->linesegs);
			}
		}

		free(trajs);
	}

	//add a new vertex to the end of the list, if it succeeds, return true
	inline bool append(Trajectory *s)
	{
		if(isFull ())
			if(!extend(100))
				return false;             //if not enough memory available, return false
		trajs[ntrajs] = s;
		//copyElem(s, polist[nporbits]);
		s->index=ntrajs;
		ntrajs++;
		return true;
	} 

	inline bool del_End() //delete the vertex at the end of the list
	{
		if(isEmpty())  return false;
		ntrajs --;
		return true;
	} 

	inline void copy_Elem(Trajectory *s, Trajectory *d)
	{
	}

	//delete the corresponding  vertex, if it succeeds, return true
	inline bool del_Node(Trajectory *s) 
	{
		if(isEmpty())  return false;

		//find the vertex, if find it, delete and move the following vertices forward
		//otherwise, return false;

		int i, pos = -1;

		for(i = 0; i < ntrajs; i++)
		{
			if(trajs[i] == s)
			{
				pos = i;
				break;
			}
		}

		if(pos == -1) return false;

		//delete it
		for(i = pos; i < ntrajs-1; i++)
		{
			//we need a copy function
			copy_Elem(trajs[i], trajs[i+1]);
		}

		ntrajs--;

		return true;
	} 

	inline bool isEmpty()  //judge whether the list is empty
		{
		if(ntrajs == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(ntrajs == curMaxNumTrajs) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		Trajectory **temp = trajs;
		trajs = (Trajectory **) malloc(sizeof(Trajectory *) * (curMaxNumTrajs + step));
		if( trajs == NULL)
		{
			//fail
			//char rout[256], var[256];
			//sprintf(rout, "%s", "TrajectoryList::extend");
			//sprintf(var, "%s", "trajs");

			//write_mem_error(rout, var, 1);
			curMaxNumTrajs = 0;
			trajs = temp;
			exit(-1);

			return false;
		}

		int i;

		for(i = 0; i < curMaxNumTrajs; i++)
			trajs[i] = temp[i];
		for(i = curMaxNumTrajs; i < curMaxNumTrajs+step; i++)
			trajs[i] = NULL;

		curMaxNumTrajs += step;

		free(temp);
		return true;
	}

	inline void reset()
	{
		ntrajs = 0;
	}

	/*we now put the separatrix calculation here */
	//void cal_startpt_sep(int triangleID, icVector3 sep_vector, double saddle_ce[3], double newpos[3]);
	//void cal_separatrices();

}; //end of TrajectoryList class



typedef struct Intersection{
	double gpos[2];  /*the global coordinates of the intersection*/
	int cellid;      /*which cell it locates at*/
	int index;       /*the index of the intersection*/
	int majorline_id, minorline_id; /*which major and minor tensor lines it locates on*/
	int majlineseg, minlineseg;    /*which line segments it locates at*/
	int *adj_edges;      /*the other intersects that directly connect with it*/
	unsigned char nadjedges;   /*the number of its connected intersections*/
	bool endpt;      /*false--intersection; true--end point*/
	/*the four visual points for road visualization*/
	double majpt1[2], majpt2[2], majpt3[2], majpt4[2];
	//unsigned char minpt_order; /*0--*/
	bool inside_region;
	bool deleted;

	/*   whether the intersection is a major road intersects a minor road
	     or two major (or minor) road intersects
		 0 -- a major road intersects a minor road
		 1 -- two major roads intersect
		 2 -- two minor roads intersect
	*/
	unsigned char intersect_type;
}Intersection;

class IntersectionList{
public:
	Intersection **intersects;
	int nelems;
	int curMaxNum;

	IntersectionList(int init_size = 100)
	{
		if(init_size == 0)
		{
			intersects = NULL;
			curMaxNum = nelems = 0;
			return;
		}

		intersects = (Intersection **)malloc(sizeof(Intersection*)*init_size);
		if(intersects == NULL)
			exit(-1);
		curMaxNum = init_size;
		for(int i=0; i<curMaxNum; i++)
		{
			intersects[i]=(Intersection *)malloc(sizeof(Intersection));
		}
		nelems = 0;
	}

	~IntersectionList()
	{
		if(intersects != NULL)
		{
			int i;
			for(i=0; i<curMaxNum; i++)
			{
				if(intersects[i] != NULL)
				{
					free(intersects[i]);
					intersects[i]=NULL;
				}
			}
			free(intersects);
			intersects=NULL;
		}
	}

	bool isfull()
	{
		if(nelems >= curMaxNum) return true;
		return false;
	}

	bool extend(int step=20)
	{
		Intersection **temp = intersects;
		intersects=(Intersection**)malloc(sizeof(Intersection*)*(curMaxNum+step));
		if(intersects == NULL)
			return false;

		int i;
		for(i=0; i<curMaxNum; i++)
			intersects[i]=temp[i];
		for(i=curMaxNum; i<curMaxNum+step; i++)
			intersects[i]=(Intersection *)malloc(sizeof(Intersection));
		curMaxNum += step;
		free(temp);
		return true;

	}

	void addNew(Intersection *intersect)
	{
		if(nelems>=curMaxNum)
		{
			if(!extend(50))
				exit(-1);
		}

		intersects[nelems] = intersect;
		intersect->index=nelems;   /*crucial: remember to set the index for each intersection*/
		nelems++;
	}

	/*other functions, such as intersection computation*/
};


/* Graph edge */
typedef struct Graph_Edge{
public:
	bool cancel;
	bool visited;
	int edge_index;                   //Edge index (unique)
	int node_index1, node_index2;     //Indices of the two nodes this edge connects
	int index;
	unsigned char roadtype;
}Graph_Edge;





class Graph_EdgeList{
public:
	Graph_Edge **edges;
	int nedges;
	int curMaxNumGedges;

	// The similar list operations as VertexList
	// The similar list operations
	Graph_EdgeList(int initsize = 1000) //construction
	{
		edges = (Graph_Edge **)malloc(sizeof(Graph_Edge *)*initsize);
		curMaxNumGedges = initsize;
		nedges = 0;

		if(edges == NULL)
		{
			exit(-1);
		}
		
		for(int i = 0; i < initsize; i++)
			edges[i] = NULL;
	
	} 

	~Graph_EdgeList()
	{
		if(edges != NULL)
		{
			int i;
			for(i = 0; i < curMaxNumGedges; i++)
			{
				if(edges[i]!=NULL)
					free(edges[i]);
			}
			free(edges);
		}
	}

	//add a new vertex to the end of the list, if it succeeds, return true
	inline bool append(Graph_Edge *s)
	{
		if(isFull ())
			if(!extend(100))
				return false;             //if not enough memory available, return false
		edges[nedges] = s;
		edges[nedges]->index = nedges;    //crucial to indexing the edge
		//copyElem(s, polist[nporbits]);
		nedges++;
		return true;
	} 

	inline bool del_End() //delete the vertex at the end of the list
	{
		if(isEmpty())  return false;
		nedges --;
		return true;
	} 

	inline void copy_Elem(Graph_Edge *s, Graph_Edge *d)
	{
	}

	//delete the corresponding  vertex, if it succeeds, return true
	inline bool del_Node(Graph_Edge *s) 
	{
		if(isEmpty())  return false;

		//find the vertex, if find it, delete and move the following vertices forward
		//otherwise, return false;

		int i, pos = -1;

		for(i = 0; i < nedges; i++)
		{
			if(edges[i] == s)
			{
				pos = i;
				break;
			}
		}

		if(pos == -1) return false;

		//delete it
		for(i = pos; i < nedges-1; i++)
		{
			//we need a copy function
			copy_Elem(edges[i], edges[i+1]);
		}

		nedges--;

		return true;
	} 

	inline bool isEmpty()  //judge whether the list is empty
		{
		if(nedges == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(nedges == curMaxNumGedges) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		Graph_Edge **temp = edges;
		edges = (Graph_Edge **) malloc(sizeof(Graph_Edge *) * (curMaxNumGedges + step));
		if( edges == NULL)
		{
			//fail
			//exit(-1);

			edges = temp;
			return false;
		}

		int i;
		for(i = 0; i < curMaxNumGedges; i++)
			edges[i] = temp[i];
		
		for(i = curMaxNumGedges; i < curMaxNumGedges+step; i++)
			edges[i] = NULL;

		curMaxNumGedges += step;

		free(temp);
		return true;
	}

	inline void reset()
	{
		nedges = 0;
	}

}; //end of Graph_EdgeList class



/* Street graph edge */
typedef struct StreetGraphEdge{
public:
	bool cancel;
	bool visited;
	//int edge_index;                   //Edge index (unique)
	int node_index1, node_index2;     //Indices of the two nodes this edge connects
	int index;

	/*We store a list of intermediate points for this edge of the street map 11/06/2007*/
	Point **inter_pts;
	int ninter_pts;
	
	unsigned char roadtype;

	/*for block operations*/
	int regionblocks[2];  //an edge can only be shared by at most two region blocks
	unsigned char nregionblocks;
	bool inversed_orient;     //to record which orientation has been used in the block construction
	                          //false-- node_index1->node_index2;  true--node_index1<-node_index2

	//bool deleted;
}StreetGraphEdge;





class StreetGraphEdge_List{
public:
	StreetGraphEdge **edges;
	int nedges;
	int curMaxNumGedges;

	// The similar list operations as VertexList
	// The similar list operations
	StreetGraphEdge_List(int initsize = 1000) //construction
	{
		edges = (StreetGraphEdge **)malloc(sizeof(StreetGraphEdge *)*initsize);
		curMaxNumGedges = initsize;
		nedges = 0;

		if(edges == NULL)
		{
			exit(-1);
		}
		
		for(int i = 0; i < initsize; i++)
			edges[i] = NULL;
	
	} 

	~StreetGraphEdge_List()
	{
		if(edges != NULL)
		{
			int i;
			for(i = 0; i < curMaxNumGedges; i++)
			{
				if(edges[i]->inter_pts!=NULL)
				{
					int j;
					for(j=0; j<edges[i]->ninter_pts; j++)
					{
						if(edges[i]->inter_pts[j]!=NULL)
						{
							free(edges[i]->inter_pts[j]);
							edges[i]->inter_pts=NULL;
						}
					}
					free(edges[i]->inter_pts);
					edges[i]->inter_pts=NULL;
				}
				if(edges[i]!=NULL)
				{
					free(edges[i]);
					edges[i]=NULL;
				}
			}
			free(edges);
			edges=NULL;
		}
	}

	//add a new vertex to the end of the list, if it succeeds, return true
	inline bool append(StreetGraphEdge *s)
	{
		if(isFull ())
			if(!extend(100))
				return false;             //if not enough memory available, return false
		edges[nedges] = s;
		edges[nedges]->index = nedges;    //crucial to indexing the edge
		//copyElem(s, polist[nporbits]);
		nedges++;
		return true;
	} 

	inline bool del_End() //delete the vertex at the end of the list
	{
		if(isEmpty())  return false;
		nedges --;
		return true;
	} 

	inline void copy_Elem(StreetGraphEdge *s, StreetGraphEdge *d)
	{
	}

	//delete the corresponding  vertex, if it succeeds, return true
	inline bool del_Node(StreetGraphEdge *s) 
	{
		if(isEmpty())  return false;

		//find the vertex, if find it, delete and move the following vertices forward
		//otherwise, return false;

		int i, pos = -1;

		for(i = 0; i < nedges; i++)
		{
			if(edges[i] == s)
			{
				pos = i;
				break;
			}
		}

		if(pos == -1) return false;

		//delete it
		for(i = pos; i < nedges-1; i++)
		{
			//we need a copy function
			copy_Elem(edges[i], edges[i+1]);
		}

		nedges--;

		return true;
	} 

	inline bool isEmpty()  //judge whether the list is empty
		{
		if(nedges == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(nedges == curMaxNumGedges) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		StreetGraphEdge **temp = edges;
		edges = (StreetGraphEdge **) malloc(sizeof(StreetGraphEdge *) * (curMaxNumGedges + step));
		if( edges == NULL)
		{
			//fail
			//exit(-1);

			edges = temp;
			return false;
		}

		int i;
		for(i = 0; i < curMaxNumGedges; i++)
			edges[i] = temp[i];
		
		for(i = curMaxNumGedges; i < curMaxNumGedges+step; i++)
			edges[i] = NULL;

		curMaxNumGedges += step;

		free(temp);
		return true;
	}

	inline void reset()
	{
		nedges = 0;
	}

}; //end of StreetGraphEdge_List class




class StreetNet{
public:
	IntersectionList *nodelist;
	//IntersectionList *danglepts;   /*the end points of tensor lines*/
	
	//Graph_EdgeList *edgelist;
	StreetGraphEdge_List *edgelist;

	StreetNet(int nnodes = 0, int nedges = 0)
	{
		if(nnodes == 0)
		{
			nodelist=NULL;
			//danglepts=NULL;
		}
		else
		{
			nodelist=new IntersectionList(nnodes);
			//danglepts=new IntersectionList((int)(nnodes/4));
		}
		if(nedges == 0)
			edgelist=NULL;
		else
		{
			//edgelist=new Graph_EdgeList(nedges);
			edgelist=new StreetGraphEdge_List(nedges);
		}
	}

	~StreetNet()
	{
		if(nodelist != NULL)
		{
			delete nodelist;
			nodelist=NULL;
		}
		if(edgelist != NULL)
		{
			delete edgelist;
			edgelist=NULL;
		}
		//if(danglepts != NULL)
		//	delete danglepts;
	}

	void reset_streetnet()
	{
		int i;
		if(nodelist->curMaxNum!=0)
		{
			for(i=0;i<nodelist->curMaxNum;i++)
			{
				if(nodelist->intersects[i]!=NULL)
				{
					free(nodelist->intersects[i]);
					nodelist->intersects[i]=NULL;
				}
			}
		}
		if(edgelist->curMaxNumGedges!=0)
		{
			for(i=0;i<edgelist->curMaxNumGedges;i++)
			{
				//if(edgelist->edges[i]!=NULL)
				//{
				//	free(edgelist->edges[i]);
				//	edgelist->edges[i]=NULL;
				//}
				if(edgelist->edges[i] != NULL)
				{
					if(edgelist->edges[i]->inter_pts!=NULL)
					{
						int j;
						for(j=0; j<edgelist->edges[i]->ninter_pts; j++)
						{
							if(edgelist->edges[i]->inter_pts[j]!=NULL)
							{
								free(edgelist->edges[i]->inter_pts[j]);
								edgelist->edges[i]->inter_pts[j]=NULL;
							}
						}
						free(edgelist->edges[i]->inter_pts);
						edgelist->edges[i]->inter_pts=NULL;
					}
					free(edgelist->edges[i]);
					edgelist->edges[i]=NULL;
				}
			}
		}
		nodelist->nelems=0;
		edgelist->nedges=0;
	}

	void add_edge_to_node(int node, int edgeindex);
	/*other graph operations, such as searching/traverse, editing, deformation
	topology computation*/
};


typedef struct IntersectionInfo{
	int intersect_id;  /*through this reference id, we can find out the corresponding
					   graph information, such as the degree of the intersection*/
	int lineseg_id;
}IntersectionInfo;

class TensorLineIntersectionInfoList{
public: 
	IntersectionInfo **infolist;
	int nelems;
	int curMaxNum;
	bool majormin;          /* false--major direction;  true--minor direction */
	int trajID;             // the corresponding trajectory index

	TensorLineIntersectionInfoList(bool majormin=false, int trajID=-1, int init_size = 20)
	{
		if(init_size == 0)
		{
			infolist = NULL;
			nelems = 0;
			curMaxNum = 0;
			this->majormin=majormin;
			this->trajID=trajID;
			return;
		}

		infolist = (IntersectionInfo **)malloc(sizeof(IntersectionInfo*)*init_size);
		if(infolist == NULL) exit(-1);
		int i;
		for(i=0; i<init_size; i++)
			infolist[i]=NULL;
		curMaxNum = init_size;
		nelems = 0;
		this->majormin=majormin;
		this->trajID=trajID;
	}

	~TensorLineIntersectionInfoList()
	{
		if(infolist != NULL)
		{
			int i;
			for(i=0; i<curMaxNum; i++)
			{
				if(infolist[i]!=NULL)
					free(infolist[i]);
			}
			free(infolist);
		}
	}

	void reset()
	{
		nelems = 0;
	}

	bool isfull()
	{
		if(nelems>=curMaxNum) return true;
		return false;
	}

	bool extend(int step=20)
	{
		IntersectionInfo **temp=infolist;
		infolist=(IntersectionInfo**)malloc(sizeof(IntersectionInfo*)*(curMaxNum+step));
		if(infolist == NULL) return false;
		int i;
		for(i=0; i<curMaxNum; i++)
			infolist[i]=temp[i];
		for(i=curMaxNum; i<curMaxNum+step; i++)
			infolist[i] = (IntersectionInfo*)malloc(sizeof(IntersectionInfo));
		curMaxNum+=step;
		free(temp);
		return true;
	}

	int cal_pos_on_same_lineseg(IntersectionInfo *newinfo, int start, int end); //implemented in evenstreamlines.cpp

	void sorted_add(IntersectionInfo *newinfo)
	{
		if(nelems>=curMaxNum)
		{
			if(!extend())
				exit(-1);
		}

		/*    sorted insert!   
		      bug: if two intersections are located at one line segment, could be wrong!
		*/
		int i, pos=0;
		int pos_end=0;
		bool found=false;
		for(i=0; i<nelems; i++)
		{
			if(!found && newinfo->lineseg_id<infolist[i]->lineseg_id)
			{
				pos = i;
				break;
			}
			else if(newinfo->lineseg_id==infolist[i]->lineseg_id)
			{
				/* we have now at least two intersections being located
				   at the same line segment
				*/
				if(!found)
				{
					pos_end=pos=i;
					found=true;
				}
				else
				{
					pos_end=i;
				}
				//break;
			}
			else
			{
				if(found && newinfo->lineseg_id<infolist[i]->lineseg_id)
				{
					break;
				}
			}
		}

		if(found)
		{
			int temppos=cal_pos_on_same_lineseg(newinfo, pos, pos_end);
			if(temppos>=0)
				pos=temppos;
		}

		if(i>=nelems) pos=nelems;
		for(i=nelems; i>pos; i--)
			infolist[i]=infolist[i-1];
		infolist[pos]=newinfo;

		nelems ++;
		return ;
	}
};




///////////////////////////////////////////////////
enum road_type{
	MINOR,
	MAJOR,
	HIGHWAY,
	BOUNDARY
	//FREEWAY
	//PLAIN,
};


typedef struct RoadLineSeg{
	double start[2], end[2];
	int cell_id;
}RoadLineSeg;

class OneRoadVis
{
public:
	RoadLineSeg **roadline1, **roadline2;
	int nroadlines1, nroadlines2;
	int curMaxNum1, curMaxNum2;

	OneRoadVis(int init_size = 100)
	{
		if(init_size == 0)
		{
			roadline1 = roadline2 = NULL;
			nroadlines1=nroadlines2=curMaxNum1=curMaxNum2 = 0;
			return;
		}

		roadline1=(RoadLineSeg**)malloc(sizeof(RoadLineSeg*)*init_size);
		if(roadline1==NULL) exit(-1);
		roadline2=(RoadLineSeg**)malloc(sizeof(RoadLineSeg*)*init_size);
		if(roadline2==NULL) exit(-1);

		int i;
		for(i=0; i<init_size; i++)
		{
			//roadline1[i]=(RoadLineSeg*)malloc(sizeof(RoadLineSeg));
			//roadline2[i]=(RoadLineSeg*)malloc(sizeof(RoadLineSeg));
			roadline1[i]=NULL;
			roadline2[i]=NULL;
		}
			
		nroadlines1=nroadlines2 = 0;
		curMaxNum1=curMaxNum2 = init_size;
	}

	~OneRoadVis()
	{
		int i;
		if(roadline1 != NULL)
		{
			for(i=0; i<curMaxNum1; i++)
			{
				if(roadline1[i] != NULL)
					free(roadline1[i]);
			}
			free(roadline1);
		}
		if(roadline2 != NULL)
		{
			for(i=0; i<curMaxNum2; i++)
			{
				if(roadline2[i] != NULL)
					free(roadline2[i]);
			}
			free(roadline2);
		}
	}

	void addNew(RoadLineSeg *newline, bool flag)
	{
		//FILE *fp;
		if(isfull(flag))
		{
			if(!extend(50, flag))
			{
				//fp=fopen("extendfail.txt", "w");
				exit(-1);
			}
		}

		if(!flag) /*to road 1*/
		{
			roadline1[nroadlines1]=newline;
			nroadlines1++;
		}
		else /*to road 2*/
		{
			roadline2[nroadlines2]=newline;
			nroadlines2++;
		}
	}

	bool isfull(bool flag)
	{
		if(!flag)
		{
			if(nroadlines1 >= curMaxNum1) return true;
			return false;
		}
		else
		{
			if(nroadlines2 >= curMaxNum2) return true;
			return false;
		}
	}

	bool extend(int step = 50, bool flag = false)
	{
		int i;
		if(!flag)
		{
			/*extend line 1*/
			RoadLineSeg **temp = roadline1;
			roadline1=(RoadLineSeg**)malloc(sizeof(RoadLineSeg*)*(curMaxNum1+step));
			if(roadline1 == NULL)
				return false;
			for(i=0; i<curMaxNum1; i++)
				roadline1[i]=temp[i];
			for(i=curMaxNum1; i<curMaxNum1+step; i++)
				roadline1[i]=NULL;
			free(temp);
			curMaxNum1+=step;
			return true;
		}
		else
		{
			/*extend line 2*/
			RoadLineSeg **temp = roadline2;
			roadline2=(RoadLineSeg**)malloc(sizeof(RoadLineSeg*)*(curMaxNum2+step));
			if(roadline2 == NULL)
				return false;
			for(i=0; i<curMaxNum2; i++)
				roadline2[i]=temp[i];
			for(i=curMaxNum2; i<curMaxNum2+step; i++)
				roadline2[i]=NULL;
			free(temp);
			curMaxNum2+=step;
			return true;
		}
	}

};


class StreetVis
{
public:
	OneRoadVis **majDirRoads, **minDirRoads;
	int nmajDirRoads, nminDirRoads;

	StreetVis(int init_maj=50, int init_min=50)
	{
		/*major*/
		int i;
		if(init_maj==0)
		{
			majDirRoads=NULL;
			nmajDirRoads = 0;
		}
		else
		{
			majDirRoads=new OneRoadVis *[init_maj];
			if(majDirRoads == NULL) exit(-1);
			nmajDirRoads = init_maj;
			for(i=0; i<init_maj; i++)
				majDirRoads[i]=NULL;
		}

		/*minor*/
		if(init_maj==0)
		{
			minDirRoads=NULL;
			nminDirRoads = 0;
		}
		else
		{
			minDirRoads=new OneRoadVis *[init_min];
			if(minDirRoads == NULL) exit(-1);
			nminDirRoads = init_min;
			for(i=0; i<init_min; i++)
				minDirRoads[i]=NULL;
		}
	}

	~StreetVis()
	{
		int i;
		if(majDirRoads != NULL)
		{
			for(i=0; i<nmajDirRoads; i++)
			{
				if(majDirRoads[i]!=NULL)
					delete majDirRoads[i];
			}
			delete [] majDirRoads;
		}
		if(minDirRoads != NULL)
		{
			for(i=0; i<nminDirRoads; i++)
			{
				if(minDirRoads[i]!=NULL)
					delete minDirRoads[i];
			}
			delete [] minDirRoads;
		}
	}

	void init_lineseglist();
	void construct_roads_vis(TrajectoryList *major, TrajectoryList *minor);
	void ave_normal(double l1[2], double l2[2], double ave_n[2]);
};




class Parcel
{
public:
	int *edges;
	int nedges;

	int *intersects;
	int nintersects;

	Point **parcel_vertices;
	int nparcel_vertices;

	Parcel()
	{
		edges = NULL;
		intersects = NULL;
		parcel_vertices = NULL;
		nedges = 0;
		nintersects = 0;
		nparcel_vertices = 0;
	}

	~Parcel()
	{
		if(edges != NULL)
			free(edges);
		if(intersects != NULL)
			free(intersects);
		if(parcel_vertices != NULL)
			free(parcel_vertices);
	}
};

class ParcelList
{
public:
	Parcel **elems;
	int nelems;
	int curMaxNumElems;

	ParcelList(int initsize = 100) //construction
	{
		elems = (Parcel **)malloc(sizeof(Parcel *)*initsize);
		curMaxNumElems = initsize;
		nelems = 0;

		if(elems == NULL)
		{
			exit(-1);
		}
		
		for(int i = 0; i < initsize; i++)
			elems[i] = NULL;
	
	} 

	~ParcelList()
	{
		if(elems != NULL)
		{
			int i;
			for(i = 0; i < curMaxNumElems; i++)
			{
				if(elems[i]!=NULL)
					delete elems[i];
			}
			free(elems);
		}
	}

	//add a new vertex to the end of the list, if it succeeds, return true
	inline bool append(Parcel *s)
	{
		if(isFull ())
			if(!extend(100))
				return false;             //if not enough memory available, return false
		elems[nelems] = s;

		nelems++;
		return true;
	} 

	inline bool del_End() //delete the vertex at the end of the list
	{
		if(isEmpty())  return false;
		nelems --;
		return true;
	} 

	inline void copy_Elem(Parcel *s, Parcel *d)
	{
	}

	//delete the corresponding  vertex, if it succeeds, return true
	inline bool del_Node(Parcel *s) 
	{
		if(isEmpty())  return false;

		//find the vertex, if find it, delete and move the following vertices forward
		//otherwise, return false;

		int i, pos = -1;

		for(i = 0; i < nelems; i++)
		{
			if(elems[i] == s)
			{
				pos = i;
				break;
			}
		}

		if(pos == -1) return false;

		//delete it
		for(i = pos; i < nelems-1; i++)
		{
			//we need a copy function
			copy_Elem(elems[i], elems[i+1]);
		}

		nelems--;

		return true;
	} 

	inline bool isEmpty()  //judge whether the list is empty
	{
		if(nelems == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(nelems == curMaxNumElems) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		Parcel **temp = elems;
		elems = (Parcel **) malloc(sizeof(Parcel *) * (curMaxNumElems + step));
		if( elems == NULL)
		{
			//fail
			//exit(-1);

			elems = temp;
			return false;
		}

		int i;
		for(i = 0; i < curMaxNumElems; i++)
			elems[i] = temp[i];
		
		for(i = curMaxNumElems; i < curMaxNumElems+step; i++)
			elems[i] = NULL;

		curMaxNumElems += step;

		free(temp);
		return true;
	}

	inline void reset()
	{
		nelems = 0;
	}
};


/*
Data structure for blocks
*/
typedef struct RegionBlock{
	int index;
	int *nodelist;
	unsigned char num_nodes;
	int *edgelist;
	unsigned char num_edges;
	bool closed_curved;
	int trajid;               /* the index of the corresponding sketch curve */
}RegionBlock;

class RegionBlockList
{
public:
	RegionBlock **blocks;
	int nelems;
	int curMaxNum;

	/*constructor*/
	RegionBlockList(int initsize=50)
	{
		if(initsize==0)
		{
			blocks=NULL;
			curMaxNum=nelems=0;
			return;
		}

		blocks=(RegionBlock**)malloc(sizeof(RegionBlock*)*initsize);
		if(blocks==NULL)
			exit(-1);
		int i;
		for(i=0;i<initsize;i++)
			blocks[i]=NULL;
		curMaxNum=initsize;
		nelems=0;
	}
	~RegionBlockList()
	{
		int i;
		if(blocks!=NULL)
		{
			for(i=0;i<curMaxNum;i++)
			{
				if(blocks[i]!=NULL)
				{
					if(blocks[i]->nodelist!=NULL)
						free(blocks[i]->nodelist);
					if(blocks[i]->edgelist!=NULL)
						free(blocks[i]->edgelist);
					free(blocks[i]);
				}
			}
			free(blocks);
		}
	}

	/*List operation*/
	//add a new vertex to the end of the list, if it succeeds, return true
	inline bool append(RegionBlock *s)
	{
		if(isFull ())
			if(!extend(100))
				return false;             //if not enough memory available, return false
		blocks[nelems] = s;
		blocks[nelems]->index=nelems;
		nelems++;
		return true;
	} 

	inline bool del_End() //delete the vertex at the end of the list
	{
		if(isEmpty())  return false;
		nelems --;
		return true;
	} 

	inline void copy_Elem(RegionBlock *s, RegionBlock *d)
	{
	//int index;
	//int *nodelist;
	//unsigned char num_nodes;
	//int *edgelist;
	//unsigned char num_edges;
	}

	//delete the corresponding  vertex, if it succeeds, return true
	inline bool del_Node(RegionBlock *s) 
	{
		if(isEmpty())  return false;

		//find the vertex, if find it, delete and move the following vertices forward
		//otherwise, return false;

		int i, pos = -1;

		for(i = 0; i < nelems; i++)
		{
			if(blocks[i] == s)
			{
				pos = i;
				break;
			}
		}

		if(pos == -1) return false;

		//delete it
		for(i = pos; i < nelems-1; i++)
		{
			//we need a copy function
			copy_Elem(blocks[i], blocks[i+1]);
		}

		nelems--;

		return true;
	} 

	inline bool isEmpty()  //judge whether the list is empty
	{
		if(nelems == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(nelems >= curMaxNum) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		RegionBlock **temp = blocks;
		blocks = (RegionBlock **) malloc(sizeof(RegionBlock *) * (curMaxNum + step));
		if( blocks == NULL)
		{
			//fail
			//char rout[256], var[256];
			//sprintf(rout, "%s", "SamplePtList::extend");
			//sprintf(var, "%s", "samples");

			//write_mem_error(rout, var, 1);
			curMaxNum  = 0;
			blocks = temp;
			exit(-1);
			return false;
		}

		int i;
		for(i = 0; i < curMaxNum; i++)
			blocks[i] = temp[i];
		for(i = curMaxNum; i < curMaxNum+step; i++)
			blocks[i] = NULL;
		curMaxNum += step;

		free(temp);
		return true;
	}

	inline void reset()
	{
		nelems = 0;
	}
};


/*
The class for a brush
*/
class Brush
{
public:
	ctr_point **brushpts;
	int nelems;
	int curMaxNum;

	/*constructor & destructor*/
	Brush(int initsize=400)
	{
		if(initsize==0)
		{
			brushpts=NULL;
			nelems=curMaxNum=0;
			return;
		}
		brushpts=(ctr_point**)malloc(sizeof(ctr_point*)*initsize);
		int i;
		for(i=0;i<initsize;i++)
			brushpts[i]=NULL/*(ctr_point*)malloc(sizeof(ctr_point))*/;
		nelems=0;
		curMaxNum=initsize;
	}

	~Brush()
	{
		int i;
		if(brushpts!=NULL)
		{
			for(i=0;i<curMaxNum;i++)
			{
				if(brushpts[i]!=NULL)
					free(brushpts[i]);
			}
			free(brushpts);
		}
	}

	/*List operations*/
	//add a new vertex to the end of the list, if it succeeds, return true
	inline bool append(ctr_point *s)
	{
		if(isFull ())
			if(!extend(100))
				return false;             //if not enough memory available, return false
		brushpts[nelems] = s;
		nelems++;
		return true;
	} 

	inline bool del_End() //delete the vertex at the end of the list
	{
		if(isEmpty())  return false;
		nelems --;
		return true;
	} 

	inline void copy_Elem(ctr_point *s, ctr_point *d)
	{
	//int index;
	//int *nodelist;
	//unsigned char num_nodes;
	//int *edgelist;
	//unsigned char num_edges;
	}


	inline bool isEmpty()  //judge whether the list is empty
	{
		if(nelems == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(nelems >= curMaxNum) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		ctr_point **temp = brushpts;
		brushpts = (ctr_point **) malloc(sizeof(ctr_point *) * (curMaxNum + step));
		if( brushpts == NULL)
		{
			//write_mem_error(rout, var, 1);
			curMaxNum  = 0;
			brushpts = temp;
			exit(-1);
			return false;
		}

		int i;
		for(i = 0; i < curMaxNum; i++)
			brushpts[i] = temp[i];
		for(i = curMaxNum; i < curMaxNum+step; i++)
			brushpts[i] = NULL;
		curMaxNum += step;

		free(temp);
		return true;
	}

	inline void reset()
	{
		nelems = 0;
	}
};

/*
Each sketch will be represented as a sequence of points
similar to the brush stroke interface
*/
class SketchList
{
public:
	Brush *brushlist;
	int nbrushes;
	int curMaxNum;

	/*constructor & destructor*/
	SketchList(int initsize=50)
	{
		if(initsize==0)
		{
			brushlist=NULL;
			nbrushes=curMaxNum=0;
			return;
		}
        brushlist=new Brush[initsize];
		nbrushes=0;
		curMaxNum=initsize;
	}

	~SketchList()
	{
		if(brushlist!=NULL)
		{
			delete [] brushlist;
		}
	}

	/*list operations*/
	inline bool isFull()
	{
		if(nbrushes >= curMaxNum) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 10)
	{
		Brush *temp = brushlist;
		brushlist = new Brush [curMaxNum + step];
		if( brushlist == NULL)
		{
			//write_mem_error(rout, var, 1);
			curMaxNum  = 0;
			brushlist = temp;
			exit(-1);
			return false;
		}

		/*we need to copy brushlist*/
		copy_brushlist(temp, curMaxNum);

		curMaxNum += step;

		delete [] temp;
		return true;
	}

	inline void copy_brushlist(Brush *s, int nb)
	{
		int i;
		Brush *onebrush;
		for(i=0;i<nb; i++)
		{
			onebrush=&s[i];
			copy_onebrush(onebrush, &this->brushlist[i]);
		}
	}

	inline void copy_onebrush(Brush *s, Brush *d)
	{
		if(d->curMaxNum<s->nelems)
		{
			if(!d->extend(s->nelems-d->curMaxNum))
				exit(-1);
		}
		int i;
		//for(i=0;i<s->nelems;i++)
		//{
		//	d->brushpts[i]=s->brushpts[i];
		//}


		/*   new solution to avoid possible crush 09/16/2008  */
		for(i=0; i<s->nelems; i++)
		{
			d->brushpts[i]         = (ctr_point*)malloc(sizeof(ctr_point));
			d->brushpts[i]->x      = s->brushpts[i]->x;
			d->brushpts[i]->y      = s->brushpts[i]->y;
			d->brushpts[i]->cellid = s->brushpts[i]->cellid;
		}

		d->nelems=s->nelems;
	}

	inline void reset()
	{
		nbrushes = 0;
	}
};


/*   For scalar field design elements 12/29/2007  */

#define SCALARMAX  40
#define SCALARMIN  -40

typedef struct ScalarSingularElem{
	int index;
	double pos[2];
	unsigned char type;
	bool deleted;

	/*  for editing  */
}ScalarSingularElem;


class ScalarSingularElemList
{
public:
	ScalarSingularElem **scalarsingularelems;
	int nelems;
	int curMaxNum;

	/*   constructor and destructor   */
	ScalarSingularElemList(int initsize=40)
	{
		if(initsize==0)
		{
			scalarsingularelems=0;
			curMaxNum=nelems=0;
			return;
		}

		scalarsingularelems=(ScalarSingularElem**)malloc(sizeof(ScalarSingularElem*)*initsize);
		nelems=0;
		curMaxNum=initsize;

		int i;
		for(i=0;i<initsize;i++)
			scalarsingularelems[i]=NULL;
	}

	~ScalarSingularElemList()
	{
		if(scalarsingularelems!=NULL)
		{
			int i;
			for(i=0;i<curMaxNum;i++)
			{
				if(scalarsingularelems[i]!=NULL)
				{
					free(scalarsingularelems[i]);
					scalarsingularelems[i]=NULL;
				}
			}
			free(scalarsingularelems);
		}
	}

	/*  List operations  */
	//add a new vertex to the end of the list, if it succeeds, return true
	inline bool append(ScalarSingularElem *s)
	{
		if(isFull ())
			if(!extend(100))
				return false;             //if not enough memory available, return false
		scalarsingularelems[nelems] = s;
		s->index=nelems;
		nelems++;
		return true;
	} 

	inline bool del_End() //delete the vertex at the end of the list
	{
		if(isEmpty())  return false;
		nelems --;
		return true;
	} 

	inline void copy_Elem(ScalarSingularElem *s, ScalarSingularElem *d)
	{
	}

	//delete the corresponding  vertex, if it succeeds, return true
	inline bool del_Node(ScalarSingularElem *s) 
	{
		if(isEmpty())  return false;

		//find the vertex, if find it, delete and move the following vertices forward
		//otherwise, return false;

		int i, pos = -1;

		for(i = 0; i < nelems; i++)
		{
			if(scalarsingularelems[i] == s)
			{
				pos = i;
				break;
			}
		}

		if(pos == -1) return false;

		//delete it
		for(i = pos; i < nelems-1; i++)
		{
			//we need a copy function
			copy_Elem(scalarsingularelems[i], scalarsingularelems[i+1]);
		}

		nelems--;

		return true;
	} 

	inline bool isEmpty()  //judge whether the list is empty
	{
		if(nelems == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(nelems == curMaxNum) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		ScalarSingularElem **temp = scalarsingularelems;
		scalarsingularelems = (ScalarSingularElem **) malloc(sizeof(ScalarSingularElem *) 
			* (curMaxNum + step));
		if( temp == NULL)
		{
			//fail
			//char rout[256], var[256];
			//sprintf(rout, "%s", "SamplePtList::extend");
			//sprintf(var, "%s", "samples");

			//write_mem_error(rout, var, 1);
			curMaxNum = 0;
			scalarsingularelems = temp;
			exit(-1);
			return false;
		}

		int i;
		for(i = 0; i < curMaxNum; i++)
			scalarsingularelems[i] = temp[i];
		for(i = curMaxNum; i < curMaxNum+step; i++)
			scalarsingularelems[i] = NULL;
		curMaxNum += step;

		free(temp);
		return true;
	}

	inline void reset()
	{
		int i;
		for(i=0;i<nelems;i++)
		{
			if(scalarsingularelems[i]!=NULL)
			{
				free(scalarsingularelems[i]);
				scalarsingularelems[i]=NULL;
			}
		}
		nelems = 0;
	}
};

/**/
bool is_repeated_elem(int *a, int b, int num);

