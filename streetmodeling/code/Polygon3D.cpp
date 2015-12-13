#include "stdafx.h"
#include ".\polygon3d.h"
#include "numerical.h"

char *elem_names[] = { /* list of the kinds of elements in the user's object */
  "vertex", "face"
};

//PlyProperty vert_props[] = { /* list of property information for a vertex */
//  {"x", Float32, Float32, offsetof(Vertex,x), 0, 0, 0, 0},
//  {"y", Float32, Float32, offsetof(Vertex,y), 0, 0, 0, 0},
//  {"z", Float32, Float32, offsetof(Vertex,z), 0, 0, 0, 0},
//  {"nx", Float32, Float32, offsetof(Vertex,nx), 0, 0, 0, 0},
//  {"ny", Float32, Float32, offsetof(Vertex,ny), 0, 0, 0, 0},
//  {"nz", Float32, Float32, offsetof(Vertex,nz), 0, 0, 0, 0},
//  {"prob_on_path", Float32, Float32, offsetof(Vertex,prob_on_path), 0, 0, 0, 0},
//};
PlyProperty vert_props[] = { /* list of property information for a vertex */
  {"x", Float64, Float64, offsetof(Vertex,x), 0, 0, 0, 0},
  {"y", Float64, Float64, offsetof(Vertex,y), 0, 0, 0, 0},
  {"z", Float64, Float64, offsetof(Vertex,z), 0, 0, 0, 0},
  {"nx", Float64, Float64, offsetof(Vertex,nx), 0, 0, 0, 0},
  {"ny", Float64, Float64, offsetof(Vertex,ny), 0, 0, 0, 0},
  {"nz", Float64, Float64, offsetof(Vertex,nz), 0, 0, 0, 0},
  {"prob_on_path", Float64, Float64, offsetof(Vertex,prob_on_path), 0, 0, 0, 0},
};

PlyProperty face_props[] = { /* list of property information for a face */
  {"vertex_indices", Int32, Int32, offsetof(Face,verts),
   1, Uint8, Uint8, offsetof(Face,nverts)},
  {"area", Float32, Float32, offsetof(Face,area), 0, 0, 0, 0},
};

/*------------------------------------------------------------*/
////Global variables, temporary solution at this moment (06/23/05)

int MaxNumSingularElems;                 //Maximum number of singular elements
int MaxNumRegularElems;                  //Maximum number of regular elements
int MaxNumSingularities;                 //Maximum number of being captured singularities
int MaxNumTrajectories;                  //Maximum number of possible trajectories
int MaxNumSeparatrices;                   //Maximum number of group of separatrices
                                         //(it should be flexible for future pen-and-ink sketch)
int MaxNumLinesegsPerTraj;               //Maximum number of line segments for each trajectory
int MaxNumLimitCycles;

SingularElement *singularelem;          //Singular elememts' list
int cur_singelem_index;

RegularElement *regularelem;            //regular elememts' list
int cur_regelem_index;
int prev_num_reg_elem;                  //for shape design generated limit cycle

Singularities *singularities;           //being captured singularites' list
int cur_singularity_index;

LimitCycle *limitcycles;                 //limit cycles data structure
int cur_limitcycle_index;

//Trajectory *trajectories2;               //new trajectories variables using new data structure

LineSeg **trajectories;                 //trajectories' list
int cur_traj_index;
int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory

Separatrices *separatrices;             //array for group of separatrices
int cur_separatrices_index;
/*------------------------------------------------------------*/

////Define global variables for underneath mesh
Polygon3D Object;

int TotalEdgesNum = 0;                  //Donot know why I have to do this, but I have to

/*------------------------------------------------------------*/

////Back up variable to store the very original field for field reflection operation (under global frame)
icVector2 *backup_field;

icVector2 *pre_field, *post_field;

/*------------------------------------------------------------*/

////Global variables for shape control
//ctr_point *pts;          // allocate our control points array
int resolution = 400;    // how many points in our output array
ctr_point *out_pts;

ctr_point *control_pts;        // allocate our control point array
int num_shapecontrol_pts;
int num_curvepts_output;
int MaxNumShapeControlPts;
int HermiteStep;

icVector2 *CtrPts_tangent;
/*------------------------------------------------------------*/

////Global variables for Conley relation graph 10/14/05
GraphNode *graphnodes = NULL;
GraphEdge *graphedges = NULL;
int cur_node_index;
int cur_graphedge_index;
/*------------------------------------------------------------*/

////Global variables for multiple repellers and attractors selection
int *repeller_nodes;     //the indices of the selected repellers in the conley graph
int NumRepellers;        //the number of the being selected repellers
int *attractor_nodes;    //the indices of the selected attractors in the conley graph
int NumAttractors;       //the number of the being selected attractors

/*------------------------------------------------------------*/
////For connected nodes searching and complex pair cancellation
int maxlength_searcharray = 50;
int *NodeSearchArray = NULL;
int cur_endposition;
int cur_search_node_pos;

int *MediaNodes = NULL;
int MaxMediaNodes = 50;
int Num_MediaNodes = 0;

/*------------------------------------------------------------*/

double triangle_approx_size = 0.;  ////variable for recording the size of the triangle, for regular mesh only!!!!
double avg_edge_length = 0.; //record the average edge length 1/8/06
extern DP htry;

extern void GetCenterofATriangle(int, double []);

/////////////////////////////////////////////////////
///////////////////////////
//////////////////////////////////////////////////////////////////



//////////////////////////////////////

CPolygon3D::CPolygon3D(void)
{
}

CPolygon3D::~CPolygon3D(void)
{
}

Polygon3D CPolygon3D::get_object()
{
	char buffer[128];

	sprintf(buffer, "%s%s.ply", temp_dir, object_name);
    read_object_file(buffer);      //read the information from specific .ply file
    calc_bounding_sphere();
 	PreprocessVertex();            //Set the range of the object between [0,1]x[0,1]
    calc_bounding_sphere();        //deciding the range and the center of the object
	triangulate();                 //ensure the meshe is triangular meshe
	calc_face_normals_and_areas(); //calculate the normals for all faces
	get_vertex_normal();           //calculate the normals for all vertices in the mesh

	for(int i=0; i < Object.nverts; i++)
	{
		//Object.vlist[i]->x = 2.14 * (Object.vlist[i]->x - 0.5) + 0.5;
		//Object.vlist[i]->y = 2.14 * (Object.vlist[i]->y - 0.5) + 0.5;
		//Object.vlist[i]->z = 2.14 * (Object.vlist[i]->z - 0.5) + 0.5;
		Object.vlist[i]->x = 2.5 * (Object.vlist[i]->x - 0.5) + 0.5;
		Object.vlist[i]->y = 2.5 * (Object.vlist[i]->y - 0.5) + 0.5;
		Object.vlist[i]->z = 2.5 * (Object.vlist[i]->z - 0.5) + 0.5;
	}

	///Initialize the number of the edges for Polyhedron and vertex
	//Object.nedges = 0;

	////Get the triangle size here (for triangular surfaces, we need to calculate them for each triangle!!!)
	////11/05/05
	icVector2 temp_vec;
	temp_vec.entry[0] = Object.vlist[Object.flist[0]->verts[0]]->x - Object.vlist[Object.flist[0]->verts[1]]->x;
	temp_vec.entry[1] = Object.vlist[Object.flist[0]->verts[0]]->y - Object.vlist[Object.flist[0]->verts[1]]->y;
	triangle_approx_size = length(temp_vec);

	double zerovec[2] = {0.};

	for(int i = 0; i < Object.nverts; i++){
		Object.vlist[i]->VertID = i;
		Object.vlist[i]->Num_edge = 0;
		Object.vlist[i]->InRegion = 0;
		Object.vlist[i]->OnBoundary = 0;

		Object.vlist[i]->vec.entry[0] = Object.vlist[i]->vec.entry[1] = 0;
	}

	for(int i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->direct_vec[0].set(zerovec);
		Object.flist[i]->direct_vec[1].set(zerovec);
		Object.flist[i]->direct_vec[2].set(zerovec);

		////Reset edge number
		Object.flist[i]->edges[0] = NULL;
		Object.flist[i]->edges[1] = NULL;
		Object.flist[i]->edges[2] = NULL;

		Object.flist[i]->inDesignCellCycle = 0;

		////2/18/06
		Object.flist[i]->trajlist = NULL;
		Object.flist[i]->num_passing_trajs = 0;

		GetCenterofATriangle(i, Object.flist[i]->center);

	}

	////Allocate space for corner table
	///We only consider triangular mesh here
	///we need to avoid to reallocation for it again and again!!!!!!!

	Object.ncorners = Object.nfaces * 3;
	Object.clist = (Corner**) malloc( sizeof(Corner*) * Object.nfaces * 3);

	////Allocate space for backup field in global space
	backup_field = (icVector2 *) malloc(sizeof(icVector2) * Object.nverts);

	pre_field = (icVector2 *) malloc(sizeof(icVector2) * Object.nverts);
	post_field = (icVector2 *) malloc(sizeof(icVector2) * Object.nverts);

	vd_n = 0.0;
	vd_f = Object.center.entry[2] + 3.0 * Object.radius;


	return Object;
}

void CPolygon3D::read_object_file(const char * filename)
{
	unsigned int i,j;
	int elem_count;
	char *elem_name;
	FILE *inFile;

  /*** Read in the original PLY object ***/

	inFile = fopen(filename, "r");
    in_ply = MyPlyLoader.read_ply (inFile);

  for (i = 0; i < in_ply->num_elem_types; i++) {

    /* prepare to read the i'th list of elements */
    elem_name = MyPlyLoader.setup_element_read_ply (in_ply, i, &elem_count);

    if (MyPlyLoader.equal_strings ("vertex", elem_name)) {

      /* create a vertex list to hold all the vertices */
      Object.vlist = (Vertex **) malloc (sizeof (Vertex *) * elem_count*1);
      Object.nverts = elem_count;


      /* set up for getting vertex elements */

      MyPlyLoader.setup_property_ply (in_ply, &vert_props[0]);
      MyPlyLoader.setup_property_ply (in_ply, &vert_props[1]);
      MyPlyLoader.setup_property_ply (in_ply, &vert_props[2]);
      MyPlyLoader.setup_property_ply (in_ply, &vert_props[3]);
      MyPlyLoader.setup_property_ply (in_ply, &vert_props[4]);
      MyPlyLoader.setup_property_ply (in_ply, &vert_props[5]);

	  //Object.vert_other = MyPlyLoader.get_other_properties_ply (in_ply, 
			//		     offsetof(Vertex,other_props));

      /* grab all the vertex elements */
      for (j = 0; j < elem_count; j++) {
        Object.vlist[j] = (Vertex *) malloc (sizeof (Vertex));
        MyPlyLoader.get_element_ply (in_ply, (void *) Object.vlist[j]);
      }
    }
    else if (MyPlyLoader.equal_strings ("face", elem_name)) {

      /* create a list to hold all the face elements */
      Object.flist = (Face **) malloc (sizeof (Face *) * elem_count*1);
      Object.nfaces = elem_count;


      /* set up for getting face elements */

      MyPlyLoader.setup_property_ply (in_ply, &face_props[0]);
      //Object.face_other = MyPlyLoader.get_other_properties_ply (in_ply, 
					 //    offsetof(Face,other_props));

      /* grab all the face elements */
      for (j = 0; j < elem_count; j++) {
        Object.flist[j] = (Face *) malloc (sizeof (Face));
        MyPlyLoader.get_element_ply (in_ply, (void *) Object.flist[j]);
      }
    }
    else
      MyPlyLoader.get_other_element_ply (in_ply);
  }

  MyPlyLoader.close_ply (in_ply);
  fclose(inFile);
}

void CPolygon3D::calc_bounding_sphere(void)
{
	unsigned int i;
	icVector3 min, max;

	for (i=0; i<Object.nverts; i++) {
		if (i==0)  {
				min.set(Object.vlist[i]->x, Object.vlist[i]->y, Object.vlist[i]->z);
				max.set(Object.vlist[i]->x, Object.vlist[i]->y, Object.vlist[i]->z);
		}
		else {
				if (Object.vlist[i]->x < min.entry[0])
				min.entry[0] = Object.vlist[i]->x;
				if (Object.vlist[i]->x > max.entry[0])
				max.entry[0] = Object.vlist[i]->x;
				if (Object.vlist[i]->y < min.entry[1])
				min.entry[1] = Object.vlist[i]->y;
				if (Object.vlist[i]->y > max.entry[1])
				max.entry[1] = Object.vlist[i]->y;
				if (Object.vlist[i]->z < min.entry[2])
				min.entry[2] = Object.vlist[i]->z;
				if (Object.vlist[i]->z > max.entry[2])
				max.entry[2] = Object.vlist[i]->z;
			}
	}
	Object.center = (min + max) * 0.5;
	Object.radius = length(Object.center - min);
	rot_center = Object.center * 1.0;
	char temp[128];

	sprintf(temp, "%s%s_info.txt", temp_dir, object_name);
}

void CPolygon3D::calc_face_normals_and_areas(void)
{
	unsigned int i;
	int j, k, l;
	icVector3 v0, v1, v2;
	Face *face;
	int *verts;
	double temp_x[100], temp_y[100];  // assume on polygon has more than 100 vertices.
	int largest_face;

	Object.area = 0.0;
	largest_face = 0;
	for (i=0; i<Object.nverts; i++)
		Object.vlist[i]->normal.set(0.0);
	for (i=0; i<Object.nfaces; i++){
		Object.flist[i]->area = 0.0;
		face = Object.flist[i];
		if (face->nverts >largest_face)
			largest_face = face->nverts;
		verts = face->verts;		
		v1.set(Object.vlist[verts[0]]->x, Object.vlist[verts[0]]->y, Object.vlist[verts[0]]->z);
		v2.set(Object.vlist[verts[1]]->x, Object.vlist[verts[1]]->y, Object.vlist[verts[1]]->z);
		v0.set(Object.vlist[verts[face->nverts-1]]->x, Object.vlist[verts[face->nverts-1]]->y, Object.vlist[verts[face->nverts-1]]->z);
		Object.flist[i]->normal = cross(v0-v1, v2-v1);
		normalize(Object.flist[i]->normal);
		v1 -= v0;
		v2 = cross(Object.flist[i]->normal, v1);
		normalize(v1);
		normalize(v2);
		for (j=0; j<face->nverts; j++){
			v0.set(Object.vlist[verts[j]]->x-Object.vlist[verts[0]]->x, Object.vlist[verts[j]]->y-Object.vlist[verts[0]]->y, Object.vlist[verts[j]]->z-Object.vlist[verts[0]]->z);
			temp_x[j] = dot(v1, v0);
			temp_y[j] = dot(v2, v0);
		}
		for (j=0; j<face->nverts; j++){
			Object.flist[i]->area += (temp_x[j]*temp_y[(j+1) % face->nverts] - temp_x[(j+1) % face->nverts]*temp_y[j]);
			Object.vlist[face->verts[j]]->normal = Object.vlist[face->verts[j]]->normal + Object.flist[i]->normal;
		}
		Object.flist[i]->area /= 2;
		Object.flist[i]->area = fabs(Object.flist[i]->area);
		Object.area += Object.flist[i]->area;
	}
	for (i=0; i<Object.nverts; i++)
		normalize(Object.vlist[i]->normal);

	double test_val;
	for (j=-2; j<=2; j++) {
		for (k=-2; k<=2; k++) { 
			for (l=-2; l<=2; l++) {
				test_val = 0.0;
				icVector3 test((double)j, (double)k, (double)l);
				test *= Object.radius*100;
				test += Object.center;
				for (i=0; i<Object.nfaces; i++){
					icVector3 cent(Object.vlist[Object.flist[i]->verts[0]]->x, Object.vlist[Object.flist[i]->verts[0]]->y, Object.vlist[Object.flist[i]->verts[0]]->z);

					test_val += dot(test-cent, Object.flist[i]->normal)*Object.flist[i]->area;
				}
				test_val /= Object.area;
				if ((j==0) && (k==0) && (l==0)){
					if (test_val<0) {
						orientation = 0;
					}
					else {
						orientation = 1;
					}
				}
			}
		}
	}
	if (orientation == 1){
		for (i=0; i<Object.nfaces; i++)
			Object.flist[i]->normal *= -1.0;
	}
}

int CPolygon3D::triangulate(void)
{
	unsigned long i, j, l;
	Face *face;
	int *verts;
	long new_faces, first_vertex, last_vertex;
	
	new_faces = Object.nfaces;
	for (i=0; i<Object.nfaces; i++) {
		face = Object.flist[i];
		if (face->nverts>3) {
	    verts = face->verts;
			Object.vlist[Object.nverts] = (Vertex *)malloc(sizeof(Vertex));
			Object.vlist[Object.nverts]->x = 0.0;
			Object.vlist[Object.nverts]->y = 0.0;
			Object.vlist[Object.nverts]->z = 0.0;
			for (j=0; j<face->nverts; j++){ 
				Object.vlist[Object.nverts]->x += Object.vlist[verts[j]]->x;
				Object.vlist[Object.nverts]->y += Object.vlist[verts[j]]->y;
				Object.vlist[Object.nverts]->z += Object.vlist[verts[j]]->z;
			}
			Object.vlist[Object.nverts]->x /= face->nverts;
			Object.vlist[Object.nverts]->y /= face->nverts;
			Object.vlist[Object.nverts]->z /= face->nverts;
			Object.nverts++;

			for (l=0; l<face->nverts-1; l++){
				Object.flist[new_faces+l] = (Face *)malloc(sizeof(Face));
				Object.flist[new_faces+l]->verts = (int *)malloc(sizeof(int) * 3);
				Object.flist[new_faces+l]->verts[0] = face->verts[l];
				Object.flist[new_faces+l]->verts[1] = face->verts[l+1];
				Object.flist[new_faces+l]->verts[2] = Object.nverts-1;
				Object.flist[new_faces+l]->nverts = 3;
			}
			new_faces += face->nverts - 1;
			first_vertex = face->verts[0];
			last_vertex = face->verts[face->nverts -1];
			free(Object.flist[i]->verts);
			Object.flist[i]->verts = (int *)malloc(sizeof(int)*3);
			Object.flist[i]->verts[0] = last_vertex;
			Object.flist[i]->verts[1] = first_vertex;
			Object.flist[i]->verts[2] = Object.nverts-1;
			Object.flist[i]->nverts = 3;
		}
	}
	Object.nfaces = new_faces;

	return 0;
}

void CPolygon3D::get_vertex_normal(void)
{
	int i, j;

	for (i=0; i<Object.nverts; i++) 
		Object.vlist[i]->normal.set(0.0);
	for (i=0; i<Object.nfaces; i++) {
		for (j=0; j<3; j++)
			Object.vlist[Object.flist[i]->verts[j]]->normal += Object.flist[i]->normal * Object.flist[i]->area;
	}
	for (i=0; i<Object.nverts; i++)
		normalize(Object.vlist[i]->normal);

}


#include "time.h"

void jitter_vertex_coords()
{
	int i, j;
	Vertex *v, *adj;
	double sum_x, sum_y;
	double rand_w;

	srand((unsigned)time(NULL));

	/*save the original*/
	for(i=0; i<Object.nverts; i++)
	{
		v = Object.vlist[i];
		v->bx = v->x;
		v->by = v->y;
	}

	for(i=0; i<Object.nverts; i++)
	{
		v = Object.vlist[i];
		sum_x = 0;
		sum_y = 0;

		if(v->Num_edge < 6)
			continue;

		for(j=0; j<v->Num_edge; j++)
		{
			adj = Object.vlist[v->edges[j]->verts[0]];
			if(v == adj)
				adj = Object.vlist[v->edges[j]->verts[1]];
			sum_x += adj->bx;
			sum_y += adj->by;
		}

		//rand_w = (double)rand()/6*RAND_MAX;

		v->x = sum_x/v->Num_edge;
		
		v->y = sum_y/v->Num_edge;
	}
}


void CPolygon3D::PreprocessVertex(void)
{


	for(int i = 0; i < Object.nverts; i++)
	{
		Object.vlist[i]->x = Object.vlist[i]->x *( 0.5/Object.radius) + 0.5;
		Object.vlist[i]->y = Object.vlist[i]->y *( 0.5/Object.radius) + 0.5;
		Object.vlist[i]->z = Object.vlist[i]->z *( 0.5/Object.radius) ;
		//Object.vlist[i]->x = Object.vlist[i]->x *( 0.3/Object.radius) + 0.5;
		//Object.vlist[i]->y = Object.vlist[i]->y *( 0.3/Object.radius) + 0.5;
		//Object.vlist[i]->z = Object.vlist[i]->z *( 0.3/Object.radius) ;
	}
}

void CPolygon3D::SetFileDirectory(char* dir, char* name)
{
	strcpy(temp_dir, dir);
	strcpy(object_name, name);
}

/*********************************************************
Initialize the local frame for each triangle for tracing
*********************************************************/
void CPolygon3D::InitLocalFrame(void)
{
	//icVector3 T, B, N, TM; //define vectors for local frame calculation

	Face *face;
	icVector2 PP;
	int i;


	////Build a local frame for each triangle
	for( i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		face->index = i;  ////get the index of each triangle

		/////using v1v0 as X axis
		face->LX.entry[0] = Object.vlist[face->verts[1]]->x - Object.vlist[face->verts[0]]->x;
		face->LX.entry[1] = Object.vlist[face->verts[1]]->y - Object.vlist[face->verts[0]]->y;

		normalize(face->LX);

		////This is a right hand coordinate system
		//face->LY = cross(face->normal,face->LX);

		////or using following way to find LY
		face->LY.entry[0] = -face->LX.entry[1];
		face->LY.entry[1] =  face->LX.entry[0];

		//normalize(face->LY);

		////At the same time, we need to calculate the coordinates of its three vertices under local frame
		face->xy[0][0] = face->xy[0][1] = 0.;    //the first vertex is the origin of the local frame

		////calculate the local coordinates for verts[1]
		PP.entry[0] = Object.vlist[face->verts[1]]->x - Object.vlist[face->verts[0]]->x;
		PP.entry[1] = Object.vlist[face->verts[1]]->y - Object.vlist[face->verts[0]]->y;

		face->xy[1][0] = dot(PP, face->LX);
		face->xy[1][1] = 0;

		////Now we need to calculate the coordinates for verts[2]
		PP.entry[0] = Object.vlist[face->verts[2]]->x - Object.vlist[face->verts[0]]->x;
		PP.entry[1] = Object.vlist[face->verts[2]]->y - Object.vlist[face->verts[0]]->y;

		face->xy[2][0] = dot(PP, face->LX);
		face->xy[2][1] = dot(PP, face->LY);
	}
}


/************************************************************************
Get the length of the edge 2/15/06
************************************************************************/
double CPolygon3D::GetEdgeLength(int v1, int v2)
{
	icVector3 len;
	len.entry[0] = Object.vlist[v2]->x - Object.vlist[v1]->x;
	len.entry[1] = Object.vlist[v2]->y - Object.vlist[v1]->x;
	return length(len);
}


/************************************************************************
Get the information of all edges in the model after file reading 1/6/05
************************************************************************/

void CPolygon3D::GetEdge()
{
	///First create and initialize the head knot of the edge link

	Object.elist = new Edge;
	Object.elist->index = -1;
	Object.elist->next = NULL;
	Cur_elink = Object.elist;

	global_edge_id = 0;
	Object.nedges = 0;

	////////////////Define variables for vertices and faces operation
	Vertex **vert = Object.vlist;
	///////////////////////////////
	int i, j, m, n;
	int Cur_vert, Next_vert;

	for( i = 0; i < Object.nfaces; i++)
	{
		for( j = 0; j < Object.flist[i]->nverts; j++)
		{
			//We need to check the neighbor vertex on that surface
			if( j == Object.flist[i]->nverts - 1){
				Cur_vert = Object.flist[i]->verts[j];
				Next_vert = Object.flist[i]->verts[0];
			}
			else{
				Cur_vert = Object.flist[i]->verts[j];          //extract the ID of the i'th face and j'th vertex
				Next_vert = Object.flist[i]->verts[j+1];       //extract the ID of the i'th face and j+1'th vertex
			}

			//check if there is any edge between them or not
			if(vert[Cur_vert]->Num_edge == 0 || vert[Next_vert]->Num_edge == 0)
			//there must be no edge between them at this moment
			{
				///Create new notes for the edge link
				Edge *new_edge = new Edge;
				new_edge->index = global_edge_id;              //first edge will be marked 0, and so on...
				global_edge_id ++;

				/////Initialize the id of the adjacent faces that share the edge
				new_edge->tris[0] = -1;
				new_edge->tris[1] = -1;

				new_edge->tris[0] = i;                        //this is the first surface sharing the edge
				new_edge->verts[0] = Cur_vert;                //Save the ids of current vertices as the terminals of the edge
				new_edge->verts[1] = Next_vert;               //Using the current orientation!!!! 1/11
				new_edge->visited = 0;                        //for my subdivision
				new_edge->next = NULL;
				Cur_elink->next = new_edge;                   //Add to the current edge link
				Cur_elink = new_edge;

				/*compute the edge length 07/21/07*/
				new_edge->length = GetEdgeLength(Cur_vert, Next_vert);


				vert[Cur_vert]->edges = 
					Extend_Elist(vert[Cur_vert]->edges, vert[Cur_vert]->Num_edge);
				vert[Next_vert]->edges = 
					Extend_Elist(vert[Next_vert]->edges, vert[Next_vert]->Num_edge);

				vert[Cur_vert]->edges[vert[Cur_vert]->Num_edge] = new_edge;
				vert[Next_vert]->edges[vert[Next_vert]->Num_edge] = new_edge;

				vert[Cur_vert]->Num_edge++;
				vert[Next_vert]->Num_edge++;
						
				/////Add to the face
				Object.flist[i]->edges[j] = new_edge;         //Link the new edge to the associated face

                ////The total number of edges add one
				Object.nedges++;
		    }
			else{
				for( m = 0; m < vert[Cur_vert]->Num_edge; m++)
					for( n = 0; n < vert[Next_vert]->Num_edge; n++)
					{
						if( vert[Cur_vert]->edges[m]->index
							== vert[Next_vert]->edges[n]->index) 
						//There already has an edge between these two vertices
						{

							vert[Cur_vert]->edges[m]->tris[1] = i;

							Object.flist[i]->edges[j] = vert[Cur_vert]->edges[m];

							goto LL;                          //if same edge ID has been found, jump out of the loop

						}
					}
				
LL:				if( m > vert[Cur_vert]->Num_edge - 1 )
				//Did not find an existing edge between these two vertices
				{
					///Create new notes for the edge link
					Edge *new_edge = new Edge;
					new_edge->index = global_edge_id;         //first edge will be marked 0, and so on...
					global_edge_id ++;
					/////Initialize the id of the adjacent faces that share the edge
					new_edge->tris[0] = -1;
					new_edge->tris[1] = -1;

					new_edge->tris[0] = i;                    //this is the first surface sharing the edge
					new_edge->verts[0] = Cur_vert;            //Save the ids of current vertices as the terminals of the edge
					new_edge->verts[1] = Next_vert;
					new_edge->visited = 0;                    //for my subdivision
					new_edge->next = NULL;
					Cur_elink->next = new_edge;               //Add to the current edge link
					Cur_elink = new_edge;
					
					/*compute the edge length 07/21/07*/
					new_edge->length = GetEdgeLength(Cur_vert, Next_vert);

					///Add the ID of the edge into corresponding vertices and faces

					vert[Cur_vert]->edges = 
						Extend_Elist(vert[Cur_vert]->edges, vert[Cur_vert]->Num_edge);
					vert[Next_vert]->edges = 
						Extend_Elist(vert[Next_vert]->edges, vert[Next_vert]->Num_edge);

					vert[Cur_vert]->edges[vert[Cur_vert]->Num_edge] = new_edge;
					vert[Next_vert]->edges[vert[Next_vert]->Num_edge] = new_edge;

					vert[Cur_vert]->Num_edge++;
					vert[Next_vert]->Num_edge++;

					/////Add to the face

				    Object.flist[i]->edges[j] = new_edge;     //Link the new edge to the associated face

					///The total number of edges add one
					Object.nedges++;
				}
			}

		}
	}

	//Object.nedges = global_edge_id;
	TotalEdgesNum = Object.nedges;

	/////////////
	Edge *temp_e = Object.elist->next;
	for(i = 0; i < Object.nedges; i++)
	{
		temp_e->OppVerts[0] = GetOppositeVertices(Object.flist[temp_e->tris[0]], temp_e->verts);
		if(temp_e->tris[1] > 0)
			temp_e->OppVerts[1] = GetOppositeVertices(Object.flist[temp_e->tris[1]], temp_e->verts);

		else
			temp_e->OppVerts[1] = -1;

		temp_e = temp_e->next;
	}

	/*construct the array style edge list for later convenience*/
	Object.edgelist=(Edge**)malloc(sizeof(Edge*)*Object.nedges);
	Edge *e = Object.elist->next;
	for(i=0; i<Object.nedges; i++)
	{
		Object.edgelist[i]= e;
		e=e->next;
	}

	//////Get the average length of the edges, we use the average length of a triangle to approximate it
	////since we employ a regular mesh in 2D plane
	double length_sum = 0;
	Face *face = Object.flist[0];
	Vertex *v1, *v2, *v3;
	icVector2 t_length;
	v1 = Object.vlist[face->verts[0]];
	v2 = Object.vlist[face->verts[1]];
	v3 = Object.vlist[face->verts[2]];

	t_length.entry[0] = v1->x - v2->x;
	t_length.entry[1] = v1->y - v2->y;
	length_sum += length(t_length);
	
	t_length.entry[0] = v3->x - v2->x;
	t_length.entry[1] = v3->y - v2->y;
	length_sum += length(t_length);
	
	t_length.entry[0] = v1->x - v3->x;
	t_length.entry[1] = v1->y - v3->y;
	length_sum += length(t_length);

	avg_edge_length = length_sum /3.;

	htry = avg_edge_length * 50;
}


/**************************************
**************************************/
	
int CPolygon3D::GetOppositeVertices(Face *face, int verts[2])
{
	for(int i = 0; i < face->nverts; i++)
	{
		if(face->verts[i] != verts[0] && face->verts[i] != verts[1])
			return (face->verts[i]);
	}
}

/***************************************************************
Extend the number of the edge incident to a specific vertex
The function extend one space each time
***************************************************************/

Edge  **CPolygon3D::Extend_Elist(Edge **edge_link, int Num_edges)
{
    Edge **temp = edge_link;
	Edge **new_edge_link = (Edge **) malloc (sizeof (Edge*)*(Num_edges+1)); //Extend the link
	if( Num_edges > 0)
	{
		for(int i = 0; i < Num_edges; i++)
			new_edge_link[i] = temp[i];
		free (temp);
	}
   
	return new_edge_link;
}

/***************************************************************
Extend the number of the edge incident to a specific vertex
The function extend one space each time
***************************************************************/

int *CPolygon3D::Extend_link(int *edge_link, int Num_edges)
{
    int *temp = edge_link;
	int *temp_edge_link = new int[Num_edges + 1];
	if( Num_edges > 0)
	{
		for(int i = 0; i < Num_edges; i++)
			temp_edge_link[i] = temp[i];
		delete temp;
	}

	return temp_edge_link;
}

/***************************************************************
Get the vectors on vertices in local frame
This routine must be called after building the local frame
and after get the global vector field on vertices!!!!
***************************************************************/

void CPolygon3D::GetLocalVector()
{
	Vertex *vert;
	Face *face;

    int i, j;

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];
			GlobalToLocal(vert->vec.entry, face->direct_vec[j], face);
		}
	}
}


/*******************************************************************
This routine is used to transfer the vector on global frame to
specific local frame
*******************************************************************/
	
void CPolygon3D::GlobalToLocal(double v[2], icVector2 &iv, Face *face)
{

	icVector2 lx, ly;

	lx.entry[0] = face->LX.entry[0];
	lx.entry[1] = face->LX.entry[1];
	ly.entry[0] = face->LY.entry[0];
	ly.entry[1] = face->LY.entry[1];

	normalize(lx);
	normalize(ly);


	iv.entry[0] = lx.entry[0]*v[0] + lx.entry[1]*v[1];
	iv.entry[1] = ly.entry[0]*v[0] + ly.entry[1]*v[1];

}

/************************************************
Clear the vector field
************************************************/

void CPolygon3D::ClearVectors()
{
	int   i, j; 
	Face *face;
	int *verts;
	   
    //Using triangle mesh to generate the vector field, here is still a 2D vector field
	for (i=0; i<Object.nfaces; i++) {
		face = Object.flist[i];
		verts = face->verts;
		for (j=0; j<face->nverts; j++) {
			Object.vlist[verts[j]]->vec.entry[0] = 0;
			Object.vlist[verts[j]]->vec.entry[1] = 0;
		}
	}	
	
	////transfer the vector into local frame
	GetLocalVector();
}

/*********************************************************************************************
This method is designed to find the opposite corner
It must be called after all the corners have been traversed and basic information (max, min) have been set up
There may be some corner withour opposite corner, because their opposite edge is the boundary of the object
But for a closed 3D object, all the corners should have an opposite corner
***************************************************************************/
void CPolygon3D::FindOpposite()
{
	int i, j, k;

	////Get the opposite corner without sorting
	Vertex **vert = Object.vlist;
	int face_index;

	for( i = 0; i < Object.nfaces; i++)
	{
		for( j = 0; j < Object.flist[i]->nverts; j++)
		{
			if(Object.clist[i*3 + j]->o < 0)
			{  
				/////Getting the opposite face ID through 'c.e' information
				if(Object.clist[i*3 + j]->e->tris[0] != i)
					face_index = Object.clist[i*3 + j]->e->tris[0];
				else
					face_index = Object.clist[i*3 + j]->e->tris[1];

				if(face_index < 0) ////the oposite triangle does not exist!
					continue;

				for( k = 0; k < Object.flist[face_index]->nverts; k++)
				{
					///Search all the corners belong to the opposite face
					///Finding the corner that c'.e has 
					if( Object.clist[face_index*3 + k]->e == Object.clist[i*3 + j]->e)
					{ // This corner is the opposite corner we want to find
						//Because the opposite corners are pairs
						//So we found one pair of opposite corner here
						Object.clist[i*3 + j]->o = face_index*3 + k;
						Object.clist[i*3 + j]->ot = face_index; //get the opposite triangle

						Object.clist[face_index*3 + k]->o = i*3 + j;
						Object.clist[face_index*3 + k]->ot = i;
					}
				}
			}
		}
	}
}

/***************************************************
Building corner table for the mesh
***************************************************/
void CPolygon3D::BuildCorner()
{
    int i, j, m;

	////////////////Define variables for corner information building
	Vertex **vert = Object.vlist;


	for( i = 0; i < Object.nfaces * 3; i++)
	{
		Object.clist[i] = (Corner *) malloc( sizeof(Corner));
		Object.clist[i]->index = 0;
		Object.clist[i]->Edge_count = 0;
		Object.clist[i]->edge[0] = NULL;
		Object.clist[i]->edge[1] = NULL;
	}

	for( i = 0; i < Object.nverts; i++)
		Object.vlist[i]->Num_corners = 0;

	if(Object.elist == NULL)
		GetEdge();

	for( i = 0; i < Object.nfaces; i++)
	{
		for( j = 0; j < Object.flist[i]->nverts; j++)
		{
			//////////////////////////////////////////////////////////////////////////////
			////////Add some information to corner structure 1/21/05
			///Get the current id of triangle to the corner
			Object.clist[i*3 + j]->index = i*3 + j;
			Object.clist[i*3 + j]->t = i;
			Object.clist[i*3 + j]->v = Object.flist[i]->verts[j]; //Get the ID of current vertex

			/////At the same time, we need to add this corner to current vertex
			AddCornertoVertex((i*3 + j), Object.vlist[Object.flist[i]->verts[j]]);

			//Get its previous corner ID
			if( (j - 1) < 0)
				Object.clist[i*3 + j]->p = i*3 + 2;
			else
				Object.clist[i*3 + j]->p = i*3 + (j-1);

			//Get its next corner ID
			Object.clist[i*3 + j]->n = i*3 + (j+1)%3;

			///Get the angle of the corner
			Object.clist[i*3 + j]->angle = GetAngle(Object.flist[i]->verts[j], 
				                                    Object.flist[i]->verts[(3+(j-1))%3], 
						                            Object.flist[i]->verts[(j+1)%3]);



			///Get the edges constitute the corner
			for( m = 0; m < 3; m++)
			{

				if( Object.flist[i]->edges[m]->verts[0] == Object.flist[i]->verts[j]
					||Object.flist[i]->edges[m]->verts[1] == Object.flist[i]->verts[j])
				{
					////If one of the end points of the edge is the current vertex
					////It means that the edge is incident to the vertex
					///Add this edge to the corner associate with the vertex
					if( Object.clist[i*3 + j]->Edge_count == 0)
					{
						Object.clist[i*3 + j]->edge[0] = Object.flist[i]->edges[m];
						Object.clist[i*3 + j]->Edge_count++;
					}
					else
					{
						Object.clist[i*3 + j]->edge[1] = Object.flist[i]->edges[m];
						Object.clist[i*3 + j]->Edge_count++;
					}
				}
				else //it means that the edge is not incident to the vertex,
					//so it must be the opposite edge of that vertex
					Object.clist[i*3 + j]->e = Object.flist[i]->edges[m];
			}

			/////////////////////////////////////////////////////////////////////////
			////Initialize the  opposite corner o to -1
			Object.clist[i*3 + j]->o = -1;
			Object.clist[i*3 + j]->ot = -1;
		}
	}
	
	FindOpposite();     //finding the opposite corner of each corner
	
	
	////Sort the edges for each corner, make edge[0] store the edge that is the opposite edge of corner c.n
	////edge[1]  store the edge that is the opposite edge of corner c.p;

	Corner *c;
	for( i = 0 ; i < Object.ncorners; i++)
	{
		c = Object.clist[Object.clist[i]->n];
		Object.clist[i]->edge[0] = c->e;

		c = Object.clist[Object.clist[i]->p];
		Object.clist[i]->edge[1] = c->e;
	}
}

/************************************************************************
This routine implements add an corner into the vertex
************************************************************************/
void CPolygon3D::AddCornertoVertex(int CornerIndex, Vertex *v)
{
	v->Corners = Extend_link(v->Corners, v->Num_corners);

	v->Corners[v->Num_corners] = CornerIndex;

	v->Num_corners++;

}

/***********************************************************************
Sort the corner list for each vertex in order  12/29/05
***********************************************************************/

void CPolygon3D::SortCornerOnVertex()
{
	int *Newlist;
	int i, j;
	Vertex *vert;
	Corner *c;
	int flag = 0;

	for(i = 0 ; i < Object.nverts; i++)
	{
		vert = Object.vlist[i];
		
		//We need to deal with the corner opposite to the boundary in plane !!!
		//if(vert->Num_corners <= 3)
		//if(vert->Num_corners <= 5) /*to load the simplified mesh!07/17/07*/
		//	continue;

		/*if one of its edge is boundary edge, this vertex is a boundary vertex, continue*/
		for(j=0; j<vert->Num_edge; j++)
		{
			if(vert->edges[j]->tris[0] == -1
				||vert->edges[j]->tris[1] == -1)
			{
				flag = 1;
				break;
			}
		}

		if(flag == 1)
		{
			flag = 0;
			continue;
		}

		Newlist = new int[vert->Num_corners];

		//Newlist[0] = vert->Corners[GetFirstCorner(vert)];
		Newlist[0] = vert->Corners[0];
		for(j = 1; j < vert->Num_corners; j++)
		{

			c = Object.clist[Object.clist[Object.clist[Newlist[j-1]]->p]->o];
			Newlist[j] = c->p;
		}

		for(j = 0; j < vert->Num_corners; j ++)
		{
			vert->Corners[j] = Newlist[j];
		}

		delete[] Newlist;
	}
}

/*******************************************************************
*******************************************************************/
void CPolygon3D::AllocateAng()
{
	int i, j;
	Vertex *vert, *adjv;
	icVector2 edgevec;
	Corner *c;
	double /*sum_ang, */cur_ang;
	double a, b;
	int orient;

	orient = 0;
	/*sum_ang = */cur_ang = 0;
	a = b = 0;



	for(i = 0; i < Object.nverts; i++)
	{
		vert = Object.vlist[i];
		//sum_ang = 0;
		vert->Anglist = new double[vert->Num_corners];

		//////////////////////////////////////////////////////////////////////

		////choose first edge of first corner to do the projection to get the initial angle
		c = Object.clist[vert->Corners[0]];
		if(c->edge[0] != NULL)
		{
			if( c->edge[0]->verts[0] != i)
				adjv = Object.vlist[c->edge[0]->verts[0]];
			else
				adjv = Object.vlist[c->edge[0]->verts[1]];
		}

		else 
			continue;

		edgevec.entry[0] = adjv->x - vert->x;
		edgevec.entry[1] = adjv->y - vert->y;

		normalize(edgevec);

		////adjust the angle to make sure that it is between 0 and 2Pi
		vert->Anglist[0] = atan2(edgevec.entry[1],edgevec.entry[0]);
		if(vert->Anglist[0] < 0)
			vert->Anglist[0] += 2*M_PI;

		cur_ang = vert->Anglist[0];  ////set it as current angle position

		//////Store it as the beginning angle for the first corner in the corners' list
		Object.clist[vert->Corners[0]]->r = 1;
		Object.clist[vert->Corners[0]]->BeginAng = vert->Anglist[0];
		Object.clist[vert->Corners[vert->Num_corners-1]]->EndAng = cur_ang;


		////we may need to first decide the orientation of the angle
		orient = GetOrientation(vert, cur_ang, c->edge[1]);
		c->orient = orient;

		////calculate the angles for other corners/triangles
		for(j = 1 ; j < vert->Num_corners; j++)
		{
			c = Object.clist[vert->Corners[j-1]];

			if(orient > 0)
			{
				vert->Anglist[j] = c->angle + cur_ang;

				if(vert->Anglist[j] >= 2*M_PI)
					vert->Anglist[j] -= 2*M_PI;
			}
			else{
				vert->Anglist[j] = cur_ang - c->angle;

				if(vert->Anglist[j] < 0)  ////modified on 04/23/05
					vert->Anglist[j] += 2*M_PI;
			}

			cur_ang = vert->Anglist[j];
		     
			//////Store it as the beginning angle for the first corner in the corners' list
			Object.clist[vert->Corners[j-1]]->EndAng = cur_ang;
			Object.clist[vert->Corners[j]]->BeginAng = cur_ang;
			Object.clist[vert->Corners[j]]->r = 1;
			Object.clist[vert->Corners[j]]->orient = orient;
		}
	}
}


/**********************************************************************************
This routine is designed to calculate the angle of the corner
**********************************************************************************/

double CPolygon3D::GetAngle(int v, int vp, int vn)
{
    icVector2 v1, v2;
	Vertex **verts = Object.vlist;

	v1.entry[0] = verts[vp]->x - verts[v]->x; 
	v1.entry[1] = verts[vp]->y - verts[v]->y; 

	v2.entry[0] = verts[vn]->x - verts[v]->x; 
	v2.entry[1] = verts[vn]->y - verts[v]->y; 
	
	///Getting the angle of this two vectors
	normalize(v1);
	normalize(v2);

	double ang = acos(dot(v1, v2));

	return ang;
}


/*********************************************************
Get the orientation of edge
*********************************************************/
int CPolygon3D::GetOrientation(Vertex *p, double cur_ang, Edge *e1)
{
	double a, b;
	icVector2 edge1, normal;
	Vertex *adjv;
	double ang1, ang_diff;
	

	if( Object.vlist[e1->verts[0]] != p)
		adjv = Object.vlist[e1->verts[0]];
	else
		adjv = Object.vlist[e1->verts[1]];

	edge1.entry[0] = adjv->x - p->x;
	edge1.entry[1] = adjv->y - p->y;


	ang1 = atan2(edge1.entry[1], edge1.entry[0]);
	if(ang1 < 0)
		ang1 += 2*M_PI;

    ang_diff = ang1 - cur_ang;

	if(ang_diff <= -M_PI)
		ang_diff += 2 * M_PI;
	if(ang_diff > M_PI)
		ang_diff -= 2 * M_PI;

	if(ang_diff >= 0)
		return 1;
	else
		return -1;
}

/**********************************************************
Select the first corner that most aligh with tangent plane
It seems that we do not need this in 2D plane
**********************************************************/

//int CPolygon3D::GetFirstCorner(Vertex *p)
//{
//	int i, firstID = 0;
//	double dot1, maxdot = -1;
//    Face *face;
//	for(i = 0; i < p->Num_corners; i++)
//	{
//		face = Object.flist[Object.clist[p->Corners[i]]->t];
//		dot1 = dot(p->normal, face->normal);
//
//		if(dot1 > maxdot)
//		{
//			firstID = i;
//			maxdot = dot1;
//		}
//	}
//
//	return firstID;
//
//}
