#pragma once
#include "lib/icVector.h"
#include "plyloader.h"

#include "VFDataStructure.h"


/* write an inline function for memory reallocation */

template <class T>
T *Extend_space(T *origin, int size, int &flag)
{
	T *temp = origin;
	flag = 0;

	if((origin = (T*)realloc((T*) origin, sizeof(T) * size)) == NULL)
	{
		MessageBox(NULL, "fail to reallocate!", "Error", MB_OK);
		origin = temp;
		flag = 1;
	}

	return origin;
}


class CPolygon3D
{
public:
	CPolygon3D(void);
	~CPolygon3D(void);


    //Polygon3D Object;
    double vd_n, vd_f;
//	char object_name[128];  //using global varialbes at this moment
//	char temp_dir[128]; 
	PlyFile *in_ply;
	icVector3 rot_center;
	int orientation;

    char object_name[128];
    char temp_dir[128]; 

	CPlyLoader MyPlyLoader;

	//Polygon3D Object;      //06/24/05

	//void get_object();
	Polygon3D get_object();
	void read_object_file(const char * filename);
	void calc_bounding_sphere(void);
	void calc_face_normals_and_areas(void);
	int triangulate(void);
	void get_vertex_normal(void);

	void PreprocessVertex(void);
	void SetFileDirectory(char* dir, char* name);

	void InitLocalFrame(void);       ////Building the local frame for each triangle

	void ClearVectors();


	//////For edge building
	int GetOppositeVertices(Face *face, int verts[2]);
	double GetEdgeLength(int v1, int v2);

	void GetEdge(); //Used to build the edge table
    Edge  **Extend_Elist(Edge **edge_link, int Num_edges);
	int global_edge_id ;
    Edge *Cur_elink;

	/////build the local vectors 
	void GetLocalVector();
	void GlobalToLocal(double v[2], icVector2 &iv, Face *face);

	////For corner building
    static int *Extend_link(int *edge_link, int Num_edges);
    void BuildCorner();
    void AddCornertoVertex(int CornerIndex, Vertex *v);
    void FindOpposite();

	void SortCornerOnVertex();
	void AllocateAng();
    double GetAngle(int v, int vp, int vn);
    int GetOrientation(Vertex *p, double cur_ang, Edge *e1);


};
