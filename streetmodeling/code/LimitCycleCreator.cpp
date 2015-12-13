/////
////LimitCycleCreator.cpp

#include "stdafx.h"
#include "LimitCycleCreator.h"

#include "VFSynthesis.h"

#include "VFDataStructure.h"
#include "LocalTracing.h"

#include "RegionSmoothing.h"
#include "HermiteCurve.h"

#include "topologyedit.h"

/////////////////////////////////////////////////
extern Polygon3D Object;

extern RegularElement *regularelem;            //regular elememts' list
extern int cur_regelem_index;
extern int MaxNumRegularElems;           //Maximum number of regular elements

extern double  dmax ;


/////variables for shape design
//extern ctr_point *pts;          // allocate our control points array
extern int resolution;    // how many points in our output array
extern ctr_point *out_pts;

extern ctr_point *control_pts;        // allocate our control point array
extern int num_shapecontrol_pts;
extern int num_curvepts_output;
extern int MaxNumShapeControlPts;
extern int HermiteStep;
extern icVector2 *CtrPts_tangent;


////Global variables for new limit cycle creation by shape design
int *DesignCurveCellCycle;
int num_triangles_designcurve;
int MaxNumTrianglesDesignCurve;


////New data structure for new limit cycle creation 5/20/06
DesignTriangleCycle myCycle;

int *Boundaryverts;
int num_boundverts = 0;
int MaxNumBoundVerts;


RegionBoundary InnerBoundary;
RegionBoundary OuterBoundary;
int MaxNumEdgesOnBoundary;


////Variable to store the points and the corresponding triangles that contain the points
LineSeg *designcurve;
int num_lineseg_designcurve;

Boundaryverts_List boundaryvertlist;  //store the boundary list
int MaxNumBoundaryVerts;

extern TriangularRegion intersectRegion;     ////The intersect region

extern Edge **Cycle_edge;
extern int num_cycleedges;
extern int MaxEdgeInCycle;

extern Vertex **regionverts;                ////mesh vertices inside user selected region
extern Edge **regionedge;                   ////mesh edges of user selected region
extern int Num_verts;                       ////number of inner vertices
extern int Num_edges;

extern InnerVertices repellerInnerverts;
extern bool TriangleSearch(int *acycle, int num_cycletriangles, int oneTriangle, int &position);
extern void AddToBoundaryEdgeList(Edge*);
extern void AddToInnerVerts(Vertex *vert, int type);

extern void CalNormalAtEdge(Edge *cur_edge, Face *face, int type);


extern void Hermitecurve(int n, ctr_point *interpts, icVector2 *T, ctr_point *output, int step, int &num_output);
extern void Hermitecurve_open(int n, ctr_point *interpts, icVector2 *T, ctr_point *output, int step, int &num_output);

extern void  InitLimitCycleDetection();
extern void InitLimitCycleStructure();

extern void InitRegionSmooth();

extern int TriangleDetect(double,double);


/////For paper images
int *InnerTriangles = NULL;
int num_innertriangles;
int limitcycletype;

////Testing
int test_num_cycleedges = 0;


//////////////////////////////////////////////////////////////////////////////
////New method to generate limit cycles through shape design

int GetOtherTriangle(Edge *cur_e)
{
	if(Object.flist[cur_e->tris[0]]->inDesignCellCycle == 1)
		return cur_e->tris[1];
	else
		return cur_e->tris[0];
}


////Judge the connecting relationship between two triangles
bool IsNeighborTriangles(int triangle1, int triangle2, int &SingleSharingVert)
{
	int i, j;
	Face *face1, *face2;
	face1 = Object.flist[triangle1];
	face2 = Object.flist[triangle2];

	int vert1, vert2;
	int num_sharingverts = 0;

	for(i = 0; i < face1->nverts; i++)
	{
		vert1 = face1->verts[i];
		for(j = 0; j < face2->nverts; j++)
		{
			vert2 = face2->verts[j];
			if(vert1 == vert2){
				SingleSharingVert = vert2;
				num_sharingverts ++;
			}
		}
	}

	if(num_sharingverts > 1 || num_sharingverts < 1) return true;
	else	return false;  ////not a neighbor triangle, but they need to share at least one vertex!!!!
}


////Connecting the two triangle with other intermedia triangles
////Note that here triangle1 is the previous triangle
////This routine may cause breaking down of the program, you make need to capture the exception here

void ConnectTwoTriangles(int triangle1, int triangle2, int SingleSharingVert)
{
	int i;
	Corner *c1, *c2, *temp_c/*, *temp_o*/;

	int *path1, *path2;
	int numtriangles_path1, numtriangles_path2;

	path1 = (int *) malloc(sizeof(int) * Object.vlist[SingleSharingVert]->Num_corners);
	path2 = (int *) malloc(sizeof(int) * Object.vlist[SingleSharingVert]->Num_corners);
	numtriangles_path1 = numtriangles_path2 = 0;

	////Get the two corners that associate with the two triangles respectively
	for(i = 0; i < 3; i++)
	{
		if(Object.clist[3*triangle1+i]->v == SingleSharingVert)
			c1 = Object.clist[3*triangle1+i];

		if(Object.clist[3*triangle2+i]->v == SingleSharingVert)
			c2 = Object.clist[3*triangle2+i];
	}

	////if these two triangle can be connected by a single intermedia triangle
	if(Object.clist[c1->p]->ot == Object.clist[c2->n]->ot)
	{
		DesignCurveCellCycle[num_triangles_designcurve] = Object.clist[c1->p]->ot;
		Object.flist[Object.clist[c1->p]->ot]->inDesignCellCycle = 1;
		num_triangles_designcurve++;
		return;
	}

	if(Object.clist[c1->n]->ot == Object.clist[c2->p]->ot)
	{
		DesignCurveCellCycle[num_triangles_designcurve] = Object.clist[c1->n]->ot;
		Object.flist[Object.clist[c1->n]->ot]->inDesignCellCycle = 1;
		num_triangles_designcurve++;
		return;
	}

	////Find the intermedia triangles to connect these two triangles
	temp_c = c1;
	while(temp_c != c2)
	{
		path1[numtriangles_path1] = Object.clist[temp_c->p]->ot;
		//Object.flist[Object.clist[temp_c->p]->ot]->inDesignCellCycle = 1;
		numtriangles_path1 ++;
		temp_c = Object.clist[Object.clist[Object.clist[temp_c->p]->o]->p];
	}

	temp_c = c1;
	while(temp_c != c2)
	{
		path2[numtriangles_path2] = Object.clist[temp_c->n]->ot;
		//Object.flist[Object.clist[temp_c->n]->ot]->inDesignCellCycle = 1;
		numtriangles_path2 ++;
		temp_c = Object.clist[Object.clist[Object.clist[temp_c->n]->o]->n];
	}

	if(numtriangles_path1 <= numtriangles_path2)
	{
		for(i = 0; i < numtriangles_path1; i++)
		{
			DesignCurveCellCycle[num_triangles_designcurve] = path1[i];
			Object.flist[path1[i]]->inDesignCellCycle = 1;
			num_triangles_designcurve++;
		}
	}
	else{
		for(i = 0; i < numtriangles_path2; i++)
		{
			DesignCurveCellCycle[num_triangles_designcurve] = path2[i];
			Object.flist[path2[i]]->inDesignCellCycle = 1;
			num_triangles_designcurve++;
		}
	}

	free(path1);
	free(path2);
}


////Used by the new cell cycle detection algorithm 5/22/06
void ConnectTwoTriangles2(int triangle1, int triangle2, int SingleSharingVert, 
						  int *DesignCurveCellCycle, int &num_triangles_designcurve)
{
	int i;
	Corner *c1, *c2, *temp_c/*, *temp_o*/;

	int *path1, *path2;
	int numtriangles_path1, numtriangles_path2;

	path1 = (int *) malloc(sizeof(int) * Object.vlist[SingleSharingVert]->Num_corners);
	path2 = (int *) malloc(sizeof(int) * Object.vlist[SingleSharingVert]->Num_corners);
	numtriangles_path1 = numtriangles_path2 = 0;

	////Get the two corners that associate with the two triangles respectively
	for(i = 0; i < 3; i++)
	{
		if(Object.clist[3*triangle1+i]->v == SingleSharingVert)
			c1 = Object.clist[3*triangle1+i];

		if(Object.clist[3*triangle2+i]->v == SingleSharingVert)
			c2 = Object.clist[3*triangle2+i];
	}

	////if these two triangle can be connected by a single intermedia triangle
	if(Object.clist[c1->p]->ot == Object.clist[c2->n]->ot)
	{
		DesignCurveCellCycle[num_triangles_designcurve] = Object.clist[c1->p]->ot;
		Object.flist[Object.clist[c1->p]->ot]->inDesignCellCycle = 1;
		num_triangles_designcurve++;
		return;
	}

	if(Object.clist[c1->n]->ot == Object.clist[c2->p]->ot)
	{
		DesignCurveCellCycle[num_triangles_designcurve] = Object.clist[c1->n]->ot;
		Object.flist[Object.clist[c1->n]->ot]->inDesignCellCycle = 1;
		num_triangles_designcurve++;
		return;
	}

	////Find the intermedia triangles to connect these two triangles
	temp_c = c1;
	while(temp_c != c2)
	{
		path1[numtriangles_path1] = Object.clist[temp_c->p]->ot;
		//Object.flist[Object.clist[temp_c->p]->ot]->inDesignCellCycle = 1;
		numtriangles_path1 ++;
		temp_c = Object.clist[Object.clist[Object.clist[temp_c->p]->o]->p];
	}

	temp_c = c1;
	while(temp_c != c2)
	{
		path2[numtriangles_path2] = Object.clist[temp_c->n]->ot;
		//Object.flist[Object.clist[temp_c->n]->ot]->inDesignCellCycle = 1;
		numtriangles_path2 ++;
		temp_c = Object.clist[Object.clist[Object.clist[temp_c->n]->o]->n];
	}

	if(numtriangles_path1 <= numtriangles_path2)
	{
		for(i = 0; i < numtriangles_path1; i++)
		{
			DesignCurveCellCycle[num_triangles_designcurve] = path1[i];
			Object.flist[path1[i]]->inDesignCellCycle = 1;
			num_triangles_designcurve++;
		}
	}
	else{
		for(i = 0; i < numtriangles_path2; i++)
		{
			DesignCurveCellCycle[num_triangles_designcurve] = path2[i];
			Object.flist[path2[i]]->inDesignCellCycle = 1;
			num_triangles_designcurve++;
		}
	}

	free(path1);
	free(path2);
}



////Get the triangular strip that contains the design hermite curve for limit cycle generation
////Note that we need to make sure the cell cycle is a manifold region
void GetCellCycleforDesignCurve()
{
	int i;
	int position;
	int cur_triangle = -1;

	int SingleSharingVert;


	////initialize
	for(i = 0; i < Object.nfaces; i++)
		Object.flist[i]->inDesignCellCycle = 0;

	num_triangles_designcurve = 0;

	////After calculating the points on the design curve, we perform triangle detection
	for(i = 0; i < num_curvepts_output; i++)
	{
		cur_triangle = TriangleDetect(out_pts[i].x, out_pts[i].y);
		
		////if the triangle is not inside current list, add it into the list
		if(!TriangleSearch(DesignCurveCellCycle, num_triangles_designcurve, cur_triangle, position)
			&& cur_triangle >= 0 && cur_triangle < Object.nfaces-1)
		{
			////Add more constrains to make sure that current being added triangle has one sharing edge
			////with previous triangle. If this condition is not satisfied, add some intermedia triangles
            
			////if new added triangle has 2 vertices sharing with previous triangle, simply add it 
			////else we need to find the intermedia triangles

			if(num_triangles_designcurve > 0
				&& !IsNeighborTriangles(DesignCurveCellCycle[num_triangles_designcurve-1], cur_triangle, SingleSharingVert))
			{
				ConnectTwoTriangles(DesignCurveCellCycle[num_triangles_designcurve-1], cur_triangle, SingleSharingVert);
			}


			DesignCurveCellCycle[num_triangles_designcurve] = cur_triangle;
			num_triangles_designcurve++;
			Object.flist[cur_triangle]->inDesignCellCycle = 1;
		}
	}

	////we still need to test the first and last triangles in the list
	if(!IsNeighborTriangles(DesignCurveCellCycle[num_triangles_designcurve-1], DesignCurveCellCycle[0], SingleSharingVert))
		ConnectTwoTriangles(DesignCurveCellCycle[num_triangles_designcurve-1], DesignCurveCellCycle[0], SingleSharingVert);
}



////Judge whether the two triangles share the specific vertex
////2/21/06
bool ShareTheSameVertex(int t1, int t2, int &v)
{
	Face *face1 = Object.flist[t1];
	Face *face2 = Object.flist[t2];

	int i, j;

	for(i = 0; i < face1->nverts; i++)
	{
		for(j = 0; j < face2->nverts; j++)
		{
			if(Object.vlist[face1->verts[i]]->VertID == v &&
				Object.vlist[face2->verts[j]]->VertID == v)
				return true;
		}
	}
	return false;
}


////new method to get design cell cycle
void GetDesignCellCycle_new()
{
	int i;
	double cur_p[2], pre_p[2];
	double alpha[3];
	icVector2 VP;
	int cur_triangle = TriangleDetect(out_pts[0].x, out_pts[0].y);
	Face *cur_face = Object.flist[cur_triangle];

	int Passvertornot;
	int position;
	Vertex *SingleSharingVert;

    ////initialize part
	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->inDesignCellCycle = 0;
		Object.flist[i]->discard = 0;
	}

	num_triangles_designcurve = 0;

	num_lineseg_designcurve = 0;


	DesignCurveCellCycle[num_triangles_designcurve] = cur_triangle;
	num_triangles_designcurve++;
	Object.flist[cur_triangle]->inDesignCellCycle = 1;

	VP.entry[0] = out_pts[0].x - Object.vlist[cur_face->verts[0]]->x;
	VP.entry[1] = out_pts[0].y - Object.vlist[cur_face->verts[0]]->y;


	pre_p[0] = dot(VP, cur_face->LX);
	pre_p[1] = dot(VP, cur_face->LY);

	//// Store to the design curve line segment for future usage 12/28/05
	//designcurve[num_lineseg_designcurve].gstart[0] = out_pts[0].x;
	//designcurve[num_lineseg_designcurve].gstart[1] = out_pts[0].y;
	//designcurve[num_lineseg_designcurve].gend[0] = out_pts[1].x;
	//designcurve[num_lineseg_designcurve].gend[1] = out_pts[1].y;
	//designcurve[num_lineseg_designcurve].start[0] = pre_p[0];
	//designcurve[num_lineseg_designcurve].start[1] = pre_p[1];
	//designcurve[num_lineseg_designcurve].Triangle_ID = cur_triangle;
	//num_lineseg_designcurve ++;
	
	
	///// Finding the triangle strip containing the design curve
	for(i = 1; i < resolution; i++)
	{
		VP.entry[0] = out_pts[i].x - Object.vlist[cur_face->verts[0]]->x;
		VP.entry[1] = out_pts[i].y - Object.vlist[cur_face->verts[0]]->y;

		cur_p[0] = dot(VP, cur_face->LX);
		cur_p[1] = dot(VP, cur_face->LY);

		////
		Get2DBarycentricFacters(cur_triangle, cur_p[0], cur_p[1], alpha);
		
		if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1 
			&& (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1 
			&& (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
		{
			pre_p[0] = cur_p[0];
			pre_p[1] = cur_p[1];
			
			/*-------------------------------------------------------------------*/
			//// Store to the design curve line segment for future usage 12/28/05
			//designcurve[num_lineseg_designcurve].gstart[0] = out_pts[i-1].x;
			//designcurve[num_lineseg_designcurve].gstart[1] = out_pts[i-1].y;
			//designcurve[num_lineseg_designcurve].gend[0] = out_pts[i].x;
			//designcurve[num_lineseg_designcurve].gend[1] = out_pts[i].y;
			//designcurve[num_lineseg_designcurve].start[0] = pre_p[0];
			//designcurve[num_lineseg_designcurve].start[1] = pre_p[1];
	  //      designcurve[num_lineseg_designcurve].Triangle_ID = cur_triangle;
			//num_lineseg_designcurve ++;
			/*-------------------------------------------------------------------*/
			continue;  ////still in the same triangle, continue to test next point

		}

		else{

LL:			GetNextTriangle_new(cur_triangle, pre_p, cur_p, Passvertornot, alpha);
			//GetNextTriangle_new2(cur_triangle, pre_p, cur_p, DesignCurveCellCycle, num_triangles_designcurve);

			if(cur_triangle == -1) ////reach boundary of the mesh
				return;

			if(Passvertornot > 0)  ////it passed a vertex
			{
				SingleSharingVert = Object.vlist[cur_face->verts[Passvertornot - 1]];

				////If the two triangles do not share the same vertex, we need to add more segments before
				////we perform following operations 2/5/06
				////because we have cell cycle extension process, so it is not very important issue
				////but if these two triangles do not share the same vertex, we do not want to call the following routine

				if(ShareTheSameVertex(DesignCurveCellCycle[num_triangles_designcurve-1], cur_triangle, SingleSharingVert->VertID))
				{
					ConnectTwoTriangles(DesignCurveCellCycle[num_triangles_designcurve-1], cur_triangle, SingleSharingVert->VertID);
				}
			}
			
			////it seems that we still need to test whether it is the triangle we've accessed before

			////if it is a new triangle in the list, add it
			if(!TriangleSearch(DesignCurveCellCycle, num_triangles_designcurve, cur_triangle, position)
				&& !(InBoundary(cur_triangle)) && Object.flist[cur_triangle]->discard != 1) //modified at 3/29/06
			{
				if(num_triangles_designcurve >= MaxNumTrianglesDesignCurve-1)
				{
					MaxNumTrianglesDesignCurve += 50;
					DesignCurveCellCycle = (int *)realloc(DesignCurveCellCycle, \
						sizeof(int)*MaxNumTrianglesDesignCurve);
				}
				DesignCurveCellCycle[num_triangles_designcurve] = cur_triangle;
				num_triangles_designcurve++;
				Object.flist[cur_triangle]->inDesignCellCycle = 1;

				cur_face = Object.flist[cur_triangle];
				cur_face->discard = 1;  //3/29/06
			}
         
			if(i < resolution -1)
			{
				pre_p[0] = cur_p[0];
				pre_p[1] = cur_p[1];

				cur_face = Object.flist[cur_triangle];
						
				VP.entry[0] = out_pts[i+1].x - Object.vlist[cur_face->verts[0]]->x;
				VP.entry[1] = out_pts[i+1].y - Object.vlist[cur_face->verts[0]]->y;

				cur_p[0] = dot(VP, cur_face->LX);
				cur_p[1] = dot(VP, cur_face->LY);

				Get2DBarycentricFacters(cur_triangle, cur_p[0], cur_p[1], alpha);

				if(!( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1 
					&& (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1 
					&& (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1))
				{
					goto LL;
				}
			}

			////Note that we do not consider the line segment overlaps with an edge
		}


		/*-------------------------------------------------------------------*/
		//// Store to the design curve line segment for future usage 12/28/05
		//designcurve[num_lineseg_designcurve].gstart[0] = out_pts[i-1].x;
		//designcurve[num_lineseg_designcurve].gstart[1] = out_pts[i-1].y;
		//designcurve[num_lineseg_designcurve].gend[0] = out_pts[i].x;
		//designcurve[num_lineseg_designcurve].gend[1] = out_pts[i].y;
		//designcurve[num_lineseg_designcurve].start[0] = pre_p[0];
		//designcurve[num_lineseg_designcurve].start[1] = pre_p[1];
	 //   designcurve[num_lineseg_designcurve].Triangle_ID = cur_triangle;
		//num_lineseg_designcurve ++;
		/*-------------------------------------------------------------------*/
	}

	/////Get the inner design cell cycle for visualization only 3/27/06
    if(InnerTriangles != NULL)
		free(InnerTriangles);
	InnerTriangles = (int*)malloc(sizeof(int)*num_triangles_designcurve);

    for(i = 0; i < num_triangles_designcurve; i++)
	{
		InnerTriangles[i] = DesignCurveCellCycle[i];
	}

	num_innertriangles = num_triangles_designcurve;
}


////new method to get design cell cycle 5/20/06
int GetOneTriangle(int &pos, double globalp[2], int &cur_triangle, int num_curvepts, 
				   int *DesignCurveCellCycle, int &num_triangles_designcurve)
{
	int Passvertornot;
	Vertex *SingleSharingVert;
	double cur_p[2], pre_p[2] = {0.};
	double alpha[3];
	icVector2 VP;
	Face *cur_face = Object.flist[cur_triangle];
	Vertex *v0 = Object.vlist[cur_face->verts[0]];

	VP.entry[0] = globalp[0] - v0->x;
	VP.entry[1] = globalp[1] - v0->y;
	
	pre_p[0] = dot(VP, cur_face->LX);
	pre_p[1] = dot(VP, cur_face->LY);

	while(pos < num_curvepts)
	{
		//Get the out_pts[pos]
		//Judge whether it is inside cur_triangle
		VP.entry[0] = out_pts[pos].x - v0->x;
		VP.entry[1] = out_pts[pos].y - v0->y;

		cur_p[0] = dot(VP, cur_face->LX);
		cur_p[1] = dot(VP, cur_face->LY);
		Get2DBarycentricFacters(cur_triangle, cur_p[0], cur_p[1], alpha);
		
		if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1 
			&& (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1 
			&& (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
		{
			pre_p[0] = cur_p[0];
			pre_p[1] = cur_p[1];

			globalp[0] = out_pts[pos].x;
			globalp[1] = out_pts[pos].y;

			pos++;

			if(pos >= num_curvepts)
				return -1;  //we finish the searching
		}

		else
			break;
	}
		

	GetNextTriangle_new(cur_triangle, pre_p, cur_p, Passvertornot, alpha);

	if(cur_triangle == -1) ////reach boundary of the mesh
		return -1;

	if(Passvertornot > 0)  ////it passed a vertex
	{
		SingleSharingVert = Object.vlist[cur_face->verts[Passvertornot - 1]];

		////If the two triangles do not share the same vertex, we need to add more segments before
		////we perform following operations 2/5/06
		////because we have cell cycle extension process, so it is not very important issue
		////but if these two triangles do not share the same vertex, we do not want to call the following routine

		if(ShareTheSameVertex(DesignCurveCellCycle[num_triangles_designcurve-1], cur_triangle, SingleSharingVert->VertID))
		{
			ConnectTwoTriangles(DesignCurveCellCycle[num_triangles_designcurve-1], cur_triangle, SingleSharingVert->VertID);
		}
	}
	
	//Update the globalp variable
	VP = cur_p[0] * cur_face->LX + cur_p[1] * cur_face->LY;
	globalp[0] = v0->x + VP.entry[0];
	globalp[1] = v0->y + VP.entry[1];

	return cur_triangle;
}


int GetOneTriangle2(int &pos, double globalp[2], int &cur_triangle, int num_curvepts, 
				   int *DesignCurveCellCycle, int &num_triangles_designcurve, int &Passvertornot)
{
	Vertex *SingleSharingVert;
	double cur_p[2], pre_p[2] = {0.};
	double alpha[3];
	icVector2 VP;
	Face *cur_face = Object.flist[cur_triangle];
	Vertex *v0 = Object.vlist[cur_face->verts[0]];

	VP.entry[0] = globalp[0] - v0->x;
	VP.entry[1] = globalp[1] - v0->y;
	
	pre_p[0] = dot(VP, cur_face->LX);
	pre_p[1] = dot(VP, cur_face->LY);

	while(pos < num_curvepts)
	{
		//Get the out_pts[pos]
		//Judge whether it is inside cur_triangle
		VP.entry[0] = out_pts[pos].x - v0->x;
		VP.entry[1] = out_pts[pos].y - v0->y;

		cur_p[0] = dot(VP, cur_face->LX);
		cur_p[1] = dot(VP, cur_face->LY);
		Get2DBarycentricFacters(cur_triangle, cur_p[0], cur_p[1], alpha);
		
		if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1 
			&& (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1 
			&& (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
		{
			pre_p[0] = cur_p[0];
			pre_p[1] = cur_p[1];

			globalp[0] = out_pts[pos].x;
			globalp[1] = out_pts[pos].y;

			pos++;

			if(pos >= num_curvepts)
				return -1;  //we finish the searching
		}

		else
			break;
	}
		

	GetNextTriangle_new(cur_triangle, pre_p, cur_p, Passvertornot, alpha);

	if(cur_triangle == -1) ////reach boundary of the mesh
		return -1;

	if(Passvertornot > 0)  ////it passed a vertex
	{
		//Do not update the global point for this case!!! 5/22/06
		//globalp[0] = Object.vlist[cur_face->verts[Passvertornot-1]]->x;
		//globalp[1] = Object.vlist[cur_face->verts[Passvertornot-1]]->y;
	}
	else
	{
		//Update the globalp variable
		VP = cur_p[0] * cur_face->LX + cur_p[1] * cur_face->LY;
		globalp[0] = v0->x + VP.entry[0];
		globalp[1] = v0->y + VP.entry[1];
	}

	return cur_triangle;
}


void GetDesignCellCycle_new3()
{
	int i;
	int cur_triangle = TriangleDetect(out_pts[0].x, out_pts[0].y);
	if(cur_triangle < 0)
	{
		MessageBox(NULL, "fail to find the first triangle!", "Error", MB_OK);
		return;
	}

	double globalp[2] = {0.};
	int Passvertornot, pre_pos;

    ////initialize part
	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->inDesignCellCycle = 0;
		Object.flist[i]->discard = 0;
	}

	num_triangles_designcurve = 0;
	DesignCurveCellCycle[num_triangles_designcurve] = cur_triangle;
	num_triangles_designcurve++;
	Object.flist[cur_triangle]->inDesignCellCycle = 1;

	globalp[0] = out_pts[0].x;
	globalp[1] = out_pts[1].y;

	////find the cell cycle (closed triangle strip) containing the design curve
	for(i = 1; i < resolution; )
	{
		//cur_triangle = GetOneTriangle(i, globalp, cur_triangle, resolution, 
		//	DesignCurveCellCycle, num_triangles_designcurve);
		Passvertornot=0;

		pre_pos = i;

		cur_triangle = GetOneTriangle2(i, globalp, cur_triangle, resolution, 
			DesignCurveCellCycle, num_triangles_designcurve, Passvertornot);
		if(cur_triangle == -1)
			return;

		if(Passvertornot>0)
		{
			///jitter the position of current point 5/22/06
			out_pts[i].x += 2*pow(-1., rand()%2)*(double)rand()/((double)RAND_MAX + 1)*Object.radius/sqrt((double)Object.nfaces);
			out_pts[i].y += 2*pow(-1., rand()%2)*(double)rand()/((double)RAND_MAX + 1)*Object.radius/sqrt((double)Object.nfaces);

			i = pre_pos;
			continue;
		}
		
		//if it is new triangle, add to the triangle list
		if(!IsRepeated(DesignCurveCellCycle, cur_triangle, num_triangles_designcurve)
			&& !(InBoundary(cur_triangle)) && Object.flist[cur_triangle]->discard != 1) 
		{
			if(num_triangles_designcurve >= MaxNumTrianglesDesignCurve-1)
			{
				MaxNumTrianglesDesignCurve += 50;
				DesignCurveCellCycle = (int *)realloc(DesignCurveCellCycle, \
					sizeof(int)*MaxNumTrianglesDesignCurve);
			}
			DesignCurveCellCycle[num_triangles_designcurve] = cur_triangle;
			num_triangles_designcurve++;
			Object.flist[cur_triangle]->inDesignCellCycle = 1;

			Object.flist[cur_triangle]->discard = 1; //avoid to include other cell belonging to other cycle
		}
	}
}



void GetDesignConnectedEdgelist()
{
	int i;
	int cur_sel_tri;
	int vert1, vert2;

	for(i = 0 ; i < num_shapecontrol_pts; i++)
	{
		cur_sel_tri = TriangleDetect(control_pts[i].x, control_pts[i].y);
		vert1 = Object.flist[cur_sel_tri]->verts[0];

		cur_sel_tri = TriangleDetect(control_pts[(i+1)%num_shapecontrol_pts].x, 
			control_pts[(i+1)%num_shapecontrol_pts].y);
		vert2 = Object.flist[cur_sel_tri]->verts[0];

		RecursiveDijistra(vert1, vert2, vert1, 0);
	}
}



////use edges to simulate the curve 2/5/06

////Get the design cell cycle according to the selected edges
////The input is the Cycle_edge list

void GetDesignCellCycle_new2()
{
	int i, j;
	Edge *cur_edge;
	int cur_triangle;
	int position = -1;

    ////initialize part
	for(i = 0; i < Object.nfaces; i++)
		Object.flist[i]->inDesignCellCycle = 0;

	////Get the connected edge list according to the recorded controlling point
	GetDesignConnectedEdgelist();

	num_triangles_designcurve = 0;

	for(i = 0; i < Num_edges; i++)
	{
		//cur_edge = Cycle_edge[i];
		cur_edge = regionedge[i];

		for(j = 0; j < 2; j++)
		{
			cur_triangle = cur_edge->tris[j];
			////if it is a new triangle in the list, add it
			if(!TriangleSearch(DesignCurveCellCycle, num_triangles_designcurve, cur_triangle, position))
			{
				if(num_triangles_designcurve >= MaxNumTrianglesDesignCurve - 1)
				{
					MaxNumTrianglesDesignCurve += 50;
					DesignCurveCellCycle = (int*)realloc(DesignCurveCellCycle, sizeof(MaxNumTrianglesDesignCurve));
				}

				DesignCurveCellCycle[num_triangles_designcurve] = cur_triangle;
				num_triangles_designcurve++;
				Object.flist[cur_triangle]->inDesignCellCycle = 1;

				//cur_face = Object.flist[cur_triangle];
			}
		}
	}
}



////Begin from the vertices on current boundaries to extend the cell cycle containing the design curve
void ExtendDesignCellCycle()
{
	int i, j;
	Edge *cur_e;
	Vertex *cur_v;
	Corner *cur_c;
	int position;
	int cur_triangle;

	for(i = 0; i < Object.nfaces; i++)
		Object.flist[i]->discard = 0;

	BuildBoundaryEdgeList(DesignCurveCellCycle, num_triangles_designcurve); ////get current boundary edges list

	////Initialize. Added at 3/9/06
	for(i = 0; i < num_cycleedges; i++)
	{
		cur_e = Cycle_edge[i];

		cur_v = Object.vlist[cur_e->verts[0]];
		cur_v->repell_flag = 0;

		cur_v = Object.vlist[cur_e->verts[1]];
		cur_v->repell_flag = 0;
	}

	////Get the vertices on the current boundaries
	repellerInnerverts.num = 0;
	for(i = 0; i < num_cycleedges; i++)
	{
		cur_e = Cycle_edge[i];

		cur_v = Object.vlist[cur_e->verts[0]];
		AddToInnerVerts(cur_v, 0);
		
		cur_v = Object.vlist[cur_e->verts[1]];
		AddToInnerVerts(cur_v, 0);
	}

	////Extend the cell cycle by adding all the adjacent triangles of the vertices into the cycle list
	////Suppose current vertices on the boundaries are stored in the variable "repellerInnerverts"
	for(i = 0; i < repellerInnerverts.num; i++)
	{
		cur_v = repellerInnerverts.vertslist[i];

		for(j = 0; j < cur_v->Num_corners; j++)
		{
			cur_c = Object.clist[cur_v->Corners[j]];
			cur_triangle = cur_c->t;
			if(!TriangleSearch(DesignCurveCellCycle, num_triangles_designcurve, cur_triangle, position)
				&& 	!(InBoundary(cur_triangle)) && Object.flist[cur_triangle]->discard != 1) //3/29/06
			{
				if(num_triangles_designcurve >= MaxNumTrianglesDesignCurve - 1)
				{
					MaxNumTrianglesDesignCurve += 50;
					DesignCurveCellCycle = (int*)realloc(DesignCurveCellCycle, sizeof(MaxNumTrianglesDesignCurve));
				}

				DesignCurveCellCycle[num_triangles_designcurve] = cur_triangle;
				Object.flist[cur_triangle]->inDesignCellCycle = 1;
				num_triangles_designcurve ++;

				Object.flist[cur_triangle]->discard = 1; //3/29/06
			}

		}
	}

	////reset
	for(i = 0; i < num_cycleedges; i++)
	{
		Cycle_edge[i]->OnBoundary = 0;
		Cycle_edge[i]->BoundaryVisited = 0;
	}
	num_cycleedges = 0;

	for(i = 0; i < repellerInnerverts.num; i++)
	{
		////reset the flag
		repellerInnerverts.vertslist[i]->repell_flag = 0;
	}
	repellerInnerverts.num = 0;
}


void GetNextTriangle_new(int &face_id, double pre_p[2], double cur_p[2], int &PassVertornot, double alpha[3])
{
	int which_edge = -1;

	int prev_face_id = face_id;

	Face *prev_face = Object.flist[face_id];

	Vertex *vert = NULL;

	PassVertornot = 0;
	
	////We should put pass vertex testing here before testing crossing edge
	int passornot = 0;
	
	CrossVert_new(face_id, cur_p, pre_p, passornot);
	
	if(passornot > 0)
	{
		PassVertornot = passornot;
		return ;
	}

	face_id = prev_face_id;  //////added on 06/08/05

	double param_t[2] = {0.};

	CrossBoundary3(pre_p, cur_p, face_id, alpha, which_edge, param_t);

	if(param_t[0] == -1 && param_t[1] == -1)
	{
		face_id = prev_face_id;   ////something wrong here
		return;
	}

	////if not passing a vertex, judge which triangle it will enter later
	PassEdge(face_id, which_edge);

	////calculate the intersection point on the boundary

}


/*--------------------------------------------------------------------------*/
extern bool IsRepeated(int *a, int b, int num);

void GetNextTriangle_new2(int &face_id, double pre_p[2], double cur_p[2], 
						  int *designcellcycle, int &num_cells)
{
	int which_edge = -1;

	int prev_face_id = face_id;

	Face *prev_face = Object.flist[face_id];

	double alpha[3] = {0.};

LL:	CrossBoundary(face_id, cur_p, which_edge);

	double t[2] = {0.};

	////Calculate the intersection, the result has been store in previous
	if(which_edge == 0)
		GetIntersection(pre_p, cur_p, prev_face->xy[1], prev_face->xy[2], t);
	else if(which_edge == 1)
		//GetIntersection(pre_p, cur_p, prev_face->xy[0], prev_face->xy[2], t);
		GetIntersection(pre_p, cur_p, prev_face->xy[2], prev_face->xy[0], t);
	else
		GetIntersection(pre_p, cur_p, prev_face->xy[0], prev_face->xy[1], t);
	

	if(t[0] == -1 && t[1] == -1)
		face_id = prev_face_id;   ////something wrong here


	////if not passing a vertex, judge which triangle it will enter later
	PassEdge(face_id, which_edge);

	////calculate the intersection point on the boundary
	pre_p[0] = pre_p[0] + t[0] * (cur_p[0] - pre_p[0]);
	pre_p[1] = pre_p[1] + t[0] * (cur_p[1] - pre_p[1]);

	Get2DBarycentricFacters(face_id, cur_p[0], cur_p[1], alpha);

	if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1 
	&& (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1 
	&& (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
	return; //if next point falls in this triangle, return

	if(!IsRepeated(designcellcycle, face_id, num_cells))
	{
		if(num_cells >= MaxNumTrianglesDesignCurve-1)
		{
			MaxNumTrianglesDesignCurve += 50;
			DesignCurveCellCycle = (int *)realloc(DesignCurveCellCycle, \
				sizeof(int)*MaxNumTrianglesDesignCurve);
		}

		designcellcycle[num_cells] = face_id;
		num_cells++;

		Object.flist[face_id]->inDesignCellCycle = 1;
	}

	goto LL;

}

/*----------------------------------------------------------------------------*/
////reuse the tracing method to get the next triangle through vertex 3/27/06

void TriangleThroughVertex_2(int vert_id, int &theone, icVector2 dir)
{
	Vertex *vert = Object.vlist[vert_id];

	////using a tricky way to get next triangle

   /*---------------------------------------------------------------------*/
	////Second method to judge the next triangle when the tracing curve
	////passes through a vertex

	int i;
	double vang;
	Face *face;
	Corner *c;
	int NewTID = -1;
    int orient;

	normalize(dir);

	vang = atan2(dir.entry[1], dir.entry[0]);

	//if(type == 1)
	//	vang += M_PI;

	if(vang < 0)
		vang += 2*M_PI;

	if(vang > 2*M_PI)
		vang -= 2*M_PI;

	////Get the orientation of the angle allocation
	orient = Object.clist[vert->Corners[0]]->orient;
	//orient = vert->oi;

	for( i = 0; i < vert->Num_corners; i++)
	{
		c = Object.clist[vert->Corners[i]];
		////first, we check which angle area the vector on the vertex will fall in
		if(orient > 0)
		{
			if(c->BeginAng > c->EndAng)
			{
				if(vang >= c->BeginAng || vang < c->EndAng)
				{
					NewTID = i;
					break;
				}
			}
			else{
				if(vang >= c->BeginAng && vang < c->EndAng)
				{
					NewTID = i;
					break;
				}
			}
		}
		else{
			if(c->BeginAng < c->EndAng)
			{
				if(vang <= c->BeginAng || vang > c->EndAng)
				{
					NewTID = i;
					break;
				}
			}
			else{
				if(vang <= c->BeginAng && vang > c->EndAng)
				{
					NewTID = i;
					break;
				}
			}
		}
	}

	if(NewTID == -1)  //reach boundary ?? 12/29/05
	{
		theone = theone;
		return;
	}

	face = Object.flist[Object.clist[vert->Corners[NewTID]]->t];

	theone = face->index;


}


void CrossVert_new(int &face_id, double cur_p[2], double pre_p[2],int &passornot)
{
	double vert[2];
	int alpha_index = 0;
	double max_alpha ;
    int newtriangleid = 0;
	int crossVert;

	Face *face = Object.flist[face_id];


	////New way to get the possible crossing vertex
	icVector2 test_dir;
	test_dir.entry[0] = cur_p[0] - face->xy[0][0];
	test_dir.entry[1] = cur_p[1] - face->xy[0][1];
	max_alpha = length(test_dir);
	alpha_index = 0;

	test_dir.entry[0] = cur_p[0] - face->xy[1][0];
	test_dir.entry[1] = cur_p[1] - face->xy[1][1];
	if(length(test_dir) < max_alpha)
	{
		max_alpha = length(test_dir);
	    alpha_index = 1;
	}

	test_dir.entry[0] = cur_p[0] - face->xy[2][0];
	test_dir.entry[1] = cur_p[1] - face->xy[2][1];
	if(length(test_dir) < max_alpha)
	{
		max_alpha = length(test_dir);
	    alpha_index = 2;
	}

	crossVert = face->verts[alpha_index];

	vert[0] = face->xy[alpha_index][0];
	vert[1] = face->xy[alpha_index][1];

	double A, B, C;
	A = pre_p[1] - cur_p[1];
	B = cur_p[0] - pre_p[0];
	C = (pre_p[0]*cur_p[1] - cur_p[0]*pre_p[1]);

	double pending = A*vert[0] + B*vert[1] + C;

	if(fabs(pending) < 1e-8) ////passing the vertex
	{
		double g_curp[2], g_prep[2];
		icVector2 gvec;
		gvec = pre_p[0]*face->LX + pre_p[1]*face->LY;
		g_prep[0] = Object.vlist[face->verts[0]]->x + gvec.entry[0];
		g_prep[1] = Object.vlist[face->verts[0]]->y + gvec.entry[1];

		gvec = cur_p[0]*face->LX + cur_p[1]*face->LY;
		g_curp[0] = Object.vlist[face->verts[0]]->x + gvec.entry[0];
		g_curp[1] = Object.vlist[face->verts[0]]->y + gvec.entry[1];
		//newtriangleid = TriangleDetect(g_curp[0], g_curp[1]);	//

		gvec.entry[0] = g_curp[0] - g_prep[0];
		gvec.entry[1] = g_curp[1] - g_prep[1];

        TriangleThroughVertex_2(crossVert, newtriangleid, gvec);
		passornot = alpha_index+1;
		return;
	}

	passornot = 0;
}



////Building the edges list consist of the cell cycle
void BuildBoundaryEdgeList(int *DesignCellCycle, int num_celltriangle)
{
	////search the direct former and later triangles of current triangle
	////find the sharing edges of these two pairs of triangle, mark the two edges of current triangle
	////the remaining edge must be on the boundary

	int i, j;
	Face *curf, *face1, *face2;
	Edge *cur_e;

	num_cycleedges = 0;

	for(i = 0; i < num_celltriangle; i++)
	{
		curf = Object.flist[DesignCellCycle[i]];

		for(j = 0; j < curf->nverts; j++)
		{
			cur_e = curf->edges[j];

			if(cur_e->tris[0] < 0 || cur_e->tris[1] < 0)
			{
				AddToBoundaryEdgeList(cur_e);
				cur_e->OnBoundary = 1;
				continue;
			}

			face1 = Object.flist[cur_e->tris[0]];
			face2 = Object.flist[cur_e->tris[1]];

			if((face1->inDesignCellCycle == 1 && face2->inDesignCellCycle != 1)\
			   ||(face1->inDesignCellCycle != 1 && face2->inDesignCellCycle == 1))
			{
				AddToBoundaryEdgeList(cur_e);
				cur_e->OnBoundary = 1;
			}
		}
	}
}

////More general routine for setting the boundary edges (not save them into another data structure!)
////4/13/06
void BuildBoundaryEdgeList(int *region, int num, int type)
{
	int i, j;
	Face *cur_f;
	Edge *cur_e;
	
	num_cycleedges = 0;

	for(i = 0; i < Object.nfaces; i++)
	{
		cur_f = Object.flist[i];

		for(j = 0; j < 3; j++)
		{
			cur_e = cur_f->edges[j];

			cur_e->BoundaryVisited = 0;
			cur_e->OnBoundary = 0;
		}
	}

	for(i = 0; i < num; i++)
	{
		cur_f = Object.flist[region[i]];

		for(j = 0; j < 3; j++)
		{
			cur_e = cur_f->edges[j];

			if(cur_e->tris[0] < 0 || cur_e->tris[1] < 0) //reach the domain/mesh boundary
			{
				//continue;
				AddToBoundaryEdgeList(cur_e);  //we need to add the edge at the mesh boundary 5/28/06
				cur_e->OnBoundary = 1;
			}

			////if both the adjacent faces of the edge are inside the region,  
			////it can't be a boundary edge

			if(IsRepeated(region, cur_e->tris[0], num) 
				&& IsRepeated(region, cur_e->tris[1], num))
				continue;

			if(cur_e->OnBoundary == 0)  //if the edge is not in the boundary list, add it
			{
				AddToBoundaryEdgeList(cur_e);
				cur_e->OnBoundary = 1;
			}
		}
	}
}

////New GetBoundary routine

void GetTwoBoundaries(Edge **edgelist, int numedges, int boundarytype)
{
	int i;
	Edge *cur_edge;
	int EndVertID;
	Vertex *cur_vert;

	if(boundarytype == 0) ////edgelist is the inner boundary
	{
		////Get inner boundary
		for(i = 0; i < numedges; i++)
		{
			InnerBoundary.edgelist[i] = edgelist[i];
		}
		InnerBoundary.num = numedges;

		////Get outer boundary
		for(i = 0; i < num_cycleedges; i++)
		{
			cur_edge = Cycle_edge[i];
			if(cur_edge->BoundaryVisited != 1)
				break;
		}
		 
		OuterBoundary.edgelist[0] = cur_edge;
		cur_edge->BoundaryVisited = 1;
		OuterBoundary.num = 1;
		EndVertID = cur_edge->verts[0];
	
		cur_vert = Object.vlist[cur_edge->verts[1]];

		////May be put into a routine
		while(cur_vert->VertID != EndVertID)   ////Not form a closed edges list
		{
			for(i = 0; i < cur_vert->Num_edge; i++)
			{
				cur_edge = cur_vert->edges[i];
				if(cur_edge->OnBoundary == 1 && cur_edge->BoundaryVisited != 1)
				{
					OuterBoundary.edgelist[OuterBoundary.num] = cur_edge;
					OuterBoundary.num ++;
					cur_edge->BoundaryVisited = 1;
					break;                         ////we just need to find one edge!
	 			}
			}
			
			////Get next testing vertex
			if(OuterBoundary.edgelist[OuterBoundary.num-1]->verts[0] != cur_vert->VertID)
				cur_vert = Object.vlist[OuterBoundary.edgelist[OuterBoundary.num-1]->verts[0]];
			else
				cur_vert = Object.vlist[OuterBoundary.edgelist[OuterBoundary.num-1]->verts[1]];
		}
	}

	else    ////edgelist is the outer boundary
	{         
		////Get outer boundary
		for(i = 0; i < numedges; i++)
		{
			OuterBoundary.edgelist[i] = edgelist[i];
		}
		OuterBoundary.num = numedges;

		////Get outer boundary
		for(i = 0; i < num_cycleedges; i++)
		{
			cur_edge = Cycle_edge[i];
			if(cur_edge->BoundaryVisited!= 1)
				break;
		}
		 
		InnerBoundary.edgelist[0] = cur_edge;
		cur_edge->BoundaryVisited = 1;
		InnerBoundary.num = 1;
		EndVertID = cur_edge->verts[0];
	
		cur_vert = Object.vlist[cur_edge->verts[1]];

		////May be put into a routine
		while(cur_vert->VertID != EndVertID)   ////Not form a closed edges list
		{
			for(i = 0; i < cur_vert->Num_edge; i++)
			{
				cur_edge = cur_vert->edges[i];
				if(cur_edge->OnBoundary == 1 && cur_edge->BoundaryVisited != 1)
				{
					InnerBoundary.edgelist[InnerBoundary.num] = cur_edge;
					InnerBoundary.num ++;
					cur_edge->BoundaryVisited = 1;
					break;                         ////we just need to find one edge!
	 			}
			}
			
			////Get next testing vertex
			if(InnerBoundary.edgelist[InnerBoundary.num-1]->verts[0] != cur_vert->VertID)
				cur_vert = Object.vlist[InnerBoundary.edgelist[InnerBoundary.num-1]->verts[0]];
			else
				cur_vert = Object.vlist[InnerBoundary.edgelist[InnerBoundary.num-1]->verts[1]];
		}
	}
}

//////////////////
void GetBoundary()
{
	////First we need to perform manifold testing and correct the original region
	int i;
	int EndVertID;
	Face *face;
	Vertex *cur_vert;
	Edge *cur_edge;

	Edge **temp_edgelist = (Edge **) malloc(sizeof(Edge *) * MaxEdgeInCycle);

	int temp_num_edges, other_num_edges;
	temp_num_edges = 0;

	////rebuild boundary again
	//for(i = 0; i < num_cycleedges; i++)
	//{
	//	Cycle_edge[i]->OnBoundary = 0;
	//}

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		for(int j = 0; j < 3; j++)
		{
			cur_edge = face->edges[j];
			cur_edge->OnBoundary = 0;
			cur_edge->BoundaryVisited = 0;
		}
	}

	BuildBoundaryEdgeList(DesignCurveCellCycle, num_triangles_designcurve);

	////Initial the flag of all the edges on current boundaries
	//for(i = 0; i < num_cycleedges; i++)
	//{
	//	Cycle_edge[i]->BoundaryVisited = 0;
	//}

	temp_num_edges = 0;
	temp_edgelist[0] = Cycle_edge[0];
	temp_num_edges++;
	EndVertID = Cycle_edge[0]->verts[0];
	temp_edgelist[0]->BoundaryVisited = 1;

	cur_vert = Object.vlist[Cycle_edge[0]->verts[1]];

	while(cur_vert->VertID != EndVertID)   ////Not form a closed edges list
	{
		for(i = 0; i < cur_vert->Num_edge; i++)
		{
			cur_edge = cur_vert->edges[i];
			if(cur_edge->OnBoundary == 1 && cur_edge->BoundaryVisited == 0)
			{
				temp_edgelist[temp_num_edges] = cur_edge;
				temp_num_edges ++;
				cur_edge->BoundaryVisited = 1;
	 		}
	    }
		
		////Get next testing vertex
		if(temp_edgelist[temp_num_edges-1]->verts[0] != cur_vert->VertID)
			cur_vert = Object.vlist[temp_edgelist[temp_num_edges-1]->verts[0]];
		else
			cur_vert = Object.vlist[temp_edgelist[temp_num_edges-1]->verts[1]];
    }

 	////Get inner and outer boundaries

    other_num_edges = num_cycleedges - temp_num_edges;

	if(other_num_edges > temp_num_edges)  ////temp_edgelist store the inner boundary
	{ 
	    GetTwoBoundaries(temp_edgelist, temp_num_edges, 0);	
	}
	else
	    GetTwoBoundaries(temp_edgelist, temp_num_edges, 1);	

	free(temp_edgelist);
}


/////////////
void GetNormalsForBoundaryEdges(int type)
{
	int i;
	Edge *cur_edge;
	Face *face;

    //////
	for(i = 0; i < num_cycleedges; i ++)
	{
		cur_edge = Cycle_edge[i];

		////Get the face that inside the region
		face = Object.flist[cur_edge->tris[0]];

		if(face->inDesignCellCycle != 1)
			face = Object.flist[cur_edge->tris[1]];

		CalNormalAtEdge(cur_edge, face, type);
	}
}


////Find the common vertex of two edges
int FindSharedVertex(Edge *edge1, Edge *edge2)
{
	int j, k;
	for(j = 0; j < 2; j++)
	{
		for(k = 0; k < 2; k++)
		{
			if(edge1->verts[j] == edge2->verts[k])
			{
				return (edge1->verts[j]);
			}
		}
	}
}

////Testing codes here 08/18/05
void GetVectorsOnBoundaries(int type)
{
	int i, j, k;
	Edge *cur_edge, *next_edge;
	Vertex *cur_vert = NULL;
	icVector2 vert_normal;
	icVector2 test1, test2;
	double theta;

	////Get the vectors on inner boundary
	for(i = 0; i < InnerBoundary.num - 1; i++)
	{
		cur_edge = InnerBoundary.edgelist[i];
		next_edge = InnerBoundary.edgelist[(i+1)%InnerBoundary.num];
		theta = (40./90.) * (M_PI / 2.);

		////find the common ending vertex of these two edges
		for(j = 0; j < 2; j++)
		{
			for(k = 0; k < 2; k++)
			{
				if(cur_edge->verts[j] == next_edge->verts[k])
				{
					cur_vert = Object.vlist[cur_edge->verts[j]];
					break;
				}
			}
		}

		////Calculate the vector on the common vertex according to the normals of the two adjacent edges
		//vert_normal = cur_edge->repell_normal + next_edge->repell_normal;
		vert_normal.entry[0] = cur_edge->normal.entry[0] + next_edge->normal.entry[0];
		vert_normal.entry[1] = cur_edge->normal.entry[1] + next_edge->normal.entry[1];
		normalize(vert_normal);
 

		//test1 = cur_edge->repell_normal;
		//test2 = next_edge->repell_normal;
		//normalize(test1);
		//normalize(test2);

		//if(dot(test1, test2) < 0.9)  ////this is a corner
		//{
		//	theta = M_PI / 2.;
		//}

		////rotate the normal counter clockwisely to get the vector
		if(type == 0)
		{
			cur_vert->vec.entry[0] = (cos(theta) * vert_normal.entry[0] - sin(theta) * vert_normal.entry[1])/200;
			cur_vert->vec.entry[1] = (sin(theta) * vert_normal.entry[0] + cos(theta) * vert_normal.entry[1])/200;
		}
		else
		{
			cur_vert->vec.entry[0] = (cos(theta) * vert_normal.entry[0] + sin(theta) * vert_normal.entry[1])/200;
			cur_vert->vec.entry[1] = (-sin(theta) * vert_normal.entry[0] + cos(theta) * vert_normal.entry[1])/200;
		}

		//////For correct output, You need to normalize the result vector

		cur_vert->OnBoundary = 1;
	}

	////Get the vectors on outer boundary
	for(i = 0; i < OuterBoundary.num - 1; i++)
	{
		cur_edge = OuterBoundary.edgelist[i];
		next_edge = OuterBoundary.edgelist[(i+1)%OuterBoundary.num];
		theta = (20./90.) * (M_PI / 2.);

		////find the common ending vertex of these two edges
		for(j = 0; j < 2; j++)
		{
			for(k = 0; k < 2; k++)
			{
				if(cur_edge->verts[j] == next_edge->verts[k])
				{
					cur_vert = Object.vlist[cur_edge->verts[j]];
					break;
				}
			}
		}

		////Calculate the vector on the common vertex according to the normals of the two adjacent edges
		//vert_normal = cur_edge->repell_normal + next_edge->repell_normal;
		//vert_normal.entry[0] = cur_edge->repell_normal.entry[0] + next_edge->repell_normal.entry[0];
		//vert_normal.entry[1] = cur_edge->repell_normal.entry[1] + next_edge->repell_normal.entry[1];
		vert_normal.entry[0] = cur_edge->normal.entry[0] + next_edge->normal.entry[0];
		vert_normal.entry[1] = cur_edge->normal.entry[1] + next_edge->normal.entry[1];
		normalize(vert_normal);
		

		//test1 = cur_edge->repell_normal;
		//test2 = next_edge->repell_normal;
		//normalize(test1);
		//normalize(test2);

		//if(dot(test1, test2) < 0.9)  ////this is a corner
		//{
		//	theta = M_PI / 2.;
		//}

		////rotate the normal clockwisely to get the vector
		if(type == 0)
		{
			cur_vert->vec.entry[0] = (cos(theta) * vert_normal.entry[0] + sin(theta) * vert_normal.entry[1])/200;
			cur_vert->vec.entry[1] = (-sin(theta) * vert_normal.entry[0] + cos(theta) * vert_normal.entry[1])/200;
		}
		else
		{
			cur_vert->vec.entry[0] = (cos(theta) * vert_normal.entry[0] - sin(theta) * vert_normal.entry[1])/200;
			cur_vert->vec.entry[1] = (sin(theta) * vert_normal.entry[0] + cos(theta) * vert_normal.entry[1])/200;
		}

		cur_vert->OnBoundary = 1;
	}
}


////Building the smoothing region by adding other vertices into the inner vertices list
////Note that we need to allocate index for those vertices
void BuildSmoothRegion()
{
	int i;
	Vertex *cur_vert;

	Num_verts = 0;

	for(i = 0; i < Object.nverts; i++)
	{
		cur_vert = Object.vlist[i];
		if(cur_vert->OnBoundary == 0)
		{
			regionverts[Num_verts] = cur_vert;
			cur_vert->RegionListID = Num_verts;
			cur_vert->InRegion = 1;
			Num_verts ++;
		}
	}
}


////Normalize the vectors on the vertices of the boundary
void NormalizeVectorsOnBoundary()
{
	int i, j;
	Edge *cur_edge;
	Vertex *cur_vert;
    double r;

	for(i = 0; i < num_cycleedges; i++)
	{
		cur_edge = Cycle_edge[i];

		for(j = 0; j < 2; j++)
		{
			cur_vert = Object.vlist[cur_edge->verts[0]];

			r = length(cur_vert->vec);
			r *= r;
						
			if (r < DistanceThreshold) 
			{
				r = DistanceThreshold;
				cur_vert->vec *= dmax/r; 
			}

			r = length(cur_vert->vec);
			r *= r;

			if (r > dmax*dmax) { 
				r  = sqrt(r); 
				cur_vert->vec *= dmax/r; 
			}
		}
	}
}



////Routines for new limit cycle shape design

////Using the hermite tangential vector calculation to get the tangential vector
icVector2 GetTangentialVec(Edge *edge1, Edge *edge2, int shared_vert)
{
	Vertex *vert1, *vert2;
	icVector2 Pstart, Pend, result;

	if(edge1->verts[0] != shared_vert)
		vert1 = Object.vlist[edge1->verts[0]];
	else
		vert1 = Object.vlist[edge1->verts[1]];

	if(edge2->verts[0] != shared_vert)
		vert2 = Object.vlist[edge2->verts[0]];
	else
		vert2 = Object.vlist[edge2->verts[1]];

	Pstart.entry[0] = vert1->x;
	Pstart.entry[1] = vert1->y;
	Pend.entry[0] = vert2->x;
	Pend.entry[1] = vert2->y;

	result = -GetTangent(Pend, Pstart);
	normalize(result);
	return result;
}


////Get the orientation of the two boundaries
void GetOrientationofBoundaries(int &same_orient, int &clockwise)
{
	Edge *edge1, *edge2;
	Vertex *cur_vert;
	icVector2 innerVec, outerVec;
	int index;
	icVector2 edge_normal, rotated_normal;

	////1. Calculate the tangential vectors on the first vertices of the two boundaries respectively

	////Inner boundary
	edge1 = InnerBoundary.edgelist[InnerBoundary.num-1];
	edge2 = InnerBoundary.edgelist[0];
	cur_vert = Object.vlist[FindSharedVertex(edge1, edge2)];
	innerVec = GetTangentialVec(edge1, edge2, cur_vert->VertID);

	////Outer boundary
	edge1 = OuterBoundary.edgelist[OuterBoundary.num-1];
	edge2 = OuterBoundary.edgelist[0];
	cur_vert = Object.vlist[FindSharedVertex(edge1, edge2)];
	outerVec = GetTangentialVec(edge1, edge2, cur_vert->VertID);

	////2. judge whether the two boundaries have the same orientation
	double dot_result = dot(innerVec, outerVec);

	if(dot_result < 0)  same_orient = 0;          ////Not the same
	else/* if(dot_result > 0) */ same_orient = 1;     ////Same orientation
	////we do not consider the dot_result == 0!! Now

	////3. Get the orientation of the inner boundary
	index = 0;
	edge1 = InnerBoundary.edgelist[index];
	edge_normal = edge1->repell_normal;
	rotated_normal.entry[0] = cos(M_PI)*edge_normal.entry[0] - sin(M_PI)*edge_normal.entry[1];
	rotated_normal.entry[1] = sin(M_PI)*edge_normal.entry[0] + cos(M_PI)*edge_normal.entry[1];

	dot_result = dot(innerVec, rotated_normal);

	if(dot_result > 0)
		clockwise = 1;                           ////clockwise orientation
	else/* if(dot_result < 0)*/
		clockwise = 0;                           ////counter clockwise orientation
	////we do not consider the dot_result == 0!! Now
}


void GetVectorsforAllBoundaryVerts(int same_orient, int clockwise, int type)
{
	////1. Calculate the vectors on inner boundary
	int i;
	Edge *cur_e, *next_e;
	Vertex *cur_v;
	icVector2 temp_vec;
	double theta;
	int sign = 1;

	////1. For inner boundary
	theta = (30./90.)*M_PI/2.;             ////Here we use constant rotation angle
	if(clockwise == 0)    
	{
		theta += M_PI;    ////rotate the tangential vector -(M_PI+theta)
	}

	for(i = 0; i < InnerBoundary.num; i++)
	{

		cur_e = InnerBoundary.edgelist[i];
		next_e = InnerBoundary.edgelist[(i+1)%InnerBoundary.num];

		cur_v = Object.vlist[FindSharedVertex(cur_e, next_e)];
		temp_vec = GetTangentialVec(cur_e, next_e, cur_v->VertID);

		if(type == 0)             ////want to get a repeller
		{
			cur_v->vec.entry[0] = (cos(theta) * temp_vec.entry[0] + sin(theta) * temp_vec.entry[1])/500;
			cur_v->vec.entry[1] = (-sin(theta) * temp_vec.entry[0] + cos(theta) * temp_vec.entry[1])/500;
		}
		else{                     ////want to get a repeller
			cur_v->vec.entry[0] = (cos(theta) * temp_vec.entry[0] - sin(theta) * temp_vec.entry[1])/500;
			cur_v->vec.entry[1] = (sin(theta) * temp_vec.entry[0] + cos(theta) * temp_vec.entry[1])/500;
		}
		cur_v->OnBoundary = 1;
	}

	////2. For outer boundary
	theta = (30./90.)*M_PI/2.;             ////Here we use constant rotation angle
	if(same_orient == 0)         ////not the same orientation
	{
		sign = -1;
	}
	
	if(clockwise == 0)    
	{
		theta += M_PI;    ////rotate the tangential vector -(M_PI+theta)
	}

	for(i = 0; i < OuterBoundary.num; i++)
	{

		cur_e = OuterBoundary.edgelist[i];
		next_e = OuterBoundary.edgelist[(i+1)%OuterBoundary.num];

		cur_v = Object.vlist[FindSharedVertex(cur_e, next_e)];
		temp_vec = GetTangentialVec(cur_e, next_e, cur_v->VertID);


		if(type == 0)             ////want to get a repeller
		{
			cur_v->vec.entry[0] = sign*(cos(theta) * temp_vec.entry[0] - sin(theta) * temp_vec.entry[1])/1000;
			cur_v->vec.entry[1] = sign*(sin(theta) * temp_vec.entry[0] + cos(theta) * temp_vec.entry[1])/1000;
		}
		else{                     ////want to get a repeller
			cur_v->vec.entry[0] = sign*(cos(theta) * temp_vec.entry[0] + sin(theta) * temp_vec.entry[1])/1000;
			cur_v->vec.entry[1] = sign*(-sin(theta) * temp_vec.entry[0] + cos(theta) * temp_vec.entry[1])/1000;
		}
		cur_v->OnBoundary = 1;
	}
}


////New method to calculate the vectors on the boundaries 08/26/05
void GetVectorsOnBoundaries_new(int type)
{
	int same_orient, clockwise;

	GetOrientationofBoundaries(same_orient, clockwise);

	GetVectorsforAllBoundaryVerts(same_orient, clockwise, type);
}




////The following two routins are used to store/set the boundary vertices and their corresponding vectors
////2/20/06
void SaveBoundaryVerts()
{
	int i;

	boundaryvertlist.num_boundarverts = 0;
	for(i = 0; i < Object.nverts; i++)
	{
		if(Object.vlist[i]->OnBoundary == 1)
		{
			if(boundaryvertlist.num_boundarverts >= MaxNumBoundaryVerts-1)
			{
				MaxNumBoundaryVerts += 100;
				boundaryvertlist.boundaryverts = (Boundaryvert*)realloc(boundaryvertlist.boundaryverts,
					sizeof(Boundaryvert)*MaxNumBoundaryVerts);
			}

			boundaryvertlist.boundaryverts[boundaryvertlist.num_boundarverts].vertID
			    = i;
			boundaryvertlist.boundaryverts[boundaryvertlist.num_boundarverts].vec
				= Object.vlist[i]->vec;
			normalize(boundaryvertlist.boundaryverts[boundaryvertlist.num_boundarverts].vec);
			boundaryvertlist.boundaryverts[boundaryvertlist.num_boundarverts].vec *= 0.005;
			boundaryvertlist.num_boundarverts++;
		}
	}
}



void SetBoundaryVerts()
{
	int i, vertID;

	for(i = 0; i < boundaryvertlist.num_boundarverts; i++)
	{

		vertID = boundaryvertlist.boundaryverts[i].vertID;
		Object.vlist[vertID]->vec = boundaryvertlist.boundaryverts[i].vec;
		Object.vlist[vertID]->OnBoundary = 1;
		Object.vlist[vertID]->RegionListID = -1;
	}
}


/*----------------------------------------------------------------------------------*/
////Routines for limit cycle shape design


void AllocShapeDesignVars()
{
	MaxNumShapeControlPts = 1000;

	control_pts = (ctr_point*)malloc(sizeof(ctr_point) * MaxNumShapeControlPts);
	
	CtrPts_tangent = (icVector2 *)malloc(sizeof(icVector2) * MaxNumShapeControlPts);

    out_pts = (ctr_point*)malloc(sizeof(ctr_point) * resolution);        // allocate our control point array

	num_shapecontrol_pts = 0;
	num_curvepts_output = 0;
	HermiteStep = 0;

    ////Temporary global variables for new limit cycle creation by shape design (modified at 3/25/06)
	MaxNumTrianglesDesignCurve = (int)Object.nfaces/5000;
    DesignCurveCellCycle = (int *) malloc(sizeof(int) * MaxNumTrianglesDesignCurve);
    num_triangles_designcurve = 0;

	////Allocate memory for the data structure for the new limit cycle creation method 5/20/06
	myCycle.MaxNumTrianglesDesignCurve = (int)Object.nfaces/4;
	myCycle.DesignCellCycle = (int *) malloc(sizeof(int) * myCycle.MaxNumTrianglesDesignCurve);
	myCycle.Appro_dir = (icVector2*) malloc(sizeof(icVector2) * myCycle.MaxNumTrianglesDesignCurve);
	myCycle.bases = (LocalPts*) malloc(sizeof(LocalPts) * myCycle.MaxNumTrianglesDesignCurve);
	/////////////////////////////////////////////////////////////////////////////////////////////

	MaxNumEdgesOnBoundary = (int)Object.nfaces/2;
	InnerBoundary.edgelist = (Edge **)malloc(sizeof(Edge *)*MaxNumEdgesOnBoundary);
	OuterBoundary.edgelist = (Edge **)malloc(sizeof(Edge *)*MaxNumEdgesOnBoundary);

    ////variable for shape design
	designcurve = (LineSeg*) malloc(sizeof(LineSeg*) * (int)(Object.nfaces/10));
    num_lineseg_designcurve = 0;

	////Allocate memory for boundary vertex list 2/20/06
	MaxNumBoundaryVerts = (int)Object.nverts/2000;
	boundaryvertlist.boundaryverts = (Boundaryvert *) malloc(sizeof(boundaryvertlist)*MaxNumBoundaryVerts);
	boundaryvertlist.num_boundarverts = 0;

	MaxNumBoundVerts = (int)Object.nverts/4;
	Boundaryverts = (int*)malloc(sizeof(int) * MaxNumBoundVerts);
	num_boundverts = 0;

}



void FinalizeShapeDesign()
{
	free(control_pts);
	free(CtrPts_tangent);
	free(out_pts);
	free(DesignCurveCellCycle);
	free(InnerBoundary.edgelist);
	free(OuterBoundary.edgelist);
	free(designcurve);

	free(Boundaryverts);
	free(boundaryvertlist.boundaryverts);

	free(myCycle.DesignCellCycle);
	free(myCycle.Appro_dir);
	free(myCycle.bases);
}


////Add the user selected point into the points list
void AddToShapeCtrPtsList(double x, double y)
{
	//int vert;

	if(num_shapecontrol_pts >= MaxNumShapeControlPts - 1)
	{
		MessageBox(NULL, "Can not add more control points!", "Error", MB_OK);
		return;
	}
	control_pts[num_shapecontrol_pts].x = x;
	control_pts[num_shapecontrol_pts].y = y;
	//control_pts[num_shapecontrol_pts].z = 0;


	num_shapecontrol_pts ++;
}


void CalHermiteCurve()
{
	if(num_shapecontrol_pts < 2)
	{
		MessageBox(NULL, "Not enough interpolated points", "Error", MB_OK);
		return;
	}

	//resolution = 400; /*10/10/2007*/

	HermiteStep = (int)(resolution / num_shapecontrol_pts);

	Hermitecurve(num_shapecontrol_pts, control_pts, CtrPts_tangent, out_pts, HermiteStep, num_curvepts_output);
	resolution = num_curvepts_output;

}

void CalOpenHermiteCurve()
{
	if(num_shapecontrol_pts < 2)
	{
		MessageBox(NULL, "Not enough interpolated points", "Error", MB_OK);
		return;
	}

	HermiteStep = (int)(resolution / num_shapecontrol_pts);

	Hermitecurve_open(num_shapecontrol_pts, control_pts, CtrPts_tangent, out_pts, HermiteStep, num_curvepts_output);
	resolution = num_curvepts_output;
}


////According to the control point to generate the limit cycle elements
////Using the tangent vectors stored for the control points
////Testing codes 08/13/05, we generate limit cycle using divergent element here
void GenerateLimitCycle(int type)
{
	////
   
    int i;

	for(i = 0; i < num_shapecontrol_pts; i++)
	{
		////1. Set basis 
		if(cur_regelem_index >= MaxNumRegularElems -1)
		{
			//Allocate new space for the element list
			MaxNumRegularElems += 50;
			regularelem = (RegularElement*)realloc(regularelem, sizeof(RegularElement) * MaxNumRegularElems);
		}
		
		regularelem[cur_regelem_index].base[0] = control_pts[i].x;
		regularelem[cur_regelem_index].base[1] = control_pts[i].y;
		regularelem[cur_regelem_index].ID = NAMEOFREGELEM + cur_regelem_index;
		regularelem[cur_regelem_index].type = type;   

		////Initialize the transformation parameters for this regular element
		regularelem[cur_regelem_index].transform_matrix.setIdentity();
		regularelem[cur_regelem_index].transposeRot.setIdentity();

		regularelem[cur_regelem_index].rotang = 0;
		regularelem[cur_regelem_index].s = 1;


		////2. Set the direction
		regularelem[cur_regelem_index].Direct = CtrPts_tangent[i];

		////Normalize the regular vector and set it to the default strength

		normalize(regularelem[cur_regelem_index].Direct);

		regularelem[cur_regelem_index].Direct.entry[0] *= RegularStrength;
		regularelem[cur_regelem_index].Direct.entry[1] *= RegularStrength;
		
		//You might need to update the transform matrix for convergent and divergent elements
		if(regularelem[cur_regelem_index].type > 0)
			SetTransform(regularelem[cur_regelem_index].base,\
				regularelem[cur_regelem_index].Direct,\
				regularelem[cur_regelem_index].transform_matrix,\
				regularelem[cur_regelem_index].transposeRot);
		
		cur_regelem_index ++;
	}
}


///////////////////////////////////////////////////////////

void LimitCycleGeneration(int type)
{
	////save the indices of boundary vertices
	SaveBoundaryVerts();

	//GetDesignCellCycle_new();
	GetDesignCellCycle_new3();

	if(CalEulerValue(DesignCurveCellCycle, num_triangles_designcurve) != 0)
	{
		MessageBox(NULL, "not form a ring-shaped region, try it again!", "Error", MB_OK);
		num_shapecontrol_pts = 0;
	    num_curvepts_output = 0;

		num_cycleedges = 0;
		num_triangles_designcurve = 0;
		num_innertriangles = 0;
		return;
	}



	////Get the vertices on current boundary and extend the cell cycle from these vertices
	ExtendDesignCellCycle();
	if(CalEulerValue(DesignCurveCellCycle, num_triangles_designcurve) != 0)
	{
		MessageBox(NULL, "not form a ring-shaped region, try it again!", "Error", MB_OK);
		num_shapecontrol_pts = 0;
	    num_curvepts_output = 0;

		num_cycleedges = 0;
		num_triangles_designcurve = 0;
		num_innertriangles = 0;
		return;
	}

    GetBoundary();
	GetNormalsForBoundaryEdges(type);

	GetVectorsOnBoundaries(type);
    //GetVectorsOnBoundaries_new(0);
	
	SetBoundaryVerts();  ////set the previous boundaries (i.e. previous limit cycles)

	BuildSmoothRegion();

	RegionSmooth();

	////You still need to normalize the vectors on the vertices of the boundary
    //NormalizeVectorsOnBoundary();
	 NormalizeField();

}

/*-------------------------------------------------------------*/
////Save the control points for one limit cycle
////3/2/06
void SaveDesignCurve()
{
	int i;
	FILE *fp = fopen("designcurve.txt", "w");
	
	////write control points
	fprintf(fp, "Control points:%d\n", num_shapecontrol_pts);
	fprintf(fp, "type:%d\n", limitcycletype);
	for(i = 0; i < num_shapecontrol_pts; i++)
	{
		fprintf(fp, "%f,%f\n", control_pts[i].x, control_pts[i].y);
	}
	fclose(fp);

}

////Load the curve from the file
void LoadDesignCurve()
{
	int i;
	FILE *fp = fopen("designcurve.txt", "r");

	if(fp == NULL)
	{
		MessageBox(NULL, "Can not open file!", "Error", MB_OK);
		return;
	}

	//int num_shapecontrol_pts;
	float x, y;

	////write control points
	fscanf(fp, "Control points:%d\n", &num_shapecontrol_pts);
	fscanf(fp, "type:%d\n", &limitcycletype);
	for(i = 0; i < num_shapecontrol_pts; i++)
	{
		fscanf(fp, "%f,%f\n", &x, &y);
		control_pts[i].x = x;
		control_pts[i].y = y;
	}
	fclose(fp);
}


/*
Get an adaptive resolution for the design curve 3/27/06
*/
int GetResolution()
{
	double edgelen;
	icVector2 edgevec;
	edgevec.entry[0] = Object.vlist[Object.flist[0]->verts[1]]->x-Object.vlist[Object.flist[0]->verts[0]]->x;
	edgevec.entry[1] = Object.vlist[Object.flist[0]->verts[1]]->y-Object.vlist[Object.flist[0]->verts[0]]->y;
	edgelen = length(edgevec);

	int resolute = (int) 4* GetWholeLengthofDesignCurve()/edgelen;

	return resolute;
}

double GetWholeLengthofDesignCurve()
{
	double whole_len = 0;
	icVector2 len_vec;

	for(int i = 0; i < num_shapecontrol_pts; i++)
	{
		len_vec.entry[0] = control_pts[(i+1)%num_shapecontrol_pts].x - control_pts[i].x;
		len_vec.entry[1] = control_pts[(i+1)%num_shapecontrol_pts].y - control_pts[i].y;

		whole_len += length(len_vec);
	}

	return whole_len;
}



void ResetCurBoundary()
{
	int i;
	Vertex *vert;

	for(i = 0; i < num_cycleedges; i++)
	{
		vert = Object.vlist[Cycle_edge[i]->verts[0]];
		vert->OnBoundary = 0;

		vert = Object.vlist[Cycle_edge[i]->verts[1]];
		vert->OnBoundary = 0;
	}
}

/*
Judge whether a triangle contains a boundary vertex
*/
bool InBoundary(int triangleid) //3/29/06
{
	int i;
	Face *face = Object.flist[triangleid];

	for(i = 0; i < face->nverts; i++)
	{
		if(Object.vlist[face->verts[i]]->OnBoundary == 1)
			return true;
	}

	return false;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
void InitLimitCycleShapeDesign()
{
	ResetLimitCycleShapeDesign();
	InitUnderneathMesh();
}



void ResetLimitCycleShapeDesign()
{
	////Initialize
	num_shapecontrol_pts = 0;
	num_curvepts_output = 0;

	num_cycleedges = 0;
	num_triangles_designcurve = 0;

	InnerBoundary.num = 0;
	OuterBoundary.num = 0;

	num_boundverts = 0;
	boundaryvertlist.num_boundarverts = 0;

	myCycle.num_triangles_designcurve = 0;
}

/*---------------------------------------------------------------------------------------*/
///////////////////////////////////////////////////////////////////////////////////////////
////// Our brand new limit cycle creation method (new method for setting boundary vertices


icVector2 GlobalToLocal_Vector(icVector2 glob_vec, int triangle)
{
	Face *face = Object.flist[triangle];

	icVector2 loc_vec;

	icVector2 lx = face->LX;
	icVector2 ly = face->LY;

	normalize(lx);
	normalize(ly);

	loc_vec.entry[0] = lx.entry[0]*glob_vec.entry[0] + lx.entry[1]*glob_vec.entry[1];
	loc_vec.entry[1] = ly.entry[0]*glob_vec.entry[0] + ly.entry[1]*glob_vec.entry[1];
	
	return loc_vec;
}


void GlobalToLocal_Pt(double gp[2], int triangle, double lp[2])
{
	Face *face = Object.flist[triangle];

	icVector2 VP;

	VP.entry[0] = gp[0] - Object.vlist[face->verts[0]]->x;
	VP.entry[1] = gp[1] - Object.vlist[face->verts[0]]->y;

	lp[0] = dot(VP, face->LX);
	lp[1] = dot(VP, face->LY);
}


void GetDesignCellCycle_new4()
{
	int i;
	int pre_triangle, pre_pos;
	int cur_triangle = TriangleDetect(out_pts[0].x, out_pts[0].y);
	double globalp[2] = {0.};
	double pre_globalp[2], loc_base[2] = {0.};
	icVector2 glob_vec, loc_vec;
	int Passvertornot;

    ////initialize part
	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->inDesignCellCycle = 0;
		Object.flist[i]->discard = 0;
	}

	for(i = 0; i < myCycle.MaxNumTrianglesDesignCurve; i++)
		myCycle.Appro_dir[i].set(0, 0);

	myCycle.num_triangles_designcurve = 0;
	myCycle.DesignCellCycle[myCycle.num_triangles_designcurve] = cur_triangle;
	myCycle.num_triangles_designcurve++;
	Object.flist[cur_triangle]->inDesignCellCycle = 1;
	myCycle.num_triangles_designcurve = 1;

	pre_globalp[0] = globalp[0] = out_pts[0].x;
	pre_globalp[1] = globalp[1] = out_pts[1].y;
	pre_triangle = cur_triangle;

	////find the cell cycle (closed triangle strip) containing the design curve
	for(i = 1; i < resolution; )
	{
		Passvertornot=0;

		pre_pos = i;

		cur_triangle = GetOneTriangle2(i, globalp, cur_triangle, resolution,
			myCycle.DesignCellCycle, myCycle.num_triangles_designcurve, Passvertornot);

		if(cur_triangle == -1)
			return;

		if(Passvertornot>0)
		{
			///jitter the position of current point 5/22/06
			out_pts[i].x += pow(-1., rand()%2)*(double)rand()/((double)RAND_MAX + 1)*Object.radius/sqrt((double)Object.nfaces);
			out_pts[i].y += pow(-1., rand()%2)*(double)rand()/((double)RAND_MAX + 1)*Object.radius/sqrt((double)Object.nfaces);

			i = pre_pos;
			continue;
		}

		
		if(cur_triangle == myCycle.DesignCellCycle[0] && i > 1) //We probably get a closed cycle!
		{
			//set the possible last triangle
			//Set the approximate vector and basis for vector calculation
			glob_vec.entry[0] = globalp[0] - pre_globalp[0];
			glob_vec.entry[1] = globalp[1] - pre_globalp[1];

			//set vector
			//loc_vec = GlobalToLocal_Vector(glob_vec, pre_triangle);
			loc_vec = GlobalToLocal_Vector(glob_vec, 
				myCycle.DesignCellCycle[myCycle.num_triangles_designcurve-1]);
			myCycle.Appro_dir[myCycle.num_triangles_designcurve-1] = loc_vec;

			//set basis
			//GlobalToLocal_Pt(pre_globalp, pre_triangle, loc_base);
			GlobalToLocal_Pt(pre_globalp, myCycle.DesignCellCycle[myCycle.num_triangles_designcurve-1], loc_base);
			myCycle.bases[myCycle.num_triangles_designcurve-1].x[0] = loc_base[0];
			myCycle.bases[myCycle.num_triangles_designcurve-1].x[1] = loc_base[1];
		}
		
		//if the design curve cross a vertex, we need to make sure that the new added triangle has
		//one edge share with previous one, otherwise, we need to add more triangle to connect them
		//if(Passvertornot>0)
		//{
		//	//Vertex *SingleSharingVert = Object.vlist[Object.flist[cur_triangle]->verts[Passvertornot - 1]];

		//	//////If the two triangles do not share the same vertex, we need to add more segments before
		//	//////we perform following operations 2/5/06
		//	//////because we have cell cycle extension process, so it is not very important issue
		//	//////but if these two triangles do not share the same vertex, we do not want to call the following routine

		//	//if(ShareTheSameVertex(DesignCurveCellCycle[num_triangles_designcurve-1], cur_triangle, SingleSharingVert->VertID))
		//	//{
		//	//	ConnectTwoTriangles2(DesignCurveCellCycle[num_triangles_designcurve-1], cur_triangle, SingleSharingVert->VertID,
		//	//		myCycle.DesignCellCycle, myCycle.num_triangles_designcurve);
		//	//}

		//	///jitter the position of current point 5/22/06
		//}


		//if it is new triangle, add to the triangle list
		if(!IsRepeated(myCycle.DesignCellCycle, cur_triangle, myCycle.num_triangles_designcurve)
			&& !(InBoundary(cur_triangle)) && Object.flist[cur_triangle]->discard != 1) 
		{
			if(myCycle.num_triangles_designcurve >= myCycle.MaxNumTrianglesDesignCurve-1)
			{
				myCycle.MaxNumTrianglesDesignCurve += 50;
				myCycle.DesignCellCycle = (int *)realloc(myCycle.DesignCellCycle, \
					sizeof(int)*myCycle.MaxNumTrianglesDesignCurve);
			}
			myCycle.DesignCellCycle[myCycle.num_triangles_designcurve] = cur_triangle;

			//Set the approximate vector and basis for vector calculation
			glob_vec.entry[0] = globalp[0] - pre_globalp[0];
			glob_vec.entry[1] = globalp[1] - pre_globalp[1];

			//set vector
			//loc_vec = GlobalToLocal_Vector(glob_vec, pre_triangle);
			loc_vec = GlobalToLocal_Vector(glob_vec, 
				myCycle.DesignCellCycle[myCycle.num_triangles_designcurve-1]);
			myCycle.Appro_dir[myCycle.num_triangles_designcurve-1] = loc_vec;

			//set basis
			//GlobalToLocal_Pt(pre_globalp, pre_triangle, loc_base);
			GlobalToLocal_Pt(pre_globalp, myCycle.DesignCellCycle[myCycle.num_triangles_designcurve-1], loc_base);
			myCycle.bases[myCycle.num_triangles_designcurve-1].x[0] = loc_base[0];
			myCycle.bases[myCycle.num_triangles_designcurve-1].x[1] = loc_base[1];

			//

			myCycle.num_triangles_designcurve++;
			Object.flist[cur_triangle]->inDesignCellCycle = 1;

			Object.flist[cur_triangle]->discard = 1; //avoid to include other cell belonging to other cycle
		}

		pre_globalp[0] = globalp[0];
		pre_globalp[1] = globalp[1];

		pre_triangle = cur_triangle;
	}
}

//////////////////
void GetBoundaryVerts(int *DesignCellCycle, int num_triangles_designcurve)
{
	////First we need to perform manifold testing and correct the original region
	int i;
	//int EndVertID;
	Face *face;
	//Vertex *cur_vert;
	Edge *cur_edge;

	Edge **temp_edgelist = (Edge **) malloc(sizeof(Edge *) * MaxEdgeInCycle);

	int temp_num_edges/*, other_num_edges*/;
	temp_num_edges = 0;

	////rebuild boundary again
	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		for(int j = 0; j < 3; j++)
		{
			cur_edge = face->edges[j];
			cur_edge->OnBoundary = 0;
			cur_edge->BoundaryVisited = 0;
		}
	}

	BuildBoundaryEdgeList(DesignCellCycle, num_triangles_designcurve);

	num_boundverts = 0;

	for(i = 0; i < num_cycleedges; i++)
	{
		cur_edge = Cycle_edge[i];
		
		if(!IsRepeated(Boundaryverts, cur_edge->verts[0], num_boundverts))
		{
			Boundaryverts[num_boundverts] = cur_edge->verts[0];
			num_boundverts ++;
		}
		
		if(!IsRepeated(Boundaryverts, cur_edge->verts[1], num_boundverts))
		{
			Boundaryverts[num_boundverts] = cur_edge->verts[1];
			num_boundverts ++;
		}
	}

}

////
void GetVectorsOnBoundaries(int *Boundaryverts, int numverts, int type)
{
	int i, j;
	Face *face;
	Vertex *cur_vert = NULL;
	Corner *c;
	int position;
	int num_vecs;
	icVector2 VP, loc_vec, glob_vec;
	double a[2] = {0.};

	for(i = 0; i < numverts; i++)
	{
		num_vecs = 0;
		glob_vec.set(0, 0);
		
		//get current vertex
		cur_vert = Object.vlist[Boundaryverts[i]];

		//Calculate the average vector for each boundary vertex
		for(j = 0; j < cur_vert->Num_corners; j++)
		{
			c = Object.clist[cur_vert->Corners[j]];

			if(TriangleSearch(myCycle.DesignCellCycle, myCycle.num_triangles_designcurve, c->t, position))
			{
				face = Object.flist[myCycle.DesignCellCycle[position]];

				VP.entry[0] = cur_vert->x - Object.vlist[face->verts[0]]->x;
				VP.entry[1] = cur_vert->y - Object.vlist[face->verts[0]]->y;

				a[0] = dot(VP, face->LX);
				a[1] = dot(VP, face->LY);

				loc_vec = GetALocalVec(a, myCycle.Appro_dir[position], myCycle.bases[position].x, type);

				VP = loc_vec.entry[0]*face->LX + loc_vec.entry[1]*face->LY;

				glob_vec = glob_vec + VP;

				num_vecs++;
			}
		}

		//Set the vector as the average of the vectors obtained above
		cur_vert->vec.entry[0] = 2*glob_vec.entry[0]/num_vecs;
		cur_vert->vec.entry[1] = 2*glob_vec.entry[1]/num_vecs;

	    cur_vert->OnBoundary = 1;

		/*save the vector values before normalization 02/21/07 */
		cur_vert->vec_J = cur_vert->vec;

	}

}



void LimitCycleGeneration2(int type) throw(...)
{
	////save the indices of boundary vertices
	SaveBoundaryVerts();

	GetDesignCellCycle_new4();

	if(CalEulerValue(myCycle.DesignCellCycle, myCycle.num_triangles_designcurve) != 0)
	{
		MessageBox(NULL, "not form a ring-shaped region, try it again!", "Error", MB_OK);
		num_shapecontrol_pts = 0;
	    num_curvepts_output = 0;

		num_cycleedges = 0;
		myCycle.num_triangles_designcurve = 0;
		num_innertriangles = 0;
		return;
	}


    GetBoundaryVerts(myCycle.DesignCellCycle, myCycle.num_triangles_designcurve);

	GetVectorsOnBoundaries(Boundaryverts, num_boundverts, type);

	//////////////Testing
	test_num_cycleedges = num_cycleedges;
	//////////////////////
	
	SetBoundaryVerts();  ////set the previous boundaries (i.e. previous limit cycles)

	BuildSmoothRegion();

	RegionSmooth();

	////You still need to normalize the vectors on the vertices of the boundary
	NormalizeField();

}

bool CreateALimitCycle(/*int *trianglelist, int num*/)
{
	CopyRegion(myCycle.DesignCellCycle,intersectRegion.trianglelist ,myCycle.num_triangles_designcurve);
	intersectRegion.num = myCycle.num_triangles_designcurve;

	if(IntersectedRegionSingCapture())
	{
		////We need to undo the 
		Undo();
		UnsetTriangleList(myCycle.DesignCellCycle, myCycle.num_triangles_designcurve);
		return false;
	}

	return true;
}


void UnsetTriangleList(int *trianglelist, int num)
{
	int i, j;
	Face *face;
	Vertex *vert;

	for(i = 0; i < num; i++)
	{
		face = Object.flist[trianglelist[i]];
		face->inDesignCellCycle = 0;

		for(j = 0; j < 3; j++)
		{
			vert = Object.vlist[face->verts[j]];
			vert->OnBoundary = 0;
		}
	}
}

double POCurl_Coeff = 40;


////This will return the vector under current local frame
icVector2 GetALocalVec(double loc_p[2], icVector2 loc_vec, double loc_base[2], int type)
{
	double dx = loc_p[0] - loc_base[0];
	double dy = loc_p[1] - loc_base[1];

	double r = dx*dx + dy*dy;

	double theta = atan2(loc_vec.entry[1], loc_vec.entry[0]);

	double t_vx = cos(theta)*dx + sin(theta)*dy;
	double t_vy = -sin(theta)*dx + cos(theta)*dy;

	if(type == 0) //repeller
	{
		//t_vx = -1./(40*r*sqrt(1+t_vy*t_vy));
		//t_vy = t_vy/(r*sqrt(1+t_vy*t_vy));
					t_vx = 1./r;
					t_vy = t_vy*POCurl_Coeff/r;
	}
	else //attractor
	{
		//t_vx = 1./(40*r*sqrt(1+t_vy*t_vy));
		//t_vy = -t_vy/(r*sqrt(1+t_vy*t_vy));
					t_vx = 1./r;
					t_vy = -t_vy*POCurl_Coeff/r;
	}

	double vx = (cos(theta)*t_vx - sin(theta)*t_vy)/r;
	double vy = (sin(theta)*t_vx + cos(theta)*t_vy)/r;

	//normalize the vector
    r = vx*vx + vy*vy;
			
	if (r > dmax*dmax) { 
		r  = sqrt(r); 
		vx *= dmax/r; 
		vy *= dmax/r; 
	}

	return icVector2(vx, vy);
}





