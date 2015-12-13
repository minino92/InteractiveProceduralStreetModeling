

/*
Taustep.cpp

created by Guoning Chen 01/10/07
*/

#include "stdafx.h"

#include "FindSCC.h"

#include "RegionSmoothing.h"

#include "LocalTracing.h"

#include "lib/icVector.h"

#include "VFDataStructure.h"

#include "Taumap.h"

extern Polygon3D Object;
extern GraphEdge *sccedges;
extern int num_sccedges;
extern int curMaxNumDirGraphEdges;
extern GraphNode2 *sccnodes;
extern int num_sccnodes;

int pre_num_sccedges;  /*record the previous number of sccedges for the adaptive framework*/



extern int cur_traj_index;
extern int *num_linesegs_curtraj;

extern Point *point;
extern int Num_SmoothRegionpoints;
extern Vertex **regionverts;                ////mesh vertices inside user selected region
extern int Num_verts;
extern int Num_edges;

extern MCGNode *mcgnodes ;
extern MCGEdge *mcgedges ;
extern int cur_mcgnode_index;
extern int num_sccomps;               //record the current number of SCC components

extern int cur_end_edgelist;                 /*the last index of the edge in the list*/


int *tri_strip = NULL;
int ntris_in_strip;
int curMaxTrisInStrip = 20;

/* new triangle T' */
icVector2 newP[3];
int en_tris[3];

double g_dt; /*global variable for storing the integration time */
double hstep;
extern int globalface; /*using tracing here*/
bool usespatialtau = false;
double trace_time;

int num_toolarge_regions = 0;


double newPx[7], newPy[7]; /*save the global coords only here! 02/26/07*/
extern int SelectTriangleID;

extern time_t g_rawtime;


/*------------------07/05/2007------------------*/
/*Important global variable for spatial tau!!!
*/
extern icVector2 g_cur_vec;
extern double g_vec_mag;


/*------------------08/01/2007------------------*/

TraceSamplePtList *_forward, *backward_spts;


extern bool IsRepeated(int *, int, int);
extern void CaptureSing();

/* Declaration of the fuctions */

void trace_Mesh_Tau(double);
void trace_Tri_Tau(int, double,int);
void trace_all_Verts(double tau, int backward);
void trace_Ver(int tri, double st[2], double en[2], int &en_tri, double, int);

/*new idea of building the $\tau$ map 03/19/07*/

void construct_tau_map(double tau);
void trace_all_Vers(double tau, int backward);
void trace_one_Ver(double tau, int backward);
void build_edges_Ver(int verid);
void build_edges_all_Vers(int backward);
void construct_tau_map(double tau);
void trace_all_edges_build_edges(double tau, int backward);
void trace_center_tris_build_edge(int tri, double tau, int backward);
void trace_all_centers_tris_build_edges(double tau, int backward);

/*note that this is only a test routine 3/20/07*/
void trace_samples_all_edges_build_edges(double tau, int backward, int nsamples);

/*adaptively sampling along the edge*/
void trace_recursive_an_edge(double v1[2], double v2[2], int &t1, int tri, int neighbor_tri,
							 double tau, int backward, int &level);
void trace_recursive_an_edge_interpolate(double v1[2], double v2[2], int &t1, int tri, int neighbor_tri,
							 double &tau1, double &tau2, int backward, int &level);

void trace_all_edges_build_di_edges_adp(double tau, int backward);
void trace_an_edge_build_di_edges_adp(double st1[2], double st2[2], int t1, int t2, int tri,
									  int neighbor_tri, double tau, int backward, int edge_id);
void trace_an_edge_build_di_edges_adp_interpolate(double st1[2], double st2[2], int t1, int t2, int tri,
									  int neighbor_tri, double tau1, double tau2, int backward);


int trace_in_triangle_tau(int &face_id, double globalp[2], int type, double tau, int &flag);


/*for testing the image evolution of a selected triangle*/
void trace_Edge(Edge *e, double);  //02/11/07
void trace_Boundary(int backward, double);

void add_Extra_Edges_Tau(int tri, double tau);
bool has_Edge_From_To(int tri1, int tri2);

void get_A_Tri_Strip(double [], int, double [], int);

void get_Tri_Strip();

void get_Next_Neighbor_Tri(int &face_id, double pre[2], double cur[2],
					 int &PassVertornot, double alpha[3], icVector2 line_dir);

void test_Cross_Vertex(int &face_id, double cur_p[2], 
					   double pre_p[2],int &passornot, icVector2 line_dir);

void get_Covered_Tris();

void get_Boundary_2(int *tri_strip, int ntris_in_strip);

void reset_Edge_Flags();

void reset_Triangle_Flags();

void init_All_Edges();
void init_all_Vers(); /*for adaptive framework*/

void trace_for_A_Strip(double p1[2], int t1, double p2[2], int t2);

int get_Type_of_Point(double alpha[3], int &which_vertex, int &which_edge);

void add_To_Triangle_Strip(int tri_id); /* We assume the default strip is tri_strip */

void trace_and_save(int tri, double tau);

Edge *find_Edge(int tri, int which_edge);
Vertex *find_Vertex(int tri, int which_vertex);
void fill_Gap_In_Ring(Edge **edgelist, int n_edges, int *triangles, int ntris);

GraphEdge *find_unused_edge(GraphEdge *oneunused, int &index, bool &flag);


bool contain_separation_feature(int tri)
{
	int i;
	Edge *e;
	Face *f = Object.flist[tri];

	for(i = 0; i < 3; i++)
	{
		e = f->edges[i];

		if(e->find_sep > 0)
			return true;
	}

	return false;
}

bool contain_attachment_feature(int tri)
{
	int i;
	Edge *e;
	Face *f = Object.flist[tri];

	for(i = 0; i < 3; i++)
	{
		e = f->edges[i];

		if(e->find_attp > 0)
			return true;
	}

	return false;
}

/* Implementation of the functions */

void trace_Ver(int tri, double st[2], double en[2], int &en_tri, double tau, int backward)
{
	int i;
	int flag = 0;
	double globalp[2];
	int pre_face, cur_face;


	pre_face = cur_face = tri;

	num_linesegs_curtraj[cur_traj_index] = 0;


	globalp[0] = st[0];
	globalp[1] = st[1];

	/* Reset the integration time */
	g_dt = 0;
	trace_time = 0;

	int steps = (int)(tau*20);

	for(i = 0; i < 5*NUMTRACINGTRIANGLE/*&& num_linesegs_curtraj[cur_traj_index]< steps */&& g_dt < tau; i++)
	{

		if(cur_face == -1)
		{
			en[0] = globalp[0];
			en[1] = globalp[1];
			en_tri = -1; /* changed at 02/25/07 */
			return;
		}

		/*if the triangle containing attachement points for backward tracing
		or separation points for _forward tracing, stop it 03/15/07
		This test has been proven not working well 03/19/07*/
		//if(backward == 0)
		//{
		//	if(contain_separation_feature(cur_face))
		//	{
		//		en[0] = globalp[0];
		//		en[1] = globalp[1];
		//		en_tri = cur_face;
		//		return;
		//	}
		//}
		//else
		//{
		//	if(contain_attachment_feature(cur_face))
		//	{
		//		en[0] = globalp[0];
		//		en[1] = globalp[1];
		//		en_tri = cur_face;
		//		return;
		//	}
		//}

		pre_face = cur_face;

		/* since we only need the last point, there is no need to save all the integrated points*/
		//cur_face = TraceInATriangle2(cur_face, globalp, backward, flag); 
		cur_face = trace_in_triangle_tau(cur_face, globalp, backward, tau, flag); 

		/*recompute the g_dt according to the intersection!!!
		we take care of it by using the real length of each line segment!*/
		//g_dt -= hstep;
		
		//if(flag == 3 || flag == 4 /*|| pre_face == cur_face*/ ) /*the accurate tracing*/
		if(flag == 3 || flag == 4 || pre_face == cur_face ) /*the approximate tracing*/
		{
			en[0] = globalp[0];
			en[1] = globalp[1];
			en_tri = cur_face;
			return;
		}
	}

	en[0] = globalp[0];
	en[1] = globalp[1];
	en_tri = cur_face;
}


/*****************************************************************************
Only trace in a triangle, not store the tracing points
For limit cycle detection, we may need to store two recently calculated points
Modified again on 02/20/07
*****************************************************************************/
int trace_in_triangle_tau(int &face_id, double globalp[2], int type, double tau, int &flag)
{
	int i;
	double alpha[3];
	double cur_point[2], pre_point[2];
	double vert0[2];
	icVector2 VP, globalv, line_length;

	if(face_id < 0)
		return -1;

	Face *face = Object.flist[face_id];

	Face *pre_f = face;
	
	////Temporary curve point array

	////initialize
	VP.entry[0] = globalp[0] - Object.vlist[face->verts[0]]->x;
	VP.entry[1] = globalp[1] - Object.vlist[face->verts[0]]->y;

	pre_point[0] = cur_point[0] = dot(VP, face->LX);
	pre_point[1] = cur_point[1] = dot(VP, face->LY);

	vert0[0] = Object.vlist[face->verts[0]]->x;   ////for update the global point
	vert0[1] = Object.vlist[face->verts[0]]->y;

	globalface = face_id;

	////////////////////////////////////////////////////
    for(i = 0; i < TRACESTEPS; i++)
	{
		////1. calculate the barycentric coordinates for current point
		Get2DBarycentricFacters(face_id, cur_point[0], cur_point[1], alpha);

		////2. if current point is inside current triangle
		if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1 
			&& (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1 
			&& (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
		{

			pre_point[0] = cur_point[0];
			pre_point[1] = cur_point[1];

			/*change to use other integration scheme 07/09/07*/
			//if(get_nextpt_RK4_adp(pre_point, cur_point, face_id, alpha, type))
			//if(get_nextpt_euler(pre_point, cur_point, face_id, alpha, type))
			if(get_nextpt_2ndeuler(pre_point, cur_point, face_id, alpha, type))
			//if(compute_next_pt(pre_point, cur_point, face_id, alpha, type))
			{
				////update the global point
				face = Object.flist[face_id];

				globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

				globalp[0] = vert0[0] + globalv.entry[0];
				globalp[1] = vert0[1] + globalv.entry[1];

				/*calculate the integration length 03/20/07*/
				line_length.entry[0] = cur_point[0]-pre_point[0];
				line_length.entry[1] = cur_point[1]-pre_point[1];

				/*spatial tau*/
				if(usespatialtau)
				{
					g_dt += length(line_length); 
					trace_time += length(line_length)/g_vec_mag;
				}

				/*temporal tau*/
				else
					g_dt += length(line_length)/g_vec_mag;
				
				//double len = length(line_length);
				//double vec_len = fabs(dot(line_length, g_cur_vec));
				//g_dt += len*len/vec_len;

				/*important to stop in time*/
				if(g_dt > tau)
				{
					flag = 3;
					//globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

					//globalp[0] = vert0[0] + globalv.entry[0];
					//globalp[1] = vert0[1] + globalv.entry[1];
					return face_id;
				}
			}

			else{  ////the curve reach a singularity
				flag = 3;

				////Store the record into global line segment array
                
				return face_id;
			}
		}

		////3. if the point is out of current triangle
		else{
			double t[2] = {0.};

			int PassVertornot = 0;
            
			GetNextTriangle(face_id, pre_point, cur_point, t, type, PassVertornot, alpha);

			////update the global point here
			if(PassVertornot > 0) /* Make a big change on 01/30/07 */
			{

				/* We need to move away from the vertex a little bit 01/30/07 */
				/* Here we can try a simple jitter, but make sure that it falls in the same triangle! */
			
				//we first need to know which vertex it is in the new triangle
				int vertid = pre_f->verts[PassVertornot-1];
				Face *cur_f = Object.flist[face_id];
				int vert_new = 0;
				for(int k = 0; k < 3; k++)
				{
					if(cur_f->verts[k] == vertid)
					{
						vert_new = k;
						break;
					}
				}

				alpha[vert_new]=1-0.0001;	
				alpha[(vert_new+1)%3]=0.00005;
				alpha[(vert_new+2)%3]=0.00005;


				/* Get the new cur_point */
				cur_point[0] = alpha[0]*cur_f->xy[0][0]+alpha[1]*cur_f->xy[1][0]+alpha[2]*cur_f->xy[2][0];
				cur_point[1] = alpha[0]*cur_f->xy[0][1]+alpha[1]*cur_f->xy[1][1]+alpha[2]*cur_f->xy[2][1];

				globalv = cur_point[0] * cur_f->LX + cur_point[1] * cur_f->LY;

				globalp[0] = Object.vlist[cur_f->verts[0]]->x + globalv.entry[0];
				globalp[1] = Object.vlist[cur_f->verts[0]]->y + globalv.entry[1];

				/*calculate the length of last line segment*/
				line_length.entry[0] = cur_point[0]-pre_point[0];
				line_length.entry[1] = cur_point[1]-pre_point[1];
				//
				if(usespatialtau)
				{
					g_dt += length(line_length);
					trace_time += length(line_length)/g_vec_mag;
				}

				/*temporal tau*/
				else
					g_dt += length(line_length)/g_vec_mag;
				
				/*spatial tau*/
				//double len = length(line_length);
				//double vec_len = fabs(dot(line_length, g_cur_vec));
				//g_dt += len*len/vec_len;
			}

			else{
				////transfer it to the global coordinates
				globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

				globalp[0] = vert0[0] + globalv.entry[0];
				globalp[1] = vert0[1] + globalv.entry[1];
				
				/*calculate the length of last line segment*/
				line_length.entry[0] = cur_point[0]-pre_point[0];
				line_length.entry[1] = cur_point[1]-pre_point[1];

				if(usespatialtau)
				{
					g_dt += length(line_length);
					trace_time += length(line_length)/g_vec_mag;
				}

				/*temporal tau*/
				else
					g_dt += length(line_length)/g_vec_mag;
				
				//double len = length(line_length);
				//double vec_len = fabs(dot(line_length, g_cur_vec));
				//g_dt += len*len/vec_len;
			}
			return face_id;
		}
	}

	return face_id;
}


/*
Instead trace per triangle, we now trace per vertex and save the tracing result
*/

void trace_all_Verts(double tau, int backward)
{
	int i;
	int tri_id, end_tris;
	Vertex *v;
	icVector2 stP, newP;

	//double lp[2]; 
	double alpha[3];
	Face *face;

	for(i = 0; i < Object.nverts; i++)
	{
		v = Object.vlist[i];

		if(length(v->vec_J) < 1e-10) continue;   //probably no vector value on it

		//tri_id = Object.clist[v->Corners[0]]->t;

		/* we can choose the triangle that the flow will lead the vertex go into */
		TriangleThroughVertex(i, tri_id, backward);

		if(tri_id < 0)
		{
			tri_id = Object.clist[v->Corners[0]]->t;
		}

		//face = Object.flist[tri_id];
		//for(int k = 0; k < 3; k++)
		//{
		//	if(i == face->verts[k])
		//		break;
		//}

		//alpha[k] = 1-0.00001;
		//alpha[(k+1)%3]=alpha[(k+2)%3]=0.000005;

		//lp[0] = alpha[0]*face->xy[0][0]+alpha[1]*face->xy[1][0]+alpha[2]*face->xy[2][0];
		//lp[1] = alpha[0]*face->xy[0][1]+alpha[1]*face->xy[1][1]+alpha[2]*face->xy[2][1];

		//icVector2 glv = lp[0]*face->LX + lp[1]*face->LY;
		//stP.entry[0] = Object.vlist[face->verts[0]]->x+glv.entry[0];
		//stP.entry[1] = Object.vlist[face->verts[0]]->y+glv.entry[1];

		stP.entry[0] = v->x;
		stP.entry[1] = v->y;

		/*associate with new method 03/19/07*/
		//double loc_tau = tau/(10*v->length[1]);
		//trace_Ver(tri_id, stP.entry, newP.entry, end_tris, loc_tau, backward);
		
		trace_Ver(tri_id, stP.entry, newP.entry, end_tris, tau, backward);

		v->nx = newP.entry[0];
		v->ny = newP.entry[1];

		v->RegionListID = end_tris;

		if(backward == 0)
			v->end_tri[0] = end_tris;
		else
			v->end_tri[1] = end_tris;

		/*save the information for this vertex*/
		if(backward == 0)
			v->tau[0] = trace_time;
		else
			v->tau[1] = trace_time;

	}
}

/*
Judge whether a point is inside the image hull
May have big error here!
*/
//bool is_in_region(double x, double y)
//{
//	//int i;
//	double /*theta,*/ sum_ang = 0;
//	icVector2 v1, v2;
//
//	////Calculate the sum of the angle
//	//for( i = 0; i < Num_SmoothRegionpoints; i++)
//	//{
//	//	v1.entry[0] = point[i][0] - x;
//	//	v1.entry[1] = point[i][1] - y;
//
//	//	v2.entry[0] = point[(i+1)%Num_SmoothRegionpoints][0] - x;
//	//	v2.entry[1] = point[(i+1)%Num_SmoothRegionpoints][1] - y;
//
//	//	normalize(v1);
//	//	normalize(v2);
//
//	//	////The following stuff may slow down the performance of smoothing
//
//	//	theta = GetDirectionalAngBetween2Vec(v1, v2);
//	//	sum_ang += theta;
//	//}
//
//	if( fabs(sum_ang) >= 2*M_PI - 1e-10)
//		return true;
//	else
//		return false;
//}
//
//
/*
Perform \tau tracing for each triangle
*/
void trace_Tri_Tau(int tri_id, double tau, int backward)
{
	Face *face = Object.flist[tri_id];

	Vertex *v;

	icVector2 stP;


	/*if the input tau is equal to zero*/
	if(tau < 1e-4)
	{
		/* we should use the orginal method to decide the directed edges */
		build_DirGraph_Tri_no_Tau(tri_id, backward);
		return;
	}

	int i;


	/* Trace the three vertices to get the corresponding mapping positions */
	for(i = 0; i < 3; i++)
	{
		v = Object.vlist[face->verts[i]];
		stP.entry[0] = v->x;
		stP.entry[1] = v->y;
	
		//How about we move away from the real vertex !! 02/04/07 (does not help much!)
		//double alpha[3], lp[2];
		//alpha[i] = 1-0.00001;
		//alpha[(i+1)%3]=alpha[(i+2)%3]=0.000005;

		//lp[0] = alpha[0]*face->xy[0][0]+alpha[1]*face->xy[1][0]+alpha[2]*face->xy[2][0];
		//lp[1] = alpha[0]*face->xy[0][1]+alpha[1]*face->xy[1][1]+alpha[2]*face->xy[2][1];

		//icVector2 glv = lp[0]*face->LX + lp[1]*face->LY;
		//stP.entry[0] = Object.vlist[face->verts[0]]->x+glv.entry[0];
		//stP.entry[1] = Object.vlist[face->verts[0]]->y+glv.entry[1];

		/*------------------------------------------------------*/
		//we try reuse the tracing result from each vertex
		newP[i].entry[0] = v->nx;
		newP[i].entry[1] = v->ny;
		en_tris[i] = v->RegionListID;

		//trace_Ver(tri_id, stP.entry, newP[i].entry, en_tris[i], tau, backward);
	}

	/******************************************************************************/
	/* !!!!! Here, we add an adaptive step here 03/09/07 !!!!!*/
	/* trace the center of the triangle */
	//double center[2];
	//int center_end_tri = tri_id;

	//{
	//	stP.entry[0] = (Object.vlist[face->verts[0]]->x
	//		+Object.vlist[face->verts[1]]->x + Object.vlist[face->verts[2]]->x)/3.;
	//	stP.entry[1] = (Object.vlist[face->verts[0]]->y
	//		+Object.vlist[face->verts[1]]->y + Object.vlist[face->verts[2]]->y)/3.;
	//	
	//	trace_Ver(tri_id, stP.entry, center, center_end_tri, tau, backward);
	//}

	//while (tau > 1)
	//{
	//	get_Covered_Tris();
	//	if(!IsRepeated(tri_strip, center_end_tri, ntris_in_strip))
	//	{
	//		/*decrease the tau and do tracing again*/
	//		tau /= 2;
	//		for(i = 0; i < 3; i++)
	//		{
	//			v = Object.vlist[face->verts[i]];
	//			stP.entry[0] = v->x;
	//			stP.entry[1] = v->y;
	//		
	//			trace_Ver(tri_id, stP.entry, newP[i].entry, en_tris[i], tau, backward);
	//		}
	//	}
	//	else
	//	{
	//		add_Extra_Edges_Tau(tri_id, tau);
	//		return;
	//	}
	//}


	/* Get the disk of triangles that covered by the new obtained triangle */
	//get_Tri_Strip();
	get_Covered_Tris();

	//FILE *fp = fopen("broketri.txt", "w");
	//fprintf(fp, "finish dealing triangle %d \n", tri_id);
	//fclose(fp);

	/* Add edges from the original triangle tri_id to all the triangles incident to the 
	obtained inner vertices*/
	add_Extra_Edges_Tau(tri_id, tau);

	//FILE *fp_n = fopen("num_edges.txt", "a");
	//fprintf(fp_n, "%d extra edges are added for triangle %d.\n", ntris_in_strip, tri_id);
	//fclose(fp_n);

	//if(backward == 0)
	//{
	//	if(tri_id == SelectTriangleID)
	//	{
	//		for(i = 0; i < 3; i++)
	//		{
	//			newPx[i] = newP[i].entry[0];
	//			newPy[i] = newP[i].entry[1];
	//		}
	//	}
	//}
}


/*
Judge whether there is an edge starting from tri1 to tri2
*/
bool has_Edge_From_To(int from, int to)
{
	int i;
	GraphEdge e;

	for(i = 0; i < sccnodes[from].nedges; i++)
	{
		e = sccedges[sccnodes[from].edges[i]];

		if(e.node_index1 != from)
			continue;

		if(e.node_index2 == to)
			return true;

	}
	return false;
}

bool is_boundary_tri(int tri)
{
	int i;
	for(i = 0; i < 3; i++)
	{
		if(Object.clist[tri*3+i]->o< 0)
			return true;
	}

	return false;
}

void add_Extra_Edges_Tau(int tri, double tau)
{
	int i, sccnodeid;
	//Vertex *v;
	//Corner *c;

	if(tri >= Object.nfaces)
	{
		//for nodes corresponding to boundary edges

		/* we first need to search the index of the node coording to the input tri */
		for(i = Object.nfaces; i < num_sccnodes; i++)
		{
			if(sccnodes[i].node_index == tri)
			{
				sccnodeid = i;
				break;
			}
		}
		
	}

	else{
		//for regular nodes corresponding to triangles
		//for(i = 0; i < ntris_in_strip; i++)
		//{
		//		if(has_Edge_From_To(tri, tri_strip[i]))
		//			continue;

		//		/* If there is no link between the two triangles, add a new edge */
		//		SCC_AddToEdge(tri, tri_strip[i], num_sccedges);

		//		/* add the edge to the nodes */
		//		SCC_AddEdgeToNode(tri, num_sccedges-1);
		//		SCC_AddEdgeToNode(tri_strip[i], num_sccedges-1);
		//}

		sccnodeid = tri;
	}

	for(i = 0; i < ntris_in_strip; i++)
	{
		/*if the node point to itself, we need to check whether it is a boundary triangle or not
		if it is, we don't add the edge! 02/25/07 */
		if(sccnodeid == tri_strip[i])
		{
			if(is_boundary_tri(sccnodeid))
				continue;
		}

		if(has_Edge_From_To(sccnodeid, tri_strip[i]))
			continue;

		/* If there is no link between the two triangles, add a new edge */
		SCC_AddToEdge(sccnodeid, tri_strip[i], num_sccedges);
		sccedges[num_sccedges].flow_length = tau;

		/* add the edge to the nodes */
		SCC_AddEdgeToNode(sccnodeid, num_sccedges-1);
		SCC_AddEdgeToNode(tri_strip[i], num_sccedges-1);
	}
}



void add_To_Triangle_Strip(int tri_id)
{
	if(tri_id < 0)
		return;

	if(ntris_in_strip >= curMaxTrisInStrip)
	{
		tri_strip = (int*)realloc(tri_strip, sizeof(int)*(curMaxTrisInStrip + 10));

		if(tri_strip == NULL)
		{
			MessageBox(NULL, "failed to realloc tri_strip!", "Error", MB_OK);
			exit(-1);
		}

		curMaxTrisInStrip += 10;
	}

	if(!IsRepeated(tri_strip, tri_id, ntris_in_strip))
	{
		tri_strip[ntris_in_strip] = tri_id;
		ntris_in_strip ++;
	}
}


Edge *find_Edge(int tri, int which_edge)
{
	Face *face = Object.flist[tri];
	Edge *e;
	int i;
	int v1, v2;
	
	if(which_edge == 1) //it is edge v1v2
	{
		v1 = face->verts[1];
		v2 = face->verts[2];

	}
	else if(which_edge == 2)  //it is edge v0v2
	{
		v1 = face->verts[0];
		v2 = face->verts[2];
	}
	else if(which_edge == 3) //it is edge v0v1
	{
		v1 = face->verts[0];
		v2 = face->verts[1];
	}

	else
		return NULL;

	for(i = 0; i < 3; i++)
	{
		e = face->edges[i];

		if( (e->verts[0] == v1 && e->verts[1] == v2)
			|| (e->verts[0] == v2 && e->verts[1] == v1))
			return e;
	}
}


Vertex *find_Vertex(int tri, int which_vertex)
{
	Face *face = Object.flist[tri];
	
	if(which_vertex >= 1 && which_vertex < 4) 
		return Object.vlist[face->verts[which_vertex-1]];

	else
		return NULL;
}


void trace_Edge(Edge *e, double tau)
{
	int i;
	Vertex *v;
	for(i = 0; i < 2; i++)
	{
		v = Object.vlist[e->verts[i]];

		//stP.entry[0] = v->x;
		//stP.entry[1] = v->y;
		//trace_Ver(tri_id, stP.entry, newP[i].entry, en_tris[i], tau, backward);

		/*--------------------------------------------*/
		newP[i].entry[0] = v->nx;
		newP[i].entry[1] = v->ny;
		en_tris[i] = v->RegionListID;
	}

	/* Get the disk of triangles that covered by the new obtained triangle */
	get_A_Tri_Strip(newP[0].entry, en_tris[0], newP[1].entry, en_tris[1]); 

	/* Add edges from the original triangle tri_id to all the triangles incident to the 
	obtained inner vertices*/
	add_Extra_Edges_Tau(Object.nfaces+e->index, tau); //need to change this routine

}

void trace_Boundary(int backward, double tau)
{
	int i, j;

	Face *face;
	Edge *e;

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		for(j = 0; j < 3; j++)
		{
			e = face->edges[j];

			if(e->tris[0] < 0 || e->tris[1] < 0)
			{
				trace_Edge(e, tau);
			}
		}
	}
}

extern double dmax;

/*
Update the directed graph with \tau mapping
*/
void trace_Mesh_Tau(double tau)
{
	int i;
	//double loc_tau;
	//double loc_ang, anisotropic;

	/* reset all the edges here 03/03/07 */
	init_All_Edges();

	trace_all_Verts(tau, 1);

	/* we first perform backward tracing */
	for(i = 0; i < Object.nfaces; i++)
	{
		/* use adaptive tau here */

		/* Here we need to test different ways to decide the tau!! 02/17/07 */
		//double min, max;
		//min = abs(Object.flist[i]->evalues[0]);
		//max = abs(Object.flist[i]->evalues[1]);
		//if( abs(Object.flist[i]->evalues[1]) < min)
		//{
		//	min = abs(Object.flist[i]->evalues[1]);
		//	max = abs(Object.flist[i]->evalues[0]);
		//}
		//min = min/max; /*perform a normalization based on the larger eigen value */
		
		//loc_tau = dmax*tau/min; //function 1, not good
		
		//loc_tau = (log(max)-log(min))*tau; //function 2, not good

		/*function 3: we consider both anisotropic and curl of the local flow.
		The idea is: the higher the anisotropicity and curl, the smaller tau we should use
		*/
		//loc_tau = 2*tau/(1+ .2*(log(max)-log(min))+tan(Object.flist[i]->length[2]));
		
		/*
		function 4: we find that the more curl the flow is, the more time(larger tau) we need
		*/
		//loc_tau = 4*tau*Object.flist[i]->length[2]/(1+ .5*(log(max)-log(min)));
		/*
		function 5:
		*/
		//loc_ang = Object.flist[i]->length[2];
		//anisotropic = log(max)-log(min);
		////if(loc_ang < 1e-12 && anisotropic < 1e-10)
		//if((loc_ang < 1e-10 || Object.flist[i]->length[0] < 1e-8)
		//  && anisotropic < 1e-3) /*for saddle like triangle, we don't trace*/
		//	loc_tau = 0;

		//else{
		//	if(loc_ang > 1.48)
		//		loc_ang = 1.48;

		//	loc_tau = tau*tan(loc_ang)/(1+ (log(max)-log(min)));
		//}

		//trace_Tri_Tau(i, loc_tau, 1);
		trace_Tri_Tau(i, tau , 1);
	}

	/* we need to get the image of all boundary edges 02/11/07*/

	///* Then, we inverse the partially obtained directed graph */
	ReverseEdges();

	/* reset all the edges here 03/03/07 */
	init_All_Edges();

	trace_all_Verts(tau, 0);
	/* Now, perform _forward tracing again to get the union graph */
	for(i = 0; i < Object.nfaces; i++)
	{
		/* use adaptive tau here */
		//double min;
		//min = abs(Object.flist[i]->evalues[0]);
		//if( abs(Object.flist[i]->evalues[1]) < min)
		//	min = abs(Object.flist[i]->evalues[1]);
		//loc_tau = 0.0001*tau/min;
		//double min, max;
		//min = abs(Object.flist[i]->evalues[0]);
		//max = abs(Object.flist[i]->evalues[1]);
		//if( abs(Object.flist[i]->evalues[1]) < min)
		//{
		//	min = abs(Object.flist[i]->evalues[1]);
		//	max = abs(Object.flist[i]->evalues[0]);
		//}
		//min = min/max; /*perform a normalization based on the larger eigen value */
		//loc_tau = dmax*tau/min;  // function 1, not good
		//loc_tau = (log(max)-log(min))*tau;  //function 2, not good 
		
		/*function 3: we consider both anisotropic and curl of the local flow.
		The idea is: the higher the anisotropicity and curl, the smaller tau we should use
		*/
		//loc_tau = 2*tau/(1+ .2*(log(max)-log(min))+tan(Object.flist[i]->length[2]));
		
		/*
		function 4: we find that the more curl the flow is, the more time(larger tau) we need
		*/
		//loc_tau = 4*tau*Object.flist[i]->length[2]/(1+ .5*(log(max)-log(min)));

		/*
		function 5:
		*/
		//loc_ang = Object.flist[i]->length[2];
		//anisotropic = log(max)-log(min);
		//if((loc_ang < 1e-10 || Object.flist[i]->length[0] < 1e-8)
		//  && anisotropic < 1e-3) /*for saddle like triangle, we don't trace*/
		//	loc_tau = 0;
		//else{
		//	if(loc_ang > 1.48)
		//		loc_ang = 1.48;

		//	loc_tau = tau*tan(loc_ang)/(1+ (log(max)-log(min)));
		//}

		//trace_Tri_Tau(i, loc_tau, 0);

		trace_Tri_Tau(i, tau, 0);
	}
	
	/* we need to get the image of all boundary edges 02/11/07*/

	//trace_all_centers_tris_build_edges(tau, 1);
	//trace_all_centers_tris_build_edges(tau, 0);

}


void get_A_Tri_Strip(double p1[2], int t1, double p2[2], int t2)
{

	/* We need to consider different cases before we really do tracing here 01/30/07 */
	/*
	Fortunetly, we can use the barycentric coordinates of each point to classify them into following
	cases:
	1) The two points are on the same edge
	   1.1) the two points are at the same vertex of the edge (not possible in linear field)
	   1.2) the two points are the two ending vertices of the edge (need no tracing)
	   1.3) some same/different points on the edge (need no tracing)
    2) The two points are on different edges
	   2.1) two edges are connected/sharing vertex (need no tracing)
	 * 2.2) two edges are not connected (perform tracing)
	      2.2.1) two points are both vertices
		  2.2.2) one point is vertex, the other is not
		  2.2.3) two points are both not vertices
	3) One point is on the edge, the other is inside a triangle
	   3.1) the edge is one edge of the triangle (need no tracing)
	 * 3.2) the edge is not any edge of the triangle (perform tracing)
	      3.2.1) one point is vertex, the other is inside a triangle
		  3.2.2) one point is on edge but not vertex, the other is inside a triangle
	4) Both two points are inside triangles
	   4.1) they are in the same triangle (need no tracing)
	 * 4.2) they are in different triangles (perform tracing)
	   (It seems that I've only implemented 4) !!)
	
	Note: boundary cases can be casted to group 1) and 2)

	-------------------------------------------------------------------------------------------
	To program more efficiently, we re-classify them as follows

	1) Both two points are in the same triangles (need no tracing)
	   1.1) The two points are on the same edge (add the neighboring triangle and current triangle)
	       1.1.1) at least one point is vertex (add the one-ring neighboring triangles)
	   1.2) The two points are on different edge (add the neighboring triangles and current triangle)
	       1.2.1) at least one point is vertex (add the one-ring neighboring triangles)
	   1.3) The two points are totally inside one triangle (just add the triangle)
    2) Two points are in different triangles
	   2.1) The two points are on the same edge (add the two triangles, since they share the edge)
	       2.1.1) at least one point is vertex (add the one-ring neighboring triangles)
	   * 2.2) The two points are on the different edges (perform tracing)
	       2.2.1) at least one point is vertex (add the one-ring neighboring triangles)
	   * 2.3) One point is on the edge, the other is not (perform tracing)
	       2.3.1) at least one point is vertex (add the one-ring neighboring triangles)
	       2.3.2) the point is not a vertex (add the triangles sharing the edge)
	   * 2.4) The two points are totally inside two different triangles (perform tracing)
	*/

	if(t1 < 0 || t2 < 0)
		return;

	/* According to the classification above, we first need to calculate the barycentric coordinates */
	double alpha1[3],alpha2[3], ta, tb;
	icVector2 VP;
	int class1, class2;
	int which_vertex1, which_vertex2, which_edge1, which_edge2;
	Edge *e1, *e2;
	Vertex *v1, *v2;

	int i/*, j*/;

	/* Calculate bary centric coordinates for point 1 */
	VP.entry[0] = p1[0]-Object.vlist[Object.flist[t1]->verts[0]]->x;
	VP.entry[1] = p1[1]-Object.vlist[Object.flist[t1]->verts[0]]->y;
	ta = dot(VP, Object.flist[t1]->LX);
	tb = dot(VP, Object.flist[t1]->LY);
	Get2DBarycentricFacters(t1, ta, tb, alpha1);
	/* According to its alpha values to classify this point */
	class1 = get_Type_of_Point(alpha1, which_vertex1, which_edge1);


	/* Calculate bary centric coordinates for point 2 */
	VP.entry[0] = p2[0]-Object.vlist[Object.flist[t2]->verts[0]]->x;
	VP.entry[1] = p2[1]-Object.vlist[Object.flist[t2]->verts[0]]->y;
	ta = dot(VP, Object.flist[t2]->LX);
	tb = dot(VP, Object.flist[t2]->LY);
	Get2DBarycentricFacters(t2, ta, tb, alpha2);
	/* According to its alpha values to classify this point */
	class2 = get_Type_of_Point(alpha2, which_vertex2, which_edge2);

	//if the two points fall in the same triangle
	if(t1 == t2)
	{
		//add t1 to the list
		if(ntris_in_strip >= curMaxTrisInStrip)
		{
			tri_strip = (int*)realloc(tri_strip, sizeof(int)*(curMaxTrisInStrip + 10));

			if(tri_strip == NULL)
			{
				MessageBox(NULL, "failed to realloc tri_strip!", "Error", MB_OK);
				exit(-1);
			}

			curMaxTrisInStrip += 10;
		}

		if(class1 == 0 && class2 == 0) /* both two points fall in the same triangle */
		{
			if(!IsRepeated(tri_strip, t1, ntris_in_strip))
			{
				tri_strip[ntris_in_strip] = t1;
				ntris_in_strip ++;
			}
			return;
		}

		else /* one or two points fall on the edge */
		{
			if(which_edge1 > 0 && which_edge2 > 0) //both points fall on the edges
			{
				/* get the two edges e1 and e2 */
				e1 = find_Edge(t1, which_edge1);
				e2 = find_Edge(t2, which_edge2);

				/* If they are the same edge */
				if(e1 == e2)
				{
L1:					if(which_vertex1 > 0 && which_vertex2 > 0) //both points are vertex
					{
						/* get the two vertices v1 and v2 */
						v1 = find_Vertex(t1, which_vertex1);
						v2 = find_Vertex(t2, which_vertex2);

						//if v1 = v2, v1/v2 must be a fixed point
						for(i = 0; v2 != NULL && i < v1->Num_corners; i++)
						{
							add_To_Triangle_Strip(Object.clist[v1->Corners[i]]->t);
						}
						
						/* if v1 != v2, add their one ring neighboring triangles */
						if(v1 != v2)
						{
							for(i = 0; v2 != NULL && i < v2->Num_corners; i++)
							{
								add_To_Triangle_Strip(Object.clist[v2->Corners[i]]->t);
							}
						}
					}
					else if(which_vertex1 > 0 || which_vertex2 > 0) //one point is vertex
					{
						if(which_vertex1 > 0)
						{
							/* get vertex v1 */
							v1 = find_Vertex(t1, which_vertex1);

							/* add its one ring neighboring triangles */
							for(i = 0; i < v1->Num_corners; i++)
							{
								add_To_Triangle_Strip(Object.clist[v1->Corners[i]]->t);
							}
						}
						
						else if(which_vertex2 > 0)
						{
							/* get vertex v2 */
						    v2 = find_Vertex(t2, which_vertex2);

							/* add its one ring neighboring triangles */
							for(i = 0; v2 != NULL && i < v2->Num_corners; i++)
							{
								add_To_Triangle_Strip(Object.clist[v2->Corners[i]]->t);
							}
						}
					}
					else
					{
						/* simply add the two triangles sharing the edge*/
						if(e1 != NULL)
						{
							add_To_Triangle_Strip(e1->tris[0]);
							add_To_Triangle_Strip(e1->tris[1]);
						}
					}
				}

				else{	/* If they are not the same edge */
					if(which_vertex1 > 0 && which_vertex2 > 0) //both points are vertex
					{
						/* get the two vertices v1 and v2 */
						v1 = find_Vertex(t1, which_vertex1);
						v2 = find_Vertex(t2, which_vertex2);

						//if v1 = v2, v1/v2 must be a fixed point
						for(i = 0; v1 != NULL && i < v1->Num_corners; i++)
						{
							add_To_Triangle_Strip(Object.clist[v1->Corners[i]]->t);
						}
						
						/* if v1 != v2, add their one ring neighboring triangles */
						if(v1 != v2)
						{
							for(i = 0; v2 != NULL && i < v2->Num_corners; i++)
							{
								add_To_Triangle_Strip(Object.clist[v2->Corners[i]]->t);
							}
						}
					}

					else if(which_vertex1 > 0 || which_vertex2 > 0) //one point is vertex
					{
						if(which_vertex1 > 0)
						{
							/* get vertex v1 */
							v1 = find_Vertex(t1, which_vertex1);

							/* add its one ring neighboring triangles */
							for(i = 0; v1 != NULL && i < v1->Num_corners; i++)
							{
								add_To_Triangle_Strip(Object.clist[v1->Corners[i]]->t);
							}
						}
						
						else if(which_vertex2 > 0)
						{
							/* get vertex v2 */
						    v2 = find_Vertex(t2, which_vertex2);

							/* add its one ring neighboring triangles */
							for(i = 0; v2 != NULL && i < v2->Num_corners; i++)
							{
								add_To_Triangle_Strip(Object.clist[v2->Corners[i]]->t);
							}
						}
					}
					else
					{
						/* simply add the three triangles sharing e1 and e2, respectively */
						if(e1 != NULL)
						{
							add_To_Triangle_Strip(e1->tris[0]);
							add_To_Triangle_Strip(e1->tris[1]);
						}

						if(e2 != NULL)
						{
							add_To_Triangle_Strip(e2->tris[0]);
							add_To_Triangle_Strip(e2->tris[1]);
						}
					}
				}
			}

			else if(which_edge1 > 0 || which_edge2 > 0) //one points on the edge
			{
				if( which_edge1 > 0 )
				{
					/* get the edge e1 */
					e1 = find_Edge(t1, which_edge1);

					if( which_vertex1 > 0)
					{
						/* get the vertex v1 */
						v1 = find_Vertex(t1, which_vertex1);

						/* add the one ring neighboring triangles of v1 */
						if(v1 != NULL)
						{
							/* add its one ring neighboring triangles */
							for(i = 0; v1 != NULL && i < v1->Num_corners; i++)
							{
								add_To_Triangle_Strip(Object.clist[v1->Corners[i]]->t);
							}
						}
					}
					else{
						/*add the neighboring triangle sharing e1 with t2*/
						if(e1 != NULL)
						{
							add_To_Triangle_Strip(e1->tris[0]);
							add_To_Triangle_Strip(e1->tris[1]);
						}
					}

				}

				else if (which_edge2 > 0)
				{
					/* get the edge e2 */
					e2 = find_Edge(t2, which_edge2);

					if( which_vertex2 > 0)
					{
						/* get the vertex v2 */
						v2 = find_Vertex(t2, which_vertex2);

						/* add the one ring neighboring triangles of v2 */
						if(v2 != NULL)
						{
							/* add its one ring neighboring triangles */
							for(i = 0; v2 != NULL && i < v2->Num_corners; i++)
							{
								add_To_Triangle_Strip(Object.clist[v2->Corners[i]]->t);
							}
						}
					}
					else{
						/*add the neighboring triangle sharing e2 with t1*/
						if(e2 != NULL)
						{
							add_To_Triangle_Strip(e2->tris[0]);
							add_To_Triangle_Strip(e2->tris[1]);
						}
					}

				}
			}
		}
	}

	
	else /* the two points are not in the same triangle */
	{
		if(which_edge1 > 0 && which_edge2 > 0) //both points fall on the edges
		{
			/* get the two edges e1 and e2 */
			e1 = find_Edge(t1, which_edge1);
			e2 = find_Edge(t2, which_edge2);

			/* If they are the same edge */
			if(e1 == e2)
			{
				goto L1;
			}

			else
			{
				if(which_vertex1 > 0 && which_vertex2 > 0) //both points are vertex
				{
					/* get the two vertices v1 and v2 */
					v1 = find_Vertex(t1, which_vertex1);
					v2 = find_Vertex(t2, which_vertex2);

					//if v1 = v2, v1/v2 must be a fixed point
					for(i = 0; v1 != NULL && i < v1->Num_corners; i++)
					{
						add_To_Triangle_Strip(Object.clist[v1->Corners[i]]->t);
					}
					
					/* if v1 != v2, add their one ring neighboring triangles */
					if(v1 != v2)
					{
						for(i = 0; v2 != NULL && i < v2->Num_corners; i++)
						{
							add_To_Triangle_Strip(Object.clist[v2->Corners[i]]->t);
						}
						trace_for_A_Strip(p1, t1, p2, t2);
					}
					return;
				}

				else if(which_vertex1 > 0 || which_vertex2 > 0) //one point is vertex
				{
					if(which_vertex1 > 0)
					{
						/* get vertex v1 */
						v1 = find_Vertex(t1, which_vertex1);

						/* add its one ring neighboring triangles */
						for(i = 0; v1 != NULL && i < v1->Num_corners; i++)
						{
							add_To_Triangle_Strip(Object.clist[v1->Corners[i]]->t);
						}
						trace_for_A_Strip(p2, t2, p1, t1);
					}
					
					else if(which_vertex2 > 0)
					{
						/* get vertex v2 */
						v2 = find_Vertex(t2, which_vertex2);

						/* add its one ring neighboring triangles */
						for(i = 0; v2 != NULL && i < v2->Num_corners; i++)
						{
							add_To_Triangle_Strip(Object.clist[v2->Corners[i]]->t);
						}
						trace_for_A_Strip(p1, t1, p2, t2);
					}
				}
				else
				{
					/* simply add the three triangles sharing e1 and e2, respectively */
					if(e1 != NULL)
					{
						add_To_Triangle_Strip(e1->tris[0]);
						add_To_Triangle_Strip(e1->tris[1]);
					}

					if(e2 != NULL)
					{
						add_To_Triangle_Strip(e2->tris[0]);
						add_To_Triangle_Strip(e2->tris[1]);
					}
				}
				trace_for_A_Strip(p1, t1, p2, t2);

			}
		}

		/* Perform tracing starting from p1 to p2 */
		else if(which_edge1 > 0 || which_edge2 > 0) //one points on the edge
		{
			if( which_edge1 > 0 )
			{
				/* get the edge e1 */
				e1 = find_Edge(t1, which_edge1);

				if( which_vertex1 > 0)
				{
					/* get the vertex v1 */
					v1 = find_Vertex(t1, which_vertex1);

					/* add the one ring neighboring triangles of v1 */
					if(v1 != NULL)
					{
						/* add its one ring neighboring triangles */
						for(i = 0; v1 != NULL && i < v1->Num_corners; i++)
						{
							add_To_Triangle_Strip(Object.clist[v1->Corners[i]]->t);
						}
					}
				}
				else{
					/*add the neighboring triangle sharing e1 with t2*/
					if(e1 != NULL)
					{
						add_To_Triangle_Strip(e1->tris[0]);
						add_To_Triangle_Strip(e1->tris[1]);
					}
				}

				trace_for_A_Strip(p2, t2, p1, t1);

			}

			else if (which_edge2 > 0)
			{
				/* get the edge e2 */
				e2 = find_Edge(t2, which_edge2);

				if( which_vertex2 > 0)
				{
					/* get the vertex v2 */
					v2 = find_Vertex(t2, which_vertex2);

					/* add the one ring neighboring triangles of v2 */
					if(v2 != NULL)
					{
						/* add its one ring neighboring triangles */
						for(i = 0; v2 != NULL && i < v2->Num_corners; i++)
						{
							add_To_Triangle_Strip(Object.clist[v2->Corners[i]]->t);
						}
					}
				}
				else{
					/*add the neighboring triangle sharing e2 with t1*/
					if(e2 != NULL)
					{
						add_To_Triangle_Strip(e2->tris[0]);
						add_To_Triangle_Strip(e2->tris[1]);
					}
				}
				trace_for_A_Strip(p1, t1, p2, t2);
			}
		}
		
		else
			trace_for_A_Strip(p1, t1, p2, t2);
	}

}


void get_Next_Neighbor_Tri(int &face_id, double pre[2], double cur[2],
					 int &PassVertornot, double alpha[3], icVector2 line_dir)
{
	int which_edge = -1;
	double param_t[2];

	int prev_face_id = face_id;

	Face *prev_face = Object.flist[face_id];

	Vertex *vert = NULL;

	PassVertornot = 0;
	
	////We should put pass vertex testing here before testing crossing edge
	test_Cross_Vertex(face_id, cur, pre, PassVertornot, line_dir);
	if(PassVertornot > 0)
	{
		return ;
	}

	face_id = prev_face_id;  //////
	CrossBoundary3(pre, cur, face_id, alpha, which_edge, param_t);

	if(param_t[0] == -1 && param_t[1] == -1)
	{
		face_id = prev_face_id;   ////something wrong here
		return;
	}

	////if not passing a vertex, judge which triangle it will enter later
	PassEdge(face_id, which_edge);
}


void test_Cross_Vertex(int &face_id, double cur_p[2], double pre_p[2], int &passornot, icVector2 line_dir)
{
	int i;
	double vert[2];
	//double max_alpha ;
    int newtriangleid = 0;
	int crossVert;
	Face *face = Object.flist[face_id];

	double A, B, C, pending;
	A = pre_p[1] - cur_p[1];
	B = cur_p[0] - pre_p[0];
	C = (pre_p[0]*cur_p[1] - cur_p[0]*pre_p[1]);

	passornot = 0;

	for(i = 0; i < 3; i++)
	{
		vert[0] = face->xy[i][0];
		vert[1] = face->xy[i][1];
		pending = A*vert[0] + B*vert[1] + C;
	    ////We also need to make sure that the vertex is between 'pre' and 'cur' points
		if(fabs(pending) == 0.00) ////passing the vertex
		{
			////Test whether the vertex is between 'pre' and 'cur' points
			double t;
			if(pre_p[0] != cur_p[0])
			{
				t = (vert[0] - pre_p[0])/(cur_p[0] - pre_p[0]);

			}
			else{
				t = (vert[1] - pre_p[1])/(cur_p[1] - pre_p[1]);
			}

			if(t < 0 || t > 1)
			{
				passornot = 0;
				continue;
			}

			crossVert = face->verts[i];

			////////////////////////////////////
			newtriangleid = face_id;

			//judge which triangle if might enter using the direction of the line 01/15/07
			icVector2 temp = Object.vlist[crossVert]->vec;
			//Object.vlist[crossVert]->vec = length(temp)*line_dir;
			Object.vlist[crossVert]->vec = 0.05*line_dir;
			TriangleThroughVertex(crossVert, newtriangleid, 0);
			Object.vlist[crossVert]->vec = temp;

			face_id = newtriangleid;	
			passornot = i+1;
			return;
		}
	}

	passornot = 0;
}


void get_Tri_Strip()
{
	if(tri_strip != NULL)
		free(tri_strip);

	curMaxTrisInStrip = 20;
	tri_strip = (int*)malloc(sizeof(int)*curMaxTrisInStrip);
	ntris_in_strip = 0;

	int i;

	for(i = 0; i < 3; i++)
	{
		get_A_Tri_Strip(newP[i].entry, en_tris[i], newP[(i+1)%3].entry, en_tris[(i+1)%3]);
	}
}
 


/****************************************************************/
/*03/19/07*/

/*new idea of building the $\tau$ map*/

void construct_tau_map(double tau)
{

	/*to visualize the relationship between spatial tau and temporal tau
	08/01/07*/
	init_samplepts_tautracing();

	/***1. first trace backward*/
	trace_all_Verts(tau, 1);

	/*1.2. build the edges according to the result*/
	build_edges_all_Vers(1);

	//trace_all_edges_build_edges(tau, 1);
	//trace_samples_all_edges_build_edges(tau, 1, 3); /*hard coded sampling*/
	trace_all_centers_tris_build_edges(tau, 1);
	trace_all_edges_build_di_edges_adp(tau, 1);

	/***2. perform _forward tracing*/
	trace_all_Verts(tau, 0);
	
	/*2.2. build the edges according to the result*/
	build_edges_all_Vers(0);

	//trace_all_edges_build_edges(tau, 0);
	//trace_samples_all_edges_build_edges(tau, 0, 3); /*hard coded sampling*/
	trace_all_centers_tris_build_edges(tau, 0);
	trace_all_edges_build_di_edges_adp(tau, 0);

	
	/*to visualize the relationship between spatial tau and temporal tau
	08/01/07*/
	//assign_color();
	assign_vertex_color();

	compute_density();

	assign_density_colors();

	//test_pro(tau, 0.1); /*test Konstantin's ideas*/

}

/*
This routine add one edge to the used element in the edge list
*/
void add_to_used_edge(int node_from, int node_to, int &num_edges, GraphEdge *edge)
{
	edge->cancelled = false;
	edge->node_index1 = node_from;
	edge->node_index2 = node_to;
	num_edges++;
}

/*
This routine should build an edge from node_from to node_to.
It will find an empty element in the edge list
*/
void build_one_edge(int node_from, int node_to)
{
	/*find an empty elment in the edge list
	if we can not find one, just build a new one*/
	GraphEdge *oneunused = NULL;
	int edge_index;
	bool flag = false;
	oneunused = find_unused_edge(oneunused, edge_index, flag);
	if(flag)
	{
		add_to_used_edge(node_from, node_to, num_sccedges, oneunused);
	}
	else
	{
		SCC_AddToEdge(node_from, node_to, num_sccedges);
		edge_index = num_sccedges - 1;
	}

	/*add the edge into the node's edge list*/
	SCC_AddEdgeToNode(node_from, edge_index);
	SCC_AddEdgeToNode(node_to, edge_index);
}


/*
build the edges based on the tracing result of one vertex
*/
void build_edges_Ver(int vertid, int backward)
{
	Vertex *v = Object.vlist[vertid];

		/*------------------------------------------------------*/
		//we try reuse the tracing result from each vertex
		//newP[i].entry[0] = v->nx;
		//newP[i].entry[1] = v->ny;
		//en_tris[i] = v->RegionListID;

	/* if it is backward tracing,
	build the edges from the end triangle to the one-ring of the vertex.
	if it is _forward tracing,
	build the edges from the one-ring of the vertex to the end triangle
	Note that the end triangle is saved in v->RegionListID now 03/19/07*/
	
	int i;

	for(i = 0; i < v->Num_corners; i++)
	{
		//if(Object.clist[v->Corners[i]]->t < 0 || v->RegionListID < 0)
			//continue;

		if(Object.clist[v->Corners[i]]->t < 0)
			continue;

		if(backward == 0) /*_forward tracing*/
		{
			//if(has_Edge_From_To(Object.clist[v->Corners[i]]->t, v->RegionListID))
			//	continue;
			//build_one_edge(Object.clist[v->Corners[i]]->t, v->RegionListID);
			
			if(v->end_tri[0] < 0)
				continue;

			if(has_Edge_From_To(Object.clist[v->Corners[i]]->t, v->end_tri[0]))
				continue;

			build_one_edge(Object.clist[v->Corners[i]]->t, v->end_tri[0]);
		}

		else /*backward tracing*/
		{
			//if(has_Edge_From_To(v->RegionListID, Object.clist[v->Corners[i]]->t))
			//	continue;
			//build_one_edge(v->RegionListID, Object.clist[v->Corners[i]]->t);
			
			if(v->end_tri[1] < 0)
				continue;

			if(has_Edge_From_To(v->end_tri[1], Object.clist[v->Corners[i]]->t))
				continue;
			build_one_edge(v->end_tri[1], Object.clist[v->Corners[i]]->t);
		}
	}
}

void build_edges_all_Vers(int backward)
{
	int i;

	for(i = 0; i < Object.nverts; i++)
	{
		cur_end_edgelist = 0;
		build_edges_Ver(i, backward);
	}
}


/*
trace from the middle point of each edge and see whether we can remove
the discontinuity 03/19/07
*/
void trace_all_edges_build_edges(double tau, int backward)
{
	int i, j;
	Face *face;
	Edge *e;
	Vertex *v1, *v2;
	double st[2];
	int endtri;

	init_All_Edges();

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		for(j = 0; j < 3; j++)
		{
			e = face->edges[j];
			if(e->visited == 1)
				continue;

			e->visited = 1;

			v1 = Object.vlist[e->verts[0]];
			v2 = Object.vlist[e->verts[1]];
			st[0] = (v1->x + v2->x)/2.;
			st[1] = (v1->y + v2->y)/2.;

			//trace_v
			trace_Ver(i, st, st, endtri, tau, backward);

			if(endtri < 0)
				continue;

			if(backward == 0)  /*_forward tracing*/
			{
				if(e->tris[0] >= 0)
				{
					if(!has_Edge_From_To(e->tris[0], endtri))
					{
						build_one_edge(e->tris[0], endtri);
					}
				}
				
				if(e->tris[1] >= 0)
				{
					if(!has_Edge_From_To(e->tris[1], endtri))
					{
						build_one_edge(e->tris[1], endtri);
					}
				}
			}
			else /*backward tracing*/
			{
				if(e->tris[0] >= 0)
				{
					if(!has_Edge_From_To(endtri, e->tris[0]))
					{
						build_one_edge(endtri, e->tris[0]);
					}
				}
				
				if(e->tris[1] >= 0)
				{
					if(!has_Edge_From_To(endtri, e->tris[1]))
					{
						build_one_edge(endtri, e->tris[1]);
					}
				}
			}
		}
	}
}




/*
In this routine, we trace several samples along the edge to find out the 
edges. This is really rough test. If it works, it means that
we do need to sample inside triangle
*/

void trace_samples_all_edges_build_edges(double tau, int backward, int nsamples)
{
	int i, j, k;
	Face *face;
	Edge *e;
	Vertex *v1, *v2;
	int endtri;

	icVector2 *st = new icVector2[nsamples];
	double interval = 1./(nsamples+1);
	double lamda;

	init_All_Edges();

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		for(j = 0; j < 3; j++)
		{
			e = face->edges[j];
			if(e->visited == 1)
				continue;

			e->visited = 1;

			v1 = Object.vlist[e->verts[0]];
			v2 = Object.vlist[e->verts[1]];

			/*get the set of samples*/
			for(k = 1; k <= nsamples; k++)
			{
				lamda = k*interval;
				st[k-1].entry[0] = lamda*v1->x-(1-lamda)*v2->x;
				st[k-1].entry[1] = lamda*v1->y-(1-lamda)*v2->y;

				trace_Ver(i, st[k-1].entry, st[k-1].entry, endtri, tau, backward);

				if(endtri < 0)
					continue;

				if(backward == 0)  /*_forward tracing*/
				{
					if(e->tris[0] >= 0)
					{
						if(!has_Edge_From_To(e->tris[0], endtri))
						{
							build_one_edge(e->tris[0], endtri);
						}
					}
					
					if(e->tris[1] >= 0)
					{
						if(!has_Edge_From_To(e->tris[1], endtri))
						{
							build_one_edge(e->tris[1], endtri);
						}
					}
				}
				else /*backward tracing*/
				{
					if(e->tris[0] >= 0)
					{
						if(!has_Edge_From_To(endtri, e->tris[0]))
						{
							build_one_edge(endtri, e->tris[0]);
						}
					}
					
					if(e->tris[1] >= 0)
					{
						if(!has_Edge_From_To(endtri, e->tris[1]))
						{
							build_one_edge(endtri, e->tris[1]);
						}
					}
				}
			}

		}
	}
}


/*
trace the center of the triangle to find the edge
*/
void trace_all_centers_tris_build_edges(double tau, int backward)
{
	int i;

	for(i = 0; i < Object.nfaces; i++)
	{
		/*Do we need adaptive \tau here ?*/
		trace_center_tris_build_edge(i, tau, backward);
	}
}


void trace_center_tris_build_edge(int tri, double tau, int backward)
{
	double center[2] = {0.};

	Face *face = Object.flist[tri];

	int i, endtri;

	for(i = 0; i < 3; i++)
	{
		center[0] += Object.vlist[face->verts[i]]->x;
		center[1] += Object.vlist[face->verts[i]]->y;
	}

	center[0] /= 3.;
	center[1] /= 3.;

	/*start tracing here*/
	trace_Ver(tri, center, center, endtri, tau, backward);

	if(endtri < 0)
		return;

	if(backward == 0)  /*_forward tracing*/
	{
		if(!has_Edge_From_To(tri, endtri))
		{
			build_one_edge(tri, endtri);
		}
	}
		
	else /*backward tracing*/
	{
		if(!has_Edge_From_To(endtri, tri))
		{
			build_one_edge(endtri, tri);
		}
	}
}



void trace_samples_in_tris_build_edges(double tau, int backward, int nsamples)
{
	int i;

	for(i = 0; i < Object.nfaces; i++)
	{
		/*Do we need adaptive \tau here ?*/
		//trace_samples_tri_build_edges(i, tau, backward);
	}
}


void trace_samples_tri_build_edges(int tri, double tau, int backward, int nsamples)
{
	double center[2] = {0.};

	Face *face = Object.flist[tri];

	int i, endtri;

	for(i = 0; i < 3; i++)
	{
		center[0] += Object.vlist[face->verts[i]]->x;
		center[1] += Object.vlist[face->verts[i]]->y;
	}

	center[0] /= 3.;
	center[1] /= 3.;

	/*start tracing here*/
	trace_Ver(tri, center, center, endtri, tau, backward);

	if(endtri < 0)
		return;

	if(backward == 0)  /*_forward tracing*/
	{
		if(!has_Edge_From_To(tri, endtri))
		{
			build_one_edge(tri, endtri);
		}
	}
		
	else /*backward tracing*/
	{
		if(!has_Edge_From_To(endtri, tri))
		{
			build_one_edge(endtri, tri);
		}
	}
}

/*
If the two triangle share the same edge, they are close neighbors
03/21/07
*/
bool are_close_neighbors(int t1, int t2)
{
	Face *f1 = Object.flist[t1];
	Face *f2 = Object.flist[t2];

	int i, j;

	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			if(f1->edges[i] == f2->edges[j])
				return true;
		}
	}
	return false;
}


/*how about if they share a vertex
If the two triangles share a vertex, they are pseudo-close neighbors
*/
bool are_pseudo_close_neighbors(int t1, int t2, int &verid)
{
	Face *f1 = Object.flist[t1];
	Face *f2 = Object.flist[t2];

	int i, j;
	verid = -1;

	/*if the two triangles share an edge*/
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			if(f1->verts[i] == f2->verts[j])
			{
				verid = f1->verts[i];
				return true;
			}
		}
	}

	return false;
}

/*
An recursive algorithm for implementation of adaptively sampling 
along an edge 03/21/07
*/
void trace_recursive_an_edge(double v1[2], double v2[2], int &t1, int tri, int neighbor_tri,
							 double tau, int backward, int &level)
{
	double stack_v[2];  /*this could only have one element in one level*/
	int top = 0;
	int t2;
	double v2_end[2];

	//int share_ver;

	level++;

	/*boundary!!*/
	if(t1 < 0)
		return;
		
	trace_Ver(tri, v2, v2_end, t2, tau, backward);

	/*boundary!!*/
	if(t2 < 0)
		return;

	//icVector3 dis;
	//dis.entry[0] = v1[0]-v2[0];
	//dis.entry[1] = v1[1]-v2[1];
	//dis.entry[2] = v1[2]-v2[2];

	//if(length(dis) < 1e-8)
	//{
	//	v1[0] = v2[0];
	//	v1[1] = v2[1];
	//	v1[2] = v2[2];
	//	t1 = t2;
	//	/*need to build edges here*/
	//	if(backward == 0) /*_forward tracing*/
	//	{
	//		if(!has_Edge_From_To(tri, t1))
	//			build_one_edge(tri, t1);
	//		if(neighbor_tri>=0&&!has_Edge_From_To(neighbor_tri, t1))
	//			build_one_edge(neighbor_tri, t1);
	//	}
	//	else /*backward tracing*/
	//	{
	//		if(!has_Edge_From_To(t1, tri))
	//			build_one_edge(t1, tri);
	//		if(neighbor_tri>=0&&!has_Edge_From_To(t1, neighbor_tri))
	//			build_one_edge(t1, neighbor_tri);
	//	}
	//	return;
	//}


	/*if the two points are too close, stop!*/
	if(level > 5)
	{
		/*trace v2, get t2 */
		v1[0] = v2[0];
		v1[1] = v2[1];
		t1 = t2;
		/*need to build edges here*/
		if(backward == 0) /*_forward tracing*/
		{
			if(!has_Edge_From_To(tri, t1))
				build_one_edge(tri, t1);
			if(neighbor_tri>=0&&!has_Edge_From_To(neighbor_tri, t1))
				build_one_edge(neighbor_tri, t1);
		}
		else /*backward tracing*/
		{
			if(!has_Edge_From_To(t1, tri))
				build_one_edge(t1, tri);
			if(neighbor_tri>=0&&!has_Edge_From_To(t1, neighbor_tri))
				build_one_edge(t1, neighbor_tri);
		}
		return;
	}

	if(t1 == t2 || are_close_neighbors(t1, t2))
	{
		v1[0] = v2[0];
		v1[1] = v2[1];
		t1 = t2;
		
		/*need to build edges here*/
		if(backward == 0) /*_forward tracing*/
		{
			if(!has_Edge_From_To(tri, t1))
				build_one_edge(tri, t1);
			if(neighbor_tri>=0&&!has_Edge_From_To(neighbor_tri, t1))
				build_one_edge(neighbor_tri, t1);
		}
		else /*backward tracing*/
		{
			if(!has_Edge_From_To(t1, tri))
				build_one_edge(t1, tri);
			if(neighbor_tri>=0&&!has_Edge_From_To(t1, neighbor_tri))
				build_one_edge(t1, neighbor_tri);
		}
		return;
	}


	/*consider further approximation of the closure of the image 07/24/07*/
	//else if(are_pseudo_close_neighbors(t1, t2, share_ver))
	//{
	//	v1[0] = v2[0];
	//	v1[1] = v2[1];
	//	t1 = t2;

	//	/*add the one-ring neighborhood of the share_ver to the edge list*/
	//	Vertex* v = Object.vlist[share_ver];
	//	Corner *c;
	//	for(int i=0; i<v->Num_corners; i++)
	//	{
	//		c = Object.clist[v->Corners[i]];
	//		/*need to build edges here*/
	//		if(backward == 0) /*_forward tracing*/
	//		{
	//			if(!has_Edge_From_To(tri, c->t))
	//				build_one_edge(tri, c->t);
	//			if(neighbor_tri>=0&&!has_Edge_From_To(neighbor_tri, c->t))
	//				build_one_edge(neighbor_tri, c->t);
	//		}
	//		else /*backward tracing*/
	//		{
	//			if(!has_Edge_From_To(c->t, tri))
	//				build_one_edge(c->t, tri);
	//			if(neighbor_tri>=0&&!has_Edge_From_To(c->t, neighbor_tri))
	//				build_one_edge(c->t, neighbor_tri);
	//		}
	//	}
	//	return;
	//}

	else
	{
		/*push v2 into the stack and indicate that the stack is not empty any more*/
		stack_v[0] = v2[0];
		stack_v[1] = v2[1];
		top++;

		/*use the middle point of this line segment as current v2*/
		v2[0] = (v1[0]+v2[0])/2.;
		v2[1] = (v1[1]+v2[1])/2.;		

		/*recursively call the routine itself*/
		trace_recursive_an_edge(v1, v2, t1, tri, neighbor_tri, tau, backward, level);
		level --;
	}

	/*if the stack is not empty*/
	if(top > 0)
	{
		/*pop up the point*/
		v2[0] = stack_v[0];
		v2[1] = stack_v[1];

		top--;
		/*recursively call the routine itself*/
		trace_recursive_an_edge(v1, v2, t1, tri, neighbor_tri, tau, backward, level);
		level--;
	}
}




/*
An recursive algorithm for implementation of adaptively sampling 
along an edge 03/21/07
*/
void nonrecursive_trace_an_edge(double v1[2], double v2[2], int &t1, int tri, int neighbor_tri,
							 double tau, int backward, int &level, int edge_id)
{
	icVector2 *stack_v = (icVector2*)malloc(sizeof(icVector2)*10);  /*this could only have one element in one level*/
	int curMaxNumStackElems = 10;
	int top = 0;
	int t2;
	double v2_end[2];

	//int share_ver;

	level = 0;

	bool usestack = false;

	do{
		/*boundary!!*/
		while(1)
		{
			level++;

			if(t1 < 0)
				break;
				
			trace_Ver(tri, v2, v2_end, t2, tau, backward);

			/*boundary!!*/
			if(t2 < 0)
				break;

			/*if the two points are too close, stop!*/
			if(level > 5)
			{
				/*trace v2, get t2 */
				v1[0] = v2[0];
				v1[1] = v2[1];
				t1 = t2;
				/*need to build edges here*/
				if(backward == 0) /*_forward tracing*/
				{
					if(!has_Edge_From_To(tri, t1))
						build_one_edge(tri, t1);
					if(neighbor_tri>=0&&!has_Edge_From_To(neighbor_tri, t1))
						build_one_edge(neighbor_tri, t1);
				}
				else /*backward tracing*/
				{
					if(!has_Edge_From_To(t1, tri))
						build_one_edge(t1, tri);
					if(neighbor_tri>=0&&!has_Edge_From_To(t1, neighbor_tri))
						build_one_edge(t1, neighbor_tri);
				}
				break;
			}

			if(t1 == t2 || are_close_neighbors(t1, t2))
			{
				v1[0] = v2[0];
				v1[1] = v2[1];
				t1 = t2;
				
				/*need to build edges here*/
				if(backward == 0) /*_forward tracing*/
				{
					if(!has_Edge_From_To(tri, t1))
						build_one_edge(tri, t1);
					if(neighbor_tri>=0&&!has_Edge_From_To(neighbor_tri, t1))
						build_one_edge(neighbor_tri, t1);
				}
				else /*backward tracing*/
				{
					if(!has_Edge_From_To(t1, tri))
						build_one_edge(t1, tri);
					if(neighbor_tri>=0&&!has_Edge_From_To(t1, neighbor_tri))
						build_one_edge(t1, neighbor_tri);
				}
				break;
			}

			else
			{
				/*push v2 into the stack and indicate that the stack is not empty any more*/
				stack_v[top].entry[0] = v2[0];
				stack_v[top].entry[1] = v2[1];

				top++;

				/*extend the list*/
				if(top >= curMaxNumStackElems)
				{
					icVector2 *temp = stack_v;
					stack_v = (icVector2*)malloc(sizeof(icVector2)*(curMaxNumStackElems+10));
					
					if(stack_v == NULL)
					{
						fprintf(stderr, "not enough memory for stack in the adaptive edge sampling\n");
						exit(-1);
					}

					for(int j=0; j<curMaxNumStackElems-10; j++)
						stack_v[j]=temp[j];

					free(temp);
				}

				/*use the middle point of this line segment as current v2*/
				v2[0] = (v1[0]+v2[0])/2.;
				v2[1] = (v1[1]+v2[1])/2.;	

				usestack = true;

			}

		}

		/*to visualize the relationship between spatial tau and temporal tau
			08/01/07*/
			/*save the results*/
			{
				TraceSamplePt *s = (TraceSamplePt*)malloc(sizeof(TraceSamplePt));
				if(!usestack)
				{
					s->pos[0] = v2[0];
					s->pos[1] = v2[1];
				}
				else
				{
					s->pos[0] = stack_v[top].entry[0];
					s->pos[1] = stack_v[top].entry[1];
				}
				s->time = trace_time;
				s->edge = edge_id;
				if(backward==0) _forward->addNew(s);
				else backward_spts->addNew(s);
			}


		level--;
		/*if the stack is not empty*/
		if(top > 0)
		{
			/*pop up the point*/
			v2[0] = stack_v[top-1].entry[0];
			v2[1] = stack_v[top-1].entry[1];

			top--;
		}

		else
		{
			free(stack_v);
			return;
		}
	}while (1);
}


/*
Here if the tau's at two vertices are different, we need to do the interpolation to get the tau
We can do this is because tau is a continuous function 03/22/07
*/
void trace_recursive_an_edge_interpolate(double v1[2], double v2[2], int &t1, int tri, int neighbor_tri,
							 double &tau1, double &tau2, int backward, int &level)
{
	double stack_v[3];  /*this could only have one element in one level*/
	int top = 0;
	int t2;
	double v2_end[2];

	level++;

	/*boundary!!*/
	if(t1 < 0)
		return;
		
	trace_Ver(tri, v2, v2_end, t2, tau2, backward);

	/*boundary!!*/
	if(t2 < 0)
		return;

	icVector3 dis;
	dis.entry[0] = v1[0]-v2[0];
	dis.entry[1] = v1[1]-v2[1];
	dis.entry[2] = v1[2]-v2[2];

	if(length(dis) < 1e-8)
	{
		v1[0] = v2[0];
		v1[1] = v2[1];
		v1[2] = v2[2];
		t1 = t2;
		/*need to build edges here*/
		if(backward == 0) /*_forward tracing*/
		{
			if(!has_Edge_From_To(tri, t1))
				build_one_edge(tri, t1);
			if(neighbor_tri>=0&&!has_Edge_From_To(neighbor_tri, t1))
				build_one_edge(neighbor_tri, t1);
		}
		else /*backward tracing*/
		{
			if(!has_Edge_From_To(t1, tri))
				build_one_edge(t1, tri);
			if(neighbor_tri>=0&&!has_Edge_From_To(t1, neighbor_tri))
				build_one_edge(t1, neighbor_tri);
		}
		return;
	}

	/*if the two points are too close, stop!*/
	if(level > 6)
	{
		/*trace v2, get t2 */
		v1[0] = v2[0];
		v1[1] = v2[1];
		t1 = t2;
		tau1 = tau2;
		/*need to build edges here*/
		if(backward == 0) /*_forward tracing*/
		{
			if(!has_Edge_From_To(tri, t1))
				build_one_edge(tri, t1);
			if(neighbor_tri>=0&&!has_Edge_From_To(neighbor_tri, t1))
				build_one_edge(neighbor_tri, t1);
		}
		else /*backward tracing*/
		{
			if(!has_Edge_From_To(t1, tri))
				build_one_edge(t1, tri);
			if(neighbor_tri>=0&&!has_Edge_From_To(t1, neighbor_tri))
				build_one_edge(t1, neighbor_tri);
		}
		return;
	}

	if(t1 == t2 || are_close_neighbors(t1, t2))
	{
		v1[0] = v2[0];
		v1[1] = v2[1];
		t1 = t2;
		tau1 = tau2;
		
		/*need to build edges here*/
		if(backward == 0) /*_forward tracing*/
		{
			if(!has_Edge_From_To(tri, t1))
				build_one_edge(tri, t1);
			if(neighbor_tri>=0&&!has_Edge_From_To(neighbor_tri, t1))
				build_one_edge(neighbor_tri, t1);
		}
		else /*backward tracing*/
		{
			if(!has_Edge_From_To(t1, tri))
				build_one_edge(t1, tri);
			if(neighbor_tri>=0&&!has_Edge_From_To(t1, neighbor_tri))
				build_one_edge(t1, neighbor_tri);
		}
		return;
	}

	else
	{
		/*push v2 into the stack and indicate that the stack is not empty any more*/
		stack_v[0] = v2[0];
		stack_v[1] = v2[1];
		stack_v[2] = tau2;   /*save the tau of the point*/
		top++;

		/*use the middle point of this line segment as current v2*/
		v2[0] = (v1[0]+v2[0])/2.;
		v2[1] = (v1[1]+v2[1])/2.;		

		double tau = (tau1+tau2)/2;

		/*recursively call the routine itself*/
		trace_recursive_an_edge_interpolate(v1, v2, t1, tri, neighbor_tri, tau1, tau, backward, level);
		level --;
	}

	/*if the stack is not empty*/
	if(top > 0)
	{
		/*pop up the point*/
		v2[0] = stack_v[0];
		v2[1] = stack_v[1];
		tau2 = stack_v[2];

		top--;
		/*recursively call the routine itself*/
		trace_recursive_an_edge_interpolate(v1, v2, t1, tri, neighbor_tri, tau1, tau2, backward, level);
		level--;
	}
}

/*
trace from all the edges using adaptive framework
03/21/07
*/
void trace_all_edges_build_di_edges_adp(double tau, int backward)
{
	int i, j;
	Face *face;
	Edge *e;
	Vertex *v1, *v2;
	double st1[2], st2[2];
	//int endtri;

	init_All_Edges();

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		for(j = 0; j < 3; j++)
		{
			e = face->edges[j];
			if(e->visited == 1)
				continue;

			e->visited = 1;

			v1 = Object.vlist[e->verts[0]];
			v2 = Object.vlist[e->verts[1]];

			/*currently, we don't deal with the boundary edges!!!!*/
			if(v1->RegionListID < 0 || v2->RegionListID < 0)
				continue;

			/*if the images of the two vertices are already continuous, need not
			consider this edge any more*/
			if(v1->RegionListID == v2->RegionListID
				|| are_close_neighbors(v1->RegionListID, v2->RegionListID))
				continue;
			
			st1[0] = v1->x;
			st1[1] = v1->y;
			st2[0] = v2->x;
			st2[1] = v2->y;

			int neighbor_tri = e->tris[0];
			if(neighbor_tri == i)
				neighbor_tri = e->tris[1];


			trace_an_edge_build_di_edges_adp(st1, st2, v1->RegionListID, v2->RegionListID, i, neighbor_tri,
				tau, backward, e->index);

		}
	}
}

/*
trace from one edge using adaptive framework 03/21/07
*/
void trace_an_edge_build_di_edges_adp(double st1[2], double st2[2], int t1, int t2, int tri,
									  int neighbor_tri, double tau, int backward, int edge_id)
{
	/*do one more recursive here*/
	int level = 1;
	double stack_st[2] = {st2[0], st2[1]}; /*save vertex v2*/

	st2[0] = (st1[0]+st2[0])/2;
	st2[1] = (st1[1]+st2[1])/2;

	double middle_p[2] = {st2[0], st2[1]};/*save the middle point*/

	/*call the recursive adaptive edge sampling for the first half of the edge*/
	//trace_recursive_an_edge(st1, st2, t1, tri, neighbor_tri, tau, backward, level, edge_id);
	nonrecursive_trace_an_edge(st1, st2, t1, tri, neighbor_tri, tau, backward, level, edge_id);

	level = 1;
	/*call the recursive adaptive edge sampling for the second half of the edge*/
	//trace_recursive_an_edge(stack_st, middle_p, t2, tri, neighbor_tri, tau, backward, level, edge_id);
	nonrecursive_trace_an_edge(stack_st, middle_p, t2, tri, neighbor_tri, tau, backward, level, edge_id);
}


void trace_an_edge_build_di_edges_adp_interpolate(double st1[2], double st2[2], int t1, int t2, int tri,
									  int neighbor_tri, double tau1, double tau2, int backward)
{
	/*do one more recursive here*/
	int level = 1;
	double stack_st[2] = {st2[0], st2[1]}; /*save vertex v2*/
    double v2_tau = tau2;                  /*save the tau for vertex v2*/

	st2[0] = (st1[0]+st2[0])/2;
	st2[1] = (st1[1]+st2[1])/2;

	double middle_p[2] = {st2[0], st2[1]}; /*save the middle point*/
	double middle_tau = (tau1+tau2)/2;     /*save the tau for the middle point*/

	/*call the recursive adaptive edge sampling for the first half of the edge*/
	trace_recursive_an_edge_interpolate(st1, st2, t1, tri, neighbor_tri, tau1, middle_tau, backward, level);

	level = 1;
	/*call the recursive adaptive edge sampling for the second half of the edge*/
	trace_recursive_an_edge_interpolate(stack_st, middle_p, t2, tri, neighbor_tri, v2_tau, middle_tau,
		backward, level);
}






/*
*/
int get_Type_of_Point(double alpha[3], int &which_vertex, int &which_edge)
{
	int type = 0;
	which_vertex = -1;
	which_edge = -1;

	if(alpha[0] == 1 || alpha[1] == 1 || alpha[2] == 1)
	{
		/* on vertex */
		if(alpha[0] == 1)
			which_vertex = 1;
		else if(alpha[1] == 1)
			which_vertex = 2;
		else
			which_vertex = 3;

		type = 1;
	}

	if(alpha[0] == 0 || alpha[1] == 0 || alpha[2] == 0)
	{
		/* on edge */
		if(alpha[0] == 0)
			which_edge = 1;
		else if(alpha[1] == 0)
			which_edge = 2;
		else
			which_edge = 3;

		type = 1;

	}

	/* inside the triangle */
	return type;

}


void trace_for_A_Strip(double p1[2], int t1, double p2[2], int t2)
{
	int cur_fid = t1;
	Face *cur_face = Object.flist[t1];

	double l_p1[2], l_p2[2];
	double alpha[3], g_p1[2];

	//int i;

	icVector2 line_dir;  //this is a global direction vector ! 01/22/07
	line_dir.entry[0] = p2[0] - p1[0];
	line_dir.entry[1] = p2[1] - p1[1];
    normalize(line_dir);

	//transfer p2 to local coordinates
	icVector2 VP;

	g_p1[0] = p1[0];
	g_p1[1] = p1[1];

	int count = 0;  // avoid dead loop 

	add_To_Triangle_Strip(t1);
	add_To_Triangle_Strip(t2);

	while(cur_fid != t2 && count < Object.nverts/10 )
	{

		if(cur_fid < 0)
			return;

		cur_face = Object.flist[cur_fid];

		VP.entry[0] = g_p1[0] - Object.vlist[cur_face->verts[0]]->x;
		VP.entry[1] = g_p1[1] - Object.vlist[cur_face->verts[0]]->y;

		l_p1[0] = dot(VP, cur_face->LX);
		l_p1[1] = dot(VP, cur_face->LY);

		VP.entry[0] = p2[0] - Object.vlist[cur_face->verts[0]]->x;
		VP.entry[1] = p2[1] - Object.vlist[cur_face->verts[0]]->y;

		l_p2[0] = dot(VP, cur_face->LX);
		l_p2[1] = dot(VP, cur_face->LY);

		//judge whether it is in current triangle
		Get2DBarycentricFacters(cur_fid, l_p2[0], l_p2[1], alpha);

		//if the point is inside current triangle
		if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1 
			&& (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1 
			&& (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
		{
			return;
		}

		else
		{
			//find the new triangle, add it to the tri_strip list, set it as current triangle
			int PassVertornot = 0;
			Face *pre_f = Object.flist[cur_fid];

			//double t[2];
			get_Next_Neighbor_Tri(cur_fid, l_p1, l_p2, PassVertornot, alpha, line_dir);

			if(PassVertornot > 0)
			{
				////we should not directly use the vertex as next point!!
				////we may move a little bit along the VF direction, but make sure it is still inside
				////the new triangle

				Vertex *PassedVert = Object.vlist[pre_f->verts[PassVertornot-1]];

				/* Add the triangles in the one-ring neighborhood of PassedVert */
				for(int k = 0; k < PassedVert->Num_corners; k++)
				{
					add_To_Triangle_Strip(Object.clist[PassedVert->Corners[k]]->t);
				}
				
				/*!!!!! Important improvement here 01/23/07 !!!!!*/
				/* We need to move off the vertex a little bit */

				//g_p1[0] = PassedVert->x + SEPARATRIXSTEP*line_dir.entry[0]/100.;
				//g_p1[1] = PassedVert->y + SEPARATRIXSTEP*line_dir.entry[1]/100.;

				cur_face = Object.flist[cur_fid];

				/* We need to make sure that the new point is still inside the triangle */
				{
					//double alpha1[3] = {0.};
           
					//VP.entry[0] = g_p1[0] - Object.vlist[cur_face->verts[0]]->x;
					//VP.entry[1] = g_p1[1] - Object.vlist[cur_face->verts[0]]->y;

					//l_p1[0] = dot(VP, cur_face->LX);
					//l_p1[1] = dot(VP, cur_face->LY);

					/* Method 1*/
					//int count2 = 2;

					//while(!FallIntheTriangle(cur_fid, g_p1) && count2 < 20)
					//{
					//	double scale = SEPARATRIXSTEP/pow(2, count2);
					//	g_p1[0] = PassedVert->x + scale*line_dir.entry[0]/20.;
					//	g_p1[1] = PassedVert->y + scale*line_dir.entry[1]/20.;
					//	count2++;
					//}

					//if(count2 > 20)
					//	return;

					/* Method 2*/
					int vertid = pre_f->verts[PassVertornot-1];
					Face *cur_f = Object.flist[cur_fid];
					int vert_new = 0;
					for(int k = 0; k < 3; k++)
					{
						if(cur_f->verts[k] == vertid)
						{
							vert_new = k;
							break;
						}
					}

					alpha[vert_new]=1-0.0001;	
					alpha[(vert_new+1)%3]=0.00005;
					alpha[(vert_new+2)%3]=0.00005;


					/* Get the new cur_point */
					l_p1[0] = alpha[0]*cur_f->xy[0][0]+alpha[1]*cur_f->xy[1][0]+alpha[2]*cur_f->xy[2][0];
					l_p1[1] = alpha[0]*cur_f->xy[0][1]+alpha[1]*cur_f->xy[1][1]+alpha[2]*cur_f->xy[2][1];

					icVector2 globalv = l_p1[0] * cur_f->LX + l_p1[1] * cur_f->LY;

					g_p1[0] = Object.vlist[cur_f->verts[0]]->x + globalv.entry[0];
					g_p1[1] = Object.vlist[cur_f->verts[0]]->y + globalv.entry[1];
				}
			}

			else{
				//set the intersection as current l_p1
				l_p1[0] = l_p2[0];
				l_p1[1] = l_p2[1];

				icVector2 globalv = l_p1[0] * pre_f->LX + l_p1[1] * pre_f->LY;
				g_p1[0] = Object.vlist[pre_f->verts[0]]->x + globalv.entry[0];
				g_p1[1] = Object.vlist[pre_f->verts[0]]->y + globalv.entry[1];
			}
			
			//add the new triangle into the list
			//if(ntris_in_strip >= curMaxTrisInStrip)
			//{
			//	tri_strip = (int*)realloc(tri_strip, sizeof(int)*(curMaxTrisInStrip + 10));

			//	if(tri_strip == NULL)
			//	{
			//		MessageBox(NULL, "failed to realloc tri_strip!", "Error", MB_OK);
			//		exit(-1);
			//	}

			//	curMaxTrisInStrip += 10;
			//}

			//if(!IsRepeated(tri_strip, cur_fid, ntris_in_strip))
			//{
			//	tri_strip[ntris_in_strip] = cur_fid;
			//	ntris_in_strip ++;
			//}
			add_To_Triangle_Strip(cur_fid);
		}
		count++;
	}

	if(count >= Object.nverts/10)
	{
		int test = 0;
	}

}



extern int CalEulerValue(int *trianglelist, int num);
extern void GetBoundary();


extern RegionBoundary InnerBoundary;
extern RegionBoundary OuterBoundary;
extern int MaxNumEdgesOnBoundary;
extern Edge **regionedge;                   ////mesh edges of user selected region
extern Edge **Cycle_edge;                   ////mesh edges of the boundary of the cell cycle
extern int MaxEdgeInCycle;
extern int num_cycleedges;                  ////number of edges consist of the cell cycle 

extern void AddToEdgeList(Edge *);
extern void GetTwoBoundaries(Edge **, int, int);	
extern void BuildBoundaryEdgeList(int *, int);

void get_Covered_Tris()
{
	//1. get the strip
	get_Tri_Strip();

	if(ntris_in_strip == 0) return;

	//2. if the strip is a disk, return the disk
	if(CalEulerValue(tri_strip, ntris_in_strip)==1)
	{
		return;
	}

	//3. else, find the outer boundary, and one inner vertex
	int i;
	for(i = 0; i < ntris_in_strip; i++)
	{
		Object.flist[tri_strip[i]]->inDesignCellCycle = 1;

	}

	get_Boundary_2(tri_strip, ntris_in_strip);

	//reset the inner and outer boundaries
	Edge *e;
	for(i = 0; i < InnerBoundary.num; i++)
	{
		e = InnerBoundary.edgelist[i];
		e->OnBoundary = 0;
	}
	//copy the edges to the region boundary edge list
	Num_edges = 0;
	for(i = 0; i < OuterBoundary.num; i++)
	{
		AddToEdgeList(OuterBoundary.edgelist[i]);
	}

	//copy the obtained new vertices (global points) of the new triangle to the Point list
	for(i = 0; i < 3; i++)
	{
		point[i].x = newP[i].entry[0]; //it seems that we should set up the boundary vertices instead of 
		point[i].y = newP[i].entry[1];
	}
	Num_SmoothRegionpoints = 3;

	//4. find out other inner vertices if any
	Num_verts = 0;
	FindOutInnerVerts();

	//5. add the one ring of each inner vertex to the triangle strip
	Vertex *v;
	Corner *c;
	int j;
	for(i = 0; i < Num_verts; i++)
	{
		v = regionverts[i];

		for(j = 0; j < v->Num_corners; j++)
		{
			c = Object.clist[v->Corners[j]];

			add_To_Triangle_Strip(c->t);

		}
	}

	reset_Triangle_Flags();
	reset_Edge_Flags();
	Num_verts = 0;
}

void reset_Triangle_Flags(/*int *tri_strip, int ntris_in_strip*/)
{
	int i;
	Face *face;

	/* We only need to reset the triangles in the previous strip 01/22/07*/
	for(i = 0; i < ntris_in_strip; i++)
	{
		face = Object.flist[tri_strip[i]];
		face->inDesignCellCycle = 0;
	}

	/* We also need to reset the flags of the vertices in the region */
	int j;
	Vertex *v;
	for(i = 0 ; i < ntris_in_strip; i++)
	{
		face = Object.flist[tri_strip[i]];
		for(j = 0; j < 3; j++)
		{
			v = Object.vlist[face->verts[j]];
			v->InRegion = 0;
			v->OnBoundary = 0;
			v->RegionListID = -1;
		}
	}
}


void reset_Edge_Flags(/*int *tri_strip, int ntris_in_strip*/)
{
	int i, j;
	Face *face;
	Edge *cur_edge;

	/* We only need to reset the triangles in the previous strip 01/22/07*/
	for(i = 0; i < ntris_in_strip; i++)
	{
		face = Object.flist[tri_strip[i]];
		//face->inDesignCellCycle = 0;

		for(j = 0; j < 3; j++)
		{
			cur_edge = face->edges[j];
			cur_edge->OnBoundary = 0;
			cur_edge->BoundaryVisited = 0;
		}
	}

}


void init_All_Edges()
{
	int i, j;
	Face *face;
	Edge *cur_edge;
	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		for(j = 0; j < 3; j++)
		{
			cur_edge = face->edges[j];
			cur_edge->OnBoundary = 0;
			cur_edge->BoundaryVisited = 0;
			cur_edge->visited = 0;
		}
	}
}

/*
This routine can not gaurantee to return the correct pointer!!
*/
int *add_To_Int_Array(int *array, int &num, int a, int &Max, int ex_step)
{
	if(num >= Max)
	{
		array = (int *)realloc(array, sizeof(int)*(Max+ex_step));
		if(array == NULL)
			exit(-1);
		/* We need to somehow write the execution error into some log */
		//err<<

		Max += ex_step;
	}
	
	if(!IsRepeated(array, a, num))
	{
		array[num] = a;
		num++;
	}

	return array;
}

/*
find out the non-manifold vertices and fix them by adding more triangles
*/
void fill_Gap_In_Ring(Edge **edgelist, int n_edges, int *triangles, int ntris)
{
	/* first, find all the non-manifold vertices */
	int i;
	Vertex *v;
	Edge *e;

	for(i = 0; i < n_edges; i++)
	{
		e = edgelist[i];

		Object.vlist[e->verts[0]]->tau_visited = 0;
		Object.vlist[e->verts[1]]->tau_visited = 0;
	}

    /* Count the number of visiting of each vertex */
	for(i = 0; i < n_edges; i++)
	{
		e = edgelist[i];

		Object.vlist[e->verts[0]]->tau_visited ++;
		Object.vlist[e->verts[1]]->tau_visited ++;
	}

	/* Fill up the gap by adding more neighboring triangles */
	for(i = 0; i < n_edges; i++)
	{
		e = edgelist[i];

		v = Object.vlist[e->verts[0]];

		if(v->tau_visited > 2) //non-manifold vertex /*could be problematic since we change the type to be bool*/
		{
			for(int j = 0; j < v->Num_corners; j++)
			{
				/* fill up the gap by adding triangles in the STAR structure */
				if(!IsRepeated(triangles, Object.clist[v->Corners[j]]->t, ntris))
				{
					/*add to the strip*/
					if(ntris_in_strip > curMaxTrisInStrip)
					{
						tri_strip = (int*)realloc(tri_strip, sizeof(int)*(curMaxTrisInStrip+20));
						if(tri_strip == NULL)
						{
							MessageBox(NULL, "fail to allocate mem for tri_strip!", "Error", MB_OK);
							exit(-1);
						}
						tri_strip[ntris_in_strip] = Object.clist[v->Corners[j]]->t;
					}
				}
			}
		}
	}
}


extern void AddToBoundaryEdgeList(Edge *cur_e);
extern int num_cycleedges;
/*
New routine for building the edge list based on the input triangle strip
Here, we do not require the 'inDesignCellCycle' being set
02/08/07
*/
void build_Boundary_Edgelist(int *triangle, int ntris)
{
	int i, j;
	Edge *cur_e;
	Face *face;

	/*we need to make sure other routines will not change the content
	in Cycle_edge first 03/03/07 */

	for(i = 0; i < num_cycleedges; i++)
	{
		cur_e = Cycle_edge[i];
		cur_e->OnBoundary = 0;
		cur_e->BoundaryVisited = 0;
	}

	num_cycleedges = 0;

	for(i = 0; i < ntris; i++)
	{
		face = Object.flist[triangle[i]];

		for(j = 0; j < 3; j++)
		{
			cur_e = face->edges[j];

			if(IsRepeated(triangle, cur_e->tris[0], ntris)
				&& IsRepeated(triangle, cur_e->tris[1], ntris))
				continue;

			AddToBoundaryEdgeList(cur_e);
			cur_e->OnBoundary = 1;  //Do we need that now? Yes we do
		}
	}
}


//////////////////
void get_Boundary_2(int *tri_strip, int ntris_in_strip)
{
	////First we need to perform manifold testing and correct the original region
	int i;
	int EndVertID;
	//Face *face;
	Vertex *cur_vert;
	Edge *cur_edge;

	Edge **temp_edgelist = (Edge **) malloc(sizeof(Edge *) * MaxEdgeInCycle);

	int temp_num_edges, other_num_edges;
	temp_num_edges = other_num_edges = 0;

	////rebuild boundary again
	//for(i = 0; i < num_cycleedges; i++)
	//{
	//	Cycle_edge[i]->OnBoundary = 0;
	//}

	//for(i = 0; i < Object.nfaces; i++)
	//{
	//	face = Object.flist[i];

	//	for(int j = 0; j < 3; j++)
	//	{
	//		cur_edge = face->edges[j];
	//		cur_edge->OnBoundary = 0;
	//		cur_edge->BoundaryVisited = 0;
	//	}
	//}

	/* We only need to reset the triangles in the previous strip 01/22/07*/
	//reset_Edge_Flags();

	//BuildBoundaryEdgeList(tri_strip, ntris_in_strip);
	build_Boundary_Edgelist(tri_strip, ntris_in_strip);
	fill_Gap_In_Ring(Cycle_edge, num_cycleedges, tri_strip, ntris_in_strip);
	build_Boundary_Edgelist(tri_strip, ntris_in_strip);

	/* we need to find out all the non-manifold vertices and make some complement */

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

	int count = 0; //avoid dead loop

	while(cur_vert->VertID != EndVertID && count < Object.nfaces)   ////Not form a closed edges list
	{
		for(i = 0; i < cur_vert->Num_edge; i++)
		{
			cur_edge = cur_vert->edges[i];
			if(cur_edge->OnBoundary == 1 && cur_edge->BoundaryVisited == 0) //this may need to change
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
		
		count++;
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


//end of Taustep.cpp


extern SCCList scclist;

/* Test code 02/04/07 */
void write_down_SCC()
{
	int i, j, k;
	int node, anode;

	int countedges = 0;

	FILE *fp = fopen("thescc.txt", "w");
	for(i = 0; i < scclist.num_sccs; i++)
	{
		if(scclist.scccomponents[i].num_nodes <= 2 && scclist.scccomponents[i].num_singularities <= 0)
			continue;

		else if(scclist.scccomponents[i].num_nodes > 2 && scclist.scccomponents[i].valid == 0)
			continue;
		
		
		fprintf(fp, "#nodes:%d\n",scclist.scccomponents[i].num_nodes);

		for(j = 0; j < scclist.scccomponents[i].num_nodes; j++)
		{
			node = scclist.scccomponents[i].nodes[j];
			countedges += sccnodes[node].nedges;
		}

		fprintf(fp, "#edges:%d\n", countedges);

		for(j = 0; j < scclist.scccomponents[i].num_nodes; j++)
		{
			node = scclist.scccomponents[i].nodes[j];
			
			int countnodeedges = 0;
			for(k = 0; k < sccnodes[node].nedges; k++)
			{
				anode = sccedges[sccnodes[node].edges[k]].node_index2;
				if(anode != node)
					countnodeedges++;
			}
			
			fprintf(fp, "the node: %d has %d edges\n", node, countnodeedges);

			for(k = 0; k < sccnodes[node].nedges; k++)
			{
				anode = sccedges[sccnodes[node].edges[k]].node_index2;
				if(anode != node)
				fprintf(fp, "%d \n", anode);
			}

			//fprintf(fp, "\n");
		}
	}


	fclose(fp);
}


/*
Test the tau distortion
*/
double middlex[3], middley[3];
void trace_and_save(int tri, double tau)
{
	if(tri < 0)
		return;

	int i;
	Face *face = Object.flist[tri];
	Vertex *v;
	double st[2];
	int temp_tri;

	for(i = 0; i < 3; i++)
	{
		v = Object.vlist[face->verts[i]];
		st[0] = v->x;
		st[1] = v->y;
	
		trace_Ver(tri, st, newP[i].entry, temp_tri, tau, 0);
	}

	for(i = 0; i < 3; i++)
	{
		newPx[2*i] = newP[i].entry[0];
		newPy[2*i] = newP[i].entry[1];
	}

	/*we trace the middle point of the edges of the triangle */
	Vertex *vn;
	for(i = 0; i < 3; i++)
	{
		v  = Object.vlist[face->verts[i]];
		vn = Object.vlist[face->verts[(i+1)%3]];

		middlex[i] = st[0] = (v->x + vn->x)/2.;
		middley[i] = st[1] = (v->y + vn->y)/2.;
		
		trace_Ver(tri, st, newP[i].entry, temp_tri, tau, 0);
	}
	
	for(i = 0; i < 3; i++)
	{
		newPx[2*i+1] = newP[i].entry[0];
		newPy[2*i+1] = newP[i].entry[1];
	}

	/* trace the center of the triangle */
	Vertex *vp;
	vp = Object.vlist[face->verts[0]];
	v  = Object.vlist[face->verts[1]];
	vn = Object.vlist[face->verts[2]];

	st[0] = (v->x + vp->x + vn->x)/3.;
	st[1] = (v->y + vp->y + vn->y)/3.;
		
	trace_Ver(tri, st, newP[0].entry, temp_tri, tau, 0);

	newPx[6] = newP[0].entry[0];
	newPy[6] = newP[0].entry[1];
}


/**********************************************************************/
/*further tuning tau 03/13/07*/


/*Judge whether a region is desired good region (Morse set)
Current criterion is the number of fixed points inside the region
!!We assume the fixed point extraction has been performed and the result has been
saved for this region!!
Good: true;   has only one fixed point
Bad:  false;  has more than one fixed point
*/

bool is_good_region(int scc_index)
{
	/*if the region is too large*/
	if(scclist.scccomponents[scc_index].num_nodes > (int)num_sccnodes/4)
		return false;

	/*if the region contains more than more fixed points*/
	if(scclist.scccomponents[scc_index].num_singularities > 1
		&& scclist.scccomponents[scc_index].num_nodes > 2)
		return false;

	return true;
}

/*
If a region is one-ring like region and it is classified as saddle-like region
we know there at least one opposite periodic orbits in there, therefore
we need to perform larger \tau tracing
(if the region has one triangle that all its three vertices are on the boundary,
then there is no hope that we can separate the two regions!!!!!)
*/
bool is_good_one_ring_region(int scc_index)
{
	if(CalEulerValue(scclist.scccomponents[scc_index].nodes, 
		scclist.scccomponents[scc_index].num_nodes)==1)
		return true;   /*it is a disk*/
	
	if(mcgnodes[scclist.scccomponents[scc_index].node_index].type == 2)
		return false;

	return true;
}


/*03/03/07*/

/*
Judge whether the directed edge "edge_index" is the outgoing edge
of the node "snode"
*/
bool is_outedge(int edge_index, int snode)
{
	if(sccedges[edge_index].node_index1 == snode) return true;
	return false;
}

/*
delete one edge from the edge list
*/

bool del_one_edge_from(int *edges, int &n, int del_e)
{
	int i, pos;

	for(i = 0; i < n; i++)
	{
		if(edges[i] == del_e)
		{
			pos = i;
			break;
		}
	}

	if(i >= n)
		return false;

	for(i = pos; i < n-1; i++)
	{
		edges[i] = edges[i+1];
	}

	n--;
	return true;
}



/*
Remove all the outgoing edges of a specified node 
*/
void remove_edges_from(int node)
{
	int i, othernode;

	for(i = 0; i < sccnodes[node].nedges; i++)
	{
		//if(is_outedge(sccnodes[node].edges[i], node))
		//{
			/*delete the edge from the edge list */
			del_one_edge_from(sccnodes[node].edges, sccnodes[node].nedges, sccnodes[node].edges[i]);

			othernode = sccedges[sccnodes[node].edges[i]].node_index1;

			if(othernode == node)
				othernode = sccedges[sccnodes[node].edges[i]].node_index2;

			/*delete the same edge from the ending node connected by the edge */
			del_one_edge_from(sccnodes[othernode].edges, sccnodes[othernode].nedges, sccnodes[node].edges[i]);

			sccedges[sccnodes[node].edges[i]].cancelled = true;
			num_sccedges--;
		//}
	}
}

/*
Remove all the edges of the nodes inside the specified scc component
*/
void remove_edges_from_scc(int scc_index)
{
	int i;
	
	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		remove_edges_from(scclist.scccomponents[scc_index].nodes[i]);
	}
}

/*
remove edges from all the bad region
*/

void remove_edges_from_toolarge_regions()
{
	int i; 

	for(i = 0; i < num_sccomps; i++)
	{
		if(scclist.scccomponents[i].bad_flag)
		{
			remove_edges_from_scc(i);
		}
	}
}


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/*new routine for removing edges from the large region
---------------------------------------------------------03/28/07
*/

void remove_edges_from_2(int node)
{
	int i/*, othernode*/;

	for(i = 0; i < sccnodes[node].nedges; i++)
	{
		/*delete the edge from the edge list */
		/*Mark all the edges to be deleted*/
		sccedges[sccnodes[node].edges[i]].cancelled = true;
		num_sccedges--;
	}
}


void remove_edges_from_scc_2(int scc_index)
{
	int i;
	
	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		remove_edges_from_2(scclist.scccomponents[scc_index].nodes[i]);
	}
}

void release_edgelist_node_2(int node)
{
	if(sccnodes[node].edges != NULL)
		free(sccnodes[node].edges);
	sccnodes[node].edges = NULL;
	sccnodes[node].nedges = 0;
}

void release_edgelists_scc_2(int scc_index)
{
	int i;
	
	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		release_edgelist_node_2(scclist.scccomponents[scc_index].nodes[i]);
	}
}



void remove_edges_from_toolarge_regions_2()
{
	int i; 

	for(i = 0; i < num_sccomps; i++)
	{
		if(scclist.scccomponents[i].bad_flag)
		{
			remove_edges_from_scc_2(i);

			/*release all the edge list of the nodes*/
			release_edgelists_scc_2(i);
		}
	}

	

}

/*
After getting all the Morse sets(regions) and detect the fixed points in 
each region, we track only the saddle-like region and see
whether it is a good region or not
*/
void find_toolarge_regions()
{
	int i;

	for(i = 0; i < num_sccomps; i++)
	{
		scclist.scccomponents[i].bad_flag = false;
		if(scclist.scccomponents[i].num_nodes<=2 && scclist.scccomponents[i].num_singularities == 0)
			continue;

		if(scclist.scccomponents[i].num_nodes>2 && scclist.scccomponents[i].valid == 0)
			continue;


		if(/*!is_good_one_ring_region(i)||*/!is_good_region(i))
		{
			scclist.scccomponents[i].bad_flag = true;
			num_toolarge_regions++;
		}

	}
}


/*
We need to save the previous regions for convergence determinance
*/
typedef struct Toolarge_Region{
	int pre_scc_index;
	int num;
}Toolarge_Region;
Toolarge_Region *num_tris_in_toolarge_regions;
int cur_num_tris_in_toolarge_regions_id;

/*
We still need to keep track of the special triangle selected from the original SCC regions !!!
*/
void make_copy_num_tri_in_toolarge_regions()
{
	num_tris_in_toolarge_regions = (Toolarge_Region*)malloc(sizeof(Toolarge_Region)*num_toolarge_regions);

	if(num_tris_in_toolarge_regions == NULL)
	{
		time ( &g_rawtime );
		FILE *fp = fopen("mem_error.txt", "a");
		fprintf(fp, "fail to allocate memory for num_tris_in_bad_regionis in routine \
			make_copy_num_tri_in_bad_regions\n");
		fclose(fp);
		fprintf(fp, "current date and time are :  %s. \n", ctime (&g_rawtime) );
		exit(-1);
	}

	int i;
	cur_num_tris_in_toolarge_regions_id = 0;

	for(i = 0; i < num_sccomps; i++)
	{
		if(scclist.scccomponents[i].bad_flag)
		{
			num_tris_in_toolarge_regions[cur_num_tris_in_toolarge_regions_id].num = 
				scclist.scccomponents[i].num_nodes;

			num_tris_in_toolarge_regions[cur_num_tris_in_toolarge_regions_id].pre_scc_index = i;
		}
	}
}


void get_boundary_vers_for_SCC(int scc_index)
{
	/*build the boundary edge list*/
	build_Boundary_Edgelist(scclist.scccomponents[scc_index].nodes,
		scclist.scccomponents[scc_index].num_nodes);

	/*mark all the boundary vertices*/
	int i;
	Edge *e;
	for(i = 0; i < num_cycleedges; i++)
	{
		e = Cycle_edge[i];

		Object.vlist[e->verts[0]]->OnBoundary = 1;
		Object.vlist[e->verts[1]]->OnBoundary = 1;

	}

}

/*
find out one unused edge from the edge list and return its pointer and index
*/
GraphEdge * find_unused_edge(GraphEdge *oneunused, int &index, bool &flag)
{
	int i;
	flag = false;

	for(i = 0; i < pre_num_sccedges; i++)
	{
		if(sccedges[i].cancelled == true)
		{
			oneunused = &sccedges[i];
			index = oneunused->edge_index;
			flag = true;
			return oneunused;
		}
	}

	oneunused = NULL;
	index = -1;
	return oneunused;
}

/*
Add the new edges to the sccedges and try to reuse the empty slot
*/
void add_new_edge_reuse(int tri, int backward)
{
	int i, sccnodeid;
	//Vertex *v;
	//Corner *c;

	if(tri >= Object.nfaces)
	{
		//for nodes corresponding to boundary edges

		/* we first need to search the index of the node according to the input tri */
		for(i = Object.nfaces; i < num_sccnodes; i++)
		{
			if(sccnodes[i].node_index == tri)
			{
				sccnodeid = i;
				break;
			}
		}
		
	}

	else{
		sccnodeid = tri;
	}

	for(i = 0; i < ntris_in_strip; i++)
	{
		/*if the node point to itself, we need to check whether it is a boundary triangle or not
		if it is, we don't add the edge! 02/25/07 */
		if(sccnodeid == tri_strip[i])
		{
			if(is_boundary_tri(sccnodeid))
				continue;
		}

		/*we don't add repeated edges*/
		if(has_Edge_From_To(sccnodeid, tri_strip[i]))
			continue;

		/* If there is no link between the two triangles, add a new edge */
		GraphEdge *oneunused = NULL;
		int edge_index = -1;
		bool flag = false;
		oneunused = find_unused_edge(oneunused, edge_index, flag);
		if(flag)
		{
			if(backward == 0) /*_forward tracing*/
				add_to_used_edge(sccnodeid, tri_strip[i], 
					num_sccedges, oneunused);
			else
				add_to_used_edge(tri_strip[i], sccnodeid,
					num_sccedges, oneunused);
		}
		else
		{
			if(backward == 0)
				SCC_AddToEdge(sccnodeid, tri_strip[i], num_sccedges);
			else
				SCC_AddToEdge(tri_strip[i], sccnodeid, num_sccedges);
			edge_index = num_sccedges - 1;
		}

		/*add the edge into the node's edge list*/
		SCC_AddEdgeToNode(sccnodeid, edge_index);
		SCC_AddEdgeToNode(tri_strip[i], edge_index);
	}
}

/*
rebuild the subgraph related to the specified triangle used new tracing result
*/
void rebuild_subgraph_tri(int tri)
{
	/* Trace the three vertices to get the corresponding mapping positions */
	Face *face = Object.flist[tri];
	Vertex *v;
	int i;
	for(i = 0; i < 3; i++)
	{
		v = Object.vlist[face->verts[i]];

		/*------------------------------------------------------*/
		//we try reuse the tracing result from each vertex
		/*use the information saved in the new member of Vertex*/
		newP[i].entry[0] = v->endp[0].entry[0];
		newP[i].entry[1] = v->endp[0].entry[1];
		en_tris[i] = v->end_tri[0];
	}

	/*get the set of image triangles*/
	get_Covered_Tris();

	/*add to the graph, we need to reuse the previous 'canceled' edges*/
}

/*
Rebuild the images of the triangles in the specified SCC
*/
void rebuild_subgraph_scc(int scc_index)
{
	int i;

	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
	}
}

/*
Adjust the \tau of the region once (for the interior vertices only!!!03/15/07)
*/
void retrace_in_region(int scc_index, int backward)
{
	/*we first find out all the interior vertices of the region*/

	get_boundary_vers_for_SCC(scc_index);

	int i, j;
	Face *face;
	Vertex *v;

	/*set the flag of all interior vertices*/
	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		face = Object.flist[scclist.scccomponents[scc_index].nodes[i]];

		for(j = 0; j < 3; j++)
		{
			Object.vlist[face->verts[j]]->visited = false;
		}
	}

	/*starting adjust \tau tracing here*/
	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		face = Object.flist[scclist.scccomponents[scc_index].nodes[i]];

		for(j = 0; j < 3; j++)
		{
			v = Object.vlist[face->verts[j]];

			if(v->visited || v->OnBoundary == 1)
				continue;  /*the vertex has been traced or on the boundary*/

			v->visited = true;

			/*otherwise, we find a mean to get a proper \tau for it*/
			//double tau = v->tau[0]+2*tan(v->length[2]); //???Can we use the curl information here to adjust that

			//if(tau > 300 && v->tau[0] > 299) /*some upper limit of the global tau*/ 
			//	continue; /*give up*/
			//else
			//	tau = 300;

			/*now we are using new idea, the tau could be relatively small 03/21/07*/
			double tau;
			
			if(backward == 0)
				tau = v->tau[0]+.1; //???Can we use the curl information here to adjust that
			else
				tau = v->tau[1]+.1;

			if(tau > 20)
			{
				if(backward == 0)
				{
					if(	v->tau[0] > 19.75) /*some upper limit of the global tau*/
						tau = 20;
					else
						continue; /*give up*/
				}
				else
				{
					if(	v->tau[1] > 19.75) /*some upper limit of the global tau*/
						tau = 20;
					else
						continue; /*give up*/
				}
			}

			/*we don't consider the decreasing of tau now
			it means that we always increase tau when necessary
			03/21/07*/
			if(backward == 0) /*_forward tracing*/
			{
				//trace_Ver(v->end_tri[0], v->endp[0].entry, v->endp[0].entry, v->end_tri[0],
				//	(tau-v->tau[0]), backward);

				/*how about we do not use previous result to do tracing*/
				int tri_id = -1;
				TriangleThroughVertex(face->verts[j], tri_id, backward);

				if(tri_id < 0)
				{
					tri_id = Object.clist[v->Corners[0]]->t;
				}

				double st[2] = {v->x, v->y};
				trace_Ver(tri_id, st, v->endp[0].entry, v->end_tri[0],
					tau, backward);

				v->tau[0] += tau;
			}
			else{
				//trace_Ver(v->end_tri[1], v->endp[1].entry, v->endp[1].entry, v->end_tri[1],
				//	(tau-v->tau[1]), backward);
				
				/*how about we do use previous result to do tracing*/
				int tri_id = -1;
				TriangleThroughVertex(face->verts[j], tri_id, backward);

				if(tri_id < 0)
				{
					tri_id = Object.clist[v->Corners[0]]->t;
				}

				double st[2] = {v->x, v->y};
				trace_Ver(tri_id, st, v->endp[1].entry, v->end_tri[1],
					tau, backward);
				v->tau[1] += tau;
			}
		}
	}

	/*reset the flags for vertices*/
	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		int j;
		face = Object.flist[scclist.scccomponents[scc_index].nodes[i]];

		for(j = 0; j < face->nverts; j++)
		{
			Object.vlist[face->verts[j]]->visited = false;
		}
	}

	/*reconstruct image only for these triangles using new idea03/21/07*/
	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		int j;
		face = Object.flist[scclist.scccomponents[scc_index].nodes[i]];

		for(j = 0; j < face->nverts; j++)
		{
			if(Object.vlist[face->verts[j]]->OnBoundary == 1
				||Object.vlist[face->verts[j]]->visited)
				continue;

			Object.vlist[face->verts[j]]->visited = true;

			build_edges_Ver(face->verts[j], backward);
		}
	}

	/*
	According to the type of tracing to add the extra edges
	*/
	Edge *e;
	double st1[2], st2[2]/*, tau1, tau2*/;
	int neighbor_tri;
	Vertex *v1, *v2;

	/*Do the adaptive sampling along inner edges here*/
	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		face = Object.flist[scclist.scccomponents[scc_index].nodes[i]];

		for(j = 0; j < 3; j++)
		{
			e = face->edges[j];

			if(e->visited == 1)
				continue;

			if(Object.vlist[e->verts[0]]->OnBoundary == 1
				&& Object.vlist[e->verts[1]]->OnBoundary == 1)
				continue;

			e->visited = 1;

			v1 = Object.vlist[e->verts[0]];
			v2 = Object.vlist[e->verts[1]];

			/*currently, we don't deal with the boundary edges!!!!*/
			if(backward == 0)
			{
				if(v1->end_tri[0] < 0 || v2->end_tri[0] < 0)
					continue;
			}
			else
			{
				if(v1->end_tri[1] < 0 || v2->end_tri[1] < 0)
					continue;
			}


			/*if the images of the two vertices are already continuous, need not
			consider this edge any more*/
			if(backward == 0)
			{
				if(v1->end_tri[0] == v2->end_tri[0]
					|| are_close_neighbors(v1->end_tri[0], v2->end_tri[0]))
				continue;
			}
			else
			{
				if(v1->end_tri[1] == v2->end_tri[1]
					|| are_close_neighbors(v1->end_tri[1], v2->end_tri[1]))
				continue;
			}
			
			st1[0] = v1->x;
			st1[1] = v1->y;
			st2[0] = v2->x;
			st2[1] = v2->y;

			neighbor_tri = e->tris[0];
			if(neighbor_tri == i)
				neighbor_tri = e->tris[1];

			if(backward == 0)
			{
				trace_an_edge_build_di_edges_adp_interpolate(st1, st2, v1->end_tri[0], v2->end_tri[0], 
					face->index, neighbor_tri, v1->tau[0], v2->tau[0], backward);
			}
			else
			{
				trace_an_edge_build_di_edges_adp_interpolate(st1, st2, v1->end_tri[1], v2->end_tri[1], 
					face->index, neighbor_tri, v1->tau[1], v2->tau[1], backward);
			}

			//trace_an_edge_build_di_edges_adp_interpolate(st1, st2, v1->RegionListID, v2->RegionListID, 
			//	i, neighbor_tri,
			//	v1->tau[0], v2->tau[0], backward);
		}
	}

	/*remove the marker of the boundary vertices*/
	for(i = 0; i < num_cycleedges; i++)
	{
		e = Cycle_edge[i];

		Object.vlist[e->verts[0]]->OnBoundary = 0;
		Object.vlist[e->verts[1]]->OnBoundary = 0;

	}
}


void retrace_toolarge_regions(int backward)
{
	int i;

	for(i = 0; i < num_sccomps; i++)
	{
		if(!scclist.scccomponents[i].bad_flag)
			continue;

		retrace_in_region(i, backward);
	}
}

/*
In this framework, we start from \tau = 0, and get the first set of Morse sets.
Then, we judge whether the corresponding region of each Morse set is good or not.
For bad regions, we mark them, and remove all the related edges of the nodes inside the region
For each bad region, we increase \tau certain amount, each region can have different scale,
this can be decided by the local feature of the flow in the region
Try this larger \tau, 
*/

void init_SCC_for_adaptive_framework()
{
	//initialize the nodes
	for(int i = 0; i < Object.nfaces; i++)
	{
		sccnodes[i].visited = 0;
		sccnodes[i].parent = i;
		sccnodes[i].levels = 0;
		//Object.flist[i]->contain_singularity = 0;
	}
	scclist.scccomponents = NULL;
	scclist.num_sccs = 0;
}

/*
save the first tracing information for a vertex
*/
void save_init_ver_tracing(int vertid)
{
	int tri_id;

	Vertex *v = Object.vlist[vertid];
	
	TriangleThroughVertex(vertid, tri_id, 0);

	if(tri_id < 0)
	{
		tri_id = Object.clist[v->Corners[0]]->t;
	}

	v->tau[0] = v->tau[1] = 0;
	v->end_tri[0] = v->end_tri[1] = tri_id;
	v->endp[0].set(v->x, v->y);
	v->endp[1].set(v->x, v->y);

	v->OnBoundary = 0;

}

/*
save the first tracing information for all vertices 
for the adaptive framework
*/
void save_init_vers()
{
	int i;

	for(i = 0; i < Object.nverts; i++)
	{
		save_init_ver_tracing(i);
	}
}


void adaptive_tau_framework()
{
	/*initialization*/
    InitForSCC();

	/*1. start \tau = 0 */
	//init_all_Vers();
	BuildConnectedGraph(0);
	CaptureSing();

	/*save the initial information for the vertex*/
	save_init_vers();

	/*Find the SCC for this \tau=0 map */
	FindSCC();
	mark_All_Valid_SCCS();

	int count = 0;

	num_toolarge_regions = 0;
	find_toolarge_regions();

	FILE *fp = fopen("adaptive.txt", "a");
	fprintf(fp, "count%d: find %d too large regions. \n", count, num_toolarge_regions);
	fprintf(fp, "count%d: current edges %d. \n", count, num_sccedges);
	fclose(fp);

	/*adaptive step here*/
	while(num_toolarge_regions > 0 && count < 5)
	{
		/*we need to record the length of the edge list*/
		//if(num_sccedges > pre_num_sccedges)
		pre_num_sccedges = cur_end_edgelist;

		/*remove the edges associated with the nodes in those "too large" region*/
		remove_edges_from_toolarge_regions();

		//remove_edges_from_toolarge_regions_2(); /*03/28/07*/

		/*retrace to find the images of the nodes in those "too large" region*/
		/*note that the previous tau should be saved for each vertex*/

		/*we haven't done the adaptive edge sampling here 03/22/07*/
		retrace_toolarge_regions(1);
		retrace_toolarge_regions(0);

		/*extract the SCCs*/
		//finalize_scclist();
		init_SCC_for_adaptive_framework();
		FindSCC();
		mark_All_Valid_SCCS();

		num_toolarge_regions = 0;
		find_toolarge_regions();

		count ++;

		fp = fopen("adaptive.txt", "a");
		fprintf(fp, "count%d: find %d too large regions. \n", count, num_toolarge_regions);
		fprintf(fp, "count%d: current edges %d. \n", count, num_sccedges);
		fclose(fp);
	}
}


void init_all_Vers()
{
	int i;
	for(i = 0; i < Object.nverts; i++)
	{
		Object.vlist[i]->tau[0] = Object.vlist[i]->tau[1] = 0;
	}
}



/*visualize the relationship between spatial tau and temporal tau*/

void init_samplepts_tautracing()
{
	_forward = new TraceSamplePtList();
	backward_spts = new TraceSamplePtList();
}

void finalize_samplepts_tautracing()
{
}


extern void  HsvRgb( float hsv[3], float rgb[3] );


void assign_vertex_color()
{
	float smallest, largest;
	int i;
	Vertex *v;
	float hsv[3] = {0, 1., 1.};
	float rgb[3] = {0.};

	/*1. consider _forward results first*/
	smallest = 1e20;
	largest = 0;

	/*process vertices first*/
	for(i=0; i<Object.nverts; i++)
	{
		v = Object.vlist[i];

		if((v->tau[0]+v->tau[1])>largest) largest = (v->tau[0]+v->tau[1]);
		if((v->tau[0]+v->tau[1])<smallest) smallest = (v->tau[0]+v->tau[1]);
	}

	/*compute the deviation*/
	double mean = (largest+smallest)/2;
	double sd = 0;
	for(i=0; i<Object.nverts; i++)
	{
		double  ttime = Object.vlist[i]->tau[0]+Object.vlist[i]->tau[1];
		sd += (ttime - mean)*(ttime - mean);
	}

	sd = sd/Object.nverts;
	sd = sqrt(sd);

	double min_time = mean - 2*sd;
	double max_time = mean + 2*sd;



	/*assign colors for those points*/
	//double diff = largest-smallest;
	double diff = max_time-min_time;
	for(i=0; i<Object.nverts; i++)
	{
		v = Object.vlist[i];

		if((v->tau[0]+v->tau[1]) > max_time) hsv[0] = 0;
		else if((v->tau[0]+v->tau[1]) < min_time) hsv[0] = 240.;
		else
			hsv[0] = 240 - 240*(v->tau[1]+v->tau[0]-min_time)/diff;
		HsvRgb(hsv, rgb);
		v->b_color[0] = rgb[0];
		v->b_color[1] = rgb[1];
		v->b_color[2] = rgb[2];
	}
}


void assign_color()
{
	float smallest, largest;
	int i;
	Vertex *v;
	float hsv[3] = {0, 1., 1.};
	float rgb[3] = {0.};

	/*1. consider _forward results first*/
	smallest = 1e20;
	largest = 0;

	/*process vertices first*/
	for(i=0; i<Object.nverts; i++)
	{
		v = Object.vlist[i];

		if(v->tau[0]>largest) largest = v->tau[0];
		if(v->tau[0]<smallest) smallest = v->tau[0];
	}

	/*process the points on edges or triangle centers*/
	for(i=0; i<_forward->nelems; i++)
	{
		if(_forward->samplepts[i]->time>largest)
			largest = _forward->samplepts[i]->time;
		if(_forward->samplepts[i]->time<smallest)
			smallest = _forward->samplepts[i]->time;
	}

	/*we let backward and _forward have the same scale*/
	///*assign colors for those points*/
	//double diff = largest-smallest;
	//for(i=0; i<Object.nverts; i++)
	//{
	//	v = Object.vlist[i];

	//	hsv[0] = 256*(v->tau[0]-smallest)/diff;
	//	HsvRgb(hsv, rgb);
	//	v->f_color[0] = rgb[0];
	//	v->f_color[1] = rgb[1];
	//	v->f_color[2] = rgb[2];
	//}

	//for(i=0;i<_forward->nelems;i++)
	//{
	//	hsv[0] = 256*(_forward->samplepts[i]->time-smallest)/diff;
	//	HsvRgb(hsv, rgb);
	//	_forward->samplepts[i]->color[0] = rgb[0];
	//	_forward->samplepts[i]->color[1] = rgb[1];
	//	_forward->samplepts[i]->color[2] = rgb[2];
	//}


	/*2. consider backward results first*/
	//smallest = 1e20;
	//largest = 0;

	/*process vertices first*/
	for(i=0; i<Object.nverts; i++)
	{
		v = Object.vlist[i];

		if(v->tau[1]>largest) largest = v->tau[1];
		if(v->tau[1]<smallest) smallest = v->tau[1];
	}

	/*process the points on edges or triangle centers*/
	for(i=0; i<backward_spts->nelems; i++)
	{
		if(backward_spts->samplepts[i]->time>largest)
			largest = backward_spts->samplepts[i]->time;
		if(backward_spts->samplepts[i]->time<smallest)
			smallest = backward_spts->samplepts[i]->time;
	}

	/*assign colors for those points*/
	double diff = largest-smallest;
	for(i=0; i<Object.nverts; i++)
	{
		v = Object.vlist[i];

		hsv[0] = 256*(v->tau[1]-smallest)/diff;
		HsvRgb(hsv, rgb);
		v->b_color[0] = rgb[0];
		v->b_color[1] = rgb[1];
		v->b_color[2] = rgb[2];
	}

	for(i=0;i<backward_spts->nelems;i++)
	{
		hsv[0] =256*(backward_spts->samplepts[i]->time-smallest)/diff;
		HsvRgb(hsv, rgb);
		backward_spts->samplepts[i]->color[0] = rgb[0];
		backward_spts->samplepts[i]->color[1] = rgb[1];
		backward_spts->samplepts[i]->color[2] = rgb[2];
	}

	/*assign color for _forward tracing samples*/
	for(i=0; i<Object.nverts; i++)
	{
		v = Object.vlist[i];

		hsv[0] = 256*(v->tau[0]-smallest)/diff;
		HsvRgb(hsv, rgb);
		v->f_color[0] = rgb[0];
		v->f_color[1] = rgb[1];
		v->f_color[2] = rgb[2];
	}

	for(i=0;i<_forward->nelems;i++)
	{
		hsv[0] = 256*(_forward->samplepts[i]->time-smallest)/diff;
		HsvRgb(hsv, rgb);
		_forward->samplepts[i]->color[0] = rgb[0];
		_forward->samplepts[i]->color[1] = rgb[1];
		_forward->samplepts[i]->color[2] = rgb[2];
	}

}


/*This routine compute the number of sample points
associated with each vertex
It will be used to visualize the density of the sample points in the domain
*/
void compute_density()
{
	/**/
	int i;
	Edge *e;
	int edge_index;
	int count = 0;
	for(i=0; i<Object.nverts; i++)
	{
		Object.vlist[i]->s_count = 2;
	}

	/*consider _forward samples*/
	for(i=0; i<_forward->nelems; i++)
	{
		edge_index = _forward->samplepts[i]->edge;
		e = Object.edgelist[edge_index];
		Object.vlist[e->verts[0]]->s_count ++;
		Object.vlist[e->verts[1]]->s_count ++;
	}

	/*consider backward samples*/
	for(i=0; i<backward_spts->nelems; i++)
	{
		edge_index = backward_spts->samplepts[i]->edge;
		e = Object.edgelist[edge_index];
		Object.vlist[e->verts[0]]->s_count ++;
		Object.vlist[e->verts[1]]->s_count ++;
	}
}

/*assign the color to each vertex according to the 
number of samples associated with it*/
void assign_density_colors()
{
	int largest, smallest;
	largest = 2;
	smallest = 987654;
	int i;
	float hsv[3] = {0, 1, 1};
	float rgb[3] = {0.};
	
	/*get the maximum number and minimum number*/
	for(i=0; i<Object.nverts; i++)
	{
		if(Object.vlist[i]->s_count > largest) largest = Object.vlist[i]->s_count;
		if(Object.vlist[i]->s_count < smallest) smallest = Object.vlist[i]->s_count;
	}

	/*compute the deviation*/
	int mean = (largest+smallest)/2;
	double sd = 0;
	for(i=0; i<Object.nverts; i++)
	{
		int  vecLen = Object.vlist[i]->s_count;
		sd += (vecLen - mean)*(vecLen - mean);
	}

	sd = sd/Object.nverts;
	sd = sqrt(sd);

	double min_count = mean - 1.5*sd;
	double max_count = mean + 0.5*sd;


	/*assign color accordingly*/
	for(i=0; i<Object.nverts; i++)
	{
		if(Object.vlist[i]->s_count > max_count) hsv[0] = 0;
		else if(Object.vlist[i]->s_count < min_count) hsv[0] = 240.;
		else
			hsv[0] = 240.-240*(Object.vlist[i]->s_count-min_count)*1./((max_count-min_count)*1.);
		HsvRgb(hsv, rgb);
		Object.vlist[i]->d_color[0] = rgb[0];
		Object.vlist[i]->d_color[1] = rgb[1];
		Object.vlist[i]->d_color[2] = rgb[2];
	}
}


/****---------------------------------------------------------****/
/****---------------------------------------------------------****/
/****---------------------------------------------------------****/
/*Konstantin's accurate tau map calculation 07/25/07*/


/*we need a data structure to store the images of the samples of
an edge in order*/

#include "GlView.h"

//AnEdgeSampleImg *edgeSampleImgList;
extern CGlView *g_pclGlView;
EdgeSampleImgList *t_AnedgeSampleImgList;
EdgeInfoLookUp *edgeSampleImgList, *lookuptable;

void accurate_trace_an_edge(int v1, int v2, int tri, double tau, double s_tau)
{
	/*create a stack using linesegment structure*/
	TauLineSeg **stack_lines = (TauLineSeg**)malloc(sizeof(TauLineSeg*)*10);
	for(int i=0; i<10; i++) /*allocate space for each element in the stack*/
	{
		stack_lines[i] = (TauLineSeg*)malloc(sizeof(TauLineSeg));
	}

	int nlines = 0;
	int curMaxNumStackLines = 10;
	int tri1, en_tri1;
	double st1[2] = {0.};
	double en1[2] = {0.};
	int tri2, en_tri2;
	double st2[2] = {0.};
	double en2[2] = {0.};
	double middle_p[2] = {0.};

	icVector2 dis;
	int i;

	double g_g_dt = 0;

	bool last_line = false;

	/*trace the two vertices of current line segment, consider _forward tracing only*/
	/*1. first step*/
	Vertex *v = Object.vlist[v1];
	st1[0]=v->x;
	st1[1]=v->y;

	v = Object.vlist[v2];
	st2[0]=v->x;
	st2[1]=v->y;

	tri1 = tri2 = tri;

	dis.entry[0] = st1[0]-st2[0];
	dis.entry[1] = st1[1]-st2[1];

	double edge_length = length(dis);

	/*2. Loop to find the more accurate image of the edge*/
	do{

		/*for one line segment*/
		while(g_g_dt < tau)
		{
			/*V1*/
			trace_Ver(tri1, st1, en1, en_tri1, s_tau, 0);
			/*V2*/
			trace_Ver(tri2, st2, en2, en_tri2, s_tau, 0);

			g_g_dt += s_tau; /*!!or g_g_dt += g_dt !!*/

			if(g_g_dt > tau)
			{
				/*we need to record the info of the two ending points*/
				/*we also need to avoid redundant information!*/

				t_AnedgeSampleImgList->sampleimgs[t_AnedgeSampleImgList->nsampleimgs] = 
					//new AnEdgeSampleImg();
					(AnEdgeSampleImg*)malloc(sizeof(AnEdgeSampleImg));
				t_AnedgeSampleImgList->sampleimgs[t_AnedgeSampleImgList->nsampleimgs]->l_pos[0]
					= en1[0];
				t_AnedgeSampleImgList->sampleimgs[t_AnedgeSampleImgList->nsampleimgs]->l_pos[1]
					= en1[1];
				t_AnedgeSampleImgList->sampleimgs[t_AnedgeSampleImgList->nsampleimgs]->trinagle
					= en_tri1;
				t_AnedgeSampleImgList->nsampleimgs++;

				if(t_AnedgeSampleImgList->nsampleimgs == 1)
				{
					int test = 0;
				}

				/*if it is the last line segment of the edge, we also need to save the second ending point*/
				if(last_line)
				{
					t_AnedgeSampleImgList->sampleimgs[t_AnedgeSampleImgList->nsampleimgs] = 
						//new AnEdgeSampleImg();
					(AnEdgeSampleImg*)malloc(sizeof(AnEdgeSampleImg));

					t_AnedgeSampleImgList->sampleimgs[t_AnedgeSampleImgList->nsampleimgs]->l_pos[0]
						= en2[0];
					t_AnedgeSampleImgList->sampleimgs[t_AnedgeSampleImgList->nsampleimgs]->l_pos[1]
						= en2[1];
					t_AnedgeSampleImgList->sampleimgs[t_AnedgeSampleImgList->nsampleimgs]->trinagle
						= en_tri2;
					t_AnedgeSampleImgList->nsampleimgs++;
					return;
				}

				if(	t_AnedgeSampleImgList->nsampleimgs >= t_AnedgeSampleImgList->curMaxNum)
				{
					/*extend*/
					if(!t_AnedgeSampleImgList->extend(10))
						exit(-1);
				}


			}

			dis.entry[0] = en1[0]-en2[0];
			dis.entry[1] = en1[1]-en2[1];

			if(length(dis)<=edge_length)
			{
				st1[0] = en1[0];
				st1[1] = en1[1];
				tri1 = en_tri1;
				st2[0] = en2[0];
				st2[1] = en2[1];
				tri2 = en_tri2;
				continue;
			}
			else /*we need to split the line segment into two parts*/
			{
				/*compute the middle point*/
				middle_p[0] = (en1[0]+en2[0])/2.;
				middle_p[1] = (en1[1]+en2[1])/2.;

				/*find the triangle contains the middle point*/
				//g_pclGlView->HitProcessforSelectUnderneathMesh(middle_p[0], middle_p[1]);

				/*we can use the two ending triangles to find out the triangle
				containing the middle point of the line segment connecting the 
				images of the two samples,
				it means that we need to do tracing from one end point to the other
				along the line segment*/
				curMaxTrisInStrip = 20;
				tri_strip = (int*)malloc(sizeof(int)*curMaxTrisInStrip);
				ntris_in_strip = 0;

				{
					get_A_Tri_Strip(en1, en_tri1, en2, en_tri2);
				}

				/*judge whether the middle point falls in one of these triangles*/
				int middle_triangle = en_tri1;
				for(int l = 0; l < ntris_in_strip; l++)
				{
					if(FallIntheTriangle(tri_strip[l], middle_p))
					{
						middle_triangle = tri_strip[l];
						break;
					}
				}

				free(tri_strip);

				/*push the second part of the original line segment into the stack*/
				TauLineSeg *oneline = (TauLineSeg*)malloc(sizeof(TauLineSeg));
				oneline->gend[0] = en2[0];
				oneline->gend[1] = en2[1];
				oneline->tri2 = en_tri2;
				oneline->gstart[0] = middle_p[0];
				oneline->gstart[1] = middle_p[1];
				//oneline->tri1 = g_pclGlView->newelementtriangle;
				oneline->tri1 = middle_triangle;

				stack_lines[nlines] = oneline;
				nlines++;

				/*set the first part of the original line segment as the current line*/

				st2[0] = middle_p[0];
				st2[1] = middle_p[1];
				tri2 = oneline->tri2;

				st1[0] = en1[0];
				st1[1] = en1[1];
				tri1 = en_tri1;

				/*extend the stack*/
				if(nlines>=curMaxNumStackLines)
				{
					TauLineSeg **temp = stack_lines;
					stack_lines = (TauLineSeg**)malloc(sizeof(TauLineSeg*)*(curMaxNumStackLines+10));

					if(stack_lines == NULL)
						exit(-1);
					
					for(i=0; i<curMaxNumStackLines; i++)
						stack_lines[i] = temp[i];

					curMaxNumStackLines+=10;
				}
			}
		}

		/*NOTE: we need to adjust g_dt to be the proper value*/

		g_g_dt -= s_tau; 

		/**/
		if(nlines > 0)
		{
			/*pop up one line segment*/
			st1[0] = stack_lines[nlines-1]->gstart[0];
			st1[1] = stack_lines[nlines-1]->gstart[1];
			tri1 = stack_lines[nlines-1]->tri1;

			st2[0] = stack_lines[nlines-1]->gend[0];
			st2[1] = stack_lines[nlines-1]->gend[1];
			tri2 = stack_lines[nlines-1]->tri2;
			
			nlines--;

			if(nlines <= 0) last_line = true;  /*set the last_line flag*/
		}
	
	}while(nlines>=0);
}

void accurate_construct(double tau, double s_tau)
{
	/*decide a smaller tau, s_tau*/

	int i, j;
	Face *face;
	Edge *e;

	/*initialization*/
	//edgeSampleImgList = new EdgeSampleImgList(3); /*we need at most 3 temporary lists for a triangle*/
	edgeSampleImgList = new EdgeInfoLookUp(3);  //3 edge sample image lists for the closure computation

	lookuptable = new EdgeInfoLookUp();

	/*we need to set up a *look up* table to store the information of the edges
	that haven't been finalized. That is, not both of the two triangles of 
	the edge have been processed.
	*/

	/*consider each triangle*/
	for(i=0; i<Object.nfaces; i++)
	{
		face = Object.flist[i];

		for(j=0; j<3; j++)
		{
			e = face->edges[j];

			/*if it is partially processed, search the look up table*/

			/*currently, we can just start tracing the edge again*/
		}
	}
}


void test_pro(double tau, double s_tau)
{
	int i;

	Edge *e = Object.edgelist[2000];

	/*need to initialize the line segment list*/
	t_AnedgeSampleImgList = new EdgeSampleImgList(10);

	//void accurate_trace_an_edge(int v1, int v2, double tau, double s_tau);

	/*compute the image of this edge*/
	accurate_trace_an_edge(e->verts[0], e->verts[1], e->tris[0], tau, s_tau);
}