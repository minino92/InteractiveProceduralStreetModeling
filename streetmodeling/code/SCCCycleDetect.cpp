/*
SCCCyclceDetect.cpp
In this module, we use strongly connected components to locate the invariant set first;
Then, we find the special points in those *valid* components,
Based on these special points, it is possible to detect all the cycles much quicker
*/

#include "stdafx.h"

#include "scccycledetect.h"

#include "FindSepAttPoints.h"

#include "FindSCC.h"
#include "VFDataStructure.h"
#include "VFAnalysis.h"

#include "LocalTracing.h"

#include "topologyedit.h"

#include "limitcycledetect.h"


extern SCCList scclist;
extern Polygon3D Object;

extern Singularities *singularities;           //being captured singularites' list
extern int cur_singularity_index;
extern LineSeg **trajectories;                 //trajectories' list
extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory
extern Separatrices *separatrices;             //array for group of separatrices


extern int MaxNumTrajectories;

extern int MaxNumLinesegsPerTraj;


////Actually, we need to find another efficient data structure to store the whole information of a limit cycle
////for future operations, such as reshape, moving
extern LimitCycle *limitcycles;
extern int cur_limitcycle_index;
extern int MaxNumLimitCycles;

extern int *cellcycle;
extern int num_celltriangle;                ////number of triangles in the found cell cycle
extern int MaxTriangleInCellCycle;

extern TriangularRegion intersectRegion;     ////The intersect region

extern bool IsMixedEdge(Edge *cur_edge);
extern void StoreCurStreamline(LineSeg *streamline, int num);

extern void UpdateListInSingularity(int singID, int limitcycle);
extern void UpdateListInLimitCycle(int limitcycle, int saddleID);
extern void BuildHandlerforLimitCycle(int cur_limitcycle_index, double fixedpoint[2], int ver);
extern double g_dt; /*02/27/07*/


double ave_length;


double sum_flow_length = 0;

/////////////////////////////////////////////////////////
/*
Judge whether a strongly connected component is valid or not
*/
bool IsValidSCC(int scc_index)
{
	if(scclist.scccomponents[scc_index].num_nodes < 2)
		return false;  //the minimum number of triangles that can contain a limit cycle is 2

	/* this calculation may have problem 07/22/06 */
	int Euler_num = CalEulerValue(scclist.scccomponents[scc_index].nodes, scclist.scccomponents[scc_index].num_nodes);
	scclist.scccomponents[scc_index].num_boundaries = 1 + (2-Euler_num)/2;

	////if the scc forms a topological disk
	if(Euler_num == 1)
	{
		//copy to the intersectRegion
		CopyRegion(scclist.scccomponents[scc_index].nodes, intersectRegion.trianglelist,
			scclist.scccomponents[scc_index].num_nodes);
		intersectRegion.num = scclist.scccomponents[scc_index].num_nodes;

		if(!IntersectedRegionSingCapture())
			return false;  //do not find any singularity inside this SCC
	}

	return true;
}


/*
Calculate the special points (separation or attachment points in those valid SCC's)
*/
void GetSpPtsForValidSCC()
{
	int i, j;
	Face *face;
	Edge *cur_e;

	////first calculate the Jacobian for all triangles and all vertices
	////Later, we can try to calculate those triangles and vertices inside the SCC
	//CalJacobianForWholeMesh(); //Move this call to the "BuildDirectedGraph()"

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		for(j = 0; j < 3; j++)
		{
			cur_e = face->edges[j];
			cur_e->visited = 0;
			cur_e->OnBoundary = 0;
			cur_e->sep.set(0, 0);
			cur_e->attp.set(0, 0);
			cur_e->find_attp = 0;
			cur_e->find_sep = 0;
		}
	}

	for(i = 0; i < scclist.num_sccs; i++)
	{
		if(!IsValidSCC(i))
		{
			scclist.scccomponents[i].valid = 0;
			continue;
		}

		scclist.scccomponents[i].valid = 1;   //this is a valid SCC

		//calculate the separation and attachment points
		CalAllSpecialPointsForASCC(i);
	}
}


/*
Judge whether the two intersections on the same edge means finding a cycle
*/
bool IsValidIntersections(Edge *cur_e, int scc_index)
{
	//if the line segments between the two intersections is small

	if(scclist.scccomponents[scc_index].num_boundaries > 1 && 
		IsMixedEdge(cur_e))
		return false;               //we probably may miss the tiny cycle in a ring shaped SCC 07/22/06

	//if the vector values on the two intersections are not almost parallel to each other

	//first, we need to get the interpolated vector values on the two intersections of the edge
	icVector2 vec1, vec2;
	icVector2 edge_length, dis;
	double ra1, ra2;

	edge_length.entry[0] = Object.vlist[cur_e->verts[1]]->x - Object.vlist[cur_e->verts[0]]->x;
	edge_length.entry[1] = Object.vlist[cur_e->verts[1]]->y - Object.vlist[cur_e->verts[0]]->y;

	dis.entry[0] = cur_e->intersections[0].entry[0] - Object.vlist[cur_e->verts[0]]->x;
	dis.entry[1] = cur_e->intersections[0].entry[1] - Object.vlist[cur_e->verts[0]]->y;
	ra1 = length(dis) / length(edge_length);
	
	dis.entry[0] = cur_e->intersections[1].entry[0] - Object.vlist[cur_e->verts[0]]->x;
	dis.entry[1] = cur_e->intersections[1].entry[1] - Object.vlist[cur_e->verts[0]]->y;
	ra2 = length(dis) / length(edge_length);

	vec1 = (1-ra1)*Object.vlist[cur_e->verts[0]]->vec + ra1*Object.vlist[cur_e->verts[1]]->vec;
	vec2 = (1-ra2)*Object.vlist[cur_e->verts[0]]->vec + ra2*Object.vlist[cur_e->verts[1]]->vec;

	double cross_r = vec1.entry[0]*vec2.entry[1] - vec1.entry[1]*vec2.entry[0];

	//if(abs(cross_r) > 1e-6)  ////modified at 08/11/06
	if(abs(cross_r) > 1e-10)  ////modified at 07/17/07 for the example field 1
		return false;

	return true;
}


/*
It will be called if the same edge has been intersected consectively, which may corresponding
to a tangent point on the edge
*/
bool ApproxTangent(Edge *cur_e, double seg_length, icVector2 new_intersect)
{
	icVector2 approx_dis;
	if(cur_e->num_intersections > 0)
	{
		approx_dis = new_intersect - cur_e->intersections[0];

		if(abs(length(approx_dis)-seg_length) < 1e-4)
			return true;
	}

	return false;
}

/*
During calculating the closed streamline of a detected limit cycle, we need to mark all its neighboring
special points as "dead" points
For attracting limit cycles, mark those close attachment points;
for repelling limit cycles, mark those close separation points.
*/
void MarkNeighborPoints(double p[2], int triangle, int scc_index, int type)
{
	int *NearbyTriangles = NULL;
	int num_triangles = 0;
	int i, j;
	Face *face;
	Edge *cur_e;

	NearbyTriangles = GetDisc(p, triangle, ave_length, 1, NearbyTriangles, num_triangles);
	for(i = 0; i < num_triangles; i++)
	{
		//Do not change the states of those special points not belonging to current SCC
		if(!IsRepeated(scclist.scccomponents[scc_index].nodes, NearbyTriangles[i],
			scclist.scccomponents[scc_index].num_nodes))
			continue;

		face = Object.flist[NearbyTriangles[i]];

		for(j = 0; j < 3; j++)
		{
			cur_e = face->edges[j];

			////No need to find out whether the edge contains a sep / att point or not
			if(type == 0)
				cur_e->sep_visit = 1;
			else
				cur_e->att_visit = 1;

		}
	}
	free(NearbyTriangles);

}


/*
Most part of this routine is similar to the local tracing,
except that we need to return the intersected edge and the intersection on the edge
Here the type still means the forward (0) or backward (1) tracing
*/
extern int globalface;

Edge *g_theedge = NULL;
Edge *g_chosenedge = NULL;

int TraceInTriangleForDetect(double g[2], int &face_id, int type, 
							 Edge *cur_e, icVector2 &intersect, double &seg_length, int &flag )
{
	//similar to the regular tracing, except that you need to mark those special points too close to the cycle
	int i;
	double alpha[3];
	double cur_point[2], pre_point[2];
	double vert0[2];
	icVector2 VP, globalv;

	if(face_id < 0)
		return -1;

	Face *face = Object.flist[face_id];

	Face *pre_f = face;
	
	////initialize
	VP.entry[0] = g[0] - Object.vlist[face->verts[0]]->x;
	VP.entry[1] = g[1] - Object.vlist[face->verts[0]]->y;

	pre_point[0] = cur_point[0] = dot(VP, face->LX);
	pre_point[1] = cur_point[1] = dot(VP, face->LY);

	vert0[0] = Object.vlist[face->verts[0]]->x;   ////for update the global point
	vert0[1] = Object.vlist[face->verts[0]]->y;

	globalface = face_id;

	icVector2 dis;
	seg_length = 0;

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
			/* Calculate the length of the previously obtained line segment */
			dis.entry[0] = cur_point[0] - pre_point[0];
			dis.entry[1] = cur_point[1] - pre_point[1];
			seg_length += length(dis);
			
			pre_point[0] = cur_point[0];
			pre_point[1] = cur_point[1];

			/*change to use other integration scheme 07/09/07*/
			//if(ToNextPoint(pre_point, cur_point, face_id, alpha, type))
			//if(get_nextpt_2ndeuler(pre_point, cur_point, face_id, alpha, type))
			if(compute_next_pt(pre_point, cur_point, face_id, alpha, type))
			{
				////update the global point
				face = Object.flist[face_id];

				globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

				g[0] = vert0[0] + globalv.entry[0];
				g[1] = vert0[1] + globalv.entry[1];

			}

			else{  ////the curve reaches a singularity
				flag = 3;

				return face_id;
			}

			////Get the length
		}

		////3. if the point is out of current triangle
		else{
			double t[2] = {0.};

			int PassVertornot = 0;

			int presave_face = face_id;
			
			CrossVertex2(face_id, cur_point, pre_point, type, PassVertornot);
            
			////update the global point here
			if(PassVertornot > 0)
			{
				////we should not directly use the vertex as next point!!
				////we may move a little bit along the VF direction, but make sure it is still inside
				////the new triangle

				Vertex *PassedVert = Object.vlist[pre_f->verts[PassVertornot-1]];
				g[0] = PassedVert->x;
				g[1] = PassedVert->y;
           

				cur_point[0] = pre_f->xy[PassVertornot-1][0];
				cur_point[1] = pre_f->xy[PassVertornot-1][1];

				cur_e = NULL;  //do not consider passing vertex cases now!!! 07/23/06
			}


			else{
				face_id = presave_face;  

				int which_edge = -1;

				CrossBoundary3(pre_point, cur_point, face_id, alpha, which_edge, t);

				PassEdge(face_id, which_edge);

				////transfer it to the global coordinates
				globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

				intersect.entry[0] = g[0] = vert0[0] + globalv.entry[0];
				intersect.entry[1] = g[1] = vert0[1] + globalv.entry[1];

				////Get the corresponding edge
				for(int k = 0; k < 3; k++)
				{
					int vertindex = face->verts[which_edge];

					cur_e = face->edges[k];
					g_theedge = cur_e;

					if(cur_e->verts[0] != vertindex && cur_e->verts[1] != vertindex)
						break;
				}
			}

			return face_id;
		}
	}
	
	return face_id;
}


/*
Added on 07/05/2007
*/
bool add_To_Cellcycle(int triangle)
{
	if(triangle < 0)
		return false;

	if(num_celltriangle >= MaxTriangleInCellCycle)
	{
		cellcycle = (int*)realloc(cellcycle, sizeof(int)*(MaxTriangleInCellCycle + 50));
		if(cellcycle == NULL)
		{
			MessageBox(NULL, "fail to allocate cellcycle!", "Error", MB_OK);
			exit(-1);
		}

		MaxTriangleInCellCycle += 50;
	}
	
	if(!IsRepeated(cellcycle, triangle, num_celltriangle))
	{
		cellcycle[num_celltriangle] = triangle;
		num_celltriangle++;
		return true;
	}

	return false;
}



bool TraceForDetect(double g[2], int &triangle, int type, int scc_index, Edge *chosen_edge,  int &flag)
{
	int i;
	flag = 0;
	double globalp[2];
	int pre_face, cur_face;

	pre_face = cur_face = triangle;

	globalp[0] = g[0];   globalp[1] = g[1];


	Edge *the_e = NULL;
	icVector2 intersect;
	chosen_edge = NULL;
	Edge *pre_e = NULL;

	double seg_length = 0;
	
	num_celltriangle = 0;
	add_To_Cellcycle(triangle);

	for(i = 0; i < 6*NUMTRACINGTRIANGLE; i++)
	{

		if(cur_face == -1)
		{
			flag = 1;       //fail, reaches boundary of the mesh, no cycle exists on this path
			return false;
		}

		pre_face = cur_face;

		////Here we need a sub routine that can return the intersection and the corresponding edge 
		////when it enters a new triangle

		pre_e = the_e;

		cur_face = TraceInTriangleForDetect(globalp, cur_face, type, the_e, intersect, seg_length, flag); 

		////Move the repeated cycle judgement here 07/25/06
		int pre_limit_index;
		//if(IsPreviousCycle_SCC(scc_index, cur_face, pre_limit_index))
		if(IsPreviousCycle_all(cur_face, pre_limit_index, 1-type))
		{
			flag = 1;
			return false;
		}

		if(flag == 3 || pre_face == cur_face) 
		{
			flag = 1;      //fail, reach singularity or other error, no cycle has been found on this path
			return false;
		}

		if(!IsRepeated(scclist.scccomponents[scc_index].nodes, cur_face, scclist.scccomponents[scc_index].num_nodes))
		{
			flag = 1;     //fail, reach the boundary of the SCC
			return false;
		}

		the_e = g_theedge;

		if(the_e == NULL)  //cross a vertex
			continue;

		/*before 07/05/2007*/

		//if(the_e->num_intersections > 1)
		//{
		//	//do not consider the tangent point now
		//	if(pre_e == the_e)
		//		if(ApproxTangent(the_e, seg_length, intersect))
		//			continue;

		//	//validate the intersections on the edge
		//	the_e->intersections[1] = intersect;        //save the current intersection to the_e->intersections[1]
		//	if(IsValidIntersections(the_e, scc_index))  //this judgement may not be very good for high curled field
		//	{
		//		g_chosenedge = chosen_edge = the_e;
		//		chosen_edge->intersections[1] = intersect;
		//		triangle = cur_face;
		//		return true;                             //we find a cycle here 07/23/06
		//	}
		//}

		//the_e->intersections[0] = intersect;  //save the intersection with the edge
		////the_e->num_intersections = 1;
		//the_e->num_intersections ++;

		/*at 07/05/2007*/
		/* we need to judge whether it forms a loop! 02/08/07 */
		if(!add_To_Cellcycle(cur_face) && num_celltriangle > 1)
		{
			////Move the repeated cycle judgement here 07/25/06
			int pre_limit_index;
			//if(IsPreviousCycle_SCC(scc_index, cur_face, pre_limit_index))
			if(IsPreviousCycle_all(cur_face, pre_limit_index, 1-type))
			{
				flag = 1;
				return false;
			}

			the_e = g_theedge;

			if(the_e == NULL)  //cross a vertex
				continue;

			if(the_e->num_intersections > 1)
			{
				//do not consider the tangent point now
				if(pre_e == the_e)
					if(ApproxTangent(the_e, seg_length, intersect))
						continue;

				//validate the intersections on the edge
				the_e->intersections[1] = intersect;        //save the current intersection to the_e->intersections[1]
				if(IsValidIntersections(the_e, scc_index))  //this judgement may not be very good for highly curled field
				{
					g_chosenedge = chosen_edge = the_e;
					chosen_edge->intersections[1] = intersect;
					triangle = cur_face;
					return true;                             //we find a cycle here 07/23/06
				}
			}

			the_e->intersections[0] = intersect;  //save the intersection with the edge
			//the_e->num_intersections = 1;
			the_e->num_intersections ++;
		}
	}

	return false;  //can not converge!

}


/*
The routine is used to get the cycle from the triangle strip
We may also need to get the streamline starting from the real cycle
New added 02/08/07
*/

bool get_Cycle(int *triangles, int &num, int cur_tri)
{
	if(num == 1)
		return false;

	/* find out cur_tri in the triangles */
	int i;
	for(i = 0; i < num; i++)
	{
		if(cur_tri == triangles[i])
			break;
	}

	if(i >= num)
		return false;

	int n_repeated = num - i;
	int *temp = (int*)malloc(sizeof(int)*(n_repeated+1));

	for(int j = i; j < num; j++)
	{
		temp[j-i] = triangles[j];
	}

	for(i = 0; i < n_repeated; i++)
	{
		triangles[i] = temp[i];
	}
	num = n_repeated;

	free(temp);
	return true;
}


/*
Calculate the closed streamline of the limit cycle according to the input edge and triangle
*/
bool GetClosedStreamline(Edge *the_e, int triangle, int type, int scc_index)
{
	int i;
	int flag = 0;
	double globalp[2];
	int pre_face, cur_face;

	icVector2 curp, dis;

	////Initialization part
	int stop = 1;

	pre_face = cur_face = triangle;

	////We use the latest intersection as the approximate fixed point to start the tracing
	globalp[0] = the_e->intersections[1].entry[0];   
	globalp[1] = the_e->intersections[1].entry[1];

	num_celltriangle = 0;

	////If current trajectories list is not enough to store all the trajectories
	////extend it!
	if(cur_traj_index >= MaxNumTrajectories - 1)
	{
		MaxNumTrajectories += 100;

		num_linesegs_curtraj = (int*)realloc(num_linesegs_curtraj, sizeof(int)*MaxNumTrajectories);
		trajectories = (LineSeg **)realloc(trajectories, sizeof(LineSeg *) * MaxNumTrajectories);

		////extend the old trajectories
		for(i = 0; i < MaxNumTrajectories-100; i++)
			trajectories[i] = (LineSeg *)realloc(trajectories[i], sizeof(LineSeg ) * MaxNumLinesegsPerTraj);

		////allocate the new trajectories
		for( ; i < MaxNumTrajectories; i++)
			trajectories[i] = (LineSeg *)malloc(sizeof(LineSeg) * MaxNumLinesegsPerTraj);

	}

	num_linesegs_curtraj[cur_traj_index] = 0;

	cellcycle[0] = cur_face;
	num_celltriangle = 1;

	/*02/27/07*/
	g_dt = 0;

	for(i = 0; i < 2*NUMTRACINGTRIANGLE/*(int)Object.nfaces/2*/; i++)
	{
		if(cur_face < 0)
			return false;

		pre_face = cur_face;
		cur_face = TraceInATriangle(cur_face, globalp, type, flag); ////0 means always forward tracing here

		if(flag == 3 || flag == 4 || pre_face == cur_face ) 
		{
			return false;
		}

		if(!IsRepeated(cellcycle, cur_face, num_celltriangle))
		{
			////we may need to extend the cell cycle here firstly
            if(num_celltriangle >= MaxTriangleInCellCycle - 1)
			{
				MaxTriangleInCellCycle += 50;
				cellcycle = (int*)realloc((int *)cellcycle, sizeof(int)*MaxTriangleInCellCycle);
			}

			cellcycle[num_celltriangle] = cur_face;
			num_celltriangle++;
		}
		
		else{ ////form a cell cycle, probably we can stop tracing now
			//we need to get the real cycle!
			if(!get_Cycle(cellcycle, num_celltriangle, cur_face))
				return false;

			if(stop >= 4)
				return true;
			stop ++;
		}

		/*
		The current tracing curve should not go outside of current SCC !
		*/
		if(!IsRepeated(scclist.scccomponents[scc_index].nodes, cur_face, scclist.scccomponents[scc_index].num_nodes))
		{
			flag = 1;     //fail, reach the boundary of the SCC
			return false;
		}

		/*
		We test whether it is previous detected cycle again here 07/25/06
		*/
		int pre_limit_index = -1;
		if(IsPreviousCycle_all(cur_face, pre_limit_index, 1-type))
		{
			//we need to build the connection later
			return false;
		}



		////Mark those near by special points having the same type as the limit cycle as "dead" point
		MarkNeighborPoints(globalp, pre_face, scc_index, 1-type);

	}
	return true;
}

/*
The following two routines build the singularity lists for all the SCC components
*/

void BuildSingularityListForASCC(int scc_index)
{
	//first, we need to know how singularity inside the component
	CopyRegion(scclist.scccomponents[scc_index].nodes, intersectRegion.trianglelist,
		scclist.scccomponents[scc_index].num_nodes);
	intersectRegion.num = scclist.scccomponents[scc_index].num_nodes;

	int totalindex = 0;
	int num_sings = IntersectedRegionSingCount(totalindex);

	//Reallocate the space for singular triangle list
	scclist.scccomponents[scc_index].singular_tri = 
		(int*) malloc(sizeof(int)*(num_sings+1));
	scclist.scccomponents[scc_index].num_singularities = 0;

	//
	int i;

	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		if(Object.flist[scclist.scccomponents[scc_index].nodes[i]]->contain_singularity == 0)
			continue;

		//add to the singular triangle list
		scclist.scccomponents[scc_index].singular_tri[scclist.scccomponents[scc_index].num_singularities]=
			scclist.scccomponents[scc_index].nodes[i];
		scclist.scccomponents[scc_index].num_singularities ++;
	}
}

void BuildSingularityListForSCCs()
{
	int i;
	
	for(i = 0; i < scclist.num_sccs; i++)
	{
		if(scclist.scccomponents[i].singular_tri != NULL)
		{
			free(scclist.scccomponents[i].singular_tri);
			scclist.scccomponents[i].singular_tri = NULL;
			scclist.scccomponents[i].num_singularities = 0;
		}

		BuildSingularityListForASCC(i);
	}
}

/*
Find a close point to the singularity but still in the same triangle
We assume the coordinates of the import point are global coordinates
*/
void FindAClosePoint(double old[2], int triangle, double out[2])
{
	Face *face = Object.flist[triangle];
	icVector2 VP;
	VP.entry[0] = old[0] - Object.vlist[face->verts[0]]->x;
	VP.entry[1] = old[1] - Object.vlist[face->verts[0]]->y;

	double a, b, alpha[3];
	a = dot(VP, face->LX);
	b = dot(VP, face->LY);

	Get2DBarycentricFacters(triangle, a, b, alpha);

	/* we first try the center of the triangle */
	out[0] = 0.333*Object.vlist[face->verts[0]]->x + 0.333*Object.vlist[face->verts[1]]->x
		+ 0.334*Object.vlist[face->verts[2]]->x;
	out[1] = 0.333*Object.vlist[face->verts[0]]->y + 0.333*Object.vlist[face->verts[1]]->y
		+ 0.334*Object.vlist[face->verts[2]]->y;

	VP.entry[0] = out[0] - old[0];
	VP.entry[1] = out[1] - old[1];

	////if it is too close to the input point, we move further to close to one of the edge
	if(length(VP) < 1e-6)
	{
		out[0] = 0.01*Object.vlist[face->verts[0]]->x + 0.495*Object.vlist[face->verts[1]]->x
			+ 0.495*Object.vlist[face->verts[2]]->x;
		out[1] = 0.01*Object.vlist[face->verts[0]]->y + 0.495*Object.vlist[face->verts[1]]->y
			+ 0.495*Object.vlist[face->verts[2]]->y;
	}
}

/*
Trace from the neighborhood of a singularity
*/
void TraceFromSingularity(int scc_index, int type)
{
	int i, singID;
	double sing_c[2] = {0.};
	double out[2] = {0.};
	int flag = -1;
	Edge *chosen_edge = NULL;
	int the_triangle;

	for(i = 0; i < scclist.scccomponents[scc_index].num_singularities; i++)
	{
		//get the singularity id
		singID = Object.flist[scclist.scccomponents[scc_index].singular_tri[i]]->singularity_index;

		//if the type of singularity is the same as the input "type", continue

		flag = -1;
		if(type == 0)  //center singularity can not be repeller
		{
			if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
				continue;
			
			ResetEdgeIntersections(scc_index); //reset the intersection information before starting tracing

			sing_c[0] = singularities[singID].gcx;
			sing_c[1] = singularities[singID].gcy;
			FindAClosePoint(sing_c, singularities[singID].Triangle_ID, out);

			//start tracing backward to locate repelling limit cycle if existed
			the_triangle = singularities[singID].Triangle_ID;
			if(TraceForDetect(out, the_triangle, 1, scc_index, chosen_edge, flag))
			{
				//we found a cycle, but we still need to return the intersection and the edge
				//for calculating the closed streamline for the limit cycle

				int pre_limit_index = -1;

				chosen_edge = g_chosenedge;
				if(GetClosedStreamline(chosen_edge, the_triangle, 1, scc_index))
				{
					//store the information for the limit cycle
					StoreCurrentCellCycle(cellcycle, num_celltriangle);  

					limitcycles[cur_limitcycle_index].node_index = -1; 

					StoreCurStreamline(trajectories[cur_traj_index], num_linesegs_curtraj[cur_traj_index]);     

					////store the type of the limit cycle
					limitcycles[cur_limitcycle_index].type = 0;

					////Initialize the connected list for the graph
					limitcycles[cur_limitcycle_index].connected_limitcycle = NULL;
					limitcycles[cur_limitcycle_index].num_connectedcycles = 0;
					limitcycles[cur_limitcycle_index].connected_saddle = NULL;
					limitcycles[cur_limitcycle_index].num_connectedsaddles = 0;
					cur_limitcycle_index++;
				}
			}

			////probably we need to mark its neighboring separation points 07/25/06
			MarkNeighborPoints(out, singularities[singID].Triangle_ID, scc_index, 0);
		}
		
		else  //center singularity can not be attractor
		{
			if(singularities[singID].type == SINK || singularities[singID].type == AFOCUS)
				continue;
			
			ResetEdgeIntersections(scc_index); //reset the intersection information before starting tracing

			sing_c[0] = singularities[singID].gcx;
			sing_c[1] = singularities[singID].gcy;
			FindAClosePoint(sing_c, singularities[singID].Triangle_ID, out);

			//start tracing forward to locate attracting limit cycle if existed
			the_triangle = singularities[singID].Triangle_ID;
			if(TraceForDetect(out, the_triangle, 0, scc_index, chosen_edge, flag))
			{
				//we found a cycle, but we still need to return the intersection and the edge
				//for calculating the closed streamline for the limit cycle

				int pre_limit_index = -1;

				chosen_edge = g_chosenedge;
				if(GetClosedStreamline(chosen_edge, the_triangle, 0, scc_index))
				{
					//store the information for the limit cycle
					StoreCurrentCellCycle(cellcycle, num_celltriangle);  

					limitcycles[cur_limitcycle_index].node_index = -1; 

					StoreCurStreamline(trajectories[cur_traj_index], num_linesegs_curtraj[cur_traj_index]);     

					////store the type of the limit cycle
					limitcycles[cur_limitcycle_index].type = 1;

					////Initialize the connected list for the graph
					limitcycles[cur_limitcycle_index].connected_limitcycle = NULL;
					limitcycles[cur_limitcycle_index].num_connectedcycles = 0;
					limitcycles[cur_limitcycle_index].connected_saddle = NULL;
					limitcycles[cur_limitcycle_index].num_connectedsaddles = 0;
					cur_limitcycle_index++;
				}
			}

			////probably we need to mark its neighboring separation points 07/25/06
			MarkNeighborPoints(out, singularities[singID].Triangle_ID, scc_index, 1);
		}
	}
}

/*
This routine use regular method to detect the limit cycle inside those SCC's without any special points
This routine may need to deal with the embedded cycle cases
*/
void RegularDetect(int scc_index, int type)
{
	//first, decide how many singularities inside the region, where are they, what are their types

	if(scclist.scccomponents[scc_index].num_boundaries == 1)
	{
		//if it is a topological disk, start from center fixed point
		/* we do not perform tracing from fixed points now 02/21/07*/
		//TraceFromSingularity(scc_index, type); //type: 0 -- repelling cycle || 1 -- attracting cycle
	}

	//if it is not a topological disk, pick a point close to boundaries (or on the boundary)
	else
	{
		//set up the "intersectRegion" to get all the boundaries of the region

		//According to the boundary feature to start tracing
		//1) if it is incoming boundary, perform forward tracing
		//2) if it is outgoing boundary, perform backward tracing
	}
}


/*
Judge whether a separation/attachment point on an edge is a valid point or not. It is valid, If it satisfies:
1) not close to the being visited separation point in the same SCC
2) not close to the previously detected limit cycle
We assume that the "cur_e" contain a separation point
Note: we need to set a threshold here. Currently, the threshold is relative to the average edge length
or the Object.radius, we need to make sure it is less than one triangle width
*/

bool IsValidSep(Edge *the_e, int triangle, int scc_index)
{
	// The following setting seems should be put in the tracing code as well

	int *NearbyTriangles = NULL;
	int num_triangles = 0;
	int i, j;
	Face *face;
	Edge *cur_e;

	NearbyTriangles = GetDisc(the_e->sep.entry, triangle, ave_length, 1, NearbyTriangles, num_triangles);

	for(i = 0; i < num_triangles; i++)
	{
		//Do not change the states of those special points not belonging to current SCC
		if(!IsRepeated(scclist.scccomponents[scc_index].nodes, NearbyTriangles[i],
			scclist.scccomponents[scc_index].num_nodes))
			continue;

		face = Object.flist[NearbyTriangles[i]];

		for(j = 0; j < 3; j++)
		{
			cur_e = face->edges[j];

			if(cur_e->valid == 0 || cur_e->find_sep == 0 || cur_e == the_e)
				continue;

			if(cur_e->sep_visit == 1)  //it is a dead point
			{
				the_e->sep_visit = 2;
				return false;
			}

		}
	}
	free(NearbyTriangles);

	return true;
}



bool IsValidAtt(Edge *the_e, int triangle, int scc_index)
{
	int *NearbyTriangles = NULL;
	int num_triangles = 0;
	int i, j;
	Face *face;
	Edge *cur_e;

	NearbyTriangles = GetDisc(the_e->attp.entry, triangle, ave_length, 1, NearbyTriangles, num_triangles);

	for(i = 0; i < num_triangles; i++)
	{
		//Do not change the states of those special points not belonging to current SCC
		if(!IsRepeated(scclist.scccomponents[scc_index].nodes, NearbyTriangles[i],
			scclist.scccomponents[scc_index].num_nodes))
			continue;

		face = Object.flist[NearbyTriangles[i]];

		for(j = 0; j < 3; j++)
		{
			cur_e = face->edges[j];

			if(cur_e->valid == 0 || cur_e->find_attp == 0 || cur_e == the_e)
				continue;

			if(cur_e->att_visit == 1)  //it is a dead point
			{
				the_e->att_visit = 2;
				return false;
			}

		}
	}
	free(NearbyTriangles);

	return true;
}

/*
Get the set of triangles fall inside the distance circle with the center being the input point p
*/
int *GetDisc(double p[3], int triangle, double dsep, double discsize,
			 int *NearbyTriangles, int &num_triangles)
{
	int cur_id;
	Vertex *vert;
	Face *face;
	icVector2 dis;
	Corner *c;

	int i, j;
	
	////Reset the visited flags
	for(i = 0; i < Object.nverts; i++)
	{
		Object.vlist[i]->visited = 0;
		Object.vlist[i]->distance = 1e48;
	}

	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->discard = 0;
	}

	////Find all the triangle
	num_triangles = 0;
	NearbyTriangles = Extend_link(NearbyTriangles, num_triangles);
	NearbyTriangles[0] = triangle;
	num_triangles++;
	cur_id = 0;

	while(cur_id < num_triangles)
	{
		face = Object.flist[NearbyTriangles[cur_id]];

		////It seems the following judgement is not necessary
		if(face->discard == 1)
		{
			cur_id ++;
			continue; 
		}

		
		for(i = 0; i < face->nverts; i++)
		{
			vert = Object.vlist[face->verts[i]];

			if(vert->visited == 1)
				continue;

			dis.entry[0] = vert->x - p[0];
			dis.entry[1] = vert->y - p[1];

			vert->distance = length(dis);

			vert->visited = 1; //set the visited flag
			if(vert->distance > discsize*dsep)
				continue;

			////if the distance between the vertex and the input point is smaller than the threshold
			////we need to add all its adjacent triangles into the triangle list
			for(j = 0; j < vert->Num_corners; j++)
			{
				c = Object.clist[vert->Corners[j]];

				if(c->t < 0) //reach the boundary!
					continue;

				if(IsRepeated(NearbyTriangles, c->t, num_triangles))
					continue;

				NearbyTriangles = Extend_link(NearbyTriangles, num_triangles);
				NearbyTriangles[num_triangles] = c->t;
				num_triangles++;
			}

		}

		face->discard = 1;
		cur_id++;
	}

	return NearbyTriangles;
}

/*
This routine use the information of special points to detect the limit cycle inside those SCC's without any special points
This routine may need to deal with the embedded cycle cases
*/
void SpecialPtsBasedDetect(int scc_index, int type)
{
	//search all the edges inside the SCC[scc_index]
	//if we find an edge containing sep or att point, validate it
	//if it is a valid point, start tracing from it
	//1)for a sep point, trace backward
	//2)for a att point, trace forward

	int i, j;
	Face *face;
	Edge *cur_e;
	int flag = -1;
	Edge *chosen_edge = NULL;
	int the_triangle;

	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		face = Object.flist[scclist.scccomponents[scc_index].nodes[i]];


		for(j = 0; j < 3; j++)
		{
			cur_e = face->edges[j];

			if(cur_e->valid == 0) //not valid edge 
				continue;

			if(type == 0)
			{
				if(cur_e->find_sep == 0)  //no separation point
					continue;

				if(cur_e->sep_visit == 1) //this is a dead point
					continue;

				if(!IsValidSep(cur_e, face->index, scc_index)) //this is not a good separation point
					continue;

				ResetEdgeIntersections(scc_index); //reset the intersection information before starting tracing
				flag = -1;

				//start tracing backward to locate repelling limit cycle if existed
				the_triangle = face->index;
				if(TraceForDetect(cur_e->sep.entry, the_triangle, 1, scc_index, chosen_edge, flag))
				{
					//we found a cycle, but we still need to return the intersection and the edge
					//for calculating the closed streamline for the limit cycle

					chosen_edge = g_chosenedge;
					if(GetClosedStreamline(chosen_edge, the_triangle, 1, scc_index))
					{
						//store the information for the limit cycle
						StoreCurrentCellCycle(cellcycle, num_celltriangle);  

						limitcycles[cur_limitcycle_index].node_index = -1; ////10/16/05

						StoreCurStreamline(trajectories[cur_traj_index], num_linesegs_curtraj[cur_traj_index]);     

						////store the type of the limit cycle
						limitcycles[cur_limitcycle_index].type = 0;

						////Initialize the connected list for the graph
						limitcycles[cur_limitcycle_index].connected_limitcycle = NULL;
						limitcycles[cur_limitcycle_index].num_connectedcycles = 0;
						limitcycles[cur_limitcycle_index].connected_saddle = NULL;
						limitcycles[cur_limitcycle_index].num_connectedsaddles = 0;
						limitcycles[cur_limitcycle_index].connected_l = 0;
						limitcycles[cur_limitcycle_index].connected_r = 0;
						BuildHandlerforLimitCycle(cur_limitcycle_index, 
							chosen_edge->intersections[1].entry, chosen_edge->verts[0]);
						cur_limitcycle_index++;
					}
				}

				////probably we need to mark its neighboring separation points 07/25/06
				MarkNeighborPoints(cur_e->sep.entry, face->index, scc_index, 0);

				//set the point as dead
				cur_e->sep_visit = 1;
			}

			else
			{
				if(cur_e->find_attp == 0)  //no attachment point
					continue;

				if(cur_e->att_visit == 1) //this is a dead point
					continue;

				if(!IsValidAtt(cur_e, face->index, scc_index)) //this is not a good attachment point
					continue;

				ResetEdgeIntersections(scc_index); //reset the intersection information before starting tracing

				//start tracing forward to locate attracting limit cycle if existed
				the_triangle = face->index;
				if(TraceForDetect(cur_e->attp.entry, the_triangle, 0, scc_index, chosen_edge, flag))
				{
					//we found a cycle, but we still need to return the intersection and the edge
					//for calculating the closed streamline for the limit cycle
					
					chosen_edge = g_chosenedge;
					if(GetClosedStreamline(chosen_edge, the_triangle, 0, scc_index))
					{
						//store the information for the limit cycle
						StoreCurrentCellCycle(cellcycle, num_celltriangle);  

						limitcycles[cur_limitcycle_index].node_index = -1; ////10/16/05

						StoreCurStreamline(trajectories[cur_traj_index], num_linesegs_curtraj[cur_traj_index]);     

						////store the type of the limit cycle
						limitcycles[cur_limitcycle_index].type = 1;

						////Initialize the connected list for the graph
						limitcycles[cur_limitcycle_index].connected_limitcycle = NULL;
						limitcycles[cur_limitcycle_index].num_connectedcycles = 0;
						limitcycles[cur_limitcycle_index].connected_saddle = NULL;
						limitcycles[cur_limitcycle_index].num_connectedsaddles = 0;
						limitcycles[cur_limitcycle_index].connected_l = 0;
						limitcycles[cur_limitcycle_index].connected_r = 0;
						BuildHandlerforLimitCycle(cur_limitcycle_index, 
							chosen_edge->intersections[1].entry, chosen_edge->verts[0]);
						cur_limitcycle_index++;
					}
				}

				////probably we need to mark its neighboring separation points 07/25/06
				MarkNeighborPoints(cur_e->sep.entry, face->index, scc_index, 1);

				//set the point as dead
				cur_e->att_visit = 1;
			}
		}
	}
}




/*
The following two judgement are very rough judgement, since two limit cycle can pass a same triangle
*/
bool IsPreviousCycle_SCC(int scc_index, int triangle, int &pre_limit_index)
{
	return false;
}

bool IsPreviousCycle_all(int triangle, int &pre_limit_index, int type)
{
	int i;
	
	for(i = 0; i < cur_limitcycle_index; i++)
	{
		//if(limitcycles[i].type == type) //can not have same type passing one triangle, can they? 07/23/06
		//	continue;

		if(limitcycles[i].type != type) 
			continue;

		if(IsRepeated(limitcycles[i].cellcycle, triangle, limitcycles[i].num_triangles))
		{
			pre_limit_index = i;  //return the limit cycle index
			return true;
		}
	}

	return false;
}


void ResetEdgeIntersections(int scc_index)
{
	int i, j;
	Face *face;
	Edge *cur_e;

	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		face = Object.flist[i];

		for(j = 0; j < 3; j++)
		{
			cur_e = face->edges[j];

			cur_e->num_intersections = 0;
			cur_e->intersections[0].set(0, 0);
			cur_e->intersections[1].set(0, 0);
			//cur_e->pre_length = 0.0;
		}
	}
}

/*
Assume we have already found out all the SCC, and get all the satisfied special points
*/

void SCCCycleDetect()
{
	int i, j;

	////Calculate the edge length here as the distance threshold
	icVector2 temp;
	temp.entry[0] = Object.vlist[Object.flist[0]->verts[0]]->x - 
		Object.vlist[Object.flist[0]->verts[1]]->x;
	temp.entry[1] = Object.vlist[Object.flist[0]->verts[0]]->y - 
		Object.vlist[Object.flist[0]->verts[1]]->y;
	ave_length = length(temp);
	/////////////////////////////////////////////////////////////

	//initialize the intersections on all edges
	Face *face;
	Edge *e;
	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];
		for(j = 0; j < 3; j++)
		{
			e = face->edges[j];
			e->num_intersections = 0;
			e->intersections[0].set(0, 0);
			e->intersections[1].set(0, 0);
		}
	}

	for(i = 0; i < scclist.num_sccs; i++)
	{
		if(scclist.scccomponents[i].valid == 0)
			continue;

		////Find the possible attracting limit cycles in this SCC
		if(scclist.scccomponents[i].num_attpts == 0)
		{
			////use regular method to trace
			RegularDetect(i, 1);   //1 means attracting
		}

		else
		{
			SpecialPtsBasedDetect(i, 1);
		}

		////Find the possible repelling limit cycles in this SCC
		if(scclist.scccomponents[i].num_seppts == 0)
		{
			////use regular method to trace
			RegularDetect(i, 0);   //0 means repelling
		}

		else
		{
			SpecialPtsBasedDetect(i, 0);
		}
	}
}

extern int num_sccnodes, num_sccedges;
double global_tau;

void DetectLimitCycle_SCC()
{
	//Find out all the SCC
	InitForSCC();
	//BuildConnectedGraph(10);
	BuildConnectedGraph(global_tau);
	
	FILE *fp = fopen("detectprocess.txt", "w");
	fprintf(fp, "finish constructing the directed graph!\n");
		fprintf(fp, "%d nodes and %d edges in the graph.\n", num_sccnodes, num_sccedges);
	fclose(fp);
	
	FindSCC();

	fp = fopen("detectprocess.txt", "a");
	fprintf(fp, "finish finding SCC!\n");
	fprintf(fp, "%d SCC's have been found!\n", scclist.num_sccs);
	fclose(fp);

	//Calculate the special points inside the SCC
	GetSpPtsForValidSCC();

	fp = fopen("detectprocess.txt", "a");
	fprintf(fp, "finish calculating special points!\n");
	//fprintf(fp, "%d points have been found!\n", scclist.num_sccs);
	fprintf(fp, "start detecting periodic orbits!\n");
	fclose(fp);

	//Detect limit cycle
	cur_limitcycle_index = 0;
	SCCCycleDetect();
	
	fp = fopen("detectprocess.txt", "a");
	fprintf(fp, "finish extracting periodic orbits!\n");
	fclose(fp);
}

/*****************************************************************/
/*
Routines for building the connections involving limit cycles
*/

int TraceInTriangleForConnection(double g[2], int &face_id, int type, int &flag)
{
	//similar to the regular tracing, except that you need to mark those special points too close to the cycle
	int i;
	double alpha[3];
	double cur_point[2], pre_point[2];
	double vert0[2];
	icVector2 VP, globalv;

	icVector2 dis;  //to calculate the distance for each line segment  08/10/06

	if(face_id < 0)
		return -1;

	Face *face = Object.flist[face_id];

	Face *pre_f = face;
	
	////initialize
	VP.entry[0] = g[0] - Object.vlist[face->verts[0]]->x;
	VP.entry[1] = g[1] - Object.vlist[face->verts[0]]->y;

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
			//if(ToNextPoint(pre_point, cur_point, face_id, alpha, type))
			//if(get_nextpt_2ndeuler(pre_point, cur_point, face_id, alpha, type))
			if(compute_next_pt(pre_point, cur_point, face_id, alpha, type))
			{
				////update the global point
				face = Object.flist[face_id];

				globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

				g[0] = vert0[0] + globalv.entry[0];
				g[1] = vert0[1] + globalv.entry[1];
			
				////08/10/06
				dis.entry[0] = cur_point[0] - pre_point[0];
				dis.entry[1] = cur_point[1] - pre_point[1];
				sum_flow_length += length(dis);

			}

			else{  ////the curve reaches a singularity
				flag = 3;
				globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

				g[0] = vert0[0] + globalv.entry[0];
				g[1] = vert0[1] + globalv.entry[1];

				return face_id;
			}

		}

		////3. if the point is out of current triangle
		else{
			double t[2] = {0.};

			int PassVertornot = 0;

			int presave_face = face_id;
			
			CrossVertex2(face_id, cur_point, pre_point, type, PassVertornot);
            
			////update the global point here
			if(PassVertornot > 0)
			{
				////we should not directly use the vertex as next point!!
				////we may move a little bit along the VF direction, but make sure it is still inside
				////the new triangle

				Vertex *PassedVert = Object.vlist[pre_f->verts[PassVertornot-1]];
				g[0] = PassedVert->x;
				g[1] = PassedVert->y;
           
				cur_point[0] = pre_f->xy[PassVertornot-1][0];
				cur_point[1] = pre_f->xy[PassVertornot-1][1];

				////08/10/06 get the length of the last line segment
				dis.entry[0] = cur_point[0] - pre_point[0];
				dis.entry[1] = cur_point[1] - pre_point[1];
				sum_flow_length += length(dis);
			}


			else{
				face_id = presave_face;  

				int which_edge = -1;

				CrossBoundary3(pre_point, cur_point, face_id, alpha, which_edge, t);

				PassEdge(face_id, which_edge);

				////transfer it to the global coordinates
				globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;
				
				g[0] = vert0[0] + globalv.entry[0];
				g[1] = vert0[1] + globalv.entry[1];
				
				////08/10/06 get the length of the last line segment
				dis.entry[0] = cur_point[0] - pre_point[0];
				dis.entry[1] = cur_point[1] - pre_point[1];
				sum_flow_length += length(dis);

			}

			return face_id;
		}
	}
	
	return face_id;
}


void TraceandBuildConnection(double s[2], int triangle, int singID, int type)
{
	int i;
	int flag = 0;
	double globalp[2];
	int pre_face, cur_face;

	pre_face = cur_face = triangle;

	globalp[0] = s[0];   globalp[1] = s[1];

	for(i = 0; i < 2*NUMTRACINGTRIANGLE; i++)
	{

		if(cur_face == -1)
		{
			return ;
		}

		pre_face = cur_face;

		////Here we need a sub routine that can return the intersection and the corresponding edge 
		////when it enters a new triangle

		cur_face = TraceInTriangleForConnection(globalp, cur_face, type, flag); 

		////Move the repeated cycle judgement here 07/25/06
		int pre_limit_index = -1;
		if(IsPreviousCycle_all(cur_face, pre_limit_index, 1-type))
		{
			//build connection here
			if(pre_limit_index < 0)
				continue;

			UpdateListInSingularity(singID, pre_limit_index);
			UpdateListInLimitCycle(pre_limit_index, singID);
			return;
		}

		if(flag == 3 || pre_face == cur_face) 
		{
			flag = 1;      //fail, reach singularity or other error, no cycle has been found on this path
			return;
		}
	}
}

/*
Trace and build the connection between nonsaddle singularities and limit cycles
*/
void TraceForConnection_nonsaddle(int singID)
{
	double sing_c[2], out[2];

	//first, find a point close to the center of the singularity
	sing_c[0] = singularities[singID].gcx;
	sing_c[1] = singularities[singID].gcy;
	FindAClosePoint(sing_c, singularities[singID].Triangle_ID, out);

	if(singularities[singID].type == SOURCE || singularities[singID].type == RFOCUS)
	{
		//we need to trace forward
		TraceandBuildConnection(out, singularities[singID].Triangle_ID, singID, 0);
	}
	else if(singularities[singID].type == SINK || singularities[singID].type == AFOCUS)
	{
		//we need to trace backward
		TraceandBuildConnection(out, singularities[singID].Triangle_ID, singID, 1);
	}
}

/*
Trace and build the connection between saddles and limit cycles
*/
void TraceForConnection_saddle(int saddle)
{
	double sing_center[2], newpos[2];
	icVector2 sep_vector;
	int start_triangle;

	sing_center[0] = singularities[saddle].gcx;
	sing_center[1] = singularities[saddle].gcy;

	//follow the four separatrices respectively

	//positive outgoing
	sep_vector = singularities[saddle].outgoing;
	start_triangle = singularities[saddle].Triangle_ID;
	FixSepBeginningPos(start_triangle, sep_vector, sing_center, newpos);
	TraceandBuildConnection(newpos, start_triangle, saddle, 0);

	//positive incoming
	sep_vector = singularities[saddle].incoming;
	start_triangle = singularities[saddle].Triangle_ID;
	FixSepBeginningPos(start_triangle, sep_vector, sing_center, newpos);
	TraceandBuildConnection(newpos, start_triangle, saddle, 1);

	//negative outgoing
	sep_vector = -singularities[saddle].outgoing;
	start_triangle = singularities[saddle].Triangle_ID;
	FixSepBeginningPos(start_triangle, sep_vector, sing_center, newpos);
	TraceandBuildConnection(newpos, start_triangle, saddle, 0);

	//negative incoming
	sep_vector = -singularities[saddle].incoming;
	start_triangle = singularities[saddle].Triangle_ID;
	FixSepBeginningPos(start_triangle, sep_vector, sing_center, newpos);
	TraceandBuildConnection(newpos, start_triangle, saddle, 1);
}



/*
Build the connections between singularities and limit cycles
*/
void ConnectionForSingCycle()
{
	int i;
	int limitcycle;

	for(i = 0; i < cur_singularity_index; i++)
	{
		if(singularities[i].type == SOURCE || singularities[i].type == RFOCUS
			||singularities[i].type == SINK || singularities[i].type == AFOCUS)
		{
			TraceForConnection_nonsaddle(i);
		}

		else if(singularities[i].type == SADDLE)
		{
			//we need to trace along its four separatrices
			TraceForConnection_saddle(i);
		}
	}
}

extern bool GetSingularyID(double x, double y, int &singular_id);

extern int GettheSide_loc(icVector2 orient, double basis[2], double p[2]);

void MarkTheLimitCycle(double globalp[2], int cur_face, int pre_limit_index)
{
	int i;
	Face *face = Object.flist[cur_face];
	double loc_p[2];
	icVector2 VP;
	icVector2 orient;
	double basis[2], head[2];
	VP.entry[0] = globalp[0] - Object.vlist[face->verts[0]]->x;
	VP.entry[1] = globalp[1] - Object.vlist[face->verts[0]]->y;
	loc_p[0] = dot(VP, face->LX);
	loc_p[1] = dot(VP, face->LY);

	int stop = 0;


	for(i = 0; i < limitcycles[pre_limit_index].num_linesegs; i++)
	{
		if(limitcycles[pre_limit_index].closed_streamline[i].Triangle_ID == cur_face && stop == 0)
		{
			//orient.entry[0] = limitcycles[pre_limit_index].closed_streamline[i].end[0]
			//	-limitcycles[pre_limit_index].closed_streamline[i].start[0];
			//orient.entry[1] = limitcycles[pre_limit_index].closed_streamline[i].end[1]
			//	-limitcycles[pre_limit_index].closed_streamline[i].start[1];

			basis[0] = limitcycles[pre_limit_index].closed_streamline[i].start[0];
			basis[1] = limitcycles[pre_limit_index].closed_streamline[i].start[1];

			stop = 1;

		}

		if(stop == 1 && limitcycles[pre_limit_index].closed_streamline[i].Triangle_ID != cur_face)
		{
			head[0] = limitcycles[pre_limit_index].closed_streamline[i-1].end[0];
			head[1] = limitcycles[pre_limit_index].closed_streamline[i-1].end[1];

			orient.entry[0] = head[0] - basis[0];
			orient.entry[1] = head[1] - basis[1];
			break;
		}
	}

	////
	normalize(orient);

	int side = GettheSide_loc(orient, basis, loc_p);

	if(side == 0)
		limitcycles[pre_limit_index].connected_l = 1;
	else if(side == 1)
		limitcycles[pre_limit_index].connected_r = 1;
}

/*
Trace and build the connection between two cycles
*/
void TraceandBuildConnection_cycle(double s[2], int triangle, int cycleID, int type)
{
	int i;
	int flag = 0;
	double globalp[2];
	int pre_face, cur_face;

	pre_face = cur_face = triangle;

	globalp[0] = s[0];   globalp[1] = s[1];

	for(i = 0; i < 2*NUMTRACINGTRIANGLE; i++)
	{

		if(cur_face == -1)
		{
			return ;
		}

		pre_face = cur_face;

		////Here we need a sub routine that can return the intersection and the corresponding edge 
		////when it enters a new triangle

		cur_face = TraceInTriangleForConnection(globalp, cur_face, type, flag); 

		////Move the repeated cycle judgement here 07/25/06
		int pre_limit_index = -1;
		if(IsPreviousCycle_all(cur_face, pre_limit_index, 1-type))
		{
			//build connection here
			if(pre_limit_index < 0)
				continue;

			UpdateCycleListInLimitCycle(cycleID, pre_limit_index);
			UpdateCycleListInLimitCycle(pre_limit_index, cycleID);

			//judge which side it is and mark this side of pre_limit_index "connected"
			MarkTheLimitCycle(globalp, cur_face, pre_limit_index);

			return;
		}

		if(flag == 3 || pre_face == cur_face) 
		{
			//if(flag == 3)
			//{
				//
				int singID = Object.flist[cur_face]->singularity_index;
				if(singID >= 0)
				{
					UpdateListInSingularity(singID, cycleID);
					UpdateListInLimitCycle(cycleID, singID);
					singularities[singID].connected = 1;
					limitcycles[cycleID].singularID = singID;  ////if we set the flag correctly, this will be fine!
				}
			//}

			//flag = 1;      //fail, reach singularity or other error, no cycle has been found on this path
			return;
		}
	}
}



//
void TraceForConnection_cycle(int cycle, int connected_side)
{
	//find out the edge contain the fixed point
	double fixed_p[2], v1[2], v2[2], out[2], alpha[3];
	int triangle = limitcycles[cycle].closed_streamline[0].Triangle_ID;
	if(triangle < 0) return;  /*07/06/2007*/
	Face *face = Object.flist[triangle];

	//get the local point to calculate the barycentric coordinates
	fixed_p[0] = limitcycles[cycle].closed_streamline[0].start[0];
	fixed_p[1] = limitcycles[cycle].closed_streamline[0].start[1];

	Get2DBarycentricFacters(triangle, fixed_p[0], fixed_p[1], alpha);

	//get the global coordinates of the fixed point
	fixed_p[0] = limitcycles[cycle].closed_streamline[0].gstart[0];
	fixed_p[1] = limitcycles[cycle].closed_streamline[0].gstart[1];

	//Do we need to consider the case that fixed point locates at a vertex ?
	if(abs(alpha[0])<1e-8) //fixed point is on the edge v1v2
	{
		v1[0] = Object.vlist[face->verts[1]]->x;
		v1[1] = Object.vlist[face->verts[1]]->y;
		
		v2[0] = Object.vlist[face->verts[2]]->x;
		v2[1] = Object.vlist[face->verts[2]]->y;
	}

	else if(abs(alpha[1])<1e-8) //fixed point is on the edge v2v0
	{
		v1[0] = Object.vlist[face->verts[2]]->x;
		v1[1] = Object.vlist[face->verts[2]]->y;
		
		v2[0] = Object.vlist[face->verts[0]]->x;
		v2[1] = Object.vlist[face->verts[0]]->y;
	}

	else{  //fixed point iss on the edge v0v1
		v1[0] = Object.vlist[face->verts[0]]->x;
		v1[1] = Object.vlist[face->verts[0]]->y;
		
		v2[0] = Object.vlist[face->verts[1]]->x;
		v2[1] = Object.vlist[face->verts[1]]->y;
	}

	//get the first point
	out[0] = (v1[0]+fixed_p[0])/2.;
	out[1] = (v1[1]+fixed_p[1])/2.;

	////Judge whether this side has been marked "connected" or not
	//first, we need to tranfer the point into local coordinates
	if(connected_side != -1)
	{
		icVector2 VP, orient;
		double p_loc[2];
		VP.entry[0] = out[0] - Object.vlist[face->verts[0]]->x;
		VP.entry[1] = out[1] - Object.vlist[face->verts[0]]->y;
		p_loc[0] = dot(VP, face->LX);
		p_loc[1] = dot(VP, face->LY);

		orient.entry[0] = limitcycles[cycle].closed_streamline[0].end[0]-
			limitcycles[cycle].closed_streamline[0].start[0];
		orient.entry[1] = limitcycles[cycle].closed_streamline[0].end[1]-
			limitcycles[cycle].closed_streamline[0].start[1];

		normalize(orient);
		int which_side = GettheSide_loc(orient, limitcycles[cycle].closed_streamline[0].start, p_loc);
		
		if(which_side != connected_side)
		{
			sum_flow_length = 0.0;                 //initialize the flow length
			TraceandBuildConnection_cycle(out, triangle, cycle, limitcycles[cycle].type);
			return;
		}

		else
		{
			//get the second point
			out[0] = (v2[0]+fixed_p[0])/2.;
			out[1] = (v2[1]+fixed_p[1])/2.;

			sum_flow_length = 0.0;                 //initialize the flow length
			TraceandBuildConnection_cycle(out, triangle, cycle, limitcycles[cycle].type);
		}
	}

	else
	{
			sum_flow_length = 0.0;                 //initialize the flow length
		TraceandBuildConnection_cycle(out, triangle, cycle, limitcycles[cycle].type);
		
		//get the second point
		out[0] = (v2[0]+fixed_p[0])/2.;
		out[1] = (v2[1]+fixed_p[1])/2.;

			sum_flow_length = 0.0;                 //initialize the flow length
		TraceandBuildConnection_cycle(out, triangle, cycle, limitcycles[cycle].type);
	}

}

/*
Build the connections between limit cycle pairs
*/

void ConnectionForCyclePair()
{
	int i;

	for(i = 0; i < cur_limitcycle_index; i++)
	{
		if(limitcycles[i].connected_l == 1 && limitcycles[i].connected_r == 1)
			continue;

		else if(limitcycles[i].connected_l == 1)
			TraceForConnection_cycle(i, 0);

		else if(limitcycles[i].connected_r == 1)
			TraceForConnection_cycle(i, 1);

		else
			TraceForConnection_cycle(i, -1);
	}
}



void BuildConnectionForCycle()
{
	ConnectionForCyclePair();
	//ConnectionForSingCycle();
}
