////LocalTracing.cpp

#include "stdafx.h"

#include "LocalTracing.h"

#include "VFDataStructure.h"

#include "Numerical.h"

#include "GL/glut.h"

#include "GlView.h"

#include "ClosedStreamlineTracing.h"

////Global variables
int globalface;

extern int MaxNumTrajectories;
extern int MaxNumLinesegsPerTraj;
extern int MaxNumSeparatrices;                   //Maximum number of group of separatrices

extern Polygon3D Object;
//extern Trajectory *trajectories2;
extern LineSeg **trajectories;                 //trajectories' list
extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory
extern Singularities *singularities;           //being captured singularites' list
extern int cur_singularity_index;
extern Separatrices *separatrices;             //array for group of separatrices
extern int cur_separatrices_index;

extern int pre_cur_traj_index;

extern void getVector(double, double, icVector2&, double&);
extern void OneBackwardRungeKutta(double, double, double &, double &, int);

extern double triangle_approx_size; 


extern double avg_edge_length; //record the average edge length 1/8/06


//////Variable for limit cycle tangent curve detection
extern double tang_pre[3], tang_cur[3];  ////They are 3D global coordinates
extern icVector2 intersect_edge;     // record the direction of the current intersected edge

extern double ave_length;

////Important global variable for adaptive stepsize of integration
DP htry = 1.;

////////Testing variables
int cur_sing;
double problemx, problemy;

/* For \tau map here 02/27/07 */
extern double g_dt;

extern double hstep;  /*for tau map*/



bool StoreToGlobalList(CurvePoints *temp, int num);

/*------------------07/05/2007------------------*/
/*Important global variable for spatial tau!!!
*/
icVector2 g_cur_vec;
double g_vec_mag;


extern bool UseNormalizedVF;

bool RK2_or_RK4 = true;
int Mag_Scheme = 0;
bool test_nonnormalized_tracing = false;

/////////////////////////////////////////////////////////////////////////
/*********************************************************************************
Routines for local tracing
*********************************************************************************/

////detect the loop in the trajectory to improve the speed
bool LoopDetect(int *triangles, int num, int oneTriangle)
{
	int i;

	for(i = 0; i < num; i++)
	{
		if(triangles[i] == oneTriangle)
			return true;
	}

	return false;
}


 //Calculate the trajectory using local frame
void CalLocalTracing(int face_id, double x, double y, int type)
{
	int i;
	int flag = 0;
	double globalp[2];
	int pre_face, cur_face;

	int loop_flag = 0;  ////for detect a closed loop from the beginning triangle
	int *triangles = (int*)malloc(sizeof(int) * NUMTRACINGTRIANGLE);
	int num_triangles = 0;

	pre_face = cur_face = face_id;

	globalp[0] = x;   globalp[1] = y;

	////If current trajectories list is not enough to store all the trajectories
	////extend it!
	if(cur_traj_index >= MaxNumTrajectories)
	{
		MaxNumTrajectories += 100;

		num_linesegs_curtraj = (int*)realloc(num_linesegs_curtraj, sizeof(int)*MaxNumTrajectories);
		trajectories = (LineSeg **)realloc(trajectories, sizeof(LineSeg *) * MaxNumTrajectories);

		////extend the old trajectories
		for(i = 0; i < MaxNumTrajectories-100; i++)
			trajectories[i] = (LineSeg *)realloc(trajectories[i], sizeof(LineSeg ) * MaxNumLinesegsPerTraj);

		////allocate the new trajectories
		for( ; i < MaxNumTrajectories; i++)
		{
			trajectories[i] = (LineSeg *)malloc(sizeof(LineSeg) * MaxNumLinesegsPerTraj);
			num_linesegs_curtraj[i] = 0;
		}

		////You need to extend the new trajectories variable here 08/29/05
	}

	num_linesegs_curtraj[cur_traj_index] = 0;


	////We may need to adjust the beginning point to make sure it fall into the triangle cur_face,
	////Modified at 1/5/06 for limit cycle detection
//LL:	int test_triangle = TriangleDetect(globalp[0], globalp[1]);
//	if(test_triangle != cur_face)
//	{
//		globalp[0] += 1e-9;
//		//OneBackwardRungeKutta(globalp[0], globalp[1], globalp[0], globalp[1], type);
//	    cur_face = test_triangle;
//		goto LL;
//	}
		
	//OneBackwardRungeKutta(globalp[0], globalp[1], globalp[0], globalp[1], type);
	//cur_face = TriangleDetect(globalp[0], globalp[1]);

	////if the beginning point is quite close to a vertex of the triangle, move it a little bit away from the vertex


	triangles[num_triangles] = cur_face;
	num_triangles++;

	////Testing codes here
	problemx = globalp[0];
	problemy = globalp[1];

	/*02/27/07*/
	g_dt = 0;


	for(i = 0; i < 10*NUMTRACINGTRIANGLE; i++)
	{

		if(cur_face == -1)
		{
			//MessageBox(NULL, "Not a legal triangle!", "Error", MB_OK);
			free(triangles);
			return;
		}

		////Set the contain_separatrix flag 1/3/06
		Object.flist[cur_face]->contain_separatrix = 1;

		pre_face = cur_face;
		cur_face = TraceInATriangle(cur_face, globalp, type, flag); ////0 means always forward tracing here

		if(flag == 3 || flag == 4 || pre_face == cur_face || loop_flag == 1) 
		{
			free(triangles);
			return;
		}

		////Not accurate to use triangle loop detection only
		////after we find a triangle loop, we need to test whether the curve constitutues a "closed" orbit!!!
		////we can use a small neighbor hood to measure the distance !

		//Do not use Loop judgement here 12/09/2005
		//if(LoopDetect(triangles, num_triangles, cur_face))
		//	loop_flag = 1;
		//else{
		//	triangles[num_triangles] = cur_face;
		//	num_triangles++;
		//}
	}

	////saved for limit cycle detection
	pre_cur_traj_index = cur_traj_index;

	free(triangles);

}

////calculate the trajectory in a single triangle
int TraceInATriangle(int &face_id, double globalp[2], int type, int &flag)
{
	int i;
	double alpha[3];
	double cur_point[2], pre_point[2];
	double vert0[2];
	icVector2 VP, globalv;

	if(face_id < 0)
		return -1;

	Face *face = Object.flist[face_id];

	Face *pre_f = face;
	
	////Temporary curve point array

	CurvePoints *temp_point_list = (CurvePoints*) malloc(sizeof(CurvePoints) * TRACESTEPS);
	int NumPoints = 0;

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
			////store the point into the temp curve points list

			temp_point_list[NumPoints].gpx = globalp[0];
			temp_point_list[NumPoints].gpy = globalp[1];
			temp_point_list[NumPoints].lpx = cur_point[0];
			temp_point_list[NumPoints].lpy = cur_point[1];
			temp_point_list[NumPoints].triangleid = face->index;  
			NumPoints++;

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

				globalp[0] = vert0[0] + globalv.entry[0];
				globalp[1] = vert0[1] + globalv.entry[1];

			}

			else{  ////the curve reach a singularity
				flag = 3;

				////Store the record into global line segment array
                
				if(!StoreToGlobalList(temp_point_list, NumPoints))
				{
					////Not enough memory
					flag = 4;
					free(temp_point_list);
					return face_id;
				}

				free(temp_point_list);

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
				////we should not directly use the vertex as next point!!
				////we may move a little bit along the VF direction, but make sure it is still inside
				////the new triangle

				//Vertex *PassedVert = Object.vlist[pre_f->verts[PassVertornot-1]];
				//globalp[0] = PassedVert->x;
				//globalp[1] = PassedVert->y;
				//
				//cur_point[0] = pre_f->xy[PassVertornot-1][0];
				//cur_point[1] = pre_f->xy[PassVertornot-1][1];

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

			}

			else{
				////transfer it to the global coordinates
				globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

				globalp[0] = vert0[0] + globalv.entry[0];
				globalp[1] = vert0[1] + globalv.entry[1];
			}

			////Add the intersection point to the temporary points' list
			temp_point_list[NumPoints].gpx = globalp[0];
			temp_point_list[NumPoints].gpy = globalp[1];
			temp_point_list[NumPoints].lpx = cur_point[0];
			temp_point_list[NumPoints].lpy = cur_point[1];
			temp_point_list[NumPoints].triangleid = face->index;  ////cause problem 05/25/05
			NumPoints++;

			if(NumPoints > 1){
 				////Store the record into global line segment array
               if(!StoreToGlobalList(temp_point_list, NumPoints))
			   {   ////Not enough memory
				   flag = 4;
				   free(temp_point_list);
				   return face_id;
			   }
			}

			free(temp_point_list);
			return face_id;
		}

	}

    StoreToGlobalList(temp_point_list, NumPoints);
	free(temp_point_list);
	return face_id;
}


/*****************************************************************************
Only trace in a triangle, not store the tracing points
For limit cycle detection, we may need to store two recently calculated points
Modified again on 02/20/07
*****************************************************************************/
int TraceInATriangle2(int &face_id, double globalp[2], int type, int &flag)
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

			//if(ToNextPoint(pre_point, cur_point, face_id, alpha, type))
			//if(get_nextpt_euler(pre_point, cur_point, face_id, alpha, type))
			//if(get_nextpt_2ndeuler(pre_point, cur_point, face_id, alpha, type))
			if(compute_next_pt(pre_point, cur_point, face_id, alpha, type))
			{
				////update the global point
				face = Object.flist[face_id];

				globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

				globalp[0] = vert0[0] + globalv.entry[0];
				globalp[1] = vert0[1] + globalv.entry[1];

				/*calculate the integration length 03/20/07*/
				line_length.entry[0] = cur_point[0]-pre_point[0];
				line_length.entry[1] = cur_point[1]-pre_point[1];

				g_dt += length(line_length);
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
				
				g_dt += length(line_length);

			}

			else{
				////transfer it to the global coordinates
				globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

				globalp[0] = vert0[0] + globalv.entry[0];
				globalp[1] = vert0[1] + globalv.entry[1];
				
				/*calculate the length of last line segment*/
				line_length.entry[0] = cur_point[0]-pre_point[0];
				line_length.entry[1] = cur_point[1]-pre_point[1];

				g_dt += length(line_length);
			}

			return face_id;
		}

	}

	return face_id;
}

/*
Trace particle for one step for particle advection
*/
int TraceParticleforOneStep(double &x, double &y, int &face_id, int type)
{
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
	VP.entry[0] = x - Object.vlist[face->verts[0]]->x;
	VP.entry[1] = y - Object.vlist[face->verts[0]]->y;

	pre_point[0] = cur_point[0] = dot(VP, face->LX);
	pre_point[1] = cur_point[1] = dot(VP, face->LY);

	vert0[0] = Object.vlist[face->verts[0]]->x;   ////for update the global point
	vert0[1] = Object.vlist[face->verts[0]]->y;

	Get2DBarycentricFacters(face_id, cur_point[0], cur_point[1], alpha);

	Euler_ToNextPoint(pre_point, cur_point, face_id, alpha, type); //call Euler method to move forward

	//// if current point is inside current triangle
	if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-8) && alpha[0] <= 1 
		&& (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-8) && alpha[1] <= 1 
		&& (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-8) && alpha[2] <= 1)
	{
		globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;
		x = vert0[0] + globalv.entry[0];
		y = vert0[1] + globalv.entry[1];
	}

	////3. if the point is out of current triangle
	else{
		double t[2] = {0.};

		int PassVertornot = 0;
        
		GetNextTriangle(face_id, pre_point, cur_point, t, type, PassVertornot, alpha);

		if(face_id < 0)
		{
			return -1;
		}

		////update the globalpoint here
		if(PassVertornot > 0)
		{     
			x = Object.vlist[Object.flist[face_id]->verts[PassVertornot-1]]->x;
			y = Object.vlist[Object.flist[face_id]->verts[PassVertornot-1]]->y;
		}

		else{
			//// transfer it to the global coordinates
			globalv = cur_point[0] * pre_f->LX + cur_point[1] * pre_f->LY;

			x = vert0[0] + globalv.entry[0];
			y = vert0[1] + globalv.entry[1];
		}
	}
	return face_id;

}


/**********************************************************************
Get the barycentric coordinates in local frame
**********************************************************************/
void Get2DBarycentricFacters(int face_id, double a, double b, double alpha[3])
{
	////first, we need to calculate the local coordinates of the three vertices
	//// Maybe we should put them into the initialtracing routine

	Face *face = Object.flist[face_id];

	alpha[2] = b/face->xy[2][1];

    if(fabs(alpha[2]) < 1e-10) alpha[2] = 0.;

	alpha[1] = (a - alpha[2]*face->xy[2][0])/ face->xy[1][0];

	if(fabs(alpha[1]) < 1e-10) alpha[1] = 0.;

	alpha[0] = 1 - alpha[1] - alpha[2];
}


/***************************************************************
Get the vector at any point inside the triangle local frame
For 2D planar vector field, we can treat each triangle contains
only piecewise linear vector field
So we can simply interpolate the vectors on three vertices
using the barycentric coordinates
***************************************************************/
icVector2 GetVectorAtPoints(int face_id, double alpha[3])
{
	Face *face = Object.flist[face_id];
	icVector2 result;

	result = alpha[0]*face->direct_vec[0] + alpha[1]*face->direct_vec[1] + alpha[2]*face->direct_vec[2];

	return result;
}


/*
Try the second order Euler method 1/9/06
*/
void SecondOrderEulerStep(double first[2], double second[2], icVector2 VecAtPoint, int &face_id, double alpha[3], int type)
{
	////Using first order Euler method to get the next point
	double temp[2] = {0.};
	
	if(type == 0)
	{
		temp[0] = first[0] + VecAtPoint.entry[0];
		temp[1] = first[1] + VecAtPoint.entry[1];
	}
	else
	{
		temp[0] = first[0] - VecAtPoint.entry[0];
		temp[1] = first[1] - VecAtPoint.entry[1];
	}

	////get the vector at next point
	Get2DBarycentricFacters(face_id, temp[0], temp[1], alpha);

	icVector2 VecAtPoint2 = GetVectorAtPoints(face_id, alpha);

	if(type == 0)
	{
		second[0] = first[0] + (VecAtPoint.entry[0] +  VecAtPoint2.entry[0])/2;
		second[1] = first[1] + (VecAtPoint.entry[1] +  VecAtPoint2.entry[1])/2;
	}
	else
	{
		second[0] = first[0] - (VecAtPoint.entry[0] +  VecAtPoint2.entry[0])/2;
		second[1] = first[1] - (VecAtPoint.entry[1] +  VecAtPoint2.entry[1])/2;
	}
}


/*
Here we use first order Euler to get next point
*/
double euler_stepsize = 0.15;

bool get_nextpt_euler(double first[2], double second[2], int &face_id, double alpha[3], int type)
{
	icVector2 VecAtPoint = GetVectorAtPoints(face_id, alpha);

	if(length(VecAtPoint) < 1e-12) return false;

	if(type == 0)
	{
		second[0] = first[0] + euler_stepsize*VecAtPoint.entry[0];
		second[1] = first[1] + euler_stepsize*VecAtPoint.entry[1];
	}

	else
	{
		second[0] = first[0] - euler_stepsize*VecAtPoint.entry[0];
		second[1] = first[1] - euler_stepsize*VecAtPoint.entry[1];
	}

	hstep = euler_stepsize*length(VecAtPoint);

	g_dt += hstep;
}

/*
second order
*/
double scaling_cons = 8;
double RK2_step;

bool get_nextpt_2ndeuler_2(double first[2], double second[2], int &face_id, double alpha[3], int type)
{
	////Using first order Euler method to get the next point
	icVector2 V1 = GetVectorAtPoints(face_id, alpha);
	if(length(V1) < 1e-14) return false;

	double temp[2] = {0.};
	icVector2 t_v;

	double shortest = get_shortestedge_tri(face_id);
	double scalar = 0.05;
	//double scalar = 10;

	if(Mag_Scheme == 0)
	{
		RK2_step = scalar*shortest*0.5;
	}
	else if(Mag_Scheme == 1)
	{
		RK2_step = scalar*shortest;
	}
	else
	{
		RK2_step = scalar*shortest*2;
	}

	t_v = V1;
	normalize(t_v);

	if(type == 0)
	{
		temp[0] = first[0] + RK2_step*t_v.entry[0];
		temp[1] = first[1] + RK2_step*t_v.entry[1];
	}
	else
	{
		temp[0] = first[0] - RK2_step*t_v.entry[0];
		temp[1] = first[1] - RK2_step*t_v.entry[1];
	}

	/*compute K2*/
	Get2DBarycentricFacters(face_id, temp[0], temp[1], alpha);
	icVector2 V2 = GetVectorAtPoints(face_id, alpha);

	icVector2 total_v = 0.5*(V1+V2);
	t_v = total_v;
	normalize(t_v);

	if(type == 0)
	{
		second[0] = first[0] + RK2_step*t_v.entry[0];
		second[1] = first[1] + RK2_step*t_v.entry[1];
	}

	else
	{
		second[0] = first[0] - RK2_step*t_v.entry[0];
		second[1] = first[1] - RK2_step*t_v.entry[1];
	}

	/*we need to accumulate the time rather than the distance */
	/*calculate the t according to the step and the magnitude of the vector,
	but we need to first project the vector onto the line segment to get the real length*/
	icVector2 line_v;
	line_v.entry[0] = second[0]-first[0];
	line_v.entry[1] = second[1]-first[1];
	normalize(line_v);
	double proj_len = dot(line_v, total_v); /*here we use the dot product in the local frame*/

	/*for temporal tau*/

	g_cur_vec = total_v;
	g_vec_mag = fabs(proj_len);
}



bool get_nextpt_2ndeuler(double first[2], double second[2], int &face_id, double alpha[3], int type)
{
	////Using first order Euler method to get the next point
	icVector2 VecAtPoint = GetVectorAtPoints(face_id, alpha);
	if(length(VecAtPoint) < 1e-20) return false;

	double temp[2] = {0.};

	icVector2 before_VecAtPt = VecAtPoint;

	double shortest = get_shortestedge_tri(face_id);


	if(type == 0)
	{
		//if(UseNormalizedVF){
		//temp[0] = first[0] + euler_stepsize*VecAtPoint.entry[0];
		//temp[1] = first[1] + euler_stepsize*VecAtPoint.entry[1];
		//}
		////temp[0] = first[0] + 0.05*euler_stepsize*VecAtPoint.entry[0];/*07/25/07*/
		////temp[1] = first[1] + 0.05*euler_stepsize*VecAtPoint.entry[1];
		//else{
		//temp[0] = first[0] + 10*euler_stepsize*VecAtPoint.entry[0];/*07/25/07*/
		//temp[1] = first[1] + 10*euler_stepsize*VecAtPoint.entry[1];
		//}

		temp[0] = first[0] + euler_stepsize*VecAtPoint.entry[0];
		temp[1] = first[1] + euler_stepsize*VecAtPoint.entry[1];
	}
	else
	{
		//if(UseNormalizedVF){
		//temp[0] = first[0] - euler_stepsize*VecAtPoint.entry[0];
		//temp[1] = first[1] - euler_stepsize*VecAtPoint.entry[1];
		//}
		//
		////temp[0] = first[0] - 0.05*euler_stepsize*VecAtPoint.entry[0]; /*07/25/07*/
		////temp[1] = first[1] - 0.05*euler_stepsize*VecAtPoint.entry[1];
		//else{
		//temp[0] = first[0] - 10*euler_stepsize*VecAtPoint.entry[0]; /*07/25/07*/
		//temp[1] = first[1] - 10*euler_stepsize*VecAtPoint.entry[1];
		//}

		temp[0] = first[0] - euler_stepsize*VecAtPoint.entry[0];
		temp[1] = first[1] - euler_stepsize*VecAtPoint.entry[1];
	}

	////get the vector at next point
	Get2DBarycentricFacters(face_id, temp[0], temp[1], alpha);

	icVector2 VecAtPoint2 = GetVectorAtPoints(face_id, alpha);
	

	/*
	we can adjust the step size here based on the difference between VecAtPoint and VecAtPoint2
	*/
	
	icVector2 total_v;

	/*The following is just mid-point average, it is not always true in every step!*/
	total_v.entry[0] = (VecAtPoint.entry[0] +  VecAtPoint2.entry[0])/2;
	total_v.entry[1] = (VecAtPoint.entry[1] +  VecAtPoint2.entry[1])/2;
	
	//normalize(total_v);/*07/30/07*/

	if(type == 0)
	{
		second[0] = first[0] + euler_stepsize*total_v.entry[0];
		second[1] = first[1] + euler_stepsize*total_v.entry[1];
	}

	else
	{
		second[0] = first[0] - euler_stepsize*total_v.entry[0];
		second[1] = first[1] - euler_stepsize*total_v.entry[1];
	}

	//hstep = euler_stepsize*length(total_v);
	
	/*********07/05/2007**********/


	/*we need to accumulate the time rather than the distance */
	/*calculate the t according to the step and the magnitude of the vector,
	but we need to first project the vector onto the line segment to get the real length*/

	{
	icVector2 line_v;
	line_v.entry[0] = second[0]-first[0];
	line_v.entry[1] = second[1]-first[1];
	normalize(line_v);
	double proj_len = dot(line_v, total_v); /*here we use the dot product in the local frame*/

	/*for temporal tau*/

	g_cur_vec = total_v;
	g_vec_mag = fabs(proj_len);
	}
}


/*
For streamline tracing.
Here we allow user select the tracing scheme he/she wants
*/

bool compute_next_pt(double first[2], double second[2], int &face_id, double alpha[3], int type)
{
	if(RK2_or_RK4 == false)
	{
		if(!test_nonnormalized_tracing)
		{
			if(get_nextpt_2ndeuler(first, second, face_id, alpha, type))
				return true;
			else
				return false;
		}
		else
		{
			if(get_nextpt_2ndeuler_2(first, second, face_id, alpha, type))
				return true;
			else
				return false;
		}
	}
	else{
		if(!test_nonnormalized_tracing)
		{
			if(ToNextPoint(first, second, face_id, alpha, type))
			//if(get_nextpt_RK4_adp(first, second, face_id, alpha, type))
			//if(get_nextpt_RK4(first, second, face_id, alpha, type))
				return true;
			else
				return false;
		}
		else
		{
			if(get_nextpt_RK4(first, second, face_id, alpha, type))
				return true;
			else
				return false;
		}
	}
}


double RK4_step;
/*
Implement the fourth order Runge-Kutta integration without using adaptive step
*/
bool get_nextpt_RK4(double first[2], double second[2], int &face_id, double alpha[3], int type)
{
	////Using 4th Runge-Kutta method to get the next point
	icVector2 t_v;
	double temp[2] = {0.};

	/*compute K1*/
	icVector2 V1 = GetVectorAtPoints(face_id, alpha);

	if(length(V1) < 1e-14) return false;

	double shortest = get_shortestedge_tri(face_id);
	double scalar = 0.06;

	if(Mag_Scheme == 0)
	{
		RK4_step = scalar*shortest*0.5;
	}
	else if(Mag_Scheme == 1)
	{
		RK4_step = scalar*shortest;
	}
	else
	{
		RK4_step = scalar*shortest*2;
	}

	t_v = V1;
	normalize(t_v);

	if(type == 0)
	{
		temp[0] = first[0] + RK4_step/2*t_v.entry[0];
		temp[1] = first[1] + RK4_step/2*t_v.entry[1];
	}
	else
	{
		temp[0] = first[0] - RK4_step/2*t_v.entry[0];
		temp[1] = first[1] - RK4_step/2*t_v.entry[1];
	}

	/*compute K2*/
	Get2DBarycentricFacters(face_id, temp[0], temp[1], alpha);
	icVector2 V2 = GetVectorAtPoints(face_id, alpha);
	t_v = V2;
	normalize(t_v);
	
	if(type == 0)
	{
		temp[0] = first[0] + RK4_step/2*t_v.entry[0];
		temp[1] = first[1] + RK4_step/2*t_v.entry[1];
	}
	else
	{
		temp[0] = first[0] - RK4_step/2*t_v.entry[0];
		temp[1] = first[1] - RK4_step/2*t_v.entry[1];
	}

	/*compute K3*/
	Get2DBarycentricFacters(face_id, temp[0], temp[1], alpha);
	icVector2 V3 = GetVectorAtPoints(face_id, alpha);
	t_v = V3;
	normalize(t_v);
	
	if(type == 0)
	{
		temp[0] = first[0] + RK4_step*t_v.entry[0];
		temp[1] = first[1] + RK4_step*t_v.entry[1];
	}
	else
	{
		temp[0] = first[0] - RK4_step*t_v.entry[0];
		temp[1] = first[1] - RK4_step*t_v.entry[1];
	}
	
	/*compute K4*/
	Get2DBarycentricFacters(face_id, temp[0], temp[1], alpha);
	icVector2 V4 = GetVectorAtPoints(face_id, alpha);

	icVector2 total_v = 1./6.*(V1+2*V2+2*V3+V4);
	t_v = total_v;
	normalize(t_v);

	if(type == 0)
	{
		second[0] = first[0] + RK4_step*t_v.entry[0];
		second[1] = first[1] + RK4_step*t_v.entry[1];
	}

	else
	{
		second[0] = first[0] - RK4_step*t_v.entry[0];
		second[1] = first[1] - RK4_step*t_v.entry[1];
	}

	//hstep = euler_stepsize*length(total_v);
	
	/*********07/05/2007**********/
	//if(time_or_patial == 0)  /*use spatial tau*/
		//g_dt += hstep;


	/*we need to accumulate the time rather than the distance */
	/*calculate the t according to the step and the magnitude of the vector,
	but we need to first project the vector onto the line segment to get the real length*/
	icVector2 line_v;
	line_v.entry[0] = second[0]-first[0];
	line_v.entry[1] = second[1]-first[1];
	normalize(line_v);
	double proj_len = dot(line_v, total_v); /*here we use the dot product in the local frame*/

	/*for temporal tau*/
	g_cur_vec = total_v;
	g_vec_mag = fabs(proj_len);
}


/*****************************************************************
Runge Kutta to get next point for local tracing
For simplicity, here first and second still use global coordinates
Otherwise, we may need a global variable to tell 
Runge Kutta integrator which triangle it is in now
*****************************************************************/

bool ToNextPoint(double first[2], double second[2], int &face_id, double alpha[3], int type)
{

	icVector2 VecAtPoint = GetVectorAtPoints(face_id, alpha);

	/*Based on Eugene's suggestion, we normalize the vector here 07/19/07*/
	//normalize(VecAtPoint);
	//VecAtPoint = DistanceThreshold*VecAtPoint;

	if(length(VecAtPoint) < 1e-20 ) return false; /*for saddle-saddle connection example only 07/11/07*/
	
	////calling Runge Kutta to get the next point
    const int N=2;
    int i,j;
    DP eps,hdid,hnext, t = 1.0;

    Vec_DP by(N),dydx(N),dysav(N),ysav(N),yscal(N);

	ysav[0] = first[0];
	ysav[1] = first[1];

	if(type == 0)
        localderive(t, ysav, dysav);
	else
		localinverse_derive(t, ysav, dysav);

    for (i=0;i<N;i++) yscal[i]=1;

	/*set the htry to be the smallest edge length  07/30/07*/

	/*--------------------------------------------------------------------*/
	////calling adaptive stepsize runge-kutta method here to get next step
	for (i=0;i<15;i++) {
		eps=exp(-DP(i+1));
		t = 1.0;
		for (j=0;j<N;j++) {
			by[j]=ysav[j];
			dydx[j]=dysav[j];
		}
	
	   if(type == 0)
		   rkqs(by,dydx,t,htry,eps,yscal,hdid,hnext,localderive);
	   else
		   rkqs(by,dydx,t,htry,eps,yscal,hdid,hnext,localinverse_derive);

	}

	////Set the new stepsize, note that this stepsize will affect the limit cycle detection
	htry = hnext;
	if(hnext >= 2.)
		htry = 2.;


	/////store some important information
	second[0] = by[0];
	second[1] = by[1];
	/*--------------------------------------------------------------------*/
	
	icVector2 line_v;
	line_v.entry[0] = second[0]-first[0];
	line_v.entry[1] = second[1]-first[1];
	normalize(line_v);
	double proj_len = dot(line_v, VecAtPoint); /*here we use the dot product in the local frame*/


	/*for temporal tau*/
	g_cur_vec = VecAtPoint;

	g_vec_mag = fabs(proj_len);
	//g_dt += hdid/g_vec_mag;

	return true;
}


double RK4_const = 10;

bool get_nextpt_RK4_adp(double first[2], double second[2], int &face_id, double alpha[3], int type)
{

	icVector2 VecAtPoint = GetVectorAtPoints(face_id, alpha);

	/*Based on Eugene's suggestion, we normalize the vector here 07/19/07*/
	//normalize(VecAtPoint);
	//VecAtPoint = DistanceThreshold*VecAtPoint;

	if(length(VecAtPoint) < 1e-15 ) return false; /*for saddle-saddle connection example only 07/11/07*/
	
	////calling Runge Kutta to get the next point
    const int N=2;
    int i,j;
    DP eps,hdid,hnext, t = 1.0;

    Vec_DP by(N),dydx(N),dysav(N),ysav(N),yscal(N);

	ysav[0] = first[0];
	ysav[1] = first[1];

	if(type == 0)
        localderive(t, ysav, dysav);
	else
		localinverse_derive(t, ysav, dysav);

    for (i=0;i<N;i++) yscal[i]=1;

	/*set the htry to be the smallest edge length  07/30/07*/

	double shortest = get_shortestedge_tri(face_id);

	if(Mag_Scheme == 0) /*use half length of the shortest edge*/
	{
		htry = shortest/2.*RK4_const;
	}

	else if(Mag_Scheme == 1) /*use full length of the shortest edge*/
	{
		htry = shortest*RK4_const;
	}

	else
	{
		htry = shortest*2*RK4_const;
	}
	
	if(htry >= 2.)
		htry = 2.;

	/*--------------------------------------------------------------------*/
	////calling adaptive stepsize runge-kutta method here to get next step
	for (i=0;i<10;i++) {
		eps=exp(-DP(i+1));
		t = 1.0;
		for (j=0;j<N;j++) {
			by[j]=ysav[j];
			dydx[j]=dysav[j];
		}
	
	   if(type == 0)
		   rkqs(by,dydx,t,htry,eps,yscal,hdid,hnext,localderive);
	   else
		   rkqs(by,dydx,t,htry,eps,yscal,hdid,hnext,localinverse_derive);

	}

	////Set the new stepsize, note that this stepsize will affect the limit cycle detection
	htry = hnext;
	if(hnext >= 2.)
		htry = 2.;


	/////store some important information
	second[0] = by[0];
	second[1] = by[1];
	/*--------------------------------------------------------------------*/
	
	icVector2 line_v;
	line_v.entry[0] = second[0]-first[0];
	line_v.entry[1] = second[1]-first[1];
	normalize(line_v);
	double proj_len = dot(line_v, VecAtPoint); /*here we use the dot product in the local frame*/


	/*for temporal tau*/
	g_cur_vec = VecAtPoint;

	g_vec_mag = fabs(proj_len);
	//g_dt += hdid/g_vec_mag;

	return true;
}


bool ToNextPoint2(double first[2], double second[2], int &face_id, double alpha[3], int type)
{

	icVector2 VecAtPoint = GetVectorAtPoints(face_id, alpha);

	if(length(VecAtPoint) < 1e-18 ) return false;
	
	////calling Runge Kutta to get the next point
    const int N=2;
    int i,j;
    DP eps,hdid,hnext, t = 1.0;

    Vec_DP by(N),dydx(N),dysav(N),ysav(N),yscal(N);

	ysav[0] = first[0];
	ysav[1] = first[1];

	if(type == 0)
        localderive(t, ysav, dysav);
	else
		localinverse_derive(t, ysav, dysav);

    for (i=0;i<N;i++) yscal[i]=1.0;
    
	////set a small stepsize for more accurate tracing for limit cycle detection 1/17/06
	htry = 0.8;

	/*--------------------------------------------------------------------*/
	////calling adaptive stepsize runge-kutta method here to get next step
	for (i=0;i<5;i++) {
		eps=exp(-DP(i+1));
		t = 1.0;
		for (j=0;j<N;j++) {
			by[j]=ysav[j];
			dydx[j]=dysav[j];
		}
	
	   if(type == 0)
		   rkqs(by,dydx,t,htry,eps,yscal,hdid,hnext,localderive);
	   else
		   rkqs(by,dydx,t,htry,eps,yscal,hdid,hnext,localinverse_derive);

	}

	////Set the new stepsize, note that this stepsize will affect the limit cycle detection
	htry = hnext;
	if(hnext >= 0.8)
		htry = 0.8;

	/////store some important information
	second[0] = by[0];
	second[1] = by[1];
	/*--------------------------------------------------------------------*/

	return true;
}


////Larger step for streamline method 2/19/06
bool ToNextPoint3(double first[2], double second[2], int &face_id, double alpha[3], int type)
{

	icVector2 VecAtPoint = GetVectorAtPoints(face_id, alpha);

	if(length(VecAtPoint) < 1e-18 ) return false;
	
	////calling Runge Kutta to get the next point
    const int N=2;
    int i,j;
    DP eps,hdid,hnext, t = 1.0;

    Vec_DP by(N),dydx(N),dysav(N),ysav(N),yscal(N);

	ysav[0] = first[0];
	ysav[1] = first[1];

	if(type == 0)
        localderive(t, ysav, dysav);
	else
		localinverse_derive(t, ysav, dysav);

    for (i=0;i<N;i++) yscal[i]=1.0;
    
	////set a small stepsize for more accurate tracing for limit cycle detection 1/17/06
	//htry = 3;

	/*--------------------------------------------------------------------*/
	////calling adaptive stepsize runge-kutta method here to get next step
	for (i=0;i<3;i++) {
		eps=exp(-DP(i+1));
		t = 1.0;
		for (j=0;j<N;j++) {
			by[j]=ysav[j];
			dydx[j]=dysav[j];
		}
	
	   if(type == 0)
		   rkqs(by,dydx,t,htry,eps,yscal,hdid,hnext,localderive);
	   else
		   rkqs(by,dydx,t,htry,eps,yscal,hdid,hnext,localinverse_derive);

	}

	////Set the new stepsize, note that this stepsize will affect the limit cycle detection
	htry = hnext;
	if(hnext >= 2)
		htry = 2;

	/////store some important information
	second[0] = by[0];
	second[1] = by[1];
	/*--------------------------------------------------------------------*/

	return true;
}



////Using euler method to perform integration for limit cycle detection
////The stepsize will definitely affect the detection of limit cycle
////Added at 1/17/06
bool Euler_ToNextPoint(double first[2], double second[2], int &face_id, double alpha[3], int type)
{

	icVector2 VecAtPoint = GetVectorAtPoints(face_id, alpha);

	if(length(VecAtPoint) < 1e-15 ) return false;
	
	////Using first order Euler method to test 12/11/05
	if(type == 0)
	{
		second[0] = first[0] + VecAtPoint.entry[0]/5;
		second[1] = first[1] + VecAtPoint.entry[1]/5;
	}
	else
	{
		second[0] = first[0] - VecAtPoint.entry[0]/5;
		second[1] = first[1] - VecAtPoint.entry[1]/5;
	}

	////How about second order Euler method? 1/9/06
	//SecondOrderEulerStep(first, second, VecAtPoint, face_id, alpha, type);

	return true;
}

/*****************************************************************
store temp curve point array to global line segment array
*****************************************************************/

bool StoreToGlobalList(CurvePoints *temp, int num)
{
	int i;
	int tempid = num_linesegs_curtraj[cur_traj_index];
	icVector3 templ;


	////if the number of the line segements over the maximum number of the line segments each trajectory can store
	////extend the space for each trajectory
	if(tempid + num - 1 >= MaxNumLinesegsPerTraj - 1)
	{
		MaxNumLinesegsPerTraj += 200;
		for(i = 0; i < MaxNumTrajectories; i++)
		{
			 trajectories[i] = (LineSeg*) realloc(trajectories[i], sizeof(LineSeg) * MaxNumLinesegsPerTraj);
			 if(trajectories[i] == NULL)
			 {
				 MessageBox(NULL, "Not enough memory!", "Error", MB_OK);
				 return false;
			 }
		}
	}

	////using new trajectory data structure
	//if(tempid + num - 1 >= trajectories2[cur_traj_index].cur_MaxNumLinesegs)
	//{
	//	trajectories2[cur_traj_index].cur_MaxNumLinesegs += 200;
	//	trajectories2[cur_traj_index].line_segs = (LineSeg *)realloc( trajectories2[cur_traj_index].line_segs, \
	//		sizeof(LineSeg) * trajectories2[cur_traj_index].cur_MaxNumLinesegs);
	//	
	//	if(trajectories2[cur_traj_index].line_segs == NULL)
	//	{
	//		MessageBox(NULL, "Not enough memory!", "Error", MB_OK);
	//		return false;
	//	}
	//}

	for( i = 0; i < num-1; i++)
	{
		////Build the line segment
		trajectories[cur_traj_index][tempid+i].gstart[0] = temp[i].gpx;
		trajectories[cur_traj_index][tempid+i].gstart[1] = temp[i].gpy;
		trajectories[cur_traj_index][tempid+i].start[0] = temp[i].lpx;
		trajectories[cur_traj_index][tempid+i].start[1] = temp[i].lpy;

 		trajectories[cur_traj_index][tempid+i].gend[0] = temp[i+1].gpx;
		trajectories[cur_traj_index][tempid+i].gend[1] = temp[i+1].gpy;
		trajectories[cur_traj_index][tempid+i].end[0] = temp[i+1].lpx;
		trajectories[cur_traj_index][tempid+i].end[1] = temp[i+1].lpy;

		//templ.entry[0] = temp[i+1].lpx - temp[i].lpx;
		//templ.entry[1] = temp[i+1].lpy - temp[i].lpy;
		templ.entry[0] = temp[i+1].gpx - temp[i].gpx;
		templ.entry[1] = temp[i+1].gpy - temp[i].gpy;
		trajectories[cur_traj_index][tempid+i].length = length(templ);

		trajectories[cur_traj_index][tempid+i].Triangle_ID = temp[i].triangleid;
	}

	num_linesegs_curtraj[cur_traj_index] = tempid + num - 1;
	return true;
}


double GetSepLength(int trajID)
{
	int i;

	double sum = 0.;

	for(i = 0; i < num_linesegs_curtraj[trajID]; i++)
	{
		sum += trajectories[trajID][i].length;
	}

	return sum;
}

/*****************************************************************
store temp curve point array to global array
Entry:
Output:
*****************************************************************/
void GetNextTriangle(int &face_id, double pre[2], double cur[2], double param_t[2], int type, 
					 int &PassVertornot, double alpha[3])
{
	int which_edge = -1;

	int prev_face_id = face_id;

	Face *prev_face = Object.flist[face_id];

	Vertex *vert = NULL;

	PassVertornot = 0;
	
	////We should put pass vertex testing here before testing crossing edge
	CrossVertex2(face_id, cur, pre, type, PassVertornot);
	if(PassVertornot > 0)
	{
		return ;
	}

	face_id = prev_face_id;  //////added on 06/08/05

	CrossBoundary3(pre, cur, face_id, alpha, which_edge, param_t);


	if(param_t[0] == -1 && param_t[1] == -1)
	{
		face_id = prev_face_id;   ////something wrong here
		return;
	}

	////if not passing a vertex, judge which triangle it will enter later
	PassEdge(face_id, which_edge);

}


////Find the index of boundary that the trajectory will cross
void  CrossBoundary(int prev_id,  double point2D[2], int &which_edge)
{
	double theta02, theta0p;
	Face *prev = Object.flist[prev_id];

	theta02 = atan2(prev->xy[2][1], prev->xy[2][0]);

	if(theta02 < 0) theta02 += 2*M_PI;

	theta0p = atan2(point2D[1], point2D[0]);

	if(theta0p < 0) theta0p += 2*M_PI;

	////Get the next triangle pointer in different cases

	if(prev->xy[2][1] > 0)
	{
		if( theta0p > 0 && theta0p < theta02)
		{
			////Next face should be the opposite triangle to vertex verts[0]
			which_edge = 0;
		}
		else if( theta0p > 0 && theta0p < M_PI)
		{
			////Next face should be the opposite triangle to vertex verts[1]
			which_edge = 1;
		}
		else
		{
			////Next face should be the opposite triangle to vertex verts[2]
			which_edge = 2;
		}
	}

	else{
		if( theta0p > theta02 && theta0p < 2*M_PI)
		{
			////Next face should be the opposite triangle to vertex verts[0]
			which_edge = 0;
		}
		else if(theta0p > M_PI)
		{
			////Next face should be the opposite triangle to vertex verts[1]
			which_edge = 1;
		}
		else
		{
			////Next face should be the opposite triangle to vertex verts[2]
			which_edge = 2;
		}
	}
}


////New added at 4/29/06
void  CrossBoundary3(double pre[2], double cur[2], int face_id, double alpha[3], int &which_edge, double t[2])
{
	Face *face = Object.flist[face_id];

	if(alpha[0] < 0 && alpha[1] < 0)
	{
		//Calculate the intersection with edge v0v2


		if(GetIntersection2(pre, cur, face->xy[0], face->xy[2], t)==1)
		{
			which_edge = 1;
			cur[0] = face->xy[0][0] + t[1] * (face->xy[2][0] - face->xy[0][0]);
			cur[1] = face->xy[0][1] + t[1] * (face->xy[2][1] - face->xy[0][1]);
			return;
		}
		
		//Calculate the intersection with edge v1v2

		if(GetIntersection2(pre, cur, face->xy[1], face->xy[2], t)==1)
		{
			which_edge = 0;
			cur[0] = face->xy[1][0] + t[1] * (face->xy[2][0] - face->xy[1][0]);
			cur[1] = face->xy[1][1] + t[1] * (face->xy[2][1] - face->xy[1][1]);
			return;
		}
	}

	else if(alpha[0] < 0 && alpha[2] < 0)
	{
		//Calculate the intersection with edge v0v1

		if(GetIntersection2(pre, cur, face->xy[0], face->xy[1], t)==1)
		{
			which_edge = 2;
			cur[0] = face->xy[0][0] + t[1] * (face->xy[1][0] - face->xy[0][0]);
			cur[1] = face->xy[0][1] + t[1] * (face->xy[1][1] - face->xy[0][1]);
			return;
		}
		
		//Calculate the intersection with edge v1v2
		
		if(GetIntersection2(pre, cur, face->xy[1], face->xy[2], t)==1)
		{
			which_edge = 0;
			cur[0] = face->xy[1][0] + t[1] * (face->xy[2][0] - face->xy[1][0]);
			cur[1] = face->xy[1][1] + t[1] * (face->xy[2][1] - face->xy[1][1]);
			return;
		}
	}

	else if(alpha[1] < 0 && alpha[2] < 0)
	{
		//Calculate the intersection with edge v0v1

		if(GetIntersection2(pre, cur, face->xy[0], face->xy[1], t)==1)
		{
			which_edge = 2;
			cur[0] = face->xy[0][0] + t[1] * (face->xy[1][0] - face->xy[0][0]);
			cur[1] = face->xy[0][1] + t[1] * (face->xy[1][1] - face->xy[0][1]);
			return;
		}
		
		
		if(GetIntersection2(pre, cur, face->xy[0], face->xy[2], t)==1)
		{
			which_edge = 1;
			cur[0] = face->xy[0][0] + t[1] * (face->xy[2][0] - face->xy[0][0]);
			cur[1] = face->xy[0][1] + t[1] * (face->xy[2][1] - face->xy[0][1]);
			return;
		}
	}

	else if(alpha[0] < 0)
	{
		which_edge = 0;
		GetIntersection2(pre, cur, face->xy[1], face->xy[2], t);
		cur[0] = face->xy[1][0] + t[1] * (face->xy[2][0] - face->xy[1][0]);
		cur[1] = face->xy[1][1] + t[1] * (face->xy[2][1] - face->xy[1][1]);
		return;
	}

	else if(alpha[1] < 0)
	{
		which_edge = 1;
		GetIntersection2(pre, cur, face->xy[2], face->xy[0], t);
		cur[0] = face->xy[2][0] + t[1] * (face->xy[0][0] - face->xy[2][0]);
		cur[1] = face->xy[2][1] + t[1] * (face->xy[0][1] - face->xy[2][1]);
		return;
	}

	else if(alpha[2] < 0)
	{
		which_edge = 2;
		GetIntersection2(pre, cur, face->xy[0], face->xy[1], t);
		cur[0] = face->xy[0][0] + t[1] * (face->xy[1][0] - face->xy[0][0]);
		cur[1] = face->xy[0][1] + t[1] * (face->xy[1][1] - face->xy[0][1]);
		return;
	}
}


////New routine for crossing vertex testing
void CrossVertex(int &face_id, double cur_p[2], double pre_p[2],int type, int &passornot)
{
//	double tx, ty;
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

	if(fabs(pending) < 1e-10) ////passing the vertex
	{
		//////Set a smaller time step for integration 1/6/06
		//htry = 0.7;

		TriangleThroughVertex(crossVert, newtriangleid, type);
		face_id = newtriangleid;	
		passornot = alpha_index+1;
		return;
	}

	passornot = 0;
}


/*
New routine for crossing vertex testing 4/30/06
*/
void CrossVertex2(int &face_id, double cur_p[2], double pre_p[2],int type, int &passornot)
{
	int i;
	double vert[2];
	double max_alpha ;
    int newtriangleid = 0;
	int crossVert;
	Face *face = Object.flist[face_id];

	double A, B, C, pending;
	A = pre_p[1] - cur_p[1];
	B = cur_p[0] - pre_p[0];
	C = (pre_p[0]*cur_p[1] - cur_p[0]*pre_p[1]);

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

			TriangleThroughVertex(crossVert, newtriangleid, type);
			face_id = newtriangleid;	
			passornot = i+1;
			return;
		}
	}

	passornot = 0;
}


////Modified on 09/06/05, still not handle all the cases that pass cross vertex
void TriangleThroughVertex(int vert_id, int &theone, int type)
{
	Vertex *vert = Object.vlist[vert_id];

	////using a tricky way to get next triangle

	////1. move the vertex along the vector direction a little bit
	//double step_v[2];
	//icVector2 vec = vert->vec;
	//normalize(vec);

	//if(type == 0)
	//{
	//	step_v[0] = vert->x + vec.entry[0]*triangle_approx_size/8.;
	//	step_v[1] = vert->y + vec.entry[1]*triangle_approx_size/8.;
	//}
	//else
	//{
	//	step_v[0] = vert->x - vec.entry[0]*triangle_approx_size/8.;
	//	step_v[1] = vert->y - vec.entry[1]*triangle_approx_size/8.;
	//}

	//theone = TriangleDetect(step_v[0], step_v[1]); ////not a stable method, we need to use geometric method

	////////Testing codes here 07/06/05
	//if(theone == -1)
	//{
	//	//if it is really outof the whole mesh
	//	theone = theone;

	//	//otherwise, we may need the help of the corner table method 12/29/05
	//}


	/*---------------------------------------------------------------------*/
	////Second method to judge the next triangle when the tracing curve
	////passes through a vertex

	int i;
	double vang;
	Face *face;
	Corner *c;
	int NewTID = -1;
	int orient;
	icVector2 vec = vert->vec;
	normalize(vec);

	vang = atan2(vec.entry[1], vec.entry[0]);

	if(type == 1) // consider the inverse of the flow
		vang += M_PI;

	if(vang < 0)
		vang += 2*M_PI;

	if(vang >= 2*M_PI)
		vang -= 2*M_PI;

	////Get the orientation of the angle allocation
	//orient = Object.clist[vert->Corners[0]]->orient;

	/* Change the orientation judgement 01/22/07 */
	if(Object.flist[Object.clist[vert->Corners[0]]->t]->xy[2][1] < 0) 
		orient = 1;
	else
		orient = 0;

	//orient = vert->oi;

	for( i = 0; i < vert->Num_corners; i++)
	{
		c = Object.clist[vert->Corners[i]];
		////first, we check which angle area the vector on the vertex will fall in
		if(orient > 0)
		{
			if(c->BeginAng > c->EndAng)
			{
				if((vang >= c->BeginAng && vang < 2*M_PI) 
					|| (vang < c->EndAng && vang >= 0 ))
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
				if((vang <= c->BeginAng && vang >= 0)|| (vang > c->EndAng && vang < 2*M_PI))
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


/********************************************************************
If we have already judge which edge the curve will cross
We can use the edge information to get next triangle
********************************************************************/
void PassEdge(int &face_id, int which_edge)
{
	////Using edge information to get next triangle

	Face *face = Object.flist[face_id];

	Edge *theone ;
	int vertindex;

	for(int i = 0; i < 3; i++)
	{
		vertindex = face->verts[which_edge];

		theone = face->edges[i];

		if(theone->verts[0] != vertindex && theone->verts[1] != vertindex)
			break;
	}

	if(theone->tris[0] != face_id)
		face_id = theone->tris[0];
	else
		face_id = theone->tris[1];

}

/***************************************************************
Get the intersection t values when tracing curve crossing edges
****************************************************************/
void GetIntersection(double p1[2], double p2[2], double q1[2], double q2[2], double t[2])
{
	double deta1 = p2[0] - p1[0];
	double deta2 = p2[1] - p1[1];
	double deta3 = q2[0] - q1[0];
	double deta4 = q2[1] - q1[1];

	if(fabs(deta1*deta4 - deta2*deta3) < 1e-40) 
	{
		t[0] = t[1] = -1;
		return;
	}

	////t1 is the parameter value of line p1p2 

	double t1 = (q1[0] * deta4 - q1[1] * deta3 - p1[0] * deta4 + p1[1] * deta3)\
		/(deta1 * deta4 - deta2 * deta3);
	
	////t2 is the parameter value of line q1q2 
	double t2 = (p1[0] * deta2 - p1[1] * deta1 - q1[0] * deta2 + q1[1] * deta1)\
		/(deta3 * deta2 - deta4 * deta1);

	////t2 may be useful in future

	t[0] = t1;
	t[1] = t2;

}



/***************************************************************
New method to calculate the intersection of two line segments
****************************************************************/

/* meaning of return value
 0----Intersection dosn't exists                                                   
 1----Intersection exists.                                                        
 2----two line segments are parallel.                                         
 3----two line segments are collinear, but not overlap.                      
 4----two line segments are collinear, and share one same end point.       
 5----two line segments are collinear, and overlap.                           
*/    

int GetIntersection2(double PointA[2], double PointB[2], double PointC[2], double PointD[2], double t[2])
{

    double delta;
    double t1,t2;
    double a,b,c,d;
    double xba,yba,xdc,ydc,xca,yca;

    xba=PointB[0]-PointA[0];    yba=PointB[1]-PointA[1];
    xdc=PointD[0]-PointC[0];    ydc=PointD[1]-PointC[1];
    xca=PointC[0]-PointA[0];    yca=PointC[1]-PointA[1];

    delta=xba*ydc-yba*xdc;
    t1=xca*ydc-yca*xdc;
    t2=xca*yba-yca*xba;

    if(delta!=0)
    {
        t[0]=t1/delta;   t[1]=t2/delta;
        /*two segments intersect (including intersect at end points)*/
        //if ( t[0]<=1 && t[0]>=0 && t[1]<=1 && t[1]>=0 ) return 1;
        if ( t[0]<=1 && (t[0]>=0 || fabs (t[0])<=1.e-8)
			&& t[1]<=1 && (t[1]>=0|| fabs (t[1])<=1.e-8)) //set threshold to allow some numerical errors
			return 1;
        else return 0; 
    }

    else
    {       
        /* AB & CD are parallel. */
        if ( (t1!=0) && (t2!=0) ) return 2;

        /* when AB & CD are collinear */

        /*if AB isn't a vertical line segment, project to x-axis */
        if(PointA[0]!=PointB[0])   
        {
            a=MIN(PointA[0],PointB[0]); b=MAX(PointA[0],PointB[0]);
            c=MIN(PointC[0],PointD[0]); d=MAX(PointC[0],PointD[0]);

            if ( (d<a) || (c>b) ) return  3;
            else if( (d==a) || (c==b) ) return 4;  
            else return 5;
        }

        else         /* if AB is a vertical line segment, project to y-axis */  
        {

            a=MIN(PointA[1],PointB[1]); b=MAX(PointA[1],PointB[1]);
            c=MIN(PointC[1],PointD[1]); d=MAX(PointC[1],PointD[1]); 

            if( (d<a) || (c>b) ) return  3;
            else if( (d==a) || (c==b) ) return 4;
            else return 5;
        }
    }
}

int cal_intersect(double PointA[2], double PointB[2], double PointC[2], double PointD[2], double t[2])
{

    double delta;
    double t1,t2;
    double a,b,c,d;
    double xba,yba,xdc,ydc,xca,yca;

    xba=PointB[0]-PointA[0];    yba=PointB[1]-PointA[1];
    xdc=PointD[0]-PointC[0];    ydc=PointD[1]-PointC[1];
    xca=PointC[0]-PointA[0];    yca=PointC[1]-PointA[1];

    delta=xba*ydc-yba*xdc;
    t1=xca*ydc-yca*xdc;
    t2=xca*yba-yca*xba;

    if(delta!=0)
    {
        t[0]=t1/delta;   t[1]=t2/delta;
        /*two segments intersect (including intersect at end points)*/
        //if ( t[0]<=1 && t[0]>=0 && t[1]<=1 && t[1]>=0 ) return 1;

		/*  to account for some numerical errors, we set a larger threshold  */
        if ( t[0]<=1+1.e-4 && (t[0]>=0 || fabs (t[0])<=1.e-4)
			&& t[1]<=1+1.e-4 && (t[1]>=0|| fabs (t[1])<=1.e-4)) //set threshold to allow some numerical errors
			return 1;
        else return 0; 
    }

    else
    {       
        /* AB & CD are parallel. */
        if ( (t1!=0) && (t2!=0) ) return 2;

        /* when AB & CD are collinear */

        /*if AB isn't a vertical line segment, project to x-axis */
        if(PointA[0]!=PointB[0])   
        {
            a=MIN(PointA[0],PointB[0]); b=MAX(PointA[0],PointB[0]);
            c=MIN(PointC[0],PointD[0]); d=MAX(PointC[0],PointD[0]);

            if ( (d<a) || (c>b) ) return  3;
            else if( (d==a) || (c==b) ) return 4;  
            else return 5;
        }

        else         /* if AB is a vertical line segment, project to y-axis */  
        {

            a=MIN(PointA[1],PointB[1]); b=MAX(PointA[1],PointB[1]);
            c=MIN(PointC[1],PointD[1]); d=MAX(PointC[1],PointD[1]); 

            if( (d<a) || (c>b) ) return  3;
            else if( (d==a) || (c==b) ) return 4;
            else return 5;
        }
    }
}


int cal_intersect_2(double PointA[2], double PointB[2], double PointC[2], double PointD[2], double t[2])
{

    double delta;
    double t1,t2;
    double a,b,c,d;
    double xba,yba,xdc,ydc,xca,yca;

    xba=PointB[0]-PointA[0];    yba=PointB[1]-PointA[1];
    xdc=PointD[0]-PointC[0];    ydc=PointD[1]-PointC[1];
    xca=PointC[0]-PointA[0];    yca=PointC[1]-PointA[1];

    delta=xba*ydc-yba*xdc;
    t1=xca*ydc-yca*xdc;
    t2=xca*yba-yca*xba;

    if(delta!=0)
    {
        t[0]=t1/delta;   t[1]=t2/delta;
        /*two segments intersect (including intersect at end points)*/
        //if ( t[0]<=1 && t[0]>=0 && t[1]<=1 && t[1]>=0 ) return 1;

		/*  to account for some numerical errors, we set a larger threshold  */
        if ( t[0]<=1+1.e-5 && (t[0]>=0 || fabs (t[0])<=1.e-5)
			&& t[1]<=1+1.e-5 && (t[1]>=0|| fabs (t[1])<=1.e-5)) //set threshold to allow some numerical errors
			return 1;
        else return 0; 
    }

    else
    {       
        /* AB & CD are parallel. */
        if ( (t1!=0) && (t2!=0) ) return 2;

        /* when AB & CD are collinear */

        /*if AB isn't a vertical line segment, project to x-axis */
        if(PointA[0]!=PointB[0])   
        {
            a=MIN(PointA[0],PointB[0]); b=MAX(PointA[0],PointB[0]);
            c=MIN(PointC[0],PointD[0]); d=MAX(PointC[0],PointD[0]);

            if ( (d<a) || (c>b) ) return  3;
            else if( (d==a) || (c==b) ) return 4;  
            else return 5;
        }

        else         /* if AB is a vertical line segment, project to y-axis */  
        {

            a=MIN(PointA[1],PointB[1]); b=MAX(PointA[1],PointB[1]);
            c=MIN(PointC[1],PointD[1]); d=MAX(PointC[1],PointD[1]); 

            if( (d<a) || (c>b) ) return  3;
            else if( (d==a) || (c==b) ) return 4;
            else return 5;
        }
    }
}


/*
   Calculate the major road intersections
*/
int cal_majIntersect(double PointA[2], double PointB[2], double PointC[2], 
					 double PointD[2], double t[2])
{

    double delta;
    double t1,t2;
    double a,b,c,d;
    double xba,yba,xdc,ydc,xca,yca;

    xba=PointB[0]-PointA[0];    yba=PointB[1]-PointA[1];
    xdc=PointD[0]-PointC[0];    ydc=PointD[1]-PointC[1];
    xca=PointC[0]-PointA[0];    yca=PointC[1]-PointA[1];

    delta=xba*ydc-yba*xdc;
    t1=xca*ydc-yca*xdc;
    t2=xca*yba-yca*xba;

    if(delta!=0)
    {
        t[0]=t1/delta;   t[1]=t2/delta;
        /*two segments intersect (including intersect at end points)*/
        //if ( t[0]<=1 && t[0]>=0 && t[1]<=1 && t[1]>=0 ) return 1;

		/*  to account for some numerical errors, we set a larger threshold  */
        if ( t[0]<=1+1.e-3 && (t[0]>=0 || fabs (t[0])<=1.e-3)
			&& t[1]<=1+1.e-3 && (t[1]>=0|| fabs (t[1])<=1.e-3)) //set threshold to allow some numerical errors
			return 1;
        else return 0; 
    }

    else
    {       
        /* AB & CD are parallel. */
        if ( (t1!=0) && (t2!=0) ) return 2;

        /* when AB & CD are collinear */

        /*if AB isn't a vertical line segment, project to x-axis */
        if(PointA[0]!=PointB[0])   
        {
            a=MIN(PointA[0],PointB[0]); b=MAX(PointA[0],PointB[0]);
            c=MIN(PointC[0],PointD[0]); d=MAX(PointC[0],PointD[0]);

            if ( (d<a) || (c>b) ) return  3;
            else if( (d==a) || (c==b) ) return 4;  
            else return 5;
        }

        else         /* if AB is a vertical line segment, project to y-axis */  
        {

            a=MIN(PointA[1],PointB[1]); b=MAX(PointA[1],PointB[1]);
            c=MIN(PointC[1],PointD[1]); d=MAX(PointC[1],PointD[1]); 

            if( (d<a) || (c>b) ) return  3;
            else if( (d==a) || (c==b) ) return 4;
            else return 5;
        }
    }
}


/*************************************************************
Runge Kutta integrator driver for local tracing
*************************************************************/
void localderive(const DP t, Vec_I_DP &y, Vec_O_DP &dydx)
{
	double alpha[3];
	Get2DBarycentricFacters(globalface, y[0], y[1], alpha);

	////first, we need to get the change of next points
	icVector2 v = GetVectorAtPoints(globalface, alpha);

	/*normalize it before RK  07/25/07*/
	//normalize(v);
	//v = 100*DistanceThreshold*v;

	dydx[0] = v.entry[0];
	dydx[1] = v.entry[1];

}
	
void localinverse_derive(const DP t, Vec_I_DP &y, Vec_O_DP &dydx)
{
	double alpha[3];
	Get2DBarycentricFacters(globalface, y[0], y[1], alpha);

	////first, we need to get the change of next points
	icVector2 v = GetVectorAtPoints(globalface, alpha);
	
	/*normalize it before RK  07/25/07*/
	//normalize(v);
	//v = 100*DistanceThreshold*v;

	dydx[0] = -v.entry[0];
	dydx[1] = -v.entry[1];
}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
////Routines and variables for limit cycle detection
////it will be put to limit cycle detector file
int TriangleDetect(double x, double y)
{
	GLuint  selectBuffer[SELECTBUFFERSIZE];
	GLint	vp[4] = {0, 0 ,1, 1};
	int hits;
	int selectedTriangle = -1;


	////Build the selection buffer here
	glSelectBuffer(SELECTBUFFERSIZE, selectBuffer);

   	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);

	gluPickMatrix(x, y, 1e-9, 1e-9, vp );
		
	glOrtho(0,1,  0,1,  0,50);

	////If one of the element has been selected for being edited

	CGlView::IBFVEffect(GL_SELECT);

	hits = glRenderMode(GL_RENDER);

	if(hits>0)
	{
		////Because this is 2D plane, we need not worry about the overlap of objects
		////It will have at most one object being selected at one click (it may not be true)
		int objname = selectBuffer[3];

		selectedTriangle = objname - NAMEOFTRIANGLE;

	}

	else{
		selectedTriangle = -1;
	}

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode( GL_MODELVIEW );

	return selectedTriangle;
}


/*---------------------------------------------------------------------*/
////Routines for separatrices calculation

/********************************************************
Judge whether the direction is incoming or outgoing from
the saddle points
********************************************************/
int InComing(double x, double y, double cx, double cy)
{
	//Using the eigenvalue to judge whether it is incoming or outgoing

	double  x1, y1, len1, len2;
	double mag;
	icVector2 vec[2];


	getVector(x, y, vec[0], mag);

	len1 = sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy));

	x1 = x + vec[0].entry[0] * INOROUTJUDGE;
	y1 = y + vec[0].entry[1] * INOROUTJUDGE;


	getVector(x1, y1, vec[1], mag);

	x1 = x + (vec[0].entry[0] + vec[1].entry[0])/2 * INOROUTJUDGE;
	y1 = y + (vec[0].entry[1] + vec[1].entry[1])/2 * INOROUTJUDGE;

	getVector(x1, y1, vec[1], mag);

	len2 = sqrt((x1-cx)*(x1-cx) + (y1-cy)*(y1-cy));

	if(len1 > len2)
		return 1;
	else
		return 0;

}


/********************************************************************
This routine calculates and store one single separatrix
********************************************************************/
	
void CalSingleSeparatrix(int Triangle_ID, double x, double y, int inout)
{
	int temp_index = 0;

	////Call the local tracing of trajectory
	CalLocalTracing(Triangle_ID, x, y, inout);

	cur_traj_index ++;
}


/********************************************************************
This routine calculates all the separatrices from captured saddle
********************************************************************/

void CalSeparatrices()
{
	int i;
	double newpos[2] = {0.};
	double sing_center[2] = {0.};
	icVector2 sep_vector;

	int start_triangle;

	////initialize the calculation of separatrices
	for(i = 0; i <= cur_traj_index; i++)
		num_linesegs_curtraj[i] = 0;

	////Initialize the separatrix flag for each triangle 1/3/06
	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->contain_separatrix = 0;
	}

	cur_traj_index = 0;
	cur_separatrices_index = 0;

	InitSampleptsList();

	/*** Calculate the average length of the edges 07/25/06 ***/
	icVector2 temp;
	temp.entry[0] = Object.vlist[Object.flist[0]->verts[0]]->x - 
		Object.vlist[Object.flist[0]->verts[1]]->x;
	temp.entry[1] = Object.vlist[Object.flist[0]->verts[0]]->y - 
		Object.vlist[Object.flist[0]->verts[1]]->y;
	ave_length = length(temp);
	/////////////////////////////////////////////////////

	/////

	int flag = -1;
	for(i = 0; i < cur_singularity_index ; i++)
	{
		if(singularities[i].type == SADDLE)
		{
			if(cur_separatrices_index >= MaxNumSeparatrices - 1)
			{
				MaxNumSeparatrices += 25;
				separatrices = (Separatrices*)realloc(separatrices, sizeof(Separatrices)*MaxNumSeparatrices);
			}

			////store the seperatrix group index to the corresponding saddle 10/13/05
			singularities[i].separtices = cur_separatrices_index;

			sing_center[0] = singularities[i].gcx;
			sing_center[1] = singularities[i].gcy;

			////Calculate the positive outgoing separatrix
			sep_vector = singularities[i].outgoing;
			start_triangle = singularities[i].Triangle_ID;
			FixSepBeginningPos(start_triangle, sep_vector, sing_center, newpos);
			htry = 0.1;
			flag = -1;
			//CalSingleSeparatrix(start_triangle, newpos[0], newpos[1], 0);
			TraceAvoidClosed(newpos, start_triangle, ave_length/4, ave_length/4, 0, flag);
			separatrices[cur_separatrices_index].sep1 = cur_traj_index-1;
			if(flag == 5)  //form a loop, search for the connected limit cycle
			{
				//Get the last point in local
				double p[2];
				p[0] = trajectories[cur_traj_index-1][num_linesegs_curtraj[cur_traj_index-1]-1].start[0];
				p[1] = trajectories[cur_traj_index-1][num_linesegs_curtraj[cur_traj_index-1]-1].start[1];
				int Cur_t = trajectories[cur_traj_index-1][num_linesegs_curtraj[cur_traj_index-1]-1].Triangle_ID;
				CloseToLimitCycle(i, p, Cur_t, 0, 1);
			}

			////Get and save the length of the separatrix  08/10/06
			separatrices[cur_separatrices_index].length1 = GetSepLength(cur_traj_index-1);
			
			////Calculate the positive incoming separatrix
			sep_vector = singularities[i].incoming;
			start_triangle = singularities[i].Triangle_ID;
			FixSepBeginningPos(start_triangle, sep_vector, sing_center, newpos);
			htry = 0.1;
			flag = -1;
			//CalSingleSeparatrix(start_triangle, newpos[0], newpos[1], 1);
			TraceAvoidClosed(newpos, start_triangle, ave_length/4, ave_length/4, 1, flag);
			separatrices[cur_separatrices_index].sep2 = cur_traj_index-1;
			if(flag == 5)  //form a loop, search for the connected limit cycle
			{
				//Get the last point in local
				double p[2];
				p[0] = trajectories[cur_traj_index-1][num_linesegs_curtraj[cur_traj_index-1]-1].start[0];
				p[1] = trajectories[cur_traj_index-1][num_linesegs_curtraj[cur_traj_index-1]-1].start[1];
				int Cur_t = trajectories[cur_traj_index-1][num_linesegs_curtraj[cur_traj_index-1]-1].Triangle_ID;
				CloseToLimitCycle(i, p, Cur_t, 1, 2);
			}
			////Get and save the length of the separatrix  08/10/06
			separatrices[cur_separatrices_index].length2 = GetSepLength(cur_traj_index-1);
			
			////Calculate the negative outgoing separatrix
			sep_vector = -singularities[i].outgoing;
			start_triangle = singularities[i].Triangle_ID;
			FixSepBeginningPos(start_triangle, sep_vector, sing_center, newpos);
			htry = 0.1;
			flag = -1;
			//CalSingleSeparatrix(start_triangle, newpos[0], newpos[1], 0);
			TraceAvoidClosed(newpos, start_triangle, ave_length/4, ave_length/4, 0, flag);
			separatrices[cur_separatrices_index].sep3 = cur_traj_index-1;
			if(flag == 5)  //form a loop, search for the connected limit cycle
			{
				//Get the last point in local
				double p[2];
				p[0] = trajectories[cur_traj_index-1][num_linesegs_curtraj[cur_traj_index-1]-1].start[0];
				p[1] = trajectories[cur_traj_index-1][num_linesegs_curtraj[cur_traj_index-1]-1].start[1];
				int Cur_t = trajectories[cur_traj_index-1][num_linesegs_curtraj[cur_traj_index-1]-1].Triangle_ID;
				CloseToLimitCycle(i, p, Cur_t, 0, 3);
			}
			////Get and save the length of the separatrix  08/10/06
			separatrices[cur_separatrices_index].length3 = GetSepLength(cur_traj_index-1);

			////Calculate the negative incoming separatrix
			sep_vector = -singularities[i].incoming;
			start_triangle = singularities[i].Triangle_ID;
			FixSepBeginningPos(start_triangle, sep_vector, sing_center, newpos);
			htry = 0.1;
			flag = -1;
		    //CalSingleSeparatrix(start_triangle, newpos[0], newpos[1], 1);
			TraceAvoidClosed(newpos, start_triangle, ave_length/4, ave_length/4, 1, flag);
			separatrices[cur_separatrices_index].sep4 = cur_traj_index-1;
			if(flag == 5)  //form a loop, search for the connected limit cycle
			{
				//Get the last point in local
				double p[2];
				p[0] = trajectories[cur_traj_index-1][num_linesegs_curtraj[cur_traj_index-1]-1].start[0];
				p[1] = trajectories[cur_traj_index-1][num_linesegs_curtraj[cur_traj_index-1]-1].start[1];
				int Cur_t = trajectories[cur_traj_index-1][num_linesegs_curtraj[cur_traj_index-1]-1].Triangle_ID;
				CloseToLimitCycle(i, p, Cur_t, 1, 4);
			}
			////Get and save the length of the separatrix  08/10/06
			separatrices[cur_separatrices_index].length4 = GetSepLength(cur_traj_index-1);
			

			cur_separatrices_index++;


		}
	}

	pre_cur_traj_index = cur_traj_index;

}



/////////////////////////////////////////////////////////////////////
/************************************************************
Judge whether a point falls into a specific triangle
using the barycentric coordinates of the point under the 
local frame of the triangle   1/1/06
*************************************************************/
bool FallIntheTriangle(int triangleID, double pos[2])
{
	double alpha[3] = {0.};
	Face *face = Object.flist[triangleID];

	icVector2 VP;

	VP.entry[0] = pos[0] - Object.vlist[face ->verts[0]]->x;
	VP.entry[1] = pos[1] - Object.vlist[face ->verts[0]]->y;

	///1. transfer to local coordinates
	double a = dot(VP, face->LX);
	double b = dot(VP, face->LY);

	///2. calculate the barycentric coordinates of the point under the local frame of the triangle
	Get2DBarycentricFacters(triangleID, a, b, alpha);

	///3. judge whether the point falls into the triangle
	if( (alpha[0] >= 0 /*|| alpha[0] > 1e-8*/) && alpha[0] <= 1 
		&& (alpha[1] >= 0 /*|| alpha[1] > 1e-8*/) && alpha[1] <= 1 
		&& (alpha[2] >= 0 /*|| alpha[2] > 1e-8*/) && alpha[2] <= 1)
		return true;
	else
		return false;

}


/*Find the smallest edge distance in the triangle*/
double get_shortestedge_tri(int tri)
{
	int i;
	Face *face = Object.flist[tri];
	double shortest = face->edges[0]->length;

	for(i=1; i<3; i++)
	{
		if(face->edges[i]->length < shortest)
			shortest = face->edges[i]->length;
	}

	return shortest;
}

/****************************************************************
Fix the beginning position of the calculation of a separatrix
1/1/06
****************************************************************/

void FixSepBeginningPos(int &triangleID, icVector2 sep_vector, double sing_center[2], double newpos[2])
{
	int count = 2;

	double shortest = get_shortestedge_tri(triangleID)/5;

	////Get the first position of the beginning point
	//newpos[0] = sing_center[0] + SEPARATRIXSTEP * sep_vector.entry[0];
	//newpos[1] = sing_center[1] + SEPARATRIXSTEP * sep_vector.entry[1];
	newpos[0] = sing_center[0] + shortest * sep_vector.entry[0];
	newpos[1] = sing_center[1] + shortest * sep_vector.entry[1];

	while(count <= 100)
	{
		if(FallIntheTriangle(triangleID, newpos))
			return;
	    
		//newpos[0] = sing_center[0] + SEPARATRIXSTEP/(2*count) * sep_vector.entry[0];
		//newpos[1] = sing_center[1] + SEPARATRIXSTEP/(2*count) * sep_vector.entry[1];
		//newpos[0] = sing_center[0] + shortest/(2*count) * sep_vector.entry[0];
		//newpos[1] = sing_center[1] + shortest/(2*count) * sep_vector.entry[1];
		newpos[0] = sing_center[0] + shortest/(count) * sep_vector.entry[0];
		newpos[1] = sing_center[1] + shortest/(count) * sep_vector.entry[1];

		count++;
	}

	////if we can not find a point inside current triangle, probably we can find a point
	////falls in the neighboring triangle

	////we can use the similar idea of tracing and passing an edge of currrent triangle 07/25/06

}


/****************************************************************
Find the beginning position of the calculation of a separatrix
07/21/07
****************************************************************/

void FixSepBeginningPos_2(int &triangleID, icVector2 sep_vector, double sing_center[2], double newpos[2])
{
	int count = 2;

	double shortest = get_shortestedge_tri(triangleID);

	////Get the first position of the beginning point
	//newpos[0] = sing_center[0] + SEPARATRIXSTEP * sep_vector.entry[0];
	//newpos[1] = sing_center[1] + SEPARATRIXSTEP * sep_vector.entry[1];
	newpos[0] = sing_center[0] + shortest * sep_vector.entry[0];
	newpos[1] = sing_center[1] + shortest * sep_vector.entry[1];

	while(count <= 500)
	{
		if(FallIntheTriangle(triangleID, newpos))
			return;
	    
		//newpos[0] = sing_center[0] + SEPARATRIXSTEP/(2*count) * sep_vector.entry[0];
		//newpos[1] = sing_center[1] + SEPARATRIXSTEP/(2*count) * sep_vector.entry[1];
		newpos[0] = sing_center[0] + shortest/(2*count) * sep_vector.entry[0];
		newpos[1] = sing_center[1] + shortest/(2*count) * sep_vector.entry[1];

		count++;
	}

	////if we can not find a point inside current triangle, probably we can find a point
	////falls in the neighboring triangle

	////we can use the similar idea of tracing and passing an edge of currrent triangle 07/25/06

}