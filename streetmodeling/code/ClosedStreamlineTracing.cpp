/*
ClosedStreamlineTracing.cpp

This file contains routines fulfilling the tracing of those closed streamlines, which means that they 
may form a loop finally.
The routines here try to avoid tracing several rounds of the same loop.
We use the ideas from the evenly placed streamlines.
*/

#include "stdafx.h"

#include "ClosedStreamlineTracing.h"
#include "VFDataStructure.h"

#include "LocalTracing.h"

#include "SCCCycledetect.h"

SamplePtsList *sampleptslist;       //each streamline will have a sample point list
int num_samplelists;
int MaxNumSampleLists;


extern int MaxNumSingularities;                 //Maximum number of being captured singularities

extern int MaxNumTrajectories;
extern int MaxNumLinesegsPerTraj;
extern int MaxNumSeparatrices;                   //Maximum number of group of separatrices
extern LimitCycle *limitcycles;                 //limit cycles data structure
extern int cur_limitcycle_index;

extern Polygon3D Object;
//extern Trajectory *trajectories2;
extern LineSeg **trajectories;                 //trajectories' list
extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory
extern Singularities *singularities;           //being captured singularites' list
extern int cur_singularity_index;
extern Separatrices *separatrices;             //array for group of separatrices
extern int cur_separatrices_index;

extern int globalface;

extern double g_dt; /*02/27/07*/

extern double sum_flow_length;

extern bool StoreToGlobalList(CurvePoints *temp, int num);
/*******************************/
//Global variables for marking the flags
int stop_triangle = -1;
double stop_p_loc[2] = {0.};

/*
Get one sample point on the regular streamline, not periodic orbit
*/
bool GetOneSamplePointfromAStreamline(int traj, int &cur_lineindex, int &movetonext,
					   double prept[2], double curpt[2], 
					   double interval, double &cur_length)
{
	int i;

	icVector2 len_vec;
	double alpha;

	int num_lines = num_linesegs_curtraj[traj];
	
	curpt[0] = curpt[1] = -1;

	if(cur_length >= interval)
	{
		alpha = (cur_length-interval)/trajectories[traj][cur_lineindex].length;
		curpt[0] = alpha*trajectories[traj][cur_lineindex].gstart[0] 
		          + (1-alpha)*trajectories[traj][cur_lineindex].gend[0];
		curpt[1] = alpha*trajectories[traj][cur_lineindex].gstart[1] 
		          + (1-alpha)*trajectories[traj][cur_lineindex].gend[1];

		cur_length -= interval;
		return true;
	}

	else{
		cur_lineindex++;
		if(cur_lineindex >= num_lines)
		{
			return false;
		}

		cur_length += trajectories[traj][cur_lineindex].length;
		return false;
	}

}



/*
Get a set of sample points during tracing a new streamline
This sample points are temporary and for preventing the closed loop of a streamline
It will keep finding all the samples from 'cur_line' line segment till the end of 
current streamline
*/
SamplePts *GetSamplePtsWhenTracing(int traj, double interval, int &cur_line, int &movetonext, double &cur_length, 
							 SamplePts *samples, int &num_samples, int &MaxNum)
{

	while(cur_line < num_linesegs_curtraj[traj] )
	{
		if(GetOneSamplePointfromAStreamline(traj, cur_line, movetonext, samples[num_samples-1].gpt,
			samples[num_samples].gpt, interval, cur_length))
		{
			samples[num_samples].triangle = trajectories[traj][cur_line].Triangle_ID;

			if(samples[num_samples].triangle < 0 || samples[num_samples].triangle > Object.nfaces)
			{
				continue;
			}
			num_samples++;

			if(num_samples > MaxNum)
			{
				SamplePts *temp;
				if((temp = (SamplePts *)realloc((SamplePts*)samples, 
					sizeof(SamplePts)* (MaxNum+200))) == NULL)
				{
					MessageBox(NULL, "Memory reallocation failed!", "Error", MB_OK);
				}
				else
				{
					samples = temp;
					MaxNum+=200;
				}
			}
		}

	}

	cur_line--;

	if(cur_line < 0)
		cur_line = 0;

	return samples;
}



/*
Judge whether a input point is too close to a set of sample points
*/
bool CloseToCurSamplePt(double p[3], int triangle, SamplePts *samples, int num_samples,
						double separate_dist, double sample_interval)
{
	int i;
	icVector2 dis;
	int stop_locate = floor(separate_dist/sample_interval)+20;

	for(i = 0 ; i < num_samples-stop_locate; i++)
	{
		dis.entry[0] = p[0] - samples[i].gpt[0];
		dis.entry[1] = p[1] - samples[i].gpt[1];

		if(length(dis) < separate_dist)
			return true;
	}
	return false;
}


/*
Avoid the separation calculation loop several rounds
*/
bool TraceAvoidClosed(double seed_p[3], int triangle, 
					double Sample_interval, double loopdsep, int type, int &flag)
{
	int i;
	flag = -1;

	int pre_face, cur_face;
	double globalp[3] = {0.};
	int cur_line = 0;
	double cur_length = 0;
	int movetonext = 0;
		
	if(triangle < 0 || triangle > Object.nfaces)
		return false;

	pre_face = cur_face = triangle;
	globalp[0] = seed_p[0];
	globalp[1] = seed_p[1];


	num_linesegs_curtraj[cur_traj_index] = 0;

	////We put the allocation of sample point list here
	if(sampleptslist[cur_traj_index].samples != NULL)
		free(sampleptslist[cur_traj_index].samples);
	 
	sampleptslist[cur_traj_index].MaxNumSamples = 1800;
	sampleptslist[cur_traj_index].samples = (SamplePts *)malloc(sizeof(SamplePts)*
		sampleptslist[cur_traj_index].MaxNumSamples);
	sampleptslist[cur_traj_index].colorindex = cur_traj_index*rand()%3;

	//////////Add the first point into the sampling point list 
	sampleptslist[cur_traj_index].num_samples = 0;

	sampleptslist[cur_traj_index].samples[0].gpt[0] = seed_p[0];
	sampleptslist[cur_traj_index].samples[0].gpt[1] = seed_p[1];
	sampleptslist[cur_traj_index].samples[0].triangle = triangle;
	sampleptslist[cur_traj_index].num_samples++;
	cur_line = 0;
	cur_length = 0;
	//////////////////////////////////////////////////////////////////////////


	//////////Forward tracing
	pre_face = cur_face = triangle;
	globalp[0] = seed_p[0];
	globalp[1] = seed_p[1];
	flag = -1;

	int stop = 0;  //08/08/06

	/*02/27/07*/
	g_dt = 0;


	for(i = 0; i < /*NUMTRACINGTRIANGLE*/2000; i++)
	//for(i = 0; i < NUMTRACINGTRIANGLE; i++)
	{

		if(cur_face < 0)
		{
			break;
		}
		
		pre_face = cur_face;
        cur_face = TraceInATriangleAvoidClosed(cur_face, globalp, type,  loopdsep,
			Sample_interval, flag);

		////if it is close to a triangle containing the limit cycle, stop it 08/08/06
		int pre_limit;
		if(IsPreviousCycle_all(cur_face, pre_limit, 1-type))
		{
			//flag = 5; //it will finally form a loop if it tends to the limit cycle :) 08/08/06
			//break;
			if(stop == 0)  //record the first stop triangle
			{
				stop_triangle = cur_face;
				Face *tf = Object.flist[cur_face];
				icVector2 VP;
				VP.entry[0] = globalp[0] - Object.vlist[tf->verts[0]]->x;
				VP.entry[1] = globalp[1] - Object.vlist[tf->verts[0]]->y;

				stop_p_loc[0] = dot(VP, tf->LX);
				stop_p_loc[1] = dot(VP, tf->LY);
			}

			stop ++;
		}

		if(stop == 4)
		{
			flag = 5;
			break;
		}

		////We need to select the sampling points from current trajectory  3/9/06
		//GetSamplePtsWhenTracing(cur_traj_index, Sample_interval, cur_line, movetonext, cur_length,
		//   sampleptslist[cur_traj_index].samples, sampleptslist[cur_traj_index].num_samples,
		//   sampleptslist[cur_traj_index].MaxNumSamples);
		

	    if(flag == 3 || flag == 4 || flag == 5 || pre_face == cur_face ) 
		{
			break;
		}
	} 

	cur_traj_index++;
}

/*
Judge the side give a normalized vector using local coordinates
Return: 0--left side;  1--right side
*/
int GettheSide_loc(icVector2 orient, double basis[2], double p[2])
{
	icVector2 temp;
	temp.entry[0] = p[0] - basis[0];
	temp.entry[1] = p[1] - basis[1];

	normalize(temp);

	//cross product
	double result = orient.entry[0]*temp.entry[1] - orient.entry[1]*temp.entry[0];

	if(result > 0) //left side
		return 0;

	else if(result < 0) //right side
		return 1;

	else               //can not tell
		return -1;
}

/*
The distance should be related to the size of the triangle in the mesh 08/06/06
*/
bool GetSingularyID(double x, double y, int &singular_id)
{
	int i;
	for( i = 0; i < cur_singularity_index; i++)
	{
		if((fabs((double)(singularities[i].gcx - x)) < 0.005)
			&&(fabs((double)(singularities[i].gcy - y)) < 0.005))
		{
			singular_id = i;
			return true;
		}
	}
	if( i >= cur_singularity_index) return false;
}


extern void UpdateListInSingularity(int singID, int limitcycle);
extern void UpdateListInLimitCycle(int limitcycle, int saddleID);
extern bool TriangleSearch(int *acycle, int num_cycletriangles, int oneTriangle, int &position);

/*
Judge whether a point is close to a limit cycle.
Here type == 0 means it is outgoing sep, which should be connected to attracting component
whearas, type == 1 means it is incoming sep, which should be connected to repelling one
*/
bool CloseToLimitCycle(int saddle, double p[2], int triangle, int type, int sep_id)
{
	int i;
	int pos = -1;
	int cycleid;
	int Rep_t;

	//////08/08/06
	triangle = stop_triangle;
	p[0] = stop_p_loc[0];
	p[1] = stop_p_loc[1];

	for(i = 0; i < cur_limitcycle_index; i++)
	{
		if(limitcycles[i].type != 1-type)
			continue;

		if(TriangleSearch(limitcycles[i].cellcycle, limitcycles[i].num_triangles, triangle, pos))
		{
			cycleid = i;
			break;
		}
	}

	if(i >= cur_limitcycle_index)  //do not close to any limit cycle
		return false;

	////Store the length of the connection/separatrix between the saddle and the limit cycle 08/10/06
	int sep = singularities[saddle].separtices;
	switch(sep_id){
		case 1:
			sum_flow_length = separatrices[sep].length1;
			break;
		case 2:
			sum_flow_length = separatrices[sep].length2;
			break;
		case 3:
			sum_flow_length = separatrices[sep].length3;
			break;
		case 4:
			sum_flow_length = separatrices[sep].length4;
			break;
	}

	//build the connection 
	UpdateListInSingularity(saddle, cycleid);
	UpdateListInLimitCycle(cycleid, saddle);


	//mark the corresponding side of the limit cycle "connected"
	//////Rep_t = limitcycles[i].cellcycle[pos];

	//get the first line segment inside the triangle and get the start point and the orientation of the linesegment
	icVector2 orient;
	double basis[2], head[2];
	int stop = 0;

	//for(i = 0; i < limitcycles[cycleid].num_linesegs; i++)
	//{
	//	if(limitcycles[cycleid].closed_streamline[i].Triangle_ID == triangle)
	//	{
	//		orient.entry[0] = limitcycles[cycleid].closed_streamline[i].end[0]
	//			-limitcycles[cycleid].closed_streamline[i].start[0];
	//		orient.entry[1] = limitcycles[cycleid].closed_streamline[i].end[1]
	//			-limitcycles[cycleid].closed_streamline[i].start[1];

	//		basis[0] = limitcycles[cycleid].closed_streamline[i].start[0];
	//		basis[1] = limitcycles[cycleid].closed_streamline[i].start[1];
	//		break;

	//	}
	//}
	
	for(i = 0; i < limitcycles[cycleid].num_linesegs; i++)
	{
		if(limitcycles[cycleid].closed_streamline[i].Triangle_ID == triangle && stop == 0)
		{
			basis[0] = limitcycles[cycleid].closed_streamline[i].start[0];
			basis[1] = limitcycles[cycleid].closed_streamline[i].start[1];
			stop = 1;
		}

		if(stop == 1 && limitcycles[cycleid].closed_streamline[i].Triangle_ID != triangle)
		{
			head[0] = limitcycles[cycleid].closed_streamline[i-1].end[0];
			head[1] = limitcycles[cycleid].closed_streamline[i-1].end[1];

			orient.entry[0] = head[0] - basis[0];
			orient.entry[1] = head[1] - basis[1];
			break;
		}
	}

	//////
	normalize(orient);

	int side = GettheSide_loc(orient, basis, p);

	if(side == 0)
		limitcycles[cycleid].connected_l = 1;
	else if(side == 1)
		limitcycles[cycleid].connected_r = 1;

	return true;
}

/*
Trace inside a triangle
*/
int TraceInATriangleAvoidClosed(int &face_id, double globalp[2], int type, 
								double loopsep, double sample_interval, int &flag)
{
	int i;
	double alpha[3];
	double cur_point[2], pre_point[2];
	double vert0[2];
	icVector2 VP, globalv;

	Face *face = Object.flist[face_id];

	Face *pre_f = face;
	
	////Temporary curve point array

	CurvePoints *temp_point_list = (CurvePoints*) malloc(sizeof(CurvePoints) * 50);
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
    for(i = 0; i < 50; i++)
	{
		////1. calculate the barycentric coordinates for current point
		Get2DBarycentricFacters(face_id, cur_point[0], cur_point[1], alpha);

		////2. if current point is inside current triangle
		if( (alpha[0] >= 0 ||fabs(alpha[0]) <= 1e-10) && alpha[0] <= 1 
			&& (alpha[1] >= 0 || fabs(alpha[1]) <= 1e-10) && alpha[1] <= 1 
			&& (alpha[2] >= 0 || fabs(alpha[2]) <= 1e-10) && alpha[2] <= 1)
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

			//if(get_nextpt_2ndeuler(pre_point, cur_point, face_id, alpha, type))
			/*change to use other integration scheme 07/09/07*/
			//if(ToNextPoint(pre_point, cur_point, face_id, alpha, type))
			if(compute_next_pt(pre_point, cur_point, face_id, alpha, type))
			{

				globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

				globalp[0] = vert0[0] + globalv.entry[0];
				globalp[1] = vert0[1] + globalv.entry[1];

			}

			else{  
				////the curve reach a singularity

				flag = 3;

				////Store the record into global line segment array
                
				if(!StoreToGlobalList(temp_point_list, NumPoints))
				{
					////Not enough memory
					flag = 4;
					free(temp_point_list);
					return face_id;
				}

				//Get the singularity id, and mark it as "connected" 08/06/06
				globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

				globalp[0] = vert0[0] + globalv.entry[0];
				globalp[1] = vert0[1] + globalv.entry[1];
				int singID = -1;

				//if(GetSingularyID(globalp[0], globalp[1], singID))
				//{
				//	singularities[singID].connected = 1;
				//}

				singID = Object.flist[face_id]->singularity_index;
					singularities[singID].connected = 1;

				free(temp_point_list);

				return face_id;
			}

			////if it is close to a limit cycle, stop and build connection between them 08/06/06
			//int limitcycleid = -1;
			//if(CloseToLimitCycle(globalp, face_id, limitcycleid))

			////We may also need to compare the current point with the sampling point on itself!
			//if(CloseToCurSamplePt(globalp, face_id, sampleptslist[cur_traj_index].samples,
			//	sampleptslist[cur_traj_index].num_samples, loopsep, sample_interval)) //scale the separate distance
			//{
			//	if(!StoreToGlobalList(temp_point_list, NumPoints))
			//	{
			//		////Not enough memory
			//		flag = 4;
			//		free(temp_point_list);
			//		return face_id;
			//	}

			//	flag = 5; //form a loop! 4/18/06
			//	free(temp_point_list);
			//	return face_id;
			//}
		}

		////3. if the point is out of current triangle
		else{
			double t[2];

			int PassVertornot = 0;
   			int which_edge = -1;
         
			GetNextTriangle(face_id, pre_point, cur_point, t, type, PassVertornot, alpha);

			////update the global point here
			if(PassVertornot > 0)
			{
				////we should not directly use the vertex as next point!!
				////we may move a little bit along the VF direction, but make sure it is still inside
				////the new triangle

				//Vertex *PassedVert = Object.vlist[pre_f->verts[PassVertornot-1]];
				//globalp[0] = PassedVert->x;
				//globalp[1] = PassedVert->y;

				//cur_point[0] = pre_f->xy[PassVertornot-1][0];
				//cur_point[1] = pre_f->xy[PassVertornot-1][1];

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

	if(NumPoints > 0)
		StoreToGlobalList(temp_point_list, NumPoints);

	free(temp_point_list);

	return face_id;
}

/*
Initialize
*/
void AllocateSampleptsList()
{
	////allocate the sample point list
	sampleptslist = (SamplePtsList *)malloc(sizeof(SamplePtsList)*MaxNumTrajectories);
	int i;
	for(i = 0; i < MaxNumTrajectories; i++)
	{
		sampleptslist[i].samples = NULL;
		sampleptslist[i].num_samples = 0;
	}

	//for(i = 0; i < Object.nfaces; i++)
	//{
	//	Object.flist[i]->MaxSampNum = 5;
	//	Object.flist[i]->samplepts = (SampleListInTriangle *)malloc(sizeof(SampleListInTriangle)*5);
	//	Object.flist[i]->num_samplepts = 0;
	//}
}


void FinalizeSampleptsList()
{
	int i;
	for(i = 0; i < MaxNumTrajectories; i++)
	{
		if(sampleptslist[i].samples != NULL)
		free(sampleptslist[i].samples);
	}

	free(sampleptslist);
}


void InitSampleptsList()
{
	int i;
	for(i = 0; i < MaxNumTrajectories; i++)
	{
		sampleptslist[i].num_samples = 0;
	}
}
