////RegionSmoothing.cpp
#include "stdafx.h"
#include "RegionSmoothing.h"
#include "VFDataStructure.h"
#include "LocalTracing.h"
#include "tensoranalysis.h"

#include "Numerical.h"

extern int stack[200000];                  /////stack for flood search inside a closed region
extern int stack_top;
extern QuadMesh *quadmesh;

extern icVector2 tenline_dir_global;

int *region_quadverts=NULL;                ////mesh vertices inside user selected region
int nregion_quadverts = 0;
int curMaxRegionQuadVerts;

int *boundarycells=NULL;
int nboundarycells=0;

int firstVertID_quad = -1;

extern Point *point;                       ////we may initial it as 50 points, over 50, we can extend it
extern int Num_SmoothRegionpoints;                     ////Number of points that user selected


////variables for numerical calculation routines
extern Vec_INT *ija_p;
extern Vec_DP *sa_p;

extern double GetDirectionalAngBetween2Vec(icVector2 v1, icVector2 v2);


void alloc_quad_regionsmooth()
{
	region_quadverts=(int *)malloc(sizeof(int)*quadmesh->nverts);
	curMaxRegionQuadVerts=quadmesh->nverts;
	nregion_quadverts = 0;

	boundarycells=(int*)malloc(sizeof(int)*quadmesh->nfaces);
	nboundarycells = 0;
}

void init_quad_regionsmooth()
{
	nregion_quadverts = 0;

	stack_top = -1;
	int i;

	for(i = 0; i < quadmesh->nverts; i++){
		quadmesh->quad_verts[i]->InRegion = false;
		quadmesh->quad_verts[i]->OnBoundary = false;
		quadmesh->quad_verts[i]->RegionListID = -1;
	}

	for(i=0; i<quadmesh->nfaces; i++)
		quadmesh->quadcells[i]->OnBoundary = false;
}

int get_firstInsideVertex()
{
	/*check the vertices of the cell strip and find out the one that is not marked as 
	"OnBoundary"*/
	int i, j;
	QuadCell *face;
	QuadVertex *v;
	for(i=0; i<nboundarycells; i++)
	{
		face = quadmesh->quadcells[boundarycells[i]];

		for(j=0; j<face->nverts; j++)
		{
			if(!quadmesh->quad_verts[face->verts[j]]->OnBoundary)
				return face->verts[j];
		}
	}

	return -1;
}


/*    The following two routines try to judge whether a specified vertex of the quad mesh
      or a random 2D point falls in a user specified region 
*/

bool is_inregion(int vert_id)
{
	int i;
	double theta, sum_ang = 0;
	icVector2 v1, v2;
	QuadVertex *p = quadmesh->quad_verts[vert_id];

	////Calculate the sum of the angle
	for( i = 0; i < Num_SmoothRegionpoints; i++)
	{
		v1.entry[0] = point[i].x - p->x;
		v1.entry[1] = point[i].y - p->y;

		v2.entry[0] = point[(i+1)%Num_SmoothRegionpoints].x - p->x;
		v2.entry[1] = point[(i+1)%Num_SmoothRegionpoints].y - p->y;

		normalize(v1);
		normalize(v2);

		////The following stuff may slow down the performance of smoothing

		theta = GetDirectionalAngBetween2Vec(v1, v2);
		sum_ang += theta;
	}

	if( fabs(sum_ang) >= 2*M_PI - 1e-8)
		return true;
	else
		return false;
}


bool is_inregion(double x, double y)
{
	int i;
	double theta, sum_ang = 0;
	icVector2 v1, v2;
	//QuadVertex *p = quadmesh->quad_verts[vert_id];

	////Calculate the sum of the angle
	for( i = 0; i < Num_SmoothRegionpoints; i++)
	{
		v1.entry[0] = point[i].x - x;
		v1.entry[1] = point[i].y - y;

		v2.entry[0] = point[(i+1)%Num_SmoothRegionpoints].x - x;
		v2.entry[1] = point[(i+1)%Num_SmoothRegionpoints].y - y;

		normalize(v1);
		normalize(v2);

		////The following stuff may slow down the performance of smoothing

		theta = GetDirectionalAngBetween2Vec(v1, v2);
		sum_ang += theta;
	}

	if( fabs(sum_ang) >= 2*M_PI - 1.e-8)
		return true;
	else
		return false;
}


bool is_inregion_appro(double x, double y, double threshold)
{
	int i;
	double theta, sum_ang = 0;
	icVector2 v1, v2;
	//QuadVertex *p = quadmesh->quad_verts[vert_id];

	////Calculate the sum of the angle
	for( i = 0; i < Num_SmoothRegionpoints; i++)
	{
		v1.entry[0] = point[i].x - x;
		v1.entry[1] = point[i].y - y;

		v2.entry[0] = point[(i+1)%Num_SmoothRegionpoints].x - x;
		v2.entry[1] = point[(i+1)%Num_SmoothRegionpoints].y - y;

		normalize(v1);
		normalize(v2);

		////The following stuff may slow down the performance of smoothing

		theta = GetDirectionalAngBetween2Vec(v1, v2);
		sum_ang += theta;
	}

	if( fabs(sum_ang) >= 2*M_PI - threshold)
		return true;
	else
		return false;
}




/*
We assume that we have found all the boundary cells
*/
void mark_boundVerts()
{
	int i, j;
	QuadCell *face;
	QuadVertex *v;
	for(i=0; i<nboundarycells; i++)
	{
		face = quadmesh->quadcells[boundarycells[i]];

		for(j=0; j<face->nverts; j++)
		{
			v = quadmesh->quad_verts[face->verts[j]];
			if(!is_inregion(face->verts[j]))
				v->OnBoundary=true;
			else
				firstVertID_quad=face->verts[j];
		}
	}
}

void search_innerVerts_nonrecursive(int seed)
{
	int i, j;

	region_quadverts[nregion_quadverts] = quadmesh->quad_verts[seed]->index;
	nregion_quadverts++;

	QuadVertex *vert; 
	QuadVertex *adj_v;

	QuadEdge *temp_e;

	////search all the inner vertices here
	i = 0;
	while(i < nregion_quadverts)
	{
		vert = quadmesh->quad_verts[region_quadverts[i]];  ////get current source vertex
		 
		////search all its adjacent vertices and judge whether it should be added to the inner vertices list or not
		for(j = 0; j < vert->Num_edge; j++)
		{
			temp_e = vert->edges[j];

			if(quadmesh->quad_verts[temp_e->verts[0]] != vert)
				adj_v = quadmesh->quad_verts[temp_e->verts[0]];
			else
				adj_v = quadmesh->quad_verts[temp_e->verts[1]];

			if(adj_v == NULL)
			{
				MessageBox(NULL, "Danggling edge is found!", "Error mesh", MB_OK);
				exit(-1);
			}

			////if the vertex is not on boundary and has not been marked, push it into the stack
			if(!adj_v->InRegion && !adj_v->OnBoundary && is_inregion(adj_v->index))
			{
				if(nregion_quadverts >= curMaxRegionQuadVerts)
				{
					curMaxRegionQuadVerts += 100;
					region_quadverts = (int*)realloc(region_quadverts, 
						sizeof(int) * curMaxRegionQuadVerts);
					
				}

				region_quadverts[nregion_quadverts] = adj_v->index;
				quadmesh->quad_verts[region_quadverts[nregion_quadverts]]->RegionListID 
					= nregion_quadverts;
				adj_v->InRegion = true;
				nregion_quadverts ++;
			}
		}
		i++;
	}
}


void find_innerVerts()
{
	//int firstVertID;
	//firstVertID = get_firstInsideVertex();

	/*We assume that firstVertID_quad >=0 now*/

	region_quadverts[0] = quadmesh->quad_verts[firstVertID_quad]->index;
	nregion_quadverts=1;

	////calling the seed search here
	search_innerVerts_nonrecursive(firstVertID_quad);

}

void find_boundarycells_oneline(double start[2], int start_cell, double end[2],
								int end_cell, int *celllist, int &ncells)
{
	double pre_p[2]={start[0], start[1]};
	double cur_p[2]={end[0], end[1]};
	int i;
	int cur_cell = start_cell;
	QuadCell *face;
	icVector2 linedir;
	linedir.set((end[0]-start[0]),(end[1]-start[1]));
	normalize(linedir);
	tenline_dir_global = linedir;

	/*add the first cell*/
	celllist[ncells]=cur_cell;
	ncells++;

	icVector2 t_major[4];

	int counter=0;

	while(cur_cell != end_cell && counter<50)
	{
		if(!is_in_cell(cur_cell, end[0], end[1]))
		{
			/*find the next cell the curve will enter*/

			/*we first store the major vectors of the vertices
			of current cell*/
			face = quadmesh->quadcells[cur_cell];
			for(i=0; i<4; i++)
				t_major[i]=quadmesh->quad_verts[face->verts[i]]->major;

			/*replace them with the line segment direction*/
			for(i=0; i<4; i++)
				quadmesh->quad_verts[face->verts[i]]->major=linedir;

			int passvertornot = 0;
			//get_next_cell(cur_cell, pre_p, cur_p, passvertornot, 0);
			get_next_cell_2(cur_cell, pre_p, cur_p, passvertornot, 0);

			if(cur_cell<0 || cur_cell>quadmesh->nfaces)
				return;

			/*add to the cell list*/
			if(!quadmesh->quadcells[cur_cell]->OnBoundary)
			{
				celllist[ncells]=cur_cell;
				ncells++;

				quadmesh->quadcells[cur_cell]->OnBoundary=true;
			}

			if(passvertornot == 0)
			{
				pre_p[0] = cur_p[0];
				pre_p[1] = cur_p[1];
			}
			cur_p[0] = end[0];
			cur_p[1] = end[1];

			/*store back the original vectors*/
			for(i=0; i<4; i++)
				quadmesh->quad_verts[face->verts[i]]->major=t_major[i];
		}

		counter++;
	}

	if(counter>=50)
	{
		int test=0;
	}
}

extern int get_cellID_givencoords(double, double);

void find_boundarycells()
{
	int i;
	double start[2], end[2];
	nboundarycells=0;
	for(i=0; i<Num_SmoothRegionpoints; i++)
	{
		start[0] = point[i].x;
		start[1] = point[i].y;
		point[i].cellid=get_cellID_givencoords(start[0], start[1]);
		end[0] = point[(i+1)%Num_SmoothRegionpoints].x;
		end[1] = point[(i+1)%Num_SmoothRegionpoints].y;
		point[(i+1)%Num_SmoothRegionpoints].cellid=get_cellID_givencoords(end[0], end[1]);

		find_boundarycells_oneline(start, point[i].cellid, end, point[(i+1)%Num_SmoothRegionpoints].cellid,
            boundarycells, nboundarycells);
	}

}





/*use the Conjugate gradient solver from Numerical Recipe to solve 
the laplacian smoothing*/

void construct_sparseSys_quad(Vec_DP &tsa, Vec_INT &tija, Mat_DP &bound_v)
{
	int i, j; 
	QuadVertex *adj_v, *cur_v = NULL;
	QuadEdge *adj_e;
	int *RegionIndex;
	int num_nonzeroarow = 0;
	int num_elements = 0;
	
	////fill the diagnal elements first
	for(i = 0; i < nregion_quadverts; i++)
		//tsa[i] = 1 + SMOOTHSTEP;
		tsa[i] = 1. ;  ////11/06/05

	num_elements = nregion_quadverts+1;

	////difficult part to store other non-zero, non-diagnal elements
	for(i = 0; i < nregion_quadverts; i++)
	{
		cur_v = quadmesh->quad_verts[region_quadverts[i]];
		RegionIndex = new int[cur_v->Num_edge];

		num_nonzeroarow = 0;

		for(j = 0; j < cur_v->Num_edge; j++)
		{
			adj_e = cur_v->edges[j];

			////get the adjacent vertex on the other side of current edge
			if(quadmesh->quad_verts[adj_e->verts[0]] != cur_v)
				adj_v = quadmesh->quad_verts[adj_e->verts[0]];
			else
				adj_v = quadmesh->quad_verts[adj_e->verts[1]];

			////if the adjacent vertex is on boundary, put it to the right side
			if(adj_v->OnBoundary || adj_v->RegionListID < 0)
			{
				RegionIndex[j] = -1;

				bound_v[i][0] += (1./cur_v->Num_edge)*adj_v->Jacobian.entry[0][0];  
				bound_v[i][1] += (1./cur_v->Num_edge)*adj_v->Jacobian.entry[0][1];
				bound_v[i][2] += (1./cur_v->Num_edge)*adj_v->Jacobian.entry[1][0];  
				bound_v[i][3] += (1./cur_v->Num_edge)*adj_v->Jacobian.entry[1][1];
			}
			else
			{
				RegionIndex[j] = adj_v->RegionListID;
				num_nonzeroarow++;   ////add one more non-zero, non-diagnal element in current row
			}
		}

		////sorting the regionindex array
		BubbleSorting(RegionIndex, cur_v->Num_edge);

		////move all non -1 element in the 'RegionIndex' to the front
		if(num_nonzeroarow < cur_v->Num_edge)
		{
			int firstnonfuyielement = cur_v->Num_edge - num_nonzeroarow;

			for(j = 0; j < num_nonzeroarow; j++)
			{
				RegionIndex[j] = RegionIndex[firstnonfuyielement + j];
			}
		}

		////Add elements to the corresponding positions of sa and ija array
		tija[i] = num_elements;
		for(j = 0; j < num_nonzeroarow; j++)
		{
			//tsa[num_elements + j] = -SMOOTHSTEP/cur_v->Num_edge;
			tsa[num_elements + j] = -1./cur_v->Num_edge;  ////11/06/05
			tija[num_elements + j] = RegionIndex[j];
		}

		num_elements += num_nonzeroarow;

		delete [] RegionIndex;
	}

	tsa[nregion_quadverts] = 0.;
	tija[nregion_quadverts] = num_elements-1;
}


/////For one smoothing of the vector field inside user-defined region
void smooth_Jac_quadregion()
{
	if(nregion_quadverts < 2)
	{
		MessageBox(NULL, "Can not find enough inner vertices", "error", MB_OK);
		return;
	}

	int i;
	QuadVertex *cur_v;
    int NMAX=10*nregion_quadverts;

    Vec_INT tempija(NMAX);
    Vec_DP tempsa(NMAX);

    DP err;
    Vec_DP b(nregion_quadverts),bcmp(nregion_quadverts),x(nregion_quadverts);

	Mat_DP bound_vert(0.0, nregion_quadverts, 4);  ////store the vectors that on the boundary vertices

    const int ITOL=1,ITMAX=200;
    const DP TOL=1.0e-10;
    int iter;

	//icVector2 tempv;

	construct_sparseSys_quad(tempsa, tempija, bound_vert);

    Vec_INT ija(tempija[nregion_quadverts]);
    Vec_DP sa(tempija[nregion_quadverts]);

	for(i = 0; i < tempija[nregion_quadverts]; i++)
	{
		ija[i] = tempija[i];
		sa[i] = tempsa[i];
	}

    ija_p = &ija;
	sa_p = &sa;
     

	////Smoothing the 00 component of the Jacobian
	for (i = 0; i < nregion_quadverts; i++) 
	{
		cur_v = quadmesh->quad_verts[region_quadverts[i]];
        x[i]=0.0;
		
		b[i]= bound_vert[i][0];
    }
    linbcg(b,x,ITOL,TOL,ITMAX,iter,err);

	////Store the results back to the vertices
	for (i = 0; i < nregion_quadverts; i++) 
	{
		cur_v = quadmesh->quad_verts[region_quadverts[i]];
		if(cur_v->inland) /*consider the geograph map*/
		cur_v->Jacobian.entry[0][0] = x[i];
    }

	////Smoothing the 01 component of the Jacobian
	for (i = 0; i < nregion_quadverts; i++) 
	{
		cur_v = quadmesh->quad_verts[region_quadverts[i]];
        x[i]=0.0;

		b[i]= bound_vert[i][1];
    }
    linbcg(b,x,ITOL,TOL,ITMAX,iter,err);

	////Store the results back to the vertices
	for (i = 0; i < nregion_quadverts; i++) 
	{
		cur_v = quadmesh->quad_verts[region_quadverts[i]];
		if(cur_v->inland) /*consider the geograph map*/
		cur_v->Jacobian.entry[0][1] = x[i];
    }

	////Smoothing the 10 component of the Jacobian
	for (i = 0; i < nregion_quadverts; i++) 
	{
		cur_v = quadmesh->quad_verts[region_quadverts[i]];
        x[i]=0.0;

		b[i]= bound_vert[i][2];
    }
    linbcg(b,x,ITOL,TOL,ITMAX,iter,err);

	////Store the results back to the vertices
	for (i = 0; i < nregion_quadverts; i++) 
	{
		cur_v = quadmesh->quad_verts[region_quadverts[i]];
		if(cur_v->inland) /*consider the geograph map*/
		cur_v->Jacobian.entry[1][0] = x[i];
    }

	////Smoothing the 11 component of the Jacobian
	for (i = 0; i < nregion_quadverts; i++) 
	{
		cur_v = quadmesh->quad_verts[region_quadverts[i]];
        x[i]=0.0;

		b[i]= bound_vert[i][3];
    }
    linbcg(b,x,ITOL,TOL,ITMAX,iter,err);

	////Store the results back to the vertices
	for (i = 0; i < nregion_quadverts; i++) 
	{
		cur_v = quadmesh->quad_verts[region_quadverts[i]];
		if(cur_v->inland) /*consider the geograph map*/
		cur_v->Jacobian.entry[1][1] = x[i];
    }
}


#include ".\glview.h"

/*test: visualize the inner vertices*/
void display_innerverts()
{
	int i;
	QuadVertex *v;
	glColor3f(1, 1, 0);
	glPointSize(3.);
	glBegin(GL_POINTS);
	for(i=0; i<nregion_quadverts; i++)
	{
		v = quadmesh->quad_verts[region_quadverts[i]];
		glVertex2f(v->x, v->y);
	}
	glEnd();
	glPointSize(1.);
}