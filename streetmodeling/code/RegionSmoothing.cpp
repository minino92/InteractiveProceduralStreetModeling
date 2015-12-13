////RegionSmoothing.cpp
#include "stdafx.h"
#include "RegionSmoothing.h"
#include "VFDataStructure.h"
//#include "LocalTracing.h"

#include "Numerical.h"




/*-------------------------------------------------------------------------------*/
////some global variables for region smoothing

int stack[200000];                  /////stack for flood search inside a closed region
int stack_top;

extern int SelectTriangleID;
extern Polygon3D Object;
extern double dmax;
extern Singularities *singularities;           //being captured singularites' list

////we need to links here to record the user selected region
////one is for visual show, the other is for smoothing
Point *point;                       ////we may initial it as 50 points, over 50, we can extend it
int Num_SmoothRegionpoints;                     ////Number of points that user selected
int MaxNumPoints;                   ////Maximum element in the point array, we can extend it if needed
int MaxNumEdges;
int MaxNumVerts;

////The other list is to store the boundary of the region based on underneath mesh
Vertex **boundaryverts;              ////mesh vertices on the boundary of user selected region
Vertex **regionverts;                ////mesh vertices inside user selected region
Edge **regionedge;                   ////mesh edges of user selected region

//icVector2 *beforesmooth, *aftersmooth;

extern icVector2 *pre_field, *post_field, *backup_field;

int Num_edges;
int Num_verts;                       ////number of inner vertices
int Num_boundaryverts;              ////number of boundary vertices

int InitialStateOn;

int Num_recursive;

////variables for numerical calculation routines
Vec_INT *ija_p;
Vec_DP *sa_p;


/*-------------------------------------------------------------------------------*/
////Implementation of routines for region smoothing
 ////judge wether a vertex is inside the region	
////Get the 2D Euler distance of two vertices
double GetLength(int ver1_id, int ver2_id)
{
	if(ver1_id == ver2_id)
		return 0.0;

	Vertex *ver1 = Object.vlist[ver1_id];
	Vertex *ver2 = Object.vlist[ver2_id];

	//icVector3 v;
	//v.entry[0] = ver1->x - ver2->x;
	//v.entry[1] = ver1->y - ver2->y;
	//v.entry[2] = ver1->z - ver2->z;

	//return (length(v));

	double v[3];
	v[0] = ver1->x - ver2->x;
	v[1] = ver1->y - ver2->y;
	v[2] = ver1->z - ver2->z;

	return (sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
}

void AddToEdgeList(Edge *e)
{
	////if the number of the edges not exceed the maximum capacity of the array, keep to add it and update the number

	if(Num_edges < MaxNumEdges)
	{
		regionedge[Num_edges] = e;

		Num_edges ++;
	}

	else{
		////we need to extend the array

		//int flag = 0;

		//regionedge = Extend_space(regionedge, (MaxNumEdges + 100), flag);

		//if(flag == 1)
		//	return;

		//boundaryverts = Extend_space(boundaryverts, (MaxNumEdges + 100), flag);

		//if(flag == 1)
		//	return;

		MaxNumEdges += 100;

		regionedge = (Edge**)realloc(regionedge, sizeof(Edge *) * MaxNumEdges);
		boundaryverts = (Vertex**) realloc(boundaryverts, sizeof(Vertex *) * (MaxNumEdges+1));

		regionedge[Num_edges] = e;
		Object.vlist[e->verts[0]]->OnBoundary = 1;
		Object.vlist[e->verts[1]]->OnBoundary = 1;

		Num_edges ++;

	}
}


bool InRegion(int vert_id)
{
	int i;
	double theta, sum_ang = 0;
	icVector2 v1, v2;
	Vertex *p = Object.vlist[vert_id];

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

	if( fabs(sum_ang) >= 2*M_PI - 1e-10)
		return true;
	else
		return false;
}

////Overload routine of InRegion to judge whether any input point is inside the region 
bool InRegion(double x, double y)
{
	int i;
	double theta, sum_ang = 0;
	icVector2 v1, v2;

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

	if( fabs(sum_ang) >= 2*M_PI - 1e-10)
		return true;
	else
		return false;
}
	
double GetDirectionalAngBetween2Vec(icVector2 v1, icVector2 v2)
{
	double ang1, ang2, dif_ang;

	ang1 = atan2(v1.entry[1], v1.entry[0]);
	ang2 = atan2(v2.entry[1], v2.entry[0]);

	dif_ang = ang2 - ang1;

	if(dif_ang < -M_PI)
		dif_ang += 2.*M_PI;

	if(dif_ang > M_PI)
		dif_ang -= 2.*M_PI;

	return dif_ang;

}

	
////Recursive search the edges on the boundary
////Recursive level maybe limit the program seriously!!!! 06/06/05

////Update Recursive search
////Some weird codes about using the variable boundaryverts[] may cause problem in future!!!!
void RecursiveDijistra(int s_id, int d_id, int cur_v_id, double dis_to_s)
{
	int i;
	double dis = 1e38;
	double curdistos, curdistod, seldistos;
	double temp_d = 0;

	Vertex *s, *d, *cur_v;
	s = Object.vlist[s_id];
	d = Object.vlist[d_id];
	cur_v = Object.vlist[cur_v_id];

	Vertex *adj_v = NULL, *selected_v = NULL;
	Edge *selected_e = NULL, *cur_edge = NULL;

	curdistos = curdistod = 0;
	
	////
	if(s_id == d_id)
		return;

	////Store the source vertex into the boundaryverts list

	boundaryverts[Num_boundaryverts] = cur_v;
	s->OnBoundary = 1;
	Num_boundaryverts ++;

	////Search all the adjacent vertices to find the one that is closest to 'd'
	for( i = 0 ; i < cur_v->Num_edge; i++)
	{
		cur_edge = cur_v->edges[i];

		if( Object.vlist[cur_edge->verts[0]] != cur_v)
			adj_v = Object.vlist[cur_edge->verts[0]];
		else
			adj_v = Object.vlist[cur_edge->verts[1]];
        
		curdistos = GetLength(adj_v->VertID, s->VertID);
		curdistod = GetLength(adj_v->VertID, d->VertID);

		temp_d = curdistos + curdistod;

		if(temp_d < dis && curdistos >= dis_to_s)
		{
			selected_v = adj_v;
			selected_e = cur_edge;
			dis = temp_d;
			seldistos = curdistos;
		}
	}
	
	////Add the selected edge to the boundary edges list
	AddToEdgeList(selected_e);

	if(selected_v == d) 
	{
		boundaryverts[Num_boundaryverts] = d;
		d->OnBoundary = 1;
		Num_boundaryverts++;
		return;
	}

	RecursiveDijistra(s_id, d_id, selected_v->VertID, seldistos);  ////if we do not reach the final vertex, we keep on searching
		
}

	


////Add the user selected point into the points list
void AddToPointList(double x, double y, int face)
{
	int vert;

	if(Num_SmoothRegionpoints >= MaxNumPoints - 1)
	{
		MaxNumPoints += 50;
		point = (Point*)realloc(point, sizeof(Point) * MaxNumPoints);
	}
	point[Num_SmoothRegionpoints].x = x;
	point[Num_SmoothRegionpoints].y = y;
	point[Num_SmoothRegionpoints].cellid = face;

	////Select the random vertex of the being chosen triangle 
	//vert = Object.flist[SelectTriangleID]->verts[0];
	//boundaryverts[Num_SmoothRegionpoints] = Object.vlist[vert];

	Num_SmoothRegionpoints ++;

}

	
void EndPointList()
{
	int i;
	
	if(Num_SmoothRegionpoints < 3)
	{
		MessageBox(NULL, "Not enough points! Please select at least 3 points to define the region!",\
			"error", MB_OK);
		Num_SmoothRegionpoints = 0;
		return;
	}

	/////call the recursive dijistra search

	for(i = 0; i < Num_SmoothRegionpoints; i++)
	{
		RecursiveDijistra(boundaryverts[i]->VertID, boundaryverts[(i+1)%Num_SmoothRegionpoints]->VertID, boundaryverts[i]->VertID, 0);
	}

	////Get the vertices inside the region
	FindOutInnerVerts();

}

int GetFirstInsideVertex()
{

	for(int i = 0; i < Num_edges; i++)
	{
		if(InRegion(Object.vlist[regionedge[i]->OppVerts[0]]->VertID))
		{
			return (Object.vlist[regionedge[i]->OppVerts[0]]->VertID);
		}

		else if(regionedge[i]->OppVerts[1] > 0 && InRegion(Object.vlist[regionedge[i]->OppVerts[1]]->VertID))
		{
			return (Object.vlist[regionedge[i]->OppVerts[1]]->VertID);
		}
	}
}

	
void FindOutInnerVerts()
{
	int firstVertID;

	firstVertID = GetFirstInsideVertex();

	regionverts[0] = Object.vlist[firstVertID];
	Num_verts++;

	////calling the seed search here
	//SeedSearch(firstVertID);
	NonRecursiveSearch(firstVertID);

}


void SeedSearch(int seed)
{
	Num_recursive++;

	Vertex *vert = Object.vlist[seed];
	Vertex *adj_v;

	Edge *temp_e;

	for(int i = 0; i < vert->Num_edge; i++)
	{
		temp_e = vert->edges[i];

		if(Object.vlist[temp_e->verts[0]] != vert)
			adj_v = Object.vlist[temp_e->verts[0]];
		else
			adj_v = Object.vlist[temp_e->verts[1]];

		////if the vertex is not on boundary and has not been marked, push it into the stack
		if(adj_v->InRegion != 1 && adj_v->OnBoundary != 1 && InRegion(adj_v->VertID))
		{
			stack_top++;
			stack[stack_top] = adj_v->VertID;
		}
	}

	if(stack_top == -1) return;

	if(stack_top >= 199999) 
	{
		MessageBox(NULL, "Out of the boundary of the stack", "error", MB_OK);
		stack_top = 0;
		return;
	}

	////Pop the vertex on the top of the stack and mark it as the vertex inside the region
	if(Num_verts >= MaxNumVerts - 1)
	{
		MaxNumVerts += 100;
		regionverts = (Vertex**)realloc(regionverts, sizeof(Vertex*) * MaxNumVerts);
	}
    regionverts[Num_verts] = Object.vlist[stack[stack_top]];
	regionverts[Num_verts]->RegionListID = Num_verts;
	regionverts[Num_verts]->InRegion = 1;
    Num_verts++;

	stack_top--;
	SeedSearch(stack[stack_top+1]);
        
}


void NonRecursiveSearch(int seed)
{
	int i, j;

	regionverts[Num_verts] = Object.vlist[seed];
	Num_verts++;

	Vertex *vert; 
	Vertex *adj_v;

	Edge *temp_e;

	////search all the inner vertices here
	i = 0;
	while(i < Num_verts)
	{
		vert = regionverts[i];  ////get current source vertex
		 
		////search all its adjacent vertices and judge whether it should be added to the inner vertices list or not
		for(j = 0; j < vert->Num_edge; j++)
		{
			temp_e = vert->edges[j];

			if(Object.vlist[temp_e->verts[0]] != vert)
				adj_v = Object.vlist[temp_e->verts[0]];
			else
				adj_v = Object.vlist[temp_e->verts[1]];

			if(adj_v == NULL)
			{
				MessageBox(NULL, "Danggling edge is found!", "Error mesh", MB_OK);
				exit(-1);
			}

			////if the vertex is not on boundary and has not been marked, push it into the stack
			if(adj_v->InRegion != 1 && adj_v->OnBoundary != 1 && InRegion(adj_v->VertID))
			{
				if(Num_verts >= MaxNumVerts - 1)
				{
					MaxNumVerts += 100;
					regionverts = (Vertex**)realloc(regionverts, sizeof(Vertex*) * MaxNumVerts);
					
				}

				regionverts[Num_verts] = adj_v;
				regionverts[Num_verts]->RegionListID = Num_verts;
				adj_v->InRegion = 1;
				Num_verts ++;
			}
		}

		i++;
	}
}

	
void InitRegionSmooth()
{
	Num_edges = 0;
	Num_SmoothRegionpoints = 0;
	Num_verts = 0;
	Num_boundaryverts = 0;

	Num_recursive = 0;

	stack_top = -1;

	for(int i = 0; i < Object.nverts; i++){
		Object.vlist[i]->InRegion = 0;
		Object.vlist[i]->OnBoundary = 0;
		Object.vlist[i]->RegionListID = -1;

	}
}


////build the sparse linear matrix in compact form
void BuildTheSparseLinearSystem(Vec_DP &tsa, Vec_INT &tija, Mat_DP &bound_v)
{
	int i, j; 
	Vertex *adj_v, *cur_v = NULL;
	Edge *adj_e;
	int *RegionIndex;
	int num_nonzeroarow = 0;
	int num_elements = 0;
	
	////fill the diagnal elements first
	for(i = 0; i < Num_verts; i++)
		//tsa[i] = 1 + SMOOTHSTEP;
		tsa[i] = 1. ;  ////11/06/05

	num_elements = Num_verts+1;

	////difficult part to store other non-zero, non-diagnal elements
	for(i = 0; i < Num_verts; i++)
	{
		cur_v = regionverts[i];
		RegionIndex = new int[cur_v->Num_edge];

		num_nonzeroarow = 0;

		for(j = 0; j < cur_v->Num_edge; j++)
		{
			adj_e = cur_v->edges[j];

			////get the adjacent vertex on the other side of current edge
			if(Object.vlist[adj_e->verts[0]] != cur_v)
				adj_v = Object.vlist[adj_e->verts[0]];
			else
				adj_v = Object.vlist[adj_e->verts[1]];

			////if the adjacent vertex is on boundary, put it to the right side
			if(adj_v->OnBoundary == 1 || adj_v->RegionListID < 0)
			{
				RegionIndex[j] = -1;

				//bound_v[i][0] += (SMOOTHSTEP/cur_v->Num_edge)*adj_v->vec.entry[0];
				//bound_v[i][1] += (SMOOTHSTEP/cur_v->Num_edge)*adj_v->vec.entry[1];
				bound_v[i][0] += (1./cur_v->Num_edge)*adj_v->vec.entry[0];  ////11/06/05
				bound_v[i][1] += (1./cur_v->Num_edge)*adj_v->vec.entry[1];
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

	tsa[Num_verts] = 0.;
	tija[Num_verts] = num_elements-1;
}

	
void BubbleSorting(int *a, int size)
{
    int i, j, temp;
    for ( i = 0; i < size; i++ )    // controls passes through the list 
    {
        for ( j = 0; j < size - 1; j++ )   // performs adjacent comparisons 
        {
            if(a[j] > a[j+1])      // determines if a swap should occur 
            {
                temp = a[j];       // swap is performed 
                a[j] = a[j+1];
                a[j+1] = temp;
            }
        }
    }
}


////Save the field before doing anything
void SavePreviousField()
{
	////Store the original field before performing smoothing operation
	for(int i = 0; i < Object.nverts; i++)
	{
		pre_field[i] = Object.vlist[i]->vec;
	}
}


////Save the field after doing something
void SavePostField()
{
	////Store the original field before performing smoothing operation
	for(int i = 0; i < Object.nverts; i++)
	{
		post_field[i] = Object.vlist[i]->vec;
	}
}


/////For one smoothing of the vector field inside user-defined region
void RegionSmooth()
{
	if(Num_verts < 2)
	{
		MessageBox(NULL, "Can not find enough inner vertices", "error", MB_OK);
		return;
	}

	int i;
	Vertex *cur_v;
    int NMAX=10*Num_verts;

    Vec_INT tempija(NMAX);
    Vec_DP tempsa(NMAX);

    DP err;
    Vec_DP b(Num_verts),bcmp(Num_verts),x(Num_verts);

	Mat_DP bound_vert(0.0, Num_verts, 2);  ////store the vectors that on the boundary vertices

    const int ITOL=1,ITMAX=100;
    const DP TOL=1.0e-10;
    int iter;

	//icVector2 tempv;

	BuildTheSparseLinearSystem(tempsa, tempija, bound_vert);

    Vec_INT ija(tempija[Num_verts]);
    Vec_DP sa(tempija[Num_verts]);

	for(i = 0; i < tempija[Num_verts]; i++)
	{
		ija[i] = tempija[i];
		sa[i] = tempsa[i];
	}

    ija_p = &ija;
	sa_p = &sa;
     

	////Smoothing the x component of the vectors
	for (i = 0; i < Num_verts; i++) 
	{
		cur_v = regionverts[i];
        x[i]=0.0;
		
		b[i]= bound_vert[i][0];
    }
    linbcg(b,x,ITOL,TOL,ITMAX,iter,err);

	////Store the results back to the vertices
	for (i = 0; i < Num_verts; i++) 
	{
		cur_v = regionverts[i];
		cur_v->vec.entry[0] = x[i];
    }

	////Smoothing the y component of the vectors
	for (i = 0; i < Num_verts; i++) 
	{
		cur_v = regionverts[i];
        x[i]=0.0;

		b[i]= bound_vert[i][1];
    }
    linbcg(b,x,ITOL,TOL,ITMAX,iter,err);

	////Store the results back to the vertices
	for (i = 0; i < Num_verts; i++) 
	{
		cur_v = regionverts[i];
		cur_v->vec.entry[1] = x[i];
		
		cur_v->vec_J = cur_v->vec;   /*we save the original vector before normalization 02/21/07*/
    }

	////For correct output, You need to normalize the result vector


    double r;

	for(i = 0; i < Num_verts; i++)
	{
		cur_v = regionverts[i];
	    r = length(cur_v->vec);
		r *= r;
					
		if (r < DistanceThreshold) 
		{
			r = DistanceThreshold;
			cur_v->vec *= dmax/r; 
		}

	    r = length(cur_v->vec);
		r *= r;

		if (r > dmax*dmax) { 
			r  = sqrt(r); 
			cur_v->vec *= dmax/r; 
		}
	}
	
	////Store the field after performing smoothing operation
	//for(i = 0; i < Num_verts; i++)
	//{
	//	aftersmooth[i] = regionverts[i]->vec;
	//}

}

void Undo()
{
	int i;

	//for(i = 0; i < Num_verts; i++)
	//{
	//	regionverts[i]->vec = beforesmooth[i];
	//}

	for(i = 0; i < Object.nverts; i++)
		Object.vlist[i]->vec = pre_field[i];
}


void Redo()
{
	int i;

	//for(i = 0; i < Num_verts; i++)
	//{
	//	regionverts[i]->vec = aftersmooth[i];
	//}

	for(i = 0; i < Object.nverts; i++)
		Object.vlist[i]->vec = post_field[i];
}


void AllocateVarforSmoothing()
{
	MaxNumPoints = 500;
	point = new Point[MaxNumPoints];

	Num_SmoothRegionpoints = 0;

	MaxNumEdges = 3;
	//MaxNumVerts = 3000;
	MaxNumVerts = /*Object.nverts+*/3;   ////set as the number of the vertices in the mesh for large region smoothing

	////The other list is to store the boundary of the region based on underneath mesh
	boundaryverts = (Vertex **)malloc(sizeof(Vertex *) * (MaxNumEdges+1));
	regionverts = (Vertex **)malloc(sizeof(Vertex *) * MaxNumVerts);
	regionedge = (Edge **)malloc(sizeof(Edge *) * MaxNumEdges);
	Num_edges = 0;
	Num_verts = 0;
	Num_boundaryverts = 0;
}

void FinalizeSmoothing()
{
	try{
		free(boundaryverts);
		free(regionverts);
		free(regionedge);

		free(backup_field);
		free(pre_field);
		free(post_field);
	}
	catch(CMemoryException *e)
	{
		char *s0;
		//e->GetErrorMessage(s0, 1);
		CString temp;
		//temp.Format("%s", s0);
		CString s = "Error: " + temp;
		const char *err = s;
		MessageBox(NULL, err, "Error", MB_OK);
		e->Delete();
		return;
	}
	catch(CException *e)
	{
	}
	//catch()
	//{
	//}

//_except_handler(
}


void ResetSmoothVars()
{
	Num_edges = 0;
	Num_verts = 0;
	Num_boundaryverts = 0;
	
	for(int i = 0; i < Object.nverts; i++){
		Object.vlist[i]->InRegion = 0;
		Object.vlist[i]->OnBoundary = 0;
		Object.vlist[i]->RegionListID = -1;

	}
}

