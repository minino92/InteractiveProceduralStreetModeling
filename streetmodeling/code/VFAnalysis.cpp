#include "stdafx.h"

#include "VFAnalysis.h"

#include "VFDataStructure.h"

int FindTriangleIndex;
int *MarkTriangleID = NULL;
	
double r1, r2, i1, i2;  //Make these Jaccobian Matrix components global !!!!!
double a, b, c, d, e, f; //it is useless here
icMatrix3x3 FieldMatrix;
double x_cp, y_cp;      //use to get the singularity coordinates 3/25/06

extern int MaxNumSingularities;                 //Maximum number of being captured singularities
extern int MaxNumTrajectories;                  //Maximum number of possible trajectories
                                         //(it should be flexible for future pen-and-ink sketch)
extern int MaxNumLinesegsPerTraj;               //Maximum number of line segments for each trajectory

extern Singularities *singularities;           //being captured singularites' list
extern int cur_singularity_index;
extern LineSeg **trajectories;                 //trajectories' list
extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory

extern Polygon3D Object;
extern bool UseNormalizedVF;

/*--------------------------------------------------------------------------*/
////Operations for vector field analysis

/**************************************************
Capture singularities using new method
**************************************************/

void CaptureSing(void)
{
	unsigned int i, j;
	Face *face;
	int *verts;
	icVector2 vec[3];      //vectors at three vertices
    double  theta[3];      //for storing the angles between two vector for Gauss circle judgement

	double  vec_ang[3];  //the angle for each vector under the polar frame of current triangle
	double  ang_sum;

	////Initialize
	FindTriangleIndex = 0;

	if(MarkTriangleID == NULL)
		MarkTriangleID = (int*) malloc(sizeof(int) * MaxNumSingularities); //default range is 200

	////These flags are used for pair cancellation
	////Initialize all the flags!!!
	for(i = 0; i < Object.nfaces; i++)
	{
		Object.flist[i]->contain_singularity = 0;
		Object.flist[i]->singularity_index = -1;
	}

	////initialize the singularities list
	cur_singularity_index = 0;

	////Calculate the Gaussian angle
	for (i = 0; i < Object.nfaces; i++) {
		face = Object.flist[i];
		verts = face->verts;

		ang_sum = 0;

		//For each triangle, calculate the vector for all vertices
		for (j=0; j < face->nverts; j++) {
			////using the vectors stored in the vertices
			vec[j] =  Object.vlist[verts[j]]->vec;
			//normalize(vec[j]);

			vec_ang[j] = atan2(vec[j].entry[1], vec[j].entry[0]);

			////
			if(vec_ang[j] < 0) vec_ang[j] += 2 * M_PI;
		}

		for(j = 0; j < face->nverts; j++)
		{
			if(j == 2)
			{
				theta[j] = vec_ang[0] - vec_ang[2];
			}
			else{
				theta[j] = vec_ang[j+1] - vec_ang[j];
			}

			if( theta[j] < -M_PI)
				theta[j] += 2 * M_PI;
			
			if( theta[j] > M_PI)
				theta[j] -= 2 * M_PI;


			ang_sum += theta[j];
		}

		if(fabs(ang_sum) >= (2 * M_PI - 0.5))
		{
			//The triangle must have singularities inside, mark it as yellow color
			//Still need to judge whether it is one of current singularities or not
			MarkTriangleID[FindTriangleIndex] = i;
			FindTriangleIndex ++;

			face->contain_singularity = 1;         ////Mark current triangle as one containing singularity
			//face->singularity_index = FindTriangleIndex-1; ////which singularity it contains

			if(FindTriangleIndex >= MaxNumSingularities - 1)
			{
				MaxNumSingularities += 50;
				MarkTriangleID = (int*) realloc(MarkTriangleID, sizeof(int) * MaxNumSingularities);
				singularities = (Singularities*) realloc(singularities, sizeof(Singularities) * MaxNumSingularities);
			}
		}
	}

	////Calculate the coordinates of being found singularities
	if(FindTriangleIndex > 0)
	{
		if(UseNormalizedVF)
		    ComputeSing();
		else
			compute_fixedpt_info();
	}
}

/***************************************************************
Compute the accurate coordinates after capture the ID of 
the triangles that contain singularities using new method 1/16
***************************************************************/
void ComputeSing(void)
{
	int i, j;
	Face *face;
	int *verts;
	double a, b, c, a1, b1; //
	double x[3], y[3]; //
	double vx[3], vy[3];//
	double vdx, vdy; //
	double vzero_x, vzero_y;
    double temp_ratio = 0;//

	////Initialize the unknown singularities link
	vzero_x = vzero_y = 0;

	for(i = 0; i < FindTriangleIndex; i++)
	{
		//For each being captured triangle, compute the coordinates of singularities inside it
		face = Object.flist[MarkTriangleID[i]];
		verts = face->verts;
			
		//For each triangle, calculate the vector for all vertices
		for (j=0; j<face->nverts; j++) {
            x[j] = Object.vlist[verts[j]]->x;
			y[j] = Object.vlist[verts[j]]->y;

			////using the vectors stored in the vertices
			
			/* Use normalized vector field*/
			vx[j] = Object.vlist[verts[j]]->vec.entry[0];
			vy[j] = Object.vlist[verts[j]]->vec.entry[1];

			/* Use non-normalized vector field*/
			//vx[j] = Object.vlist[verts[j]]->vec_J.entry[0];
			//vy[j] = Object.vlist[verts[j]]->vec_J.entry[1];

			////Testing using local vectors here
			//vx[j] = face->direct_vec[j].entry[0];
			//vy[j] = face->direct_vec[j].entry[1];

			if(vx[j] == 0 && vy[j] == 0)
			{
				//add this vertex into the unknow singularities link
				////Add the point into the unkown singularities link
				singularities[cur_singularity_index].gcx = x[j];
				singularities[cur_singularity_index].gcy = y[j];

				singularities[cur_singularity_index].type = GetSingType(MarkTriangleID[i]); //Judge the type of this singularity

				////Store the Jacobian matrix of the singularity
		
				//double M[][2] = {FieldMatrix.entry[0][0],FieldMatrix.entry[0][1],
				//	FieldMatrix.entry[1][0],FieldMatrix.entry[1][1]};
				singularities[cur_singularity_index].Jacobian.set(FieldMatrix); //store the field matrix
				
				singularities[cur_singularity_index].Triangle_ID = MarkTriangleID[i];

				//singularities[cur_singularity_index].node_index = -1;           //10/16/05
			//face->singularity_index = FindTriangleIndex-1; ////which singularity it contains
				Object.flist[MarkTriangleID[i]]->singularity_index = cur_singularity_index;

				if(singularities[cur_singularity_index].type == SADDLE)
				{
					////Calculate its major direction here 1/27/05
					////Here, r1 is always larger than r2
					singularities[cur_singularity_index].outgoing.entry[0]\
						= singularities[cur_singularity_index].Jacobian.entry[0][1]; //outgoing direction
					singularities[cur_singularity_index].outgoing.entry[1]\
						= r1 - singularities[cur_singularity_index].Jacobian.entry[0][0];

					singularities[cur_singularity_index].incoming.entry[0]\
						= singularities[cur_singularity_index].Jacobian.entry[0][1]; //outgoing direction
					singularities[cur_singularity_index].incoming.entry[1]\
						= r2 - singularities[cur_singularity_index].Jacobian.entry[0][0];
					//singularities[cur_singularity_index].incoming.entry[0]\
					//	= r2 - singularities[cur_singularity_index].Jacobian.entry[1][1];
					//singularities[cur_singularity_index].incoming.entry[1]\
					//	= singularities[cur_singularity_index].Jacobian.entry[1][0]; //incoming direction

					normalize(singularities[cur_singularity_index].outgoing);
					normalize(singularities[cur_singularity_index].incoming);

				}

				singularities[cur_singularity_index].connected_limitcycles = NULL;
				singularities[cur_singularity_index].num_connected_limitcycles = 0;

				singularities[cur_singularity_index].connected = 0;

				cur_singularity_index++;

			}
		}

		/////Calculate the coordinates of singularities here

		/////The following codes using bilinear method to get the center of new critical points
		//if(vx[2]!=0){
  //          temp_ratio = vy[2]/vx[2];
		//}
		//else{
		//	///Switch the index of vectors to make sure it is not divided by zero
		//}

		//////Calculate the a1 first according to the following equations
		//// vdx = a1*vx[0] + b1* vx[1];
		//// vdy = a1*vy[0] + b1* vy[1];
		//// a1 + b1 = 1;
		//// vdy/vdx = temp_ratio;

  //      a1 = (vy[1] - temp_ratio * vx[1])/
		//	(vy[1] - vy[0] + temp_ratio*vx[0] - temp_ratio*vx[1]);
		//b1 = 1 - a1;

		/////Getting vdx, vdy using a1 and b1
		//vdx = a1 * vx[0] + b1 * vx[1];
		//vdy = a1 * vy[0] + b1 * vy[1];

		/////Calculate c using following equation
		////(1 - c)*vdx + c*vx[2] = 0
  //      c = vdx / (vdx - vx[2]); //need to judge whether vdx == vx[2] or not!!!!

		/////Get the a, b
		//a = (1 - c)*a1;
		//b = (1 - c)*b1;

		/////Get the coordinates of this singularities
		//if((a>=0)&&(a<=1)&&(b>=0)&&(b<=1)&&(c>=0)&&(c<=1)){
		//	vzero_x = a * x[0] + b * x[1] + c * x[2];  ////using the barycentric coefficients to calculate 
		//	vzero_y = a * y[0] + b * y[1] + c * y[2];  ////the center of the singularity
		//}

		///Add to the being captured singularities list
		singularities[cur_singularity_index].type = GetSingType(MarkTriangleID[i]); //Judge the type of this singularity

		////3/25/06
		//if(singularities[cur_singularity_index].type == SINK)
		//{
			singularities[cur_singularity_index].gcx = x_cp; ////set the center coordinates
			singularities[cur_singularity_index].gcy = y_cp;
		//}

		//else
		//{
		//singularities[cur_singularity_index].gcx = vzero_x; ////set the center coordinates
		//singularities[cur_singularity_index].gcy = vzero_y;
		//}

		singularities[cur_singularity_index].Jacobian.set(FieldMatrix); //store the field matrix
			
		singularities[cur_singularity_index].Triangle_ID = MarkTriangleID[i];

		//singularities[cur_singularity_index].node_index = -1;           //10/16/05
		Object.flist[MarkTriangleID[i]]->singularity_index = cur_singularity_index;
		
		if(singularities[cur_singularity_index].type == SADDLE)
		{
			////Here, r1 is always larger than r2
			
			//singularities[cur_singularity_index].outgoing.entry[0]\
			//	= singularities[cur_singularity_index].Jacobian.entry[0][1]; //outgoing direction
			//singularities[cur_singularity_index].outgoing.entry[1]\
			//	= r1 - singularities[cur_singularity_index].Jacobian.entry[0][0];

			//singularities[cur_singularity_index].incoming.entry[0]\
			//	= r2 - singularities[cur_singularity_index].Jacobian.entry[1][1];
			//singularities[cur_singularity_index].incoming.entry[1]\
			//	= singularities[cur_singularity_index].Jacobian.entry[1][0]; //incoming direction

			
			singularities[cur_singularity_index].outgoing.entry[0]\
				= -singularities[cur_singularity_index].Jacobian.entry[0][1]; //outgoing direction
			singularities[cur_singularity_index].outgoing.entry[1]\
				= singularities[cur_singularity_index].Jacobian.entry[0][0]- r1;

			singularities[cur_singularity_index].incoming.entry[0]\
				= -singularities[cur_singularity_index].Jacobian.entry[0][1]; //incoming direction
			singularities[cur_singularity_index].incoming.entry[1]\
				= singularities[cur_singularity_index].Jacobian.entry[0][0]- r2;
			//singularities[cur_singularity_index].incoming.entry[0]\
			//	= r2 - singularities[cur_singularity_index].Jacobian.entry[1][1];
			//singularities[cur_singularity_index].incoming.entry[1]\
			//	= singularities[cur_singularity_index].Jacobian.entry[1][0]; //incoming direction



			normalize(singularities[cur_singularity_index].outgoing);
			normalize(singularities[cur_singularity_index].incoming);
		}
		singularities[cur_singularity_index].connected_limitcycles = NULL;
		singularities[cur_singularity_index].num_connected_limitcycles = 0;
		singularities[cur_singularity_index].connected = 0;
		cur_singularity_index++;

	}
}

/*
This routine uses normalized vector field to find the type of the 
fixed point
*/
int GetSingType(int MarkTriangleID)
{
	Face *face;
	int *verts;

	double x[3], y[3]; //
	double vx[3], vy[3];//

	int i, j;

	///These information has been obtained before
	///Can we reuse them and enhance the performance 
    face = Object.flist[MarkTriangleID];
	verts = face->verts;

	//For each triangle, calculate the vector for all vertices
	for (j=0; j<face->nverts; j++) {
        x[j] = Object.vlist[verts[j]]->x;
		y[j] = Object.vlist[verts[j]]->y;

		/*Use normalized vector field*/
		vx[j] = Object.vlist[verts[j]]->vec.entry[0];
		vy[j] = Object.vlist[verts[j]]->vec.entry[1];

		/*Use non-normalized vector field 07/19/07*/
		//vx[j] = Object.vlist[verts[j]]->vec_J.entry[0];
		//vy[j] = Object.vlist[verts[j]]->vec_J.entry[1];
	}


	//////////////////////////////////////////////////////////////////////////////
	/////Second way to calculate the a b c d e f values

	double coord[3][3],  *inver_coord ;  //inver_coord[3][3];

	for(i = 0; i < 3; i++)
	{
		coord[0][i] = x[i];
		coord[1][i] = y[i];
		coord[2][i] = 1.;
	}

	inver_coord = MatrixOpp((double*)coord, 3, 3);

	icMatrix3x3 result, rightM;
	result.set(vx[0],vx[1],vx[2],  vy[0],vy[1],vy[2],  1,1,1);
	rightM.set(inver_coord[0], inver_coord[1], inver_coord[2],
				inver_coord[3], inver_coord[4], inver_coord[5],
				inver_coord[6], inver_coord[7], inver_coord[8]);


	result.rightMultiply(rightM);

	FieldMatrix.set(result);

	a = result.entry[0][0];
	b = result.entry[0][1];
	c = result.entry[0][2];
	d = result.entry[1][0];
	e = result.entry[1][1];
	f = result.entry[1][2];

	//need to store it as a part of the information of the elements or unknow singularities

	//use to calculate the coordinates of the singularity 3/25/06
	x_cp = (f*b - c*e)/(a*e - d*b);
	y_cp = (c*d - a*f)/(a*e - b*d);

	//Getting the eigenvalue here

	double A, B, C, delta;

	A = 1;
	B = -(a + e);
	C = (a * e - b * d);

	delta =  B*B - 4*A*C;

	if( delta >= 0)
	{
		i1 = i2 = 0;
		r1 = (-B + sqrt(delta))/2;
		r2 = (-B - sqrt(delta))/2;
	}
	else
	{
		r1 = r2 = -B/2.;
		i1 = sqrt(-delta)/2;
		i2 = -sqrt(-delta)/2;
	}

    ////
	//if(fabs(r1) < 1e-4) r1 = 0;
	//if(fabs(r2) < 1e-4) r2 = 0;

	////////////////////////////////////////////////////////////////
	if( r1 > 0 && r2 > 0 && i1 == 0 && i2 == 0)
		return 1; //it is a source

	else if( r1 < 0 && r2 < 0 && i1 == 0 && i2 == 0)
		return 2; //it is a sink

	else if( r1 * r2 < 0 && i1 == 0 && i2 == 0)
		return 3; //it is a saddle

	else if( r1 == 0 && r2 == 0 && i1 != 0 && i2 != 0)
	//else if( fabs(r1) <= 1e-4 && fabs(r2) <= 1e-4 && i1 != 0 && i2 != 0)
		return 4; //it is a center

	else if( r1 < 0 && r2 < 0 && i1 != 0 && i2 != 0)
		return 6; //it is an attracting focus

	else if( r1 > 0 && r2 > 0 && i1 != 0 && i2 != 0)
		return 7; //it is a repelling focus

	else
		return 0; //Unknow, there must be some error here!

}


/*----------------------------------------------------------------------*/
////Routines for matrix operations

/****************************************************
Maybe better for 2X2 and 3X3 matrix inversion
****************************************************/

double *MatrixOpp(double A[],int m,int n) //inverse 
{ 
    int i,j,x,y,k; 
    double *SP=NULL,*AB=NULL,*B=NULL,XX,*C; 
    SP=(double *)malloc(m*n*sizeof(double)); 
    AB=(double *)malloc(m*n*sizeof(double)); 
    B=(double *)malloc(m*n*sizeof(double)); 
    
    XX=Surplus(A,m,n); 
    XX=1/XX; 
    
    for(i=0;i<m;i++) 
	for(j=0;j<n;j++) 
	{
		for(k=0;k<m*n;k++) 
			B[k]=A[k]; 
			{
				for(x=0;x<n;x++) 
				B[i*n+x]=0; 
				for(y=0;y<m;y++) 
				B[m*y+j]=0; 
				B[i*n+j]=1; 
				SP[i*n+j]=Surplus(B,m,n); 
				AB[i*n+j]=XX*SP[i*n+j]; 
			} 
	} 

    C=MatrixInver(AB,m,n); 

	free(SP);
	free(AB);
	free(B);
    
    return C; 
} 
    
double * MatrixInver(double A[],int m,int n) //zhuanzhi
{ 
    int i,j; 
    double *B=NULL; 
    B=(double *)malloc(m*n*sizeof(double)); 
    
    for(i=0;i<n;i++) 
	for(j=0;j<m;j++) 
		B[i*m+j]=A[j*n+i]; 
    
    return B; 
} 
    
double Surplus(double A[],int m,int n) //hanglieshi
{ 
    
    int i,j,k,p,r; 
    double XX,temp=1,temp1=1,s=0,s1=0; 
    
    if(n==2) 
    {for(i=0;i<m;i++) 
    for(j=0;j<n;j++) 
    if((i+j)%2) temp1*=A[i*n+j]; 
    else temp*=A[i*n+j]; 
    XX=temp-temp1;} 
    else{ 
    for(k=0;k<n;k++) 
    {for(i=0,j=k;i<m,j<n;i++,j++) 
    temp*=A[i*n+j]; 
    if(m-i) 
    {for(p=m-i,r=m-1;p>0;p--,r--) 
    temp*=A[r*n+p-1];} 
    s+=temp; 
    temp=1; 
    } 
    
    for(k=n-1;k>=0;k--) 
    {for(i=0,j=k;i<m,j>=0;i++,j--) 
    temp1*=A[i*n+j]; 
    if(m-i) 
    {for(p=m-1,r=i;r<m;p--,r++) 
    temp1*=A[r*n+p];} 
    s1+=temp1; 
    temp1=1; 
    } 
    
    XX=s-s1;} 
    return XX; 
} 



#include "CalDeformation.h"


int get_singularity_type(int MarkTriangleID)
{
	double jacobian[2][2];
	//double c, f;
	get_vf_formula_tri(MarkTriangleID, jacobian, c, f);


	a = jacobian[0][0];
	b = jacobian[0][1];
	d = jacobian[1][0];
	e = jacobian[1][1];

	//need to store it as a part of the information of the elements or unknow singularities

	//use to calculate the coordinates of the singularity 3/25/06
	x_cp = (f*b - c*e)/(a*e - d*b);
	y_cp = (c*d - a*f)/(a*e - b*d);

	//Getting the eigenvalue here

	double A, B, C, delta;

	A = 1;
	B = -(a + e);
	C = (a * e - b * d);

	delta =  B*B - 4*A*C;

	if( delta >= 0)
	{
		i1 = i2 = 0;
		r1 = (-B + sqrt(delta))/2;
		r2 = (-B - sqrt(delta))/2;
	}
	else
	{
		r1 = r2 = -B/2.;
		i1 = sqrt(-delta)/2;
		i2 = -sqrt(-delta)/2;
	}

	////////////////////////////////////////////////////////////////
	if( r1 > 0 && r2 > 0 && i1 == 0 && i2 == 0)
		return 1; //it is a source

	else if( r1 < 0 && r2 < 0 && i1 == 0 && i2 == 0)
		return 2; //it is a sink

	else if( r1 * r2 < 0 && i1 == 0 && i2 == 0)
		return 3; //it is a saddle

	else if( r1 == 0 && r2 == 0 && i1 != 0 && i2 != 0)
	//else if( fabs(r1) <= 1e-4 && fabs(r2) <= 1e-4 && i1 != 0 && i2 != 0)
		return 4; //it is a center

	else if( r1 < 0 && r2 < 0 && i1 != 0 && i2 != 0)
		return 6; //it is an attracting focus

	else if( r1 > 0 && r2 > 0 && i1 != 0 && i2 != 0)
		return 7; //it is a repelling focus

	else
		return 0; //Unknow, there must be some error here!
}


void get_loc_Jacobian(int tri, double Jacobian[2][2], double &c, double &f)
{
	double x[3], y[3], vx[3], vy[3];
	int i, j;
	Face *face;
	int *verts;

	///These information has been obtained before
	///Can we reuse them and enhance the performance 
    face = Object.flist[tri];
	verts = face->verts;

	//For each triangle, calculate the vector for all vertices
	for (j=0; j<face->nverts; j++) {
            x[j] = face->xy[j][0];
            y[j] = face->xy[j][1];

			/* Use non-normalized vector field*/
			vx[j] = face->direct_vec[j].entry[0];
			vy[j] = face->direct_vec[j].entry[1];
	}

	/*define local variable, all the calculations are done locally
	and need to transfer back to global frame*/

	double la, lb, lc, ld, le, lf;

	/*for v0 (0, 0)*/
	lc = vx[0];
	lf = vy[0];

	/*for v1(x1, 0)*/
	la = (vx[1]-lc)/x[1];
	ld = (vy[1]-lf)/x[1];

	/*from v2(x2, y2)*/
	lb = (vx[2]-la*x[2]-lc)/y[2];
	le = (vy[2]-ld*x[2]-lf)/y[2];

	Jacobian[0][0] = la;
	Jacobian[0][1] = lb;
	Jacobian[1][0] = ld;
	Jacobian[1][1] = le;

	c = lc;
	f = lf;
}

extern 	void Get2DBarycentricFacters(int , double , double,  double[]);

/*we use the Jacobian obtained through local coordinate calculation
to judge the type of the fixed point*/
int get_singularity_type_loc(int tri, double Jacobian[2][2], double x_loc, double y_loc, double alpha[3])
{
	a = Jacobian[0][0];
	b = Jacobian[0][1];
	d = Jacobian[1][0];
	e = Jacobian[1][1];

	//need to store it as a part of the information of the elements or unknow singularities

	//use to calculate the coordinates of the fixed point
	/*note that the coordinates here is local coordinates!*/
	
	//double x_cp_loc = (f*b - c*e)/(a*e - d*b);
	//double y_cp_loc = (c*d - a*f)/(a*e - b*d);
	x_loc = (f*b - c*e)/(a*e - d*b);
	y_loc = (c*d - a*f)/(a*e - b*d);

	/*get the bary centri coordinates of this point*/
	//Get2DBarycentricFacters(tri, x_cp_loc, y_cp_loc, alpha);
	Get2DBarycentricFacters(tri, x_loc, y_loc, alpha);

	//Face *face = Object.flist[tri];
	//x_cp=alpha[0]*Object.vlist[face->verts[0]]->x
	//	+alpha[1]*Object.vlist[face->verts[1]]->x
	//	+alpha[2]*Object.vlist[face->verts[2]]->x;
	//y_cp=alpha[0]*Object.vlist[face->verts[0]]->y
	//	+alpha[1]*Object.vlist[face->verts[1]]->y
	//	+alpha[2]*Object.vlist[face->verts[2]]->y;

	//Getting the eigenvalue here
	
	/*NOTE that the eigen vector now is also local*/

	double A, B, C, delta;

	A = 1;
	B = -(a + e);
	C = (a * e - b * d);

	delta =  B*B - 4*A*C;

	if( delta >= 0)
	{
		i1 = i2 = 0;
		r1 = (-B + sqrt(delta))/2;
		r2 = (-B - sqrt(delta))/2;
	}
	else
	{
		r1 = r2 = -B/2.;
		i1 = sqrt(-delta)/2;
		i2 = -sqrt(-delta)/2;
	}

	////////////////////////////////////////////////////////////////
	if( r1 > 0 && r2 > 0 && i1 == 0 && i2 == 0)
		return 1; //it is a source

	else if( r1 < 0 && r2 < 0 && i1 == 0 && i2 == 0)
		return 2; //it is a sink

	else if( r1 * r2 < 0 && i1 == 0 && i2 == 0)
		return 3; //it is a saddle

	else if( r1 == 0 && r2 == 0 && i1 != 0 && i2 != 0)
	//else if( fabs(r1) <= 1e-4 && fabs(r2) <= 1e-4 && i1 != 0 && i2 != 0)
		return 4; //it is a center

	else if( r1 < 0 && r2 < 0 && i1 != 0 && i2 != 0)
		return 6; //it is an attracting focus

	else if( r1 > 0 && r2 > 0 && i1 != 0 && i2 != 0)
		return 7; //it is a repelling focus

	else
		return 0; //Unknow, there must be some error here!

}


/***************************************************************
Compute the accurate coordinates after capture the ID of 
the triangles that contain fixed points.
In this new routine, we ignore the case that 
the fixed points locate at vertices
***************************************************************/
void compute_fixedpt_info(void)
{
	int i, j;
	Face *face;
	int *verts;
	double x_loc, y_loc; 
	icVector2 glob;
	double Jacobian[2][2];
	double alpha[3];

	////Initialize the unknown singularities link

	for(i = 0; i < FindTriangleIndex; i++)
	{
		//For each being captured triangle, compute the coordinates of singularities inside it
		face = Object.flist[MarkTriangleID[i]];
		verts = face->verts;

		get_loc_Jacobian(MarkTriangleID[i], Jacobian, c, f);

		/////Calculate the coordinates of singularities here

		///Add to the being captured singularities list
		//singularities[cur_singularity_index].type = GetSingType(MarkTriangleID[i]); //Judge the type of this singularity
		//singularities[cur_singularity_index].type = get_singularity_type(MarkTriangleID[i]); //Judge the type of this singularity
		singularities[cur_singularity_index].type = get_singularity_type_loc(MarkTriangleID[i], Jacobian,
			x_loc, y_loc, alpha); //Judge the type of this singularity

		/*compute the global coordinates of the fixed point*/
		x_cp=alpha[0]*Object.vlist[face->verts[0]]->x
			+alpha[1]*Object.vlist[face->verts[1]]->x
			+alpha[2]*Object.vlist[face->verts[2]]->x;
		y_cp=alpha[0]*Object.vlist[face->verts[0]]->y
			+alpha[1]*Object.vlist[face->verts[1]]->y
			+alpha[2]*Object.vlist[face->verts[2]]->y;

		singularities[cur_singularity_index].gcx = x_cp; ////set the center coordinates
		singularities[cur_singularity_index].gcy = y_cp;

		//singularities[cur_singularity_index].Jacobian.set(FieldMatrix); //store the field matrix
		double M[][3]={{Jacobian[0][0], Jacobian[0][1], 0},
			{Jacobian[1][0], Jacobian[1][1], 0},
			{0, 0, 1}
		};
		singularities[cur_singularity_index].Jacobian.set(M);
			
		singularities[cur_singularity_index].Triangle_ID = MarkTriangleID[i];

		//singularities[cur_singularity_index].node_index = -1;           //10/16/05
		Object.flist[MarkTriangleID[i]]->singularity_index = cur_singularity_index;
		
		if(singularities[cur_singularity_index].type == SADDLE)
		{
			
			singularities[cur_singularity_index].outgoing.entry[0]\
				= -singularities[cur_singularity_index].Jacobian.entry[0][1]; //outgoing direction
			singularities[cur_singularity_index].outgoing.entry[1]\
				= singularities[cur_singularity_index].Jacobian.entry[0][0]- r1;

			/*--------------------------------------------------*/
			/*need to map to global vector07/21/07*/
			glob = singularities[cur_singularity_index].outgoing.entry[0]*face->LX+
				singularities[cur_singularity_index].outgoing.entry[1]*face->LY;
			singularities[cur_singularity_index].outgoing = glob;
			/*--------------------------------------------------*/

			singularities[cur_singularity_index].incoming.entry[0]\
				= -singularities[cur_singularity_index].Jacobian.entry[0][1]; //incoming direction
			singularities[cur_singularity_index].incoming.entry[1]\
				= singularities[cur_singularity_index].Jacobian.entry[0][0]- r2;

			/*--------------------------------------------------*/
			/*need to map to global vector07/21/07*/
			glob = singularities[cur_singularity_index].incoming.entry[0]*face->LX+
				singularities[cur_singularity_index].incoming.entry[1]*face->LY;
			singularities[cur_singularity_index].incoming = glob;
			/*--------------------------------------------------*/

			normalize(singularities[cur_singularity_index].outgoing);
			normalize(singularities[cur_singularity_index].incoming);
		}
		singularities[cur_singularity_index].connected_limitcycles = NULL;
		singularities[cur_singularity_index].num_connected_limitcycles = 0;
		singularities[cur_singularity_index].connected = 0;
		cur_singularity_index++;

	}
}
