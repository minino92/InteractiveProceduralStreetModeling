/*
This file contains routines computing the deformation tensor of a given triangle

created and modified by Guoning Chen
Copyright c 2007
*/

#include "stdafx.h"

#include "VFDataStructure.h"

void get_deformation_tensor(icVector2 x[3], icVector2 x_img[3], icMatrix2x2 mat);
void get_deformation_tensor(double vx[3], double vy[3], 
							 double vx_img[3], double vy_img[3], 
							 double mat[2][2]);
void get_vf_formula_tri(int tri);
int get_type_deformation(double mat[2][2]);
int get_type_deformation(icMatrix2x2 mat);


extern Polygon3D Object;

/*
calculate the inverse matrix of a 2x2 matrix
*/
double cal_determinant2x2(double a[][2])
{
	return (a[0][0]*a[1][1]-a[0][1]*a[1][0]);
}


double cal_determinant(icMatrix2x2 mat)
{
	return (mat.entry[0][0]*mat.entry[1][1]-mat.entry[0][1]*mat.entry[1][0]);
}


/*
compute the inverse matrix of a 2X2 matrix
*/
bool cal_inverse2x2(double a[][2], double inverse_a[][2])
{
	double det;
	if((det=cal_determinant2x2(a)) == 0)  //degenerate case
		return false;
	inverse_a[0][0] = a[1][1]/det;
	inverse_a[1][1] = a[0][0]/det;
	inverse_a[0][1] = -a[0][1]/det;
	inverse_a[1][0] = -a[1][0]/det;
	return true;
}

/*
compute the inverse matrix of a 2X2 matrix
*/
bool cal_inverse2x2(icMatrix2x2 a, icMatrix2x2 inverse_a)
{
	double det;
	if((det=cal_determinant(a)) == 0)  //degenerate case
		return false;
	inverse_a.entry[0][0] = a.entry[1][1]/det;
	inverse_a.entry[1][1] = a.entry[0][0]/det;
	inverse_a.entry[0][1] = -a.entry[0][1]/det;
	inverse_a.entry[1][0] = -a.entry[1][0]/det;
	return true;
}



/*
calculate the eigen values
*/
int cal_eigenval_matrix2x2(double mat[2][2], double evalues[2])
{
	double la, lb, lc, ld, A, B, C, delta;

	la = mat[0][0];
	lb = mat[0][1];
	lc = mat[1][0];
	ld = mat[1][1];

	evalues[0] = evalues[1] = 0.;

	A = 1;
	B = -(la + ld);
	C = (la * ld - lb * lc);                

	delta =  B*B - 4*A*C;

	if(delta >= 0)
	{
		evalues[0] = (-B + sqrt(delta))/2;
		evalues[1] = (-B - sqrt(delta))/2;

		if(delta == 0)
			return 1;

		return 2;
	}

	else{  /*contain imaginary eigen values*/

		/*just get the real part of the eigen values*/
		evalues[0] = evalues[1] = -B/2.;

		if(B==0)
			return -1;

		return -2;  //two complex eigen values
	}
}

/*
calculate the eigen values
*/
int cal_eigenval_matrix2x2(icMatrix2x2 mat, double evalues[2])
{
	double la, lb, lc, ld, A, B, C, delta;

	la = mat.entry[0][0];
	lb = mat.entry[0][1];
	lc = mat.entry[1][0];
	ld = mat.entry[1][1];

	evalues[0] = evalues[1] = 0.;

	A = 1;
	B = -(la + ld);
	C = (la * ld - lb * lc);                

	delta =  B*B - 4*A*C;

	if(delta >= 0)
	{
		evalues[0] = (-B + sqrt(delta))/2;
		evalues[1] = (-B - sqrt(delta))/2;

		if(delta == 0)
			return 1;

		return 2;
	}

	else{  /*contain imaginary eigen values*/

		/*just get the real part of the eigen values*/
		evalues[0] = evalues[1] = -B/2.;

		if(B==0)
			return -1;

		return -2;  //two complex eigen values
	}
}


void cal_eigenvectors(icMatrix2x2 mat, double evalues[2], icVector2 ev[2])
{
	//first calculate the dominant of the matrix, if it is zero, return
	if(abs(cal_determinant(mat)) < 1e-30)
		return;

	ev[0].entry[0] = mat.entry[1][1] - evalues[0];
	ev[0].entry[1] = - mat.entry[1][0];
	
	ev[1].entry[0] = mat.entry[1][1] - evalues[1];
	ev[1].entry[1] = - mat.entry[1][0];

}


/* The routine of calculate eigen vectors for a given 2x2 matrix */

void cal_eigenvec_for_Matrix2x2(icMatrix2x2 mat, icVector2 evec[2])
{
	double la, lb, lc, ld, A, B, C, delta;

	la = mat.entry[0][0];
	lb = mat.entry[0][1];
	lc = mat.entry[1][0];
	ld = mat.entry[1][1];

	double evalues[2] = {0.};

	A = 1;
	B = -(la + ld);
	C = (la * ld - lb * lc);                

	delta =  B*B - 4*A*C;

	if(delta >= 0)
	{
		evalues[0] = (-B + sqrt(delta))/2;
		evalues[1] = (-B - sqrt(delta))/2;

		////find the largest (absolute) eigen values, and store it to the first element
		//if(abs(evalues[0]) < abs(evalues[1]))
		//{
			double temp = evalues[0];
			evalues[0] = evalues[1];
			evalues[1] = temp;
		//}

		//for real eigen values, we calculate the eigen vectors of it
		cal_eigenvectors(mat, evalues, evec);

		normalize(evec[0]);
		normalize(evec[1]);
	}

	else
	{
		evec[0].set(0, 0);
		evec[1].set(0, 0);
	}
}


/*
Calculate the deformation tensor according to the original triangle 
and the approximate image of the triangle
*/
void get_deformation_tensor(icVector2 x[3], icVector2 x_img[3], icMatrix2x2 mat)
{
	icVector2 n_x[3], n_x_img[3];
	int i;

	/*translate the (x1, y1) to be the origin*/
	n_x[0].entry[0] = 0;
	n_x[1].entry[1] = 0;
	for(i=1; i<3; i++)
		n_x[i] = x[i]-x[0];
	for(i=0; i<3; i++)
		n_x_img[i] = x_img[i]-x[0];

	/*solve for the deformation tensor*/
	double a, b, c, d, e, f;
	c = n_x_img[0].entry[0];
	f = n_x_img[0].entry[1];

	icMatrix2x2 M, inverse_M;
	double M_entries[][2]={{n_x[1].entry[0], n_x[1].entry[1]},{n_x[2].entry[0], n_x[2].entry[1]}};
	M.set(M_entries);

	/*calculate the inverse matrix of M*/
	cal_inverse2x2(M, inverse_M);

	icVector2 B;

	/*solve for a and b*/
	B.set(n_x_img[1].entry[0]-c, n_x_img[2].entry[0]-c);
	a = inverse_M.entry[0][0]*B.entry[0]+inverse_M.entry[0][1]*B.entry[1];
	b = inverse_M.entry[1][0]*B.entry[0]+inverse_M.entry[1][1]*B.entry[1];


	/*solve for d and e*/
	B.set(n_x_img[1].entry[1]-f, n_x_img[2].entry[1]-f);
	d = inverse_M.entry[0][0]*B.entry[0]+inverse_M.entry[0][1]*B.entry[1];
	e = inverse_M.entry[1][0]*B.entry[0]+inverse_M.entry[1][1]*B.entry[1];

	double deform_entries[2][2]={{a, b},{d, e}};
	mat.set(deform_entries);
}


/*
Calculate the deformation tensor according to the original triangle 
and the approximate image of the triangle
*/
void get_deformation_tensor(double vx[3], double vy[3], 
							 double vx_img[3], double vy_img[3], 
							 double mat[2][2])
{
	double n_vx[3], n_vy[3], n_vx_img[3], n_vy_img[3];
	int i;

	/*translate the (x1, y1) to be the origin*/
	n_vx[0] = n_vy[0] = 0;
	for(i=1; i<3; i++)
	{
		n_vx[i] = vx[i]-vx[0];
		n_vy[i] = vy[i]-vy[0];
	}
	for(i=0; i<3; i++)
	{
		n_vx_img[i] = vx_img[i]-vx[0];
		n_vy_img[i] = vy_img[i]-vy[0];
	}

	/*solve for the deformation tensor*/
	double a, b, c, d, e, f;
	c = n_vx_img[0];
	f = n_vy_img[0];

	double M_entries[][2]={{n_vx[1], n_vy[1]},{n_vx[2], n_vy[2]}};
	double inverse_M[2][2]={0.};

	/*calculate the inverse matrix of M*/
	cal_inverse2x2(M_entries, inverse_M);

	double B[2]={n_vx_img[1]-c, n_vx_img[2]-c};

	/*solve for a and b*/
	a = inverse_M[0][0]*B[0]+inverse_M[0][1]*B[1];
	b = inverse_M[1][0]*B[0]+inverse_M[1][1]*B[1];


	/*solve for d and e*/
	B[0] = n_vy_img[1]-f; B[1] = n_vy_img[2]-f;
	d = inverse_M[0][0]*B[0]+inverse_M[0][1]*B[1];
	e = inverse_M[1][0]*B[0]+inverse_M[1][1]*B[1];

	mat[0][0] = a;
	mat[0][1] = b;
	mat[1][0] = d;
	mat[1][1] = e;
}

void get_vf_formula_tri(int tri, double Jacobian[2][2], double &c, double &f)
{
	double vx[3], vy[3], vx_img[3], vy_img[3];
	double n_vx[3], n_vy[3], n_vx_img[3], n_vy_img[3];
	int i, j;
	Face *face;
	int *verts;

	///These information has been obtained before
	///Can we reuse them and enhance the performance 
    face = Object.flist[tri];
	verts = face->verts;

	//For each triangle, calculate the vector for all vertices
	for (j=0; j<face->nverts; j++) {

            vx[j] = Object.vlist[verts[j]]->x;
			vy[j] = Object.vlist[verts[j]]->y;

			/* Use non-normalized vector field*/
			vx_img[j] = Object.vlist[verts[j]]->vec_J.entry[0];
			vy_img[j] = Object.vlist[verts[j]]->vec_J.entry[1];
	}

	/*translate the (x1, y1) to be the origin*/
	n_vx[0] = n_vy[0] = 0;
	for(i=1; i<3; i++)
	{
		n_vx[i] = vx[i]-vx[0];
		n_vy[i] = vy[i]-vy[0];
	}
	for(i=0; i<3; i++)
	{
		n_vx_img[i] = vx_img[i]-vx[0];
		n_vy_img[i] = vy_img[i]-vy[0];
	}

	/*solve for the deformation tensor*/
	double a, b, d, e;
	c = n_vx_img[0];
	f = n_vy_img[0];

	double M_entries[][2]={{n_vx[1], n_vy[1]},{n_vx[2], n_vy[2]}};
	double inverse_M[2][2]={0.};

	/*calculate the inverse matrix of M*/
	cal_inverse2x2(M_entries, inverse_M);

	double B[2]={n_vx_img[1]-c, n_vx_img[2]-c};

	/*solve for a and b*/
	a = inverse_M[0][0]*B[0]+inverse_M[0][1]*B[1];
	b = inverse_M[1][0]*B[0]+inverse_M[1][1]*B[1];


	/*solve for d and e*/
	B[0] = n_vy_img[1]-f; B[1] = n_vy_img[2]-f;
	d = inverse_M[0][0]*B[0]+inverse_M[0][1]*B[1];
	e = inverse_M[1][0]*B[0]+inverse_M[1][1]*B[1];

	Jacobian[0][0] = a;
	Jacobian[0][1] = b;
	Jacobian[1][0] = d;
	Jacobian[1][1] = e;
}



int get_type_deformation(double mat[2][2])
{
	double evalues[2] = {0.};
	cal_eigenval_matrix2x2(mat, evalues);
	if(evalues[0] > 0 && evalues[1] > 0)
		return 0;   //two positive eigen values
	if(evalues[0] < 0 && evalues[0] < 0)
		return 1;
	return 2;
}


int get_type_deformation(icMatrix2x2 mat)
{
	double evalues[2] = {0.};
	cal_eigenval_matrix2x2(mat, evalues);
	if(evalues[0] > 0 && evalues[1] > 0)
		return 0;   //two positive eigen values
	if(evalues[0] < 0 && evalues[0] < 0)
		return 1;
	return 2;
}

/*
Calculate the deformation tensor and its type for specified triangle 
*/
void cal_deform_tri(int tri, double t)
{
	Face *face = Object.flist[tri];
	double vx[3], vy[3], vx_img[3], vy_img[3];
	int i;

	/*trace the three vertices with fixed time t to obtain the images of them*/
	for(i=0; i<3; i++)
	{
		Vertex *v = Object.vlist[face->verts[i]];
		vx_img[i] = vx[i] = v->x; 
		vy_img[i] = vy[i] = v->y;

		/*trace and update vx_img[i], vy_img[i]*/
	}

	/*we use global coordinates to compute the deformation tensor of the triangle*/
	double deform_tensor[2][2] = {0.};
	get_deformation_tensor(vx, vy, vx_img, vy_img, deform_tensor);

	int type = get_type_deformation(deform_tensor);

	/*color the triangle according to the type
	0:  green (divergent); 1: red (convergent); 2:  blue (stretching)
	*/
}