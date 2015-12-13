/*
This file implements the tensor field analysis
Especially degenerate point/singularity detection
NOTE: the data structure has combined with current data structure
if the data structure is changed, according changes are necessary.
09/18/2007
*/
#include "stdafx.h"

#include "VFDataStructure.h"
#include "caldeformation.h"

#include "tensorvis.h"

#include "VFAnalysis.h"

#include "tensoranalysis.h"

int ndegenerate_tris = 0;
int *degenerate_tris = NULL;
int MaxNumDegeneratePts = 200;

extern Polygon3D Object;
extern QuadMesh *quadmesh;
extern double  ten_tmax ;
extern double  ten_dmax ;
extern double euler_stepsize;


double ep = 1.e-8;  //error
double RK45_hstep; //the step size for current integration step
double RK45_hnext; //the predict step size for next step

int g_face, g_type;
double predict_stepsize;/* = quadmesh->xinterval;*/



extern bool StoreToGlobalList(CurvePoints *temp, int num);
extern int get_cellID_picking(double x, double y);
extern int get_cellID_givencoords(double x, double y);

/*we will have a tensor degenerate point class later
class DegeneratePtList
{
public:
      DegeneratePt *degpts;
	  int ndegpts;
	  int curMaxNum;
}
*/

DegeneratePt *degpts = NULL;
int ndegpts = 0;
//int curMaxNumDegPts;

/*
Initialize the degenerate point list
*/
void init_degpts()
{
	//curMaxNumDegPts = 200;
	if(degpts == NULL)
	{
		degpts = (DegeneratePt*)malloc(sizeof(DegeneratePt)*MaxNumDegeneratePts);
		if(degpts == NULL)
			exit(-1);
	}
	ndegpts = 0;

	int i;
	for(i=0; i<quadmesh->nfaces; i++)
	{
		quadmesh->quadcells[i]->degpt_index = -1;
	}
}


/*
extract the triangles that contain degenerate points/singularities using tensor index
*/
void detect_degeneratepts()
{
	unsigned int i, j;
	Face *face;
	int *verts;
	icVector2 vec[3];      //vectors at three vertices
    double  theta[3];      //for storing the angles between two vector for Gauss circle judgement

	double  vec_ang[3];  //the angle for each vector under the polar frame of current triangle
	double  ang_sum;

	////Initialize
	ndegenerate_tris = 0;

	if(degenerate_tris == NULL)
		degenerate_tris = (int*) malloc(sizeof(int) * MaxNumDegeneratePts); //default range is 200

	////These flags are used for pair cancellation
	////Initialize all the flags!!!

	//for(i = 0; i < Object.nfaces; i++)
	//{
	//	Object.flist[i]->contain_singularity = 0;
	//	Object.flist[i]->singularity_index = -1;
	//}

	////initialize the singularities list

	////Calculate the Poincare index
	for (i = 0; i < Object.nfaces; i++) {
		face = Object.flist[i];
		verts = face->verts;

		ang_sum = 0;

		for(j=0; j<3; j++)
			vec_ang[j] = Object.vlist[verts[j]]->tensor_major_ang;

		for(j = 0; j < face->nverts; j++)
		{
			double temp1, temp2;
			temp1 = vec_ang[(j+1)%3] - vec_ang[j];
			//temp2 = vec_ang[(j+1)%3] - vec_ang[j] - M_PI;

			//if(fabs(temp1)<fabs(temp2))
			//	theta[j] = temp1;
			//else
			//	theta[j] = temp2;

			if(fabs(temp1) > M_PI/2.)
			{
				if(temp1 > 0)
					theta[j] = -(M_PI-temp1);
				if(temp1 < 0)
					theta[j] = -(-M_PI-temp1);
			}
			else
				theta[j] = temp1;

			ang_sum += theta[j];
		}

		double index = ang_sum/(2*M_PI);


		if(fabs(index) >= 0.5)
		{
			//The triangle must have singularities inside, mark it as yellow color
			//Still need to judge whether it is one of current singularities or not
			degenerate_tris[ndegenerate_tris] = i;
			ndegenerate_tris ++;

			if(fabs(index-0.5)<1e-8)  /*it is a wedge*/
				Object.flist[i]->degenerate_type = 0;
			else if(fabs(index+0.5)<1e-8)  /*it is a trisector*/
				Object.flist[i]->degenerate_type = 1;
			else if(fabs(index-1)<1e-8)    /*it is a node/center*/
				Object.flist[i]->degenerate_type = 2;
			else if(fabs(index+1)<1e-8)    /*it is a saddle*/
				Object.flist[i]->degenerate_type = 3;

			if(ndegenerate_tris >= MaxNumDegeneratePts - 1)
			{
				MaxNumDegeneratePts += 50;
				degenerate_tris = (int*) realloc(degenerate_tris, sizeof(int) * MaxNumDegeneratePts);
				degpts = (DegeneratePt*)realloc(degpts, sizeof(DegeneratePt)*MaxNumDegeneratePts);
			}
		}
	}

	if(ndegenerate_tris > 0)
		locate_degeneratepts_tranvec();
}

/*
The following routine uses the obtained vector field 
from the original (multiply major eigenvector by 2)
tensor field to extract the degenerate points.
*/
void detect_degeneratepts_tranvec()
{
	unsigned int i, j;
	Face *face;
	Vertex *v;
	icVector2 vec[3];      //vectors at three vertices
    double  theta[3];      //for storing the angles between two vector for Gauss circle judgement

	double  vec_ang[3];  //the angle for each vector under the polar frame of current triangle
	double  ang_sum;

	////Initialize
	ndegenerate_tris = 0;

	if(degenerate_tris == NULL)
		degenerate_tris = (int*) malloc(sizeof(int) * MaxNumDegeneratePts); //default range is 200



	////Calculate the Poincare index
	for (i = 0; i < Object.nfaces; i++) {
		face = Object.flist[i];

		ang_sum = 0;

		for(j=0; j<3; j++)
		{
			v = Object.vlist[face->verts[j]];
			vec_ang[j] = atan2(v->tran_vec.entry[1], v->tran_vec.entry[0]);
			if(vec_ang[j] < 0) vec_ang[j] += 2 * M_PI;
		}

		for(j = 0; j < face->nverts; j++)
		{
			theta[j] = vec_ang[(j+1)%3] - vec_ang[j];

			if( theta[j] < -M_PI)
				theta[j] += 2 * M_PI;
			
			if( theta[j] > M_PI)
				theta[j] -= 2 * M_PI;

			ang_sum += theta[j];
		}

		double index = ang_sum/(2*M_PI);


		if(fabs(index) >= 1.)
		{
			//The triangle must have singularities inside, mark it as yellow color
			//Still need to judge whether it is one of current singularities or not
			degenerate_tris[ndegenerate_tris] = i;
			ndegenerate_tris ++;

			if(fabs(index-1)<1e-8)  /*it is a wedge*/
				Object.flist[i]->degenerate_type = 0;
			else if(fabs(index+1)<1e-8)  /*it is a trisector*/
				Object.flist[i]->degenerate_type = 1;
			else if(fabs(index-2)<1e-8)    /*it is a node/center*/
				Object.flist[i]->degenerate_type = 2;
			else if(fabs(index+2)<1e-8)    /*it is a saddle*/
				Object.flist[i]->degenerate_type = 3;

			if(ndegenerate_tris >= MaxNumDegeneratePts - 1)
			{
				MaxNumDegeneratePts += 50;
				degenerate_tris = (int*) realloc(degenerate_tris, sizeof(int) * MaxNumDegeneratePts);
				degpts = (DegeneratePt*)realloc(degpts, sizeof(DegeneratePt)*MaxNumDegeneratePts);
			}
		}
	}

	if(ndegenerate_tris>0) /*we find some degenerate triangles*/
		locate_degeneratepts_tranvec();
}

void locate_degeneratepts()
{
}

/***************************************************************
Compute the accurate coordinates after capture the ID of 
the triangles that contain singularities using new method 1/16
***************************************************************/
void locate_degeneratepts_tranvec(void)
{
	int i, j, k;
	Face *face;
	int *verts;
	double a, b, c, d, e, f; //
	double x[3], y[3]; //
	double vx[3], vy[3];//
	double vdx, vdy; //
	double x_cp, y_cp;

	////Initialize the unknown singularities link
	ndegpts = 0;

	for(i = 0; i < ndegenerate_tris; i++)
	{
		//For each being captured triangle, compute the coordinates of singularities inside it
		face = Object.flist[degenerate_tris[i]];
		verts = face->verts;
			
		//For each triangle, calculate the vector for all vertices
		for (j=0; j<face->nverts; j++) {
            x[j] = Object.vlist[verts[j]]->x;
			y[j] = Object.vlist[verts[j]]->y;

			////using the vectors stored in the vertices
			
			/* Use normalized vector field*/
			vx[j] = Object.vlist[verts[j]]->tran_vec.entry[0];
			vy[j] = Object.vlist[verts[j]]->tran_vec.entry[1];

		}

		/////Calculate the coordinates of singularities here
		double coord[3][3],  *inver_coord ;  //inver_coord[3][3];

		for(k = 0; k < 3; k++)
		{
			coord[0][k] = x[k];
			coord[1][k] = y[k];
			coord[2][k] = 1.;
		}

		inver_coord = MatrixOpp((double*)coord, 3, 3);

		icMatrix3x3 result, rightM;
		result.set(vx[0],vx[1],vx[2],  vy[0],vy[1],vy[2],  1,1,1);
		rightM.set(inver_coord[0], inver_coord[1], inver_coord[2],
					inver_coord[3], inver_coord[4], inver_coord[5],
					inver_coord[6], inver_coord[7], inver_coord[8]);


		result.rightMultiply(rightM);

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

		degpts[ndegpts].gcx = x_cp;
		degpts[ndegpts].gcy = y_cp;
		degpts[ndegpts].degpt_index = ndegpts;
		degpts[ndegpts].type = face->degenerate_type;
		degpts[ndegpts].Triangle_ID = degenerate_tris[i];

		/*compute the separatrices*/
		degpts[ndegpts].nseps = 0;
		compute_separatrixVec_degpt(ndegpts);

		ndegpts++;
	}
}


void compute_tensor_at(int tri, double alpha[3], icMatrix2x2 &ten)
{
	Face *face=Object.flist[tri];
	Vertex *v;
	int i;
	ten.set(0.);
	for(i=0; i<3; i++)
	{
		v = Object.vlist[face->verts[i]];
		ten.entry[0][0] += alpha[i]*v->Jacobian.entry[0][0];
		ten.entry[0][1] += alpha[i]*v->Jacobian.entry[0][1];
		ten.entry[1][0] += alpha[i]*v->Jacobian.entry[1][0];
		ten.entry[1][1] += alpha[i]*v->Jacobian.entry[1][1];
	}
}

/*
given a degenerate point, compute its separatrices if exist
*/
void compute_separatrixVec_degpt(int degpt_index)
{
	Face *face = Object.flist[degpts[degpt_index].Triangle_ID];
	Vertex *v;
	icMatrix2x2 dv[3]; // the deviators of the tensors at the three vertices

	int i;

	/*obtain the deviators*/
	for(i=0; i<3; i++)
	{
		dv[i].set(0.);
		v=Object.vlist[face->verts[i]];
		double half_trace = 0.5*(v->Jacobian.entry[0][0]+v->Jacobian.entry[1][1]);
		dv[i].entry[0][0] = v->Jacobian.entry[0][0]-half_trace;
		dv[i].entry[1][1] = v->Jacobian.entry[1][1]-half_trace;
		dv[i].entry[0][1] = v->Jacobian.entry[0][1];
		dv[i].entry[1][0] = v->Jacobian.entry[1][0];
	}

	/*translate the tensor system according to the center of the degenerate point*/
	//for(i=0; i<3; i++)
	//{
	//}

	double a, b, c, d, e, f; //
	double x[3], y[3]; //
	double vx[3], vy[3];//

	for (i=0; i<face->nverts; i++) {
		v=Object.vlist[face->verts[i]];
        x[i] = v->x;
		y[i] = v->y;

		////using the vectors stored in the vertices
		
		/* Use normalized vector field*/
		vx[i] = dv[i].entry[0][0];
		vy[i] = dv[i].entry[0][1];

	}

	/////Calculate the jacobian of this linear system
	double coord[3][3],  *inver_coord ;  //inver_coord[3][3];
	for(int k = 0; k < 3; k++)
	{
		coord[0][k] = x[k];
		coord[1][k] = y[k];
		coord[2][k] = 1.;
	}

	inver_coord = MatrixOpp((double*)coord, 3, 3);

	icMatrix3x3 result, rightM;
	result.set(vx[0],vx[1],vx[2],  vy[0],vy[1],vy[2],  1,1,1);
	rightM.set(inver_coord[0], inver_coord[1], inver_coord[2],
				inver_coord[3], inver_coord[4], inver_coord[5],
				inver_coord[6], inver_coord[7], inver_coord[8]);


	result.rightMultiply(rightM);

	a = result.entry[0][0];
	b = result.entry[0][1];
	c = result.entry[0][2];
	d = result.entry[1][0];
	e = result.entry[1][1];
	f = result.entry[1][2];

	/*compute the solution for a cubic equations*/
	double solutions[4] = {0.};
	//int nroots = solve_ten_cubic(e, (d+2*b), (2*a-e), -d, solutions);
	int nroots = solve_ten_cubic_3(e, (d+2*b), (2*a-e), -d, solutions);

	if(nroots == 0)
		return;
	else if(nroots == 1 /*|| nroots == 4*/)
	{
		/*it is the separatrix of a wedge*/
		degpts[degpt_index].nseps = 1;
		degpts[degpt_index].s1_ang = atan(solutions[0]);
		degpts[degpt_index].s[0].entry[0] = cos(degpts[degpt_index].s1_ang);
		degpts[degpt_index].s[0].entry[1] = sin(degpts[degpt_index].s1_ang);
	}
	else if(nroots == 2|| nroots == 4)
	{
		/*they are the separatrices of a wedge*/
		degpts[degpt_index].nseps = 2;
		degpts[degpt_index].s1_ang = atan(solutions[0]);
		degpts[degpt_index].s[0].entry[0] = cos(degpts[degpt_index].s1_ang);
		degpts[degpt_index].s[0].entry[1] = sin(degpts[degpt_index].s1_ang);

		degpts[degpt_index].s2_ang = atan(solutions[1]);
		degpts[degpt_index].s[1].entry[0] = cos(degpts[degpt_index].s2_ang);
		degpts[degpt_index].s[1].entry[1] = sin(degpts[degpt_index].s2_ang);
	}
	else if(nroots == 3) /*they are the separatrices of a trisector*/
	{
		degpts[degpt_index].nseps = 3;
		degpts[degpt_index].s1_ang = atan(solutions[0]);
		degpts[degpt_index].s[0].entry[0] = cos(degpts[degpt_index].s1_ang);
		degpts[degpt_index].s[0].entry[1] = sin(degpts[degpt_index].s1_ang);

		degpts[degpt_index].s2_ang = atan(solutions[1]);
		degpts[degpt_index].s[1].entry[0] = cos(degpts[degpt_index].s2_ang);
		degpts[degpt_index].s[1].entry[1] = sin(degpts[degpt_index].s2_ang);

		degpts[degpt_index].s3_ang = atan(solutions[2]);
		degpts[degpt_index].s[2].entry[0] = cos(degpts[degpt_index].s3_ang);
		degpts[degpt_index].s[2].entry[1] = sin(degpts[degpt_index].s3_ang);
	}

}

int get_sign(double x)
{
	if(x>=0.0) return 1;
	else return -1;
}

int solve_ten_cubic(double a, double b, double c, double d, double solutions[4])
{
	if(fabs(a) < 1e-6)
		return solve_ten_quadratic(b, c, d, solutions);
	b /= a;
	c /= a;
	d /= a;
	double Q = (9.*b*c-27.*d-2.*b*b*b)/54.;
	double D = (3.*c-b*b)/9.;
	double R = D*D*D + Q*Q;
	double S, T;

	if(R>=0)
	{
		double s_R = sqrt(R);
		S = get_sign(Q+s_R)*pow(fabs(Q+s_R), 1./3);
		T = get_sign(Q-s_R)*pow(fabs(Q-s_R), 1./3);

		solutions[0] = -b/3.+(S+T);
		solutions[1] = -b/3.-(S+T)/2.;
		solutions[2] = -b/3.-(S+T)/2.;
		solutions[3] = fabs(sqrt(3.)/2.*(S-T));
		if(solutions[3] == 0) /*no complex roots*/
		{
			if(solutions[1] != solutions[0] )
				return 2;  /*two real roots*/
			else 
				return 1;
		}
		else /*contains two complex roots*/
		{
			return 4;
		}
	}

	else   /*we have distinct real roots*/
	{
		double th = acos(Q/sqrt(-D*D*D));
		solutions[0] = 2*sqrt(-D)*cos(th/3.)-b/3;
		solutions[1] = 2*sqrt(-D)*cos((th+2*M_PI)/3.)-b/3;
		solutions[2] = 2*sqrt(-D)*cos((th+4*M_PI)/3.)-b/3.;
		return 3;
	}
	return 0;
}


int solve_ten_cubic_3(double a, double b, double c, double d, double solutions[4])
{
	if(fabs(a) < 1e-6)
		return solve_ten_quadratic(b, c, d, solutions);
	double Q = (9*a*b*c-27*a*a*d-2*b*b*b)/(54*a*a*a);
	double D = (3*a*c-b*b)/(9*a*a);
	double R = D*D*D + Q*Q;
	double S, T;

	if(R>=0)
	{
		double s_R = sqrt(R);
		S = get_sign(Q+s_R)*pow(fabs(Q+s_R), 1./3);
		T = get_sign(Q-s_R)*pow(fabs(Q-s_R), 1./3);

		solutions[0] = -b/(3*a)+(S+T);
		solutions[1] = -b/(3*a)-(S+T)/2.;
		solutions[2] = -b/(3*a)-(S+T)/2.;
		solutions[3] = fabs(sqrt(3.)/2.*(S-T));
		if(solutions[3] == 0) /*no complex roots*/
		{
			if(solutions[1] != solutions[0] )
				return 2;  /*two real roots*/
			else 
				return 1;
		}
		else /*contains two complex roots*/
		{
			return 4;
		}
	}

	else   /*we have distinct real roots*/
	{
		double th = acos(Q/sqrt(-D*D*D));
		solutions[0] = 2*sqrt(-D)*cos(th/3.)-b/(3*a);
		solutions[1] = 2*sqrt(-D)*cos((th+2*M_PI)/3.)-b/(3*a);
		solutions[2] = 2*sqrt(-D)*cos((th+4*M_PI)/3.)-b/(3*a);
		return 3;
	}
	return 0;
}

/*   a*x3 + b*x2 + c*x + d = 0

Formula: 
  Step 1: Calculate p and q
          p = ( 3*c/a - (b/a)2 ) / 3
          q = ( 2*(b/a)3 - 9*b*c/a/a + 27*d/a ) / 27
  Step 2: Calculate discriminant D
          D = (p/3)3 + (q/2)2
  Step 3: Depending on the sign of D, you follow different strategy.
          If D<0, three distinct real roots.
          If D=0, three real roots of which at least two are equal.
          If D>0, one real and two complex roots.
  Step 3a: For D>0 and D=0
          Calculate u and v
          u = cubic_root(-q/2 + sqrt(D))
          v = cubic_root(-q/2 - sqrt(D))
          Find the three transformed roots
          y1 = u + v
          y2 = -(u+v)/2 + i (u-v)*sqrt(3)/2
          y3 = -(u+v)/2 - i (u-v)*sqrt(3)/2
  Step 3b: Alternately, for D<0, a trigonometric formulation is more convenient
          y1 =  2 * sqrt(|p|/3) * cos(phi/3)
          y2 = -2 * sqrt(|p|/3) * cos((phi+pi)/3)
          y3 = -2 * sqrt(|p|/3) * cos((phi-pi)/3)
          where phi = acos(-q/2/sqrt(|p|3/27))
                pi  = 3.141592654...
  Step 4  Finally, find the three roots
          x = y - b/a/3

*/
int solve_ten_cubic_2(double a, double b, double c, double d, double solutions[4])
{
    double p = (3*c/a - (b/a)*(b/a))/3.;
    double q = (2*pow((b/a),3) - 9*b*c/a/a + 27*d/a)/27.;
	double D = pow((p/3),3) + (q/2)*(q/2);

	if(D < 0 )/*3 distinct real roots*/
	{
		double phi = acos(-q/2/sqrt(pow(fabs(p),3)/27));
        solutions[0] =  2 * sqrt(fabs(p)/3) * cos(phi/3);
        solutions[1] = -2 * sqrt(fabs(p)/3) * cos((phi+M_PI)/3);
        solutions[2] = -2 * sqrt(fabs(p)/3) * cos((phi-M_PI)/3);
	}
	else if(D>=0)
	{
        double u = pow((-q/2 + sqrt(D)), 1./3.);
        double v = pow((-q/2 - sqrt(D)), 1./3.);
		solutions[0] = u + v;

		if(D == 0)
		{
			solutions[1] = -(u+v)/2;
			if(solutions[0] == solutions[1])
				return 1;
			else
				return 2;
		}

		return 1;
	}
}


int solve_ten_quadratic(double a, double b, double c, double solutions[2])
{
	if(a == 0)
	{
		solutions[0] = -c/b;
		return 1;
	}

	if(a == 0 && b == 0)
		return 0;

	double D = b*b-4*a*c;
	if(D>=0)
	{
		if(D==0)
		{
			solutions[0] = -b/(2*a); return 1;
		}
		else
		{
			D=sqrt(D);
			//solutions[0] = -b/(2*a)+D;
			//solutions[1] = -b/(2*a)-D;
			solutions[0] = (-b+D)/(2*a);
			solutions[1] = (-b-D)/(2*a);
			return 2;
		}
	}
	return 0;  /*there is no real root*/
}




/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/
/*******************************************************************************************/

/*******************************************************************************************/

/*Quad mesh
The following codes compute the tensor field using the Jacobian 
of the design vector field
*/
void cal_all_eigenvecs_quad()
{
	int i;

	for(i=0; i<quadmesh->nverts; i++)
	{
		if(quadmesh->quad_verts[i]->inland)
			cal_eigenvecs_onevert_quad(i);
	}

	/*normalize the major and minor field*/
	normalized_tensorfield_quad();
}

extern void cal_eigen_vector_sym(icMatrix2x2 ten, icVector2 ev[]);

void cal_eigenvecs_onevert_quad(int ver)
{
	QuadVertex *v=quadmesh->quad_verts[ver];
	double evalues[2];
	icVector2 ev[2];

	if(fabs(v->Jacobian.entry[0][0])<=1.e-7
		&&fabs(v->Jacobian.entry[0][1])<=1.e-7
		&&fabs(v->Jacobian.entry[1][0])<=1.e-7
		&&fabs(v->Jacobian.entry[1][1])<=1.e-7)
	{
		v->major.set(0,0);
		v->minor.set(0,0);
		v->major_ang=0;
		v->minor_ang=0;
		return;
	}

	/*calculate eigen values*/
	//cal_eigenval_matrix2x2(v->Jacobian, evalues);

	///*calculate eigen vectors*/
	//cal_eigenvectors(v->Jacobian, evalues, ev);

	/*calculate the eigen vectors using symmetric tensor eigen vector calculation 09/26/2007*/
	cal_eigen_vector_sym(v->Jacobian, ev);

	v->major = ev[0];
	v->minor = ev[1];

	/*compute the angle*/
	v->major_ang = atan2(v->major.entry[1], v->major.entry[0]);
	v->minor_ang = atan2(v->minor.entry[1], v->minor.entry[0]);

	/*save the angle of major field for degenerate points/singularity detection*/
	//if(v->major_ang<0)
	//	v->tensor_major_ang = M_PI+v->major_ang;
	//else
		v->tensor_major_ang = v->major_ang;
		
	/*obtain the multiplied vector, we will use this obtained
	vector field to extract singularities*/
	v->tran_vec.set(cos(2*v->tensor_major_ang), sin(2*v->tensor_major_ang));


	/*transfer to cos^2*/
	double major_ang_cos = cos(v->major_ang);
	double major_ang_sin = sin(v->major_ang);

	v->major_ang = major_ang_cos*major_ang_cos;
	if(major_ang_cos<0)
		v->major_cos = true;
	else
		v->major_cos = false;
	if(major_ang_sin<0)
		v->major_sin = true;
	else
		v->major_sin = false;

	double minor_ang_cos = cos(v->minor_ang);
	double minor_ang_sin = sin(v->minor_ang);

	v->minor_ang = minor_ang_cos*minor_ang_cos;
	//v->minor_ang = minor_ang_sin*minor_ang_sin;
	if(minor_ang_cos<0)
		v->minor_cos = true;
	else
		v->minor_cos = false;
	if(minor_ang_sin<0)
		v->minor_sin = true;
	else
		v->minor_sin = false;

}


/*
It actually normalizes the major and minor vector fields
*/
void normalized_tensorfield_quad()
{
	/*normalize the major and minor field*/
	int i;
    double r;
	QuadVertex *cur_v;

	for(i = 0; i < quadmesh->nverts; i++)
	{
		cur_v = quadmesh->quad_verts[i];

		/*normalize major field*/
		r = length(cur_v->major);
		r *= r;
		if (r < DistanceThreshold) 
		{
			r = DistanceThreshold;
			cur_v->major *= ten_dmax/r; 
		}
		r = length(cur_v->major);
		r *= r;
		if (r > ten_dmax*ten_dmax) { 
			r  = sqrt(r); 
			cur_v->major *= ten_dmax/r; 
		}
		
		/*normalize minor field*/
		r = length(cur_v->minor);
		r *= r;
		if (r < DistanceThreshold) 
		{
			r = DistanceThreshold;
			cur_v->minor *= ten_dmax/r; 
		}
		r = length(cur_v->minor);
		r *= r;
		if (r > ten_dmax*ten_dmax) { 
			r  = sqrt(r); 
			cur_v->minor *= ten_dmax/r; 
		}
	}
}




/*we use quad mesh to locate the degenerate points*/
void locate_degpts_cells_tranvec_quad(void)
{
	unsigned int i, j;
	QuadCell *face;
	QuadVertex *v;
	icVector2 vec[4];      //vectors at four vertices
    double  theta[4];      //for storing the angles between two vector for Gauss circle judgement

	double  vec_ang[4];  //the angle for each vector under the polar frame of current triangle
	double  ang_sum;

	////Initialize
	ndegenerate_tris = 0;

	if(degenerate_tris == NULL)
		degenerate_tris = (int*) malloc(sizeof(int) * MaxNumDegeneratePts); //default range is 200


	////Calculate the Poincare index
	for (i = 0; i < quadmesh->nfaces; i++) {
		face = quadmesh->quadcells[i];

		ang_sum = 0;

		for(j=0; j<face->nverts; j++)
		{
			v = quadmesh->quad_verts[face->verts[j]];
			vec_ang[j] = atan2(v->tran_vec.entry[1], v->tran_vec.entry[0]);
			if(vec_ang[j] < 0) vec_ang[j] += 2 * M_PI;
		}

		for(j = 0; j < face->nverts; j++)
		{
			theta[j] = vec_ang[(j+1)%face->nverts] - vec_ang[j];

			if( theta[j] < -M_PI)
				theta[j] += 2 * M_PI;
			
			if( theta[j] > M_PI)
				theta[j] -= 2 * M_PI;

			ang_sum += theta[j];
		}

		double index = ang_sum/(2*M_PI);


		/*here we need to allow some numerical errors 09/26/2007*/
		if(fabs(index) >= 1.- 1.e-1)
		{
			//The triangle must have singularities inside, mark it as yellow color
			//Still need to judge whether it is one of current singularities or not
			degenerate_tris[ndegenerate_tris] = i;
			ndegenerate_tris ++;

			if(fabs(index-1)<1e-2)  /*it is a wedge*/
				quadmesh->quadcells[i]->degenerate_type = 0;
			else if(fabs(index+1)<1e-2)  /*it is a trisector*/
				quadmesh->quadcells[i]->degenerate_type = 1;
			else if(fabs(index-2)<1e-2)    /*it is a node/center*/
				quadmesh->quadcells[i]->degenerate_type = 2;
			else if(fabs(index+2)<1e-2)    /*it is a saddle*/
				quadmesh->quadcells[i]->degenerate_type = 3;

			if(ndegenerate_tris >= MaxNumDegeneratePts - 1)
			{
				MaxNumDegeneratePts += 50;
				degenerate_tris = (int*) realloc(degenerate_tris, sizeof(int) * MaxNumDegeneratePts);
				degpts = (DegeneratePt*)realloc(degpts, sizeof(DegeneratePt)*MaxNumDegeneratePts);
			}
		}
	}

	if(ndegenerate_tris>0) /*we find some degenerate triangles*/
		compute_degpts_pos_tranvec_quad();
}

/*Quad mesh
Here we use the Gaussian circle and the angles of the corresponding vectors
of the vertices to find out the a.
In this routine, we linearly interpolate v0 and v1, v2 and v3, respectively,
and try to find an a such that, v0v1 and v2v3 have opposite directions.
*/
void compute_a_alongx_degptlocate(double &a, icVector2 v0, icVector2 v1, icVector2 v2, icVector2 v3,
								  icVector2 &v0v1, icVector2 &v2v3)
{
	/*use binary search*/
	/*initialization*/
	double a0, a1;
	bool orient = false;  //false -- CCW, true -- CW
	a0=0.; a1=1.;
	a = 0.5;
	v0v1 = 0.5*(v0+v1);
	v2v3 = 0.5*(v3+v2);
	double theta1, theta2, theta;
	double theta_v0, theta_v1, theta_v2, theta_v3, theta_v03;
	theta_v0 = atan2(v0.entry[1], v0.entry[0]);
	//if(theta_v0<0) theta_v0 += 2*M_PI;
	//theta_v1 = atan2(v1.entry[1], v1.entry[0]);
	//if(theta_v1<0) theta_v1 += 2*M_PI;
	//theta_v2 = atan2(v2.entry[1], v2.entry[0]);
	//if(theta_v2<0) theta_v2 += 2*M_PI;
	theta_v3 = atan2(v3.entry[1], v3.entry[0]);
	//if(theta_v3<0) theta_v3 += 2*M_PI;

	theta_v03 = theta_v3-theta_v0;
	if(theta_v03>=0)
		orient = false; //CCW
	else
		orient = true;
	if(theta_v03<-M_PI) orient = false; //CCW
	if(theta_v03>M_PI) orient = true;   //CW

	/*NOTE: we interpolate from v0 to v1
	                       from v3 to v2 */
	//normalize(v0v1);
	//normalize(v2v3);
	theta1 = atan2(v0v1.entry[1], v0v1.entry[0]); //obtain angle for v0v1;
	//if(theta1<0) theta1 += 2*M_PI;
	theta2 = atan2(v2v3.entry[1], v2v3.entry[0]); //obtain angle for v2v3;
	//if(theta2<0) theta2 += 2*M_PI;
	/*subtract the two angles*/
	theta = theta1-theta2;

	bool s_orient = false;
	while(fabs(fabs(theta)-M_PI)>1e-8 && fabs(a0-a1)>1e-9) /*if they are not opposite to each other*/
	{
		if(theta>=0)
			s_orient = false;                   //CCW
		else
			s_orient = true;                    //CW

		if(theta>M_PI) s_orient = true;         //CW
		if(theta<-M_PI) s_orient = false;       //CCW

		/*if they have the same orientation*/
		if((orient&&s_orient) || (!orient&&!s_orient))
		{
			/*we need to increase a*/
			a1 = a;
			a = (a0+a1)/2.;
		}
		else
		{
			/*we need to decrease a*/
			a0 = a;
			a = (a0+a1)/2.;
		}

		/*recalculate v0v1 and v2v3*/
		v0v1 = (1-a)*v0+a*v1;
		v2v3 = (1-a)*v3+a*v2;

		/*recompute the angles of v0v1 and v2v3*/
		theta1 = atan2(v0v1.entry[1], v0v1.entry[0]); //obtain angle for v0v1;
		if(theta1<0) theta1 += 2*M_PI;
		theta2 = atan2(v2v3.entry[1], v2v3.entry[0]); //obtain angle for v2v3;
		if(theta2<0) theta2 += 2*M_PI;
		/*subtract the two angles*/
		theta = theta1-theta2;

	}
	
	v0v1 = (1-a)*v0+a*v1;
	v2v3 = (1-a)*v3+a*v2;
}

/*Quad mesh
compute the exact position of the degenerate points 09/26/2007
*/
void compute_degpts_pos_tranvec_quad()
{
	int i;
	double x_cp, y_cp;
	ndegpts = 0;

	for(i=0; i<ndegenerate_tris; i++)
	{
		compute_onedegpt_pos_tranvec_quad(degenerate_tris[i], x_cp, y_cp);

		/*save the information to the degenerate point list*/
		degpts[ndegpts].gcx = x_cp;
		degpts[ndegpts].gcy = y_cp;
		degpts[ndegpts].degpt_index = ndegpts;
		degpts[ndegpts].type = quadmesh->quadcells[degenerate_tris[i]]->degenerate_type;
		degpts[ndegpts].Triangle_ID = degenerate_tris[i];
		quadmesh->quadcells[degenerate_tris[i]]->degpt_index = ndegpts;

		/*compute the separatrices*/
		degpts[ndegpts].nseps = 0;

		ndegpts++;
	}
}



void compute_onedegpt_pos_tranvec_quad(int cellid, double &x, double &y)
{
	/*first, we find an a along x direction, such that with this coefficient
	the interpolated vectors between v0v1 and v2v3 have opposite direction*/

	/*second, on this a, we find an b along y direction, such that the
	magnitude of v0v1 equal the magnitude of v2v3*/
	QuadCell *qc = quadmesh->quadcells[cellid];

	//QuadVertex *v00 = quadmesh->quad_verts[qc->verts[0]];
	//QuadVertex *v01 = quadmesh->quad_verts[qc->verts[1]];
	//QuadVertex *v10 = quadmesh->quad_verts[qc->verts[2]];
	//QuadVertex *v11 = quadmesh->quad_verts[qc->verts[3]];
	QuadVertex *v00 = quadmesh->quad_verts[qc->verts[0]];
	QuadVertex *v01 = quadmesh->quad_verts[qc->verts[3]];
	QuadVertex *v10 = quadmesh->quad_verts[qc->verts[1]];
	QuadVertex *v11 = quadmesh->quad_verts[qc->verts[2]];

	icVector2 v0v1, v2v3;
	double a, b;

	/*get a: the most difficult step*/
	compute_a_alongx_degptlocate(a, v00->tran_vec, v10->tran_vec, v11->tran_vec, v01->tran_vec,
		v0v1, v2v3);

	/*get b. after first step, v0v1 and v2v3 are opposite to each other*/
	if(fabs(v0v1.entry[0])>1e-8) /*use x direction*/
	{
		b = (v0v1.entry[0])/(v0v1.entry[0]-v2v3.entry[0]);
	}
	else /*use y direction*/
	{
		b = (v0v1.entry[1])/(v0v1.entry[1]-v2v3.entry[1]);
	}

	/*obtain the position*/
	x = bilinear_interpolate(a, b, v00->x, v01->x, v10->x, v11->x);
	y = bilinear_interpolate(a, b, v00->y, v01->y, v10->y, v11->y);
}

/*Quad mesh
given a quad cell index, return the x and y ranges of this regular cell
*/
void get_x_y_ranges(int id, double &xstart, double &xend, 
					double &ystart, double &yend)
{
	QuadCell *qc = quadmesh->quadcells[id];
	int i;
	QuadVertex *qv = quadmesh->quad_verts[qc->verts[0]];
	xstart = xend = qv->x;
	ystart = yend = qv->y;

	for(i=1; i<qc->nverts; i++)
	{
		if(xstart > qv->x)
			xstart = qv->x;
		if(xend < qv->x)
			xend = qv->x;
		if(ystart > qv->y)
			ystart = qv->y;
		if(yend < qv->y)
			yend = qv->y;
	}
}


bool is_in_cell(int id, double x, double y)
{
	if(id<0) return false;

	QuadCell *q = quadmesh->quadcells[id];
	double xleft=q->x_start_coord-1.e-8;
	double xright=q->x_start_coord+quadmesh->xinterval+1.e-8;
	double ybuttom=q->y_start_coord-1.e-8;
	double yupper=q->y_start_coord+quadmesh->yinterval+1.e-8;
	if((x>=xleft && x<=xright)
		&&(y>=ybuttom && y<=yupper))
		return true;
	return false;

}

bool is_in_reg_cell(int id, double x, double y)
{
	if(id<0) return false;
	int i=(x-quadmesh->xstart)/quadmesh->xinterval;
	int j=(y-quadmesh->ystart)/quadmesh->yinterval;

	if(id==(j*(quadmesh->XDIM-1)+i))
		return true;
	return false;
}

double bilinear_interpolate(double a, double b, 
						  double f00, double f01, double f10, double f11)
{
	return (f00*(1-a)*(1-b)+f10*a*(1-b)+f01*(1-a)*b+f11*a*b);
}
/*Quad mesh
Obtain the tensor in a quad cell using bilinear interpolation
NOTE: we assume that the given point always falls in this cell
09/26/2007
*/
void compute_tensor_at_quad(int face, double x, double y, icMatrix2x2 &ten)
{
	//double xstart, xend, ystart, yend;
	//get_x_y_ranges(face, xstart, xend, ystart, yend);

	QuadCell *qc = quadmesh->quadcells[face];

	/*get the x coeff and y coeff*/
	double a = (x-qc->x_start_coord)/quadmesh->xinterval;
	double b = (y-qc->y_start_coord)/quadmesh->yinterval;

	if(fabs(a)<1e-6)
		a = 0;
	if(fabs(b)<1e-6)
		b = 0;

	//if(a<0 ||b<0)
	//{
	//	int test = 0;
	//}

	/*obtain the vertices of this cell in the order of
	  v00 v01
	  v10 v11
	*/
	QuadVertex *v00 = quadmesh->quad_verts[qc->verts[0]];
	QuadVertex *v01 = quadmesh->quad_verts[qc->verts[3]];
	QuadVertex *v10 = quadmesh->quad_verts[qc->verts[1]];
	QuadVertex *v11 = quadmesh->quad_verts[qc->verts[2]];

	/*the all the components of the interpolated tensor, respectively*/
	ten.entry[0][0] = bilinear_interpolate(a, b, v00->Jacobian.entry[0][0], v01->Jacobian.entry[0][0],
		v10->Jacobian.entry[0][0], v11->Jacobian.entry[0][0]);
	ten.entry[0][1] = bilinear_interpolate(a, b, v00->Jacobian.entry[0][1], v01->Jacobian.entry[0][1],
		v10->Jacobian.entry[0][1], v11->Jacobian.entry[0][1]);
	ten.entry[1][0] = bilinear_interpolate(a, b, v00->Jacobian.entry[1][0], v01->Jacobian.entry[1][0],
		v10->Jacobian.entry[1][0], v11->Jacobian.entry[1][0]);
	ten.entry[1][1] = bilinear_interpolate(a, b, v00->Jacobian.entry[1][1], v01->Jacobian.entry[1][1],
		v10->Jacobian.entry[1][1], v11->Jacobian.entry[1][1]);

}

/*the following, we try to implement tensor field tracing 09/20/2007*/

#include "LocalTracing.h"
#include "Numerical.h"

extern int MaxNumTrajectories;
extern int MaxNumLinesegsPerTraj;

//extern Trajectory *trajectories2;
extern LineSeg **trajectories;                 //trajectories' list
extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory
extern int globalface;
extern DP htry;

icVector2 tenline_dir_loc;       //record the local direction of current tensor line 09/20/2007
icVector2 tenline_dir_global;    //record the global direction of current tensor line 09/20/2007
icVector2 tenline_dir_global_p;  //record the previous global direction of current tensor line 09/20/2007

void ten_trace(int face_id, double x, double y)
{
	int i;
	int flag = 0;
	double globalp[2];
	int pre_face, cur_face;

	//bool loop_flag = false;  ////for detect a closed loop from the beginning triangle
	//int *triangles = (int*)malloc(sizeof(int) * NUMTRACINGTRIANGLE);
	//int num_triangles = 0;

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

	/*we need to pick a direction for it first*/
	double alpha[3];
	Get2DBarycentricFacters(face_id, x, y, alpha);

	icMatrix2x2 ten;
	compute_tensor_at(face_id, alpha, ten);

	double evalues[2] = {0.};
	icVector2 ev[2];
	cal_eigenval_matrix2x2(ten, evalues);
	cal_eigenvectors(ten, evalues, ev);
	tenline_dir_global = ev[0];
	tenline_dir_loc.entry[0] = dot(tenline_dir_global, Object.flist[face_id]->LX);
	tenline_dir_loc.entry[1] = dot(tenline_dir_global, Object.flist[face_id]->LY);


	/*along ev[0] direction*/

	for(i = 0; i < 10*NUMTRACINGTRIANGLE; i++)
	{

		if(cur_face == -1) /*reach boundary*/
		{
			return;
		}

		pre_face = cur_face;
		cur_face = trace_ten_in_one_tri(cur_face, globalp, flag); 

		if(flag == 1 || flag == 2 || pre_face == cur_face/* || loop_flag */) 
		{
			break;
		}
	}

	/*along negative ev[0] direction*/
	pre_face = cur_face = face_id;
	globalp[0] = x;   globalp[1] = y;
	flag = 0;
	tenline_dir_global = -ev[0];
	tenline_dir_loc.entry[0] = dot(tenline_dir_global, Object.flist[face_id]->LX);
	tenline_dir_loc.entry[1] = dot(tenline_dir_global, Object.flist[face_id]->LY);

	for(i = 0; i < 10*NUMTRACINGTRIANGLE; i++)
	{

		if(cur_face == -1) /*reach boundary*/
		{
			return;
		}

		pre_face = cur_face;
		cur_face = trace_ten_in_one_tri(cur_face, globalp, flag); 

		if(flag == 1 || flag == 2 || pre_face == cur_face/* || loop_flag */) 
		{
			return;
		}

	}
}


////calculate the trajectory in a single triangle
int trace_ten_in_one_tri(int &face_id, double globalp[2], int &flag)
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

	CurvePoints *temp_point_list = (CurvePoints*) malloc(sizeof(CurvePoints) * 100);
	int NumPoints = 0;

	////initialize
	VP.entry[0] = globalp[0] - Object.vlist[face->verts[0]]->x;
	VP.entry[1] = globalp[1] - Object.vlist[face->verts[0]]->y;

	pre_point[0] = cur_point[0] = dot(VP, face->LX);
	pre_point[1] = cur_point[1] = dot(VP, face->LY);

	vert0[0] = Object.vlist[face->verts[0]]->x;   ////for update the global point
	vert0[1] = Object.vlist[face->verts[0]]->y;

	globalface = face_id;

	/*we need to obtain the direction of current tensor line.
	The good thing is that we can reuse the previous obtained direction knowing
	there may be some errors.
	*/

	///*obtain the global direction of current tensor line 09/20/2007*/
	
	//tenline_dir_global.entry[0] = globalp[0] - temp_point_list[NumPoints].gpx;
	//tenline_dir_global.entry[1] = globalp[1] - temp_point_list[NumPoints].gpy;


	////////////////////////////////////////////////////
    for(i = 0; i < 100; i++)
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
			if(compute_next_pt_tensor(pre_point, cur_point, face_id, alpha))

			{
				////update the global point
				face = Object.flist[face_id];

				globalv = cur_point[0] * face->LX + cur_point[1] * face->LY;

				globalp[0] = vert0[0] + globalv.entry[0];
				globalp[1] = vert0[1] + globalv.entry[1];

				/*obtain the local direction of current tensor line 09/20/2007*/
				tenline_dir_loc.entry[0] = cur_point[0]-pre_point[0];
				tenline_dir_loc.entry[1] = cur_point[1]-pre_point[1];

				tenline_dir_global_p = tenline_dir_global;

				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = globalp[0] - temp_point_list[NumPoints-1].gpx;
				tenline_dir_global.entry[1] = globalp[1] - temp_point_list[NumPoints-1].gpy;

			}

			else{  ////the curve reach a singularity
				flag = 1;

				////Store the record into global line segment array
                
				if(!StoreToGlobalList(temp_point_list, NumPoints))
				{
					////Not enough memory
					flag = 2;
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
            
			get_next_tri(face_id, pre_point, cur_point, t, PassVertornot, alpha);

			////update the global point here
			if(PassVertornot > 0) /* Make a big change on 01/30/07 */
			{
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

			/*obtain the local direction of current tensor line 09/20/2007*/
			tenline_dir_loc.entry[0] = cur_point[0]-pre_point[0];
			tenline_dir_loc.entry[1] = cur_point[1]-pre_point[1];

			/*obtain the global direction of current tensor line 09/20/2007*/
			tenline_dir_global.entry[0] = globalp[0] - temp_point_list[NumPoints-1].gpx;
			tenline_dir_global.entry[1] = globalp[1] - temp_point_list[NumPoints-1].gpy;
			
			tenline_dir_loc.entry[0] = dot(tenline_dir_global, Object.flist[face_id]->LX);
			tenline_dir_loc.entry[1] = dot(tenline_dir_global, Object.flist[face_id]->LY);

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
				   flag = 2;
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

/*************************************************************
Runge Kutta integrator driver for local tracing
*************************************************************/
void localderive_tensor(const DP t, Vec_I_DP &y, Vec_O_DP &dydx)
{
	double alpha[3];
	Get2DBarycentricFacters(globalface, y[0], y[1], alpha);

	/*first, we need to get the tensor at current location*/
	icMatrix2x2 ten;
	compute_tensor_at(globalface, alpha, ten);

	/*second, we compute its major eigen vector*/
	double evalues[2] = {0.};
	icVector2 ev[2];
	cal_eigenval_matrix2x2(ten, evalues);
	cal_eigenvectors(ten, evalues, ev);

	/*we also need to judge whether the direction of the obtained vector and
	the direction of the current tensor line are the same using dot product*/
	double re = dot(tenline_dir_global, ev[0]);
	if(re<0) 
		ev[0] = -ev[0];

	/*transfer back to local vector*/
	icVector2 loc_ev;
	Face *face = Object.flist[globalface];
	loc_ev.entry[0] = dot(ev[0], face->LX);
	loc_ev.entry[1] = dot(ev[0], face->LY);

	dydx[0] = loc_ev.entry[0];
	dydx[1] = loc_ev.entry[1];

	//if(dot(tenline_dir_loc, loc_ev)<0)
	//	dydx[0] = -loc_ev.entry[0];
	//	dydx[1] = -loc_ev.entry[1];

}


bool compute_next_pt_tensor(double first[2], double second[2], int &face_id, double alpha[3])
{

	/*we need to interpolate the tensor and calculate the major eigen vector here*/

	//icVector2 VecAtPoint = GetVectorAtPoints(face_id, alpha);

	/*the following can be optimized 09/20/2007*/
	icMatrix2x2 ten;
	compute_tensor_at(face_id, alpha, ten);

	double evalues[2] = {0.};
	icVector2 ev[2];
	cal_eigenval_matrix2x2(ten, evalues);
	cal_eigenvectors(ten, evalues, ev);

	if(length(ev[0]) < 1e-18 ) 
		return false; /*reach singularity*/
	
	////calling Runge Kutta to get the next point
    const int N=2;
    int i,j;
    DP eps,hdid,hnext, t = 1.0;

    Vec_DP by(N),dydx(N),dysav(N),ysav(N),yscal(N);

	ysav[0] = first[0];
	ysav[1] = first[1];

    localderive_tensor(t, ysav, dysav);

    for (i=0;i<N;i++) yscal[i]=1;

	/*--------------------------------------------------------------------*/
	////calling adaptive stepsize runge-kutta method here to get next step
	for (i=0;i<10;i++) {
		eps=exp(-DP(i+1));
		t = 1.0;
		for (j=0;j<N;j++) {
			by[j]=ysav[j];
			dydx[j]=dysav[j];
		}
	
		rkqs(by,dydx,t,htry,eps,yscal,hdid,hnext,localderive_tensor);

	}

	////Set the new stepsize, note that this stepsize will affect the limit cycle detection
	htry = hnext;
	if(hnext >= 2.)
		htry = 2.;


	/////store some important information
	second[0] = by[0];
	second[1] = by[1];

	return true;
}


/*
get next triangle when tracing
*/
void get_next_tri(int &face_id, double pre[2], double cur[2], double param_t[2], 
					 int &PassVertornot, double alpha[3])
{
	int which_edge = -1;

	int prev_face_id = face_id;

	Face *prev_face = Object.flist[face_id];

	Vertex *vert = NULL;

	PassVertornot = 0;
	
	////We should put pass vertex testing here before testing crossing edge
	//double re = dot(vert

	cross_AVertex(face_id, cur, pre, PassVertornot);

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


/*
New routine for crossing vertex testing 4/30/06
*/
void cross_AVertex(int &face_id, double cur_p[2], double pre_p[2], int &passornot)
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
		if(fabs(pending) < 1e-8) ////passing the vertex
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

			double re = dot(Object.vlist[crossVert]->major, 
				tenline_dir_global);
			if(re>=0)
                TriangleThroughVertex(crossVert, newtriangleid, 0);
			else
                TriangleThroughVertex(crossVert, newtriangleid, 1);

			face_id = newtriangleid;	
			passornot = i+1;
			return;
		}
	}

	passornot = 0;
}



/*********************************************************************************/


/*the followings implement the tensor field tracing under the quad mesh
NOTE: currently, we consider the local tracing under the global
09/26/2007*/
extern void get_tensor(double, double, double []);

void ten_trace_quad(int cell_id, double x, double y, int type)
{
	int i;
	int flag = 0;
	double globalp[2];
	int pre_face, cur_face;
	double t[4];

	pre_face = cur_face = cell_id;
	globalp[0] = x;   globalp[1] = y;
	euler_stepsize = quadmesh->xinterval/8.;

	predict_stepsize = quadmesh->xinterval/2.;

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

	/*we need to pick a direction for it first*/

	icMatrix2x2 ten;
	//compute_tensor_at_quad(cell_id, x, y, ten);
	get_tensor(x, y, t);
	ten.entry[0][0]=t[0];
	ten.entry[0][1]=t[1];
	ten.entry[1][0]=t[2];
	ten.entry[1][1]=t[3];

	double evalues[2] = {0.};
	icVector2 ev[2], startdir;
	cal_eigen_vector_sym(ten, ev);

	if(type == 1) /*use minor eigen vector field*/
		ev[0] = ev[1];

	tenline_dir_global_p = tenline_dir_global = startdir = ev[0];  /*obtain the major eigen vector*/

	/*along positive ev[0] direction*/

	int NUMTRACINGCELLS = (int)(quadmesh->nfaces/2);
	for(i = 0; i < NUMTRACINGCELLS; i++)
	{

		if(cur_face == -1 || cur_face>=quadmesh->nfaces) /*reach boundary*/
		{
			break;
		}

		/*it reaches any boundary, we should stop as well*/
		if(globalp[0]>quadmesh->xend-1.e-8||globalp[0]<quadmesh->xstart+1.e-8
			||globalp[1]>quadmesh->yend-1.e-8||globalp[1]<quadmesh->ystart+1.e-8)
			break;

		pre_face = cur_face;
		cur_face = trace_ten_in_one_tri_quad(cur_face, globalp, flag, type); 

		if(flag == 1 || flag == 2 || pre_face == cur_face/* || loop_flag */) 
		{
			break;
		}
	}

	/*along negative ev[0] direction*/
	pre_face = cur_face = cell_id;
	globalp[0] = x;   globalp[1] = y;
	flag = 0;
	tenline_dir_global_p = tenline_dir_global = -startdir;


	for(i = 0; i < NUMTRACINGCELLS; i++)
	{

		if(cur_face == -1|| cur_face>=quadmesh->nfaces) /*reach boundary*/
		{
			return;
		}

		pre_face = cur_face;
		cur_face = trace_ten_in_one_tri_quad(cur_face, globalp, flag, type); 

		if(flag == 1 || flag == 2 || pre_face == cur_face/* || loop_flag */) 
		{
			return;
		}

	}
}


////calculate the trajectory in a single triangle
int trace_ten_in_one_tri_quad(int &face_id, double globalp[2], int &flag, int type)
{
	int i;
	double pre_point[2];
	if(face_id < 0)
		return -1;

	QuadCell *face = quadmesh->quadcells[face_id];

	QuadCell *pre_f = face;
	
	////Temporary curve point array

	CurvePoints *temp_point_list = (CurvePoints*) malloc(sizeof(CurvePoints) * 150);
	int NumPoints = 0;

	////initialize

	/*the tracing will be performed under the global frame*/
	globalface = face_id;

	pre_point[0] = globalp[0];
	pre_point[1] = globalp[1];

	///*obtain the global direction of current tensor line 09/20/2007*/

	////////////////////////////////////////////////////
    for(i = 0; i < 150; i++)
	{
		////1. judge whether current point still falls in current cell

		////2. if current point is inside current triangle
		if( is_in_cell(face_id, globalp[0], globalp[1]))
		{
			////store the point into the temp curve points list

			temp_point_list[NumPoints].gpx = globalp[0];
			temp_point_list[NumPoints].gpy = globalp[1];
			temp_point_list[NumPoints].triangleid = face->index;  
			NumPoints++;

			pre_point[0] = globalp[0];
			pre_point[1] = globalp[1];
			
			/*change to use other integration scheme 07/09/07*/
			//if(compute_next_pt_tensor_quad_global(pre_point, globalp, face_id))
			//if(get_nextpt_2ndeuler_ten_quad(pre_point, globalp, face_id, type))
			if(get_nextpt_RK23_ten_quad(pre_point, globalp, face_id, type))
			//if(get_nextpt_RK45_ten_quad(pre_point, globalp, face_id, type))
			{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = globalp[0] - pre_point[0];
				tenline_dir_global.entry[1] = globalp[1] - pre_point[1];
			}

			else{  ////the curve reach a singularity/degenerate point
				flag = 1;

				////Store the record into global line segment array
                
				if(!StoreToGlobalList(temp_point_list, NumPoints))
				{
					////Not enough memory
					flag = 2;
					free(temp_point_list);
					return face_id;
				}

				free(temp_point_list);

				return face_id;
			}
		}

		////3. if the point is out of current cell
		else{
			/*!!!!!!need to judge which cell it will enter!!!!!*/
			int PassVertornot = 0;
			//get_next_cell(face_id, pre_point, globalp, PassVertornot, type);
			get_next_cell_2(face_id, pre_point, globalp, PassVertornot, type);

			if(PassVertornot>0)  /*cross a vertex*/
			{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global.entry[0] = pre_point[0] - globalp[0];
				tenline_dir_global.entry[1] = pre_point[1] - globalp[1];

				/**/
				temp_point_list[NumPoints].gpx = globalp[0];
				temp_point_list[NumPoints].gpy = globalp[1];
				temp_point_list[NumPoints].triangleid = face->index;  ////cause problem 05/25/05
				NumPoints++;

				temp_point_list[NumPoints].gpx = pre_point[0];
				temp_point_list[NumPoints].gpy = pre_point[1];
				temp_point_list[NumPoints].triangleid = face_id;  ////cause problem 05/25/05
				NumPoints++;
			}
			else{
				/*obtain the global direction of current tensor line 09/20/2007*/
				tenline_dir_global_p = tenline_dir_global;

				tenline_dir_global.entry[0] = globalp[0] - pre_point[0];
				tenline_dir_global.entry[1] = globalp[1] - pre_point[1];

				////Add the intersection point to the temporary points' list
				temp_point_list[NumPoints].gpx = globalp[0];
				temp_point_list[NumPoints].gpy = globalp[1];
				temp_point_list[NumPoints].triangleid = face->index;  ////cause problem 05/25/05
				NumPoints++;
			}

			if(NumPoints > 1){
 				////Store the record into global line segment array
               if(!StoreToGlobalList(temp_point_list, NumPoints))
			   {   ////Not enough memory
				   flag = 2;
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

/*************************************************************
Runge Kutta integrator driver for global tracing in quad mesh
*************************************************************/
void derive_tensor_quad_global(const DP t, Vec_I_DP &y, Vec_O_DP &dydx)
{
	/*first, we need to get the tensor at current location*/
	icMatrix2x2 ten;
	//compute_tensor_at_quad(globalface, y[0], y[1], ten);
	double tt[4];
	get_tensor(y[0], y[1], tt);
	ten.entry[0][0]=tt[0];
	ten.entry[0][1]=tt[1];
	ten.entry[1][0]=tt[2];
	ten.entry[1][1]=tt[3];

	/*second, we compute its major eigen vector*/
	icVector2 ev[2];
	cal_eigen_vector_sym(ten, ev);

	/*we also need to judge whether the direction of the obtained vector and
	the direction of the current tensor line are the same using dot product*/
	double re = dot(tenline_dir_global, ev[0]);
	if(re<0) 
		ev[0] = -ev[0];

	dydx[0] = ev[0].entry[0];
	dydx[1] = ev[1].entry[1];
}


bool compute_next_pt_tensor_quad_global(double first[2], double second[2], int &face_id)
{

	/*we need to interpolate the tensor and calculate the major eigen vector here*/

	/*the following can be optimized 09/20/2007*/
	icMatrix2x2 ten;
	compute_tensor_at_quad(face_id, first[0], first[1], ten);

	icVector2 ev[2];
	cal_eigen_vector_sym(ten, ev);

	if(length(ev[0]) < 1e-18 ) 
		return false; /*reach singularity*/
	
	////calling Runge Kutta to get the next point
    const int N=2;
    int i,j;
    DP eps,hdid,hnext, t = 1.0;

    Vec_DP by(N),dydx(N),dysav(N),ysav(N),yscal(N);

	ysav[0] = first[0];
	ysav[1] = first[1];

    derive_tensor_quad_global(t, ysav, dysav);

    for (i=0;i<N;i++) yscal[i]=1;

	/*--------------------------------------------------------------------*/
	////calling adaptive stepsize runge-kutta method here to get next step
	for (i=0;i<15;i++) {
		eps=exp(-DP(i+1));
		t = 1.0;
		for (j=0;j<N;j++) {
			by[j]=ysav[j];
			dydx[j]=dysav[j];
		}
	
		rkqs(by,dydx,t,htry,eps,yscal,hdid,hnext,derive_tensor_quad_global);

	}

	////Set the new stepsize, note that this stepsize will affect the limit cycle detection
	htry = hnext;
	if(hnext >= 1.1)
		htry = 1.1;


	/////store some important information
	second[0] = by[0];
	second[1] = by[1];

	return true;
}


/*
The following routines help determine the next quad cell during trace
*/

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

extern int GetIntersection2(double PointA[2], double PointB[2], double PointC[2], double PointD[2], double t[2]);

/*
New routine for crossing vertex testing 
This is suitable for the tracing under global frame 
09/27/2007
*/
bool cross_vertex_ten_quad(int &face_id, double cur_p[2], double pre_p[2], int &passornot, int type)
{
	int i;
	double vert[2];
	double max_alpha ;
    int newtriangleid = 0;
	int crossVert;
	QuadCell *face = quadmesh->quadcells[face_id];
	QuadVertex *v;

	double A, B, C, pending;
	A = pre_p[1] - cur_p[1];
	B = cur_p[0] - pre_p[0];
	C = (pre_p[0]*cur_p[1] - cur_p[0]*pre_p[1]);

	for(i = 0; i < face->nverts; i++)
	{
		v = quadmesh->quad_verts[face->verts[i]];
		vert[0] = v->x;
		vert[1] = v->y;
		pending = A*vert[0] + B*vert[1] + C;
	    ////We also need to make sure that the vertex is between 'pre' and 'cur' points
		//if(fabs(pending) == 0.0) ////passing the vertex
		if(fabs(pending) <=1.e-8) ////passing the vertex
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

			get_cell_through_ver(crossVert, newtriangleid, type);
			if(newtriangleid <0 || newtriangleid>=quadmesh->nfaces) 
				return false;
			face_id = newtriangleid;
			passornot = i+1;
			cur_p[0] = quadmesh->quad_verts[crossVert]->x;
			cur_p[1] = quadmesh->quad_verts[crossVert]->y;
			//pre_p[0] = quadmesh->quad_verts[crossVert]->x;
			//pre_p[1] = quadmesh->quad_verts[crossVert]->y;

			/*make it off the vertex a little bit 10/02/2007*/
			icVector2 tvec;
			QuadVertex *tv = quadmesh->quad_verts[crossVert];
			if(type==0)
				tvec=tv->major;
			else
				tvec=tv->minor;
			normalize(tvec);
			double re=dot(tvec, tenline_dir_global);
			if(re<0) tvec=-tvec;

			/*  try to avoid the vertex  */
			pre_p[0] = tv->x+0.01*quadmesh->xinterval*tvec.entry[0];
			pre_p[1] = tv->y+0.01*quadmesh->xinterval*tvec.entry[1];

			int test_cell=get_cellID_givencoords(pre_p[0], pre_p[1]);

			//if(test_cell != newtriangleid)
			//{
			//	int test=0;
			//}

			return true;
		}
	}

	passornot = 0;
	return false;
}



/*get the cell id that the vector will enter
09/27/2006
*/
void get_cell_through_ver(int vertid, int &cell, int type)
{
	QuadVertex *v = quadmesh->quad_verts[vertid];
	icVector2 vec; 
	if(type == 0)
	    vec = v->major;
	else
		vec = v->minor;

	double re=dot(vec, tenline_dir_global);
	if(re<0) vec=-vec;
	normalize(vec);

	//double x, y;

	double x = v->x+0.02*quadmesh->xinterval*vec.entry[0];
	double y = v->y+0.02*quadmesh->xinterval*vec.entry[1];

	/*judge which cell contains (x, y)*/
	//int i;
	//for(i=0; i<v->ncells; i++)
	//{
	//	if(v->cells[i]->index == cell)
	//		continue;
	//	if(is_in_cell(v->cells[i]->index, x, y))
	//	{
	//		cell = v->cells[i]->index;
	//		return;
	//	}
	//}
	//cell = get_cellID_picking(x, y);
	cell = get_cellID_givencoords(x, y);
}


/*get the next cell the tensor line will enter next 09/27/2007*/
void get_next_cell(int &face_id, double pre[2], double cur[2], 
					 int &PassVertornot, int type)
{
	/*for horizontal cases*/
	if(fabs(cur[1]-pre[1])<1e-9)
	{
		if(cur[0]>pre[0] && cur[0]<=quadmesh->xend-1.e-9) /*move to the right cell*/
		{
			cur[0]=quadmesh->quadcells[face_id]->x_start_coord+quadmesh->xinterval+1.1e-8;
			face_id++;
			return;
		}
		else if(cur[0]>pre[0] && cur[0]>quadmesh->xend-1.e-9) /*out of mesh*/
		{
			face_id=-1;
			return;
		}

		else if(cur[0]<pre[0] && cur[0]>=quadmesh->xstart+1.e-9)/*move to the left cell*/
		{
			cur[0]=quadmesh->quadcells[face_id]->x_start_coord-1.1e-8;
			face_id--;
			return;
		}
		else if(cur[0]<pre[0] && cur[0]<quadmesh->xstart+1.e-9)/*out of mesh*/
		{
			face_id=-1;
			return;
		}
	}

	/*for vertical cases*/
	if(fabs(cur[0]-pre[0])<1e-9)
	{
		if(cur[1]>pre[1] && cur[1]<=quadmesh->yend-1.e-9) /*move to the upper cell*/
		{
			cur[1]=quadmesh->quadcells[face_id]->y_start_coord+quadmesh->yinterval+1.1e-8;
			face_id+=(quadmesh->XDIM-1);
			return;
		}

		else if(cur[1]>pre[1] && cur[1]>quadmesh->yend-1.e-9)
		{
			face_id = -1;
			return;
		}

		else if(cur[1]<pre[1] && cur[1]>=quadmesh->ystart+1.e-9) /*move to the lower cell*/
		{
			cur[1]=quadmesh->quadcells[face_id]->y_start_coord-1.1e-8;
			face_id-=(quadmesh->XDIM-1);
			return;
		}
		
		else if(cur[1]<pre[1] && cur[1]<quadmesh->ystart+1.e-9)
		{
			face_id = -1;
			return;
		}
	}

	int pre_f = face_id;
	if(cross_vertex_ten_quad(face_id, cur, pre, PassVertornot, type))
		return;

	face_id = pre_f;
	/*get the angle of current line L*/
	icVector2 dir;
	dir.entry[0] = cur[0]-pre[0];
	dir.entry[1] = cur[1]-pre[1];
	double theta = atan2(dir.entry[1], dir.entry[0]);

	double t[2];
	double C[2], D[2];
	QuadCell *f = quadmesh->quadcells[face_id];

	if(theta>=-M_PI/2 && theta<=M_PI/2.)
	{
		/*calculate the intersection between edge e1 and L*/
		C[0] = quadmesh->quad_verts[f->edges[1]->verts[0]]->x;
		C[1] = quadmesh->quad_verts[f->edges[1]->verts[0]]->y;
		D[0] = quadmesh->quad_verts[f->edges[1]->verts[1]]->x;
		D[1] = quadmesh->quad_verts[f->edges[1]->verts[1]]->y;
		GetIntersection2(pre, cur, C, D, t);

		if(t[1]>0 && t[1]<1)
		{
			face_id = face_id+1;
			/*calculate the intersection*/
			cur[0]=C[0]+t[1]*(D[0]-C[0]);
			cur[1]=C[1]+t[1]*(D[1]-C[1]);
		}
		else
		{
			if(theta>0)
			{
				/*calculate the intersection between edge e2 and L*/
				C[0] = quadmesh->quad_verts[f->edges[2]->verts[0]]->x;
				C[1] = quadmesh->quad_verts[f->edges[2]->verts[0]]->y;
				D[0] = quadmesh->quad_verts[f->edges[2]->verts[1]]->x;
				D[1] = quadmesh->quad_verts[f->edges[2]->verts[1]]->y;
				GetIntersection2(pre, cur, C, D, t);
				
				face_id = face_id+(quadmesh->XDIM-1);
				/*calculate the intersection*/
				cur[0]=C[0]+t[1]*(D[0]-C[0]);
				cur[1]=C[1]+t[1]*(D[1]-C[1]);
			}
			else
			{
				/*calculate the intersection between edge e0 and L*/
				C[0] = quadmesh->quad_verts[f->edges[0]->verts[0]]->x;
				C[1] = quadmesh->quad_verts[f->edges[0]->verts[0]]->y;
				D[0] = quadmesh->quad_verts[f->edges[0]->verts[1]]->x;
				D[1] = quadmesh->quad_verts[f->edges[0]->verts[1]]->y;
				GetIntersection2(pre, cur, C, D, t);
				
				face_id = face_id-(quadmesh->XDIM-1);
				/*calculate the intersection*/
				cur[0]=C[0]+t[1]*(D[0]-C[0]);
				cur[1]=C[1]+t[1]*(D[1]-C[1]);
			}
		}
	}
	else
	{
		/*calculate the intersection between edge e3 and L*/
		C[0] = quadmesh->quad_verts[f->edges[3]->verts[0]]->x;
		C[1] = quadmesh->quad_verts[f->edges[3]->verts[0]]->y;
		D[0] = quadmesh->quad_verts[f->edges[3]->verts[1]]->x;
		D[1] = quadmesh->quad_verts[f->edges[3]->verts[1]]->y;
		GetIntersection2(pre, cur, C, D, t);
		
		if(t[1]>0 && t[1]<1)
		{
			face_id = face_id-1;
			/*calculate the intersection*/
			cur[0]=C[0]+t[1]*(D[0]-C[0]);
			cur[1]=C[1]+t[1]*(D[1]-C[1]);
		}
		else
		{
			if(theta>0)
			{
				/*calculate the intersection between edge e2 and L*/
				C[0] = quadmesh->quad_verts[f->edges[2]->verts[0]]->x;
				C[1] = quadmesh->quad_verts[f->edges[2]->verts[0]]->y;
				D[0] = quadmesh->quad_verts[f->edges[2]->verts[1]]->x;
				D[1] = quadmesh->quad_verts[f->edges[2]->verts[1]]->y;
				GetIntersection2(pre, cur, C, D, t);
				
				face_id = face_id+(quadmesh->XDIM-1);
				/*calculate the intersection*/
				cur[0]=C[0]+t[1]*(D[0]-C[0]);
				cur[1]=C[1]+t[1]*(D[1]-C[1]);
			}
			else
			{
				/*calculate the intersection between edge e2 and L*/
				C[0] = quadmesh->quad_verts[f->edges[0]->verts[0]]->x;
				C[1] = quadmesh->quad_verts[f->edges[0]->verts[0]]->y;
				D[0] = quadmesh->quad_verts[f->edges[0]->verts[1]]->x;
				D[1] = quadmesh->quad_verts[f->edges[0]->verts[1]]->y;
				GetIntersection2(pre, cur, C, D, t);
				
				face_id = face_id-(quadmesh->XDIM-1);
				/*calculate the intersection*/
				cur[0]=C[0]+t[1]*(D[0]-C[0]);
				cur[1]=C[1]+t[1]*(D[1]-C[1]);
			}
		}
	}
}


/*get the next cell the tensor line will enter next 09/27/2007*/
void get_next_cell_2(int &face_id, double pre[2], double cur[2], 
					 int &PassVertornot, int type)
{

	/*for horizontal cases*/
	if(fabs(cur[1]-pre[1])<1e-8)
	{
		if(cur[0]>pre[0] && cur[0]<=quadmesh->xend-1.e-8) /*move to the right cell*/
		{
			cur[0]=quadmesh->quadcells[face_id]->x_start_coord+quadmesh->xinterval+1.1e-8;
			face_id++;
			return;
		}

		else if(cur[0]>pre[0] && cur[0]>quadmesh->xend-1.e-8) /*out of mesh*/
		{
			face_id=-1;
			return;
		}

		else if(cur[0]<pre[0] && cur[0]>=quadmesh->xstart+1.e-8)/*move to the left cell*/
		{
			cur[0]=quadmesh->quadcells[face_id]->x_start_coord-1.1e-8;
			face_id--;
			return;
		}

		else if(cur[0]<pre[0] && cur[0]<quadmesh->xstart+1.e-8)/*out of mesh*/
		{
			face_id=-1;
			return;
		}
	}

	/*for vertical cases*/
	if(fabs(cur[0]-pre[0])<1e-8)
	{
		if(cur[1]>pre[1] && cur[1]<=quadmesh->yend-1.e-8) /*move to the upper cell*/
		{
			cur[1]=quadmesh->quadcells[face_id]->y_start_coord+quadmesh->yinterval+1.1e-8;
			face_id+=(quadmesh->XDIM-1);
			return;
		}

		else if(cur[1]>pre[1] && cur[1]>quadmesh->yend-1.e-8)
		{
			face_id = -1;
			return;
		}

		else if(cur[1]<pre[1] && cur[1]>=quadmesh->ystart+1.e-8) /*move to the lower cell*/
		{
			cur[1]=quadmesh->quadcells[face_id]->y_start_coord-1.1e-8;
			face_id-=(quadmesh->XDIM-1);
			return;
		}
		
		else if(cur[1]<pre[1] && cur[1]<quadmesh->ystart+1.e-8)
		{
			face_id = -1;
			return;
		}
	}

	int pre_f = face_id;
	if(cross_vertex_ten_quad(face_id, cur, pre, PassVertornot, type))
		return;

	face_id = pre_f;

	double t[2];
	double C[2], D[2];
	QuadCell *f = quadmesh->quadcells[face_id];

	if(cur[0]<f->x_start_coord) /*left */
	{
		/*calculate the intersection between edge e3 and L*/
		C[0] = quadmesh->quad_verts[f->verts[0]]->x;
		C[1] = quadmesh->quad_verts[f->verts[0]]->y;
		D[0] = quadmesh->quad_verts[f->verts[3]]->x;
		D[1] = quadmesh->quad_verts[f->verts[3]]->y;
		if(GetIntersection2(pre, cur, C, D, t)==1)
		{
			face_id = face_id-1;
			/*calculate the intersection*/
			cur[0]=C[0]+t[1]*(D[0]-C[0]);
			//cur[0]=f->x_start_coord-1.e-8;
			cur[1]=C[1]+t[1]*(D[1]-C[1]);

			return;
		}
	}
	else if(cur[0]>f->x_start_coord+quadmesh->xinterval)  /*right*/
	{
		/*calculate the intersection between edge e2 and L*/
		C[0] = quadmesh->quad_verts[f->verts[1]]->x;
		C[1] = quadmesh->quad_verts[f->verts[1]]->y;
		D[0] = quadmesh->quad_verts[f->verts[2]]->x;
		D[1] = quadmesh->quad_verts[f->verts[2]]->y;
		
		if(GetIntersection2(pre, cur, C, D, t)==1)
		{
			face_id = face_id+1;
			/*calculate the intersection*/
			cur[0]=C[0]+t[1]*(D[0]-C[0]);
			//cur[0]=f->x_start_coord+quadmesh->xinterval+1.e-8;
			cur[1]=C[1]+t[1]*(D[1]-C[1]);

			return;
		}
	}


	if(cur[1]<f->y_start_coord) /*bottom */
	{
		/*calculate the intersection between edge e0 and L*/
		C[0] = quadmesh->quad_verts[f->verts[0]]->x;
		C[1] = quadmesh->quad_verts[f->verts[0]]->y;
		D[0] = quadmesh->quad_verts[f->verts[1]]->x;
		D[1] = quadmesh->quad_verts[f->verts[1]]->y;
		
		if(GetIntersection2(pre, cur, C, D, t)==1)
		{
			face_id = face_id-quadmesh->XDIM+1;
			/*calculate the intersection*/
			cur[0]=C[0]+t[1]*(D[0]-C[0]);
			cur[1]=C[1]+t[1]*(D[1]-C[1]);
			//cur[1]=f->y_start_coord-1.e-8;

			return;
		}
	}
	else if(cur[1]>f->y_start_coord+quadmesh->yinterval)  /*upper*/
	{
		/*calculate the intersection between edge e2 and L*/
		C[0] = quadmesh->quad_verts[f->verts[2]]->x;
		C[1] = quadmesh->quad_verts[f->verts[2]]->y;
		D[0] = quadmesh->quad_verts[f->verts[3]]->x;
		D[1] = quadmesh->quad_verts[f->verts[3]]->y;
		
		if(GetIntersection2(pre, cur, C, D, t)==1)
		{
			face_id = face_id+quadmesh->XDIM-1;
			/*calculate the intersection*/
			cur[0]=C[0]+t[1]*(D[0]-C[0]);
			cur[1]=C[1]+t[1]*(D[1]-C[1]);
			//cur[1]=f->y_start_coord+quadmesh->yinterval+1.e-8;

			return;
		}
	}

	face_id = -1;
}



bool get_nextpt_2ndeuler_ten_quad(double first[2], double second[2], int &face_id, int type)
{
	////Using first order Euler method to get the next point
	double t[4];
	icMatrix2x2 ten;
	compute_tensor_at_quad(face_id, first[0], first[1], ten);
	/*old before 10/08/2007*/
	//get_tensor(first[0], first[1], t);
	//ten.entry[0][0]=t[0];
	//ten.entry[0][1]=t[1];
	//ten.entry[1][0]=t[2];
	//ten.entry[1][1]=t[3];

	/*second, we compute its major eigen vector*/
	icVector2 ev[2];
	cal_eigen_vector_sym(ten, ev);

	if(type == 1) /*use minor eigen vector field*/
		ev[0] = ev[1];

	double re = dot(tenline_dir_global, ev[0]);
	if(re<0) 
		ev[0] = -ev[0];

	icVector2 VecAtPoint=ev[0];

	if(length(VecAtPoint) < 1e-10) return false;

	double temp[2] = {0.};

	icVector2 before_VecAtPt = VecAtPoint;


	temp[0] = first[0] + euler_stepsize*VecAtPoint.entry[0];
	temp[1] = first[1] + euler_stepsize*VecAtPoint.entry[1];


	compute_tensor_at_quad(face_id, temp[0], temp[1], ten);
	/*old before 10/08/2007*/
	//get_tensor(temp[0], temp[1], t);
	//ten.entry[0][0]=t[0];
	//ten.entry[0][1]=t[1];
	//ten.entry[1][0]=t[2];
	//ten.entry[1][1]=t[3];

	/*second, we compute its major eigen vector*/
	cal_eigen_vector_sym(ten, ev);

	if(type == 1) /*use minor eigen vector field*/
		ev[0] = ev[1];

	re = dot(tenline_dir_global, ev[0]);
	if(re<0) 
		ev[0] = -ev[0];

	icVector2 VecAtPoint2=ev[0];
	
	/*
	we can adjust the step size here based on the difference between VecAtPoint and VecAtPoint2
	*/

	icVector2 total_v;

	/*The following is just mid-point average, it is not always true in every step!*/
	total_v.entry[0] = (VecAtPoint.entry[0] +  VecAtPoint2.entry[0])/2;
	total_v.entry[1] = (VecAtPoint.entry[1] +  VecAtPoint2.entry[1])/2;
	

	second[0] = first[0] + euler_stepsize*total_v.entry[0];
	second[1] = first[1] + euler_stepsize*total_v.entry[1];
}


/*  
    Considering the scalar field for the asymmetric tensor field tracing
*/
#include "shareinterfacevars.h"
extern SharedInterfaceVars sharedvars;
extern void cal_eigen_vector_asym(icMatrix2x2, double, icVector2 []);
extern void compute_phi_in_quad(int face, double x, double y, double &phi);


void get_tenvec_quad(double cur_p[2], double vec[2])
{
	double t[4];
	icMatrix2x2 ten;
	icVector2 ev[2];
	double re;

	compute_tensor_at_quad(g_face, cur_p[0], cur_p[1], ten);

	/*old before 10/08/2007*/
	//get_tensor(cur_p[0], cur_p[1], t);
	//ten.entry[0][0]=t[0];
	//ten.entry[0][1]=t[1];
	//ten.entry[1][0]=t[2];
	//ten.entry[1][1]=t[3];
	/*compute its major eigen vector*/

	if(!sharedvars.ApplyAsymFldOn/*sharedvars.ShowScalarFieldOn*/)
		cal_eigen_vector_sym(ten, ev);
	else
	{
		/*  need to obtain the interpolated phi  */
		double phi=0;
		compute_phi_in_quad(g_face, cur_p[0], cur_p[1], phi);
		cal_eigen_vector_asym(ten, phi, ev);
	}

	if(g_type == 1) /*use minor eigen vector field*/
		ev[0] = ev[1];

	re = dot(tenline_dir_global, ev[0]);
	if(re<0) 
		ev[0] = -ev[0];
	//normalize(ev[0]);

	vec[0] = ev[0].entry[0];
	vec[1] = ev[0].entry[1];
}


/*we try to implement the RK45 as numerical recipe*/
void RK45_2d(double pre_p[2], double next_p[2], double &hstep, double &hnext,
			 double eps, double &eps_did, void deriv(double cur_p[2], double vec[2]))
{
	double dp0[2], dp1[2], dp2[2], dp3[2];
	double temp[2] = {0.};
	double t_vec[2];
	
	/*compute dp0*/
	deriv(pre_p, t_vec);
	dp0[0] = hstep*t_vec[0];
	dp0[1] = hstep*t_vec[1];

	/*compute dp1*/
	temp[0]=pre_p[0]+dp0[0]/2;
	temp[1]=pre_p[1]+dp0[1]/2;
	deriv(temp, t_vec);
	dp1[0] = hstep*t_vec[0];
	dp1[1] = hstep*t_vec[1];
	
	/*compute dp2*/
	temp[0]=pre_p[0]+dp1[0]/2;
	temp[1]=pre_p[1]+dp1[1]/2;
	deriv(temp, t_vec);
	dp2[0] = hstep*t_vec[0];
	dp2[1] = hstep*t_vec[1];

	/*compute dp3*/
	temp[0]=pre_p[0]+dp2[0];
	temp[1]=pre_p[1]+dp2[1];
	deriv(temp, t_vec);
	dp3[0] = hstep*t_vec[0];
	dp3[1] = hstep*t_vec[1];

	/*compute the next position using dp0, dp1, dp2 and dp3*/
	next_p[0]=pre_p[0]+dp0[0]/6+dp1[0]/3+dp2[0]/3+dp3[0]/6;
	next_p[1]=pre_p[1]+dp0[1]/6+dp1[1]/3+dp2[1]/3+dp3[1]/6;

	/*evaluate the error*/
	deriv(next_p, t_vec);
	
	icVector2 ep;
	ep.entry[0]=(dp3[0]-hstep*t_vec[0])/6;
	ep.entry[1]=(dp3[1]-hstep*t_vec[1])/6;

	double error = length(ep);

	/*adjust the step size accordingly*/
	hnext = hstep;
	if(error<eps/10.) hnext = 2*hstep;
	if(error>eps*5) hnext = hstep/2;

	eps_did = error;
}

/*
Use RK45 to do tensor line computation
*/
bool get_nextpt_RK45_ten_quad(double first[2], double second[2], int &face_id, int type)
{
	g_face = face_id;
	g_type = type;

	double t_vec[2] = {0.};
	get_tenvec_quad(first, t_vec);
	icVector2 vec;
	vec.entry[0] = t_vec[0];
	vec.entry[1] = t_vec[1];
	if(length(vec) < 1.e-10) return false;

	double eps_did, eps = 1.e-9;
	int i;
	double hstep, hnext;

	hstep = predict_stepsize;
	if(hstep > quadmesh->xinterval)
	{
		hstep = quadmesh->xinterval;
		predict_stepsize = hstep;
	}


	for(i=0; i<10; i++)
	{
		hstep = predict_stepsize;
		RK45_2d(first, second, hstep, hnext, eps, eps_did, get_tenvec_quad);
		predict_stepsize = hnext;
		if(eps_did<eps)
			break;
	}
	return true;
}


/*we try to implement the RK45 as numerical recipe*/
void RK23_2d(double pre_p[2], double next_p[2], double &hstep_loc, double &hnext,
			 double eps, double &eps_did, void deriv(double cur_p[2], double vec[2]))
{
	double dp0[2], dp1[2];
	double temp[2] = {0.};
	double t_vec[2];
	
	/*compute dp0*/
	deriv(pre_p, t_vec);
	dp0[0] = hstep_loc*t_vec[0];
	dp0[1] = hstep_loc*t_vec[1];

	/*compute dp1*/
	temp[0]=pre_p[0]+dp0[0];
	temp[1]=pre_p[1]+dp0[1];
	deriv(temp, t_vec);
	dp1[0] = hstep_loc*t_vec[0];
	dp1[1] = hstep_loc*t_vec[1];
	

	/*compute the next position using dp0, dp1, dp2 and dp3*/
	next_p[0]=pre_p[0]+dp0[0]/2+dp1[0]/2;
	next_p[1]=pre_p[1]+dp0[1]/2+dp1[1]/2;

	/*evaluate the error*/
	deriv(next_p, t_vec);
	
	icVector2 ep;
	ep.entry[0]=(dp1[0]-hstep_loc*t_vec[0])/3;
	ep.entry[1]=(dp1[1]-hstep_loc*t_vec[1])/3;

	double error = length(ep);

	/*adjust the step size accordingly*/
	hnext = hstep_loc;
	if(error<eps/10.) hnext = 2*hstep_loc;
	if(error>eps*5) hnext = hstep_loc/2;

	eps_did = error;
}

/*
Use RK23 to do tensor line computation
*/
bool get_nextpt_RK23_ten_quad(double first[2], double second[2], int &face_id, int type)
{
	g_face = face_id;
	g_type = type;

	double t_vec[2] = {0.};
	get_tenvec_quad(first, t_vec);
	icVector2 vec;
	vec.entry[0] = t_vec[0];
	vec.entry[1] = t_vec[1];
	if(length(vec) < 1.e-10) return false;

	double eps_did, eps = 1.e-9;
	int i;
	double hstep_loc, hnext;

	hstep_loc = predict_stepsize;
	if(hstep_loc > quadmesh->xinterval/2.)
	{
		hstep_loc = quadmesh->xinterval/2.;
		predict_stepsize = hstep_loc;
	}


	for(i=0; i<5; i++)
	{
		hstep_loc = predict_stepsize;
		RK23_2d(first, second, hstep_loc, hnext, eps, eps_did, get_tenvec_quad);
		predict_stepsize = hnext;
		if(eps_did<eps)
			break;
	}
	return true;
}


void global_test_trace(double x, double y)
{
	double t[4];
	icMatrix2x2 ten;
	icVector2 ev[2], startdir;
	double re;
	double p[2]={x, y};
	double prep[2]={x, y};
	int nlines = 0;
	num_linesegs_curtraj[cur_traj_index] = 0;

	get_tensor(p[0], p[1], t);
	ten.entry[0][0]=t[0];
	ten.entry[0][1]=t[1];
	ten.entry[1][0]=t[2];
	ten.entry[1][1]=t[3];

	/*second, we compute its major eigen vector*/
	cal_eigen_vector_sym(ten, ev);
	tenline_dir_global=startdir=ev[0];

	euler_stepsize = quadmesh->xinterval/5;

	/*start from the positive direction*/
	while(p[0]>-0.1&&p[0]<=1.1&&p[1]>=-0.1&&p[1]<=1.1 && nlines<1600)
	{
		/**/
		prep[0]=p[0];
		prep[1]=p[1];

		/*first order Euler*/
		//get_tensor(p[0], p[1], t);
		//ten.entry[0][0]=t[0];
		//ten.entry[0][1]=t[1];
		//ten.entry[1][0]=t[2];
		//ten.entry[1][1]=t[3];

		///*second, we compute its major eigen vector*/
		//cal_eigen_vector_sym(ten, ev);
		//re = dot(tenline_dir_global, ev[0]);
		//if(re<0) 
		//	ev[0] = -ev[0];
		//
		////if(get_nextpt_2ndeuler_ten_quad(p, p, 0))
		//if(length(ev[0])<1e-18) break;
		//normalize(ev[0]);
		//ev[0] = quadmesh->xinterval/10*ev[0];
		//p[0] += ev[0].entry[0];
		//p[1] += ev[0].entry[1];
		get_nextpt_2ndeuler_ten_quad_gl(prep, p);

		/*save to the trajectory data structure*/
		trajectories[cur_traj_index][nlines].gend[0]=p[0];
		trajectories[cur_traj_index][nlines].gend[1]=p[1];
		trajectories[cur_traj_index][nlines].gstart[0]=prep[0];
		trajectories[cur_traj_index][nlines].gstart[1]=prep[1];
		nlines++;

		tenline_dir_global.entry[0] = p[0]-prep[0];
		tenline_dir_global.entry[1] = p[1]-prep[1];
	}

	/*negative direction*/
	tenline_dir_global=-startdir;
	p[0]=prep[0]=x;
	p[1]=prep[1]=y;

	/*start from the positive direction*/
	while(p[0]>-0.1&&p[0]<=1.1&&p[1]>=-0.1&&p[1]<=1.1 && nlines<1600)
	{
		/**/
		prep[0]=p[0];
		prep[1]=p[1];

		/*first order Euler*/
		//get_tensor(p[0], p[1], t);
		//ten.entry[0][0]=t[0];
		//ten.entry[0][1]=t[1];
		//ten.entry[1][0]=t[2];
		//ten.entry[1][1]=t[3];

		///*second, we compute its major eigen vector*/
		//cal_eigen_vector_sym(ten, ev);
		//re = dot(tenline_dir_global, ev[0]);
		//if(re<0) 
		//	ev[0] = -ev[0];
		//
		////if(get_nextpt_2ndeuler_ten_quad(p, p, 0))
		//if(length(ev[0])<1e-18) break;
		//normalize(ev[0]);
		//ev[0] = quadmesh->xinterval/10*ev[0];
		//p[0] += ev[0].entry[0];
		//p[1] += ev[0].entry[1];

		get_nextpt_2ndeuler_ten_quad_gl(prep, p);

		/*save to the trajectory data structure*/
		trajectories[cur_traj_index][nlines].gend[0]=p[0];
		trajectories[cur_traj_index][nlines].gend[1]=p[1];
		trajectories[cur_traj_index][nlines].gstart[0]=prep[0];
		trajectories[cur_traj_index][nlines].gstart[1]=prep[1];
		nlines++;

		tenline_dir_global.entry[0] = p[0]-prep[0];
		tenline_dir_global.entry[1] = p[1]-prep[1];
	}

	num_linesegs_curtraj[cur_traj_index] = nlines;
}


/*
2nd integrator using global coordinates only
*/
bool get_nextpt_2ndeuler_ten_quad_gl(double first[2], double second[2])
{
	////Using first order Euler method to get the next point
	icMatrix2x2 ten;
	double t[4] = {0.};
	get_tensor(first[0], first[1], t);
	ten.entry[0][0]=t[0];
	ten.entry[0][1]=t[1];
	ten.entry[1][0]=t[2];
	ten.entry[1][1]=t[3];

	/*second, we compute its major eigen vector*/
	icVector2 ev[2];
	cal_eigen_vector_sym(ten, ev);
	double re = dot(tenline_dir_global, ev[0]);
	if(re<0) 
		ev[0] = -ev[0];

	icVector2 VecAtPoint=ev[0];

	if(length(VecAtPoint) < 1e-20) return false;
	normalize(VecAtPoint);

	double temp[2] = {0.};

	temp[0] = first[0] + euler_stepsize*VecAtPoint.entry[0];
	temp[1] = first[1] + euler_stepsize*VecAtPoint.entry[1];


	get_tensor(temp[0], temp[1], t);
	ten.entry[0][0]=t[0];
	ten.entry[0][1]=t[1];
	ten.entry[1][0]=t[2];
	ten.entry[1][1]=t[3];

	/*second, we compute its major eigen vector*/
	cal_eigen_vector_sym(ten, ev);
	re = dot(tenline_dir_global, ev[0]);
	if(re<0) 
		ev[0] = -ev[0];

	icVector2 VecAtPoint2=ev[0];
	normalize(VecAtPoint2);
	

	/*
	we can adjust the step size here based on the difference between VecAtPoint and VecAtPoint2
	*/
	
	icVector2 total_v;

	/*The following is just mid-point average, it is not always true in every step!*/
	total_v.entry[0] = (VecAtPoint.entry[0] +  VecAtPoint2.entry[0])/2;
	total_v.entry[1] = (VecAtPoint.entry[1] +  VecAtPoint2.entry[1])/2;
	

	second[0] = first[0] + euler_stepsize*total_v.entry[0];
	second[1] = first[1] + euler_stepsize*total_v.entry[1];
}
