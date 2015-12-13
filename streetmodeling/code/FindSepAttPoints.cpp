/*
FindSepAttPoints.cpp
*/


#include "stdafx.h"

#include "FindSepAttPoints.h"

#include "FindSCC.h"
#include "VFDataStructure.h"
#include "VFAnalysis.h"


extern SCCList scclist;
extern Polygon3D Object;

extern bool UseNormalizedVF;

double  R3Filter = 0;

int SymmetricOn = 0;

void CalEigenForEdges(Edge *cur_e);


/* The routine of calculate eigen vectors for a given triangle */

void CalEigenForTriangle(int triangle)
{
	Face *face;
	int *verts;

	double x[3], y[3]; //
	double vx[3], vy[3];//

	double la, lb, lc, ld, le, lf;

	int i, j;

	///These information has been obtained before
	///Can we reuse them and enhance the performance 
    face = Object.flist[triangle];
	verts = face->verts;

	//For each triangle, calculate the vector for all vertices
	for (j=0; j<face->nverts; j++) {
        x[j] = Object.vlist[verts[j]]->x;
		y[j] = Object.vlist[verts[j]]->y;

		if(UseNormalizedVF)
		{
   			vx[j] = Object.vlist[verts[j]]->vec.entry[0];
			vy[j] = Object.vlist[verts[j]]->vec.entry[1];
		}
		else
		{
			vx[j] = Object.vlist[verts[j]]->vec_J.entry[0]; //07/19/06
			vy[j] = Object.vlist[verts[j]]->vec_J.entry[1];
		}
	}


	double ta, tb, tc, td, in_ta, in_tb, in_tc, in_td, det;

	//solve la, lb, lc
	ta = x[1]-x[0];
	tb = y[1]-y[0];
	tc = x[2]-x[0];
	td = y[2]-y[0];
	det = ta*td-tb*tc;

	in_ta = td/det;
	in_tb = -tb/det;
	in_tc = -tc/det;
	in_td = ta/det;

	la = in_ta*(vx[1]-vx[0])+in_tb*(vx[2]-vx[0]);
	lb = in_tc*(vx[1]-vx[0])+in_td*(vx[2]-vx[0]);

	//solve ld, le, lf

	ld = in_ta*(vy[1]-vy[0])+in_tb*(vy[2]-vy[0]);
	le = in_tc*(vy[1]-vy[0])+in_td*(vy[2]-vy[0]);


	/**------Calculate the eigen values and vectors for the given matrix ------**/
	double temp1[][2]= {{la, lb}, {ld, le}};
	double trace = (la + le)/2.;
	//double temp2[][2]= {{la, (lb+ld)/2.}, {(ld+lb)/2., le}};  //use symmetric matrix
	double temp2[][2]= {{la-trace, (lb+ld)/2.}, {(ld+lb)/2., le-trace}};  //use symmetric matrix
    icMatrix2x2 mat;

	if(SymmetricOn == 0)
	{
		mat.set(temp1);
	}
	else{
		mat.set(temp2);
	}

	face->Jacobian.set(mat); //save the Jacobian of the triangle  07/20/06

	//Calculate the eigen values
	double A, B, C, delta;

	double eigen[2] = {0.};
	icVector2 evec[2];


	A = 1;

	if(SymmetricOn == 0)
	{
		B = -(la + le);
		C = (la * le - lb * ld);                //use non-symmetric
	}
	else
	{
		//B = -(la + le);
		//C = (la * le - (lb+ld)* (lb+ld)/4.);    //use symmetric matrix
		B = -(la + le - 2*trace);
		C = ((la-trace) * (le-trace) - (lb+ld)* (lb+ld)/4.);    //use symmetric matrix
	}

	delta =  B*B - 4*A*C;

	if( delta >= 0)
	{
		eigen[0] = (-B + sqrt(delta))/2;
		eigen[1] = (-B - sqrt(delta))/2;

		////find the largest (absolute) eigen values, and store it to the first element
		//if(abs(eigen[0]) < abs(eigen[1]))
		//{
		//	double temp = eigen[0];
		//	eigen[0] = eigen[1];
		//	eigen[1] = temp;
		//}

		//for real eigen values, we calculate the eigen vectors of it
		GetEigenVectors(mat, eigen, evec);

		face->evalues[0] = eigen[0];
		face->evalues[1] = eigen[1];
		face->eigen[0] = evec[0];
		face->eigen[1] = evec[1];

		normalize(face->eigen[0]);
		normalize(face->eigen[1]);

		//////Modified at 07/20/06
		//face->eigen[1].entry[0] = -face->eigen[0].entry[1];
		//face->eigen[1].entry[1] = face->eigen[0].entry[0];
	}

	else
	{
		face->eigen[0].set(0, 0);
		face->eigen[1].set(0, 0);
	}

	//double temp3[][2]= {{la-trace, (lb+ld)/2.}, {(ld+lb)/2., le-trace}};  //use symmetric matrix

	//double theta = atan(((lb+ld)/2.)/(la-trace))/2.;

	//face->eigen[0].entry[0] = cos(theta);
	//face->eigen[0].entry[1] = sin(theta);

	//face->eigen[1].entry[0] = cos(theta + M_PI/2.);
	//face->eigen[1].entry[1] = sin(theta + M_PI/2.);

	//normalize(face->eigen[0]);
	//normalize(face->eigen[1]);

}


void CalJacobianForTriangle(int triangle)
{
	Face *face;
	int *verts;

	double x[3], y[3]; //
	double vx[3], vy[3];//

	double la, lb, lc, ld, le, lf;

	int i, j;

	///These information has been obtained before
	///Can we reuse them and enhance the performance 
    face = Object.flist[triangle];
	verts = face->verts;

	//For each triangle, calculate the vector for all vertices
	for (j=0; j<face->nverts; j++) {
        x[j] = Object.vlist[verts[j]]->x;
		y[j] = Object.vlist[verts[j]]->y;

		if(UseNormalizedVF)
		{
			vx[j] = Object.vlist[verts[j]]->vec.entry[0];
			vy[j] = Object.vlist[verts[j]]->vec.entry[1];
		}

		else{
			/* we use the vector before normalization to calculate the Jacobian of the triangle
			02/21/07*/
			vx[j] = Object.vlist[verts[j]]->vec_J.entry[0]; //07/19/06
			vy[j] = Object.vlist[verts[j]]->vec_J.entry[1];
		}
	}


	//new method to get the coefficients 07/19/06
	double ta, tb, tc, td, in_ta, in_tb, in_tc, in_td, det;

	//solve la, lb, lc
	ta = x[1]-x[0];
	tb = y[1]-y[0];
	tc = x[2]-x[0];
	td = y[2]-y[0];
	det = ta*td-tb*tc;

	in_ta = td/det;
	in_tb = -tb/det;
	in_tc = -tc/det;
	in_td = ta/det;

	la = in_ta*(vx[1]-vx[0])+in_tb*(vx[2]-vx[0]);
	lb = in_tc*(vx[1]-vx[0])+in_td*(vx[2]-vx[0]);

	////solve ld, le, lf

	ld = in_ta*(vy[1]-vy[0])+in_tb*(vy[2]-vy[0]);
	le = in_tc*(vy[1]-vy[0])+in_td*(vy[2]-vy[0]);



	double temp1[][2]= {{la, lb}, {ld, le}};

	//set the general tensor for the triangle
	face->Jacobian.set(temp1);

	/* Calculate the eigen values of the Jacobian */

	double A, B, C, evalues[2] = {0.};

	A = 1;
	B = -(la + le);
	C = (la * le - lb * ld);                

	double delta =  B*B - 4*A*C;

	if(delta >= 0)
	{
		evalues[0] = (-B + sqrt(delta))/2;
		evalues[1] = (-B - sqrt(delta))/2;

		face->evalues[0] = evalues[0];
		face->evalues[1] = evalues[1];

		/* can we some how scale them down to (-1, 1)? */
	}

	else
	{
		/* get the real part of the eigen values */
		face->evalues[0] = -B/2;
		face->evalues[1] = -B/2;

		/* second thought: get the magnitude of the complex eigen values */

		/* third thought: get the imaginary part of the complex eigen values */
		//face->evalues[0] = sqrt(abs(delta))/2;
		//face->evalues[1] = sqrt(abs(delta))/2;
	}
}


void CalEigenForVertex(int vertID)
{
	int i;
	Face *face;
	icMatrix2x2 mat;
	Vertex *v = Object.vlist[vertID];

	mat.entry[0][0] = mat.entry[0][1] = mat.entry[1][0] = mat.entry[1][1] = 0.;

	for(i = 0; i < v->Num_corners; i++)
	{
		if(v->Corners[i] < 0)
			continue;

		face = Object.flist[Object.clist[v->Corners[i]]->t];
		mat = mat + face->Jacobian;
	}

	double la, lb, lc, ld, A, B, C, delta;

	la = mat.entry[0][0]/v->Num_corners;
	lb = mat.entry[0][1]/v->Num_corners;
	lc = mat.entry[1][0]/v->Num_corners;
	ld = mat.entry[1][1]/v->Num_corners;
	//la = mat.entry[0][0];
	//lb = mat.entry[0][1];
	//lc = mat.entry[1][0];
	//ld = mat.entry[1][1];

	if(SymmetricOn == 1)
	{
		//Calculate the symmetric matrix
		double ta, tb, tc, td;
		double trace = la+ld;

		ta = la-trace;
		td = ld-trace;
		tc = tb = (lb+lc)/2.;

		mat.entry[0][0] = ta;
		mat.entry[0][1] = tb;
		mat.entry[1][0] = tc;
		mat.entry[1][1] = td;

		la = ta;
		lb = tb;
		lc = tc;
		ld = td;
	}

	////save the Jacobian  07/20/06
	v->Jacobian.set(mat);  //this is still a general tensor

	double evalues[2] = {0.};

	A = 1;
	B = -(la + ld);
	C = (la * ld - lb * lc);                

	delta =  B*B - 4*A*C;

	if( delta >= 0)
	{
		evalues[0] = (-B + sqrt(delta))/2;
		evalues[1] = (-B - sqrt(delta))/2;

		////find the largest (absolute) eigen values, and store it to the first element
		//if(abs(evalues[0]) < abs(evalues[1]))
		//{
		//	double temp = evalues[0];
		//	evalues[0] = evalues[1];
		//	evalues[1] = temp;
		//}

		//for real eigen values, we calculate the eigen vectors of it
		GetEigenVectors(mat, evalues, v->evec);

		normalize(v->evec[0]);
		normalize(v->evec[1]);
		//v->evec[1].entry[0] = -v->evec[0].entry[1];
		//v->evec[1].entry[1] = v->evec[0].entry[0];

	}

	else
	{
		v->evec[0].set(0, 0);
		v->evec[1].set(0, 0);
	}
}



void CalEigenForEdge(Edge *cur_e)
{
	Vertex *v1, *v2;
	icMatrix2x2 mat;

	v1 = Object.vlist[cur_e->verts[0]];
	v2 = Object.vlist[cur_e->verts[1]];

	mat = v1->Jacobian + v2->Jacobian;
	double la, lb, lc, ld, A, B, C, delta;

	la = mat.entry[0][0];
	lb = mat.entry[0][1];
	lc = mat.entry[1][0];
	ld = mat.entry[1][1];

	////save the Jacobian of the edge
	cur_e->Jacobian.set(mat);

	if(SymmetricOn == 1)
	{
		//Calculate the symmetric matrix
		double ta, tb, tc, td;
		double trace = la+ld;

		ta = la-trace;
		td = ld-trace;
		tc = tb = (lb+lc)/2.;

		mat.entry[0][0] = ta;
		mat.entry[0][1] = tb;
		mat.entry[1][0] = tc;
		mat.entry[1][1] = td;

		la = ta;
		lb = tb;
		lc = tc;
		ld = td;
	}

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
		GetEigenVectors(mat, evalues, cur_e->evec);

		normalize(cur_e->evec[0]);
		normalize(cur_e->evec[1]);
	}

	else
	{
		cur_e->evec[0].set(0, 0);
		cur_e->evec[1].set(0, 0);
	}
}


void CalEigenForAllEdges()
{
	int i, j;
	Face *face;
	Edge *cur_e;

	//initialization
	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		for(j = 0; j < 3; j++)
		{
			cur_e = face->edges[j];
			cur_e->visited = 0;
			cur_e->evec[0].set(0, 0);
			cur_e->evec[1].set(0, 0);
		}
	}

	//calculate the eigen vectors
	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		for(j = 0; j < 3; j++)
		{
			cur_e = face->edges[j];

			if(cur_e->visited == 1)
				continue;

			CalEigenForEdge(cur_e);
		}
	}
}



void GetEigenVectors(icMatrix2x2 mat, double evalues[2], icVector2 ev[2])
{
	//first calculate the dominant of the matrix, if it is zero, return
	if(abs(CalDeterminant(mat)) < 1e-30)
		return;

	ev[0].entry[0] = mat.entry[1][1] - evalues[0];
	ev[0].entry[1] = - mat.entry[1][0];
	
	ev[1].entry[0] = mat.entry[1][1] - evalues[1];
	ev[1].entry[1] = - mat.entry[1][0];

	//ev[0].entry[0] = mat.entry[0][1];
	//ev[0].entry[1] = evalues[0] - mat.entry[0][0];
	//
	//ev[1].entry[0] = mat.entry[0][1];
	//ev[1].entry[1] = evalues[1] - mat.entry[0][0];

	//ev[1].entry[0] = evalues[1] - mat.entry[1][1];
	//ev[1].entry[1] = mat.entry[1][0];
}

double CalDeterminant(icMatrix2x2 mat)
{
	return (mat.entry[0][0]*mat.entry[1][1]-mat.entry[0][1]*mat.entry[1][0]);
}


void CalEigenVecsForSCC()
{
	int i, j;

	for(i = 0; i < scclist.num_sccs; i++)
	{
		if(scclist.scccomponents[i].num_nodes < 2)
			continue;

		for(j = 0; j < scclist.scccomponents[i].num_nodes; j++)
		{
			CalEigenForTriangle(scclist.scccomponents[i].nodes[j]);
			//CalJacobianForTriangle(scclist.scccomponents[i].nodes[j]);
		}
	}

	/*   Calculate the eigen vectors on vertices instead of faces 07/20/06 */
	for(i = 0; i < Object.nverts; i++)
	{
		Object.vlist[i]->visited = 0;
	}
	
	Face *face;

	for(i = 0; i < scclist.num_sccs; i++)
	{
		if(scclist.scccomponents[i].num_nodes < 2)
			continue;

		for(j = 0; j < scclist.scccomponents[i].num_nodes; j++)
		{
			face = Object.flist[scclist.scccomponents[i].nodes[j]];

			for(int k = 0; k < face->nverts; k++)
			{
				if(Object.vlist[face->verts[k]]->visited == 1)
					continue;

				CalEigenForVertex(face->verts[k]);
				Object.vlist[face->verts[k]]->visited = 1;
			}
		}
	}
}


/* To calculate one separation point 07/17/06 */
bool ExistSeparationPoint(Edge *cur_e, int triangle, double sep[2])
{
	Vertex *v1, *v2;
	v1 = Object.vlist[cur_e->verts[0]];
	v2 = Object.vlist[cur_e->verts[1]];

	////Another filter to remove some edges 07/21/06
	double cr = cur_e->evec[0].entry[0]*cur_e->evec[1].entry[1]\
		-cur_e->evec[1].entry[0]*cur_e->evec[0].entry[1];

	if(abs(cr) < 1e-8)
		return false;

	////if the R3 component of the Jacobian of the edge is smaller than some threshold, ignore that
	double r1, r2;
	
	icVector3 r_coef;
	r_coef.entry[0] = (cur_e->Jacobian.entry[0][0]+cur_e->Jacobian.entry[1][1])/2.;
	r_coef.entry[1] = (cur_e->Jacobian.entry[0][1]-cur_e->Jacobian.entry[1][0])/2.;

	r1 = (cur_e->Jacobian.entry[0][0]-cur_e->Jacobian.entry[1][1])/2.;
	r2 = (cur_e->Jacobian.entry[0][1]+cur_e->Jacobian.entry[1][0])/2.;

	r_coef.entry[2] = sqrt(r1*r1+r2*r2); 

	////Remove those edges that have small R1, R2 and R3 magnitude
	if(length(r_coef) < 0.01)
		return false;

	normalize(r_coef);

	////Remove those edges that have weak R3 component
	if(r_coef.entry[2] < R3Filter)
		return false;

	double re1 = cur_e->evec[0].entry[0]*v1->vec.entry[1] - cur_e->evec[0].entry[1]*v1->vec.entry[0];
	double re2 = cur_e->evec[0].entry[0]*v2->vec.entry[1] - cur_e->evec[0].entry[1]*v2->vec.entry[0];

	if(re1*re2 >= 0)
		return false;

	if(abs(re1) <= 1e-40)
	{
		sep[0] = v1->x;
		sep[1] = v1->y;
		//return true;
		return false;
	}

	if(abs(re2) <= 1e-40)
	{
		sep[0] = v2->x;
		sep[1] = v2->y;
		//return true;
		return false;
	}

	/* The following is not an accurate estimation */
	double ratio = abs(re1)/(abs(re2)+abs(re1));
	//double ratio;
	//ratio = (vec_sep.entry[0]*v1->vec.entry[1]-v1->vec.entry[0]*vec_sep.entry[1])/
	//	((v2->vec.entry[0]-v1->vec.entry[0])*vec_sep.entry[1]-(v2->vec.entry[1]-v1->vec.entry[1])*vec_sep.entry[0]);

	sep[0] = (1-ratio)*v1->x + ratio*v2->x;
	sep[1] = (1-ratio)*v1->y + ratio*v2->y;
	return true;
}


/* To calculate one attachment point */
bool ExistAttachmentPoint(Edge *cur_e, int triangle, double attp[2])
{
	Vertex *v1, *v2;
	v1 = Object.vlist[cur_e->verts[0]];
	v2 = Object.vlist[cur_e->verts[1]];

	////Another filter to remove some edges 07/21/06
	double cr = cur_e->evec[0].entry[0]*cur_e->evec[1].entry[1]\
		-cur_e->evec[1].entry[0]*cur_e->evec[0].entry[1];

	if(abs(cr) < 1e-8)
		return false;

	////if the R3 component of the Jacobian of the edge is smaller than some threshold, ignore that
	icVector3 r_coef;
	r_coef.entry[0] = (cur_e->Jacobian.entry[0][0]+cur_e->Jacobian.entry[1][1])/2.;
	r_coef.entry[1] = (cur_e->Jacobian.entry[0][1]-cur_e->Jacobian.entry[1][0])/2.;

	double r1, r2;
	r1 = (cur_e->Jacobian.entry[0][0]-cur_e->Jacobian.entry[1][1])/2.;
	r2 = (cur_e->Jacobian.entry[0][1]+cur_e->Jacobian.entry[1][0])/2.;
	
	r_coef.entry[2] = sqrt(r1*r1+r2*r2); 

	////Remove those edges that have small R1, R2 and R3 magnitude
	if(length(r_coef) < 0.01)
		return false;

	normalize(r_coef);

	////Remove those edges that have weak R3 component
	if(r_coef.entry[2] < R3Filter)
		return false;
	
	double re1 = cur_e->evec[1].entry[0]*v1->vec.entry[1] - cur_e->evec[1].entry[1]*v1->vec.entry[0];
	double re2 = cur_e->evec[1].entry[0]*v2->vec.entry[1] - cur_e->evec[1].entry[1]*v2->vec.entry[0];

	if(re1 * re2 >= 0)
		return false;

	if(abs(re1) <= 1e-40)
	{
		attp[0] = v1->x;
		attp[1] = v1->y;
		//return true;
		return false;
	}

	if(abs(re2) <= 1e-40)
	{
		attp[0] = v2->x;
		attp[1] = v2->y;
		//return true;
		return false;
	}

	/* The following is not an accurate estimation */
	double ratio = abs(re1)/(abs(re2)+abs(re1));
	//double ratio;
	//ratio = (vec_att.entry[0]*v1->vec.entry[1]-v1->vec.entry[0]*vec_att.entry[1])/
	//	((v2->vec.entry[0]-v1->vec.entry[0])*vec_att.entry[1]-(v2->vec.entry[1]-v1->vec.entry[1])*vec_att.entry[0]);

	attp[0] = (1-ratio)*v1->x + ratio*v2->x;
	attp[1] = (1-ratio)*v1->y + ratio*v2->y;
	return true;
}


void CalAllSpecialPoints()
{
	int i, j, k;
	Face *face;
	Edge *cur_e;
	double p[2] = {0.};

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
		if(scclist.scccomponents[i].num_nodes < 2)
			continue;

		for(j = 0; j < scclist.scccomponents[i].num_nodes; j++)
		{
			////First, you need to mark all the boundary edges, and don't calculate the special points on them!
			////Now, just let me ignore that step and experiment a little bit further

			face = Object.flist[scclist.scccomponents[i].nodes[j]];

			for(k = 0; k < 3; k++)
			{
				cur_e = face->edges[k];

				if(cur_e->visited != 0)
					continue;

				////calculate the separation point
				if(ExistSeparationPoint(cur_e, face->index,p))
				{
					cur_e->sep.entry[0] = p[0];
					cur_e->sep.entry[1] = p[1];
					cur_e->find_sep = 1;
				}
				
				////calculate the attachment point
				if(ExistAttachmentPoint(cur_e, face->index,p))
				{
					cur_e->attp.entry[0] = p[0];
					cur_e->attp.entry[1] = p[1];
					cur_e->find_attp = 1;
				}
			}
		}
	}
}


void CalEigenVecsForMesh()
{
	int i, j;

	for(i = 0; i < Object.nfaces; i++)
	{
		//CalEigenForTriangle(i);
		CalJacobianForTriangle(i);
	}

	/*   Calculate the eigen vectors on vertices instead of faces 07/20/06 */
	for(i = 0; i < Object.nverts; i++)
	{
		CalEigenForVertex(i);
	}

	CalEigenForAllEdges();
}



void CalAllSpecialPoints_alltraingles()
{
	int i, j, k;
	Face *face;
	Edge *cur_e;
	double p[2] = {0.};

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

	for(j = 0; j < Object.nfaces; j++)
	{
		////First, you need to mark all the boundary edges, and don't calculate the special points on them!
		////Now, just let me ignore that step and experiment a little bit further

		face = Object.flist[j];

		for(k = 0; k < 3; k++)
		{
			cur_e = face->edges[k];

			if(cur_e->visited != 0)
				continue;

			////calculate the separation point
			if(ExistSeparationPoint(cur_e, face->index,p))
			{
				cur_e->sep.entry[0] = p[0];
				cur_e->sep.entry[1] = p[1];
				cur_e->find_sep = 1;
			}
			
			////calculate the attachment point
			if(ExistAttachmentPoint(cur_e, face->index,p))
			{
				cur_e->attp.entry[0] = p[0];
				cur_e->attp.entry[1] = p[1];
				cur_e->find_attp = 1;
			}

			cur_e->visited = 1;
		}
	}
}

void InitFindSepAttPoints()
{
}


/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*************************************************************************/

/*  the following routines are used for the limit cycle detection 07/22/06  */

void CalEigenForMatrix2x2(icMatrix2x2 mat, icVector2 evec[2])
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
		GetEigenVectors(mat, evalues, evec);

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
Calculate the eigen vectors of the edges inside the specific SCC using its symmetric Jacobian
*/
void CalEigenVecsForEdgesInSCC(int scc_index)
{
	int i, j;
	Face *face;
	Edge *cur_e;

	icMatrix2x2 mat;

	double ta, tb, tc, td, trace;

	for(i = 0; i < scclist.scccomponents[scc_index].num_nodes; i++)
	{
		face = Object.flist[scclist.scccomponents[scc_index].nodes[i]];

		for(j = 0; j < 3; j++)
		{
			cur_e = face->edges[j];

			if(cur_e->visited == 1)
				continue;

			////the orginal matrix saved in each edge is not a symmetric one,
			////so you need to get the symmetric one first
			trace = cur_e->Jacobian.entry[0][0]+cur_e->Jacobian.entry[1][1];

			ta = cur_e->Jacobian.entry[0][0]-trace;
			td = cur_e->Jacobian.entry[1][1]-trace;
			tc = tb = (cur_e->Jacobian.entry[0][1]+cur_e->Jacobian.entry[1][0])/2.;

			mat.entry[0][0] = ta;
			mat.entry[0][1] = tb;
			mat.entry[1][0] = tc;
			mat.entry[1][1] = td;
			 
			CalEigenForMatrix2x2(mat, cur_e->evec);
			
			cur_e->visited = 1;
		}
	}
}

/* To calculate one separation point 07/17/06 */
bool ExistSeparationPoint_SCC(Edge *cur_e, int triangle, double sep[2])
{
	Vertex *v1, *v2;
	v1 = Object.vlist[cur_e->verts[0]];
	v2 = Object.vlist[cur_e->verts[1]];

	icMatrix2x2 mat = v1->Jacobian + v2->Jacobian;

	cur_e->Jacobian.set(mat);

	////if the R3 component of the Jacobian of the edge is smaller than some threshold, ignore that
	double r1, r2;
	
	icVector3 r_coef;
	r_coef.entry[0] = (mat.entry[0][0]+mat.entry[1][1])/2.;
	r_coef.entry[1] = (mat.entry[0][1]-mat.entry[1][0])/2.;

	r1 = (mat.entry[0][0]-mat.entry[1][1])/2.;
	r2 = (mat.entry[0][1]+mat.entry[1][0])/2.;

	r_coef.entry[2] = sqrt(r1*r1+r2*r2); 

	////Remove those edges that have small R1, R2 and R3 magnitude
	//we calculate as many points as possible (07/09/07)
	//if(length(r_coef) < 0.04)
	//{
	//	return false;
	//}

	normalize(r_coef);

	////Remove those edges that have weak R3 component
	//if(r_coef.entry[2] < R3Filter)

	//if(r_coef.entry[2] < 0.2)
	//{
	//	return false;
	//}

	////Calculate the eigen vector of the symmetric matrix of the Jacobian of the edge
	double trace = (mat.entry[0][0] + mat.entry[1][1])/2.;
	double antidiag = (mat.entry[0][1] + mat.entry[1][0])/2.;
	mat.entry[0][0] = mat.entry[0][0] - trace;
	mat.entry[1][1] = mat.entry[1][1] - trace;
	mat.entry[0][1] = antidiag;
	mat.entry[1][0] = antidiag;
	CalEigenForMatrix2x2(mat, cur_e->evec);  //use asymmetric here 07/28/06

	////Another filter to remove some edges 07/21/06
	double cr = cur_e->evec[0].entry[0]*cur_e->evec[1].entry[1]\
		-cur_e->evec[1].entry[0]*cur_e->evec[0].entry[1];

	if(abs(cr) < 1e-10)
	{
		return false;
	}


	double re1 = cur_e->evec[0].entry[0]*v1->vec.entry[1] - cur_e->evec[0].entry[1]*v1->vec.entry[0];
	double re2 = cur_e->evec[0].entry[0]*v2->vec.entry[1] - cur_e->evec[0].entry[1]*v2->vec.entry[0];

	if(re1*re2 >= 0)
	{
		return false;
	}

	if(abs(re1) <= 1e-40)
	{
		sep[0] = v1->x;
		sep[1] = v1->y;
		return false;
	}

	if(abs(re2) <= 1e-40)
	{
		sep[0] = v2->x;
		sep[1] = v2->y;
		return false;
	}

	/* The following is not an accurate estimation */
	double ratio = abs(re1)/(abs(re2)+abs(re1));

	sep[0] = (1-ratio)*v1->x + ratio*v2->x;
	sep[1] = (1-ratio)*v1->y + ratio*v2->y;
	return true;
}


/* To calculate one attachment point */
bool ExistAttachmentPoint_SCC(Edge *cur_e, int triangle, double attp[2])
{
	Vertex *v1, *v2;
	v1 = Object.vlist[cur_e->verts[0]];
	v2 = Object.vlist[cur_e->verts[1]];


	//icMatrix2x2 mat = v1->Jacobian + v2->Jacobian;

	//cur_e->Jacobian.set(mat);


	////if the R3 component of the Jacobian of the edge is smaller than some threshold, ignore that
	icVector3 r_coef;
	r_coef.entry[0] = (cur_e->Jacobian.entry[0][0]+cur_e->Jacobian.entry[1][1])/2.;
	r_coef.entry[1] = (cur_e->Jacobian.entry[0][1]-cur_e->Jacobian.entry[1][0])/2.;

	double r1, r2;
	r1 = (cur_e->Jacobian.entry[0][0]-cur_e->Jacobian.entry[1][1])/2.;
	r2 = (cur_e->Jacobian.entry[0][1]+cur_e->Jacobian.entry[1][0])/2.;
	
	r_coef.entry[2] = sqrt(r1*r1+r2*r2); 

	////Remove those edges that have small R1, R2 and R3 magnitude

	//for the paper examples, we calculate as many points as possible
	//if(length(r_coef) < 0.04)
	//{
	//	cur_e->valid = 0;
	//	return false;
	//}

	normalize(r_coef);

	////Remove those edges that have weak R3 component

	//if(r_coef.entry[2] < 0.2)
	//{
	//	return false;
	//}

	////Calculate the eigen vector of the symmetric matrix of the Jacobian of the edge
	//double trace = (mat.entry[0][0] + mat.entry[1][1])/2.;
	//double antidiag = (mat.entry[0][1] + mat.entry[1][0])/2.;
	//mat.entry[0][0] = mat.entry[0][0] - trace;
	//mat.entry[1][1] = mat.entry[1][1] - trace;
	//mat.entry[0][1] = antidiag;
	//mat.entry[1][0] = antidiag;
	//CalEigenForMatrix2x2(mat, cur_e->evec);

	////Another filter to remove some edges 07/21/06
	double cr = cur_e->evec[0].entry[0]*cur_e->evec[1].entry[1]\
		-cur_e->evec[1].entry[0]*cur_e->evec[0].entry[1];

	if(abs(cr) < 1e-10)
	{
		return false;
	}

	double re1 = cur_e->evec[1].entry[0]*v1->vec.entry[1] - cur_e->evec[1].entry[1]*v1->vec.entry[0];
	double re2 = cur_e->evec[1].entry[0]*v2->vec.entry[1] - cur_e->evec[1].entry[1]*v2->vec.entry[0];

	if(re1 * re2 >= 0)
	{
		return false;
	}

	if(abs(re1) <= 1e-40)
	{
		attp[0] = v1->x;
		attp[1] = v1->y;
		return false;
	}

	if(abs(re2) <= 1e-40)
	{
		attp[0] = v2->x;
		attp[1] = v2->y;
		return false;
	}

	/* The following is not an accurate estimation */
	double ratio = abs(re1)/(abs(re2)+abs(re1));

	attp[0] = (1-ratio)*v1->x + ratio*v2->x;
	attp[1] = (1-ratio)*v1->y + ratio*v2->y;
	return true;
}




/*
Calculate the Jacobian for one vertex
*/
void CalJacobianForVertex(int vertID)
{
	int i;
	Face *face;
	icMatrix2x2 mat;
	Vertex *v = Object.vlist[vertID];

	mat.entry[0][0] = mat.entry[0][1] = mat.entry[1][0] = mat.entry[1][1] = 0.;

	for(i = 0; i < v->Num_corners; i++)
	{
		if(v->Corners[i] < 0)
			continue;

		face = Object.flist[Object.clist[v->Corners[i]]->t];
		mat = mat + face->Jacobian;
	}

	double la, lb, lc, ld, A, B, C, delta;

	la = mat.entry[0][0]/v->Num_corners;
	lb = mat.entry[0][1]/v->Num_corners;
	lc = mat.entry[1][0]/v->Num_corners;
	ld = mat.entry[1][1]/v->Num_corners;

	////save the Jacobian  07/20/06
	v->Jacobian.set(mat);  //this is still a general tensor

}


extern void cal_Curl_For_Mesh();
extern void cal_Div_For_Mesh();
extern void cal_Ang_Curl_Div();

/*
For real detection, we don't need to really calculate the e-vectors for triangles and vertices
*/
void CalJacobianForWholeMesh()
{	
	int i, j;

	for(i = 0; i < Object.nfaces; i++)
	{
		CalJacobianForTriangle(i);
	}

	/*   Calculate the eigen vectors on vertices instead of faces 07/20/06 */
	for(i = 0; i < Object.nverts; i++)
	{
		CalJacobianForVertex(i);
	}


	/* 02/12/07 */
    cal_Curl_For_Mesh(); /* Calculate the curl of the mesh */
    cal_Div_For_Mesh();  /* Calculate the divergence of the whole mesh */
	cal_Ang_Curl_Div();
}


void CalAllSpecialPointsForASCC(int scc_index)
{
	int j, k;
	Face *face;
	Edge *cur_e;
	double p[2] = {0.};

	scclist.scccomponents[scc_index].num_seppts = 0;
	scclist.scccomponents[scc_index].num_attpts = 0;

	for(j = 0; j < scclist.scccomponents[scc_index].num_nodes; j++)
	{
		////First, you need to mark all the boundary edges, and don't calculate the special points on them!
		////Now, just let me ignore that step and experiment a little bit further

		face = Object.flist[scclist.scccomponents[scc_index].nodes[j]];

		for(k = 0; k < 3; k++)
		{
			cur_e = face->edges[k];

			if(cur_e->visited == 1)
				continue;

			cur_e->valid = 0;  //the default setting is invalide

			////calculate the separation point
			if(ExistSeparationPoint_SCC(cur_e, face->index,p))
			{
				cur_e->sep.entry[0] = p[0];
				cur_e->sep.entry[1] = p[1];
				cur_e->find_sep = 1;
				cur_e->valid = 1;
				cur_e->sep_visit = 0;

				scclist.scccomponents[scc_index].num_seppts ++;
			}
			
			////calculate the attachment point
			if(ExistAttachmentPoint_SCC(cur_e, face->index,p))
			{
				cur_e->attp.entry[0] = p[0];
				cur_e->attp.entry[1] = p[1];
				cur_e->find_attp = 1;
				cur_e->valid = 1;
				cur_e->att_visit = 0;

				scclist.scccomponents[scc_index].num_attpts ++;
			}

			cur_e->visited = 1;
		}
	}
}



/*************************************************************************/

