

/*
Calculate the divergence and curl of the vector field everywhere of the mesh (triangle)
*/

#include "stdafx.h"

#include "FindSepAttPoints.h"

#include "VFDataStructure.h"

extern Polygon3D Object;

/* Global maximum and minimum curl and divergence */
double max_div, min_div, max_curl, min_curl;
double max_curl_ang, min_curl_ang;

void cal_Curl_For_Mesh();
void cal_Div_For_Mesh();

void cal_Appro_Curl_Ver(int vertid);
void cal_Jacobian_Ver(int vertid);
void cal_Jacobian_All_Vers();

void decomp_Jac_ver();
void decomp_Jac_all_verts();


void cal_Curl_For_Mesh()
{
	int i;
	Face *face;
	max_curl = 0;
	min_curl = 1e40;

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];
		face->length[0] = fabs(face->Jacobian.entry[1][0]-face->Jacobian.entry[0][1]);
		if(face->length[0] > max_curl) max_curl = face->length[0];
		if(face->length[0] < min_curl) min_curl = face->length[0];
	}
}

void cal_Div_For_Mesh()
{
	int i;
	Face *face;
	max_div = 0;
	min_div = 1e40;

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];
		face->length[1] = fabs(face->Jacobian.entry[0][0]+face->Jacobian.entry[1][1]);

		if(face->length[1] > max_div) max_div = face->length[1];
		if(face->length[1] < min_div) min_div = face->length[1];
	}
}

void cal_Ang_Curl_Div()
{
	int i;
	Face *face;
	max_curl_ang = 0;
	min_curl_ang = 90;

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		if(face->length[0] < 1e-12) /*for Jacobian without curl term, set the angle to be 0 */
			face->length[2] = 0;
		else
			face->length[2] = atan2(face->length[0], face->length[1]);

		if(face->length[2] > max_curl_ang) max_curl_ang = face->length[2];
		if(face->length[2] < min_curl_ang) min_curl_ang = face->length[2];
	}
}

/*
Calculate the Jacobian for one specified vertex according to its one ring neighbors
This routine should be called after calculating the Jacobi for all triangles
*/
void cal_Jacobian_Ver(int vertid)
{
	int i;
	Vertex *v = Object.vlist[vertid];
	double a[4] = {0.};

	for(i = 0; i < v->Num_corners; i++)
	{
		a[0] += Object.flist[Object.clist[v->Corners[i]]->t]->Jacobian.entry[0][0];
		a[1] += Object.flist[Object.clist[v->Corners[i]]->t]->Jacobian.entry[0][1];
		a[2] += Object.flist[Object.clist[v->Corners[i]]->t]->Jacobian.entry[1][0];
		a[3] += Object.flist[Object.clist[v->Corners[i]]->t]->Jacobian.entry[1][1];
	}

	/*we have to use area of each triangle as the weight for irregular mesh*/
	v->Jacobian.entry[0][0] = a[0]/v->Num_corners;
	v->Jacobian.entry[0][1] = a[1]/v->Num_corners;
	v->Jacobian.entry[1][0] = a[2]/v->Num_corners;
	v->Jacobian.entry[1][1] = a[3]/v->Num_corners;
}

/*
This routine calculates the Jacobian of all the vertices of the Mesh
*/
void cal_Jacobian_All_Vers()
{
	max_curl_ang = 0;
	min_curl_ang = 90;
	max_curl = 0;
	min_curl = 1e40;
	max_div = 0;
	min_div = 1e40;

	int i;

	for(i = 0; i < Object.nverts; i++)
	{
		cal_Jacobian_Ver(i);

		cal_Appro_Curl_Ver(i);
	}
}


/*
This routine will calculate the exact curl on the specified vertex
*/
void cal_Exact_Curl_Ver(int vertid)
{
	Vertex *v = Object.vlist[vertid];
}


/*
This routine calculate the divergence and curl on a vertex according to its Jacobian
This routine should be called after calculating the Jacobian on the vertex
*/
void cal_Appro_Curl_Ver(int vertid)
{
	Vertex *v = Object.vlist[vertid];

	v->length[0] = fabs(v->Jacobian.entry[1][0] - v->Jacobian.entry[0][1]);
	v->length[1] = fabs(v->Jacobian.entry[0][0] + v->Jacobian.entry[1][1]);
	v->length[2] = atan2(v->length[0], v->length[1]);

	if(v->length[0] > max_curl) max_curl = v->length[0];
	if(v->length[0] < min_curl) min_curl = v->length[0];

	if(v->length[1] > max_div) max_div = v->length[1];
	if(v->length[1] < min_div) min_div = v->length[1];

	if(v->length[2] > max_curl_ang) max_curl_ang = v->length[2];
	if(v->length[2] < min_curl_ang) min_curl_ang = v->length[2];
}


/*We decompose the Jacobian matrix for each Vertex*/

void decomp_Jac_ver(int vertid)
{
	Vertex *v = Object.vlist[vertid];

	//icMatrix2x2 D, S, R;

	double gama_d, gama_s, gama_r;

	double la, lb, lc, ld;
	la = v->Jacobian.entry[0][0];
	lb = v->Jacobian.entry[0][1];
	lc = v->Jacobian.entry[1][0];
	ld = v->Jacobian.entry[1][1];

	gama_d = (la+ld)/2.;
	gama_r = (lc-lb)/2.;
	gama_s = sqrt((la-ld)*(la-ld)+(lb+lc)*(lb+lc))/2.;

	icVector3 temp;
	temp.entry[0] = gama_d;
	temp.entry[1] = gama_r;
	temp.entry[2] = gama_s;

	normalize(temp);
	
	v->gama_d = temp.entry[0];
	v->gama_r = temp.entry[1];
	v->gama_s = temp.entry[2];

	double theta = atan2(v->gama_d, v->gama_r);
	theta += M_PI;

	v->hue = theta*180/M_PI;
	v->saturation = 1-v->gama_s*gama_s;
}


void decomp_Jac_all_verts()
{
	int i;

	for(i = 0; i < Object.nverts; i++)
	{
		decomp_Jac_ver(i);
	}
}