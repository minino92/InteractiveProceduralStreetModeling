/*
This file implements the scalar field creation.
This scalar field will be used to obtain the asymmetric tensor field.

NOTE: the data structure has combined with current data structure
if the data structure is changed, corresponding changes are necessary.
*/
#include "stdafx.h"

#include <gl/glut.h> 

#include "VFDataStructure.h"

#include "tensorvis.h"

#include "tensoranalysis.h"

#include "tensordesign.h"

#include "scalardesign.h"

ScalarSingularElemList *scalar_singular_elems=NULL;

extern QuadMesh *quadmesh;
extern void  HsvRgb( float hsv[3], float rgb[3] );

/*
    compute the eigenvectors of all vertices of the quad mesh
	with the asymmetric terms being considered
*/

/*
	using new way to compute the eigen vectors of a given 2x2 symmetric tensor
	and the corresponding \phi
*/
void cal_eigen_vector_asym(icMatrix2x2 m, double phi, icVector2 ev[2])
{
	/*first get the deviator of m*/
	icMatrix2x2 dev;
	double half_trace = 0.5*(m.entry[0][0]+m.entry[1][1]);
	dev.entry[0][0] = m.entry[0][0]-half_trace;
	dev.entry[1][1] = m.entry[1][1]-half_trace;
	dev.entry[0][1] = m.entry[0][1];
	dev.entry[1][0] = m.entry[1][0];

	/*compute the eigen vectors*/
	double theta = atan2(dev.entry[0][1], dev.entry[0][0]);

	/*  according to \phi value, compute the relative direction of the 
	    major and minor eigen vectors
	*/

	phi = phi/180.*M_PI;

	double term1=sqrt(sin(phi+M_PI/4.));
	double term2=sqrt(cos(phi+M_PI/4.));
	icVector2 temp_ev[2];

	temp_ev[0].entry[0]=term1+term2;
	temp_ev[0].entry[1]=term1-term2;

	temp_ev[1].entry[0]=term1-term2;
	temp_ev[1].entry[1]=term1+term2;

	/*   rotate them with theta/2 degree  */
	double sin_theta_half=sin(theta/2.);
	double cos_theta_half=cos(theta/2.);

	ev[0].entry[0]=cos_theta_half*temp_ev[0].entry[0]-sin_theta_half*temp_ev[0].entry[1];
	ev[0].entry[1]=sin_theta_half*temp_ev[0].entry[0]+cos_theta_half*temp_ev[0].entry[1];

	ev[1].entry[0]=cos_theta_half*temp_ev[1].entry[0]-sin_theta_half*temp_ev[1].entry[1];
	ev[1].entry[1]=sin_theta_half*temp_ev[1].entry[0]+cos_theta_half*temp_ev[1].entry[1];

	normalize(ev[0]);
	normalize(ev[1]);

}


void cal_eigenvecs_onevert_quad_asym(int ver)
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

	cal_eigen_vector_asym(v->Jacobian, v->phi, ev);

	v->major = ev[0];
	v->minor = ev[1];

	/*compute the angle*/
	v->major_ang = atan2(v->major.entry[1], v->major.entry[0]);
	v->minor_ang = atan2(v->minor.entry[1], v->minor.entry[0]);

	/*save the angle of major field for degenerate points/singularity detection*/
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

void cal_all_eigenvecs_asym()
{
	int i;

	for(i=0;i<quadmesh->nverts;i++)
	{
		cal_eigenvecs_onevert_quad_asym(i);
	}
}

/*
    visualize the scalar field using color coding 
*/
void vis_scalarfield()
{
	int i,j;
	QuadCell *face;
	QuadVertex *v;

	float hsv[3]={0, 1.,1.};
	float rgb[3]={0.};

	//float max=SCALARMAX/180.*M_PI;
	//float min=SCALARMIN/180.*M_PI;
	
	float max=SCALARMAX;
	float min=SCALARMIN;


	glShadeModel(GL_SMOOTH);
	for(i=0;i<quadmesh->nfaces;i++)
	{
		face=quadmesh->quadcells[i];

		glBegin(GL_POLYGON);
		for(j=0;j<face->nverts;j++)
		{
			v=quadmesh->quad_verts[face->verts[j]];

			if(!v->inland)
				continue;

			/*   obtain the color according the phi value of the vertex   */
			if(v->phi>max)
			{
				hsv[0] = 0.;
			}
			else if(v->phi<min)
			{
				hsv[0] = 240;
			}
			else
			{
				double temp=6 * (v->phi-min)/(max-min);

				double adj=temp-floor(temp);

				hsv[0] = 240.- 240*adj;
			}

			HsvRgb(hsv, rgb);
			glColor3fv(rgb);
			glVertex2f(v->x, v->y);
		}
		glEnd();
	}
}


/*
    compute the \phi values in the mesh using the idea of the summation of basis fields
*/

double cal_phi_at(double x, double y)
{
	double dx, dy, r;
	int i;

	double phi=0;

	for(i=0;i<scalar_singular_elems->nelems;i++)
	{
		if(scalar_singular_elems->scalarsingularelems[i]->deleted)
			continue;

		dx=x-scalar_singular_elems->scalarsingularelems[i]->pos[0];
		dy=y-scalar_singular_elems->scalarsingularelems[i]->pos[1];

		r= dx*dx+dy*dy;
    	
		if (r < DistanceThreshold)   r = DistanceThreshold;

		/*  we should let user decide the decaying factor  
		    currently, we are using 25
		*/

		if(scalar_singular_elems->scalarsingularelems[i]->type==0)/*Max*/
			phi += exp(-25*r)*SCALARMAX; 
			//phi += SCALARMAX/(r*sqrt(r)); 
		   
		else
			phi += exp(-25*r)*SCALARMIN;
			//phi += SCALARMIN/(r*sqrt(r));
	}

	return phi;
}

void cal_phi_allverts()
{
	int i;
	QuadVertex *v;

	for(i=0;i<quadmesh->nverts;i++)
	{
		v=quadmesh->quad_verts[i];

		v->phi=cal_phi_at(v->x, v->y);
	}
}


void compute_phi_in_quad(int face, double x, double y, double &phi)
{

	QuadCell *qc = quadmesh->quadcells[face];

	/*get the x coeff and y coeff*/
	double a = (x-qc->x_start_coord)/quadmesh->xinterval;
	double b = (y-qc->y_start_coord)/quadmesh->yinterval;

	if(fabs(a)<1e-6)
		a = 0;
	if(fabs(b)<1e-6)
		b = 0;

	/*obtain the vertices of this cell in the order of
	  v00 v01
	  v10 v11
	*/
	QuadVertex *v00 = quadmesh->quad_verts[qc->verts[0]];
	QuadVertex *v01 = quadmesh->quad_verts[qc->verts[3]];
	QuadVertex *v10 = quadmesh->quad_verts[qc->verts[1]];
	QuadVertex *v11 = quadmesh->quad_verts[qc->verts[2]];

	/*the all the components of the interpolated tensor, respectively*/
	phi = bilinear_interpolate(a, b, v00->phi, v01->phi,
		v10->phi, v11->phi);

}



/*
    propagate the \phi values in the mesh using the constrained optimization
	NOTE: we set the vertices of the cells that contain design elements as the
	boundary vertices !
*/

void pro_phi_allverts()
{
	/*  initialization  */

	/*  obtain the inner vertices  */

	/*  construct the sparse linear system  */

	/*  solve the system  */
}



/*
    normalize the scalar field
*/

void normalize_scalar_phi()
{
	int i;
	QuadVertex *v;
	double max_phi, min_phi;

	/*  obtain the current maximum and minimum value  */
	max_phi=min_phi=quadmesh->quad_verts[0]->phi;
	for(i=1;i<quadmesh->nverts;i++)
	{
		v=quadmesh->quad_verts[i];

		if(max_phi<v->phi) max_phi=v->phi;
		if(min_phi>v->phi) min_phi=v->phi;
	}

	/*   normalize   */
	for(i=0;i<quadmesh->nverts;i++)
	{
		v=quadmesh->quad_verts[i];
		v->phi=(v->phi-min_phi)/(max_phi-min_phi)*(SCALARMAX-SCALARMIN)+SCALARMIN;
	}
}


/*
    initializing the design element list
*/

void init_scalar_singular_elemlist()
{
	if(scalar_singular_elems == NULL)
	{
		scalar_singular_elems=new ScalarSingularElemList();
	}

	scalar_singular_elems->reset();
}


void add_new_scalar_singular_elem(double x, double y, unsigned char type)
{
	ScalarSingularElem *elem=(ScalarSingularElem *)malloc(sizeof(ScalarSingularElem));
	elem->pos[0]=x;
	elem->pos[1]=y;
	elem->type=type;
	elem->deleted=false;
	scalar_singular_elems->append(elem);
}


/*
    A testing routine
*/

void test_scalar_field()
{
	init_scalar_singular_elemlist();

	/* add two new elements*/

	add_new_scalar_singular_elem(0.5, 0.5, 0);
	add_new_scalar_singular_elem(.1, .4, 1);

	/*  calculate the phi values  */
	cal_phi_allverts();

	normalize_scalar_phi();

	/*  calculate the eigen vectors */
	cal_all_eigenvecs_asym();
}