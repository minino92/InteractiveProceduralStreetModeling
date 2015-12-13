/*
This file implements the tensor field visualization
Here, we assume that the tensor field has been obtained and 
assigned at each vertex.
NOTE: the data structure has combined with current data structure
if the data structure is changed, according changes are necessary.
09/17/2007
*/
#include "stdafx.h"

#include ".\glview.h"

#include "VFDataStructure.h"
#include "caldeformation.h"

#include "tensorvis.h"

extern void cal_Jacobian_Ver(int vertid);/*approximate Jacobian at each vertex*/
extern Polygon3D Object;
#define NPIX  512
#define	NPN 128
extern int     iframe ; 
extern int     Npat   ;
extern int     alpha  ;
//extern double  tmax   ;
//extern double  dmax   ;
extern int REALWINSIZE;

#define SCALE 3.
double  ten_tmax   = NPIX/(SCALE*NPN);
double  ten_dmax   = SCALE/NPIX;


int major_iframe = 0;
int minor_iframe = 0;

GLuint tentextnames[16];
GLubyte major_tex1[NPIX][NPIX][3], major_tex2[NPIX][NPIX][3], major_tex[NPIX][NPIX][3],
		minor_tex1[NPIX][NPIX][3], minor_tex2[NPIX][NPIX][3], minor_tex[NPIX][NPIX][3],
		major_alpha_map[NPIX][NPIX][3], minor_alpha_map[NPIX][NPIX][3];
	
GLubyte major_temp[NPIX][NPIX][4], minor_temp[NPIX][NPIX][4];

/*the quad mesh object for the tensor field design*/
QuadMesh *quadmesh = NULL;

bool flag_loadmap = false;

extern unsigned char *fittedmap1;
extern unsigned char *displaymap;

extern double zoom_factor;
extern double trans_x;
extern double trans_y;

extern bool is_on_local_editing;


#include "shareinterfacevars.h"
extern SharedInterfaceVars sharedvars;

#include "imgboundaryextract.h"

extern MapBoundaryList *mapboundarylist;

unsigned char final_tex[800][800][3];

/*
The following codes compute the tensor field using the Jacobian 
of the design vector field
*/
void cal_eigenvector_vertices()
{
	int i;

	for(i=0; i<Object.nverts; i++)
	{
		cal_eigenvector_perver(i);
	}

	/*normalize the major and minor field*/
    double r;
	Vertex *cur_v;

	for(i = 0; i < Object.nverts; i++)
	{
		cur_v = Object.vlist[i];

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


void cal_eigenvector_perver(int ver)
{
	Vertex *v=Object.vlist[ver];
	double evalues[2];
	icVector2 ev[2];

	/*calculate eigen values*/
	cal_eigenval_matrix2x2(v->Jacobian, evalues);

	/*calculate eigen vectors*/
	cal_eigenvectors(v->Jacobian, evalues, ev);

	v->major = ev[0];
	v->minor = ev[1];

	/*compute the angle*/
	v->major_ang = atan2(v->major.entry[1], v->major.entry[0]);
	v->minor_ang = atan2(v->minor.entry[1], v->minor.entry[0]);

	/*save the angle of major field for degenerate points/singularity detection*/

	//if(v->major_ang<0)
	//	v->tensor_major_ang = M_PI+v->major_ang;
	//else
	//	v->tensor_major_ang = v->major_ang;
		
	/*obtain the multiplied vector, we will use this obtained
	vector field to extract singularities*/
	v->tran_vec.set(cos(2*v->tensor_major_ang), sin(2*v->tensor_major_ang));
	//v->tran_vec = evalues[0]*v->tran_vec;  /*scaled by the larger eigen value*/


	/*transfer to cos^2*/
	double major_ang_cos = cos(v->major_ang);
	double major_ang_sin = sin(v->major_ang);

	v->major_ang = major_ang_cos*major_ang_cos;
	//v->major_ang = major_ang_sin*major_ang_sin;
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
void normalized_tensorfield()
{
	/*normalize the major and minor field*/
	int i;
    double r;
	Vertex *cur_v;

	for(i = 0; i < Object.nverts; i++)
	{
		cur_v = Object.vlist[i];

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


/*
The followings are codes for symmetric tensor field visualization
*/

/*
NOTE: we should also be able to transfer local tensor to a global tensor
*/

/*--------------------------------------------------------------------------------*/
//create pattern for texture mapping to create flow effect visualization
//the size of the pattern is 64x64
/*--------------------------------------------------------------------------------*/
void make_tens_Patterns(void)
{
   int lut[256];
   int phase[NPN][NPN];
   GLubyte pat[NPN][NPN][4];
   GLubyte spat[NPN][NPN][4];
   int i, j, k, t;
    
   for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
   for (i = 0; i < NPN; i++)
   for (j = 0; j < NPN; j++) phase[i][j] = rand() % 256; 

   for (k = 200; k < Npat+200; k++) {
      //t = k*256/Npat;                           //t is used to control the animation of the image
      for (i = 0; i < NPN; i++) 
      for (j = 0; j < NPN; j++) {

		  /*major*/
          //pat[i][j][0] = (GLubyte)lut[ phase[i][j] % 255];
          //pat[i][j][1] = 20;
          //pat[i][j][2] = 20;
          //pat[i][j][3] = alpha;
          pat[i][j][0] = 
          pat[i][j][1] = 
          pat[i][j][2] = (GLubyte)lut[ phase[i][j] % 255];
          pat[i][j][3] = alpha+5;
            
		  /*minor*/
		  //spat[i][j][0] = 20;
    //      spat[i][j][1] = 20;
    //      spat[i][j][2] = (GLubyte)lut[ phase[i][j] % 255];
    //      spat[i][j][3] = alpha;
		  
		  spat[i][j][0] = 
          spat[i][j][1] = 
          spat[i][j][2] = (GLubyte)lut[ phase[i][j] % 255];
          spat[i][j][3] = alpha+5;
      }

	  glNewList(k + 1, GL_COMPILE);  //major texture
      glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, 
                   GL_RGBA, GL_UNSIGNED_BYTE, pat);
      glEndList();


	  glNewList(k + 1 + 100, GL_COMPILE);       //This is for minor image
      glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, 
                   GL_RGBA, GL_UNSIGNED_BYTE, spat);
      glEndList();   
   }
}


void render_ibfv_tens(bool major_minor, bool x_y)
{
	int i, j;
	Face *face;
	int *verts;
	double px, py;

	/*major*/
	if(!major_minor) /*major field*/
	{
		if(!x_y) /*positive x*/
		{
			for (i=0; i<Object.nfaces; i++) {
				face = Object.flist[i];
				verts = face->verts;
				glBegin(GL_POLYGON);
				for (j=0; j<face->nverts; j++) {
					glTexCoord2f(Object.vlist[verts[j]]->x, Object.vlist[verts[j]]->y);

					if(Object.vlist[verts[j]]->major_cos)
					{
						px = Object.vlist[verts[j]]->x - Object.vlist[verts[j]]->major.entry[0];
						py = Object.vlist[verts[j]]->y - Object.vlist[verts[j]]->major.entry[1];
					}
					else
					{
						px = Object.vlist[verts[j]]->x + Object.vlist[verts[j]]->major.entry[0];
						py = Object.vlist[verts[j]]->y + Object.vlist[verts[j]]->major.entry[1];
					}

					glVertex2f(px, py);
				}
				glEnd();
			}
		}
		else  /*positive y*/
		{
			for (i=0; i<Object.nfaces; i++) {
				face = Object.flist[i];
				verts = face->verts;
				glBegin(GL_POLYGON);
				for (j=0; j<face->nverts; j++) {
					glTexCoord2f(Object.vlist[verts[j]]->x, Object.vlist[verts[j]]->y);
					//px = Object.vlist[verts[j]]->x + Object.vlist[verts[j]]->major.entry[0];
					//py = Object.vlist[verts[j]]->y + fabs(Object.vlist[verts[j]]->major.entry[1]);
					if(Object.vlist[verts[j]]->major_sin)
					{
						px = Object.vlist[verts[j]]->x - Object.vlist[verts[j]]->major.entry[0];
						py = Object.vlist[verts[j]]->y - Object.vlist[verts[j]]->major.entry[1];
					}
					else
					{
						px = Object.vlist[verts[j]]->x + Object.vlist[verts[j]]->major.entry[0];
						py = Object.vlist[verts[j]]->y + Object.vlist[verts[j]]->major.entry[1];
					}
					glVertex2f(px, py);
				}
				glEnd();
			}
		}
	}
	else  /*minor field*/
	{
		if(!x_y) /*positive x*/
		{
			for (i=0; i<Object.nfaces; i++) {
				face = Object.flist[i];
				verts = face->verts;
				glBegin(GL_POLYGON);
				for (j=0; j<face->nverts; j++) {
					glTexCoord2f(Object.vlist[verts[j]]->x, Object.vlist[verts[j]]->y);
					//px = Object.vlist[verts[j]]->x + fabs(Object.vlist[verts[j]]->minor.entry[0]);
					//py = Object.vlist[verts[j]]->y + Object.vlist[verts[j]]->minor.entry[1];
					if(Object.vlist[verts[j]]->minor_cos)
					{
						px = Object.vlist[verts[j]]->x - Object.vlist[verts[j]]->minor.entry[0];
						py = Object.vlist[verts[j]]->y - Object.vlist[verts[j]]->minor.entry[1];
					}
					else
					{
						px = Object.vlist[verts[j]]->x + Object.vlist[verts[j]]->minor.entry[0];
						py = Object.vlist[verts[j]]->y + Object.vlist[verts[j]]->minor.entry[1];
					}
					glVertex2f(px, py);
				}
				glEnd();
			}
		}
		else  /*positive y*/
		{
			for (i=0; i<Object.nfaces; i++) {
				face = Object.flist[i];
				verts = face->verts;
				glBegin(GL_POLYGON);
				for (j=0; j<face->nverts; j++) {
					glTexCoord2f(Object.vlist[verts[j]]->x, Object.vlist[verts[j]]->y);
					//px = Object.vlist[verts[j]]->x + Object.vlist[verts[j]]->minor.entry[0];
					//py = Object.vlist[verts[j]]->y + fabs(Object.vlist[verts[j]]->minor.entry[1]);
					if(Object.vlist[verts[j]]->minor_sin)
					{
						px = Object.vlist[verts[j]]->x - Object.vlist[verts[j]]->minor.entry[0];
						py = Object.vlist[verts[j]]->y - Object.vlist[verts[j]]->minor.entry[1];
					}
					else
					{
						px = Object.vlist[verts[j]]->x + Object.vlist[verts[j]]->minor.entry[0];
						py = Object.vlist[verts[j]]->y + Object.vlist[verts[j]]->minor.entry[1];
					}
					glVertex2f(px, py);
				}
				glEnd();
			}
		}
	}
	
	if(!major_minor)
		major_iframe ++;
	else
		minor_iframe ++;

	glEnable(GL_BLEND); 

	/*double check the following codes*/
	if(!major_minor)
		glCallList(major_iframe % Npat + 1 + 200);
	else
		glCallList(minor_iframe % Npat + 1 + 200 + 100);

    glBegin(GL_QUAD_STRIP);
		glTexCoord2f(0.0,  0.0);  glVertex2f(0.0, 0.0);
		glTexCoord2f(0.0,  ten_tmax); glVertex2f(0.0, 1.0);
		glTexCoord2f(ten_tmax, 0.0);  glVertex2f(1.0, 0.0);
		glTexCoord2f(ten_tmax, ten_tmax); glVertex2f(1.0, 1.0);
	glEnd();
}


void render_ibfv_tens_quad(bool major_minor, bool x_y)
{
	int i, j;
	QuadCell *face;
	QuadVertex *vert;
	double px, py;

	//glViewport(0, 0, (GLsizei) NPIX, (GLsizei) NPIX);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity(); 
	gluOrtho2D(0, 1, 0, 1);


	/*major*/
	if(!major_minor) /*major field*/
	{
		if(!x_y) /*positive x*/
		{
			for (i=0; i<quadmesh->nfaces; i++) {
				face = quadmesh->quadcells[i];
				glBegin(GL_POLYGON);
				for (j=0; j<face->nverts; j++) {
					vert = quadmesh->quad_verts[face->verts[j]];
					glTexCoord2f(vert->x, vert->y);

					if(vert->major_cos)
					{
						px = vert->x - vert->major.entry[0];
						py = vert->y - vert->major.entry[1];
					}
					else
					{
						px = vert->x + vert->major.entry[0];
						py = vert->y + vert->major.entry[1];
					}

					glVertex2f(px, py);
				}
				glEnd();
			}
		}
		else  /*positive y*/
		{
			for (i=0; i<quadmesh->nfaces; i++) {
				face = quadmesh->quadcells[i];
				glBegin(GL_POLYGON);
				for (j=0; j<face->nverts; j++) {
					vert = quadmesh->quad_verts[face->verts[j]];
					glTexCoord2f(vert->x, vert->y);
					if(vert->major_sin)
					{
						px = vert->x - vert->major.entry[0];
						py = vert->y - vert->major.entry[1];
					}
					else
					{
						px = vert->x + vert->major.entry[0];
						py = vert->y + vert->major.entry[1];
					}
					glVertex2f(px, py);
				}
				glEnd();
			}
		}
	}
	else  /*minor field*/
	{
		if(!x_y) /*positive x*/
		{
			for (i=0; i<quadmesh->nfaces; i++) {
				face = quadmesh->quadcells[i];
				glBegin(GL_POLYGON);
				for (j=0; j<face->nverts; j++) {
					vert = quadmesh->quad_verts[face->verts[j]];
					glTexCoord2f(vert->x, vert->y);
					if(vert->minor_cos)
					{
						px = vert->x - vert->minor.entry[0];
						py = vert->y - vert->minor.entry[1];
					}
					else
					{
						px = vert->x + vert->minor.entry[0];
						py = vert->y + vert->minor.entry[1];
					}
					glVertex2f(px, py);
				}
				glEnd();
			}
		}
		else  /*positive y*/
		{
			for (i=0; i<quadmesh->nfaces; i++) {
				face = quadmesh->quadcells[i];
				glBegin(GL_POLYGON);
				for (j=0; j<face->nverts; j++) {
					vert = quadmesh->quad_verts[face->verts[j]];
					glTexCoord2f(vert->x, vert->y);
					if(vert->minor_sin)
					{
						px = vert->x - vert->minor.entry[0];
						py = vert->y - vert->minor.entry[1];
					}
					else
					{
						px = vert->x + vert->minor.entry[0];
						py = vert->y + vert->minor.entry[1];
					}
					glVertex2f(px, py);
				}
				glEnd();
			}
		}
	}
	
	if(!major_minor)
		major_iframe ++;
	else
		minor_iframe ++;

	glEnable(GL_BLEND); 

	/*double check the following codes*/
	if(!major_minor)
		glCallList(major_iframe % Npat + 1 + 200);
	else
		glCallList(minor_iframe % Npat + 1 + 200 + 100);

    glBegin(GL_QUAD_STRIP);
		glTexCoord2f(0.0,  0.0);  glVertex2f(0.0, 0.0);
		glTexCoord2f(0.0,  ten_tmax); glVertex2f(0.0, 1.0);
		glTexCoord2f(ten_tmax, 0.0);  glVertex2f(1.0, 0.0);
		glTexCoord2f(ten_tmax, ten_tmax); glVertex2f(1.0, 1.0);
	glEnd();
	
	glPopMatrix();

}


extern bool is_not_inland(int);

/*for later blending*/
void render_tensor_blend()
{
	/*use the mesh to display the texture instead 09/21/2007*/
	//glClearColor(1, 1, 1, 1);
	//glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	
	glBegin(GL_QUAD_STRIP);
		glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
		glTexCoord2f(0.0, 1.0);  glVertex2f(0.0, 1.0);
		glTexCoord2f(1.0, 0.0);  glVertex2f(1.0, 0.0);
		glTexCoord2f(1.0, 1.0);  glVertex2f(1.0, 1.0);
	glEnd();

	//int i, j;
	//QuadCell *face;
	//QuadVertex *v;
	//for(i=0; i<quadmesh->nfaces; i++)
	//{
	//	face=quadmesh->quadcells[i];

	//	glBegin(GL_POLYGON);
	//	for(j=0; j<face->nverts; j++)
	//	{
	//		v=quadmesh->quad_verts[face->verts[j]];
	//		//if(v->inland)
	//			glTexCoord2f(v->x, v->y);
	//		glVertex2f(v->x, v->y);
	//	}
	//	glEnd();
	//}
}


void render_tensor_final(unsigned char majororminor)
{
	/*use the mesh to display the texture instead 09/21/2007*/
	glClearColor(1, 1, 1, 1);
	//glClearColor(0.6, .7, 0.9, 1.);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	
	if(majororminor == 0)
	{
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
			GL_RGB, GL_UNSIGNED_BYTE, major_tex);
	}
	else if(majororminor == 1)
	{
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
			GL_RGB, GL_UNSIGNED_BYTE, minor_tex);
	}

	glBegin(GL_QUAD_STRIP);
		glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
		glTexCoord2f(0.0, 1.0);  glVertex2f(0.0, 1.0);
		glTexCoord2f(1.0, 0.0);  glVertex2f(1.0, 0.0);
		glTexCoord2f(1.0, 1.0);  glVertex2f(1.0, 1.0);
	glEnd();
}




void render_tensor_blend(unsigned char majororminor)
{
	/*use the mesh to display the texture instead 09/21/2007*/
	glDrawBuffer(GL_BACK);
	glClearColor(1, 1, 1, 1.);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	

	if(majororminor == 0)
	{
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
			GL_RGB, GL_UNSIGNED_BYTE, major_tex);
	}
	else if(majororminor == 1)
	{
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
			GL_RGB, GL_UNSIGNED_BYTE, minor_tex);
	}

	int i, j;
	QuadCell *face;
	QuadVertex *v;
	for(i=0; i<quadmesh->nfaces; i++)
	{
		face=quadmesh->quadcells[i];

		glBegin(GL_POLYGON);
		for(j=0; j<face->nverts; j++)
		{
			v=quadmesh->quad_verts[face->verts[j]];
			//if(!v->inland) continue;

			glTexCoord2f(v->x, v->y);
			glVertex2f(v->x, v->y);
		}
		glEnd();
	}

	/*How about we blend the map here 11/10/2007*/
	//////glEnable(GL_BLEND);
	//////glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
	//////	GL_RGB, GL_UNSIGNED_BYTE, fittedmap1);
	//////glBegin(GL_QUAD_STRIP);
	//////	glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
	//////	glTexCoord2f(0.0, 1.0);  glVertex2f(0.0, 1.0);
	//////	glTexCoord2f(1.0, 0.0);  glVertex2f(1.0, 0.0);
	//////	glTexCoord2f(1.0, 1.0);  glVertex2f(1.0, 1.0);
	//////glEnd();
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_COLOR_MATERIAL);
	glShadeModel(GL_SMOOTH);

	for(i=0; i<quadmesh->nfaces; i++)
	{
		face=quadmesh->quadcells[i];

		glBegin(GL_POLYGON);
		for(j=0; j<face->nverts; j++)
		{
			v=quadmesh->quad_verts[face->verts[j]];
			if(v->inland)
				glColor4f(1.,1.,1., 0.);
			else
				glColor4f(0.6, .7, 0.9, 1.);
			glVertex2f(v->x, v->y);
		}
		glEnd();
	}
	glShadeModel(GL_FLAT);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_BLEND);
	glEnable(GL_TEXTURE_2D);

}


/**/
void render_tensor_blend_localEdit(unsigned char majororminor)
{
	/*use the mesh to display the texture instead 09/21/2007*/
	glDrawBuffer(GL_BACK);
	glClearColor(1, 1, 1, 1.);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	

	if(majororminor == 0)
	{
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
			GL_RGB, GL_UNSIGNED_BYTE, major_tex);
	}
	else if(majororminor == 1)
	{
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
			GL_RGB, GL_UNSIGNED_BYTE, minor_tex);
	}

	int i, j;
	QuadCell *face;
	QuadVertex *v;
	for(i=0; i<quadmesh->nfaces; i++)
	{
		face=quadmesh->quadcells[i];

		glBegin(GL_POLYGON);
		for(j=0; j<face->nverts; j++)
		{
			v=quadmesh->quad_verts[face->verts[j]];
			//if(!v->inland) continue;

			glTexCoord2f(v->x, v->y);
			glVertex2f(v->x, v->y);
		}
		glEnd();
	}

	/*How about we blend the map here 11/10/2007*/
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_COLOR_MATERIAL);
	glShadeModel(GL_SMOOTH);

	for(i=0; i<quadmesh->nfaces; i++)
	{
		face=quadmesh->quadcells[i];

		glBegin(GL_POLYGON);
		for(j=0; j<face->nverts; j++)
		{
			v=quadmesh->quad_verts[face->verts[j]];
			if(v->InRegion)
				glColor4f(1.,1.,1., 0.);
			else
				glColor4f(0.82, .81, 0.8, 1.);
			glVertex2f(v->x, v->y);
		}
		glEnd();
	}
	glShadeModel(GL_FLAT);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_BLEND);
	glEnable(GL_TEXTURE_2D);

}


/*Here we assume the size of the map is always 512x512*/
void render_a_map(unsigned char *map)
{
	glClearColor(1, 1, 1, 1);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glEnable(GL_TEXTURE_2D);

	if(map==NULL) return;

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
			GL_RGB, GL_UNSIGNED_BYTE, map);

	glBegin(GL_QUAD_STRIP);
		glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
		glTexCoord2f(0.0, 1.0);  glVertex2f(0.0, 1.0);
		glTexCoord2f(1.0, 0.0);  glVertex2f(1.0, 0.0);
		glTexCoord2f(1.0, 1.0);  glVertex2f(1.0, 1.0);
	glEnd();

	if(sharedvars.ShowExtractedBoundsOn && mapboundarylist!=NULL)
	{
		glDisable(GL_TEXTURE_2D);
		glEnable(GL_COLOR_MATERIAL);
		//glColor3f(.8,.8,.8);
		glLineWidth(1.5);
		mapboundarylist->display_boundaries(GL_RENDER);
		glLineWidth(1.);
	}
}


/*
In this routine, we render the alpha map into another texture
For simplicity, we render only the alpha map for the 
positive x direction.
NOTE: this will need to render only once when a tensor field is created
but it needs to be updated after the tensor field is changed
major_minor: false--major, true--minor
*/
void render_alpha_map(bool major_minor)
{
	/**/
	int i, j;
	Face *face;
	int *verts;

	//glEnable(GL_BLEND);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
	glShadeModel(GL_SMOOTH);
	if(!major_minor)
	{
		for (i=0; i<Object.nfaces; i++) {
			face = Object.flist[i];
			verts = face->verts;
			glBegin(GL_POLYGON);
			for (j=0; j<face->nverts; j++) {
				//glColor4f(Object.vlist[verts[j]]->major_ang,
				//	Object.vlist[verts[j]]->major_ang,
				//	Object.vlist[verts[j]]->major_ang,
				//	Object.vlist[verts[j]]->major_ang);
				glColor3f(Object.vlist[verts[j]]->major_ang,
					Object.vlist[verts[j]]->major_ang,
					Object.vlist[verts[j]]->major_ang);
				
				//GLubyte color=(GLubyte)(Object.vlist[verts[j]]->major_ang*255);
				//glColor3b(color, color, color);
				glVertex2f(Object.vlist[verts[j]]->x, Object.vlist[verts[j]]->y);
			}
			glEnd();
		}

		/*copy to the major_alpha_map*/
		glReadBuffer(GL_BACK);
		glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, major_alpha_map);
	}
	else
	{
		for (i=0; i<Object.nfaces; i++) {
			face = Object.flist[i];
			verts = face->verts;
			glBegin(GL_POLYGON);
			for (j=0; j<face->nverts; j++) {
				//glColor4f(Object.vlist[verts[j]]->minor_ang,
				//	Object.vlist[verts[j]]->minor_ang,
				//	Object.vlist[verts[j]]->minor_ang,
				//	Object.vlist[verts[j]]->minor_ang);
				glColor3f(Object.vlist[verts[j]]->minor_ang,
					Object.vlist[verts[j]]->minor_ang,
					Object.vlist[verts[j]]->minor_ang);
				glVertex2f(Object.vlist[verts[j]]->x, Object.vlist[verts[j]]->y);
			}
			glEnd();
		}
		
		/*copy to the minor_alpha_map*/
		glReadBuffer(GL_BACK);
		glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, minor_alpha_map);
	}
	glDisable(GL_BLEND);
}


void render_alpha_map_quad(bool major_minor)
{
	/**/
	int i, j;
	QuadCell *face;
	QuadVertex *vert;

	glViewport(0, 0, (GLsizei)NPIX, (GLsizei)NPIX);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_BLEND);
	glShadeModel(GL_SMOOTH);
	if(!major_minor)
	{
		for (i=0; i<quadmesh->nfaces; i++) {
			face = quadmesh->quadcells[i];
			glBegin(GL_POLYGON);
			for (j=0; j<face->nverts; j++) {
				vert = quadmesh->quad_verts[face->verts[j]];
				glColor3f(vert->major_ang,
					vert->major_ang,
					vert->major_ang);
				
				glVertex2f(vert->x, vert->y);
			}
			glEnd();
		}

		/*copy to the major_alpha_map*/
		glReadBuffer(GL_BACK);
		glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, major_alpha_map);
	}
	else
	{
		for (i=0; i<quadmesh->nfaces; i++) {
			face = quadmesh->quadcells[i];
			glBegin(GL_POLYGON);
			for (j=0; j<face->nverts; j++) {
				vert = quadmesh->quad_verts[face->verts[j]];
				glColor3f(vert->minor_ang,
					vert->minor_ang,
					vert->minor_ang);
				glVertex2f(vert->x, vert->y);
			}
			glEnd();
		}
		
		/*copy to the minor_alpha_map*/
		glReadBuffer(GL_BACK);
		glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, minor_alpha_map);
	}
	glDisable(GL_BLEND);
	glViewport(0, 0, (GLsizei)REALWINSIZE, (GLsizei)REALWINSIZE);
	//glViewport(0, 0, (GLsizei)(REALWINSIZE*zoom_factor), (GLsizei)(REALWINSIZE*zoom_factor));
}

/*
generate texture objects
*/
void init_texture_objs()
{
	/*generate texture names*/
	glGenTextures(8, tentextnames);

	/*bind the texture maps*/
	glBindTexture(GL_TEXTURE_2D, tentextnames[0]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, major_tex1);

	glBindTexture(GL_TEXTURE_2D, tentextnames[1]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, major_tex2);


	glBindTexture(GL_TEXTURE_2D, tentextnames[2]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, major_tex);

	glBindTexture(GL_TEXTURE_2D, tentextnames[3]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, minor_tex1);

	glBindTexture(GL_TEXTURE_2D, tentextnames[4]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, minor_tex2);


	glBindTexture(GL_TEXTURE_2D, tentextnames[5]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, minor_tex);
}

/*initialize the textures*/
void tensor_init_tex()
{
	glDrawBuffer(GL_BACK);

	glCallList(201);

    glBegin(GL_QUAD_STRIP);
		glTexCoord2f(0.0,  0.0);  glVertex2f(0.0, 0.0);
		glTexCoord2f(0.0,  ten_dmax); glVertex2f(0.0, 1.0);
		glTexCoord2f(ten_dmax, 0.0);  glVertex2f(1.0, 0.0);
		glTexCoord2f(ten_dmax, ten_dmax); glVertex2f(1.0, 1.0);
	glEnd();
	
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, major_tex1);

	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, major_tex2);
	
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, major_tex);


	glCallList(301);

    glBegin(GL_QUAD_STRIP);
		glTexCoord2f(0.0,  0.0);  glVertex2f(0.0, 0.0);
		glTexCoord2f(0.0,  ten_dmax); glVertex2f(0.0, 1.0);
		glTexCoord2f(ten_dmax, 0.0);  glVertex2f(1.0, 0.0);
		glTexCoord2f(ten_dmax, ten_dmax); glVertex2f(1.0, 1.0);
	glEnd();
	
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, minor_tex1);

	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, minor_tex2);
	
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, minor_tex);
}

void clear_all_tens()
{
	int i;
	for(i=0; i<quadmesh->nverts; i++)
	{
		quadmesh->quad_verts[i]->Jacobian.set(0.);
		quadmesh->quad_verts[i]->major.set(0.,0.);
		quadmesh->quad_verts[i]->minor.set(0.,0.);
	}
}

void reset_inland()
{
	int i;
	for(i=0; i<quadmesh->nverts; i++)
	{
		quadmesh->quad_verts[i]->inland=true;
	}
}

void vis_alpha_map()
{
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGB, GL_UNSIGNED_BYTE, major_alpha_map);
	render_tensor_blend();
}

/*
we calculate the texture for major field
*/

void major_vis()
{
	/*rendering the positive x direction*/
	glBindTexture(GL_TEXTURE_2D, tentextnames[0]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGB, GL_UNSIGNED_BYTE, major_tex1);

	render_ibfv_tens(false, false);

	/*save image*/
	glDisable(GL_BLEND);
	//glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 
	//				0, 0, NPIX, NPIX, 0);
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, major_tex1);

	/*assign the alpha values for each texels according to the alpha map*/
	//int i, j;
	//for(i=0; i<NPIX; i++)
	//	for(j=0; j<NPIX; j++)
	//	{
	//		major_tex1[i][j][3] = major_alpha_map[i][j][0];
	//	}
	/***************************************************************/

	/*rendering the positive y direction*/
	major_iframe--;
	glBindTexture(GL_TEXTURE_2D, tentextnames[1]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGB, GL_UNSIGNED_BYTE, major_tex2);

	render_ibfv_tens(false, true);
	glDisable(GL_BLEND);
	//glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 
	//				0, 0, NPIX, NPIX, 0);
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, major_tex2);
	//for(i=0; i<NPIX; i++)
	//	for(j=0; j<NPIX; j++)
	//	{
	//		major_tex2[i][j][3] = 255-major_alpha_map[i][j][0];
	//	}
	/***************************************************************/

	/*blend them*/

	int i, j;
	for(i=0; i<NPIX; i++) /*y direction*/
		for(j=0; j<NPIX; j++)
		{
			major_temp[i][j][0] = major_tex1[i][j][0];
			major_temp[i][j][1] = major_tex1[i][j][1];
			major_temp[i][j][2] = major_tex1[i][j][2];
			major_temp[i][j][3] = major_alpha_map[i][j][0];
		}

	for(i=0; i<NPIX; i++) /*x direction*/
		for(j=0; j<NPIX; j++)
		{
			minor_temp[i][j][0] = major_tex2[i][j][0];
			minor_temp[i][j][1] = major_tex2[i][j][1];
			minor_temp[i][j][2] = major_tex2[i][j][2];
			minor_temp[i][j][3] = 255-major_alpha_map[i][j][0];
		}

	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
	//	GL_RGB, GL_UNSIGNED_BYTE, major_tex1);
		
	glShadeModel(GL_SMOOTH);

	glEnable(GL_BLEND);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, major_temp);
	render_tensor_blend();
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, minor_temp);
	render_tensor_blend();

	////glEnable(GL_BLEND);
	////glBindTexture(GL_TEXTURE_2D, tentextnames[0]);
	////glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NPIX, NPIX, 0,
	////	GL_RGBA, GL_UNSIGNED_BYTE, major_tex1);
	////render_tensor_blend();
	////glBindTexture(GL_TEXTURE_2D, tentextnames[1]);
	////glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NPIX, NPIX, 0,
	////	GL_RGBA, GL_UNSIGNED_BYTE, major_tex2);
	////render_tensor_blend();
	////glBindTexture(GL_TEXTURE_2D, tentextnames[2]);
	////glDisable(GL_BLEND);

	//////glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 
	//////				0, 0, NPIX, NPIX, 0);
	////

	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, major_tex);



	/*Let's use hard code*/

	//for(i=0; i<NPIX; i++)
	//{
	//	for(j=0; j<NPIX; j++)
	//	{
	//		major_tex[i][j][0]=
	//		major_tex[i][j][1]=
	//		major_tex[i][j][2]=(unsigned char)((255.-major_alpha_map[i][j][0])/255.*major_tex1[i][j][0]+
	//			major_alpha_map[i][j][0]/255.*major_tex2[i][j][0]);
	//	}
	//}

	

	/*try to improve quality of the texture.
	Work for gray texture only now
	09/18/2007*/
	for(i=0; i<NPIX; i++)
	{
		for(j=0; j<NPIX; j++)
		{
			unsigned int sum = major_tex[i][j][0]+major_tex[i][j][1]+major_tex[i][j][2];
			if(sum<125)
			{
				double coeff = (2+(125.-sum)/125.);
				major_tex[i][j][0]=
				major_tex[i][j][1]=
				major_tex[i][j][2]=(unsigned char)major_tex[i][j][0]/coeff;
			}
		}
	}
}


void major_vis_quad()
{
	/*reset the view point here 11/09/2007*/
	glViewport(0, 0, (GLsizei)NPIX, (GLsizei)NPIX);

	/*rendering the positive x direction*/
	glBindTexture(GL_TEXTURE_2D, tentextnames[0]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGB, GL_UNSIGNED_BYTE, major_tex1);

	render_ibfv_tens_quad(false, false);

	/*save image*/
	glDisable(GL_BLEND);
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, major_tex1);

	/***************************************************************/

	/*rendering the positive y direction*/
	major_iframe--;
	glBindTexture(GL_TEXTURE_2D, tentextnames[1]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGB, GL_UNSIGNED_BYTE, major_tex2);

	render_ibfv_tens_quad(false, true);

	glDisable(GL_BLEND);
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, major_tex2);
	/***************************************************************/

	/*blend them*/

	int i, j;
	for(i=0; i<NPIX; i++) /*y direction*/
		for(j=0; j<NPIX; j++)
		{
			major_temp[i][j][0] = major_tex1[i][j][0];
			major_temp[i][j][1] = major_tex1[i][j][1];
			major_temp[i][j][2] = major_tex1[i][j][2];
			major_temp[i][j][3] = major_alpha_map[i][j][0];
		}

	for(i=0; i<NPIX; i++) /*x direction*/
		for(j=0; j<NPIX; j++)
		{
			minor_temp[i][j][0] = major_tex2[i][j][0];
			minor_temp[i][j][1] = major_tex2[i][j][1];
			minor_temp[i][j][2] = major_tex2[i][j][2];
			minor_temp[i][j][3] = 255-major_alpha_map[i][j][0];
		}

		
	//glShadeModel(GL_SMOOTH);

	glEnable(GL_BLEND);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, major_temp);
	render_tensor_blend();
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, minor_temp);
	render_tensor_blend();

	////

	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, major_tex);
		
	glClearColor(1, 1, 1, 1);
	glClear(GL_COLOR_BUFFER_BIT);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGBA, GL_UNSIGNED_BYTE, major_temp);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGBA, GL_UNSIGNED_BYTE, minor_temp);

	glViewport(0, 0, (GLsizei)REALWINSIZE, (GLsizei)REALWINSIZE);
	//glViewport(0, 0, (GLsizei)(REALWINSIZE*zoom_factor), (GLsizei)(REALWINSIZE*zoom_factor));
}


void render_majorfield()
{
	major_vis();
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGB, GL_UNSIGNED_BYTE, major_tex);
	render_tensor_blend();
}



void render_majorfield_quad()
{

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity(); 
	gluOrtho2D(0, 1, 0, 1);
	
	major_vis_quad();

	glPopMatrix();

	//glDisable(GL_BLEND);


	if(!flag_loadmap)
	{
		if(!is_on_local_editing)
			render_tensor_final(0);
		else
			render_tensor_blend_localEdit(0);

	}
	else
	{
		if(sharedvars.ShowTheMapOn)
			render_a_map(displaymap);
		else
			render_tensor_blend(0);
	}
	
	//glReadBuffer(GL_BACK);
	//glReadPixels(0, 0, 800, 800, GL_RGB, GL_UNSIGNED_BYTE, final_tex);

	//glClearColor(1, 1, 1, 1);
	//glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	//glEnable(GL_TEXTURE_2D);
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 800, 800, 0,
	//		GL_RGB, GL_UNSIGNED_BYTE, final_tex);

	//glBegin(GL_QUAD_STRIP);
	//	glTexCoord2f(0.0, 0.0);  glVertex2f(0.0, 0.0);
	//	glTexCoord2f(0.0, 1.0);  glVertex2f(0.0, 1.0);
	//	glTexCoord2f(1.0, 0.0);  glVertex2f(1.0, 0.0);
	//	glTexCoord2f(1.0, 1.0);  glVertex2f(1.0, 1.0);
	//glEnd();


	glDisable(GL_TEXTURE_2D);
}


void render_minorfield()
{
	minor_vis();
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGB, GL_UNSIGNED_BYTE, minor_tex);
	render_tensor_blend();
}


void render_minorfield_quad()
{
	minor_vis_quad();
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGB, GL_UNSIGNED_BYTE, minor_tex);
	render_tensor_blend();
}

/*
we calculate the texture for minor field
*/

void minor_vis()
{
	/*rendering the positive x direction*/
	glBindTexture(GL_TEXTURE_2D, tentextnames[3]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGB, GL_UNSIGNED_BYTE, minor_tex1);

	render_ibfv_tens(true, false);

	/*save image*/
	glDisable(GL_BLEND);
	//glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 
	//				0, 0, NPIX, NPIX, 0);
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, minor_tex1);

	/***************************************************************/

	/*rendering the positive y direction*/
	minor_iframe--;
	glBindTexture(GL_TEXTURE_2D, tentextnames[4]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGB, GL_UNSIGNED_BYTE, minor_tex2);

	render_ibfv_tens(true, true);
	glDisable(GL_BLEND);
	//glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 
	//				0, 0, NPIX, NPIX, 0);
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, minor_tex2);
	/***************************************************************/

	/*blend them*/

	////int i, j;
	////for(i=0; i<NPIX; i++)
	////	for(j=0; j<NPIX; j++)
	////	{
	////		major_temp[i][j][0] = minor_tex1[i][j][0];
	////		major_temp[i][j][1] = minor_tex1[i][j][1];
	////		major_temp[i][j][2] = minor_tex1[i][j][2];
	////		major_temp[i][j][3] = minor_alpha_map[i][j][0];
	////	}

	////glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
	////	GL_RGB, GL_UNSIGNED_BYTE, minor_tex2);
	////render_tensor_blend();
	////glEnable(GL_BLEND);
	////glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NPIX, NPIX, 0,
	////	GL_RGBA, GL_UNSIGNED_BYTE, major_temp);
	////render_tensor_blend();
	
	int i, j;
	for(i=0; i<NPIX; i++) /*x direction*/
		for(j=0; j<NPIX; j++)
		{
			major_temp[i][j][0] = minor_tex1[i][j][0];
			major_temp[i][j][1] = minor_tex1[i][j][1];
			major_temp[i][j][2] = minor_tex1[i][j][2];
			major_temp[i][j][3] = minor_alpha_map[i][j][0];
		}

	for(i=0; i<NPIX; i++) /*y direction*/
		for(j=0; j<NPIX; j++)
		{
			minor_temp[i][j][0] = minor_tex2[i][j][0];
			minor_temp[i][j][1] = minor_tex2[i][j][1];
			minor_temp[i][j][2] = minor_tex2[i][j][2];
			minor_temp[i][j][3] = 255-minor_alpha_map[i][j][0];
		}


	glEnable(GL_BLEND);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, major_temp);
	render_tensor_blend();
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, minor_temp);
	render_tensor_blend();

	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, minor_tex);
	
	/*try to improve quality of the texture.
	Work for gray texture only now
	09/18/2007*/
	for(i=0; i<NPIX; i++)
	{
		for(j=0; j<NPIX; j++)
		{
			unsigned int sum = minor_tex[i][j][0]+minor_tex[i][j][1]+minor_tex[i][j][2];
			if(sum<125)
			{
				double coeff = (2+(125.-sum)/125.);
				minor_tex[i][j][0]=
				minor_tex[i][j][1]=
				minor_tex[i][j][2]=(unsigned char)minor_tex[i][j][0]/coeff;
			}
		}
	}
}

void minor_vis_quad()
{
	glViewport(0, 0, (GLsizei)NPIX, (GLsizei)NPIX);
	
	/*rendering the positive x direction*/
	glBindTexture(GL_TEXTURE_2D, tentextnames[3]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGB, GL_UNSIGNED_BYTE, minor_tex1);

	render_ibfv_tens_quad(true, false);

	/*save image*/
	glDisable(GL_BLEND);
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, minor_tex1);

	/***************************************************************/

	/*rendering the positive y direction*/
	minor_iframe--;
	glBindTexture(GL_TEXTURE_2D, tentextnames[4]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGB, GL_UNSIGNED_BYTE, minor_tex2);

	render_ibfv_tens_quad(true, true);
	glDisable(GL_BLEND);
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, minor_tex2);
	/***************************************************************/

	
	int i, j;
	for(i=0; i<NPIX; i++) /*x direction*/
		for(j=0; j<NPIX; j++)
		{
			major_temp[i][j][0] = minor_tex1[i][j][0];
			major_temp[i][j][1] = minor_tex1[i][j][1];
			major_temp[i][j][2] = minor_tex1[i][j][2];
			major_temp[i][j][3] = minor_alpha_map[i][j][0];
		}

	for(i=0; i<NPIX; i++) /*y direction*/
		for(j=0; j<NPIX; j++)
		{
			minor_temp[i][j][0] = minor_tex2[i][j][0];
			minor_temp[i][j][1] = minor_tex2[i][j][1];
			minor_temp[i][j][2] = minor_tex2[i][j][2];
			minor_temp[i][j][3] = 255-minor_alpha_map[i][j][0];
		}


	glEnable(GL_BLEND);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, major_temp);
	render_tensor_blend();
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, minor_temp);
	render_tensor_blend();

	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, NPIX, NPIX, GL_RGB, GL_UNSIGNED_BYTE, minor_tex);
	
	glViewport(0, 0, (GLsizei)REALWINSIZE, (GLsizei)REALWINSIZE);
}

/*
this will show both major and minor field
*/

void mix_vis()
{
	/*obtain major*/
	major_vis();

	/*obtain minor*/
	minor_vis();

	/*blend them*/
	int i, j;
	for(i=0; i<NPIX; i++)
		for(j=0; j<NPIX; j++)
		{
			major_temp[i][j][0] = minor_tex[i][j][0];
			major_temp[i][j][1] = minor_tex[i][j][1];
			major_temp[i][j][2] = minor_tex[i][j][2];
			major_temp[i][j][3] = 130;
		}
	
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGB, GL_UNSIGNED_BYTE, major_tex);
	render_tensor_blend();
	
	glEnable(GL_BLEND);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, major_temp);
	render_tensor_blend();
}





void mix_vis_quad()
{
	/*obtain major*/
	major_vis_quad();

	/*obtain minor*/
	minor_vis_quad();

	/*blend them*/
	int i, j;
	for(i=0; i<NPIX; i++)
		for(j=0; j<NPIX; j++)
		{
			major_temp[i][j][0] = minor_tex[i][j][0];
			major_temp[i][j][1] = minor_tex[i][j][1];
			major_temp[i][j][2] = minor_tex[i][j][2];
			major_temp[i][j][3] = 130;
		}
	
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, NPIX, NPIX, 0,
		GL_RGB, GL_UNSIGNED_BYTE, major_tex);
	render_tensor_blend();
	
	glEnable(GL_BLEND);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, NPIX, NPIX, 0,
		GL_RGBA, GL_UNSIGNED_BYTE, major_temp);
	render_tensor_blend();
}





/*display all the quad mesh*/
void display_quadmesh(GLenum mode)
{
	QuadCell *face;
	QuadVertex *vert;
	int i, j;

	/*display all the quad mesh first*/
	glEnable(GL_BLEND);
	glColor4f(.8, 0, 0.6, 0.3);
	for(i=0;i<quadmesh->nfaces;i++)
	{
		face = quadmesh->quadcells[i];
		if(!face->OnBoundary)
			continue;

		glBegin(GL_POLYGON);
		for(j = 0; j < face->nverts; j++)
		{
			vert = quadmesh->quad_verts[face->verts[j]];
			glVertex2f(vert->x, vert->y);
		}
		glEnd();
	}

	glColor3f(.8,.6, 0);
	for(i=0;i<quadmesh->nfaces;i++)
	{
		if(mode == GL_SELECT)
			glLoadName(NAMEOFTRIANGLE+i);

		face = quadmesh->quadcells[i];
		glBegin(GL_LINE_LOOP);
		for(j = 0; j < face->nverts; j++)
		{
			vert = quadmesh->quad_verts[face->verts[j]];
			glVertex2f(vert->x, vert->y);
		}
		glEnd();
	}

	glDisable(GL_BLEND);

}

void display_quadmesh_select(GLenum mode)
{
	QuadCell *face;
	QuadVertex *vert;
	int i, j;

	/*display all the triangle mesh first*/
	glColor3f(.8,.6, 0);
	for(i=0;i<quadmesh->nfaces;i++)
	{
		if(mode == GL_SELECT)
			glLoadName(NAMEOFTRIANGLE+i);

		face = quadmesh->quadcells[i];
		glBegin(GL_POLYGON);
		for(j = 0; j < face->nverts; j++)
		{
			vert = quadmesh->quad_verts[face->verts[j]];
			glVertex2f(vert->x, vert->y);
		}
		glEnd();
	}
}


/*********************************************************************************/

/*The following will implement the generation of a regular quadratic mesh 09/25/2007*/
/*the following code we generate the regular quadractic mesh 09/25/2007
NOTE: the result will be saved into the "Object" variable directly.
*/
//typedef struct QuadVer{
//	int id;      /*id = i*ydim+j*/
//	double x, y;
//}


QuadMesh::QuadMesh()
{
	quad_verts = NULL;
	quadcells = NULL;
	elist = NULL;
	edgelist = NULL;

	nverts = nfaces = nedges = 0;   //number of vertices, faces, edges respectively
}


QuadMesh::QuadMesh(int xdim, int ydim, double xstart, double xend, 
				   double ystart, double yend)
{
	quad_verts = NULL;
	quadcells = NULL;
	elist = NULL;
	edgelist = NULL;

	nverts = nfaces = nedges = 0;   //number of vertices, faces, edges respectively

	gen_regquad_mesh(xdim, ydim, xstart, xend, ystart, yend);

	/*construct edge list*/
	construct_edges();
	orient_edges_cells();

	XDIM = xdim;
	YDIM = ydim;
	this->xstart = xstart;
	this->xend = xend;
	this->ystart = ystart;
	this->yend = yend;
}


void QuadMesh::finalize_quad_verts()
{
	int i;
	if(quad_verts != NULL)
	{
		for(i=0; i<nverts; i++)
		{
			if(quad_verts[i]!=NULL)
			{
				if(quad_verts[i]->edges != NULL)
					free(quad_verts[i]->edges);
				free(quad_verts[i]);
			}
		}
		free(quad_verts);
	}
}

/*generate a regular quad mesh*/
bool QuadMesh::gen_regquad_mesh(int xdim, int ydim, double xstart, double xend, 
				   double ystart, double yend)
{
	/*generate the vertex list first*/
	if(gen_regquad_vertices(xdim, ydim, xstart, xend, ystart, yend))
	{
		/*construct the faces to obtain the corresponding connectivity information*/
		gen_regquad_faces(xdim, ydim);

		return true;
	}
	else
		return false;
}


/*generate the vertex list of the quad mesh*/
bool QuadMesh::gen_regquad_vertices(int xdim, int ydim, double xstart, double xend, 
				   double ystart, double yend)
{
	if(xstart > xend || yend < ystart) return false;
	int i, j;
	double d_x, d_y, x_coord, y_coord;
	xinterval = d_x = (xend-xstart)/(xdim-1);
	yinterval = d_y = (yend-ystart)/(ydim-1);
	int quad_vert_id;

	/*first, allocate the memory for the quad mesh*/
	if(quad_verts != NULL)
		finalize_quad_verts();
	nverts = xdim*ydim;
	quad_verts = (QuadVertex**)malloc(sizeof(QuadVertex*)*nverts);
	for(i=0; i<nverts; i++)
		quad_verts[i] = (QuadVertex*)malloc(sizeof(QuadVertex));

	for(j=0; j<ydim; j++)
	{
		y_coord = ystart+j*d_y; /*obtain y coordinates for jth row*/
		for(i=0; i<xdim; i++)
		{
			quad_vert_id = j*xdim+i;
			x_coord = xstart+i*d_x; /*obtain x coordinates for ith column*/
			quad_verts[quad_vert_id]->index = quad_vert_id;
			quad_verts[quad_vert_id]->x = x_coord;
			quad_verts[quad_vert_id]->y = y_coord;
		}
	}

	init_vertices();

	return true;
}


/*intialize the vertex list in the very beginning*/
void QuadMesh::init_vertices()
{
	int i;
	QuadVertex *qv;
	for(i=0; i<nverts; i++)
	{
		qv = quad_verts[i];
		qv->edges = NULL;
		qv->Num_edge = 0;
		qv->InRegion = false;
		qv->distance = 1.e50;
		qv->OnBoundary = false;
		qv->RegionListID = -1;
		//qv->repell_flag = qv->attract_flag = 0;
		qv->visited = false;
		qv->Jacobian.set(0.);
		qv->major.set(0.);
		qv->minor.set(0.);
		qv->ncells = 0;
		qv->cells = NULL;

		qv->inland = true;
		qv->inveg = false;

		qv->which_region=0;
		qv->inbrushregion=false;

		qv->phi=0.;

		qv->density=1.;
	}
}


/*finalize the cell list
*/
void QuadMesh::finalize_quad_cells()
{
	int i;
	if(quadcells != NULL)
	{
		for(i=0; i<nfaces; i++)
		{
			if(quadcells[i] != NULL)
			{
				free(quadcells[i]->verts);
				free(quadcells[i]);
			}
		}
		free(quadcells);
		quadcells=NULL;
	}
}


/*compute the faces/connectivities of the regular quad mesh
NOTE: this routine should be called after generating the vertex list
*/
void QuadMesh::gen_regquad_faces(int xdim, int ydim)
{
	int i, j;
	QuadCell *f;

	/*allocate space for the quad cell list*/
	if(quadcells != NULL)
	{
		finalize_quad_cells();
	}

	nfaces = (xdim-1)*(ydim-1);
	quadcells = (QuadCell**)malloc(sizeof(QuadCell*)*nfaces);
	for(i=0; i<nfaces; i++)
	{
		quadcells[i] = (QuadCell*)malloc(sizeof(QuadCell));
		quadcells[i]->verts = (int*)malloc(sizeof(int)*4);
	}

	int cell_index = 0;

	for(j=0; j<ydim-1; j++)
	{
		for(i=0; i<xdim-1; i++)
		{
			cell_index = j*(xdim-1)+i;
			quadcells[cell_index]->verts[0] = j*xdim+i;
			quadcells[cell_index]->verts[1] = j*xdim+(i+1);
			quadcells[cell_index]->verts[2] = (j+1)*xdim+(i+1);
			quadcells[cell_index]->verts[3] = (j+1)*xdim+i;
			quadcells[cell_index]->nverts=4;
			quadcells[cell_index]->index = cell_index;

			/**/
			quadcells[cell_index]->xstart = i;
			quadcells[cell_index]->xend = i+1;
			quadcells[cell_index]->ystart = j;
			quadcells[cell_index]->yend = j+1;
			quadcells[cell_index]->x_start_coord = quad_verts[quadcells[cell_index]->verts[0]]->x;
			quadcells[cell_index]->y_start_coord = quad_verts[quadcells[cell_index]->verts[0]]->y;

			quadcells[cell_index]->degpt_index = -1;

			quadcells[cell_index]->OnBoundary=false;
			quadcells[cell_index]->visited=false;
			//quadcells[cell_index]->repell_inregion=false;
			//quadcells[cell_index]->attract_inregion=false;

			/*for even tensor line placement*/
			//quadcells[cell_index]->samplepts=NULL;
			//quadcells[cell_index]->num_samplepts = 0;
			quadcells[cell_index]->maj_samplepts=NULL;
			quadcells[cell_index]->maj_nsamplepts = 0;
			quadcells[cell_index]->MAJMaxSampNum = 0;
			quadcells[cell_index]->min_samplepts=NULL;
			quadcells[cell_index]->min_nsamplepts = 0;
			quadcells[cell_index]->MINMaxSampNum = 0;

			/*for the computation of the intersections 10/02/2007*/
			quadcells[cell_index]->majorlines = NULL;
			quadcells[cell_index]->hasmajor = false;
			quadcells[cell_index]->minorlines = NULL;
			quadcells[cell_index]->hasminor = false;

			quadcells[cell_index]->sketchlines = NULL;

			quadcells[cell_index]->intersectlist = NULL;
			quadcells[cell_index]->nintersects = 0;

			quadcells[cell_index]->streetgraphedgelist = NULL;
			quadcells[cell_index]->nstreetgraphedges = 0;

			quadcells[cell_index]->which_region = 0;
			quadcells[cell_index]->in_region=false;

			quadcells[cell_index]->mapbounds=NULL;
			quadcells[cell_index]->nmapbounds=0;

			//quadcells[cell_index]->in_veg=false;

			/*add it to the corresponding vertices*/
			/*v0*/
			QuadVertex *v;
			v = quad_verts[quadcells[cell_index]->verts[0]];
			v->cells= extend_celllist_ver(v->cells, v->ncells);
			v->cells[v->ncells] = quadcells[cell_index];
			v->ncells ++;
			/*v1*/
			v = quad_verts[quadcells[cell_index]->verts[1]];
			v->cells= extend_celllist_ver(v->cells, v->ncells);
			v->cells[v->ncells] = quadcells[cell_index];
			v->ncells ++;
			/*v2*/
			v = quad_verts[quadcells[cell_index]->verts[2]];
			v->cells= extend_celllist_ver(v->cells, v->ncells);
			v->cells[v->ncells] = quadcells[cell_index];
			v->ncells ++;
			/*v3*/
			v = quad_verts[quadcells[cell_index]->verts[3]];
			v->cells= extend_celllist_ver(v->cells, v->ncells);
			v->cells[v->ncells] = quadcells[cell_index];
			v->ncells ++;
		}
	}
}

void QuadMesh::construct_edges()
{
	///First create and initialize the head knot of the edge link
	QuadEdge *Cur_elink;

	elist = new QuadEdge();
	elist->index = -1;
	elist->next = NULL;
	Cur_elink = elist;

	int edge_id = 0;
	nedges = 0;

	////////////////Define variables for vertices and faces operation
	QuadVertex **vert = quad_verts;

	///////////////////////////////
	int i, j, m, n;
	int Cur_vert, Next_vert;

	for( i = 0; i < nfaces; i++)
	{
		for( j = 0; j < quadcells[i]->nverts; j++)
		{
			//We need to check the neighbor vertex on that surface
			if( j == quadcells[i]->nverts - 1){
				Cur_vert = quadcells[i]->verts[j];
				Next_vert = quadcells[i]->verts[0];
			}
			else{
				Cur_vert = quadcells[i]->verts[j];          //extract the ID of the i'th face and j'th vertex
				Next_vert = quadcells[i]->verts[j+1];       //extract the ID of the i'th face and j+1'th vertex
			}

			//check if there is any edge between them or not
			if(vert[Cur_vert]->Num_edge == 0 || vert[Next_vert]->Num_edge == 0)
			//there must be no edge between them at this moment
			{
				///Create new notes for the edge link
				QuadEdge *new_edge = new QuadEdge;
				new_edge->index = edge_id;              //first edge will be marked 0, and so on...
				edge_id ++;

				/////Initialize the id of the adjacent faces that share the edge
				new_edge->tris[0] = -1;
				new_edge->tris[1] = -1;

				new_edge->tris[0] = i;                        //this is the first surface sharing the edge
				new_edge->verts[0] = Cur_vert;                //Save the ids of current vertices as the terminals of the edge
				new_edge->verts[1] = Next_vert;               //Using the current orientation!!!! 1/11
				new_edge->visited = false;                        //for my subdivision
				new_edge->next = NULL;
				Cur_elink->next = new_edge;                   //Add to the current edge link
				Cur_elink = new_edge;

				/*compute the edge length 07/21/07*/
				//new_edge->length = GetEdgeLength(Cur_vert, Next_vert);


				vert[Cur_vert]->edges = 
					Extend_Elist(vert[Cur_vert]->edges, vert[Cur_vert]->Num_edge);
				vert[Next_vert]->edges = 
					Extend_Elist(vert[Next_vert]->edges, vert[Next_vert]->Num_edge);

				vert[Cur_vert]->edges[vert[Cur_vert]->Num_edge] = new_edge;
				vert[Next_vert]->edges[vert[Next_vert]->Num_edge] = new_edge;

				vert[Cur_vert]->Num_edge++;
				vert[Next_vert]->Num_edge++;
						
				/////Add to the face
				quadcells[i]->edges[j] = new_edge;         //Link the new edge to the associated face

                ////The total number of edges add one
				nedges++;
		    }
			else{
				for( m = 0; m < vert[Cur_vert]->Num_edge; m++)
					for( n = 0; n < vert[Next_vert]->Num_edge; n++)
					{
						if( vert[Cur_vert]->edges[m]->index
							== vert[Next_vert]->edges[n]->index) 
						//There already has an edge between these two vertices
						{

							vert[Cur_vert]->edges[m]->tris[1] = i;

							quadcells[i]->edges[j] = vert[Cur_vert]->edges[m];

							goto LL;                          //if same edge ID has been found, jump out of the loop

						}
					}
				
LL:				if( m > vert[Cur_vert]->Num_edge - 1 )
				//Did not find an existing edge between these two vertices
				{
					///Create new notes for the edge link
					QuadEdge *new_edge = new QuadEdge;
					new_edge->index = edge_id;         //first edge will be marked 0, and so on...
					edge_id ++;
					/////Initialize the id of the adjacent faces that share the edge
					new_edge->tris[0] = -1;
					new_edge->tris[1] = -1;

					new_edge->tris[0] = i;                    //this is the first surface sharing the edge
					new_edge->verts[0] = Cur_vert;            //Save the ids of current vertices as the terminals of the edge
					new_edge->verts[1] = Next_vert;
					new_edge->visited = 0;                    //for my subdivision
					new_edge->next = NULL;
					Cur_elink->next = new_edge;               //Add to the current edge link
					Cur_elink = new_edge;
					
					/*compute the edge length 07/21/07*/
					//new_edge->length = GetEdgeLength(Cur_vert, Next_vert);

					///Add the ID of the edge into corresponding vertices and faces

					vert[Cur_vert]->edges = 
						Extend_Elist(vert[Cur_vert]->edges, vert[Cur_vert]->Num_edge);
					vert[Next_vert]->edges = 
						Extend_Elist(vert[Next_vert]->edges, vert[Next_vert]->Num_edge);

					vert[Cur_vert]->edges[vert[Cur_vert]->Num_edge] = new_edge;
					vert[Next_vert]->edges[vert[Next_vert]->Num_edge] = new_edge;

					vert[Cur_vert]->Num_edge++;
					vert[Next_vert]->Num_edge++;

					/////Add to the face

				    quadcells[i]->edges[j] = new_edge;     //Link the new edge to the associated face

					///The total number of edges add one
					nedges++;
				}
			}

		}
	}

	//Object.nedges = global_edge_id;
	//TotalEdgesNum = nedges;

	/*construct the array style edge list for later convenience*/
	edgelist=(QuadEdge**)malloc(sizeof(QuadEdge*)*nedges);
	QuadEdge *e = elist->next;
	for(i=0; i<nedges; i++)
	{
		edgelist[i]= e;
		e=e->next;
	}
}


QuadEdge  **QuadMesh::Extend_Elist(QuadEdge **edge_link, int Num_edges)
{
    QuadEdge **temp = edge_link;
	QuadEdge **new_edge_link = (QuadEdge **) malloc (sizeof (QuadEdge*)*(Num_edges+1)); //Extend the link
	if( Num_edges > 0)
	{
		for(int i = 0; i < Num_edges; i++)
			new_edge_link[i] = temp[i];
		free (temp);
	}
   
	return new_edge_link;
}

/*extend the cell list for vertex*/
QuadCell  **QuadMesh::extend_celllist_ver(QuadCell **cells, int ncells)
{
    QuadCell **temp = cells;
	QuadCell **new_cells = (QuadCell **) malloc (sizeof (QuadCell*)*(ncells+1)); //Extend the link
	if( ncells > 0)
	{
		for(int i = 0; i < ncells; i++)
			new_cells[i] = temp[i];
		free (temp);
	}
    cells = new_cells;
	return new_cells;
}


/*
orient the edge list in each cell in the following order
   e2
e3    e1
   e0
*/
void QuadMesh::orient_edges_cells()
{
	int i, j, k;
	QuadCell *face;
	QuadEdge *e[4];
	
	for(i=0; i<nfaces; i++)
	{
		face = quadcells[i];
		for(j=0; j<face->nverts; j++)
		{
			if((face->edges[j]->verts[0]==face->verts[0]&&face->edges[j]->verts[1]==face->verts[1])
				||(face->edges[j]->verts[0]==face->verts[1]&&face->edges[j]->verts[1]==face->verts[0]))
				e[0] = face->edges[j];
			else if((face->edges[j]->verts[0]==face->verts[1]&&face->edges[j]->verts[1]==face->verts[2])
				||(face->edges[j]->verts[0]==face->verts[2]&&face->edges[j]->verts[1]==face->verts[1]))
				e[1] = face->edges[j];
			else if((face->edges[j]->verts[0]==face->verts[2]&&face->edges[j]->verts[1]==face->verts[3])
				||(face->edges[j]->verts[0]==face->verts[3]&&face->edges[j]->verts[1]==face->verts[2]))
				e[2] = face->edges[j];
			else
				e[3] = face->edges[j];
		}

		for(j=0; j<face->nverts; j++)
			face->edges[j]=e[j];
	}
}
