/*
This file implements the tensor field creation
NOTE: the data structure has combined with current data structure
if the data structure is changed, corresponding changes are necessary.
09/19/2007
*/
#include "stdafx.h"

#include "VFDataStructure.h"

#include "tensorvis.h"

#include "tensoranalysis.h"

#include "tensordesign.h"

#include "sketchdesign.h"

Degenerate_Design *ten_designelems = NULL;
int ntenelems = 0;
int curMaxNumTenDesignElems = 100;

TenRegularElem *ten_regularelems = NULL;
int nten_regelems = 0;
int curMaxNumTenRegElems =100;

/*   A global variable to record the previous (or current) tensor field   
     The variable is useful when one to add new tensor on current(or previous) field
*/
icMatrix2x2 *pre_ten=NULL;
int ten_blend_factor = .5;  /*we use linear interpolation to blend two tensor fields*/

#include "SketchDesign.h"
extern int NDesignRegions;
extern int cur_chosen_region; /*for multi-region design*/


extern Polygon3D Object;
extern int chosen_tenelem_ID;
extern QuadMesh *quadmesh;

extern unsigned char cur_max_reg_index;

extern int get_cellID_givencoords(double x, double y);


//#include "shareinterfacevars.h"
//extern SharedInterfaceVars sharedvars;

bool please_comb_prefield;

icMatrix2x2 *pre_tenfield=NULL;

/*  for local editing, we keep the global field for a while */
//int cur_sel_region = 0;

/*    intialize the pre_ten variable to store the tensor field    */
void init_pre_ten()
{
	if(pre_ten!=NULL)
		delete [] pre_ten;

	pre_ten=new icMatrix2x2[quadmesh->nverts];
}

void store_cur_ten()
{
	int i;
	for(i=0;i<quadmesh->nverts;i++)
		/*pre_ten*/pre_tenfield[i].set(quadmesh->quad_verts[i]->Jacobian);
}


/*using new way to compute the eigen vectors of a given 2x2 symmetric tensor
09/24/2007*/
void cal_eigen_vector_sym(icMatrix2x2 m, icVector2 ev[2])
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

	//if(theta < 0) theta += 2*M_PI;

	/*major eigen vector*/
	ev[0].entry[0] = cos(theta/2.);
	ev[0].entry[1] = sin(theta/2.);

	//ev[0] = half_trace*ev[0];

	/*minor eigen vector*/
	ev[1].entry[0] = cos((theta+M_PI)/2.);
	ev[1].entry[1] = sin((theta+M_PI)/2.);

	//ev[1] = half_trace*ev[1];
}


/*
calculate the eigen vector systems of the symmetric tensor fields
over the mesh  09/24/2007
*/

extern void normalized_tensorfield();

void cal_sym_eigenvecs_vertices()
{
	int i, j;
	Vertex *v;
	icVector2 ev[2];
	for(i=0; i<Object.nverts; i++)
	{
		v = Object.vlist[i];
		cal_eigen_vector_sym(v->Jacobian, ev);
		v->major = ev[0];
		v->minor = ev[1];
	}

	/*normalize the major and minor field*/
	normalized_tensorfield();
}


/*
intialize the design element lists
*/
void init_ten_designelems()
{
	/*initialize singular design element*/
	if(ten_designelems == NULL)
	{
		ten_designelems=(Degenerate_Design *)malloc(sizeof(Degenerate_Design)*curMaxNumTenDesignElems);
		if(ten_designelems == NULL)
			exit(-1);
	}
	ntenelems = 0;

	/*initialize regular design element*/
	//if(ten_regularelems == NULL)
	//{
	//	ten_regularelems=(TenRegularElem *)malloc(sizeof(TenRegularElem)*curMaxNumTenRegElems);
	//	if(ten_regularelems == NULL)
	//		exit(-1);
	//}
	//nten_regelems = 0;
	init_regular_ten_designelems();
}

void init_regular_ten_designelems()
{
	/*initialize regular design element*/
	if(ten_regularelems == NULL)
	{
		ten_regularelems=(TenRegularElem *)malloc(sizeof(TenRegularElem)*curMaxNumTenRegElems);
		if(ten_regularelems == NULL)
			exit(-1);
	}
	nten_regelems = 0;
}


/*Add an elememt into the list
*/
void addto(double x, double y, int triangle, int type)
{
	/*extend space if not enough*/
	if(ntenelems >= curMaxNumTenDesignElems)
	{
		ten_designelems = (Degenerate_Design*)realloc(ten_designelems, sizeof(Degenerate_Design)*
			(curMaxNumTenDesignElems+20));
		if(ten_designelems == NULL)
			exit(-1);
		curMaxNumTenDesignElems += 20;
	}

	ten_designelems[ntenelems].centerx = x;
	ten_designelems[ntenelems].centery = y;
	ten_designelems[ntenelems].ID = ntenelems;
	ten_designelems[ntenelems].Triangle_ID = triangle;
	ten_designelems[ntenelems].type = type;
	ten_designelems[ntenelems].transform_matrix.setIdentity();
	
	ten_designelems[ntenelems].rotang = 0;
	ten_designelems[ntenelems].s = 
	ten_designelems[ntenelems].sx =
	ten_designelems[ntenelems].sy = 1;
	ten_designelems[ntenelems].deleted = false;

	/*mark which region the design element belongs to 11/21/2007*/
	ten_designelems[ntenelems].which_region=get_region_id(x, y);

	/*  may be a bug 1/18/2008 */
	cur_chosen_region = ten_designelems[ntenelems].which_region;
			
	//Add the editbox for properties editing
	init_tenelem_EditBox(ntenelems, x, y);

	ntenelems++;
}

void init_tenelem_EditBox(int index, double x, double y)
{
	////initial the original edit box
	ten_designelems[index].editbox.p1.entry[0] = x - EDITBOXSIZE;  ////low left point
	ten_designelems[index].editbox.p1.entry[1] = y - EDITBOXSIZE;

	ten_designelems[index].editbox.p2.entry[0] = x - EDITBOXSIZE;  ////upper left point
	ten_designelems[index].editbox.p2.entry[1] = y + EDITBOXSIZE;

	ten_designelems[index].editbox.p3.entry[0] = x + EDITBOXSIZE;  ////upper right point
	ten_designelems[index].editbox.p3.entry[1] = y + EDITBOXSIZE;

	ten_designelems[index].editbox.p4.entry[0] = x + EDITBOXSIZE;  ////low right point
	ten_designelems[index].editbox.p4.entry[1] = y - EDITBOXSIZE;

	ten_designelems[index].editbox.Up.entry[0] = x;  ////rotation controling point
	ten_designelems[index].editbox.Up.entry[1] = y + 2*EDITBOXSIZE;


	////initial current editbox
	ten_designelems[index].cur_editbox.p1.entry[0] = x - EDITBOXSIZE;  ////low left point
	ten_designelems[index].cur_editbox.p1.entry[1] = y - EDITBOXSIZE;

	ten_designelems[index].cur_editbox.p2.entry[0] = x - EDITBOXSIZE;  ////upper left point
	ten_designelems[index].cur_editbox.p2.entry[1] = y + EDITBOXSIZE;

	ten_designelems[index].cur_editbox.p3.entry[0] = x + EDITBOXSIZE;  ////upper right point
	ten_designelems[index].cur_editbox.p3.entry[1] = y + EDITBOXSIZE;

	ten_designelems[index].cur_editbox.p4.entry[0] = x + EDITBOXSIZE;  ////low right point
	ten_designelems[index].cur_editbox.p4.entry[1] = y - EDITBOXSIZE;
	
	ten_designelems[index].cur_editbox.Up.entry[0] = x;  ////rotation controling point
	ten_designelems[index].cur_editbox.Up.entry[1] = y + 2*EDITBOXSIZE;
}

#define LOWER 1

/*
sum up the basis field
*/
//void get_tensor(double x, double y, double t[4])
//{
//   int i;
//   double  dx, dy, vx, vy, t00, t01, t10, t11, r=0.;
//   double d;
//
//   double tt00, tt01, tt10, tt11;
//
//   icMatrix3x3 tempJacobian, transposerot, transpose_tran;
//   double ang;
//   //double newdirect[2] = {0, 0};
//
//   //double sim_magnitude = 0;  ////
//
//   vx = vy = 0.;
//
//   t[0]=t[1]=t[2]=t[3]=0.;
//
//   ///Combine all the degenerate elements 
//   for(i = 0; i < ntenelems; i++)
//   {
//	   if(ten_designelems[i].ID >= 0 && !ten_designelems[i].deleted){
//			dx = x - ten_designelems[i].centerx;
//			dy = y - ten_designelems[i].centery;
//
//			r  = dx*dx + dy*dy; 
//			d = exp(-1000*r);
//
//    		if (r < DistanceThreshold)   r = DistanceThreshold;
//
//			tempJacobian.set(ten_designelems[i].transform_matrix);
//
//			ang = ten_designelems[i].rotang;
//
//			//transposerot.set(cos(ang), sin(ang), 0,
//			//                 -sin(ang),  cos(ang), 0,
//			//				 0,0,1);
//			transposerot.set(cos(ang/2), sin(ang/2), 0,
//			                 -sin(ang/2),  cos(ang/2), 0,
//							 0,0,1);
//			tempJacobian.rightMultiply(transposerot);
//
//			transpose_tran.set(tempJacobian.entry[0][0], tempJacobian.entry[1][0], tempJacobian.entry[2][0],
//				tempJacobian.entry[0][1], tempJacobian.entry[1][1], tempJacobian.entry[2][1],
//				tempJacobian.entry[0][2], tempJacobian.entry[1][2], tempJacobian.entry[2][2]);
//
//
//			if(ten_designelems[i].type == 0) /*wedge*/
//			{
//				/*one order*/
//			
//				/*second order*/
//				tt00 = (dx * transpose_tran.entry[0][0] + dy * transpose_tran.entry[0][1]);  
//				tt01 = (dy * transpose_tran.entry[0][0] - dx * transpose_tran.entry[0][1]);
//				tt10 = (dx * transpose_tran.entry[1][0] + dy * transpose_tran.entry[1][1]);
//				tt11 = (dy * transpose_tran.entry[1][0] - dx * transpose_tran.entry[1][1]);
//
//			}
//			else if(ten_designelems[i].type == 1) /*trisector*/
//			{
//				tt00 = (dx * transpose_tran.entry[0][0] - dy * transpose_tran.entry[0][1]);  
//				tt01 = (-dy * transpose_tran.entry[0][0] - dx * transpose_tran.entry[0][1]);
//				tt10 = (dx * transpose_tran.entry[1][0] - dy * transpose_tran.entry[1][1]);
//				tt11 = (-dy * transpose_tran.entry[1][0] - dx * transpose_tran.entry[1][1]);
//				
//			}
//			else if(ten_designelems[i].type == 2) /*node*/
//			{
//				tt00 = ((dx*dx-dy*dy)*transpose_tran.entry[0][0]+2*dx*dy*transpose_tran.entry[0][1])/(sqrt(r));
//				tt01 = (2*dx*dy*transpose_tran.entry[0][0]-(dx*dx-dy*dy)*transpose_tran.entry[0][1])/(sqrt(r));
//				tt10 = ((dx*dx-dy*dy)*transpose_tran.entry[1][0]+2*dx*dy*transpose_tran.entry[1][1])/(sqrt(r));
//				tt11 = (2*dx*dy*transpose_tran.entry[1][0]-(dx*dx-dy*dy)*transpose_tran.entry[1][1])/(sqrt(r));
//				
//			}
//			else if(ten_designelems[i].type == 3) /*center*/
//			{
//				tt00 = ((dy*dy-dx*dx)*transpose_tran.entry[0][0]-2*dx*dy*transpose_tran.entry[0][1])/(sqrt(r));
//				tt01 = (-2*dx*dy*transpose_tran.entry[0][0]+(dx*dx-dy*dy)*transpose_tran.entry[0][1])/(sqrt(r));
//				tt10 = ((dy*dy-dx*dx)*transpose_tran.entry[1][0]-2*dx*dy*transpose_tran.entry[1][1])/(sqrt(r));
//				tt11 = (-2*dx*dy*transpose_tran.entry[1][0]+(dx*dx-dy*dy)*transpose_tran.entry[1][1])/(sqrt(r));
//			}
//			else if(ten_designelems[i].type == 4) /*saddle*/
//			{
//				tt00 = ((dx*dx-dy*dy)*transpose_tran.entry[0][0]-2*dx*dy*transpose_tran.entry[0][1])/(sqrt(r));
//				tt01 = (-2*dx*dy*transpose_tran.entry[0][0]-(dx*dx-dy*dy)*transpose_tran.entry[0][1])/(sqrt(r));
//				tt10 = ((dx*dx-dy*dy)*transpose_tran.entry[1][0]-2*dx*dy*transpose_tran.entry[1][1])/(sqrt(r));
//				tt11 = (-2*dx*dy*transpose_tran.entry[1][0]-(dx*dx-dy*dy)*transpose_tran.entry[1][1])/(sqrt(r));
//			}
//
//			/*you may also need to multiply the transpose of the transformation matrix*/
//				t00 = (tt00 * tempJacobian.entry[0][0] + tt01 * tempJacobian.entry[1][0])/r;  
//				t01 = (tt00 * tempJacobian.entry[0][1] + tt01 * tempJacobian.entry[1][1])/r;
//				t10 = (tt10 * tempJacobian.entry[0][0] + tt11 * tempJacobian.entry[1][0])/r;
//				t11 = (tt10 * tempJacobian.entry[0][1] + tt11 * tempJacobian.entry[1][1])/r;
//
//
//			t[0] += t00/sqrt(r)*LOWER;
//			t[1] += t01/sqrt(r)*LOWER;
//			t[2] += t10/sqrt(r)*LOWER;
//			t[3] += t11/sqrt(r)*LOWER;
//
//	   }
//	
//   }
//
//   /*the following we combine the regular element*/
//   icMatrix3x3 regten, temp;
//   for(i = 0; i < nten_regelems; i++)
//   {
//	   if(ten_regularelems[i].ID > 0 && !ten_regularelems[i].deleted)
//	   {
//			dx = x - ten_regularelems[i].base[0];
//			dy = y - ten_regularelems[i].base[1];
//
//			r  = dx*dx + dy*dy; 
//
//			if (r < DistanceThreshold)   r = DistanceThreshold;
//
//			if(ten_regularelems[i].type == 0) ////regular element
//			{
//				double strength = length(ten_regularelems[i].Direct);
//
//				t[0] += strength*cos(2*ten_regularelems[i].rotang)/(r*sqrt(r));
//				t[1] += strength*sin(2*ten_regularelems[i].rotang)/(r*sqrt(r));
//				t[2] += strength*sin(2*ten_regularelems[i].rotang)/(r*sqrt(r));
//				t[3] += -strength*cos(2*ten_regularelems[i].rotang)/(r*sqrt(r));
//			}
//	   }
//   }
//}
//
//
//
//
//
//
//
void get_tensor(double x, double y, double t[4])
{
   int i;
   double  dx, dy, vx, vy, t00, t01, t10, t11, r=0.;
   double d;

   icMatrix3x3 tempJacobian, transposerot;
   double ang;
   //double newdirect[2] = {0, 0};

   //double sim_magnitude = 0;  ////

   vx = vy = 0.;

   t[0]=t[1]=t[2]=t[3]=0.;

   ///Combine all the degenerate elements 
   for(i = 0; i < ntenelems; i++)
   {
	   if(ten_designelems[i].ID >= 0 && !ten_designelems[i].deleted){
			dx = x - ten_designelems[i].centerx;
			dy = y - ten_designelems[i].centery;

			r  = dx*dx + dy*dy; 
			d = exp(-1000*r);

    		if (r < DistanceThreshold)   r = DistanceThreshold;

			tempJacobian.set(ten_designelems[i].transform_matrix);

			ang = ten_designelems[i].rotang;

			transposerot.set(cos(ang), sin(ang), 0,
			                 -sin(ang),  cos(ang), 0,
							 0,0,1);

			tempJacobian.rightMultiply(transposerot);

			if(ten_designelems[i].type == 0) /*wedge*/
			{
				t00 = (dx * tempJacobian.entry[0][0] + dy * tempJacobian.entry[0][1])/r;  
				t01 = (dy * tempJacobian.entry[0][0] - dx * tempJacobian.entry[0][1])/r;
				t10 = (dx * tempJacobian.entry[1][0] + dy * tempJacobian.entry[1][1])/r;
				t11 = (dy * tempJacobian.entry[1][0] - dx * tempJacobian.entry[1][1])/r;
			}
			else if(ten_designelems[i].type == 1) /*trisector*/
			{
				t00 = (dx * tempJacobian.entry[0][0] - dy * tempJacobian.entry[0][1])/r;  
				t01 = (-dy * tempJacobian.entry[0][0] - dx * tempJacobian.entry[0][1])/r;
				t10 = (dx * tempJacobian.entry[1][0] - dy * tempJacobian.entry[1][1])/r;
				t11 = (-dy * tempJacobian.entry[1][0] - dx * tempJacobian.entry[1][1])/r;
			}
			else if(ten_designelems[i].type == 2) /*node*/
			{
				t00 = ((dx*dx-dy*dy)*tempJacobian.entry[0][0]+2*dx*dy*tempJacobian.entry[0][1])/(r*sqrt(r));
				t01 = (2*dx*dy*tempJacobian.entry[0][0]-(dx*dx-dy*dy)*tempJacobian.entry[0][1])/(r*sqrt(r));
				t10 = ((dx*dx-dy*dy)*tempJacobian.entry[1][0]+2*dx*dy*tempJacobian.entry[1][1])/(r*sqrt(r));
				t11 = (2*dx*dy*tempJacobian.entry[1][0]-(dx*dx-dy*dy)*tempJacobian.entry[1][1])/(r*sqrt(r));
				
			}
			else if(ten_designelems[i].type == 3) /*center*/
			{
				t00 = ((dy*dy-dx*dx)*tempJacobian.entry[0][0]-2*dx*dy*tempJacobian.entry[0][1])/(r*sqrt(r));
				t01 = (-2*dx*dy*tempJacobian.entry[0][0]+(dx*dx-dy*dy)*tempJacobian.entry[0][1])/(r*sqrt(r));
				t10 = ((dy*dy-dx*dx)*tempJacobian.entry[1][0]-2*dx*dy*tempJacobian.entry[1][1])/(r*sqrt(r));
				t11 = (-2*dx*dy*tempJacobian.entry[1][0]+(dx*dx-dy*dy)*tempJacobian.entry[1][1])/(r*sqrt(r));
			}
			else if(ten_designelems[i].type == 4) /*saddle*/
			{
				t00 = ((dx*dx-dy*dy)*tempJacobian.entry[0][0]-2*dx*dy*tempJacobian.entry[0][1])/(r*sqrt(r));
				t01 = (-2*dx*dy*tempJacobian.entry[0][0]-(dx*dx-dy*dy)*tempJacobian.entry[0][1])/(r*sqrt(r));
				t10 = ((dx*dx-dy*dy)*tempJacobian.entry[1][0]-2*dx*dy*tempJacobian.entry[1][1])/(r*sqrt(r));
				t11 = (-2*dx*dy*tempJacobian.entry[1][0]-(dx*dx-dy*dy)*tempJacobian.entry[1][1])/(r*sqrt(r));
			}

			/*you may also need to multiply the transpose of the transformation matrix*/
			

			t[0] += t00/sqrt(r)*LOWER;
			t[1] += t01/sqrt(r)*LOWER;
			t[2] += t10/sqrt(r)*LOWER;
			t[3] += t11/sqrt(r)*LOWER;

	   }
	
   }

   /*the following we combine the regular element*/
   icMatrix3x3 regten, temp;
   for(i = 0; i < nten_regelems; i++)
   {
	   if(ten_regularelems[i].ID > 0 && !ten_regularelems[i].deleted)
	   {
			dx = x - ten_regularelems[i].base[0];
			dy = y - ten_regularelems[i].base[1];

			r  = dx*dx + dy*dy; 

			if (r < DistanceThreshold)   r = DistanceThreshold;

			if(ten_regularelems[i].type == 0) ////regular element
			{
				/*the creation of regular element before 09/26/2007*/
				double strength = length(ten_regularelems[i].Direct);

				t[0] += strength*cos(2*ten_regularelems[i].rotang)/(r*sqrt(r));
				t[1] += strength*sin(2*ten_regularelems[i].rotang)/(r*sqrt(r));
				t[2] += strength*sin(2*ten_regularelems[i].rotang)/(r*sqrt(r));
				t[3] += -strength*cos(2*ten_regularelems[i].rotang)/(r*sqrt(r));

			}

	   }
   }

}

/*
calculate the tensor based on the design elements
*/
void cal_alltensor()
{
	int i;
	Vertex *v;
	double t[4]={0.};

	for(i=0; i<Object.nverts; i++)
	{
		v = Object.vlist[i];
		get_tensor(v->x, v->y, t);
		v->Jacobian.entry[0][0] = t[0];
		v->Jacobian.entry[0][1] = t[1];
		v->Jacobian.entry[1][0] = t[2];
		v->Jacobian.entry[1][1] = t[3];
	}
}

/*Quad mesh
calculate the tensor based on the design elements
*/
void cal_alltensor_quad()
{
	int i;
	QuadVertex *v;
	double t[4]={0.};

	for(i=0; i<quadmesh->nverts; i++)
	{
		v = quadmesh->quad_verts[i];

		if(!v->inland)
			continue;

		get_tensor(v->x, v->y, t);
		v->Jacobian.entry[0][0] = t[0];
		v->Jacobian.entry[0][1] = t[1];
		v->Jacobian.entry[1][0] = t[2];
		v->Jacobian.entry[1][1] = t[3];
	}
}


/***********************************************************************/
/*
Region based tensor field design
*/
/*
sum up the basis field
*/
//void get_tensor_inReg(double x, double y, double t[4], int regionid)
//{
//   int i;
//   double  dx, dy, vx, vy, t00, t01, t10, t11, r=0.;
//   double d;
//
//   double tt00, tt01, tt10, tt11;
//
//   icMatrix3x3 tempJacobian, transposerot, transpose_tran;
//   double ang;
//
//   vx = vy = 0.;
//
//   t[0]=t[1]=t[2]=t[3]=0.;
//
//   ///Combine all the degenerate elements 
//   for(i = 0; i < ntenelems; i++)
//   {
//	   if(ten_designelems[i].ID >= 0 && !ten_designelems[i].deleted
//		   && ten_designelems[i].which_region==regionid)
//	   {
//			dx = x - ten_designelems[i].centerx;
//			dy = y - ten_designelems[i].centery;
//
//			r  = dx*dx + dy*dy; 
//			d = exp(-1000*r);
//
//    		if (r < DistanceThreshold)   r = DistanceThreshold;
//
//			tempJacobian.set(ten_designelems[i].transform_matrix);
//
//			ang = ten_designelems[i].rotang;
//
//			transposerot.set(cos(ang/2), sin(ang/2), 0,
//			                 -sin(ang/2),  cos(ang/2), 0,
//							 0,0,1);
//			tempJacobian.rightMultiply(transposerot);
//
//			transpose_tran.set(tempJacobian.entry[0][0], tempJacobian.entry[1][0], tempJacobian.entry[2][0],
//				tempJacobian.entry[0][1], tempJacobian.entry[1][1], tempJacobian.entry[2][1],
//				tempJacobian.entry[0][2], tempJacobian.entry[1][2], tempJacobian.entry[2][2]);
//
//
//			if(ten_designelems[i].type == 0) /*wedge*/
//			{
//				tt00 = (dx * transpose_tran.entry[0][0] + dy * transpose_tran.entry[0][1]);  
//				tt01 = (dy * transpose_tran.entry[0][0] - dx * transpose_tran.entry[0][1]);
//				tt10 = (dx * transpose_tran.entry[1][0] + dy * transpose_tran.entry[1][1]);
//				tt11 = (dy * transpose_tran.entry[1][0] - dx * transpose_tran.entry[1][1]);
//			}
//			else if(ten_designelems[i].type == 1) /*trisector*/
//			{
//				tt00 = (dx * transpose_tran.entry[0][0] - dy * transpose_tran.entry[0][1]);  
//				tt01 = (-dy * transpose_tran.entry[0][0] - dx * transpose_tran.entry[0][1]);
//				tt10 = (dx * transpose_tran.entry[1][0] - dy * transpose_tran.entry[1][1]);
//				tt11 = (-dy * transpose_tran.entry[1][0] - dx * transpose_tran.entry[1][1]);
//				
//			}
//			else if(ten_designelems[i].type == 2) /*node*/
//			{
//				tt00 = ((dx*dx-dy*dy)*transpose_tran.entry[0][0]+2*dx*dy*transpose_tran.entry[0][1])/(sqrt(r));
//				tt01 = (2*dx*dy*transpose_tran.entry[0][0]-(dx*dx-dy*dy)*transpose_tran.entry[0][1])/(sqrt(r));
//				tt10 = ((dx*dx-dy*dy)*transpose_tran.entry[1][0]+2*dx*dy*transpose_tran.entry[1][1])/(sqrt(r));
//				tt11 = (2*dx*dy*transpose_tran.entry[1][0]-(dx*dx-dy*dy)*transpose_tran.entry[1][1])/(sqrt(r));
//				
//			}
//			else if(ten_designelems[i].type == 3) /*center*/
//			{
//			
//				tt00 = ((dy*dy-dx*dx)*transpose_tran.entry[0][0]-2*dx*dy*transpose_tran.entry[0][1])/(sqrt(r));
//				tt01 = (-2*dx*dy*transpose_tran.entry[0][0]+(dx*dx-dy*dy)*transpose_tran.entry[0][1])/(sqrt(r));
//				tt10 = ((dy*dy-dx*dx)*transpose_tran.entry[1][0]-2*dx*dy*transpose_tran.entry[1][1])/(sqrt(r));
//				tt11 = (-2*dx*dy*transpose_tran.entry[1][0]+(dx*dx-dy*dy)*transpose_tran.entry[1][1])/(sqrt(r));
//			}
//			else if(ten_designelems[i].type == 4) /*saddle*/
//			{
//
//				tt00 = ((dx*dx-dy*dy)*transpose_tran.entry[0][0]-2*dx*dy*transpose_tran.entry[0][1])/(sqrt(r));
//				tt01 = (-2*dx*dy*transpose_tran.entry[0][0]-(dx*dx-dy*dy)*transpose_tran.entry[0][1])/(sqrt(r));
//				tt10 = ((dx*dx-dy*dy)*transpose_tran.entry[1][0]-2*dx*dy*transpose_tran.entry[1][1])/(sqrt(r));
//				tt11 = (-2*dx*dy*transpose_tran.entry[1][0]-(dx*dx-dy*dy)*transpose_tran.entry[1][1])/(sqrt(r));
//			}
//
//			/*you may also need to multiply the transpose of the transformation matrix*/
//			t00 = (tt00 * tempJacobian.entry[0][0] + tt01 * tempJacobian.entry[1][0])/r;  
//			t01 = (tt00 * tempJacobian.entry[0][1] + tt01 * tempJacobian.entry[1][1])/r;
//			t10 = (tt10 * tempJacobian.entry[0][0] + tt11 * tempJacobian.entry[1][0])/r;
//			t11 = (tt10 * tempJacobian.entry[0][1] + tt11 * tempJacobian.entry[1][1])/r;
//
//
//			t[0] += t00/sqrt(r)*LOWER;
//			t[1] += t01/sqrt(r)*LOWER;
//			t[2] += t10/sqrt(r)*LOWER;
//			t[3] += t11/sqrt(r)*LOWER;
//
//	   }
//	
//   }
//
//   /*the following we combine the regular element*/
//   icMatrix3x3 regten, temp;
//   for(i = 0; i < nten_regelems; i++)
//   {
//	   if(ten_regularelems[i].ID > 0 && !ten_regularelems[i].deleted
//		   && ten_regularelems[i].which_region == regionid)
//	   {
//			dx = x - ten_regularelems[i].base[0];
//			dy = y - ten_regularelems[i].base[1];
//
//			r  = dx*dx + dy*dy; 
//
//			if (r < DistanceThreshold)   r = DistanceThreshold;
//
//			if(ten_regularelems[i].type == 0) ////regular element
//			{
//				double strength = length(ten_regularelems[i].Direct);
//
//				t[0] += strength*cos(2*ten_regularelems[i].rotang)/(r*sqrt(r));
//				t[1] += strength*sin(2*ten_regularelems[i].rotang)/(r*sqrt(r));
//				t[2] += strength*sin(2*ten_regularelems[i].rotang)/(r*sqrt(r));
//				t[3] += -strength*cos(2*ten_regularelems[i].rotang)/(r*sqrt(r));
//
//			}
//
//	   }
//   }
//
//}



void get_tensor_inReg(double x, double y, double t[4], int regionid, bool inbrush)
{
   int i;
   double  dx, dy, vx, vy, t00, t01, t10, t11, r=0.;
   double d;

   icMatrix3x3 tempJacobian, transposerot;
   double ang;

   vx = vy = 0.;

   t[0]=t[1]=t[2]=t[3]=0.;

   ///Combine all the degenerate elements 
   for(i = 0; i < ntenelems; i++)
   {
	   if(ten_designelems[i].ID >= 0 && !ten_designelems[i].deleted
		   && ten_designelems[i].which_region==regionid)
	   {
			dx = x - ten_designelems[i].centerx;
			dy = y - ten_designelems[i].centery;

			r  = dx*dx + dy*dy; 
			d = exp(-1000*r);

    		if (r < DistanceThreshold)   r = DistanceThreshold;

			tempJacobian.set(ten_designelems[i].transform_matrix);

			ang = ten_designelems[i].rotang;

			transposerot.set(cos(ang), sin(ang), 0,
			                 -sin(ang),  cos(ang), 0,
							 0,0,1);

			tempJacobian.rightMultiply(transposerot);

			if(ten_designelems[i].type == 0) /*wedge*/
			{
				t00 = (dx * tempJacobian.entry[0][0] + dy * tempJacobian.entry[0][1])/r;  
				t01 = (dy * tempJacobian.entry[0][0] - dx * tempJacobian.entry[0][1])/r;
				t10 = (dx * tempJacobian.entry[1][0] + dy * tempJacobian.entry[1][1])/r;
				t11 = (dy * tempJacobian.entry[1][0] - dx * tempJacobian.entry[1][1])/r;
			}
			else if(ten_designelems[i].type == 1) /*trisector*/
			{
				t00 = (dx * tempJacobian.entry[0][0] - dy * tempJacobian.entry[0][1])/r;  
				t01 = (-dy * tempJacobian.entry[0][0] - dx * tempJacobian.entry[0][1])/r;
				t10 = (dx * tempJacobian.entry[1][0] - dy * tempJacobian.entry[1][1])/r;
				t11 = (-dy * tempJacobian.entry[1][0] - dx * tempJacobian.entry[1][1])/r;
			}
			else if(ten_designelems[i].type == 2) /*node*/
			{
				t00 = ((dx*dx-dy*dy)*tempJacobian.entry[0][0]+2*dx*dy*tempJacobian.entry[0][1])/(r*sqrt(r));
				t01 = (2*dx*dy*tempJacobian.entry[0][0]-(dx*dx-dy*dy)*tempJacobian.entry[0][1])/(r*sqrt(r));
				t10 = ((dx*dx-dy*dy)*tempJacobian.entry[1][0]+2*dx*dy*tempJacobian.entry[1][1])/(r*sqrt(r));
				t11 = (2*dx*dy*tempJacobian.entry[1][0]-(dx*dx-dy*dy)*tempJacobian.entry[1][1])/(r*sqrt(r));
				
			}
			else if(ten_designelems[i].type == 3) /*center*/
			{
				t00 = ((dy*dy-dx*dx)*tempJacobian.entry[0][0]-2*dx*dy*tempJacobian.entry[0][1])/(r*sqrt(r));
				t01 = (-2*dx*dy*tempJacobian.entry[0][0]+(dx*dx-dy*dy)*tempJacobian.entry[0][1])/(r*sqrt(r));
				t10 = ((dy*dy-dx*dx)*tempJacobian.entry[1][0]-2*dx*dy*tempJacobian.entry[1][1])/(r*sqrt(r));
				t11 = (-2*dx*dy*tempJacobian.entry[1][0]+(dx*dx-dy*dy)*tempJacobian.entry[1][1])/(r*sqrt(r));
			}
			else if(ten_designelems[i].type == 4) /*saddle*/
			{
				t00 = ((dx*dx-dy*dy)*tempJacobian.entry[0][0]-2*dx*dy*tempJacobian.entry[0][1])/(r*sqrt(r));
				t01 = (-2*dx*dy*tempJacobian.entry[0][0]-(dx*dx-dy*dy)*tempJacobian.entry[0][1])/(r*sqrt(r));
				t10 = ((dx*dx-dy*dy)*tempJacobian.entry[1][0]-2*dx*dy*tempJacobian.entry[1][1])/(r*sqrt(r));
				t11 = (-2*dx*dy*tempJacobian.entry[1][0]-(dx*dx-dy*dy)*tempJacobian.entry[1][1])/(r*sqrt(r));
			}

			/*you may also need to multiply the transpose of the transformation matrix*/
			

			t[0] += t00/sqrt(r)*LOWER;
			t[1] += t01/sqrt(r)*LOWER;
			t[2] += t10/sqrt(r)*LOWER;
			t[3] += t11/sqrt(r)*LOWER;

	   }
	
   }

   /*the following we combine the regular element*/
   icMatrix3x3 regten, temp;
   for(i = 0; i < nten_regelems; i++)
   {
	   if(ten_regularelems[i].ID > 0 && !ten_regularelems[i].deleted
		    && ten_regularelems[i].which_region == regionid)
	   {
			dx = x - ten_regularelems[i].base[0];
			dy = y - ten_regularelems[i].base[1];

			r  = dx*dx + dy*dy; 

			if (r < DistanceThreshold)   r = DistanceThreshold;

			if(ten_regularelems[i].type == 0) ////regular element
			{
				/*the creation of regular element before 09/26/2007*/
				double strength = length(ten_regularelems[i].Direct);

				t[0] += strength*cos(2*ten_regularelems[i].rotang)/(r*sqrt(r));
				t[1] += strength*sin(2*ten_regularelems[i].rotang)/(r*sqrt(r));
				t[2] += strength*sin(2*ten_regularelems[i].rotang)/(r*sqrt(r));
				t[3] += -strength*cos(2*ten_regularelems[i].rotang)/(r*sqrt(r));

			}

	   }

	   else if(/*sharedvars.GenTenFromSketchesOn &&*/ ten_regularelems[i].ID > 0 && 
		   !ten_regularelems[i].deleted && ten_regularelems[i].which_region == 0
		   && inbrush)
	   {
			dx = x - ten_regularelems[i].base[0];
			dy = y - ten_regularelems[i].base[1];

			r  = dx*dx + dy*dy; 

			if (r < DistanceThreshold)   r = DistanceThreshold;

			if(ten_regularelems[i].type == 0) ////regular element
			{
				/*the creation of regular element before 09/26/2007*/
				double strength = length(ten_regularelems[i].Direct);

				t[0] += strength*cos(2*ten_regularelems[i].rotang)/(r*sqrt(r));
				t[1] += strength*sin(2*ten_regularelems[i].rotang)/(r*sqrt(r));
				t[2] += strength*sin(2*ten_regularelems[i].rotang)/(r*sqrt(r));
				t[3] += -strength*cos(2*ten_regularelems[i].rotang)/(r*sqrt(r));

			}
	   }
   }

}


/*Quad mesh
calculate the tensor based on the design elements and user specified region
*/
void cal_tensorvals_quad_inReg(/*int regionid*/)
{
	int i;
	QuadVertex *v;
	double t[4]={0.};

	for(i=0; i<quadmesh->nverts; i++)
	{
		v = quadmesh->quad_verts[i];

		if(!v->inland)
			continue;

		get_tensor_inReg(v->x, v->y, t, v->which_region, v->inbrushregion);
		if(please_comb_prefield/*sharedvars.SelStreetRegToEditOn*/)
		{
			if(v->which_region == cur_max_reg_index)
			{
				v->Jacobian.entry[0][0] = t[0];
				v->Jacobian.entry[0][1] = t[1];
				v->Jacobian.entry[1][0] = t[2];
				v->Jacobian.entry[1][1] = t[3];
			}
			else{
				v->Jacobian.entry[0][0] = 0.1*pre_tenfield[i].entry[0][0]+t[0];
				v->Jacobian.entry[0][1] = 0.1*pre_tenfield[i].entry[0][1]+t[1];
				v->Jacobian.entry[1][0] = 0.1*pre_tenfield[i].entry[1][0]+t[2];
				v->Jacobian.entry[1][1] = 0.1*pre_tenfield[i].entry[1][1]+t[3];
			}
		}
		else{
			//get_tensor_inReg(v->x, v->y, t, v->which_region, v->inbrushregion);
			v->Jacobian.entry[0][0] = t[0];
			v->Jacobian.entry[0][1] = t[1];
			v->Jacobian.entry[1][0] = t[2];
			v->Jacobian.entry[1][1] = t[3];
		}

		//if(/*combine with previous 12/2/2007*/)
		//{
		//	v->Jacobian.entry[0][0]=(1-ten_blend_factor)*t[0]+pre_ten[i].entry[0][0];
		//	v->Jacobian.entry[0][1]=(1-ten_blend_factor)*t[1]+pre_ten[i].entry[0][1];
		//	v->Jacobian.entry[1][0]=(1-ten_blend_factor)*t[2]+pre_ten[i].entry[1][0];
		//	v->Jacobian.entry[1][1]=(1-ten_blend_factor)*t[3]+pre_ten[i].entry[1][1];
		//}
	}
}


bool exist_Elems_inReg(int regionid)
{
	int i;
   ///Combine all the degenerate elements 
   for(i = 0; i < ntenelems; i++)
   {
	   if(ten_designelems[i].ID >= 0 && !ten_designelems[i].deleted
		   && ten_designelems[i].which_region==regionid)
	   {
		   return true;
	   }
	
   }

   for(i = 0; i < nten_regelems; i++)
   {
	   if(ten_regularelems[i].ID > 0 && !ten_regularelems[i].deleted
		    && ten_regularelems[i].which_region == regionid)
	   {
		   return true;
	   }

	   //else if(/*sharedvars.GenTenFromSketchesOn &&*/ ten_regularelems[i].ID > 0 && 
		  // !ten_regularelems[i].deleted && ten_regularelems[i].which_region == 0
		  // && inbrush)
	   //{
		  // return true;
	   //}
   }

   return false;
}

/*  overlap function allow user to change the local tensor field only  */
void cal_tensorvals_quad_inReg(int regionid)
{
	int i;
	QuadVertex *v;
	double t[4]={0.};

	/*  see whether we have elements fall inside this region  */
	/*  if it is not, we use the stored global field  */
	bool existElem=exist_Elems_inReg(regionid);

	if(!existElem)
	{
		/* use stored global field  */

		for(i=0; i<quadmesh->nverts; i++)
		{
			v = quadmesh->quad_verts[i];

			if(!v->inland)
				continue;

			if(v->which_region != regionid)
				continue;

			v->Jacobian.entry[0][0] = pre_tenfield[i].entry[0][0];
			v->Jacobian.entry[0][1] = pre_tenfield[i].entry[0][1];
			v->Jacobian.entry[1][0] = pre_tenfield[i].entry[1][0];
			v->Jacobian.entry[1][1] = pre_tenfield[i].entry[1][1];
		}
		return;
	}


	for(i=0; i<quadmesh->nverts; i++)
	{
		v = quadmesh->quad_verts[i];

		if(!v->inland)
			continue;

		if(v->which_region != regionid)
			continue;

		//get_tensor_inReg(v->x, v->y, t, v->which_region, v->inbrushregion);
		get_tensor_inReg(v->x, v->y, t, regionid, v->inbrushregion);
		if(please_comb_prefield/*sharedvars.SelStreetRegToEditOn*/)
		{
			if(v->which_region == cur_max_reg_index)
			{
				v->Jacobian.entry[0][0] = t[0];
				v->Jacobian.entry[0][1] = t[1];
				v->Jacobian.entry[1][0] = t[2];
				v->Jacobian.entry[1][1] = t[3];
			}
			else{
				v->Jacobian.entry[0][0] = 0.1*pre_tenfield[i].entry[0][0]+t[0];
				v->Jacobian.entry[0][1] = 0.1*pre_tenfield[i].entry[0][1]+t[1];
				v->Jacobian.entry[1][0] = 0.1*pre_tenfield[i].entry[1][0]+t[2];
				v->Jacobian.entry[1][1] = 0.1*pre_tenfield[i].entry[1][1]+t[3];
			}
		}
		else{
			//get_tensor_inReg(v->x, v->y, t, v->which_region, v->inbrushregion);
			v->Jacobian.entry[0][0] = t[0];
			v->Jacobian.entry[0][1] = t[1];
			v->Jacobian.entry[1][0] = t[2];
			v->Jacobian.entry[1][1] = t[3];
		}

		//if(/*combine with previous 12/2/2007*/)
		//{
		//	v->Jacobian.entry[0][0]=(1-ten_blend_factor)*t[0]+pre_ten[i].entry[0][0];
		//	v->Jacobian.entry[0][1]=(1-ten_blend_factor)*t[1]+pre_ten[i].entry[0][1];
		//	v->Jacobian.entry[1][0]=(1-ten_blend_factor)*t[2]+pre_ten[i].entry[1][0];
		//	v->Jacobian.entry[1][1]=(1-ten_blend_factor)*t[3]+pre_ten[i].entry[1][1];
		//}
	}
}



void cal_eigenvecs_quad_inReg(int regionid)
{
	int i;

	for(i=0; i<quadmesh->nverts; i++)
	{
		if(quadmesh->quad_verts[i]->inland
			&& quadmesh->quad_verts[i]->which_region==regionid)
			cal_eigenvecs_onevert_quad(i);
	}

	/*normalize the major and minor field*/
	normalized_tensorfield_quad();
}


/**************************************************************************************/


/*Update the transform matrix of sepecific singular element
*/
extern double sx, sy, uniforms;
extern double rotateAng;

void update_tenSingularElem_transform(int TransformType, int elem_id)
{
	icMatrix3x3 trans, intermedia, inversetranslate;

	////Note that we need to translate to the center of the element
	////then perform the specified transformation
	////and translate back
	ten_designelems[elem_id].transform_matrix.set(\
		1, 0, ten_designelems[elem_id].centerx,\
		0, 1, ten_designelems[elem_id].centery,\
		0, 0, 1);

	inversetranslate.set(1, 0, -ten_designelems[elem_id].centerx,
		0, 1, -ten_designelems[elem_id].centery,
		0, 0, 1);

	if(TransformType == 1) ////If we scale the element uniformly
	{
		ten_designelems[elem_id].s += uniforms;

		////if the scalar less than 0, it will change the type of singular element!!!
		if(ten_designelems[elem_id].s <= 0)
			ten_designelems[elem_id].s = 0.1;
	}

	else if(TransformType == 2)
	{
		ten_designelems[elem_id].sx += sx;
		ten_designelems[elem_id].sy += sy;

		////if the scalar less than 0, it will change the type of singular element!!!
		if(ten_designelems[elem_id].sx <= 0)
			ten_designelems[elem_id].sx = 0.1;
		if(ten_designelems[elem_id].sy <= 0)
			ten_designelems[elem_id].sy = 0.1;
	}
	else{  ////it is rotation
		ten_designelems[elem_id].rotang += rotateAng;
	}
		

	////set the rotation matrix
	//intermedia.set(cos(ten_designelems[elem_id].rotang), -sin(ten_designelems[elem_id].rotang), 0,
	//				sin(ten_designelems[elem_id].rotang), cos(ten_designelems[elem_id].rotang),  0,
	//				0,0,1);
	
	intermedia.set(1, 0, 0,
					0, 1,  0,
					0,0,1);

	ten_designelems[elem_id].transform_matrix.rightMultiply(intermedia);

	////set the uniform scaling matrix
	intermedia.set(ten_designelems[elem_id].s, 0, 0,
					0, ten_designelems[elem_id].s, 0,
					0, 0, 1);
	ten_designelems[elem_id].transform_matrix.rightMultiply(intermedia);

	////set the non uniform scaling matrix
	intermedia.set(ten_designelems[elem_id].sx, 0, 0,
					0, ten_designelems[elem_id].sy, 0,
					0, 0, 1);

	ten_designelems[elem_id].transform_matrix.rightMultiply(intermedia);
	ten_designelems[elem_id].transform_matrix.rightMultiply(inversetranslate);

	update_tenElem_EditBox(elem_id);
}

void update_tenElem_EditBox(int elem_id)
{
	double origin[3];

	icMatrix3x3 trans, transposerot, inversetranslate, intermedia;

	trans.set(1, 0, ten_designelems[elem_id].centerx,
		0, 1, ten_designelems[elem_id].centery,
		0, 0, 1);
	inversetranslate.set(1, 0, -ten_designelems[elem_id].centerx,
		0, 1, -ten_designelems[elem_id].centery,
		0, 0, 1);

	////set the rotation matrix
	intermedia.set(cos(ten_designelems[elem_id].rotang), -sin(ten_designelems[elem_id].rotang), 0,
					sin(ten_designelems[elem_id].rotang), cos(ten_designelems[elem_id].rotang),  0,
					0,0,1);

	trans.rightMultiply(intermedia);

	////set the uniform scaling matrix
	intermedia.set(ten_designelems[elem_id].s, 0, 0,
					0, ten_designelems[elem_id].s, 0,
					0, 0, 1);
	trans.rightMultiply(intermedia);

	////set the non-uniform scaling matrix
	intermedia.set(ten_designelems[elem_id].sx, 0, 0,
					0, ten_designelems[elem_id].sy, 0,
					0, 0, 1);
	trans.rightMultiply(intermedia);

	trans.rightMultiply(inversetranslate);


	////Update p1
	origin[0] = ten_designelems[elem_id].editbox.p1.entry[0];
	origin[1] = ten_designelems[elem_id].editbox.p1.entry[1];
	origin[2] = 1;

		ten_designelems[elem_id].cur_editbox.p1.entry[0] = trans.entry[0][0] * origin[0]\
		+ trans.entry[0][1] * origin[1]\
		+ trans.entry[0][2] * origin[2];
	
	ten_designelems[elem_id].cur_editbox.p1.entry[1] = trans.entry[1][0] * origin[0]\
		+ trans.entry[1][1] * origin[1]\
		+ trans.entry[1][2] * origin[2];

	////Update p2
	origin[0] = ten_designelems[elem_id].editbox.p2.entry[0];
	origin[1] = ten_designelems[elem_id].editbox.p2.entry[1];
	origin[2] = 1;

	ten_designelems[elem_id].cur_editbox.p2.entry[0] = trans.entry[0][0] * origin[0]\
		+ trans.entry[0][1] * origin[1]\
		+ trans.entry[0][2] * origin[2];
	
	ten_designelems[elem_id].cur_editbox.p2.entry[1] = trans.entry[1][0] * origin[0]\
		+ trans.entry[1][1] * origin[1]\
		+ trans.entry[1][2] * origin[2];

	////Update p3
	origin[0] = ten_designelems[elem_id].editbox.p3.entry[0];
	origin[1] = ten_designelems[elem_id].editbox.p3.entry[1];
	origin[2] = 1;

	ten_designelems[elem_id].cur_editbox.p3.entry[0] = trans.entry[0][0] * origin[0]\
		+ trans.entry[0][1] * origin[1]\
		+ trans.entry[0][2] * origin[2];
	
	ten_designelems[elem_id].cur_editbox.p3.entry[1] = trans.entry[1][0] * origin[0]\
		+ trans.entry[1][1] * origin[1]\
		+ trans.entry[1][2] * origin[2];
	
	////Update p4
	origin[0] = ten_designelems[elem_id].editbox.p4.entry[0];
	origin[1] = ten_designelems[elem_id].editbox.p4.entry[1];
	origin[2] = 1;

	ten_designelems[elem_id].cur_editbox.p4.entry[0] = trans.entry[0][0] * origin[0]\
		+ trans.entry[0][1] * origin[1]\
		+ trans.entry[0][2] * origin[2];
	
	ten_designelems[elem_id].cur_editbox.p4.entry[1] = trans.entry[1][0] * origin[0]\
		+ trans.entry[1][1] * origin[1]\
		+ trans.entry[1][2] * origin[2];

	////Update Up
	origin[0] = ten_designelems[elem_id].editbox.Up.entry[0];
	origin[1] = ten_designelems[elem_id].editbox.Up.entry[1];
	origin[2] = 1;

	ten_designelems[elem_id].cur_editbox.Up.entry[0] = trans.entry[0][0] * origin[0]\
		+ trans.entry[0][1] * origin[1]\
		+ trans.entry[0][2] * origin[2];
	
	ten_designelems[elem_id].cur_editbox.Up.entry[1] = trans.entry[1][0] * origin[0]\
		+ trans.entry[1][1] * origin[1]\
		+ trans.entry[1][2] * origin[2];


	//ten_designelems[elem_id].cur_editbox.p1.entry[0] = ten_designelems[elem_id].transform_matrix.entry[0][0] * origin[0]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[0][1] * origin[1]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[0][2] * origin[2];
	//
	//ten_designelems[elem_id].cur_editbox.p1.entry[1] = ten_designelems[elem_id].transform_matrix.entry[1][0] * origin[0]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[1][1] * origin[1]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[1][2] * origin[2];

	//////Update p2
	//origin[0] = ten_designelems[elem_id].editbox.p2.entry[0];
	//origin[1] = ten_designelems[elem_id].editbox.p2.entry[1];
	//origin[2] = 1;

	//ten_designelems[elem_id].cur_editbox.p2.entry[0] = ten_designelems[elem_id].transform_matrix.entry[0][0] * origin[0]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[0][1] * origin[1]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[0][2] * origin[2];
	//
	//ten_designelems[elem_id].cur_editbox.p2.entry[1] = ten_designelems[elem_id].transform_matrix.entry[1][0] * origin[0]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[1][1] * origin[1]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[1][2] * origin[2];

	//////Update p3
	//origin[0] = ten_designelems[elem_id].editbox.p3.entry[0];
	//origin[1] = ten_designelems[elem_id].editbox.p3.entry[1];
	//origin[2] = 1;

	//ten_designelems[elem_id].cur_editbox.p3.entry[0] = ten_designelems[elem_id].transform_matrix.entry[0][0] * origin[0]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[0][1] * origin[1]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[0][2] * origin[2];
	//
	//ten_designelems[elem_id].cur_editbox.p3.entry[1] = ten_designelems[elem_id].transform_matrix.entry[1][0] * origin[0]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[1][1] * origin[1]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[1][2] * origin[2];
	//
	//////Update p4
	//origin[0] = ten_designelems[elem_id].editbox.p4.entry[0];
	//origin[1] = ten_designelems[elem_id].editbox.p4.entry[1];
	//origin[2] = 1;

	//ten_designelems[elem_id].cur_editbox.p4.entry[0] = ten_designelems[elem_id].transform_matrix.entry[0][0] * origin[0]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[0][1] * origin[1]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[0][2] * origin[2];
	//
	//ten_designelems[elem_id].cur_editbox.p4.entry[1] = ten_designelems[elem_id].transform_matrix.entry[1][0] * origin[0]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[1][1] * origin[1]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[1][2] * origin[2];

	//////Update Up
	//origin[0] = ten_designelems[elem_id].editbox.Up.entry[0];
	//origin[1] = ten_designelems[elem_id].editbox.Up.entry[1];
	//origin[2] = 1;

	//ten_designelems[elem_id].cur_editbox.Up.entry[0] = ten_designelems[elem_id].transform_matrix.entry[0][0] * origin[0]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[0][1] * origin[1]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[0][2] * origin[2];
	//
	//ten_designelems[elem_id].cur_editbox.Up.entry[1] = ten_designelems[elem_id].transform_matrix.entry[1][0] * origin[0]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[1][1] * origin[1]\
	//	+ ten_designelems[elem_id].transform_matrix.entry[1][2] * origin[2];
}

/*
Update the transformation matrix of a specific design elements
*/
void update_tenElem_Transform(int TransformType, int choose_ID)
{
	icMatrix3x3 trans, intermedia, inversetranslate;

	////Set the transformation
	////Suppose that we perform scale first

	if(TransformType > 0)
	{
		////Update the Jacobian matrix of selected element
		if(chosen_tenelem_ID >= 0 && chosen_tenelem_ID < NAMEOFREGELEM)
		{
             update_tenSingularElem_transform(TransformType, chosen_tenelem_ID);
		}
	
		//else if(chosen_tenelem_ID >= NAMEOFREGELEM && chosen_tenelem_ID < NAMEOFSINGCONTROL )
		//{ ////we choose a regular, convergent or divergent element
  //           update_tenRegElem_transform(TransformType, chosen_tenelem_ID - NAMEOFREGELEM);
		//}

	}
}


/*set the basis position of the regular design element*/
void set_ten_regBasis(double x, double y, int type)
{
	////Add to the regular elements list
	if(nten_regelems >= curMaxNumTenRegElems-1 )
	{
		//Allocate new space for the element list
		curMaxNumTenRegElems += 50;
		ten_regularelems = (TenRegularElem*)realloc(ten_regularelems, sizeof(TenRegularElem) * curMaxNumTenRegElems);
		if(ten_regularelems == NULL)
			exit(-1);

	}
	
	ten_regularelems[nten_regelems].end[0] = ten_regularelems[nten_regelems].base[0] = x;
	ten_regularelems[nten_regelems].end[1] = ten_regularelems[nten_regelems].base[1] = y;
	ten_regularelems[nten_regelems].ID = NAMEOFREGELEM + nten_regelems+1;
	ten_regularelems[nten_regelems].type = type;
	ten_regularelems[nten_regelems].Direct.set(0, 0);

	////Initialize the transformation parameters for this regular element
	ten_regularelems[nten_regelems].transform_matrix.setIdentity();
	//ten_regularelems[nten_regelems].transposeRot.setIdentity();

	ten_regularelems[nten_regelems].rotang = 0;
	ten_regularelems[nten_regelems].s = 1;

	ten_regularelems[nten_regelems].deleted=false;

	/*mark which region the design element belongs to 11/21/2007*/
	ten_regularelems[nten_regelems].which_region=get_region_id(x, y);
	
	/*  may be a bug 1/18/2008 */
	cur_chosen_region = ten_regularelems[nten_regelems].which_region;

	nten_regelems ++;
}

	
/*set direction function for regular element before 09/26/2007*/

void set_ten_regDir(double x, double y)
{
	////Update the Cur_regularID here

	ten_regularelems[nten_regelems - 1].end[0] = x ;
	ten_regularelems[nten_regelems - 1].end[1] = y ;

	ten_regularelems[nten_regelems - 1].Direct.entry[0] = x - ten_regularelems[nten_regelems - 1].base[0];
	ten_regularelems[nten_regelems - 1].Direct.entry[1] = y - ten_regularelems[nten_regelems - 1].base[1];

	ten_regularelems[nten_regelems - 1].rotang = 
		atan2(ten_regularelems[nten_regelems - 1].Direct.entry[1],
			ten_regularelems[nten_regelems - 1].Direct.entry[0]);


}

void update_ten_regDir(int id, double x, double y)
{
	////Update the Cur_regularID here

	ten_regularelems[id].end[0] = x ;
	ten_regularelems[id].end[1] = y ;

	ten_regularelems[id].Direct.entry[0] = x - ten_regularelems[id].base[0];
	ten_regularelems[id].Direct.entry[1] = y - ten_regularelems[id].base[1];

	ten_regularelems[id].rotang = 
		atan2(ten_regularelems[id].Direct.entry[1],
			ten_regularelems[id].Direct.entry[0]);

}


/*
save the tensor field per vertex
*/

void save_tenField_perVer(char *filename)
{
	FILE *vfp = fopen(filename, "w");
	if(vfp == NULL)
	{
		MessageBox(NULL, "Can't open the file!", "Error", MB_OK);
		return;
	}

	int i; 
	for(i = 0; i < quadmesh->nverts; i++)
	{
		fprintf(vfp,"%f,%f,%f,%f\n", quadmesh->quad_verts[i]->Jacobian.entry[0][0], 
			quadmesh->quad_verts[i]->Jacobian.entry[0][1],
			quadmesh->quad_verts[i]->Jacobian.entry[1][0],
			quadmesh->quad_verts[i]->Jacobian.entry[1][1]);
	}
	fclose(vfp);
}


/*
load the pre-stored tensor field
*/
void load_tenField_perVer(char *filename)
{
    float t00, t01, t10, t11;

	FILE *vfp = fopen(filename,"r");

	if(vfp == NULL)
	{
		MessageBox(NULL, "Can't open the file!", "Error", MB_OK);
		return;
	}

	for(int i = 0; i < quadmesh->nverts; i++)
	{
		fscanf(vfp, "%f,%f,%f,%f\n", &t00, &t01, &t10, &t11);
		quadmesh->quad_verts[i]->Jacobian.entry[0][0] = t00;
		quadmesh->quad_verts[i]->Jacobian.entry[0][1] = t01;
		quadmesh->quad_verts[i]->Jacobian.entry[1][0] = t10;
		quadmesh->quad_verts[i]->Jacobian.entry[1][1] = t11;

	}
	fclose(vfp);

	save_cur_field();
}

/*
   save the design elements of current tensor fields
*/

bool save_tenField_elems(char *filename)
{
	FILE *fp=fopen(filename, "w");

	if(fp == NULL)
		return false;

	int counter=0;
	int i;

	/*  save singular elements  */
	for(i=0;i<ntenelems;i++)
	{
		if(ten_designelems[i].deleted)
			continue;
		counter++;
	}

	fprintf(fp, "#singular elements: %d\n", counter);

	counter=0;
	for(i=0;i<ntenelems; i++)
	{
		if(ten_designelems[i].deleted)
			continue;

		fprintf(fp,"elem %d: \n",counter);
		fprintf(fp,"type: %d\n", ten_designelems[i].type);
		fprintf(fp,"pos: %f, %f\n", ten_designelems[i].centerx, ten_designelems[i].centery);
		/*  we may also need to store the user editing results  */
		fprintf(fp,"rotation: %f\n", ten_designelems[i].rotang);
		fprintf(fp,"scale: %f, %f, %f\n", ten_designelems[i].s, 
			ten_designelems[i].sx, ten_designelems[i].sy);
		fprintf(fp,"tranform matrix: %f, %f, %f, %f\n", 
			ten_designelems[i].transform_matrix.entry[0][0],
			ten_designelems[i].transform_matrix.entry[0][1],
			ten_designelems[i].transform_matrix.entry[1][0],
			ten_designelems[i].transform_matrix.entry[1][1]);

		counter++;
	}

	/*  save regular elements  */
	counter=0;
	for(i=0;i<nten_regelems;i++)
	{
		if(ten_regularelems[i].deleted)
			continue;
		counter++;
	}

	fprintf(fp, "#regular elements: %d\n", counter);

	counter=0;
	for(i=0;i<nten_regelems;i++)
	{
		if(ten_regularelems[i].deleted)
			continue;

		fprintf(fp,"elem %d: \n", counter);
		fprintf(fp,"basis: %f, %f\n", ten_regularelems[i].base[0],ten_regularelems[i].base[1]);
		fprintf(fp,"end: %f, %f\n", ten_regularelems[i].end[0],ten_regularelems[i].end[1]);
		fprintf(fp,"Direction: %f, %f\n", ten_regularelems[i].Direct.entry[0],
			ten_regularelems[i].Direct.entry[1]);
		fprintf(fp,"rotation: %f\n", ten_regularelems[i].rotang);
		fprintf(fp,"scale: %f\n", ten_regularelems[i].s);
		counter++;
	}

	fclose(fp);
	return true;
}

bool load_tenField_elems(char *filename)
{
	FILE *fp=fopen(filename, "r");

	if(fp == NULL)
		return false;

	int counter=0;
	int i;
	float x, y, a0, a1, a2, a3;
	float s, sx, sy, rotang;
	int type;
	int elemindex;

	/*  load singular elements  */

	fscanf(fp, "#singular elements: %d\n", &counter);

	/*   enlarge the space if not enough   */
	if(counter >= curMaxNumTenDesignElems)
	{
		ten_designelems = (Degenerate_Design*)realloc(ten_designelems, sizeof(Degenerate_Design)*
			(curMaxNumTenDesignElems+counter));
		if(ten_designelems == NULL)
			exit(-1);
		curMaxNumTenDesignElems += counter;
	}

	ntenelems=0;
	for(i=0;i<counter; i++)
	{
		fscanf(fp,"elem %d: \n",&elemindex);
		fscanf(fp,"type: %d\n", &type);
		fscanf(fp,"pos: %f, %f\n", &x, &y);
		fscanf(fp,"rotation: %f\n", &rotang);
		fscanf(fp,"scale: %f, %f, %f\n", &s, &sx, &sy);
		fscanf(fp,"tranform matrix: %f, %f, %f, %f\n", 
			&a0, &a1, &a2, &a3);

		ten_designelems[ntenelems].centerx = x;
		ten_designelems[ntenelems].centery = y;
		ten_designelems[ntenelems].ID = ntenelems;
		ten_designelems[ntenelems].Triangle_ID = get_cellID_givencoords(x,y);
		ten_designelems[ntenelems].type = type;
		ten_designelems[ntenelems].transform_matrix.set(a0,a1,0,   a2,a3,0,   0,0,1);
		
		ten_designelems[ntenelems].rotang = rotang;
		ten_designelems[ntenelems].s = s;
		ten_designelems[ntenelems].sx =sx;
		ten_designelems[ntenelems].sy =sy;
		ten_designelems[ntenelems].deleted = false;

		/*mark which region the design element belongs to 11/21/2007*/
		ten_designelems[ntenelems].which_region=get_region_id(x, y);
				
		//Add the editbox for properties editing
		init_tenelem_EditBox(ntenelems, x, y);

		update_tenElem_EditBox(ntenelems);

		ntenelems++;
	}

	/*  load regular elements  */

	fscanf(fp, "#regular elements: %d\n", &counter);
	if(counter >= curMaxNumTenRegElems-1 )
	{
		//Allocate new space for the element list
		curMaxNumTenRegElems += counter;
		ten_regularelems = (TenRegularElem*)realloc(ten_regularelems, sizeof(TenRegularElem) * curMaxNumTenRegElems);
		if(ten_regularelems == NULL)
			exit(-1);

	}

	nten_regelems=0;
	for(i=0;i<counter;i++)
	{
		fscanf(fp,"elem %d: \n", &elemindex);
		fscanf(fp,"basis: %f, %f\n", &x, &y);
		fscanf(fp,"end: %f, %f\n", &a0,&a1);
		fscanf(fp,"Direction: %f, %f\n", &a2, &a3);
		fscanf(fp,"rotation: %f\n", &rotang);
		fscanf(fp,"scale: %f\n", &s);

		ten_regularelems[nten_regelems].base[0] = x;
		ten_regularelems[nten_regelems].base[1] = y;
		ten_regularelems[nten_regelems].end[0] = a0;
		ten_regularelems[nten_regelems].end[1] = a1;
		ten_regularelems[nten_regelems].ID = NAMEOFREGELEM + nten_regelems+1;
		ten_regularelems[nten_regelems].type = 0;
		ten_regularelems[nten_regelems].Direct.set(a2, a3);

		////Initialize the transformation parameters for this regular element
		ten_regularelems[nten_regelems].transform_matrix.setIdentity();
		//ten_regularelems[nten_regelems].transposeRot.setIdentity();

		ten_regularelems[nten_regelems].rotang = rotang;
		ten_regularelems[nten_regelems].s = s;

		ten_regularelems[nten_regelems].deleted=false;

		/*mark which region the design element belongs to 11/21/2007*/
		ten_regularelems[nten_regelems].which_region=get_region_id(x, y);

		nten_regelems++;
	}

	fclose(fp);
	return true;
}

/*
    store current tensor field into the global variable "pre_field"
*/
void save_cur_field()
{
	int i;
	for(i=0;i<quadmesh->nverts;i++)
	{
		pre_tenfield[i].set(quadmesh->quad_verts[i]->Jacobian);
	}
}
