
////VFSynthesis.cpp
#include "stdafx.h"
#include "VFSynthesis.h"
#include "VFDataStructure.h"
#include "LocalTracing.h"

#include "RegionSmoothing.h"
//#include "HermiteCurve.h"

#include "LimitCycleCreator.h"

#include "GL/glut.h"


//////////////////////////////////////////////////
////Part of the variables used by this module

extern int NewAddedElementType;         ////The original variable is defined in IBFVDlg.cpp

extern int MaxNumSingularElems;
extern int MaxNumRegularElems;           //Maximum number of regular elements
extern int MaxNumTrajectories;           //Maximum number of possible trajectories
                                         //(it should be flexible for future pen-and-ink sketch)
extern int MaxNumLinesegsPerTraj;               //Maximum number of line segments for each trajectory

extern SingularElement *singularelem;          //Singular elememts' list
extern int cur_singelem_index;

extern RegularElement *regularelem;            //regular elememts' list
extern int cur_regelem_index;

extern Singularities *singularities;           //being captured singularites' list
extern int cur_singularity_index;

extern LineSeg **trajectories;                 //trajectories' list
extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory

extern int cur_separatrices_index;

extern Polygon3D Object;
extern icVector2 *backup_field;

extern double  dmax ;

////Transformation variables for element edit
//extern double rotmatrix[3][3];
extern double sx, sy, uniforms;
extern double rotateAng;

extern int prev_num_reg_elem;


extern Vertex **regionverts;                ////mesh vertices inside user selected region
extern int Num_verts;                       ////number of inner vertices



extern void  InitLimitCycleDetection();
extern void InitLimitCycleStructure();

double max_mag, min_mag;

bool UseNormalizedVF = true;


/////////////////////////////////////////////////////////
//Add basic singular element into the singular element list
void AddElement(int type, double x, double y)
{

	if(cur_singelem_index >= MaxNumSingularElems -1 )
	{
		//Allocate new space for the element list
		MaxNumSingularElems += 50;
		singularelem = (SingularElement*)realloc(singularelem, sizeof(SingularElement) * MaxNumSingularElems);
	}
		

	//if still have place, then add to the element list
	singularelem[cur_singelem_index].centerx = x;
	singularelem[cur_singelem_index].centery = y;
	
	//singularelem[cur_singelem_index].Triangle_ID = TriangleDetect(x, y);

	singularelem[cur_singelem_index].type = type;
	singularelem[cur_singelem_index].ID = cur_singelem_index+1;

	////Initialize the Jacobian matrix for the element
	InitialJacobian(singularelem[cur_singelem_index].Jacobian, type);
	
	////Initialize the transform matrix of the element
	singularelem[cur_singelem_index].transform_matrix.setIdentity();

	singularelem[cur_singelem_index].rotang = 0;
	singularelem[cur_singelem_index].s = 
	singularelem[cur_singelem_index].sx =
	singularelem[cur_singelem_index].sy = 1;
	singularelem[cur_singelem_index].deleted = false;
			
	//Add the editbox for properties editing
	InitEditBox(cur_singelem_index, x, y);

	cur_singelem_index ++;
}


//the overloaded function for add a new element
void AddElement(int type, int triangle)
{
	if(cur_singelem_index >= MaxNumSingularElems-1 )
	{
		//Allocate new space for the element list
		MaxNumSingularElems += 50;
		singularelem = (SingularElement*)realloc(singularelem, sizeof(SingularElement) * MaxNumSingularElems);
	}
	
	///////////////////////////////////////////////////
	///We pick the center of the triangle is the location of the singularity
	//singularelem[cur_singelem_index].gcenter ;
	
	singularelem[cur_singelem_index].Triangle_ID = triangle;

	singularelem[cur_singelem_index].type = type;
	singularelem[cur_singelem_index].ID = cur_singelem_index+1;

	////Initialize the Jacobian matrix for the element
	InitialJacobian(singularelem[cur_singelem_index].Jacobian, type);
	
	////Initialize the transform matrix of the element
	singularelem[cur_singelem_index].transform_matrix.setIdentity();

	singularelem[cur_singelem_index].rotang = 0;
	singularelem[cur_singelem_index].s = 
	singularelem[cur_singelem_index].sx =
	singularelem[cur_singelem_index].sy = 1;
	singularelem[cur_singelem_index].deleted = false;
			
	//Add the editbox for properties editing
	//InitEditBox(cur_singelem_index, x, y);

	cur_singelem_index ++;
}



/************************************************************
Set the three vectors at the 3 vectices of the selected triangle
according to the type of the singularity
************************************************************/

void SetVectorsAtTriangleforNewElem(int type, int triangle)
{
	int i;
	icVector2 localv;
	icVector2 globalv;
	Face *face;

	if(triangle < 0)
		return;

	face = Object.flist[triangle];


    ////New method to set the three vectors, 
	////this method only works well for regular triangle
	if(type == SOURCE)
	{
		face->direct_vec[0].entry[0]=
			face->xy[0][0] - (face->xy[1][0]+face->xy[2][0])/2;
		face->direct_vec[0].entry[1]=
			face->xy[0][1] - (face->xy[1][1]+face->xy[2][1])/2;
		normalize(face->direct_vec[0]);
		
		face->direct_vec[1].entry[0]=
			face->xy[1][0] - (face->xy[0][0]+face->xy[2][0])/2;
		face->direct_vec[1].entry[1]=
			face->xy[1][1] - (face->xy[0][1]+face->xy[2][1])/2;
		normalize(face->direct_vec[1]);

		face->direct_vec[2].entry[0]=
			face->xy[2][0] - (face->xy[1][0]+face->xy[0][0])/2;
		face->direct_vec[2].entry[1]=
			face->xy[2][1] - (face->xy[1][1]+face->xy[0][1])/2;
		normalize(face->direct_vec[2]);
	}

	else if(type == SINK)
	{
		face->direct_vec[0].entry[0]=
			(face->xy[1][0]+face->xy[2][0])/2 - face->xy[0][0];
		face->direct_vec[0].entry[1]=
			(face->xy[1][1]+face->xy[2][1])/2 - face->xy[0][1];
		normalize(face->direct_vec[0]);
		
		face->direct_vec[1].entry[0]=
			(face->xy[0][0]+face->xy[2][0])/2 - face->xy[1][0];
		face->direct_vec[1].entry[1]=
			(face->xy[0][1]+face->xy[2][1])/2 - face->xy[1][1];
		normalize(face->direct_vec[1]);

		face->direct_vec[2].entry[0]=
			(face->xy[1][0]+face->xy[0][0])/2 - face->xy[2][0];
		face->direct_vec[2].entry[1]=
			(face->xy[1][1]+face->xy[0][1])/2 - face->xy[2][1];
		normalize(face->direct_vec[2]);
	}

	//we may adopt other method to set the vectors on the three vertices to create a saddle!! 12/01/05
	else if(type == SADDLE)  ////May have some exception here
	{
		//double cos30 = sqrt(3.)/2;
		double cos30 = sqrt(2.)/2;

		face->direct_vec[0].entry[0]=
			(face->xy[1][0]+face->xy[2][0])/2 - face->xy[0][0];
		face->direct_vec[0].entry[1]=
			(face->xy[1][1]+face->xy[2][1])/2 - face->xy[0][1];
		normalize(face->direct_vec[0]);
		

		face->direct_vec[2].entry[0]=
			-cos30 * face->direct_vec[0].entry[0] - cos30/*0.5*/ * face->direct_vec[0].entry[1];
		face->direct_vec[2].entry[1]=
			cos30/*0.5*/ * face->direct_vec[0].entry[0] - cos30 * face->direct_vec[0].entry[1];
		normalize(face->direct_vec[1]);
		
		face->direct_vec[1].entry[0]=
			-cos30 * face->direct_vec[0].entry[0] + cos30/*0.5*/ * face->direct_vec[0].entry[1];
		face->direct_vec[1].entry[1]=
			-cos30/*0.5*/ * face->direct_vec[0].entry[0] - cos30 * face->direct_vec[0].entry[1];
		normalize(face->direct_vec[2]);
	}

	else if(type == CWCENTER)
	{
		face->direct_vec[0].entry[0]=
			-((face->xy[1][1]+face->xy[2][1])/2 - face->xy[0][1]);
		face->direct_vec[0].entry[1]=
			(face->xy[1][0]+face->xy[2][0])/2 - face->xy[0][0];
		normalize(face->direct_vec[0]);
		
		face->direct_vec[1].entry[0]=
			-((face->xy[0][1]+face->xy[2][1])/2 - face->xy[1][1]);
		face->direct_vec[1].entry[1]=
			(face->xy[0][0]+face->xy[2][0])/2 - face->xy[1][0];
		normalize(face->direct_vec[1]);

		face->direct_vec[2].entry[0]=
			-((face->xy[1][1]+face->xy[0][1])/2 - face->xy[2][1]);
		face->direct_vec[2].entry[1]=
			(face->xy[1][0]+face->xy[0][0])/2 - face->xy[2][0];
		normalize(face->direct_vec[2]);
	}

	else if(type == CCWCENTER)
	{
		face->direct_vec[0].entry[0]=
			(face->xy[1][1]+face->xy[2][1])/2 - face->xy[0][1];
		face->direct_vec[0].entry[1]=
			-((face->xy[1][0]+face->xy[2][0])/2 - face->xy[0][0]);
		normalize(face->direct_vec[0]);
		
		face->direct_vec[1].entry[0]=
			(face->xy[0][1]+face->xy[2][1])/2 - face->xy[1][1];
		face->direct_vec[1].entry[1]=
			-((face->xy[0][0]+face->xy[2][0])/2 - face->xy[1][0]);
		normalize(face->direct_vec[1]);

		face->direct_vec[2].entry[0]=
			(face->xy[1][1]+face->xy[0][1])/2 - face->xy[2][1];
		face->direct_vec[2].entry[1]=
			-((face->xy[1][0]+face->xy[0][0])/2 - face->xy[2][0]);
		normalize(face->direct_vec[2]);
	}

	for(i = 0; i < 3; i++)
	{
		globalv = face->direct_vec[i].entry[0]*face->LX + face->direct_vec[i].entry[1]*face->LY;
		Object.vlist[face->verts[i]]->vec = 0.10 * globalv;
	}
}


////Set the vectors for the regular element 1/2/06
void SetVectorsForNewRegularElem(int type, int triangle, icVector2 global_vec)
{
	int i;
	Face *face;
	face = Object.flist[triangle];
	icVector2 local_vec, globalv;

	if (type == 0) ////for regular element
	{
		////first, transfer to local direction
		local_vec.entry[0] = dot(global_vec, face->LX);
		local_vec.entry[1] = dot(global_vec, face->LY);

		////then, set the vectors on the three vertices of the triangle 
		face->direct_vec[0] = face->direct_vec[1] = face->direct_vec[2] = local_vec;
	}
	
	for(i = 0; i < 3; i++)
	{
		globalv = face->direct_vec[i].entry[0]*face->LX + face->direct_vec[i].entry[1]*face->LY;
		Object.vlist[face->verts[i]]->vec = 0.15 * globalv;
	}
}

/*------------------------------------------*/
//Normalize the whole field
void NormalizeField()
{
	int i;
    double r;
	Vertex *cur_v;

	for(i = 0; i < Object.nverts; i++)
	{
		cur_v = Object.vlist[i];
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
}


/*---------------------------------------------------------------------------------*/
//create the field through region smoothing
void BuildtheField()
{
	Face *face;
	int i, j;

	for(i = 0; i < cur_singelem_index; i++)
	{
		face = Object.flist[singularelem[i].Triangle_ID];
		////Set the three vectors at the selected triangle
		SetVectorsAtTriangleforNewElem(singularelem[i].type, singularelem[i].Triangle_ID);

		for(j = 0; j < face->nverts; j++)
		{
			Object.vlist[face->verts[j]]->OnBoundary = 1;
		}
	}

	////Initial the regular element
	for(i = 0; i < cur_regelem_index; i++)
	{
		face = Object.flist[regularelem[i].basis_triangle];

		////Set the three vectors at the selected triangle
		SetVectorsForNewRegularElem(regularelem[i].type, regularelem[i].basis_triangle, regularelem[i].Direct);
		
		for(j = 0; j < face->nverts; j++)
		{
			Object.vlist[face->verts[j]]->OnBoundary = 1;
		}
	}

	///Build the smoothing region
    BuildRegionforAddElement();

	////smooth the region to get the vectors on all other vertices
	RegionSmooth();

	NormalizeField();

	//// Calculate the local vector field
	GetLocalVector();

}


////Save the field after smoothing process
////since after smoothing, we remove some fixed points including those
////element points, so we need to use current detected fixed points as current element points
////to let user keep adding new elements 4/06/06, not finished yet!
void SaveFieldAfterSmoothing()
{
	int i;

	cur_singelem_index = 0;

	for(i = 0; i < cur_singularity_index; i++)
	{
		singularelem[i].Triangle_ID = singularities[i].Triangle_ID;
		singularelem[i].centerx = singularities[i].gcx;
		singularelem[i].centery = singularities[i].gcy;
		singularelem[i].Jacobian.set(singularities[i].Jacobian);
		singularelem[i].type = singularities[i].type;
	}

	cur_singelem_index = cur_singularity_index;
}





/*-----------------------------------------------------------*/

void BuildRegionforAddElement()
{
	int i;
	Vertex *cur_vert;

	Num_verts = 0;

	for(i = 0; i < Object.nverts; i++)
	{
		cur_vert = Object.vlist[i];
		if(cur_vert->OnBoundary == 0)
		{
			regionverts[Num_verts] = cur_vert;
			cur_vert->RegionListID = Num_verts;
			cur_vert->InRegion = 1;
			Num_verts ++;
		}
	}
}



//////////////////////////////////////////////////////////////////
////Add regular, convergent or divergent element
void SetBasis(CPoint point, int type)
{
	GLdouble position[3];
	ScreenToWorld(position, point);

	////Add to the regular elements list
	if(cur_regelem_index >= MaxNumRegularElems-1 )
	{
		//Allocate new space for the element list
		MaxNumRegularElems += 50;
		regularelem = (RegularElement*)realloc(regularelem, sizeof(RegularElement) * MaxNumRegularElems);

	}
	
	regularelem[cur_regelem_index].base[0] = position[0];
	regularelem[cur_regelem_index].base[1] = position[1];
	regularelem[cur_regelem_index].ID = NAMEOFREGELEM + cur_regelem_index;
	regularelem[cur_regelem_index].type = type;

	////Initialize the transformation parameters for this regular element
	regularelem[cur_regelem_index].transform_matrix.setIdentity();
	regularelem[cur_regelem_index].transposeRot.setIdentity();

	regularelem[cur_regelem_index].rotang = 0;
	regularelem[cur_regelem_index].s = 1;

	cur_regelem_index ++;

}

////The overload method for generating vector field through optimization 
void SetBasis(int triangle, CPoint point, int type)
{
	GLdouble position[3];
	ScreenToWorld(position, point);

	////Add to the regular elements list
	if(cur_regelem_index >= MaxNumRegularElems-1 )
	{
		//Allocate new space for the element list
		MaxNumRegularElems += 50;
		regularelem = (RegularElement*)realloc(regularelem, sizeof(RegularElement) * MaxNumRegularElems);

	}
	
	regularelem[cur_regelem_index].base[0] = position[0];
	regularelem[cur_regelem_index].base[1] = position[1];
	regularelem[cur_regelem_index].ID = NAMEOFREGELEM + cur_regelem_index;
	regularelem[cur_regelem_index].basis_triangle = triangle;  //store the triangle containing the basis
	regularelem[cur_regelem_index].type = type;

	////Initialize the transformation parameters for this regular element
	regularelem[cur_regelem_index].transform_matrix.setIdentity();
	regularelem[cur_regelem_index].transposeRot.setIdentity();

	regularelem[cur_regelem_index].rotang = 0;
	regularelem[cur_regelem_index].s = 1;

	cur_regelem_index ++;

}


void SetDirect(CPoint point)
{
	////Update the Cur_regularID here
	GLdouble position[3];
	ScreenToWorld(position, point);

	regularelem[cur_regelem_index - 1].Direct.entry[0] = position[0] - regularelem[cur_regelem_index - 1].base[0];
	regularelem[cur_regelem_index - 1].Direct.entry[1] = position[1] - regularelem[cur_regelem_index - 1].base[1];

	////Normalize the regular vector and set it to the default strength

	normalize(regularelem[cur_regelem_index - 1].Direct);

	regularelem[cur_regelem_index - 1].Direct.entry[0] *= RegularStrength;
	regularelem[cur_regelem_index - 1].Direct.entry[1] *= RegularStrength;
	
	//You might need to update the transform matrix for convergent and divergent elements
	if(regularelem[cur_regelem_index - 1].type > 0)
		SetTransform(regularelem[cur_regelem_index - 1].base,\
		     regularelem[cur_regelem_index - 1].Direct,\
		     regularelem[cur_regelem_index - 1].transform_matrix,\
		     regularelem[cur_regelem_index - 1].transposeRot);
}


////Implement the coordinates transformation from screen space to opengl space
void ScreenToWorld(double position[3], CPoint center)
{
    int viewport[4];
	GLdouble modelMatrix[16];
    GLdouble projMatrix[16];

	glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
    glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
    glGetIntegerv(GL_VIEWPORT,viewport);
    
	gluUnProject(
        center.x,
        center.y,
        0,
        modelMatrix,
        projMatrix,
        viewport,
        //the next 3 parameters are the pointers to the final object
        //coordinates. Notice that these MUST be double's
        &position[0], //-> pointer to your own position (optional)
        &position[1], // id
        &position[2] // id
    );

	//position[0] = position[0]-0.02;
	//position[1] = 1 - position[1] + 0.02;
	position[0] = position[0]-0.01;
	position[1] = 1 - position[1]+0.01;
}


////Initialize the Jacobian matrix for new added element according to its type
void InitialJacobian(icMatrix3x3 &Jacobian, int type)
{
	if(type == SOURCE)
	{
		Jacobian.set(1, 0, 0,   0, 1, 0,   0, 0, 1);
	}
	else if(type == SINK)
	{
		Jacobian.set(-1, 0, 0,   0, -1, 0,   0, 0, 1);
	}
	else if(type == SADDLE)
	{
		Jacobian.set(1, 0, 0,   0, -1, 0,   0, 0, 1);
	}
	else if(type == CWCENTER)
	{
		Jacobian.set(0, 1, 0,   -1, 0, 0,   0, 0, 1);
	}
	else{
		Jacobian.set(0, -1, 0,   1, 0, 0,   0, 0, 1);
	}
}

////Initialize the edit box for specific singular element
////according to the center of that element
////The coordinates of the 4 points will changed according to user operation
void InitEditBox(int index, double x, double y)
{
	////initial the original edit box
	singularelem[index].editbox.p1.entry[0] = x - EDITBOXSIZE;  ////low left point
	singularelem[index].editbox.p1.entry[1] = y - EDITBOXSIZE;

	singularelem[index].editbox.p2.entry[0] = x - EDITBOXSIZE;  ////upper left point
	singularelem[index].editbox.p2.entry[1] = y + EDITBOXSIZE;

	singularelem[index].editbox.p3.entry[0] = x + EDITBOXSIZE;  ////upper right point
	singularelem[index].editbox.p3.entry[1] = y + EDITBOXSIZE;

	singularelem[index].editbox.p4.entry[0] = x + EDITBOXSIZE;  ////low right point
	singularelem[index].editbox.p4.entry[1] = y - EDITBOXSIZE;

	singularelem[index].editbox.Up.entry[0] = x;  ////rotation controling point
	singularelem[index].editbox.Up.entry[1] = y + 2*EDITBOXSIZE;


	////initial current editbox
	singularelem[index].cur_editbox.p1.entry[0] = x - EDITBOXSIZE;  ////low left point
	singularelem[index].cur_editbox.p1.entry[1] = y - EDITBOXSIZE;

	singularelem[index].cur_editbox.p2.entry[0] = x - EDITBOXSIZE;  ////upper left point
	singularelem[index].cur_editbox.p2.entry[1] = y + EDITBOXSIZE;

	singularelem[index].cur_editbox.p3.entry[0] = x + EDITBOXSIZE;  ////upper right point
	singularelem[index].cur_editbox.p3.entry[1] = y + EDITBOXSIZE;

	singularelem[index].cur_editbox.p4.entry[0] = x + EDITBOXSIZE;  ////low right point
	singularelem[index].cur_editbox.p4.entry[1] = y - EDITBOXSIZE;
	
	singularelem[index].cur_editbox.Up.entry[0] = x;  ////rotation controling point
	singularelem[index].cur_editbox.Up.entry[1] = y + 2*EDITBOXSIZE;
}


////Set the transformation matrix for convergent/divergent elements
void SetTransform(
	double base[2],
	icVector2 Direct,
	icMatrix3x3 &transform_matrix, 
	icMatrix3x3 &transposeRot
)
{
	icMatrix3x3 trans, intermedia, inversetranslate;
	icVector2 tempd;
	tempd = Direct;

	////Normalize the direction vector
	normalize(tempd);

	double theta = atan2(tempd.entry[1], tempd.entry[0]);

	////Set the initialize transformation
	trans.setIdentity();

	trans.set(1, 0, base[0],
		      0, 1, base[1],
			  0, 0, 1);

	intermedia.set(cos(theta), -sin(theta), 0,
		           sin(theta),  cos(theta), 0,
				   0, 0, 1);

	inversetranslate.set(1, 0, -base[0],
		                 0, 1, -base[1],
						 0, 0, 1);

	trans.rightMultiply(intermedia);
	trans.rightMultiply(inversetranslate);

	transform_matrix.set(intermedia);

	////store the transpose of the rotation matrix
	transposeRot.set(cos(theta),  sin(theta), 0,
		             -sin(theta), cos(theta), 0,
				     0, 0, 1);

}


/*****************************************************************
Routine for vector field combination
*****************************************************************/
void getVector(double x, double y, icVector2 &vec, double &mag)
{
   int i;
   double  dx, dy, vx, vy, tvx, tvy, r=0.;

   icMatrix3x3 tempJacobian, transposerot;
   double ang;
   double newdirect[2] = {0, 0};

   double sim_magnitude = 0;  ////

   vx = vy = 0.;

   ///Combine all the singular elements 
   for(i = 0; i < MaxNumSingularElems; i++)
   {
	   if(singularelem[i].ID > 0 && !singularelem[i].deleted){
			dx = x - singularelem[i].centerx;
			dy = y - singularelem[i].centery;

			r  = dx*dx + dy*dy; 

    		if (r < DistanceThreshold)   r = DistanceThreshold;

			tempJacobian.set(singularelem[i].transform_matrix);
			tempJacobian.rightMultiply(singularelem[i].Jacobian);

			ang = singularelem[i].rotang;

			transposerot.set(cos(ang), sin(ang), 0,
			                 -sin(ang),  cos(ang), 0,
							 0,0,1);

			tempJacobian.rightMultiply(transposerot);

			tvx = (dx * tempJacobian.entry[0][0] + dy * tempJacobian.entry[0][1])/r;  
			tvy = (dx * tempJacobian.entry[1][0] + dy * tempJacobian.entry[1][1])/r;
			vx += tvx/sqrt(r);
			vy += tvy/sqrt(r);
			//tvx = (dx * tempJacobian.entry[0][0] + dy * tempJacobian.entry[0][1]);  
			//tvy = (dx * tempJacobian.entry[1][0] + dy * tempJacobian.entry[1][1]);
			//vx += tvx;
			//vy += tvy;


			if(singularelem[i].type == SOURCE)
			{
			    sim_magnitude += 1./sqrt(sqrt(r));
			}
			
			else if(singularelem[i].type == SINK)
			{
			    sim_magnitude -= 1./sqrt(sqrt(r));
			}

	   }
	
   }


   ////Here we need to combine the regular element into the whole field
   for(i = 0; i < MaxNumRegularElems; i++)
   {
	   if(regularelem[i].ID > 0)
	   {
			dx = x - regularelem[i].base[0];
			dy = y - regularelem[i].base[1];

			r  = dx*dx + dy*dy; 

			if (r < DistanceThreshold)   r = DistanceThreshold;

			if(regularelem[i].type == 0) ////regular element
			{
				newdirect[0] = regularelem[i].Direct.entry[0] * regularelem[i].transform_matrix.entry[0][0]\
					+ regularelem[i].Direct.entry[1] * regularelem[i].transform_matrix.entry[0][1];
				newdirect[1] = regularelem[i].Direct.entry[0] * regularelem[i].transform_matrix.entry[1][0]\
					+ regularelem[i].Direct.entry[1] * regularelem[i].transform_matrix.entry[1][1];

				tvx = newdirect[0]/r;
				tvy = newdirect[1]/r;

			}

			else{
				////Note that we haven't considered the user tranformations performed on this special element 07/07/05

				double t_vx = regularelem[i].transposeRot.entry[0][0] * dx + regularelem[i].transposeRot.entry[0][1] * dy;
				double t_vy = regularelem[i].transposeRot.entry[1][0] * dx + regularelem[i].transposeRot.entry[1][1] * dy;

				if(regularelem[i].type == 1) ////convergent element
				{
					//t_vx = 1./(40*r*sqrt(1+t_vy*t_vy));
					//t_vy = -t_vy/(r*sqrt(1+t_vy*t_vy));
					t_vx = 1./r;
					t_vy = -t_vy*8./r;
				}

				else{ ////divergent element
					//t_vx = -1./(40*r*sqrt(1+t_vy*t_vy));
					//t_vy = t_vy/(r*sqrt(1+t_vy*t_vy));
					t_vx = 1./r;
					t_vy = t_vy*8./r;
				}
				tvx = (regularelem[i].transform_matrix.entry[0][0] * t_vx + regularelem[i].transform_matrix.entry[0][1] * t_vy)/r;
				tvy = (regularelem[i].transform_matrix.entry[1][0] * t_vx + regularelem[i].transform_matrix.entry[1][1] * t_vy)/r;
			}

			vx += tvx;
			vy += tvy;
	   }
   }

   ////before normalizing the vector, we may need to store the current magnitude of the speed for color plots
   //mag = sqrt(vx*vx + vy*vy);
     mag = sim_magnitude;


   ////This "normalization" just for visual results, it will not affect the feature of the field
 //  r = vx*vx + vy*vy;
	//		
	//if (r > dmax*dmax) { 
	//	r  = sqrt(r); 
	//	vx *= dmax/r; 
	//	vy *= dmax/r; 
	//}

	////Maybe combine with other basic field here 07/06/05
	////For example the field read from the file that store the vector field per vertex

	vec.entry[0] = vx;
	vec.entry[1] = vy;
}

////Get the vector on each vertex
void CalVectors()
{
	int   i/*, j*/; 
	icVector2 vect;
	double mag;
	
	max_mag = 0; min_mag = 100;

    //Using triangle mesh to generate the vector field under global frame
	for (i = 0; i < Object.nverts; i++) {
		getVector(Object.vlist[i]->x, Object.vlist[i]->y, vect, mag);

		/***********************************************************/
		/*save the vector value before normalization 02/21/07 */
		Object.vlist[i]->vec_J = 0.0001*vect;  //07/19/06
	  
		double r = vect.entry[0]*vect.entry[0] + vect.entry[1]*vect.entry[1];
				
		if (r > dmax*dmax) { 
			r  = sqrt(r); 
			vect.entry[0] *= dmax/r; 
			vect.entry[1] *= dmax/r; 
		}

		Object.vlist[i]->vec = vect;
		/***********************************************************/

		Object.vlist[i]->mag_speed = mag;

		////store to the backup field variable
		backup_field[i] = vect;
		
		////find the maximun and minmun magnitude of speed
		if(Object.vlist[i]->mag_speed > max_mag) max_mag = Object.vlist[i]->mag_speed;
		if(Object.vlist[i]->mag_speed < min_mag) min_mag = Object.vlist[i]->mag_speed;
	}	
	

	////transfer the vector into local frame
	GetLocalVector();
}

/*
For vortex field 04/05/07
u(x, y, t) = -cos(k_x*x)sin(k_y*y)*exp(-u*(kx*kx+ky*ky)*t)
v(x, y, t) = -sin(k_x*x)cos(k_y*y)*exp(-u*(kx*kx+ky*ky)*t)
*/
void getformular_4(double x, double y, double t, icVector2 &vec, icVector2 &non_norm, double &p)
{
	double vx, vy, r = 0;
	double kx = 2*M_PI/0.5;
	double ky = 2*M_PI/0.5;
	double miu = 0.2;
	double coeff = miu*(kx*kx+ky*ky)*t;

	vx = -cos(kx*x)*sin(ky*y)*exp(-coeff);
	vy = sin(kx*x)*cos(ky*y)*exp(-coeff);
	
	//vx = -cos(kx*(x-t*0.3))*sin(ky*(y-t*1.2))*exp(-coeff);
	//vy = sin(kx*(x+t*0.5))*cos(ky*(y-t*1.))*exp(-coeff);
	
	//vx = -cos(kx*(x-t*0.3))*sin(ky*(y+t*0.6))*exp(-coeff);
	//vy = sin(kx*(x+t*0.5))*cos(ky*(y-t*1.))*exp(-coeff);

	//vx = -cos(kx*(x+t*0.3))*sin(ky*(y-t*0.6))*exp(-coeff);
	//vy = sin(kx*(x+t*0.5))*cos(ky*(y-t*1.))*exp(-coeff);

	//vx = -cos(kx*(x-t*0.3))*sin(ky*(y-t*0.6))*exp(-coeff);
	//vy = sin(kx*(x+t*0.8))*cos(ky*(y-t*0.5))*exp(-coeff);

	//vx = -cos(kx*(x-t*0.5))*sin(ky*(y+t*0.5))*exp(-coeff);
	//vy = sin(kx*(x-t*0.5))*cos(ky*(y-t*0.5))*exp(-coeff);

	///*we can apply the shearing to the flow according to the time*/
	//icMatrix2x2 tran;
	////double M[][2] = {20*t, 20*(1-t), 40*(1-t), 60*t};
	//double M[][2] = {20*t, 20*(1-t), 20*(1-t), 20*t};
	////double M[][2] = {20, 0, 0, 20};
	////double M[][2] = {0, 20, 20, 0};
	//tran.set(M);
	//double vx_ = tran.entry[0][0]*vx + tran.entry[0][1]*vy;
	//double vy_ = tran.entry[1][0]*vx + tran.entry[1][1]*vy;
	//
	//vx = vx_;
	//vy = vy_;

	non_norm.entry[0] = vx;
	non_norm.entry[1] = vy;
	
	////normalize

	r = vx*vx + vy*vy;
			
	if (r > dmax*dmax) { 
		r  = sqrt(r); 
		vx *= dmax/r; 
		vy *= dmax/r; 
	}
	vec.entry[0] = vx;
	vec.entry[1] = vy;

	/*calculate the P*/
	p = -.25*(cos(2*kx*x)+cos(2*ky*y))*exp(-2*miu*(kx*kx+ky*ky)*t);

}

/*
Calculate the gradient of the pressure
p(x, y, t) = -0.25*[cos(2k_x*x)+cos(2k_y*y)]*exp(-2u*(kx*kx+ky*ky)*t)

u(x, y, t) = 0.25*2*k_x*sin(2*k_x*x)*exp(-2u*(kx*kx+ky*ky)*t)
v(x, y, t) = 0.25*2*k_y*sin(2*k_y*y)*exp(-2u*(kx*kx+ky*ky)*t)
*/

void getformular_5(double x, double y, double t, icVector2 &vec, icVector2 &non_norm, double &p)
{
	double vx, vy, r = 0;
	double kx = 2*M_PI/0.5;
	double ky = 2*M_PI/0.5;
	double miu = 0.2;
	double coeff = 2*miu*(kx*kx+ky*ky)*t;

	vx = 0.25*2*kx*sin(2*kx*x)*exp(-coeff);
	vy = 0.25*2*ky*sin(2*ky*y)*exp(-coeff);

	/*we can apply the shearing to the flow according to the time*/
	//icMatrix2x2 tran;
	////double M[][2] = {20*t, 20*(1-t), 40*(1-t), 60*t};
	//double M[][2] = {20*t, 20*(1-t), 20*(1-t), 20*t};
	////double M[][2] = {20, 0, 0, 20};
	////double M[][2] = {0, 20, 20, 0};
	//tran.set(M);
	//double vx_ = tran.entry[0][0]*vx + tran.entry[0][1]*vy;
	//double vy_ = tran.entry[1][0]*vx + tran.entry[1][1]*vy;
	//
	//vx = vx_;
	//vy = vy_;

	non_norm.entry[0] = vx;
	non_norm.entry[1] = vy;
	
	////normalize

	r = vx*vx + vy*vy;
			
	if (r > dmax*dmax) { 
		r  = sqrt(r); 
		vx *= dmax/r; 
		vy *= dmax/r; 
	}
	vec.entry[0] = vx;
	vec.entry[1] = vy;

	/*calculate the P*/
	p = -.25*(cos(2*kx*x)+cos(2*ky*y))*exp(-2*miu*(kx*kx+ky*ky)*t);

}


/*Get the vector on each vertex
04/05/07
**/
void CalVectors_formula(double t)
{
	int   i/*, j*/; 
	icVector2 vect;
	double mag;
	
    //Using triangle mesh to generate the vector field under global frame
	for (i = 0; i < Object.nverts; i++) {
		//getVector_formula1(Object.vlist[i]->x, Object.vlist[i]->y, vect, mag);
		//getformular_4(Object.vlist[i]->x, Object.vlist[i]->y, t, vect, Object.vlist[i]->vec_J,
		//	Object.vlist[i]->mag_speed);
		getformular_5(Object.vlist[i]->x, Object.vlist[i]->y, t, vect, Object.vlist[i]->vec_J,
			Object.vlist[i]->mag_speed);
		
		////store to the backup field variable
		backup_field[i] = vect;
		Object.vlist[i]->vec = vect;
		Object.vlist[i]->vec_J = vect;
		
		////find the maximun and minmun magnitude of speed 04/10/07
		if(Object.vlist[i]->mag_speed > max_mag) max_mag = Object.vlist[i]->mag_speed;
		if(Object.vlist[i]->mag_speed < min_mag) min_mag = Object.vlist[i]->mag_speed;
	}	

	////transfer the vector into local frame
	NormalizeField();
	GetLocalVector();
}



/*
For the MCG paper! 07/05/07
We create new example for showing the uncertainty along the flow visualization pipeline
*/

/*
pitchfork bifurcation for periodic orbit, in polar coordinates
r' = r( (r-1)(k-(r-1)^2))
theta' = 1
*/
void formula1(double x, double y, icVector2 &vec, icVector2 &non_norm)
{
	double vx, vy, dx, dy;

	dx = 4*(x-0.501);
	dy = 4*(y-0.501);
	if(dx == 0 && dy == 0)
	{
		vx = 0;
		vy = 0;
		vec.set(0, 0);
		non_norm.set(0,0);
		return;
	}

	double cons_k = 0.05;
	double ss2 = dx*dx+dy*dy;
	double ss1 = sqrt(ss2);
	double A = (ss1-1)*(cons_k-ss2+2*ss1-1);

	vx = dx*A-dy;
	vy = dy*A+dx;

	vec.set(vx, vy);
	non_norm.set(vx, vy);
}

/*
saddle saddle connection
x' = -xy
y' = y(y^2-1)
*/
void formula2(double x, double y, icVector2 &vec, icVector2 &non_norm)
{
	double vx, vy, dx, dy;

	dx = 4.*(x-0.5001);
	dy = 4.*(y-0.5001);

	vx = -dx*dy;
	vy = /*dy**/(dy*dy-1);
	//vx = dy;
	//vy = 0.5*vx;
	vec.set(vx, vy);
	non_norm.set(vx, vy);
}

/*
existence or not existence of homoclinic orbit
x' = y
y' = ky - x(1-x)
*/
void formula3(double x, double y, icVector2 &vec, icVector2 &non_norm)
{
	double vx, vy, dx, dy;

	/*
	When k=0.01 and we zoom out the region 4 times, the separatrix can cross through periodic orbit!!
	07/05/2007
	*/
	//dx = 3*(x-0.53);
	//dy = 3*(y-0.51);
	dx = 3*(x-0.501);
	dy = 3*(y-0.501);

	double cons_k = 0.00250;

	vx = -dy;
	vy = cons_k*dy-dx*(1-dx);
	vec.set(vx, vy);
	non_norm.set(vx, vy);
}


/*
Monkey saddle
x'=x^2-y^2
y'=-2xy
*/

void formula4(double x, double y, icVector2 &vec, icVector2 &non_norm)
{
	double vx, vy, dx, dy;

	//dx = 3*(x-0.50);
	//dy = 3*(y-0.50);
	dx = 3*(x-0.501);
	dy = 3*(y-0.501);


	vx = dx*dx-dy*dy;
	vy = -2*dx*dy;
	vec.set(vx, vy);
	non_norm.set(vx, vy);
}


/*
Monkey saddle
x'=x^2-y^2
y'=2xy
*/

void formula5(double x, double y, icVector2 &vec, icVector2 &non_norm)
{
	double vx, vy, dx, dy;


	//dx = 3*(x-0.50);
	//dy = 3*(y-0.50);
	dx = 3*(x-0.501);
	dy = 3*(y-0.501);


	vx = dx*dx-dy*dy;
	vy = 2*dx*dy;
	vec.set(vx, vy);
	non_norm.set(vx, vy);
}


void create_ex1()
{
	int   i/*, j*/; 
	icVector2 vect;
	double max_mag = 10000;
	
    //Using triangle mesh to generate the vector field under global frame
	for (i = 0; i < Object.nverts; i++) {
		
		formula2(Object.vlist[i]->x, Object.vlist[i]->y, vect, vect);
		
		////store to the backup field variable
		//backup_field[i] = vect;
		Object.vlist[i]->vec = vect;
		Object.vlist[i]->vec_J = vect;
		
		//if(max_mag < length(vect))
		//	max_mag = length(vect);
		if(max_mag > length(vect))
			max_mag = length(vect);
	}	

	////transfer the vector into local frame
	
	NormalizeField();
	
	/*scale the vector field according to the maximum magnitude
	07/19/07*/
	for(i=0; i<Object.nverts; i++)
	{
		Object.vlist[i]->vec_J.entry[0] = Object.vlist[i]->vec_J.entry[0]/max_mag*0.06/**DistanceThreshold*/;
		Object.vlist[i]->vec_J.entry[1] = Object.vlist[i]->vec_J.entry[1]/max_mag*0.06/**DistanceThreshold*/;
	}

	GetLocalVector();
}



/*
We perturbate the vector values with certain percentage
07/09/07
*/
void perturbate_vecs()
{
	int i;
	Vertex *v;
	double vx, vy;
	icVector2 perturb;
	double theta;
	int sign;
	double percentage = 0.3;

	for(i=0; i<Object.nverts; i++)
	{
		v = Object.vlist[i];

		/*randomly generate two values*/
		//vx = (double)rand()/RAND_MAX;
		//vy = (double)rand()/RAND_MAX;

		///*perturbation scheme 1: add a random vector value*/
		//perturb.set(vx, vy);
		//normalize(perturb);
		//perturb = length(v->vec)*perturb;
		//perturb = percentage*perturb;

		//v->vec = v->vec + perturb;

		/*perturbation scheme 2: add a random rotation*/
		if(rand()%2 == 0)
			sign = 1;
		else
			sign = -1;

		theta = (double)rand()/(10*RAND_MAX);

		vx = v->vec.entry[0]*cos(theta)-v->vec.entry[1]*sin(theta);
		vy = v->vec.entry[0]*sin(theta)+v->vec.entry[1]*cos(theta);

		v->vec.set(vx, vy);
	}

	////transfer the vector into local frame
	NormalizeField();
	GetLocalVector();
}







/***************************************************************
Get the vectors on vertices in local frame
This routine must be called after building the local frame
and after get the global vector field on vertices!!!!
***************************************************************/

void GetLocalVector()
{
	Vertex *vert;
	Face *face;

    int i, j;

	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];

			/*Use normalized vector field*/
			if(UseNormalizedVF)
				GlobalToLocal(vert->vec.entry, face->direct_vec[j], i);
			/*Use non-normalized vector field*/
			else
				GlobalToLocal(vert->vec_J.entry, face->direct_vec[j], i);
		}
	}
}

/*******************************************************************
This routine is used to transfer the vector on global frame to
specific local frame
*******************************************************************/
void GlobalToLocal(double v[2], icVector2 &iv, int face_id)
{
	if(face_id < 0)
		return;

	Face *face = Object.flist[face_id];

	icVector2 lx, ly;

	lx = face->LX;
	ly = face->LY;

	normalize(lx);
	normalize(ly);


	iv.entry[0] = lx.entry[0]*v[0] + lx.entry[1]*v[1];
	iv.entry[1] = ly.entry[0]*v[0] + ly.entry[1]*v[1];

}


void InitUnderneathMesh()
{
	//double zerovec[2] = {0.};

	for(int i = 0; i < Object.nverts; i++){
		Object.vlist[i]->InRegion = 0;
		Object.vlist[i]->OnBoundary = 0;

		//Object.vlist[i]->vec.set(zerovec);
	}

	for(int i = 0; i < Object.nfaces; i++)
	{
		//Object.flist[i]->direct_vec[0].set(zerovec);
		//Object.flist[i]->direct_vec[1].set(zerovec);
		//Object.flist[i]->direct_vec[2].set(zerovec);

		Object.flist[i]->inDesignCellCycle = 0;
	}

}

void ClearAllField()
{
	int i = 0;

	////Clear all the elements list, singularities' list, trajectories' list
    cur_singelem_index = 0;
	for(i = 0; i < MaxNumSingularElems; i++)
	{
		singularelem[i].ID = -1;
		singularelem[i].type = -1;
	}
    
    cur_regelem_index = 0;
	for(i = 0; i < MaxNumRegularElems; i++)
	{
		regularelem[i].ID = -1;
		regularelem[i].type = -1;
	}
    
	////Singularities
    cur_singularity_index = 0;

	////Trajectories and separatrices
    cur_traj_index = 0;
	cur_separatrices_index = 0;

	for(i = 0; i < MaxNumTrajectories; i ++)
	{
		num_linesegs_curtraj[i] = 0;
	}

	////Reset all the flag of the mesh object
	InitUnderneathMesh();
	
	////Reset the variables of limit cycle detection
    InitLimitCycleDetection();

	InitLimitCycleStructure();

	////Reset the variables of region smoothing
    InitRegionSmooth();

	////Reset the field stored on each vertex
	CalVectors();

	////Reset rotation and reflection parameters here
//	RotateDegreeofField = 0;

	////Reset all the related flag here
}

void UpdateTransform(int TransformType, int choose_ID)
{
	icMatrix3x3 trans, intermedia, inversetranslate;

	////Set the transformation
	////Suppose that we perform scale first

	if(TransformType > 0)
	{
		////Update the Jacobian matrix of selected element
		if(choose_ID > 0 && choose_ID < NAMEOFREGELEM)
		{
             UpdateSinElemTransform(TransformType, choose_ID - 1);
		}
	
		else if(choose_ID >= NAMEOFREGELEM && choose_ID < NAMEOFSINGCONTROL )
		{ ////we choose a regular, convergent or divergent element
             UpdateRegElemTransform(TransformType, choose_ID - NAMEOFREGELEM);
		}

	}
}

////Update the transform matrix of sepecific singular element
void UpdateSinElemTransform(int TransformType, int elem_id)
{
	icMatrix3x3 trans, intermedia, inversetranslate;

	////Note that we need to translate to the center of the element
	////then perform the specified transformation
	////and translate back
	singularelem[elem_id].transform_matrix.set(\
		1, 0, singularelem[elem_id].centerx,\
		0, 1, singularelem[elem_id].centery,\
		0, 0, 1);

	inversetranslate.set(1, 0, -singularelem[elem_id].centerx,
		0, 1, -singularelem[elem_id].centery,
		0, 0, 1);

	if(TransformType == 1) ////If we scale the element uniformly
	{
		singularelem[elem_id].s += uniforms;

		////if the scalar less than 0, it will change the type of singular element!!!
		if(singularelem[elem_id].s <= 0)
			singularelem[elem_id].s = 0.1;
	}

	else if(TransformType == 2)
	{
		singularelem[elem_id].sx += sx;
		singularelem[elem_id].sy += sy;

		////if the scalar less than 0, it will change the type of singular element!!!
		if(singularelem[elem_id].sx <= 0)
			singularelem[elem_id].sx = 0.1;
		if(singularelem[elem_id].sy <= 0)
			singularelem[elem_id].sy = 0.1;
	}
	else{  ////it is rotation
		singularelem[elem_id].rotang += rotateAng;
	}
		

	////set the rotation matrix
	intermedia.set(cos(singularelem[elem_id].rotang), -sin(singularelem[elem_id].rotang), 0,
					sin(singularelem[elem_id].rotang), cos(singularelem[elem_id].rotang),  0,
					0,0,1);

	singularelem[elem_id].transform_matrix.rightMultiply(intermedia);

	////set the uniform scaling matrix
	intermedia.set(singularelem[elem_id].s, 0, 0,
					0, singularelem[elem_id].s, 0,
					0, 0, 1);
	singularelem[elem_id].transform_matrix.rightMultiply(intermedia);

	////set the non uniform scaling matrix
	intermedia.set(singularelem[elem_id].sx, 0, 0,
					0, singularelem[elem_id].sy, 0,
					0, 0, 1);

	singularelem[elem_id].transform_matrix.rightMultiply(intermedia);
	singularelem[elem_id].transform_matrix.rightMultiply(inversetranslate);

	UpdateEditBox(elem_id);
}


////Update the transform matrix of specific regular element
void UpdateRegElemTransform(int TransformType, int elem_id)
{
	icMatrix3x3 trans, intermedia, inversetranslate;

	regularelem[elem_id].transform_matrix.set(\
		1, 0, regularelem[elem_id].base[0],\
		0, 1, regularelem[elem_id].base[1],\
		0, 0, 1);

	inversetranslate.set(1, 0, -regularelem[elem_id].base[0],
		0, 1, -regularelem[elem_id].base[1],
		0, 0, 1);

	if(TransformType == 1) ////If we scale the field
	{
		regularelem[elem_id].s += uniforms;

		if(regularelem[elem_id].s < 0)
			regularelem[elem_id].s = 0.1;
	}
	else{  ////it is rotation
		regularelem[elem_id].rotang += rotateAng;
	}
		
	////set the rotation matrix
	intermedia.set(cos(regularelem[elem_id].rotang), -sin(regularelem[elem_id].rotang), 0,
					sin(regularelem[elem_id].rotang), cos(regularelem[elem_id].rotang),  0,
					0,0,1);

	regularelem[elem_id].transform_matrix.rightMultiply(intermedia);
	
	////set the scale matrix
	intermedia.set(regularelem[elem_id].s, 0, 0,
					0, regularelem[elem_id].s, 0,
					0, 0, 1);

	regularelem[elem_id].transform_matrix.rightMultiply(intermedia);

	regularelem[elem_id].transform_matrix.rightMultiply(inversetranslate);
}


void UpdateEditBox(int elem_id)
{
	double origin[3];

	////Update p1
	origin[0] = singularelem[elem_id].editbox.p1.entry[0];
	origin[1] = singularelem[elem_id].editbox.p1.entry[1];
	origin[2] = 1;

	singularelem[elem_id].cur_editbox.p1.entry[0] = singularelem[elem_id].transform_matrix.entry[0][0] * origin[0]\
		+ singularelem[elem_id].transform_matrix.entry[0][1] * origin[1]\
		+ singularelem[elem_id].transform_matrix.entry[0][2] * origin[2];
	
	singularelem[elem_id].cur_editbox.p1.entry[1] = singularelem[elem_id].transform_matrix.entry[1][0] * origin[0]\
		+ singularelem[elem_id].transform_matrix.entry[1][1] * origin[1]\
		+ singularelem[elem_id].transform_matrix.entry[1][2] * origin[2];

	////Update p2
	origin[0] = singularelem[elem_id].editbox.p2.entry[0];
	origin[1] = singularelem[elem_id].editbox.p2.entry[1];
	origin[2] = 1;

	singularelem[elem_id].cur_editbox.p2.entry[0] = singularelem[elem_id].transform_matrix.entry[0][0] * origin[0]\
		+ singularelem[elem_id].transform_matrix.entry[0][1] * origin[1]\
		+ singularelem[elem_id].transform_matrix.entry[0][2] * origin[2];
	
	singularelem[elem_id].cur_editbox.p2.entry[1] = singularelem[elem_id].transform_matrix.entry[1][0] * origin[0]\
		+ singularelem[elem_id].transform_matrix.entry[1][1] * origin[1]\
		+ singularelem[elem_id].transform_matrix.entry[1][2] * origin[2];

	////Update p3
	origin[0] = singularelem[elem_id].editbox.p3.entry[0];
	origin[1] = singularelem[elem_id].editbox.p3.entry[1];
	origin[2] = 1;

	singularelem[elem_id].cur_editbox.p3.entry[0] = singularelem[elem_id].transform_matrix.entry[0][0] * origin[0]\
		+ singularelem[elem_id].transform_matrix.entry[0][1] * origin[1]\
		+ singularelem[elem_id].transform_matrix.entry[0][2] * origin[2];
	
	singularelem[elem_id].cur_editbox.p3.entry[1] = singularelem[elem_id].transform_matrix.entry[1][0] * origin[0]\
		+ singularelem[elem_id].transform_matrix.entry[1][1] * origin[1]\
		+ singularelem[elem_id].transform_matrix.entry[1][2] * origin[2];
	
	////Update p4
	origin[0] = singularelem[elem_id].editbox.p4.entry[0];
	origin[1] = singularelem[elem_id].editbox.p4.entry[1];
	origin[2] = 1;

	singularelem[elem_id].cur_editbox.p4.entry[0] = singularelem[elem_id].transform_matrix.entry[0][0] * origin[0]\
		+ singularelem[elem_id].transform_matrix.entry[0][1] * origin[1]\
		+ singularelem[elem_id].transform_matrix.entry[0][2] * origin[2];
	
	singularelem[elem_id].cur_editbox.p4.entry[1] = singularelem[elem_id].transform_matrix.entry[1][0] * origin[0]\
		+ singularelem[elem_id].transform_matrix.entry[1][1] * origin[1]\
		+ singularelem[elem_id].transform_matrix.entry[1][2] * origin[2];

	////Update Up
	origin[0] = singularelem[elem_id].editbox.Up.entry[0];
	origin[1] = singularelem[elem_id].editbox.Up.entry[1];
	origin[2] = 1;

	singularelem[elem_id].cur_editbox.Up.entry[0] = singularelem[elem_id].transform_matrix.entry[0][0] * origin[0]\
		+ singularelem[elem_id].transform_matrix.entry[0][1] * origin[1]\
		+ singularelem[elem_id].transform_matrix.entry[0][2] * origin[2];
	
	singularelem[elem_id].cur_editbox.Up.entry[1] = singularelem[elem_id].transform_matrix.entry[1][0] * origin[0]\
		+ singularelem[elem_id].transform_matrix.entry[1][1] * origin[1]\
		+ singularelem[elem_id].transform_matrix.entry[1][2] * origin[2];
}

/*------------------------------------------------------------------------------------------*/
////Saved and reload current field on each vertex into / from a file 

void SaveFieldPerVertex(const char *filename)
{
	//FILE *vfp = fopen("savedfieldpervertex.txt", "w");
	FILE *vfp = fopen(filename, "w");
	if(vfp == NULL)
	{
		MessageBox(NULL, "Can't open the file!", "Error", MB_OK);
		return;
	}

	for(int i = 0; i < Object.nverts; i++)
	{
		fprintf(vfp,"%f,%f\n", Object.vlist[i]->vec.entry[0], Object.vlist[i]->vec.entry[1]);
	}
	fclose(vfp);
}

void ReloadFieldPerVertex(const char *filename)
{
    float vx, vy;

	//FILE *vfp = fopen("savedfieldpervertex.txt","r");
	FILE *vfp = fopen(filename,"r");

	if(vfp == NULL)
	{
		MessageBox(NULL, "Can't open the file!", "Error", MB_OK);
		return;
	}

	for(int i = 0; i < Object.nverts; i++)
	{
		fscanf(vfp, "%f,%f\n", &vx, &vy);
		Object.vlist[i]->vec.entry[0] = vx;
		Object.vlist[i]->vec.entry[1] = vy;

		/*save the vector before normalization 02/21/07*/
		Object.vlist[i]->vec_J = Object.vlist[i]->vec;
	}
	fclose(vfp);
}


////////
////save the field through saving its elements
void SaveFieldElements(const char *filename)
{
	int i;
	//FILE *vfp = fopen("savedfieldelements.txt", "w");
		
	FILE *vfp = fopen(filename, "w");

	if(vfp == NULL)
	{
		MessageBox(NULL, "Can't open the file!", "Error", MB_OK);
		return;
	}

	////write singular elements first
	fprintf(vfp, "SingularElem:%d\n", cur_singelem_index);
	for(i = 0; i < cur_singelem_index; i++)
	{
		WriteSingularElem(vfp, i);
	}

	////write regular elements then
	fprintf(vfp, "RegularElem:%d\n", cur_regelem_index);
	for(i = 0; i < cur_regelem_index; i++)
	{
		WriteRegularElem(vfp, i);
	}

	fclose(vfp);
}


////save to a temporary file before doing other advanced calculation
void  SaveBeforeCal()
{
	int i;
	FILE *vfp = fopen("temp_savedfield.txt", "w");

	////write singular elements first
	fprintf(vfp, "SingularElem:%d\n", cur_singelem_index);
	for(i = 0; i < cur_singelem_index; i++)
	{
		WriteSingularElem(vfp, i);
	}

	////write regular elements then
	fprintf(vfp, "RegularElem:%d\n", cur_regelem_index);
	for(i = 0; i < cur_regelem_index; i++)
	{
		WriteRegularElem(vfp, i);
	}

	fclose(vfp);
}


void ReloadFieldElements(const char *filename)
{
	int i;
	int StoreNumSingularElems, StoreNumRegularElems;
	//FILE *vfp = fopen("savedfieldelements.txt", "r");
	FILE *vfp = fopen(filename, "r");

	if(vfp == NULL)
	{
		MessageBox(NULL, "Can't open the file!", "Error", MB_OK);
		return;
	}

	////write singular elements first
	fscanf(vfp, "SingularElem:%d\n", &StoreNumSingularElems);
	for(i = 0; i < StoreNumSingularElems; i++)
	{
		ReadSingularElem(vfp, i);
	}
	cur_singelem_index = StoreNumSingularElems;

	////write regular elements then
	fscanf(vfp, "RegularElem:%d\n", &StoreNumRegularElems);
	for(i = 0; i < StoreNumRegularElems; i++)
	{
		ReadRegularElem(vfp, i);
	}
	cur_regelem_index = StoreNumRegularElems;

	fclose(vfp);
}


////Write or Read a single element according to the input index

void WriteSingularElem(FILE *fp, int index)
{
	fprintf(fp, "%d\n", singularelem[index].ID);
	fprintf(fp, "%f,%f\n", singularelem[index].centerx, singularelem[index].centery);
	fprintf(fp, "%d\n", singularelem[index].type);
	
	////write transformation matrix
	fprintf(fp, "%f,%f,%f\n", singularelem[index].transform_matrix.entry[0][0],\
		 singularelem[index].transform_matrix.entry[0][1], singularelem[index].transform_matrix.entry[0][2]);
	fprintf(fp, "%f,%f,%f\n", singularelem[index].transform_matrix.entry[1][0],\
		 singularelem[index].transform_matrix.entry[1][1], singularelem[index].transform_matrix.entry[1][2]);
	fprintf(fp, "%f,%f,%f\n", singularelem[index].transform_matrix.entry[2][0],\
		 singularelem[index].transform_matrix.entry[2][1], singularelem[index].transform_matrix.entry[2][2]);
	
	////write Jacobian matrix
	fprintf(fp, "%f,%f,%f\n", singularelem[index].Jacobian.entry[0][0],\
		 singularelem[index].Jacobian.entry[0][1], singularelem[index].Jacobian.entry[0][2]);
	fprintf(fp, "%f,%f,%f\n", singularelem[index].Jacobian.entry[1][0],\
		 singularelem[index].Jacobian.entry[1][1], singularelem[index].Jacobian.entry[1][2]);
	fprintf(fp, "%f,%f,%f\n", singularelem[index].Jacobian.entry[2][0],\
		 singularelem[index].Jacobian.entry[2][1], singularelem[index].Jacobian.entry[2][2]);

	////we need not write the edit box, it can be reconstructed after reading the previous information
}

void ReadSingularElem(FILE *fp, int index)
{
	float cx, cy, col0, col1, col2;
	fscanf(fp, "%d\n", &singularelem[index].ID);
    fscanf(fp, "%f,%f\n", &cx, &cy);
	singularelem[index].centerx = cx;
	singularelem[index].centery = cy;
	fscanf(fp, "%d\n", &singularelem[index].type);
	
	////read transformation matrix
	fscanf(fp, "%f,%f,%f\n", &col0, &col1, &col2);
	singularelem[index].transform_matrix.entry[0][0] = col0;
	singularelem[index].transform_matrix.entry[0][1] = col1;
	singularelem[index].transform_matrix.entry[0][2] = col2;

	fscanf(fp, "%f,%f,%f\n", &col0, &col1, &col2);
	singularelem[index].transform_matrix.entry[1][0] = col0;
	singularelem[index].transform_matrix.entry[1][1] = col1;
	singularelem[index].transform_matrix.entry[1][2] = col2;

	fscanf(fp, "%f,%f,%f\n", &col0, &col1, &col2);
	singularelem[index].transform_matrix.entry[2][0] = col0;
	singularelem[index].transform_matrix.entry[2][1] = col1;
	singularelem[index].transform_matrix.entry[2][2] = col2;
	
	////read Jacobian matrix
	fscanf(fp, "%f,%f,%f\n", &col0, &col1, &col2);
	singularelem[index].Jacobian.entry[0][0] = col0;
	singularelem[index].Jacobian.entry[0][1] = col1;
	singularelem[index].Jacobian.entry[0][2] = col2;

	fscanf(fp, "%f,%f,%f\n", &col0, &col1, &col2);
	singularelem[index].Jacobian.entry[1][0] = col0;
	singularelem[index].Jacobian.entry[1][1] = col1;
	singularelem[index].Jacobian.entry[1][2] = col2;

	fscanf(fp, "%f,%f,%f\n", &col0, &col1, &col2);
	singularelem[index].Jacobian.entry[2][0] = col0;
	singularelem[index].Jacobian.entry[2][1] = col1;
	singularelem[index].Jacobian.entry[2][2] = col2;

	////Initialize the transform matrix and other variables
	//singularelem[index].transform_matrix.setIdentity();
	//singularelem[index].s = singularelem[index].sx = singularelem[index].sy = 1;
	singularelem[index].rotang = 0;

	////Initialize the edit box for this singular element
	InitEditBox(index, cx, cy);
	
	//singularelem[index].Triangle_ID = TriangleDetect(cx, cy);
}


void WriteRegularElem(FILE *fp, int index)
{
	fprintf(fp, "%d\n", regularelem[index].ID);
	fprintf(fp, "%f,%f\n", regularelem[index].base[0], regularelem[index].base[1]);
	fprintf(fp, "%f,%f\n", regularelem[index].Direct.entry[0], regularelem[index].Direct.entry[1]);
	fprintf(fp, "%d\n", regularelem[index].type);
	
	////write transformation matrix
	fprintf(fp, "%f,%f,%f\n", regularelem[index].transform_matrix.entry[0][0],\
		 regularelem[index].transform_matrix.entry[0][1], regularelem[index].transform_matrix.entry[0][2]);
	fprintf(fp, "%f,%f,%f\n", regularelem[index].transform_matrix.entry[1][0],\
		 regularelem[index].transform_matrix.entry[1][1], regularelem[index].transform_matrix.entry[1][2]);
	fprintf(fp, "%f,%f,%f\n", regularelem[index].transform_matrix.entry[2][0],\
		 regularelem[index].transform_matrix.entry[2][1], regularelem[index].transform_matrix.entry[2][2]);
	
	////write transpose matrix
	fprintf(fp, "%f,%f,%f\n", regularelem[index].transposeRot.entry[0][0],\
		 regularelem[index].transposeRot.entry[0][1], regularelem[index].transposeRot.entry[0][2]);
	fprintf(fp, "%f,%f,%f\n", regularelem[index].transposeRot.entry[1][0],\
		 regularelem[index].transposeRot.entry[1][1], regularelem[index].transposeRot.entry[1][2]);
	fprintf(fp, "%f,%f,%f\n", regularelem[index].transposeRot.entry[2][0],\
		 regularelem[index].transposeRot.entry[2][1], regularelem[index].transposeRot.entry[2][2]);
}


void ReadRegularElem(FILE *fp, int index)
{
	float x, y, col0, col1, col2;
	fscanf(fp, "%d\n", &regularelem[index].ID);
	fscanf(fp, "%f,%f\n", &x, &y);
	regularelem[index].base[0] = x;
	regularelem[index].base[1] = y;
	fscanf(fp, "%f,%f\n", &x, &y);
	regularelem[index].Direct.entry[0] = x;
	regularelem[index].Direct.entry[1] = y;
	fscanf(fp, "%d\n", &regularelem[index].type);
	
	////read transformation matrix
	fscanf(fp, "%f,%f,%f\n", &col0, &col1, &col2);
	regularelem[index].transform_matrix.entry[0][0] = col0;
	regularelem[index].transform_matrix.entry[0][1] = col1;
	regularelem[index].transform_matrix.entry[0][2] = col2;

	fscanf(fp, "%f,%f,%f\n", &col0, &col1, &col2);
	regularelem[index].transform_matrix.entry[1][0] = col0;
	regularelem[index].transform_matrix.entry[1][1] = col1;
	regularelem[index].transform_matrix.entry[1][2] = col2;

	fscanf(fp, "%f,%f,%f\n", &col0, &col1, &col2);
	regularelem[index].transform_matrix.entry[2][0] = col0;
	regularelem[index].transform_matrix.entry[2][1] = col1;
	regularelem[index].transform_matrix.entry[2][2] = col2;

	
	////Read transpose rotation matrix
	fscanf(fp, "%f,%f,%f\n", &col0, &col1, &col2);
	regularelem[index].transposeRot.entry[0][0] = col0;
	regularelem[index].transposeRot.entry[0][1] = col1;
	regularelem[index].transposeRot.entry[0][2] = col2;


	fscanf(fp, "%f,%f,%f\n", &col0, &col1, &col2);
	regularelem[index].transposeRot.entry[1][0] = col0;
	regularelem[index].transposeRot.entry[1][1] = col1;
	regularelem[index].transposeRot.entry[1][2] = col2;


	fscanf(fp, "%f,%f,%f\n", &col0, &col1, &col2);
	regularelem[index].transposeRot.entry[2][0] = col0;
	regularelem[index].transposeRot.entry[2][1] = col1;
	regularelem[index].transposeRot.entry[2][2] = col2;

	regularelem[index].s = 1;
	regularelem[index].rotang = 0;

}


////Rotate or reflect the whole field
void RotateTheWholeField(double rotAng)
{
	icMatrix3x3 rotMatrix;
	int   i; 
	//Face *face;
	//int *verts;
	double s[3] = {0.};

	rotMatrix.set(cos(rotAng), -sin(rotAng), 0,\
		sin(rotAng), cos(rotAng), 0,\
		0, 0, 1);

	for(i = 0; i < Object.nverts; i++)
	{
		s[0] = Object.vlist[i]->vec.entry[0];
		s[1] = Object.vlist[i]->vec.entry[1];


		Object.vlist[i]->vec.entry[0] = rotMatrix.entry[0][0] * s[0] \
			+ rotMatrix.entry[0][1] * s[1];
		Object.vlist[i]->vec.entry[1] = rotMatrix.entry[1][0] * s[0] \
			+ rotMatrix.entry[1][1] * s[1];
	}

    GetLocalVector();
}


void ReflectTheWholeField(double rotAng)
{
	icMatrix3x3 rotMatrix;
	int   i; 
	//Face *face;
	//int *verts;
	double s[3] = {0.};

	rotMatrix.set(cos(rotAng), -sin(rotAng), 0,\
		-sin(rotAng), -cos(rotAng), 0,\
		0, 0, 1);


	for(i = 0; i < Object.nverts; i++)
	{
		s[0] = Object.vlist[i]->vec.entry[0];
		s[1] = Object.vlist[i]->vec.entry[1];


		Object.vlist[i]->vec.entry[0] = rotMatrix.entry[0][0] * s[0] \
			+ rotMatrix.entry[0][1] * s[1];
		Object.vlist[i]->vec.entry[1] = rotMatrix.entry[1][0] * s[0] \
			+ rotMatrix.entry[1][1] * s[1];
	}

    GetLocalVector();
}





/* The following routine implements the Jacobian based design 
   Note that we consider only one singular design element in the field only
   And we use global coordinates here.
   For surfaces, we can calculate the three vector values in the local frame of the triangle
   containing the design element, then perform constrained optimization to propagate
   the vector field to the whole surface
*/

void create_Jac_Field_Single(double a, double b, double c, 
							 double d, double x0, double y0)
{
	/* Here, we calculate the vector values on each vertices */
	int i;
	Vertex *v;

	for(i = 0; i < Object.nverts; i++)
	{
		v = Object.vlist[i];

		v->vec.entry[0] = a*(v->x - x0) + b*(v->y - y0);
		v->vec.entry[1] = c*(v->x - x0) + d*(v->y - y0);

		v->vec_J = v->vec;   /*we save the original vector before normalization*/
	}

	NormalizeField();  /*normalize the obtained vector field */
	GetLocalVector();
}



void TestRoutine()
{
	Face *face = Object.flist[369];
	Vertex *vert;
	icVector2 globalv;
	int i;
    
	////New method to set the three vectors
	face->direct_vec[0].entry[0]=
		face->xy[0][0] - (face->xy[1][0]+face->xy[2][0])/2;
	face->direct_vec[0].entry[1]=
		face->xy[0][1] - (face->xy[1][1]+face->xy[2][1])/2;
	normalize(face->direct_vec[0]);
	
	face->direct_vec[1].entry[0]=
		face->xy[1][0] - (face->xy[0][0]+face->xy[2][0])/2;
	face->direct_vec[1].entry[1]=
		face->xy[1][1] - (face->xy[0][1]+face->xy[2][1])/2;
	normalize(face->direct_vec[1]);

	face->direct_vec[2].entry[0]=
		face->xy[2][0] - (face->xy[1][0]+face->xy[0][0])/2;
	face->direct_vec[2].entry[1]=
		face->xy[2][1] - (face->xy[1][1]+face->xy[0][1])/2;
	normalize(face->direct_vec[2]);

	for(i = 0; i < 3; i++)
	{
		globalv = face->direct_vec[i].entry[0]*face->LX + face->direct_vec[i].entry[1]*face->LY;
		Object.vlist[face->verts[i]]->vec = 0.12 * globalv;
	}

	for(i = 0; i < face->nverts; i++)
	{
		vert = Object.vlist[face->verts[i]];
		vert->OnBoundary = 1;
	}

	BuildSmoothRegion();

	RegionSmooth();

	GetLocalVector();
	//CalVectors();
}


