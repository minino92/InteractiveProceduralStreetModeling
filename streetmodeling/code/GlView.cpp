// GlView.cpp : implementation file
//

#include "stdafx.h"
#include "IBFV2D.h"
#include "GlView.h"
#include ".\glview.h"
#include <math.h>
#include "VFDataStructure.h"

#include "RegionSmoothing.h"
#include "topologyedit.h"
#include "VFSynthesis.h"
#include "LimitCycleCreator.h"

#include "ClosedStreamlineTracing.h"
#include "NewLimitCycleDetect.h"

#include "LocalTracing.h"

#include "FindSCC.h"

#include "Taumap.h"

#include "scalardesign.h"

/*------------------------------------------------------------*/
//some constants for ibfv, some may be useless in future!!!06/23/05
#define	NPN 64
#define NMESH  100
//#define DM  ((double) (1.0/(NMESH-1.0)))
#define NPIX  512
#define SCALE 1.5

extern int REALWINSIZE;
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
//some global variables for ibfv
int     iframe = 0; 
int     Npat   = 32;
int     alpha  = (0.08*255);
double  tmax   = NPIX/(SCALE*NPN);
double  dmax   = SCALE/REALWINSIZE/*NPIX*/;
/*------------------------------------------------------------*/

/////////
int choose_ID, chosen_tenelem_ID = -1;
int SelectTriangleID;
int Separatrix_saddle;  //for separatrix modification

////Transformation variables for element edit
//double rotmatrix[3][3];
double sx, sy, uniforms;
double rotateAng;

////Global variables for whole field rotation and reflection
double RotateDegreeofField;

//get enhanced IBFV using two images
//GLuint Textures[2];
GLubyte f_tex[NPIX][NPIX][3], b_tex[NPIX][NPIX][3], applied_tex[NPIX][NPIX][3];


extern int MaxNumSingularElems;                 //Maximum number of singular elements
extern int MaxNumRegularElems;                  //Maximum number of regular elements
extern int MaxNumSingularities;                 //Maximum number of being captured singularities
extern int MaxNumTrajectories;                  //Maximum number of possible trajectories
                                                //(it should be flexible for future pen-and-ink sketch)
extern int MaxNumSeparatrices;                   //Maximum number of group of separatrices
extern int MaxNumLinesegsPerTraj;               //Maximum number of line segments for each trajectory

extern SingularElement *singularelem;          //Singular elememts' list
extern int cur_singelem_index;
extern RegularElement *regularelem;            //regular elememts' list
extern int cur_regelem_index;
extern Singularities *singularities;           //being captured singularites' list
extern int cur_singularity_index;
//extern Trajectory *trajectories2;
extern LineSeg **trajectories;                 //trajectories' list
extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory
extern Separatrices *separatrices;             //array for group of separatrices
extern int cur_separatrices_index;
extern LimitCycle *limitcycles;                 //limit cycles data structure
extern int cur_limitcycle_index;
extern int MaxNumLimitCycles;

extern Polygon3D Object;

////Global variables for limit cycle detection
extern int **CellCycleList;
extern int Cur_CellCycleIndex;                    //current cell cycle index
extern int *NumTriangleInEachCellCycle;           //store the number of triangles in each cell cycle

extern void AllocVarforLimitCycleDetect();
extern void finalizeLimitCycleDetect();           ////release memory allocate to variables of limit cycle detection

extern void AllocBoundaryList();

extern int Num_SmoothRegionpoints;
extern Point *point;
extern void AddToPointList(double, double, int);      ////add to region points list for smoothing
//extern void AllocateVarforSmoothing();


/////variables for shape design
//extern ctr_point *pts;          // allocate our control points array
extern int resolution;    // how many points in our output array
extern ctr_point *out_pts;

extern ctr_point *control_pts;        // allocate our control point array
extern int num_shapecontrol_pts;
extern int num_curvepts_output;
extern int MaxNumShapeControlPts;
extern int HermiteStep;

extern int *InnerTriangles ;
extern int num_innertriangles;


/////Testing variables 07/25/05
extern double MarkNextBx[2], MarkNextBy[2];
extern int *TriangleList;
extern int *cellcycle;
extern int num_celltriangle;
extern int num_trianglelist;                ////number of triangles in the triangle list
extern Edge *OneSharedEdge;
extern double problemx, problemy;
extern int test_numcelledges;




extern int *DesignCurveCellCycle;
extern int num_triangles_designcurve;


extern Edge **Cycle_edge;
extern int num_cycleedges;

extern int test_num_cycleedges;


extern RegionBoundary InnerBoundary;
extern RegionBoundary OuterBoundary;
////////////////////////////////////////////

extern DesignTriangleCycle myCycle;


/////////////////////////////////////////////
////variables for the mouse pick up singularity or limit cycle for Conley relation graph displaying
////10/16/05
int picked_sing, picked_limitcycle;
extern int picked_node;
extern int num_related_edges;
extern int *related_edges;        ////The corresponding edges incident to the node above
int *related_traj;                ////related trajectories for highlighting

extern double max_mag, min_mag;

extern int *repeller_nodes;     //the indices of the selected repellers in the conley graph
extern int NumRepellers;        //the number of the being selected repellers
extern int *attractor_nodes;    //the indices of the selected attractors in the conley graph
extern int NumAttractors;       //the number of the being selected attractors
extern int 	Num_MediaNodes ;

///////
extern GraphNode *graphnodes;
extern GraphEdge *graphedges;
extern int cur_node_index;
extern int cur_graphedge_index;

extern GraphNode2 *sccnodes;
extern GraphEdge  *sccedges;

extern MCGNode *mcgnodes;
extern MCGEdge *mcgedges;

////variables for highlighting the corresponding separatrices
extern int *related_trajs;
extern int num_related_trajs;

/* --- For strongly connected component --- */
extern SCCList scclist;

extern SamplePtsList *sampleptslist;       //each streamline will have a sample point list
extern int num_samplelists;





/***** 01/15/07 ******/
//Test the finding of being covered triangles
extern int *tri_strip;
extern int ntris_in_strip;
extern int cur_selectTris;
extern int en_tris[3];
extern icVector2 newP[3];
extern double newPx[7], newPy[7]; /*save the global coords only here! 02/26/07*/


/*Test the path between two triangles 03/13/07*/
extern int glob_t1, glob_t2;
extern int *path;
extern int num_tris_inpath;


/*  the major road network  */
extern StreetNet *majRoadnet;



extern unsigned char *popdensitymap_disp;
extern unsigned char *vegmap_disp;

struct jitter_struct{
         double x;
         double y;
} jitter_para;

jitter_struct ji1[1] = {{0.0, 0.0}};
jitter_struct ji16[16] = {{0.125, 0.125}, {0.375, 0.125}, {0.625, 0.125}, 
{0.875, 0.125},{0.125, 0.375}, {0.375, 0.375}, {0.625, 0.375}, {0.875, 0.375},
 {0.125, 0.625}, {0.375, 0.625}, {0.625, 0.625}, {0.875, 0.625},{0.125, 0.875}, 
{0.375, 0.875}, {0.625, 0.875}, {0.875, 0.875}};





/* For visualizing the curl 02/14/07 */
double curl_percent;
double curl_angle;



/*08/01/07*/
extern TraceSamplePtList *_forward, *backward_spts;
extern int back__forward;


/*tensor field visualization 09/20/2007*/
extern int ndegenerate_tris ;
extern int *degenerate_tris ;
extern Degenerate_Design *ten_designelems ;
extern int ntenelems ;
extern TenRegularElem *ten_regularelems;
extern int nten_regelems;
extern int curMaxNumTenRegElems;

/*for graph level editing of the street net 11/06/2007*/
int selectedIntersect=-1;
double zoom_factor = 1.;
double trans_x=0;
double trans_y=0;

double VisMinRoadWidth=4.;  //using opengl line width
double VisMajRoadWidth=4.;  //using opengl line width
double VisHighWayWidth=4.;  //using opengl line width


/*The quad mesh for tensor field design 09/25/2007*/
extern QuadMesh *quadmesh;
extern void display_quadmesh(GLenum mode);
extern void display_quadmesh_select(GLenum mode);
extern void display_innerverts();
extern void alloc_quad_regionsmooth();

#include "EvenlyStreamlines.h"
extern EvenStreamlinePlace *major, *minor;

extern EvenStreamlinePlace *major_level1, *minor_level1 ;
extern SeedList *seedsalongbounds;

extern RegionBlockList *sketchblocklist;
extern StreetNet *sketchnet;
extern void vis_regionblocks(RegionBlockList *, StreetNet *);


extern StreetNet *streetnet;

extern bool is_on_local_editing;


///***************************************************/
///*                                                 */
///*        The variables for the loaded maps       */
//
///*   Geographical map:  we try to load a geograph map from file*/
//unsigned char *map1=NULL;
//unsigned char *fittedmap1=NULL;
//unsigned char *displaymap=NULL;
//unsigned char *streetmapbackground=NULL;
//int boundindex1, boundindex2;
//bool ylargerthanx=false; /*false: width>=height;  true: width<height*/
//int map1_w, map1_h;
//
///*       Population density map      */
//unsigned char *popdensitymap=NULL;
//int popdensitymap_w, popdensitymap_h;
//unsigned char *popdensitymap_fit=NULL;
//unsigned char *popdensitymap_dis=NULL;
//
///*       Hieght field       */
//unsigned char *heightfield=NULL;
//int heightfield_w, heightfield_h;
//unsigned char *heightfield_fit=NULL;
//unsigned char *heightfield_dis=NULL;

#include "loadmaps.h"

/***************************************************/
/*                                                 */
/*        The variables for the loaded maps       */

/*   Geographical map:  we try to load a geograph map from file*/
extern unsigned char *map1;
extern unsigned char *fittedmap1;
extern unsigned char *displaymap;
extern unsigned char *streetmapbackground;
extern int boundindex1, boundindex2;
extern bool ylargerthanx; /*false: width>=height;  true: width<height*/
extern int map1_w, map1_h;

/*       Population density map      */
extern unsigned char *popdensitymap;
extern int popdensitymap_w, popdensitymap_h;
extern unsigned char *popdensitymap_fit;
extern unsigned char *popdensitymap_dis;

/*       Hieght field       */
extern unsigned char *heightfield;
extern int heightfield_w, heightfield_h;
extern unsigned char *heightfield_fit;
extern unsigned char *heightfield_dis;


extern bool flag_loadmap;

extern double majorDensity ;
extern double minorDensity ;

extern int cur_chosen_region;

//double designgrid_xinterval, designgrid_yinterval;

#include "shareinterfacevars.h"
extern SharedInterfaceVars sharedvars;


float inverse_tran[16];


/*    it is a tranformation function here 12/06/2007    */
void transform_fun()
{
	glTranslatef(trans_x, trans_y, 0);
	glTranslatef(0.5, 0.5, 0);
	glScalef(zoom_factor, zoom_factor, zoom_factor);
	glTranslatef(-.5,-.5, 0);
}


void cal_inverse_transform()
{
	//wglMakeCurrent(m_hDC, m_hglRC);
	glMatrixMode(GL_MODELVIEW_MATRIX);
	glPushMatrix();
	glLoadIdentity();
	glTranslatef(-trans_x, -trans_y, 0);
	glTranslatef(0.5, 0.5, 0);
	glScalef(1./zoom_factor, 1./zoom_factor, 1./zoom_factor);
	glTranslatef(-.5,-.5, 0);

	glGetFloatv(GL_MODELVIEW_MATRIX, inverse_tran);

	glPopMatrix();
}

void transform_point3D(float p[3], float rot_mat[16])
{
	double tmp[3] = {0.};

	tmp[0] = rot_mat[0] * p[0] + rot_mat[4] * p[1] + rot_mat[8]  * p[2] + rot_mat[12];
	tmp[1] = rot_mat[1] * p[0] + rot_mat[5] * p[1] + rot_mat[9]  * p[2] + rot_mat[13];
	tmp[2] = rot_mat[2] * p[0] + rot_mat[6] * p[1] + rot_mat[10] * p[2] + rot_mat[14];

	p[0] = tmp[0];
	p[1] = tmp[1];
	p[2] = tmp[2];
}


// CGlView

////Not a good way !!!! But no choice at present
int CGlView::MoveOrStop = 1;

IMPLEMENT_DYNAMIC(CGlView, CWnd)

CGlView::CGlView(CWnd *pclWnd)
{

    m_pclWnd = pclWnd;
    m_hWnd   = pclWnd->m_hWnd;
    m_hDC    = ::GetDC(m_hWnd);

	////Initialize possible used vector field variables here
	InitVFVariables();
	
	////Initialize flags
    InitFlag();

	MoveOrStop = 1;
	EditModeOn = 0;
	SingularitiesOn = 1;
	RegularElemOn = 1;
	
	EditModeOn = 0;
	TrajectoryOn = 0;
	SeparatricesOn = 0;
	LimitCycleOn = 0;
	SmoothOn = 0;
	DisplaySmoothRegionOn = 0;
	PickPointOn = 0;

	choose_ID = -1;
	SelectTriangleID = -1;
	TraceBeginTriangleID = -1;

	zoom_factor = 1.;
	trans_x=0;
	trans_y=0;
}

CGlView::CGlView()
{
	MoveOrStop = 1;
	EditModeOn = 0;
	SingularitiesOn = 1;
	RegularElemOn = 1;
	
	EditModeOn = 0;
	TrajectoryOn = 0;
	SeparatricesOn = 0;
	LimitCycleOn = 0;
	SmoothOn = 0;
	DisplaySmoothRegionOn = 0;
	PickPointOn = 0;

	choose_ID = -1;
	SelectTriangleID = -1;
	TraceBeginTriangleID = -1;
}

CGlView::~CGlView()
{
	finalize();
}


BEGIN_MESSAGE_MAP(CGlView, CWnd)
	ON_WM_DESTROY()
	ON_WM_ERASEBKGND()
END_MESSAGE_MAP()



// CGlView message handlers
int CGlView::OnCreate() 
{

    m_hDC = ::GetDC(this->m_hWnd);

    if(!SetPixelformat(m_hDC))
    {
	::MessageBox(::GetFocus(),"SetPixelformat Failed!","Error",MB_OK);
	return -1;
    }

    m_hglRC = wglCreateContext(m_hDC);
    int i= wglMakeCurrent(m_hDC,m_hglRC);

	InitGL();	

	return 0;
}

void CGlView::OnDestroy()
{
	CWnd::OnDestroy();

	// TODO: Add your message handler code here
    wglMakeCurrent(NULL,NULL);
    wglDeleteContext(m_hglRC);	
}

BOOL CGlView::OnEraseBkgnd(CDC* pDC)
{
	// TODO: Add your message handler code here and/or call default

	return CWnd::OnEraseBkgnd(pDC);
}

/*--------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------*/
// Other extended operations for OpenGl window setting

BOOL CGlView::SetPixelformat(HDC hdc)
{

    PIXELFORMATDESCRIPTOR *ppfd; 
    int pixelformat; 
 
    PIXELFORMATDESCRIPTOR pfd = { 
    sizeof(PIXELFORMATDESCRIPTOR),  //  size of this pfd 
    1,                     // version number 
    PFD_DRAW_TO_WINDOW |   // support window 
    PFD_SUPPORT_OPENGL |   // support OpenGL 
    PFD_GENERIC_FORMAT |
    PFD_DOUBLEBUFFER,      // double buffered 
    PFD_TYPE_RGBA,         // RGBA type 
    32,                    // 24-bit color depth 
    0, 0, 0, 0, 0, 0,      // color bits ignored 
    8,                     // no alpha buffer 
    0,                     // shift bit ignored 
    8,                     // no accumulation buffer 
    0, 0, 0, 0,            // accum bits ignored 
    32,                    // 32-bit z-buffer	 
    1,                     // no stencil buffer 
    8,                     // no auxiliary buffer 
    PFD_MAIN_PLANE,        // main layer 
    0,                     // reserved 
    0, 0, 0                // layer masks ignored 
    }; 

   
    ppfd = &pfd;

 
    if ( (pixelformat = ChoosePixelFormat(hdc, ppfd)) == 0 ) 
    { 
        ::MessageBox(NULL, "ChoosePixelFormat failed", "Error", MB_OK); 
        return FALSE; 
    } 
 
    if (SetPixelFormat(hdc, pixelformat, ppfd) == FALSE) 
    { 
        ::MessageBox(NULL, "SetPixelFormat failed", "Error", MB_OK); 
        return FALSE; 
    } 
 
    return TRUE; 

}


int CGlView::InitGL(GLvoid)								// All Setup For OpenGL Goes Here
{
	//glViewport(0, 0, (GLsizei) NPIX, (GLsizei) NPIX);
	glViewport(0, 0, (GLsizei)REALWINSIZE, (GLsizei)REALWINSIZE);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	gluOrtho2D(0, 1, 0, 1);
	//glTranslatef(-1.0, -1.0, -1.0); 
	//glScalef(2.0, 2.0, 1.0);
	glTexParameteri(GL_TEXTURE_2D, 
					GL_TEXTURE_WRAP_S, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, 
					GL_TEXTURE_WRAP_T, GL_REPEAT); 
	glTexParameteri(GL_TEXTURE_2D, 
					GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, 
					GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, 
					GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClear(GL_COLOR_BUFFER_BIT);

	glDisable(GL_STENCIL_TEST);

	return TRUE;										// Initialization Went OK
}
/*--------------------------------------------------------------------------------*/

void  CGlView::HsvRgb( float hsv[3], float rgb[3] )
{
	float h, s, v;			// hue, sat, value
	float r, g, b;			// red, green, blue
	float i, f, p, q, t;		// interim values

	// guarantee valid input:
	h = hsv[0] / 60.;
	while( h >= 6. )	h -= 6.;
	while( h <  0. ) 	h += 6.;

	s = hsv[1];
	if( s < 0. )
		s = 0.;
	if( s > 1. )
		s = 1.;

	v = hsv[2];
	if( v < 0. )
		v = 0.;
	if( v > 1. )
		v = 1.;

	// if sat==0, then is a gray:
	if( s == 0.0 )
	{
		rgb[0] = rgb[1] = rgb[2] = v;
		return;
	}

	// get an rgb from the hue itself:
	i = floor( h );
	f = h - i;
	p = v * ( 1. - s );
	q = v * ( 1. - s*f );
	t = v * ( 1. - ( s * (1.-f) ) );

	switch( (int) i )
	{
		case 0:
			r = v;	g = t;	b = p;
			break;
	
		case 1:
			r = q;	g = v;	b = p;
			break;
	
		case 2:
			r = p;	g = v;	b = t;
			break;
	
		case 3:
			r = p;	g = q;	b = v;
			break;
	
		case 4:
			r = t;	g = p;	b = v;
			break;
	
		case 5:
			r = v;	g = p;	b = q;
			break;
	}

	rgb[0] = r;
	rgb[1] = g;
	rgb[2] = b;
}

/*--------------------------------------------------------------------------------*/
//Display routine
/*--------------------------------------------------------------------------------*/
#include "tensorvis.h"
//int CGlView::DrawGLScene(GLenum mode)					// Here's Where We Do All The Drawing
//{
//
//	//glClearColor(1, 1, 1, 1);
//	//glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
//	//glEnable(GL_TEXTURE_2D);
//
//	if(displayRoadNetOn)
//	{
//		glLineWidth(1.0);
//		glClear(GL_ACCUM_BUFFER_BIT);
//		for(int i = 0; i < 16; i++)
//		{
//			glPushMatrix ();
//			glTranslatef (ji16[i].x*1.0/**quadmesh->radius*//512, ji16[i].y*1.0/**quadmesh->radius*//512, 0.0);
//
//			/*we should put the transformation here 11/06/2007*/
//			glTranslatef(trans_x, trans_y, 0);
//			glTranslatef(0.5, 0.5, 0);
//			glScalef(zoom_factor, zoom_factor, zoom_factor);
//			glTranslatef(-.5, -.5, 0);
//
//
//			//display_roads(mode);
//			/*directly use the obtained tensor lines to visualize the street network*/
//			display_roads_width(mode);
//
//			if(sharedvars.ShowRegionBlocksOn)
//				vis_regionblocks();
//			if(showStreetGraphOn)
//				display_streetnet(mode);
//
//			if(showTensorLineOn)
//			{
//				/*09/30/2007 draw evenly placed tensor lines here*/
//				display_major_tenlines(mode);
//				display_minor_tenlines(mode);
//			}
//			
//			if(sharedvars.DesignGridOn)
//				display_design_grid();
//
//			if(sharedvars.EnableSketchBasedDesign||	sharedvars.ShowSketchesOn)
//			{
//				display_brush_sketch(GL_RENDER);
//				if(sharedvars.ShowRegionBlocksOn)
//					vis_regionblocks(sketchblocklist, sketchnet);
//			}
//
//			glPopMatrix ();
//			glAccum(GL_ACCUM, 1.0/16);
//		}
//		glAccum (GL_RETURN, 1.0);
//		glReadBuffer(GL_BACK);
//		glLineWidth(1.);
//		SwapBuffers(m_hDC);
//		return TRUE;
//	}
//
////	if(displayDisMapOn == 1)
////	{
////	glClearColor(1, 1, 1, 1);
////	glClear(GL_COLOR_BUFFER_BIT);
////	glDisable(GL_TEXTURE_2D);
////	glEnable(GL_COLOR_MATERIAL);
////extern void vis_distance();
////		vis_distance();
////extern void display_zero_level();
////		display_zero_level();
////
////			if(sharedvars.DesignGridOn)
////				display_design_grid();
////
////		SwapBuffers(m_hDC);
////		return TRUE;
////	}
//
//	else if(sharedvars.ShowTheMapOn)
//	{
//			/*we should put the transformation here 11/06/2007*/
//		glLineWidth(1.0);
//		glClear(GL_ACCUM_BUFFER_BIT);
//		for(int i = 0; i < 16; i++)
//		{
//			glPushMatrix ();
//			glTranslatef (ji16[i].x*1.0/**quadmesh->radius*//512, ji16[i].y*1.0/**quadmesh->radius*//512, 0.0);
//			glTranslatef(trans_x, trans_y, 0);
//			glTranslatef(0.5, 0.5, 0);
//			glScalef(zoom_factor, zoom_factor, zoom_factor);
//			glTranslatef(-.5, -.5, 0);
//
//		render_a_map(streetmapbackground);
//			glPopMatrix ();
//			glAccum(GL_ACCUM, 1.0/16);
//		}
//		glAccum (GL_RETURN, 1.0);
//		glReadBuffer(GL_BACK);
//		glLineWidth(1.);
//	}
//
//	else
//	{
//		glDisable(GL_COLOR_MATERIAL);
//
//
//		if(StreamlineBasedOn == 0 && showTensorOn == 0)
//		{
//			IBFVEffect(mode);
//		   
//			glDisable(GL_TEXTURE_2D);
//
//			if(ColorPlotOn == 1)
//			{
//				////Draw some color here
//				glEnable(GL_BLEND); 
//				DisplayColorPlots();
//				glDisable(GL_BLEND); 
//			}
//		}
//
//		else if(sharedvars.ShowIBFVOn && showTensorOn == 1) /*visualize tensor field*/
//		{
//			glEnable(GL_TEXTURE_2D);
//			glShadeModel(GL_SMOOTH);
//
//			/*using triangle mesh*/
//			////render_majorfield();
//			//////render_minorfield();
//			//////mix_vis();
//			//////vis_alpha_map();
//
//			/*using quad mesh 09/25/2007*/
//			render_majorfield_quad();
//
//			glDisable(GL_TEXTURE_2D);
//
//			glEnable(GL_COLOR_MATERIAL);
//			
//		}
//
//		else
//		{
//			glClearColor(1, 1, 1, 1);
//			glClear(GL_COLOR_BUFFER_BIT);
//			glDisable(GL_TEXTURE_2D);
//		}
//	}
//			
//	glDisable(GL_TEXTURE_2D);
//
//    glEnable(GL_COLOR_MATERIAL);
//
//	if(EditModeOn == 1)
//		DisplayEditBox(mode);
//    
//	if((EditModeOn == 1 || sharedvars.MoveElemOn ||	sharedvars.RemoveElemOn) && showTensorOn == 1)
//		display_tenElem_EditBox(mode);
//
//
//	/*
//	Tensor element editing
//	*/
//	if(chosen_tenelem_ID >= 0 && EditModeOn == 1 && showTensorOn == 1)
//	{
//		////if it is singular element being selected
//		if(chosen_tenelem_ID < NAMEOFREGELEM)
//		    display_singularElem_controlPts_tensor(mode, chosen_tenelem_ID);
//	}
//
//	if(RegularElemOn == 1 && showTensorOn == 1)
//	    display_tenRegElem(mode);
//
//
//	if(DisplaySmoothRegionOn == 1)
//		DisplaySmoothRegion();
//	
//	////Using antialiasing 07/27/06
//	glLineWidth(2.0);
//	glClear(GL_ACCUM_BUFFER_BIT);
//	for(int i = 0; i < 16; i++)
//	{
//		glPushMatrix ();
//		glTranslatef (ji16[i].x*1.0/512, ji16[i].y*1.0/512, 0.0);
//
//		if(showTensorLineOn)
//		{
//			/*09/30/2007 draw evenly placed tensor lines here*/
//			display_major_tenlines(mode);
//			display_minor_tenlines(mode);
//
//			//display_major_samps();
//			//display_minor_samps();
//		}
//
//		if(	showStreetGraphOn)
//		{
//			display_streetnet(mode);
//		}
//
//		if(sharedvars.DesignGridOn)
//		    display_design_grid();
//
//
//
//	if(SingularitiesOn == 1)
//	{
//		display_degenerate_pts(mode);
//	}
//		
//
//		glPopMatrix ();
//		glAccum(GL_ACCUM, 1.0/16);
//	}
//	glAccum (GL_RETURN, 1.0);
//	glReadBuffer(GL_BACK);
//	glLineWidth(1.);
//	
//
//	if(ShowSamplePtsOn == 1)
//		DisplaySamplingPts();
//
//
//	if(RegularElemOn == 1)
//	    DisplayRegularElemt(mode);
//
//
//	////Display the shape control
//	if(ShapeControlPtsOn == 1 )
//	{
//		DisplayDesignCurve();
//		DisplayShapeControlPts(GL_RENDER);
//
//		//if(showTensorOn == 1)
//		display_brush();
//
//	}
//
//	if(sharedvars.EnableSketchBasedDesign||	sharedvars.ShowSketchesOn)
//	{
//		display_brush_sketch(GL_RENDER);
//		if(sharedvars.ShowRegionBlocksOn)
//			vis_regionblocks(sketchblocklist, sketchnet);
//	}
//
//	//////////////////////////////////////////////////////////
//	/////Testing displaying codes
//	if(ShowCancelSmoothRegion == 1)
//	{
//		if(!showTensorOn)
//			TestDisplayRegion(mode);
//		else
//		{
//			display_quadmesh(mode);
//			display_innerverts();
//		}
//		ShowAllBoundaries();
//	}
//
//
//	/*display the tracing samples with colors*/
//	//if(showTraceSamplingPts == 1)
//	//{
//	//	display_trace_SampltPts(back__forward);
//	//}
//
//
//	//if(showEdgeImageOn == 1)
//	//{
//	//	display_image_edge();
//	//}
//
//
//
//	glDisable(GL_COLOR_MATERIAL);
//  	glEnable(GL_TEXTURE_2D);
//
//	SwapBuffers(m_hDC);
//	return TRUE;										// Keep Going
//}


extern int *innerintersections ;
				extern int ninnerintersections ;
				extern int *region_quadverts;
				extern int nregion_quadverts;
//extern int *innercells;
//extern int ninnercells;
extern Trajectory *regionboundary;
extern int *contour_cells;
extern int ncontour_cells;


int CGlView::with_anti_aliasing(GLenum mode)
{
	if(displayRoadNetOn)
	{
		glLineWidth(1.0);
		glClear(GL_ACCUM_BUFFER_BIT);
		//glEnable(GL_COLOR_MATERIAL);
		for(int i = 0; i < 16; i++)
		{
			glPushMatrix ();
			glTranslatef (ji16[i].x*1./**quadmesh->radius*//512, ji16[i].y*1./**quadmesh->radius*//512, 0.0);

			/*we should put the transformation here 11/06/2007*/
			transform_fun();


			//////display_roads(mode);

			/*directly use the obtained tensor lines to visualize the street network*/
			if(sharedvars.ShowRoadMapOn)
			{
				if(!sharedvars.ShowLineStyleStreetsOn)
					display_roads_width(mode);
				else
				{
					display_road_graph_linestyle();
				}
			}

			else if(sharedvars.ShowStreetUseNetworkOn)
				display_road_use_network();

			if(sharedvars.ShowRegionBlocksOn)
				vis_regionblocks();
			if(showStreetGraphOn)
				display_streetnet(mode);

			if(sharedvars.ShowMajRoadNetworkOn)
				display_majRoadnet(mode);

			if(showTensorLineOn)
			{
				/*09/30/2007 draw evenly placed tensor lines here*/
				display_major_tenlines(mode);
				display_minor_tenlines(mode);

				if(sharedvars.ShowMajRoadsOn)
				{
					if(sharedvars.ShowMajRoadGoogleStyleOn)
						display_majRoad_googlestyle(mode);
					else
						display_level1_tenlines(mode);
				}
			}

			if(sharedvars.ShowInitSeedsOn)
				display_init_seeds(mode);
			
			if(sharedvars.DesignGridOn)
				display_design_grid();

			if(sharedvars.EnableSketchBasedDesign||	sharedvars.ShowSketchesOn)
			{
				display_brush_sketch_thin(GL_RENDER);
				display_brush_sketch_google(GL_RENDER);
				if(sharedvars.ShowRegionBlocksOn)
					vis_regionblocks(sketchblocklist, sketchnet);
			}
			
			//   Display the shape control
			if(ShapeControlPtsOn == 1 )
			{
				DisplayDesignCurve();
				DisplayShapeControlPts(GL_RENDER);

				//if(showTensorOn == 1)
				display_brush();
			}

			if(DisplaySmoothRegionOn == 1 )
			{
				DisplaySmoothRegion();
				if(regionboundary != NULL)
					{

						/*show all the marked cells */
						//for(int j=0;j<quadmesh->nfaces;j++)
						//{
						//	QuadCell *face=quadmesh->quadcells[j];
						//	if(!face->is_contour_cell) continue;
						//		glColor3f(0, 1, 0);
						//	glBegin(GL_LINE_LOOP);
						//	for(int k=0;k<face->nverts;k++)
						//	{
						//		QuadVertex *v=quadmesh->quad_verts[face->verts[k]];
						//		glVertex2f(v->x,v->y);
						//	}
						//	glEnd();

						//}

						/*  show all the contour cells  */
						//for(int j=0;j<ncontour_cells;j++)
						//{
						//	QuadCell *face=quadmesh->quadcells[contour_cells[j]];
						//	if(!face->is_contour_cell) continue;
						//		glColor3f(1, 1, 0);
						//	glBegin(GL_LINE_LOOP);
						//	for(int k=0;k<face->nverts;k++)
						//	{
						//		QuadVertex *v=quadmesh->quad_verts[face->verts[k]];
						//		glVertex2f(v->x,v->y);
						//	}
						//	glEnd();
						//}

						for(int j=0;j<regionboundary->nlinesegs;j++)
						{
							if(regionboundary->linesegs[j].Triangle_ID<0) continue;

							if(i%2==0)
								glColor3f(0., 1., 1.);
							else
								glColor3f(1, 0, 1);

							glBegin(GL_LINES);
								glVertex2f(regionboundary->linesegs[j].gstart[0],
									regionboundary->linesegs[j].gstart[1]);
								glVertex2f(regionboundary->linesegs[j].gend[0],
									regionboundary->linesegs[j].gend[1]);
							glEnd();

							//QuadCell *face=quadmesh->quadcells[regionboundary->linesegs[i].Triangle_ID];
							//	glColor3f(0, 1, 0);
							//glBegin(GL_LINE_LOOP);
							//for(int k=0;k<face->nverts;k++)
							//{
							//	QuadVertex *v=quadmesh->quad_verts[face->verts[k]];
							//	glVertex2f(v->x,v->y);
							//}
							//glEnd();
						}

					}

				if(innerintersections!=NULL && region_quadverts != NULL
					&& sharedvars.SelStreetRegToEditOn)
				{
				//
				///*  display inner vertices and inner intersections */
					int i;
				//	glColor3f(1,0,0);
				//	glPointSize(3.);
				//	glBegin(GL_POINTS);
				//	for(i=0;i<ninnerintersections;i++)
				//	{
				//		glVertex2f(streetnet->nodelist->intersects[innerintersections[i]]->gpos[0],
				//			streetnet->nodelist->intersects[innerintersections[i]]->gpos[1]);
				//	}
				//	glEnd();
				//	
				//	glColor3f(0,1,0);
				//	glBegin(GL_POINTS);
				//	for(i=0;i<nregion_quadverts;i++)
				//	{
				//		glVertex2f(quadmesh->quad_verts[region_quadverts[i]]->x,
				//			quadmesh->quad_verts[region_quadverts[i]]->y);
				//	}
				//	glEnd();
				//	glPointSize(1.);

				//	glColor3f(0,0,1);
				//	for(i=0;i<ninnercells;i++)
				//	{
				//		QuadCell *face=quadmesh->quadcells[innercells[i]];
				//		glBegin(GL_LINE_LOOP);
				//		for(int j=0;j<face->nverts;j++)
				//		{
				//			QuadVertex *v=quadmesh->quad_verts[face->verts[j]];
				//			glVertex2f(v->x,
				//				v->y);
				//		}
				//		glEnd();
				//	}
				//	glPointSize(1.);

					//glColor3f(0,1,0);
					//for(i=0;i<streetnet->edgelist->nedges;i++)
					//{
					//	StreetGraphEdge *edge=streetnet->edgelist->edges[i];
					//	if(edge->cancel || !edge->visited)
					//		continue;

					//	Intersection *intersect1=streetnet->nodelist->intersects[edge->node_index1];
					//	Intersection *intersect2=streetnet->nodelist->intersects[edge->node_index2];
					//	glBegin(GL_LINES);
					//	glVertex2f(intersect1->gpos[0], intersect1->gpos[1]);
					//	glVertex2f(intersect2->gpos[0], intersect2->gpos[1]);
					//	glEnd();
					//}

					}
			}



			glPopMatrix ();
			glAccum(GL_ACCUM, 1.0/16);
		}
		glAccum (GL_RETURN, 1.0);
		glReadBuffer(GL_BACK);
		glLineWidth(1.);
		SwapBuffers(m_hDC);
		return TRUE;
	}

	else if(sharedvars.ShowPopDensityMapOn)
	{
		glPushMatrix ();
		transform_fun();
		if(popdensitymap_disp!=NULL)
		render_a_map(popdensitymap_disp);
		glPopMatrix();
		SwapBuffers(m_hDC);
		return TRUE;
	}

	else if(sharedvars.ShowVegMapOn)
	{
		glPushMatrix ();
		transform_fun();
		if(vegmap_disp!=NULL)
		render_a_map(vegmap_disp);
		glPopMatrix();
		SwapBuffers(m_hDC);
		return TRUE;
	}

	else if(sharedvars.ShowScalarFieldOn)
	{
		glClearColor(1, 1, 1, 1);
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	
		glDisable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);
		/*  remember to switch to the smooth mode  */
		glEnable(GL_COLOR_MATERIAL);
		glShadeModel(GL_SMOOTH);

		glPushMatrix ();

		transform_fun();
		vis_scalarfield();

		glPopMatrix();

		SwapBuffers(m_hDC);
		return TRUE;
	}
	else
	{
		////Using antialiasing 07/27/06

		//glDrawBuffer(GL_FRONT_AND_BACK);
		glClearColor(0.8, 0.8, 0.9, 1);
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
			
		glPushMatrix ();

		transform_fun();


			if(sharedvars.ShowTheMapOn)
			{
				render_a_map(displaymap);
			}

			else
			{
				//glDisable(GL_COLOR_MATERIAL);

			    if(sharedvars.ShowIBFVOn && showTensorOn == 1) /*visualize tensor field*/
				{
					glEnable(GL_TEXTURE_2D);
					glShadeModel(GL_FLAT);

					/*using quad mesh 09/25/2007*/
					render_majorfield_quad();  /*  showing major field only  */

					if(is_on_local_editing)
					{
						if(streetnet!=NULL)
						{
							/* rend the back of the roads */
							display_road_graph_back(0, MINOR);
							display_road_graph_back(0, MAJOR);
							display_road_graph_back(0, HIGHWAY);

							/* rend the front of the roads */
							display_road_graph_front(0, MINOR);
							display_road_graph_front(0, MAJOR);
							display_road_graph_front(0, HIGHWAY);
						}
					}

				}

				else
				{
					glClearColor(1, 1, 1, 1);
					glClear(GL_COLOR_BUFFER_BIT);
					glDisable(GL_TEXTURE_2D);
				}
			}
					
			glDisable(GL_TEXTURE_2D);

			glEnable(GL_COLOR_MATERIAL);

		glLineWidth(2.0);
		glClear(GL_ACCUM_BUFFER_BIT);
		for(int i = 0; i < 16; i++)
		{
			glPushMatrix ();
			glTranslatef (ji16[i].x*1.0/512, ji16[i].y*1.0/512, 0.0);

			if(EditModeOn == 1)
				DisplayEditBox(mode);
		    
			if((EditModeOn == 1 || sharedvars.MoveElemOn ||	sharedvars.RemoveElemOn) 
				&& showTensorOn == 1)
				display_tenElem_EditBox(mode);


			/*
			Tensor element editing
			*/
			if(chosen_tenelem_ID >= 0 && EditModeOn == 1 && showTensorOn == 1)
			{
				////if it is singular element being selected
				if(chosen_tenelem_ID < NAMEOFREGELEM)
					display_singularElem_controlPts_tensor(mode, chosen_tenelem_ID);
			}

			if(RegularElemOn == 1 && showTensorOn == 1)
				display_tenRegElem(mode);


			if(DisplaySmoothRegionOn == 1)
				DisplaySmoothRegion();
			

			if(showTensorLineOn)
			{
				/*09/30/2007 draw evenly placed tensor lines here*/
				display_major_tenlines(mode);
				display_minor_tenlines(mode);

				if(sharedvars.ShowMajRoadsOn)
				{
					if(sharedvars.ShowMajRoadGoogleStyleOn)
						display_majRoad_googlestyle(mode);
					else
						display_level1_tenlines(mode);
				}
			}

			if(sharedvars.ShowInitSeedsOn)
				display_init_seeds(mode);

			if(	showStreetGraphOn)
			{
				display_streetnet(mode);
			}

			if(sharedvars.ShowMajRoadNetworkOn)
				display_majRoadnet(mode);

			if(sharedvars.DesignGridOn)
				display_design_grid();

			if(SingularitiesOn == 1)
			{
				display_degenerate_pts(mode);
			}
				
			if(sharedvars.EnableSketchBasedDesign||	sharedvars.ShowSketchesOn)
			{
				display_brush_sketch_thin(GL_RENDER);
				if(sharedvars.ShowRegionBlocksOn)
					vis_regionblocks(sketchblocklist, sketchnet);
			}

				glPopMatrix ();
				glAccum(GL_ACCUM, 1.0/16);
			}
			glAccum (GL_RETURN, 1.0);
			glReadBuffer(GL_BACK);
			glLineWidth(1.);

			//if(ShowSamplePtsOn == 1)
			//	DisplaySamplingPts();


			//if(RegularElemOn == 1)
			//	DisplayRegularElemt(mode);


			////Display the shape control
			if(ShapeControlPtsOn == 1 )
			{
				DisplayDesignCurve();
				DisplayShapeControlPts(GL_RENDER);

				//if(showTensorOn == 1)
				display_brush();
			}


			//////////////////////////////////////////////////////////
			/////Testing displaying codes
			if(ShowCancelSmoothRegion == 1)
			{
				if(!showTensorOn)
					TestDisplayRegion(mode);
				else
				{
					display_quadmesh(mode);
					display_innerverts();
				}
				ShowAllBoundaries();
			}

		//	glPopMatrix ();
		//	glAccum(GL_ACCUM, 1.0/16);
		//}
		//glAccum (GL_RETURN, 1.0);
		//glReadBuffer(GL_BACK);
		//glLineWidth(1.);

		glPopMatrix ();

		glDisable(GL_COLOR_MATERIAL);
  		glEnable(GL_TEXTURE_2D);

		SwapBuffers(m_hDC);
		return TRUE;										// Keep Going
	}
}
int CGlView::without_anti_aliasing(GLenum mode)
{
	if(displayRoadNetOn)
	{

			/*we should put the transformation here 11/06/2007*/
		glPushMatrix();

			transform_fun();
			//////display_roads(mode);

			/*directly use the obtained tensor lines to visualize the street network*/
			if(sharedvars.ShowRoadMapOn)
			{
				if(!sharedvars.ShowLineStyleStreetsOn)
					display_roads_width(mode);
				else
				{
					display_road_graph_linestyle();
				}
			}

			else if(sharedvars.ShowStreetUseNetworkOn)
				display_road_use_network();

			if(sharedvars.ShowRegionBlocksOn)
				vis_regionblocks();
			if(showStreetGraphOn)
				display_streetnet(mode);

			if(sharedvars.ShowMajRoadNetworkOn)
				display_majRoadnet(mode);

			if(showTensorLineOn)
			{
				/*09/30/2007 draw evenly placed tensor lines here*/
				display_major_tenlines(mode);
				display_minor_tenlines(mode);

				if(sharedvars.ShowMajRoadsOn)
				{
					if(sharedvars.ShowMajRoadGoogleStyleOn)
						display_majRoad_googlestyle(mode);
					else
						display_level1_tenlines(mode);
				}
			}

			if(sharedvars.ShowInitSeedsOn)
				display_init_seeds(mode);
			
			if(sharedvars.DesignGridOn)
				display_design_grid();

			if(sharedvars.EnableSketchBasedDesign||	sharedvars.ShowSketchesOn)
			{
				display_brush_sketch_thin(GL_RENDER);
				display_brush_sketch_google(GL_RENDER);
				if(sharedvars.ShowRegionBlocksOn)
					vis_regionblocks(sketchblocklist, sketchnet);
			}
			
			//   Display the shape control
			if(ShapeControlPtsOn == 1 )
			{
				DisplayDesignCurve();
				DisplayShapeControlPts(GL_RENDER);

				display_brush();
			}

			if(DisplaySmoothRegionOn == 1 )
			{
				DisplaySmoothRegion();
				if(regionboundary != NULL)
					{

						for(int j=0;j<regionboundary->nlinesegs;j++)
						{
							if(regionboundary->linesegs[j].Triangle_ID<0) continue;

							//if(j%2==0)
							//	glColor3f(0., 1., 1.);
							//else
							//	glColor3f(1, 0, 1);

							glBegin(GL_LINES);
								glVertex2f(regionboundary->linesegs[j].gstart[0],
									regionboundary->linesegs[j].gstart[1]);
								glVertex2f(regionboundary->linesegs[j].gend[0],
									regionboundary->linesegs[j].gend[1]);
							glEnd();

						}

					}

				if(innerintersections!=NULL && region_quadverts != NULL
					&& sharedvars.SelStreetRegToEditOn)
				{

					}
			}

		glPopMatrix();

		glLineWidth(1.);
		SwapBuffers(m_hDC);
		return TRUE;
	}

	else if(sharedvars.ShowPopDensityMapOn)
	{
		glPushMatrix ();
		transform_fun();
		if(popdensitymap_disp!=NULL)
		render_a_map(popdensitymap_disp);
		glPopMatrix();
		SwapBuffers(m_hDC);
		return TRUE;
	}

	else if(sharedvars.ShowVegMapOn)
	{
		glPushMatrix ();
		transform_fun();
		if(vegmap_disp!=NULL)
		render_a_map(vegmap_disp);
		glPopMatrix();
		SwapBuffers(m_hDC);
		return TRUE;
	}

	else if(sharedvars.ShowScalarFieldOn)
	{
		glClearColor(1, 1, 1, 1);
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	
		glDisable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);
		/*  remember to switch to the smooth mode  */
		glEnable(GL_COLOR_MATERIAL);
		glShadeModel(GL_SMOOTH);

		glPushMatrix ();

		transform_fun();
		vis_scalarfield();

		glPopMatrix();

		SwapBuffers(m_hDC);
		return TRUE;
	}
	else
	{
		////Using antialiasing 07/27/06

		//glDrawBuffer(GL_FRONT_AND_BACK);
		//237 234 226 
		glClearColor(0.93, 0.93, 0.87, 1);
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
			
		glPushMatrix ();

		transform_fun();


			if(sharedvars.ShowTheMapOn)
			{
				render_a_map(displaymap);
			}

			else
			{
				//glDisable(GL_COLOR_MATERIAL);

			    if(sharedvars.ShowIBFVOn && showTensorOn == 1) /*visualize tensor field*/
				{
					glEnable(GL_TEXTURE_2D);
					glShadeModel(GL_FLAT);

					/*using quad mesh 09/25/2007*/
					render_majorfield_quad();  /*  showing major field only  */
					if(is_on_local_editing)
					{
						if(streetnet!=NULL)
						{
							/* rend the back of the roads */
							display_road_graph_back(0, MINOR);
							display_road_graph_back(0, MAJOR);
							display_road_graph_back(0, HIGHWAY);

							/* rend the front of the roads */
							display_road_graph_front(0, MINOR);
							display_road_graph_front(0, MAJOR);
							display_road_graph_front(0, HIGHWAY);
						}
					}

				}

				else
				{
					glClearColor(1, 1, 1, 1);
					glClear(GL_COLOR_BUFFER_BIT);
					glDisable(GL_TEXTURE_2D);
				}
			}
					
			glDisable(GL_TEXTURE_2D);

			glEnable(GL_COLOR_MATERIAL);

		glLineWidth(2.0);

			if(EditModeOn == 1)
				DisplayEditBox(mode);
		    
			if((EditModeOn == 1 || sharedvars.MoveElemOn ||	sharedvars.RemoveElemOn) 
				&& showTensorOn == 1)
				display_tenElem_EditBox(mode);


			/*
			Tensor element editing
			*/
			if(chosen_tenelem_ID >= 0 && EditModeOn == 1 && showTensorOn == 1)
			{
				////if it is singular element being selected
				if(chosen_tenelem_ID < NAMEOFREGELEM)
					display_singularElem_controlPts_tensor(mode, chosen_tenelem_ID);
			}

			if(RegularElemOn == 1 && showTensorOn == 1)
				display_tenRegElem(mode);


			if(DisplaySmoothRegionOn == 1)
				DisplaySmoothRegion();
			

			if(showTensorLineOn)
			{
				/*09/30/2007 draw evenly placed tensor lines here*/
				display_major_tenlines(mode);
				display_minor_tenlines(mode);

				if(sharedvars.ShowMajRoadsOn)
				{
					if(sharedvars.ShowMajRoadGoogleStyleOn)
						display_majRoad_googlestyle(mode);
					else
						display_level1_tenlines(mode);
				}
			}

			if(sharedvars.ShowInitSeedsOn)
				display_init_seeds(mode);

			if(	showStreetGraphOn)
			{
				display_streetnet(mode);
			}

			if(sharedvars.ShowMajRoadNetworkOn)
				display_majRoadnet(mode);

			if(sharedvars.DesignGridOn)
				display_design_grid();

			if(SingularitiesOn == 1)
			{
				display_degenerate_pts(mode);
			}
				
			if(sharedvars.EnableSketchBasedDesign||	sharedvars.ShowSketchesOn)
			{
				display_brush_sketch_thin(GL_RENDER);
				if(sharedvars.ShowRegionBlocksOn)
					vis_regionblocks(sketchblocklist, sketchnet);
			}


			////Display the shape control
			if(ShapeControlPtsOn == 1 )
			{
				DisplayDesignCurve();
				DisplayShapeControlPts(GL_RENDER);

				//if(showTensorOn == 1)
				display_brush();
			}


			//////////////////////////////////////////////////////////
			/////Testing displaying codes
			if(ShowCancelSmoothRegion == 1)
			{
				if(!showTensorOn)
					TestDisplayRegion(mode);
				else
				{
					display_quadmesh(mode);
					display_innerverts();
				}
				ShowAllBoundaries();
			}


		glPopMatrix ();

		glDisable(GL_COLOR_MATERIAL);
  		glEnable(GL_TEXTURE_2D);

		SwapBuffers(m_hDC);
		return TRUE;										// Keep Going
	}
}


int CGlView::DrawGLScene(GLenum mode)					// Here's Where We Do All The Drawing
{
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	glDisable(GL_LIGHTING);



	if(sharedvars.AntiAliasingOn)
		return with_anti_aliasing(mode);
	else
		return without_anti_aliasing(mode);
}


void CGlView::DisplayCoveredTris()
{
	Face *face;
	Vertex *v;

	int i, j;
	for(i = 0; i < ntris_in_strip; i++)
	{
		face = Object.flist[tri_strip[i]];

		glEnable(GL_BLEND);
		glColor4f(1, 0, 1, 0.4);
		glBegin(GL_TRIANGLES);
		for(j = 0; j < 3; j++)
		{
			v = Object.vlist[face->verts[j]];

			glVertex2f(v->x, v->y);
		}
		glEnd();
		glDisable(GL_BLEND);
		
		glColor3f(0.6, 0.6, 0.6);
		glBegin(GL_LINE_LOOP);
		for(j = 0; j < 3; j++)
		{
			v = Object.vlist[face->verts[j]];

			glVertex2f(v->x, v->y);
		}
		glEnd();
	}

	glEnable(GL_BLEND);
	for(i = 0; i < cur_selectTris; i++)
	{
		face = Object.flist[en_tris[i]];

		glColor4f(1, 0, 1, 0.4);
		glBegin(GL_TRIANGLES);
		for(j = 0; j < 3; j++)
		{
			v = Object.vlist[face->verts[j]];

			glVertex2f(v->x, v->y);
		}
		glEnd();
	}
	glDisable(GL_BLEND);

	glBegin(GL_LINES);
	for(i = 0; i < cur_selectTris; i++)
	{
		glColor3f(1, 1, 1);
		glVertex2f(newP[i].entry[0], newP[i].entry[1]);
		glVertex2f(newP[(i+1)%3].entry[0], newP[(i+1)%3].entry[1]);
	}
	glEnd();
}

//unsigned char color_vf[512][512][3];
void save_to_color_map()
{
	int i;
	float largest, smallest;

	/*for R*/
	smallest = largest = Object.vlist[0]->vec.entry[0];
	for(i=0; i<Object.nverts; i++)
	{
		if(Object.vlist[i]->vec.entry[0]<smallest) smallest = Object.vlist[i]->vec.entry[0];
		if(Object.vlist[i]->vec.entry[0]>largest) largest = Object.vlist[i]->vec.entry[0];
	}

	/*normalize the vector to get the red color*/
	for(i=0; i<Object.nverts; i++)
	{
		Object.vlist[i]->vf_r = 
			(Object.vlist[i]->vec.entry[0]-smallest)/(largest-smallest);
	}


	/*for G*/
	smallest = largest = Object.vlist[0]->vec.entry[1];
	for(i=0; i<Object.nverts; i++)
	{
		if(Object.vlist[i]->vec.entry[0]<smallest) smallest = Object.vlist[i]->vec.entry[1];
		if(Object.vlist[i]->vec.entry[0]>largest) largest = Object.vlist[i]->vec.entry[1];
	}

	/*normalize the vector to get the red color*/
	for(i=0; i<Object.nverts; i++)
	{
		Object.vlist[i]->vf_g = 
			(Object.vlist[i]->vec.entry[1]-smallest)/(largest-smallest);
	}

	/*for B*/
	for(i=0; i<Object.nverts; i++)
		Object.vlist[i]->vf_b = .5;
}




void CGlView::VF_to_ColorMap(GLenum mode)
{
}

////Using IBFV method to generate flow effect
void CGlView::IBFVEffect(GLenum mode)
{
	int   i, j; 
	double px, py;
	
	Face *face;
	int *verts;

	/*--------------- The 2IBFV -------------*/
	glDrawBuffer(GL_BACK);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 512, 512, 0,
		GL_RGB, GL_UNSIGNED_BYTE, f_tex);
	/*--------------- The 2IBFV -------------*/
	
	for (i=0; i<Object.nfaces; i++) {
		face = Object.flist[i];
		verts = face->verts;

		if(mode == GL_SELECT)
			glLoadName(NAMEOFTRIANGLE + i);

		glBegin(GL_POLYGON);
		for (j=0; j<face->nverts; j++) {
			glTexCoord2f(Object.vlist[verts[j]]->x, Object.vlist[verts[j]]->y);
			px = Object.vlist[verts[j]]->x + Object.vlist[verts[j]]->vec.entry[0];
			py = Object.vlist[verts[j]]->y + Object.vlist[verts[j]]->vec.entry[1];
			glVertex2f(px, py);
		}
		glEnd();
	}	
	
	iframe = iframe + 1;

	glEnable(GL_BLEND); 
	if( MoveOrStop == 0) //moving image
		glCallList(iframe % Npat + 1);
	else   //static image
		glCallList(iframe % Npat + 1 + 100);

    glBegin(GL_QUAD_STRIP);
		glTexCoord2f(0.0,  0.0);  glVertex2f(0.0, 0.0);
		glTexCoord2f(0.0,  tmax); glVertex2f(0.0, 1.0);
		glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
		glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
	glEnd();
	
	//glDisable(GL_BLEND);
	//glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 
	//				0, 0, NPIX, NPIX, 0);

	/*--------------- The 2IBFV -------------*/
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, 512, 512, GL_RGB, GL_UNSIGNED_BYTE, f_tex);

	if(MoveOrStop == 1)
	{
		/* Calculate backward texture */
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 512, 512, 0,
			GL_RGB, GL_UNSIGNED_BYTE, b_tex);
		
		for (i=0; i<Object.nfaces; i++) {
			face = Object.flist[i];
			verts = face->verts;

			if(mode == GL_SELECT)
				glLoadName(NAMEOFTRIANGLE + i);

			glBegin(GL_POLYGON);
			for (j=0; j<face->nverts; j++) {
				glTexCoord2f(Object.vlist[verts[j]]->x, Object.vlist[verts[j]]->y);
				px = Object.vlist[verts[j]]->x - Object.vlist[verts[j]]->vec.entry[0];
				py = Object.vlist[verts[j]]->y - Object.vlist[verts[j]]->vec.entry[1];
				glVertex2f(px, py);
			}
			glEnd();
		}	
		iframe = iframe + 1;
		
		glEnable(GL_BLEND); 
		if( MoveOrStop == 0) //moving image
			glCallList(iframe % Npat + 1);
		else   //static image
			glCallList(iframe % Npat + 1 + 100);

		glBegin(GL_QUAD_STRIP);
			glTexCoord2f(0.0,  0.0);  glVertex2f(0.0, 0.0);
			glTexCoord2f(0.0,  tmax); glVertex2f(0.0, 1.0);
			glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
			glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
		glEnd();
		
		glReadBuffer(GL_BACK);
		glReadPixels(0, 0, 512, 512, GL_RGB, GL_UNSIGNED_BYTE, b_tex);

		//blend two images

		for(int x = 0; x < 512; x++)
		{
			for(int y = 0; y < 512; y++)
			{
				applied_tex[x][y][0] = (int)(f_tex[x][y][0] + b_tex[x][y][0])/2.;
				applied_tex[x][y][1] = (int)(f_tex[x][y][1] + b_tex[x][y][1])/2.;
				applied_tex[x][y][2] = (int)(f_tex[x][y][2] + b_tex[x][y][2])/2.;
			}
		}
	   
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 512, 512, 0,
			GL_RGB, GL_UNSIGNED_BYTE, applied_tex);
	}

	else{
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 512, 512, 0,
			GL_RGB, GL_UNSIGNED_BYTE, f_tex);
	}

	glBegin(GL_QUAD_STRIP);
		glTexCoord2f(0.0,  0.0);  glVertex2f(0.0, 0.0);
		glTexCoord2f(0.0,  1.); glVertex2f(0.0, 1.0);
		glTexCoord2f(1., 0.0);  glVertex2f(1.0, 0.0);
		glTexCoord2f(1., 1.); glVertex2f(1.0, 1.0);
	glEnd();
}


////Select the color for the visual icons according to the type of the singularity
/*current color scheme  08/25/05
source: pure green (0, 1, 0)
repeller: light green (0, 1, 0.5)
sink:   pure red (1, 0, 0)
attractor:  orange (1, 0.5, 0)
saddle: pure blue (0, 0, 1)
center: light red (1, 0, 1)
*/

void CGlView::SetColorByType(int type)
{
	if(type == SOURCE){
		glColor4f(0.,1., 0.,1);
	}

	else if(type == SINK){
		glColor4f(1.,0., 0.,1);
	}

	else if(type == SADDLE){
		glColor4f(0.,0.2,1., 1);
	}

	else if(type == CWCENTER){
		glColor4f(1.,0., 1.,1);
	}

	else if(type == CCWCENTER){
		glColor4f(0.3, 1., 1.,1);
	}

	else if(type == AFOCUS){
		//glColor4f(1., 1., 0., 1);
		//glColor4f(1.,0.5, 0, 1);
		glColor4f(1.,0, 0, 1);
	}

	else if(type == RFOCUS){
		//glColor4f(0, 1, 0.7, 1);
		glColor4f(0.,1., 0, 1);
	}
}
	
void CGlView::DisplayCapturedSin(GLenum mode)
{
	int i;

	int singular_id = -1;

	for(i = 0; i < cur_singularity_index; i++)
	{
		SetColorByType(singularities[i].type);

		if(IsInCenter(singularities[i].gcx, singularities[i].gcy, singular_id))
		{
			if(PairCancelOn == 1 || LimitSingCancelOn == 1){
				if(mode == GL_SELECT)
					glLoadName(NAMEOFSINGULARITY + i);////Assign name for singularities
			}
			else{
				if(mode == GL_SELECT)
					glLoadName(singular_id);  ////using the ID of the element as its name (just the same as NAMEOFSINGELEM)
			}

		    DrawSolidCircle(singularities[i].gcx, singularities[i].gcy);
			glColor3f(0, 0, 0);
			draw_hollow_circle(singularities[i].gcx, singularities[i].gcy);
		}

		else{
			if((PairCancelOn == 1 || LimitSingCancelOn == 1) && mode == GL_SELECT)
				glLoadName(NAMEOFSINGULARITY + i);

            //DrawMarkTriangle(singularities[i].gcx, singularities[i].gcy);
		    DrawSolidCircle(singularities[i].gcx, singularities[i].gcy);

			glColor3f(0, 0, 0);
			draw_hollow_circle(singularities[i].gcx, singularities[i].gcy);
		}
	}
}




void CGlView::DisplayAllCapturedSin(GLenum mode)
{
	int i;

	for(i = 0; i < cur_singularity_index; i++)
	{
		SetColorByType(singularities[i].type);

		if(mode == GL_SELECT)
			glLoadName(NAMEOFSINGULARITY + i);

        DrawMarkTriangle(singularities[i].gcx, singularities[i].gcy);
	}
}

////Display the 9 controlling points for selected singular element
void CGlView::DisplaySingularControlPoints(GLenum mode, int SelectedElem)
{
	glColor3f(1, 1, 1);
	////The first four control points locate at the vertices of the editbox
	////They are used to perform uniform scaling
	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+1);  ////name for control point on p1(low left point)
	DrawSolidCircle(singularelem[SelectedElem].cur_editbox.p1.entry[0],\
		singularelem[SelectedElem].cur_editbox.p1.entry[1]);

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+2);  ////name for control point on p2(upper left point)
	DrawSolidCircle(singularelem[SelectedElem].cur_editbox.p2.entry[0],\
		singularelem[SelectedElem].cur_editbox.p2.entry[1]);

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+3);  ////name for control point on p3(upper right point)
	DrawSolidCircle(singularelem[SelectedElem].cur_editbox.p3.entry[0],\
		singularelem[SelectedElem].cur_editbox.p3.entry[1]);

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+4);  ////name for control point on p4(low right point)
	DrawSolidCircle(singularelem[SelectedElem].cur_editbox.p4.entry[0],\
		singularelem[SelectedElem].cur_editbox.p4.entry[1]);

	////The following 4 points locate at the edge of edit box
	////They are used to perform non-uniform scaling
	double x, y;

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+5);  ////name for control point on p1p2 (left)
	x = (singularelem[SelectedElem].cur_editbox.p1.entry[0]+singularelem[SelectedElem].cur_editbox.p2.entry[0])/2;
	y = (singularelem[SelectedElem].cur_editbox.p1.entry[1]+singularelem[SelectedElem].cur_editbox.p2.entry[1])/2;
	DrawSolidCircle(x, y);

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+6);  ////name for control point on p2p3 (upper)
	x = (singularelem[SelectedElem].cur_editbox.p2.entry[0]+singularelem[SelectedElem].cur_editbox.p3.entry[0])/2;
	y = (singularelem[SelectedElem].cur_editbox.p2.entry[1]+singularelem[SelectedElem].cur_editbox.p3.entry[1])/2;
	DrawSolidCircle(x, y);

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+7);  ////name for control point on p3p4 (right)
	x = (singularelem[SelectedElem].cur_editbox.p3.entry[0]+singularelem[SelectedElem].cur_editbox.p4.entry[0])/2;
	y = (singularelem[SelectedElem].cur_editbox.p3.entry[1]+singularelem[SelectedElem].cur_editbox.p4.entry[1])/2;
	DrawSolidCircle(x, y);

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+8);  ////name for control point on p4p1 (bottom)
	x = (singularelem[SelectedElem].cur_editbox.p4.entry[0]+singularelem[SelectedElem].cur_editbox.p1.entry[0])/2;
	y = (singularelem[SelectedElem].cur_editbox.p4.entry[1]+singularelem[SelectedElem].cur_editbox.p1.entry[1])/2;
	DrawSolidCircle(x, y);

	////The following point is for rotation control

	////draw a small line segment between rotation controling point and upper edge controling point
	x = (singularelem[SelectedElem].cur_editbox.p2.entry[0]+singularelem[SelectedElem].cur_editbox.p3.entry[0])/2;
	y = (singularelem[SelectedElem].cur_editbox.p2.entry[1]+singularelem[SelectedElem].cur_editbox.p3.entry[1])/2;
	double ux = singularelem[SelectedElem].cur_editbox.Up.entry[0];
	double uy = singularelem[SelectedElem].cur_editbox.Up.entry[1];
	glBegin(GL_LINES);
	glVertex2f(x, y);
	glVertex2f(ux, uy);
	glEnd();

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+9);
	DrawSolidCircle(ux, uy);
}


////Display the 9 controlling points for selected singular element
void CGlView::display_singularElem_controlPts_tensor(GLenum mode, int SelectedElem)
{
	//glColor3f(1, 1, 1);
	////The first four control points locate at the vertices of the editbox
	////They are used to perform uniform scaling

	if(ten_designelems[SelectedElem].type == 0) /*wedge*/
		glColor3f(1, 0, 0);
	else if(ten_designelems[SelectedElem].type == 1) /*trisector*/
		glColor3f(0, 1, 0);
	else if(ten_designelems[SelectedElem].type == 2) /*node*/
		glColor3f(1, 1, 0);
	else if(ten_designelems[SelectedElem].type == 3) /*center*/
		glColor3f(1, 0, 1);
	else if(ten_designelems[SelectedElem].type == 4) /*saddle*/
		glColor3f(0, 0, 1);

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+1);  ////name for control point on p1(low left point)
	DrawSolidCircle_size(ten_designelems[SelectedElem].cur_editbox.p1.entry[0],\
		ten_designelems[SelectedElem].cur_editbox.p1.entry[1], 0.006/zoom_factor);

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+2);  ////name for control point on p2(upper left point)
	DrawSolidCircle_size(ten_designelems[SelectedElem].cur_editbox.p2.entry[0],\
		ten_designelems[SelectedElem].cur_editbox.p2.entry[1], 0.006/zoom_factor);

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+3);  ////name for control point on p3(upper right point)
	DrawSolidCircle_size(ten_designelems[SelectedElem].cur_editbox.p3.entry[0],\
		ten_designelems[SelectedElem].cur_editbox.p3.entry[1], 0.006/zoom_factor);

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+4);  ////name for control point on p4(low right point)
	DrawSolidCircle_size(ten_designelems[SelectedElem].cur_editbox.p4.entry[0],\
		ten_designelems[SelectedElem].cur_editbox.p4.entry[1], 0.006/zoom_factor);

	////The following 4 points locate at the edge of edit box
	////They are used to perform non-uniform scaling
	double x, y;

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+5);  ////name for control point on p1p2 (left)
	x = (ten_designelems[SelectedElem].cur_editbox.p1.entry[0]
		+ten_designelems[SelectedElem].cur_editbox.p2.entry[0])/2;
	y = (ten_designelems[SelectedElem].cur_editbox.p1.entry[1]
		+ten_designelems[SelectedElem].cur_editbox.p2.entry[1])/2;
	DrawSolidCircle_size(x, y, 0.007/zoom_factor);

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+6);  ////name for control point on p2p3 (upper)
	x = (ten_designelems[SelectedElem].cur_editbox.p2.entry[0]
		+ten_designelems[SelectedElem].cur_editbox.p3.entry[0])/2;
	y = (ten_designelems[SelectedElem].cur_editbox.p2.entry[1]
		+ten_designelems[SelectedElem].cur_editbox.p3.entry[1])/2;
	DrawSolidCircle_size(x, y, 0.006/zoom_factor);

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+7);  ////name for control point on p3p4 (right)
	x = (ten_designelems[SelectedElem].cur_editbox.p3.entry[0]
		+ten_designelems[SelectedElem].cur_editbox.p4.entry[0])/2;
	y = (ten_designelems[SelectedElem].cur_editbox.p3.entry[1]
		+ten_designelems[SelectedElem].cur_editbox.p4.entry[1])/2;
	DrawSolidCircle_size(x, y, 0.006/zoom_factor);

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+8);  ////name for control point on p4p1 (bottom)
	x = (ten_designelems[SelectedElem].cur_editbox.p4.entry[0]
		+ten_designelems[SelectedElem].cur_editbox.p1.entry[0])/2;
	y = (ten_designelems[SelectedElem].cur_editbox.p4.entry[1]
		+ten_designelems[SelectedElem].cur_editbox.p1.entry[1])/2;
	DrawSolidCircle_size(x, y, 0.006/zoom_factor);

	////The following point is for rotation control

	////draw a small line segment between rotation controling point and upper edge controling point
	x = (ten_designelems[SelectedElem].cur_editbox.p2.entry[0]
		+ten_designelems[SelectedElem].cur_editbox.p3.entry[0])/2;
	y = (ten_designelems[SelectedElem].cur_editbox.p2.entry[1]
		+ten_designelems[SelectedElem].cur_editbox.p3.entry[1])/2;
	double ux = ten_designelems[SelectedElem].cur_editbox.Up.entry[0];
	double uy = ten_designelems[SelectedElem].cur_editbox.Up.entry[1];
	glBegin(GL_LINES);
	glVertex2f(x, y);
	glVertex2f(ux, uy);
	glEnd();

	if(mode == GL_SELECT)
		glLoadName(NAMEOFSINGCONTROL+9);
	DrawSolidCircle_size(ux, uy, 0.007/zoom_factor);
}


////Display the 2 controlling points for selected regular element
////The controling points should be tranformed with the arrow
void CGlView::DisplayRegularControlPoints(GLenum mode, int SelectedElem)
{
	glColor3f(1, 1, 1);

	if(mode == GL_SELECT)
		glLoadName(NAMEOFREGCONTROL+1);   ////control point at the base
	DrawSolidCircle(regularelem[SelectedElem].base[0], regularelem[SelectedElem].base[1]);

    ////Transform the rotation controling point according to the arrow transformation matrix
	double tempdirx, tempdiry, newdirx, newdiry;
	double rot_ang = regularelem[SelectedElem].rotang;

	tempdirx = regularelem[SelectedElem].s * regularelem[SelectedElem].Direct.entry[0];
	tempdiry = regularelem[SelectedElem].s * regularelem[SelectedElem].Direct.entry[1];

	newdirx = cos(rot_ang)*tempdirx - sin(rot_ang)*tempdiry;
	newdiry = sin(rot_ang)*tempdirx + cos(rot_ang)*tempdiry;

	if(mode == GL_SELECT)
		glLoadName(NAMEOFREGCONTROL+2);   ////control point at the head
	//DrawSolidCircle(regularelem[SelectedElem].base[0]+ARROWSCALE*regularelem[SelectedElem].Direct.entry[0]/RegularStrength, \
	//	regularelem[SelectedElem].base[1]+ARROWSCALE*regularelem[SelectedElem].Direct.entry[1]/RegularStrength);
	DrawSolidCircle(regularelem[SelectedElem].base[0]+ARROWSCALE*newdirx/RegularStrength, \
		regularelem[SelectedElem].base[1]+ARROWSCALE*newdiry/RegularStrength);

}


////display the regular elements using arrows
void CGlView::DisplayRegularElemt(GLenum mode)
{
	for(int i = 0; i < MaxNumRegularElems; i++)
	{
	   if(regularelem[i].ID>0)
	   {
		   ////Display the arrow
		   if(mode == GL_SELECT)
			   glLoadName(regularelem[i].ID);

		   	////perform the user defined transformation for editing
		   glPushMatrix();
           glTranslatef(regularelem[i].base[0], regularelem[i].base[1], 0);
		   glRotatef(360*regularelem[i].rotang/(2*M_PI), 0, 0, 1);
		   glScalef(regularelem[i].s, regularelem[i].s, 1);
		   glTranslatef(-regularelem[i].base[0], -regularelem[i].base[1], 0);

		   if(regularelem[i].type == 0) ////basic regular element
		       glColor3f(0, 1, 1);
		   else if(regularelem[i].type == 1) ////convergent element
		       //glColor3f(1, 0.5, 0);
		       glColor3f(1, 1, 0);       //for the paper
		   else                         ////divergent element
		       glColor3f(0, 1, 0.5);

		   DrawArrow(regularelem[i].base, regularelem[i].Direct.entry);

		   glPopMatrix();
	   }
	}
}


////display the edit box for all elememts when under editing mode
////display the control points for the chosen element
void CGlView::DisplayEditBox(GLenum mode)
{
	int i;

	for(i = 0; i < MaxNumSingularElems; i++)
	{
		if(singularelem[i].ID >= 0 && !singularelem[i].deleted)
		{
			SetColorByType(singularelem[i].type);

			DrawEditBox(singularelem[i].cur_editbox.p1,\
				singularelem[i].cur_editbox.p2,\
				singularelem[i].cur_editbox.p3,\
				singularelem[i].cur_editbox.p4);

			DrawEditBox_back(singularelem[i].cur_editbox.p1,\
				singularelem[i].cur_editbox.p2,\
				singularelem[i].cur_editbox.p3,\
				singularelem[i].cur_editbox.p4);
		}
	}

}


////Display the trajectory
void CGlView::DisplayTrajectory()
{
	int i;

	for(i = 0; i <= cur_traj_index; i++)
	{
		DrawSingleTrajectory(i);
	}

	////Just display current trajectory
		
	//DrawSingleTrajectory(cur_traj_index);

}


void CGlView::DisplaySmoothRegion()
{
	int i;
	//glColor3f(1, 1, 1);
		glColor3f(0.5, 0, 0.8);

	//if(sharedvars.SelStreetRegToEditOn)
	//	glColor3f(1, 0, 0.8);

	glLineWidth(1.5);

	for(i = 0; i < Num_SmoothRegionpoints; i++)
	{
		glBegin(GL_LINES);
		glVertex2f(point[i].x, point[i].y);
		glVertex2f(point[(i+1)%(Num_SmoothRegionpoints)].x, point[(i+1)%(Num_SmoothRegionpoints)].y);
		glEnd();
	}
	glLineWidth(1.);
}


void CGlView::DisplayShapeControlPts(GLenum mode)
{
	int i;

	glColor3f(1, 1, 1);
	////Draw the control hull

	
	if(num_shapecontrol_pts > 1) 
	{
		glBegin(GL_LINE_STRIP);
		for(i = 0; i < num_shapecontrol_pts; i++)
		{
			glVertex2f(control_pts[i].x, control_pts[i].y);
		}
		glEnd();
	}

	//Draw the line segment between the last and first control points
	//glBegin(GL_LINE);
	//glVertex2f(control_pts[num_shapecontrol_pts-1].x, control_pts[num_shapecontrol_pts-1].y);
	//glVertex2f(control_pts[0].x, control_pts[0].y);
	//glEnd();

	////Draw control points
	if(sharedvars.BrushInterfaceOn && num_shapecontrol_pts<=1) return;

	glColor3f(1, 0, 0);
	glPointSize(5.);
	for(i = 0; i < num_shapecontrol_pts; i++)
	{
		if(mode == GL_SELECT)
			glLoadName(NAMEOFSHAPECONTROL+i);

		//DrawControlRectangle(control_pts[i].x, control_pts[i].y);

		glBegin(GL_POINTS);
		glVertex2f(control_pts[i].x, control_pts[i].y);
		glEnd();
	}
	glEnd();
}


void CGlView::DisplayOldDesignCellCycle()
{
	int i, j;
	Face *face;
	Vertex *cur_v;

	glColor3f(1, 0, 1);

    glLineWidth(1.);
    for(i = 0; i < num_innertriangles; i++)
	{
		face = Object.flist[InnerTriangles[i]];
		 
		glBegin(GL_LINE_LOOP);
		for(j = 0; j < 3; j++)
		{
			cur_v = Object.vlist[face->verts[j]];
			glVertex2f(cur_v->x, cur_v->y);
		}
		glEnd();
	}
	glLineWidth(1.);
}



extern Edge **regionedge;                   ////mesh edges of user selected region
extern int Num_edges;

void CGlView::DisplayDesignCurve()
{
	int i;

	glColor3f(1, 0, 0);
	glBegin(GL_LINE_STRIP);
		for(i = 0; i < num_curvepts_output; i++)
			glVertex2f(out_pts[i].x, out_pts[i].y);
	glEnd();

	//glColor3f(1, 0, 1);
	//glPointSize(5.);
	//glBegin(GL_POINTS);
	//for(i = 0; i < resolution; i++)
	//	glVertex2f(out_pts[i].x, out_pts[i].y);
	//glEnd();
	//glPointSize(1.);

	//Edge *cur_e;
	//Vertex *v;
	//glBegin(GL_LINES);
	//for(i = 0 ; i < Num_edges; i++)
	//{
	//	cur_e = regionedge[i];
	//	v = Object.vlist[cur_e->verts[0]];
	//	glVertex2f(v->x, v->y);
	//	
	//	v = Object.vlist[cur_e->verts[1]];
	//	glVertex2f(v->x, v->y);
	//}
	//glEnd();
}


////Display a particle
void CGlView::DisplayParticle()
{
	/*--------------------------------------------------*/
	GLfloat ambient[] = { 1.0, 0.6, 0, 1.0 };
	GLfloat diffuse[] = { 1.0, 1.0, 0.0, 1.0 };
	GLfloat specular[] = { 1.0, 1.0, 1.0, 1.0 };
    int shiny = 300;

	GLfloat light_ambient0[] = { 1., 0.4, 0, 1.0 };
	GLfloat light_diffuse0[] = { 1., 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_diffuse2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);

	light_position[0] = 2.1;
	light_position[1] = 0.0;
	light_position[2] = 1.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -1.1;
	light_position[1] = 0.5;
	light_position[2] = 1.0;
	glLightfv(GL_LIGHT1, GL_POSITION, light_position);

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
	glMaterialf(GL_FRONT, GL_SHININESS, (GLfloat)shiny);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 

	glOrtho(0, 1,  0, 1,  0, 50);

	glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	glShadeModel(GL_SMOOTH);

	/*--------------------------------------------------*/

	icVector2 cur_vec;
	double mag;
	//getVector(px, py, cur_vec, mag);
	//px += cur_vec.entry[0];
	//py += cur_vec.entry[1];

	//Modified at 5/22/06
	particle_triangle = TraceParticleforOneStep(px, py, particle_triangle, 0);

 	glColor3f(0.8, 0.8, 0.8);
   
	glPushMatrix();
	glTranslatef(px, py, 0.);
	glutSolidSphere(0.01, 40, 40 );
	glPopMatrix();

	/*--------------------------------------------------*/
    //InitGL();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity(); 
	glTranslatef(-1.0, -1.0, -1.0); 
	glScalef(2.0, 2.0, 1.0);
	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_DEPTH_TEST);
	/*--------------------------------------------------*/
}

void CGlView::DisplayColorPlots()
{
	int   i, j; 
	
	Face *face;
	int *verts;

	hsv[1] = hsv[2] = 1.;
				
	glShadeModel (GL_SMOOTH);

	for (i=0; i<Object.nfaces; i++) {
		face = Object.flist[i];
		verts = face->verts;

		glBegin(GL_POLYGON);
		for (j=0; j<face->nverts; j++) {
			hsv[0] = 240 - 240 * (Object.vlist[verts[j]]->mag_speed - min_mag)/(max_mag - min_mag);

			HsvRgb(hsv, rgb);
			glColor4f(rgb[0], rgb[1], rgb[2], 0.4);
			glVertex2f(Object.vlist[verts[j]]->x, Object.vlist[verts[j]]->y);
		}
		glEnd();
	}
}


/////////

////////
///Testing codes 12/27/05
void CGlView::DisplayNewTriangleStrip()
{
	int i;
	Edge *cur_e;
	Vertex *cur_v;

	glColor3f(1, 1, 1);
	glLineWidth(1.5);
    glBegin(GL_LINES);
	for(i = 0; i < num_cycleedges; i++)
	{
		cur_e = Cycle_edge[i];
		cur_v = Object.vlist[cur_e->verts[0]];
		glVertex2f(cur_v->x, cur_v->y);

		cur_v = Object.vlist[cur_e->verts[1]];
		glVertex2f(cur_v->x, cur_v->y);
	}
	glEnd();
	glLineWidth(1.0);

	////Display the normals on the boundary
	Vertex *vert = NULL;
	glColor3f(0, 1, 1);
	//for(i = 0; i < num_cycleedges; i++)
	//{
	//	cur_e = Cycle_edge[i];

	//	vert = Object.vlist[cur_e->verts[0]];
	//	glPushMatrix();
	//	glTranslatef(vert->x, vert->y, 0);
	//	glRotatef(atan2(cur_e->normal.entry[1],cur_e->normal.entry[0])*360/(2*M_PI), 0, 0, 1);
	//	glScalef(ARROWSCALE/4, ARROWSCALE/4, 1);
	//		glBegin(GL_LINES);
	//		glVertex2f(0, 0);
	//		glVertex2f(1, 0);
	//		glEnd();

	//		////Draw the wings of the arrow
	//		glBegin(GL_LINES);
	//		glVertex2f(1, 0);
	//		glVertex2f(0.8, 0.16);

	//		glVertex2f(1, 0);
	//		glVertex2f(0.8, -0.16);
	//		glEnd();
	//	glPopMatrix();
	//}

	////Display the new vectors of the vertices on the boundary
	for(i = 0; i < num_cycleedges; i++)
	{
		cur_e = Cycle_edge[i];

		vert = Object.vlist[cur_e->verts[0]];
		glPushMatrix();
		glTranslatef(vert->x, vert->y, 0);
		glRotatef(atan2(vert->vec.entry[1],vert->vec.entry[0])*360/(2*M_PI), 0, 0, 1);
		glScalef(ARROWSCALE/4, ARROWSCALE/4, 1);
			glBegin(GL_LINES);
			glVertex2f(0, 0);
			glVertex2f(1, 0);
			glEnd();

			////Draw the wings of the arrow
			glBegin(GL_LINES);
			glVertex2f(1, 0);
			glVertex2f(0.8, 0.16);

			glVertex2f(1, 0);
			glVertex2f(0.8, -0.16);
			glEnd();
		glPopMatrix();
		
		
		vert = Object.vlist[cur_e->verts[1]];
		glPushMatrix();
		glTranslatef(vert->x, vert->y, 0);
		glRotatef(atan2(vert->vec.entry[1],vert->vec.entry[0])*360/(2*M_PI), 0, 0, 1);
		glScalef(ARROWSCALE/4, ARROWSCALE/4, 1);
			glBegin(GL_LINES);
			glVertex2f(0, 0);
			glVertex2f(1, 0);
			glEnd();

			////Draw the wings of the arrow
			glBegin(GL_LINES);
			glVertex2f(1, 0);
			glVertex2f(0.8, 0.16);

			glVertex2f(1, 0);
			glVertex2f(0.8, -0.16);
			glEnd();
		glPopMatrix();
	}
}




/*--------------------------------------------------------------*/
//Visual icons drawing

////Judge whether the captured singularity is an element
////Then tell program using different icons to mark element
bool CGlView::IsInCenter(double x, double y, int &singular_id)
{
	int i;
	for( i = 0; i < MaxNumSingularElems; i++)
	{
		if(singularelem[i].ID > 0 && !singularelem[i].deleted){
			if((fabs((double)(singularelem[i].centerx - x)) < 0.015)
				&&(fabs((double)(singularelem[i].centery - y)) < 0.015))
			{
				singular_id = singularelem[i].ID;
				return true;
			}
		}
	}
	if( i >= MaxNumSingularElems) return false;
}

/*****************************************************************
Draw an small triangle for captured singularities
*****************************************************************/	
void CGlView::DrawMarkTriangle(double cx, double cy)
{
	glBegin(GL_TRIANGLES);
		glVertex2f(cx-0.01, cy-0.01);
		glVertex2f(cx, cy+0.01);
		glVertex2f(cx+0.01, cy-0.01);
	glEnd();
}

void CGlView::DrawControlRectangle(double cx, double cy)
{
	glBegin(GL_POLYGON);
		glVertex2f(cx-0.01/zoom_factor, cy-0.01/zoom_factor);
		glVertex2f(cx-0.01/zoom_factor, cy+0.01/zoom_factor);
		glVertex2f(cx+0.01/zoom_factor, cy+0.01/zoom_factor);
		glVertex2f(cx+0.01/zoom_factor, cy-0.01/zoom_factor);
	glEnd();
}


/*****************************************************************
Draw a solid circle in 2D plane for center type singularities
*****************************************************************/
void CGlView::DrawSolidCircle(double cx, double cy)
{
	int i;
	//double R = 0.0001;
	double R = 0.008/zoom_factor;
	double theta, deta ;
	deta = 2 * M_PI/49.;
	double x, y;
	theta = 0.;
	glBegin(GL_POLYGON);
	for(i = 0; i < 50; i++, theta += deta)
	{
		x =  cx + R * cos(theta);
		y =  cy + R * sin(theta);

		glVertex2f(x, y);
	}
	glEnd();
}

/*****************************************************************
Draw a hollow circle in 2D plane for center type singularities
*****************************************************************/
void CGlView::draw_hollow_circle(double cx, double cy)
{
	int i;
	//double R = 0.0001;
	double R = 0.0085;
	double theta, deta ;
	deta = 2 * M_PI/49.;
	double x, y;
	theta = 0.;
	glBegin(GL_LINE_LOOP);
	for(i = 0; i < 50; i++, theta += deta)
	{
		x =  cx + R * cos(theta);
		y =  cy + R * sin(theta);

		glVertex2f(x, y);
	}
	glEnd();
}

	
void CGlView::draw_hollow_circle_size(double cx, double cy, double R)
{
	int i;
	double theta, deta ;
	deta = 2 * M_PI/49.;
	double x, y;
	theta = 0.;
	glBegin(GL_LINE_LOOP);
	for(i = 0; i < 50; i++, theta += deta)
	{
		x =  cx + R * cos(theta);
		y =  cy + R * sin(theta);

		glVertex2f(x, y);
	}
	glEnd();
}


void CGlView::DrawSolidCircle_size(double cx, double cy, double R)
{
	int i;
	double theta, deta ;
	deta = 2 * M_PI/49.;
	double x, y;
	theta = 0.;
	glBegin(GL_POLYGON);
	for(i = 0; i < 50; i++, theta += deta)
	{
		x =  cx + R * cos(theta);
		y =  cy + R * sin(theta);

		glVertex2f(x, y);
	}
	glEnd();
}

void CGlView::DrawSolidRect_size(double cx, double cy, double R)
{
	glBegin(GL_POLYGON);
	glVertex2f(cx-R, cy-R);
	glVertex2f(cx-R, cy+R);
	glVertex2f(cx+R, cy+R);
	glVertex2f(cx+R, cy-R);
	glEnd();
}


void CGlView::draw_hollow_rect_size(double cx, double cy, double R)
{
	glBegin(GL_LINE_LOOP);
	glVertex2f(cx-R, cy-R);
	glVertex2f(cx-R, cy+R);
	glVertex2f(cx+R, cy+R);
	glVertex2f(cx+R, cy-R);
	glEnd();
}

void CGlView::draw_rect_size(double cx, double cy, double R)
{
	glBegin(GL_POLYGON);
	glVertex2f(cx-R, cy-R);
	glVertex2f(cx-R, cy+R);
	glVertex2f(cx+R, cy+R);
	glVertex2f(cx+R, cy-R);
	glEnd();
}


void CGlView::DrawAParticle(double cx, double cy)
{
	int i;
	double R = 0.008;
	double theta, deta ;
	deta = 2 * M_PI/49.;
	double x, y;
	theta = 0.;
	glBegin(GL_POLYGON);
	for(i = 0; i < 50; i++, theta += deta)
	{
		x =  cx + R * cos(theta);
		y =  cy + R * sin(theta);

		glVertex2f(x, y);
	}
	glEnd();
}

/*****************************************************************
Draw an arrow in 2D plane to plot the regular element
*****************************************************************/	
void CGlView::DrawArrow(double base[2], double Direc[2])
{
    glPushMatrix();
	glTranslatef(base[0], base[1], 0);
	glRotatef(atan2(Direc[1],Direc[0])*360/(2*M_PI), 0, 0, 1);
	glScalef(ARROWSCALE, ARROWSCALE, 1);
	DrawUnitArrow();
	glPopMatrix();
}

/*****************************************************************
Draw an unit arrow (1, 0) in 2D plane to plot the regular element
*****************************************************************/	
void CGlView::DrawUnitArrow()
{
	glLineWidth(3.0);
	glBegin(GL_LINES);
	glVertex2f(0, 0);
	glVertex2f(1, 0);
	glEnd();

	////Draw the wings of the arrow
	glBegin(GL_LINES);
	glVertex2f(1, 0);
	glVertex2f(0.8, 0.16);

	glVertex2f(1, 0);
	glVertex2f(0.8, -0.16);
	glEnd();

	glLineWidth(1.);
}

/******************************************************************
Draw the basic edit box for element
******************************************************************/
void CGlView::DrawEditBox(icVector2 p1, icVector2 p2, icVector2 p3, icVector2 p4)
{
	glBegin(GL_LINE_LOOP);
		glVertex2f(p1.entry[0], p1.entry[1]);
		glVertex2f(p2.entry[0], p2.entry[1]);
		glVertex2f(p3.entry[0], p3.entry[1]);
		glVertex2f(p4.entry[0], p4.entry[1]);
    glEnd();
}

/******************************************************************
Draw the basic edit box for element editing
******************************************************************/
void CGlView::DrawEditBox_back(icVector2 p1, icVector2 p2, icVector2 p3, icVector2 p4)
{
	glEnable(GL_BLEND);
	glColor4f(1,1,1,0);
	glBegin(GL_POLYGON);
		glVertex2f(p1.entry[0], p1.entry[1]);
		glVertex2f(p2.entry[0], p2.entry[1]);
		glVertex2f(p3.entry[0], p3.entry[1]);
		glVertex2f(p4.entry[0], p4.entry[1]);
    glEnd();
	glDisable(GL_BLEND);
}


/******************************************************************
Draw specific trajectory according to its index
******************************************************************/

void CGlView::DrawSingleTrajectory(int index)
{ 
	int numlinesegs = num_linesegs_curtraj[index];
	
	//glDepthFunc(GL_LEQUAL);
	for(int i = 0; i < numlinesegs; i++)
	{
		glBegin(GL_LINES);
		glVertex2f(trajectories[index][i].gstart[0],trajectories[index][i].gstart[1]);
		glVertex2f(trajectories[index][i].gend[0],trajectories[index][i].gend[1]);
		glEnd();
	}
	//glDepthFunc(GL_LESS);

	//int numlinesegs = trajectories2[index].num_linesegs;
	//
	//for(int i = 0; i < numlinesegs; i++)
	//{
	//	glBegin(GL_LINES);
	//	glVertex2f(trajectories2[index].line_segs[i].gstart[0],trajectories2[index].line_segs[i].gstart[1]);
	//	glVertex2f(trajectories2[index].line_segs[i].gend[0],trajectories2[index].line_segs[i].gend[1]);
	//	glEnd();
	//}

	////Testing codes 12/29/05
	////Draw the last triangle for each trajectory
	//if(numlinesegs == 0)
	//	return;

	//int i;
	//Face *cur_f = Object.flist[trajectories[index][numlinesegs-1].Triangle_ID];

	//glColor3f(1, 0, 0);
	//glBegin(GL_LINE_LOOP);
	//for(i = 0; i < 3; i++)
	//{
	//	glVertex2f(Object.vlist[cur_f->verts[i]]->x, Object.vlist[cur_f->verts[i]]->y);
	//}
	//glEnd();
}


////Draw all the separatrices in the field
void CGlView::DrawSeparatrices()
{
	//Draw the separatrices begin from each saddle and long their major direction

	int i;

	////For highlighting the corresponding separatrices to the selected edges in Conley graph
	////1/28/06

	glLineWidth(5.);
	glColor3f(1, 0.5, 0);
	for(i = 0; i < num_related_trajs; i++)
	{
		DrawSingleTrajectory(related_trajs[i]);
	}


	glLineWidth(3.);
	//for(i = 1; i < cur_separatrices_index && i<2; i++)/*for vis2007 saddle-saddle connection*/
	for(i = 0; i < cur_separatrices_index /*&& i<2*/; i++)
	{
		//glColor3f(1, 1, 0);
		glColor3f(1, 0, 0);
		DrawSingleTrajectory(separatrices[i].sep1);

		glColor3f(0, 1, 0);
		DrawSingleTrajectory(separatrices[i].sep2);

		//glColor3f(1, 1, 0);
		glColor3f(1, 0, 0);
		DrawSingleTrajectory(separatrices[i].sep3);
		
		glColor3f(0, 1, 0);
		DrawSingleTrajectory(separatrices[i].sep4);

	}
	glLineWidth(2.);
}


void CGlView::DrawLimitCycleLegend(double cx, double cy, double bx, double by, int type)
{
	////1. Draw the base->center connected line
	if(type == 0) ////repeller
		glColor3f(0, 1., 0.);
	else
		glColor3f(1, 0., 0.);

	glBegin(GL_LINES);
	glVertex2f(cx, cy);
	glVertex2f(bx, by);
	glEnd();

	////2. Draw the outer rectangle
	glColor3f(1, 1, 0);
	glPushMatrix();
	//glScalef(1.2, 1.2, 0);
	DrawControlRectangle(cx, cy);

	////3. Draw the inner circle
	if(type == 0) ////repeller
		glColor3f(0, 1, 0.7);
	else
		glColor3f(1, 0.5, 0);
	DrawSolidCircle(cx, cy);
	glPopMatrix();
}


void CGlView::DisplayEigenVectorForSaddle()
{
	int i;
	for(i = 0; i < cur_singularity_index; i++)
	{
		if(singularities[i].type == SADDLE)
		{
			glColor3f(1, 0, 0);
			DrawEigenVector(singularities[i].outgoing, singularities[i].gcx, singularities[i].gcy);
			
			glColor3f(0, 1, 0);
			DrawEigenVector(singularities[i].incoming, singularities[i].gcx, singularities[i].gcy);
		}
	}
}


	
void CGlView::DrawEigenVector(icVector2 vec, double cx, double cy)
{
	glBegin(GL_LINES);
	glVertex2f(cx + 0.01*vec.entry[0],  cy + 0.01*vec.entry[1]);
	glVertex2f(cx - 0.01*vec.entry[0],  cy - 0.01*vec.entry[1]);
	glEnd();
}


void CGlView::DrawHighLighted()
{
}


void CGlView::BuildHighLighted()
{
	//int i;
	//int node1, node2;

	if(num_related_edges == 0)
	{
		return;
	}

	//for(i = 0; i < num_related_edges; i++)
	//{
	//	node1 = graphedges[related_edges[i]].node_index1;

	//	if(graphnodes[node1].type == 2)  ////It is a saddle
	//	{
	//		////we need to test the separatrices of it
	//	}
	//}
}


void CGlView::DisplayAllCellCycle()
{
	int cn, i, j;
	Face *face;
	Vertex *vert;
	
	//for(cn = 0; cn < Cur_CellCycleIndex; cn++)
	//{
	//	for(i = 0; i < NumTriangleInEachCellCycle[cn]; i++)
	//	{
	//		face = Object.flist[CellCycleList[cn][i]];
	//		glBegin(GL_LINE_LOOP);
	//		for(j = 0; j < face->nverts; j++)
	//		{
	//			vert = Object.vlist[face->verts[j]];
	//			glVertex2f(vert->x, vert->y);
	//		}
	//		glEnd();
	//	}
	//}

	for(cn = 0; cn < cur_limitcycle_index; cn++)
	{
		if(limitcycles[cn].type == 0)
			glColor3f(0, 1, 0.7);
		else
			glColor3f(1, 0.5, 0);

		for(i = 0; i < limitcycles[cn].num_triangles; i++)
		{
			face = Object.flist[limitcycles[cn].cellcycle[i]];
			glBegin(GL_LINE_LOOP);
			for(j = 0; j < face->nverts; j++)
			{
				vert = Object.vlist[face->verts[j]];
				glVertex2f(vert->x, vert->y);
			}
			glEnd();
		}
	}
}



////A testing variable here 1/11/06
extern int test_aboundtriangle, second_begintriangle;

Edge *test_sel_edge = NULL;

/******************************************************************************
This routine is used to display the closed streamline of all the limit cycles
******************************************************************************/
void CGlView::DisplayLimitCycles()
{
	int i, j;
	//int numlinesegs = num_linesegs_curtraj[index];

	////We may draw the highlight
	if(picked_node > 0 && graphnodes[picked_node].LimitCycleID >=0)
	{
		int picked_limit = graphnodes[picked_node].LimitCycleID;
		glLineWidth(2.);
		glColor3f(1, 0.5, 0);
		for(i = 0; i < limitcycles[picked_limit].num_linesegs; i++)
		{
			glBegin(GL_LINES);
			glVertex2f(limitcycles[picked_limit].closed_streamline[i].gstart[0], 
				limitcycles[picked_limit].closed_streamline[i].gstart[1]);
			glVertex2f(limitcycles[picked_limit].closed_streamline[i].gend[0], 
				limitcycles[picked_limit].closed_streamline[i].gend[1]);
			glEnd();
		}
	}

	
	glLineWidth(5.);
	for(i = 0; i < cur_limitcycle_index; i++)
	{

		if(limitcycles[i].type == 0)
		    //glColor3f(0, 1, 0.7);
		    glColor3f(0, 1, 0.);
		else
			//glColor3f(1, 0.5, 0);
		    glColor3f(1, 0, 0.);

		for(j = 0; j < limitcycles[i].num_linesegs; j++)
		{
			glBegin(GL_LINES);
			glVertex2f(limitcycles[i].closed_streamline[j].gstart[0], limitcycles[i].closed_streamline[j].gstart[1]);
			glVertex2f(limitcycles[i].closed_streamline[j].gend[0], limitcycles[i].closed_streamline[j].gend[1]);
			glEnd();
		}
		
		////Highlight the fixed point
		//glColor3f(0, 0, 0);
		//glPointSize(5.);
		//glBegin(GL_POINTS);
		//	glVertex2f(limitcycles[i].closed_streamline[0].gstart[0], limitcycles[i].closed_streamline[0].gstart[1]);
		//glEnd();
		//glPointSize(1.);
	}
	glLineWidth(1.);

	/*------------------------------------------------*/
	////Testing codes 1/5/06
	//Face *face = Object.flist[1370];
	//glBegin(GL_LINE_LOOP);
	//for(i = 0; i < 3; i++)
	//{
	//	glVertex2f(Object.vlist[face->verts[i]]->x, Object.vlist[face->verts[i]]->y);
	//}
	//glEnd();
	
	//face = Object.flist[325];
	//glColor3f(1, 0, 0);
	//glBegin(GL_LINE_LOOP);
	//for(i = 0; i < 3; i++)
	//{
	//	glVertex2f(Object.vlist[face->verts[i]]->x, Object.vlist[face->verts[i]]->y);
	//}
	//glEnd();

	//glColor3f(1, 1, 0);
	//glPointSize(3.);
	//glBegin(GL_POINTS);
	//glVertex2f(MarkNextBx[0], MarkNextBx[1]);
	//glEnd();

	//if(test_aboundtriangle >= 0)
	//{
	//	Face *face = Object.flist[test_aboundtriangle];
	//	glColor3f(1, 0, 0);
	//	glBegin(GL_LINE_LOOP);
	//	for(i = 0; i < 3; i++)
	//	{
	//		glVertex2f(Object.vlist[face->verts[i]]->x, Object.vlist[face->verts[i]]->y);
	//	}
	//	glEnd();
	//	
	//}

	//if(second_begintriangle >= 0)
	//{
	//	Face *face = Object.flist[second_begintriangle];
	//	glColor3f(1, 1, 0);
	//	glBegin(GL_LINE_LOOP);
	//	for(i = 0; i < 3; i++)
	//	{
	//		glVertex2f(Object.vlist[face->verts[i]]->x, Object.vlist[face->verts[i]]->y);
	//	}
	//	glEnd();
	//}

	//if(num_celltriangle > 0) //display the cell cycle
	//{
	//	Face *face;

	//	glColor3f(1, 0, 1);
	//	for(i = 0; i < num_celltriangle; i++)
	//	{
	//		face = Object.flist[cellcycle[i]];
	//	
	//		glBegin(GL_LINE_LOOP);
	//		for(j = 0; j < 3; j++)
	//		{
	//			glVertex2f(Object.vlist[face->verts[j]]->x, Object.vlist[face->verts[j]]->y);
	//		}
	//		glEnd();
	//	}
	//}

	//if(test_sel_edge != NULL)
	//{
	//	Vertex *v;
	//	glColor3f(0, 1, 0);
	//	glBegin(GL_LINES);
	//	for(i = 0; i < 2; i++)
	//	{
	//		v = Object.vlist[test_sel_edge->verts[i]];
	//		glVertex2f(v->x, v->y);
	//	}
	//	glEnd();

	//	for(j = 0; j < 2; j++)
	//	{
	//		v = Object.vlist[test_sel_edge->verts[j]];
	//		glPushMatrix();
	//		glTranslatef(v->x, v->y, 0);
	//		glRotatef(atan2(v->vec.entry[1],v->vec.entry[0])*360/(2*M_PI), 0, 0, 1);
	//		glScalef(ARROWSCALE/2, ARROWSCALE/2, 1);
	//			glBegin(GL_LINES);
	//			glVertex2f(0, 0);
	//			glVertex2f(1, 0);
	//			glEnd();

	//			////Draw the wings of the arrow
	//			glBegin(GL_LINES);
	//			glVertex2f(1, 0);
	//			glVertex2f(0.8, 0.16);

	//			glVertex2f(1, 0);
	//			glVertex2f(0.8, -0.16);
	//			glEnd();
	//		glPopMatrix();
	//	}
	//}

}


void CGlView::DisplayLimitCycleLegends(GLenum mode)
{
	int i;

	for(i = 0; i < cur_limitcycle_index; i++)
	{
		if(mode == GL_SELECT)
			glLoadName(NAMEOFLIMITCYCLES + i);

		DrawLimitCycleLegend(limitcycles[i].legend_center[0], limitcycles[i].legend_center[1], \
			limitcycles[i].legend_base[0], limitcycles[i].legend_base[1], limitcycles[i].type);
	}
}


///********************************************************
//This routine is used to draw a single separatrix
//********************************************************/
//
//void CGlView::DrawSingSeparatrix(int SeparatrixID)
//{
//	int i;
//
//	glBegin(GL_LINE_STRIP);
//	for(i = 0; i < Cur_point_index[CurveID]-1; i++)
//	{
//		glVertex2f(curve_point_list[CurveID][i].gpx, curve_point_list[CurveID][i].gpy);
//	}
//	glEnd();
//}

/*--------------------------------------------------------------------------------*/
//create pattern for texture mapping to create flow effect visualization
//the size of the pattern is 64x64
/*--------------------------------------------------------------------------------*/
void CGlView::makePatterns(void)
{
   int lut[256];
   int phase[NPN][NPN];
   GLubyte pat[NPN][NPN][4];
   GLubyte spat[NPN][NPN][4];
   int i, j, k, t;
    
   for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
   for (i = 0; i < NPN; i++)
   for (j = 0; j < NPN; j++) phase[i][j] = rand() % 256; 

   for (k = 0; k < Npat; k++) {
      t = k*256/Npat;                           //t is used to control the animation of the image
      for (i = 0; i < NPN; i++) 
      for (j = 0; j < NPN; j++) {
          pat[i][j][0] =
          pat[i][j][1] =
          pat[i][j][2] = lut[(t + phase[i][j]) % 255];
          pat[i][j][3] = alpha;
            
		  spat[i][j][0] = 
          spat[i][j][1] = 
          spat[i][j][2] = lut[ phase[i][j] % 255];
          spat[i][j][3] = alpha;
      }
      glNewList(k + 1, GL_COMPILE);
      glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, 
                   GL_RGBA, GL_UNSIGNED_BYTE, pat);
      glEndList();

	  glNewList(k + 1 + 100, GL_COMPILE);       //This is for static image
      glTexImage2D(GL_TEXTURE_2D, 0, 4, NPN, NPN, 0, 
                   GL_RGBA, GL_UNSIGNED_BYTE, spat);
      glEndList();   
   }
}

void CGlView::getInitTex(void)
{
	glDrawBuffer(GL_BACK);

	glCallList(0);

    glBegin(GL_QUAD_STRIP);
		glTexCoord2f(0.0,  0.0);  glVertex2f(0.0, 0.0);
		glTexCoord2f(0.0,  tmax); glVertex2f(0.0, 1.0);
		glTexCoord2f(tmax, 0.0);  glVertex2f(1.0, 0.0);
		glTexCoord2f(tmax, tmax); glVertex2f(1.0, 1.0);
	glEnd();
	
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, 512, 512, GL_RGB, GL_UNSIGNED_BYTE, f_tex);

	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, 512, 512, GL_RGB, GL_UNSIGNED_BYTE, b_tex);
}


/*--------------------------------------------------------------------------------*/
////Initialize the vector field related global variables here!!!
void CGlView::InitVFVariables()
{
	int i;

    MaxNumSingularElems = 5;                  //Maximum number of singular elements
    MaxNumRegularElems = 1;                  //Maximum number of regular elements
    MaxNumSingularities = 1;                 //Maximum number of being captured singularities
    MaxNumTrajectories = 4;                  //Maximum number of possible trajectories
                                               //(it should be flexible for future pen-and-ink sketch)
	MaxNumSeparatrices = 1;                   //Maximum number of group of separatrices
    MaxNumLinesegsPerTraj = 1600;               //Maximum number of line segments for each trajectory

	MaxNumLimitCycles = 2;

	////Elements
    singularelem = (SingularElement*) malloc(sizeof(SingularElement)* MaxNumSingularElems);  //Singular elememts' list
    cur_singelem_index = 0;
	for(i = 0; i < MaxNumSingularElems; i++)
	{
		singularelem[i].ID = -1;
		singularelem[i].type = -1;
	}
    
	regularelem = (RegularElement*) malloc(sizeof(RegularElement) * MaxNumRegularElems);     //regular elememts' list
    cur_regelem_index = 0;
	for(i = 0; i < MaxNumRegularElems; i++)
	{
		regularelem[i].ID = -1;
		regularelem[i].type = -1;
	}
    
	////Singularities
	singularities = (Singularities*) malloc(sizeof(Singularities) * MaxNumSingularities);    //being captured singularites' list
    cur_singularity_index = 0;

	////Limit cycles
	limitcycles = (LimitCycle *)malloc(sizeof(LimitCycle) * MaxNumLimitCycles);
	cur_limitcycle_index = 0;
		
	////There are some memeory issues here if I set the pointers as NULL
 	for(i = 0; i < MaxNumLimitCycles; i++)
	{
		//limitcycles[i].cellcycle = NULL;
		limitcycles[i].num_triangles = 0;

		//limitcycles[i].closed_streamline = NULL;
		limitcycles[i].num_linesegs = 0;
		limitcycles[i].singularID = -1;
		limitcycles[i].type = -1;

		//limitcycles[i].connected_limitcycle = NULL;
		limitcycles[i].num_connectedcycles = 0;
		
		//limitcycles[i].connected_saddle = NULL;
		limitcycles[i].num_connectedsaddles = 0;
	}


	////Allocate of new trajectories data structure
	//trajectories2 = (Trajectory *)malloc(sizeof(Trajectory) * MaxNumTrajectories);
	//for(i = 0; i < MaxNumTrajectories; i++)
	//{
	//	trajectories2[i].cur_MaxNumLinesegs = 500;
	//	trajectories2[i].num_linesegs = 0;
	//	trajectories2[i].line_segs = (LineSeg *)malloc(sizeof(LineSeg) * trajectories2[i].cur_MaxNumLinesegs);
	//}

	////Trajectories and separatrices
	trajectories = (LineSeg**) malloc(sizeof(LineSeg*) * MaxNumTrajectories);                //trajectories' list
	for(i = 0 ; i < MaxNumTrajectories; i++)
	{
		trajectories[i] = (LineSeg*) malloc(sizeof(LineSeg) * MaxNumLinesegsPerTraj);
	}
    cur_traj_index = 0;

	num_linesegs_curtraj = (int*) malloc(sizeof(int)*MaxNumTrajectories);  //an array stores the number of line segments for corresponding trajectory
	for(i = 0; i < MaxNumTrajectories; i ++)
	{
		num_linesegs_curtraj[i] = 0;
	}

	separatrices = (Separatrices*) malloc(sizeof(Separatrices) * MaxNumSeparatrices);
	cur_separatrices_index = 0;


	////Allocate the space for the variables needed by limit cycle detection
	AllocVarforLimitCycleDetect();
	AllocBoundaryList();

	////Allocate the space for the variables needed by region smoothing
	AllocateVarforSmoothing();

	////Allocate the space for the variables needed by topology editing
	AllocateVarforTopologyEdit();
	AllocateVarforMultiRegion();

	////Allocate the space for the shape editing
	AllocShapeDesignVars();

	////For streamline based method 2/18/06
	//AllocateVarForStreamlinePlace();

	AllocForSCC();

	AllocateSampleptsList();

	alloc_quad_regionsmooth();

	////Initialize variables for transformation
	sx = sy = uniforms = 1;
    rotateAng = 0;

	RotateDegreeofField = 0;


}


/////////////////
////Initialize flags
void CGlView::InitFlag()
{
	MoveOrStop = 1;
	EditModeOn = 0;
	MoveElemOn = 0;
	SingularitiesOn = 1;
	RegularElemOn = 1;
	
	EditModeOn = 0;
	TrajectoryOn = 0;
	SeparatricesOn = 1; //turn on the separatrices
	LimitCycleOn = 1;  //show limit cycle
	SmoothOn = 0;
	DisplaySmoothRegionOn = 0;
	PickPointOn = 0;

	//ControlPointsOn = 0;
	PairCancelOn = 0;
	sing1 = sing2 = -1;                       ////Two singularities for pair cancelation
	pair_counter = 0;

	SingularityMoveOn = 0;
	source_triangle = target_triangle = oldsingularityID = -1;
	triangle_counter = 0;
	newposx = newposy = 0;

	LimitSingCancelOn = 0;
	pairsing = pairlimit = -1;
	limitsing_counter = 0;

	LimitPairCancelOn = 0;
	pairlimit_1 = pairlimit_2 = -1;
	pairlimit_counter = 0;

	LimitCycleRelocateOn = 0;
	LimitCycleDeformOn = 0;

	////Initial shape editing flags
	NewCurveOn = 0;
	ShapeControlPtsOn = 0;
	ShapeEditOn = 0;
    num_shapecontrol_pts = 0;
    num_curvepts_output = 0;
	which_shapectrpt = -1;
	FinisheCtrptsel = 0;


	bx = by = 0;
	TraceBeginTriangleID = 4068;

	px = py = 0;
	particle_triangle = 0;
	ReleaseParticleOn = 0;

	ColorPlotOn = 0;

    SelectRepellerOrAttractor = 0;  ////select repeller or attractor?
    ClearRepellandAttractList();

    ////Initialize the variables for element editing
	TransformType = 0;                      ////1 scale, 2 rotation
	which_control_point = -1;               ////which control point has been selected
	choose_ID = -1;
	SelectTriangleID = -1;
	//TraceBeginTriangleID = -1;

	////Testing variables
	ShowCancelSmoothRegion = 0;

	////Separatrix Editing
	SeparatrixEditOn = 0;

	////streamline based visualization
	StreamlineBasedOn = 0;

	ShowSCCOn = 0;

	ShowEigenVecOn = 0;

	ShowSamplePtsOn = 0;


	ShowCoveredTris = 0;

	ShowHighCurl = 0;

	ShowTestPathOn = 0;

	ShowAllEdgesOn = 0;

	ShowHighDiv = 0;

	ShowMagOn = 0; //04/10/07

	ShowDecompOn = 0; //04/12/07

	showTraceSamplingPts = 0; /*08/01/07*/

	showSampDensityOn = 0;
	showEdgeImageOn = 0;

	showMCGConnectionOn = 0;

	//showTensorOn = 0; /*tensor visualization*/

	displayRoadNetOn = false; /*display the obtain road network 10/04/2007*/

	displayDisMapOn = 0;  /*display the distance map using color coding for fast marching 10/08/2007*/

	streetNetEditOn = false;

	showTensorLineOn=false;

	showStreetGraphOn=false;

	deleteElemOn = false;

	zoom_factor = 1;
	trans_x=trans_y=0;
}


////////release the memory allocated for the object
void CGlView::finalize()
{
	int i;

	free(singularelem);
	free(regularelem);
	free(singularities);

	for(i = 0; i < MaxNumTrajectories; i++)
	{
		if(trajectories[i] != NULL)
		    free(trajectories[i]);
	}
	free(trajectories);

	free(num_linesegs_curtraj);
	free(separatrices);


	finalizeLimitCycleDetect();

	FinalizeSmoothing();
	FinalizeTopologyEdit();
	//FianlizeMultiRegion();
	FinalizeShapeDesign();
	//FinalizeStreamlinePlace();

	////we need more finalize process here 2/18/06
}

/***************************************************************
for mouse selection process
***************************************************************/
	
void CGlView::HitProcess(double ss, double st)
{
	GLuint  selectBuffer[SELECTBUFFERSIZE];
	GLint	vp[4] = {0, 0 ,1, 1};
	int hits;

	////Build the selection buffer here
	glSelectBuffer(SELECTBUFFERSIZE, selectBuffer);

   	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);
	
	if(SmoothOn == 1 || TrajectoryOn == 1)
		gluPickMatrix(ss, st, 1e-8, 1e-8, vp );  ////set smaller pick window for triangle selection

	else
		gluPickMatrix(ss, st, 0.005, 0.005, vp );  ////set a larger pick window for element selection
		
	glOrtho(0, 1,  0, 1,  0, 50);


	////If one of the element has been selected for being edited
	TransformType = 0;
	which_control_point = -1;
	which_shapectrpt = -1;

	if(SmoothOn == 1 || TrajectoryOn == 1 || SingularityMoveOn == 1)
		IBFVEffect(GL_SELECT);    ////Draw underneath mesh, allocate name for each triangle

	else if(ShapeEditOn == 1)
		DisplayShapeControlPts(GL_SELECT);
	
	else if(PairCancelOn == 1 )  //for singularity pair cancellation only
	{
		DisplayCapturedSin(GL_SELECT);
	}

	else if(LimitSingCancelOn == 1)  //for limit cycle & singularity pair cancellation
	{
		DisplayCapturedSin(GL_SELECT);
		DisplayLimitCycleLegends(GL_SELECT);
	}

	else if(LimitPairCancelOn == 1 || LimitCycleRelocateOn == 1 )  //for limit cycle pair cancellation
		DisplayLimitCycleLegends(GL_SELECT);

	else if(LimitCycleDeformOn == 1)
	{
		//DisplayLimitCycleLegends(GL_SELECT);
		DisplayShapeControlPts(GL_SELECT);
	}

	////Testing codes for separatrix shape editing
	else if(SeparatrixEditOn == 1)
	{
		DisplayCapturedSin(GL_SELECT); ////we need to find another way to allocate name for element and its control points
		DisplayShapeControlPts(GL_SELECT);
	}

	else{
		
		DisplayCapturedSin(GL_SELECT); ////we need to find another way to allocate name for element and its control points

		DisplayRegularElemt(GL_SELECT);

		if(choose_ID>0)
		{
			DisplaySingularControlPoints(GL_SELECT, choose_ID-1);
			DisplayRegularControlPoints(GL_SELECT, choose_ID-NAMEOFREGELEM);
		}

	}

	hits = glRenderMode(GL_RENDER);


	if(hits>0)
	{
		////Because this is 2D plane, we need not worry about the overlap of objects
		////It will have at most one object being selected at one click
		//int objname = selectBuffer[3];

		////New hit processing here 06/29/05
		if(hits == 1)
		{
			if(selectBuffer[3] > 0 && selectBuffer[3] < NAMEOFSINGCONTROL )
			{
				////select an element
				choose_ID = selectBuffer[3];
					
				/////12/28/05
				if(SeparatrixEditOn == 1)
				{
					Separatrix_saddle = choose_ID;  ////save the selected singularity
				}
			}

			else if(selectBuffer[3] >= NAMEOFSINGCONTROL && selectBuffer[3] < NAMEOFSINGULARITY)
			{    ////select a control point
				if(selectBuffer[3] < NAMEOFREGCONTROL) ////control point for singular element
				{
					ControlPointForSinElem(selectBuffer[3]-NAMEOFSINGCONTROL);
				}
				else{                                  ////control point for regular element
                    ControlPointForRegElem(selectBuffer[3]-NAMEOFREGCONTROL);
				}
			}

			else if(selectBuffer[3] >= NAMEOFSINGULARITY && selectBuffer[3] < NAMEOFLIMITCYCLES/*NAMEOFTRIANGLE*/)
			{
				////select a singularity
				//if(pair_counter == 0){
				//	sing1 = selectBuffer[3] - NAMEOFSINGULARITY;
				//	pair_counter++;
				//}
				//else if(pair_counter == 1){
				//	sing2 = selectBuffer[3] - NAMEOFSINGULARITY;
				//	pair_counter++;
				//}

				////More general multiple repellers and attractors cancellation
				if(SelectRepellerOrAttractor == 0 && PairCancelOn == 1)
				{
					sing1 = selectBuffer[3] - NAMEOFSINGULARITY;

					glMatrixMode(GL_PROJECTION);
					glPopMatrix();

					glMatrixMode( GL_MODELVIEW );
					AddToRepellerList(sing1, 0);
				}
				else if(SelectRepellerOrAttractor == 1 && PairCancelOn == 1)
				{
					sing2 = selectBuffer[3] - NAMEOFSINGULARITY;

					glMatrixMode(GL_PROJECTION);
					glPopMatrix();

					glMatrixMode( GL_MODELVIEW );
					AddToAttractorList(sing2, 0);
				}
				
				//// 2/06/06 for limit cycle and singularity pair cancellation
				else if(LimitSingCancelOn == 1) //user select the singularity
				{
					if(limitsing_counter == 0)
					{
						pairsing = selectBuffer[3] - NAMEOFSINGULARITY;
						limitsing_counter = 1;
					}

					else{
						////we have already select one node
						if(pairsing >= 0 )//we have already selected a singularity!
						{
							MessageBox("Please select a limit cycle!", "Error", MB_OK);
							goto LL;
						}
						else
						{
							pairsing = selectBuffer[3] - NAMEOFSINGULARITY;
							limitsing_counter = 2;
						}
					}
				}
				
				/////12/28/05
				if(SeparatrixEditOn == 1)
				{
					Separatrix_saddle = choose_ID;  ////save the selected singularity
				}


				/////12/29/05  if it is the limit cycle and singularity pair cancellation
				if(LimitSingCancelOn == 1)
				{
					if(pair_counter == 0)  // we force the first select one should be singularity
					{
						sing1 = selectBuffer[3] - NAMEOFSINGULARITY;
						pair_counter++;
					}
				}
				///////////////////////
			}

			else if(selectBuffer[3] >= NAMEOFLIMITCYCLES && selectBuffer[3] < NAMEOFSHAPECONTROL)
			{
				////SELECT A LIMIT CYCLE 08/31/05

				////if it is limit cycle pair cancellation operation
				if(LimitPairCancelOn == 1)
				{
					if(pairlimit_counter == 0){
						pairlimit_1 = selectBuffer[3] - NAMEOFLIMITCYCLES;
						pairlimit_counter++;
					}
					else if(pairlimit_counter == 1){
						pairlimit_2 = selectBuffer[3] - NAMEOFLIMITCYCLES;
						pairlimit_counter++;
					}
				}

				//// 2/06/06 for limit cycle and singularity pair cancellation
				if(LimitSingCancelOn == 1)
				{
					if(limitsing_counter == 0)
					{
						pairlimit = selectBuffer[3] - NAMEOFLIMITCYCLES;
						limitsing_counter = 1;
					}

					else{
						////we have already select one node
						if(pairlimit >= 0 )//we have already selected a singularity!
						{
							MessageBox("Please select a singularity!", "Error", MB_OK);
							goto LL;
						}
						else
						{
							pairlimit = selectBuffer[3] - NAMEOFLIMITCYCLES;
							limitsing_counter = 2;
						}
					}
				}

				if(LimitCycleRelocateOn == 1 || LimitCycleDeformOn == 1)
				{
					pairlimit = selectBuffer[3] - NAMEOFLIMITCYCLES;
				}
			}

			////Get the id of the shape control point
			else if(selectBuffer[3] >= NAMEOFSHAPECONTROL && selectBuffer[3] < NAMEOFTRIANGLE
				&& FinisheCtrptsel == 1)  ////After we finish all the selection of control points
			{ 
				if(selectBuffer[3] >= NAMEOFSHAPECONTROL)
					which_shapectrpt = selectBuffer[3] - NAMEOFSHAPECONTROL;
			}

			////Get the id of the shape control point for limit cycle deformation 3/10/06
			else if(selectBuffer[3] >= NAMEOFSHAPECONTROL && selectBuffer[3] < NAMEOFTRIANGLE
				&& LimitCycleDeformOn == 1)  ////After we finish all the selection of control points
			{ 
				if(selectBuffer[3] >= NAMEOFSHAPECONTROL)
					which_shapectrpt = selectBuffer[3] - NAMEOFSHAPECONTROL;
			}

			else{  ////select a triangle
				if(SmoothOn == 1 && PickPointOn == 1)  ////under select smoothing region mode
				{
					SelectTriangleID = selectBuffer[3] - NAMEOFTRIANGLE;
					if(showTensorOn==0)
						AddToPointList(ss, st, SelectTriangleID);
				}
				else if(TrajectoryOn == 1)             ////under trajectory drawing mode
					SelectTriangleID = selectBuffer[3] - NAMEOFTRIANGLE;

				else if(SingularityMoveOn == 1)        ////under singularity movement operation
				{
					if(triangle_counter == 0)
					{
						source_triangle = selectBuffer[3] - NAMEOFTRIANGLE;
						if(Object.flist[source_triangle]->singularity_index < 0){
							glMatrixMode(GL_PROJECTION);
							glPopMatrix();

							glMatrixMode( GL_MODELVIEW );
							
							MessageBox("Not select a legal singularity, try again!", "Error");
							
							return;
						}
						else{
							oldsingularityID = Object.flist[source_triangle]->singularity_index;
						    triangle_counter ++;
						}
					}
					else if(triangle_counter == 1){
						if(selectBuffer[3] - NAMEOFTRIANGLE > Object.nfaces
							||selectBuffer[3] - NAMEOFTRIANGLE < 0 )
						{
							glMatrixMode(GL_PROJECTION);
							glPopMatrix();

							glMatrixMode( GL_MODELVIEW );
							
							MessageBox("Not select a legal singularity, try again!", "Error");
							return;
						}
						else{
							target_triangle = selectBuffer[3] - NAMEOFTRIANGLE;
							triangle_counter ++;
							
							////record the global coordinates of the new position
							newposx = ss;
							newposy = st;
						}
					}
				}
			}
		}

		else{ ////hits>1, normally, hits == 2
			////This situation should only happen under editing mode
			int large_objname, small_objname;
			if(selectBuffer[3] > selectBuffer[7])
			{
				large_objname = selectBuffer[3];
				small_objname = selectBuffer[7];
			}
			else{
				large_objname = selectBuffer[7];
				small_objname = selectBuffer[3];
			}

			///////////////////////////////////////////////
			if(EditModeOn == 1) ////Priory to control point
			{
				if(large_objname < NAMEOFREGCONTROL) ////control point for singular element
				{
					ControlPointForSinElem(large_objname-NAMEOFSINGCONTROL);
				}
				else{                                  ////control point for regular element
                    ControlPointForRegElem(large_objname-NAMEOFREGCONTROL);
				}
			}

			else if(SmoothOn == 1 && PickPointOn == 1)  ////under select smoothing region mode
			{
				SelectTriangleID = small_objname - NAMEOFTRIANGLE;
					if(showTensorOn==0)
						AddToPointList(ss, st, SelectTriangleID);
			}

			else if(NewCurveOn == 0 && ShapeEditOn == 1)
			{
				which_shapectrpt = large_objname - NAMEOFSHAPECONTROL;
			}

			else if(TrajectoryOn == 1)             ////under trajectory drawing mode
				SelectTriangleID = small_objname - NAMEOFTRIANGLE;
			
			else if(SingularityMoveOn == 1)        ////under singularity movement operation
			{
				if(triangle_counter == 0)
				{
					source_triangle = small_objname - NAMEOFTRIANGLE;
					if(Object.flist[source_triangle]->singularity_index < 0){
						glMatrixMode(GL_PROJECTION);
						glPopMatrix();

						glMatrixMode( GL_MODELVIEW );
						
						MessageBox("Not select a legal singularity, try again!", "Error");
						
						return;
					}
					else{
						oldsingularityID = Object.flist[source_triangle]->singularity_index;
						triangle_counter ++;
					}
				}
				else if(triangle_counter == 1){
					if(small_objname - NAMEOFTRIANGLE > Object.nfaces || small_objname - NAMEOFTRIANGLE < 0){
						glMatrixMode(GL_PROJECTION);
						glPopMatrix();

						glMatrixMode( GL_MODELVIEW );
						MessageBox("Not select a legal singularity, try again!", "Error");
						return;
					}
					else{
						target_triangle = small_objname - NAMEOFTRIANGLE;
						triangle_counter ++;

						////record the global coordinates of the new position
						newposx = ss;
						newposy = st;
					}
				}
			}
		}
	}

	else{
		choose_ID = -1;
		SelectTriangleID = -1;
	}



LL:	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode( GL_MODELVIEW );

}



////Hit process for periodic orbit deform 3/30/06 
void CGlView::HitProcessforLimitDeform(double ss, double st)
{
	GLuint  selectBuffer[SELECTBUFFERSIZE];
	GLint	vp[4] = {0, 0 ,1, 1};
	int hits;

	////Build the selection buffer here
	glSelectBuffer(SELECTBUFFERSIZE, selectBuffer);

   	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);

	gluPickMatrix(ss, st, 0.005, 0.005, vp );  ////set a larger pick window for element selection
	glOrtho(0, 1,  0, 1,  0, 50);
		
	if(LimitCycleDeformOn == 1)
		DisplayShapeControlPts(GL_SELECT);

	which_shapectrpt = -1;

	hits = glRenderMode(GL_RENDER);

	if(hits>0)
	{
		if(selectBuffer[3] >= NAMEOFSHAPECONTROL && selectBuffer[3] < NAMEOFTRIANGLE )
		{
			////select a singularity
			which_shapectrpt = selectBuffer[3]-NAMEOFSHAPECONTROL;
		}

		else if(hits > 1 && selectBuffer[7] >= NAMEOFSHAPECONTROL && selectBuffer[7] < NAMEOFTRIANGLE)
		{
			which_shapectrpt = selectBuffer[7]-NAMEOFSHAPECONTROL;
		}
	}

	else{
		which_shapectrpt = -1;
	}
	
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode( GL_MODELVIEW );
}


////Hit process for limit cycle selection 3/30/06 
void CGlView::HitProcessforLimitCycleSelect(double ss, double st)
{
	GLuint  selectBuffer[SELECTBUFFERSIZE];
	GLint	vp[4] = {0, 0 ,1, 1};
	int hits;

	////Build the selection buffer here
	glSelectBuffer(SELECTBUFFERSIZE, selectBuffer);

   	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);

	gluPickMatrix(ss, st, 0.005, 0.005, vp );  ////set a larger pick window for element selection
	glOrtho(0, 1,  0, 1,  0, 50);
		
	if(LimitCycleDeformOn == 1)
		DisplayLimitCycleLegends(GL_SELECT);

	pairlimit = -1;

	hits = glRenderMode(GL_RENDER);

	if(hits>0)
	{
		if(selectBuffer[3] >= NAMEOFLIMITCYCLES && selectBuffer[3] < NAMEOFSHAPECONTROL )
		{
			////select a limit cycle
			pairlimit = selectBuffer[3]-NAMEOFLIMITCYCLES;
		}

		else if(hits > 1 && selectBuffer[7] >= NAMEOFLIMITCYCLES && selectBuffer[7] < NAMEOFSHAPECONTROL)
		{
			pairlimit = selectBuffer[7]-NAMEOFLIMITCYCLES;
		}
	}

	else{
		pairlimit = -1;
	}
	
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode( GL_MODELVIEW );
}

////////
////Hit process for Conley graph 
void CGlView::HitProcessforGraph(double ss, double st)
{
	GLuint  selectBuffer[SELECTBUFFERSIZE];
	GLint	vp[4] = {0, 0 ,1, 1};
	int hits;

	////Build the selection buffer here
	glSelectBuffer(SELECTBUFFERSIZE, selectBuffer);

   	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);

	gluPickMatrix(ss, st, 0.005, 0.005, vp );  ////set a larger pick window for element selection
	glOrtho(0, 1,  0, 1,  0, 50);
		
	DisplayAllCapturedSin(GL_SELECT);             ////
	DisplayLimitCycleLegends(GL_SELECT);

	picked_sing = -1;
	picked_limitcycle = -1;
	picked_node = -1;

	hits = glRenderMode(GL_RENDER);

	if(hits>0)
	{
		if(selectBuffer[3] >= NAMEOFSINGULARITY && selectBuffer[3] < NAMEOFLIMITCYCLES )
		{
			////select a singularity
			picked_sing = selectBuffer[3]-NAMEOFSINGULARITY;
			picked_node = singularities[picked_sing].node_index;
		}
		else if(selectBuffer[3] >= NAMEOFLIMITCYCLES && selectBuffer[3] < NAMEOFSHAPECONTROL)
		{
			////select a limit cycle
			picked_limitcycle = selectBuffer[3] - NAMEOFLIMITCYCLES;
			picked_node = limitcycles[picked_limitcycle].node_index;
		}
	}

	else{
		picked_sing = -1;
		picked_limitcycle = -1;
		picked_node = -1;
	}
	
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode( GL_MODELVIEW );
}


////Hit processing for adding new elements by region smoothing 11/30/05
void CGlView::HitProcessforSelectUnderneathMesh(double ss, double st)
{
	GLuint  selectBuffer[SELECTBUFFERSIZE];
	GLint	vp[4] = {0, 0 ,1, 1};
	int hits;

	////Build the selection buffer here
	glSelectBuffer(SELECTBUFFERSIZE, selectBuffer);

   	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);
	
	gluPickMatrix(ss, st, 1e-8, 1e-8, vp );  ////set smaller pick window for triangle selection

	glOrtho(0, 1,  0, 1,  0, 50);

	IBFVEffect(GL_SELECT);    ////Draw underneath mesh, allocate name for each triangle

	hits = glRenderMode(GL_RENDER);

	if(hits>0)
	{
		if(selectBuffer[3] >= NAMEOFTRIANGLE )
		{
			////select a singularity
			newelementtriangle = selectBuffer[3]-NAMEOFTRIANGLE;
		}
	}

	else{
		newelementtriangle = -1;
	}
	
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode( GL_MODELVIEW );
}



////Routines for multiple repellers or attractors selection
void CGlView::AddToRepellerList(int index, int singorlimitcycle)
{
	////Deal with the selection of a saddle
	if(((NumRepellers > 0)
		||(NumAttractors == 0 && NumRepellers == 1 && singularities[graphnodes[repeller_nodes[0]].singularityID].type == SADDLE))
        && singularities[index].type == SADDLE)
	{
		MessageBox("Can not select a saddle as repeller now!", "Error", MB_OK);
		return;
	}

	if(NumRepellers == 1 && singularities[graphnodes[repeller_nodes[0]].singularityID].type == SADDLE)
	{
		MessageBox("You have selected a saddle before!", "Error", MB_OK);
		return;
	}

	////Deal with other cases
	if((singularities[index].type != SOURCE && singularities[index].type != SADDLE
		&& singularities[index].type != RFOCUS && singorlimitcycle == 0)
		|| (singorlimitcycle == 1 && limitcycles[index].type != 0))
	{
		MessageBox("This is not a repeller!", "Error", MB_OK);
		return;
	}
	
	////check repeatation
	int i;
	int nodewanttoadd;

	if(singorlimitcycle == 0)
		nodewanttoadd = singularities[index].node_index;
	else
		nodewanttoadd = limitcycles[index].node_index;

	for(i = 0; i < NumRepellers; i++)
	{
		if(nodewanttoadd == repeller_nodes[i])
		{
			MessageBox("The node has been selected!", "Error", MB_OK);
			return;
		}
	}

	repeller_nodes = Extend_link(repeller_nodes, NumRepellers);
	    
	repeller_nodes[NumRepellers] = nodewanttoadd;

	NumRepellers++;
}


///We need to put the following two routines to the "topologyedit" files
void CGlView::AddToAttractorList(int index, int singorlimitcycle)
{
	////Deal with the selection of a saddle
	if(((NumAttractors > 0) 
		||(NumAttractors == 0 && NumRepellers == 1 && singularities[graphnodes[repeller_nodes[0]].singularityID].type == SADDLE))
		&& singularities[index].type == SADDLE)
	{
		MessageBox("Can not select a saddle as attractor now!", "Error", MB_OK);
		return;
	}

	if(NumAttractors == 1 && singularities[graphnodes[attractor_nodes[0]].singularityID].type == SADDLE)
	{
		MessageBox("You have selected a saddle before!", "Error", MB_OK);
		return;
	}

	////Deal with other cases
	if((singularities[index].type != SINK && singularities[index].type != SADDLE
		&& singularities[index].type != AFOCUS && singorlimitcycle == 0)
		||(singorlimitcycle == 1 && limitcycles[index].type != 1))
	{
		MessageBox("This is not a attractor!", "Error", MB_OK);
		return;
	}

	////check repeatation
	int i;
	int nodewanttoadd;

	if(singorlimitcycle == 0)
		nodewanttoadd = singularities[index].node_index;
	else
		nodewanttoadd = limitcycles[index].node_index;

	for(i = 0; i < NumAttractors; i++)
	{
		if(nodewanttoadd == attractor_nodes[i])
		{
			MessageBox("The node has been selected!", "Error", MB_OK);
			return;
		}
	}

	attractor_nodes = Extend_link(attractor_nodes, NumAttractors);

	attractor_nodes[NumAttractors] = nodewanttoadd;

	NumAttractors++;
}


void CGlView::ClearRepellandAttractList()
{
	if(repeller_nodes != NULL)
	{
		free(repeller_nodes);
		repeller_nodes = NULL;
	}
	if(attractor_nodes != NULL)
	{
		free(attractor_nodes);
		attractor_nodes = NULL;
	}

	NumRepellers = 0;
	NumAttractors = 0;
	Num_MediaNodes = 0;
}


void CGlView::ControlPointForSinElem(int Point_id)
{
	if(Point_id == LOWLEFT || Point_id == UPPERLEFT ||Point_id == UPPERRIGHT || Point_id == LOWRIGHT)
	{
		TransformType = 1; ////Uniform scaling
	}
	else if(Point_id == LEFT || Point_id == RIGHT || Point_id == UPPER || Point_id == BUTTOM)
	{
		TransformType = 2; ////non-Uniform scaling
	}
	else
	{
		TransformType = 3; ////Rotation
	}
	which_control_point = Point_id;
}

	
void CGlView::ControlPointForRegElem(int Point_id)
{
	if(Point_id + ARROWBASE - 1 == ARROWBASE)
	{
		TransformType = 1; ////scale
	}
	else
		TransformType = 3; ////Rotation

	which_control_point = Point_id + ARROWBASE - 1;
}



//////Testing functions
void CGlView::DisplayElementTriangles()
{
	int i, j/*, face_id*/;
    Face *face;
	Vertex *vert;

	glColor3f(0,1,0);
	for(i = 0; i < cur_singelem_index; i++)
	{
		face = Object.flist[singularelem[i].Triangle_ID];
		glBegin(GL_LINE_LOOP);
		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];
			glVertex2f(vert->x, vert->y);
		}
		glEnd();
	}

	////Display the litter arrows on the vertices
		//glColor3f(1, 1, 0);
		//for(i = 0; i < cur_singelem_index; i++)
		//{
		//	face = Object.flist[singularelem[i].Triangle_ID];

		//	for(j = 0; j < face->nverts; j++)
		//	{
		//		vert = Object.vlist[face->verts[j]];
		//		glPushMatrix();
		//		glTranslatef(vert->x, vert->y, 0);
		//		glRotatef(atan2(vert->vec.entry[1],vert->vec.entry[0])*360/(2*M_PI), 0, 0, 1);
		//		glScalef(ARROWSCALE/5, ARROWSCALE/5, 1);
		//			glBegin(GL_LINES);
		//			glVertex2f(0, 0);
		//			glVertex2f(1, 0);
		//			glEnd();

		//			////Draw the wings of the arrow
		//			glBegin(GL_LINES);
		//			glVertex2f(1, 0);
		//			glVertex2f(0.8, 0.16);

		//			glVertex2f(1, 0);
		//			glVertex2f(0.8, -0.16);
		//			glEnd();
		//		glPopMatrix();
		//	}
		//}
}


void CGlView::DisplaySingularTriangles()
{
	int i, j;
    Face *face;
	Vertex *vert;

	glColor3f(1,0,0);
	for(i = 0; i < cur_singularity_index; i++)
	{
		face = Object.flist[singularities[i].Triangle_ID];

		glBegin(GL_LINE_LOOP);
		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];
			glVertex2f(vert->x, vert->y);
		}
		glEnd();
	}
}

	
void CGlView::DisplayWholeMeshWithFieldArrows()
{
	int i, j;
    Face *face;
	Vertex *vert;

	glColor3f(0.8,0.8,0.8);
	for(i = 0; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		glBegin(GL_LINE_LOOP);
		for(j = 0; j < face->nverts; j++)
		{
			vert = Object.vlist[face->verts[j]];
			glVertex2f(vert->x, vert->y);
		}
		glEnd();
	}
	
	////Display the litter arrows on the vertices
		glColor3f(1, 1, 0);
		for(i = 0; i < Object.nfaces; i++)
		{
			face = Object.flist[i];

			for(j = 0; j < face->nverts; j++)
			{
				vert = Object.vlist[face->verts[j]];
				glPushMatrix();
				glTranslatef(vert->x, vert->y, 0);
				glRotatef(atan2(vert->vec.entry[1],vert->vec.entry[0])*360/(2*M_PI), 0, 0, 1);
				glScalef(ARROWSCALE/5, ARROWSCALE/5, 1);
					glBegin(GL_LINES);
					glVertex2f(0, 0);
					glVertex2f(1, 0);
					glEnd();

					////Draw the wings of the arrow
					glBegin(GL_LINES);
					glVertex2f(1, 0);
					glVertex2f(0.8, 0.16);

					glVertex2f(1, 0);
					glVertex2f(0.8, -0.16);
					glEnd();
				glPopMatrix();
			}
		}
}


/*
*/
void  HsvRgb( float hsv[3], float rgb[3] )
{
	float h, s, v;			// hue, sat, value
	float r, g, b;			// red, green, blue
	float i, f, p, q, t;		// interim values

	// guarantee valid input:
	h = hsv[0] / 60.;
	while( h >= 6. )	h -= 6.;
	while( h <  0. ) 	h += 6.;

	s = hsv[1];
	if( s < 0. )
		s = 0.;
	if( s > 1. )
		s = 1.;

	v = hsv[2];
	if( v < 0. )
		v = 0.;
	if( v > 1. )
		v = 1.;

	// if sat==0, then is a gray:
	if( s == 0.0 )
	{
		rgb[0] = rgb[1] = rgb[2] = v;
		return;
	}

	// get an rgb from the hue itself:
	i = floor( h );
	f = h - i;
	p = v * ( 1. - s );
	q = v * ( 1. - s*f );
	t = v * ( 1. - ( s * (1.-f) ) );

	switch( (int) i )
	{
		case 0:
			r = v;	g = t;	b = p;
			break;
	
		case 1:
			r = q;	g = v;	b = p;
			break;
	
		case 2:
			r = p;	g = v;	b = t;
			break;
	
		case 3:
			r = p;	g = q;	b = v;
			break;
	
		case 4:
			r = t;	g = p;	b = v;
			break;
	
		case 5:
			r = v;	g = p;	b = q;
			break;
	}

	rgb[0] = r;
	rgb[1] = g;
	rgb[2] = b;
}


void GetColor20(int num, float rgb[3])
{
	float hsv[3] = {0, 1., 1.};

	hsv[0] = 360 * (double)num/20.;

	HsvRgb(hsv,rgb);
}

void GetColor8(int num, float rgb[3])
{
	float hsv[3] = {0, 1., 1.};

	hsv[0] = 360 * (double)num/8.;

	HsvRgb(hsv,rgb);
}

void getcolor_scc_frac(int num, int frac, float rgb[3])
{
	float hsv[3] = {0, 1, 1};

	hsv[0] = 360 * (double)num/frac;

	HsvRgb(hsv,rgb);
}



void CGlView::DisplaySCC()
{
	int i, j, k;
	Face *face;
	Vertex *vert;

	int countdisplay = 0;
	for(i = 0; i < scclist.num_sccs; i++)
	{
		if(scclist.scccomponents[i].num_nodes <= 2 && scclist.scccomponents[i].num_singularities <= 0)
			//||(scclist.scccomponents[i].num_nodes > 2 && scclist.scccomponents[i].valid == 0))
			continue;

		else if(scclist.scccomponents[i].num_nodes > 2 && scclist.scccomponents[i].valid == 0)
			continue;


		else
		{
			//if(i%7 == 0)
			//	glColor4f(1, 0, 0, 0.4);
			//else if(i%7 == 1)
			//	glColor4f(0, 1, 0, 0.4);
			//else if(i%7 == 2)
			//	glColor4f(0, 0, 1, 0.4);
			//else if(i%7 == 3)
			//	glColor4f(1, 1, 0, 0.4);
			//else if(i%7 == 4)
			//	glColor4f(1, 0, 1, 0.4);
			//else if(i%7 == 5)
			//	glColor4f(0, 1, 1, 0.4);
			//else
			//	glColor4f(1, 1, 1, 0.4);

			float rgb[3] = {0.};
			//GetColor20(i%30, rgb);
			//GetColor8((countdisplay+0)%8, rgb);
			getcolor_scc_frac(countdisplay%7, 7, rgb);
			countdisplay++;
			glColor4f(rgb[0], rgb[1], rgb[2], 0.4);

			/*for the TVCG paper 08/26/2007*/
			//if(graphnodes[scclist.scccomponents[i].node_index].type == 1
			//	&& graphnodes[scclist.scccomponents[i].node_index].labelindex == 1) 
			//glColor4f(0, 1, 0.5, 0.4);

			glBegin(GL_TRIANGLES);
			for(j = 0; j < scclist.scccomponents[i].num_nodes; j++)
			{
				face = Object.flist[scclist.scccomponents[i].nodes[j]];

				for(k = 0; k < face->nverts; k++)
				{
					vert = Object.vlist[face->verts[k]];
					glVertex2f(vert->x, vert->y);
				}
			}
			glEnd();

			//glColor3f(1, 1, 1);
			//for(j = 0; j < scclist.scccomponents[i].num_nodes; j++)
			//{
			//	face = Object.flist[scclist.scccomponents[i].nodes[j]];

			//glBegin(GL_LINE_LOOP);
			//	for(k = 0; k < face->nverts; k++)
			//	{
			//		vert = Object.vlist[face->verts[k]];
			//		glVertex2f(vert->x, vert->y);
			//	}
			//glEnd();
			//}

			/**  Display mixed edges and the vectors on the two ending vertices 07/18/06 **/
			//Edge *cur_e;
			//glColor3f(0, 0, 0);
			//glBegin(GL_LINES);
			//for(j = 0; j < scclist.scccomponents[i].num_nodes; j++)
			//{
			//	face = Object.flist[scclist.scccomponents[i].nodes[j]];

			//	for(k = 0; k < 3; k++)
			//	{
			//		cur_e = face->edges[k];

			//		if(cur_e->mixed == 0 && cur_e->OnBoundary == 0)
			//			continue;

			//		if(cur_e->mixed == 2 && cur_e->OnBoundary == 0)
			//			continue;

			//		if(cur_e->mixed == 1)
			//			glColor3f(0, 0, 0);

			//		else if(cur_e->mixed == 2 && cur_e->OnBoundary == 1)
			//			glColor3f(0, 1, 0);

			//		else 
			//			glColor3f(1, 1, 0);

			//		vert = Object.vlist[cur_e->verts[0]];
			//		glVertex2f(vert->x, vert->y);

			//		vert = Object.vlist[cur_e->verts[1]];
			//		glVertex2f(vert->x, vert->y);
			//	}
			//}
			//glEnd();

			//glColor3f(1, 0, 0);
			//for(j = 0; j < scclist.scccomponents[i].num_nodes; j++)
			//{
			//	face = Object.flist[scclist.scccomponents[i].nodes[j]];

			//	for(k = 0; k < 3; k++)
			//	{
			//		cur_e = face->edges[k];

			//		if(cur_e->mixed == 0)
			//			continue;
			//		
			//		if(cur_e->mixed == 2 && cur_e->OnBoundary == 0)
			//			continue;

			//		vert = Object.vlist[cur_e->verts[0]];

			//		glPushMatrix();
			//		glTranslatef(vert->x, vert->y, 0);
			//		glRotatef(atan2(vert->vec.entry[1], vert->vec.entry[0])*180./M_PI, 0, 0, 1);
			//		glScalef(ARROWSCALE/4., ARROWSCALE/4., 1);
			//			glBegin(GL_LINES);
			//			glVertex2f(0, 0);
			//			glVertex2f(1, 0);
			//			glEnd();

			//			////Draw the wings of the arrow
			//			glBegin(GL_LINES);
			//			glVertex2f(1, 0);
			//			glVertex2f(0.8, 0.16);

			//			glVertex2f(1, 0);
			//			glVertex2f(0.8, -0.16);
			//			glEnd();
			//		glPopMatrix();
			//		
			//		vert = Object.vlist[cur_e->verts[1]];
			//		
			//		glPushMatrix();
			//		glTranslatef(vert->x, vert->y, 0);
			//		glRotatef(atan2(vert->vec.entry[1], vert->vec.entry[0])*180./M_PI, 0, 0, 1);
			//		glScalef(ARROWSCALE/4., ARROWSCALE/4., 1);
			//			glBegin(GL_LINES);
			//			glVertex2f(0, 0);
			//			glVertex2f(1, 0);
			//			glEnd();

			//			////Draw the wings of the arrow
			//			glBegin(GL_LINES);
			//			glVertex2f(1, 0);
			//			glVertex2f(0.8, 0.16);

			//			glVertex2f(1, 0);
			//			glVertex2f(0.8, -0.16);
			//			glEnd();
			//		glPopMatrix();
			//	}
			//}

		}
	}

	/*Display the highlighted SCC*/
	if(/*ShowMCGOn == 1 && */picked_node >= 0)
	{
		int *tris = scclist.scccomponents[mcgnodes[picked_node].scc_index].nodes;
		for(i = 0; i < scclist.scccomponents[mcgnodes[picked_node].scc_index].num_nodes; i++)
		{
			glColor3f(1, 1, 1);

			if(tris[i] < 0)
				continue;

			face = Object.flist[tris[i]];

			glBegin(GL_LINE_LOOP);
			for(j = 0; j < 3; j++)
			{
				glVertex2f(Object.vlist[face->verts[j]]->x, Object.vlist[face->verts[j]]->y);
			}
			glEnd();
		}
	}
}


extern double max_div, min_div, max_curl, min_curl;
extern double max_curl_ang, min_curl_ang;

/* 02/12/07 */
void CGlView::DisplayHighlyCurl()
{
	double curl_threshold; 

	int i, j;
	Face *face;

	//glColor4f(1, 0, 1, 0.3);
	for(i = 0 ; i < Object.nfaces; i++)
	{
		face = Object.flist[i];

		if(face->length[2] < curl_angle)
			continue;
		
		curl_threshold = (face->length[0]-min_curl)/(max_curl-min_curl);

		if(curl_threshold < curl_percent)
			continue;



		glBegin(GL_TRIANGLES);
		//glBegin(GL_POLYGON);
		for(j = 0; j < face->nverts; j++)
		{
		float hsv[3] = {0, 1., 1.};
		float rgb[3] = {0.};
			hsv[0] = 256 * (max_curl_ang - Object.vlist[face->verts[j]]->length[2])/(max_curl_ang-min_curl_ang);
			HsvRgb(hsv, rgb);

			glColor4f(rgb[0], rgb[1], rgb[2], 0.3);
			glVertex2f(Object.vlist[face->verts[j]]->x, Object.vlist[face->verts[j]]->y);
		}
		glEnd();
	}
}


/*
Display divergence
*/
void CGlView::DisplayHighlyDiv()
{
	double curl_threshold; 

	int i, j;
	Face *face;

	//glColor4f(1, 0, 1, 0.3);
	for(i = 0 ; i < Object.nfaces; i++)
	{
		face = Object.flist[i];


		glBegin(GL_TRIANGLES);
		//glBegin(GL_POLYGON);
		for(j = 0; j < face->nverts; j++)
		{
		float hsv[3] = {0, 1., 1.};
		float rgb[3] = {0.};
			hsv[0] = 256 * (max_div - Object.vlist[face->verts[j]]->length[1])/(max_div-min_div);
			HsvRgb(hsv, rgb);

			glColor4f(rgb[0], rgb[1], rgb[2], 0.3);
			glVertex2f(Object.vlist[face->verts[j]]->x, Object.vlist[face->verts[j]]->y);
		}
		glEnd();
	}
}


/*
Display divergence
*/
void CGlView::DisplayMag()
{

	int i, j;
	Face *face;

	for(i = 0 ; i < Object.nfaces; i++)
	{
		face = Object.flist[i];


		glBegin(GL_TRIANGLES);
		for(j = 0; j < face->nverts; j++)
		{
		float hsv[3] = {0, 1., 1.};
		float rgb[3] = {0.};
		hsv[0] = 256 * (max_mag - Object.vlist[face->verts[j]]->mag_speed)/(max_mag-min_mag);
		HsvRgb(hsv, rgb);

			glColor4f(rgb[0], rgb[1], rgb[2], 0.3);
			glVertex2f(Object.vlist[face->verts[j]]->x, Object.vlist[face->verts[j]]->y);
		}
		glEnd();
	}
}



/*
Display the decomposition regions using asymmetry method
04/12/07
*/
void CGlView::DisplayDecomp_regions()
{

	int i, j;
	Face *face;
	Vertex *v;

	for(i = 0 ; i < Object.nfaces; i++)
	{
		face = Object.flist[i];


		glBegin(GL_TRIANGLES);
		for(j = 0; j < face->nverts; j++)
		{
			v = Object.vlist[face->verts[j]];
			float hsv[3] = {v->hue, v->saturation, 1.};
			float rgb[3] = {0.};
			HsvRgb(hsv, rgb);

			glColor4f(rgb[0], rgb[1], rgb[2], 0.3);
			glVertex2f(Object.vlist[face->verts[j]]->x, Object.vlist[face->verts[j]]->y);
		}
		glEnd();
	}
}


/*
Display the detected path between two nodes
*/
void CGlView::DisplayTestPath()
{
	int i, j;
	Face *face;
	float center1[2], center2[2];

	glColor3f(1, 1, 0);
	for(i = 0; i < num_tris_inpath; i++)
	{
		face = Object.flist[path[i]];
		center1[0] = center1[1] = 0;
		glBegin(GL_LINE_LOOP);
		for(j = 0; j < face->nverts; j++)
		{
			glVertex2f(Object.vlist[face->verts[j]]->x, Object.vlist[face->verts[j]]->y);
			center1[0] += Object.vlist[face->verts[j]]->x;
			center1[1] += Object.vlist[face->verts[j]]->y;
		}
		glEnd();

		if(i < num_tris_inpath-1)
		{
			face=Object.flist[path[i+1]];
			center2[0] = center2[1] = 0;
			for(j = 0; j < face->nverts; j++)
			{
				center2[0] += Object.vlist[face->verts[j]]->x;
				center2[1] += Object.vlist[face->verts[j]]->y;
			}
			glBegin(GL_LINES);
			glVertex2f(center1[0]/3, center1[1]/3);
			glVertex2f(center2[0]/3, center2[1]/3);
			glEnd();
		}
	}

	if(glob_t1 >= 0)
	{
		/*STAR*/
		glColor4f(0, 1, 0, 0.4);
		face = Object.flist[glob_t1];
		glBegin(GL_TRIANGLES);
		for(j = 0; j < face->nverts; j++)
		{
			glVertex2f(Object.vlist[face->verts[j]]->x, Object.vlist[face->verts[j]]->y);
		}
		glEnd();
	}

	if(glob_t2 >= 0)
	{
		/*END*/
		glColor4f(1, 0, 0, 0.4);
		face = Object.flist[glob_t2];
		glBegin(GL_TRIANGLES);
		for(j = 0; j < face->nverts; j++)
		{
			glVertex2f(Object.vlist[face->verts[j]]->x, Object.vlist[face->verts[j]]->y);
		}
		glEnd();
	}
}


/*
Display all the edges incident to the node of the selected triangle
*/

void CGlView::DisplayAllEdges(int node)
{
	int i, j, other_node;
	double center1[2], center2[2];
	int out_flag;

	Face *face = Object.flist[node];
	center1[0] = center1[1] = 0.;
	glColor3f(0, 1, 0);
	glBegin(GL_TRIANGLES);
	for(j = 0; j < face->nverts; j++)
	{
		glVertex2f(Object.vlist[face->verts[j]]->x, Object.vlist[face->verts[j]]->y);
		center1[0] += Object.vlist[face->verts[j]]->x;
		center1[1] += Object.vlist[face->verts[j]]->y;
	}
	glEnd();

	center1[0] /= 3;
	center1[1] /= 3;

	for (i = 0; i < sccnodes[node].nedges; i++)
	{
		other_node = sccedges[sccnodes[node].edges[i]].node_index2;
		out_flag = 1;
		if(other_node == node)
		{
			other_node = sccedges[sccnodes[node].edges[i]].node_index1;
			out_flag = 0;
		}

		/*Display the triangle it can reach*/
		face = Object.flist[other_node];
		glColor3f(1, 1, 0);
		glBegin(GL_LINE_LOOP);
		center2[0] = center2[1] = 0;
		for(j = 0; j < face->nverts; j++)
		{
			glVertex2f(Object.vlist[face->verts[j]]->x, Object.vlist[face->verts[j]]->y);
			center2[0] += Object.vlist[face->verts[j]]->x;
			center2[1] += Object.vlist[face->verts[j]]->y;
		}
		glEnd();

		center2[0] /= 3;
		center2[1] /= 3;

		if(out_flag == 1)
			glColor3f(0, 1, 0);
		else
			glColor3f(1, 0, 0);

		glBegin(GL_LINES);
			glVertex2f(center1[0], center1[1]);
			glVertex2f(center2[0], center2[1]);
		glEnd();
	}
}


void CGlView::DisplayEigenVecForSCC()
{
	int i, j;
	Face *face;
	Vertex *vert;
	Edge *cur_e;
	
	glColor3f(0.6, 0.6, 0.6);
	//glBegin(GL_LINES);
	//for(i = 0; i < scclist.num_sccs; i++)
	//{
	//	if(scclist.scccomponents[i].num_nodes < 2)
	//		continue;

	//	for(j = 0; j < scclist.scccomponents[i].num_nodes; j++)
	//	{
	//		face = Object.flist[scclist.scccomponents[i].nodes[j]];

	//		for(int k = 0; k < 3; k++)
	//		{
	//			cur_e = face->edges[k];

	//			glVertex2f(Object.vlist[cur_e->verts[0]]->x, Object.vlist[cur_e->verts[0]]->y);
	//			
	//			glVertex2f(Object.vlist[cur_e->verts[1]]->x, Object.vlist[cur_e->verts[1]]->y);
	//		}

	//	}
	//}
	//glEnd();

	//for(i = 0; i < scclist.num_sccs; i++)
	//{
	//	if(scclist.scccomponents[i].num_nodes < 2)
	//		continue;

	//	for(j = 0; j < scclist.scccomponents[i].num_nodes; j++)
	//	{
	//		face = Object.flist[scclist.scccomponents[i].nodes[j]];

	//		glColor3f(1, 0, 0); //major eigen vector
	//		glPushMatrix();
	//		glTranslatef(face->center[0], face->center[1], 0);
	//		glRotatef(atan2(face->eigen[0].entry[1], face->eigen[0].entry[0])*180./M_PI, 0, 0, 1);
	//		glScalef(ARROWSCALE/4., ARROWSCALE/4., 1);
	//			glBegin(GL_LINES);
	//			glVertex2f(0, 0);
	//			glVertex2f(1, 0);
	//			glEnd();

	//			////Draw the wings of the arrow
	//			glBegin(GL_LINES);
	//			glVertex2f(1, 0);
	//			glVertex2f(0.8, 0.16);

	//			glVertex2f(1, 0);
	//			glVertex2f(0.8, -0.16);
	//			glEnd();
	//		glPopMatrix();
	//		
	//		glColor3f(0, 1, 0); //minor eigen vector
	//		glPushMatrix();
	//		glTranslatef(face->center[0], face->center[1], 0);
	//		glRotatef(atan2(face->eigen[1].entry[1], face->eigen[1].entry[0])*180./M_PI, 0, 0, 1);
	//		glScalef(ARROWSCALE/4., ARROWSCALE/4., 1);
	//			glBegin(GL_LINES);
	//			glVertex2f(0, 0);
	//			glVertex2f(1, 0);
	//			glEnd();

	//			////Draw the wings of the arrow
	//			glBegin(GL_LINES);
	//			glVertex2f(1, 0);
	//			glVertex2f(0.8, 0.16);

	//			glVertex2f(1, 0);
	//			glVertex2f(0.8, -0.16);
	//			glEnd();
	//		glPopMatrix();
	//	}
	//}


	//glPointSize(3.);
	//glBegin(GL_POINTS);
	//for(i = 0; i < scclist.num_sccs; i++)
	//{
	//	if(scclist.scccomponents[i].num_nodes < 2)
	//		continue;

	//	for(j = 0; j < scclist.scccomponents[i].num_nodes; j++)
	//	{
	//		face = Object.flist[scclist.scccomponents[i].nodes[j]];

	//		for(int k = 0; k < 3; k++)
	//		{
	//			cur_e = face->edges[k];
	//			glColor3f(0, 0, 0); //sep points

	//			if(cur_e->find_sep == 1)
	//			glVertex2f(cur_e->sep.entry[0], cur_e->sep.entry[1]);
	//			
	//			glColor3f(0.8, 0, 1); //att points
	//			if(cur_e->find_attp == 1)
	//			glVertex2f(cur_e->attp.entry[0], cur_e->attp.entry[1]);
	//		}

	//	}
	//}
	//glEnd();
	//glPointSize(1.);


	//glBegin(GL_LINES);
	//	for(j = 0; j < Object.nfaces; j++)
	//	{
	//		face = Object.flist[j];

	//		for(int k = 0; k < 3; k++)
	//		{
	//			cur_e = face->edges[k];

	//			glVertex2f(Object.vlist[cur_e->verts[0]]->x, Object.vlist[cur_e->verts[0]]->y);
	//			
	//			glVertex2f(Object.vlist[cur_e->verts[1]]->x, Object.vlist[cur_e->verts[1]]->y);
	//		}

	//	}
	//glEnd();

	//for(j = 0; j < Object.nfaces; j++)
	//{
	//	face = Object.flist[j];

	//	glColor3f(1, 0, 0); //major eigen vector
	//	glPushMatrix();
	//	glTranslatef(face->center[0], face->center[1], 0);
	//	glRotatef(atan2(face->eigen[0].entry[1], face->eigen[0].entry[0])*180./M_PI, 0, 0, 1);
	//	glScalef(ARROWSCALE/4., ARROWSCALE/4., 1);
	//		glBegin(GL_LINES);
	//		glVertex2f(0, 0);
	//		glVertex2f(1, 0);
	//		glEnd();

	//		////Draw the wings of the arrow
	//		glBegin(GL_LINES);
	//		glVertex2f(1, 0);
	//		glVertex2f(0.8, 0.16);

	//		glVertex2f(1, 0);
	//		glVertex2f(0.8, -0.16);
	//		glEnd();
	//	glPopMatrix();
	//	
	//	glColor3f(0, 1, 0); //minor eigen vector
	//	glPushMatrix();
	//	glTranslatef(face->center[0], face->center[1], 0);
	//	glRotatef(atan2(face->eigen[1].entry[1], face->eigen[1].entry[0])*180./M_PI, 0, 0, 1);
	//	glScalef(ARROWSCALE/4., ARROWSCALE/4., 1);
	//		glBegin(GL_LINES);
	//		glVertex2f(0, 0);
	//		glVertex2f(1, 0);
	//		glEnd();

	//		////Draw the wings of the arrow
	//		glBegin(GL_LINES);
	//		glVertex2f(1, 0);
	//		glVertex2f(0.8, 0.16);

	//		glVertex2f(1, 0);
	//		glVertex2f(0.8, -0.16);
	//		glEnd();
	//	glPopMatrix();
	//}

	//for(j = 0; j < Object.nverts; j++)
	//{
	//	vert = Object.vlist[j];

	//	glColor3f(1, 0, 0); //major eigen vector
	//	glPushMatrix();
	//	glTranslatef(vert->x, vert->y, 0);
	//	glRotatef(atan2(vert->evec[0].entry[1], vert->evec[0].entry[0])*180./M_PI, 0, 0, 1);
	//	glScalef(ARROWSCALE/6., ARROWSCALE/6., 1);
	//		glBegin(GL_LINES);
	//		glVertex2f(0, 0);
	//		glVertex2f(1, 0);
	//		glEnd();

	//		glBegin(GL_LINES);
	//		glVertex2f(0, 0);
	//		glVertex2f(-1, 0);
	//		glEnd();

	//		////Draw the wings of the arrow
	//		//glBegin(GL_LINES);
	//		//glVertex2f(1, 0);
	//		//glVertex2f(0.8, 0.16);

	//		//glVertex2f(1, 0);
	//		//glVertex2f(0.8, -0.16);
	//		//glEnd();
	//	glPopMatrix();
	//	
	//	glColor3f(0, 1, 0); //minor eigen vector
	//	glPushMatrix();
	//	glTranslatef(vert->x, vert->y, 0);
	//	glRotatef(atan2(vert->evec[1].entry[1], vert->evec[1].entry[0])*180./M_PI, 0, 0, 1);
	//	glScalef(ARROWSCALE/6., ARROWSCALE/6., 1);
	//		glBegin(GL_LINES);
	//		glVertex2f(0, 0);
	//		glVertex2f(1, 0);
	//		glEnd();

	//		glBegin(GL_LINES);
	//		glVertex2f(0, 0);
	//		glVertex2f(-1, 0);
	//		glEnd();

	//		////Draw the wings of the arrow
	//		//glBegin(GL_LINES);
	//		//glVertex2f(1, 0);
	//		//glVertex2f(0.8, 0.16);

	//		//glVertex2f(1, 0);
	//		//glVertex2f(0.8, -0.16);
	//		//glEnd();
	//	glPopMatrix();
	//}


	glPointSize(3.);
	glBegin(GL_POINTS);

	for(j = 0; j < Object.nfaces; j++)
	{
		face = Object.flist[j];

		for(int k = 0; k < 3; k++)
		{
			cur_e = face->edges[k];
			//glColor3f(0, 0, 1);   //sep points

			glColor3f(1, 0, 1); //att points
			if(cur_e->find_sep == 1)
			{
				//if(cur_e->visited == 1)
				//	glColor3f(1, 1, 1);
				glVertex2f(cur_e->sep.entry[0], cur_e->sep.entry[1]);
			}
			
			//glColor3f(0.8, 1, 0); //att points
			glColor3f(0, 1, 1);   //sep points
			if(cur_e->find_attp == 1)
			{
				//if(cur_e->visited == 1)
				//	glColor3f(1, 1, 1);
				glVertex2f(cur_e->attp.entry[0], cur_e->attp.entry[1]);
			}
		}

	}

	glEnd();
	glPointSize(1.);
}


void CGlView::DisplaySamplingPts()
{
	int i, j;

	glColor3f(0, 1, 0);
	glPointSize(3.);
	glBegin(GL_POINTS);

	for(i = 0; i < cur_traj_index; i++)
	{
		for(j = 0; j < sampleptslist[i].num_samples; j++)
		{
			glVertex2f(sampleptslist[i].samples[j].gpt[0], sampleptslist[i].samples[j].gpt[1]);
		}
	}

	glEnd();
	glPointSize(1.);
	glColor3f(0, 0, 0);
}




void CGlView::display_trace_SampltPts(int back__forward)
{
	int i;
	Vertex *v;

	glPointSize(3.);
	glBegin(GL_POINTS);

	/*display vertices*/
	for(i=0; i<Object.nverts; i++)
	{
		v=Object.vlist[i];
		if(back__forward == 0)           /*_forward color*/
			glColor3fv(v->f_color);
		else if(back__forward == 1)      /*backward color*/
			glColor3fv(v->b_color);

		glVertex2f(v->x, v->y);
	}

	/*display sample points on edges and inside triangles*/
	TraceSamplePtList *temp;

	if(back__forward < 2)
	{
		
		if(back__forward ==0) temp = _forward;
		else if(back__forward == 1) temp=backward_spts;
		for(i=0; i<temp->nelems; i++)
		{
			glColor3fv(temp->samplepts[i]->color);
			glVertex2fv(temp->samplepts[i]->pos);
		}
	}
	else
	{
		temp = _forward;
		for(i=0; i<temp->nelems; i++)
		{
			glColor3fv(temp->samplepts[i]->color);
			glVertex2fv(temp->samplepts[i]->pos);
		}
		
		temp = backward_spts;
		for(i=0; i<temp->nelems; i++)
		{
			glColor3fv(temp->samplepts[i]->color);
			glVertex2fv(temp->samplepts[i]->pos);
		}
	}

	glEnd();
}


	
void CGlView::display_sample_density()
{
	int i, j;
	Vertex *v;
	Face *face;
	
	for(i=0; i<Object.nfaces; i++)
	{
		face = Object.flist[i];

		glBegin(GL_TRIANGLES);
		for(j=0; j<3; j++)
		{
			v = Object.vlist[face->verts[j]];

			glColor4f(v->d_color[0], v->d_color[1], v->d_color[2], 0.4);

			glVertex2f(v->x, v->y);
		}
		glEnd();
	}
}


void getcolor_mcgconnector_frac(int num, int frac, float rgb[3])
{
	float hsv[3] = {0, 0.5, 0.5};

	hsv[0] = 360 * (double)num/frac;

	HsvRgb(hsv,rgb);
}
			
	
extern int cur_mcgedge_index;

void CGlView::display_MCG_connections()
{
	int i, j, k;

	int node1, node2;
	Vertex *v;

	for(i=0; i<cur_mcgedge_index; i++)
	{
		if(mcgedges[i].cancel)
			continue;

		node1 = mcgedges[i].node_index1;
		node2 = mcgedges[i].node_index2;

		for(j=0; j<mcgedges[i].ntris; j++)
		{
			/**/
			//float rgb[3] = {0.};
			//glColor4f(0.3, 0.3, 0.3, 0.8);

			//if(mcgnodes[node1].type == 0)
			//	glColor4f(0, 0.35, 0, 0.8);
			//else if(mcgnodes[node2].type == 1)
			//	glColor4f(0.35, 0, 0, 0.8);
			//else
			//	glColor4f(0.0, 0.0, 0.4, 0.8);
			
			float rgb[3] = {0.};
			glColor4f(1, 1, 1, 0.3);

			Face *face = Object.flist[mcgedges[i].triangles[j]];

			glBegin(GL_TRIANGLES);
			for(k=0; k<3; k++)
			{
				glVertex2f(Object.vlist[face->verts[k]]->x,
					Object.vlist[face->verts[k]]->y);
			}
			glEnd();

			double center[2] = {0.};
			center[0] = (Object.vlist[face->verts[0]]->x+Object.vlist[face->verts[1]]->x
				+Object.vlist[face->verts[2]]->x)/3.;
			center[1] = (Object.vlist[face->verts[0]]->y+Object.vlist[face->verts[1]]->y
				+Object.vlist[face->verts[2]]->y)/3.;
			
			if(mcgnodes[node1].type == 0)
				glColor4f(0, 0.6, 0, 0.8);
			else if(mcgnodes[node2].type == 1)
				glColor4f(0.9, 0, 0, 0.8);
			else
				glColor4f(0.0, 0.0, 0.9, 0.8);

			DrawSolidCircle_size(center[0], center[1], 0.007);
		}
	}
}



//#include "taumap.h"
extern EdgeSampleImgList *t_AnedgeSampleImgList;
	
void CGlView::display_image_edge()
{
	int i;

	glBegin(GL_LINES);
	for(i=0; i<t_AnedgeSampleImgList->nsampleimgs; i++)
	{
		glVertex2f(t_AnedgeSampleImgList->sampleimgs[i]->l_pos[0],
			t_AnedgeSampleImgList->sampleimgs[i]->l_pos[1]);
	}
	glEnd();
}



void CGlView::display_ver_time()
{
	int i, j;
	Face *face;
	Vertex *v;
	
	for(i=0; i<Object.nfaces; i++)
	{
		face = Object.flist[i];

		glBegin(GL_TRIANGLES);
		for(j=0; j<3; j++)
		{
			v = Object.vlist[face->verts[j]];
			glColor4f(v->b_color[0], v->b_color[1], v->b_color[2], 0.4);
			glVertex2f(v->x, v->y);
		}
		glEnd();
	}
}


extern DegeneratePt *degpts ;
extern int ndegpts ;

/**------ Display the detected degenerate triangles ------**/
void CGlView::display_degenerate_tris()
{
	int i, j;
	Face *face;
	Vertex *v;

	for(i=0; i<ndegenerate_tris; i++)
	{
		face = Object.flist[degenerate_tris[i]];

		if(face->degenerate_type == 0)
			glColor3f(1, 0, 0);
		else if(face->degenerate_type == 1)
			glColor3f(0, 1, 0);
		else if(face->degenerate_type == 2)
			glColor3f(1, 1, 0);
		else if(face->degenerate_type == 3)
			glColor3f(0, 0, 1);
		else
			glColor3f(1, 1, 1);

		glBegin(GL_LINE_LOOP);
		for(j=0; j<3; j++)
		{
			v = Object.vlist[face->verts[j]];
			glVertex2f(v->x, v->y);
		}
		glEnd();
	}
}


/**------ Display the detected degenerate quads/cells ------**/
void CGlView::display_degenerate_cells()
{
	int i, j;
	QuadCell *face;
	QuadVertex *v;

	for(i=0; i<ndegenerate_tris; i++)
	{
		face = quadmesh->quadcells[degenerate_tris[i]];

		if(face->degenerate_type == 0)
			glColor3f(1, 0, 0);
		else if(face->degenerate_type == 1)
			glColor3f(0, 1, 0);
		else if(face->degenerate_type == 2)
			glColor3f(1, 1, 0);
		else if(face->degenerate_type == 3)
			glColor3f(0, 0, 1);
		else
			glColor3f(1, 1, 1);

		glBegin(GL_LINE_LOOP);
		for(j=0; j<face->nverts; j++)
		{
			v = quadmesh->quad_verts[face->verts[j]];
			glVertex2f(v->x, v->y);
		}
		glEnd();
	}
}

	
bool CGlView::is_ten_designElem(double x, double y, int &degptindex)
{
	int i;
	for( i = 0; i < ntenelems; i++)
	{
		if(ten_designelems[i].ID >= 0 && !ten_designelems[i].deleted){
			if((fabs((ten_designelems[i].centerx - x)) < 0.01)
				&&(fabs((ten_designelems[i].centery - y)) < 0.01))
			{
				degptindex = ten_designelems[i].ID;
				return true;
			}
		}
	}
	if( i >= ntenelems) return false;
}

/*visualize the extracted degenerate points 09/19/2007*/



void CGlView::display_degenerate_pts(GLenum mode)
{
	int i, singular_id = 0;

	for(i = 0; i < ndegpts; i++)
	{
		if(degpts[i].type == 0) /*wedge*/
			glColor3f(1, 0, 0);
		else if(degpts[i].type == 1) /*trisector*/
			glColor3f(0, 1, 0);
		else if(degpts[i].type == 2) /*node*/
			glColor3f(1, 1, 0);
		else if(degpts[i].type == 3) /*center*/
			glColor3f(1, 0, 1);
		else if(degpts[i].type == 4) /*saddle*/
			glColor3f(0, 0, 1);
		else
			glColor3f(1, 1, 1);

		if(is_ten_designElem(degpts[i].gcx, degpts[i].gcy, singular_id))
		{
			if(PairCancelOn == 1 /*|| LimitSingCancelOn == 1*/){
				if(mode == GL_SELECT)
					glLoadName(NAMEOFSINGULARITY + i);////Assign name for design element
					//glLoadName(NAMEOFSINGULARITY + singular_id);////Assign name for design element
			}
			else{
				//this is for design element editing
				if(mode == GL_SELECT)
					glLoadName(singular_id+1);  ////using the ID of the element as its name (just the same as NAMEOFSINGELEM)
			}
		}

		else{
			if((PairCancelOn == 1 || LimitSingCancelOn == 1) && mode == GL_SELECT)
				glLoadName(NAMEOFSINGULARITY + i);
		}


		DrawSolidCircle_size(degpts[i].gcx, degpts[i].gcy, 0.006/zoom_factor);

		glColor3f(0, 0, 0);
		glLineWidth(1.4);
		draw_hollow_circle_size(degpts[i].gcx, degpts[i].gcy, 0.0065/zoom_factor);
	}
}

void CGlView::display_seps_directions()
{
	int i, j;

	for(i=0; i<ndegpts; i++)
	{
		int nseps = degpts[i].nseps;
		glLineWidth(1.5);
		for(j=0; j<nseps; j++)
		{
			if(j==0)
				glColor3f(1, 0.7, 0);
			else if(j==1)
				glColor3f(0, 0.8, 0.8);
			else
				glColor3f(0.8, 0, 0.8);

			glBegin(GL_LINES);
			glVertex2f(degpts[i].gcx, degpts[i].gcy);
			glVertex2f(degpts[i].gcx + 0.03*degpts[i].s[j].entry[0],
				degpts[i].gcy + 0.03*degpts[i].s[j].entry[1]);
			glEnd();
		}
	}
}


void CGlView::display_tenElem_EditBox(GLenum mode)
{
	int i;

	glLineWidth(4);

	for(i = 0; i < ntenelems; i++)
	{
		if(ten_designelems[i].ID >= 0 && !ten_designelems[i].deleted)
		{
			if(mode == GL_SELECT)
			{
				glLoadName(i+1);  ////using the ID of the element as its name (just the same as NAMEOFSINGELEM)
			}
			if(ten_designelems[i].type == 0) /*wedge*/
				glColor3f(1, 0, 0);
			else if(ten_designelems[i].type == 1) /*trisector*/
				glColor3f(0, 1, 0);
			else if(ten_designelems[i].type == 2) /*node*/
				glColor3f(1, 1, 0);
			else if(ten_designelems[i].type == 3) /*center*/
				glColor3f(1, 0, 1);
			else if(ten_designelems[i].type == 4) /*saddle*/
				glColor3f(0, 0, 1);

			DrawEditBox(ten_designelems[i].cur_editbox.p1,\
				ten_designelems[i].cur_editbox.p2,\
				ten_designelems[i].cur_editbox.p3,\
				ten_designelems[i].cur_editbox.p4);

		}
	}
	glLineWidth(1.);
}


void CGlView::ten_hitProcess(double ss, double st)
{
	GLuint  selectBuffer[SELECTBUFFERSIZE];
	GLint	vp[4] = {0, 0 ,1, 1};
	int hits;

	choose_ID = -1;
	selectedIntersect = -1;

	////Build the selection buffer here
	glSelectBuffer(SELECTBUFFERSIZE, selectBuffer);

   	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);
	

	if(sharedvars.RegionSmoothOn || sharedvars.SelStreetRegToEditOn)
		gluPickMatrix(ss, st, 1e-8, 1e-8, vp );  ////set smaller pick window for triangle selection

	else
		gluPickMatrix(ss, st, 0.005, 0.005, vp );  ////set a larger pick window for element selection
		
	glOrtho(0, 1,  0, 1,  0, 50);


	////If one of the element has been selected for being edited
	TransformType = 0;
	which_control_point = -1;
	which_shapectrpt = -1;


	/*using quad mesh 09/26/2007*/

	/*we apply the transformation here 11/08/2007*/
	transform_fun();


	//else if(ShapeEditOn == 1)
	//	DisplayShapeControlPts(GL_SELECT);
	
	//else if(PairCancelOn == 1 )  //for singularity pair cancellation only
	//{
	//	display_degenerate_pts(GL_SELECT);
	//}

	if(deleteElemOn ||sharedvars.EditElementOn || sharedvars.MoveElemOn)
	{
		display_degenerate_pts(GL_SELECT); ////we need to find another way to allocate name for element and its control points
		display_tenElem_EditBox(GL_SELECT);
		display_tenRegElem(GL_SELECT);
		if(chosen_tenelem_ID>=0 && chosen_tenelem_ID < NAMEOFREGELEM)
			display_singularElem_controlPts_tensor(GL_SELECT, chosen_tenelem_ID);
	}

	else if(sharedvars.RegionSmoothOn || sharedvars.SelStreetRegToEditOn)
		display_quadmesh_select(GL_SELECT);


	else{
		
		display_degenerate_pts(GL_SELECT); ////we need to find another way to allocate name for element and its control points
		display_tenElem_EditBox(GL_SELECT);

		if(chosen_tenelem_ID>=0 && chosen_tenelem_ID < NAMEOFREGELEM)
			display_singularElem_controlPts_tensor(GL_SELECT, chosen_tenelem_ID);
		
		display_tenRegElem(GL_SELECT);

		//DisplayRegularControlPoints(GL_SELECT, choose_ID-NAMEOFREGELEM);

		if(/*showTensorOn==1 && */streetNetEditOn)
			display_streetnet(GL_SELECT);
	}

	hits = glRenderMode(GL_RENDER);


	if(hits>0)
	{
		////Because this is 2D plane, we need not worry about the overlap of objects
		////It will have at most one object being selected at one click
		//int objname = selectBuffer[3];

		////New hit processing here 06/29/05
		if(hits == 1)
		{
			if(selectBuffer[3] > 0 && selectBuffer[3] < NAMEOFSINGCONTROL)
			{
				////select an element
				chosen_tenelem_ID = selectBuffer[3]-1;					
				
				if(deleteElemOn)
				{
					if(chosen_tenelem_ID < NAMEOFREGELEM)
					{
						/*set the corresponding elements as deleted*/
						ten_designelems[chosen_tenelem_ID].deleted=true;
						cur_chosen_region=ten_designelems[chosen_tenelem_ID].which_region;
					}

					else 
					{
						/*set the corresponding elements as deleted*/
						ten_regularelems[selectBuffer[3]-NAMEOFREGELEM-1].deleted=true;
						cur_chosen_region=
							ten_regularelems[selectBuffer[3]-NAMEOFREGELEM-1].which_region;
					}
				}

				if(chosen_tenelem_ID < NAMEOFREGELEM)
				{
					/*set the corresponding elements as deleted*/
					cur_chosen_region=ten_designelems[chosen_tenelem_ID].which_region;
				}

				else 
				{
					/*set the corresponding elements as deleted*/
					cur_chosen_region=
						ten_regularelems[selectBuffer[3]-NAMEOFREGELEM-1].which_region;
				}
			}

			else if(selectBuffer[3] >= NAMEOFSINGCONTROL && selectBuffer[3] < NAMEOFSINGULARITY)
			{    ////select a control point
				if(selectBuffer[3] < NAMEOFREGCONTROL) ////control point for singular element
				{
					ControlPointForSinElem(selectBuffer[3]-NAMEOFSINGCONTROL);
				}
				else{                                  ////control point for regular element
                    ControlPointForRegElem(selectBuffer[3]-NAMEOFREGCONTROL);
				}
			}
			
			else if(TrajectoryOn == 1 && selectBuffer[3]>=NAMEOFTRIANGLE) /*select the tracing triangle*/
			{
				SelectTriangleID = selectBuffer[3]-NAMEOFTRIANGLE;
			}
			

			else if(streetNetEditOn && selectBuffer[3]>=NAMEOFINTERSECTS)
			{
				selectedIntersect = selectBuffer[3]-NAMEOFINTERSECTS;
			}

			else if(sharedvars.RegionSmoothOn||sharedvars.SelStreetRegToEditOn)  ////under select smoothing region mode
			{
				SelectTriangleID = selectBuffer[3] - NAMEOFTRIANGLE;
				/*  We need to change the s and t value here 12/11/2007*/
				ss=(ss-.5-trans_x)/zoom_factor+.5;
				st=(st-.5-trans_y)/zoom_factor+.5;
				AddToPointList(ss, st, SelectTriangleID);
			}
		}

		else{ ////hits>1, normally, hits == 2
			////This situation should only happen under editing mode
			int large_objname, small_objname;
			if(selectBuffer[3] > selectBuffer[7])
			{
				large_objname = selectBuffer[3];
				small_objname = selectBuffer[7];
			}
			else{
				large_objname = selectBuffer[7];
				small_objname = selectBuffer[3];
			}

			///////////////////////////////////////////////
			if(EditModeOn == 1) ////Priory to control point
			{
				if(large_objname < NAMEOFREGCONTROL) ////control point for singular element
				{
					ControlPointForSinElem(large_objname-NAMEOFSINGCONTROL);
				}
				else{                                  ////control point for regular element
                    ControlPointForRegElem(large_objname-NAMEOFREGCONTROL);
				}

				if(chosen_tenelem_ID < NAMEOFREGELEM)
				{
					/*set the corresponding elements as deleted*/
					cur_chosen_region=ten_designelems[chosen_tenelem_ID].which_region;
				}

				else 
				{
					/*set the corresponding elements as deleted*/
					cur_chosen_region=
						ten_regularelems[selectBuffer[3]-NAMEOFREGELEM-1].which_region;
				}
			}

			else if(streetNetEditOn && large_objname>=NAMEOFINTERSECTS)
			{
				selectedIntersect = large_objname-NAMEOFINTERSECTS;
			}
			
			else if(deleteElemOn)
			{
				if(small_objname < NAMEOFREGELEM)
				{
					chosen_tenelem_ID = small_objname-1;	
					/*set the corresponding elements as deleted*/
					degpts[chosen_tenelem_ID].deleted=true;
				}

				else if(small_objname>=NAMEOFREGELEM&&small_objname<NAMEOFSINGCONTROL)
				{
					chosen_tenelem_ID = small_objname-NAMEOFREGELEM;	
					/*set the corresponding elements as deleted*/
					ten_regularelems[chosen_tenelem_ID].deleted=true;
				}
				if(chosen_tenelem_ID < NAMEOFREGELEM)
				{
					/*set the corresponding elements as deleted*/
					cur_chosen_region=ten_designelems[chosen_tenelem_ID].which_region;
				}

				else 
				{
					/*set the corresponding elements as deleted*/
					cur_chosen_region=
						ten_regularelems[selectBuffer[3]-NAMEOFREGELEM-1].which_region;
				}
			}
			
			else if(sharedvars.RegionSmoothOn||sharedvars.SelStreetRegToEditOn)  ////under select smoothing region mode
			{
				SelectTriangleID = small_objname - NAMEOFTRIANGLE;
				/*  We need to change the s and t value here 12/11/2007*/
				ss=(ss-.5-trans_x)/zoom_factor+.5;
				st=(st-.5-trans_y)/zoom_factor+.5;
				AddToPointList(ss, st, SelectTriangleID);
			}

		}
	}

	else{
		chosen_tenelem_ID = -1;
		SelectTriangleID = -1;
	}


	//glMatrixMode( GL_MODELVIEW );
	//glPopMatrix();


LL:	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	//glMatrixMode( GL_MODELVIEW );
	//glLoadIdentity();

}





////display the regular elements using arrows
void CGlView::display_tenRegElem(GLenum mode)
{
	glLineWidth(2.);
	for(int i = 0; i < nten_regelems; i++)
	{
	   if(ten_regularelems[i].ID>0)
	   {
		   ////Display the arrow
		   if(mode == GL_SELECT)
			   glLoadName(ten_regularelems[i].ID);

		   if(ten_regularelems[i].deleted)
			   continue;

		   	////perform the user defined transformation for editing

		   if(ten_regularelems[i].type == 0) ////basic regular element
		       glColor3f(0, 1, 1);
		   else if(ten_regularelems[i].type == 1) ////convergent element
		       glColor3f(1, 0.5, 0);
		       //glColor3f(1, 1, 0);       //for the paper
		   else                         ////divergent element
		       glColor3f(0, 1, 0.5);

		   if(chosen_tenelem_ID >= NAMEOFREGELEM && chosen_tenelem_ID < NAMEOFSINGCONTROL
			   && chosen_tenelem_ID-NAMEOFREGELEM == i)
			   glColor3f(1, 1, 1);


		   /*we draw two line segments*/
		   glBegin(GL_LINES);
		   glVertex2f(ten_regularelems[i].base[0], ten_regularelems[i].base[1]);
		   glVertex2f(ten_regularelems[i].end[0], ten_regularelems[i].end[1]);
		   glEnd();

		   glBegin(GL_LINES);
		   glVertex2f(ten_regularelems[i].base[0], ten_regularelems[i].base[1]);
		   glVertex2f(ten_regularelems[i].base[0]-ten_regularelems[i].Direct.entry[0], 
			   ten_regularelems[i].base[1]-ten_regularelems[i].Direct.entry[1]);
		   glEnd();

		   /*draw three control points*/

		   //if(mode == GL_SELECT)
			  // glLoadName(NAMEOFREGCONTROL+1);   ////control point at the base
		   DrawSolidRect_size(ten_regularelems[i].base[0], ten_regularelems[i].base[1], 0.006/zoom_factor);

		   //if(mode == GL_SELECT)
			  // glLoadName(NAMEOFREGCONTROL+2);   ////control point at the base
		   DrawSolidCircle_size(ten_regularelems[i].base[0]-ten_regularelems[i].Direct.entry[0], 
			   ten_regularelems[i].base[1]-ten_regularelems[i].Direct.entry[1], 0.006/zoom_factor);
		   
		   //if(mode == GL_SELECT)
			  // glLoadName(NAMEOFREGCONTROL+3);   ////control point at the end
		   DrawSolidCircle_size(ten_regularelems[i].end[0], ten_regularelems[i].end[1], 0.006/zoom_factor);

	   }
	}
}



/*display the major tensor line placement*/
void CGlView::display_major_tenlines(GLenum mode)
{
	if(major == NULL)
		return;
	int i, j, k;
	Trajectory *cur_traj;
	glColor3f(1, 0, 0);
	glLineWidth(1.5);
	//glClear(GL_ACCUM_BUFFER_BIT);
	//for(int i = 0; i < 16; i++)
	//{
	//	glPushMatrix ();
	//	glTranslatef (ji16[i].x*2.0*quadmesh->radius/512, ji16[i].y*2.0*quadmesh->radius/512, 0.0);
		
		for(j=0; j<major->evenstreamlines->ntrajs; j++)
		{
			cur_traj = major->evenstreamlines->trajs[j];

			for(k=0; k<cur_traj->nlinesegs; k++)
			{
				glBegin(GL_LINES);
				glVertex2f(cur_traj->linesegs[k].gstart[0], cur_traj->linesegs[k].gstart[1]);
				glVertex2f(cur_traj->linesegs[k].gend[0], cur_traj->linesegs[k].gend[1]);
				glEnd();
			}
		}
		
	//	/*draw the tensor line*/
	//	glPopMatrix ();
	//	glAccum(GL_ACCUM, 1.0/16);
	//}
	//glAccum (GL_RETURN, 1.0);
	//glReadBuffer(GL_BACK);
	//glLineWidth(1.);
}

/*display the major tensor line placement*/
void CGlView::display_minor_tenlines(GLenum mode)
{
	if(minor == NULL)
		return;
	int i, j, k;
	Trajectory *cur_traj;
	glLineWidth(1.5);
	glColor3f(0, 1, 0);
	//glClear(GL_ACCUM_BUFFER_BIT);
	//for(int i = 0; i < 16; i++)
	//{
	//	glPushMatrix ();
	//	glTranslatef (ji16[i].x*2.0*quadmesh->radius/512, ji16[i].y*2.0*quadmesh->radius/512, 0.0);
		
		/*draw the tensor line*/
		for(j=0; j<minor->evenstreamlines->ntrajs; j++)
		{
			cur_traj = minor->evenstreamlines->trajs[j];

			for(k=0; k<cur_traj->nlinesegs; k++)
			{
				glBegin(GL_LINES);
				glVertex2f(cur_traj->linesegs[k].gstart[0], cur_traj->linesegs[k].gstart[1]);
				glVertex2f(cur_traj->linesegs[k].gend[0], cur_traj->linesegs[k].gend[1]);
				glEnd();
			}
		}

	//	glPopMatrix ();
	//	glAccum(GL_ACCUM, 1.0/16);
	//}
	//glAccum (GL_RETURN, 1.0);
	//glReadBuffer(GL_BACK);
	//glLineWidth(1.);
}


/*    display the major roads   */
void CGlView::display_level1_tenlines(GLenum mode)
{
	if(major_level1 == NULL||minor_level1 == NULL)
		return;

	int /*i,*/ j, k;
	Trajectory *cur_traj;
	glColor3f(.95, .6, 0);
	glLineWidth(2.5);
		
	for(j=0; j<major_level1->evenstreamlines->ntrajs; j++)
	{
		cur_traj = major_level1->evenstreamlines->trajs[j];

		if(mode==GL_SELECT)
			glLoadName(j+NAMEOFMAJROADS);

		for(k=0; k<cur_traj->nlinesegs; k++)
		{
			glBegin(GL_LINES);
			glVertex2f(cur_traj->linesegs[k].gstart[0], cur_traj->linesegs[k].gstart[1]);
			glVertex2f(cur_traj->linesegs[k].gend[0], cur_traj->linesegs[k].gend[1]);
			glEnd();
		}
	}
	
	glColor3f(.65, .9, 0);
	for(j=0; j<minor_level1->evenstreamlines->ntrajs; j++)
	{
		cur_traj = minor_level1->evenstreamlines->trajs[j];
		
		if(mode==GL_SELECT)
			glLoadName(j+NAMEOFMAJROADS+major_level1->evenstreamlines->ntrajs);

		for(k=0; k<cur_traj->nlinesegs; k++)
		{
			glBegin(GL_LINES);
			glVertex2f(cur_traj->linesegs[k].gstart[0], cur_traj->linesegs[k].gstart[1]);
			glVertex2f(cur_traj->linesegs[k].gend[0], cur_traj->linesegs[k].gend[1]);
			glEnd();
		}
	}
}


/*    display the obtained seeds  
      NOTE: this is very important to let user know where are the seeds
*/
void CGlView::display_init_seeds(GLenum mode)
{

	int i;
	float hsv[3]={0.,1.,1.};
	float rgb[3]={0.};

	/*   display the boundary seeds    */
	if(seedsalongbounds !=NULL)
	{
		double rang=seedsalongbounds->max_weight-seedsalongbounds->min_weight;

		if(rang==0) return;

		//glColor3f(1, 0, 0);
		glPointSize(5.);
		for(i=0;i<seedsalongbounds->nseeds;i++)
		{
			if(mode == GL_SELECT)
				glLoadName(NAMEOFSEEDS+1);  ////

			//hsv[0]=250.*((double)i/(seedsalongbounds->nseeds-1));
			hsv[0]=250-250*(seedsalongbounds->seeds[i]->weight-seedsalongbounds->min_weight)
				/rang;

			if(fabs(hsv[0])<1.e-8) hsv[0]=0.;

			HsvRgb(hsv, rgb);
			glColor3fv(rgb);

			glBegin(GL_POINTS);
				glVertex2f(seedsalongbounds->seeds[i]->pos[0],
					seedsalongbounds->seeds[i]->pos[1]);
			glEnd();
		}
		glPointSize(1.);
	}
}



/*display the major tensor line placement*/
void CGlView::display_major_samps()
{
	if(major == NULL)
		return;
	int i, j, k;
	Trajectory *cur_traj;
	glLineWidth(1.5);
	glColor3f(1, 1, 0);

	SamplePtList *cur_samplist;

	glColor3f(0, 1, 0);
	glPointSize(5.);
	/*draw the tensor line*/

	for(j=0; j<major->evenstreamlines->ntrajs; j++)
	{
		cur_samplist = major->samplepts[j];

		if(cur_samplist==NULL) continue;

		glBegin(GL_POINTS);
		for(k=0; k<cur_samplist->nsamples; k++)
		{
			glVertex2f(cur_samplist->samples[k]->gpt[0], cur_samplist->samples[k]->gpt[1]);
		}
		glEnd();
	}
	glPointSize(1.);
}


/*display the major tensor line placement*/
void CGlView::display_minor_samps()
{
	if(minor == NULL)
		return;
	int i, j, k;
	Trajectory *cur_traj;
	glLineWidth(1.5);
	glColor3f(1, 1, 0);

	SamplePtList *cur_samplist;

	glColor3f(0, 1, 0);
	glPointSize(5.);
	/*draw the tensor line*/

	for(j=0; j<minor->evenstreamlines->ntrajs; j++)
	{
		cur_samplist = minor->samplepts[j];

		if(cur_samplist==NULL) continue;

		glBegin(GL_POINTS);
		for(k=0; k<cur_samplist->nsamples; k++)
		{
			glVertex2f(cur_samplist->samples[k]->gpt[0], cur_samplist->samples[k]->gpt[1]);
		}
		glEnd();
	}
	glPointSize(1.);

}


/*     Display the street network    */


void CGlView::display_intersections(GLenum mode)
{
	if(streetnet==NULL ||streetnet->nodelist==NULL) return;
	int i;
	Intersection *curin;

	/*display the end points*/
	//for(i=0; i<streetnet->danglepts->nelems; i++)
	//{
	//	curin=streetnet->danglepts->intersects[i];
	//	glColor3f(0, 1, 0.6);
	//	DrawSolidRect_size(curin->gpos[0], curin->gpos[1], 0.005);
	//	glColor3f(0, 0, 0);
	//	draw_hollow_rect_size(curin->gpos[0], curin->gpos[1], 0.0055);
	//}

	/*display the intersections*/
	for(i=0; i<streetnet->nodelist->nelems; i++)
	{
		curin=streetnet->nodelist->intersects[i];
		if(!curin->endpt)
		{

			if(mode == GL_SELECT)
			{
				glLoadName(NAMEOFINTERSECTS + i);////Assign name for intersection
			}
			glColor3f(1, 1, 1);
			DrawSolidCircle_size(curin->gpos[0], curin->gpos[1], 0.006);
			glColor3f(0, 0, 0);
			draw_hollow_circle_size(curin->gpos[0], curin->gpos[1], 0.0065);
		}
		else
		{ 
			//This is end point
			if(mode == GL_SELECT)
			{
				glLoadName(NAMEOFINTERSECTS + i);////Assign name for intersection
			}

			glColor3f(0, 1, 0.6);
			DrawSolidRect_size(curin->gpos[0], curin->gpos[1], 0.004);
			glColor3f(0, 0, 0);
			draw_hollow_rect_size(curin->gpos[0], curin->gpos[1], 0.0045);
		}
	}
}


/*display the street graph (intersections and their connections)*/
void CGlView::display_streetnet(GLenum mode)
{
	if(streetnet==NULL ||streetnet->nodelist==NULL) return;
	int i;
	Intersection *curin;
	//Graph_Edge *edge;
	StreetGraphEdge *edge;

	/*display the edges first*/
	glLineWidth(1.);
	glColor3f(0, 0, 0);
	for(i=0; i<streetnet->edgelist->nedges; i++)
	{
		edge = streetnet->edgelist->edges[i];

		if(edge->cancel)
			continue;

		glBegin(GL_LINES);
			curin = streetnet->nodelist->intersects[edge->node_index1];
			glVertex2f(curin->gpos[0], curin->gpos[1]);
			curin = streetnet->nodelist->intersects[edge->node_index2];
			glVertex2f(curin->gpos[0], curin->gpos[1]);
		glEnd();
	}

	/*display the intersections*/
	if(!sharedvars.ShowIntersectsOn)
		return;

	for(i=0; i<streetnet->nodelist->nelems; i++)
	{
		curin=streetnet->nodelist->intersects[i];
		if(curin->deleted)
			continue;

		if(!curin->endpt)
		{
			if(mode == GL_SELECT)
			{
				glLoadName(NAMEOFINTERSECTS + i);////Assign name for intersection
			}
			glColor3f(1, 1, 1);
			//DrawSolidCircle_size(curin->gpos[0], curin->gpos[1], 0.006);
			DrawSolidCircle_size(curin->gpos[0], curin->gpos[1], 0.0045/zoom_factor);
			glColor3f(0, 0, 0);
			//draw_hollow_circle_size(curin->gpos[0], curin->gpos[1], 0.0065);
			draw_hollow_circle_size(curin->gpos[0], curin->gpos[1], 0.005/zoom_factor);
		}
		else if(sharedvars.ShowStreetGraphEndPointsOn)
		{
			if(mode == GL_SELECT)
			{
				glLoadName(NAMEOFINTERSECTS + i);////Assign name for intersection
			}
			glColor3f(0, 1, 0.6);
			//DrawSolidRect_size(curin->gpos[0], curin->gpos[1], 0.004);
			DrawSolidRect_size(curin->gpos[0], curin->gpos[1], 0.003/zoom_factor);
			glColor3f(0, 0, 0);
			//draw_hollow_rect_size(curin->gpos[0], curin->gpos[1], 0.0045);
			draw_hollow_rect_size(curin->gpos[0], curin->gpos[1], 0.0035/zoom_factor);
		}

		/*display the four points*/
		//if(!curin->endpt)
		//{
		//glPointSize(3.);
		//glBegin(GL_POINTS);
		//glColor3f(1, 0, 0);
		//glVertex2f(curin->majpt1[0], curin->majpt1[1]);
		//glColor3f(0, 1, 0);
		//glVertex2f(curin->majpt2[0], curin->majpt2[1]);
		//glColor3f(0, 0, 1);
		//glVertex2f(curin->majpt3[0], curin->majpt3[1]);
		//glColor3f(0, 0, 0);
		//glVertex2f(curin->majpt4[0], curin->majpt4[1]);
		//glEnd();
		//glPointSize(1.);
		//}
	}
}

/*
    display the major road network /graph (intersections and their connections)
*/
void CGlView::display_majRoadnet(GLenum mode)
{
	if(majRoadnet==NULL ||majRoadnet->nodelist==NULL) return;
	int i;
	Intersection *curin;
	//Graph_Edge *edge;
	StreetGraphEdge *edge;

	/*display the edges first*/
	glLineWidth(1.);
	glColor3f(0, 0, 0);
	for(i=0; i<majRoadnet->edgelist->nedges; i++)
	{
		edge = majRoadnet->edgelist->edges[i];

		if(edge->cancel)
			continue;

		glBegin(GL_LINES);
			curin = majRoadnet->nodelist->intersects[edge->node_index1];
			glVertex2f(curin->gpos[0], curin->gpos[1]);
			curin = majRoadnet->nodelist->intersects[edge->node_index2];
			glVertex2f(curin->gpos[0], curin->gpos[1]);
		glEnd();
	}

	/*display the intersections*/
	if(!sharedvars.ShowIntersectsOn)
		return;

	for(i=0; i<majRoadnet->nodelist->nelems; i++)
	{
		curin=majRoadnet->nodelist->intersects[i];
		if(curin->deleted)
			continue;

		if(!curin->endpt)
		{
			if(mode == GL_SELECT)
			{
				glLoadName(NAMEOFINTERSECTS + i);////Assign name for intersection
			}
			glColor3f(1, 1, 1);
			DrawSolidCircle_size(curin->gpos[0], curin->gpos[1], 0.0045/zoom_factor);
			glColor3f(0, 0, 0);
			draw_hollow_circle_size(curin->gpos[0], curin->gpos[1], 0.005/zoom_factor);
		}
		else if(sharedvars.ShowStreetGraphEndPointsOn)
		{
			if(mode == GL_SELECT)
			{
				glLoadName(NAMEOFINTERSECTS + i);////Assign name for intersection
			}
			glColor3f(0, 1, 0.6);
			DrawSolidRect_size(curin->gpos[0], curin->gpos[1], 0.003/zoom_factor);
			glColor3f(0, 0, 0);
			draw_hollow_rect_size(curin->gpos[0], curin->gpos[1], 0.0035/zoom_factor);
		}

	}
}
extern StreetVis *roadnetvis;

void CGlView::display_roads(GLenum mode)
{
	if(roadnetvis == NULL)
		return;

	/*clear the back ground*/
	glDisable(GL_TEXTURE_2D);
	glEnable(GL_COLOR_MATERIAL);
	//glClearColor(0.8, .8, 0.9, 1.);
	
	//153 179 204 	
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(0.6, .7, 0.8, 1.);

	int i, j;
	RoadLineSeg *curlineseg;
	OneRoadVis *curroad;

	/**/
	QuadCell *face;
	QuadVertex *v;
	glShadeModel(GL_SMOOTH);
	for(i=0; i<quadmesh->nfaces; i++)
	{
		face = quadmesh->quadcells[i];

		//if(is_not_inland(i))
		//	continue;
		
		glBegin(GL_POLYGON);
		for(j=0; j<face->nverts; j++)
		{
			v=quadmesh->quad_verts[face->verts[j]];
			if(v->inland)
				glColor3f(0.93, 0.93, 0.87);
			else
				glColor3f(0.6, .7, 0.8);
			glVertex2f(v->x, v->y);
		}
		glEnd();

	}

	/*fill up the colors*/
	for(i=0; i<roadnetvis->nmajDirRoads; i++)
	{
		curroad = roadnetvis->majDirRoads[i];

		//if(major->evenstreamlines->trajs[i]->roadtype != HIGHWAY
		//	&& major->evenstreamlines->trajs[i]->roadtype != FREEWAY)
		//	continue;

		if(major->evenstreamlines->trajs[i]->roadtype == HIGHWAY)
			 //255 250 115 
			glColor3f(1, .98, 0.45);
		else if(major->evenstreamlines->trajs[i]->roadtype == MAJOR)
			//242 191 36 
			glColor3f(.95, .75, 0.14);
		else
			glColor3f(1, 1, 1);
		
		for(j=0; j<curroad->nroadlines1-1; j++)
		{
			glBegin(GL_POLYGON);
			curlineseg=curroad->roadline1[j];
			glVertex2f(curlineseg->start[0], curlineseg->start[1]);
			curlineseg=curroad->roadline1[j+1];
			glVertex2f(curlineseg->start[0], curlineseg->start[1]);

			curlineseg=curroad->roadline2[j+1];
			glVertex2f(curlineseg->end[0], curlineseg->end[1]);
			curlineseg=curroad->roadline2[j];
			glVertex2f(curlineseg->end[0], curlineseg->end[1]);
			glEnd();
		}
	}

	for(i=0; i<roadnetvis->nminDirRoads; i++)
	{
		curroad = roadnetvis->minDirRoads[i];

		//if(minor->evenstreamlines->trajs[i]->roadtype != HIGHWAY
		//	&& minor->evenstreamlines->trajs[i]->roadtype != FREEWAY)
		//	continue;

		if(minor->evenstreamlines->trajs[i]->roadtype == HIGHWAY)
			 //255 250 115 
			glColor3f(1, .98, 0.45);
		else if(minor->evenstreamlines->trajs[i]->roadtype == MAJOR)
			//242 191 36 
			glColor3f(.95, .75, 0.14);
		else
			glColor3f(1, 1, 1);
		
		//glBegin(GL_POLYGON);
		//for(j=0; j<curroad->nroadlines1; j++)
		//{
		//	curlineseg=curroad->roadline1[j];
		//	glVertex2f(curlineseg->start[0], curlineseg->start[1]);
		//	glVertex2f(curlineseg->end[0], curlineseg->end[1]);

		//	curlineseg=curroad->roadline2[j];
		//	glVertex2f(curlineseg->end[0], curlineseg->end[1]);
		//	glVertex2f(curlineseg->start[0], curlineseg->start[1]);
		//}
		//glEnd();

		for(j=0; j<curroad->nroadlines1-1; j++)
		{
			glBegin(GL_POLYGON);
			curlineseg=curroad->roadline1[j];
			glVertex2f(curlineseg->start[0], curlineseg->start[1]);
			curlineseg=curroad->roadline1[j+1];
			glVertex2f(curlineseg->start[0], curlineseg->start[1]);

			curlineseg=curroad->roadline2[j+1];
			glVertex2f(curlineseg->end[0], curlineseg->end[1]);
			curlineseg=curroad->roadline2[j];
			glVertex2f(curlineseg->end[0], curlineseg->end[1]);
			glEnd();
		}
	}


	/*draw the enclosure of the roads*/
	glColor3f(0, 0, 0);

	/*major roads*/
	for(i=0; i<roadnetvis->nmajDirRoads; i++)
	{
		curroad = roadnetvis->majDirRoads[i];
		if(major->evenstreamlines->trajs[i]->roadtype == MINOR)
			glLineWidth(1.0);
		//else if(major->evenstreamlines->trajs[i]->roadtype == PLAIN)
		//	glLineWidth(1.2);
		else if(major->evenstreamlines->trajs[i]->roadtype == MAJOR)
			glLineWidth(1.5);
		else
			glLineWidth(2.);

		/*display road 1*/
		for(j=0; j<curroad->nroadlines1; j++)
		{
			curlineseg = curroad->roadline1[j];
			glBegin(GL_LINES);
			glVertex2f(curlineseg->start[0], curlineseg->start[1]);
			glVertex2f(curlineseg->end[0], curlineseg->end[1]);
			glEnd();
		}
		/*display road 2*/
		for(j=0; j<curroad->nroadlines2; j++)
		{
			curlineseg = curroad->roadline2[j];
			glBegin(GL_LINES);
			glVertex2f(curlineseg->start[0], curlineseg->start[1]);
			glVertex2f(curlineseg->end[0], curlineseg->end[1]);
			glEnd();
		}
	}

	/*minor roads*/
	for(i=0; i<roadnetvis->nminDirRoads; i++)
	{
		curroad = roadnetvis->minDirRoads[i];

		if(minor->evenstreamlines->trajs[i]->roadtype == MINOR)
			glLineWidth(1.0);
		//else if(minor->evenstreamlines->trajs[i]->roadtype == PLAIN)
		//	glLineWidth(1.2);
		else if(minor->evenstreamlines->trajs[i]->roadtype == MAJOR)
			glLineWidth(1.5);
		else
			glLineWidth(2.);

		/*display road 1*/
		for(j=0; j<curroad->nroadlines1; j++)
		{
			curlineseg = curroad->roadline1[j];
			glBegin(GL_LINES);
			glVertex2f(curlineseg->start[0], curlineseg->start[1]);
			glVertex2f(curlineseg->end[0], curlineseg->end[1]);
			glEnd();
		}
		/*display road 2*/
		for(j=0; j<curroad->nroadlines2; j++)
		{
			curlineseg = curroad->roadline2[j];
			glBegin(GL_LINES);
			glVertex2f(curlineseg->start[0], curlineseg->start[1]);
			glVertex2f(curlineseg->end[0], curlineseg->end[1]);
			glEnd();
		}
	}
}


void CGlView::display_roads_width(GLenum mode)
{
	/*clear the back ground*/
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
	//glEnable(GL_COLOR_MATERIAL);
	//glClearColor(0.8, .8, 0.9, 1.);
	
	//153 179 204 	
	glClearColor(0.6, .7, 0.8, 1.);
	glClear(GL_COLOR_BUFFER_BIT);

	int i, j;
	RoadLineSeg *curlineseg;
	OneRoadVis *curroad;

	/**/


	if(flag_loadmap)
		render_a_map(streetmapbackground);
	else
	{
		QuadCell *face;
		QuadVertex *v;
		glShadeModel(GL_SMOOTH);
		for(i=0; i<quadmesh->nfaces; i++)
		{
			face = quadmesh->quadcells[i];

			glBegin(GL_POLYGON);
			for(j=0; j<face->nverts; j++)
			{
				v=quadmesh->quad_verts[face->verts[j]];
				if(v->inland)
					//glColor3f(0.83, 0.82, 0.8);
					//237 234 226 
					glColor3f(0.93, 0.93, 0.87);
				else
					glColor3f(0.6, .7, 0.8);
				glVertex2f(v->x, v->y);
			}
			glEnd();
		}
	}

	glDisable(GL_TEXTURE_2D);

	/*display the back*/
	if(major==NULL||minor==NULL)
		return;

	//for(i=0; i<major->evenstreamlines->ntrajs; i++)
	//	display_road_width_back(i, 0, false);
	//
	//for(i=0; i<minor->evenstreamlines->ntrajs; i++)
	//	display_road_width_back(i, 0, true);
	//
	//for(i=0; i<minor->evenstreamlines->ntrajs; i++)
	//	display_road_width_front(i, 0, true);

	//for(i=major->evenstreamlines->ntrajs-1; i>=0; i--)
	//	display_road_width_front(i, 0, false);
	//


	/* use the graph edge to visualize the obtained roads */
	if(streetnet==NULL)
		return;

	/* rend the back of the roads */
	display_road_graph_back(0, MINOR);
	display_road_graph_back(0, MAJOR);
	display_road_graph_back(0, HIGHWAY);

	/* rend the front of the roads */
	display_road_graph_front(0, MINOR);
	display_road_graph_front(0, MAJOR);
	display_road_graph_front(0, HIGHWAY);
	
	//glEnable(GL_DEPTH_TEST);
}


/*
Currently, we still use the obtained tensor lines to visualize the road network
Later, we should consider to use the new sampled points to visualize the roads
*/
void CGlView::display_road_width_back(int id, double width, bool majormin)
{
	int i;
	Trajectory *cur_traj;
		
	if(!majormin)
		cur_traj = major->evenstreamlines->trajs[id];
	else
		cur_traj = minor->evenstreamlines->trajs[id];

	//glLineWidth(4);
	if(cur_traj->roadtype==MINOR)
		glLineWidth(VisMinRoadWidth);
	else if(cur_traj->roadtype==HIGHWAY)
		glLineWidth(VisMinRoadWidth+2);
	else if(cur_traj->roadtype==MAJOR)
		glLineWidth(VisMinRoadWidth+1);

	glColor3f(0, 0, 0);
	for(i=0; i<cur_traj->nlinesegs; i++)
	{
		glBegin(GL_LINES);
		glVertex2f(cur_traj->linesegs[i].gstart[0], cur_traj->linesegs[i].gstart[1]);
		glVertex2f(cur_traj->linesegs[i].gend[0], cur_traj->linesegs[i].gend[1]);
		glEnd();
	}
}

void CGlView::display_road_width_front(int id, double width, bool majormin)
{
	int i;
	Trajectory *cur_traj;
		
	if(!majormin)
		cur_traj = major->evenstreamlines->trajs[id];
	else
		cur_traj = minor->evenstreamlines->trajs[id];


	//glLineWidth(2.);
	if(cur_traj->roadtype==MINOR)
		glLineWidth(VisMinRoadWidth-2);
	else if(cur_traj->roadtype==HIGHWAY)
		glLineWidth(VisMinRoadWidth);
	else if(cur_traj->roadtype==MAJOR)
		glLineWidth(VisMinRoadWidth-1);

	if(cur_traj->roadtype==HIGHWAY)
		glColor3f(.95, .75, 0.14);
	else if(cur_traj->roadtype==MAJOR)
		glColor3f(1, .98, 0.45);
	else
		glColor3f(1, 1, 1);

	for(i=0; i<cur_traj->nlinesegs; i++)
	{
		glBegin(GL_LINES);
		glVertex2f(cur_traj->linesegs[i].gstart[0], cur_traj->linesegs[i].gstart[1]);
		glVertex2f(cur_traj->linesegs[i].gend[0], cur_traj->linesegs[i].gend[1]);
		glEnd();
	}
}



/*
    We visualize the road network using the samples along the tensor lines
*/

void CGlView::display_road_graph_back(double width, unsigned char roadtype)
{
	int i,j;
	StreetGraphEdge *edge;
		
	if(roadtype==MINOR)
	{
		glLineWidth(VisMinRoadWidth);
		glPointSize(VisMinRoadWidth-1);
	}
	else if(roadtype==HIGHWAY)
	{
		glLineWidth(VisMinRoadWidth+3);
		glPointSize(VisMinRoadWidth+1);
	}
	else if(roadtype==MAJOR)
	{
		glLineWidth(VisMinRoadWidth+2);
		glPointSize(VisMinRoadWidth);
	}

	glColor3f(0, 0, 0);
	for(i=0; i<streetnet->edgelist->nedges; i++)
	{
		edge=streetnet->edgelist->edges[i];
		if(edge->roadtype!=roadtype)
			continue;

		if(edge->cancel)
			continue;

		//glBegin(GL_LINE_STRIP);
		//for(j=0;j<edge->ninter_pts;j++)
		//{
		//	glVertex2f(edge->inter_pts[j]->x, edge->inter_pts[j]->y);
		//}
		//glEnd();

		//glBegin(GL_LINES);
		//for(j=0;j<edge->ninter_pts-1;j++)
		//{
		//	glVertex2f(edge->inter_pts[j]->x, edge->inter_pts[j]->y);
		//	glVertex2f(edge->inter_pts[j+1]->x, edge->inter_pts[j+1]->y);
		//}
		//glEnd();
		
		double size;
		if(roadtype==MINOR)
			size=(VisMinRoadWidth)/1200.;
		else if(roadtype==HIGHWAY)
			size=(VisMinRoadWidth+2)/1200.;
		else if(roadtype==MAJOR)
			size=(VisMinRoadWidth+1)/1100.;
		
		size=size/zoom_factor;

		icVector2 line_dir, norm;
		double pre_up[2], pre_low[2];

		line_dir.entry[0]=edge->inter_pts[1]->x-edge->inter_pts[0]->x;
		line_dir.entry[1]=edge->inter_pts[1]->y-edge->inter_pts[0]->y;
		normalize(line_dir);
		norm.entry[0]=-line_dir.entry[1];
		norm.entry[1]=line_dir.entry[0];
		norm = size/1.3* norm;

		pre_up[0]=edge->inter_pts[0]->x+norm.entry[0];
		pre_up[1]=edge->inter_pts[0]->y+norm.entry[1];
		pre_low[0]=edge->inter_pts[0]->x-norm.entry[0];
		pre_low[1]=edge->inter_pts[0]->y-norm.entry[1];

		for(j=1;j<edge->ninter_pts;j++)
		{
			if(j<edge->ninter_pts-1)
			{
				line_dir.entry[0]=edge->inter_pts[j+1]->x-edge->inter_pts[j-1]->x;
				line_dir.entry[1]=edge->inter_pts[j+1]->y-edge->inter_pts[j-1]->y;
			}
			else{
				line_dir.entry[0]=edge->inter_pts[j]->x-edge->inter_pts[j-1]->x;
				line_dir.entry[1]=edge->inter_pts[j]->y-edge->inter_pts[j-1]->y;
			}
			normalize(line_dir);
			norm.entry[0]=-line_dir.entry[1];
			norm.entry[1]=line_dir.entry[0];

			norm = size/1.3* norm;

			glBegin(GL_POLYGON);
			glVertex2f(pre_up[0], pre_up[1]);
			glVertex2f(edge->inter_pts[j]->x+norm.entry[0],
				edge->inter_pts[j]->y+norm.entry[1]);
			glVertex2f(edge->inter_pts[j]->x-norm.entry[0],
				edge->inter_pts[j]->y-norm.entry[1]);
			glVertex2f(pre_low[0], pre_low[1]);
			glEnd();

			pre_up[0]=edge->inter_pts[j]->x+norm.entry[0];
			pre_up[1]=edge->inter_pts[j]->y+norm.entry[1];
			pre_low[0]=edge->inter_pts[j]->x-norm.entry[0];
			pre_low[1]=edge->inter_pts[j]->y-norm.entry[1];
		}

		/*  display the end point */
		if(roadtype==MINOR)
			size=(VisMinRoadWidth)/1520.;
		else if(roadtype==HIGHWAY)
			size=(VisMinRoadWidth+3)/1590.;
		else if(roadtype==MAJOR)
			size=(VisMinRoadWidth+2)/1600.;

		Intersection *intersect1=streetnet->nodelist->intersects[edge->node_index1];
		Intersection *intersect2=streetnet->nodelist->intersects[edge->node_index2];
		

		if(intersect1->endpt)
		{
			//glBegin(GL_POINTS);
			//glVertex2f(intersect1->gpos[0],intersect1->gpos[1]);
			//glEnd();
			DrawSolidCircle_size(intersect1->gpos[0], intersect1->gpos[1], size/zoom_factor);
		}

		if(intersect2->endpt)
		{
			//glBegin(GL_POINTS);
			//glVertex2f(intersect2->gpos[0],intersect2->gpos[1]);
			//glEnd();
			DrawSolidCircle_size(intersect2->gpos[0], intersect2->gpos[1], size/zoom_factor);
		}
	}

}

void CGlView::display_road_graph_front(double width, unsigned char roadtype)
{
	int i,j;
	StreetGraphEdge *edge;

	if(roadtype==MINOR)
	{
		glLineWidth(VisMinRoadWidth-2);
		glPointSize(VisMinRoadWidth-2.5);
	}
	else if(roadtype==HIGHWAY)
	{
		glLineWidth(VisMinRoadWidth);
		glPointSize(VisMinRoadWidth-0.5);
	}
	else if(roadtype==MAJOR)
	{
		glLineWidth(VisMinRoadWidth-0.5);
		glPointSize(VisMinRoadWidth-1.5);
	}

	if(roadtype==HIGHWAY)
		glColor3f(.95, .75, 0.14);
	else if(roadtype==MAJOR)
		glColor3f(1, .98, 0.45);
	else
		glColor3f(1, 1, 1);

	for(i=0; i<streetnet->edgelist->nedges; i++)
	{
		edge=streetnet->edgelist->edges[i];
		if(edge->roadtype!=roadtype)
			continue;
		if(edge->cancel)
			continue;


		//glBegin(GL_LINES);
		//for(j=0;j<edge->ninter_pts-1;j++)
		//{
		//	glVertex2f(edge->inter_pts[j]->x, edge->inter_pts[j]->y);
		//	glVertex2f(edge->inter_pts[j+1]->x, edge->inter_pts[j+1]->y);
		//}
		//glEnd();

		double size;
		if(roadtype==MINOR)
			size=(VisMinRoadWidth-2)/1250.;
		else if(roadtype==HIGHWAY)
			size=(VisMinRoadWidth)/1250.;
		else if(roadtype==MAJOR)
			size=(VisMinRoadWidth-1)/1250.;
		
		size=size/zoom_factor;

		icVector2 line_dir, norm;
		double pre_up[2], pre_low[2];

		line_dir.entry[0]=edge->inter_pts[1]->x-edge->inter_pts[0]->x;
		line_dir.entry[1]=edge->inter_pts[1]->y-edge->inter_pts[0]->y;
		normalize(line_dir);
		norm.entry[0]=-line_dir.entry[1];
		norm.entry[1]=line_dir.entry[0];
		norm = size/1.3* norm;

		pre_up[0]=edge->inter_pts[0]->x+norm.entry[0];
		pre_up[1]=edge->inter_pts[0]->y+norm.entry[1];
		pre_low[0]=edge->inter_pts[0]->x-norm.entry[0];
		pre_low[1]=edge->inter_pts[0]->y-norm.entry[1];

		for(j=1;j<edge->ninter_pts;j++)
		{
			if(j<edge->ninter_pts-1)
			{
				line_dir.entry[0]=edge->inter_pts[j+1]->x-edge->inter_pts[j-1]->x;
				line_dir.entry[1]=edge->inter_pts[j+1]->y-edge->inter_pts[j-1]->y;
			}
			else{
				line_dir.entry[0]=edge->inter_pts[j]->x-edge->inter_pts[j-1]->x;
				line_dir.entry[1]=edge->inter_pts[j]->y-edge->inter_pts[j-1]->y;
			}

			normalize(line_dir);
			norm.entry[0]=-line_dir.entry[1];
			norm.entry[1]=line_dir.entry[0];

			norm = size/1.3* norm;

			glBegin(GL_POLYGON);
			glVertex2f(pre_up[0], pre_up[1]);
			glVertex2f(edge->inter_pts[j]->x+norm.entry[0],
				edge->inter_pts[j]->y+norm.entry[1]);
			glVertex2f(edge->inter_pts[j]->x-norm.entry[0],
				edge->inter_pts[j]->y-norm.entry[1]);
			glVertex2f(pre_low[0], pre_low[1]);
			glEnd();

			pre_up[0]=edge->inter_pts[j]->x+norm.entry[0];
			pre_up[1]=edge->inter_pts[j]->y+norm.entry[1];
			pre_low[0]=edge->inter_pts[j]->x-norm.entry[0];
			pre_low[1]=edge->inter_pts[j]->y-norm.entry[1];
		}


		/*  display the end point */
		if(roadtype==MINOR)
			size=(VisMinRoadWidth-2)/1530.;
		else if(roadtype==HIGHWAY)
			size=(VisMinRoadWidth)/1520.;
		else if(roadtype==MAJOR)
			size=(VisMinRoadWidth-1)/1500.;

		Intersection *intersect1=streetnet->nodelist->intersects[edge->node_index1];
		Intersection *intersect2=streetnet->nodelist->intersects[edge->node_index2];
		

		if(intersect1->endpt)
		{
			//glBegin(GL_POINTS);
			//glVertex2f(intersect1->gpos[0],intersect1->gpos[1]);
			//glEnd();
			DrawSolidCircle_size(intersect1->gpos[0], intersect1->gpos[1], size/zoom_factor);
		}

		if(intersect2->endpt)
		{
			//glBegin(GL_POINTS);
			//glVertex2f(intersect2->gpos[0],intersect2->gpos[1]);
			//glEnd();
			DrawSolidCircle_size(intersect2->gpos[0], intersect2->gpos[1], size/zoom_factor);
		}
	}

}


/*
 bug:
 This is only a temporary function for major road 'visualization' !!!
 1/13/2008
*/
void CGlView::display_majRoad_googlestyle_back_temp(GLenum mode, bool majormin)
{
	int i,j,k;

	if(major_level1==NULL || minor_level1==NULL)
		return;
		
	glLineWidth(VisMinRoadWidth+2);

	Trajectory *curtraj;
	LineSeg *curline;
	EvenStreamlinePlace *curplace;

	if(!majormin) /*  major direction  */
	{
		curplace=major_level1;
	}
	else
		curplace=minor_level1;

	double size=(VisMinRoadWidth+1)/1000.;
	size=size/zoom_factor;

	glColor3f(0, 0, 0);

	for(i=0; i<curplace->evenstreamlines->ntrajs; i++)
	{
		curtraj=curplace->evenstreamlines->trajs[i];

		size=(VisMinRoadWidth+1)/1000.;
		size=size/zoom_factor;

		icVector2 line_dir, norm;
		double pre_up[2], pre_low[2];

		curline=&curtraj->linesegs[0];

		line_dir.entry[0]=curline->gend[0]-curline->gstart[0];
		line_dir.entry[1]=curline->gend[1]-curline->gstart[1];
		normalize(line_dir);
		norm.entry[0]=-line_dir.entry[1];
		norm.entry[1]=line_dir.entry[0];
		norm = size/1.3* norm;

		pre_up[0]= curline->gstart[0]+norm.entry[0];
		pre_up[1]= curline->gstart[1]+norm.entry[1];
		pre_low[0]= curline->gstart[0]-norm.entry[0];
		pre_low[1]= curline->gstart[1]-norm.entry[1];

		for(j=1;j<curtraj->nlinesegs;j++)
		{
			if(j<curtraj->nlinesegs-1)
			{
				line_dir.entry[0]=curtraj->linesegs[j+1].gend[0]-
					curtraj->linesegs[j-1].gstart[0];
				line_dir.entry[1]=curtraj->linesegs[j+1].gend[1]-
					curtraj->linesegs[j-1].gstart[1];
			}
			else{
				line_dir.entry[0]=curtraj->linesegs[j].gend[0]-
					curtraj->linesegs[j-1].gstart[0];
				line_dir.entry[1]=curtraj->linesegs[j].gend[1]-
					curtraj->linesegs[j-1].gstart[1];
			}
			normalize(line_dir);
			norm.entry[0]=-line_dir.entry[1];
			norm.entry[1]=line_dir.entry[0];

			norm = size/1.3* norm;

			curline=&curtraj->linesegs[j];

			glBegin(GL_POLYGON);
			glVertex2f(pre_up[0], pre_up[1]);
			glVertex2f(curline->gstart[0]+norm.entry[0],
				curline->gstart[1]+norm.entry[1]);
			glVertex2f(curline->gstart[0]-norm.entry[0],
				curline->gstart[1]-norm.entry[1]);
			glVertex2f(pre_low[0], pre_low[1]);
			glEnd();

			pre_up[0]=curline->gstart[0]+norm.entry[0];
			pre_up[1]=curline->gstart[1]+norm.entry[1];
			pre_low[0]=curline->gstart[0]-norm.entry[0];
			pre_low[1]=curline->gstart[1]-norm.entry[1];
		}

		/*   consider last line segment  */
		glBegin(GL_POLYGON);
		glVertex2f(pre_up[0], pre_up[1]);
		glVertex2f(curline->gend[0]+norm.entry[0],
			curline->gend[1]+norm.entry[1]);
		glVertex2f(curline->gend[0]-norm.entry[0],
			curline->gend[1]-norm.entry[1]);
		glVertex2f(pre_low[0], pre_low[1]);
		glEnd();

		
		
		/*  display the end point */
		size=(VisMinRoadWidth+2)/1610.;

		curline=&curtraj->linesegs[0];
		DrawSolidCircle_size(curline->gstart[0], curline->gstart[1], size/zoom_factor);

		curline=&curtraj->linesegs[curtraj->nlinesegs-1];
		DrawSolidCircle_size(curline->gend[0], curline->gend[1], size/zoom_factor);
	}
}

void CGlView::display_majRoad_googlestyle_front_temp(GLenum mode, bool majormin)
{
	int i,j,k;

	if(major_level1==NULL || minor_level1==NULL)
		return;
		
	glLineWidth(VisMinRoadWidth+2);

	Trajectory *curtraj;
	LineSeg *curline;
	EvenStreamlinePlace *curplace;

	if(!majormin) /*  major direction  */
	{
		curplace=major_level1;
	}
	else
		curplace=minor_level1;

	glColor3f(1, .98, 0.45);

	double size=(VisMinRoadWidth+1)/1400;
	size=size/zoom_factor;

	for(i=0; i<curplace->evenstreamlines->ntrajs; i++)
	{
		curtraj=curplace->evenstreamlines->trajs[i];

		size=(VisMinRoadWidth+1)/1400.;
		size=size/zoom_factor;

		icVector2 line_dir, norm;
		double pre_up[2], pre_low[2];

		curline=&curtraj->linesegs[0];

		line_dir.entry[0]=curline->gend[0]-curline->gstart[0];
		line_dir.entry[1]=curline->gend[1]-curline->gstart[1];
		normalize(line_dir);
		norm.entry[0]=-line_dir.entry[1];
		norm.entry[1]=line_dir.entry[0];
		norm = size/1.3* norm;

		pre_up[0]= curline->gstart[0]+norm.entry[0];
		pre_up[1]= curline->gstart[1]+norm.entry[1];
		pre_low[0]= curline->gstart[0]-norm.entry[0];
		pre_low[1]= curline->gstart[1]-norm.entry[1];

		for(j=1;j<curtraj->nlinesegs;j++)
		{
			if(j<curtraj->nlinesegs-1)
			{
				line_dir.entry[0]=curtraj->linesegs[j+1].gend[0]-
					curtraj->linesegs[j-1].gstart[0];
				line_dir.entry[1]=curtraj->linesegs[j+1].gend[1]-
					curtraj->linesegs[j-1].gstart[1];
			}
			else{
				line_dir.entry[0]=curtraj->linesegs[j].gend[0]-
					curtraj->linesegs[j-1].gstart[0];
				line_dir.entry[1]=curtraj->linesegs[j].gend[1]-
					curtraj->linesegs[j-1].gstart[1];
			}
			normalize(line_dir);
			norm.entry[0]=-line_dir.entry[1];
			norm.entry[1]=line_dir.entry[0];

			norm = size/1.3* norm;

			curline=&curtraj->linesegs[j];

			glBegin(GL_POLYGON);
			glVertex2f(pre_up[0], pre_up[1]);
			glVertex2f(curline->gstart[0]+norm.entry[0],
				curline->gstart[1]+norm.entry[1]);
			glVertex2f(curline->gstart[0]-norm.entry[0],
				curline->gstart[1]-norm.entry[1]);
			glVertex2f(pre_low[0], pre_low[1]);
			glEnd();

			pre_up[0]=curline->gstart[0]+norm.entry[0];
			pre_up[1]=curline->gstart[1]+norm.entry[1];
			pre_low[0]=curline->gstart[0]-norm.entry[0];
			pre_low[1]=curline->gstart[1]-norm.entry[1];
		}

		/*   consider last line segment  */
		glBegin(GL_POLYGON);
		glVertex2f(pre_up[0], pre_up[1]);
		glVertex2f(curline->gend[0]+norm.entry[0],
			curline->gend[1]+norm.entry[1]);
		glVertex2f(curline->gend[0]-norm.entry[0],
			curline->gend[1]-norm.entry[1]);
		glVertex2f(pre_low[0], pre_low[1]);
		glEnd();

		
		
		/*  display the end point */
		size=(VisMinRoadWidth-1)/1500.;

		curline=&curtraj->linesegs[0];
		DrawSolidCircle_size(curline->gstart[0], curline->gstart[1], size/zoom_factor);

		curline=&curtraj->linesegs[curtraj->nlinesegs-1];
		DrawSolidCircle_size(curline->gend[0], curline->gend[1], size/zoom_factor);
	}
}


void CGlView::display_majRoad_googlestyle(GLenum mode)
{
	/*clear the back ground*/
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
	
	//153 179 204 	
	glClearColor(0.6, .7, 0.8, 1.);
	glClear(GL_COLOR_BUFFER_BIT);

	int i, j;

	if(flag_loadmap)
		render_a_map(streetmapbackground);
	else
	{
		QuadCell *face;
		QuadVertex *v;
		glShadeModel(GL_SMOOTH);
		for(i=0; i<quadmesh->nfaces; i++)
		{
			face = quadmesh->quadcells[i];

			glBegin(GL_POLYGON);
			for(j=0; j<face->nverts; j++)
			{
				v=quadmesh->quad_verts[face->verts[j]];
				if(v->inland)
					//glColor3f(0.8, 0.8, 0.8);
					//237 234 226 
					glColor3f(0.93, 0.93, 0.87);
				else
					glColor3f(0.6, .7, 0.8);
				glVertex2f(v->x, v->y);
			}
			glEnd();
		}
	}

	glDisable(GL_TEXTURE_2D);

	display_majRoad_googlestyle_back_temp(mode, false);
	display_majRoad_googlestyle_back_temp(mode, true);
	display_majRoad_googlestyle_front_temp(mode, false);
	display_majRoad_googlestyle_front_temp(mode, true);
}


/*
    We display the streets use the obtained network
*/
void CGlView::display_road_use_network()
{
	int i, j;
	StreetGraphEdge *edge;
	Intersection *intersect1, *intersect2;
	/*clear the back ground*/
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	//glEnable(GL_COLOR_MATERIAL);
	//glClearColor(0.8, .8, 0.9, 1.);
	
	//153 179 204 	
	glClearColor(0.6, .7, 0.8, 1.);
	glClear(GL_COLOR_BUFFER_BIT);

	/**/


	if(flag_loadmap)
		render_a_map(streetmapbackground);
	else
	{
		QuadCell *face;
		QuadVertex *v;
		glShadeModel(GL_SMOOTH);
		for(i=0; i<quadmesh->nfaces; i++)
		{
			face = quadmesh->quadcells[i];

			glBegin(GL_POLYGON);
			for(j=0; j<face->nverts; j++)
			{
				v=quadmesh->quad_verts[face->verts[j]];
				if(v->inland)
					//glColor3f(0.8, 0.8, 0.8);
					//237 234 226 
					glColor3f(0.93, 0.93, 0.87);
				else
					glColor3f(0.6, .7, 0.8);
				glVertex2f(v->x, v->y);
			}
			glEnd();
		}
	}

	glDisable(GL_TEXTURE_2D);

	if(streetnet==NULL) return;

	glColor3f(0, 0, 0);
	for(i=0; i<streetnet->edgelist->nedges; i++)
	{
		edge=streetnet->edgelist->edges[i];

		if(edge->cancel)
			continue;

		intersect1=streetnet->nodelist->intersects[edge->node_index1];
		intersect2=streetnet->nodelist->intersects[edge->node_index2];

		if(edge->roadtype==MINOR)
			glLineWidth(VisMinRoadWidth-2);
		else if(edge->roadtype==HIGHWAY)
			glLineWidth(VisMinRoadWidth);
		else if(edge->roadtype==MAJOR)
			glLineWidth(VisMinRoadWidth-1);


		glBegin(GL_LINES);
		glVertex2f(intersect1->gpos[0],intersect1->gpos[1]);
		glVertex2f(intersect2->gpos[0],intersect2->gpos[1]);
		glEnd();
	}
}

void CGlView::display_road_graph_linestyle()
{
	//int i,j;
	//StreetGraphEdge *edge;

	//if(roadtype==MINOR)
	//{
	//	glLineWidth(VisMinRoadWidth-2);
	//}
	//else if(roadtype==HIGHWAY)
	//{
	//	glLineWidth(VisMinRoadWidth);
	//}
	//else if(roadtype==MAJOR)
	//{
	//	glLineWidth(VisMinRoadWidth-1);
	//}

	//glColor3f(0,0,0);
	//for(i=0; i<streetnet->edgelist->nedges; i++)
	//{
	//	edge=streetnet->edgelist->edges[i];
	//	if(edge->roadtype!=roadtype)
	//		continue;
	//	if(edge->cancel)
	//		continue;

	//	glBegin(GL_LINES);
	//	for(j=0;j<edge->ninter_pts-1;j++)
	//	{
	//		glVertex2f(edge->inter_pts[j]->x, edge->inter_pts[j]->y);
	//		glVertex2f(edge->inter_pts[j+1]->x, edge->inter_pts[j+1]->y);
	//	}
	//	glEnd();

	//}

	int i, j;
	StreetGraphEdge *edge;
	/*clear the back ground*/
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	
	//153 179 204 	
	glClearColor(0.6, .7, 0.8, 1.);
	glClear(GL_COLOR_BUFFER_BIT);

	/**/


	if(flag_loadmap)
		render_a_map(streetmapbackground);
	else
	{
		QuadCell *face;
		QuadVertex *v;
		glShadeModel(GL_SMOOTH);
		for(i=0; i<quadmesh->nfaces; i++)
		{
			face = quadmesh->quadcells[i];

			glBegin(GL_POLYGON);
			for(j=0; j<face->nverts; j++)
			{
				v=quadmesh->quad_verts[face->verts[j]];
				if(v->inland)
					//glColor3f(0.8, 0.8, 0.8);
					//237 234 226 
					glColor3f(0.93, 0.93, 0.87);
				else
					glColor3f(0.6, .7, 0.8);
				glVertex2f(v->x, v->y);
			}
			glEnd();
		}
	}

	glDisable(GL_TEXTURE_2D);

	if(streetnet==NULL) return;

	glColor3f(0, 0, 0);
	for(i=0; i<streetnet->edgelist->nedges; i++)
	{
		edge=streetnet->edgelist->edges[i];

		if(edge->cancel)
			continue;

		if(edge->roadtype==MINOR)
			glLineWidth(VisMinRoadWidth-2);
		else if(edge->roadtype==HIGHWAY)
			glLineWidth(VisMinRoadWidth);
		else if(edge->roadtype==MAJOR)
			glLineWidth(VisMinRoadWidth-1);

		glBegin(GL_LINE_STRIP);
		for(j=0;j<edge->ninter_pts;j++)
		{
			glVertex2f(edge->inter_pts[j]->x, edge->inter_pts[j]->y);
		}
		glEnd();
	}

}





/*for the brush stroke editing of tensor field design*/
ctr_point **brushpts = NULL;
int nbrushpts = 0;
int curMaxNumBrushPts = 1500;

/*visualize the brush stroke when under the brush editing mode*/

void CGlView::display_brush()
{
	if(brushpts == NULL)
		return;

	int i;
	glColor3f(1, 0, 0);
	glLineWidth(3.);
	glBegin(GL_LINE_STRIP);
	for(i=0; i<nbrushpts; i++)
	{
		glVertex2f(brushpts[i]->x, brushpts[i]->y);
	}
	glEnd();
	glLineWidth(1.);
}


extern SketchList *brushes;
extern TrajectoryList *sketchlines;
extern StreetNet *sketchnet;

/*visualize the brush-based sketch stroke when under the brush editing mode*/

void CGlView::display_brush_sketch_thin(GLenum mode)
{
	if(brushes == NULL)
		return;

	int i, j;

	/*display the obtained trajectories  11/22/2007*/
	for(j=0;j<sketchlines->ntrajs;j++)
	{
		//if(mode==GL_SELECT)
		//	glLoadName(NAMEOFBRUSHES+j);

		glLineWidth(2.);
		glColor3f(1, 1, 1);
		glBegin(GL_LINES);
		for(i=0; i<sketchlines->trajs[j]->nlinesegs; i++)
		{
			glVertex2f(sketchlines->trajs[j]->linesegs[i].gstart[0], 
				sketchlines->trajs[j]->linesegs[i].gstart[1]);
			glVertex2f(sketchlines->trajs[j]->linesegs[i].gend[0], 
				sketchlines->trajs[j]->linesegs[i].gend[1]);
		}
		glEnd();
		
	}
	glLineWidth(1.);

	for(j=0;j<brushes->nbrushes+1;j++)
	{
		if(mode==GL_SELECT)
			glLoadName(NAMEOFBRUSHES+j);

		//glLineWidth(7.);
		glLineWidth(1.);
		glColor3f(0, 0, 0);
		glBegin(GL_LINE_STRIP);
		for(i=0; i<brushes->brushlist[j].nelems; i++)
		{
			glVertex2f(brushes->brushlist[j].brushpts[i]->x, 
				brushes->brushlist[j].brushpts[i]->y);
		}
		glEnd();
		
		//glLineWidth(5.);
		glLineWidth(1.);
		glColor3f(.95, .75, 0.14);
		glBegin(GL_LINE_STRIP);
		for(i=0; i<brushes->brushlist[j].nelems; i++)
		{
			glVertex2f(brushes->brushlist[j].brushpts[i]->x, 
				brushes->brushlist[j].brushpts[i]->y);
		}
		glEnd();
	}

	if(sketchlines==NULL)
		return;

	for(j=0;j<sketchlines->ntrajs;j++)
	{
		//if(mode==GL_SELECT)
		//	glLoadName(NAMEOFBRUSHES+j);

		//glLineWidth(7.);
		glLineWidth(1.);
		glColor3f(0, 0, 0);
		glBegin(GL_LINES);
		for(i=0; i<sketchlines->trajs[j]->nlinesegs; i++)
		{
			glVertex2f(sketchlines->trajs[j]->linesegs[i].gstart[0], 
				sketchlines->trajs[j]->linesegs[i].gstart[1]);
			glVertex2f(sketchlines->trajs[j]->linesegs[i].gend[0], 
				sketchlines->trajs[j]->linesegs[i].gend[1]);
		}
		glEnd();
		
		//glLineWidth(5.);
		glLineWidth(1.);
		glColor3f(.95, .75, 0.14);
		glBegin(GL_LINES);
		for(i=0; i<sketchlines->trajs[j]->nlinesegs; i++)
		{
			if(sketchlines->trajs[j]->linesegs[i].Triangle_ID==5101)
				glColor3f(1,0,0);
			else
			{
				if(i%2==0)
					glColor3f(.95, .75, 0.14);
				else
					glColor3f(.14, .75, 0.14);
			}

			glVertex2f(sketchlines->trajs[j]->linesegs[i].gstart[0], 
				sketchlines->trajs[j]->linesegs[i].gstart[1]);
			glVertex2f(sketchlines->trajs[j]->linesegs[i].gend[0], 
				sketchlines->trajs[j]->linesegs[i].gend[1]);
		}
		glEnd();		
	}
	glLineWidth(1.);

	///* display the obtain graph if not empty*/
	if(sharedvars.ShowSegmentGraphOn)
	{
		if(sketchnet==NULL||sketchnet->nodelist==NULL||sketchnet->nodelist->nelems==0) return;
		Intersection *curin;
		StreetGraphEdge *edge;

		///*display the edges first*/
		glLineWidth(1.);
		glColor3f(0, 0, 0);
		for(i=0; i<sketchnet->edgelist->nedges; i++)
		{
			edge = sketchnet->edgelist->edges[i];
			glBegin(GL_LINES);
				curin = sketchnet->nodelist->intersects[edge->node_index1];
				glVertex2f(curin->gpos[0], curin->gpos[1]);
				curin = sketchnet->nodelist->intersects[edge->node_index2];
				glVertex2f(curin->gpos[0], curin->gpos[1]);
			glEnd();
		}

		/*display the intersections*/
		for(i=0; i<sketchnet->nodelist->nelems; i++)
		{
			curin=sketchnet->nodelist->intersects[i];
			if(!curin->endpt)
			{
				if(mode == GL_SELECT)
				{
					glLoadName(NAMEOFINTERSECTS + i);////Assign name for intersection
				}
				glColor3f(1, 1, 1);
				DrawSolidCircle_size(curin->gpos[0], curin->gpos[1], 0.0045/zoom_factor);
				glColor3f(0, 0, 0);
				draw_hollow_circle_size(curin->gpos[0], curin->gpos[1], 0.005/zoom_factor);
			}
			else if(sharedvars.ShowStreetGraphEndPointsOn)
			{
				if(mode == GL_SELECT)
				{
					glLoadName(NAMEOFINTERSECTS + i);////Assign name for intersection
				}
				glColor3f(0, 1, 0.6);
				DrawSolidRect_size(curin->gpos[0], curin->gpos[1], 0.003/zoom_factor);
				glColor3f(0, 0, 0);
				draw_hollow_rect_size(curin->gpos[0], curin->gpos[1], 0.0035/zoom_factor);
			}
		}
	}

//extern int * boundarycells;
//extern int nboundarycells;
//    glColor3f(1,1,0);
//    for(i=1;i<nboundarycells;i++)
//	{
//		QuadCell *qc=quadmesh->quadcells[boundarycells[i]];
//		glBegin(GL_LINE_LOOP);
//		for(j=0;j<qc->nverts;j++)
//		{
//			QuadVertex *cv=quadmesh->quad_verts[qc->verts[j]];
//			glVertex2f(cv->x, cv->y);
//		}
//		glEnd();
//	}
//
//
//extern int *boundvertlist;
//extern int nboundverts;
//		/*  display the vertices in the obtain narrow band  */
//    glBegin(GL_POINTS);
//	glPointSize(5.);
//	glColor3f(1,1,1);
//	for(i=0;i<nboundverts;i++)
//	{
//		QuadVertex *cv=quadmesh->quad_verts[boundvertlist[i]];
//		glVertex2f(cv->x, cv->y);
//	}
//	glEnd();

}


/*   display the sketch curve as google style for the layed editing visualization  */
void CGlView::display_brush_sketch_google(GLenum mode)
{
	if(brushes == NULL)
		return;

	int i, j;

	//for(i=0;i<brushes->nbrushes;i++)
	//{
	//	if(mode==GL_SELECT)
	//		glLoadName(NAMEOFBRUSHES+i);

	//	glColor3f(0, 0, 0);

	//	//glLineWidth(7.);
	//	//glBegin(GL_LINE_STRIP);
	//	//for(j=0; j<brushes->brushlist[i].nelems; j++)
	//	//{
	//	//	glVertex2f(brushes->brushlist[i].brushpts[j]->x, 
	//	//		brushes->brushlist[i].brushpts[j]->y);
	//	//}
	//	//glEnd();

	//	double size=(VisMinRoadWidth+2)/1000;

	//	icVector2 line_dir, norm;
	//	double pre_up[2], pre_low[2];

	//	Brush *curbrush=&brushes->brushlist[i];

	//	line_dir.entry[0]=curbrush->brushpts[1]->x-curbrush->brushpts[0]->x;
	//	line_dir.entry[1]=curbrush->brushpts[1]->y-curbrush->brushpts[0]->y;
	//	normalize(line_dir);
	//	norm.entry[0]=-line_dir.entry[1];
	//	norm.entry[1]=line_dir.entry[0];
	//	norm = size/1.3* norm;

	//	pre_up[0]= curbrush->brushpts[0]->x+norm.entry[0];
	//	pre_up[1]= curbrush->brushpts[0]->y+norm.entry[1];
	//	pre_low[0]= curbrush->brushpts[0]->x-norm.entry[0];
	//	pre_low[1]= curbrush->brushpts[0]->y-norm.entry[1];

	//	for(j=1;j<brushes->brushlist[i].nelems;j++)
	//	{
	//		if(j<brushes->brushlist[i].nelems-1)
	//		{
	//			line_dir.entry[0]=curbrush->brushpts[j+1]->x-
	//				curbrush->brushpts[j-1]->x;
	//			line_dir.entry[1]=curbrush->brushpts[j+1]->y-
	//				curbrush->brushpts[j-1]->y;
	//		}
	//		else{
	//			line_dir.entry[0]=curbrush->brushpts[j]->x-
	//				curbrush->brushpts[j-1]->x;
	//			line_dir.entry[1]=curbrush->brushpts[j]->y-
	//				curbrush->brushpts[j-1]->y;
	//		}
	//		normalize(line_dir);
	//		norm.entry[0]=-line_dir.entry[1];
	//		norm.entry[1]=line_dir.entry[0];

	//		norm = size/1.3* norm;


	//		glBegin(GL_POLYGON);
	//		glVertex2f(pre_up[0], pre_up[1]);
	//		glVertex2f(curbrush->brushpts[j]->x+norm.entry[0],
	//			curbrush->brushpts[j]->y+norm.entry[1]);
	//		glVertex2f(curbrush->brushpts[j]->x-norm.entry[0],
	//			curbrush->brushpts[j]->y-norm.entry[1]);
	//		glVertex2f(pre_low[0], pre_low[1]);
	//		glEnd();

	//		pre_up[0]=curbrush->brushpts[j]->x+norm.entry[0];
	//		pre_up[1]=curbrush->brushpts[j]->y+norm.entry[1];
	//		pre_low[0]=curbrush->brushpts[j]->x-norm.entry[0];
	//		pre_low[1]=curbrush->brushpts[j]->y-norm.entry[1];
	//	}
	//	
	//	/*  display the end point */
	//	size=(VisMinRoadWidth+2)/1200.;

	//	DrawSolidCircle_size(curbrush->brushpts[0]->x, curbrush->brushpts[0]->y, size/zoom_factor);
	//	DrawSolidCircle_size(curbrush->brushpts[j-1]->x, curbrush->brushpts[j-1]->y, size/zoom_factor);



	//	size=(VisMinRoadWidth+2)/1400.;

	//	glColor3f(.95, .75, 0.14);
	//	//glLineWidth(5.);
	//	//glBegin(GL_LINE_STRIP);
	//	//for(j=0; j<brushes->brushlist[i].nelems; j++)
	//	//{
	//	//	glVertex2f(brushes->brushlist[i].brushpts[j]->x, 
	//	//		brushes->brushlist[i].brushpts[j]->y);
	//	//}
	//	//glEnd();

	//	line_dir.entry[0]=curbrush->brushpts[1]->x-curbrush->brushpts[0]->x;
	//	line_dir.entry[1]=curbrush->brushpts[1]->y-curbrush->brushpts[0]->y;
	//	normalize(line_dir);
	//	norm.entry[0]=-line_dir.entry[1];
	//	norm.entry[1]=line_dir.entry[0];
	//	norm = size/1.3* norm;

	//	pre_up[0]= curbrush->brushpts[0]->x+norm.entry[0];
	//	pre_up[1]= curbrush->brushpts[0]->y+norm.entry[1];
	//	pre_low[0]= curbrush->brushpts[0]->x-norm.entry[0];
	//	pre_low[1]= curbrush->brushpts[0]->y-norm.entry[1];

	//	for(j=1;j<brushes->brushlist[i].nelems;j++)
	//	{
	//		if(j<brushes->brushlist[i].nelems-1)
	//		{
	//			line_dir.entry[0]=curbrush->brushpts[j+1]->x-
	//				curbrush->brushpts[j-1]->x;
	//			line_dir.entry[1]=curbrush->brushpts[j+1]->y-
	//				curbrush->brushpts[j-1]->y;
	//		}
	//		else{
	//			line_dir.entry[0]=curbrush->brushpts[j]->x-
	//				curbrush->brushpts[j-1]->x;
	//			line_dir.entry[1]=curbrush->brushpts[j]->y-
	//				curbrush->brushpts[j-1]->y;
	//		}
	//		normalize(line_dir);
	//		norm.entry[0]=-line_dir.entry[1];
	//		norm.entry[1]=line_dir.entry[0];

	//		norm = size/1.3* norm;


	//		glBegin(GL_POLYGON);
	//		glVertex2f(pre_up[0], pre_up[1]);
	//		glVertex2f(curbrush->brushpts[j]->x+norm.entry[0],
	//			curbrush->brushpts[j]->y+norm.entry[1]);
	//		glVertex2f(curbrush->brushpts[j]->x-norm.entry[0],
	//			curbrush->brushpts[j]->y-norm.entry[1]);
	//		glVertex2f(pre_low[0], pre_low[1]);
	//		glEnd();

	//		pre_up[0]=curbrush->brushpts[j]->x+norm.entry[0];
	//		pre_up[1]=curbrush->brushpts[j]->y+norm.entry[1];
	//		pre_low[0]=curbrush->brushpts[j]->x-norm.entry[0];
	//		pre_low[1]=curbrush->brushpts[j]->y-norm.entry[1];
	//	}

	//	/*  display the end point */
	//	size=(VisMinRoadWidth+1)/1500.;

	//	DrawSolidCircle_size(curbrush->brushpts[0]->x, curbrush->brushpts[0]->y, size/zoom_factor);
	//	DrawSolidCircle_size(curbrush->brushpts[j-1]->x, curbrush->brushpts[j-1]->y, size/zoom_factor);
	//}


	if(sketchlines==NULL)
		return;

	for(i=0;i<sketchlines->ntrajs;i++)
	{
		//if(mode==GL_SELECT)
		//	glLoadName(NAMEOFBRUSHES+i);

		glColor3f(0, 0, 0);
		double size=(VisMinRoadWidth+2)/1000;

		size=size/zoom_factor;

		icVector2 line_dir, norm;
		double pre_up[2], pre_low[2];

		Trajectory *curtraj=sketchlines->trajs[i];

		LineSeg *curline=&curtraj->linesegs[0];

		line_dir.entry[0]=curline->gend[0]-curline->gstart[0];
		line_dir.entry[1]=curline->gend[1]-curline->gstart[1];
		normalize(line_dir);
		norm.entry[0]=-line_dir.entry[1];
		norm.entry[1]=line_dir.entry[0];
		norm = size/1.3* norm;

		pre_up[0]= curline->gstart[0]+norm.entry[0];
		pre_up[1]= curline->gstart[1]+norm.entry[1];
		pre_low[0]= curline->gstart[0]-norm.entry[0];
		pre_low[1]= curline->gstart[1]-norm.entry[1];

		for(j=1;j<curtraj->nlinesegs;j++)
		{
			if(j<curtraj->nlinesegs-1)
			{
				line_dir.entry[0]=curtraj->linesegs[j+1].gend[0]-
					curtraj->linesegs[j-1].gstart[0];
				line_dir.entry[1]=curtraj->linesegs[j+1].gend[1]-
					curtraj->linesegs[j-1].gstart[1];
			}
			else{
				line_dir.entry[0]=curtraj->linesegs[j].gend[0]-
					curtraj->linesegs[j-1].gstart[0];
				line_dir.entry[1]=curtraj->linesegs[j].gend[1]-
					curtraj->linesegs[j-1].gstart[1];
			}
			normalize(line_dir);
			norm.entry[0]=-line_dir.entry[1];
			norm.entry[1]=line_dir.entry[0];

			norm = size/1.3* norm;


			glBegin(GL_POLYGON);
			glVertex2f(pre_up[0], pre_up[1]);
			glVertex2f(curtraj->linesegs[j].gstart[0]+norm.entry[0],
				curtraj->linesegs[j].gstart[1]+norm.entry[1]);
			glVertex2f(curtraj->linesegs[j].gstart[0]-norm.entry[0],
				curtraj->linesegs[j].gstart[1]-norm.entry[1]);
			glVertex2f(pre_low[0], pre_low[1]);
			glEnd();

			pre_up[0]=curtraj->linesegs[j].gstart[0]+norm.entry[0];
			pre_up[1]=curtraj->linesegs[j].gstart[1]+norm.entry[1];
			pre_low[0]=curtraj->linesegs[j].gstart[0]-norm.entry[0];
			pre_low[1]=curtraj->linesegs[j].gstart[1]-norm.entry[1];
		}
		/*   consider last line segment  */
		curline=&curtraj->linesegs[curtraj->nlinesegs-1];
		glBegin(GL_POLYGON);
		glVertex2f(pre_up[0], pre_up[1]);
		glVertex2f(curline->gend[0]+norm.entry[0],
			curline->gend[1]+norm.entry[1]);
		glVertex2f(curline->gend[0]-norm.entry[0],
			curline->gend[1]-norm.entry[1]);
		glVertex2f(pre_low[0], pre_low[1]);
		glEnd();
	
		/*  display the end point */
		size=(VisMinRoadWidth+2)/1200.;

		curline=&curtraj->linesegs[0];
		DrawSolidCircle_size(curline->gstart[0], curline->gstart[1], size/zoom_factor);
		curline=&curtraj->linesegs[curtraj->nlinesegs-1];
		DrawSolidCircle_size(curline->gend[0], curline->gend[1], size/zoom_factor);


		/*   visualize front  */

		size=(VisMinRoadWidth+2)/1400.;
		size=size/zoom_factor;

		glColor3f(.95, .75, 0.14);
		curline=&curtraj->linesegs[0];
		line_dir.entry[0]=curline->gend[0]-curline->gstart[0];
		line_dir.entry[1]=curline->gend[1]-curline->gstart[1];
		normalize(line_dir);
		norm.entry[0]=-line_dir.entry[1];
		norm.entry[1]=line_dir.entry[0];
		norm = size/1.3* norm;

		pre_up[0]= curline->gstart[0]+norm.entry[0];
		pre_up[1]= curline->gstart[1]+norm.entry[1];
		pre_low[0]= curline->gstart[0]-norm.entry[0];
		pre_low[1]= curline->gstart[1]-norm.entry[1];

		for(j=1;j<curtraj->nlinesegs;j++)
		{
			if(j<curtraj->nlinesegs-1)
			{
				line_dir.entry[0]=curtraj->linesegs[j+1].gend[0]-
					curtraj->linesegs[j-1].gstart[0];
				line_dir.entry[1]=curtraj->linesegs[j+1].gend[1]-
					curtraj->linesegs[j-1].gstart[1];
			}
			else{
				line_dir.entry[0]=curtraj->linesegs[j].gend[0]-
					curtraj->linesegs[j-1].gstart[0];
				line_dir.entry[1]=curtraj->linesegs[j].gend[1]-
					curtraj->linesegs[j-1].gstart[1];
			}
			normalize(line_dir);
			norm.entry[0]=-line_dir.entry[1];
			norm.entry[1]=line_dir.entry[0];

			norm = size/1.3* norm;


			glBegin(GL_POLYGON);
			glVertex2f(pre_up[0], pre_up[1]);
			glVertex2f(curtraj->linesegs[j].gstart[0]+norm.entry[0],
				curtraj->linesegs[j].gstart[1]+norm.entry[1]);
			glVertex2f(curtraj->linesegs[j].gstart[0]-norm.entry[0],
				curtraj->linesegs[j].gstart[1]-norm.entry[1]);
			glVertex2f(pre_low[0], pre_low[1]);
			glEnd();

			pre_up[0]=curtraj->linesegs[j].gstart[0]+norm.entry[0];
			pre_up[1]=curtraj->linesegs[j].gstart[1]+norm.entry[1];
			pre_low[0]=curtraj->linesegs[j].gstart[0]-norm.entry[0];
			pre_low[1]=curtraj->linesegs[j].gstart[1]-norm.entry[1];
		}
		/*   consider last line segment  */
		curline=&curtraj->linesegs[curtraj->nlinesegs-1];
		glBegin(GL_POLYGON);
		glVertex2f(pre_up[0], pre_up[1]);
		glVertex2f(curline->gend[0]+norm.entry[0],
			curline->gend[1]+norm.entry[1]);
		glVertex2f(curline->gend[0]-norm.entry[0],
			curline->gend[1]-norm.entry[1]);
		glVertex2f(pre_low[0], pre_low[1]);
		glEnd();


		/*  display the end point */
		size=(VisMinRoadWidth+2)/1500.;

		curline=&curtraj->linesegs[0];
		DrawSolidCircle_size(curline->gstart[0], curline->gstart[1], size/zoom_factor);
		curline=&curtraj->linesegs[curtraj->nlinesegs-1];
		DrawSolidCircle_size(curline->gend[0], curline->gend[1], size/zoom_factor);
		
	}


	///* display the obtain graph if not empty*/
	if(sharedvars.ShowSegmentGraphOn)
	{
		if(sketchnet==NULL||sketchnet->nodelist==NULL||sketchnet->nodelist->nelems==0) return;
		Intersection *curin;
		StreetGraphEdge *edge;

		///*display the edges first*/
		glLineWidth(1.);
		glColor3f(0, 0, 0);
		for(i=0; i<sketchnet->edgelist->nedges; i++)
		{
			edge = sketchnet->edgelist->edges[i];
			glBegin(GL_LINES);
				curin = sketchnet->nodelist->intersects[edge->node_index1];
				glVertex2f(curin->gpos[0], curin->gpos[1]);
				curin = sketchnet->nodelist->intersects[edge->node_index2];
				glVertex2f(curin->gpos[0], curin->gpos[1]);
			glEnd();
		}

		/*display the intersections*/
		for(i=0; i<sketchnet->nodelist->nelems; i++)
		{
			curin=sketchnet->nodelist->intersects[i];
			if(!curin->endpt)
			{
				if(mode == GL_SELECT)
				{
					glLoadName(NAMEOFINTERSECTS + i);////Assign name for intersection
				}
				glColor3f(1, 1, 1);
				DrawSolidCircle_size(curin->gpos[0], curin->gpos[1], 0.0045/zoom_factor);
				glColor3f(0, 0, 0);
				draw_hollow_circle_size(curin->gpos[0], curin->gpos[1], 0.005/zoom_factor);
			}
			else if(sharedvars.ShowStreetGraphEndPointsOn)
			{
				if(mode == GL_SELECT)
				{
					glLoadName(NAMEOFINTERSECTS + i);////Assign name for intersection
				}
				glColor3f(0, 1, 0.6);
				DrawSolidRect_size(curin->gpos[0], curin->gpos[1], 0.003/zoom_factor);
				glColor3f(0, 0, 0);
				draw_hollow_rect_size(curin->gpos[0], curin->gpos[1], 0.0035/zoom_factor);
			}
		}
	}

}

/*
Display the design grid to give user a feeling about how dense the tensor line could be
*/
void CGlView::display_design_grid()
{
	double design_yinterval=(majorDensity*quadmesh->yinterval);
	double design_xinterval=(minorDensity*quadmesh->xinterval);
	int i;
	int nlines_ydir=1./design_yinterval+1;
	int nlines_xdir=1./design_xinterval+1;

	glColor3f(0.8, 0.8, 0.5);
	glLineWidth(1.);
	//glBegin(GL_LINES);
	//for(i=0;i<nlines_ydir;i++)
	//{
	//	if(i*design_yinterval<1.1)
	//	glVertex2f(0, i*design_yinterval);
	//	glVertex2f(1, i*design_yinterval);
	//}
	//glEnd();

	//glBegin(GL_LINES);
	//for(i=0;i<nlines_xdir;i++)
	//{
	//	if(i*design_xinterval<1.1)
	//	glVertex2f(i*design_xinterval, 0);
	//	glVertex2f(i*design_xinterval, 1);
	//}
	//glEnd();

	QuadCell* face;
	QuadVertex *v;
	int j;
	for(i=0;i<quadmesh->nfaces;i++)
	{
		face=quadmesh->quadcells[i];

		glBegin(GL_LINE_LOOP);
		for(j=0;j<face->nverts;j++)
		{
			v=quadmesh->quad_verts[face->verts[j]];
			glVertex2f(v->x, v->y);
		}
		glEnd();
	}

	//glColor3f(1,0,0);
	//face=quadmesh->quadcells[668];
	//	glBegin(GL_LINE_LOOP);
	//	for(j=0;j<face->nverts;j++)
	//	{
	//		v=quadmesh->quad_verts[face->verts[j]];
	//		glVertex2f(v->x, v->y);
	//	}
	//	glEnd();
	//face=quadmesh->quadcells[1997];
	//	glBegin(GL_LINE_LOOP);
	//	for(j=0;j<face->nverts;j++)
	//	{
	//		v=quadmesh->quad_verts[face->verts[j]];
	//		glVertex2f(v->x, v->y);
	//	}
	//	glEnd();
	//face=quadmesh->quadcells[5872];
	//	glBegin(GL_LINE_LOOP);
	//	for(j=0;j<face->nverts;j++)
	//	{
	//		v=quadmesh->quad_verts[face->verts[j]];
	//		glVertex2f(v->x, v->y);
	//	}
	//	glEnd();
}
