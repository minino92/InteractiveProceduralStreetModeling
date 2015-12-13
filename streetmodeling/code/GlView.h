#pragma once

#include <gl/glut.h> 
#include "lib/icVector.h"


// CGlView

class CGlView : public CWnd
{
	DECLARE_DYNAMIC(CGlView)

// Construction
public:
	CGlView(CWnd *pclWnd);

	CGlView();
	virtual ~CGlView();

// Attributes
public:

	HDC  m_hDC;		        // GDI Device Context 
    HGLRC	m_hglRC;		// Rendering Context

    CWnd *m_pclWnd;
    HWND m_hWnd;

	/***----------------------------------------------------------***/
	////member variables for visual control
	static int MoveOrStop;
	int EditModeOn;
	int MoveElemOn;
	int SingularitiesOn;
	int RegularElemOn;
	int TrajectoryOn;
	int SeparatricesOn;
	int LimitCycleOn;

	/***----------------------------------------------------------***/
	int SmoothOn;                           ////flag for region smooth
	int DisplaySmoothRegionOn;
	int PickPointOn;                        ////flag for region pick up
	
	/***----------------------------------------------------------***/
	////Singularity pair cancellation
	int PairCancelOn;                       ////for singularities pair cancelation
	int sing1, sing2;                       ////Two singularities for pair cancelation
	int pair_counter;

	int SelectRepellerOrAttractor;
	
	/***----------------------------------------------------------***/
    ////Singularity movement
	int SingularityMoveOn;
	int source_triangle, target_triangle, oldsingularityID;
	int triangle_counter;
	double newposx, newposy;

	/***----------------------------------------------------------***/
	////limit cycle + singularity pair cancellation
	int LimitSingCancelOn;
	int pairsing, pairlimit;
	int limitsing_counter;
    
	////limit cycle pair cancellation
	int LimitPairCancelOn;
	int pairlimit_1, pairlimit_2;
	int pairlimit_counter;
	
	
	/***----------------------------------------------------------***/
	////For limit cycle relocate
	int LimitCycleRelocateOn;
	int LimitCycleDeformOn;


	////Member variables and flags for shape editing and control
	int NewCurveOn;
	int ShapeControlPtsOn;
	int ShapeEditOn;
	int which_shapectrpt;
	int FinisheCtrptsel;

	////Separatrix editing
	int SeparatrixEditOn;
	
	/***----------------------------------------------------------***/

	////Streamline based method
	int StreamlineBasedOn;  //2/18/06

	/***----------------------------------------------------------***/
	////Member variables for mouse selection
	int TraceBeginTriangleID;

	int TransformType;                     ////1 scale, 2 rotation
	int which_control_point;               ////which control point has been selected

	/***----------------------------------------------------------***/
	int newelementtriangle;

	////memeber variables for beginning point of a trajectory
	double bx, by;

	////member variables for current position of a particle
	double px, py;
	int particle_triangle;
	int ReleaseParticleOn;

	float hsv[3], rgb[3];
	int ColorPlotOn;

	/**Test finding SCC**/
	int ShowSCCOn;

	/** Test calculating eigen vectors for given triangle **/
	int ShowEigenVecOn;

	int ShowSamplePtsOn;

	/////////////////////////
	////Testing variables
	int ShowCancelSmoothRegion;



// Operations for OpenGl window setting and displaying
public:
	BOOL SetPixelformat(HDC hdc);  //set the pixel format for the OpenGL window
	int InitGL(GLvoid);	           //initialize the OpenGl envrionment
	int DrawGLScene(GLenum mode);       //The same function as display routine to show the visual effect in the opengl window
	static void IBFVEffect(GLenum mode);
	int OnCreate();
	afx_msg void OnDestroy();
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);


protected:
	DECLARE_MESSAGE_MAP()

/*--------------------------------------------------------------------------*/
//Operations for ibfv tool visual effect displaying
public:
    void makePatterns(void);  ////To initialize the texture pattern for flow effect visualization
	void getInitTex(void);

	void InitVFVariables();   ////initialize the vector field related global variables here
	void InitFlag();          ////initialize some flags here

	void finalize();          ////release the memory

//Operations for visual icons and trajectories
	void DisplayEditBox(GLenum mode);          ////display the edit box for specified elememt
	void DisplaySingularControlPoints(GLenum mode, int SelectedElem); ////we just need to display 9 points
	void DisplayRegularControlPoints(GLenum mode, int SelectedElem); ////we only need to display 2 points
	void DisplayElem(GLenum mode);                ////display the center legends for elements
	void DisplayRegularElemt(GLenum mode);        ////display the regular elements using arrows
	
	void DisplayCapturedSin(GLenum mode);         ////display the captured singularities according to current mode
    void DisplayAllCapturedSin(GLenum mode);       ////display all the captured singularities

	void DisplayLimitCycleLegends(GLenum mode);

	void DisplayTrajectory();
	void DisplaySeparatrices();
	void DisplaySmoothRegion();

	void DisplayShapeControlPts(GLenum mode);                      ////Display the control point for shape control
	void DisplayDesignCurve();

	bool IsInCenter(double x, double y, int &singular_id);

	void SetColorByType(int);
	void DrawMarkTriangle(double cx, double cy);        ////draw a triangular legend
	void DrawSolidCircle(double cx, double cy);         ////draw a solid circle legend
	void draw_hollow_circle(double cx, double cy);

	void DrawControlRectangle(double cx, double cy);    ////draw a rectangular legend
    void DrawUnitArrow();                               ////draw an arrow legend
	void DrawArrow(double base[2], double Direc[2]);
    void DrawEditBox(icVector2 p1, icVector2 p2, icVector2 p3, icVector2 p4);
	void DrawEditBox_back(icVector2 p1, icVector2 p2, icVector2 p3, icVector2 p4);
	void DrawSingleTrajectory(int);
	void DrawSeparatrices();
	void DisplayEigenVectorForSaddle();
	//void DrawAGroupSeparatrix(int SeparatrixID);
	void DrawLimitCycleLegend(double cx, double cy, double bx, double by, int type);

	void DrawEigenVector(icVector2 vec, double cx, double cy);
	
	void BuildHighLighted();
	void DrawHighLighted();

	/**-------------------------------------------------------------------------**/
	////For multiple repellers and attractors selection  11/20/05
	void AddToRepellerList(int, int);
	void AddToAttractorList(int, int);
	void ClearRepellandAttractList();

	/*---------------------------05/30/05----------------------------------------*/
	////Display cell cycle detection 
	void DisplayAllCellCycle();
	void DisplayLimitCycles();

	/**Show SCC components having more than 2 triangles **/
	void DisplaySCC();

	/**Show the eigen vectors of the triangles in the SCC components having more than 2 triangles **/
	void DisplayEigenVecForSCC();

	void DisplaySamplingPts();

    //////functions for mouse screen selection
	void HitProcess(double ss, double st);
    void HitProcessforGraph(double ss, double st);
	void HitProcessforSelectUnderneathMesh(double ss, double st);
	void HitProcessforLimitDeform(double ss, double st);
	void HitProcessforLimitCycleSelect(double ss, double st);

	void ControlPointForSinElem(int );
	void ControlPointForRegElem(int );


	//////Testing functions
	void DisplayElementTriangles();
	void DisplaySingularTriangles();

	void DisplayWholeMeshWithFieldArrows();

	////Display the boundary of the new triangle strip 12/27/05
	void DisplayNewTriangleStrip();

	void DisplayParticle();
	void DrawAParticle(double cx, double cy);

	void DisplayColorPlots();
    void HsvRgb(float hsv[3], float rgb[3]);

	////For testing the cell cycle containing the design curve
	void DisplayOldDesignCellCycle();

	/* 01/15/07 Test the being covered triangles */
	int ShowCoveredTris;
	void DisplayCoveredTris();

	/* 02/12/07 */
	/** Show highly curled region */
	int ShowHighCurl;
	void DisplayHighlyCurl();

	int ShowHighDiv;
	void DisplayHighlyDiv();

	/**Display the path between two triangle**/
	int ShowTestPathOn;
	void DisplayTestPath();
	void DisplayAllEdges(int);
	int ShowAllEdgesOn;

	/*to display the P component for Jim's data 04/10/07*/
	void DisplayMag();
	int ShowMagOn;

	/*Display the decomposition of the regions using asymmetry method 04/12/07*/
	void DisplayDecomp_regions();
	int ShowDecompOn;

	/*just for fun*/
	void VF_to_ColorMap(GLenum mode);

	/*display the sample points using different colors*/
	int showTraceSamplingPts;
	void display_trace_SampltPts(int);

	int showSampDensityOn;
	void display_sample_density();

	int showEdgeImageOn;
	void display_image_edge();
	void display_ver_time();


	int showMCGConnectionOn;
	void display_MCG_connections();
	void DrawSolidCircle_size(double cx, double cy, double R);
	void DrawSolidRect_size(double cx, double cy, double R);
	void draw_hollow_circle_size(double cx, double cy, double R);
	void draw_hollow_rect_size(double cx, double cy, double R);
	void draw_rect_size(double cx, double cy, double R);

	int showTensorOn;
	void display_degenerate_tris();
	void display_degenerate_pts(GLenum mode);
	void display_seps_directions();
	void display_singularElem_controlPts_tensor(GLenum mode, int SelectedElem);

	bool is_ten_designElem(double x, double y, int &degptindex);
	void display_tenElem_EditBox(GLenum mode);
	void ten_hitProcess(double ss, double st);
	void display_tenRegElem(GLenum mode);
	void display_degenerate_cells();
	
	void display_major_tenlines(GLenum mode);
	void display_minor_tenlines(GLenum mode);

	void display_level1_tenlines(GLenum mode);
	void display_init_seeds(GLenum mode);

	/*   Testing the sample points 12/04/2007   */
	void display_major_samps();
	void display_minor_samps();

	void display_intersections(GLenum mode);
	void display_streetnet(GLenum mode);

	void display_majRoadnet(GLenum mode);

	bool displayRoadNetOn;
	void display_roads(GLenum mode);

	void display_brush();

	void display_road_width_back(int id, double width, bool majormin);
	void display_road_width_front(int id, double width, bool majormin);
	void display_roads_width(GLenum mode);
	void display_road_graph_back(double width, unsigned char roadtype);
	void display_road_graph_front(double width, unsigned char roadtype);

	void display_road_use_network();
	void display_road_graph_linestyle();


	int displayDisMapOn;

	bool streetNetEditOn;

	bool showTensorLineOn;
	bool showStreetGraphOn;

	bool deleteElemOn;

	void display_design_grid();

	void display_brush_sketch_thin(GLenum mode);
	void display_brush_sketch_google(GLenum mode);

	void display_majRoad_googlestyle(GLenum mode);
	void display_majRoad_googlestyle_back_temp(GLenum mode, bool majormin);
	void display_majRoad_googlestyle_front_temp(GLenum mode, bool majormin);

	/*  for anti-aliasing visualization  */
	int with_anti_aliasing(GLenum mode);
	int without_anti_aliasing(GLenum mode);


};


void cal_inverse_transform();
void transform_point3D(float p[3], float rot_mat[16]);
