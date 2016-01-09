// IBFVDlg.cpp : implementation file
//

#include "stdafx.h"
#include "IBFV2D.h"
#include "IBFVDlg.h"
#include ".\ibfvdlg.h"


#include "VFSynthesis.h"
#include "LimitCycleCreator.h"
#include "VFAnalysis.h"
#include "LocalTracing.h"
#include "LimitCycleDetect.h"
#include "RegionSmoothing.h"
#include "topologyedit.h"
#include "LimitCycleEdit.h"
#include "ConleyRelationGraph.h"
#include "SeparatrixEdit.h"
#include "NewLimitCycleDetect.h"

#include "BmpProcess.h"
#include "FindSCC.h"
#include "FindSepAttPoints.h"

#include "SCCCycleDetect.h"

#include "shareinterfacevars.h"

#include "ImgBoundaryExtract.h"
#include "BoundaryBasedTenGen.h"

#include "SketchDesign.h"

#include "loadmaps.h"

#include "scalardesign.h"

#include "SaveProjectSetting.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


/*Here we define a global variable to encapsulate the interface variables 
that may be shared by multiple dialogs 11/08/2007*/

SharedInterfaceVars sharedvars;

void init_sharedvars()
{
	sharedvars.BrushInterfaceOn=false;
	sharedvars.CombineWithOtherPartsOn=true;
	sharedvars.EditElementOn=false;
	sharedvars.MoveElemOn=false;
	sharedvars.RemoveElemOn=false;
	sharedvars.ShowExtractedBoundsOn=true;
	sharedvars.ShowRegElemOn=true;
	sharedvars.ShowRoadMapOn=false;
	sharedvars.ShowSingularitiesOn=true;
	sharedvars.ShowStreetGraphOn=false;
	sharedvars.ShowTensorLinesOn=false;
	sharedvars.ShowTheMapOn=false;
	sharedvars.TensorDesignOn=false;
	sharedvars.MeanFilterOn=true;
	sharedvars.GenTenFromExtractedBoundsOn=true;
	sharedvars.UseAllBoundsOn=true;
	sharedvars.DesignGridOn=false;
	sharedvars.ConnectDeadEndsOn=false;
	sharedvars.ShowRegionBlocksOn=false;
	sharedvars.ShowStreetGraphEndPointsOn=true;
	sharedvars.EnableSketchBasedDesign=false;
	//sharedvars.CombinePopDensityOn=false;
	sharedvars.CombinePopDensityOn=true;  /*  modify on 1/19/2008*/
	sharedvars.GenTenFromSketchesOn=true;
	sharedvars.ShowSketchesOn=false;
	sharedvars.CloseLoopOn=true;
	sharedvars.ShowPopDensityMapOn=false;
	sharedvars.ShowIBFVOn=true;
	sharedvars.EditStreetNetOn=false;
	sharedvars.ShowIntersectsOn=true;
	sharedvars.ShowStreetUseNetworkOn=false;
	sharedvars.SelStreetRegToEditOn=false;
	sharedvars.ShowLineStyleStreetsOn=false;
	sharedvars.UseBoundsAsSketchesOn=true;
	sharedvars.UseBoundsAsRoadsOn=false;
	sharedvars.UseMajRoadsAsSketchesOn=false;
	sharedvars.ShowScalarFieldOn=false;
	sharedvars.EnableHeighfieldDesignOn=false;
	sharedvars.ShowMajRoadsOn=false;
	sharedvars.ShowInitSeedsOn=false;
	sharedvars.ShowSegmentGraphOn=false;
	sharedvars.ShowVegMapOn=false;
	sharedvars.ApplyAsymFldOn=false;
	sharedvars.RemoveMajDeadEndsOn=true;
	
	sharedvars.ShowMajRoadGoogleStyleOn=false;

	sharedvars.AntiAliasingOn=false;
	sharedvars.BrushLikeRegSelOn=false;


	/*  in the major road setting dialog  */
	sharedvars.AllowMajRoadCrossRiverOn=false;
	sharedvars.AllowMajRoadFollowBoundaryOn=false;
	sharedvars.ShowMajRoadNetworkOn=false;
	sharedvars.AllowCrossSingularitiesOn=false;
	sharedvars.JobardMethodOn=false;

	/*  in the minor road setting dialog  */
	sharedvars.AllowMinCloseToMajOn=true;
	sharedvars.RemDeadEndsTraceOn=false;
	sharedvars.UseNewMedRemDeadEndOn=false;

	/*  in the save project setting dialog  */
	sharedvars.rdSaveProjDesignElem=0;
	sharedvars.rdSaveProjSketches=0;
	sharedvars.rdSaveProjBrushes=0;
	sharedvars.rdSaveProjMajRoadNetwork=0;
	sharedvars.rdSaveProjOtherSetting=0;
	sharedvars.rdSaveProjStreetNetwork=0;

	/* for the radio buttons*/
	sharedvars.rdSketchMajorHighways=false;
	sharedvars.rdScalarSingularElem=false;
}

extern Polygon3D Object;
extern QuadMesh *quadmesh;


CGlView *g_pclGlView = NULL;

MyTabCtrl *g_m_tbCtrl = NULL;

int cur_selectTris = 0;

extern int NDesignRegions;
extern int cur_chosen_region; 

/**/
extern double zoom_factor ;
extern double trans_x;
extern double trans_y;
extern double street_sample_interval;
extern double sx, sy, uniforms;
extern double rotateAng;
extern float inverse_tran[16];


/*for tensor field vis/analysis/design*/
#include "tensorvis.h"
#include "caldeformation.h"
#include "tensordesign.h"
#include "tensoranalysis.h"
#include "computeroadvis.h"
#include "evenlystreamlines.h"
#include "regionsmooth_quad.h"

extern void detect_degeneratepts();
extern void detect_degeneratepts_tranvec();

extern ctr_point *control_pts;        // allocate our control point array
extern int num_shapecontrol_pts;
extern int num_curvepts_output;
extern int resolution;


extern int chosen_tenelem_ID;
extern Degenerate_Design *ten_designelems;
extern TenRegularElem *ten_regularelems ;
extern int nten_regelems ;
extern StreetVis *roadnetvis ;
extern EvenStreamlinePlace *major, *minor;
extern bool brushinterfaceOn;
extern bool closedbrush;
extern int *boundarycells;
extern int nboundarycells;
extern StreetNet *streetnet;
extern int selectedIntersect;

extern ctr_point **brushpts;
extern int nbrushpts;
extern int curMaxNumBrushPts;

int REALWINSIZE=800;

//#include "ppm.h"
extern bool flag_loadmap;
extern unsigned char *map1, *fittedmap1, *displaymap, *streetmapbackground;
extern int map1_w, map1_h;
extern int boundindex1, boundindex2;
extern bool ylargerthanx;
extern bool flag_loadvegmap;

/*       vegetation map      */
extern unsigned char *vegmap_fit;
extern unsigned char *vegmap_disp;

/*       Population density map      */
extern unsigned char *popdensitymap_fit;
extern unsigned char *popdensitymap_disp;

/*       Hieght field       */
extern unsigned char *heightfield_fit;
extern unsigned char *heightfield_dis;

extern int MeanFilterSize;

//extern bool *map_boundary_mask;


ImgProcess *imgprocess=NULL;
ImgBoundary *imgboundaries=NULL;
extern MapBoundaryList *mapboundarylist;
extern SeedList *seedsalongbounds;

extern void obtain_smooth_region_multibounds(double widtheachbound);
extern void obtain_field_basis();

extern double BoundRegionWidth;

extern void add_to_current_brushPts(double, double);

extern SketchList *brushes;

/*   for globally updating the number of design brushes    */
extern int *g_m_edNumBrushes;

extern bool highwayexisted;
extern bool majorroadsexisted;

extern int gen_regElem_sketches;


#include "MyDlg1.h"
#include "MyDlg2.h"
#include "StreetNetPanel.h"

extern MyDlg1 *g_mydlg1;
extern MyDlg2 *g_mydlg2;
extern StreetNetPanel *g_streetnetpanel;


extern icMatrix2x2 *pre_tenfield;

extern bool please_comb_prefield;

double map_xrang = 10000;  // default is 10km
double map_yrang = 10000;  // default is 10km

extern bool upstreaming_edit;
extern icMatrix2x2 *pre_ten;

//extern int cur_sel_region;
extern bool is_on_local_editing;



/*
Keep track of the mouse movement to obtain the brush path
*/
void add_to_brushEdPts(double x, double y)
{
	int i;
	if(nbrushpts >= curMaxNumBrushPts)
	{
		ctr_point **temp=brushpts;
		brushpts=(ctr_point**)malloc(sizeof(ctr_point*)*(curMaxNumBrushPts+100));
		if(brushpts == NULL)
			exit(-1);
		for(i=0; i<curMaxNumBrushPts; i++)
			brushpts[i]=temp[i];
		for(i=curMaxNumBrushPts; i<=curMaxNumBrushPts+100; i++)
			brushpts[i]=(ctr_point*)malloc(sizeof(ctr_point));
		curMaxNumBrushPts += 100;
	}

	brushpts[nbrushpts]->x=x;
	brushpts[nbrushpts]->y=y;
	nbrushpts++;
}


extern int MaxNumPoints;
extern Point *point;                       ////we may initial it as 50 points, over 50, we can extend it
extern int Num_SmoothRegionpoints;                     ////Number of points that user selected
extern ctr_point *control_pts;        // allocate our control point array
extern int num_shapecontrol_pts;

void compute_brush_Reg()
{
	/**/
	if(nbrushpts <= 2)
		return;

	int i;
	icVector2 brushDir, norm;

	Num_SmoothRegionpoints=0;

	for(i=0;i<num_shapecontrol_pts-1;i++)
	{
		brushDir.entry[0]=control_pts[i+1].x-control_pts[i].x;
		brushDir.entry[1]=control_pts[i+1].y-control_pts[i].y;
		normalize(brushDir);
		norm.entry[0]=-brushDir.entry[1];
		norm.entry[1]=brushDir.entry[0];

		/* */
		if(Num_SmoothRegionpoints>=MaxNumPoints)
		{
			MaxNumPoints += 50;
			point = (Point*)realloc(point, sizeof(Point) * MaxNumPoints);
		}

		if(i==0)
		{
			point[Num_SmoothRegionpoints].x=
				brushpts[i]->x+BoundRegionWidth*quadmesh->xinterval*norm.entry[0];
			point[Num_SmoothRegionpoints].y=
				brushpts[i]->y+BoundRegionWidth*quadmesh->xinterval*norm.entry[1];
		}
		else
		{
			point[Num_SmoothRegionpoints].x=
				control_pts[i].x+BoundRegionWidth*quadmesh->xinterval*norm.entry[0];
			point[Num_SmoothRegionpoints].y=
				control_pts[i].y+BoundRegionWidth*quadmesh->xinterval*norm.entry[1];
		}

		Num_SmoothRegionpoints++;
	}

	/* add the last point */
	point[Num_SmoothRegionpoints].x=
		control_pts[num_shapecontrol_pts-1].x+BoundRegionWidth*quadmesh->xinterval*norm.entry[0];
	point[Num_SmoothRegionpoints].y=
		control_pts[num_shapecontrol_pts-1].y+BoundRegionWidth*quadmesh->xinterval*norm.entry[1];
	Num_SmoothRegionpoints++;

	/*inversed*/
	/* add the first point (i.e. the last point in previous iteration)  */
	point[Num_SmoothRegionpoints].x=
		control_pts[num_shapecontrol_pts-1].x-BoundRegionWidth*quadmesh->xinterval*norm.entry[0];
	point[Num_SmoothRegionpoints].y=
		control_pts[num_shapecontrol_pts-1].y-BoundRegionWidth*quadmesh->xinterval*norm.entry[1];
	Num_SmoothRegionpoints++;

	for(i=num_shapecontrol_pts-2;i>=0;i--)
	{
		brushDir.entry[0]=control_pts[i+1].x-control_pts[i].x;
		brushDir.entry[1]=control_pts[i+1].y-control_pts[i].y;

		normalize(brushDir);
		norm.entry[0]=brushDir.entry[1];
		norm.entry[1]=-brushDir.entry[0];
		
		if(Num_SmoothRegionpoints>=MaxNumPoints)
		{
			MaxNumPoints += 50;
			point = (Point*)realloc(point, sizeof(Point) * MaxNumPoints);
		}
		if(i==0)
		{
			point[Num_SmoothRegionpoints].x=
				brushpts[i]->x+BoundRegionWidth*quadmesh->xinterval*norm.entry[0];
			point[Num_SmoothRegionpoints].y=
				brushpts[i]->y+BoundRegionWidth*quadmesh->xinterval*norm.entry[1];
		}
		else
		{
			point[Num_SmoothRegionpoints].x=
				control_pts[i].x+BoundRegionWidth*quadmesh->xinterval*norm.entry[0];
			point[Num_SmoothRegionpoints].y=
				control_pts[i].y+BoundRegionWidth*quadmesh->xinterval*norm.entry[1];
		}

		Num_SmoothRegionpoints++;
	}
}


//void compute_brush_Reg()
//{
//	/**/
//	if(nbrushpts <= 2)
//		return;
//
//	int i;
//	icVector2 brushDir, norm;
//
//	Num_SmoothRegionpoints=0;
//
//	double pre_p[2]={brushpts[0]->x, brushpts[0]->y};
//	for(i=0;i<nbrushpts-1;i++)
//	{
//		if(i%10!=0)
//			continue;
//
//		if(i==0)
//		{
//			brushDir.entry[0]=brushpts[i+1]->x-brushpts[i]->x;
//			brushDir.entry[1]=brushpts[i+1]->y-brushpts[i]->y;
//		}
//		else
//		{
//			brushDir.entry[0]=brushpts[i]->x-pre_p[0];
//			brushDir.entry[1]=brushpts[i]->y-pre_p[1];
//		}
//		pre_p[0]=brushpts[i]->x;
//		pre_p[1]=brushpts[i]->y;
//
//		normalize(brushDir);
//		norm.entry[0]=-brushDir.entry[1];
//		norm.entry[1]=brushDir.entry[0];
//
//		/* */
//		if(Num_SmoothRegionpoints>=MaxNumPoints)
//		{
//			MaxNumPoints += 50;
//			point = (Point*)realloc(point, sizeof(Point) * MaxNumPoints);
//		}
//
//		point[Num_SmoothRegionpoints].x=
//			brushpts[i]->x+BoundRegionWidth*quadmesh->xinterval*norm.entry[0];
//		point[Num_SmoothRegionpoints].y=
//			brushpts[i]->y+BoundRegionWidth*quadmesh->xinterval*norm.entry[1];
//		point[Num_SmoothRegionpoints].cellid=get_cellID_givencoords(point[Num_SmoothRegionpoints].x,
//			point[Num_SmoothRegionpoints].y);
//		Num_SmoothRegionpoints++;
//	}
//
//	/* add the last point */
//	point[Num_SmoothRegionpoints].x=
//		brushpts[nbrushpts-1]->x+BoundRegionWidth*quadmesh->xinterval*norm.entry[0];
//	point[Num_SmoothRegionpoints].y=
//		brushpts[nbrushpts-1]->y+BoundRegionWidth*quadmesh->xinterval*norm.entry[1];
//	Num_SmoothRegionpoints++;
//	point[Num_SmoothRegionpoints].cellid=get_cellID_givencoords(point[Num_SmoothRegionpoints].x,
//		point[Num_SmoothRegionpoints].y);
//
//	/*inversed*/
//	/* add the first point (i.e. the last point in previous iteration)  */
//	point[Num_SmoothRegionpoints].x=
//		brushpts[nbrushpts-1]->x-BoundRegionWidth*quadmesh->xinterval*norm.entry[0];
//	point[Num_SmoothRegionpoints].y=
//		brushpts[nbrushpts-1]->y-BoundRegionWidth*quadmesh->xinterval*norm.entry[1];
//	Num_SmoothRegionpoints++;
//	point[Num_SmoothRegionpoints].cellid=get_cellID_givencoords(point[Num_SmoothRegionpoints].x,
//		point[Num_SmoothRegionpoints].y);
//	pre_p[0]=brushpts[nbrushpts-1]->x;
//	pre_p[1]=brushpts[nbrushpts-1]->y;
//
//	for(i=nbrushpts-2;i>=0;i--)
//	{
//		if(i%10!=0)
//			continue;
//
//		if(i==0)
//		{
//			brushDir.entry[0]=brushpts[i]->x-brushpts[i+1]->x;
//			brushDir.entry[1]=brushpts[i]->y-brushpts[i+1]->y;
//		}
//		else
//		{
//			brushDir.entry[0]=brushpts[i]->x-pre_p[0];
//			brushDir.entry[1]=brushpts[i]->y-pre_p[1];
//		}
//		pre_p[0]=brushpts[i]->x;
//		pre_p[1]=brushpts[i]->y;
//
//		normalize(brushDir);
//		norm.entry[0]=-brushDir.entry[1];
//		norm.entry[1]=brushDir.entry[0];
//		
//		if(Num_SmoothRegionpoints>=MaxNumPoints)
//		{
//			MaxNumPoints += 50;
//			point = (Point*)realloc(point, sizeof(Point) * MaxNumPoints);
//		}
//		point[Num_SmoothRegionpoints].x=
//			brushpts[i]->x+BoundRegionWidth*quadmesh->xinterval*norm.entry[0];
//		point[Num_SmoothRegionpoints].y=
//			brushpts[i]->y+BoundRegionWidth*quadmesh->xinterval*norm.entry[1];
//		point[Num_SmoothRegionpoints].cellid=get_cellID_givencoords(point[Num_SmoothRegionpoints].x,
//			point[Num_SmoothRegionpoints].y);
//
//		Num_SmoothRegionpoints++;
//	}
//}
//
//
/*we would like to do the evenly sampling along the brush curve
please refer to the even streamline placement*/
void sample_brush(double samp_interval)
{
	int i;
	double cur_length = 0;
	icVector2 linedir;
	double t;
	double len;
	double samp[2]={0.};

	for(i=0; i<nbrushpts-1; i++)
	{
		linedir.entry[0]=brushpts[i+1]->x-brushpts[i]->x;
		linedir.entry[1]=brushpts[i+1]->y-brushpts[i]->y;

		len=length(linedir);
		cur_length+=len;

LL:		if(cur_length>samp_interval)
		{
			t=(cur_length-samp_interval)/len;

			/*we then obtain one sample points*/
			samp[0]=t*brushpts[i]->x+(1-t)*brushpts[i+1]->x;
			samp[1]=t*brushpts[i]->y+(1-t)*brushpts[i+1]->y;

			/*add to the control point list*/
			add_to_shapeCtrPtsList(samp[0], samp[1], 
				get_cellID_givencoords(samp[0], samp[1]));

			cur_length -= samp_interval;
			goto LL;

		}
	}

	/*we also need to add the last point*/
	add_to_shapeCtrPtsList(brushpts[nbrushpts-1]->x, brushpts[nbrushpts-1]->y, 
		get_cellID_givencoords(brushpts[nbrushpts-1]->x, brushpts[nbrushpts-1]->y));

}

//void sample_brush_weak(double samp_interval)
//{
//	int i;
//	double cur_length = 0;
//	icVector2 linedir;
//	double t;
//	double len;
//	double samp[2]={0.};
//
//	for(i=0; i<nbrushpts-1; i++)
//	{
//		linedir.entry[0]=brushpts[i+1]->x-brushpts[i]->x;
//		linedir.entry[1]=brushpts[i+1]->y-brushpts[i]->y;
//
//		len=length(linedir);
//		cur_length+=len;
//
//LL:		if(cur_length>samp_interval)
//		{
//			t=(cur_length-samp_interval)/len;
//
//			/*we then obtain one sample points*/
//			add_to_shapeCtrPtsList(brushpts[i]->x, brushpts[i]->y, 
//				get_cellID_givencoords(brushpts[i]->x, brushpts[i]->y));
//
//			cur_length -= samp_interval;
//			goto LL;
//
//		}
//
//	}
//
//	/*we also need to add the last point*/
//	add_to_shapeCtrPtsList(brushpts[nbrushpts-1]->x, brushpts[nbrushpts-1]->y, 
//		get_cellID_givencoords(brushpts[nbrushpts-1]->x, brushpts[nbrushpts-1]->y));
//
//}
//
//
//
////temparory coordinates transferring for second window (Conley graph window)
//11/06/05
void ScreenToSecondWin(
   int px, int py,
   int screen_leftx, int screen_bottomy,
   int win_screen_sizex, int win_screen_sizey,
   double world_leftx, double world_bottomy,
   double win_world_sizex, double win_world_sizey,
   double &s, double &t)
{
	double ratiox = (double)(px - screen_leftx) / (double)win_screen_sizex;
	double ratioy = (double)(screen_bottomy - py)/(double)win_screen_sizey;
    
	s = (double)world_leftx + ratiox * win_world_sizex;
	t = (double)world_bottomy + ratioy * win_world_sizey;
}


// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()


// CIBFVDlg dialog



CIBFVDlg::CIBFVDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CIBFVDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);

	saveprojsetDlgID=IDD_DIALOG_SAVEPROJSETTING;
	saveprojsetDlg=new SaveProjSettingDlg();
}

void CIBFVDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_MYTAB, m_tbCtrl);
	DDX_Control(pDX, IDC_OPENGLWIN, m_ctrlOpenGlWin);
	DDX_Control(pDX, IDCANCEL, trace);
}

BEGIN_MESSAGE_MAP(CIBFVDlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	ON_WM_TIMER()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_WM_MOUSEWHEEL()
	ON_WM_RBUTTONDOWN()
	ON_WM_MBUTTONDOWN()
	ON_WM_MBUTTONUP()
	ON_BN_CLICKED(IDC_BUTTON_CLEARALL, OnBnClickedButtonClearall)
	ON_COMMAND(ID_FILE_NEW32771, OnFileNew32771)
	ON_COMMAND(ID_OPEN_LOADWATERMAP, OnOpenLoadwatermap)
	ON_COMMAND(ID_FILE_EXIT, OnFileExit)
	ON_BN_CLICKED(IDC_BUTTON_SAVEPROJECT, OnBnClickedButtonSaveproject)
	ON_BN_CLICKED(IDC_BUTTON_LOADPROJECT, OnBnClickedButtonLoadproject)
	ON_NOTIFY(TCN_SELCHANGE, IDC_MYTAB, &CIBFVDlg::OnTcnSelchangeMytab)
	ON_BN_CLICKED(IDCANCEL, &CIBFVDlg::OnBnClickedCancel)
	ON_BN_CLICKED(IDC_BUTTON_PLACETENSORLINES, &CIBFVDlg::OnBnClickedButtonPlacetensorlines)
	ON_BN_CLICKED(IDC_CHECK_TENSORDESIGNON, &CIBFVDlg::OnBnClickedCheckTensordesignon)
	ON_BN_CLICKED(IDC_RADIO_ADDAREGULAR, &CIBFVDlg::OnBnClickedRadioAddaregular)
	ON_BN_CLICKED(IDC_RADIO_ADDACENTER, &CIBFVDlg::OnBnClickedRadioAddacenter)
END_MESSAGE_MAP()


// CIBFVDlg message handlers

BOOL CIBFVDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	//MainMenu *mmenu=new MainMenu();
	//MenuItem *item = mmenu->MenuItems->Add("File");
	//item->MenuItems->Add("Open...",new EventHandler(this,&NForm::OnOpen));
	//CCeCommandBar *pCommandBar = (CCeCommandBar*)m_pWndEmptyCB;
	//pCommandBar->InsertMenuBar(IDR_MENUBAR);

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here
		


	sharedvars.ShowTheMapOn=true;
		
	fittedmap1=(unsigned char*)malloc(sizeof(unsigned char)*3*512*512);
	popdensitymap_fit=(unsigned char*)malloc(sizeof(unsigned char)*3*512*512);
	vegmap_fit=(unsigned char*)malloc(sizeof(unsigned char)*3*512*512);
		displaymap=(unsigned char*)malloc(sizeof(unsigned char)*3*512*512);
		streetmapbackground=(unsigned char*)malloc(sizeof(unsigned char)*3*512*512);



	/*create a quad mesh if not exist 10/04/2007*/

	quadmesh = new QuadMesh(100, 100, -0.01, 1.01, -0.01, 1.01);
	
	//pre_ten=(icMatrix2x2*)malloc(sizeof(icMatrix2x2)*quadmesh->nverts);

	init_verts_all();
	reset_inland();
	init_ten_designelems();
	tensor_init_tex();
	init_degpts();

	init_scalar_singular_elemlist();

	//init_streetnet();
	reset_roadnetvis();
	init_evenplace_ten();
	init_evenplace_ten();
	brushinterfaceOn=false;
	flag_loadmap = false;
	flag_loadvegmap=false;


	zoom_factor = 1;
	trans_x=trans_y=0;
	cal_inverse_transform();

	init_sharedvars();

	if(brushes == NULL)
	{
		init_brushlist();
	}
	else
	{
		delete brushes;
		brushes=NULL;
		init_brushlist();
		init_sketchnet();
	}

	num_shapecontrol_pts = 0;
	num_curvepts_output = 0;
	nbrushpts = 0;

	if(seedsalongbounds!=NULL)
	{
		delete seedsalongbounds;
		seedsalongbounds=NULL;
	}

	release_domain_boundaries();
	init_domain_boundaries();
	highwayexisted=false;
	majorroadsexisted=false;
	gen_regElem_sketches=0;
	please_comb_prefield=false;

	street_sample_interval = quadmesh->xinterval/4;

		if(pre_tenfield!=NULL)
			delete [] pre_tenfield;
		pre_tenfield=new icMatrix2x2[quadmesh->nverts];



	////Get the size of the whole dialog
	RECT rect;


	////intialize pcGlView object to show the flow visualization
	CStatic *pclStatic = (CStatic *)GetDlgItem(IDC_OPENGLWIN);
	//pclStatic->MoveWindow(10, 10, 512, 512);
	pclStatic->MoveWindow(10, 10, REALWINSIZE, REALWINSIZE);
	pclGlView = new CGlView(pclStatic);

	init_ten_designelems();
	init_scalar_singular_elemlist(); /* initialize the scalar field design element list */

	pclGlView->showTensorOn = 1;
	
	g_pclGlView = pclGlView;

	init_sharedvars();

    
	LPPOINT ptSomePoint = &firstwin_leftbottom;
	
	/*------------------------------------------*/
	////Get the size and position of the window
	m_ctrlOpenGlWin.GetClientRect(&rect);
	firstwin_leftx = rect.left;
	firstwin_rightx = rect.right;
	firstwin_bottomy = rect.bottom;
	firstwin_topy = rect.top;
	
	firstwin_leftbottom.x = rect.left;
	firstwin_leftbottom.y = rect.bottom;

	m_ctrlOpenGlWin.ScreenToClient(ptSomePoint);
	/*------------------------------------------*/

	wglMakeCurrent(pclGlView->m_hDC, pclGlView->m_hglRC);
	pclGlView->OnCreate();
	pclGlView->makePatterns();
	pclGlView->getInitTex();
	pclGlView->DrawGLScene(GL_RENDER);

	/*initialize tensor field visualization*/
	make_tens_Patterns();
	tensor_init_tex();
	
	m_MiddleButtonDown=m_LeftButtonDown=m_RightButtonDown=FALSE;
	upstreaming_edit=false;

	m_tbCtrl.InitDialogs();

	m_tbCtrl.InsertItem(0,"Tensor field editing");
	m_tbCtrl.InsertItem(1,"Street network panel");
	m_tbCtrl.InsertItem(2,"Visualization");
	//m_tbCtrl.InsertItem(3,"Street edit panel");

	m_tbCtrl.ActivateTabDialogs();

	g_m_tbCtrl=&m_tbCtrl;

	saveprojsetDlg->Create(saveprojsetDlgID, GetParent());

	SetTimer(1,30,NULL);

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CIBFVDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CIBFVDlg::OnPaint() 
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
	
	pclGlView->DrawGLScene(GL_RENDER);
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CIBFVDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

void CIBFVDlg::OnTimer(UINT nIDEvent)
{
	// TODO: Add your message handler code here and/or call default

	CDialog::OnTimer(nIDEvent);
	OnPaint();
}

void CIBFVDlg::OnLButtonDown(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	double s, t;

	GLdouble position[3];

	ScreenToWorld(position, point);
	s = position[0];
	t = position[1];


	//if( ((s < 0) || (s > 1) || (t < 0) || (t > 1))
	//	&&((point.x < secondwin_leftx+secondwin_leftbottom.x)
	//			||(point.x > secondwin_rightx+secondwin_leftbottom.x)
	//			||(point.y > secondwin_leftbottom.y)
	//			||(point.y < secondwin_leftbottom.y - secondwin_bottomy)))
	//{
	//	return;
	//}

	if( ((s < 0) || (s > 1) || (t < 0) || (t > 1)))
	{
		return;
	}

	m_LeftButtonDown = TRUE;

	if((s >= 0) && (s <= 1) && (t >= 0) && (t <= 1))
	{
		////If it is under adding element mode
		save_cur_field();


		//if it is tensor field design
		if(sharedvars.TensorDesignOn && sharedvars.tenElemType >= 0 && sharedvars.tenElemType  < 5)
		{

			s=(s-.5-trans_x)/zoom_factor+.5;
			t=(t-.5-trans_y)/zoom_factor+.5;

			addto(s, t,	get_cellID_givencoords(s,t), sharedvars.tenElemType);

			/*using quad mesh 09/25/2007*/
			//cal_alltensor_quad();
			if(sharedvars.EnableSketchBasedDesign||sharedvars.ShowSketchesOn)
			{
				/*  if we perform local tensor field design  */
				cal_tensorvals_quad_inReg(cur_chosen_region);
				cal_eigenvecs_quad_inReg(cur_chosen_region);
			}
			else{
				cal_tensorvals_quad_inReg(/*cur_chosen_region*/);
				cal_all_eigenvecs_quad();
			}

			init_degpts();
			/*render alpha map*/
			render_alpha_map_quad(false);
			render_alpha_map_quad(true);

			locate_degpts_cells_tranvec_quad();

			save_tenField_perVer("temp_tenfield_perver.txt");
		}

		else if(sharedvars.TensorDesignOn && sharedvars.tenElemType >= 5) /*it is regular element*/
		{
			s=(s-.5-trans_x)/zoom_factor+.5;
			t=(t-.5-trans_y)/zoom_factor+.5;

			m_LeftDownPos = point;
			set_ten_regBasis(s, t, 0);  /*set the basis of the regular element*/
			m_LeftButtonDown = TRUE;
		}

		else if(sharedvars.EnableHeighfieldDesignOn)
		{
			s=(s-.5-trans_x)/zoom_factor+.5;
			t=(t-.5-trans_y)/zoom_factor+.5;

			if(!sharedvars.rdScalarSingularElem) /*maximum*/
			{
				add_new_scalar_singular_elem(s,t, 0);
			}
			else  /*minimum*/
				add_new_scalar_singular_elem(s,t, 1);

			/*  calculate the phi values  */
			cal_phi_allverts();

			normalize_scalar_phi();

			/*  calculate the eigen vectors */
			//cal_all_eigenvecs_asym();

			//render_alpha_map_quad(false);
			//render_alpha_map_quad(true);
		}

		/*mouse picking mode*/
		else
		{
			////Element deletion
			if(sharedvars.RemoveElemOn)
			{
				pclGlView->ten_hitProcess(s, t);

				//if(chosen_tenelem_ID < NAMEOFREGELEM && 
				//	chosen_tenelem_ID >= 0)
				//{
				//	cur_chosen_region = ten_designelems[chosen_tenelem_ID].which_region;
				//}
				//else if (chosen_tenelem_ID >= NAMEOFREGELEM)
				//{
				//	cur_chosen_region = 
				//		ten_regularelems[chosen_tenelem_ID-NAMEOFREGELEM-1].which_region;
				//}
				
				/*using quad mesh 09/25/2007*/
				//////cal_tensorvals_quad_inReg(/*cur_chosen_region*/);
				//////cal_all_eigenvecs_quad();
				////////cal_eigenvecs_quad_inReg(cur_chosen_region);
				if(sharedvars.EnableSketchBasedDesign||sharedvars.ShowSketchesOn)
				{
					/*  if we perform local tensor field design  */
					cal_tensorvals_quad_inReg(cur_chosen_region);
					cal_eigenvecs_quad_inReg(cur_chosen_region);
				}
				else{
					cal_tensorvals_quad_inReg(/*cur_chosen_region*/);
					cal_all_eigenvecs_quad();
				}

				init_degpts();
				/*render alpha map*/
				render_alpha_map_quad(false);
				render_alpha_map_quad(true);

				locate_degpts_cells_tranvec_quad();
			}

			/*tensor field: brush interface 10/10/2007*/
			else if(sharedvars.BrushInterfaceOn || sharedvars.BrushLikeRegSelOn)
			{
				if(created == 1)
				{
					num_shapecontrol_pts = 0;
					num_curvepts_output = 0;
		
					created = 0;
				}

				if(pclGlView->FinisheCtrptsel == 1)
				{
					nbrushpts=0;
					num_shapecontrol_pts=0;
					pclGlView->FinisheCtrptsel = 0;
				}

				
				add_to_shapeCtrPtsList(s, t, get_cellID_givencoords(s,t)
					/*get_cellID_picking(s, t)*/);

				if(num_shapecontrol_pts>1)
				{
					resolution = get_Resolution_adp();

					CalOpenHermiteCurve();
				}
			}

			/*   sketch based design   */
			else if(sharedvars.EnableSketchBasedDesign)
			{
				/*do some initialization here*/
			}

			/*other editing (say graph level editing)*/
			else //if(sharedvars.)
			{
				/*do some initialization here*/
				pclGlView->ten_hitProcess(s, t);
				
				if(chosen_tenelem_ID < NAMEOFREGELEM && 
					chosen_tenelem_ID >= 0)
				{
					cur_chosen_region = ten_designelems[chosen_tenelem_ID].which_region;
				}
				else if (chosen_tenelem_ID >= NAMEOFREGELEM)
				{
					cur_chosen_region = 
						ten_regularelems[chosen_tenelem_ID-NAMEOFREGELEM-1].which_region;
				}
			}
			
			if(pclGlView->showTensorOn == TRUE)/*for tensor field editing 09/20/2007*/
				pclGlView->ten_hitProcess(s, t);
			
			else if(pclGlView->SmoothOn == 1 && m_RightButtonDown == TRUE)
			{
				pclGlView->PickPointOn = 1;
				InitRegionSmooth();

				init_quad_regionsmooth();

				m_RightButtonDown = FALSE;
			}
		}

	}


	CDialog::OnLButtonDown(nFlags, point);
}

void CIBFVDlg::OnLButtonUp(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	double  s, t, s_down, t_down;
	GLdouble position[3];

	ScreenToWorld(position, m_LeftDownPos);
	s_down = position[0];
	t_down = position[1];

	ScreenToWorld(position, point);
	s = position[0];
	t = position[1];

	if ((s == s_down) && (t == t_down))
	{
		m_LeftButtonDown = FALSE;
		return;
	}


	if( (s < 0) || (s > 1) || (t < 0) || (t > 1))
	{
		//if( s < 0) s = 0;
		//if( s > 1) s = 1;
		//if( t < 0) t = 0;
		//if( t > 1) t = 1;
		m_LeftButtonDown = FALSE;
		return;
	}

	/*for brush stroke interface 10/20/2007*/
	if(m_LeftButtonDown == TRUE && sharedvars.BrushInterfaceOn)
	{
		/*   add the last point   */
		if (fabs(s - s_old)>1.e-7 && fabs(t - t_old)>1.e-7)
		{
			add_to_brushEdPts(s, t);
		}

		num_shapecontrol_pts = 0;
		sample_brush(5*quadmesh->xinterval);
		
		resolution = get_Resolution_adp();

		int returnMsg=MessageBox("Want a closed loop?", 0, MB_YESNOCANCEL/*MB_OKCANCEL*/);

		if(returnMsg==IDYES)
		{
			CalHermiteCurve();
			closedbrush = true;
		}
		else if(returnMsg==IDNO)
		{
			CalOpenHermiteCurve();
			closedbrush=false;
		}
		else
		{
			/*  we cancel the obtained stroke  */
			num_shapecontrol_pts = 0;
				num_curvepts_output = 0;

			/*reset the brush point list*/
			nbrushpts = 0;
				Num_SmoothRegionpoints=0;
		}

		save_brush_curve("temp_brushcurve.bru");
		
		/*calculate the cell strip that contains the curve*/
		//init_dis_verts();
		//init_dis_cells();
		reset_vert_dis();
		get_cellstrip_curve_quad(boundarycells, nboundarycells);



		UpdateData(FALSE);
    
		//finish = clock();
		//m_timespent = (double)(finish - start)/CLOCKS_PER_SEC;

		pclGlView->FinisheCtrptsel = 1;
	}

	/*   for brush stroke interface for the local editing 1/20/2008   */
	else if(m_LeftButtonDown == TRUE && sharedvars.BrushInterfaceOn)
	{
		if (fabs(s - s_old)>1.e-7 && fabs(t - t_old)>1.e-7)
		{
			/*   add the last point   */
			add_to_brushEdPts(s, t);
		}

		num_shapecontrol_pts = 0;
		sample_brush(5*quadmesh->xinterval);
		
		resolution = get_Resolution_adp();

		CalOpenHermiteCurve();
		closedbrush=false;
		save_brush_curve("temp_brushcurve.bru");
		
		/*calculate the cell strip that contains the curve*/
		reset_vert_dis();
		get_cellstrip_curve_quad(boundarycells, nboundarycells);

		UpdateData(FALSE);

		pclGlView->FinisheCtrptsel = 1;
	}

	else if( m_LeftButtonDown == TRUE && sharedvars.EnableSketchBasedDesign && 
		!sharedvars.TensorDesignOn)
	{
		if (fabs(s - s_old)<1.e-7 && fabs(t - t_old)<1.e-7)
		{
			/*   add the last point   */
			int cell=get_cellID_givencoords(s,t);
			if(cell>=0&&cell<quadmesh->nfaces)
				add_to_current_brushPts(s, t);
		}

		/*we need to terminate one curve*/
		brushes->nbrushes++;
		//m_edNumBrushes=brushes->nbrushes;
		(*g_m_edNumBrushes)=brushes->nbrushes;

		if(brushes->nbrushes>=brushes->curMaxNum)
		{
			if(!brushes->extend())
			{
				MessageBox("Not enough memory for new brush!");
				m_LeftButtonDown = FALSE;
				return;
			}
		}

		if(brushes->nbrushes>0)
			convert_one_brush_to_atraj(brushes->nbrushes-1);
		
		//init_sketchnet();
		//update_sketchcurveinfo();
		//cal_sketchlines_intersects();
		//search_for_sketchbased_graph();

		g_mydlg1->UpdateData(FALSE);
		UpdateData(FALSE);
	}



    m_LeftButtonDown = FALSE;

	CDialog::OnLButtonUp(nFlags, point);
}

void CIBFVDlg::OnMouseMove(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	double s, t;

	GLdouble position[3];

	ScreenToWorld(position, point);
	s = position[0];
	t = position[1];

	if(pan_s_old==s && t== pan_t_old && m_MiddleButtonDown == TRUE)
		return;
	
	if( (s < 0) || (s > 1) || (t < 0) || (t > 1))
	{
		if( s < 0) s = 0;
		if( s > 1) s = 1;
		if( t < 0) t = 0;
		if( t > 1) t = 1;
	}
	
	double pan_dx=s-pan_s_old;
	double pan_dy=t-pan_t_old;
	pan_s_old=s;
	pan_t_old=t;


	s=(s-.5-trans_x)/zoom_factor+.5;
	t=(t-.5-trans_y)/zoom_factor+.5;

	//if ((s == s_old) && (t == t_old) && m_LeftButtonDown == TRUE)
	//	return;
	if ((fabs(s - s_old)<1.e-7 && fabs(t - t_old)<1.e-7) && m_LeftButtonDown == TRUE)
		return;
	
	double dx = s - s_old, dy = t - t_old;

			//dx=(dx-.5-trans_x)/zoom_factor+.5;
			//dy=(dy-.5-trans_y)/zoom_factor+.5;
	
	save_cur_field();

	if ( m_LeftButtonDown == TRUE && sharedvars.TensorDesignOn
		&& sharedvars.tenElemType >= 5 /*&& pclGlView->MoveElemOn == 0*/) /*it is tensor regular element design*/
	{
		/*record the direction and the end point for the direction and the strength of the element*/
		set_ten_regDir(s, t);  /*set the end point and direction of the regular element*/
		
		/*using quad mesh 09/25/2007*/
		//////////cal_alltensor_quad();
		////////cal_tensorvals_quad_inReg(/*cur_chosen_region*/);
		////////cal_all_eigenvecs_quad();
		//////////cal_eigenvecs_quad_inReg(cur_chosen_region);
			if(sharedvars.EnableSketchBasedDesign||sharedvars.ShowSketchesOn)
			{
				/*  if we perform local tensor field design  */
				cal_tensorvals_quad_inReg(cur_chosen_region);
				cal_eigenvecs_quad_inReg(cur_chosen_region);
			}
			else{
				cal_tensorvals_quad_inReg(/*cur_chosen_region*/);
				cal_all_eigenvecs_quad();
			}

		init_degpts();

		/*calculate the alpha map here*/
		render_alpha_map_quad(false);
		render_alpha_map_quad(true);
			
		locate_degpts_cells_tranvec_quad();

	}

	else if( m_LeftButtonDown == TRUE && chosen_tenelem_ID >= 0 && pclGlView->MoveElemOn == 1)
	{  ////This is for dragging the center of selected element
		//update the center of the dragged element
		if( chosen_tenelem_ID >= 0 &&  chosen_tenelem_ID < NAMEOFREGELEM)
		{
			ten_designelems[chosen_tenelem_ID].centerx += dx;
			ten_designelems[chosen_tenelem_ID].centery += dy;

			////we also need to update the editbox associated with this element 04/26/05
			init_tenelem_EditBox(chosen_tenelem_ID, ten_designelems[chosen_tenelem_ID].centerx, 
				ten_designelems[chosen_tenelem_ID].centery);
			update_tenElem_EditBox(chosen_tenelem_ID);

			cur_chosen_region=ten_designelems[chosen_tenelem_ID].which_region
				=get_region_id(ten_designelems[chosen_tenelem_ID].centerx,
				ten_designelems[chosen_tenelem_ID].centery);

		}

		else{
			/*it is regular element*/

			if(chosen_tenelem_ID <= NAMEOFSINGCONTROL && chosen_tenelem_ID >= NAMEOFREGELEM)
			{
				ten_regularelems[chosen_tenelem_ID-NAMEOFREGELEM].base[0] += dx;
				ten_regularelems[chosen_tenelem_ID-NAMEOFREGELEM].base[1] += dy;
				ten_regularelems[chosen_tenelem_ID-NAMEOFREGELEM].end[0] += dx;
				ten_regularelems[chosen_tenelem_ID-NAMEOFREGELEM].end[1] += dy;

				cur_chosen_region=ten_regularelems[chosen_tenelem_ID-NAMEOFREGELEM].which_region
					=get_region_id(ten_regularelems[chosen_tenelem_ID-NAMEOFREGELEM].base[0],
					ten_regularelems[chosen_tenelem_ID-NAMEOFREGELEM].base[1]);

			}
		}

		rotateAng = 0;

		/*using quad mesh 09/25/2007*/
		////////cal_alltensor_quad();
		//////cal_tensorvals_quad_inReg(/*cur_chosen_region*/);
		//////cal_all_eigenvecs_quad();
		////////cal_eigenvecs_quad_inReg(cur_chosen_region);
			if(chosen_tenelem_ID < NAMEOFREGELEM && 
				chosen_tenelem_ID >= 0)
			{
				cur_chosen_region = ten_designelems[chosen_tenelem_ID].which_region;
			}
			else if (chosen_tenelem_ID >= NAMEOFREGELEM)
			{
				cur_chosen_region = 
					ten_regularelems[chosen_tenelem_ID-NAMEOFREGELEM].which_region;
			}
			if(sharedvars.EnableSketchBasedDesign||sharedvars.ShowSketchesOn)
			{
				/*  if we perform local tensor field design  */
				cal_tensorvals_quad_inReg(cur_chosen_region);
				cal_eigenvecs_quad_inReg(cur_chosen_region);
			}
			else{
				cal_tensorvals_quad_inReg(/*cur_chosen_region*/);
				cal_all_eigenvecs_quad();
			}

		init_degpts();

		/*calculate the alpha map here*/
		render_alpha_map_quad(false);
		render_alpha_map_quad(true);
			
		locate_degpts_cells_tranvec_quad();
		
		OnPaint(); 
	}

	else if( m_LeftButtonDown == TRUE && chosen_tenelem_ID >= 0 && pclGlView->EditModeOn == 1
		&& pclGlView->showTensorOn == 1)
	{ 
		////This is for editing the selected element
		////Using tracing ball simulation to get the rotation matrix

		////Perform coordinates transformation to transform the s_old, s, t_old and t to the
		////new local frame system of current being edited element 07/07/05

		////1. we need to get the rotation information of current edited element
		//icMatrix3x3 inverse_rot;
		double cur_ang;

		if(chosen_tenelem_ID < NAMEOFREGELEM)
			cur_ang = ten_designelems[chosen_tenelem_ID].rotang;            //// singular element
		else{
			cur_ang = ten_regularelems[chosen_tenelem_ID - NAMEOFREGELEM].rotang /*+ ini_ang*/; ////regular element
		}


		double ldx, ldy;
		ldx = cos(cur_ang)*dx + sin(cur_ang)*dy;
		ldy = -sin(cur_ang)*dx + cos(cur_ang)*dy;

		if(chosen_tenelem_ID - NAMEOFREGELEM>=0) /*we select a regular element*/
		{
			update_ten_regDir(chosen_tenelem_ID-NAMEOFREGELEM, s, t);
		}

		////3. use the local coordinates to get the offset value of transformation
		else if(pclGlView->which_control_point == UPPERROTATE || pclGlView->which_control_point == ARROWHEAD )
		{
			double rotang = 0;

			////Get the rotation angle
			if(fabs(ldx) >= fabs(ldy))
			{
				rotateAng = -ldx*ROTATESTEP;
			}
			else
			{
				rotateAng = -ldy*ROTATESTEP;
			}

			uniforms = sx = sy = 0;
		}

		////Set the scalar here
		else if(pclGlView->which_control_point == LOWLEFT ||
			pclGlView->which_control_point == UPPERLEFT ||
			pclGlView->which_control_point == UPPERRIGHT ||
			pclGlView->which_control_point == LOWRIGHT)
		{
			////uniform scaling here
			if(pclGlView->which_control_point == LOWLEFT)
			{
				if(ldx < 0 && ldy < 0)
					uniforms = sqrt(ldx * ldx + ldy * ldy) * SCALESTEP/3;
				else if(ldx > 0 && ldy > 0 )
					uniforms = -sqrt(ldx * ldx + ldy * ldy) * SCALESTEP/3;
			}
			else if(pclGlView->which_control_point == UPPERLEFT)
			{
				if(ldx < 0 && ldy > 0)
					uniforms = sqrt(ldx * ldx + ldy * ldy) * SCALESTEP/3;
				else if(ldx > 0 && ldy < 0 )
					uniforms = -sqrt(ldx * ldx + ldy * ldy) * SCALESTEP/3;
			}
			else if(pclGlView->which_control_point == UPPERRIGHT)
			{
				if(ldx < 0 && ldy < 0)
					uniforms = -sqrt(ldx * ldx + ldy * ldy) * SCALESTEP/3;
				else if(ldx > 0 && ldy > 0 )
					uniforms = sqrt(ldx * ldx + ldy * ldy) * SCALESTEP/3;
			}
			else{
				if(ldx < 0 && ldy > 0)
					uniforms = -sqrt(ldx * ldx + ldy * ldy) * SCALESTEP/3;
				else if(ldx > 0 && ldy < 0 )
					uniforms = sqrt(ldx * ldx + ldy * ldy) * SCALESTEP/3;
			}

			sy = sx = 0.;
			rotateAng = 0;
		}

		else if(pclGlView->which_control_point == LEFT || pclGlView->which_control_point == RIGHT)
		{
			if(pclGlView->which_control_point == LEFT)
				sx = - ldx*SCALESTEP;
			else
				sx = + ldx*SCALESTEP;

			uniforms = sy = 0.;
			rotateAng = 0;

		}
		else if(pclGlView->which_control_point == UPPER || pclGlView->which_control_point == BUTTOM)
		{
			if(pclGlView->which_control_point == UPPER)
				sy = ldy*SCALESTEP;
			else
				sy = -ldy*SCALESTEP;
			
			uniforms = sx = 0.;
			rotateAng = 0;
		}

		else if(pclGlView->which_control_point == ARROWBASE)
		{
			uniforms = ldx * SCALESTEP;
			sy = sx = 0.;
			rotateAng = 0;
	    }
		
		update_tenElem_Transform(pclGlView->TransformType, chosen_tenelem_ID);

		/*using quad mesh 09/25/2007*/
		////////cal_alltensor_quad();
		//////cal_tensorvals_quad_inReg(/*cur_chosen_region*/);
		//////cal_all_eigenvecs_quad();
			if(chosen_tenelem_ID < NAMEOFREGELEM && 
				chosen_tenelem_ID >= 0)
			{
				cur_chosen_region = ten_designelems[chosen_tenelem_ID].which_region;
			}
			else if (chosen_tenelem_ID >= NAMEOFREGELEM)
			{
				cur_chosen_region = 
					ten_regularelems[chosen_tenelem_ID-NAMEOFREGELEM].which_region;
			}
			if(sharedvars.EnableSketchBasedDesign||sharedvars.ShowSketchesOn)
			{
				/*  if we perform local tensor field design  */
				cal_tensorvals_quad_inReg(cur_chosen_region);
				cal_eigenvecs_quad_inReg(cur_chosen_region);
			}
			else{
				cal_tensorvals_quad_inReg(/*cur_chosen_region*/);
				cal_all_eigenvecs_quad();
			}

		init_degpts();

		/*calculate the alpha map here*/
		render_alpha_map_quad(false);
		render_alpha_map_quad(true);
			
		locate_degpts_cells_tranvec_quad();

		OnPaint(); 
		
	}


	/*if it is under brush stroke editing mode 10/20/2007*/
	else if( m_LeftButtonDown == TRUE && (sharedvars.BrushInterfaceOn || sharedvars.BrushLikeRegSelOn)
		/*&& m_chkActiveTensorDesign == TRUE*/)
	{
		/*we may need to record the mouse movement into a point 10/20/2007*/
		//int test = 0;

		/*we first record all the points that the cursor move to*/
		add_to_brushEdPts(s, t);

		/*  compute a closed region that includes the brush */
		num_shapecontrol_pts = 0;
		sample_brush(5*quadmesh->xinterval);

		compute_brush_Reg();

		/*when user release the left button of the mouse,
		we sample those points and use them to get the Hermite curve*/
	}

	/*sketch based design 11/20/2007*/
	else if( m_LeftButtonDown == TRUE && sharedvars.EnableSketchBasedDesign)
	{
		/*  important: if (s,t) out of the design plane, don't add it!  */
		int cell=get_cellID_givencoords(s,t);
		if(cell>=0&&cell<quadmesh->nfaces)
			add_to_current_brushPts(s, t);
	}

	/*graph level editing*/
	else if(m_LeftButtonDown == TRUE && selectedIntersect >= 0
		&& sharedvars.EditStreetNetOn)
	{
		streetnet->nodelist->intersects[selectedIntersect]->gpos[0]+=dx;
		streetnet->nodelist->intersects[selectedIntersect]->gpos[1]+=dy;

	}

	if(m_MiddleButtonDown == TRUE) /*IT IS TRANSLATION*/
	{	
			//dx=(dx-.5-trans_x)/zoom_factor+.5;
			//dy=(dy-.5-trans_y)/zoom_factor+.5;
		trans_x += pan_dx;
		trans_y += pan_dy;
	}
	
	s_old = s;
	t_old = t;


	CDialog::OnMouseMove(nFlags, point);
}

BOOL CIBFVDlg::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt)
{
	// TODO: Add your message handler code here and/or call default
	double s, t;

	GLdouble position[3];

	ScreenToWorld(position, pt);
	s = position[0];
	t = position[1];

	if( ((s < 0) || (s > 1.2) || (t < -0.1) || (t > 0.9)))
	{
		return CDialog::OnMouseWheel(nFlags, zDelta, pt);
	}

	if(zDelta>0)  /*zoom in*/
	{
		zoom_factor +=0.05;
	}

	else  /*zoom out*/
	{
		if(zoom_factor>0.96)
			zoom_factor -=0.05;
	}

	mouse_zDelta=zDelta;

	return CDialog::OnMouseWheel(nFlags, zDelta, pt);
}

void CIBFVDlg::OnRButtonDown(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	if(pclGlView->SmoothOn == 1)
	{
		m_RightButtonDown = TRUE;
	}

	CDialog::OnRButtonDown(nFlags, point);
}

void CIBFVDlg::OnMButtonDown(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	m_MiddleDownPos = point;
	
	double s, t;

	GLdouble position[3];

	ScreenToWorld(position, point);
	s = position[0];
	t = position[1];


	//if( ((s < 0) || (s > 1) || (t < 0) || (t > 1))
	//	&&((point.x < secondwin_leftx+secondwin_leftbottom.x)
	//			||(point.x > secondwin_rightx+secondwin_leftbottom.x)
	//			||(point.y > secondwin_leftbottom.y)
	//			||(point.y < secondwin_leftbottom.y - secondwin_bottomy)))
	if( ((s < 0) || (s > 1) || (t < 0) || (t > 1)))
	{
		return;
	}


	//if((s >= 0) && (s <= 1) && (t >= 0) && (t <= 1))
	//{
	//}
	
	s_old = s;
	t_old = t;
	m_MiddleButtonDown = TRUE;

	CDialog::OnMButtonDown(nFlags, point);
}

void CIBFVDlg::OnMButtonUp(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default

	m_MiddleButtonDown = FALSE;
	CDialog::OnMButtonUp(nFlags, point);
}


void clean()
{
	g_pclGlView->InitFlag();
	
	////Initialize variables for transformation
	sx = sy = uniforms = 1;
    rotateAng = 0;
	
	clear_all_tens();
	init_verts_all();
	init_ten_designelems();
	tensor_init_tex();
	init_degpts();

	init_scalar_singular_elemlist();

	//init_streetnet();
	reset_roadnetvis();
	init_evenplace_ten();
	//init_roadnetvis();
	brushinterfaceOn=false;
	//flag_loadmap = false;
	//if(seedsalongbounds!=NULL)
	//{
	//	delete seedsalongbounds;
	//	seedsalongbounds=NULL;
	//}
	reset_inland();
	flag_loadvegmap=false;


	zoom_factor = 1;
	trans_x=trans_y=0;
	cal_inverse_transform();


	init_sharedvars();

	if(brushes == NULL)
	{
		init_brushlist();
	}
	else
	{
		delete brushes;
		brushes=NULL;
		init_brushlist();
		init_sketchnet();
	}

	num_shapecontrol_pts = 0;
	num_curvepts_output = 0;
	nbrushpts = 0;


	release_domain_boundaries();
	init_domain_boundaries();
	highwayexisted=false;
	majorroadsexisted=false;
	gen_regElem_sketches=0;
	please_comb_prefield=false;
	upstreaming_edit=false;
	is_on_local_editing=false;

}

void CIBFVDlg::OnBnClickedButtonClearall()
{
	// TODO: Add your control notification handler code here
	init_streetnetwork_info_in_cells();
	quadmesh->finalize_quad_verts();
	quadmesh->finalize_quad_cells();
	if(quadmesh!=NULL)
		delete quadmesh;
	//quadmesh = new QuadMesh(60, 60, -0.01, 1.01, -0.01, 1.01);
	quadmesh = new QuadMesh(100, 100, -0.01, 1.01, -0.01, 1.01);

	pclGlView->InitFlag();
	
	////Initialize variables for transformation
	sx = sy = uniforms = 1;
    rotateAng = 0;
	
	clear_all_tens();
	init_verts_all();
	reset_inland();
	init_ten_designelems();
	tensor_init_tex();
	init_degpts();

	init_scalar_singular_elemlist();

	//init_streetnet();
	reset_roadnetvis();
	init_evenplace_ten();
	//init_roadnetvis();
	brushinterfaceOn=false;
	flag_loadmap = false;
	flag_loadvegmap=false;


	zoom_factor = 1;
	trans_x=trans_y=0;
	cal_inverse_transform();


	m_tbCtrl.Reset();
	m_tbCtrl.InitDialogs();
	m_tbCtrl.ActivateTabDialogs();

	init_sharedvars();

	if(brushes == NULL)
	{
		init_brushlist();
	}
	else
	{
		delete brushes;
		brushes=NULL;
		init_brushlist();
		init_sketchnet();
	}

	num_shapecontrol_pts = 0;
	num_curvepts_output = 0;
	nbrushpts = 0;

	if(seedsalongbounds!=NULL)
	{
		delete seedsalongbounds;
		seedsalongbounds=NULL;
	}

	release_domain_boundaries();
	init_domain_boundaries();
	highwayexisted=false;
	majorroadsexisted=false;
	gen_regElem_sketches=0;
	please_comb_prefield=false;
	upstreaming_edit=false;
    is_on_local_editing=false;
}

void CIBFVDlg::OnFileNew32771()
{
	// TODO: Add your command handler code here
	OnBnClickedButtonClearall();
}

void CIBFVDlg::OnOpenLoadwatermap()
{
	// TODO: Add your command handler code here
	CFile f;

	char strFilter[] = { "Bitmap Files (*.bmp)|*.bmp|All Files (*.*)|*.*||" };

	CFileDialog FileDlg(TRUE, ".bmp", NULL, 0, strFilter);
	if( FileDlg.DoModal() == IDOK )
	{
		CString f2 = FileDlg.GetFileName();
		CString f3 = f2.Left(f2.Find('.'))+".bmp";

		const char *filename = f3;
		char filename1[255];
		strcpy(filename1, filename);
		
		map1 = BmpToTexture(filename1, &map1_w, &map1_h);	

		//if(fittedmap1!=NULL)
		//	free(fittedmap1);

		//fittedmap1=(unsigned char*)malloc(sizeof(unsigned char)*3*512*512);

		/*obtain an image that fit in the 512x512 windows*/
		get_fitted_map(map1, map1_w, map1_h, fittedmap1, 512, 512);

		/*regenerate the underlying mesh*/
		if(quadmesh!=NULL)
			delete quadmesh;

		if(map1_w>=map1_h)
		{
			double ratio = (double)map1_w/map1_h;
			double yrang = 1./ratio;
			//quadmesh = new QuadMesh(100, 100/ratio, -0.0, 1.0, 
			//	0.5-yrang/2, 0.5+yrang/2.);
			//get_mask_map(-0.0, 1.0, 0.5-yrang/2, 0.5+yrang/2, map1, map1_w, map1_h);
			quadmesh = new QuadMesh(100, 100, -0.0, 1.0, 
				0, 1);
			//quadmesh = new QuadMesh(153, 153, -0.0, 1.0, 
			//	0, 1);
			get_mask_map(-0.0, 1.0, 0., 1, fittedmap1, 512, 512);
		}
		else
		{
			double ratio = (double)map1_h/map1_w;
			double xrang = 1./ratio;
			//quadmesh = new QuadMesh(100/ratio, 100, 0.5-xrang/2, 0.5+xrang/2.,
			//	-0.0, 1.0);
			//get_mask_map(0.5-xrang/2, 0.5+xrang/2, -0.02, 1.02, map1, map1_w, map1_h);
			quadmesh = new QuadMesh(100, 100, 0., 1.,
				-0.0, 1.0);
			//quadmesh = new QuadMesh(153, 153, -0.0, 1.0, 
			//	0, 1);
			get_mask_map(-0.0, 1.0, 0., 1, fittedmap1, 512, 512);
		}

		if(pre_tenfield!=NULL)
			delete [] pre_tenfield;
		pre_tenfield=new icMatrix2x2[quadmesh->nverts];

		flag_loadmap = true;
		sharedvars.ShowTheMapOn=true;
		g_mydlg2->m_chkShowTheMapOn=TRUE;
		g_mydlg2->UpdateData(FALSE);


		init_degpts();
		init_ten_designelems();
		tensor_init_tex();


		imgprocess=new ImgProcess(fittedmap1);

		g_streetnetpanel->UpdateData();
		if(sharedvars.MeanFilterOn)
			imgprocess->mean_filter(g_streetnetpanel->m_edMeanFilterSize/*3*/);
			//imgprocess->mean_filter_usesame(9);

		imgprocess->label_different_regions();
		imgprocess->extract_regions();
		imgprocess->label_boundary_pixels();

		imgboundaries=new ImgBoundary();
		//imgboundaries->get_bound_pixels(imgprocess);

		imgboundaries->get_allboundary_pixels(imgprocess);
		imgboundaries->naive_find_sortedboundaries(imgprocess);
		
		imgboundaries->render_result(imgprocess);

		////SSS_SavePPM(512,512, imgboundaries->output, "test_output.ppm");

		//displaymap=(unsigned char*)malloc(sizeof(unsigned char)*3*512*512);
		//streetmapbackground=(unsigned char*)malloc(sizeof(unsigned char)*3*512*512);

		for(int i=0;i<3*512*512;i++)
		{
			displaymap[i]=imgboundaries->output[i];
			streetmapbackground[i]=imgboundaries->streetmap_background[i];
		}


		///*get the geometry boundary*/
		mapboundarylist=new MapBoundaryList();
		mapboundarylist->get_worldcoords_unitintervals(0.,1.,0.,1.,512,512);

		mapboundarylist->downsamp_imgboundaries(imgboundaries, 5);
		mapboundarylist->index_allboundarypts();

		cal_seeds_basedon_multibounds();

		if(sharedvars.GenTenFromExtractedBoundsOn)
		{
			//obtain_smooth_region_multibounds(BoundRegionWidth);
			obtain_field_basis();
			init_degpts();
			/*render alpha map*/
			render_alpha_map_quad(false);
			render_alpha_map_quad(true);

			locate_degpts_cells_tranvec_quad();

			/*  disable the visualization of regular elements  */
			g_mydlg2->m_chkShowRegElemOn=FALSE;
			g_pclGlView->RegularElemOn=0;
			sharedvars.ShowRegElemOn=false;
		}

	
		delete imgprocess;
		imgprocess=NULL;
		delete imgboundaries;
		imgboundaries=NULL;


		sharedvars.ShowPopDensityMapOn=false;
		g_mydlg2->m_chkShowPopDensityMapOn=FALSE;
		
		g_mydlg2->m_chkShowVegMapOn=FALSE;
		sharedvars.ShowVegMapOn=false;

		g_mydlg2->m_chkShowIBFVOn=FALSE;
		g_mydlg2->m_chkShowScalarFieldOn=FALSE;
		g_mydlg2->m_chkShowExtractedBoundsOn=FALSE;
		g_mydlg2->m_chkShowRoadMapOn=FALSE;

		sharedvars.ShowIBFVOn=false;
		sharedvars.ShowScalarFieldOn=false;
		sharedvars.ShowExtractedBoundsOn=false;
		sharedvars.ShowRoadMapOn=false;
	}
}

void CIBFVDlg::OnFileExit()
{
	// TODO: Add your command handler code here
	exit(0);
}

void CIBFVDlg::OnBnClickedButtonSaveproject()
{
	// TODO: Add your control notification handler code here

	/*  we will activate a saving project setting dialog to guide the saving process  */
	if(saveprojsetDlg == NULL)
		return;

	CRect l_rectClient;
	CRect l_rectWnd;

	GetClientRect(l_rectClient);
	GetWindowRect(l_rectWnd);

	saveprojsetDlg->GetWindowRect(l_rectClient);

	saveprojsetDlg->SetWindowPos(&wndTop,l_rectWnd.left+840,l_rectWnd.top+100,l_rectClient.Width(),
		l_rectClient.Height(),SWP_SHOWWINDOW);
	saveprojsetDlg->ShowWindow(SW_SHOW);

}

void CIBFVDlg::OnBnClickedButtonLoadproject()
{
	CFile f;
	char filename[256];

	char strFilter[] = { "Text Files (*.spj)|*.spj"};

	CFileDialog FileDlg(TRUE, ".spj", NULL, 0, strFilter);
	if( FileDlg.DoModal() == IDOK )
	{

		CString f2 = FileDlg.GetFileName();
		const char *tempf = f2;
		strcpy(filename, f2);
	
		load_a_project(filename);

		/*   generate the tensor field  */
		cal_all_eigenvecs_quad();
		init_degpts();
		/*render alpha map*/
		render_alpha_map_quad(false);
		render_alpha_map_quad(true);

		locate_degpts_cells_tranvec_quad();
	}
}


void CIBFVDlg::OnTcnSelchangeMytab(NMHDR *pNMHDR, LRESULT *pResult)
{
	// TODO: Add your control notification handler code here
	*pResult = 0;
}


void CIBFVDlg::OnBnClickedCancel()
{
	// TODO: Add your control notification handler code here
	CDialog::OnCancel();
	//StreetNetPanel::OnBnClickedButtonPlacetensorlines();
}


void CIBFVDlg::OnBnClickedButtonPlacetensorlines()
{
	// TODO: Add your control notification handler code here
	g_streetnetpanel->OnBnClickedButtonPlacetensorlines();
}


void CIBFVDlg::OnBnClickedCheckTensordesignon()
{
	// TODO: Add your control notification handler code here
	g_mydlg1->OnBnClickedCheckTensordesignon();
}


void CIBFVDlg::OnBnClickedRadioAddaregular()
{
	// TODO: Add your control notification handler code here
	g_mydlg1->OnBnClickedRadioAddaregular();
}


void CIBFVDlg::OnBnClickedRadioAddacenter()
{
	// TODO: Add your control notification handler code here
	g_mydlg1->OnBnClickedRadioAddacenter();
}
