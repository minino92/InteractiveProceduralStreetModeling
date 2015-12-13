// StreetNetPanel.cpp : implementation file
//

#include "stdafx.h"
#include "IBFV2D.h"
#include "StreetNetPanel.h"
#include ".\streetnetpanel.h"

#include "VFDataStructure.h"

#include "GlView.h"
#include ".\glview.h"
#include "shareinterfacevars.h"

#include "tensorvis.h"
#include "caldeformation.h"
#include "tensordesign.h"
#include "tensoranalysis.h"
#include "computeroadvis.h"
#include "evenlystreamlines.h"
#include "regionsmooth_quad.h"

#include "ImgBoundaryExtract.h"
#include "BmpProcess.h"

#include "BoundaryBasedTenGen.h"

#include "Loadmaps.h"

#include "Sketchdesign.h"
extern int NDesignRegions;
extern int cur_chosen_region; 

#include "MyTabCtrl.h"
extern MyTabCtrl *g_m_tbCtrl;

#include "Mydlg1.h"
#include "MyDlg2.h"
#include "StreetEditPanel.h"
#include "LocalRegionEdit.h"

extern MyDlg1 *g_mydlg1;
extern MyDlg2 *g_mydlg2;
extern StreetNetPanel *g_streetnetpanel;
extern StreetEditPanel *g_streeteditpanel;

//#include "ppm.h"

extern SharedInterfaceVars sharedvars;
extern CGlView *g_pclGlView;
extern StreetVis *roadnetvis;
extern EvenStreamlinePlace *major, *minor;
extern int *boundarycells;
extern int nboundarycells;
extern bool closedbrush;
extern QuadMesh *quadmesh;

extern int REALWINSIZE;


extern bool flag_loadmap;
extern bool flag_loadvegmap;

extern unsigned char *map1, *fittedmap1, *displaymap, *streetmapbackground;
extern int map1_w, map1_h;


/*       vegetation map      */
extern unsigned char *vegmap;
extern int vegmap_w, vegmap_h;
extern unsigned char *vegmap_fit;
extern unsigned char *vegmap_disp;

/*       Population density map      */
extern unsigned char *popdensitymap;
extern int popdensitymap_w, popdensitymap_h;
extern unsigned char *popdensitymap_fit;
extern unsigned char *popdensitymap_disp;

/*       Hieght field       */
extern unsigned char *heightfield;
extern int heightfield_w, heightfield_h;
extern unsigned char *heightfield_fit;
extern unsigned char *heightfield_dis;

extern ImgProcess *imgprocess;
extern ImgBoundary *imgboundaries;
extern MapBoundaryList *mapboundarylist;

extern void get_mask_map(double xstart, double xend, double ystart, double yend,
				  unsigned char *map, int width, int height);
extern void get_fitted_map(unsigned char *map, int width, int height,
					unsigned char *fittedmap, int fitwidth, int fitheight);
extern void set_vegflags_verts(double xstart, double xend, double ystart, double yend,
				  unsigned char *map, int width, int height);


extern double majorDensity ;
extern double minorDensity ;
extern double mintenline_length;
extern int which_level;

CEdit *g_m_EdCtrlSetStreetLevels;
int *g_m_edStreetLevel;

extern double minBoundaryLen;

/*     construct the sketch based network and its blocks    */
extern RegionBlockList *sketchblocklist;
extern StreetNet *sketchnet;
extern void construct_regionblocks_edgewise(StreetNet *net, RegionBlockList *aregionblocklist);

extern double BoundRegionWidth;
extern bool majorroadsexisted;

bool upstreaming_edit=false;
extern int num_shapecontrol_pts;
extern int num_curvepts_output;


int MeanFilterSize;

extern bool highwayexisted;

extern icMatrix2x2 *pre_tenfield;


extern double map_xrang;  // default is 10km
extern double map_yrang;  // default is 10km

extern ctr_point **brushpts ;
extern int nbrushpts ;
extern int curMaxNumBrushPts ;


extern void obtain_smooth_region_multibounds(double widtheachbound);
extern void obtain_field_basis();


//extern void reset_roadnetvis();
//extern void init_evenplace_ten();
	
extern void place();

/*  noise control global variable  */
NoiseCtrlPanel *g_noise_ctrlDlg;

extern int Num_SmoothRegionpoints;

extern double BoundRegionWidth;
extern bool is_on_local_editing;

extern bool please_comb_prefield;
extern void jitter_roads_inReg(double curScale);
extern void jitter_tens_inReg(double curScale, double curAmp);
extern void jitter_all_majRoads(double scale, double curAmp, bool majormin);

extern ctr_point *control_pts;        // allocate our control point array
extern int num_shapecontrol_pts;
extern unsigned char cur_max_reg_index;

void get_brush_approDir()
{
	int i;
	for(i=0; i<num_shapecontrol_pts-1; i++)
	{
		set_ten_regBasis(control_pts[i].x, control_pts[i].y, 0);
		set_ten_regDir(control_pts[i+1].x, control_pts[i+1].y);
	}
	
	//for(i=0; i<resolution-1; i++)
	//{
	//	set_ten_regBasis(out_pts[i].x, out_pts[i].y, 0);
	//	set_ten_regDir(out_pts[i+1].x, out_pts[i+1].y);
	//}
}

// StreetNetPanel dialog

IMPLEMENT_DYNAMIC(StreetNetPanel, CDialog)
StreetNetPanel::StreetNetPanel(CWnd* pParent /*=NULL*/)
	: CDialog(StreetNetPanel::IDD, pParent)
	//, m_edMinTenLineDensity(200.)
	//, m_edMajTenLineDensity(200.)
	, m_edMinTenLineDensity(1200)
	, m_edMajTenLineDensity(1400)
	, m_chkEditStreetNetOn(FALSE)
	, m_chkMeanFilterOn(TRUE)
	, m_chkGenTenFromExtractedBoundsOn(TRUE)
	, m_chkUseAllBoundsOn(TRUE)
	, m_edMinTensorLineLength(600.)
	, m_chkDesignGridOn(FALSE)
	//, m_edStreetLevel(2)
	, m_edStreetLevel(1)
	, m_chkConnectDeadEndsOn(FALSE)
	, m_chkShowRegionBlocksOn(FALSE)
	//, m_chkCombinePopDensityOn(FALSE)
	, m_chkCombinePopDensityOn(TRUE)
	, m_chkGenTenFromSketchesOn(TRUE)
	, m_chkCloseLoopOn(TRUE)
	, m_chkUseBoundsAsSketchesOn(TRUE)
	, m_edEffectWidthBoundary(0)
	, m_edMeanFilterSize(3)
	, m_chkShowInitSeedsOn(FALSE)
	, m_chkApplyAsymFldOn(FALSE)
	, m_chkRemoveMajDeadEndsOn(TRUE)
	, m_edMinBoundaryLen(600.)
	, m_chkSelRegToEditOn(FALSE)
	, m_chkBrushLikeRegSelOn(FALSE)
	, m_rdSubConnectMethod(FALSE)
{
	g_m_EdCtrlSetStreetLevels=&m_EdCtrlSetStreetLevels;
	g_m_edStreetLevel=&m_edStreetLevel;


	maj_road_settingDlgID=IDD_DIALOG_MAJROADSETTING;
	maj_road_settingDlg=new MajRoadSettingDlg();

	seeds_ctrlDlgID=IDD_DIALOG_SEEDPOINTCTRLDLG;
	seeds_ctrlDlg = new SeedPtCtrlDlg();


	noise_ctrlDlgID=IDD_DIALOG_NOISECTRLPANEL;
	g_noise_ctrlDlg=noise_ctrlDlg=new NoiseCtrlPanel();

	min_road_setDlgID=IDD_DIALOG_MINROADSETTING;
	min_road_setDlg=new MinRoadSetDlg();
}

StreetNetPanel::~StreetNetPanel()
{
	//if(maj_road_settingDlg!=NULL)
	//	delete maj_road_settingDlg;
}

void StreetNetPanel::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT_MINORDENSITY, m_edMinTenLineDensity);
	DDX_Text(pDX, IDC_EDIT_MAJORDENSITY, m_edMajTenLineDensity);
	DDX_Check(pDX, IDC_CHECK_EDITSTREETGRAPHON, m_chkEditStreetNetOn);
	DDX_Check(pDX, IDC_CHECK_MEANFILTERON, m_chkMeanFilterOn);
	DDX_Check(pDX, IDC_CHECK_GENTENSORBASEDONEXTRACTBOUNDARIES, m_chkGenTenFromExtractedBoundsOn);
	DDX_Check(pDX, IDC_CHECK_USEALLBOUNDARIES, m_chkUseAllBoundsOn);
	DDX_Control(pDX, IDC_SLIDER_BOUNDARYWIDTH, m_ctrlSliderBoundaryRegionWidth);
	DDX_Text(pDX, IDC_EDIT_MINTENSORLINELENGTH, m_edMinTensorLineLength);
	DDX_Check(pDX, IDC_CHECK_GRIDON, m_chkDesignGridOn);
	DDX_Text(pDX, IDC_EDIT_STREETLEVEL, m_edStreetLevel);
	DDX_Check(pDX, IDC_CHECK_CONNECTDEADENDON, m_chkConnectDeadEndsOn);
	DDX_Control(pDX, IDC_SPIN_TRACINGLEVEL, m_ctrlSpinTracingLevel);
	DDX_Check(pDX, IDC_CHECK_SHOWREGIONBLOCKS, m_chkShowRegionBlocksOn);
	DDX_Check(pDX, IDC_CHECK_COMBINEDENSITYMAP, m_chkCombinePopDensityOn);
	DDX_Check(pDX, IDC_CHECK_GENTEN_SKETCHBASED, m_chkGenTenFromSketchesOn);
	DDX_Check(pDX, IDC_CHECK_CLOSELOOP, m_chkCloseLoopOn);
	DDX_Control(pDX, IDC_EDIT_STREETLEVEL, m_EdCtrlSetStreetLevels);
	DDX_Check(pDX, IDC_CHECK_USEBOUNDARIESASSKETCHES, m_chkUseBoundsAsSketchesOn);
	DDX_Text(pDX, IDC_EDIT_WIDTHBOUNDARY, m_edEffectWidthBoundary);
	DDX_Text(pDX, IDC_EDIT2, m_edMeanFilterSize);
	DDX_Check(pDX, IDC_CHECK_SHOWINITSEEDS, m_chkShowInitSeedsOn);
	DDX_Control(pDX, IDC_BUTTON_MAJROADSETTING, m_ctrlMajRoadSettingButton);
	DDX_Control(pDX, IDC_BUTTON_SEEDSETTING, m_ctrlSeedPtSetting);
	DDX_Check(pDX, IDC_CHECK_APPLYASYMFLD, m_chkApplyAsymFldOn);
	DDX_Check(pDX, IDC_CHECK_REMOVEMAJDEADEDNON, m_chkRemoveMajDeadEndsOn);
	DDX_Text(pDX, IDC_EDIT_MINBOUNDARYLEN, m_edMinBoundaryLen);
	DDX_Check(pDX, IDC_CHECK_SELREGTOEDITON, m_chkSelRegToEditOn);
	DDX_Check(pDX, IDC_CHECK_USEBRUSHLIKEREGSEL, m_chkBrushLikeRegSelOn);
	DDX_Radio(pDX, IDC_RADIO_CONMETHOD1, m_rdSubConnectMethod);
}


void StreetNetPanel::Reset()
{
	m_edMinTenLineDensity=1.5;
	m_edMajTenLineDensity=1.5;
	m_edMinTensorLineLength=6;
}


BEGIN_MESSAGE_MAP(StreetNetPanel, CDialog)
	ON_BN_CLICKED(IDC_BUTTON_PLACETENSORLINES, OnBnClickedButtonPlacetensorlines)
	ON_BN_CLICKED(IDC_BUTTON_COMPUTESTREETNETWORK, OnBnClickedButtonComputestreetnetwork)
	ON_BN_CLICKED(IDC_CHECK_EDITSTREETGRAPHON, OnBnClickedCheckEditstreetgraphon)
	ON_BN_CLICKED(IDC_CHECK_MEANFILTERON, OnBnClickedCheckMeanfilteron)
	ON_BN_CLICKED(IDC_BUTTON_LOADAWATERMAP, OnBnClickedButtonLoadawatermap)
	ON_BN_CLICKED(IDC_CHECK_GENTENSORBASEDONEXTRACTBOUNDARIES, OnBnClickedCheckGentensorbasedonextractboundaries)
	ON_BN_CLICKED(IDC_CHECK_USEALLBOUNDARIES, OnBnClickedCheckUseallboundaries)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_BOUNDARYWIDTH, OnNMCustomdrawSliderBoundarywidth)
	ON_BN_CLICKED(IDC_CHECK_GRIDON, OnBnClickedCheckGridon)
	ON_BN_CLICKED(IDC_CHECK_CONNECTDEADENDON, OnBnClickedCheckConnectdeadendon)
	ON_NOTIFY(UDN_DELTAPOS, IDC_SPIN_TRACINGLEVEL, OnDeltaposSpinTracinglevel)
	ON_BN_CLICKED(IDC_BUTTON_CONNECTDEADEND, OnBnClickedButtonConnectdeadend)
	ON_BN_CLICKED(IDC_BUTTON_COMPUTEREGIONBLOCKS, OnBnClickedButtonComputeregionblocks)
	ON_BN_CLICKED(IDC_CHECK_SHOWREGIONBLOCKS, OnBnClickedCheckShowregionblocks)
	ON_BN_CLICKED(IDC_BUTTON_TESTDIFFREGIONDESIGN, OnBnClickedButtonTestdiffregiondesign)
	ON_COMMAND(ID_LOAD_LOADDENSITYMAP, OnLoadLoaddensitymap)
	ON_BN_CLICKED(IDC_BUTTON_LOADPOPDENSITYMAP, OnBnClickedButtonLoadpopdensitymap)
	ON_BN_CLICKED(IDC_CHECK_COMBINEDENSITYMAP, OnBnClickedCheckCombinedensitymap)
	ON_BN_CLICKED(IDC_CHECK_GENTEN_SKETCHBASED, OnBnClickedCheckGentenSketchbased)
	ON_BN_CLICKED(IDC_CHECK_CLOSELOOP, OnBnClickedCheckCloseloop)
	ON_COMMAND(ID_FILE_SAVEAS, OnFileSaveas)
	ON_BN_CLICKED(IDC_BUTTON_EXPORTSTREETNETWORK, OnBnClickedButtonExportstreetnetwork)
	ON_EN_CHANGE(IDC_EDIT_STREETLEVEL, OnEnChangeEditStreetlevel)
	ON_EN_CHANGE(IDC_EDIT_MAJORDENSITY, OnEnChangeEditMajordensity)
	ON_EN_CHANGE(IDC_EDIT_MINORDENSITY, OnEnChangeEditMinordensity)
	ON_EN_CHANGE(IDC_EDIT_MINTENSORLINELENGTH, OnEnChangeEditMintensorlinelength)
	ON_BN_CLICKED(IDC_BUTTON_ADDNOISETOROADS, OnBnClickedButtonAddnoisetoroads)
	ON_BN_CLICKED(IDC_CHECK_USEBOUNDARIESASSKETCHES, OnBnClickedCheckUseboundariesassketches)
	ON_EN_CHANGE(IDC_EDIT_WIDTHBOUNDARY, OnEnChangeEditWidthboundary)
	ON_EN_CHANGE(IDC_EDIT2, OnEnChangeEdit2)
	ON_BN_CLICKED(IDC_CHECK_SHOWINITSEEDS, OnBnClickedCheckShowinitseeds)
	ON_BN_CLICKED(IDC_BUTTON_MAJROADSETTING, OnBnClickedButtonMajroadsetting)
	ON_BN_CLICKED(IDC_BUTTON_SEEDSETTING, OnBnClickedButtonSeedsetting)
	ON_BN_CLICKED(IDC_BUTTON_LOADVEGMAP, OnBnClickedButtonLoadvegmap)
	ON_BN_CLICKED(IDC_CHECK_APPLYASYMFLD, OnBnClickedCheckApplyasymfld)
	ON_BN_CLICKED(IDC_BUTTON_ADDNOISETOFIELD, OnBnClickedButtonAddnoisetofield)
	ON_BN_CLICKED(IDC_BUTTON_NOISESETTING, OnBnClickedButtonNoisesetting)
	ON_BN_CLICKED(IDC_BUTTON_MINROADSETTING, OnBnClickedButtonMinroadsetting)
	ON_BN_CLICKED(IDC_CHECK_REMOVEMAJDEADEDNON, OnBnClickedCheckRemovemajdeadednon)
	ON_BN_CLICKED(IDC_BUTTON_IMPORTSTREETNETWORK, OnBnClickedButtonImportstreetnetwork)
	ON_EN_CHANGE(IDC_EDIT_MINBOUNDARYLEN, OnEnChangeEditMinboundarylen)
	ON_BN_CLICKED(IDC_CHECK_SELREGTOEDITON, OnBnClickedCheckSelregtoediton)
	ON_BN_CLICKED(IDC_BUTTON_REMOVESTREETSINREG, OnBnClickedButtonRemovestreetsinreg)
	ON_BN_CLICKED(IDC_BUTTON_REPLACESTREETSINSIDE, OnBnClickedButtonReplacestreetsinside)
	ON_BN_CLICKED(IDC_BUTTON_COMPSTREETNETINSIDE, OnBnClickedButtonCompstreetnetinside)
	ON_BN_CLICKED(IDC_BUTTON_ADDNOISETOLOCALFIELD, OnBnClickedButtonAddnoisetolocalfield)
	ON_BN_CLICKED(IDC_BUTTON_ADDNOISELOCALROADS, OnBnClickedButtonAddnoiselocalroads)
	ON_BN_CLICKED(IDC_CHECK_USEBRUSHLIKEREGSEL, OnBnClickedCheckUsebrushlikeregsel)
	ON_BN_CLICKED(IDC_RADIO_CONMETHOD1, OnBnClickedRadioConmethod1)
	ON_BN_CLICKED(IDC_RADIO_CONMETHOD2, OnBnClickedRadioConmethod2)
END_MESSAGE_MAP()


// StreetNetPanel message handlers

void StreetNetPanel::OnBnClickedButtonPlacetensorlines()
{
	// TODO: Add your control notification handler code here
	UpdateData(TRUE);
	
	if(upstreaming_edit)
	{
		MessageBox("This functionality is not available for this editing!");
		return;
	}
	//majorDensity = m_edMajTenLineDensity;
	//minorDensity = m_edMinTenLineDensity;
	//mintenline_length = m_edMinTensorLineLength;
	minorDensity = (m_edMinTenLineDensity/map_xrang*(quadmesh->xend-quadmesh->xstart))/
		quadmesh->xinterval;
	majorDensity = (m_edMajTenLineDensity/map_xrang*(quadmesh->xend-quadmesh->xstart))/
		quadmesh->xinterval;
	mintenline_length = (m_edMinTensorLineLength/map_xrang*(quadmesh->xend-quadmesh->xstart))
		/quadmesh->xinterval;

	which_level = m_edStreetLevel;


	reset_roadnetvis();
	init_evenplace_ten();
	
	if(flag_loadmap/*&&sharedvars.GenTenFromExtractedBoundsOn*/)
	{
		if(m_edStreetLevel == 1)
		{
			place_alternative_level1_loadmap();
			majorroadsexisted=true;
		}
		else
			place_tensorlines_consider_bounds();
	}
	else
	{
		if(m_edStreetLevel == 1)
		{
			place_alternative_level1_makeup_seeds();
			majorroadsexisted=true;
		}
		else
			place();
	}
}

void StreetNetPanel::OnBnClickedButtonComputestreetnetwork()
{
	// TODO: Add your control notification handler code here

	if(major==NULL || minor == NULL) return;

	if(upstreaming_edit)
	{
		MessageBox("This functionality is not available for this editing!");
		return;
	}

	init_streetnet();
	//reset_roadnetvis();
	compute_intersects();
	search_for_connection();

	/*  compute the removing dead ends as well 
	1/18/2008
	*/
	remove_all_deadends();
	
	init_streetnet();
	//reset_roadnetvis();
	compute_intersects();
	search_for_connection();


	/*compute the roads for visualization*/
	//init_roadnetvis();
	//compute_visual_points_intersects();
	//roadnetvis->construct_roads_vis(major->evenstreamlines, minor->evenstreamlines);

		g_mydlg2->m_chkShowRoadMapOn=TRUE;
		g_mydlg2->m_chkShowStreetUseNetworkOn=FALSE;
		g_mydlg2->m_chkShowIBFVOn=FALSE;

		g_mydlg1->m_chkBrushInterfaceOn=FALSE;
		sharedvars.BrushInterfaceOn=false;

		g_mydlg1->m_chkEnableSketchBasedDesign=FALSE;
		sharedvars.EnableSketchBasedDesign=false;
		g_mydlg1->m_chkShowSketchesOn=FALSE;
		sharedvars.ShowSketchesOn=false;

		sharedvars.ShowRoadMapOn=true;
		g_pclGlView->displayRoadNetOn=true;
		g_pclGlView->ShapeControlPtsOn=0;
		sharedvars.ShowIBFVOn=false;
		sharedvars.ShowPopDensityMapOn=false;
		sharedvars.ShowStreetUseNetworkOn=false;

		/*  Turn off the tensor line visualization here ! */
		g_mydlg2->m_chkShowTensorLinesOn=FALSE;
		sharedvars.ShowTensorLinesOn=false;
		g_pclGlView->showTensorLineOn=false;

		/*  could be a possible bug here! 1/19/2008 */

	
}

void StreetNetPanel::OnBnClickedCheckEditstreetgraphon()
{
	// TODO: Add your control notification handler code here
	if(m_chkEditStreetNetOn==FALSE)
	{
		m_chkEditStreetNetOn=TRUE;
		g_pclGlView->streetNetEditOn = true;
		g_pclGlView->TrajectoryOn = 0;
		sharedvars.EditStreetNetOn=true;

		g_mydlg1->m_chkEditElementOn = FALSE;
		g_mydlg1->m_chkMoveElemOn = FALSE;
		g_mydlg1->m_chkRemoveElemOn = FALSE;
		g_mydlg1->m_chkTensorDesignOn = FALSE;
		g_mydlg1->m_chkRegionSmoothOn=FALSE;
		g_mydlg1->m_chkBrushInterfaceOn=FALSE;
		g_mydlg1->m_chkEnableHeighfieldDesignOn=FALSE;
		g_streeteditpanel->m_chkSelStreetRegToEditOn=FALSE;

		sharedvars.EditElementOn = false;
		sharedvars.RemoveElemOn = false;
		sharedvars.MoveElemOn = false;
		sharedvars.TensorDesignOn = false;
		sharedvars.tenElemType = -1;
		sharedvars.RegionSmoothOn=false;
		sharedvars.BrushInterfaceOn=false;
		sharedvars.EnableSketchBasedDesign=false;
		sharedvars.EnableHeighfieldDesignOn=false;
		sharedvars.SelStreetRegToEditOn=false;
	}
	else
	{
		m_chkEditStreetNetOn=FALSE;
		g_pclGlView->streetNetEditOn=false;
		sharedvars.EditStreetNetOn=false;
	}
	UpdateData(FALSE);
}

void StreetNetPanel::OnBnClickedCheckMeanfilteron()
{
	// TODO: Add your control notification handler code here
	if(m_chkMeanFilterOn==FALSE)
	{
		m_chkMeanFilterOn=TRUE;
		sharedvars.MeanFilterOn=true;
	}
	else
	{
		m_chkMeanFilterOn=FALSE;
		sharedvars.MeanFilterOn=false;
	}
}

void StreetNetPanel::OnBnClickedButtonLoadawatermap()
{
	// TODO: Add your control notification handler code here
	CFile f;
	int i, j;

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

		if(sharedvars.MeanFilterOn)
			imgprocess->mean_filter(m_edMeanFilterSize/*3*/);
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

		for(i=0;i<3*512*512;i++)
		{
			displaymap[i]=imgboundaries->output[i];
			streetmapbackground[i]=imgboundaries->streetmap_background[i];
		}

		//for(i=0;i<512;i++)
		//{
		//	for(j=0;j<512;j++)
		//	{
		//		fittedmap1[3*(i*512+j)]  =imgprocess->imgs[i][j][0];
		//		fittedmap1[3*(i*512+j)+1]=imgprocess->imgs[i][j][1];
		//		fittedmap1[3*(i*512+j)+2]=imgprocess->imgs[i][j][2];
		//	}
		//}
		
		
		if(flag_loadvegmap)
		{

			/* set the vegetation flag for each vertex of the mesh */
			set_vegflags_verts(quadmesh->xstart, quadmesh->xend, 
				quadmesh->ystart, quadmesh->yend, vegmap_fit, 512,512);

			/*  we need to update the background image as well  */
			for(i=0;i<512*512;i++)
			{
				if((vegmap_fit[3*i]+vegmap_fit[3*i+1]+vegmap[3*i+2])>200)
				{
					streetmapbackground[3*i]  =vegmap_fit[3*i];
					streetmapbackground[3*i+1]=vegmap_fit[3*i+1];
					streetmapbackground[3*i+2]=vegmap_fit[3*i+2];
				}
			}
		}



		///*get the geometry boundary*/
		mapboundarylist=new MapBoundaryList();
		mapboundarylist->get_worldcoords_unitintervals(0.,1.,0.,1.,512,512);

		mapboundarylist->downsamp_imgboundaries(imgboundaries, 8);
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

void StreetNetPanel::OnBnClickedCheckGentensorbasedonextractboundaries()
{
	// TODO: Add your control notification handler code here
	if(m_chkGenTenFromExtractedBoundsOn==FALSE)
	{
		m_chkGenTenFromExtractedBoundsOn=TRUE;
		sharedvars.GenTenFromExtractedBoundsOn=true;
	}
	else
	{
		m_chkGenTenFromExtractedBoundsOn=FALSE;
		sharedvars.GenTenFromExtractedBoundsOn=false;
	}
	UpdateData(FALSE);
}

void StreetNetPanel::OnBnClickedCheckUseallboundaries()
{
	// TODO: Add your control notification handler code here
	if(m_chkUseAllBoundsOn==FALSE)
	{
		m_chkUseAllBoundsOn=TRUE;
		sharedvars.UseAllBoundsOn=true;
	}
	else
	{
		m_chkUseAllBoundsOn=FALSE;
		sharedvars.UseAllBoundsOn=false;
	}
}

BOOL StreetNetPanel::OnInitDialog()
{
	CDialog::OnInitDialog();

	// TODO:  Add extra initialization here
	m_ctrlSliderBoundaryRegionWidth.SetRange(1, 300);
	m_ctrlSliderBoundaryRegionWidth.SetPos(40);
	m_ctrlSliderBoundaryRegionWidth.SetTicFreq(1);
	m_ctrlSliderBoundaryRegionWidth.SetTic(0);

	m_edEffectWidthBoundary=(double)m_ctrlSliderBoundaryRegionWidth.GetPos();

	UpdateData(FALSE);

	m_ctrlSpinTracingLevel.SetRange(1,3);
	m_ctrlSpinTracingLevel.SetPos(1);

	//m_ctrlMajRoadSettingButton.EnableWindow(FALSE);
	m_ctrlMajRoadSettingButton.EnableWindow(TRUE);

	maj_road_settingDlg->Create(maj_road_settingDlgID, GetParent());

	seeds_ctrlDlg->Create(seeds_ctrlDlgID, GetParent());

	noise_ctrlDlg->Create(noise_ctrlDlgID, GetParent());

	min_road_setDlg->Create(min_road_setDlgID, GetParent());

	return TRUE;  // return TRUE unless you set the focus to a control
	// EXCEPTION: OCX Property Pages should return FALSE
}

void StreetNetPanel::OnNMCustomdrawSliderBoundarywidth(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);

	// TODO: Add your control notification handler code here
	BoundRegionWidth=(double)m_ctrlSliderBoundaryRegionWidth.GetPos();
	m_edEffectWidthBoundary=BoundRegionWidth;
	UpdateData(FALSE);

	*pResult = 0;
}

void StreetNetPanel::OnBnClickedCheckGridon()
{
	// TODO: Add your control notification handler code here
	if(m_chkDesignGridOn==FALSE)
	{
		m_chkDesignGridOn=TRUE;
		sharedvars.DesignGridOn=true;
		UpdateData(TRUE);
		//designgrid_xinterval=;
		//designgrid_yinterval;
		majorDensity = m_edMajTenLineDensity;
		minorDensity = m_edMinTenLineDensity;

	}
	else
	{
		m_chkDesignGridOn=FALSE;
		sharedvars.DesignGridOn=false;
	}
		UpdateData(FALSE);
}

void StreetNetPanel::OnBnClickedCheckConnectdeadendon()
{
	// TODO: Add your control notification handler code here
	if(m_chkConnectDeadEndsOn==FALSE)
	{
		m_chkConnectDeadEndsOn=TRUE;
		sharedvars.ConnectDeadEndsOn=true;
	}
	else
	{
		m_chkConnectDeadEndsOn=FALSE;
		sharedvars.ConnectDeadEndsOn=false;
	}
		
	UpdateData(FALSE);
}

void StreetNetPanel::OnDeltaposSpinTracinglevel(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMUPDOWN pNMUpDown = reinterpret_cast<LPNMUPDOWN>(pNMHDR);
	// TODO: Add your control notification handler code here

	m_edStreetLevel=m_ctrlSpinTracingLevel.GetPos();
	UpdateData(FALSE);

	*pResult = 0;
}

void StreetNetPanel::OnBnClickedButtonConnectdeadend()
{
	// TODO: Add your control notification handler code here
	if(upstreaming_edit)
	{
		MessageBox("This functionality is not available for this editing!");
		return;
	}

	remove_all_deadends();
}

void StreetNetPanel::OnBnClickedButtonComputeregionblocks()
{
	// TODO: Add your control notification handler code here
	construct_regionblocks_edgewise();
}

void StreetNetPanel::OnBnClickedCheckShowregionblocks()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowRegionBlocksOn==FALSE)
	{
		m_chkShowRegionBlocksOn=TRUE;
		sharedvars.ShowRegionBlocksOn=true;
	}
	else
	{
		m_chkShowRegionBlocksOn=FALSE;
		sharedvars.ShowRegionBlocksOn=false;
	}
}


void StreetNetPanel::OnBnClickedButtonTestdiffregiondesign()
{
	// TODO: Add your control notification handler code here

	//load_the_test_img("test_diff_regions.bmp");
	//g_mydlg1->m_ctrlcomboxDesignRegionIndex.AddString("1");


	/* WE use this button to test the generation of sketched based design*/

	FILE *fp;
	//fp=fopen("sketchbased_error.txt", "w");
	//fprintf(fp, "initializing skecthes\n");
	//fclose(fp);

	release_domain_boundaries();
	init_domain_boundaries();

	init_sketchnet();
	init_sketchlineinfo_incells();
	update_sketchcurveinfo();
	cal_sketchlines_intersects();
	search_for_sketchbased_graph();
	
	//fp=fopen("sketchbased_error.txt", "w");
	//fprintf(fp, "start extending the skecthes\n");
	//fclose(fp);

	extend_sketchcurves();
	init_sketchnet();
	init_sketchlineinfo_incells();
	update_sketchcurveinfo();
	cal_sketchlines_intersects();
	search_for_sketchbased_graph();

	//fp=fopen("sketchbased_error.txt", "w");
	//fprintf(fp, "start extracting the regions determined by the sketches\n");
	//fclose(fp);

	/*      extract the blocks        */
	init_sketchblocklist();
	construct_regionblocks_edgewise(sketchnet, sketchblocklist);
	
	//fp=fopen("sketchbased_error.txt", "w");
	//fprintf(fp, "start marking the subregions\n");
	//fclose(fp);

	/*      mark the vertices of the mesh       */
extern void color_subregionblocks(RegionBlockList *aregionblocklist, StreetNet *net);
	color_subregionblocks(sketchblocklist, sketchnet);

	highwayexisted=true;
	
	//fp=fopen("sketchbased_error.txt", "w");
	//fprintf(fp, "start generating the tensor field\n");
	//fclose(fp);

	if(sharedvars.GenTenFromSketchesOn)
	{
		gen_sketch_based_tensorfield(BoundRegionWidth);

		//fp=fopen("sketchbased_error.txt", "w");
		//fprintf(fp, "finish generating the tensor field\n");
		//fclose(fp);

		render_alpha_map_quad(false);
		render_alpha_map_quad(true);
		
		//fp=fopen("sketchbased_error.txt", "w");
		//fprintf(fp, "finish rendering the alpha map.\n");
		//fclose(fp);

	}

	/*  Disable the sketch based design  */
	/*sharedvars.EnableSketchBasedDesign=false*/;
}



/*load a population density map here */
void StreetNetPanel::OnLoadLoaddensitymap()
{
	// TODO: Add your command handler code here
	OnBnClickedButtonLoadpopdensitymap();
}

void StreetNetPanel::OnBnClickedButtonLoadpopdensitymap()
{
	// TODO: Add your control notification handler code here
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
		
		popdensitymap = BmpToTexture(filename1, &popdensitymap_w, &popdensitymap_h);	

		//if(popdensitymap_fit!=NULL)
		//	free(popdensitymap_fit);

		//popdensitymap_fit=(unsigned char*)malloc(sizeof(unsigned char)*3*512*512);

		/*obtain an image that fit in the 512x512 windows*/
		get_fitted_map(popdensitymap, popdensitymap_w, popdensitymap_h, 
			popdensitymap_fit, 512, 512);

		get_density_value(quadmesh->xstart, quadmesh->xend, quadmesh->ystart, quadmesh->yend,
			popdensitymap, popdensitymap_w, popdensitymap_h);


		/**                                       **/
		/**          Testing code here            **/
		
		/*displaymap=*/popdensitymap_disp=popdensitymap_fit;

		/**                                       **/

		//flag_loadmap = true;
		//sharedvars.ShowTheMapOn=true;


		//init_degpts();
		//init_ten_designelems();
		//tensor_init_tex();

		sharedvars.ShowPopDensityMapOn=true;
		g_mydlg2->m_chkShowPopDensityMapOn=TRUE;
		
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

void StreetNetPanel::OnBnClickedCheckCombinedensitymap()
{
	// TODO: Add your control notification handler code here
	if(m_chkCombinePopDensityOn==FALSE)
	{
		m_chkCombinePopDensityOn=TRUE;
		sharedvars.CombinePopDensityOn=true;
	}
	else
	{
		m_chkCombinePopDensityOn=FALSE;
		sharedvars.CombinePopDensityOn=false;
	}
	UpdateData(FALSE);
}

void StreetNetPanel::OnBnClickedCheckGentenSketchbased()
{
	// TODO: Add your control notification handler code here
	if(m_chkGenTenFromSketchesOn==FALSE)
	{
		m_chkGenTenFromSketchesOn=TRUE;
		sharedvars.GenTenFromSketchesOn=true;
	}
	else
	{
		m_chkGenTenFromSketchesOn=FALSE;
		sharedvars.GenTenFromSketchesOn=false;
	}
	UpdateData(FALSE);
}

void StreetNetPanel::OnBnClickedCheckCloseloop()
{
	// TODO: Add your control notification handler code here
	if(m_chkCloseLoopOn==FALSE)
	{
		m_chkCloseLoopOn=TRUE;
		sharedvars.CloseLoopOn=true;
	}
	else
	{
		m_chkCloseLoopOn=FALSE;
		sharedvars.CloseLoopOn=false;
	}
	UpdateData(FALSE);
}


/*export the obtained street network */
void StreetNetPanel::OnFileSaveas()
{
	// TODO: Add your command handler code here
	OnBnClickedButtonExportstreetnetwork();
}

void StreetNetPanel::OnBnClickedButtonExportstreetnetwork()
{
	// TODO: Add your control notification handler code here
	CFile f;
	char filename[256];

	char strFilter[] = { "Text Files (*.gra)|*.gra|All Files (*.*)|*.*||" };

	CFileDialog FileDlg(FALSE, ".gra", NULL, 0, strFilter);
	if( FileDlg.DoModal() == IDOK )
	{

		CString f2 = FileDlg.GetFileName();
		const char *tempf = f2;
		strcpy(filename, f2);
	
		save_cur_street_network(filename);
	}
}

void StreetNetPanel::OnEnChangeEditStreetlevel()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	//which_level = m_edStreetLevel;
	UpdateData();
	which_level = m_edStreetLevel;

	if(which_level<1 || which_level>2)
	{
		MessageBox("Please set the tracing level between 1 and 2", MB_OK);
		which_level = m_edStreetLevel = 2;
		UpdateData(FALSE);
		return;
	}

	if(m_edStreetLevel==1)
	{
	    m_edMinTenLineDensity *=10;
		m_edMajTenLineDensity *=10;
	}
	else if(m_edStreetLevel==2)
	{
		if(m_edMinTenLineDensity>20)
			m_edMinTenLineDensity /=10;
		if(m_edMajTenLineDensity>20)
			m_edMajTenLineDensity /=10;
	}

	if(which_level==1)
		m_ctrlMajRoadSettingButton.EnableWindow(TRUE);
	else
		m_ctrlMajRoadSettingButton.EnableWindow(FALSE);

	UpdateData(FALSE);
}

void StreetNetPanel::OnEnChangeEditMajordensity()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData();
	majorDensity = (m_edMajTenLineDensity/map_xrang*(quadmesh->xend-quadmesh->xstart))/
		quadmesh->xinterval;
}

void StreetNetPanel::OnEnChangeEditMinordensity()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData();
	//minorDensity = m_edMinTenLineDensity;
	minorDensity = (m_edMinTenLineDensity/map_xrang*(quadmesh->xend-quadmesh->xstart))/
		quadmesh->xinterval;
}

void StreetNetPanel::OnEnChangeEditMintensorlinelength()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData();
	//mintenline_length = m_edMinTensorLineLength;
	mintenline_length = (m_edMinTensorLineLength/map_xrang*(quadmesh->xend-quadmesh->xstart))
		/quadmesh->xinterval;
}

extern void jitter_all_roads(double scale, double curAmp);
extern void jitter_all_tens(double scale, double curAmp);
extern double NoiseFreq;
extern double NoiseAmp;


void StreetNetPanel::OnBnClickedButtonAddnoisetoroads()
{
	// TODO: Add your control notification handler code here
	/*   Add Perlin noise to the obtained roads  */
	//jitter_all_roads(500., 0.1);
	if(which_level==1 && majorroadsexisted)
	{
		jitter_all_majRoads(NoiseFreq, NoiseAmp, false);
		jitter_all_majRoads(NoiseFreq, NoiseAmp, true);
	}
	else
		jitter_all_roads(NoiseFreq, NoiseAmp);
}

void StreetNetPanel::OnBnClickedButtonAddnoisetofield()
{
	// TODO: Add your control notification handler code here
	//jitter_all_tens(50., 1);
	jitter_all_tens(NoiseFreq, NoiseAmp);
	cal_all_eigenvecs_quad();

	init_degpts();
	/*render alpha map*/
	render_alpha_map_quad(false);
	render_alpha_map_quad(true);

	locate_degpts_cells_tranvec_quad();
}


void StreetNetPanel::OnBnClickedCheckUseboundariesassketches()
{
	// TODO: Add your control notification handler code here
	if(m_chkUseBoundsAsSketchesOn==FALSE)
	{
		m_chkUseBoundsAsSketchesOn=TRUE;
		sharedvars.UseBoundsAsSketchesOn=true;
	}
	else
	{
		m_chkUseBoundsAsSketchesOn=FALSE;
		sharedvars.UseBoundsAsSketchesOn=false;
	}
	UpdateData(FALSE);
}

void StreetNetPanel::OnEnChangeEditWidthboundary()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData();
	BoundRegionWidth=m_edEffectWidthBoundary;
		
	m_ctrlSliderBoundaryRegionWidth.SetPos(floor(BoundRegionWidth));
}

void StreetNetPanel::OnEnChangeEdit2()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData();
	MeanFilterSize=m_edMeanFilterSize;
}

void StreetNetPanel::OnBnClickedCheckShowinitseeds()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowInitSeedsOn==FALSE)
	{
		m_chkShowInitSeedsOn=TRUE;
		sharedvars.ShowInitSeedsOn=true;
	}
	else
	{
		m_chkShowInitSeedsOn=FALSE;
		sharedvars.ShowInitSeedsOn=false;
	}
	UpdateData(FALSE);
}

void StreetNetPanel::OnBnClickedButtonMajroadsetting()
{
	//// TODO: Add your control notification handler code here
	if(maj_road_settingDlg == NULL)
		return;

	CRect l_rectClient;
	CRect l_rectWnd;

	GetClientRect(l_rectClient);
	GetWindowRect(l_rectWnd);
	GetParent()->ScreenToClient(l_rectWnd);

	maj_road_settingDlg->GetWindowRect(l_rectClient);

	maj_road_settingDlg->SetWindowPos(&wndTop,l_rectWnd.left+50,l_rectWnd.top+100,l_rectClient.Width(),
		l_rectClient.Height(),SWP_SHOWWINDOW);
	//maj_road_settingDlg->SetFocus();
	maj_road_settingDlg->ShowWindow(SW_SHOW);

}


/*  seed point control panel */
void StreetNetPanel::OnBnClickedButtonSeedsetting()
{
	// TODO: Add your control notification handler code here
	if(seeds_ctrlDlg == NULL)
		return;

	CRect l_rectClient;
	CRect l_rectWnd;

	GetClientRect(l_rectClient);
	GetWindowRect(l_rectWnd);
	GetParent()->ScreenToClient(l_rectWnd);

	seeds_ctrlDlg->GetWindowRect(l_rectClient);

	seeds_ctrlDlg->SetWindowPos(&wndTop,l_rectWnd.left+50,l_rectWnd.top+300,l_rectClient.Width(),
		l_rectClient.Height(),SWP_SHOWWINDOW);
	seeds_ctrlDlg->ShowWindow(SW_SHOW);
}

void StreetNetPanel::OnBnClickedButtonLoadvegmap()
{
	// TODO: Add your control notification handler code here
	// TODO: Add your control notification handler code here
	CFile f;
	int i;

	char strFilter[] = { "Bitmap Files (*.bmp)|*.bmp|All Files (*.*)|*.*||" };

	CFileDialog FileDlg(TRUE, ".bmp", NULL, 0, strFilter);
	if( FileDlg.DoModal() == IDOK )
	{
		CString f2 = FileDlg.GetFileName();
		CString f3 = f2.Left(f2.Find('.'))+".bmp";

		const char *filename = f3;
		char filename1[255];
		strcpy(filename1, filename);
		
		vegmap = BmpToTexture(filename1, &vegmap_w, &vegmap_h);	

		/*obtain an image that fit in the 512x512 windows*/
		get_fitted_map(vegmap, vegmap_w, vegmap_h, 
			vegmap_fit, 512, 512);

		/* set the vegetation flag for each vertex of the mesh */
		set_vegflags_verts(quadmesh->xstart, quadmesh->xend, 
			quadmesh->ystart, quadmesh->yend, vegmap_fit, 512,512);
	
		vegmap_disp=vegmap_fit;

		g_mydlg2->m_chkShowVegMapOn=TRUE;
		sharedvars.ShowVegMapOn=true;

		sharedvars.ShowPopDensityMapOn=false;
		g_mydlg2->m_chkShowPopDensityMapOn=FALSE;

		g_mydlg2->m_chkShowIBFVOn=FALSE;
		g_mydlg2->m_chkShowScalarFieldOn=FALSE;
		g_mydlg2->m_chkShowExtractedBoundsOn=FALSE;
		g_mydlg2->m_chkShowRoadMapOn=FALSE;

		sharedvars.ShowIBFVOn=false;
		sharedvars.ShowScalarFieldOn=false;
		sharedvars.ShowExtractedBoundsOn=false;
		sharedvars.ShowRoadMapOn=false;

		g_mydlg2->UpdateData(FALSE);

		flag_loadvegmap=true;

		/*  we need to update the background image as well  */
		for(i=0;i<512*512;i++)
		{
			if((vegmap_fit[3*i]+vegmap_fit[3*i+1]+vegmap[3*i+2])>200)
			{
				streetmapbackground[3*i]  =vegmap_fit[3*i];
				streetmapbackground[3*i+1]=vegmap_fit[3*i+1];
				streetmapbackground[3*i+2]=vegmap_fit[3*i+2];
			}
		}

	}
}

void StreetNetPanel::OnBnClickedCheckApplyasymfld()
{
	// TODO: Add your control notification handler code here
	if(m_chkApplyAsymFldOn==FALSE)
	{
		m_chkApplyAsymFldOn=TRUE;
		sharedvars.ApplyAsymFldOn=true;
	}
	else
	{
		m_chkApplyAsymFldOn=FALSE;
		sharedvars.ApplyAsymFldOn=false;
	}
	UpdateData(FALSE);
}


void StreetNetPanel::OnBnClickedButtonNoisesetting()
{
	// TODO: Add your control notification handler code here
	if(noise_ctrlDlg == NULL)
		return;

	CRect l_rectClient;
	CRect l_rectWnd;

	GetClientRect(l_rectClient);
	GetWindowRect(l_rectWnd);
	GetParent()->ScreenToClient(l_rectWnd);

	noise_ctrlDlg->GetWindowRect(l_rectClient);

	noise_ctrlDlg->SetWindowPos(&wndTop,l_rectWnd.left+50,l_rectWnd.top+500,l_rectClient.Width(),
		l_rectClient.Height(),SWP_SHOWWINDOW);
	noise_ctrlDlg->ShowWindow(SW_SHOW);
}

void StreetNetPanel::OnBnClickedButtonMinroadsetting()
{
	// TODO: Add your control notification handler code here
	if(maj_road_settingDlg == NULL)
		return;

	CRect l_rectClient;
	CRect l_rectWnd;

	GetClientRect(l_rectClient);
	GetWindowRect(l_rectWnd);
	GetParent()->ScreenToClient(l_rectWnd);

	min_road_setDlg->GetWindowRect(l_rectClient);

	min_road_setDlg->SetWindowPos(&wndTop,l_rectWnd.left+50,l_rectWnd.top+120,l_rectClient.Width(),
		l_rectClient.Height(),SWP_SHOWWINDOW);
	//maj_road_settingDlg->SetFocus();
	min_road_setDlg->ShowWindow(SW_SHOW);
}

void StreetNetPanel::OnBnClickedCheckRemovemajdeadednon()
{
	// TODO: Add your control notification handler code here
	if(m_chkRemoveMajDeadEndsOn==FALSE)
	{
		m_chkRemoveMajDeadEndsOn=TRUE;
		sharedvars.RemoveMajDeadEndsOn=true;
	}
	else
	{
		m_chkRemoveMajDeadEndsOn=FALSE;
		sharedvars.RemoveMajDeadEndsOn=false;
	}
	UpdateData(FALSE);
}

void StreetNetPanel::OnBnClickedButtonImportstreetnetwork()
{
	// TODO: Add your control notification handler code here
	CFile f;
	char filename[256];

	char strFilter[] = { "Text Files (*.gra)|*.gra" };

	CFileDialog FileDlg(TRUE, ".gra", NULL, 0, strFilter);
	if( FileDlg.DoModal() == IDOK )
	{

		CString f2 = FileDlg.GetFileName();
		const char *tempf = f2;
		strcpy(filename, f2);
	
		if(!load_a_street_network(filename))
		{
			MessageBox("Can not open the file!");
			return;
		}

		g_pclGlView->displayRoadNetOn=true;
		sharedvars.ShowRoadMapOn=true;
		g_mydlg2->m_chkShowRoadMapOn=TRUE;

		upstreaming_edit=true;
	}
}

void StreetNetPanel::OnEnChangeEditMinboundarylen()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData();
	minBoundaryLen=m_edMinBoundaryLen;
	UpdateData(FALSE);
}

void StreetNetPanel::OnBnClickedCheckSelregtoediton()
{
	// TODO: Add your control notification handler code here
	if(m_chkSelRegToEditOn==FALSE)
	{
		m_chkSelRegToEditOn = TRUE;
		sharedvars.SelStreetRegToEditOn=true;

		please_comb_prefield=true;

		g_mydlg1->m_chkRegionSmoothOn=TRUE;
		g_mydlg1->m_chkBrushInterfaceOn=FALSE;
		g_mydlg1->m_chkRemoveElemOn = FALSE;
		g_mydlg1->m_chkEditElementOn = FALSE;
		g_mydlg1->m_chkMoveElemOn = FALSE;
		g_mydlg1->m_chkTensorDesignOn = FALSE;
		g_mydlg1->m_chkEnableHeighfieldDesignOn = FALSE;

		sharedvars.RegionSmoothOn=true;
		sharedvars.BrushInterfaceOn = false;
		sharedvars.RemoveElemOn = false;
		sharedvars.EditElementOn = false;
		sharedvars.MoveElemOn = false;
		sharedvars.TensorDesignOn = false;
		sharedvars.EnableSketchBasedDesign=false;
		sharedvars.EnableHeighfieldDesignOn=false;

		g_mydlg1->m_rdAddVFSingularElem = -1;  
		g_pclGlView->SmoothOn = 1;              ////SELECTION MODE
		g_pclGlView->PickPointOn = 1;
		g_pclGlView->DisplaySmoothRegionOn = 1;
		Num_SmoothRegionpoints = 0;
		g_pclGlView->EditModeOn = 0;
		g_pclGlView->MoveElemOn = 0;
		g_pclGlView->ShapeControlPtsOn = 0;
		g_pclGlView->deleteElemOn=false;
	}
	else
	{
		m_chkSelRegToEditOn = FALSE;
		sharedvars.SelStreetRegToEditOn=false;
		g_pclGlView->SmoothOn = 0;              ////SELECTION MODE
		g_pclGlView->PickPointOn = 0;
		g_pclGlView->DisplaySmoothRegionOn = 0;
		sharedvars.RegionSmoothOn=false;
		please_comb_prefield=false;

		is_on_local_editing=false;

	}
	UpdateData(FALSE);
}

void StreetNetPanel::OnBnClickedButtonRemovestreetsinreg()
{
	// TODO: Add your control notification handler code here
	if(sharedvars.BrushInterfaceOn)
	{
		/*  we need to determine the region according to the level set distance  */
		compute_level_set_contour(BoundRegionWidth-10*quadmesh->xinterval);
	}
	else if(sharedvars.SelStreetRegToEditOn||sharedvars.BrushLikeRegSelOn)
	{
		remove_street_in_reg();
		mark_inner_verts_with_new_RegIndex();
		update_RegIndex_inner_and_boundary_cells();

		if(sharedvars.BrushLikeRegSelOn)
		{
			/*  generate a tensor field follow the design curve inside the region  */
			get_brush_approDir();

			cal_tensorvals_quad_inReg(cur_max_reg_index);
			cal_eigenvecs_quad_inReg(cur_max_reg_index);

			init_degpts();
			/*render alpha map*/
			render_alpha_map_quad(false);
			render_alpha_map_quad(true);

			locate_degpts_cells_tranvec_quad();
		}
	}
}

void StreetNetPanel::OnBnClickedButtonReplacestreetsinside()
{
	// TODO: Add your control notification handler code here
	g_streetnetpanel->UpdateData();

	minorDensity = (m_edMinTenLineDensity/map_xrang*(quadmesh->xend-quadmesh->xstart))/
		quadmesh->xinterval;
	majorDensity = (m_edMajTenLineDensity/map_xrang*(quadmesh->xend-quadmesh->xstart))/
		quadmesh->xinterval;
	mintenline_length = (m_edMinTensorLineLength/map_xrang*(quadmesh->xend-quadmesh->xstart))
		/quadmesh->xinterval;
	replace_inReg();
}

void StreetNetPanel::OnBnClickedButtonCompstreetnetinside()
{
	// TODO: Add your control notification handler code here
	construct_sub_graph_inReg();

	/* Add a testing  */
	int i;
	extern StreetNet *streetnet;

	/*  merge subgraph to the original street network still has problem 1/19/2008 */

	//if(m_rdSubConnectMethod==0)
	//{
		merge_subgraph_to_streetnet_3(false,false);
		merge_subgraph_to_streetnet_3(false,true);
		merge_subgraph_to_streetnet_3(true,false);
		merge_subgraph_to_streetnet_3(true,true);
	//}

	//else
	//{
	//	merge_subgraph_to_streetnet_4(false,false);
	//	merge_subgraph_to_streetnet_4(false,true);
	//	merge_subgraph_to_streetnet_4(true,false);
	//	merge_subgraph_to_streetnet_4(true,true);
	//}

	/**/
	connect_outer_deadends_Reg();

	for(i=0;i<streetnet->nodelist->nelems;i++)
	{
		if(streetnet->nodelist->intersects[i]->endpt)
			continue;

		if(streetnet->nodelist->intersects[i]->nadjedges==1)
		{
			streetnet->nodelist->intersects[i]->endpt=true;
		}
	}

	for(i=0;i<streetnet->edgelist->nedges;i++)
	{
		if(streetnet->edgelist->edges[i]->cancel)
			continue;
		if(streetnet->nodelist->intersects[streetnet->edgelist->edges[i]->node_index1]->deleted
			||streetnet->nodelist->intersects[streetnet->edgelist->edges[i]->node_index2]->deleted)
		{
			//int test=0;
			streetnet->edgelist->edges[i]->cancel=true;
		}
	}

	/*  we also need to merge the outer dead end to the new subgraph!  
         1/19/2008
	*/
	update_street_network();

	/**/
	//int i;
	//for(i=0;i<streetnet->edgelist->nedges;i++)
	//{
	//	if(streetnet->edgelist->edges[i]->node_index1>=streetnet->nodelist->nelems
	//		||streetnet->edgelist->edges[i]->node_index2>=streetnet->nodelist->nelems)
	//	{
	//		int test=0;
	//	}
	//}
}

void StreetNetPanel::OnBnClickedButtonAddnoisetolocalfield()
{
	// TODO: Add your control notification handler code here
	init_quad_regionsmooth();
	find_boundarycells();
	mark_boundVerts();
	find_innerVerts();
	find_inner_cells();

	jitter_tens_inReg(NoiseFreq, NoiseAmp);
	cal_all_eigenvecs_quad();

	init_degpts();
	/*render alpha map*/
	render_alpha_map_quad(false);
	render_alpha_map_quad(true);

	locate_degpts_cells_tranvec_quad();
}

void StreetNetPanel::OnBnClickedButtonAddnoiselocalroads()
{
	// TODO: Add your control notification handler code here
	init_quad_regionsmooth();
	find_boundarycells();
	mark_boundVerts();
	find_innerVerts();
	find_inner_cells();

	jitter_roads_inReg(NoiseFreq/*500.*/);
}



void StreetNetPanel::OnBnClickedCheckUsebrushlikeregsel()
{
	// TODO: Add your control notification handler code here
	if(m_chkBrushLikeRegSelOn==FALSE)
	{
		m_chkBrushLikeRegSelOn=TRUE;
		sharedvars.BrushLikeRegSelOn=true;

		m_chkSelRegToEditOn = FALSE;
		sharedvars.SelStreetRegToEditOn=false;

		please_comb_prefield=true;

		g_mydlg1->m_chkRegionSmoothOn=TRUE;
		g_mydlg1->m_chkBrushInterfaceOn=FALSE;
		g_mydlg1->m_chkRemoveElemOn = FALSE;
		g_mydlg1->m_chkEditElementOn = FALSE;
		g_mydlg1->m_chkMoveElemOn = FALSE;
		g_mydlg1->m_chkTensorDesignOn = FALSE;
		g_mydlg1->m_chkEnableHeighfieldDesignOn = FALSE;

		sharedvars.RegionSmoothOn=true;
		sharedvars.BrushInterfaceOn = false;
		sharedvars.RemoveElemOn = false;
		sharedvars.EditElementOn = false;
		sharedvars.MoveElemOn = false;
		sharedvars.TensorDesignOn = false;
		sharedvars.EnableSketchBasedDesign=false;
		sharedvars.EnableHeighfieldDesignOn=false;

		g_mydlg1->m_rdAddVFSingularElem = -1;  
		//g_pclGlView->SmoothOn = 1;              ////SELECTION MODE
		//g_pclGlView->PickPointOn = 1;
		g_pclGlView->DisplaySmoothRegionOn = 1;
		Num_SmoothRegionpoints = 0;
		g_pclGlView->EditModeOn = 0;
		g_pclGlView->MoveElemOn = 0;
		g_pclGlView->ShapeControlPtsOn = 1;
		g_pclGlView->deleteElemOn=false;
		num_shapecontrol_pts = 0;
	    num_curvepts_output = 0;
		/*initialize the brush point list*/
		if(brushpts == NULL)
		{
			brushpts=(ctr_point**)malloc(sizeof(ctr_point*)*curMaxNumBrushPts);
			if(brushpts==NULL)
				exit(-1);
			for(int i=0; i<curMaxNumBrushPts; i++)
				brushpts[i]=(ctr_point*)malloc(sizeof(ctr_point));
		}
		nbrushpts = 0;
	}
	else
	{
		m_chkBrushLikeRegSelOn=FALSE;
		sharedvars.BrushLikeRegSelOn=false;
		g_pclGlView->SmoothOn = 0;              ////SELECTION MODE
		g_pclGlView->PickPointOn = 0;
		g_pclGlView->DisplaySmoothRegionOn = 0;
		sharedvars.RegionSmoothOn=false;
		please_comb_prefield=false;

		is_on_local_editing=false;
	}
	UpdateData(FALSE);
}

void StreetNetPanel::OnBnClickedRadioConmethod1()
{
	// TODO: Add your control notification handler code here
	m_rdSubConnectMethod=0;
	UpdateData(FALSE);
}

void StreetNetPanel::OnBnClickedRadioConmethod2()
{
	// TODO: Add your control notification handler code here
	m_rdSubConnectMethod=1;
	UpdateData(FALSE);
}
