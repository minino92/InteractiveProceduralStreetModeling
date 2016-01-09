// MyDlg1.cpp : implementation file
//

#include "stdafx.h"
#include "IBFV2D.h"
#include "MyDlg1.h"
#include ".\mydlg1.h"
#include "shareinterfacevars.h"

#include "VFDataStructure.h"

#include "GlView.h"
#include ".\glview.h"


#include "tensorvis.h"
#include "caldeformation.h"
#include "tensordesign.h"
#include "tensoranalysis.h"
#include "computeroadvis.h"
#include "evenlystreamlines.h"
#include "regionsmooth_quad.h"

#include "ImgBoundaryExtract.h"

#include "SketchDesign.h"

extern SharedInterfaceVars sharedvars;
extern CGlView *g_pclGlView;

extern int *boundarycells;
extern int nboundarycells;
extern bool closedbrush;
extern bool brushinterfaceOn;
extern int num_shapecontrol_pts;
extern int num_curvepts_output;
extern ctr_point **brushpts;
extern int nbrushpts;
extern int curMaxNumBrushPts;

extern double brush_width;

extern int Num_SmoothRegionpoints;


/*for multiple region design*/
extern TrajectoryList *sketchlines;
extern SketchList *brushes;
extern void init_brushlist();
extern int NDesignRegions;
extern int cur_chosen_region; 

extern bool highwayexisted;

extern MapBoundaryList *mapboundarylist;

/*     construct the sketch based network and its blocks    */
extern RegionBlockList *sketchblocklist;
extern StreetNet *sketchnet;
extern void construct_regionblocks_edgewise(StreetNet *net, RegionBlockList *aregionblocklist);

extern double BoundRegionWidth;
extern int gen_regElem_sketches;

#include "StreetNetPanel.h"
#include "StreetEditPanel.h"

extern StreetNetPanel *g_streetnetpanel;
extern StreetEditPanel *g_streeteditpanel;


extern CEdit *g_m_EdCtrlSetStreetLevels;
extern int *g_m_edStreetLevel;
extern int which_level;

extern bool majorroadsexisted;

/********************/
/*   for globally updating the number of design brushes    */
int *g_m_edNumBrushes;


// MyDlg1 dialog

IMPLEMENT_DYNAMIC(MyDlg1, CDialog)
MyDlg1::MyDlg1(CWnd* pParent /*=NULL*/)
	: CDialog(MyDlg1::IDD, pParent)
	, m_chkTensorDesignOn(FALSE)
	, m_chkMoveElemOn(FALSE)
	, m_chkEditElementOn(FALSE)
	, m_chkRemoveElemOn(FALSE)
	, m_chkBrushInterfaceOn(FALSE)
	, m_chkCombineWithOtherPartsOn(TRUE)
	, m_edBrushWidth(5.)
	, m_rdAddVFSingularElem(FALSE)
	, m_chkRegionSmoothOn(FALSE)
	, m_chkEnableSketchBasedDesign(FALSE)
	, m_edNumBrushes(0)
	, m_selRegion(_T(""))
	, m_chkShowSketchesOn(FALSE)
	, m_rdSketchMajorHighways(FALSE)
	, m_chkUseBoundsAsRoadsOn(FALSE)
	, m_chkUseMajRoadsAsSketchesOn(FALSE)
	, m_rdScalarSingularElem(FALSE)
	, m_chkEnableHeighfieldDesignOn(FALSE)
	, m_chkShowSegmentGraphOn(FALSE)
{
	g_m_edNumBrushes=&m_edNumBrushes;

	m_rdAddVFSingularElem=5;
	
}

MyDlg1::~MyDlg1()
{
}

void MyDlg1::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Check(pDX, IDC_CHECK_TENSORDESIGNON, m_chkTensorDesignOn);
	DDX_Check(pDX, IDC_CHECK_MOVEELEM, m_chkMoveElemOn);
	DDX_Check(pDX, IDC_CHECK_EDITANELEM, m_chkEditElementOn);
	DDX_Check(pDX, IDC_CHECK_REMOVEELEM, m_chkRemoveElemOn);
	DDX_Check(pDX, IDC_CHECK_ACTIVATEBRUSH, m_chkBrushInterfaceOn);
	DDX_Check(pDX, IDC_CHECK_COMBINEWITHPREVIOUSTENSOR, m_chkCombineWithOtherPartsOn);
	DDX_Text(pDX, IDC_EDIT1, m_edBrushWidth);
	DDX_Radio(pDX, IDC_RADIO_ADDAWEDGE, m_rdAddVFSingularElem);
	DDX_Check(pDX, IDC_CHECK_REGIONSMOOTHON, m_chkRegionSmoothOn);
	DDX_Check(pDX, IDC_CHECK_USEDRAWSKETCHON, m_chkEnableSketchBasedDesign);
	DDX_Text(pDX, IDC_EDIT_CURRENTNUMSKETCHES, m_edNumBrushes);
	DDX_Control(pDX, IDC_COMBO_DIFFREGIONS, m_ctrlcomboxDesignRegionIndex);
	DDX_CBString(pDX, IDC_COMBO_DIFFREGIONS, m_selRegion);
	DDX_Check(pDX, IDC_CHECK_DISPLAYSKETCHES, m_chkShowSketchesOn);
	DDX_Control(pDX, IDC_SLIDER1, m_ctrlSliderBrushRegionWidth);
	DDX_Radio(pDX, IDC_RADIO_SKETCHMAJORROADS, m_rdSketchMajorHighways);
	DDX_Check(pDX, IDC_CHECK_USEBOUNDARYASROADS, m_chkUseBoundsAsRoadsOn);
	DDX_Check(pDX, IDC_CHECK_USEMAJROADSASSKETCHES, m_chkUseMajRoadsAsSketchesOn);
	DDX_Radio(pDX, IDC_RADIO_SCALARMAXIMUM, m_rdScalarSingularElem);
	DDX_Check(pDX, IDC_CHECK_DESIGNSCALARFIELD, m_chkEnableHeighfieldDesignOn);
	DDX_Check(pDX, IDC_CHECK_SHOWSEGMENTGRAPH, m_chkShowSegmentGraphOn);
	DDX_Control(pDX, IDC_BUTTON_CLEARSKETCHES, Trace);
}


BEGIN_MESSAGE_MAP(MyDlg1, CDialog)
	ON_BN_CLICKED(IDC_CHECK_TENSORDESIGNON, OnBnClickedCheckTensordesignon)
	ON_BN_CLICKED(IDC_CHECK_MOVEELEM, OnBnClickedCheckMoveelem)
	ON_BN_CLICKED(IDC_CHECK_EDITANELEM, OnBnClickedCheckEditanelem)
	ON_BN_CLICKED(IDC_CHECK_REMOVEELEM, OnBnClickedCheckRemoveelem)
	ON_BN_CLICKED(IDC_CHECK_ACTIVATEBRUSH, OnBnClickedCheckActivatebrush)
	ON_BN_CLICKED(IDC_CHECK_COMBINEWITHPREVIOUSTENSOR, OnBnClickedCheckCombinewithprevioustensor)
	ON_BN_CLICKED(IDC_BUTTON_GENTENFROMBRUSH, OnBnClickedButtonGentenfrombrush)
	ON_BN_CLICKED(IDC_BUTTON_SAVEBRUSH, OnBnClickedButtonSavebrush)
	ON_BN_CLICKED(IDC_BUTTON_LOADBRUSH, OnBnClickedButtonLoadbrush)
	ON_BN_CLICKED(IDC_RADIO_ADDAWEDGE, OnBnClickedRadioAddawedge)
	ON_BN_CLICKED(IDC_RADIO_ADDATRISECTOR, OnBnClickedRadioAddatrisector)
	ON_BN_CLICKED(IDC_RADIO_ADDANODE, OnBnClickedRadioAddanode)
	ON_BN_CLICKED(IDC_RADIO_ADDACENTER, OnBnClickedRadioAddacenter)
	ON_BN_CLICKED(IDC_RADIO_ADDASADDLE, OnBnClickedRadioAddasaddle)
	ON_BN_CLICKED(IDC_RADIO_ADDAREGULAR, OnBnClickedRadioAddaregular)
	ON_BN_CLICKED(IDC_CHECK_REGIONSMOOTHON, OnBnClickedCheckRegionsmoothon)
	ON_BN_CLICKED(IDC_BUTTON_SMOOTHONCE, OnBnClickedButtonSmoothonce)
	ON_BN_CLICKED(IDC_BUTTON_SAVEFIELDPERVER, OnBnClickedButtonSavefieldperver)
	ON_BN_CLICKED(IDC_BUTTON_LOADFIELDPERVERT, OnBnClickedButtonLoadfieldpervert)
	ON_BN_CLICKED(IDC_CHECK_USEDRAWSKETCHON, OnBnClickedCheckUsedrawsketchon)
	ON_CBN_SELCHANGE(IDC_COMBO_DIFFREGIONS, OnCbnSelchangeComboDiffregions)
	ON_BN_CLICKED(IDC_CHECK_DISPLAYSKETCHES, OnBnClickedCheckDisplaysketches)
	ON_BN_CLICKED(IDC_BUTTON_SKETCHBASEDTENGEN, OnBnClickedButtonSketchbasedtengen)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER1, OnNMCustomdrawSlider1)
	ON_BN_CLICKED(IDC_BUTTON_CLEARSKETCHES, OnBnClickedButtonClearsketches)
	ON_EN_CHANGE(IDC_EDIT1, OnEnChangeEdit1)
	ON_BN_CLICKED(IDC_BUTTON_SAVESKETCHES, OnBnClickedButtonSavesketches)
	ON_BN_CLICKED(IDC_BUTTON_LOADSKETCHES, OnBnClickedButtonLoadsketches)
	ON_BN_CLICKED(IDC_BUTTON_SEGMENTREGION, OnBnClickedButtonSegmentregion)
	ON_BN_CLICKED(IDC_RADIO_SKETCHMAJORROADS, OnBnClickedRadioSketchmajorroads)
	ON_BN_CLICKED(IDC_RADIO_SKETCHHIGHWAY, OnBnClickedRadioSketchhighway)
	ON_BN_CLICKED(IDC_CHECK_USEBOUNDARYASROADS, OnBnClickedCheckUseboundaryasroads)
	ON_BN_CLICKED(IDC_BUTTON_CLEARTHEFIELD, OnBnClickedButtonClearthefield)
	ON_BN_CLICKED(IDC_CHECK_USEMAJROADSASSKETCHES, OnBnClickedCheckUsemajroadsassketches)
	ON_BN_CLICKED(IDC_RADIO_SCALARMAXIMUM, OnBnClickedRadioScalarmaximum)
	ON_BN_CLICKED(IDC_RADIO_SCALARMINIMUM, OnBnClickedRadioScalarminimum)
	ON_BN_CLICKED(IDC_CHECK_DESIGNSCALARFIELD, OnBnClickedCheckDesignscalarfield)
	ON_BN_CLICKED(IDC_CHECK_SHOWSEGMENTGRAPH, OnBnClickedCheckShowsegmentgraph)
	ON_BN_CLICKED(IDC_BUTTON_SAVEDESIGNELEMS, OnBnClickedButtonSavedesignelems)
	ON_BN_CLICKED(IDC_BUTTON_LOADDESIGNELEMS, OnBnClickedButtonLoaddesignelems)
END_MESSAGE_MAP()


// MyDlg1 message handlers

void MyDlg1::Reset()
{
	//m_rdAddVFSingularElem=FALSE;
	m_chkTensorDesignOn=FALSE;
	m_chkMoveElemOn=FALSE;
	m_chkEditElementOn=FALSE;
	m_chkRemoveElemOn=FALSE;
	m_chkBrushInterfaceOn=FALSE;
	m_chkCombineWithOtherPartsOn=FALSE;
	m_edBrushWidth=6.;
}

void MyDlg1::OnBnClickedCheckTensordesignon()
{
	// TODO: Add your control notification handler code here
	if(m_chkTensorDesignOn == FALSE)
	{
		m_chkTensorDesignOn = TRUE;
		m_chkEditElementOn = FALSE;
		m_chkMoveElemOn = FALSE;
		m_chkRemoveElemOn = FALSE;
		m_chkRegionSmoothOn=FALSE;
		m_chkBrushInterfaceOn=FALSE;
		m_chkEnableHeighfieldDesignOn=FALSE;
		g_streeteditpanel->m_chkSelStreetRegToEditOn=FALSE;
		g_streeteditpanel->UpdateData(FALSE);


		sharedvars.TensorDesignOn = true;
		sharedvars.MoveElemOn=false;
		sharedvars.RemoveElemOn=false;
		sharedvars.EditElementOn=false;
		sharedvars.RegionSmoothOn=false;
		sharedvars.BrushInterfaceOn=false;
		sharedvars.EnableSketchBasedDesign=false;

		sharedvars.SelStreetRegToEditOn=false;

		sharedvars.EnableHeighfieldDesignOn=false;

		sharedvars.SelStreetRegToEditOn=false;
		
		g_pclGlView->deleteElemOn=false;
		g_pclGlView->EditModeOn = 0;
		g_pclGlView->MoveElemOn = 0;
		
	}
	else
	{
		m_chkTensorDesignOn = FALSE;
		sharedvars.TensorDesignOn = false;
	}
	UpdateData(FALSE);
}


void MyDlg1::OnBnClickedCheckMoveelem()
{
	// TODO: Add your control notification handler code here
	if(m_chkMoveElemOn == FALSE)
	{
		m_chkMoveElemOn = TRUE;
		m_chkEditElementOn = FALSE;
		m_chkTensorDesignOn = FALSE;
		m_chkRegionSmoothOn=FALSE;
		m_chkBrushInterfaceOn=FALSE;
		m_chkRemoveElemOn = FALSE;
		
		m_chkEnableHeighfieldDesignOn=FALSE;

		sharedvars.MoveElemOn = true;
		sharedvars.TensorDesignOn = false;
		sharedvars.tenElemType = -1;
		sharedvars.EditElementOn = false;
		sharedvars.RegionSmoothOn=false;
		sharedvars.BrushInterfaceOn=false;
		sharedvars.EnableSketchBasedDesign=false;
		sharedvars.RemoveElemOn=false;
		sharedvars.EnableHeighfieldDesignOn=false;

		sharedvars.SelStreetRegToEditOn=false;
		
		g_pclGlView->deleteElemOn=false;
		g_pclGlView->MoveElemOn = 1;
		g_pclGlView->EditModeOn = 0;
	}
	else
	{
		m_chkMoveElemOn = FALSE;
		sharedvars.MoveElemOn = false;
		g_pclGlView->MoveElemOn = 0;
	}
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedCheckEditanelem()
{
	// TODO: Add your control notification handler code here
	if(m_chkEditElementOn == FALSE)
	{
		m_chkEditElementOn = TRUE;
		m_chkMoveElemOn = FALSE;
		m_chkRemoveElemOn = FALSE;
		m_chkTensorDesignOn = FALSE;
		m_chkRegionSmoothOn=FALSE;
		m_chkBrushInterfaceOn=FALSE;
		m_chkEnableHeighfieldDesignOn=FALSE;

		sharedvars.EditElementOn = true;
		sharedvars.RemoveElemOn = false;
		sharedvars.MoveElemOn = false;
		sharedvars.TensorDesignOn = false;
		sharedvars.tenElemType = -1;
		sharedvars.RegionSmoothOn=false;
		sharedvars.BrushInterfaceOn=false;
		sharedvars.EnableSketchBasedDesign=false;

		sharedvars.SelStreetRegToEditOn=false;
		sharedvars.EnableHeighfieldDesignOn=false;

		g_pclGlView->EditModeOn = 1;
		g_pclGlView->MoveElemOn = 0;
		g_pclGlView->deleteElemOn=false;
	}
	else
	{
		m_chkEditElementOn = FALSE;
		sharedvars.EditElementOn = false;
		g_pclGlView->EditModeOn = 0;
	}
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedCheckRemoveelem()
{
	// TODO: Add your control notification handler code here
	if(m_chkRemoveElemOn == FALSE)
	{
		m_chkRemoveElemOn = TRUE;
		m_chkEditElementOn = FALSE;
		m_chkMoveElemOn = FALSE;
		m_chkTensorDesignOn = FALSE;
		m_chkRegionSmoothOn=FALSE;
		m_chkBrushInterfaceOn=FALSE;
		m_chkEnableHeighfieldDesignOn=FALSE;


		sharedvars.RemoveElemOn = true;
		sharedvars.EditElementOn = false;
		sharedvars.MoveElemOn = false;
		sharedvars.TensorDesignOn = false;
		sharedvars.tenElemType = -1;
		sharedvars.RegionSmoothOn=false;
		sharedvars.BrushInterfaceOn=false;
		sharedvars.EnableSketchBasedDesign=false;

		sharedvars.SelStreetRegToEditOn=false;
		sharedvars.EnableHeighfieldDesignOn=false;

		g_pclGlView->deleteElemOn=true;
		g_pclGlView->EditModeOn = 0;
		g_pclGlView->MoveElemOn = 0;
	}
	else
	{
		m_chkRemoveElemOn = FALSE;
		sharedvars.RemoveElemOn = false;
		g_pclGlView->deleteElemOn=false;
	}

	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedCheckActivatebrush()
{
	// TODO: Add your control notification handler code here
	if(m_chkBrushInterfaceOn==FALSE)
	{
		m_chkBrushInterfaceOn=TRUE;
		m_chkRemoveElemOn = FALSE;
		m_chkEditElementOn = FALSE;
		m_chkMoveElemOn = FALSE;
		m_chkTensorDesignOn = FALSE;
		g_streeteditpanel->m_chkSelStreetRegToEditOn=FALSE;
		g_streeteditpanel->UpdateData(FALSE);

		sharedvars.BrushInterfaceOn = true;
		sharedvars.RemoveElemOn = false;
		sharedvars.EditElementOn = false;
		sharedvars.MoveElemOn = false;
		sharedvars.TensorDesignOn = false;
		sharedvars.tenElemType = -1;
		sharedvars.SelStreetRegToEditOn=false;
		sharedvars.EnableHeighfieldDesignOn=false;
		
		m_chkRegionSmoothOn=FALSE;

		sharedvars.RegionSmoothOn=false;

		g_pclGlView->EditModeOn = 0;
		g_pclGlView->MoveElemOn = 0;
		g_pclGlView->ShapeControlPtsOn = 1;

		g_pclGlView->FinisheCtrptsel = 0;  ////Do not let user move control point now
		g_pclGlView->DisplaySmoothRegionOn = 1;
	    
		////Initialize
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
		m_chkBrushInterfaceOn=FALSE;
		sharedvars.BrushInterfaceOn = false;
		g_pclGlView->ShapeControlPtsOn = 0;
	}
	UpdateData(FALSE);
}


void MyDlg1::OnBnClickedCheckCombinewithprevioustensor()
{
	// TODO: Add your control notification handler code here
	if(m_chkCombineWithOtherPartsOn==FALSE)
	{
		m_chkCombineWithOtherPartsOn=TRUE;
		sharedvars.CombineWithOtherPartsOn=true;
		brushinterfaceOn=true;
	}
	else
	{
		if(MessageBox("This will remove all the previous design result, do you want to continue",
			"Warning", MB_OKCANCEL)==IDOK)
		{
			m_chkCombineWithOtherPartsOn=FALSE;
			sharedvars.CombineWithOtherPartsOn=false;
			brushinterfaceOn=false;
		}
	}
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedButtonGentenfrombrush()
{
	// TODO: Add your control notification handler code here
	UpdateData(TRUE);
		
	init_dis_cells();

	brush_width = m_edBrushWidth;
	fast_marching_quad_withDis(m_edBrushWidth);

	//get_brushbuffer_quad_withDis(brush_width+brush_width*1.5);
		
	init_regular_ten_designelems();

	get_bound_approDir();

	/*we further combine it with the current existing tensor field*/
	if(m_chkCombineWithOtherPartsOn == FALSE)
	{
		/*initialize all the tensor field*/
		clear_all_tens();
		//init_verts_all();
		//reset_inland();
	}
	get_brushinterior_ten();

		//smooth_Jac_quadregion();
		cal_all_eigenvecs_quad();

		init_degpts();
		/*render alpha map*/
		render_alpha_map_quad(false);
		render_alpha_map_quad(true);

		/*we still need to remove the flags inside the brush*/
		brushinterfaceOn = false;
		//sharedvars.BrushInterfaceOn=false;

	/*we need other condition flag to separate the brush stroke and boundary-based  
	generation from the image*/
	//else if(flag_loadmap = true && m_chkCombineWithOtherPartsOn ==TRUE)
	//{
	//	get_imgbound_approDir();
	//}

	locate_degpts_cells_tranvec_quad();

	UpdateData(FALSE);
}


void MyDlg1::OnBnClickedButtonSavebrush()
{
	// TODO: Add your control notification handler code here
	if(m_chkBrushInterfaceOn == TRUE)
	{
		CFile f;
		char filename[256];

		char strFilter[] = { "Text Files (*.bru)|*.bru|All Files (*.*)|*.*||" };

		CFileDialog FileDlg(FALSE, ".bru", NULL, 0, strFilter);
		if( FileDlg.DoModal() == IDOK )
		{

			CString f2 = FileDlg.GetFileName();
			CString f3 = f2.Left(f2.Find('.'))+".bru";
			const char *tempf = f3;
			strcpy(filename, f3);
	 
			save_brush_curve(filename);
		}
	}
}

void MyDlg1::OnBnClickedButtonLoadbrush()
{
	// TODO: Add your control notification handler code here
	if(m_chkBrushInterfaceOn == TRUE)
	{
		CFile f;
		char filename[256];

		char strFilter[] = { "Text Files (*.bru)|*.bru" };

		CFileDialog FileDlg(TRUE, ".bru", NULL, 0, strFilter);
		if( FileDlg.DoModal() == IDOK )
		{

			CString f2 = FileDlg.GetFileName();
			const char *tempf = f2;
			strcpy(filename, f2);
	 
			load_brush_curve(filename);

			//init_dis_verts();
			init_verts_for_brush();
			init_dis_cells();
			get_cellstrip_curve_quad(boundarycells, nboundarycells);
		}
	}
}


void MyDlg1::OnBnClickedRadioAddawedge()
{
	// TODO: Add your control notification handler code here
	m_rdAddVFSingularElem = 0;
	sharedvars.tenElemType = 0;
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedRadioAddatrisector()
{
	// TODO: Add your control notification handler code here
	m_rdAddVFSingularElem = 1;
	sharedvars.tenElemType = 1;
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedRadioAddanode()
{
	// TODO: Add your control notification handler code here
	m_rdAddVFSingularElem = 2;
	sharedvars.tenElemType = 2;
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedRadioAddacenter()
{
	// TODO: Add your control notification handler code here
	m_rdAddVFSingularElem = 3;
	sharedvars.tenElemType = 3;
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedRadioAddasaddle()
{
	// TODO: Add your control notification handler code here
	m_rdAddVFSingularElem = 4;
	sharedvars.tenElemType = 4;
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedRadioAddaregular()
{
	// TODO: Add your control notification handler code here
	m_rdAddVFSingularElem = 5;
	sharedvars.tenElemType = 5;
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedCheckRegionsmoothon()
{
	// TODO: Add your control notification handler code here
	if(m_chkRegionSmoothOn==FALSE)
	{
		m_chkRegionSmoothOn=TRUE;
		m_chkBrushInterfaceOn=FALSE;
		m_chkRemoveElemOn = FALSE;
		m_chkEditElementOn = FALSE;
		m_chkMoveElemOn = FALSE;
		m_chkTensorDesignOn = FALSE;
		m_chkEnableHeighfieldDesignOn=FALSE;

		sharedvars.RegionSmoothOn=true;
		sharedvars.BrushInterfaceOn = false;
		sharedvars.RemoveElemOn = false;
		sharedvars.EditElementOn = false;
		sharedvars.MoveElemOn = false;
		sharedvars.TensorDesignOn = false;
		sharedvars.EnableSketchBasedDesign=false;
		sharedvars.EnableHeighfieldDesignOn=false;

		m_rdAddVFSingularElem = -1;  
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
		m_chkRegionSmoothOn=FALSE;
		g_pclGlView->SmoothOn = 0;              ////SELECTION MODE
		g_pclGlView->PickPointOn = 0;
		g_pclGlView->DisplaySmoothRegionOn = 0;
		sharedvars.RegionSmoothOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedButtonSmoothonce()
{
	// TODO: Add your control notification handler code here
	if(g_pclGlView->SmoothOn == 1)
	{
		if(Num_SmoothRegionpoints>2)
		{
			init_quad_regionsmooth();
			find_boundarycells();
			mark_boundVerts();
			find_innerVerts();

			smooth_Jac_quadregion();
			cal_all_eigenvecs_quad();

			init_degpts();
			/*render alpha map*/
			render_alpha_map_quad(false);
			render_alpha_map_quad(true);

			locate_degpts_cells_tranvec_quad();
		
			Num_SmoothRegionpoints = 0;
		}
	}
}

void MyDlg1::OnBnClickedButtonSavefieldperver()
{
	// TODO: Add your control notification handler code here
	CFile f;

	char strFilter[] = { "Text Files (*.txt)|*.txt|All Files (*.*)|*.*||" };

	CFileDialog FileDlg(FALSE, ".txt", NULL, 0, strFilter);
	if( FileDlg.DoModal() == IDOK )
	{
		CString f2 = FileDlg.GetFileName();
		const char *filename = f2;
 
		char filename2[256];
		strcpy(filename2, filename);
		save_tenField_perVer(filename2);
	}
}

void MyDlg1::OnBnClickedButtonLoadfieldpervert()
{
	// TODO: Add your control notification handler code here
	CFile f;

	char strFilter[] = { "Text Files (*.txt)|*.txt|All Files (*.*)|*.*||" };

	CFileDialog FileDlg(TRUE, ".txt", NULL, 0, strFilter);
	if( FileDlg.DoModal() == IDOK )
	{

		clock_t start, finish;
		start = clock();

		CString f2 = FileDlg.GetFileName();

		const char *filename = f2;
		
		char filename2[256];
		strcpy(filename2, filename);
		load_tenField_perVer(filename2);
		
		cal_all_eigenvecs_quad();

		init_degpts();
		/*render alpha map*/
		render_alpha_map_quad(false);
		render_alpha_map_quad(true);

		locate_degpts_cells_tranvec_quad();

		num_shapecontrol_pts = 0;
	    num_curvepts_output = 0;

		//num_cycleedges = 0;
		//num_triangles_designcurve = 0;
		//num_innertriangles = 0;

		//InnerBoundary.num = 0;
		//OuterBoundary.num = 0;

		/////
		
		finish = clock();
		
		UpdateData(FALSE);
	}
}

void MyDlg1::OnBnClickedCheckUsedrawsketchon()
{
	// TODO: Add your control notification handler code here
	if(m_chkEnableSketchBasedDesign==FALSE)
	{
		m_chkEnableSketchBasedDesign=TRUE;
		m_chkShowSketchesOn=TRUE;
		m_chkBrushInterfaceOn=FALSE;
		m_chkRemoveElemOn = FALSE;
		m_chkEditElementOn = FALSE;
		m_chkMoveElemOn = FALSE;
		m_chkTensorDesignOn = FALSE;
		m_chkEnableHeighfieldDesignOn=FALSE;

		sharedvars.ShowSketchesOn=true;
		sharedvars.EnableSketchBasedDesign=true;
		sharedvars.BrushInterfaceOn = false;
		sharedvars.RemoveElemOn = false;
		sharedvars.EditElementOn = false;
		sharedvars.MoveElemOn = false;
		sharedvars.TensorDesignOn = false;
		sharedvars.tenElemType = -1;

		g_pclGlView->EditModeOn = 0;
		g_pclGlView->MoveElemOn = 0;
		g_pclGlView->ShapeControlPtsOn = 1;
		g_pclGlView->deleteElemOn = false;

		sharedvars.SelStreetRegToEditOn=false;
		sharedvars.EnableHeighfieldDesignOn=false;


		g_pclGlView->FinisheCtrptsel = 0;  ////Do not let user move control point now
	    
		////Initialize
		num_shapecontrol_pts = 0;
	    num_curvepts_output = 0;

		reset_vert_dis();

		/*initialize the brush point list*/
		if(brushes == NULL)
		{
			init_brushlist();
		}


		//else
		//{
		//	delete brushes;
		//	brushes=NULL;
		//	init_brushlist();
		//	init_sketchnet();
		//}

		//release_domain_boundaries();
		//init_domain_boundaries();
		//highwayexisted=false;

		if(m_rdSketchMajorHighways==0)
		{
			g_m_EdCtrlSetStreetLevels->EnableWindow(FALSE);
			which_level=(*g_m_edStreetLevel)=2;
			g_streetnetpanel->UpdateData(FALSE);
			//g_streetnetpanel->OnEnChangeEditStreetlevel();
			//init_level1_placement();
		}

	}
	else
	{
		m_chkEnableSketchBasedDesign=FALSE;
		sharedvars.EnableSketchBasedDesign=false;
		//highwayexisted=false;
		g_m_EdCtrlSetStreetLevels->EnableWindow(TRUE);
	}
	UpdateData(FALSE);
}

BOOL MyDlg1::OnInitDialog()
{
	CDialog::OnInitDialog();

	// TODO:  Add extra initialization here

	m_ctrlcomboxDesignRegionIndex.AddString("0");

	m_ctrlSliderBrushRegionWidth.SetRange(1, 300);
	m_ctrlSliderBrushRegionWidth.SetPos(5);
	m_ctrlSliderBrushRegionWidth.SetTicFreq(1);
	m_ctrlSliderBrushRegionWidth.SetTic(0);
	BoundRegionWidth=(double)m_ctrlSliderBrushRegionWidth.GetPos();
	m_edBrushWidth=BoundRegionWidth;

	m_rdAddVFSingularElem=5;
	sharedvars.tenElemType=5;


	UpdateData(FALSE);

	return TRUE;  // return TRUE unless you set the focus to a control
	// EXCEPTION: OCX Property Pages should return FALSE
}

void MyDlg1::OnCbnSelchangeComboDiffregions()
{
	// TODO: Add your control notification handler code here
	m_ctrlcomboxDesignRegionIndex.GetLBText(m_ctrlcomboxDesignRegionIndex.GetCurSel(), m_selRegion);
	cur_chosen_region=atoi(m_selRegion);
}

void MyDlg1::OnBnClickedCheckDisplaysketches()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowSketchesOn==FALSE)
	{
		m_chkShowSketchesOn=TRUE;
		sharedvars.ShowSketchesOn=true;
	}
	else
	{
		m_chkShowSketchesOn=FALSE;
		sharedvars.ShowSketchesOn=false;
	}
	UpdateData(FALSE);
}




void MyDlg1::OnBnClickedButtonSketchbasedtengen()
{
	// TODO: Add your control notification handler code here
	FILE *fp;
	//fp=fopen("sketchbased_error.txt", "w");
	//fprintf(fp, "initializing skecthes\n");
	//fclose(fp);
		
	release_domain_boundaries();
	init_domain_boundaries();
	
	/*   consider to use the boundaries of the map to generate the sketches  */
	if(sharedvars.UseBoundsAsSketchesOn&&mapboundarylist!=NULL)
	{
		convert_mapbounds_to_trajs();
	}

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

	if(m_rdSketchMajorHighways==0)
	{
		/*   copy to the major_level1   */
		init_level1_placement();
		copy_sketchlines_to_major_level1();

		/*   set the major road existing flag   */
		majorroadsexisted=true;
	}
	else
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

			init_degpts();
			/*render alpha map*/
			render_alpha_map_quad(false);
			render_alpha_map_quad(true);

			locate_degpts_cells_tranvec_quad();
		
		//fp=fopen("sketchbased_error.txt", "w");
		//fprintf(fp, "finish rendering the alpha map.\n");
		//fclose(fp);

	}

	/*  Disable the sketch based design  */
	//sharedvars.EnableSketchBasedDesign=false;
}

extern void compute_brush_Reg();

void MyDlg1::OnNMCustomdrawSlider1(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	// TODO: Add your control notification handler code here
	BoundRegionWidth=(double)m_ctrlSliderBrushRegionWidth.GetPos();
	m_edBrushWidth=BoundRegionWidth;
	UpdateData(FALSE);

	//compute_brush_Reg();


	*pResult = 0;
}

void MyDlg1::OnBnClickedButtonClearsketches()
{
	// TODO: Add your control notification handler code here
		/*initialize the brush point list*/
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

		release_domain_boundaries();
		init_domain_boundaries();
		highwayexisted=false;
		gen_regElem_sketches=0;
}

void MyDlg1::OnEnChangeEdit1()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData();
	BoundRegionWidth=m_edBrushWidth;
	m_ctrlSliderBrushRegionWidth.SetPos(floor(BoundRegionWidth));
}

void MyDlg1::OnBnClickedButtonSavesketches()
{
	// TODO: Add your control notification handler code here
	if(m_chkEnableSketchBasedDesign== TRUE)
	{
		CFile f;
		char filename[256];

		char strFilter[] = { "Text Files (*.skt)|*.skt|All Files (*.*)|*.*||" };

		CFileDialog FileDlg(FALSE, ".skt", NULL, 0, strFilter);
		if( FileDlg.DoModal() == IDOK )
		{

			CString f2 = FileDlg.GetFileName();
			CString f3 = f2.Left(f2.Find('.'))+".skt";
			const char *tempf = f3;
			strcpy(filename, f3);
	 
			save_sketch_brushes(filename);
		}
	}
}

void MyDlg1::OnBnClickedButtonLoadsketches()
{
	// TODO: Add your control notification handler code here
	if(m_chkEnableSketchBasedDesign == TRUE)
	{
		CFile f;
		char filename[256];

		char strFilter[] = { "Text Files (*.skt)|*.skt" };

		CFileDialog FileDlg(TRUE, ".skt", NULL, 0, strFilter);
		if( FileDlg.DoModal() == IDOK )
		{

			CString f2 = FileDlg.GetFileName();
			const char *tempf = f2;
			strcpy(filename, f2);
	 
			load_sketch_brushes(filename);
		}
	}
}

void MyDlg1::OnBnClickedButtonSegmentregion()
{
	// TODO: Add your control notification handler code here

	if(brushes==NULL) return;
	if(sketchlines == NULL) return;

	release_domain_boundaries();
	init_domain_boundaries();

	if(sharedvars.UseBoundsAsSketchesOn&&mapboundarylist!=NULL)
	{
		convert_mapbounds_to_trajs();
	}

	if(sharedvars.UseMajRoadsAsSketchesOn /*&& major_level1!=NULL&&minor_level1!=NULL*/)
	{
		convert_majRoads_to_sketches();
	}

	FILE *fp;

	init_sketchnet();
	init_sketchlineinfo_incells();
	update_sketchcurveinfo();
	cal_sketchlines_intersects();
	search_for_sketchbased_graph();
	
	if(!sharedvars.UseMajRoadsAsSketchesOn)
	{
		/*  if we use the tensor field generated major road network
		    don't do the extension of the roads 1/7/2008
		*/
		extend_sketchcurves();

		init_sketchnet();
		init_sketchlineinfo_incells();
		update_sketchcurveinfo();
		cal_sketchlines_intersects();
		search_for_sketchbased_graph();
	}
	
	/*      extract the blocks        */
	init_sketchblocklist();
	construct_regionblocks_edgewise(sketchnet, sketchblocklist);
	
	/*      mark the vertices of the mesh       */
extern void color_subregionblocks(RegionBlockList *aregionblocklist, StreetNet *net);
	color_subregionblocks(sketchblocklist, sketchnet);
	
	/*  reassign the indices for the existing elements (could contain bug) 1/12/2008  */
	gen_regElem_sketches=0;
	reset_region_indexes_all_designElems();

	if(m_rdSketchMajorHighways==0)
	{
		/*   copy to the major_level1   */
		if(sharedvars.UseMajRoadsAsSketchesOn && majorroadsexisted)
		{
			/*  this can be problematic when consider the map boundaries 12/28/2007  */
		}
		else
		{
			init_level1_placement();
			copy_sketchlines_to_major_level1();

			/*   set the major road existing flag   */
			majorroadsexisted=true;
		}
	}
	else
		highwayexisted=true;

	/*  please save the global field now 
	    possible bug  1/18/2008
	*/
	store_cur_ten();

	if(sharedvars.UseMajRoadsAsSketchesOn)
	{
		sharedvars.ShowSketchesOn=true;
		m_chkShowSketchesOn=TRUE;
	}
	
}

void MyDlg1::OnBnClickedRadioSketchmajorroads()
{
	// TODO: Add your control notification handler code here
	m_rdSketchMajorHighways=0;
	sharedvars.rdSketchMajorHighways=false;

	if(m_chkEnableSketchBasedDesign==TRUE)
	{
		g_m_EdCtrlSetStreetLevels->EnableWindow(FALSE);
		which_level=(*g_m_edStreetLevel)=2;
		g_streetnetpanel->UpdateData(FALSE);
		g_streetnetpanel->m_ctrlMajRoadSettingButton.EnableWindow(FALSE);
		//g_streetnetpanel->OnEnChangeEditStreetlevel();
	}
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedRadioSketchhighway()
{
	// TODO: Add your control notification handler code here
	m_rdSketchMajorHighways=1;
	sharedvars.rdSketchMajorHighways=true;
	if(m_chkEnableSketchBasedDesign==TRUE)
	{
		g_m_EdCtrlSetStreetLevels->EnableWindow(TRUE);
	}
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedCheckUseboundaryasroads()
{
	// TODO: Add your control notification handler code here
	if(m_chkUseBoundsAsRoadsOn==FALSE)
	{
		m_chkUseBoundsAsRoadsOn=TRUE;
		sharedvars.UseBoundsAsRoadsOn=true;
	}
	else
	{
		m_chkUseBoundsAsRoadsOn=FALSE;
		sharedvars.UseBoundsAsRoadsOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedButtonClearthefield()
{
	// TODO: Add your control notification handler code here
	clear_all_tens();
	init_ten_designelems();
	tensor_init_tex();
	init_degpts();
}

void MyDlg1::OnBnClickedCheckUsemajroadsassketches()
{
	// TODO: Add your control notification handler code here
	if(m_chkUseMajRoadsAsSketchesOn==FALSE)
	{
		m_chkUseMajRoadsAsSketchesOn=TRUE;
		sharedvars.UseMajRoadsAsSketchesOn=true;

		/*  copy the major_level1 and minor_level1 to the sketchlines!  */
	}
	else
	{
		m_chkUseMajRoadsAsSketchesOn=FALSE;
		sharedvars.UseMajRoadsAsSketchesOn=false;
	}
}

void MyDlg1::OnBnClickedRadioScalarmaximum()
{
	// TODO: Add your control notification handler code here
	m_rdScalarSingularElem=0;
	sharedvars.rdScalarSingularElem=false;
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedRadioScalarminimum()
{
	// TODO: Add your control notification handler code here
	m_rdScalarSingularElem=1;
	sharedvars.rdScalarSingularElem=true;
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedCheckDesignscalarfield()
{
	// TODO: Add your control notification handler code here
	if(m_chkEnableHeighfieldDesignOn==FALSE)
	{
		m_chkEnableHeighfieldDesignOn=TRUE;
		sharedvars.EnableHeighfieldDesignOn=true;

		m_chkTensorDesignOn = FALSE;
		m_chkEditElementOn = FALSE;
		m_chkMoveElemOn = FALSE;
		m_chkRemoveElemOn = FALSE;
		m_chkRegionSmoothOn=FALSE;
		m_chkBrushInterfaceOn=FALSE;
		g_streeteditpanel->m_chkSelStreetRegToEditOn=FALSE;
		g_streeteditpanel->UpdateData(FALSE);


		sharedvars.TensorDesignOn = false;
		sharedvars.MoveElemOn=false;
		sharedvars.RemoveElemOn=false;
		sharedvars.EditElementOn=false;
		sharedvars.RegionSmoothOn=false;
		sharedvars.BrushInterfaceOn=false;
		sharedvars.EnableSketchBasedDesign=false;

		sharedvars.EditStreetNetOn=false;

		sharedvars.SelStreetRegToEditOn=false;


		sharedvars.SelStreetRegToEditOn=false;
		g_streetnetpanel->m_chkEditStreetNetOn=FALSE;
		sharedvars.EditStreetNetOn=false;
		
		g_pclGlView->deleteElemOn=false;
		g_pclGlView->EditModeOn = 0;
		g_pclGlView->MoveElemOn = 0;
	}
	else
	{
		m_chkEnableHeighfieldDesignOn=FALSE;
		sharedvars.EnableHeighfieldDesignOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedCheckShowsegmentgraph()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowSegmentGraphOn==FALSE)
	{
		m_chkShowSegmentGraphOn=TRUE;
		sharedvars.ShowSegmentGraphOn=true;
	}
	else
	{
		m_chkShowSegmentGraphOn=FALSE;
		sharedvars.ShowSegmentGraphOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg1::OnBnClickedButtonSavedesignelems()
{
	// TODO: Add your control notification handler code here
	CFile f;

	char strFilter[] = { "Text Files (*.txt)|*.txt|All Files (*.*)|*.*||" };

	CFileDialog FileDlg(FALSE, ".txt", NULL, 0, strFilter);
	if( FileDlg.DoModal() == IDOK )
	{
		CString f2 = FileDlg.GetFileName();
		const char *filename = f2;
 
		char filename2[256];
		strcpy(filename2, filename);
		save_tenField_elems(filename2);
	}
}

void MyDlg1::OnBnClickedButtonLoaddesignelems()
{
	// TODO: Add your control notification handler code here
	CFile f;

	char strFilter[] = { "Text Files (*.txt)|*.txt|All Files (*.*)|*.*||" };

	CFileDialog FileDlg(TRUE, ".txt", NULL, 0, strFilter);
	if( FileDlg.DoModal() == IDOK )
	{

		clock_t start, finish;
		start = clock();

		CString f2 = FileDlg.GetFileName();

		const char *filename = f2;
		
		char filename2[256];
		strcpy(filename2, filename);
		load_tenField_elems(filename2);
		
		cal_tensorvals_quad_inReg();
		cal_all_eigenvecs_quad();

		init_degpts();
		/*render alpha map*/
		render_alpha_map_quad(false);
		render_alpha_map_quad(true);

		locate_degpts_cells_tranvec_quad();

		num_shapecontrol_pts = 0;
	    num_curvepts_output = 0;

		//num_cycleedges = 0;
		//num_triangles_designcurve = 0;
		//num_innertriangles = 0;

		//InnerBoundary.num = 0;
		//OuterBoundary.num = 0;

		/////
		
		finish = clock();
		
		UpdateData(FALSE);
	}
}
