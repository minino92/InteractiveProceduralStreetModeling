// StreetEditPanel.cpp : implementation file
//

#include "stdafx.h"
#include "IBFV2D.h"
#include "StreetEditPanel.h"
#include ".\streeteditpanel.h"

#include "MyDlg1.h"
#include "MyDlg2.h"
#include "StreetNetPanel.h"

#include "regionsmooth_quad.h"
#include "GlView.h"
#include ".\glview.h"
#include "shareinterfacevars.h"
#include "VFDataStructure.h"
#include "tensoranalysis.h"
#include "tensorvis.h"

#include "LocalRegionEdit.h"

extern MyDlg1 *g_mydlg1;
extern MyDlg2 *g_mydlg2;
extern StreetNetPanel *g_streetnetpanel;
extern CGlView *g_pclGlView;

extern QuadMesh *quadmesh;

extern SharedInterfaceVars sharedvars;

extern int Num_SmoothRegionpoints;

extern double BoundRegionWidth;

extern bool please_comb_prefield;


extern void jitter_roads_inReg(double curScale);
extern void jitter_tens_inReg(double curScale, double curAmp);

/*  noise control global variable  */
#include "NoiseCtrlPanel.h"
extern NoiseCtrlPanel *g_noise_ctrlDlg;

extern double NoiseFreq;
extern double NoiseAmp;

extern double majorDensity ;
extern double minorDensity ;
extern double mintenline_length;

extern double map_xrang ;  // default is 10km
extern double map_yrang ;  // default is 10km

// StreetEditPanel dialog

IMPLEMENT_DYNAMIC(StreetEditPanel, CDialog)
StreetEditPanel::StreetEditPanel(CWnd* pParent /*=NULL*/)
	: CDialog(StreetEditPanel::IDD, pParent)
	, m_chkSelStreetRegToEditOn(FALSE)
{
}

StreetEditPanel::~StreetEditPanel()
{
}

void StreetEditPanel::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Check(pDX, IDC_CHECK_SELSTREETREGTOEDIT, m_chkSelStreetRegToEditOn);
}


BEGIN_MESSAGE_MAP(StreetEditPanel, CDialog)
	ON_BN_CLICKED(IDC_CHECK_SELSTREETREGTOEDIT, OnBnClickedCheckSelstreetregtoedit)
	ON_BN_CLICKED(IDC_BUTTON_REMOVESELREGION, OnBnClickedButtonRemoveselregion)
	ON_BN_CLICKED(IDC_BUTTON_TRACEINSIDEREGION, OnBnClickedButtonTraceinsideregion)
	ON_BN_CLICKED(IDC_BUTTON_COMSTREETNETINREG, OnBnClickedButtonComstreetnetinreg)
	ON_BN_CLICKED(IDC_BUTTON_ADDNOISETOROADSINREG, OnBnClickedButtonAddnoisetoroadsinreg)
	ON_BN_CLICKED(IDC_BUTTON_NOISESETTING, OnBnClickedButtonNoisesetting)
	ON_BN_CLICKED(IDC_BUTTON_ADDNOISETOFIELDINREG, OnBnClickedButtonAddnoisetofieldinreg)
END_MESSAGE_MAP()


// StreetEditPanel message handlers

void StreetEditPanel::OnBnClickedCheckSelstreetregtoedit()
{
	// TODO: Add your control notification handler code here
	if(m_chkSelStreetRegToEditOn==FALSE)
	{
		m_chkSelStreetRegToEditOn=TRUE;
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
		m_chkSelStreetRegToEditOn=FALSE;
		sharedvars.SelStreetRegToEditOn=false;
		g_pclGlView->SmoothOn = 0;              ////SELECTION MODE
		g_pclGlView->PickPointOn = 0;
		g_pclGlView->DisplaySmoothRegionOn = 0;
		sharedvars.RegionSmoothOn=false;
		please_comb_prefield=false;
	}
}

void StreetEditPanel::OnBnClickedButtonRemoveselregion()
{
	// TODO: Add your control notification handler code here
	if(sharedvars.BrushInterfaceOn)
	{
		/*  we need to determine the region according to the level set distance  */
		compute_level_set_contour(BoundRegionWidth-10*quadmesh->xinterval);
	}
	else if(sharedvars.SelStreetRegToEditOn)
	{
		remove_street_in_reg();
		mark_inner_verts_with_new_RegIndex();
		update_RegIndex_inner_and_boundary_cells();
	}
}

void StreetEditPanel::OnBnClickedButtonTraceinsideregion()
{
	// TODO: Add your control notification handler code here
	g_streetnetpanel->UpdateData();

	minorDensity = (g_streetnetpanel->m_edMinTenLineDensity/map_xrang*(quadmesh->xend-quadmesh->xstart))/
		quadmesh->xinterval;
	majorDensity = (g_streetnetpanel->m_edMajTenLineDensity/map_xrang*(quadmesh->xend-quadmesh->xstart))/
		quadmesh->xinterval;
	mintenline_length = (g_streetnetpanel->m_edMinTensorLineLength/map_xrang*(quadmesh->xend-quadmesh->xstart))
		/quadmesh->xinterval;
	replace_inReg();
}

void StreetEditPanel::OnBnClickedButtonComstreetnetinreg()
{
	// TODO: Add your control notification handler code here
	construct_sub_graph_inReg();

	//merge_subgraph_to_streetnet(false,false);
	//merge_subgraph_to_streetnet(false,true);
	//merge_subgraph_to_streetnet(true,false);
	//merge_subgraph_to_streetnet(true,true);

	/*  merge subgraph to the original street network still has problem 1/19/2008 */
	//merge_subgraph_to_streetnet_2(false,false);
	//merge_subgraph_to_streetnet_2(false,true);
	//merge_subgraph_to_streetnet_2(true,false);
	//merge_subgraph_to_streetnet_2(true,true);
	update_street_network();

	/**/
	int i;
	extern StreetNet *streetnet;
	for(i=0;i<streetnet->edgelist->nedges;i++)
	{
		if(streetnet->edgelist->edges[i]->node_index1>=streetnet->nodelist->nelems
			||streetnet->edgelist->edges[i]->node_index2>=streetnet->nodelist->nelems)
		{
			int test=0;
		}
	}
}

void StreetEditPanel::OnBnClickedButtonAddnoisetoroadsinreg()
{
	// TODO: Add your control notification handler code here
	init_quad_regionsmooth();
	find_boundarycells();
	mark_boundVerts();
	find_innerVerts();
	find_inner_cells();

	jitter_roads_inReg(NoiseFreq/*500.*/);
}

void StreetEditPanel::OnBnClickedButtonNoisesetting()
{
	// TODO: Add your control notification handler code here
	if(g_noise_ctrlDlg == NULL)
		return;

	CRect l_rectClient;
	CRect l_rectWnd;

	GetClientRect(l_rectClient);
	GetWindowRect(l_rectWnd);
	GetParent()->ScreenToClient(l_rectWnd);

	g_noise_ctrlDlg->GetWindowRect(l_rectClient);

	g_noise_ctrlDlg->SetWindowPos(&wndTop,l_rectWnd.left+50,l_rectWnd.top+500,l_rectClient.Width(),
		l_rectClient.Height(),SWP_SHOWWINDOW);
	g_noise_ctrlDlg->ShowWindow(SW_SHOW);
}

void StreetEditPanel::OnBnClickedButtonAddnoisetofieldinreg()
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
