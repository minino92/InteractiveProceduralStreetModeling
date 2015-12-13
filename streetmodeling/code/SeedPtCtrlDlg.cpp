// SeedPtCtrlDlg.cpp : implementation file
//

#include "stdafx.h"
#include "IBFV2D.h"
#include "SeedPtCtrlDlg.h"
#include ".\seedptctrldlg.h"

#include "VFDataStructure.h"

#include "shareinterfacevars.h"
extern SharedInterfaceVars sharedvars;


extern double map_xrang;  // default is 10km
extern double map_yrang;  // default is 10km
extern QuadMesh *quadmesh;

// SeedPtCtrlDlg dialog
double DistAwayBoundaries=1.;

IMPLEMENT_DYNAMIC(SeedPtCtrlDlg, CDialog)
SeedPtCtrlDlg::SeedPtCtrlDlg(CWnd* pParent /*=NULL*/)
	: CDialog(SeedPtCtrlDlg::IDD, pParent)
	//, m_edDisAwayFromBoundaries(100.)
	, m_edSeedDistAwayBoundary(100.)
{
}

SeedPtCtrlDlg::~SeedPtCtrlDlg()
{
}

void SeedPtCtrlDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//DDX_Text(pDX, IDC_EDIT1, m_edDisAwayFromBoundaries);
	DDX_Text(pDX, IDC_EDIT1_SEEDDISTAWAYBOUNDARY, m_edSeedDistAwayBoundary);
}


BEGIN_MESSAGE_MAP(SeedPtCtrlDlg, CDialog)
	ON_EN_CHANGE(IDC_EDIT1_SEEDDISTAWAYBOUNDARY, OnEnChangeEdit1Seeddistawayboundary)
END_MESSAGE_MAP()


// SeedPtCtrlDlg message handlers


void SeedPtCtrlDlg::OnEnChangeEdit1Seeddistawayboundary()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData();
	DistAwayBoundaries=(m_edSeedDistAwayBoundary/map_xrang*(quadmesh->xend-quadmesh->xstart))/
		quadmesh->xinterval;
}
