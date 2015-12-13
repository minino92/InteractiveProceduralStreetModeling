// MajRoadSettingDlg.cpp : implementation file
//

#include "stdafx.h"
#include "IBFV2D.h"
#include "MajRoadSettingDlg.h"
#include ".\majroadsettingdlg.h"

#include "VFDataStructure.h"

#include "EvenlyStreamlines.h"

#include "shareinterfacevars.h"
extern SharedInterfaceVars sharedvars;


//double MinDistCrossRiver = 150;
double MinDistCrossRiver = 500;
double MinAngCrossRiver = 30;
double MaxDistFollowBoundary = 600;

extern double map_xrang;  // default is 10km
extern double map_yrang;  // default is 10km
extern QuadMesh *quadmesh;

// MajRoadSettingDlg dialog

IMPLEMENT_DYNAMIC(MajRoadSettingDlg, CDialog)
MajRoadSettingDlg::MajRoadSettingDlg(CWnd* pParent /*=NULL*/)
	: CDialog(MajRoadSettingDlg::IDD, pParent)
	//, m_edMinDistCrossRiver(150)
	, m_edMinDistCrossRiver(500)
	, m_chkAllowMajRoadCrossRiverOn(FALSE)
	, m_edMinAngCrossRiver(30.)
	, m_edMaxDistFollowBoundary(600)
	, m_chkAllowMajRoadFollowBoundaryOn(FALSE)
	, m_chkShowMajRoadNetworkOn(FALSE)
	, m_chkAllowCrossSingularitiesOn(FALSE)
	, m_chkJobardMethodOn(FALSE)
{
}

MajRoadSettingDlg::~MajRoadSettingDlg()
{
}

void MajRoadSettingDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT_DISCROSSRIVERS, m_edMinDistCrossRiver);
	DDX_Check(pDX, IDC_CHECK_CROSSRIVER, m_chkAllowMajRoadCrossRiverOn);
	DDX_Text(pDX, IDC_EDIT_ANGCROSSINGRIVERS, m_edMinAngCrossRiver);
	DDX_Text(pDX, IDC_EDIT_DISFOLLOWBOUNDARIES, m_edMaxDistFollowBoundary);
	DDX_Check(pDX, IDC_CHECK_MAJROADFOLLOWBOUNDARIES, m_chkAllowMajRoadFollowBoundaryOn);
	DDX_Check(pDX, IDC_CHECK_SHOWMAJROADNETWORK, m_chkShowMajRoadNetworkOn);
	DDX_Check(pDX, IDC_CHECK_ALLOWCROSSDEGPTS, m_chkAllowCrossSingularitiesOn);
	DDX_Check(pDX, IDC_CHECK_JOBARDMETHOD, m_chkJobardMethodOn);
}


BEGIN_MESSAGE_MAP(MajRoadSettingDlg, CDialog)
	ON_BN_CLICKED(IDC_CHECK_CROSSRIVER, OnBnClickedCheckCrossriver)
	ON_EN_CHANGE(IDC_EDIT_DISCROSSRIVERS, OnEnChangeEditDiscrossrivers)
	ON_EN_CHANGE(IDC_EDIT_DISFOLLOWBOUNDARIES, OnEnChangeEditDisfollowboundaries)
	ON_EN_CHANGE(IDC_EDIT_ANGCROSSINGRIVERS, OnEnChangeEditAngcrossingrivers)
	ON_BN_CLICKED(IDC_CHECK_MAJROADFOLLOWBOUNDARIES, OnBnClickedCheckMajroadfollowboundaries)
	ON_BN_CLICKED(IDC_CHECK_SHOWMAJROADNETWORK, OnBnClickedCheckShowmajroadnetwork)
	ON_BN_CLICKED(IDC_CHECK_ALLOWCROSSDEGPTS, OnBnClickedCheckAllowcrossdegpts)
	ON_BN_CLICKED(IDC_BUTTON_CLEANMAJRAODDEADENDS, OnBnClickedButtonCleanmajraoddeadends)
	ON_BN_CLICKED(IDC_CHECK_JOBARDMETHOD, OnBnClickedCheckJobardmethod)
END_MESSAGE_MAP()


// MajRoadSettingDlg message handlers

void MajRoadSettingDlg::OnBnClickedCheckCrossriver()
{
	// TODO: Add your control notification handler code here
	if(m_chkAllowMajRoadCrossRiverOn==FALSE)
	{
		m_chkAllowMajRoadCrossRiverOn=TRUE;
		sharedvars.AllowMajRoadCrossRiverOn=true;
		UpdateData();
		MinDistCrossRiver=m_edMinDistCrossRiver;
		MinDistCrossRiver=(MinDistCrossRiver/map_xrang*(quadmesh->xend-quadmesh->xstart))/
			quadmesh->xinterval;
		MinAngCrossRiver=m_edMinAngCrossRiver;
		MinAngCrossRiver=MinAngCrossRiver/180*M_PI;
	}
	else
	{
		m_chkAllowMajRoadCrossRiverOn=FALSE;
		sharedvars.AllowMajRoadCrossRiverOn=false;
	}
	UpdateData(FALSE);
}

void MajRoadSettingDlg::OnEnChangeEditDiscrossrivers()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData();
	MinDistCrossRiver=m_edMinDistCrossRiver;
	MinDistCrossRiver=(MinDistCrossRiver/map_xrang*(quadmesh->xend-quadmesh->xstart))/
		quadmesh->xinterval;
}

void MajRoadSettingDlg::OnEnChangeEditDisfollowboundaries()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData();
	MaxDistFollowBoundary=m_edMaxDistFollowBoundary;
	MaxDistFollowBoundary=(MaxDistFollowBoundary/map_xrang*(quadmesh->xend-quadmesh->xstart))/
		quadmesh->xinterval;
}

void MajRoadSettingDlg::OnEnChangeEditAngcrossingrivers()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData();
	MinAngCrossRiver=m_edMinAngCrossRiver;
	if(MinAngCrossRiver<0||MinAngCrossRiver>90.)
	{
		MessageBox("Please set the angle between 0 and 90 degree", MB_OK);
		MinAngCrossRiver=m_edMinAngCrossRiver=30.;
		UpdateData(FALSE);
	}
	MinAngCrossRiver=MinAngCrossRiver/180*M_PI;
}

void MajRoadSettingDlg::OnBnClickedCheckMajroadfollowboundaries()
{
	// TODO: Add your control notification handler code here
	if(m_chkAllowMajRoadFollowBoundaryOn==FALSE)
	{
		m_chkAllowMajRoadFollowBoundaryOn=TRUE;
		sharedvars.AllowMajRoadFollowBoundaryOn=true;
		UpdateData();
		MaxDistFollowBoundary=m_edMaxDistFollowBoundary;
		MaxDistFollowBoundary=(MaxDistFollowBoundary/map_xrang*(quadmesh->xend-quadmesh->xstart))/
			quadmesh->xinterval;
	}
	else
	{
		m_chkAllowMajRoadFollowBoundaryOn=FALSE;
		sharedvars.AllowMajRoadFollowBoundaryOn=false;
	}
	UpdateData(FALSE);
}

void MajRoadSettingDlg::OnBnClickedCheckShowmajroadnetwork()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowMajRoadNetworkOn==FALSE)
	{
		m_chkShowMajRoadNetworkOn=TRUE;
		sharedvars.ShowMajRoadNetworkOn=true;
	}
	else
	{
		m_chkShowMajRoadNetworkOn=FALSE;
		sharedvars.ShowMajRoadNetworkOn=false;
	}
	UpdateData(FALSE);
}

void MajRoadSettingDlg::OnBnClickedCheckAllowcrossdegpts()
{
	// TODO: Add your control notification handler code here
	if(m_chkAllowCrossSingularitiesOn==FALSE)
	{
		m_chkAllowCrossSingularitiesOn=TRUE;
		sharedvars.AllowCrossSingularitiesOn=true;
	}
	else
	{
		m_chkAllowCrossSingularitiesOn=FALSE;
		sharedvars.AllowCrossSingularitiesOn=false;
	}
	UpdateData(FALSE);
}

void MajRoadSettingDlg::OnBnClickedButtonCleanmajraoddeadends()
{
	// TODO: Add your control notification handler code here
	connect_majRoads_postproc();
}

void MajRoadSettingDlg::OnBnClickedCheckJobardmethod()
{
	// TODO: Add your control notification handler code here
	if(m_chkJobardMethodOn==FALSE)
	{
		m_chkJobardMethodOn=TRUE;
		sharedvars.JobardMethodOn=true;
	}
	else
	{
		m_chkJobardMethodOn=FALSE;
		sharedvars.JobardMethodOn=false;
	}
	UpdateData(FALSE);
}
