// MinRoadSetDlg.cpp : implementation file
//

#include "stdafx.h"
#include "IBFV2D.h"
#include "MinRoadSetDlg.h"
#include ".\minroadsetdlg.h"

#include "shareinterfacevars.h"
extern SharedInterfaceVars sharedvars;

// MinRoadSetDlg dialog

IMPLEMENT_DYNAMIC(MinRoadSetDlg, CDialog)
MinRoadSetDlg::MinRoadSetDlg(CWnd* pParent /*=NULL*/)
	: CDialog(MinRoadSetDlg::IDD, pParent)
	, m_chkAllowMinCloseToMajOn(TRUE)
	, m_chkRemDeadEndsTraceOn(FALSE)
	, m_chkUseNewMedRemDeadEndOn(FALSE)
{
}

MinRoadSetDlg::~MinRoadSetDlg()
{
}

void MinRoadSetDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Check(pDX, IDC_CHECK_MINCLOSETOMAJ, m_chkAllowMinCloseToMajOn);
	DDX_Check(pDX, IDC_CHECK_AUTOREMOVEDEADENDS, m_chkRemDeadEndsTraceOn);
	DDX_Check(pDX, IDC_CHECK_USENEWMEDREMDEADENDS, m_chkUseNewMedRemDeadEndOn);
}


BEGIN_MESSAGE_MAP(MinRoadSetDlg, CDialog)
	ON_BN_CLICKED(IDC_CHECK_MINCLOSETOMAJ, OnBnClickedCheckMinclosetomaj)
	ON_BN_CLICKED(IDC_CHECK_AUTOREMOVEDEADENDS, OnBnClickedCheckAutoremovedeadends)
	ON_BN_CLICKED(IDC_CHECK_USENEWMEDREMDEADENDS, OnBnClickedCheckUsenewmedremdeadends)
END_MESSAGE_MAP()


// MinRoadSetDlg message handlers

void MinRoadSetDlg::OnBnClickedCheckMinclosetomaj()
{
	// TODO: Add your control notification handler code here
	if(m_chkAllowMinCloseToMajOn==FALSE)
	{
		m_chkAllowMinCloseToMajOn=TRUE;
		sharedvars.AllowMinCloseToMajOn=true;
	}
	else
	{
		m_chkAllowMinCloseToMajOn=FALSE;
		sharedvars.AllowMinCloseToMajOn=false;
	}
	UpdateData(FALSE);
}

void MinRoadSetDlg::OnBnClickedCheckAutoremovedeadends()
{
	// TODO: Add your control notification handler code here
	if(m_chkRemDeadEndsTraceOn==FALSE)
	{
		m_chkRemDeadEndsTraceOn=TRUE;
		sharedvars.RemDeadEndsTraceOn=true;
	}
	else
	{
		m_chkRemDeadEndsTraceOn=FALSE;
		sharedvars.RemDeadEndsTraceOn=false;
	}
	UpdateData(FALSE);
}

void MinRoadSetDlg::OnBnClickedCheckUsenewmedremdeadends()
{
	// TODO: Add your control notification handler code here
	if(m_chkUseNewMedRemDeadEndOn==FALSE)
	{
		m_chkUseNewMedRemDeadEndOn=TRUE;
		sharedvars.UseNewMedRemDeadEndOn=true;
	}
	else
	{
		m_chkUseNewMedRemDeadEndOn=FALSE;
		sharedvars.UseNewMedRemDeadEndOn=false;
	}
	UpdateData(FALSE);
}
