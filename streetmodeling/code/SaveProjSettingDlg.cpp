// SaveProjSettingDlg.cpp : implementation file
//

#include "stdafx.h"
#include "IBFV2D.h"
#include "SaveProjSettingDlg.h"
#include ".\saveprojsettingdlg.h"

#include "SaveProjectSetting.h"
#include "shareinterfacevars.h"

extern SharedInterfaceVars sharedvars;

// SaveProjSettingDlg dialog

IMPLEMENT_DYNAMIC(SaveProjSettingDlg, CDialog)
SaveProjSettingDlg::SaveProjSettingDlg(CWnd* pParent /*=NULL*/)
	: CDialog(SaveProjSettingDlg::IDD, pParent)
	, m_rdSaveProjDesignElem(FALSE)
	, m_rdSaveProjSketches(FALSE)
	, m_rdSaveProjBrushes(FALSE)
	, m_rdSaveProjMajRoadNetwork(FALSE)
	, m_rdSaveProjOtherSetting(FALSE)
	, m_rdSaveProjStreetNetwork(FALSE)
{
}

SaveProjSettingDlg::~SaveProjSettingDlg()
{
}

void SaveProjSettingDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Radio(pDX, IDC_RADIO_SAVEDESIGNELEMNONE, m_rdSaveProjDesignElem);
	DDX_Radio(pDX, IDC_RADIO_NOTSAVESKETCHES, m_rdSaveProjSketches);
	DDX_Radio(pDX, IDC_RADIO_NOTSAVEBRUSHES, m_rdSaveProjBrushes);
	DDX_Radio(pDX, IDC_RADIO_NOTSAVEMAJROADS, m_rdSaveProjMajRoadNetwork);
	DDX_Radio(pDX, IDC_RADIO_NOTSAVEOTHERSETTING, m_rdSaveProjOtherSetting);
	DDX_Radio(pDX, IDC_RADIO_NOTSAVESTREETNETWORK, m_rdSaveProjStreetNetwork);
}


BEGIN_MESSAGE_MAP(SaveProjSettingDlg, CDialog)
	ON_BN_CLICKED(IDC_BUTTON_SAVEPROJECT, OnBnClickedButtonSaveproject)
	ON_BN_CLICKED(IDC_RADIO_SAVEDESIGNELEMNONE, OnBnClickedRadioSavedesignelemnone)
	ON_BN_CLICKED(IDC_RADIO_SAVEASELEMS, OnBnClickedRadioSaveaselems)
	ON_BN_CLICKED(IDC_RADIO_SAVEASPERVERFIELD, OnBnClickedRadioSaveasperverfield)
	ON_BN_CLICKED(IDC_RADIO_NOTSAVESKETCHES, OnBnClickedRadioNotsavesketches)
	ON_BN_CLICKED(IDC_RADIO_SAVESKETCHES, OnBnClickedRadioSavesketches)
	ON_BN_CLICKED(IDC_RADIO_NOTSAVEBRUSHES, OnBnClickedRadioNotsavebrushes)
	ON_BN_CLICKED(IDC_RADIO_SAVEBRUSHES, OnBnClickedRadioSavebrushes)
	ON_BN_CLICKED(IDC_RADIO_NOTSAVEMAJROADS, OnBnClickedRadioNotsavemajroads)
	ON_BN_CLICKED(IDC_RADIO_SAVEMAJROADSASTENLINES, OnBnClickedRadioSavemajroadsastenlines)
	ON_BN_CLICKED(IDC_RADIO_SAVEMAJROADSASGRAPH, OnBnClickedRadioSavemajroadsasgraph)
	ON_BN_CLICKED(IDC_RADIO_NOTSAVEOTHERSETTING, OnBnClickedRadioNotsaveothersetting)
	ON_BN_CLICKED(IDC_RADIO_SAVEOTHERSETTING, OnBnClickedRadioSaveothersetting)
	ON_BN_CLICKED(IDC_RADIO_NOTSAVESTREETNETWORK, OnBnClickedRadioNotsavestreetnetwork)
	ON_BN_CLICKED(IDC_RADIO_SAVESTREETNETWORK, OnBnClickedRadioSavestreetnetwork)
END_MESSAGE_MAP()


// SaveProjSettingDlg message handlers

void SaveProjSettingDlg::OnBnClickedButtonSaveproject()
{
	// TODO: Add your control notification handler code here
	CFile f;
	char filename[256];

	char strFilter[] = { "Text Files (*.spj)|*.spj"};

	CFileDialog FileDlg(FALSE, ".spj", NULL, 0, strFilter);
	if( FileDlg.DoModal() == IDOK )
	{

		CString f2 = FileDlg.GetFileName();
		const char *tempf = f2;
		strcpy(filename, f2);
	
		save_current_project(filename);
	}

}

void SaveProjSettingDlg::OnBnClickedRadioSavedesignelemnone()
{
	// TODO: Add your control notification handler code here
	m_rdSaveProjDesignElem=0;
	sharedvars.rdSaveProjDesignElem=0;
	UpdateData(FALSE);
}

void SaveProjSettingDlg::OnBnClickedRadioSaveaselems()
{
	// TODO: Add your control notification handler code here
	m_rdSaveProjDesignElem=1;
	sharedvars.rdSaveProjDesignElem=1;
	UpdateData(FALSE);
}

void SaveProjSettingDlg::OnBnClickedRadioSaveasperverfield()
{
	// TODO: Add your control notification handler code here
	m_rdSaveProjDesignElem=2;
	sharedvars.rdSaveProjDesignElem=2;
	UpdateData(FALSE);
}

void SaveProjSettingDlg::OnBnClickedRadioNotsavesketches()
{
	// TODO: Add your control notification handler code here
	m_rdSaveProjSketches=0;
	sharedvars.rdSaveProjSketches=0;
	UpdateData(FALSE);
}

void SaveProjSettingDlg::OnBnClickedRadioSavesketches()
{
	// TODO: Add your control notification handler code here
	m_rdSaveProjSketches=1;
	sharedvars.rdSaveProjSketches=1;
	UpdateData(FALSE);
}

void SaveProjSettingDlg::OnBnClickedRadioNotsavebrushes()
{
	// TODO: Add your control notification handler code here
	m_rdSaveProjBrushes=0;
	sharedvars.rdSaveProjBrushes=0;
	UpdateData(FALSE);
}

void SaveProjSettingDlg::OnBnClickedRadioSavebrushes()
{
	// TODO: Add your control notification handler code here
	m_rdSaveProjBrushes=1;
	sharedvars.rdSaveProjBrushes=1;
	UpdateData(FALSE);
}

void SaveProjSettingDlg::OnBnClickedRadioNotsavemajroads()
{
	// TODO: Add your control notification handler code here
	m_rdSaveProjMajRoadNetwork=0;
	sharedvars.rdSaveProjMajRoadNetwork=0;
	UpdateData(FALSE);
}

void SaveProjSettingDlg::OnBnClickedRadioSavemajroadsastenlines()
{
	// TODO: Add your control notification handler code here
	m_rdSaveProjMajRoadNetwork=1;
	sharedvars.rdSaveProjMajRoadNetwork=1;
	UpdateData(FALSE);
}

void SaveProjSettingDlg::OnBnClickedRadioSavemajroadsasgraph()
{
	// TODO: Add your control notification handler code here
	m_rdSaveProjMajRoadNetwork=2;
	sharedvars.rdSaveProjMajRoadNetwork=2;
	UpdateData(FALSE);
}

void SaveProjSettingDlg::OnBnClickedRadioNotsaveothersetting()
{
	// TODO: Add your control notification handler code here
	m_rdSaveProjOtherSetting=0;
	sharedvars.rdSaveProjOtherSetting=0;
	UpdateData(FALSE);
}

void SaveProjSettingDlg::OnBnClickedRadioSaveothersetting()
{
	// TODO: Add your control notification handler code here
	m_rdSaveProjOtherSetting=1;
	sharedvars.rdSaveProjOtherSetting=0;
	UpdateData(FALSE);
}

void SaveProjSettingDlg::OnBnClickedRadioNotsavestreetnetwork()
{
	// TODO: Add your control notification handler code here
	m_rdSaveProjStreetNetwork=0;
	UpdateData(FALSE);
}

void SaveProjSettingDlg::OnBnClickedRadioSavestreetnetwork()
{
	// TODO: Add your control notification handler code here
	m_rdSaveProjStreetNetwork=1;
	UpdateData(FALSE);
}
