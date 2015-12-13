// NoiseCtrlPanel.cpp : implementation file
//

#include "stdafx.h"
#include "IBFV2D.h"
#include "NoiseCtrlPanel.h"
#include ".\noisectrlpanel.h"


// NoiseCtrlPanel dialog

double NoiseFreq=10.;
double NoiseAmp=1.;

IMPLEMENT_DYNAMIC(NoiseCtrlPanel, CDialog)
NoiseCtrlPanel::NoiseCtrlPanel(CWnd* pParent /*=NULL*/)
	: CDialog(NoiseCtrlPanel::IDD, pParent)
	, m_edNoiseFreq(10.)
	, m_edNoiseAmp(1.)
{
}

NoiseCtrlPanel::~NoiseCtrlPanel()
{
}

void NoiseCtrlPanel::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_SLIDER_NOISEFREQ, m_ctrlNoiseFreq);
	DDX_Control(pDX, IDC_SLIDER_NOISEAMP, m_ctrlNoiseAmp);
	DDX_Text(pDX, IDC_EDIT_NOISEFREQ, m_edNoiseFreq);
	DDX_Text(pDX, IDC_EDIT_NOISEAMP, m_edNoiseAmp);
}


BEGIN_MESSAGE_MAP(NoiseCtrlPanel, CDialog)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_NOISEFREQ, OnNMCustomdrawSliderNoisefreq)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER_NOISEAMP, OnNMCustomdrawSliderNoiseamp)
	ON_EN_CHANGE(IDC_EDIT_NOISEFREQ, OnEnChangeEditNoisefreq)
	ON_EN_CHANGE(IDC_EDIT_NOISEAMP, OnEnChangeEditNoiseamp)
END_MESSAGE_MAP()


// NoiseCtrlPanel message handlers

void NoiseCtrlPanel::OnNMCustomdrawSliderNoisefreq(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	// TODO: Add your control notification handler code here
	NoiseFreq=(double)m_ctrlNoiseFreq.GetPos();
	m_edNoiseFreq=NoiseFreq;
	UpdateData(FALSE);
	*pResult = 0;
}

void NoiseCtrlPanel::OnNMCustomdrawSliderNoiseamp(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	// TODO: Add your control notification handler code here
	NoiseAmp=(double)m_ctrlNoiseAmp.GetPos()/100.;
	m_edNoiseAmp=NoiseAmp;
	UpdateData(FALSE);
	*pResult = 0;
}

void NoiseCtrlPanel::OnEnChangeEditNoisefreq()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData();
	NoiseFreq=m_edNoiseFreq;
	m_ctrlNoiseFreq.SetPos(int(NoiseFreq));
	UpdateData(FALSE);
}

void NoiseCtrlPanel::OnEnChangeEditNoiseamp()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData();
	NoiseAmp=m_edNoiseAmp;
	m_ctrlNoiseAmp.SetPos(int(NoiseAmp*100.));
	UpdateData(FALSE);
}

BOOL NoiseCtrlPanel::OnInitDialog()
{
	CDialog::OnInitDialog();

	// TODO:  Add extra initialization here

	m_ctrlNoiseFreq.SetRange(1, 500);
	m_ctrlNoiseFreq.SetPos(10);
	m_ctrlNoiseFreq.SetTicFreq(1);
	m_ctrlNoiseFreq.SetTic(0);

	m_ctrlNoiseAmp.SetRange(1, 500);
	m_ctrlNoiseAmp.SetPos(100);
	m_ctrlNoiseAmp.SetTicFreq(1);
	m_ctrlNoiseAmp.SetTic(0);

	UpdateData(FALSE);
	return TRUE;  // return TRUE unless you set the focus to a control
	// EXCEPTION: OCX Property Pages should return FALSE
}
