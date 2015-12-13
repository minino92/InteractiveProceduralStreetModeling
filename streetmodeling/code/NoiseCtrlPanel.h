#pragma once
#include "afxcmn.h"


// NoiseCtrlPanel dialog

class NoiseCtrlPanel : public CDialog
{
	DECLARE_DYNAMIC(NoiseCtrlPanel)

public:
	NoiseCtrlPanel(CWnd* pParent = NULL);   // standard constructor
	virtual ~NoiseCtrlPanel();

// Dialog Data
	enum { IDD = IDD_DIALOG_NOISECTRLPANEL };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	CSliderCtrl m_ctrlNoiseFreq;
	CSliderCtrl m_ctrlNoiseAmp;
	double m_edNoiseFreq;
	double m_edNoiseAmp;
	afx_msg void OnNMCustomdrawSliderNoisefreq(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnNMCustomdrawSliderNoiseamp(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnEnChangeEditNoisefreq();
	afx_msg void OnEnChangeEditNoiseamp();
	virtual BOOL OnInitDialog();
};
