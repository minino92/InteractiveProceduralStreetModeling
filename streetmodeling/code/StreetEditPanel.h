#pragma once


// StreetEditPanel dialog
	//NoiseCtrlPanel *noise_ctrlDlg;
	//int noise_ctrlDlgID;

class StreetEditPanel : public CDialog
{
	DECLARE_DYNAMIC(StreetEditPanel)

public:
	StreetEditPanel(CWnd* pParent = NULL);   // standard constructor
	virtual ~StreetEditPanel();

// Dialog Data
	enum { IDD = IDD_DIALOG_STREETEDITPANEL };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	BOOL m_chkSelStreetRegToEditOn;
	afx_msg void OnBnClickedCheckSelstreetregtoedit();
	afx_msg void OnBnClickedButtonRemoveselregion();
	afx_msg void OnBnClickedButtonTraceinsideregion();
	afx_msg void OnBnClickedButtonComstreetnetinreg();
	afx_msg void OnBnClickedButtonAddnoisetoroadsinreg();
	afx_msg void OnBnClickedButtonNoisesetting();
	afx_msg void OnBnClickedButtonAddnoisetofieldinreg();
};
