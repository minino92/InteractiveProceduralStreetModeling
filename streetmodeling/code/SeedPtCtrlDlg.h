#pragma once


// SeedPtCtrlDlg dialog

class SeedPtCtrlDlg : public CDialog
{
	DECLARE_DYNAMIC(SeedPtCtrlDlg)

public:
	SeedPtCtrlDlg(CWnd* pParent = NULL);   // standard constructor
	virtual ~SeedPtCtrlDlg();

// Dialog Data
	enum { IDD = IDD_DIALOG_SEEDPOINTCTRLDLG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	//double m_edDisAwayFromBoundaries;
	double m_edSeedDistAwayBoundary;
	afx_msg void OnEnChangeEdit1Seeddistawayboundary();
};
