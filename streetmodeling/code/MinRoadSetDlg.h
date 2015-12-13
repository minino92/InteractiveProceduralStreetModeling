#pragma once


// MinRoadSetDlg dialog

class MinRoadSetDlg : public CDialog
{
	DECLARE_DYNAMIC(MinRoadSetDlg)

public:
	MinRoadSetDlg(CWnd* pParent = NULL);   // standard constructor
	virtual ~MinRoadSetDlg();

// Dialog Data
	enum { IDD = IDD_DIALOG_MINROADSETTING };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	BOOL m_chkAllowMinCloseToMajOn;
	afx_msg void OnBnClickedCheckMinclosetomaj();
	BOOL m_chkRemDeadEndsTraceOn;
	afx_msg void OnBnClickedCheckAutoremovedeadends();
	BOOL m_chkUseNewMedRemDeadEndOn;
	afx_msg void OnBnClickedCheckUsenewmedremdeadends();
};
