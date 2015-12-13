#pragma once


// MajRoadSettingDlg dialog

class MajRoadSettingDlg : public CDialog
{
	DECLARE_DYNAMIC(MajRoadSettingDlg)

public:
	MajRoadSettingDlg(CWnd* pParent = NULL);   // standard constructor
	virtual ~MajRoadSettingDlg();

// Dialog Data
	enum { IDD = IDD_DIALOG_MAJROADSETTING };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	double m_edMinDistCrossRiver;
	BOOL m_chkAllowMajRoadCrossRiverOn;
	afx_msg void OnBnClickedCheckCrossriver();
	double m_edMinAngCrossRiver;
	double m_edMaxDistFollowBoundary;
	afx_msg void OnEnChangeEditDiscrossrivers();
	afx_msg void OnEnChangeEditDisfollowboundaries();
	afx_msg void OnEnChangeEditAngcrossingrivers();
	BOOL m_chkAllowMajRoadFollowBoundaryOn;
	afx_msg void OnBnClickedCheckMajroadfollowboundaries();
	BOOL m_chkShowMajRoadNetworkOn;
	afx_msg void OnBnClickedCheckShowmajroadnetwork();
	BOOL m_chkAllowCrossSingularitiesOn;
	afx_msg void OnBnClickedCheckAllowcrossdegpts();
	afx_msg void OnBnClickedButtonCleanmajraoddeadends();
	BOOL m_chkJobardMethodOn;
	afx_msg void OnBnClickedCheckJobardmethod();
};
