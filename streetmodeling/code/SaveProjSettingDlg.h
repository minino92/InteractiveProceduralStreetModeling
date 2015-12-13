#pragma once


// SaveProjSettingDlg dialog

class SaveProjSettingDlg : public CDialog
{
	DECLARE_DYNAMIC(SaveProjSettingDlg)

public:
	SaveProjSettingDlg(CWnd* pParent = NULL);   // standard constructor
	virtual ~SaveProjSettingDlg();

// Dialog Data
	enum { IDD = IDD_DIALOG_SAVEPROJSETTING };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButtonSaveproject();
	BOOL m_rdSaveProjDesignElem;
	BOOL m_rdSaveProjSketches;
	BOOL m_rdSaveProjBrushes;
	BOOL m_rdSaveProjMajRoadNetwork;
	BOOL m_rdSaveProjOtherSetting;
	afx_msg void OnBnClickedRadioSavedesignelemnone();
	afx_msg void OnBnClickedRadioSaveaselems();
	afx_msg void OnBnClickedRadioSaveasperverfield();
	afx_msg void OnBnClickedRadioNotsavesketches();
	afx_msg void OnBnClickedRadioSavesketches();
	afx_msg void OnBnClickedRadioNotsavebrushes();
	afx_msg void OnBnClickedRadioSavebrushes();
	afx_msg void OnBnClickedRadioNotsavemajroads();
	afx_msg void OnBnClickedRadioSavemajroadsastenlines();
	afx_msg void OnBnClickedRadioSavemajroadsasgraph();
	afx_msg void OnBnClickedRadioNotsaveothersetting();
	afx_msg void OnBnClickedRadioSaveothersetting();
	BOOL m_rdSaveProjStreetNetwork;
	afx_msg void OnBnClickedRadioNotsavestreetnetwork();
	afx_msg void OnBnClickedRadioSavestreetnetwork();
};
