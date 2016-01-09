#pragma once
#include "afxwin.h"
#include "afxcmn.h"


// MyDlg1 dialog

class MyDlg1 : public CDialog
{
	DECLARE_DYNAMIC(MyDlg1)

public:
	MyDlg1(CWnd* pParent = NULL);   // standard constructor
	virtual ~MyDlg1();
	
	void Reset();

// Dialog Data
	enum { IDD = IDD_DIALOG1 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	//BOOL m_rdAddVFSingularElem;
	BOOL m_chkTensorDesignOn;
	afx_msg void OnBnClickedCheckTensordesignon();
	BOOL m_chkMoveElemOn;
	afx_msg void OnBnClickedCheckMoveelem();
	BOOL m_chkEditElementOn;
	afx_msg void OnBnClickedCheckEditanelem();
	BOOL m_chkRemoveElemOn;
	afx_msg void OnBnClickedCheckRemoveelem();
	BOOL m_chkBrushInterfaceOn;
	afx_msg void OnBnClickedCheckActivatebrush();
	BOOL m_chkCombineWithOtherPartsOn;
	afx_msg void OnBnClickedCheckCombinewithprevioustensor();
	afx_msg void OnBnClickedButtonGentenfrombrush();
	double m_edBrushWidth;
	afx_msg void OnBnClickedButtonSavebrush();
	afx_msg void OnBnClickedButtonLoadbrush();
	BOOL m_rdAddVFSingularElem;
	afx_msg void OnBnClickedRadioAddawedge();
	afx_msg void OnBnClickedRadioAddatrisector();
	afx_msg void OnBnClickedRadioAddanode();
	afx_msg void OnBnClickedRadioAddacenter();
	afx_msg void OnBnClickedRadioAddasaddle();
	afx_msg void OnBnClickedRadioAddaregular();
	BOOL m_chkRegionSmoothOn;
	afx_msg void OnBnClickedCheckRegionsmoothon();
	afx_msg void OnBnClickedButtonSmoothonce();
	afx_msg void OnBnClickedButtonSavefieldperver();
	afx_msg void OnBnClickedButtonLoadfieldpervert();
	BOOL m_chkEnableSketchBasedDesign;
	afx_msg void OnBnClickedCheckUsedrawsketchon();
	int m_edNumBrushes;
	CComboBox m_ctrlcomboxDesignRegionIndex;
	virtual BOOL OnInitDialog();
	CString m_selRegion;
	afx_msg void OnCbnSelchangeComboDiffregions();
	BOOL m_chkShowSketchesOn;
	afx_msg void OnBnClickedCheckDisplaysketches();
	afx_msg void OnBnClickedButtonSketchbasedtengen();
	CSliderCtrl m_ctrlSliderBrushRegionWidth;
	afx_msg void OnNMCustomdrawSlider1(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnBnClickedButtonClearsketches();
	afx_msg void OnEnChangeEdit1();
	afx_msg void OnBnClickedButtonSavesketches();
	afx_msg void OnBnClickedButtonLoadsketches();
	afx_msg void OnBnClickedButtonSegmentregion();
	BOOL m_rdSketchMajorHighways;
	afx_msg void OnBnClickedRadioSketchmajorroads();
	afx_msg void OnBnClickedRadioSketchhighway();
	BOOL m_chkUseBoundsAsRoadsOn;
	afx_msg void OnBnClickedCheckUseboundaryasroads();
	afx_msg void OnBnClickedButtonClearthefield();
	BOOL m_chkUseMajRoadsAsSketchesOn;
	afx_msg void OnBnClickedCheckUsemajroadsassketches();
	BOOL m_rdScalarSingularElem;
	afx_msg void OnBnClickedRadioScalarmaximum();
	afx_msg void OnBnClickedRadioScalarminimum();
	BOOL m_chkEnableHeighfieldDesignOn;
	afx_msg void OnBnClickedCheckDesignscalarfield();
	BOOL m_chkShowSegmentGraphOn;
	afx_msg void OnBnClickedCheckShowsegmentgraph();
	afx_msg void OnBnClickedButtonSavedesignelems();
	afx_msg void OnBnClickedButtonLoaddesignelems();
	CButton Trace;
};
