#pragma once
#include "afxcmn.h"


// MyDlg2 dialog

class MyDlg2 : public CDialog
{
	DECLARE_DYNAMIC(MyDlg2)

public:
	MyDlg2(CWnd* pParent = NULL);   // standard constructor
	virtual ~MyDlg2();

	void Reset();

// Dialog Data
	enum { IDD = IDD_DIALOG2 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	BOOL m_chkShowSingularitiesOn;
	afx_msg void OnBnClickedCheckShowsingularity();
	BOOL m_chkShowRegElemOn;
	afx_msg void OnBnClickedCheckShowregularelem();
	BOOL m_chkShowTensorLinesOn;
	afx_msg void OnBnClickedCheckShowtensorlines();
	BOOL m_chkShowStreetGraphOn;
	afx_msg void OnBnClickedCheckShowstreetgraph();
	BOOL m_chkShowRoadMapOn;
	afx_msg void OnBnClickedCheckShowraodmap();
	afx_msg void OnBnClickedButtonSavetobmp();
	CSliderCtrl m_sliderctrlMinRoadWidth;
	virtual BOOL OnInitDialog();
	afx_msg void OnNMCustomdrawSlider1(NMHDR *pNMHDR, LRESULT *pResult);
	BOOL m_chkShowTheMapOn;
	afx_msg void OnBnClickedCheckShowthemap();
	BOOL m_chkShowExtractedBoundsOn;
	afx_msg void OnBnClickedCheckShowextractedboundaries();
	BOOL m_chkShowStreetGraphEndPointsOn;
	afx_msg void OnBnClickedCheckShowendpoints();
	BOOL m_chkShowPopDensityMapOn;
	afx_msg void OnBnClickedCheckShowpopdensitymap();
	BOOL m_chkShowIBFVOn;
	afx_msg void OnBnClickedCheckShowibfv();
	BOOL m_chkShowIntersectsOn;
	afx_msg void OnBnClickedCheckShowintersectionson();
	BOOL m_chkShowStreetUseNetworkOn;
	afx_msg void OnBnClickedCheckShowstreetusenetwork();
	BOOL m_chkShowLineStyleStreetsOn;
	afx_msg void OnBnClickedCheckShowlinestylestreeson();
	BOOL m_chkShowScalarFieldOn;
	afx_msg void OnBnClickedCheckShowscalarfield();
	afx_msg void OnBnClickedButtonTestscalarfield();
	BOOL m_chkShowMajRoadsOn;
	afx_msg void OnBnClickedCheckShowmajroads();
	double m_edRoadWidthValue;
	afx_msg void OnEnChangeEditRoadwidth();
	BOOL m_chkShowVegMapOn;
	afx_msg void OnBnClickedCheckShowvegmap();
	BOOL m_chkShowMajRoadGoogleStyleOn;
	afx_msg void OnBnClickedCheckUsegooglestyle();
	BOOL m_chkAntiAliasingOn;
	afx_msg void OnBnClickedCheckAntialiasingon();
};
