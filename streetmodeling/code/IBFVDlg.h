// IBFVDlg.h : header file
//

#pragma once

#include "GlView.h"
#include "afxcmn.h"
#include "afxwin.h"
#include "MyTabCtrl.h"
#include "Traceball.h"
#include "Polygon3D.h"

#include "SaveProjSettingDlg.h"


// CIBFVDlg dialog
class CIBFVDlg : public CDialog
{
// Construction
public:
	CIBFVDlg(CWnd* pParent = NULL);	// standard constructor
    CGlView *pclGlView;

// Dialog Data
	enum { IDD = IDD_IBFV2D_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support

	CPolygon3D mypolygon;
    BOOL m_LeftButtonDown, m_RightButtonDown, m_MiddleButtonDown;
	CPoint m_LeftDownPos;
	CPoint m_MiddleDownPos;
	double s_old, t_old;
	double s2_old, t2_old;   /*for the movement in second window*/
	double RotateDegree_old;
	double last_x, last_y;
	double pan_s_old, pan_t_old; 
	CTraceBall traceball;
	Quaternion rvec;

	short mouse_zDelta;

	double MixedEdgeRatio_old;
	double R3Filter_old;
	int TauStep_old;

	/*--------------------------------------3/27/06--*/
	int created;
	/*-----------------------------------------------*/

	//int dialog_leftx, dialog_rightx, dialog_buttomy, dialog_topy;
	CPoint firstwin_leftbottom;
	int firstwin_leftx, firstwin_bottomy, firstwin_rightx, firstwin_topy;

	CPoint secondwin_leftbottom;
	int secondwin_leftx, secondwin_bottomy, secondwin_rightx, secondwin_topy;

// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnTimer(UINT nIDEvent);
	MyTabCtrl m_tbCtrl;
	CStatic m_ctrlOpenGlWin;
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnMButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnMButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnBnClickedButtonClearall();
	afx_msg void OnFileNew32771();
	afx_msg void OnOpenLoadwatermap();
	afx_msg void OnFileExit();
	afx_msg void OnBnClickedButtonSaveproject();
	afx_msg void OnBnClickedButtonLoadproject();

	SaveProjSettingDlg *saveprojsetDlg;
	int saveprojsetDlgID;
	afx_msg void OnTcnSelchangeMytab(NMHDR *pNMHDR, LRESULT *pResult);
	CButton trace;
	afx_msg void OnBnClickedCancel();
	afx_msg void OnBnClickedButtonPlacetensorlines();
	afx_msg void OnBnClickedCheckTensordesignon();
	afx_msg void OnBnClickedRadioAddaregular();
	afx_msg void OnBnClickedRadioAddacenter();
};
