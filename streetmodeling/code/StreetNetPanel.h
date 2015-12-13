#pragma once
#include "afxcmn.h"
#include "afxwin.h"

#include "MajRoadSettingDlg.h"

#include "SeedPtCtrlDlg.h"

#include "NoiseCtrlPanel.h"

#include "MinRoadSetDlg.h"

// StreetNetPanel dialog

class StreetNetPanel : public CDialog
{
	DECLARE_DYNAMIC(StreetNetPanel)

public:
	StreetNetPanel(CWnd* pParent = NULL);   // standard constructor
	virtual ~StreetNetPanel();
	
	void Reset();

// Dialog Data
	enum { IDD = IDD_DIALOG_STREETNETPANEL };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedButtonPlacetensorlines();
	double m_edMinTenLineDensity;
	double m_edMajTenLineDensity;
	afx_msg void OnBnClickedButtonComputestreetnetwork();
	BOOL m_chkEditStreetNetOn;
	afx_msg void OnBnClickedCheckEditstreetgraphon();
	BOOL m_chkMeanFilterOn;
	afx_msg void OnBnClickedCheckMeanfilteron();
	afx_msg void OnBnClickedButtonLoadawatermap();
	BOOL m_chkGenTenFromExtractedBoundsOn;
	afx_msg void OnBnClickedCheckGentensorbasedonextractboundaries();
	BOOL m_chkUseAllBoundsOn;
	afx_msg void OnBnClickedCheckUseallboundaries();
	CSliderCtrl m_ctrlSliderBoundaryRegionWidth;
	virtual BOOL OnInitDialog();
	afx_msg void OnNMCustomdrawSliderBoundarywidth(NMHDR *pNMHDR, LRESULT *pResult);
	double m_edMinTensorLineLength;
	BOOL m_chkDesignGridOn;
	afx_msg void OnBnClickedCheckGridon();
	int m_edStreetLevel;
	BOOL m_chkConnectDeadEndsOn;
	afx_msg void OnBnClickedCheckConnectdeadendon();
	CSpinButtonCtrl m_ctrlSpinTracingLevel;
	afx_msg void OnDeltaposSpinTracinglevel(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnBnClickedButtonConnectdeadend();
	afx_msg void OnBnClickedButtonComputeregionblocks();
	BOOL m_chkShowRegionBlocksOn;
	afx_msg void OnBnClickedCheckShowregionblocks();
	afx_msg void OnBnClickedButtonTestdiffregiondesign();
	afx_msg void OnLoadLoaddensitymap();
	afx_msg void OnBnClickedButtonLoadpopdensitymap();
	BOOL m_chkCombinePopDensityOn;
	afx_msg void OnBnClickedCheckCombinedensitymap();
	BOOL m_chkGenTenFromSketchesOn;
	afx_msg void OnBnClickedCheckGentenSketchbased();
	BOOL m_chkCloseLoopOn;
	afx_msg void OnBnClickedCheckCloseloop();
	afx_msg void OnFileSaveas();
	afx_msg void OnBnClickedButtonExportstreetnetwork();
	afx_msg void OnEnChangeEditStreetlevel();
	afx_msg void OnEnChangeEditMajordensity();
	afx_msg void OnEnChangeEditMinordensity();
	afx_msg void OnEnChangeEditMintensorlinelength();
	afx_msg void OnBnClickedButtonAddnoisetoroads();
	CEdit m_EdCtrlSetStreetLevels;
	BOOL m_chkUseBoundsAsSketchesOn;
	afx_msg void OnBnClickedCheckUseboundariesassketches();
	double m_edEffectWidthBoundary;
	afx_msg void OnEnChangeEditWidthboundary();
	int m_edMeanFilterSize;
	afx_msg void OnEnChangeEdit2();
	BOOL m_chkShowInitSeedsOn;
	afx_msg void OnBnClickedCheckShowinitseeds();
	CButton m_ctrlMajRoadSettingButton;
	afx_msg void OnBnClickedButtonMajroadsetting();
	CButton m_ctrlSeedPtSetting;

    MajRoadSettingDlg *maj_road_settingDlg;
	int maj_road_settingDlgID;
	afx_msg void OnBnClickedButtonSeedsetting();

	SeedPtCtrlDlg *seeds_ctrlDlg;
	int seeds_ctrlDlgID;

	NoiseCtrlPanel *noise_ctrlDlg;
	int noise_ctrlDlgID;

	MinRoadSetDlg *min_road_setDlg;
	int min_road_setDlgID;

	afx_msg void OnBnClickedButtonLoadvegmap();
	BOOL m_chkApplyAsymFldOn;
	afx_msg void OnBnClickedCheckApplyasymfld();
	afx_msg void OnBnClickedButtonAddnoisetofield();
	afx_msg void OnBnClickedButtonNoisesetting();
	afx_msg void OnBnClickedButtonMinroadsetting();
	BOOL m_chkRemoveMajDeadEndsOn;
	afx_msg void OnBnClickedCheckRemovemajdeadednon();
	afx_msg void OnBnClickedButtonImportstreetnetwork();
	double m_edMinBoundaryLen;
	afx_msg void OnEnChangeEditMinboundarylen();
	BOOL m_chkSelRegToEditOn;
	afx_msg void OnBnClickedCheckSelregtoediton();
	afx_msg void OnBnClickedButtonRemovestreetsinreg();
	afx_msg void OnBnClickedButtonReplacestreetsinside();
	afx_msg void OnBnClickedButtonCompstreetnetinside();
	afx_msg void OnBnClickedButtonAddnoisetolocalfield();
	afx_msg void OnBnClickedButtonAddnoiselocalroads();
	BOOL m_chkBrushLikeRegSelOn;
	afx_msg void OnBnClickedCheckUsebrushlikeregsel();
	BOOL m_rdSubConnectMethod;
	afx_msg void OnBnClickedRadioConmethod1();
	afx_msg void OnBnClickedRadioConmethod2();
};
