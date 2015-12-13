// MyDlg2.cpp : implementation file
//

#include "stdafx.h"
#include "IBFV2D.h"
#include "MyDlg2.h"
#include ".\mydlg2.h"
#include "shareinterfacevars.h"

#include "GlView.h"
#include ".\glview.h"

#include "BmpProcess.h"

#include "scalardesign.h"

#include "MyDlg1.h"
#include "StreetNetPanel.h"

extern MyDlg1 *g_mydlg1;
extern MyDlg2 *g_mydlg2;
extern StreetNetPanel *g_streetnetpanel;


extern SharedInterfaceVars sharedvars;
extern CGlView *g_pclGlView;

extern double VisMinRoadWidth;  //using opengl line width
extern double VisMajRoadWidth;  //using opengl line width
extern double VisHighWayWidth;  //using opengl line width

// MyDlg2 dialog

IMPLEMENT_DYNAMIC(MyDlg2, CDialog)
MyDlg2::MyDlg2(CWnd* pParent /*=NULL*/)
	: CDialog(MyDlg2::IDD, pParent)
	, m_chkShowSingularitiesOn(TRUE)
	, m_chkShowRegElemOn(TRUE)
	, m_chkShowTensorLinesOn(FALSE)
	, m_chkShowStreetGraphOn(FALSE)
	, m_chkShowRoadMapOn(FALSE)
	, m_chkShowTheMapOn(FALSE)
	, m_chkShowExtractedBoundsOn(TRUE)
	, m_chkShowStreetGraphEndPointsOn(TRUE)
	, m_chkShowPopDensityMapOn(FALSE)
	, m_chkShowIBFVOn(FALSE)
	, m_chkShowIntersectsOn(TRUE)
	, m_chkShowStreetUseNetworkOn(FALSE)
	, m_chkShowLineStyleStreetsOn(FALSE)
	, m_chkShowScalarFieldOn(FALSE)
	, m_chkShowMajRoadsOn(FALSE)
	, m_edRoadWidthValue(0)
	, m_chkShowVegMapOn(FALSE)
	, m_chkShowMajRoadGoogleStyleOn(FALSE)
	, m_chkAntiAliasingOn(FALSE)
{
}

MyDlg2::~MyDlg2()
{
}

void MyDlg2::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Check(pDX, IDC_CHECK_SHOWSINGULARITY, m_chkShowSingularitiesOn);
	DDX_Check(pDX, IDC_CHECK_SHOWREGULARELEM, m_chkShowRegElemOn);
	DDX_Check(pDX, IDC_CHECK_SHOWTENSORLINES, m_chkShowTensorLinesOn);
	DDX_Check(pDX, IDC_CHECK_SHOWSTREETGRAPH, m_chkShowStreetGraphOn);
	DDX_Check(pDX, IDC_CHECK_SHOWRAODMAP, m_chkShowRoadMapOn);
	DDX_Control(pDX, IDC_SLIDER1, m_sliderctrlMinRoadWidth);
	DDX_Check(pDX, IDC_CHECK_SHOWTHEMAP, m_chkShowTheMapOn);
	DDX_Check(pDX, IDC_CHECK_SHOWEXTRACTEDBOUNDARIES, m_chkShowExtractedBoundsOn);
	DDX_Check(pDX, IDC_CHECK_SHOWENDPOINTS, m_chkShowStreetGraphEndPointsOn);
	DDX_Check(pDX, IDC_CHECK_SHOWPOPDENSITYMAP, m_chkShowPopDensityMapOn);
	DDX_Check(pDX, IDC_CHECK_SHOWIBFV, m_chkShowIBFVOn);
	DDX_Check(pDX, IDC_CHECK_SHOWINTERSECTIONSON, m_chkShowIntersectsOn);
	DDX_Check(pDX, IDC_CHECK_SHOWSTREETUSENETWORK, m_chkShowStreetUseNetworkOn);
	DDX_Check(pDX, IDC_CHECK_SHOWLINESTYLESTREESON, m_chkShowLineStyleStreetsOn);
	DDX_Check(pDX, IDC_CHECK_SHOWSCALARFIELD, m_chkShowScalarFieldOn);
	DDX_Check(pDX, IDC_CHECK_SHOWMAJROADS, m_chkShowMajRoadsOn);
	DDX_Text(pDX, IDC_EDIT_ROADWIDTH, m_edRoadWidthValue);
	DDX_Check(pDX, IDC_CHECK_SHOWVEGMAP, m_chkShowVegMapOn);
	DDX_Check(pDX, IDC_CHECK_USEGOOGLESTYLE, m_chkShowMajRoadGoogleStyleOn);
	DDX_Check(pDX, IDC_CHECK_ANTIALIASINGON, m_chkAntiAliasingOn);
}


BEGIN_MESSAGE_MAP(MyDlg2, CDialog)
	ON_BN_CLICKED(IDC_CHECK_SHOWSINGULARITY, OnBnClickedCheckShowsingularity)
	ON_BN_CLICKED(IDC_CHECK_SHOWREGULARELEM, OnBnClickedCheckShowregularelem)
	ON_BN_CLICKED(IDC_CHECK_SHOWTENSORLINES, OnBnClickedCheckShowtensorlines)
	ON_BN_CLICKED(IDC_CHECK_SHOWSTREETGRAPH, OnBnClickedCheckShowstreetgraph)
	ON_BN_CLICKED(IDC_CHECK_SHOWRAODMAP, OnBnClickedCheckShowraodmap)
	ON_BN_CLICKED(IDC_BUTTON_SAVETOBMP, OnBnClickedButtonSavetobmp)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER1, OnNMCustomdrawSlider1)
	ON_BN_CLICKED(IDC_CHECK_SHOWTHEMAP, OnBnClickedCheckShowthemap)
	ON_BN_CLICKED(IDC_CHECK_SHOWEXTRACTEDBOUNDARIES, OnBnClickedCheckShowextractedboundaries)
	ON_BN_CLICKED(IDC_CHECK_SHOWENDPOINTS, OnBnClickedCheckShowendpoints)
	ON_BN_CLICKED(IDC_CHECK_SHOWPOPDENSITYMAP, OnBnClickedCheckShowpopdensitymap)
	ON_BN_CLICKED(IDC_CHECK_SHOWIBFV, OnBnClickedCheckShowibfv)
	ON_BN_CLICKED(IDC_CHECK_SHOWINTERSECTIONSON, OnBnClickedCheckShowintersectionson)
	ON_BN_CLICKED(IDC_CHECK_SHOWSTREETUSENETWORK, OnBnClickedCheckShowstreetusenetwork)
	ON_BN_CLICKED(IDC_CHECK_SHOWLINESTYLESTREESON, OnBnClickedCheckShowlinestylestreeson)
	ON_BN_CLICKED(IDC_CHECK_SHOWSCALARFIELD, OnBnClickedCheckShowscalarfield)
	ON_BN_CLICKED(IDC_BUTTON_TESTSCALARFIELD, OnBnClickedButtonTestscalarfield)
	ON_BN_CLICKED(IDC_CHECK_SHOWMAJROADS, OnBnClickedCheckShowmajroads)
	ON_EN_CHANGE(IDC_EDIT_ROADWIDTH, OnEnChangeEditRoadwidth)
	ON_BN_CLICKED(IDC_CHECK_SHOWVEGMAP, OnBnClickedCheckShowvegmap)
	ON_BN_CLICKED(IDC_CHECK_USEGOOGLESTYLE, OnBnClickedCheckUsegooglestyle)
	ON_BN_CLICKED(IDC_CHECK_ANTIALIASINGON, OnBnClickedCheckAntialiasingon)
END_MESSAGE_MAP()


void MyDlg2::Reset()
{
	m_chkShowSingularitiesOn=TRUE;
	m_chkShowRegElemOn=TRUE;
	m_chkShowTensorLinesOn=FALSE;
	m_chkShowStreetGraphOn=FALSE;
	m_chkShowRoadMapOn=FALSE;
}


// MyDlg2 message handlers

void MyDlg2::OnBnClickedCheckShowsingularity()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowSingularitiesOn == FALSE)
	{
		m_chkShowSingularitiesOn=TRUE;
		sharedvars.ShowSingularitiesOn=true;
		g_pclGlView->SingularitiesOn = 1;
	}
	else
	{
		m_chkShowSingularitiesOn=FALSE;
		sharedvars.ShowSingularitiesOn=false;
		g_pclGlView->SingularitiesOn = 0;
	}
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedCheckShowregularelem()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowRegElemOn==FALSE)
	{
		m_chkShowRegElemOn=TRUE;
		g_pclGlView->RegularElemOn=1;
		sharedvars.ShowRegElemOn=true;
	}
	else
	{
		m_chkShowRegElemOn=FALSE;
		g_pclGlView->RegularElemOn=0;
		sharedvars.ShowRegElemOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedCheckShowtensorlines()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowTensorLinesOn==FALSE)
	{
		m_chkShowTensorLinesOn=TRUE;
		sharedvars.ShowTensorLinesOn=true;
		g_pclGlView->showTensorLineOn=true;
	}
	else
	{
		m_chkShowTensorLinesOn=FALSE;
		sharedvars.ShowTensorLinesOn=false;
		g_pclGlView->showTensorLineOn=false;
	}
}

void MyDlg2::OnBnClickedCheckShowstreetgraph()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowStreetGraphOn==FALSE)
	{
		m_chkShowStreetGraphOn=TRUE;
		sharedvars.ShowStreetGraphOn=true;
		g_pclGlView->showStreetGraphOn=true;
	}
	else
	{
		m_chkShowStreetGraphOn=FALSE;
		sharedvars.ShowStreetGraphOn=false;
		g_pclGlView->showStreetGraphOn=false;
	}
}

void MyDlg2::OnBnClickedCheckShowraodmap()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowRoadMapOn==FALSE)
	{
		m_chkShowRoadMapOn=TRUE;
		m_chkShowStreetUseNetworkOn=FALSE;
		m_chkShowIBFVOn=FALSE;

		g_mydlg1->m_chkBrushInterfaceOn=FALSE;
		sharedvars.BrushInterfaceOn=false;

		sharedvars.ShowRoadMapOn=true;
		g_pclGlView->displayRoadNetOn=true;
		g_pclGlView->ShapeControlPtsOn=0;
		sharedvars.ShowIBFVOn=false;
		sharedvars.ShowPopDensityMapOn=false;
		sharedvars.ShowStreetUseNetworkOn=false;
	}
	else
	{
		m_chkShowRoadMapOn=FALSE;
		sharedvars.ShowRoadMapOn=false;
		g_pclGlView->displayRoadNetOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedButtonSavetobmp()
{
	// TODO: Add your control notification handler code here
	CFile f;

	char strFilter[] = { "Bitmap Files (*.bmp)|*.bmp|All Files (*.*)|*.*||" };

	CFileDialog FileDlg(FALSE, ".bmp", NULL, 0, strFilter);
	if( FileDlg.DoModal() == IDOK )
	{
		CString f2 = FileDlg.GetFileName();
		CString f3 = f2.Left(f2.Find('.'))+".bmp";

		const char *filename = f3;

		WriteBmp(filename);
	
	}
}

BOOL MyDlg2::OnInitDialog()
{
	CDialog::OnInitDialog();

	// TODO:  Add extra initialization here
	m_sliderctrlMinRoadWidth.SetRange(2, 20);
	m_sliderctrlMinRoadWidth.SetPos(4);
	m_sliderctrlMinRoadWidth.SetTicFreq(1);
	m_sliderctrlMinRoadWidth.SetTic(0);

	m_edRoadWidthValue=(double)m_sliderctrlMinRoadWidth.GetPos();

		m_chkShowTensorLinesOn=TRUE;
		sharedvars.ShowTensorLinesOn=true;
		g_pclGlView->showTensorLineOn=true;
		
		m_chkShowMajRoadsOn=TRUE;
		sharedvars.ShowMajRoadsOn=true;

	UpdateData(FALSE);


	return TRUE;  // return TRUE unless you set the focus to a control
	// EXCEPTION: OCX Property Pages should return FALSE
}

void MyDlg2::OnNMCustomdrawSlider1(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	// TODO: Add your control notification handler code here
	VisMinRoadWidth=(double)m_sliderctrlMinRoadWidth.GetPos();
	m_edRoadWidthValue=VisMinRoadWidth;
	*pResult = 0;
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedCheckShowthemap()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowTheMapOn==FALSE)
	{
		m_chkShowTheMapOn=TRUE;
		m_chkShowIBFVOn=FALSE;
		m_chkShowPopDensityMapOn=FALSE;
		m_chkShowRoadMapOn=FALSE;

		sharedvars.ShowTheMapOn=true;
		sharedvars.ShowIBFVOn=false;
		sharedvars.ShowPopDensityMapOn=false;
		sharedvars.ShowRoadMapOn=false;
		
	}
	else
	{
		m_chkShowTheMapOn=FALSE;
		sharedvars.ShowTheMapOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedCheckShowextractedboundaries()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowExtractedBoundsOn==FALSE)
	{
		m_chkShowExtractedBoundsOn=TRUE;
		sharedvars.ShowExtractedBoundsOn=true;
	}
	else
	{
		m_chkShowExtractedBoundsOn=FALSE;
		sharedvars.ShowExtractedBoundsOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedCheckShowendpoints()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowStreetGraphEndPointsOn==FALSE)
	{
		m_chkShowStreetGraphEndPointsOn=TRUE;
		sharedvars.ShowStreetGraphEndPointsOn=true;
	}
	else
	{
		m_chkShowStreetGraphEndPointsOn=FALSE;
		sharedvars.ShowStreetGraphEndPointsOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedCheckShowpopdensitymap()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowPopDensityMapOn==FALSE)
	{
		m_chkShowPopDensityMapOn=TRUE;
		m_chkShowIBFVOn=FALSE;
		m_chkShowTheMapOn=FALSE;
		m_chkShowRoadMapOn=FALSE;

		sharedvars.ShowPopDensityMapOn=true;
		sharedvars.ShowRoadMapOn=false;
		g_pclGlView->displayRoadNetOn=false;
		sharedvars.ShowIBFVOn=false;
		sharedvars.ShowTheMapOn=false;
	}
	else
	{
		m_chkShowPopDensityMapOn=FALSE;
		sharedvars.ShowPopDensityMapOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedCheckShowibfv()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowIBFVOn==FALSE)
	{
		m_chkShowIBFVOn=TRUE;
		m_chkShowPopDensityMapOn=FALSE;
		m_chkShowTheMapOn=FALSE;
		m_chkShowRoadMapOn=FALSE;
		m_chkShowStreetUseNetworkOn=FALSE;
				
		m_chkShowScalarFieldOn=FALSE;
		sharedvars.ShowScalarFieldOn=false;


		sharedvars.ShowIBFVOn=true;
		sharedvars.ShowRoadMapOn=false;
		g_pclGlView->displayRoadNetOn=false;
		sharedvars.ShowPopDensityMapOn=false;
		sharedvars.ShowTheMapOn=false;
		sharedvars.ShowStreetUseNetworkOn=false;
	}
	else
	{
		m_chkShowIBFVOn=FALSE;
		sharedvars.ShowIBFVOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedCheckShowintersectionson()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowIntersectsOn==FALSE)
	{
		m_chkShowIntersectsOn=TRUE;
		sharedvars.ShowIntersectsOn=true;
	}
	else
	{
		m_chkShowIntersectsOn=FALSE;
		sharedvars.ShowIntersectsOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedCheckShowstreetusenetwork()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowStreetUseNetworkOn==FALSE)
	{
		m_chkShowStreetUseNetworkOn=TRUE;
		m_chkShowIBFVOn=FALSE;
		m_chkShowRoadMapOn=FALSE;
		g_pclGlView->displayRoadNetOn=true;

		sharedvars.ShowStreetUseNetworkOn=true;
		sharedvars.ShowRoadMapOn=false;
		sharedvars.ShowIBFVOn=false;
		sharedvars.ShowPopDensityMapOn=false;
	}
	else
	{
		m_chkShowStreetUseNetworkOn=FALSE;
		sharedvars.ShowStreetUseNetworkOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedCheckShowlinestylestreeson()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowLineStyleStreetsOn==FALSE)
	{
		m_chkShowLineStyleStreetsOn=TRUE;
		sharedvars.ShowLineStyleStreetsOn=true;
	}
	else
	{
		m_chkShowLineStyleStreetsOn=FALSE;
		sharedvars.ShowLineStyleStreetsOn=false;
	}
}

void MyDlg2::OnBnClickedCheckShowscalarfield()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowScalarFieldOn==FALSE)
	{
		m_chkShowScalarFieldOn=TRUE;
		sharedvars.ShowScalarFieldOn=true;

		m_chkShowIBFVOn=FALSE;
		sharedvars.ShowIBFVOn=false;
	}
	else
	{
		m_chkShowScalarFieldOn=FALSE;
		sharedvars.ShowScalarFieldOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedButtonTestscalarfield()
{
	// TODO: Add your control notification handler code here
	test_scalar_field();
}

void MyDlg2::OnBnClickedCheckShowmajroads()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowMajRoadsOn==FALSE)
	{
		m_chkShowMajRoadsOn=TRUE;
		sharedvars.ShowMajRoadsOn=true;
	}
	else
	{
		m_chkShowMajRoadsOn=FALSE;
		sharedvars.ShowMajRoadsOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg2::OnEnChangeEditRoadwidth()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	UpdateData();
	//VisMinRoadWidth=(double)m_sliderctrlMinRoadWidth.GetPos();
	VisMinRoadWidth=m_edRoadWidthValue;
	m_sliderctrlMinRoadWidth.SetPos(floor(VisMinRoadWidth));
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedCheckShowvegmap()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowVegMapOn==FALSE)
	{
		m_chkShowVegMapOn=TRUE;
	}
	else
	{
		m_chkShowVegMapOn=FALSE;
	}
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedCheckUsegooglestyle()
{
	// TODO: Add your control notification handler code here
	if(m_chkShowMajRoadGoogleStyleOn==FALSE)
	{
		m_chkShowMajRoadGoogleStyleOn=TRUE;
		sharedvars.ShowMajRoadGoogleStyleOn=true;
	}
	else
	{
		m_chkShowMajRoadGoogleStyleOn=FALSE;
		sharedvars.ShowMajRoadGoogleStyleOn=false;
	}
	UpdateData(FALSE);
}

void MyDlg2::OnBnClickedCheckAntialiasingon()
{
	// TODO: Add your control notification handler code here
	if(m_chkAntiAliasingOn==FALSE)
	{
		m_chkAntiAliasingOn=TRUE;
		sharedvars.AntiAliasingOn=true;
	}
	else
	{
		m_chkAntiAliasingOn=FALSE;
		sharedvars.AntiAliasingOn=false;
	}
	UpdateData(FALSE);
}
