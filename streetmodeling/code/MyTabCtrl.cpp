// MyTabCtrl.cpp : implementation file
//

#include "stdafx.h"
#include "IBFV2D.h"
#include "MyTabCtrl.h"
#include "MyDlg1.h"
#include "MyDlg2.h"
#include "StreetNetPanel.h"
#include "StreetEditPanel.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

MyDlg1 *g_mydlg1;
MyDlg2 *g_mydlg2;
StreetNetPanel *g_streetnetpanel;
StreetEditPanel *g_streeteditpanel;

/////////////////////////////////////////////////////////////////////////////
// MyTabCtrl

MyTabCtrl::MyTabCtrl()
{
	m_DialogID[0] =IDD_DIALOG1;
	m_DialogID[1] =IDD_DIALOG_STREETNETPANEL;
	m_DialogID[3] =IDD_DIALOG_STREETEDITPANEL;
	m_DialogID[2] =IDD_DIALOG2;

	m_Dialog[0] = g_mydlg1 = new MyDlg1();
	m_Dialog[1] = g_streetnetpanel = new StreetNetPanel();
	m_Dialog[3] = g_streeteditpanel = new StreetEditPanel();
	m_Dialog[2] = g_mydlg2 = new MyDlg2();

	m_nPageCount = 4;

}

MyTabCtrl::~MyTabCtrl()
{
}

void MyTabCtrl::InitDialogs()
{
	m_Dialog[0]->Create(m_DialogID[0],GetParent());
	m_Dialog[1]->Create(m_DialogID[1],GetParent());
	m_Dialog[2]->Create(m_DialogID[2],GetParent());
	m_Dialog[3]->Create(m_DialogID[3],GetParent());
}


void MyTabCtrl::Reset()
{
	delete m_Dialog[0];
	delete m_Dialog[1];
	delete m_Dialog[2];
	delete m_Dialog[3];
	
	//m_Dialog[0] = new MyDlg1();
	//m_Dialog[1] = new StreetNetPanel();
	//m_Dialog[2] = new StreetEditPanel();
	//m_Dialog[2] = new MyDlg2();
	m_Dialog[0] = g_mydlg1 = new MyDlg1();
	m_Dialog[1] = g_streetnetpanel = new StreetNetPanel();
	m_Dialog[3] = g_streeteditpanel = new StreetEditPanel();
	m_Dialog[2] = g_mydlg2 = new MyDlg2();
}


BEGIN_MESSAGE_MAP(MyTabCtrl, CTabCtrl)
	//{{AFX_MSG_MAP(MyTabCtrl)
	ON_NOTIFY_REFLECT(TCN_SELCHANGE, OnSelchange)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// MyTabCtrl message handlers

void MyTabCtrl::OnSelchange(NMHDR* pNMHDR, LRESULT* pResult) 
{
	// TODO: Add your control notification handler code here
	ActivateTabDialogs();
	*pResult = 0;
}

void MyTabCtrl::ActivateTabDialogs()
{

	int nSel = GetCurSel();
	if(m_Dialog[nSel]->m_hWnd)
		m_Dialog[nSel]->ShowWindow(SW_HIDE);

	CRect l_rectClient;
	CRect l_rectWnd;

	GetClientRect(l_rectClient);
	AdjustRect(FALSE,l_rectClient);
	GetWindowRect(l_rectWnd);
	GetParent()->ScreenToClient(l_rectWnd);
	l_rectClient.OffsetRect(l_rectWnd.left,l_rectWnd.top);
	for(int nCount=0; nCount < m_nPageCount; nCount++){
		m_Dialog[nCount]->SetWindowPos(&wndTop, l_rectClient.left,l_rectClient.top,l_rectClient.Width(),l_rectClient.Height(),SWP_HIDEWINDOW);
	}
	m_Dialog[nSel]->SetWindowPos(&wndTop,l_rectClient.left,l_rectClient.top,l_rectClient.Width(),l_rectClient.Height(),SWP_SHOWWINDOW);

	m_Dialog[nSel]->ShowWindow(SW_SHOW);
}
