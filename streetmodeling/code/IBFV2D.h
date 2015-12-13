// IBFV2D.h : main header file for the PROJECT_NAME application
//

#pragma once

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// main symbols


// CIBFVApp:
// See IBFV2D.cpp for the implementation of this class
//

class CIBFVApp : public CWinApp
{
public:
	CIBFVApp();

// Overrides
	public:
	virtual BOOL InitInstance();

// Implementation

	DECLARE_MESSAGE_MAP()
};

extern CIBFVApp theApp;