#pragma once
#include "afxwin.h"
#include <GL/glut.h>
#include "lib/icVector.h"

class CConleyGraphView :
	public CWnd
{
public:
	CConleyGraphView(void);
	CConleyGraphView(CWnd *pclWnd);

	virtual ~CConleyGraphView(void);

	////attributes
public:
	HDC  m_hDC;		        // GDI Device Context 
    HGLRC	m_hglRC;		// Rendering Context

    CWnd *m_pclWnd;
    HWND m_hWnd;

	////For singularities pair cancellation
	int PairCancelOn;
	int SelectRepellerOrAttractor;

	int cancel_node1, cancel_node2;
	int *cancel_nodes;
	int num_cancel_nodes;

	int PairCounter;

	int ShowMCGOn;

	////methods
public:
	BOOL SetPixelformat(HDC hdc);
	//GLvoid ReSizeGLScene(GLsizei width, GLsizei height);
	int InitGL(GLvoid);	
	int DrawGLScene(GLenum mode);
	int OnCreate();

	////Draw the components of the graph
	void DrawNodes(GLenum mode);
	void DrawEdges();
	void DrawWings(double head[2], icVector2 direct);

    void DrawSolidCircle(double cx, double cy, double R);
	void DrawSolidSquare(double cx, double cy);
	void DrawCircle(double cx, double cy);
    void DrawLimitCircle(double cx, double cy, int type);
	void DrawHighLights();
	void SetColor(int nodetype);

	void DrawLimitCycleLegend(double cx, double cy);

	void Displaylabel(int x, int y, char *string); ////display the labels, testing at 3/14/06
	void DrawLabels();
	
	
	void SaveTheImagesTemporary();

	void HitProcessforGraph(double ss, double tt);
	void HitProcessforCancel(double ss, double tt);

	void GetTheHighLightNodeandEdges();
	void GetHighLightSeparatrices();
	int FindConnectedSep(int saddle, int othersing, int type);

    void InitFlags();


	////For multiple repellers and attractors selection  11/20/05
	////Note that these routines are different from the routines in class "CGlView"
	void AddToRepellerList(int);
	void AddToAttractorList(int);
	void ClearRepellandAttractList();

	int show3DOn;

	DECLARE_MESSAGE_MAP()
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
};

void display(void);

void display_roads_3d(GLenum mode);
