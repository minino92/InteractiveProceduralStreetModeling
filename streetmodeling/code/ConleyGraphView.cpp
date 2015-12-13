#include "StdAfx.h"
#include ".\conleygraphview.h"

#include "VFDataStructure.h"




///////
extern GraphNode *graphnodes;
extern GraphEdge *graphedges;
extern int cur_node_index;
extern int cur_graphedge_index;

extern MCGNode *mcgnodes;
extern MCGEdge *mcgedges;
extern int cur_mcgnode_index;
extern int cur_mcgedge_index;


/*variables for MCG*/

/////////////////////////////////////////////////// 10/16/05
int picked_node = -1;                  ////The picked up node when mouse moves to it
int *related_edges = NULL;        ////The corresponding edges incident to the node above
int num_related_edges = 0;

////the list for the highlighted singularities and separatrices 1/28/06
int *related_trajs = NULL;
int num_related_trajs = 0;

////
extern int *repeller_nodes;     //the indices of the selected repellers in the conley graph
extern int NumRepellers;        //the number of the being selected repellers
extern int *attractor_nodes;    //the indices of the selected attractors in the conley graph
extern int NumAttractors;       //the number of the being selected attractors


/*----------------------------------------------------------*/
///1/28/06
extern Singularities *singularities;           //being captured singularites' list
//extern int cur_singularity_index;
//extern Trajectory *trajectories2;
extern LineSeg **trajectories;                 //trajectories' list
//extern int cur_traj_index;
extern int *num_linesegs_curtraj;              //an array stores the number of line segments for corresponding trajectory
extern Separatrices *separatrices;             //array for group of separatrices

extern LimitCycle *limitcycles;

extern Polygon3D Object;

float rotmat[4][4];


/*----------------------------------------------------------*/

extern int *Extend_link(int *, int);

/////////////////////////////////////////////////// 11/06/05
////variables for pair cancellation through selection from the graph
//int cancel_node1, cancel_node2;
//int *cancel_nodes;
//int num_cancel_nodes;


CConleyGraphView::CConleyGraphView(CWnd *pclWnd)
{

    m_pclWnd = pclWnd;
    m_hWnd   = pclWnd->m_hWnd;
    m_hDC    = ::GetDC(m_hWnd);

	cancel_nodes = NULL;

	//mat_ident(rotmat);
	int i;

	for (i = 0; i <= 3; i++) {
		rotmat[i][0] = 0.0;
		rotmat[i][1] = 0.0;
		rotmat[i][2] = 0.0;
		rotmat[i][3] = 0.0;
		rotmat[i][i] = 1.0;
	}

	InitFlags();

}

CConleyGraphView::CConleyGraphView(void)
{
}

CConleyGraphView::~CConleyGraphView(void)
{
}


int CConleyGraphView::OnCreate() 
{

    m_hDC = ::GetDC(this->m_hWnd);

    if(!SetPixelformat(m_hDC))
    {
	::MessageBox(::GetFocus(),"SetPixelformat Failed!","Error",MB_OK);
	return -1;
    }

    m_hglRC = wglCreateContext(m_hDC);
    int i= wglMakeCurrent(m_hDC,m_hglRC);

	InitGL();	

	return 0;
}

/*--------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------*/
// Other extended operations for OpenGl window setting

BOOL CConleyGraphView::SetPixelformat(HDC hdc)
{

    PIXELFORMATDESCRIPTOR *ppfd; 
    int pixelformat; 
 
    PIXELFORMATDESCRIPTOR pfd = { 
    sizeof(PIXELFORMATDESCRIPTOR),  //  size of this pfd 
    1,                     // version number 
    PFD_DRAW_TO_WINDOW |   // support window 
    PFD_SUPPORT_OPENGL |   // support OpenGL 
    PFD_GENERIC_FORMAT |
    PFD_DOUBLEBUFFER,      // double buffered 
    PFD_TYPE_RGBA,         // RGBA type 
    32,                    // 24-bit color depth 
    0, 0, 0, 0, 0, 0,      // color bits ignored 
    8,                     // no alpha buffer 
    0,                     // shift bit ignored 
    8,                     // no accumulation buffer 
    0, 0, 0, 0,            // accum bits ignored 
    64,                    // 32-bit z-buffer	 
    8,                     // no stencil buffer 
    8,                     // no auxiliary buffer 
    PFD_MAIN_PLANE,        // main layer 
    0,                     // reserved 
    0, 0, 0                // layer masks ignored 
    }; 
 
   
    ppfd = &pfd;

 
    if ( (pixelformat = ChoosePixelFormat(hdc, ppfd)) == 0 ) 
    { 
        ::MessageBox(NULL, "ChoosePixelFormat failed", "Error", MB_OK); 
        return FALSE; 
    } 
 
    if (SetPixelFormat(hdc, pixelformat, ppfd) == FALSE) 
    { 
        ::MessageBox(NULL, "SetPixelFormat failed", "Error", MB_OK); 
        return FALSE; 
    } 
 
    return TRUE; 

   
    ppfd = &pfd;

 
    if ( (pixelformat = ChoosePixelFormat(hdc, ppfd)) == 0 ) 
    { 
        ::MessageBox(NULL, "ChoosePixelFormat failed", "Error", MB_OK); 
        return FALSE; 
    } 
 
    if (SetPixelFormat(hdc, pixelformat, ppfd) == FALSE) 
    { 
        ::MessageBox(NULL, "SetPixelFormat failed", "Error", MB_OK); 
        return FALSE; 
    } 
 
    return TRUE; 

}

/* --------------------------------------------------
*  Initialize the opengl environment
---------------------------------------------------*/
int CConleyGraphView::InitGL(GLvoid)								// All Setup For OpenGL Goes Here
{
	glMatrixMode(GL_PROJECTION);		// Select The Projection Matrix
	glLoadIdentity();					// Reset The Projection Matrix

	// Calculate The Aspect Ratio Of The Window
	glOrtho(0.0f,4.0f,0.0f,2.0f,0.0f,50.0f);

	glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
	glLoadIdentity();									// Reset The Modelview Matrix
	
	glClearColor(1.0f, 1.0f, 1.0f, 1.f);				// Black Background

	return TRUE;										// Initialization Went OK
}

/*--------------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------------*/
//Display routine
/*--------------------------------------------------------------------------------*/
int CConleyGraphView::DrawGLScene(GLenum mode)					// Here's Where We Do All The Drawing
{
	////Draw the conley relation graph according to the graphnodes and graphedges structure
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer

	if(show3DOn == 0)
	{
		if(num_related_edges > 0 && picked_node >= 0)
			DrawHighLights();     //Draw highlighted objects if exist

		DrawEdges();          //Draw edges first
	 
   		DrawNodes(mode);          //Draw nodes

		DrawLabels();

		//SaveTheImagesTemporary();
	}

	else
	{
		display();
		display_roads_3d(mode);
	}

	SwapBuffers(m_hDC);
	return TRUE;										// Keep Going
}


////Draw the components of the graph

/* --------------------------------------------------
*  Draw the nodes of the graph
---------------------------------------------------*/
void CConleyGraphView::DrawNodes(GLenum mode)
{
	int i;

	if(ShowMCGOn == 0)
	{
		for(i = 0; i < cur_node_index; i++)
		{
			if(graphnodes[i].cancelled == 1)
				continue;

			if(mode == GL_SELECT)
				glLoadName(i+1);
			
			SetColor(graphnodes[i].type);

			if(graphnodes[i].LimitCycleID >= 0)
			{
				//DrawSolidSquare(graphnodes[i].pos_x, graphnodes[i].pos_y);
				int type = limitcycles[graphnodes[i].LimitCycleID].type;
				DrawLimitCircle(graphnodes[i].pos_x, graphnodes[i].pos_y, type);
				DrawSolidCircle(graphnodes[i].pos_x, graphnodes[i].pos_y, 0.04);
			}

			else{

				DrawSolidCircle(graphnodes[i].pos_x, graphnodes[i].pos_y, 0.04);
			}
		}
	}

	else{
		for(i = 0; i < cur_mcgnode_index; i++)
		{
			if(mcgnodes[i].cancelled == 1)
				continue;

			if(mode == GL_SELECT)
				glLoadName(i+1);
			
			SetColor(mcgnodes[i].type);

			DrawSolidCircle(mcgnodes[i].pos_x, mcgnodes[i].pos_y, 0.04);
		}
	}
}

/* --------------------------------------------------
*  Draw the solid circle for the nodes in the graph
---------------------------------------------------*/
void CConleyGraphView::DrawSolidCircle(double cx, double cy, double R)
{
	int i;
	//double R = 0.06;
	double theta, deta ;
	deta = 2 * M_PI/80.;
	double x, y;
	theta = 0.;
	glBegin(GL_POLYGON);
	for(i = 0; i < 80; i++, theta += deta)
	{
		x =  cx + R * cos(theta);
		y =  cy + R * sin(theta);

		glVertex2f(x, y);
	}
	glEnd();
}


/* --------------------------------------------------
*  Draw a circle for the nodes in the graph
---------------------------------------------------*/
void CConleyGraphView::DrawLimitCircle(double cx, double cy, int type)
{
	int i;
	double R = 0.085;
	double theta, deta ;
	deta = 2 * M_PI/80.;
	double x, y;
	theta = 0.;
	glLineWidth(1.8);
	if(type == 0) //repeller
		glColor3f(0, 1, 0);
	else 
		glColor3f(1, 0, 0);

	glBegin(GL_LINE_LOOP);
	for(i = 0; i < 80; i++, theta += deta)
	{
		x =  cx + R * cos(theta);
		y =  cy + R * sin(theta);

		glVertex2f(x, y);
	}
	glEnd();
	glLineWidth(1.);
}

/* --------------------------------------------------
*  Draw the solid square for the nodes representing
*  limit cycles in the graph
---------------------------------------------------*/
void CConleyGraphView::DrawSolidSquare(double cx, double cy)
{
	glColor3f(0.8, 1, 0);
	glBegin(GL_POLYGON);
	glVertex2f(cx-0.04, cy-0.04);
	glVertex2f(cx-0.04, cy+0.04);
	glVertex2f(cx+0.04, cy+0.04);
	glVertex2f(cx+0.04, cy-0.04);
	glEnd();
}

extern int WriteBmp(const char*);
void CConleyGraphView::SaveTheImagesTemporary()
{
	int viewport[4] = {0};

	glGetIntegerv(GL_VIEWPORT, viewport);

	int dx = viewport[2] - viewport[0];
	int dy = viewport[3] - viewport[1];

	//array = new unsigned char[3*dx*dy];
	
	glPixelStorei( GL_PACK_ALIGNMENT, 4 );
	glPixelStorei( GL_PACK_ROW_LENGTH, dx );
	//glReadBuffer( GL_FRONT );
	glReadBuffer( GL_BACK );
	//glReadPixels( 0, 0, dx, dy, GL_RGB, GL_UNSIGNED_BYTE, array );

	//WriteBmp("C_graphTemp.bmp");

}


/* --------------------------------------------------
*  Draw the legend of the nodes representing 
*  limit cyclein the graph
---------------------------------------------------*/
//void DrawLimitCycleLegend(double cx, double cy)
//{
//	DrawSolidSquare(cx, cy);
//
//	Draw
//}

////Draw the wings of an arrow according to its head and direction 1/31/06
void CConleyGraphView::DrawWings(double head[2], icVector2 direct)
{
	glPushMatrix();
	glTranslatef(head[0], head[1], 0);
	glRotatef(atan2(direct.entry[1], direct.entry[0])*360/(2*M_PI), 0, 0, 1);
	glScalef(0.33, 0.33, 1);

	////Draw the wings of the arrow
	//glBegin(GL_LINES);
	//glVertex2f(0, 0);
	//glVertex2f(-0.2, 0.16);

	//glVertex2f(0, 0);
	//glVertex2f(-0.2, -0.16);
	//glEnd();
	glBegin(GL_TRIANGLES);
	glVertex2f(0, 0);
	glVertex2f(-0.35, 0.12);
    glVertex2f(-0.35, -0.12);
	glEnd();

	glPopMatrix();
}



/* --------------------------------------------------
*  Draw the edges of the graph
---------------------------------------------------*/
void CConleyGraphView::DrawEdges()
{
	int i;
	int node1, node2;
	icVector2 direct;
	double head[2];

	glLineWidth(1.5);

	if(ShowMCGOn == 0)
	{
		for(i = 0; i < cur_graphedge_index; i++)
		{
			node1 = graphedges[i].node_index1;
			node2 = graphedges[i].node_index2;

			if(	graphnodes[node1].cancelled == 1 || graphnodes[node2].cancelled == 1)
				continue;


			glBegin(GL_LINES);
			SetColor(graphnodes[node1].type);
			glVertex2f(graphnodes[node1].pos_x, graphnodes[node1].pos_y);
			
			SetColor(graphnodes[node2].type);
			glVertex2f(graphnodes[node2].pos_x , graphnodes[node2].pos_y + 0.02);
			glEnd();

			////Draw the wings of the arrow
			direct.entry[0] = graphnodes[node2].pos_x - graphnodes[node1].pos_x;
			direct.entry[1] = graphnodes[node2].pos_y - graphnodes[node1].pos_y;
			normalize(direct);

			head[0] = graphnodes[node2].pos_x ;
			head[1] = graphnodes[node2].pos_y + 0.02;
			DrawWings(head, direct);
		}
	}

	else
	{
		for(i = 0; i < cur_mcgedge_index; i++)
		{
			node1 = mcgedges[i].node_index1;
			node2 = mcgedges[i].node_index2;

			if(mcgedges[i].cancel == true)
				continue;

			if(	mcgnodes[node1].cancelled == 1 || mcgnodes[node2].cancelled == 1)
				continue;


			glBegin(GL_LINES);
			SetColor(mcgnodes[node1].type);
			glVertex2f(mcgnodes[node1].pos_x, mcgnodes[node1].pos_y);
			
			SetColor(mcgnodes[node2].type);
			glVertex2f(mcgnodes[node2].pos_x , mcgnodes[node2].pos_y + 0.02);
			glEnd();

			////Draw the wings of the arrow
			direct.entry[0] = mcgnodes[node2].pos_x - mcgnodes[node1].pos_x;
			direct.entry[1] = mcgnodes[node2].pos_y - mcgnodes[node1].pos_y;
			normalize(direct);

			head[0] = mcgnodes[node2].pos_x ;
			head[1] = mcgnodes[node2].pos_y + 0.02;
			DrawWings(head, direct);
		}
	}

}

	
void CConleyGraphView::Displaylabel(int x, int y, char *string)
{
  int len, i;
  glRasterPos2f(x, y);
  len = (int) strlen(string);
  for (i = 0; i < len; i++)
  {
    //glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, string[i]);
     glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, string[i]);
 }
}


////Display all the labels
void CConleyGraphView::DrawLabels()
{
	int i;
	int posx, posy;
	char strings[5];
	//int r_counter, a_counter, s_counter;

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	glMatrixMode(GL_PROJECTION);		// Select The Projection Matrix
	glLoadIdentity();					// Reset The Projection Matrix


	// Calculate The Aspect Ratio Of The Window
	glOrtho(viewport[0],viewport[2],viewport[1],viewport[3],0.0f,50.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//r_counter = a_counter = s_counter = 0;
	if(ShowMCGOn == 0)
	{
		for(i = 0; i < cur_node_index; i++)
		{
			if(graphnodes[i].cancelled == 1)
				continue;

			glPushMatrix();

			posx = (int)((graphnodes[i].pos_x/4.)*(float)viewport[2])-34;
			posy = (int)((graphnodes[i].pos_y/2.)*(float)viewport[3])-5;

			glColor3f(0,0,0);
			if(graphnodes[i].type == 0) // a repeller
			{
				strings[0] = 'R';
			}

			else if(graphnodes[i].type == 1) // an attractor
			{
				strings[0] = 'A';
			}

			else   //saddle
			{
				strings[0] = 'S';
			}
				
			if(graphnodes[i].labelindex+1 >= 100)
			{
				int hundred = floor((graphnodes[i].labelindex+1)/100.);
				int ten = floor((graphnodes[i].labelindex + 1 - hundred * 100)/10.);
				int number = (graphnodes[i].labelindex+1) % 10;
				strings[1] = '000'+hundred;
				strings[2] = '000'+ten;
				strings[3] = '000'+number;
				strings[4] = '\0';
			}

			else if(graphnodes[i].labelindex+1 >= 10)
			{
				int ten = floor((graphnodes[i].labelindex+1)/10.);
				int number = (graphnodes[i].labelindex+1) % 10;
				strings[1] = '000'+ten;
				strings[2] = '000'+number;
				strings[3] = '\0';
			}

			else{
				strings[1] = '001'+graphnodes[i].labelindex;
				strings[2] = '\0';
			}

			Displaylabel(posx, posy, strings);

			glPopMatrix();
		}
	}

	else{

		for(i = 0; i < cur_mcgnode_index; i++)
		{
			if(mcgnodes[i].cancelled == 1)
				continue;

			glPushMatrix();

			posx = (int)((mcgnodes[i].pos_x/4.)*(float)viewport[2])-34;
			posy = (int)((mcgnodes[i].pos_y/2.)*(float)viewport[3])-5;

			glColor3f(0,0,0);
			if(mcgnodes[i].type == 0) // a repeller
			{
				strings[0] = 'R';
			}

			else if(mcgnodes[i].type == 1) // an attractor
			{
				strings[0] = 'A';
			}

			else   //saddle
			{
				strings[0] = 'S';
			}
				
			if(mcgnodes[i].labelindex+1 >= 100)
			{
				int hundred = floor((mcgnodes[i].labelindex+1)/100.);
				int ten = floor((mcgnodes[i].labelindex + 1 - hundred * 100)/10.);
				int number = (mcgnodes[i].labelindex+1) % 10;
				strings[1] = '000'+hundred;
				strings[2] = '000'+ten;
				strings[3] = '000'+number;
				strings[4] = '\0';
			}

			else if(mcgnodes[i].labelindex+1 >= 10)
			{
				int ten = floor((mcgnodes[i].labelindex+1)/10.);
				int number = (mcgnodes[i].labelindex+1) % 10;
				strings[1] = '000'+ten;
				strings[2] = '000'+number;
				strings[3] = '\0';
			}

			else{
				strings[1] = '001'+mcgnodes[i].labelindex;
				strings[2] = '\0';
			}

			Displaylabel(posx, posy, strings);

			glPopMatrix();
		}
	}

	////Reset the viewport back to original setting
	glMatrixMode(GL_PROJECTION);		// Select The Projection Matrix
	glLoadIdentity();					// Reset The Projection Matrix

	// Calculate The Aspect Ratio Of The Window
	glOrtho(0.0f,4.0f,0.0f,2.0f,0.0f,50.0f);
}



/* --------------------------------------------------
*  Select colors for the nodes according to their 
*  types in the graph
---------------------------------------------------*/
void CConleyGraphView::SetColor(int nodetype)
{
	if(nodetype == 0)
		glColor3f(0, 1, 0);
	else if (nodetype == 1)
		glColor3f(1, 0, 0);
	else
		glColor3f(0, 0, 1);
}

/* --------------------------------------------------
*  Draw a hollow circle for highlighting the nodes 
*  that selected by user in the graph
---------------------------------------------------*/
void CConleyGraphView::DrawCircle(double cx, double cy)
{
	int i;
	double R = 0.033;
	double theta, deta ;
	deta = 2 * M_PI/49.;
	double x, y;
	theta = 0.;
	glBegin(GL_LINE_LOOP);
	for(i = 0; i < 50; i++, theta += deta)
	{
		x =  cx + R * cos(theta);
		y =  cy + R * sin(theta);

		glVertex2f(x, y);
	}
	glEnd();
}

	
void CConleyGraphView::DrawHighLights()
{
	int i;
	int node1, node2;

	glColor3f(1, 0.6, 0.3);   ////set yellow as the highlighted color
	glLineWidth(2.8);

	if(ShowMCGOn == 0)
	{
		////1. Draw the highlighted edges first
		glBegin(GL_LINES);
		for(i = 0; i < num_related_edges; i++)
		{
			node1 = graphedges[related_edges[i]].node_index1;
			node2 = graphedges[related_edges[i]].node_index2;

			////10/28/05 temporary
			if(node1 < 0 || node2 < 0)
				return;

			if(graphnodes[node1].cancelled == 1 || graphnodes[node2].cancelled == 1)
				continue;

			////Line1
			glVertex2f(graphnodes[node1].pos_x, graphnodes[node1].pos_y-0.01);
			glVertex2f(graphnodes[node2].pos_x, graphnodes[node2].pos_y-0.01);
			
			////Line2
			glVertex2f(graphnodes[node1].pos_x, graphnodes[node1].pos_y+0.01);
			glVertex2f(graphnodes[node2].pos_x, graphnodes[node2].pos_y+0.01);
		}
		glEnd();

		////2. Draw the highlighted nodes
		DrawCircle(graphnodes[picked_node].pos_x, graphnodes[picked_node].pos_y);

		for(i = 0; i < num_related_edges; i++)
		{
			node1 = graphedges[related_edges[i]].node_index1;

			if(graphnodes[node1].cancelled == 1)
				continue;

			if(node1 == picked_node)
				node1 = graphedges[related_edges[i]].node_index2;

			DrawCircle(graphnodes[node1].pos_x, graphnodes[node1].pos_y);
		}
	}

	else
	{
		////2. Draw the highlighted nodes
		DrawCircle(mcgnodes[picked_node].pos_x, mcgnodes[picked_node].pos_y);

		//for(i = 0; i < num_related_edges; i++)
		//{
		//	node1 = mcgedges[related_edges[i]].node_index1;

		//	if(mcgnodes[node1].cancelled == 1)
		//		continue;

		//	if(node1 == picked_node)
		//		node1 = mcgedges[related_edges[i]].node_index2;

		//	DrawCircle(mcgnodes[node1].pos_x, mcgnodes[node1].pos_y);
		//}
	}
}


void CConleyGraphView::GetTheHighLightNodeandEdges()
{
	int i;
	int MaxNumRelatedEdges = 20;

	////free previous allocation
	if(related_edges != NULL)
	{
		free(related_edges);
		related_edges = NULL;
	}
	
	num_related_edges = 0;

	if(picked_node < 0)
		return;
	 
	related_edges = (int *)malloc(sizeof(int)*MaxNumRelatedEdges);

	if(ShowMCGOn == 0)
	{
		////get all the related edges that incident to the node
		for(i = 0; i < cur_graphedge_index; i++)
		{
			if(graphedges[i].node_index1 == picked_node || graphedges[i].node_index2 == picked_node)
			{
				if(num_related_edges >= MaxNumRelatedEdges - 1)
				{
					MaxNumRelatedEdges += 20;
					related_edges = (int *)realloc(related_edges, sizeof(int) * MaxNumRelatedEdges);
				}

				related_edges[num_related_edges] = i;
				num_related_edges ++;
			}
		}
	}
	else
	{
	}
}


int CConleyGraphView::FindConnectedSep(int saddle, int othersing, int type)
{
	int i;
	int sep;
	int end_triangle;
	
	if(type == 0)  //othersing is a repeller
	{
		sep = separatrices[singularities[saddle].separtices].sep2;
		end_triangle = trajectories[sep][num_linesegs_curtraj[sep]-1].Triangle_ID;
		if(end_triangle >= 0 && Object.flist[end_triangle]->contain_singularity == 1)
		{
			if(Object.flist[end_triangle]->singularity_index == othersing)
				return sep;
		}
		////Maybe its neighboring triangle contains singularity
		//if(end_triangle != -1)
		//{
		//	Corner *c;

		//	for(i = 0; i < 3; i++)
		//	{
		//		c = Object.clist[end_triangle*3+i];
		//		if(Object.flist[c->ot]->contain_singularity == 1
		//			&& Object.flist[c->ot]->singularity_index == othersing)
		//			return sep;
		//	}
		//}

		sep = separatrices[singularities[saddle].separtices].sep4;
		end_triangle = trajectories[sep][num_linesegs_curtraj[sep]-1].Triangle_ID;
		if(end_triangle >= 0 && Object.flist[end_triangle]->contain_singularity == 1)
		{
			if(Object.flist[end_triangle]->singularity_index == othersing)
				return sep;
		}

		////Maybe its neighboring triangle contains singularity
		//if(end_triangle != -1)
		//{
		//	Corner *c;

		//	for(i = 0; i < 3; i++)
		//	{
		//		c = Object.clist[end_triangle*3+i];
		//		if(Object.flist[c->ot]->contain_singularity == 1
		//			&& Object.flist[c->ot]->singularity_index == othersing)
		//			return sep;
		//	}
		//}

		return -1;
	}

	else  //other singularity is an attractor
	{
		sep = separatrices[singularities[saddle].separtices].sep1;
		end_triangle = trajectories[sep][num_linesegs_curtraj[sep]-1].Triangle_ID;
		if(end_triangle >= 0 && Object.flist[end_triangle]->contain_singularity == 1)
		{
			if(Object.flist[end_triangle]->singularity_index == othersing)
				return sep;
		}
		////Maybe its neighboring triangle contains singularity
		//if(end_triangle != -1)
		//{
		//	Corner *c;

		//	for(i = 0; i < 3; i++)
		//	{
		//		c = Object.clist[end_triangle*3+i];
		//		if(Object.flist[c->ot]->contain_singularity == 1
		//			&& Object.flist[c->ot]->singularity_index == othersing)
		//			return sep;
		//	}
		//}

		sep = separatrices[singularities[saddle].separtices].sep3;
		end_triangle = trajectories[sep][num_linesegs_curtraj[sep]-1].Triangle_ID;
		if(end_triangle >= 0 && Object.flist[end_triangle]->contain_singularity == 1)
		{
			if(Object.flist[end_triangle]->singularity_index == othersing)
				return sep;
		}
		
		////Maybe its neighboring triangle contains singularity
		//if(end_triangle != -1)
		//{
		//	Corner *c;

		//	for(i = 0; i < 3; i++)
		//	{
		//		c = Object.clist[end_triangle*3+i];
		//		if(Object.flist[c->ot]->contain_singularity == 1
		//			&& Object.flist[c->ot]->singularity_index == othersing)
		//			return sep;
		//	}
		//}

		return -1;
	}

}


////Get the corresponding highlighting separatrices 1/28/06 Happy chinese new year!
void CConleyGraphView::GetHighLightSeparatrices()
{
	int i;
	int MaxNumRelatedEdges = 20;
	int node;
	int findonesep;

	////free previous allocation
	if(related_trajs != NULL)
	{
		free(related_trajs);
		related_trajs = NULL;
	}
	
	num_related_trajs = 0;

	if(num_related_edges == 0)
		return;
	 
	related_trajs = (int *)malloc(sizeof(int)*MaxNumRelatedEdges);

	////
	int singularityID = graphnodes[picked_node].singularityID;
	int type = graphnodes[picked_node].type;

	if(type == 2) //saddle
	{
		////just highlight all its separatrices
		related_trajs[0] = separatrices[singularities[singularityID].separtices].sep1;
		related_trajs[1] = separatrices[singularities[singularityID].separtices].sep2;
		related_trajs[2] = separatrices[singularities[singularityID].separtices].sep3;
		related_trajs[3] = separatrices[singularities[singularityID].separtices].sep4;

		num_related_trajs = 4;
	}

	else //repeller or attractor
	{
		for(i = 0; i < num_related_edges; i++)
		{
			if(graphnodes[graphedges[related_edges[i]].node_index1].LimitCycleID >= 0
				|| graphnodes[graphedges[related_edges[i]].node_index2].LimitCycleID >= 0)
				continue;

			if(type == 0)
				node = graphedges[related_edges[i]].node_index2; //because the edge has direction
			else
				node = graphedges[related_edges[i]].node_index1;

			findonesep = FindConnectedSep(graphnodes[node].singularityID, singularityID, type);

			if(findonesep >= 0)
			{
			
				if(num_related_trajs >= MaxNumRelatedEdges - 1)
				{
					MaxNumRelatedEdges += 20;
					related_trajs = (int *)realloc(related_trajs, sizeof(int) * MaxNumRelatedEdges);
				}

				related_trajs[num_related_trajs] = findonesep;
				num_related_trajs ++;
			}
		}
	}
}


/*--------------------------------------------------------------------------------*/
/*Mouse selection routine
/*--------------------------------------------------------------------------------*/
void CConleyGraphView::HitProcessforGraph(double ss, double tt)
{
	GLuint  selectBuffer[SELECTBUFFERSIZE];
	GLint	vp[4] = {0, 0 ,4, 2};
	int hits;

	////Build the selection buffer here
	glSelectBuffer(SELECTBUFFERSIZE, selectBuffer);

   	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);

	gluPickMatrix(ss, tt, 0.005, 0.005, vp );  ////set a larger pick window for element selection
	glOrtho(0, 4,  0, 2,  0, 50);

	DrawNodes(GL_SELECT);
    //picked_node = -1;

	hits = glRenderMode(GL_RENDER);

	if(hits>0)
	{
		picked_node = selectBuffer[3]-1;
	}
	else{
		picked_node = -1;
	}
	
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode( GL_MODELVIEW );
}

 
void CConleyGraphView::HitProcessforCancel(double ss, double tt)
{
	GLuint  selectBuffer[SELECTBUFFERSIZE];
	GLint	vp[4] = {0, 0 ,4, 2};
	int hits;

	////Build the selection buffer here
	glSelectBuffer(SELECTBUFFERSIZE, selectBuffer);

   	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);

	gluPickMatrix(ss, tt, 0.005, 0.005, vp );  ////set a larger pick window for element selection
	glOrtho(0, 4,  0, 2,  0, 50);

	DrawNodes(GL_SELECT);
    picked_node = -1;

	hits = glRenderMode(GL_RENDER);

	if(hits>0)
	{
		if(PairCancelOn == 1)
		{

			////More general cancellation
			if(SelectRepellerOrAttractor == 0 )
			{
				AddToRepellerList(selectBuffer[3] - 1);
			}
			else 
			{
				AddToAttractorList(selectBuffer[3] - 1);
			}

			////we do not include limit cycle cancellation now 1/2/06
		}

		picked_node = selectBuffer[3]-1;
	}
	else{
		picked_node = -1;
	}
	
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode( GL_MODELVIEW );
}


void CConleyGraphView::InitFlags()
{
	PairCancelOn = 0;
	SelectRepellerOrAttractor = 0;

	PairCounter = 0;
    cancel_node1 = cancel_node2 = -1;
	if(cancel_nodes != NULL)
		free(cancel_nodes);
    cancel_nodes = NULL;

    num_cancel_nodes = 0;

	picked_node = -1;                 
    if(related_edges != NULL)       
		free(related_edges);
    num_related_edges = 0;

	ShowMCGOn = 0;

	show3DOn = 0;

}


////Routines for multiple repeller or attractor selection in the Conley graph

////Here index is the index of the node in the Conley graph
void CConleyGraphView::AddToRepellerList(int index)
{
	if(graphnodes[index].type != 0)
	{
		MessageBox("This is not a repeller!", "Error", MB_OK);
		return;
	}
	
	////check repeatation
	int i;

	for(i = 0; i < NumRepellers; i++)
	{
		if(index == repeller_nodes[i])
		{
			MessageBox("The node has been selected!", "Error", MB_OK);
			return;
		}
	}

	repeller_nodes = Extend_link(repeller_nodes, NumRepellers);
	    
	repeller_nodes[NumRepellers] = index;

	NumRepellers++;
}
	
void CConleyGraphView::AddToAttractorList(int index)
{
	if(graphnodes[index].type != 1)
	{
		MessageBox("This is not a attractor!", "Error", MB_OK);
		return;
	}

	////check repeatation
	int i;

	for(i = 0; i < NumAttractors; i++)
	{
		if(index == attractor_nodes[i])
		{
			MessageBox("The node has been selected!", "Error", MB_OK);
			return;
		}
	}

	attractor_nodes = Extend_link(attractor_nodes, NumAttractors);

	attractor_nodes[NumAttractors] = index;

	NumAttractors++;
}


////Somewhat repeat the work in class CGlView
void CConleyGraphView::ClearRepellandAttractList()
{
	if(repeller_nodes != NULL)
		free(repeller_nodes);
	if(attractor_nodes != NULL)
		free(attractor_nodes);

	NumRepellers = 0;
	NumAttractors = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
/////   3D display
#include "traceball.h"
int view_mode=0;  // 0 = othogonal, 1=perspective
int ACSIZE = 1; // for antialiasing
double min_dot;
int which_menu=1;  
int mouse_mode = -2;  // -2=no action, -1 = down, 0 = zoom, 1 = rotate x, 2 = rotate y, 3 = tranlate x, 4 = translate y, 5 = cull near 6 = cull far
int mouse_button = -1; // -1=no button, 0=left, 1=middle, 2=right
int last_x, last_y;
int orientation;

double this_depth;
float *depths;
double modelview_matrix[16], projection_matrix[16];
GLint viewport[4];
extern double zoom_factor ;
icVector3 rot_center;

int cull_mode=0;  // 0 = cull back, 1= cull front, 2=cull none
double radius_factor=0.9, translate_x=0.0, translate_y=0.0;
double vd_n, vd_f;
double spin_longitude=0.0, spin_latitude=0.0;
int red_bits, green_bits, blue_bits, alpha_bits, index_bits, depth_bits;

FILE *this_file;

float s_old, t_old;
static Quaternion rvec;

extern QuadMesh *quadmesh;


//struct jitter_struct{
//	double x;
//	double y;
//} jitter_para;

//jitter_struct ji1[1] = {{0.0, 0.0}};
//jitter_struct ji16[16] = {{0.125, 0.125}, {0.375, 0.125}, {0.625, 0.125}, {0.875, 0.125}, 
//						  {0.125, 0.375}, {0.375, 0.375}, {0.625, 0.375}, {0.875, 0.375}, 
//						  {0.125, 0.625}, {0.375, 0.625}, {0.625, 0.625}, {0.875, 0.625}, 
//
//
//
//{0.125, 0.875}, {0.375, 0.875}, {0.625, 0.875}, {0.875, 0.875}, };


void multmatrix(const Matrix m)
{ 
  int i,j, index = 0;

  GLfloat mat[16];

  for ( i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      mat[index++] = m[i][j];

  glMultMatrixf (mat);
}


//void mat_ident( Matrix m)
//{
//	int i;
//
//	for (i = 0; i <= 3; i++) {
//		m[i][0] = 0.0;
//		m[i][1] = 0.0;
//		m[i][2] = 0.0;
//		m[i][3] = 0.0;
//		m[i][i] = 1.0;
//	}
//}

void set_view(GLenum mode)
{
	icVector3 up, ray, view;
	GLfloat light_ambient0[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.9, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_diffuse2[] = { 0.8, 0.8, 1.0, 1.0 };
	GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 0.0, 0.0, 2.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);
	glLightfv(GL_LIGHT2, GL_AMBIENT, light_ambient2);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, light_diffuse2);
	glLightfv(GL_LIGHT2, GL_SPECULAR, light_specular2);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (view_mode == 0)
		glOrtho(-radius_factor*zoom_factor, radius_factor*zoom_factor, 
			-radius_factor*zoom_factor, radius_factor*zoom_factor, 0.0, 40.0);
//		glOrtho(-1.01*zoom_factor, 1.01*zoom_factor, -1.01*zoom_factor, 1.01*zoom_factor, 0.0, 40.0);
	else
		gluPerspective(45.0, 1.0, 0.1, 40.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 2.5;
	light_position[1] = 2.0;
	light_position[2] = 3.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 1.0;
	glLightfv(GL_LIGHT1, GL_POSITION, light_position);
	light_position[0] = -2.1;
	light_position[1] = -1.0;
	light_position[2] = 1.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);

}

void set_scene(GLenum mode)
{
	glTranslatef(0.0, 0.0, -3.5);
	glTranslatef(rot_center.entry[0]-Object.center.entry[0], rot_center.entry[1]-Object.center.entry[1], rot_center.entry[2]-Object.center.entry[2]);
	multmatrix( rotmat );

	glTranslatef(-rot_center.entry[0]+Object.center.entry[0], -rot_center.entry[1]+Object.center.entry[1], -rot_center.entry[2]+Object.center.entry[2]);
	glScalef(1.0/Object.radius, 1.0/Object.radius, 1.0/Object.radius);
	glRotatef(-60, 1, 0, 0);
	//glRotatef(180, 0, 1, 0);
	glTranslatef(-Object.center.entry[0], -Object.center.entry[1], -Object.center.entry[2]);
}


void display_object(GLenum mode)
{
	unsigned int i, j;
	QuadCell *face;
	int *verts;
	icVector3 up, ray, view;

	switch(which_menu) {
		case 0:
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glShadeModel (GL_FLAT);
			glDisable(GL_LIGHTING);
			glDisable(GL_LIGHT0);
			for (i=0; i<quadmesh->nfaces; i++) {
				if (mode == GL_SELECT)
					glLoadName(i+1);
				glColor3ub(0, 0, 255);
				face = quadmesh->quadcells[i];
				verts = face->verts;
				//glBegin(GL_POLYGON);
				glBegin(GL_LINE_LOOP);
				for (j=0; j<face->nverts; j++) {
					glVertex3d(quadmesh->quad_verts[verts[j]]->x, 
						quadmesh->quad_verts[verts[j]]->y, 0);
				}
				glEnd();
			}

			//display the vertices only
			//glColor3f(1, 0, 0);
			//glPointSize(3.);
			//glBegin(GL_POINTS);
			//for(i = 0; i < Object.nverts; i++)
			//{
			//	//if(i == Object.nverts-1 || i == Object.nverts - 2
			//	//	||i == Object.nverts - 1 - (16 + 2))
			//	glVertex3f(Object.vlist[i]->x, Object.vlist[i]->y, Object.vlist[i]->z);
			//}
			//glEnd();
			break;

		case 1:
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glShadeModel(GL_SMOOTH);
			glEnable(GL_LIGHTING);
			glEnable(GL_LIGHT0);
			glEnable(GL_LIGHT1);
			glEnable(GL_LIGHT2);
			for (i=0; i<quadmesh->nfaces; i++) {

				face = quadmesh->quadcells[i];
				verts = face->verts;
				glBegin(GL_POLYGON);
				for (j=0; j<face->nverts; j++) {
					glNormal3d(0, 0, 1.);
					glVertex3d(quadmesh->quad_verts[verts[j]]->x, 
						quadmesh->quad_verts[verts[j]]->y, 0);
				}
				glEnd();
			}
				 
			break;
	}
}

void display(void)
{
  int jitter;

//  glClearColor (0.0, 0.0, 0.0, 1.0);  // background in my original code
//  glClearColor (0.2, 0.2, 0.8, 1.0);  // background for rendering color coding
  glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting

  glGetIntegerv (GL_VIEWPORT, viewport);
 
//  glClear(GL_ACCUM_BUFFER_BIT);
 // for (jitter = 0; jitter < ACSIZE; jitter++) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	  set_view(GL_RENDER);
   /* glPushMatrix ();
		switch(ACSIZE){
		case 1:
			glTranslatef (ji1[jitter].x*2.0*radius_factor/viewport[2], ji1[jitter].y*2.0*radius_factor/viewport[3], 0.0);
			break;

		case 16:
			glTranslatef (ji16[jitter].x*2.0*radius_factor/viewport[2], ji16[jitter].y*2.0*radius_factor/viewport[3], 0.0);
			break;

		default:
			glTranslatef (ji1[jitter].x*2.0*radius_factor/viewport[2], ji1[jitter].y*2.0*radius_factor/viewport[3], 0.0);
			break;
		}
	*/
	  			

	set_scene(GL_RENDER);
    display_object (GL_RENDER);
/*    glPopMatrix ();
	  glutSwapBuffers();
    glAccum(GL_ACCUM, 1.0/ACSIZE);
  }
  glAccum (GL_RETURN, 1.0);
	glReadBuffer(GL_BACK);
	glReadPixels(0, 0, win_width, win_height, GL_DEPTH_COMPONENT, GL_FLOAT, depths);
  */
	//glutSwapBuffers();
 //	glFinish();

}

#include "EvenlyStreamlines.h"
extern EvenStreamlinePlace *major, *minor;
extern StreetVis *roadnetvis;

void display_roads_3d(GLenum mode)
{
	if(roadnetvis == NULL)
		return;

	/*clear the back ground*/
	//glDisable(GL_TEXTURE_2D);
	//glEnable(GL_COLOR_MATERIAL);
	//glClearColor(0.8, .8, 0.9, 1.);
	//glClear(GL_COLOR_BUFFER_BIT);

	int i, j;
	RoadLineSeg *curlineseg;
	OneRoadVis *curroad;

	/*fill up the colors*/
	for(i=0; i<roadnetvis->nmajDirRoads; i++)
	{
		curroad = roadnetvis->majDirRoads[i];

		if(major->evenstreamlines->trajs[i]->roadtype != HIGHWAY
			&& major->evenstreamlines->trajs[i]->roadtype != MAJOR)
			continue;

		if(major->evenstreamlines->trajs[i]->roadtype == HIGHWAY)
			glColor3f(1, 1, 0);
		else
			glColor3f(1, .5, 0);
		
		//for(j=0; j<curroad->nroadlines1; j++)
		//{
		//	glBegin(GL_POLYGON);
		//	curlineseg=curroad->roadline1[j];
		//	glVertex3f(curlineseg->start[0], curlineseg->start[1]);
		//	glVertex3f(curlineseg->end[0], curlineseg->end[1]);

		//	curlineseg=curroad->roadline2[j];
		//	glVertex3f(curlineseg->end[0], curlineseg->end[1]);
		//	glVertex3f(curlineseg->start[0], curlineseg->start[1]);
		//	glEnd();
		//}

		for(j=0; j<curroad->nroadlines1-1; j++)
		{
			glBegin(GL_POLYGON);
			curlineseg=curroad->roadline1[j];
			glVertex3f(curlineseg->start[0], curlineseg->start[1], 0.01);
			curlineseg=curroad->roadline1[j+1];
			glVertex3f(curlineseg->start[0], curlineseg->start[1], 0.01);

			curlineseg=curroad->roadline2[j+1];
			glVertex3f(curlineseg->end[0], curlineseg->end[1], 0.01);
			curlineseg=curroad->roadline2[j];
			glVertex3f(curlineseg->end[0], curlineseg->end[1], 0.01);
			glEnd();
		}
	}

	for(i=0; i<roadnetvis->nminDirRoads; i++)
	{
		curroad = roadnetvis->minDirRoads[i];

		if(minor->evenstreamlines->trajs[i]->roadtype != HIGHWAY
			&& minor->evenstreamlines->trajs[i]->roadtype != MAJOR)
			continue;

		if(minor->evenstreamlines->trajs[i]->roadtype == HIGHWAY)
			glColor3f(1, 1, 0);
		else
			glColor3f(1, .5, 0);
		
		glBegin(GL_POLYGON);
		for(j=0; j<curroad->nroadlines1; j++)
		{
			curlineseg=curroad->roadline1[j];
			glVertex3f(curlineseg->start[0], curlineseg->start[1], 0.01);
			glVertex3f(curlineseg->end[0], curlineseg->end[1], 0.01);

			curlineseg=curroad->roadline2[j];
			glVertex3f(curlineseg->end[0], curlineseg->end[1], 0.01);
			glVertex3f(curlineseg->start[0], curlineseg->start[1], 0.01);
		}
		glEnd();
	}


	/*draw the enclosure of the roads*/
	glDisable(GL_LIGHTING);
	glColor3f(0, 0, 0);

	/*major roads*/
	for(i=0; i<roadnetvis->nmajDirRoads; i++)
	{
		curroad = roadnetvis->majDirRoads[i];
		if(major->evenstreamlines->trajs[i]->roadtype == MINOR)
			glLineWidth(1.0);
		//else if(major->evenstreamlines->trajs[i]->roadtype == PLAIN)
		//	glLineWidth(1.2);
		else if(major->evenstreamlines->trajs[i]->roadtype == MAJOR)
			glLineWidth(1.5);
		else
			glLineWidth(2.);

		/*display road 1*/
		for(j=0; j<curroad->nroadlines1; j++)
		{
			curlineseg = curroad->roadline1[j];
			glBegin(GL_LINES);
			glVertex3f(curlineseg->start[0], curlineseg->start[1], 0.01);
			glVertex3f(curlineseg->end[0], curlineseg->end[1], 0.01);
			glEnd();
		}
		/*display road 2*/
		for(j=0; j<curroad->nroadlines2; j++)
		{
			curlineseg = curroad->roadline2[j];
			glBegin(GL_LINES);
			glVertex3f(curlineseg->start[0], curlineseg->start[1], 0.01);
			glVertex3f(curlineseg->end[0], curlineseg->end[1], 0.01);
			glEnd();
		}
	}

	/*minor roads*/
	for(i=0; i<roadnetvis->nminDirRoads; i++)
	{
		curroad = roadnetvis->minDirRoads[i];

		if(minor->evenstreamlines->trajs[i]->roadtype == MINOR)
			glLineWidth(1.0);
		//else if(minor->evenstreamlines->trajs[i]->roadtype == PLAIN)
		//	glLineWidth(1.2);
		else if(minor->evenstreamlines->trajs[i]->roadtype == MAJOR)
			glLineWidth(1.5);
		else
			glLineWidth(2.);

		/*display road 1*/
		for(j=0; j<curroad->nroadlines1; j++)
		{
			curlineseg = curroad->roadline1[j];
			glBegin(GL_LINES);
			glVertex3f(curlineseg->start[0], curlineseg->start[1], 0.01);
			glVertex3f(curlineseg->end[0], curlineseg->end[1], 0.01);
			glEnd();
		}
		/*display road 2*/
		for(j=0; j<curroad->nroadlines2; j++)
		{
			curlineseg = curroad->roadline2[j];
			glBegin(GL_LINES);
			glVertex3f(curlineseg->start[0], curlineseg->start[1], 0.01);
			glVertex3f(curlineseg->end[0], curlineseg->end[1], 0.01);
			glEnd();
		}
	}
}




BEGIN_MESSAGE_MAP(CConleyGraphView, CWnd)
	ON_WM_MOUSEMOVE()
END_MESSAGE_MAP()

void CConleyGraphView::OnMouseMove(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default

	CWnd::OnMouseMove(nFlags, point);
}
