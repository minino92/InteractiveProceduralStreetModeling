/*
   This file is typically implementing the routines for saving or loading the project setting
*/

#include "stdafx.h"

//#include "VFAnalysis.h"
//#include "LocalTracing.h"
//#include "topologyedit.h"
//
//#include "BmpProcess.h"
//
//
#include "VFDataStructure.h"
#include "shareinterfacevars.h"
////
//#include "ImgBoundaryExtract.h"
//#include "BoundaryBasedTenGen.h"
//
#include "SketchDesign.h"
//
//#include "loadmaps.h"
//
#include "scalardesign.h"
//
#include "tensorvis.h"

#include "tensoranalysis.h"
#include "evenlystreamlines.h"
#include "tensordesign.h"

#include "LocalRegionEdit.h"
#include "regionsmooth_quad.h"

//#include "tensorvis.h"

#include "SaveProjectSetting.h"

extern SharedInterfaceVars sharedvars;
extern QuadMesh *quadmesh;
extern bool majorroadsexisted;
extern unsigned char *fittedmap1;
extern bool flag_loadmap;
extern void get_mask_map(double xstart, double xend, double ystart, double yend,
				  unsigned char *map, int width, int height);


extern void clean();

/*
    What kind of information we need to save for a project:
	1. design elements if exist (or a tensor field per vertex, remember the resolution)
	2. the sketches if exist
	3. brush strokes if exist
	4. obtained major road network if exists (store the tensor lines directly; or
	store the graph; seed points)
	5. some setting (parameters: tracing density, crossing river parameters,
	filter parameter, seed point parameters, brush region parameters, etc.)
	6. other information

	The format of the file will be:
	#comments: file_name
	#mesh resolution: XDIM x YDIM
	#design elements: e | v | n (e--element based, v--vertex based, n--none)
	   (NOTE: e-- file_name_elem.txt;
	          v-- file_name_perver.txt;
			  n-- NULL;
		)
	#sketches: y | n (y-- exist, n--none)
	   (NOTE: y -- file_name_sketch.skt;
	          n -- NULL;
	   )
    #brush strokes: y | n (y--exist, n--none)
	   (NOTE: y -- file_name_brushes.bru;
	          n -- NULL;
	   )
    #major network: t | g | n 
	   (NOTE: t -- store as tensor lines, file_name_majTenLines.mjt;
	          g -- store as graph,        file_name_majNetwork.gra;
			  n -- NULL;
	   )

	#other setting: y | n
	   (NOTE: y -- exists: file_name_setting.tst
	          n -- NULL;
	   )

   The format of other setting (.tst file)

    #tracing density:
	1: maj_density (float), min_density (float) (major roads)
	2: maj_density (float), min_density (float) (minor roads)
	#density map: y | n
	   (NOTE: y -- density_map.bmp
	          n -- NULL
	   )
    #cross river: y (float) | n

	#brush region width: (float)

	#filter size: (int)

	#noise setting: freq (float), amp (float)
*/

bool save_current_project(char *filename)
{
	FILE *fp=fopen(filename, "w");

	if(fp==NULL)
		return false;

	fprintf(fp, "#comments:%s\n", filename);

	/*  obtain the name except for the extension  */
	int i, pos;
	for(i=strlen(filename)-1; i>=0; i--)
	{
		if(filename[i] == '.')
		{
			pos=i;
			break;
		}
	}

	char str[255];
	char *commonfilename=(char*)malloc(sizeof(char)*(pos+1));
	for(i=0;i<pos;i++)
	{
		commonfilename[i]=filename[i];
	}
	commonfilename[i]='\0';
	
	/*  #mesh resolution: XDIM x YDIM  */
	fprintf(fp,"#mesh resolution: %d, %d, %f, %f, %f, %f\n", quadmesh->XDIM, quadmesh->YDIM,
		quadmesh->xstart, quadmesh->xend, quadmesh->ystart, quadmesh->yend);

	/* #design elements: e | v | n (e--element based, v--vertex based, n--none) */

	char design_elem_char = 'n'; //this value will be decided by usre selection
	if(sharedvars.rdSaveProjDesignElem==1)
		design_elem_char='e';
	else if(sharedvars.rdSaveProjDesignElem==2)
		design_elem_char='v';

	fprintf(fp, "#design elements: %c\n", design_elem_char);

	if(design_elem_char == 'e')
	{
		fprintf(fp, "%s_elem.txt\n", commonfilename);

		/*  call the save element function to save the elements  */
		sprintf(str, "%s_elem.txt", commonfilename);
		
		save_tenField_elems(str);
	}
	else if(design_elem_char == 'v')
	{
		fprintf(fp, "%s_perver.txt\n", commonfilename);

		/*  call the save tensor field per vertex routine  */
		sprintf(str, "%s_perver.txt", commonfilename);

		save_tenField_perVer(str);
	}

	//fprintf(fp, "\n");


	/* #sketches: y | n (y-- exist, n--none)
	   (NOTE: y -- file_name_sketch.skt;
	          n -- NULL;
	   )*/
	char sketch_char = 'n';  //decided by the user selection
	if(sharedvars.rdSaveProjSketches==1)
		sketch_char='y';

	fprintf(fp, "#sketches: %c\n", sketch_char);

	if(sketch_char == 'y')
	{
		fprintf(fp, "%s_sketch.skt\n", commonfilename);

		/*  call the save sketch routine */
		sprintf(str, "%s_sketch.skt", commonfilename);

		save_sketch_brushes(str);
	}

	//fprintf(fp, "\n");

	/*
    #brush strokes: y | n (y--exist, n--none)
	   (NOTE: y -- file_name_brushes.bru;
	          n -- NULL;
	   )
	*/

	char brush_char = 'n';
	if(sharedvars.rdSaveProjBrushes==1)
		brush_char='y';

	fprintf(fp, "#brush strokes: %c\n", brush_char);

	if(brush_char == 'y')
	{
		fprintf(fp, "%s_brushes.bru\n", commonfilename);

		/*  call the save brushes routine  */
		sprintf(str, "%s_brushes.bru", commonfilename);
	}

	//fprintf(fp, "\n");

	/*
    #major network: t | g | n 
	   (NOTE: t -- store as tensor lines, file_name_majTenLines.mjt;
	          g -- store as graph,        file_name_majNetwork.gra;
			  n -- NULL;
	   )
	*/

	char major_net_char = 'n';
	if(sharedvars.rdSaveProjMajRoadNetwork==1)
		major_net_char='t';
	else if(sharedvars.rdSaveProjMajRoadNetwork==2)
		major_net_char='g';

	fprintf(fp, "#major network: %c\n", major_net_char);

	if(major_net_char == 't')
	{
		fprintf(fp, "%s_majTenLines.mjt\n", commonfilename);

		/*  call the save major roads as tensor line routine  */
		sprintf(str, "%s_majTenLines.mjt", commonfilename);

		save_obtained_majRoads_tenLines(str);
	}
	else if(major_net_char == 'g')
	{
		fprintf(fp, "%s_majNetwork.gra\n", commonfilename);

		/*  call the save major roads as street network routine */
		sprintf(str, "%s_majNetwork.gra", commonfilename);
	}

	//fprintf(fp, "\n");
	
	/*
    #street network: y | n 
	   (NOTE: y --  file_name_streets.gra;
			  n -- NULL;
	   )
	*/

	char street_net_char = 'n';
	if(sharedvars.rdSaveProjStreetNetwork==1)
		street_net_char='y';

	fprintf(fp, "#street network: %c\n", street_net_char);

	if(street_net_char == 'y')
	{
		fprintf(fp, "%s_street.gra\n", commonfilename);

		/*  call the save major roads as tensor line routine  */
		sprintf(str, "%s_street.gra", commonfilename);

		save_cur_street_network(str);
	}

	//fprintf(fp, "\n");

	/*
	#other setting: y | n
	   (NOTE: y -- exists: file_name_setting.tst
	          n -- NULL;
	   )
	*/

	char otherset_char = 'n';
	if(sharedvars.rdSaveProjOtherSetting==1)
		otherset_char='y';

	fprintf(fp, "#other setting: %c\n", otherset_char);

	if(otherset_char == 'y')
	{
		fprintf(fp, "%s_setting.tst\n", commonfilename);

		/*  call the save other setting routine  */
		sprintf(str, "%s_setting.tst", commonfilename);
	}

	fclose(fp);

	return true;
}


/*
    Load a previously saved project
*/

bool load_a_project(char *filename)
{
	FILE *fp=fopen(filename, "r");

	if(fp == NULL)
		return false;

	char projfilename[255];
	fscanf(fp, "#comments:%s\n", projfilename);

	/*  obtain the name except for the extension  */
	int i, pos;
	for(i=strlen(projfilename)-1; i>=0; i--)
	{
		if(projfilename[i] == '.')
		{
			pos=i;
			break;
		}
	}

	char str[255];
	char *commonfilename=(char*)malloc(sizeof(char)*(pos+1));
	for(i=0;i<pos;i++)
	{
		commonfilename[i]=projfilename[i];
	}
	commonfilename[i]='\0';
	
	/*  #mesh resolution: XDIM x YDIM  */
	int xdim, ydim;
	float xstart, xend, ystart, yend;
	fscanf(fp,"#mesh resolution: %d, %d, %f, %f, %f, %f\n", &xdim, &ydim,
		&xstart, &xend, &ystart, &yend);

	/*  create a mesh */
	if(quadmesh!=NULL)
	{
		quadmesh->finalize_quad_cells();
		quadmesh->finalize_quad_verts();
		delete quadmesh;
		quadmesh=NULL;
	}

	quadmesh=new QuadMesh(xdim, ydim, xstart, xend, ystart, yend);
	clean();

	if(flag_loadmap && fittedmap1!=NULL)
	{
		get_mask_map(xstart, xend, ystart, yend, fittedmap1, 512, 512);
	}

	char design_elem_char; //this value will be decided by usre selection
	fscanf(fp, "#design elements: %c\n", &design_elem_char);

	if(design_elem_char == 'e')
	{
		fscanf(fp, "%s\n", str);

		/*  call the load element function to load the elements  */
		load_tenField_elems(str);
		cal_tensorvals_quad_inReg();
	}
	else if(design_elem_char == 'v')
	{
		fscanf(fp, "%s\n", str);

		/*  call the load tensor field per vertex routine  */

		load_tenField_perVer(str);
	}

	//fscanf(fp, "\n");


	/* #sketches: y | n (y-- exist, n--none)
	   (NOTE: y -- file_name_sketch.skt;
	          n -- NULL;
	   )*/
	char sketch_char;  //decided by the user selection
	fscanf(fp, "#sketches: %c\n", &sketch_char);

	if(sketch_char == 'y')
	{
		fscanf(fp, "%s\n", str);

		/*  call the load sketch routine */
		load_sketch_brushes(str);
	}


	/*
    #brush strokes: y | n (y--exist, n--none)
	   (NOTE: y -- file_name_brushes.bru;
	          n -- NULL;
	   )
	*/

	char brush_char ;
	fscanf(fp, "#brush strokes: %c\n", &brush_char);

	if(brush_char == 'y')
	{
		fscanf(fp, "%s\n", str);

		/*  call the load brushes routine  */
	}

	/*
    #major network: t | g | n 
	   (NOTE: t -- store as tensor lines, file_name_majTenLines.mjt;
	          g -- store as graph,        file_name_majNetwork.gra;
			  n -- NULL;
	   )
	*/

	char major_net_char ;
	fscanf(fp, "#major network: %c\n", &major_net_char);

	if(major_net_char == 't')
	{
		fscanf(fp, "%s\n", str);

		/*  call the load major roads as tensor line routine  */
		load_majRoads_tenLines(str);

		majorroadsexisted=true;
	}
	else if(major_net_char == 'g')
	{
		fscanf(fp, "%s\n", str);

		/*  call the save major roads as street network routine */
	}

	/*
    #street network: y | n 
	   (NOTE: y --  file_name_streets.gra;
			  n -- NULL;
	   )
	*/

	char street_net_char;
	fscanf(fp, "#street network: %c\n", &street_net_char);

	if(street_net_char == 'y')
	{
		fprintf(fp, "%s\n", str);

		/*  call the save major roads as tensor line routine  */

		load_a_street_network(str);
	}


	/*
	#other setting: y | n
	   (NOTE: y -- exists: file_name_setting.tst
	          n -- NULL;
	   )
	*/

	char otherset_char ;
	fscanf(fp, "#other setting: %c\n", &otherset_char);

	if(otherset_char == 'y')
	{
		fscanf(fp, "%s\n", str);

		/*  call the load other setting routine  */
		sprintf(str, "%s_setting.tst", commonfilename);
	}

	fclose(fp);
	return true;
}