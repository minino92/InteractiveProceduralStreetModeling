
/******************************************************/
/*                                                    */
/*                     LoadMaps.cpp                   */


#include "stdafx.h"
#include "loadmaps.h"

#include "VFDataStructure.h"

//#include "BmpProcess.h"


/***************************************************/
/*                                                 */
/*        The variables for the loaded maps       */

/*   Geographical map:  we try to load a geograph map from file*/
unsigned char *map1=NULL;
unsigned char *fittedmap1=NULL;
unsigned char *displaymap=NULL;
unsigned char *streetmapbackground=NULL;
int boundindex1, boundindex2;
bool ylargerthanx=false; /*false: width>=height;  true: width<height*/
int map1_w, map1_h;

/*      vegetation map     */
unsigned char *vegmap=NULL;
int vegmap_w, vegmap_h;
unsigned char *vegmap_fit=NULL;
unsigned char *vegmap_disp=NULL;
bool flag_loadvegmap=false;

/*       Population density map      */
unsigned char *popdensitymap=NULL;
int popdensitymap_w, popdensitymap_h;
unsigned char *popdensitymap_fit=NULL;
unsigned char *popdensitymap_disp=NULL;

/*       Hieght field       */
unsigned char *heightfield=NULL;
int heightfield_w, heightfield_h;
unsigned char *heightfield_fit=NULL;
unsigned char *heightfield_disp=NULL;


/*       A variable for density table retrieval      */
//typedef struct DensityTableOneRow{
//	double 
//}


extern QuadMesh *quadmesh;


/*
based on the loaded map, we decide which vertices are in land, which are in
the water
*/
void get_mask_map(double xstart, double xend, double ystart, double yend,
				  unsigned char *map, int width, int height)
{
	double xrang = xend-xstart;
	double yrang = yend-ystart;
	double dx=xrang/(width-1);
	double dy=yrang/(height-1);
	int i;

	for(i=0; i<quadmesh->nverts; i++)
	{
		int c=(quadmesh->quad_verts[i]->x-xstart)/dx;
		int r=(quadmesh->quad_verts[i]->y-ystart)/dy;
		//int c=(quadmesh->quad_verts[i]->x-xstart)/xrang*(width-1);
		//int r=(quadmesh->quad_verts[i]->y-ystart)/yrang*(height-1);

		if(r>=height) r=height-1;
		if(c>=width) c=width-1;

		if(r<0) r=0;
		if(c<0) c=0;

		int id=(r*(width)+c);

		//if((map[3*id]+map[3*id+1]+map[3*id+2])<255)
		if(map[3*id]<100)
			quadmesh->quad_verts[i]->inland = true;
		else
			quadmesh->quad_verts[i]->inland = false;
	}
}


/*we use the simple nearest neighborhood approximation
NOTE: we also want to maintain the ratio of the original image
*/
void get_fitted_map(unsigned char *map, int width, int height,
					unsigned char *fittedmap, int fitwidth, int fitheight)
{
	int i,j;
	int approxi, approxj;
	double scalex=(double)width/fitwidth;
	double scaley=(double)height/fitheight;

	double scale;
	bool xyflag=false;  //default is x
	
	int offset_w, offset_h;

	offset_w=offset_h=0;

	if(scalex>scaley)
	{
		scale=scalex;
		int temp_height=(int)(fitheight/(scalex/scaley));
		offset_h=(int)(fitheight-temp_height)/2.-1;
		ylargerthanx=false;
		//boundindex1=offset_h;
		//boundindex2=fitheight-offset_h-1;
		boundindex1=0;
		boundindex2=fitheight-1;
	}
	else
	{
		scale=scaley;
		xyflag=true;    //width<height
		ylargerthanx=true;
		int temp_width=(int)(fitwidth/(scaley/scalex));
		offset_w=(int)(fitwidth-temp_width)/2.-1;
		//boundindex1=offset_w;
		//boundindex2=fitwidth-offset_w-1;
		boundindex1=0;
		boundindex2=fitwidth-1;
	}

	/*initialize all the pixels as white*/
	for(i=0; i<fitheight; i++)
	{
		for(j=0; j<fitwidth; j++)
		{
			//fittedmap[4*(i*fitwidth+j)]  =
			//fittedmap[4*(i*fitwidth+j)+1]=
			//fittedmap[4*(i*fitwidth+j)+2]=255;
			//fittedmap[4*(i*fitwidth+j)+3]=20;
			fittedmap[3*(i*fitwidth+j)]  =
			fittedmap[3*(i*fitwidth+j)+1]=
			fittedmap[3*(i*fitwidth+j)+2]=255;
		}
	}


	//for(i=offset_h; i<fitheight-offset_h; i++)
	for(i=0; i<fitheight; i++)
	{
		//approxi=(int)((i-offset_h)*scale);
		approxi=(int)(i*scale);
		//for(j=offset_w; j<fitwidth-offset_w; j++)
		for(j=0; j<fitwidth; j++)
		{
			//approxj=(int)((j-offset_w)*scale);
			approxj=(int)(j*scale);

			if(3*(approxi*width+approxj)<3*height*width)
			{
				//fittedmap[4*(i*fitwidth+j)]  =map[3*(approxi*width+approxj)];
				//fittedmap[4*(i*fitwidth+j)+1]=map[3*(approxi*width+approxj)+1];
				//fittedmap[4*(i*fitwidth+j)+2]=map[3*(approxi*width+approxj)+2];
				//fittedmap[4*(i*fitwidth+j)+3]=20;
				fittedmap[3*(i*fitwidth+j)]  =map[3*(approxi*width+approxj)];
				fittedmap[3*(i*fitwidth+j)+1]=map[3*(approxi*width+approxj)+1];
				fittedmap[3*(i*fitwidth+j)+2]=map[3*(approxi*width+approxj)+2];
			}
		}
	}

	/*we also need to record the indexes of the two boundaries of the map 11/13/2007*/
}



/*
	based on the loaded map, we compute the density value for each vertex
	NOTE: we assume that the input image is a gray image
*/
void get_density_value(double xstart, double xend, double ystart, double yend,
				  unsigned char *map, int width, int height)
{
	double xrang = xend-xstart;
	double yrang = yend-ystart;
	double dx=xrang/(width-1);
	double dy=yrang/(height-1);
	int i;

	/*        We first find the largest and smallest value       */
	unsigned char smallest, largest;
	smallest=255;
	largest=0;
	for(i=0;i<width*height*3;i+=3)
	{
		if(map[i]>largest)  largest=map[i];
		if(map[i]<smallest) smallest=map[i];
	}

	/*       Second, we compute a density value for each vertex     */
	unsigned char val_range=largest-smallest;

	for(i=0; i<quadmesh->nverts; i++)
	{
		int c=(quadmesh->quad_verts[i]->x-xstart)/dx;
		int r=(quadmesh->quad_verts[i]->y-ystart)/dy;

		int id=(r*(width)+c);

		/*  !!!! we may make use of more intelligent function here 11/25/2007 */

		quadmesh->quad_verts[i]->density=
			((double)(map[3*id]-smallest)/val_range);

		quadmesh->quad_verts[i]->density+=1;
	}
}



/*
*/
void set_vegflags_verts(double xstart, double xend, double ystart, double yend,
				  unsigned char *map, int width, int height)
{
	double xrang = xend-xstart;
	double yrang = yend-ystart;
	double dx=xrang/(width-1);
	double dy=yrang/(height-1);
	int i;

	for(i=0; i<quadmesh->nverts; i++)
	{
		int c=(quadmesh->quad_verts[i]->x-xstart)/dx;
		int r=(quadmesh->quad_verts[i]->y-ystart)/dy;

		if(r>=height) r=height-1;
		if(c>=width) c=width-1;

		if(r<0) r=0;
		if(c<0) c=0;

		int id=(r*(width)+c);

		if((map[3*id]+map[3*id+1]+map[3*id+2])> 200)
			quadmesh->quad_verts[i]->inveg = true;
		else
			quadmesh->quad_verts[i]->inveg = false;
	}
}