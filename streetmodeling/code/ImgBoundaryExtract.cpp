#include "stdafx.h"

#include "VFDataStructure.h"

#include "ImgBoundaryExtract.h"

#include "EvenlyStreamlines.h"

#include "sketchdesign.h"

extern int boundindex1, boundindex2;
extern bool ylargerthanx;
extern double mintenline_length;

extern TrajectoryList *sketchlines;
extern QuadMesh *quadmesh;

MapBoundaryList *mapboundarylist=NULL;
extern void get_linesegs_twopts_on_asketchline(ctr_point *p1, int cell1, 
										ctr_point *p2, int cell2,
										int &nlines, Trajectory *traj);

double minBoundaryLen=6.;

/*Implementation of the member functions of class ImgProcess*/

void ImgProcess::label_different_regions()
{
	int i, j;
	/*initialization*/
	for(i=0;i<512;i++)
	{
		for(j=0;j<512;j++)
			labels[i][j]=-2;
	}

	curMaxLabels = 200;
	parents=(int*)malloc(sizeof(int)*curMaxLabels);
	parents[0]=-1;
	label_count = 0;

	/*
	start labeling
	NOTE: here we assume that the input image is black/white color
	*/
	for(i=0; i<512; i++)
	{
		for(j=0; j<512; j++)
		{
			if(imgs[i][j][0]>=100) /*we meet a white pixel*/
			{
				if(i==0 && j==0) /*if it is the first pixel*/
				{
					labels[i][j]=label_count; /*label it*/
				}

				else if(i==0 && j>0)/*the first row of the pixels*/
				{
					/*check its previous pixel*/
					if(labels[i][j-1]>=0)
					{
						labels[i][j]=labels[i][j-1];
					}
					else
					{
						label_count++;
						labels[i][j]=label_count;

						if(label_count>=curMaxLabels)
						{
							parents=(int*)realloc(parents, sizeof(int)*(curMaxLabels+100));
							if(parents==NULL)
								exit(-1);   //memory reallocation failed
							curMaxLabels+=100;
						}
						parents[label_count]=-1;
					}
				}

				else if(j==0 && i>0)/*the first column*/
				{
					/*check the pixel in the previous level with the same x index*/
					if(labels[i-1][j]>=0)
					{
						labels[i][j]=labels[i-1][j];
					}
					else
					{
						label_count++;
						labels[i][j]=label_count;
						
						if(label_count>=curMaxLabels)
						{
							parents=(int*)realloc(parents, sizeof(int)*(curMaxLabels+100));
							if(parents==NULL)
								exit(-1);   //memeory reallocation failed
							curMaxLabels+=100;
						}
						parents[label_count]=-1;
					}
				}

				else
				{
					/*need to check its previous pixel and the pixel
					in the previous level with the same x index*/
					if(labels[i-1][j]<0 && labels[i][j-1]<0) /*if both of them are not labeled*/
					{
						label_count++;
						labels[i][j]=label_count;
						if(label_count>=curMaxLabels)
						{
							parents=(int*)realloc(parents, sizeof(int)*(curMaxLabels+100));
							if(parents==NULL)
								exit(-1);   //memeory reallocation failed
							curMaxLabels+=100;
						}

						parents[label_count]=-1;
					}

					else if(labels[i-1][j]<0&& labels[i][j-1]>=0)
						labels[i][j]=labels[i][j-1];
					else if(labels[i-1][j]>=0 && labels[i][j-1]<0)
						labels[i][j]=labels[i-1][j];

					else if(labels[i-1][j]>=0 && labels[i][j-1]>=0)
					{ /*if they are both labeled*/
						int label_previous_level=labels[i-1][j];
						int label_previous=labels[i][j-1];

						if(label_previous_level==label_previous)
						{
							labels[i][j]=label_previous;
						}
						else if(label_previous<label_previous_level)
						{
							/*if the label index of previous pixel is smaller than the label index
							of the pixel in previous level*/
							labels[i][j]=label_previous; /*we choose the smaller one*/
							/*set the parent of the labeling for the previous level pixel
							as the index of the previous pixel in the same level*/
							parents[label_previous_level]=label_previous; 
						}
						else
						{
							/*if the label index of previous pixel is smaller than the label index
							of the pixel in previous level*/
							labels[i][j]=label_previous_level; /*we choose the smaller one*/
							/*set the parent of the labeling for the previous level pixel
							as the index of the previous pixel in the same level*/
							parents[label_previous]=label_previous_level; 
						}
					}
				}
			}
		}
	}

	//parents[0]=-1;

    /*relabeling*/
	for(i=0; i<label_count+1; i++)
	{
		if(parents[i]<0)
			continue;

		j=i;
		int *temp_index=(int*)malloc(sizeof(int));
		int nelems = 0;
		temp_index[0]=j;
		nelems=1;

		while(1)
		{
			if(parents[j]<0)
				break;

			j=parents[j];

			temp_index=(int*)realloc(temp_index, sizeof(int)*(nelems+1));
			temp_index[nelems]=j;
			nelems++;
		}

		for(int k=0; k<nelems-1; k++)
			parents[temp_index[k]]=temp_index[nelems-1];

		free(temp_index);
	}

	/*  Testing codes  */
	//int pos;
	//for(i=label_count-1;i>=0;i--)
	//{
	//	if(parents[i]<0)
	//	{
	//		pos=i;
	//		break;
	//	}
	//}

	/*reset the labels of those white regions*/
	for(i=0; i<512; i++)
	{
		for(j=0; j<512; j++)
		{
			if(labels[i][j]<0)
				continue;

			if(parents[labels[i][j]]<0)
				continue;

			//if(i==510 && j==88)
			//{
			//	int test=0;
			//}

			labels[i][j]=parents[labels[i][j]];

			//if(labels[i][j]>pos)
			//{
			//	int test=0;
			//}
		}
	}
}

/*we perform a mean filer to the original image*/
void ImgProcess::mean_filter()
{
	int i,j;
	int blackcount, whitecount;
	//unsigned char *temp_img=(unsigned char*)malloc(sizeof(unsigned char)*3*512*512);

	//for(i=0;i<512;i++)
	//{
	//	for(j=0;j<512;j++)
	//	{
	//		temp_img[3*(i*512+j)]  =imgs[i][j][0];
	//		temp_img[3*(i*512+j)+1]=imgs[i][j][1];
	//		temp_img[3*(i*512+j)+2]=imgs[i][j][2];
	//	}
	//}

	for(i=1;i<511;i++)
	{
		for(j=1;j<511;j++)
		{
			blackcount=whitecount=0;
			/*count the numbers of black pixels and white pixels in the 3x3 region with current
			pixel as the center*/
			for(int k=-1; k<=1; k++)
			{
				for(int l=-1; l<=1; l++)
				{
					if(/*temp_img[3*((i+k)*512+(j+l))]*/
						origin_img[i+k][j+l][0]<100)
						blackcount++;
					else
						whitecount++;
				}
			}
			if(blackcount>whitecount)
			{
				imgs[i][j][0]=
				imgs[i][j][1]=
				imgs[i][j][2]=0;
			}
			else
			{
				imgs[i][j][0]=
				imgs[i][j][1]=
				imgs[i][j][2]=254;
			}
		}
	}

	//free(temp_img);
}


/*we perform a mean filer to the original image*/
void ImgProcess::mean_filter(int windowsize)
{
	int i,j;
	int blackcount, whitecount;
	unsigned char *temp_img=(unsigned char*)malloc(sizeof(unsigned char)*3*512*512);
	int mid_size=(int)(windowsize/2.);

	for(i=0;i<512;i++)
	{
		for(j=0;j<512;j++)
		{
			temp_img[3*(i*512+j)]  =imgs[i][j][0];
			temp_img[3*(i*512+j)+1]=imgs[i][j][1];
			temp_img[3*(i*512+j)+2]=imgs[i][j][2];
		}
	}

	for(i=mid_size;i<512-mid_size;i++)
	{
		for(j=mid_size;j<512-mid_size;j++)
		{
			blackcount=whitecount=0;
			/*count the numbers of black pixels and white pixels in the 3x3 region with current
			pixel as the center*/
			for(int k=-mid_size; k<=mid_size; k++)
			{
				for(int l=-mid_size; l<=mid_size; l++)
				{
					if(/*temp_img[3*((i+k)*512+(j+l))]<100*/
						origin_img[i+k][j+l][0]<100)
						blackcount++;
					else
						whitecount++;
				}
			}
			if(blackcount>whitecount)
			{
				imgs[i][j][0]=
				imgs[i][j][1]=
				imgs[i][j][2]=0;
			}
			else
			{
				imgs[i][j][0]=
				imgs[i][j][1]=
				imgs[i][j][2]=254;
			}
		}
	}

	free(temp_img);
}

void ImgProcess::mean_filter_usesame(int windowsize)
{
	int i,j;
	int blackcount, whitecount;
	int mid_size=(int)(windowsize/2.);

	for(i=mid_size;i<512-mid_size;i++)
	{
		for(j=mid_size;j<512-mid_size;j++)
		{
			blackcount=whitecount=0;
			/*count the numbers of black pixels and white pixels in the 3x3 region with current
			pixel as the center*/
			for(int k=-mid_size; k<=mid_size; k++)
			{
				for(int l=-mid_size; l<=mid_size; l++)
				{
					if(imgs[i+k][j+l][0]<100)
						blackcount++;
					else
						whitecount++;
				}
			}
			if(blackcount>whitecount)
			{
				imgs[i][j][0]=
				imgs[i][j][1]=
				imgs[i][j][2]=0;
			}
			else
			{
				imgs[i][j][0]=
				imgs[i][j][1]=
				imgs[i][j][2]=254;
			}
		}
	}
}



/*after labeling the white regions, we now find out the boundaries of 
those regions, respectively*/
void ImgProcess::label_boundary_pixels()
{
	int i, j, k, l;

	for(i=1; i<511; i++)
	{
		for(j=1; j<511; j++)
		{
			if(imgs[i][j][0]<100)
				continue;

			/*for a white pixel, propagate to its four direct neighbors*/
			labels[i][j-1]=labels[i][j];
			labels[i][j+1]=labels[i][j];
			labels[i-1][j]=labels[i][j];
			labels[i+1][j]=labels[i][j];
			//for(k=-1;k<=1;k++)
			//{
			//	for(l=-1;l<=1;l++)
			//		labels[i+k][j+l]=labels[i][j];
			//}
		}
	}
}


extern bool is_repeated_elem(int *, int, int);
/*
After labeling the regions, extract how many pixels inside each regions,
and how many regions are obtained
*/
void ImgProcess::extract_regions()
{
	int i,j;

	if(npixels_regs != NULL)
		free(npixels_regs);

	npixels_regs=(int*)malloc(sizeof(int)*(label_count+1));
	int *temp=(int*)malloc(sizeof(int)*(label_count+1));
	int nelems_temp=0;

	nregs=0;

	for(i=0;i<label_count+1; i++)
	{
		npixels_regs[i]=0;
		temp[i]=-1;
	}

	for(i=0;i<512;i++)
	{
		for(j=0;j<512;j++)
		{
			if(labels[i][j]<0)
				continue;

			if(!is_repeated_elem(temp, labels[i][j], nelems_temp))
			{
				temp[nelems_temp]=labels[i][j];
				nelems_temp++;
				nregs++;
			}

			npixels_regs[labels[i][j]]++;
		}
	}

	free(temp);

	if(region_label!=NULL)
	{
		free(region_label);
		region_label=NULL;
	}

	region_label=(int*)malloc(sizeof(int)*nregs);
	int temp_index=0;
	for(i=0;i<label_count+1; i++)
	{
		if(npixels_regs[i]>0){
			region_label[temp_index]=i;
			temp_index++;
		}
	}
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Implementation of the member functions of class ImgBoundary*/

ImgBoundary::~ImgBoundary()
{
	if(output!=NULL)
		free(output);

	if(boundaries!=NULL)
		delete [] boundaries;

	if(allboundpixels!=NULL)
		delete allboundpixels;

	if(streetmap_background!=NULL)
		free(streetmap_background);

	if(boundarylist!=NULL)
		delete [] boundarylist;

}


void ImgBoundary::get_bound_pixels(ImgProcess *imgprocess)
{
	/*allocate memory according to the marked pixel*/

	boundaries=new PixelBoundary[imgprocess->nregs];
	nboundaries = imgprocess->nregs;

	int i,j;

	for(i=0;i<512;i++)
	{
		for(j=0;j<512;j++)
		{
			if(imgprocess->imgs[i][j][0]<120 && imgprocess->labels[i][j]>=0 
				/*&& imgprocess->labels[i][j]<imgprocess->label_count*/	)
			{
				//boundaries[imgprocess->labels[i][j]].addNew(j,i);

				int k=find_boundary_label(imgprocess, i, j);
				boundaries[k].addNew(j,i);
			}
		}
	}
}


void ImgBoundary::get_allboundary_pixels(ImgProcess *imgprocess)
{
	allboundpixels=new PixelBoundary(5000);

	//pos_in_boundpixels

	int i,j;
	for(i=0;i<512;i++)
	{
		for(j=0;j<512; j++)
		{
			index_in_boundarypixel_list[i][j]=-1;
			if(imgprocess->imgs[i][j][0]<100 && imgprocess->labels[i][j]>=0 
				/*&& !is_boundarypixel(j, i)*/)
			{
				index_in_boundarypixel_list[i][j]=allboundpixels->nboundpts;
				allboundpixels->addNew(j,i); /*i is y, j is x*/
			}
		}
	}
}


void ImgBoundary::naive_find_sortedboundaries(ImgProcess *imgprocess)
{
	/*initialize the visited state of the boundary pixels*/
	int i;
	for(i=0;i<allboundpixels->nboundpts;i++)
	{
		allboundpixels->boundpts[i].visited=false;
	}

	//curMaxNum = imgprocess->nregs;
	curMaxNum=20;
	boundarylist=new PixelBoundary[curMaxNum];

	nboundaries=0;

	/*start searching all the boundaries*/
	for(i=0;i<allboundpixels->nboundpts;i++)
	{
		if(allboundpixels->boundpts[i].visited)
			continue;
		allboundpixels->boundpts[i].visited=true;

		/*start searching from this pixel*/
		if(nboundaries>=curMaxNum)
		{
			/*extend the space*/
			PixelBoundary *temp=boundarylist;
			boundarylist=new PixelBoundary[curMaxNum+10];
			if(boundarylist==NULL)
				exit(-1);//memory reallocation failed

			/*copy previous elements*/
			for(int j=0;j<curMaxNum;j++)
			{
				boundarylist[j].copy_boundary(&temp[j]);
			}

			curMaxNum+=10;

			delete [] temp;
		}

		//boundarylist[nboundaries]=new PixelBoundary();
		//boundarylist
		boundarylist[nboundaries].nboundpts=0;
		naive_search_one_boundary(imgprocess, allboundpixels->boundpts[i].x, 
			allboundpixels->boundpts[i].y, &boundarylist[nboundaries]);
		nboundaries++;
	}
}


/*
judge whether a pixel is on the boundary of the map
*/
bool is_boundarypixel(int x, int y)
{
	if(!ylargerthanx) /*width>height*/
	{
		if(y==boundindex1 || y==boundindex2
			||y==boundindex1-1 || y==boundindex1+1 /*|| y==boundindex1+2 || y==boundindex1-2*/
			||y==boundindex2-1 || y==boundindex2+1 || y==boundindex2+2 || y==boundindex2-2)
			return true;
	}
	else
	{
		if(x==boundindex1 || x==boundindex2			
			||x==boundindex1-1 || x==boundindex1+1 || x==boundindex1+2 || x==boundindex1-2
			||x==boundindex2-1 || x==boundindex2+1 || x==boundindex2+2 || x==boundindex2-2)
			return true;
	}
	return false;
}

/*
searching one boundary starting from pixel(starti,startj)
For this naive searching algorithm, we assume that the boundaries are
physically separated and will not have self-intersection!!
11/11/2007
*/	
void ImgBoundary::naive_search_one_boundary(ImgProcess *imgprocess, int starti, int startj, 
									  PixelBoundary *curboundary)
{
	int cur_x, cur_y;
	/*first, add the first pixel*/
	curboundary->addNew(starti, startj);

	cur_x=starti;
	cur_y=startj;

	while(1/*neithor boundary nor form a loop*/)
	{
		/*if it is already a boundary point, break*/
		//if(cur_x==0||cur_x==511||cur_y==0||cur_y==511)
		//	break;

		if(cur_x==curboundary->boundpts[0].x && cur_y==curboundary->boundpts[0].y
			&& curboundary->nboundpts>10)
			break;

		/*consider the eight neighors of current point*/

		if(cur_x<511&&index_in_boundarypixel_list[cur_y][cur_x+1]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x+1]].visited
			/*&& !is_boundarypixel(cur_x+1, cur_y)*/) /*cur_x+1, cur_y*/ //right
		{
			cur_x=cur_x+1;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}
		else if(cur_y<511&&index_in_boundarypixel_list[cur_y+1][cur_x]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y+1][cur_x]].visited
			/*&& !is_boundarypixel(cur_x, cur_y+1)*/)/*cur_x, cur_y+1*/ //up
		{
			cur_y=cur_y+1;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}
		else if(cur_x<511 && cur_y<511 && index_in_boundarypixel_list[cur_y+1][cur_x+1]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y+1][cur_x+1]].visited
			/*&& !is_boundarypixel(cur_x+1, cur_y+1)*/)/*cur_x+1, cur_y+1*/ //uppperright
		{
			cur_x=cur_x+1;
			cur_y=cur_y+1;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}
		else if(cur_y>0&&index_in_boundarypixel_list[cur_y-1][cur_x]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y-1][cur_x]].visited
			  /*&& !is_boundarypixel(cur_x, cur_y-1)*/)/*cur_x, cur_y-1*/ //down
		{
			cur_y=cur_y-1;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}
		else if( cur_x<511 && cur_y>0&&index_in_boundarypixel_list[cur_y-1][cur_x+1]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y-1][cur_x+1]].visited
			/*&& !is_boundarypixel(cur_x+1, cur_y-1)*/)/*cur_x+1, cur_y-1*/ //lowerright
		{
			cur_x=cur_x+1;
			cur_y--;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}
		else if(cur_x>0 && cur_y>0&&index_in_boundarypixel_list[cur_y-1][cur_x-1]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y-1][cur_x-1]].visited
			  /*&& !is_boundarypixel(cur_x-1, cur_y-1)*/)/*cur_x-1, cur_y-1*/ //lower left
		{
			cur_x=cur_x-1;
			cur_y--;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}
		else if(cur_x>0&&index_in_boundarypixel_list[cur_y][cur_x-1]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x-1]].visited
			 /* && !is_boundarypixel(cur_x-1, cur_y)*/)/*cur_x-1, cur_y*/ //left
		{
			cur_x=cur_x-1;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}
		else if(cur_x>0 && cur_y<511&&index_in_boundarypixel_list[cur_y+1][cur_x-1]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y+1][cur_x-1]].visited
			 /* && !is_boundarypixel(cur_x-1, cur_y+1)*/)/*cur_x-1, cur_y+1*/ //upper left
		{
			cur_x=cur_x-1;
			cur_y=cur_y+1;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}
		else
			break; // this is an end pixel of the boundary
	}

	/*reverse the order of the obtained list of pixels*/
	curboundary->inverse_list();

	/*start searching from the last pixel again!*/
	cur_x=curboundary->boundpts[curboundary->nboundpts-1].x;
	cur_y=curboundary->boundpts[curboundary->nboundpts-1].y;
	while(1/*neithor boundary nor form a loop*/)
	{

		if(cur_x==curboundary->boundpts[0].x && cur_y==curboundary->boundpts[0].y
			&& curboundary->nboundpts>10)
			break;

		/*consider the eight neighors of current point*/

		if(cur_x>0&&index_in_boundarypixel_list[cur_y][cur_x-1]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x-1]].visited
			  /*&& !is_boundarypixel(cur_x-1, cur_y)*/)/*cur_x-1, cur_y*/ //left
		{
			cur_x=cur_x-1;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}
		else if(cur_y>0&&index_in_boundarypixel_list[cur_y-1][cur_x]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y-1][cur_x]].visited
			  /*&& !is_boundarypixel(cur_x, cur_y-1)*/)/*cur_x, cur_y-1*/ //down
		{
			cur_y=cur_y-1;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}
		else if(cur_x>0 && cur_y<511&&index_in_boundarypixel_list[cur_y+1][cur_x-1]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y+1][cur_x-1]].visited
			 /*&& !is_boundarypixel(cur_x-1, cur_y+1)*/)/*cur_x-1, cur_y+1*/ //upper left
		{
			cur_x=cur_x-1;
			cur_y=cur_y+1;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}
		else if(cur_x>0 && cur_y>0&&index_in_boundarypixel_list[cur_y-1][cur_x-1]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y-1][cur_x-1]].visited
			  /*&& !is_boundarypixel(cur_x-1, cur_y-1)*/)/*cur_x-1, cur_y-1*/ //lower left
		{
			cur_x=cur_x-1;
			cur_y--;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}

		else if(cur_y<511&&index_in_boundarypixel_list[cur_y+1][cur_x]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y+1][cur_x]].visited
			 /*&& !is_boundarypixel(cur_x, cur_y+1)*/)/*cur_x, cur_y+1*/ //up
		{
			cur_y=cur_y+1;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}
		else if(cur_x<511&&index_in_boundarypixel_list[cur_y][cur_x+1]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x+1]].visited
			 /*&& !is_boundarypixel(cur_x+1, cur_y)*/) /*cur_x+1, cur_y*/ //right
		{
			cur_x=cur_x+1;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}
		else if(cur_x<511 && cur_y<511&& index_in_boundarypixel_list[cur_y+1][cur_x+1]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y+1][cur_x+1]].visited
			/*&& !is_boundarypixel(cur_x+1, cur_y+1)*/)/*cur_x+1, cur_y+1*/ //upper right
		{
			cur_x=cur_x+1;
			cur_y=cur_y+1;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}
		else if(cur_x<511 && cur_y>0&&index_in_boundarypixel_list[cur_y-1][cur_x+1]>=0
			&& !allboundpixels->boundpts[index_in_boundarypixel_list[cur_y-1][cur_x+1]].visited
			  /*&& !is_boundarypixel(cur_x+1, cur_y-1)*/)/*cur_x+1, cur_y-1*/ //lower right
		{
			cur_x=cur_x+1;
			cur_y--;
			curboundary->addNew(cur_x, cur_y);
			allboundpixels->boundpts[index_in_boundarypixel_list[cur_y][cur_x]].visited=true;
		}

		else
			break; // this is an end pixel of the boundary
	}

}


int ImgBoundary::find_boundary_label(ImgProcess *imgprocess, int i, int j)
{
	int k;
	for(k=0; k<imgprocess->nregs; k++)
	{
		if(imgprocess->labels[i][j]==imgprocess->region_label[k])
			return k;
	}
}

/*sort all the boundaries*/
void ImgBoundary::sort_boundaries()
{
	int i;
	for(i=0;i<nboundaries;i++)
	{
		sort_aboundary(i);
	}
}

	
void ImgBoundary::sort_aboundary(int boundid)
{
}

/*
We render the result to an array
*/
extern void getcolor_scc_frac(int num, int frac, float rgb[3]);

void ImgBoundary::render_result(ImgProcess *imgprocess)
{
	if(output==NULL)
		output=(unsigned char*)malloc(sizeof(unsigned char)*3*512*512);

	//float rgb[3];
	int i,j;
	//for(i=0;i<512;i++)
	//{
	//	for(j=0;j<512;j++)
	//	{
	//		if(imgprocess->imgs[i][j][0]<120)
	//		{
	//			output[3*(i*512+j)]=
	//			output[3*(i*512+j)+1]=
	//			output[3*(i*512+j)+2]=imgprocess->imgs[i][j][0];
	//		}
	//		else
	//		{
	//			//getcolor_scc_frac(imgprocess->labels[i][j], nboundaries, rgb);
	//			//output[3*(i*512+j)]=(unsigned char)(255*rgb[0]);
	//			//output[3*(i*512+j)+1]=(unsigned char)(255*rgb[1]);
	//			//output[3*(i*512+j)+2]=(unsigned char)(255*rgb[2]);
	//			output[3*(i*512+j)]=
	//			output[3*(i*512+j)+1]=
	//			output[3*(i*512+j)+2]=220;
	//		}
	//	}
	//}

	///*render the boundaries*/
	////for(i=0;i<nboundaries;i++)
	////{
	////	for(j=0;j<boundaries[i].nboundpts;j++)
	////	{
	////		int x=boundaries[i].boundpts[j].i;
	////		int y=boundaries[i].boundpts[j].j;
	////		output[3*(y*512+x)]=
	////		output[3*(y*512+x)+1]=
	////		output[3*(y*512+x)+2]=255;
	////	}
	////}


	///*render different boundaries with different colors*/
	//for(i=0;i<allboundpixels->nboundpts;i++)
	//{
	//	int x=allboundpixels->boundpts[i].x;
	//	int y=allboundpixels->boundpts[i].y;
	//	output[3*(y*512+x)]=
	//	output[3*(y*512+x)+1]=
	//	output[3*(y*512+x)+2]=255;
	//}

	//for(i=0;i<nboundaries; i++)
	//{
	//	PixelBoundary *oneb=&boundarylist[i];
	//	getcolor_scc_frac(i, nboundaries, rgb);
	//	for(j=0;j<oneb->nboundpts;j++)
	//	{
	//		int x=oneb->boundpts[j].x;
	//		int y=oneb->boundpts[j].y;
	//		output[3*(y*512+x)]=(unsigned char)(255*rgb[0]);
	//		output[3*(y*512+x)+1]=(unsigned char)(255*rgb[1]);
	//		output[3*(y*512+x)+2]=(unsigned char)(255*rgb[2]);
	//	}
	//}
	
	
	/*render the street map background*/

	if(streetmap_background==NULL)
		streetmap_background=(unsigned char*)malloc(sizeof(unsigned char)*3*512*512);

	for(i=0;i<512;i++)
	{
		for(j=0;j<512;j++)
		{
			//if(imgprocess->imgs[i][j][0]<120)
			if(imgprocess->origin_img[i][j][0]<120)
			{
				//streetmap_background[3*(i*512+j)]=
				//streetmap_background[3*(i*512+j)+1]=
				//streetmap_background[3*(i*512+j)+2]=204;

				streetmap_background[3*(i*512+j)]=235;
				streetmap_background[3*(i*512+j)+1]=230;
				streetmap_background[3*(i*512+j)+2]=220;
			}
			else
			{	//153 179 204
				streetmap_background[3*(i*512+j)]=153;
				streetmap_background[3*(i*512+j)+1]=179;
				streetmap_background[3*(i*512+j)+2]=204;
			}
		}
	}


	for(i=0;i<512;i++)
	{
		for(j=0;j<512;j++)
		{
			if(imgprocess->imgs[i][j][0]<120)
			//if(imgprocess->origin_img[i][j][0]<120)
			{
				//output[3*(i*512+j)]=
				//output[3*(i*512+j)+1]=
				//output[3*(i*512+j)+2]=204;
				output[3*(i*512+j)]=235;
				output[3*(i*512+j)+1]=230;
				output[3*(i*512+j)+2]=220;
			}
			else
			{	//153 179 204
				output[3*(i*512+j)]=153;
				output[3*(i*512+j)+1]=179;
				output[3*(i*512+j)+2]=204;
			}
		}
	}
}





/************************************************************************************/
/***********************************************************************************/
/*  Implementation of the member function of class MapBoundaryList*/
void MapBoundaryList::downsamp_imgboundaries(ImgBoundary *imgboundary, int rate)
{
	int i,j;
	nmapboundaries = 0;
	PixelBoundary *curboundary;

	/*   initialization   */
	QuadCell *face;
	for(i=0;i<quadmesh->nfaces;i++)
	{
		face=quadmesh->quadcells[i];
		if(face->mapbounds!=NULL)
		{
			free(face->mapbounds);
			face->mapbounds=NULL;
			face->nmapbounds=0;
		}
	}

	for(i=0;i<imgboundary->nboundaries;i++)
	{
		curboundary=&imgboundary->boundarylist[i];

		if(nmapboundaries>=curMaxNum)
		{
			if(!extend(10))
				exit(-1);
		}

		int boundpts_counter=0;
		icVector2 tempDist;
		double boundLen=0;
		for(j=0;j<curboundary->nboundpts;j++)
		{
			if(j%rate==0)
			{
				boundpts_counter++;

				/*  accumulate the distance of the boundary  */
				if(j<curboundary->nboundpts-1)
				{
					tempDist.entry[0]=curboundary->boundpts[j+1].x-
						curboundary->boundpts[j].x;
					tempDist.entry[1]=curboundary->boundpts[j+1].y-
						curboundary->boundpts[j].y;
					boundLen+=length(tempDist);
				}
			}
		}
		if(boundpts_counter<=3)
			continue;

		if(boundLen <= minBoundaryLen)
			continue;

		mapboundarylist[nmapboundaries].nelems=0;

		for(j=0;j<curboundary->nboundpts;j++)
		{

			/*   add to the boundary list of each cell  */
			face=quadmesh->quadcells[get_cellID_givencoords(curboundary->boundpts[j].x*dx, 
				curboundary->boundpts[j].y*dy)];
			if(face->mapbounds==NULL)
			{
				face->mapbounds=(int*)malloc(sizeof(int));
				face->mapbounds[0]=nmapboundaries;
				face->nmapbounds=1;
			}
			else
			{
				if(!is_repeated_elem(face->mapbounds, nmapboundaries, face->nmapbounds))
				{
					face->mapbounds=(int*)realloc(face->mapbounds, sizeof(int)*(face->nmapbounds+1));
					face->mapbounds[face->nmapbounds]=nmapboundaries;
					face->nmapbounds++;
				}
			}

			/*   obtain a new sample on the boundary here   */
			if(j%rate==0)
			{
				PtOnBoundary *newpt=(PtOnBoundary*)malloc(sizeof(PtOnBoundary));
				newpt->x=curboundary->boundpts[j].x*dx;
				newpt->y=curboundary->boundpts[j].y*dy;
				mapboundarylist[nmapboundaries].addNew(newpt);
			}

		}

		/*	Recall the boundaries that are closed 12/26/2007 */
		icVector2 dist;
		dist.entry[0]=curboundary->boundpts[0].x*dx-
			curboundary->boundpts[curboundary->nboundpts-1].x*dx;
		dist.entry[1]=curboundary->boundpts[0].y*dy-
			curboundary->boundpts[curboundary->nboundpts-1].y*dy;

		double dist_len=length(dist);

		if(dist_len<5.e-3)    /*  it is closed!  */
		{
			PtOnBoundary *newpt=(PtOnBoundary*)malloc(sizeof(PtOnBoundary));
			newpt->x=curboundary->boundpts[0].x*dx;
			newpt->y=curboundary->boundpts[0].y*dy;
			mapboundarylist[nmapboundaries].addNew(newpt);
		}

	//	if(mapboundarylist[nmapboundaries].nelems>3)
	//		nmapboundaries++;
	//	else
	//	{
	//		/*  we need to release the previous assign memory*/
	//		for(j=0;j<mapboundarylist[nmapboundaries].nelems;j++)
	//		{
	//			if(mapboundarylist[nmapboundaries].pts[j]!=NULL)
	//			{
	//				free(mapboundarylist[nmapboundaries].pts[j]);
	//				mapboundarylist[nmapboundaries].pts[j]=NULL;
	//			}
	//		}
	//	}
		nmapboundaries++;
	}
}

	
/*
We index all the boundary points for future editing process
*/
void MapBoundaryList::index_allboundarypts()
{
	int g_id=0;
	int i,j;
	MapBoundary *curboundary;

	for(i=0;i<nmapboundaries;i++)
	{
		curboundary=&mapboundarylist[i];
		for(j=0;j<curboundary->nelems;j++)
		{
			curboundary->pts[j]->index=g_id;
			g_id++;
		}
	}
}

//#include <gl/glut.h> 

void MapBoundaryList::display_boundaries(GLenum mode)
{
	if(mapboundarylist==NULL) return;

	float rgb[3]={0.};

	int i,j;
	MapBoundary *curboundary;
	for(i=0;i<nmapboundaries;i++)
	{
		getcolor_scc_frac(i, nmapboundaries, rgb);
		glColor3fv(rgb);
		curboundary=&mapboundarylist[i];
		glBegin(GL_LINE_STRIP);
		for(j=0;j<curboundary->nelems;j++)
			glVertex2f(curboundary->pts[j]->x, curboundary->pts[j]->y);
		glEnd();
	}

	//glPointSize(5.);
	//glColor3f(1,0,0);

	//for(i=0;i<nmapboundaries;i++)
	//{
	//	curboundary=&mapboundarylist[i];
	//	glBegin(GL_POINTS);
	//	for(j=0;j<curboundary->nelems;j++)
	//	{
	//		if(mode == GL_SELECT)
	//			glLoadName(NAMEOFSHAPECONTROL+curboundary->pts[j]->index);
	//		glVertex2f(curboundary->pts[j]->x, curboundary->pts[j]->y);
	//	}
	//	glEnd();
	//}
}


/*   we need to convert the boundaries of the map into trajectories for region
     segmentation of the local editing  12/20/2007
*/


void convert_one_bound_to_one_traj(int id)
{
	MapBoundary *abound=&mapboundarylist->mapboundarylist[id];
	Trajectory *traj=new Trajectory(id, abound->nelems);
	//Trajectory *traj=new Trajectory(id);

	int i;
	int lineid=0;
	int cell1, cell2;
	icVector2 line_dir;
	for(i=0;i<abound->nelems-1;i++)
	{
		if(abound->pts[i]->x<quadmesh->xstart+1.e-8 || abound->pts[i]->x>quadmesh->xend-1.e-8
			||abound->pts[i]->y<quadmesh->ystart+1.e-8 || abound->pts[i]->x>quadmesh->yend-1.e-8
			||abound->pts[i+1]->x<quadmesh->xstart+1.e-8 || abound->pts[i+1]->x>quadmesh->xend-1.e-8
			||abound->pts[i+1]->y<quadmesh->ystart+1.e-8 || abound->pts[i+1]->x>quadmesh->yend-1.e-8)
			break;

		cell1=get_cellID_givencoords(abound->pts[i]->x, abound->pts[i]->y);
		cell2=get_cellID_givencoords(abound->pts[i+1]->x, abound->pts[i+1]->y);

		if(cell1<0||cell1>=quadmesh->nfaces
			||cell2<0||cell2>=quadmesh->nfaces)
			break;
		
		if(traj->nlinesegs>=traj->curMaxNumLinesegs)
		{
			traj->extend_line_segments(50);
		}

		if(cell1==cell2)
		{
			traj->linesegs[lineid].gstart[0]=abound->pts[i]->x;
			traj->linesegs[lineid].gstart[1]=abound->pts[i]->y;
			traj->linesegs[lineid].gend[0]=abound->pts[i+1]->x;
			traj->linesegs[lineid].gend[1]=abound->pts[i+1]->y;
			traj->linesegs[lineid].Triangle_ID=cell1;
			line_dir.entry[0]=abound->pts[i+1]->x-abound->pts[i]->x;
			line_dir.entry[1]=abound->pts[i+1]->y-abound->pts[i]->y;
			traj->linesegs[lineid].length=length(line_dir);
			lineid++;
			traj->nlinesegs=lineid;
		}
		else
		{
			ctr_point p1,p2;
			p1.x=abound->pts[i]->x;
			p1.y=abound->pts[i]->y;
			p1.cellid=cell1;
			p2.x=abound->pts[i+1]->x;
			p2.y=abound->pts[i+1]->y;
			p2.cellid=cell2;
			get_linesegs_twopts_on_asketchline(&p1, cell1,
				&p2, cell2, lineid, traj);
			//double p1[2],p2[2];
			//p1[0]=abound->pts[i]->x;
			//p1[1]=abound->pts[i]->y;
			//p2[0]=abound->pts[i+1]->x;
			//p2[1]=abound->pts[i+1]->y;

			//get_linesegs_anytwopts(p1, cell1, 
			//			p2, cell2, traj, 0, 10);
			lineid=traj->nlinesegs;
		}
	}

	/*  we need to judge whether it is a closed tensor line  */
	line_dir.entry[0]=traj->linesegs[0].gstart[0]-
		traj->linesegs[traj->nlinesegs-1].gend[0];
	line_dir.entry[1]=traj->linesegs[0].gstart[1]-
		traj->linesegs[traj->nlinesegs-1].gend[1];

	traj->closed=false;

	if(length(line_dir)<1.e-2)
	{
		if(traj->nlinesegs>=traj->curMaxNumLinesegs)
		{
			traj->extend_line_segments(1);
		}
		traj->linesegs[traj->nlinesegs].gstart[0]=traj->linesegs[traj->nlinesegs-1].gend[0];
		traj->linesegs[traj->nlinesegs].gstart[1]=traj->linesegs[traj->nlinesegs-1].gend[1];
		traj->linesegs[traj->nlinesegs].gend[0]=traj->linesegs[0].gstart[0];
		traj->linesegs[traj->nlinesegs].gend[1]=traj->linesegs[0].gstart[1];
		traj->linesegs[traj->nlinesegs].Triangle_ID=traj->linesegs[traj->nlinesegs-1].Triangle_ID;
		line_dir.entry[0]=traj->linesegs[traj->nlinesegs].gstart[0]-
			traj->linesegs[traj->nlinesegs].gend[0];
		line_dir.entry[1]=traj->linesegs[traj->nlinesegs].gstart[1]-
			traj->linesegs[traj->nlinesegs].gend[1];
		traj->linesegs[traj->nlinesegs].length=length(line_dir);

		//if(traj->get_length()>mintenline_length)
		traj->nlinesegs++;

		traj->closed=true;
	}

	traj->traj_len=traj->get_length();
	traj->is_mapboundary=true;

	/*  we don't consider the tiny boundaries now 12/27/2007  */
	if(traj->traj_len<4*quadmesh->xinterval)
	{
		delete traj;
		return;
	}

	sketchlines->append(traj);
}

void convert_mapbounds_to_trajs()
{
	int i;
	for(i=0;i<mapboundarylist->nmapboundaries;i++)
	{
		convert_one_bound_to_one_traj(i);
	}
	
}
