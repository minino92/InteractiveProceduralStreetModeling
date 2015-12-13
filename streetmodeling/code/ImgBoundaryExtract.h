
/*
ImgBoundaryExtract.h
*/
#include <gl/glut.h> 

typedef struct PtOnBoundary {
  double x;
  double y;
  int index;
}PtOnBoundary;

typedef struct Pixel{
	int x, y;
	//int listindex;
	bool visited;
}Pixel;

typedef struct OneInterval{
	//int starti, startj, endi, endj;
	Pixel start, end;
	PtOnBoundary start_coord, end_coord;
}OneInterval;


/*
since a region may have several consecutive intervals at one scan line
due to the complex shape, we define this class to record the intervals 
that belong to one region.
The operations include split one previous interval into multiple ones
merge multiple intervals into one.
All the intervals in one scan line are ordered according to their 
x direciton index.
*/
class OneScanLine
{
public:
	OneInterval **intervals;
	int nintervals;
	int curMaxNum;

	OneScanLine(int initsize=0)
	{
		if(initsize==0)
		{
			intervals=NULL;
			nintervals=curMaxNum=0;
		}

		intervals=(OneInterval**)malloc(sizeof(OneInterval*)*initsize);
		if(intervals == NULL)
			exit(-1);  //fail

		curMaxNum=initsize;
		nintervals=0;
		int i;
		for(i=0; i<curMaxNum; i++)
			intervals[i]=NULL;
	}

	~OneScanLine()
	{
		if(intervals == NULL)
			return;

		int i;
		for(i=0; i<curMaxNum; i++)
		{
			if(intervals[i]!=NULL)
			{
				free(intervals[i]);
				intervals[i]=NULL;
			}
		}
		free(intervals);
	}
	
	//List operations
	bool isfull()
	{
		if(nintervals >= curMaxNum) return true;
		return false;
	}

	bool extend(int step=2)
	{
		OneInterval **temp = intervals;
		intervals=(OneInterval**)malloc(sizeof(OneInterval*)*(curMaxNum+step));
		if(intervals == NULL)
			return false;

		int i;
		for(i=0; i<curMaxNum; i++)
			intervals[i]=temp[i];
		for(i=curMaxNum; i<curMaxNum+step; i++)
			intervals[i]=NULL;
		curMaxNum += step;
		free(temp);
		return true;

	}

	void addNew(OneInterval *interval)
	{
		if(nintervals>=curMaxNum)
		{
			if(!extend(2))
				exit(-1);
		}

		intervals[nintervals] = interval;
		nintervals++;
	}
};

/*
Recall that on one scan line, we could have multiple intervals.
Therefore, the operations include split one previous interval into multiple ones
merge multiple intervals into one.
All the intervals in one scan line are ordered according to their 
x direciton index.
*/
//class PixelBoundary
//{
//public:
//	OneScanLine *scanlines;
//	int nscanlines;
//	int curMaxNum;
//
//	Pixel *sortedboundpts;
//
//	PtOnBoundary *boundarypts;
//	int nboundarypts;
//
//	double dx, dy;
//
//	bool activate;
//
//	//constructor
//	PixelBoundary(int initsize=512) /*we now assume the image size is 512x512 for our applications*/
//	{
//		if(initsize == 0)
//		{
//			scanlines = NULL;
//			nscanlines=curMaxNum=0;
//			return;
//		}
//
//		scanlines=new OneScanLine[initsize];
//		nscanlines = 0;
//		curMaxNum=initsize;
//
//		sortedboundpts = NULL;
//		boundarypts = NULL;
//		nboundarypts = 0;
//		activate = true;
//	}
//
//	~PixelBoundary()
//	{
//		if(scanlines == NULL)
//			return;
//		delete [] scanlines;
//
//		if(sortedboundpts != NULL)
//			free(sortedboundpts);
//		if(boundarypts != NULL)
//			free(boundarypts);
//	}
//
//	/*merge with other boundary
//	merge according to the scan lines
//	*/
//	void merge_with_other_boundary(PixelBoundary *bound)
//	{
//		int i, j;
//		for(i=0; i<bound->nscanlines; i++)
//		{
//			if(bound->scanlines[i].intervals==NULL || bound->scanlines[i].nintervals==0)
//				continue;
//
//			/*extend memory if not enough space is available*/
//			if(this->scanlines[i].intervals==NULL)
//			{
//				this->scanlines[i].intervals=(OneInterval**)malloc(sizeof(OneInterval*)
//					*bound->scanlines[i].nintervals);
//				this->scanlines[i].curMaxNum=this->scanlines[i].nintervals
//					=bound->scanlines[i].nintervals;
//			}
//
//			if(this->scanlines[i].nintervals+bound->scanlines[i].nintervals
//				>=curMaxNum)
//			{
//				if(!this->scanlines[i].extend(bound->scanlines[i].nintervals))
//					exit(-1);
//			}
//
//			/*copy the content from other boundary to current boundary*/
//			for(j=0; j<bound->scanlines[i].nintervals; j++)
//			{
//				this->scanlines[i].intervals[this->scanlines[i].nintervals]=
//					bound->scanlines[i].intervals[j];
//				this->scanlines[i].nintervals++;
//			}
//		}
//
//		bound->activate=false;
//	}
//
//
//	//other important functions
//	//bool isBelongToThisBoundary(int nstarti, int nendi, int nextj)
//	//{
//	//	if(nextj!=intervals[nintervals-1]->start.j+1)
//	//		return false;
//	//	if((nstarti>intervals[nintervals-1]->end.i)
//	//		||(nendi<intervals[nintervals-1]->start.i))
//	//		return false;
//	//	return true;
//	//}
//
//	///*sort the boundary pixels to form a closed loop
//	//We need to reimplement this 
//	//*/
//	//void sort_boundarypixels()
//	//{
//	//	if(sortedboundpts != NULL)
//	//		free(sortedboundpts);
//
//	//	sortedboundpts = (Pixel*)malloc(sizeof(Pixel)*2*nintervals);
//	//	int i, index = 0;
//	//	/*add the left end point of each interval*/
//	//	for(i=0; i<nintervals; i++)
//	//	{
//	//		sortedboundpts[index].i=intervals[i]->start.i;
//	//		sortedboundpts[index].j=intervals[i]->start.j;
//	//		index++;
//	//	}
//	//	/*add the right end point of each interval in the reversed order*/
//	//	for(i=nintervals-1; i>=0; i--)
//	//	{
//	//		sortedboundpts[index].i=intervals[i]->end.i;
//	//		sortedboundpts[index].j=intervals[i]->end.j;
//	//		index++;
//	//	}
//	//}
//
//	/*get the world coordinates for each pixel on the boundary
//	Here we assume a square windows is used
//	*/
//	//void get_worldcoords_intervals(double winstartx, double winendx, 
//	//	double winstarty, double winendy, int xsize, int ysize)
//	//{
//	//	get_worldcoords_unitintervals(winstartx, winendx, winstarty, winendy, xsize, ysize);
//
//	//	int i;
//
//	//	for(i=0; i<nintervals; i++)
//	//	{
//	//		intervals[i]->start_coord.x=(intervals[i]->start.i)*dx;
//	//		intervals[i]->start_coord.y=(intervals[i]->start.j)*dy;
//	//	}
//	//}
//
//	void get_worldcoords_unitintervals(double winstartx, double winendx, 
//		double winstarty, double winendy, int xsize, int ysize)
//	{
//		double xrang=winendx-winstartx;
//		double yrang=winendy-winstarty;
//
//		dx=xrang/(xsize-1);
//		dy=yrang/(ysize-1);
//	}
//
//	/*
//	We probably don't need such a dense boundary, 
//	we can down sample the obtain boundary
//	*/
//
//	//void downsamp(int rate)
//	//{
//	//	if(boundarypts!=NULL)
//	//		free(boundarypts);
//	//	boundarypts=(PtOnBoundary*)malloc(sizeof(PtOnBoundary)*2*nintervals);
//	//	nboundarypts = 0;
//
//	//	int i;
//	//	for(i=0; i<2*nintervals; i++)
//	//	{
//	//		if(i%rate!=0)
//	//			continue;
//
//	//		boundarypts[nboundarypts].x=sortedboundpts[i].i*dx;
//	//		boundarypts[nboundarypts].y=sortedboundpts[i].j*dx;
//	//		nboundarypts++;
//	//	}
//
//	//	/*add the last point anyway*/
//	//	if((nintervals-1)%rate != 0)
//	//	{
//	//		boundarypts[nboundarypts].x=sortedboundpts[nintervals-1].i*dx;
//	//		boundarypts[nboundarypts].y=sortedboundpts[nintervals-1].j*dx;
//	//		nboundarypts++;
//	//	}
//	//}
//};
//
//
/*we use image processing algorithms to detect the boundaries 
and order the boundary pixels*/
/*for the convenience of implementation, we set the image size is 
always 512x512*/
class ImgProcess{
public:
	unsigned char imgs[512][512][3], origin_img[512][512][3];
	int labels[512][512];
	int *parents;
	int curMaxLabels;
	int label_count;

	int *npixels_regs;  //how many pixels in each region
	int nregs;          //how many simply connected regions

	int *region_label;

	//constructor; for the initialization of the image
	ImgProcess()
	{
		int i,j;
		for(i=0; i<512; i++)
		{
			for(j=0; j<512; j++)
			{
				imgs[i][j][0]=
				imgs[i][j][1]=
				imgs[i][j][2]=255.;
			}
		}

		npixels_regs=NULL;
		parents=NULL;
		region_label=NULL;
	}

	ImgProcess(unsigned char img[512][512][3])
	{
		int i, j;
		for(i=0; i<512; i++)
		{
			for(j=0; j<512; j++)
			{
				imgs[i][j][0]=img[i][j][0];
				imgs[i][j][1]=img[i][j][1];
				imgs[i][j][2]=img[i][j][2];
			}
		}
		npixels_regs=NULL;
		parents=NULL;
		region_label=NULL;
	}
	
	ImgProcess(unsigned char *img)
	{
		int i, j;
		for(i=0; i<512; i++)
		{
			for(j=0; j<512; j++)
			{
				origin_img[i][j][0]=imgs[i][j][0]=img[3*(i*512+j)];
				origin_img[i][j][1]=imgs[i][j][1]=img[3*(i*512+j)+1];
				origin_img[i][j][2]=imgs[i][j][2]=img[3*(i*512+j)+2];
			}
		}
		npixels_regs=NULL;
		parents=NULL;
		region_label=NULL;
	}

	~ImgProcess()
	{
		if(parents!=NULL)
			free(parents);

		if(region_label!=NULL)
			free(region_label);
		
		//if(npixels_regs!=NULL)
		//	free(npixels_regs);
	}

	//member functions

	void mean_filter();
	void mean_filter(int windowsize);
	void mean_filter_usesame(int windowsize);

	void label_different_regions();

	void label_boundary_pixels();

	void extract_regions();

	friend class ImgBoundary;

};

class PixelBoundary
{
public:
	Pixel *boundpts;
	int nboundpts;
	int curMaxNum;


	PixelBoundary(int initsize=500)
	{
		if(initsize==0)
		{
			boundpts=NULL;
			nboundpts=0;

			return;
		}

		boundpts=(Pixel *)malloc(sizeof(Pixel)*initsize);
		nboundpts=0;
		curMaxNum=initsize;

	}

	~PixelBoundary()
	{
		if(boundpts!=NULL)
		{
			free(boundpts);
		}

	}
	
	void copy_boundary(PixelBoundary *otherb)
	{
		this->nboundpts=otherb->nboundpts;
		this->curMaxNum=otherb->curMaxNum;

		if(boundpts != NULL)
			free(boundpts);

		boundpts=(Pixel*)malloc(sizeof(Pixel)*otherb->nboundpts);
		int i;
		for(i=0;i<nboundpts;i++)
		{
			this->boundpts[i].x=otherb->boundpts[i].x;
			this->boundpts[i].y=otherb->boundpts[i].y;
			this->boundpts[i].visited=otherb->boundpts[i].visited;
		}
	}

	//List operation

	void addNew(int x, int y)
	{
		if(isFull())
		{
			if(!extend())
				exit(-1);
		}

		boundpts[nboundpts].x=x;
		boundpts[nboundpts].y=y;
		//boundpts[nboundpts].listindex=nboundpts;
		boundpts[nboundpts].visited=false;
		nboundpts++;
	}

	bool isFull()
	{
		if(nboundpts>=curMaxNum) return true;
		return false;
	}

	bool extend(int step=1000)
	{
		boundpts=(Pixel*)realloc(boundpts, sizeof(Pixel)*(curMaxNum+step));
		if(boundpts==NULL)
			return false;// reallocation failed

		curMaxNum+=step;
		return true;
	}

	void inverse_list()
	{
		Pixel *temp;
		temp=(Pixel*)malloc(sizeof(Pixel)*nboundpts);
		int i;
		int index=0;
		for(i=nboundpts-1;i>=0; i--)
		{
			temp[index].x=boundpts[i].x;
			temp[index].y=boundpts[i].y;
			temp[index].visited=boundpts[i].visited;
			index++;
		}

		for(i=0;i<nboundpts;i++)
		{
			boundpts[i].x=temp[i].x;
			boundpts[i].y=temp[i].y;
			boundpts[i].visited=temp[i].visited;
		}

		free(temp);
	}
};

class ImgBoundary
{
public:
	PixelBoundary *boundaries;
	int nboundaries;
	int curMaxNum;

	PixelBoundary *allboundpixels; /*store all the boundary pixels*/
	PixelBoundary *boundarylist;
	//int *pos_in_boundpixels;
	int index_in_boundarypixel_list[512][512];

	unsigned char *output;
	unsigned char *streetmap_background;

	ImgBoundary()
	{
		boundaries=NULL;
		output=NULL;
		allboundpixels=NULL;
		boundarylist=NULL;
		//pos_in_boundpixels=NULL;

		streetmap_background=NULL;
	}
	~ImgBoundary();
	//~ImgBoundary()
	//{
	//	if(output!=NULL)
	//		free(output);

	//	if(boundaries!=NULL)
	//		delete [] boundaries;

	//	if(allboundpixels!=NULL)
	//		delete allboundpixels;

	//	if(streetmap_background!=NULL)
	//		free(streetmap_background);

	//}

	//List operation

	//find out the boundaries
	void get_bound_pixels(ImgProcess *imgprocess);
	int find_boundary_label(ImgProcess *imgprocess, int i, int j);

	void get_allboundary_pixels(ImgProcess *imgprocess);
	void naive_find_sortedboundaries(ImgProcess *imgprocess);
	void naive_search_one_boundary(ImgProcess *imgprocess, int i, int j, PixelBoundary *curboundary);

	//sort the boundaries
	void sort_boundaries();
	void sort_aboundary(int);

	void render_result(ImgProcess *imgprocess);

	friend class MapBoundaryList;
};

/*
The data structure to store the boundary information under the 
world coordinate system
*/
class MapBoundary{
public:
	PtOnBoundary **pts;
	int nelems;
	int curMaxNum;

	MapBoundary(int initsize=100)
	{
		if(initsize==0)
		{
			pts=NULL;
			nelems = 0;
			return;
		}

		pts=(PtOnBoundary**)malloc(sizeof(PtOnBoundary*)*initsize);
		int i;
		for(i=0;i<initsize;i++)
			pts[i]=NULL;
		nelems=0;
		curMaxNum=initsize;
	}

	~MapBoundary()
	{
		if(pts!=NULL)
		{
			int i;
			for(i=0;i<nelems;i++)
			{
				if(pts[i]!=NULL)
				{
					free(pts[i]);
					pts[i]=NULL;
				}
			}
			free(pts);
		}
	}

	bool isFull()
	{
		if(nelems>=curMaxNum) return true;
		return false;
	}

	bool extend(int step=100)
	{
		PtOnBoundary **temp=pts;
		pts=(PtOnBoundary**)malloc(sizeof(PtOnBoundary*)*(curMaxNum+step));
		if(pts==NULL)
			return false;

		int i;
		for(i=0;i<curMaxNum;i++)
			pts[i]=temp[i];

		for(i=curMaxNum;i<curMaxNum+step;i++)
			pts[i]=NULL;

		curMaxNum+=step;
		free(temp);
		return true;
	}

	void addNew(PtOnBoundary *pt)
	{
		if(isFull())
		{
			if(!extend())
				exit(-1);
		}

		pts[nelems]=pt;
		nelems++;
	}

	void copy_boundary(MapBoundary *otherb)
	{
		this->nelems=otherb->nelems;
		this->curMaxNum=otherb->curMaxNum;

		if(pts != NULL)
			free(pts);

		pts=(PtOnBoundary**)malloc(sizeof(PtOnBoundary*)*otherb->nelems);
		int i;
		for(i=0;i<nelems;i++)
		{
			this->pts[i]=(PtOnBoundary*)malloc(sizeof(PtOnBoundary));
			this->pts[i]->x=otherb->pts[i]->x;
			this->pts[i]->y=otherb->pts[i]->y;
		}
	}
};


/*The list class to maintain all the obtained boundaries of the map
from the map image*/
class MapBoundaryList
{
public:
	MapBoundary *mapboundarylist;
	int nmapboundaries;
	int curMaxNum;
	double dx, dy;

	MapBoundaryList(int initsize=10)
	{
		if(initsize==0)
		{
			mapboundarylist=NULL;
			nmapboundaries=0;
			return;
		}
		mapboundarylist=new MapBoundary[initsize];
		nmapboundaries=0;
		curMaxNum=initsize;
	}

	~MapBoundaryList()
	{
		if(mapboundarylist!=NULL)
			delete [] mapboundarylist;
	}

	bool extend(int step=10)
	{
		/*extend the space*/
		MapBoundary *temp=mapboundarylist;
		mapboundarylist=new MapBoundary[curMaxNum+step];
		if(mapboundarylist==NULL)
			return false;//memory reallocation failed

		/*copy previous elements*/
		for(int j=0;j<curMaxNum;j++)
		{
			mapboundarylist[j].copy_boundary(&temp[j]);
		}

		curMaxNum+=10;

		delete [] temp;

		return true;
	}

	//
	void get_worldcoords_unitintervals(double winstartx, double winendx, 
		double winstarty, double winendy, int xsize, int ysize)
	{
		double xrang=winendx-winstartx;
		double yrang=winendy-winstarty;

		dx=xrang/(xsize-1);
		dy=yrang/(ysize-1);
	}

	void downsamp_imgboundaries(ImgBoundary *imgboundary, int rate);
	void display_boundaries(GLenum mode);
	void index_allboundarypts();
};

bool is_boundarypixel(int x, int y);
void get_imgbound_approDir();
void cal_init_multibound_verts_onebound(int );
void convert_one_bound_to_one_traj(int id);
void convert_mapbounds_to_trajs();

