/*
Taumap.h
*/

#include "lib/icVector.h"
#include "lib/nr.h"


typedef struct EdgePt{
	double pos[3];
	double time;
	int traj_id;
	int tri_id;
}EdgePt;




class EdgePtList
{
public:
	EdgePt **edgepts;
	int nedgepts;
	int curMaxNum;

	EdgePtList(int size = 10)
	{
		edgepts = (EdgePt **)malloc(sizeof(EdgePt *)*size);
		curMaxNum = size;
		nedgepts = 0;

		if(edgepts == NULL)
		{
			curMaxNum = 0;
			exit(-1);
		}

		for(int i = 0; i < size; i++)
			edgepts[i] = NULL;
	}

	~EdgePtList()
	{
		if(edgepts != NULL)
		{
			int i;
			for(i = 0; i < curMaxNum; i++)
			{
				if(edgepts[i] != NULL)
					free(edgepts[i]);
			}

			free(edgepts);
		}
	}

	inline bool isEmpty()  //judge whether the list is empty
		{
		if(nedgepts == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(nedgepts == curMaxNum) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		EdgePt **temp = edgepts;
		edgepts = (EdgePt **) malloc(sizeof(EdgePt *) * (curMaxNum + step));
		if( temp == NULL)
		{
			//fail
			curMaxNum = 0;
			exit(-1);

			edgepts = temp;
			return false;
		}

		int i;
		for(i = 0; i < curMaxNum; i++)
			edgepts[i] = temp[i];
		for(i = curMaxNum; i < curMaxNum+step; i++)
			edgepts[i] = NULL;
		curMaxNum += step;

		free(temp);
		return true;
	}
};



/*Konstantin's accurate tau map calculation 07/25/07*/

typedef struct TauLineSeg{
	int tri1, tri2;                    //which triangle this line segment locates in
	double start[2], end[2];            //local coordinates for start and end points
	double gstart[3], gend[3];          //global coordinates for start and end points
} TauLineSeg;


/*we need a data structure to store the images of the samples of
an edge in order*/

typedef struct AnEdgeSampleImg{
	double l_pos[2];
	int trinagle;
}AnEdgeSampleImg;

class EdgeSampleImgList{
public:
	AnEdgeSampleImg **sampleimgs;
	int nsampleimgs;
	int curMaxNum;
	int edge_index;                     //the corresponding edge index

	EdgeSampleImgList(int initsize = 3)
	{
		sampleimgs=(AnEdgeSampleImg **)malloc(sizeof(AnEdgeSampleImg*)*initsize);
		nsampleimgs = 0;
		curMaxNum = initsize;

		if(sampleimgs == NULL) //fail
			exit(-1);

		for(int i=0; i<curMaxNum; i++)
			sampleimgs[i]=NULL;
	}

	virtual ~EdgeSampleImgList()
	{
		if(sampleimgs != NULL)
		{
			int i;
			for(i = 0; i < curMaxNum; i++)
			{
				if(sampleimgs[i] != NULL)
					free(sampleimgs[i]);
			}

			free(sampleimgs);
		}
	}
	
	inline bool isEmpty()  //judge whether the list is empty
		{
		if(nsampleimgs == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(nsampleimgs == curMaxNum) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		AnEdgeSampleImg **temp = sampleimgs;
		sampleimgs = (AnEdgeSampleImg **) malloc(sizeof(AnEdgeSampleImg *) * (curMaxNum + step));
		if( sampleimgs == NULL)
		{
			//fail
			curMaxNum = 0;
			exit(-1);

		}

		for(int i = 0; i < curMaxNum; i++)
			sampleimgs[i] = temp[i];
		curMaxNum += step;

		//free(temp);
		return true;
	}
};



class EdgeInfoLookUp{
public:
	EdgeSampleImgList **edgelookuptable;
	int nelems;
	int curMaxNum;

	EdgeInfoLookUp(int initsize = 30)
	{
		edgelookuptable = new EdgeSampleImgList*[initsize];
		nelems = 0;
		curMaxNum = initsize;

		for(int i=0; i<initsize; i++)
			edgelookuptable[i] = new EdgeSampleImgList(10);

	}

	virtual ~EdgeInfoLookUp()
	{
	}

	/*other functionality*/

	EdgeSampleImgList *search_by_edgeindex(int edgeindex)
	{
		int i;
		for(i=0; i<nelems; i++)
		{
			if(edgeindex == edgelookuptable[i]->edge_index)
				return edgelookuptable[i];
		}

		return NULL;
	}

	/*remove the information of an edge from the list*/
	bool remove_by_edgeindex(int edgeindex) 
	{
		int i, pos = 0;
		for(i=0; i<nelems; i++)
		{
			if(edgeindex == edgelookuptable[i]->edge_index)
			{
				pos = i;
				break;
			}
		}

		if(pos >= nelems)
			return false;

		/*remove the element*/
		free(edgelookuptable[i]);
		for(i=pos; i<nelems-1; i++)
		{
			edgelookuptable[i]=edgelookuptable[i+1];
		}
		nelems--;
	}


	inline bool isEmpty()  //judge whether the list is empty
		{
		if(nelems == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(nelems == curMaxNum) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		EdgeSampleImgList **temp = edgelookuptable;
		edgelookuptable = new EdgeSampleImgList * [curMaxNum + step];
		if( temp == NULL)
		{
			//fail
			curMaxNum = 0;
			exit(-1);

		}

		for(int i = 0; i < curMaxNum; i++)
			edgelookuptable[i] = temp[i];
		curMaxNum += step;

		delete [] temp;
		return true;
	}
};


typedef struct TraceSamplePt{
	float pos[2];
	float time;
	float color[3]/*, b_color[3]*/;
	int tri;
	int edge;
}TraceSamplePt;


class TraceSamplePtList
{
public:
	TraceSamplePt **samplepts;
	int nelems;
	int curMaxNum;

	TraceSamplePtList(int initsize = 6000)
	{
		samplepts = (TraceSamplePt**)malloc(sizeof(TraceSamplePt*)*initsize);
		if(samplepts == NULL)
			exit(-1);
		nelems = 0;
		curMaxNum = initsize;
	}

	inline bool addNew(TraceSamplePt *s)
	{
		if(isFull())
		{
			if(!extend())
				exit(-1); //fail to reallocate
		}

		/**/
		samplepts[nelems] = s;
		nelems++;
	}

	inline bool isEmpty()  //judge whether the list is empty
		{
		if(nelems == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(nelems == curMaxNum) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		TraceSamplePt **temp = samplepts;
		samplepts = (TraceSamplePt**)malloc(sizeof(TraceSamplePt*)*(curMaxNum + step));
		if( temp == NULL)
		{
			//fail
			curMaxNum = 0;
			exit(-1);

		}

		for(int i = 0; i < curMaxNum; i++)
			samplepts[i] = temp[i];
		curMaxNum += step;

		free(temp);
		return true;
	}
};



void init_samplepts_tautracing();
void assign_color();
void assign_vertex_color();

void compute_density();
void assign_density_colors();
void test_pro(double tau, double s_tau);

