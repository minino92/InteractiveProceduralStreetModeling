/*
This file contains the declaration of class EvenStreamlinePlace

Created and modified by Guoning Chen
copyright @2007
*/

//#include "VField.h"
//#include "common_routines.h"

/** Trajectory **/

//#include "vfdatastructure.h"
//#include "Localtracing.h"



class DynList_Int{
public:
	//member variables
	int *elems;
	int  nelems;
	int  curMaxNum;
	int  extendstep;

	//member functions
	//constructions and destruction
	DynList_Int(int MaxNum = 500) 
	{
		elems = new int [MaxNum];
		curMaxNum = MaxNum;
		nelems = 0;
	}
	~DynList_Int()
	{delete [] elems; }

	inline bool extend(int step = 100) 
	{
		int *temp = elems;
		elems = new int[curMaxNum + step];
		if(elems == NULL)
		{
			elems = temp;
			return false;
		}

		for(int i = 0; i < curMaxNum; i++)
			elems[i] = temp[i];

		curMaxNum += step;

		delete [] temp;
		return true;
	}

	inline bool add_New(int t)
	{
		if(is_Full())
		{
			if(!extend())
			{
				return false;
			}
		}

		if(!is_repeated(t))
		{
			elems[nelems] = t;
			nelems++;
			return true;
		}

		return false;
	}
	
	inline bool add_New_2(int t)
	{
		if(nelems>=curMaxNum)
		{
			if(!extend())
			{
				return false;
			}
		}

		//if(!is_repeated_elem(elems, t, nelems))
		//{
		//	elems[nelems] = t;
		//	nelems++;
		//	return true;
		//}
			elems[nelems] = t;
			nelems++;
			return true;

		return false;
	}

	inline bool is_Full()
	{
		return nelems>=curMaxNum;
	}

	inline bool  del_Elem(int t)
	{
		int i, pos = 0;
		for(i = 0; i < nelems; i++)
		{
			if(elems[i] == t)
			{
				pos = i;
				break;
			}
		}
		if(pos >= nelems)
			return false;

		for(i = pos; i < nelems-1; i++)
			elems[i] = elems[i+1];

		nelems--;
		return true;
	}


	inline bool  del_Last()
	{
		if(nelems <= 0) return false;
		nelems--;
		return true;
	}

	inline bool  is_repeated(int newelem)
	{
		int i;
		for(i=0; i<nelems; i++)
		{
			if(elems[i] == newelem)
				return true;
		}
		return false;
	}
};







typedef struct Seed{
	//double pos[3];                //the coordinates of the seed point
	double pos[2];                  //the coordinates of the seed point
	int triangle;                   //the triangle contains the seed point
	//int which_traj;                 //record which trajectory it locates on
	unsigned char state;            //0—active;1—inactive; 2—can not growing
	double weight;
}Seed; // end of Seed structure

class SeedList{
public:
	Seed **seeds;
	int nseeds;
	int curMaxNumSeeds;
	double max_weight, min_weight;

	// The similar list operations
	SeedList(int initsize = 3000) //construction
	{
		seeds = (Seed **)malloc(sizeof(Seed *)*initsize);
		curMaxNumSeeds = initsize;
		nseeds = 0;

		if(seeds == NULL)
		{
			//char rout[256], var[256];
			//sprintf(rout, "%s", "SeedList Constructor");
			//sprintf(var, "%s", "seeds");

			//write_mem_error(rout, var, 0);
			curMaxNumSeeds = 0;
			exit(-1);

		}
		int i;
		for(i = 0; i < initsize; i++)
			seeds[i] = NULL;
	} 

	~SeedList()
	{
		int i;
		for(i = 0; i < curMaxNumSeeds; i++)
		{
			if(seeds[i] != NULL)
				free(seeds[i]);
		}
		free(seeds);
	}

	inline void copy_otherseedList(SeedList *otherseeds)
	{
		int i;
		if(otherseeds->nseeds>curMaxNumSeeds)
		{
			if(!extend(otherseeds->nseeds-curMaxNumSeeds))
				exit(-1);
		}

		/*copy element by element*/
		for(i=0;i<otherseeds->nseeds;i++)
		{
			if(seeds[i]==NULL)
				seeds[i]=(Seed *) malloc(sizeof(Seed));
			
			seeds[i]->pos[0]=otherseeds->seeds[i]->pos[0];
			seeds[i]->pos[1]=otherseeds->seeds[i]->pos[1];

			seeds[i]->triangle=otherseeds->seeds[i]->triangle;
			seeds[i]->state=otherseeds->seeds[i]->state;
			seeds[i]->weight=otherseeds->seeds[i]->weight;
		}
		nseeds=otherseeds->nseeds;
	}

	inline void copyandappend_otherseedList(SeedList *otherseeds)
	{
		int i;
		if(otherseeds->nseeds+nseeds>curMaxNumSeeds)
		{
			if(!extend(otherseeds->nseeds+nseeds-curMaxNumSeeds))
				exit(-1);
		}

		/*copy element by element*/
		for(i=nseeds;i<nseeds+otherseeds->nseeds;i++)
		{
			if(seeds[i]==NULL)
				seeds[i]=(Seed *) malloc(sizeof(Seed));
			
			seeds[i]->pos[0]=otherseeds->seeds[i-nseeds]->pos[0];
			seeds[i]->pos[1]=otherseeds->seeds[i-nseeds]->pos[1];

			seeds[i]->triangle=otherseeds->seeds[i-nseeds]->triangle;
			seeds[i]->state=otherseeds->seeds[i-nseeds]->state;
		}
		nseeds+=otherseeds->nseeds;
	}

	//add a new vertex to the end of the list, if it succeeds, return true
	inline bool append(Seed *s)
	{
		if(isFull ())
			if(!extend(100))
				return false;             //if not enough memory available, return false
		seeds[nseeds] = s;
		//copyElem(s, polist[nporbits]);
		nseeds++;
		return true;
	} 

	/*  insert a new seed to the proper position of the list according to its weight  */
	inline bool sorted_add(Seed *s)
	{
		if(nseeds>=curMaxNumSeeds)
		{
			if(!extend())
				//exit(-1);
				return false;
		}

		/*sorted insert!*/
		int i, pos=0;
		for(i=0; i<nseeds; i++)
		{
			if(s->weight>seeds[i]->weight)
			{
				pos = i;
				break;
			}
		}
		if(i>=nseeds) pos=nseeds;
		for(i=nseeds; i>pos; i--)
			seeds[i]=seeds[i-1];
		seeds[pos]=s;

		nseeds ++;
		return true;
	}

	inline void update_max_min_weights()
	{
		int i;
		if(nseeds==0) return;
		max_weight=min_weight=seeds[0]->weight;
		for(i=1;i<nseeds;i++)
		{
			if(seeds[i]->weight>max_weight) max_weight=seeds[i]->weight;
			if(seeds[i]->weight<min_weight) min_weight=seeds[i]->weight;
		}
	}

	inline bool del_End() //delete the vertex at the end of the list
	{
		if(isEmpty())  return false;
		nseeds --;
		return true;
	} 

	inline void copy_Elem(Seed *to, Seed *from)
	{
		to->pos[0]=from->pos[0];
		to->pos[1]=from->pos[1];
		to->triangle=from->triangle;
		to->state=from->state;
	}

	//delete the corresponding  vertex, if it succeeds, return true
	inline bool del_Node(Seed *s) 
	{
		if(isEmpty())  return false;

		//find the vertex, if find it, delete and move the following vertices forward
		//otherwise, return false;

		int i, pos = -1;

		for(i = 0; i < nseeds; i++)
		{
			if(seeds[i] == s)
			{
				pos = i;
				break;
			}
		}

		if(pos == -1) return false;

		//delete it
		for(i = pos; i < nseeds-1; i++)
		{
			//we need a copy function
			copy_Elem(seeds[i], seeds[i+1]);
		}

		nseeds--;

		return true;
	} 
	
	inline bool del_Node_byindex(int pos) 
	{
		int i;
		if(isEmpty())  return false;

		//find the vertex, if find it, delete and move the following vertices forward
		//otherwise, return false;

		if(pos == -1) return false;

		//delete it
		for(i = pos; i < nseeds-1; i++)
		{
			//we need a copy function
			copy_Elem(seeds[i], seeds[i+1]);
		}

		nseeds--;

		return true;
	} 

	inline bool isEmpty()  //judge whether the list is empty
		{
		if(nseeds == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(nseeds == curMaxNumSeeds) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		Seed **temp = seeds;
		seeds = (Seed **) malloc(sizeof(Seed *) * (curMaxNumSeeds + step));
		if( temp == NULL)
		{
			//fail
			//char rout[256], var[256];
			//sprintf(rout, "%s", "SeedList::extend");
			//sprintf(var, "%s", "seeds");

			//write_mem_error(rout, var, 1);
			curMaxNumSeeds = 0;
			seeds = temp;
			exit(-1);
			return false;
		}

		int i;
		for(i = 0; i < curMaxNumSeeds; i++)
			seeds[i] = temp[i];

		for(i=curMaxNumSeeds; i<curMaxNumSeeds+step; i++)
			seeds[i]=NULL;

		curMaxNumSeeds += step;

		free(temp);
		return true;
	}

	inline void reset()
	{
		nseeds = 0;
	}
}; //end of SeedList class





/** Sample point structure **/
typedef struct SamplePt{
	//double gpt[3];
	double gpt[2];
	int triangle;
	int traj;                                //which trajectory this sample falls on
}SamplePt; //end of SamplePt class





class SamplePtList{
public:
	SamplePt **samples;
	int nsamples;
	int curMaxNumSamplePts;

	SamplePtList(int initsize = 1000) //construction
	{
		samples = (SamplePt **)malloc(sizeof(SamplePt *)*initsize);
		curMaxNumSamplePts = initsize;
		nsamples = 0;

		if(samples == NULL)
		{
			//char rout[256], var[256];
			//sprintf(rout, "%s", "SamplePtList Constructor");
			//sprintf(var, "%s", "samples");

			//write_mem_error(rout, var, 0);
			curMaxNumSamplePts = 0;
			exit(-1);
		}
		
		for(int i = 0; i < initsize; i++)
			samples[i] = NULL;
	} 

	~SamplePtList()
	{
		if(samples!= NULL)
		{
			for(int i = 0; i < curMaxNumSamplePts; i++)
			{
				if(samples[i] != NULL)
					free(samples[i]);
			}
			free(samples);
		}
	}

	//add a new vertex to the end of the list, if it succeeds, return true
	inline bool append(SamplePt *s)
	{
		if(isFull ())
			if(!extend(100))
				return false;             //if not enough memory available, return false
		samples[nsamples] = s;
		//copyElem(s, polist[nporbits]);
		nsamples++;
		return true;
	} 

	inline bool del_End() //delete the vertex at the end of the list
	{
		if(isEmpty())  return false;
		nsamples --;
		return true;
	} 

	inline void copy_Elem(SamplePt *s, SamplePt *d)
	{
	}

	//delete the corresponding  vertex, if it succeeds, return true
	inline bool del_Node(SamplePt *s) 
	{
		if(isEmpty())  return false;

		//find the vertex, if find it, delete and move the following vertices forward
		//otherwise, return false;

		int i, pos = -1;

		for(i = 0; i < nsamples; i++)
		{
			if(samples[i] == s)
			{
				pos = i;
				break;
			}
		}

		if(pos == -1) return false;

		//delete it
		for(i = pos; i < nsamples-1; i++)
		{
			//we need a copy function
			copy_Elem(samples[i], samples[i+1]);
		}

		nsamples--;

		return true;
	} 

	inline bool isEmpty()  //judge whether the list is empty
		{
		if(nsamples == 0)   return true;
		return false;
	}

	inline bool isFull()
	{
		if(nsamples == curMaxNumSamplePts) return true;
		return false;
	}

	//extend the original list, if it succeeds, return true
	inline bool extend(int step = 100)
	{
		SamplePt **temp = samples;
		samples = (SamplePt **) malloc(sizeof(SamplePt *) * (curMaxNumSamplePts + step));
		if( temp == NULL)
		{
			//fail
			//char rout[256], var[256];
			//sprintf(rout, "%s", "SamplePtList::extend");
			//sprintf(var, "%s", "samples");

			//write_mem_error(rout, var, 1);
			curMaxNumSamplePts = 0;
			samples = temp;
			exit(-1);
			return false;
		}

		int i;
		for(i = 0; i < curMaxNumSamplePts; i++)
			samples[i] = temp[i];
		for(i = curMaxNumSamplePts; i < curMaxNumSamplePts+step; i++)
			samples[i] = NULL;
		curMaxNumSamplePts += step;

		free(temp);
		return true;
	}

	inline void reset()
	{
		nsamples = 0;
	}
}; //end of SamplePtList class





//extern struct CurvePoints;

class EvenStreamlinePlace
{
public:
	TrajectoryList *evenstreamlines;
	SeedList *seedpts;
	SamplePtList **samplepts;
	int cur_traj; /*the index of current streamline*/
	CurvePoints *tracing_points /*= (CurvePoints*) malloc(sizeof(CurvePoints) * 805)*/;
    int num_tracingpoints ;
	DynList_Int *trianglelist;
	//EdgeList *geo_boundary;
	//icVector2 tenline_dir_global;

	/*important parameters to control the even placement of the streamlines*/
	double streamlinelength;
	double dsep;
	double percentage_dsep;
	double discsize;
	double sample_interval;
	int every_nsample;
	double loopdsep;
	double dist2sing;
	double seeddist;
	double minstartdist;
	bool majororminor;

	/*to record the sample that stops the tracing of a tensor line*/
	int which_triangle;
	double samp[2];

	/*to record current computing position in the seed list
	and in the trajectory list.
	This is really useful in alternative tracing*/
	int cur_seed_pos/*, cur_traj_pos*/;

	/*   for recording the number of new computed tensor lines in the smaller region   */
	int nnewlines_inReg;

	/*record the position of the seed point for each tensor line*/
	int *seedposition_ineachtensorline;

	EvenStreamlinePlace(bool type, int initsize=300)
	{
		majororminor = type;
		evenstreamlines = new TrajectoryList(initsize);
		samplepts = (SamplePtList **)malloc(sizeof(SamplePtList*)*evenstreamlines->curMaxNumTrajs);

		seedposition_ineachtensorline=(int*)malloc(sizeof(int)*evenstreamlines->curMaxNumTrajs);

		if(samplepts == NULL || seedposition_ineachtensorline==NULL)
		{
			exit(-1);
		}
		
		for(int i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
		{
			samplepts[i] = new SamplePtList(100);

			if(samplepts[i] == NULL)
			{
				exit(-1);
			}

			seedposition_ineachtensorline[i]=0;
		}

		trianglelist = NULL;
		seedpts = NULL;
	    //geo_boundary = NULL;
	}

	~EvenStreamlinePlace()
	{
		delete evenstreamlines;

		if(samplepts != NULL)
		{
			for(int i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
			{
				if(samplepts[i] != NULL)
					free(samplepts[i]);
			}
			free(samplepts);
		}

		if(seedposition_ineachtensorline!=NULL)
			free(seedposition_ineachtensorline);
	}
	
	/*initialize the element inside the lists*/
	inline void init()
	{
		int i, j;
		for(i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
		{
			//evenstreamlines->trajs[i] = (Trajectory*)malloc(sizeof(Trajectory));
			evenstreamlines->trajs[i] = new Trajectory(i);

			if(evenstreamlines->trajs[i] == NULL)
			{
				//char rout[256], var[256];
				//sprintf(rout, "%s", "EvenStreamlinePlace::init");
				//sprintf(var, "%s", "evenstreamlines->trajs[i]");

				//write_mem_error(rout, var, 0);
				exit(-1);
			}
		}

		for(i = 0; i < evenstreamlines->curMaxNumTrajs; i++)
		{
			for(j = 0; j < samplepts[i]->curMaxNumSamplePts; j++)
			{
				//samplepts[i]->samples[j] = new SamplePt[1];
				samplepts[i]->samples[j] = (SamplePt *)malloc(sizeof(SamplePt));
				if(samplepts[i]->samples[j] == NULL)
				{
					//char rout[256], var[256];
					//sprintf(rout, "%s", "EvenStreamlinePlacement::init");
					//sprintf(var, "%s", "samplepts[i]->samples[j]");

					//write_mem_error(rout, var, 0);
					exit(-1);
				}
			}
		}
	}

	void reset_placement_quad();

	void set_default_parameters(bool fieldtype);

	//void place_streamlines(int num_initial, double dsep, double percentage_dsep, 
	//			  double discsize, double sample_interval, int every_nsample, 
	//			  double loopdsep, double dist2sing, double streamlinelength, 
	//			  double seeddist, double minstartdist, int flag);
	void place_streamlines(int type, bool brushon);
	void place_streamlines(int type, bool brushon, SeedList *preseeds);
	void place_tensorlines_alternatively(int type, bool brushon, SeedList *preseeds);

	void place_tensorlines_inReg(int type, SeedList *ini_seeds);
	bool grow_a_tensorline_inReg(double seed_p[2], int triangle, double dtest, 
											double discsize, double Sample_interval, 
											double loopdsep, double dist2sing, 
											double streamlinelength, int type, bool brushon);
	int trace_in_quad_inReg(int &face_id, double globalp[2], int type, 
						double dtest, double loopsep, double dist2sing, 
						double sample_interval, double discsize, int &flag);

	void cal_init_streamlines(int num, double streamlinelength, double dtest, double discsize,
					double Sample_interval, double loopdsep, double dist2sing, int type, 
					bool brushon);
	void cal_init_streamlines(double streamlinelength, double dtest, double discsize,
				double Sample_interval, double loopdsep, double dist2sing, int type, bool brushon,
				SeedList *userspecifiedseeds);

	void copy_sketchlines();
	
	void copy_majorRoads(bool);
	void sample_majorRoads_seeds(SeedList *, int);

	/*
	compute the initial set of streamlines considering separatrices and periodic orbits
	*/
	//void cal_init_streamlines_enhanced(double streamlinelength, double dtest, double discsize,
	//					double sample_interval, double dist2sing);

	bool is_valid_loop(int trajid, double sep);

	/*
	Trace one streamline and judgement the density during tracing
	*/
	bool grow_a_tensorline(double seed_p[2], int triangle, double dtest, 
											double discsize, double Sample_interval, 
											double loopdsep, double dist2sing, 
											double streamlinelength, int type, bool brushon);
	bool grow_a_tensorline_withoutinversing(double seed_p[2], int triangle, double dtest, 
											double discsize, double Sample_interval, double loopdsep, 
											double dist2sing, double streamlinelength, 
											int type, bool brushon);

	bool grow_a_majRoad(double seed_p[2], int triangle, double dtest, 
											double discsize, double Sample_interval, 
											double loopdsep, double dist2sing, 
											double streamlinelength, 
											int type, bool brushon);

	void comp_smallest_dist_to_one_mapBound(double pt[2], int cell, int boundID, 
											double &dist, double intersect[2],
											int &pos, icVector2 &appro_bound_dir);
	bool find_closest_bound(double pt[2], int cell, int &boundID, 
							double intersect[2], int &pos, icVector2 &appro_bound_dir);
	void trace_water_region(double pt[2], int &cell, double Sample_interval,
						double loopdsep, double min_waterwidth, double cross_angle,
						double dist_follow_bound, int type);

	void rebuild_lineinfo_for_a_tensorline(int id);
	void rebuild_all_lineinfo();



	void cal_seeds(int traj, double dsep, int every_nsample, int type, bool brushon);
	bool get_a_seed(double sample_p[3], int begin_triangle, int cur_traj, int type,
			  double end_p[3], int &end_triangle, double dsep, bool brushon);
	bool get_a_seed(double sample_p[2], icVector2 line_dir, 
			double end_p[2], int &end_triangle, double dsep, bool brushon);
	/*
	Trace inside a triangle for the even streamline placement
	*/
	int trace_in_quad_even(int &face_id, double globalp[2], int type, 
						double dtest, double loopsep, double dist2sing, 
						double sample_interval, double discsize, int &flag);
	int trace_majRoad_in_quad(int &face_id, double globalp[2], int type, 
					double dtest, double loopsep, double dist2sing, 
					double sample_interval, double discsize, int &flag);

	/*
	Get a set of sample points during tracing a new streamline
	This sample points are temporary and for preventing the closed loop of a streamline
	It will keep finding all the samples from 'cur_line' line segment till the end of 
	current streamline
	Note that we use local coordinates to get one sample point each time,
	after that, we need to transform it back to 3D global coordinates
	*/
	SamplePt **cal_samplepts_when_tracing(int traj, double interval, int &cur_line, int &movetonext, double &cur_length, 
							 SamplePt **samples, int &num_samples);

	/*
	Reverse the order of line segments of a specific streamline
	*/
	void reverse_streamline(int streamlineid);

	/*
	Make sure that all the line segments of the streamline have
	the same direction
	*/
	void reorder_streamline(int streamlineid);

	/*
	Update the smaple point list of all the triangles associated with the input sample list
	*/
	void update_samples_in_cell(int traj, SamplePt **samples, int num_samples);

	/*initial the samplept list for each cell*/
	void init_samplelist_in_cell(bool fieldtype);

	/*
	Perform a short local tracing from each sampling point to find a seed point,
	we still need to store the temporaylist to a global line segment list
	*/
	int trace_in_triangle_seed_quad(int &face_id, double globalp[3], int type, 
								int &flag, double pre_p[3], double cur_p[3], 
								double dsep, double &cur_length, int cur_traj);

	/*
	During seed point finding, we use local tracing with some constant step size
	so, we need to get the seed point having the exact distance to the streamline
	*/
	void cal_exact_seed(double cur_p[3], double pre_p[3], int triangle, 
						  double extra_length, double exact_seed[3]);


	/*distance judgement routines*/

	/* obtain the density information at a given point */
	double cal_approx_density_at(double p[2], int cell);
	/*
	Judge whether a given point is too close to existing fixed points except for the singid
	*/
	bool close_to_degPt_except(double p[3], int triangle, int singid, double threshold, double discsize);
	/*
	Judge whether current point is too close to existing streamline except for the set traj
	*/
	bool close_to_cur_streamline(double p[3], int triangle, int *Except_trajs, int num_trajs, 
							double separate_dist, double discsize, int Init);
							//int &which_triangle, double sample[2]);
	bool close_to_cur_streamline(double p[2], int triangle, icVector2 dir,
						int *Except_trajs, int num_trajs, 
						double separate_dist, double discsize/*, int Init*/);
	bool close_to_cur_streamline_maj(double p[3], int triangle, 
						int *Except_trajs, int num_trajs, 
						double separate_dist, double discsize, int &which_traj, int &which_samp);

	bool cal_approx_perpendicular_pt(int traj, int &cellid, double p[2], double pt[2], bool type);

	/*Using the similar streamline placement for 2D planar vector field*/
	bool close_to_cur_samplePt(double p[3], int triangle, SamplePt **samples, int num_samples,
							double separate_dist, double discsize, double sample_interval);

	/*
	Get one sample point on the regular streamline, not periodic orbit
	Different from planar case, we use local coordinates here
	*/
	bool cal_a_sample_of_streamline(int traj, int &cur_lineindex, int &movetonext,
						double curpt[2], double interval, double &cur_length);

	/*
	Add one sample to the particular triangle
	*/
	void add_sample_to_cell(int triangle, int which_traj, int which_sample, bool fieldtype);
	SampleListInTriangle *extend_cell_samplelist(SampleListInTriangle *samplepts, int nsamples);

	/*
	Initialize the information of the lines for the cells they cross
	*/
	void init_major_minor_line_info(bool type);

	/*
	save the placement results into a file
	*/
	void save_to_a_file(const char* filename);
	void load_from_a_file(const char *filename);
	
	void save_to_a_file(const char* filename, int flag);
	void load_from_a_file_enhanced(const char *filename);


	/*The following routines compute the Geodesic (Euclidean in 2D) distance*/

	int *cal_euclidean_dist(int triangle, double globalp[3], double dsep, double discsize, 
					   int *trianglelist, int &num_triangles);
	void cal_euclidean_dist_2(int triangle, double p[3], double dsep, double discsize, 
					DynList_Int *trianglelist);
	void reset_dist(int *NearbyTriangles, int num);
	void reset_dist(DynList_Int *trianglelist);
	int get_thirdVer_of_triangle(int v1, int v2, int triangle);

	void update_one_vertex(int verta, int vertb, int vertc, int triangle);


	void display_sampts_seeds();  /*test the seeds and samples*/

	bool extend_trajList();

};

void init_evenplace_ten();
void init_level1_placement();
void place();

void place_alternative_level1_makeup_seeds();
void place_alternative_level1(SeedList *inputseeds);
bool place_onetensorline_level1(bool type, SeedList *inputseeds, SeedList *outputseeds);
bool get_a_proper_seed(EvenStreamlinePlace *curplace, double samp[2], int &cell, int fieldtype);

void place_alternative_level1_loadmap();

void init_majRoad_intersectionlists();
void init_majRoadnet();

bool compute_majRoad_intersect_between_twolines(int majtraj, int majstart, int majend,
										int mintraj, int minstart, int minend,
										double intersect[2], 
										int &majlinesegid, int &minlinesegid);
bool has_edge_between_minRoad(int intersect1, int intersect2);
void compute_majRoad_intersects_in_cell(int cellid);
void compute_majRoad_intersects();
void search_for_connection_majRoad();
void add_edge_to_majRoad_intersectnode(int node, int edgeindex);


bool save_obtained_majRoads_tenLines(char *filename);
bool load_majRoads_tenLines(char *filename);
bool save_obtained_majRoads_network(char *filename);
bool load_majRoads_network(char *filename);








/*  routines for constructing the street network graph  */
void init_streetnet();
void init_streetnetwork_info_in_cells();

void compute_intersects();
bool compute_intersect_between_twolines(int majtraj, int majstart, int majend,
										int mintraj, int minstart, int minend,
										double intersect[2], 
										int &majlineid, int &minlineid);
void compute_intersects_in_cell(int cellid);

bool compute_intersect_between_twolines_sametype(int majtraj, int majstart, int majend,
										int mintraj, int minstart, int minend,
										double intersect[2], 
										int &majlinesegid, int &minlinesegid, bool majormin);
void compute_intersects_in_cell_sametype(int cellid, bool majormin);

void add_to_cell_intersectlist(int cellid, int intersectid, bool endpt);
void init_tensorline_intersectionlists();
void search_for_connection();
void add_edge_to_intersectnode(int node, int edgeindex);

void save_cur_street_network(char *filename);
bool load_a_street_network(char *filename);

/*remove dead ends*/
bool find_cell_contain_bothmajmin(double start[2], int &startcell, double dsep, 
								  icVector2 line_dir, int fieldtype);
void compute_extended_intersect(double p[2], int cellid, double intersect[2], 
								EvenStreamlinePlace *curplace, int fieldtype);
void adj_deadend_onetensorline(int fieldtype, int id, EvenStreamlinePlace *curplace);
void remove_all_deadends();
void remove_one_deadend_minRoad(int id, bool majormin, bool startorend,
								double search_dist);


/*  remove dead ends of the major roads  */
void connect_majRoads_postproc();

/*functions for fast marching*/
void fast_marching_quad();
void vis_distance();
void display_zero_level();
void init_dis_verts();
void init_verts_all();
void init_dis_cells();
void reset_vert_dis();
void init_verts_for_brush();

void DistanceFromLine(double cx, double cy, double ax, double ay ,
					  double bx, double by, double &distanceSegment,
					  double &distanceLine);
void fast_marching_quad_withDis(double disthred);
void init_boundvertlist();
void reset_narrowband();
void init_narrow_band();
void update_neighbor_Dis(int cur_v);


/*implement the brush interface*/
int get_Resolution_adp();
void get_cellstrip_curve_quad(int *celllist, int &ncells);
int get_cellID_picking(double x, double y);
int get_cellID_givencoords(double x, double y);

void add_to_shapeCtrPtsList(double x, double y, int cellid);
bool find_pos_curve(int cellid, int startid, int &pos);
void get_bound_approDir();
void get_brushinterior_ten();
void save_brush_curve(char *filename);
void load_brush_curve(char *filename);
void save_cur_brushverts();
void get_brushbuffer_quad_withDis(double disthred);

bool is_not_inland(int cellid);  /*is not in the land*/
bool is_not_inland_cell_weak(int cellid);

bool is_inveg(int cellid);
bool is_inveg_cell_weak(int cellid);
/*
Routines for finding the region blocks
*/
void init_regionblocklist();
void init_all_streetgraph_edges();
bool select_a_valid_edge(int nodeid, icVector2 normal, int &nextedge, int exceptedge);
void construct_a_regionblock(int *nodes, int nnodes, int *edges, int nedges);
bool form_a_block(int edgeid, bool orient);
int construct_blocks_for_an_edge(int edgeid);
void construct_regionblocks_edgewise();
void vis_regionblocks();
