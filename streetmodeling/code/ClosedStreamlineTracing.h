/*
ClosedStreamlineTracing.h

The head file of the close streamline tracing.
*/

#include "lib/icVector.h"

typedef struct SamplePts{
	double gpt[3];
	int triangle;
}SamplePts;

typedef struct SamplePtsList{
	SamplePts *samples;
	int num_samples;
	int MaxNumSamples;
	int colorindex;
}SamplePtsList;


bool TraceAvoidClosed(double seed_p[3], int triangle, 
					double Sample_interval, double loopdsep, int type, int &flag);

bool CloseToCurSamplePt(double p[3], int triangle, SamplePts *samples, int num_samples,
						double separate_dist, double sample_interval);

bool GetOneSamplePointfromAStreamline(int traj, int &cur_lineindex, int &movetonext,
					   double prept[2], double curpt[2], 
					   double interval, double &cur_length);

SamplePts *GetSamplePtsWhenTracing(int traj, double interval, int &cur_line, int &movetonext, double &cur_length, 
							 SamplePts *samples, int &num_samples, int &MaxNum);
int TraceInATriangleAvoidClosed(int &face_id, double globalp[2], int type, 
										 double loopsep, double sample_interval, int &flag);

void AllocateSampleptsList();

void FinalizeSampleptsList();

void InitSampleptsList();

bool GetSingularyID(double x, double y, int &singular_id);

bool CloseToLimitCycle(int saddle, double p[2], int triangle, int type, int sep_id);

int GettheSide_loc(icVector2 orient, double basis[2], double p[2]);

