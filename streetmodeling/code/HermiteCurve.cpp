////HermiteCurve.cpp

#include "stdafx.h"

#include "HermiteCurve.h"

#include "VFDataStructure.h"

typedef double Base4[4];

/**********************************************************
n :  number of interpolated points
*interpts: array of interpolate points
*output: array of output points on the curve
step: number of segments between two interpolate points
num_output: number of total output points in the output array
**********************************************************/
void Hermitecurve(int n, ctr_point *interpts, icVector2 *T, ctr_point *output, int step, int &num_output)
{
	int i, j;
	Base4 *Base = (Base4 *) malloc(sizeof(Base4) * step);
	double s/*, h1, h2, h3, h4*/;

	//icVector2 *T = (icVector2 *) malloc(sizeof(icVector2) * n);
	icVector2 Pstart, Pend, T1, T2;

	////calculate the blending functions once and store them into an array
	for(j = 0; j < step; j++)
	{
		s = (double) j / (double) step;
		Base[j][0] = 2 * s * s * s - 3 * s * s + 1;
		Base[j][1] = -2 * s * s * s + 3 * s * s;
		Base[j][2] = s * s * s - 2 * s * s + s;
		Base[j][3] = s * s * s - s * s;
	}

	////We may first store the tagent vectors on the interpolate points into an array
	for(i = 0; i < n; i++)
	{
		if( i == 0)
		{
			Pstart.entry[0] = interpts[n-1].x;
			Pstart.entry[1] = interpts[n-1].y;

			Pend.entry[0] = interpts[1].x;
			Pend.entry[1] = interpts[1].y;
		}
		else if( i == n-1 )
		{
			Pstart.entry[0] = interpts[n-2].x;
			Pstart.entry[1] = interpts[n-2].y;

			Pend.entry[0] = interpts[0].x;
			Pend.entry[1] = interpts[0].y;
		}
		else{
			Pstart.entry[0] = interpts[i-1].x;
			Pstart.entry[1] = interpts[i-1].y;

			Pend.entry[0] = interpts[i+1].x;
			Pend.entry[1] = interpts[i+1].y;
		}

		T[i] = -GetTangent(Pstart, Pend);
	}

	////Then, we calculate the points on the curve
	num_output = 0;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < step; j++)
		{
			output[num_output].x = interpts[i].x * Base[j][0] + interpts[(i+1)%n].x * Base[j][1]\
				+ T[i].entry[0] * Base[j][2] + T[(i+1)%n].entry[0] * Base[j][3];
			
			output[num_output].y = interpts[i].y * Base[j][0] + interpts[(i+1)%n].y * Base[j][1]\
				+ T[i].entry[1] * Base[j][2] + T[(i+1)%n].entry[1] * Base[j][3];

			num_output++;
		}
	}

	free(Base);

}


void Hermitecurve_open(int n, ctr_point *interpts, icVector2 *T, ctr_point *output, int step, int &num_output)
{
	int i, j;
	Base4 *Base = (Base4 *) malloc(sizeof(Base4) * step);
	double s/*, h1, h2, h3, h4*/;

	//icVector2 *T = (icVector2 *) malloc(sizeof(icVector2) * n);
	icVector2 Pstart, Pend, T1, T2;

	////calculate the blending functions once and store them into an array
	for(j = 0; j < step; j++)
	{
		s = (double) j / (double) step;
		Base[j][0] = 2 * s * s * s - 3 * s * s + 1;
		Base[j][1] = -2 * s * s * s + 3 * s * s;
		Base[j][2] = s * s * s - 2 * s * s + s;
		Base[j][3] = s * s * s - s * s;
	}

	////We may first store the tagent vectors on the interpolate points into an array
	for(i = 0; i < n; i++)
	{
		if( i == 0)
		{
			Pstart.entry[0] = interpts[0].x;
			Pstart.entry[1] = interpts[0].y;

			Pend.entry[0] = interpts[1].x;
			Pend.entry[1] = interpts[1].y;
		}
		else if( i == n-1 )
		{
			Pstart.entry[0] = interpts[n-2].x;
			Pstart.entry[1] = interpts[n-2].y;

			Pend.entry[0] = interpts[n-1].x;
			Pend.entry[1] = interpts[n-1].y;
		}
		else{
			Pstart.entry[0] = interpts[i-1].x;
			Pstart.entry[1] = interpts[i-1].y;

			Pend.entry[0] = interpts[i+1].x;
			Pend.entry[1] = interpts[i+1].y;
		}

		T[i] = -GetTangent(Pstart, Pend);
	}

	////Then, we calculate the points on the curve
	num_output = 0;
	for(i = 0; i < n-1; i++)
	{
		for(j = 0; j < step; j++)
		{
			output[num_output].x = interpts[i].x * Base[j][0] + interpts[(i+1)%n].x * Base[j][1]\
				+ T[i].entry[0] * Base[j][2] + T[(i+1)%n].entry[0] * Base[j][3];
			
			output[num_output].y = interpts[i].y * Base[j][0] + interpts[(i+1)%n].y * Base[j][1]\
				+ T[i].entry[1] * Base[j][2] + T[(i+1)%n].entry[1] * Base[j][3];

			num_output++;
		}
	}
	free(Base);
}




icVector2 GetTangent(icVector2 Pend, icVector2 Pstart)
{
	icVector2 result;
	result = Pend - Pstart;

	result = 0.5 * result;
	return result;
}