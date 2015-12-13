/*  AddNoise.cpp */


#include "stdafx.h"
//#include "RegionSmoothing.h"
#include "VFDataStructure.h"
#include "LocalTracing.h"
#include "tensoranalysis.h"
#include "evenlystreamlines.h"

#include "LocalRegionEdit.h"
#include "regionsmooth_quad.h"
#include "SketchDesign.h"

#include "Numerical.h"

#include "PerlinNoise.h"

extern StreetNet *streetnet;
extern QuadMesh *quadmesh;

extern double majorDensity ;
extern double minorDensity ;

extern int ninnercells, *innercells;

extern int *region_quadverts;                ////mesh vertices inside user selected region
extern int nregion_quadverts;

extern EvenStreamlinePlace *major_level1 ;
extern EvenStreamlinePlace *minor_level1 ;


//void GrayImage::GenerateNoiseTexture( Uint32 width, Float32 scale ) {//, bool useSum ) {
//    Reset();
//
//    mWidth = width;
//    mHeight = width;
//    mBuffer = new Float32[ width * width ];
//    ImprovedNoise noise;
//
//    Uint32 cntx, cnty;
//    Float32 posX, posY;
//
//    Float32 max = 0, min = 0;
//    Float32 curAmp, curScale;
//
//    //if( useSum ) {
//        for( Uint32 cnty = 0; cnty < width; ++cnty ) {
//            for( Uint32 cntx = 0; cntx < width; ++cntx ) {
//
//                curAmp = 0.5f;
//                //curAmp = amp;
//                curScale = scale * 0.5f;
//
//                mBuffer[ cntx + cnty * width ] = 0;
//                posX = static_cast< Float32 >( cntx ) / width;
//                posY = static_cast< Float32 >( cnty ) / width;
//                for( Uint32 band = 0; band < 4; ++band ) {
//                    mBuffer[ cntx + cnty * width ] += curAmp * noise.noise( curScale * 
//posX, curScale * posY, 11.5 );
//                    curScale *= 2.0f;
//                    curAmp *= 0.5f;
//                }
//
//                if( mBuffer[ cntx + cnty * width ] > max ) {
//                    max = mBuffer[ cntx + cnty * width ];
//                }
//
//                if( mBuffer[ cntx + cnty * width ] < min ) {
//                    min = mBuffer[ cntx + cnty * width ];
//                }
//            }
//        }
//    //}
//    //else {
//    //    for( Uint32 cnty = 0; cnty < width; ++cnty ) {
//    //        for( Uint32 cntx = 0; cntx < width; ++cntx ) {
//    //            posX = static_cast< Float32 >( cntx ) / width;
//    //            posY = static_cast< Float32 >( cnty ) / width;
//    //            mBuffer[ cntx + cnty * width ] = noise.noise( scale * posX, scale * 
//posY, 11.5 );
//    //
//    //            if( mBuffer[ cntx + cnty * width ] > max ) {
//    //                max = mBuffer[ cntx + cnty * width ];
//    //            }
//
//    //            if( mBuffer[ cntx + cnty * width ] < min ) {
//    //                min = mBuffer[ cntx + cnty * width ];
//    //            }
//    //        }
//    //    }
//    //}
//
//    // scale all values to be between 0 and 1
//    for( cnty = 0; cnty < width; ++cnty ) {
//        for( cntx = 0; cntx < width; ++cntx ) {
//            mBuffer[ cntx + cnty * width ] = ( mBuffer[ cntx + cnty * width ] - min ) / ( 
//max - min );
//        }
//    }
//
//    LoadTextureToVideoMemory( false );
//}



//function Noise1(integer x, integer y)
//    n = x + y * 57
//    n = (n<<13) ^ n;
//    return ( 1.0 - ( (n * (n * n * 15731 + 789221) + 1376312589) & 7fffffff) / 1073741824.0);    
//  end function
//
//  function SmoothNoise_1(float x, float y)
//    corners = ( Noise(x-1, y-1)+Noise(x+1, y-1)+Noise(x-1, y+1)+Noise(x+1, y+1) ) / 16
//    sides   = ( Noise(x-1, y)  +Noise(x+1, y)  +Noise(x, y-1)  +Noise(x, y+1) ) /  8
//    center  =  Noise(x, y) / 4
//    return corners + sides + center
//  end function
//
//  function InterpolatedNoise_1(float x, float y)
//
//      integer_X    = int(x)
//      fractional_X = x - integer_X
//
//      integer_Y    = int(y)
//      fractional_Y = y - integer_Y
//
//      v1 = SmoothedNoise1(integer_X,     integer_Y)
//      v2 = SmoothedNoise1(integer_X + 1, integer_Y)
//      v3 = SmoothedNoise1(integer_X,     integer_Y + 1)
//      v4 = SmoothedNoise1(integer_X + 1, integer_Y + 1)
//
//      i1 = Interpolate(v1 , v2 , fractional_X)
//      i2 = Interpolate(v3 , v4 , fractional_X)
//
//      return Interpolate(i1 , i2 , fractional_Y)
//
//  end function
//
//
//  function PerlinNoise_2D(float x, float y)
//
//      total = 0
//      p = persistence
//      n = Number_Of_Octaves - 1
//
//      loop i from 0 to n
//
//          frequency = 2i
//          amplitude = pi
//
//          total = total + InterpolatedNoisei(x * frequency, y * frequency) * amplitude
//
//      end of i loop
//
//      return total
//
//  end function



ImprovedNoise noise;

/*  jitter all the obtained roads  */
void jitter_all_roads(double scale, double curAmp)
{
	int i, j/*, k*/;

    curAmp = 0.01f;
    //curAmp = amp;
	double curScale;
	double posX, posY;

	if(streetnet==NULL) return;

	icVector2 ran_dir, orient, jitter_dir;

	StreetGraphEdge *edge;
	double ran_val=0;
	double dis_accum=0;

	for(i=0;i<streetnet->edgelist->nedges;i++)
	{
		edge=streetnet->edgelist->edges[i];
		dis_accum=0.0;
		ran_dir.set(0.);
		double sub_scale=.8;

		for(j=1;j<edge->ninter_pts-1;j++)
		{
			ran_val=0;
			curScale = scale * 0.5;
			//curAmp = 0.1;
			curAmp = 0.5;

			orient.entry[0]=edge->inter_pts[j+1]->x-edge->inter_pts[j]->x;
			orient.entry[1]=edge->inter_pts[j+1]->y-edge->inter_pts[j]->y;

			dis_accum+=length(orient);

			if(j!=1 && j!= edge->ninter_pts-1 && dis_accum<quadmesh->xinterval/2.5)
			{
				ran_dir = sub_scale * ran_dir;
				sub_scale *=.8;
				edge->inter_pts[j]->x += sub_scale*jitter_dir.entry[0];
				edge->inter_pts[j]->y += sub_scale*jitter_dir.entry[1];
				edge->inter_pts[j]->cellid=get_cellID_givencoords(edge->inter_pts[j]->x,
					edge->inter_pts[j]->y);
				continue;
			}

			dis_accum=0;
			sub_scale=.8;

			posX=edge->inter_pts[j]->x;
			posY=edge->inter_pts[j]->y;

			//for(k=0;k<4;k++)
			//{
   //             ran_val += curAmp * noise.noise(curScale*posX, curScale*posY, 11.5);
   //             curScale *= 2.0;
   //             curAmp *= 0.5;
			//}
            
			ran_val = curAmp * noise.noise(curScale*posX, curScale*posY, 11.5);

			if(fabs(ran_val) > (quadmesh->xinterval/6./**majorDensity/2.*/))
			{
				if(ran_val>0)
					ran_val=quadmesh->xinterval/6./**majorDensity/2.*/;
				else
					ran_val=-quadmesh->xinterval/6./**majorDensity/2.*/;
			}

			/*   generate a random direction   */
			ran_dir.entry[0]=((double)rand()/RAND_MAX);
			ran_dir.entry[1]=((double)rand()/RAND_MAX);

			normalize(ran_dir);
			ran_dir= ran_val*ran_dir;

			/*   project to the perpendicular direction of the road  */

			orient.entry[0]=edge->inter_pts[j+1]->x-edge->inter_pts[j-1]->x;
			orient.entry[1]=edge->inter_pts[j+1]->y-edge->inter_pts[j-1]->y;

			normalize(orient);

			jitter_dir=ran_dir-dot(ran_dir,orient)*orient;

			edge->inter_pts[j]->x += jitter_dir.entry[0];
			edge->inter_pts[j]->y += jitter_dir.entry[1];
			edge->inter_pts[j]->cellid=get_cellID_givencoords(edge->inter_pts[j]->x,
				edge->inter_pts[j]->y);
		}
	}
}


/*  jitter only the major roads   */
void jitter_all_majRoads(double scale, double curAmp, bool majormin)
{
	int i, j/*, k*/;

    curAmp = 0.01f;
    //curAmp = amp;
	double curScale;
	double posX, posY;
	EvenStreamlinePlace *curPlace;

	if(!majormin)
		curPlace=major_level1;
	else
		curPlace=minor_level1;

	if(curPlace==NULL)
	{
		return;
	}


	icVector2 ran_dir, orient, jitter_dir;

	//StreetGraphEdge *edge;
	Trajectory *curtraj;
	LineSeg *curline;
	double ran_val=0;
	double dis_accum=0;

	for(i=0;i<curPlace->evenstreamlines->ntrajs;i++)
	{
		curtraj=curPlace->evenstreamlines->trajs[i];

		dis_accum=0.0;
		ran_dir.set(0.);
		double sub_scale=.8;

		for(j=0;j<curtraj->nlinesegs;j++)
		{
			curline = &curtraj->linesegs[j];
			ran_val=0;
			curScale = scale * 0.5;
			
			curAmp = 0.5;

			//orient.entry[0]=curline->gend[0]-curline->gstart[0];
			//orient.entry[1]=curline->gend[1]-curline->gstart[1];

			dis_accum+=curline->length;

			int haha=rand()%3+1;

			if(j!=0 && j!= curtraj->nlinesegs-1 && dis_accum</*2*/haha*quadmesh->xinterval)
			{
				ran_dir = sub_scale * ran_dir;
				sub_scale *=.98;
				curline->gstart[0] += sub_scale*jitter_dir.entry[0];
				curline->gstart[1] += sub_scale*jitter_dir.entry[1];
				curline->Triangle_ID=get_cellID_givencoords(curline->gstart[0],
					curline->gstart[1]);
				if(j>0)
				{
					curtraj->linesegs[j-1].gend[0]=curline->gstart[0];
					curtraj->linesegs[j-1].gend[1]=curline->gstart[1];
				}
				continue;
			}

			dis_accum=0;
			sub_scale=.98;

			posX=curline->gstart[0];
			posY=curline->gstart[1];

          
			ran_val = curAmp * noise.noise(curScale*posX, curScale*posY, 11.5);

			if(fabs(ran_val) > (quadmesh->xinterval/6.))
			{
				if(ran_val>0)
					ran_val=quadmesh->xinterval/6.;
				else
					ran_val=-quadmesh->xinterval/6.;
			}

			/*   generate a random direction   */
			ran_dir.entry[0]=((double)rand()/RAND_MAX);
			ran_dir.entry[1]=((double)rand()/RAND_MAX);

			normalize(ran_dir);
			ran_dir= ran_val*ran_dir;

			/*   project to the perpendicular direction of the road  */

			orient.entry[0]=curline->gend[0]-curline->gstart[0];
			orient.entry[1]=curline->gend[1]-curline->gstart[1];
			//orient.entry[0]=curtraj->linesegs[j+1].gend[0]-curtraj->linesegs[j-1].gstart[0];
			//orient.entry[1]=curtraj->linesegs[j+1].gend[1]-curtraj->linesegs[j-1].gstart[1];

			normalize(orient);

			jitter_dir=ran_dir-dot(ran_dir,orient)*orient;

			curline->gstart[0] += jitter_dir.entry[0];
			curline->gstart[1] += jitter_dir.entry[1];
			curline->Triangle_ID=get_cellID_givencoords(curline->gstart[0],
				curline->gstart[0]);

			if(j>0)
			{
				curtraj->linesegs[j-1].gend[0]=curline->gstart[0];
				curtraj->linesegs[j-1].gend[1]=curline->gstart[1];
			}
		}
	}
}



/*  jitter the roads in a specific region  */
void jitter_roads_inReg(double curScale)
{

	convert_to_tensorline();
	if(ninnercells==0 || innercells==NULL) return;

	if(streetnet==NULL) return;

	/*  we just need to find out those edges inside the region and jitter them  */
	int i, j;
	StreetGraphEdge *edge;
	icVector2 ran_dir, orient, jitter_dir;
	double posX, posY;

	double ran_val=0;
	double dis_accum=0;

	double con_interval=10./180.*M_PI;

	double curAmp;

	for(i=0;i<streetnet->edgelist->nedges;i++)
	{
		edge=streetnet->edgelist->edges[i];
		edge->visited=false;
	}

	for(i=0;i<ninnercells;i++)
	{
		QuadCell *face=quadmesh->quadcells[innercells[i]];

		if(face->nstreetgraphedges==0 || face->streetgraphedgelist==NULL) 
			continue;

		for(j=0;j<face->nstreetgraphedges;j++)
		{
			edge=streetnet->edgelist->edges[face->streetgraphedgelist[j]];

			dis_accum=0.0;
			ran_dir.set(0.);
			con_interval=10./180.*M_PI;
			double sub_scale=sin(M_PI/2.-con_interval);

			if(edge->visited)
				continue;
			edge->visited=true;

			for(int k=1; k<edge->ninter_pts-1; k++)
			{
				/*  jitter the edge line segment using Perlin Noise  */
				ran_val=0;
				//curAmp = 0.1;
				curAmp = 0.5;

				orient.entry[0]=edge->inter_pts[k+1]->x-edge->inter_pts[k]->x;
				orient.entry[1]=edge->inter_pts[k+1]->y-edge->inter_pts[k]->y;

				dis_accum+=length(orient);

				if(!is_inregion(edge->inter_pts[k]->x, edge->inter_pts[k]->y))
					continue;

				if(k!=1 && k!= edge->ninter_pts-1 && dis_accum<quadmesh->xinterval/2.5)
				{
					ran_dir = sub_scale * ran_dir;
					con_interval+=sin(M_PI/2.-con_interval);
					sub_scale = sin(M_PI/2.-con_interval);;
					edge->inter_pts[k]->x += sub_scale*jitter_dir.entry[0];
					edge->inter_pts[k]->y += sub_scale*jitter_dir.entry[1];
					edge->inter_pts[k]->cellid=get_cellID_givencoords(edge->inter_pts[k]->x,
						edge->inter_pts[k]->y);
					continue;
				}

				dis_accum=0;
				sub_scale=.8;

				posX=edge->inter_pts[k]->x;
				posY=edge->inter_pts[k]->y;
	            
				ran_val = curAmp * noise.noise(curScale*posX, curScale*posY, 11.5);

				if(fabs(ran_val) > (quadmesh->xinterval/6./**majorDensity/2.*/))
				{
					if(ran_val>0)
						ran_val=quadmesh->xinterval/6./**majorDensity/2.*/;
					else
						ran_val=-quadmesh->xinterval/6./**majorDensity/2.*/;
				}

				/*   generate a random direction   */
				ran_dir.entry[0]=((double)rand()/RAND_MAX);
				ran_dir.entry[1]=((double)rand()/RAND_MAX);

				normalize(ran_dir);
				ran_dir= ran_val*ran_dir;

				/*   project to the perpendicular direction of the road  */

				orient.entry[0]=edge->inter_pts[k+1]->x-edge->inter_pts[k-1]->x;
				orient.entry[1]=edge->inter_pts[k+1]->y-edge->inter_pts[k-1]->y;

				normalize(orient);

				jitter_dir=ran_dir-dot(ran_dir,orient)*orient;

				edge->inter_pts[k]->x += jitter_dir.entry[0];
				edge->inter_pts[k]->y += jitter_dir.entry[1];
				edge->inter_pts[k]->cellid=get_cellID_givencoords(edge->inter_pts[k]->x,
					edge->inter_pts[k]->y);
			}
		}
	}
}

/*
    We jitter the tensor field everywhere
*/

void rot_ten(icMatrix2x2 in, icMatrix2x2 &out, double ang)
{
	double cosang=cos(ang);
	double sinang=sin(ang);

	icMatrix2x2 rot;
	double M[2][2]={{cosang, -sinang}, {sinang, cosang}};
	//rot.set(M);

	out.setIdentity();
	out.set(M);
	
	M[0][1]=sinang; M[1][0]=-sinang;
	rot.set(M);

	out.rightMultiply(in);
	out.rightMultiply(rot);

}

void jitter_ten_one_vert(int vertid, double curScale, double curAmp)
{
	QuadVertex *v=quadmesh->quad_verts[vertid];
	double posX, posY;

	double ran_val=0;
	double dis_accum=0;

	double con_interval=10./180.*M_PI;

	//double curAmp=.5;

	/* save the old tensor */
	v->origin_ten=v->Jacobian;
	
	/* jitter the tensor at the vertex  */
	posX=v->x;
	posY=v->y;
	for(int i=0;i<4;i++)
	{
		curScale *= 2.0;
		curAmp *= .5;
		ran_val += curAmp * noise.noise(curScale*posX, curScale*posY, 11.5);
	}

	/*  rotate the tensor using ran_val as the angle  */
	rot_ten(v->origin_ten, v->Jacobian, ran_val * M_PI);
}

void jitter_all_tens(double scale, double curAmp)
{
	int i;
	for(i=0;i<quadmesh->nverts;i++)
	{
		jitter_ten_one_vert(i, scale, curAmp);
	}
}

/*
    We jitter the tensor field inside a sub-region
*/
void jitter_tens_inReg(double curScale, double curAmp)
{
	int i;
	for(i=0;i<nregion_quadverts; i++)
	{
		jitter_ten_one_vert(region_quadverts[i], curScale, curAmp);
	}
}