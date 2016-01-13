using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;

namespace IPSM
{
    public enum RoadType
    {
        Principal, Secondary
    }
    public class Road
    {
        public RoadType type;
        public List<PointF> segments;
        public int nodeID1;
        public int nodeID2;
    }
    public class Node
    {
        public PointF position;
        public List<int> connectedNodeIDs;
        public List<int> connectedRoadIDs;
    }
    class StreetGraph
    {
        private TensorField mtf;
        private Dictionary<int, Node> mNodes;
        private Dictionary<int, Road> mRoads;
        //the number of random points to start drawing hyperstreamlines
        public List<PointF> mSeeds;
        private PointF mBottomLeft;
        private PointF mTopRight;
        private int mLastNodeID;
        private int mLastRoadID;
        private List<PointF[]> totalRoads;
        //Distance for road density
        private double mDistSeparation;
        private int mSeedInitMethod;
        private Size mRegionSize;
        public int Scale;
        private int recursif = 0;
        private Dictionary<PointF, bool> positionToHandle;

        public StreetGraph(PointF bottomLeft,PointF topRight,TensorField tf,double distSeparation)
        {
            mtf = tf;
            mBottomLeft = bottomLeft;
            mTopRight = topRight;
            mDistSeparation = distSeparation;
            mLastNodeID = 0;
            mLastRoadID = 0;
            mSeedInitMethod = 0;
            mRegionSize = new Size(Noise.size,Noise.size);
            mSeeds = new List<PointF>();
            Scale = Convert.ToInt32(Noise.size / mtf.NumberOfTensorsToDisplay);
            totalRoads=new List<PointF[]>();
            positionToHandle = new Dictionary<PointF, bool>();
        }
        public void createRandomSeedList(int numberOfSeeds, bool append)
        {
            if (mSeeds!=null)
            {
                mSeeds.Clear();
            }
            Random rd = new Random();
            for (int i = 0; i < numberOfSeeds; i++)
            {
                float randX = rd.Next(0, mRegionSize.Width);
                float randY = rd.Next(0, mRegionSize.Height);
                mSeeds.Add(new PointF(mBottomLeft.X+randX,mBottomLeft.Y-randY));
            }
            
        }
        public void createDensityConstrainedSeedList(int numberOfSeeds, bool append)
        {
        }
        bool pointRespectSeedSeparationDistance(PointF point, float separationDistance)
        {
            return false;
        }
        public void computeMajorHyperStreamLinesWithDist(Bitmap b, Graphics g, float dist,PointF next)
        {
            if (recursif < mSeeds.Count)
            {
                try
                {
                    PointF currentPoint = next;
                    //Major vectors
                    EigenVector ev = mtf.matrixEigenVectors[Convert.ToInt32(Math.Floor(currentPoint.X)), Convert.ToInt32(Math.Floor(currentPoint.Y))];
                    //first direction
                    PointF fin = drawMajor(currentPoint, ev, currentPoint);
                    g.DrawRectangle(new Pen(Color.Red), currentPoint.X, currentPoint.Y, 2f, 2f);
                    g.DrawLine(new Pen(Color.Gold), currentPoint, fin);
                    //second direction
                    fin = drawMajor(currentPoint, ev, currentPoint, false);
                    g.DrawLine(new Pen(Color.Gold), currentPoint, fin);
                    //-----------------------------
                    //Minor vectors
                    //first direction
                    fin = drawMinor(currentPoint, ev, currentPoint);
                    g.DrawLine(new Pen(Color.Black), currentPoint, fin);
                    //second direction
                    fin = drawMinor(currentPoint, ev, currentPoint, false);
                    g.DrawLine(new Pen(Color.Black), currentPoint, fin);
                    recursif++;
                    PointF temp = currentPoint;
                    currentPoint = new PointF(temp.X + dist, temp.Y + dist);
                    computeMajorHyperStreamLinesWithDist(b, g, dist, currentPoint);
                }
                catch (Exception e)
                {

                }                
            }            
        }
        public void computeMajorHyperStreamLinesNew(Bitmap b, Graphics g,float dist,TensorField t)
        {
            drawSomething(mSeeds[0], dist, b,g,t.matrixTensors);
        }
        public void computeMajorHyperstreamlines(Bitmap bmp,Graphics g)
        {
            createRandomSeedList(10,false);
            //int scaleI = Convert.ToInt32(Noise.size / mtf.NumberOfTensorsToDisplay);//size between tensors to display
            //int scaleJ = scaleI;
            PointF currentSeed = mSeeds[0];
            PointF previousSeed = mSeeds[0];
            for (int i = 0; i < mSeeds.Count; i++)
            {
                try
                {
                    PointF[] temp = new PointF[5];
                    temp[0] = mSeeds[i];
                    //Major vectors
                    EigenVector ev = mtf.matrixEigenVectors[Convert.ToInt32(Math.Floor(mSeeds[i].X)), Convert.ToInt32(Math.Floor(mSeeds[i].Y))];
                        //first direction
                    PointF fin = drawMajor(mSeeds[i], ev, mSeeds[i]);
                    g.DrawRectangle(new Pen(Color.Red), mSeeds[i].X, mSeeds[i].Y, 2f, 2f);
                    g.DrawLine(new Pen(Color.Gold), mSeeds[i], fin);
                    temp[1] = fin;
                        //second direction
                    fin = drawMajor(mSeeds[i], ev, mSeeds[i],false);
                    g.DrawLine(new Pen(Color.Gold), mSeeds[i], fin);
                    temp[2] = fin;
                    //-----------------------------
                    //Minor vectors
                        //first direction
                    fin = drawMinor(mSeeds[i], ev, mSeeds[i]);
                    g.DrawLine(new Pen(Color.Black), mSeeds[i], fin);
                    temp[3] = fin;
                        //second direction
                    fin = drawMinor(mSeeds[i], ev, mSeeds[i],false);
                    g.DrawLine(new Pen(Color.Black), mSeeds[i], fin);
                    temp[4] = fin;
                    totalRoads.Add(temp);
                }
                catch (Exception e)
                {
                    //there are some errors that we have to handle but we do not
                    //that is why sometime the number of the seeds points are not correct.
                }                
            }
            //we have to sort every road now with the distance
        }
       
        private bool distance2SeedsOK(PointF a, PointF b,float dist)
        {
            System.Windows.Vector v = new System.Windows.Vector(b.X - a.X, b.Y - a.Y);
            if (Math.Sqrt(v.X * v.X + v.Y * v.Y) < dist) return false;
            return true;
        }
        private PointF drawMajor(PointF point,EigenVector prev,PointF prevP,bool other=true)
        {           
            if (point.X > Noise.size || point.X<0 || point.Y<0 || point.Y > Noise.size)
            {
                return prevP;
            }
            //Prendre Ev sur point
            EigenVector ev = mtf.matrixEigenVectors[Convert.ToInt32(Math.Floor(point.X)), Convert.ToInt32(Math.Floor(point.Y))];
            if ((prev.X != ev.X && prev.Y != ev.Y))
            {
                return point;
            }
            PointF temp = new PointF(point.X, point.Y);
            if (other)
            {
                point.X = temp.X + (float)ev.X * Scale;
                point.Y = temp.Y + (float)ev.Y * Scale;
            }
            else
            {
                point.X = temp.X + (float)ev.X * Scale*(-1);
                point.Y = temp.Y + (float)ev.Y * Scale*(-1);
            }
            prev = ev;
            prevP = temp;
            PointF fin=drawMajor(point, prev,prevP,other);
            return fin;
        }
        private PointF drawMinor(PointF point, EigenVector prev, PointF prevP,bool other=true)
        {
            if (point.X > Noise.size || point.X < 0 || point.Y < 0 || point.Y > Noise.size)
            {
                return prevP;
            }
            //Prendre Ev sur point
            EigenVector ev = mtf.matrixEigenVectors[Convert.ToInt32(Math.Floor(point.X)), Convert.ToInt32(Math.Floor(point.Y))];
            if ((prev.Z != ev.Z&& prev.W != ev.W))
            {
                return point;
            }
            PointF temp = new PointF(point.X, point.Y);
            if (other)
            {
                point.X = temp.X + (float)ev.Z * Scale;
                point.Y = temp.Y + (float)ev.W * Scale;
            }
            else
            {
                point.X = temp.X + (float)ev.Z * Scale*(-1);
                point.Y = temp.Y + (float)ev.W * Scale * (-1);
            }           
            prev = ev;
            prevP = temp;
            PointF fin = drawMinor(point, prev, prevP,other);
            return fin;
        }
        private void drawSomething(PointF position,float dist,Bitmap b,Graphics g,Tensor[,] mtx)
        {
            if ((position.X >= 0 && position.X < Noise.size) && (position.Y >= 0 && position.Y < Noise.size))
            {
                if (!positionToHandle.ContainsKey(position))
                {
                    try
                    {
                        positionToHandle.Add(position, true);
                        //Draw
                        PointF currentPoint = position;
                        //Major vectors
                        EigenVector ev = mtf.matrixEigenVectors[Convert.ToInt32(Math.Floor(currentPoint.X)), Convert.ToInt32(Math.Floor(currentPoint.Y))];
                        //first direction
                        PointF fin = drawMajor(currentPoint, ev, currentPoint);
                        g.DrawRectangle(new Pen(Color.Red), currentPoint.X, currentPoint.Y, 2f, 2f);
                        g.DrawLine(new Pen(Color.Gold), currentPoint, fin);
                        //second direction
                        fin = drawMajor(currentPoint, ev, currentPoint, false);
                        g.DrawLine(new Pen(Color.Gold), currentPoint, fin);
                        //-----------------------------
                        //Minor vectors
                        //first direction
                        fin = drawMinor(currentPoint, ev, currentPoint);
                        g.DrawLine(new Pen(Color.Black), currentPoint, fin);
                        //second direction
                        fin = drawMinor(currentPoint, ev, currentPoint, false);
                        g.DrawLine(new Pen(Color.Black), currentPoint, fin);
                        //call this function for each neigbour
                        Tensor ts = mtx[Convert.ToInt32(Math.Floor(position.X)), Convert.ToInt32(Math.Floor(position.Y))];
                        PointF[] nb = new PointF[4];
                        //nb[0] = new PointF(position.X + dist * (float)Math.Cos(ts.theta), position.Y + dist * (float)Math.Sin(ts.theta));
                        //nb[1] = new PointF(position.X + dist * (float)Math.Cos(ts.theta + Math.PI/2), position.Y + dist * (float)Math.Sin(ts.theta + Math.PI/2));
                        //nb[2] = new PointF(position.X + dist * (float)Math.Cos(ts.theta + Math.PI), position.Y + dist * (float)Math.Sin(ts.theta + Math.PI));
                        //nb[3] = new PointF(position.X + dist * (float)Math.Cos(ts.theta + 1.5*Math.PI), position.Y + dist * (float)Math.Sin(ts.theta + 1.5*Math.PI));
                        nb[0] = new PointF(position.X + dist * (float)Math.Cos(ts.theta), position.Y + dist * (float)Math.Sin(ts.theta));
                        nb[1] = new PointF(position.X + dist * (float)Math.Cos(ts.theta ), position.Y - dist * (float)Math.Sin(ts.theta ));
                        nb[2] = new PointF(position.X - dist * (float)Math.Cos(ts.theta ), position.Y + dist * (float)Math.Sin(ts.theta ));
                        nb[3] = new PointF(position.X - dist * (float)Math.Cos(ts.theta ), position.Y - dist * (float)Math.Sin(ts.theta ));
                        drawSomething(nb[0], dist, b, g, mtx);
                        drawSomething(nb[1], dist, b, g, mtx);
                        drawSomething(nb[2], dist, b, g, mtx);
                        drawSomething(nb[3], dist, b, g, mtx);
                    }
                    catch (Exception e) 
                    {
                    }                    
                }                
            }            
        }
    }
}
