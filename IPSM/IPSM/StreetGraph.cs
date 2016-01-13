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
        private List<PointF> mSeeds;
        private PointF mBottomLeft;
        private PointF mTopRight;
        private int mLastNodeID;
        private int mLastRoadID;
        //Distance for road density
        private float mDistSeparation;
        private int mSeedInitMethod;
        private Size mRegionSize;
        public int Scale;
        public StreetGraph(PointF bottomLeft,PointF topRight,TensorField tf,float distSeparation)
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
        public void computeMajorHyperstreamlines(Bitmap bmp,Graphics g)
        {
            createRandomSeedList(40,false);
            //int scaleI = Convert.ToInt32(Noise.size / mtf.NumberOfTensorsToDisplay);//size between tensors to display
            //int scaleJ = scaleI;
            PointF currentSeed = mSeeds[0];
            for (int i = 0; i < mSeeds.Count; i++)
            {
                try
                {
                    //Major vectors
                    EigenVector ev = mtf.matrixEigenVectors[Convert.ToInt32(Math.Floor(mSeeds[i].X)), Convert.ToInt32(Math.Floor(mSeeds[i].Y))];
                        //first direction
                    PointF fin = drawMajor(mSeeds[i], ev, mSeeds[i]);
                    g.DrawRectangle(new Pen(Color.Red), mSeeds[i].X, mSeeds[i].Y, 2f, 2f);
                    g.DrawLine(new Pen(Color.Red), mSeeds[i], fin);
                        //second direction
                    fin = drawMajor(mSeeds[i], ev, mSeeds[i],false);
                    g.DrawLine(new Pen(Color.Gold), mSeeds[i], fin);
                    //-----------------------------
                    //Minor vectors
                        //first direction
                    fin = drawMinor(mSeeds[i], ev, mSeeds[i]);
                    g.DrawLine(new Pen(Color.Black), mSeeds[i], fin);
                        //second direction
                    fin = drawMinor(mSeeds[i], ev, mSeeds[i],false);
                    g.DrawLine(new Pen(Color.Gold), mSeeds[i], fin);
                }
                catch (Exception e)
                {
                    //there are some errors that we have to handle but we do not
                    //that is why sometime the number of the seeds points are not correct.
                }                
            }
        }
        public void computeStreetGraph(bool clearStorage)
        {
            if (clearStorage) clearStoredStreetGraph();
            if (mtf == null) throw new NullReferenceException("TensorField reference null");
            generateSeedListWithUIMethod();//here we create the seeds
            bool majorGrowth = true;
            for (int k = 0; k < mSeeds.Count; k++)
            {
                Node node1 = mNodes[++mLastNodeID];
                node1.position = mSeeds[k];

                Road road = mRoads[++mLastRoadID];
                node1.connectedRoadIDs.Add(mLastRoadID);
                road.type = RoadType.Principal;
                road.nodeID1 = mLastNodeID;

                Road road2 = mRoads[++mLastRoadID];
                node1.connectedRoadIDs.Add(mLastRoadID);
                road2.type = RoadType.Principal;
                road2.nodeID1 = mLastNodeID;

                growRoad(road, node1, majorGrowth, false, false);
                growRoad(road2, node1, majorGrowth, true, false);

                majorGrowth = !majorGrowth;
            }
        }
        public void generateStreetGraph()
        {
            computeStreetGraph(true);
            drawStreetGraph(true, false);
        }

        public void drawStreetGraph(bool p1, bool p2)
        {
            throw new NotImplementedException();
        }
        public void clearStoredStreetGraph()
        {
            mNodes.Clear();
            mRoads.Clear();
            mLastNodeID = 0;
            mLastRoadID = 0;
        }
        public void generateSeedListWithUIMethod()
        {
            //there are 3 types which can be used
            //we only consider one case
            createRandomSeedList(100, false);
        }
        public void growRoad(Road road, Node startNode, bool growInMajorDirection, bool growInOppositeDirection, bool useExceedLenStopCond)
        {

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
    }
}
