using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;

namespace IPSM
{
    enum RoadType
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
        private List<PointF> mSeeds;
        private PointF mBottomLeft;
        private PointF mTopRight;
        private int mLastNodeID;
        private int mLastRoadID;
        //Distance for road density
        private float mDistSeparation;
        private int mSeedInitMethod;
        private Size mRegionSize;
        
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
        }
        public void createRandomSeedList(int numberOfSeeds, bool append)
        {
            if (!append)
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
        public void computeMajorHyperstreamlines(bool clearStorage)
        {

        }
        public void computeStreetGraph(bool clearStorage)
        {
            if (clearStorage) clearStoredStreetGraph();
            if (mtf == null) throw new NullReferenceException("TensorField reference null");
            generateSeedListWithUIMethod();
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
    }
}
