using System;
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
        RoadType type;
        List<PointF> segments;
        int nodeID1;
        int nodeID2;
    }
    public class Node
    {
        PointF position;
        List<int> connectedNodeIDs;
        List<int> connectedRoadIDs;
    }
    class StreetGraph
    {
        public StreetGraph(PointF bottomLeft,PointF topRight,TensorField tf,float distSeparation)
        {

        }
        public void createRandomSeedList(int numberOfSeeds, bool append)
        {

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

        }
    }
}
