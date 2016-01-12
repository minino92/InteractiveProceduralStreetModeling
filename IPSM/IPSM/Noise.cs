using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IPSM
{
    class Noise
    {

        public int[,] noiseTable;
        public static int size = 800;

        public Noise()
        {
            noiseTable = new int[size, size];
            Random random = new Random();

            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    noiseTable[i, j] = random.Next(0,255);
                }
            }
        }

    }
}
