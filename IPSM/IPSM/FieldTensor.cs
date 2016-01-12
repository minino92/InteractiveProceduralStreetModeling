using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace IPSM
{
    public class FieldTensor
    {
        public System.Drawing.Point position;
        public System.Drawing.Point finalPosition;
        public Vector[,] tensorField;
        public FieldTensor()
        {
            tensorField = new Vector[Noise.size, Noise.size];
        }
        public void generateBaseField()
        {
            
        }
    }
}
