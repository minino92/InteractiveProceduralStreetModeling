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
        public double angle;
        public FieldTensor()
        {
        }
        public void CalculateAngle()
        {
            Vector direction = new Vector(finalPosition.X - position.X,finalPosition.Y - position.Y);
            direction.Normalize();
            Vector initDirection = new Vector(-1, 0);
            initDirection.Normalize();
            angle = Math.Acos(Vector.Multiply(initDirection, direction))*180/Math.PI -Math.PI/2;
            if (direction.Y > 0)
            {
                angle = 360-angle;
            }
        }
    }
}
