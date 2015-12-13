using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Forms;
using System.Threading;
using System.Threading.Tasks;

namespace IPSM
{
    public partial class IPSM : Form
    {
        Thread th = null;
        public IPSM()
        {
            InitializeComponent();
        }

        private void IPSM_Load(object sender, EventArgs e)
        {

        }

        private void IPSM_Paint(object sender, PaintEventArgs e)
        {
            // Define the size of our viewport using arbitary world coordinates
            var viewportSize = new SizeF(10, 10);
            th = new Thread(() =>
            {
                // Create a new bitmap image that is 500 by 500 pixels
                using (var bmp = new Bitmap(pictureZone.Width, pictureZone.Height, PixelFormat.Format32bppPArgb))
                {
                    // Create graphics object to draw on the bitmap
                    using (var g = pictureZone.CreateGraphics())
                    {
                        // Set up transformation so that drawing calls automatically convert world coordinates into bitmap coordinates
                        //g.TranslateTransform(0, bmp.Height * 0.5f - 1);
                        //g.ScaleTransform(bmp.Width / viewportSize.Width, -bmp.Height / viewportSize.Height);
                        //g.TranslateTransform(0, -viewportSize.Height * 0.5f);

                        // Create pen object for drawing with
                        using (var redPen = new Pen(Color.Red)) // Note that line thickness is in world coordinates!
                        {
                            // Randomization
                            var rand = new Random();

                            // Draw a 10x10 grid of vectors
                            var a = new Vector();
                            for (a.X = 0.5; a.X < bmp.Width; a.X += 1.0)
                            {
                                for (a.Y = 0.5; a.Y < bmp.Height; a.Y += 1.0)
                                {
                                    // Connect the center of this cell to a random point inside the cell
                                    var offset = new Vector(rand.NextDouble() - 0.5, rand.NextDouble() - 0.5);
                                    var b = a + offset;
                                    g.DrawLine(redPen, new PointF((float)a.X, (float)a.Y), new PointF((float)b.X, (float)b.Y));
                                }
                            }
                        }
                    }
                    //end graphics

                }
                //end bitmap
            });
            th.Start();             
        }

        private void IPSM_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (th != null)
            {
                th.Abort();
                th = null;
            }
        }        
    }
}
