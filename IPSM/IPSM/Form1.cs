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
        Noise noise;
        public IPSM()
        {
            InitializeComponent();
            noise = new Noise();

        }

        private void IPSM_Load(object sender, EventArgs e)
        {

        }

        private void IPSM_Paint(object sender, PaintEventArgs e)
        {
            pictureZone.Width = Noise.size;
            pictureZone.Height = Noise.size;
            th = new Thread(() =>
            {
                var g = pictureZone.CreateGraphics();
                using (var bmp = new Bitmap(Noise.size, Noise.size, g)) 
                {                    
                    using (var redPen = new Pen(Color.Red)) 
                    {
                        for (int i = 0; i < Noise.size; i++ )
                        {
                            for (int j = 0; j < Noise.size; j++ )
                            {
                                
                                bmp.SetPixel(i, j, Color.FromArgb(noise.noiseTable[i,j],noise.noiseTable[i,j],noise.noiseTable[i,j]));  
                            }
                        }
                        g.DrawImage(bmp, new PointF(0, 0));
                    }
                }
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
