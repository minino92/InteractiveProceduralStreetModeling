using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace IPSM
{
    class ibfv
    {

        int NPN;
        int NMESH = 100;
        //double DM;
        int NPIX = 512;
        double SCALE = 4.0;
        int iframe = 0;
        int Npat = 32;
        double alpha = (0.12 * 255);
        double sa;
        double tmax;
        double dmax;
        double M_PI = 3.14;
        byte[, ,] pat;

        public ibfv()
        {
            NPN = Noise.size;
            NMESH = Noise.size;
            //DM = (1.0 / (NMESH - 1.0));
            tmax = NPIX / (SCALE * NPN);
            dmax = SCALE / NPIX;
            pat = new byte[NPN, NPN, 4];
        }


        public void makePatterns()
        {
            int[] lut = new int[256];
            int[,] phase = new int[NPN, NPN];
             
            Random random = new Random();
            int i, j, k, t;
            for (i = 0; i < 256; i++)
            {
                lut[i] = i < 127 ? 0 : 255;
            }
            for (i = 0; i < NPN; i++)
            {
                for (j = 0; j < NPN; j++)
                {
                    phase[i, j] = random.Next(0, 255);
                };
            }
            for (k = 0; k < Npat; k++)
            {
                t = k * 256 / Npat;
                for (i = 0; i < NPN; i++)
                {
                    for (j = 0; j < NPN; j++)
                    {
                        pat[i, j, 0] =
                        pat[i, j, 1] =
                        pat[i, j, 2] = (byte)lut[(t + phase[i, j]) % 256];
                        pat[i, j, 3] = (byte)alpha;
                    }
                }
            }
        }

        /*void getDP(double x, double y, double px, double py)
        {
            double dx, dy, vx, vy, r;
            dx = x - 0.5;
            dy = y - 0.5;
            r = dx * dx + dy * dy;
            if (r < 0.0001) r = 0.0001;
            vx = sa * dx / r + 0.02;
            vy = sa * dy / r;
            r = vx * vx + vy * vy;
            if (r > dmax * dmax)
            {
                r = Math.Sqrt(r);
                vx *= dmax / r;
                vy *= dmax / r;
            }
            px = x + vx;
            py = y + vy;
        }*/

        public void display(Bitmap bmp, Graphics g, EigenVector[,] mEigenVector)
        {
            int i, j;
            int px, py;
            px = 0;
            py = 0;
            sa = 0.010 * Math.Cos(iframe * 2.0 * M_PI / 200.0);
            for (int k = 0; k < 5; k++)
            {
                for (i = 0; i < Noise.size; i++)
                {
                    for (j = 0; j < Noise.size; j++)
                    {
                        //getDP(i, j, px, py);
                        Vector dir = new Vector(mEigenVector[i, j].X, mEigenVector[i, j].Y) * 2;
                        //dir.Normalize();
                        px = (int)(i + dir.X);
                        py = (int)(j + dir.Y);
                        if (px < Noise.size && py < Noise.size && px >= 0 && py >= 0)
                        {
                            pat[px, py, 0] = (byte)((pat[i, j, 0] + pat[px, py, 0]) / 2);
                            pat[px, py, 1] = (byte)((pat[i, j, 1] + pat[px, py, 1]) / 2);
                            pat[px, py, 2] = (byte)((pat[i, j, 2] + pat[px, py, 2]) / 2);
                            bmp.SetPixel(i, j, Color.FromArgb(pat[i, j, 0], pat[i, j, 1], pat[i, j, 2]));
                        }
                    }
                }
            }
            iframe = iframe + 1;
        }

    }
}
