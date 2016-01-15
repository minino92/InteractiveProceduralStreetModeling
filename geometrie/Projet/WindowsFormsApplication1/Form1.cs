using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;

namespace WindowsFormsApplication1
{
    public partial class Form1 : Form
    {
        private Graphics g;
        private float[] centreP;
        private string path = @"Datas\coeur_30.txt";
        private List<PointF> points;
        private List<PointF> points_octants;
        public Form1()
        {
            InitializeComponent();
            centreP = new float[] { Width / 2, Height / 2 };
            points = new List<PointF>();
            points_octants = new List<PointF>();
        }

        private void Form1_Paint(object sender, PaintEventArgs e)
        {
            try
            {
                using (g = CreateGraphics())
                {
                    Bitmap map = new Bitmap(Width,Height,g);
                    g.DrawPolygon(new Pen(new SolidBrush(Color.Red), 10), new PointF[] { new PointF(0, 0), new PointF(0, 0) });
                    StreamReader reader = new StreamReader(path);
                    string pattern = @"((?<x>-?\d+),(?<y>-?\d+))";
                    while (!reader.EndOfStream)
                    {
                        string line = reader.ReadLine();
                        Regex rg = new Regex(pattern);
                        MatchCollection m=rg.Matches(line);
                        foreach (Match ms in m)
                        {
                            float[] position = ChangeBase(float.Parse(ms.Groups["x"].Value), float.Parse(ms.Groups["y"].Value));
                            if (float.Parse(ms.Groups["x"].Value) >= 0 && float.Parse(ms.Groups["y"].Value) >= 0)
                            {
                                map.SetPixel((int)position[0], (int)position[1], Color.Red);
                                points_octants.Add(new PointF(float.Parse(ms.Groups["x"].Value), float.Parse(ms.Groups["y"].Value)));
                            }
                            else
                                map.SetPixel((int)position[0],(int) position[1], Color.Black);
                                
                            points.Add(new PointF(float.Parse(ms.Groups["x"].Value), float.Parse(ms.Groups["y"].Value)));
                        }
                        
                    }
                    segmentNaiveLine(new PointF(0, 0), new PointF(10, 10), 3, 5,0, map);
                    //premierOctant(points_octants.ToArray(),map);
                    g.DrawImage(map, new Point(0, 0));
                }
            }
            catch (Exception ex)
            {

            }
        }
        /// <summary>
        /// Un répère de l'écran (0,0) top left
        /// Un répère de l'image (width/2,height/2)
        /// </summary>
        /// <param name="x">Position x dans le répère P</param>
        /// <param name="y">Position y dans le répère P</param>
        /// <returns></returns>
        private float[] ChangeBase(float x,float y)
        {            
            float[] res = new float[2];
            //u=(5,0,0) et v=(0,-5,0) => 1px -> x unité
            float[] u =new float[]{5,0,0};
            float[] v=new float[]{0,-5,0};
            float[] w=new float[]{0,0,0};//car on est en 2D
            res[0] = u[0] * x + v[0] * y + centreP[0];
            res[1] = u[1] * x + v[1] * y + centreP[1];
            return res;
        }
        /// <summary>
        /// A segment of a naive line
        /// </summary>
        private void segmentNaiveLine(PointF x0,PointF x1,int a,int b,int mu,Bitmap m)
        {
            PointF current = x0;
            int r = (int)(a * x0.X - b * x0.Y -mu);
            while (current.X < x1.X)
            {
                if (r >= (b - a))
                {
                    current.Y += 1;
                    r = r - b + a;
                }
                else
                {
                    r += a;
                }
                current.X += 1;
                float[] position = ChangeBase(current.X, (int)current.Y);
                //difference leaning points
                int xx = (int)position[0];
                int yy = (int)position[1];
                //points d'appui cf thèse partie 3.2
                if (a * current.X - b * current.Y == mu)
                {
                    m.SetPixel(xx, yy, Color.Black);
                }
                else if (a * current.X - b * current.Y == (mu + b - 1))
                {
                    m.SetPixel(xx, yy, Color.Blue);
                }
                else m.SetPixel(xx, yy, Color.Red);
            }
        }
        private void premierOctant(PointF[] E,Bitmap map)
        {
            //premier octant
            int a = 0, b = 1, mu = 0,i=1;
            PointF m = E[0];//1er points de E
            PointF s0 = m, i0 = m;//points d'appui initial sup et inf
            PointF p = new PointF();
            PointF dPoint = new PointF();
            bool segment = false;
            int r;
            while (i < E.Length)
            {
                m = E[i];
                //r =(int) (a * m.X - b * m.Y);
                //segment = r >= mu - 1 && r <= mu+b;
                //notion de distance vient de l'algo en anglais
                dPoint = new PointF(E[i].X - E[i - 1].X, E[i].Y - E[i - 1].Y);
                m.X += dPoint.X;
                m.Y += dPoint.Y;
                r = (int)(a * m.X - b * m.Y);
                if (r == mu - 1 || r == mu + b)
                {
                    if (r == mu - 1)
                    {
                        p.X = s0.X;
                        p.Y = s0.Y;
                        i0.X = m.X - b;
                        i0.Y = m.Y - a - 1;
                    }
                    else if (r == mu + b)
                    {
                        p.X = i0.X;
                        p.Y = i0.Y;
                        s0.X = m.X - b;
                        s0.Y = m.Y - a + 1;
                    }
                    b =(int)( m.X - p.X);
                    a=(int)(m.Y-p.Y);
                    mu =(int)( a * s0.X - b * s0.Y);
                    //segmentNaiveLine(s0, p, a, b, mu, map);
                    //float[] position = ChangeBase(p.X, p.Y);
                    //map.SetPixel((int)position[0],(int) position[1], Color.Black);
                }
                i++;
            }
        }
        void ExchangeAxis(int alpha, int beta, int gama, int delta, int phi, int psi)
        {
            int temp;
            temp = alpha; alpha = gama; gama = temp;
            temp = beta; beta = delta; delta = temp;
            temp = phi; phi = psi; psi = temp;
        }
        void ChangeVerticalAxis(int alpha, int beta, int gama, int delta, int phi, int psi)
        {
            if (delta != 0)
            {
                delta *= -1;
                psi *= -1;
            }
            else if (beta != 0)
            {
                beta *= -1;
                phi *= -1;
            }
        }
        void ChangeHorizontalAxis(int alpha, int beta, int gama, int delta, int phi, int psi)
        {
            if (alpha != 0)
            {
                alpha *= -1;
                phi *= -1;
            }
            else if (gama != 0)
            {
                gama *= -1;
                psi *= -1;
            }
        }
        void AdjustAxis(int dx, int dy, int ddx, int ddy,int alpha,int beta,int gama,int delta,int phi,int psi)
        {
            if (dy == 0 && ddx==dx && ddy==-dx)
            {
                ChangeVerticalAxis(alpha, beta, gama, delta, phi, psi);
            }
            else if (dx == 0 && ddx == dy && ddy == dy)
            {
                ChangeHorizontalAxis(alpha, beta, gama, delta, phi, psi);
            }
            else
            {
                if ((dx == dy && ddx == 0 && ddy == dy) || (dx == -dy && ddy == 0 && ddx == dy))
                {
                    ExchangeAxis(alpha, beta, gama, delta, phi, psi);
                }
            }
        }
        void Rotation(int xh, int yh, int dx, int dy, int alpha, int beta, int gama, int delta, int phi, int psi)
        {
            if(dx==1 && dy>=0)
            {
                alpha=1;beta=0;phi=-xh;
                gama=0;delta=1;psi=-yh;
            }
            else if(dx<=0 && dy==1)
            {
                alpha=0;beta=1;phi=-yh;
                gama=-1;delta=0;psi=xh;
            }
            else if (dx == -1 && dy < 0)
            {
                alpha = -1; beta = 0; phi = xh;
                gama = 0; delta = -1; psi = yh;
            }
            else if (dx >= 0 && dy == -1)
            {
                alpha = 0; beta = -1; phi = yh;
                gama = 1; delta = 0; psi = -xh;
            }
        }
        void MovingFrame(int k, int xk, int yk, int alpha, int beta, int gama, int delta, int phi, int psi, int a, int b, int mu, int ux, int uy, int lx, int ly)
        {
            int dx, dy, ddy, ddx;
        }
    }
}
