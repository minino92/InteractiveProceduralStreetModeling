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
        private string path = @"Datas\cercle_5.txt";
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
                    string pattern = @"((?<x>-?\d),(?<y>-?\d))";
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
                            }
                            else
                                map.SetPixel((int)position[0],(int) position[1], Color.Black);
                            //points.Add(new PointF(position[0], position[1]));
                            if (float.Parse(ms.Groups["x"].Value) >= 0 && float.Parse(ms.Groups["y"].Value) >= 0)
                            {
                                points_octants.Add(new PointF(float.Parse(ms.Groups["x"].Value), float.Parse(ms.Groups["y"].Value)));
                            }      
                            points.Add(new PointF(float.Parse(ms.Groups["x"].Value), float.Parse(ms.Groups["y"].Value)));
                        }
                        
                    }
                    //segmentNaiveLine(new PointF(0, 0), new PointF(25, 25), 3, 8,-4, map);
                    premierOctant(points_octants.ToArray(),map);
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
            float[] u =new float[]{10,0,0};
            float[] v=new float[]{0,-10,0};
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
            bool segment = false;
            int r;
            while (i < E.Length)
            {
                m = E[i];
                r =(int) (a * m.X - b * m.Y);
                segment = r >= mu - 1 && r <= mu+b;
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
                    float[] position = ChangeBase(p.X, p.Y);
                    map.SetPixel((int)position[0],(int) position[1], Color.Black);
                }
                i++;
            }
        }
    }
}
