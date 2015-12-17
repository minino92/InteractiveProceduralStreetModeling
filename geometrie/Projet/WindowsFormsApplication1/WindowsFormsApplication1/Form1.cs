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
        private string path = @"Datas\coeur_5.txt";
        public Form1()
        {
            InitializeComponent();
            centreP = new float[] { Width / 2, Height / 2 };
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
                            map.SetPixel((int)position[0],(int) position[1], Color.Red);
                        }
                        
                    }
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
        private void segmentNaiveLine(PointF x0,PointF x1,int a,int b,Bitmap m)
        {
            PointF current = x0;
            int r = (int)(a * x0.X - b * x0.Y );
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
                m.SetPixel((int)position[0], (int)position[1], Color.Black);
            }
        }
    }
}
