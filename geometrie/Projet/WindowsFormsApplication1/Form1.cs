using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Printing;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows;
using System.IO;

namespace WindowsFormsApplication1
{
    public partial class Form1 : Form
    {
        List<PointF> points;
        private Graphics g;
        private float[] centerScreen;
        private int showObjetNum = 0;
        private int N;
        private Vector dM;
        private string path = @"Datas\coeur_30.txt";
        private List<PointF> continu;

        public Form1()
        {
            InitializeComponent();
            centerScreen = new float[] { Width / 2, Height / 2 };
            points = new List<PointF>();
            continu = new List<PointF>();
        }

        public void setPoints(List<PointF> points)
        {
            this.points = points;
        }


        private void Form1_Load(object sender, EventArgs e)
        {

        }

        private void Form1_Paint(object sender, PaintEventArgs e)
        {
            try
            {
                using (g = CreateGraphics())
                {
                    Bitmap map = new Bitmap(Width, Height, g);
                    g.DrawPolygon(new Pen(new SolidBrush(Color.Red), 10), new PointF[] { new PointF(0, 0), new PointF(0, 0) });
                    StreamReader reader = new StreamReader(path);
                    string pattern = @"((?<x>-?\d+),(?<y>-?\d+))";
                    while (!reader.EndOfStream)
                    {
                        string line = reader.ReadLine();
                        Regex rg = new Regex(pattern);
                        MatchCollection m = rg.Matches(line);
                        foreach (Match ms in m)
                        {
                            int[] position = ChangeToCenterScreen(float.Parse(ms.Groups["x"].Value), float.Parse(ms.Groups["y"].Value));
                            float x = float.Parse(ms.Groups["x"].Value);
                            float y = float.Parse(ms.Groups["y"].Value);
                            if (x >= 0 && y >= 0 && y < x)
                            {
                                //map.SetPixel((int)position[0], (int)position[1], Color.Red);
                                g.DrawRectangle(new Pen(Color.Black), (int)position[0], (int)position[1], 2, 2);
                            }
                            else
                                g.DrawRectangle(new Pen(Color.Black), (int)position[0], (int)position[1], 2, 2);
                            //map.SetPixel((int)position[0], (int)position[1], Color.Black);

                            points.Add(new PointF(float.Parse(ms.Groups["x"].Value), float.Parse(ms.Groups["y"].Value)));
                        }

                    }
                    g.DrawPolygon(new Pen(new SolidBrush(Color.Red), 10), new PointF[] { new PointF(0, 0), new PointF(0, 0) });

                    N = points.Count;
                    CalculateDM(0);
                    Segmentation((int)points[0].X, (int)points[0].Y, map,g);
                    //for (int i = 0; i < points.Count; i++)
                    //{
                    //    int[] positionedPixel = ChangeToCenterScreen(points[i].X, points[i].Y);
                    //    map.SetPixel(positionedPixel[0], positionedPixel[1], Color.Red);
                    //}

                    //NaiveLine(new PointF(0, 0), new PointF(25, 25), 3, 8, 0, map, Color.Red);
                    //g.DrawPolygon(new Pen(Color.Red), continu.ToArray());
                    g.DrawImage(map, new PointF(0, 0));
                }
            }
            catch (Exception)
            {
            }
        }

        private void CalculateDM(int v)
        {
            if (v >= N - 1)
            {
                dM = new Vector(points[0].X - points[N - 1].X, points[0].Y - points[N - 1].Y);
            }
            else
            {
                dM = new Vector(points[v + 1].X - points[v].X, points[v + 1].Y - points[v].Y);
            }
        }

        private int[] ChangeToCenterScreen(float x, float y)
        {
            int[] res = new int[2];
            float[] u = new float[] { 5, 0, 0 };
            float[] v = new float[] { 0, -5, 0 };
            float[] w = new float[] { 0, 0, 0 };
            res[0] = (int)(u[0] * x + v[0] * y + centerScreen[0]);
            res[1] = (int)(u[1] * x + v[1] * y + centerScreen[1]);
            return res;
        }

        private void NaiveLine(PointF x0, PointF x1, int a, int b, int mu, Bitmap m, Color color,Graphics g,bool type=true)
        {
            PointF currentPoint = x0;
            int r = (int)(a * x0.X - b * x0.Y - mu);
            while (currentPoint.X < x1.X)
            {
                if (r >= (b - a))
                {
                    currentPoint.Y++;
                    r = r - b + a;
                }
                else
                {
                    r += a;
                }
                currentPoint.X++;
                int[] position = ChangeToCenterScreen(currentPoint.X, currentPoint.Y);

                if (type) m.SetPixel(position[0], position[1], color);
                else g.DrawRectangle(new Pen(color), position[0], position[1], 2, 2);
            }
        }

        private void ExchangeAxis(int alfa, int beta, int gamma, int delta, int phi, int psi)
        {
            int temp;
            temp = alfa; alfa = gamma; gamma = temp;
            temp = beta; beta = delta; delta = temp;
            temp = phi; phi = psi; psi = temp;
        }

        private void ChangeVerticalAxis(int alfa, int beta, int gamma, int delta, int phi, int psi)
        {
            if (delta != 0)
            {
                delta = -delta;
                psi = -psi;
            }
            else
            {
                if (beta != 0)
                {
                    beta = -beta;
                    phi = -phi;
                }
            }
        }

        private void ChangeHorizontalAxis(int alfa, int beta, int gamma, int delta, int phi, int psi)
        {
            if (alfa != 0)
            {
                alfa = -alfa;
                phi = -phi;
            }
            else
            {
                if (gamma != 0)
                {
                    gamma = -gamma;
                    psi = -psi;
                }
            }
        }

        private void Rotation(int xh, int yh, int dx, int dy,
            out int alfa, out int beta, out int gamma, out int delta, out int phi, out int psi)
        {

            alfa = 0; beta = 0; phi = 0;
            gamma = 0; delta = 0; psi = 0;
            if (dx == 1 && dy >= 0)
            {
                alfa = 1; beta = 0; phi = -xh;
                gamma = 0; delta = 1; psi = -yh;
            }
            else
            {
                if (dx <= 0 && dy == 1)
                {
                    alfa = 0; beta = 1; phi = -yh;
                    gamma = -1; delta = 0; psi = xh;
                }
                else
                {
                    if (dx == -1 && dy <= 0)
                    {
                        alfa = -1; beta = 0; phi = xh;
                        gamma = 0; delta = -1; psi = yh;
                    }
                    else
                    {
                        if (dx >= 0 && dy == -1)
                        {
                            alfa = 0; beta = -1; phi = yh;
                            gamma = 1; delta = 0; psi = -xh;
                        }
                    }
                }
            }
        }

        private void AdjustAxis(int dx, int dy, int ddx, int ddy,
            int alfa, int beta, int gamma, int delta, int phi, int psi)
        {
            if (dy == 0 && ddx == dx && ddy == -dx)
            {
                ChangeVerticalAxis(alfa, beta, gamma, delta, phi, psi);
            }
            else
            {
                if (dx == 0 && ddx == dy && ddy == dy)
                {
                    ChangeHorizontalAxis(alfa, beta, gamma, delta, phi, psi);
                }
                else
                {
                    if (dx == dy && ddx == 0 && ddy == dy || dx == -dy && ddx == dx && ddy == 0)
                    {
                        ExchangeAxis(alfa, beta, gamma, delta, phi, psi);
                    }
                }
            }
        }

        private void MovingFrame(out int k, out int xk, out int yk,
             out int a, out int b, out int mu, int h, int xh, int yh,
             out int alfa, out  int beta, out int gamma, out int delta, out int phi, out int psi,
             out int Ux, out int Uy, out int Lx, out int Ly)
        {

            int dx, dy;
            Ux = 0;
            Uy = 0;
            Lx = 0;
            Ly = 0;
            a = 0;
            b = 0;
            mu = 0;
            k = h; xk = xh; yk = yh;
            CalculateDM(h + 1);
            dx = (int)dM.X;
            dy = (int)dM.Y;
            while (k < N && (int)dM.Y == dy && (int)dM.X == dx)
            {
                k++;
                CalculateDM(k);
                xk = (int)(xk + dM.X);
                yk = (int)(yk + dM.Y);
            }
            Rotation(xh, yh, dx, dy, out alfa, out beta, out gamma, out delta, out phi, out psi);
            if (k < N)
            {
                CalculateDM(k + 1);
                AdjustAxis(dx, dy, (int)dM.X, (int)dM.Y, alfa, beta, gamma, delta, phi, psi);
                if (dx == 0 || dy == 0)
                {
                    Lx = k - h;
                }
                else
                {
                    Lx = 0;
                }
                a = 1;
                b = Lx + 1;
                mu = 0;
                Ly = 0;
                Ux = 0;
                Uy = 0;
            }
        }

        private void RecognizeSegment(int k, int xk, int yk, int a, int b, int mu,
            int alfa, int beta, int gamma, int delta, int phi, int psi,
            int Ux, int Uy, int Lx, int Ly)
        {
            int r, ddx, ddy, ddx_loc, ddy_loc, xx, yy, x_loc, y_loc;

            while (k < N)
            {
                CalculateDM(k + 1);
                ddx = (int)dM.X;
                ddy = (int)dM.Y;
                ddx_loc = alfa * ddx + beta * ddy;
                ddy_loc = gamma * ddx + delta * ddy;
                if (ddx_loc <= 0 || ddy_loc < 0)
                {
                    return;
                }
                xx = xk + ddx;
                yy = yk + ddy;
                x_loc = alfa * xx + beta * yy + phi;
                y_loc = gamma * xx + delta * yy + psi;
                r = a * x_loc - b * y_loc;
                if (r < mu - 1 || r > mu + b)
                {
                    return;
                }
                k = k + 1;
                xk = xx;
                yk = yy;
                if (r == mu + b)
                {
                    Ux = x_loc - b;
                    Uy = y_loc + 1 - a;
                    a = y_loc - Ly;
                    b = x_loc - Lx;
                    mu = a * Lx - b * Ly - b + 1;
                }
                else
                {
                    if (r == mu - 1)
                    {
                        Ux = x_loc - b;
                        Uy = y_loc - 1 - a;
                        a = y_loc - Uy;
                        b = x_loc - Ux;
                        mu = a * Ux - b * Uy;
                    }
                }
            }
        }

        private void Segmentation(int x0, int y0, Bitmap map,Graphics g)
        {
            int h, k, xk, yk, xh, yh, a, b, mu, Ux, Uy, Lx, Ly;
            int alfa, beta, gamma, delta, phi, psi;
            h = 0;
            xh = x0;
            yh = y0;
            //used for drawing continuous line
            PointF temp = new PointF();
            temp.X = xh;
            temp.Y = yh;
            while (h < N)
            {
                MovingFrame(out k, out xk, out yk, out a, out b, out mu, h, xh, yh, out alfa, out beta, out gamma, out delta, out phi, out psi, out Ux, out Uy, out Lx, out Ly);
                RecognizeSegment(k, xk, yk, a, b, mu, alfa, beta, gamma, delta, phi, psi, Ux, Uy, Lx, Ly);
                NaiveLine(new PointF(xh, yh), new PointF(xk, yk), a, b, mu, map, Color.Blue,g,false);
                int[] pp = ChangeToCenterScreen(temp.X, temp.Y);
                int[] pppp = ChangeToCenterScreen(xk, yk);
                g.DrawLine(new Pen(Color.Red), new PointF(pp[0], pp[1]), new PointF(pppp[0], pppp[1]));
                temp.X = xk;
                temp.Y = yk;
                //continu.Add(new PointF(pp[0], pp[1]));
                //continu.Add(new PointF(pppp[0], pppp[1]));
                if (k == N)
                {
                    h = N;
                }
                else
                {
                    h = k + 1;
                    CalculateDM(h);
                    xh = (int)(xk + dM.X);
                    yh = (int)(yk + dM.Y);
                }
            }

        }       
    }
}

//using System;
//using System.Collections.Generic;
//using System.ComponentModel;
//using System.Data;
//using System.Drawing;
//using System.Drawing.Drawing2D;
//using System.Linq;
//using System.Text;
//using System.Text.RegularExpressions;
//using System.Threading.Tasks;
//using System.Windows.Forms;
//using System.IO;

//namespace WindowsFormsApplication1
//{
//    public partial class Form1 : Form
//    {
//        private Graphics g;
//        private float[] centreP;
//        private string path = @"Datas\cercle_30.txt";
//        private List<PointF> points;
//        private List<PointF> points_octants;
//        //variable globales
//        private PointF X0;
//        private PointF dM;
//        private int N;//nombre de pixel dans E
//        public Form1()
//        {
//            InitializeComponent();
//            centreP = new float[] { Width / 2, Height / 2 };
//            points = new List<PointF>();
//            points_octants = new List<PointF>();
//        }

//        private void Form1_Paint(object sender, PaintEventArgs e)
//        {
//            try
//            {
//                using (g = CreateGraphics())
//                {
//                    Bitmap map = new Bitmap(Width,Height,g);
//                    g.DrawPolygon(new Pen(new SolidBrush(Color.Red), 10), new PointF[] { new PointF(0, 0), new PointF(0, 0) });
//                    StreamReader reader = new StreamReader(path);
//                    string pattern = @"((?<x>-?\d+),(?<y>-?\d+))";
//                    while (!reader.EndOfStream)
//                    {
//                        string line = reader.ReadLine();
//                        Regex rg = new Regex(pattern);
//                        MatchCollection m=rg.Matches(line);
//                        foreach (Match ms in m)
//                        {
//                            float[] position = ChangeBase(float.Parse(ms.Groups["x"].Value), float.Parse(ms.Groups["y"].Value));
//                            float x=float.Parse(ms.Groups["x"].Value);
//                            float y=float.Parse(ms.Groups["y"].Value);
//                            if ( x>= 0 &&  y>= 0 && y<x)
//                            {
//                                //map.SetPixel((int)position[0], (int)position[1], Color.Red);
//                                g.DrawRectangle(new Pen(Color.Black), (int)position[0], (int)position[1], 2, 2);
//                                points_octants.Add(new PointF(float.Parse(ms.Groups["x"].Value), float.Parse(ms.Groups["y"].Value)));
//                            }
//                            else
//                                g.DrawRectangle(new Pen(Color.Black), (int)position[0], (int)position[1], 2, 2); 
//                                //map.SetPixel((int)position[0], (int)position[1], Color.Black);
                                
//                            points.Add(new PointF(float.Parse(ms.Groups["x"].Value), float.Parse(ms.Groups["y"].Value)));
//                        }
                        
//                    }
//                    //initialisation du variable globale
//                    N = points.Count;
//                    X0 = points[0];
//                    //dM = new PointF(points[1].X - X0.X, points[1].Y - X0.Y);
//                    CalucateDM(0);
//                    Console.WriteLine("nombre de points:" + N);
//                    float[] position_ = ChangeBase(points[0].X, points[0].Y);
//                    g.DrawRectangle(new Pen(Color.DodgerBlue), (int)position_[0], (int)position_[1], 2, 2);
//                    Segmentation(X0, dM,map,g);
//                    //segmentNaiveLine(new PointF(0, 0), new PointF(10, 10), 3, 5,0, map,g);
//                    g.DrawImage(map, new Point(0, 0));
//                }
//            }
//            catch (Exception ex)
//            {

//            }
//        }
//        /// <summary>
//        /// Un répère de l'écran (0,0) top left
//        /// Un répère de l'image (width/2,height/2)
//        /// </summary>
//        /// <param name="x">Position x dans le répère P</param>
//        /// <param name="y">Position y dans le répère P</param>
//        /// <returns></returns>
//        private float[] ChangeBase(float x,float y)
//        {            
//            float[] res = new float[2];
//            //u=(5,0,0) et v=(0,-5,0) => 1px -> x unité
//            float[] u =new float[]{5,0,0};
//            float[] v=new float[]{0,-5,0};
//            float[] w=new float[]{0,0,0};//car on est en 2D
//            res[0] = u[0] * x + v[0] * y + centreP[0];
//            res[1] = u[1] * x + v[1] * y + centreP[1];
//            return res;
//        }
//        /// <summary>
//        /// A segment of a naive line
//        /// </summary>
//        private void segmentNaiveLine(PointF x0,PointF x1,int a,int b,int mu,Bitmap m,Graphics g,bool type=true)
//        {
//            PointF current = x0;
//            int r = (int)(a * x0.X - b * x0.Y -mu);
//            while (current.X < x1.X)
//            {
//                if (r >= (b - a))
//                {
//                    current.Y ++;
//                    r = r - b + a;
//                }
//                else
//                {
//                    r += a;
//                }
//                current.X ++;
//                float[] position = ChangeBase(current.X, (int)current.Y);
//                //difference leaning points
//                int xx = (int)position[0];
//                int yy = (int)position[1];
//                if (type)
//                {
//                    //points d'appui cf thèse partie 3.2
//                    if (a * current.X - b * current.Y == mu)
//                    {
//                        //m.SetPixel(xx, yy, Color.Black);
//                        g.DrawRectangle(new Pen(Color.Black), xx, yy, 2, 2);
//                    }
//                    else if (a * current.X - b * current.Y == (mu + b - 1))
//                    {
//                        //m.SetPixel(xx, yy, Color.Blue);
//                        g.DrawRectangle(new Pen(Color.Blue), xx, yy, 2, 2);
//                    }
//                    else g.DrawRectangle(new Pen(Color.Red), xx, yy, 2, 2);//m.SetPixel(xx, yy, Color.Red);
//                }
//                else g.DrawRectangle(new Pen(Color.Green), xx, yy, 2, 2);
                
//            }
//        }        
//        void ExchangeAxis(int alpha, int beta, int gama, int delta, int phi, int psi)
//        {
//            int temp;
//            temp = alpha; alpha = gama; gama = temp;
//            temp = beta; beta = delta; delta = temp;
//            temp = phi; phi = psi; psi = temp;
//        }
//        void ChangeVerticalAxis(int alpha, int beta, int gama, int delta, int phi, int psi)
//        {
//            if (delta != 0)
//            {
//                delta = -delta;
//                psi = -psi;
//            }
//            else
//            {
//                if (beta != 0)
//                {
//                    beta = -beta;
//                    phi = -phi;
//                }
//            }
//        }
//        void ChangeHorizontalAxis(int alpha, int beta, int gama, int delta, int phi, int psi)
//        {
//            if (alpha != 0)
//            {
//                alpha = -alpha;
//                phi = -phi;
//            }
//            else
//            {
//                if (gama != 0)
//                {
//                    gama = -gama;
//                    psi = -psi;
//                }
//            }                
//        }
//        void AdjustAxis(int dx, int dy, int ddx, int ddy,int alpha,int beta,int gama,int delta,int phi,int psi)
//        {
//            if (dy == 0 && ddx == dx && ddy == -dx)
//            {
//                ChangeVerticalAxis(alpha, beta, gama, delta, phi, psi);
//            }
//            else
//            {
//                if (dx == 0 && ddx == dy && ddy == dy)
//                {
//                    ChangeHorizontalAxis(alpha, beta, gama, delta, phi, psi);
//                }
//                else
//                {
//                    if ((dx == dy && ddx == 0 && ddy == dy) || (dx == -dy && ddy == 0 && ddx == dy))
//                    {
//                        ExchangeAxis(alpha, beta, gama, delta, phi, psi);
//                    }
//                }
//            }
//        }
//        void Rotation(PointF xh, int dx, int dy, out int alpha, out int beta, out int gama, out int delta, out int phi, out int psi)
//        {
//            alpha = beta = gama = delta = phi = psi = 0;
//            if (dx == 1 && dy >= 0)
//            {
//                alpha = 1; beta = 0; phi = -(int)xh.X;
//                gama = 0; delta = 1; psi = -(int)xh.Y;
//            }
//            else
//            {
//                if (dx <= 0 && dy == 1)
//                {
//                    alpha = 0; beta = 1; phi = -(int)xh.Y;
//                    gama = -1; delta = 0; psi = (int)xh.X;
//                }
//                else 
//                {
//                    if (dx == -1 && dy <= 0)
//                    {
//                        alpha = -1; beta = 0; phi = (int)xh.X;
//                        gama = 0; delta = -1; psi = (int)xh.Y;
//                    }
//                    else
//                    {
//                        if (dx >= 0 && dy == -1)
//                        {
//                            alpha = 0; beta = -1; phi = (int)xh.Y;
//                            gama = 1; delta = 0; psi = -(int)xh.X;
//                        }
//                    }
//                }                
//            } 
//        }
//        /// <summary>
//        /// This procedure returns index k et Mk coordinates
//        /// </summary>
//        /// <param name="k"></param>
//        /// <param name="xk"></param>
//        /// <param name="alpha"></param>
//        /// <param name="beta"></param>
//        /// <param name="gama"></param>
//        /// <param name="delta"></param>
//        /// <param name="phi"></param>
//        /// <param name="psi"></param>
//        /// <param name="a"></param>
//        /// <param name="b"></param>
//        /// <param name="mu"></param>
//        /// <param name="U"></param>
//        /// <param name="L"></param>
//        void MovingFrame(out int k,out PointF xk, out int alpha, out int beta, out int gama, out int delta, out int phi, out int psi, out int a, out int b, out int mu, out PointF U, out PointF L,int h,PointF xh)
//        {
//            int dx, dy, ddx, ddy;
//            L = new PointF(0, 0);
//            U = new PointF(0, 0);
//            a = 0;
//            b = 0;
//            mu = 0;
//            k=h;
//            xk=xh;            
//            CalucateDM(h + 1);
//            dx = (int)dM.X;
//            dy = (int)dM.Y;            
//            while (k < N && (int)dM.X == dx && (int)dM.Y == dy)
//            {               
//                k++;
//                CalucateDM(k);
//                xk.X += dM.X;
//                xk.Y += dM.Y;
//            }
//            Rotation(xh, dx, dy, out alpha, out beta, out gama, out delta, out phi, out psi);
//            if (k < N)
//            {
//                //here is the ddx and ddy
//                CalucateDM(k + 1);
//                AdjustAxis(dx, dy, (int)dM.X, (int)dM.Y, alpha, beta, gama, delta, phi, psi);                
//                if (dx == 0 || dy == 0)
//                {
//                    L.X = k - h;
//                }
//                else
//                {
//                    L.X = 0;
//                }
//                L.Y = 0;
//                a = 1;
//                b = (int)L.X + 1;
//                mu = 0;
//            }            
//        }
//        void CalucateDM(int index)
//        {            
//            if (index >= N - 1)
//            {
//                //what happen here
//                dM = new PointF(points[0].X - points[N - 1].X, points[0].Y - points[N - 1].Y);
//            }
//            else
//            {
//                dM = new PointF(points[index + 1].X - points[index].X, points[index + 1].Y - points[index].Y);
//            }
//        }
//        /// <summary>
//        /// 
//        /// </summary>
//        /// <param name="x0">provient du variable globale X0</param>
//        /// <param name="dm">la direction du courbe</param>
//        void Segmentation(PointF x0, PointF dm,Bitmap map,Graphics g)
//        {
//            int h = 0;
//            int k, a, b, mu;
//            int alpha, beta, gama, delta, phi, psi;
//            PointF xh=x0;//current point
//            PointF xk;
//            PointF U, L;
//            while (h < N)
//            {               
//                MovingFrame(out k,out xk, out alpha, out beta, out gama, out delta, out phi, out psi, out a, out b, out mu, out U, out L,h,xh);
//                RecognizeSegment(k,ref xk,ref a,ref b,ref mu, alpha, beta, gama, delta, phi, psi,ref U,ref L);
//                segmentNaiveLine(xh, xk, a, b, mu, map, g,false);
//                Console.WriteLine("h={0}/k={1}/mu={2}/a={3}/b={4}",h, k, mu, a, b);
//                h = k;
//                xh = xk;
//            }
//        }
//        /// <summary>
//        /// Give back vertex Mk
//        /// </summary>
//        /// <param name="k"></param>
//        /// <param name="xk"></param>
//        /// <param name="a"></param>
//        /// <param name="b"></param>
//        /// <param name="mu"></param>
//        private void RecognizeSegment(int k,ref PointF xk,ref int a,ref int b,ref int mu,int alpha,int beta,int gama,int delta,int phi,int psi,ref PointF U,ref PointF L)
//        {
//            int r, ddx, ddy, ddx_loc, ddy_loc, xx, yy, x_loc, y_loc;
//            while (k < N)
//            {
//                CalucateDM(k + 1);
//                ddx = (int)dM.X;
//                ddy = (int)dM.Y;
//                ddx_loc = alpha * ddx + beta * ddy;
//                ddy_loc = gama * ddx + delta * ddy;
//                if (ddx_loc <= 0 || ddy_loc < 0) return;
//                xx = (int)xk.X + ddx;
//                yy = (int)xk.Y + ddy;
//                x_loc = alpha * xx+beta*yy+phi;
//                y_loc = gama * xx + delta * yy + psi;
//                r = a * x_loc - b * y_loc;
//                if (r < mu - 1 || r > mu + b) return;
//                k++;
//                xk.X = xx;
//                xk.Y = yy;
//                if (r == mu + b)
//                {
//                    U.X = x_loc - b;
//                    U.Y = y_loc + 1 - a;
//                    a = y_loc - (int)L.Y;
//                    b = x_loc - (int)L.X;
//                    mu = a * (int)L.X - b * (int)L.Y + 1;
//                }
//                else
//                {
//                    if (r == mu - 1)
//                    {
//                        L.X = x_loc - b;
//                        L.Y = y_loc - 1 - a;
//                        a = y_loc - (int)U.Y;
//                        b = x_loc - (int)U.X;
//                        mu = a * (int)U.X - b * (int)U.Y;
//                    }
//                }
//            }
//        }
//    }
//}
