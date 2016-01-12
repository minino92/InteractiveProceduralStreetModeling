using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Drawing.Drawing2D;
using System.Drawing;

namespace IPSM
{
    public class Tensor
    {
        public double X{get;set;}
        public double Y{get;set;}
        public double Z{get;set;}
        public double W{get;set;}
        public float theta
        {
            get;
            set;
        }
        public static Tensor operator*(Tensor t,int l)
        {
            t.X *= l;
            t.Y *= l;
            t.Z *= l;
            t.W *= l;
            return t;
        }
    }
    public class EigenVector
    {
        public double X { get; set; }
        public double Y { get; set; }
        public double Z { get; set; }
        public double W { get; set; }
    }
    public class EigenValue
    {
        //pour une matrice symétrique carré lambda=1ou-1
        public double X { get; set; }
        public double Y { get; set; }       
    }
    class TensorField
    {
        private Tensor[,] mData;
        private EigenVector[,] mEigenVectors;
        private EigenValue[,] mEigenValues;
        private bool mFieldIsFilled;
        private bool mEigenIsComputed;

        public TensorField(int size)
        {
            mData = new Tensor[Noise.size, Noise.size];
            mFieldIsFilled = false;
            mEigenIsComputed = false;
            mEigenVectors = null;
            mEigenValues = null;
        }
        /// <summary>
        /// get tensor at position (i,j)
        /// </summary>
        /// <param name="i">ligne</param>
        /// <param name="j">colonne</param>
        /// <returns></returns>
        public Tensor getTensor(int i, int j)
        {
            return mData[i, j];
        }
        /// <summary>
        /// set tensor at position (i,j) to tensor
        /// </summary>
        /// <param name="i"></param>
        /// <param name="j"></param>
        /// <param name="tensor"></param>
        public void setTensor(int i, int j,Tensor tensor)
        {
            mData[i, j] = tensor;
        }
        /// <summary>
        /// Génération du tensorfield
        /// </summary>
        /// <param name="theta"></param>
        /// <param name="l"></param>
        public void generateGridTensorField(Bitmap bmp,Graphics g,float theta)
        {            
            fillGridBasisField(theta, 1);
            computeTensorsEigenDecomposition();
            exportEigenVectorsImage(bmp,g);
        }
        /// <summary>
        /// La formule est R*(X Y Z W)
        /// on ajoute un tenseur à chaque pixel
        /// </summary>
        /// <param name="theta"></param>
        /// <param name="l"></param>
        public void fillGridBasisField(float theta,int l)
        {
            for (int i = 0; i < Noise.size; i++)
            {
                for (int j = 0; j < Noise.size; j++)
                {
                    Tensor t = new Tensor();
                    t.X = Math.Cos(2 * theta);
                    t.Y = Math.Sin(2 * theta);
                    t.Z = Math.Sin(2 * theta);
                    t.W = -Math.Cos(2 * theta);
                    t *= l;
                    t.theta = theta;
                    mData[i, j] = t;
                }
            }
            mFieldIsFilled = true;
        }
        public void exportEigenVectorsImage(Bitmap bmp,Graphics g)
        {
            if (bmp == null) throw new NullReferenceException("Bitmap null dans exportEigenVectorsImage");
            if (!mFieldIsFilled)
            {
                throw new NotImplementedException("Tenseur field n'est pas encore prêt");
            }
            PointF origine = new PointF(0.5f, 0.5f);
            int numberOfTensorsToDisplay = 16;
            //longueur des vecteurs
            int scaleI =(Noise.size/numberOfTensorsToDisplay)*3;
            int scaleJ = (Noise.size/numberOfTensorsToDisplay)*3;
            //we increment i and j according to the scaleI and scaleJ so that we can visualize 
            //one tensor each (i+scaleI,j+scaleJ) position
            for(int i=0; i<Noise.size ; i=i+scaleI)
            {
                for(int j=0; j<Noise.size; j=j+scaleJ)
                {
                    PointF temp=new PointF(j*0.5f, i*0.5f);
                        PointF b =new PointF(origine.X+temp.X,origine.Y+temp.Y);
                        
                        EigenVector eigenMajorVector = mEigenVectors[i,j];
                        eigenMajorVector.X = eigenMajorVector.X * 0.5 / 2.0f * scaleI * 0.8;
                        eigenMajorVector.Y = eigenMajorVector.Y * 0.5 / 2.0f * scaleJ * 0.8;                    
                        PointF tip =new PointF(b.X+(float)eigenMajorVector.X,b.Y+(float)eigenMajorVector.Y);
                        temp.X=b.X;
                        temp.Y=b.Y;
                        b =new PointF(temp.X-(float) eigenMajorVector.X,temp.Y-(float)eigenMajorVector.Y);
                        g.DrawLine(new Pen(Color.Red),b,tip);
                        //roundVector2D(base);
                        //roundVector2D(tip);
                        //// Flip the y axis because the painter system has its origin
                        //// on the top left corner and the y axis points down
                        //painter.drawLine(base.x(),imageSize - base.y(),tip.x(),imageSize - tip.y());
                        PointF temp_ = new PointF(j * 0.5f, i * 0.5f);
                        PointF b_ = new PointF(origine.X + temp_.X, origine.Y + temp_.Y);

                        EigenVector eigenMinorVector = mEigenVectors[i, j];
                        eigenMajorVector.X = eigenMinorVector.Z * 0.5 / 2.0f * scaleI * 0.8;
                        eigenMajorVector.Y = eigenMinorVector.W * 0.5 / 2.0f * scaleJ * 0.8;
                        PointF tip_ = new PointF(b_.X + (float)eigenMinorVector.Z, b_.Y + (float)eigenMinorVector.W);
                        temp_.X = b_.X;
                        temp_.Y = b_.Y;
                        b_ = new PointF(temp.X - (float)eigenMajorVector.X, temp.Y - (float)eigenMajorVector.Y);
                        g.DrawLine(new Pen(Color.Black), b_, tip_);
                    //roundVector2D(base);
                    //roundVector2D(tip);
                    //// Flip the y axis because the painter system has its origin
                    //// on the top left corner and the y axis points down
                    //painter.drawLine(base.x(),imageSize - base.y(),tip.x(),imageSize - tip.y());
                    //if(drawVector2)
                    //{
                    //    painter.setPen(pen2);
                    //    QVector2D base = origin + QVector2D(j*dv, i*du);
                    //    QVector2D eigenVector = getTensorMinorEigenVector(mData[i][j]);
                    //    eigenVector.setX(eigenVector.x()*du/2.0f*scaleI*0.8);
                    //    eigenVector.setY(eigenVector.y()*dv/2.0f*scaleJ*0.8);
                    //    QVector2D tip = base + eigenVector;
                    //    base -= eigenVector;
                    //    roundVector2D(base);
                    //    roundVector2D(tip);
                    //    // Flip the y axis because the painter system has its origin
                    //    // on the top left corner and the y axis points down
                    //    painter.drawLine(base.x(),imageSize - base.y(),tip.x(),imageSize - tip.y());
                    //}
                }
            }

            //à traiter            
        }
        public int computeTensorsEigenDecomposition()
        {
            if (!mFieldIsFilled)
            {
                throw new NotImplementedException("Erreur dans le calcul des points degénérés");
            }
            //initialisation
            if (mEigenVectors == null || mEigenValues == null)
            {
                mEigenVectors = new EigenVector[Noise.size, Noise.size];
                mEigenValues = new EigenValue[Noise.size, Noise.size];
            }
            int numberOfDegeneratePoints = 0;
            for (int i = 0; i < Noise.size; i++)
            {
                for (int j = 0; j < Noise.size; j++)
                {
                    EigenVector v = new EigenVector { X = Math.Cos(mData[i, j].theta), Y = Math.Sin(mData[i, j].theta), Z = Math.Cos(mData[i, j].theta + Math.PI / 2), W = Math.Sin(mData[i, j].theta + Math.PI / 2) };
                    mEigenVectors[i, j] = v;
                    mEigenValues[i, j] = new EigenValue { X = 1, Y = -1 };
                    if (isDegenerate(mEigenVectors[i,j]))
                    {
                        numberOfDegeneratePoints++;
                    }
                }
            }
            mEigenIsComputed = true;
            return numberOfDegeneratePoints;
        }       
        public static bool isDegenerate(EigenVector ev)
        {
            double epsilon=1e-5;
            return ev.X < epsilon && ev.Y < epsilon && ev.Z < epsilon && ev.W < epsilon;
        }
    }
}
