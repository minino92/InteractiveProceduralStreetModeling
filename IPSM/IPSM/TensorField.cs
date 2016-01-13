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
        public static EigenVector operator *(EigenVector ev, int l)
        {
            EigenVector res = new EigenVector();
            res.X= ev.X * -1;
            res.Y=ev.Y * -1;
            res.Z=ev.Z * -1;
            res.W=ev.W * -1;
            return res;
        }
    }
    public class EigenValue
    {
        //pour une matrice symétrique carré lambda=1ou-1
        public double X { get; set; }
        public double Y { get; set; }       
    }
    class TensorField
    {
        public Tensor[,] matrixTensors;
        public EigenVector[,] matrixEigenVectors;
        private EigenValue[,] matrixEigenValues;
        private bool matrixFieldIsFilled;
        private bool matrixEigenIsComputed;
        private int numberOfTensorsToDisplay;

        public int NumberOfTensorsToDisplay
        {
            get { return numberOfTensorsToDisplay; }
            set { numberOfTensorsToDisplay = value; }
        }

        public TensorField(int size)
        {
            matrixTensors = new Tensor[Noise.size, Noise.size];
            matrixFieldIsFilled = false;
            matrixEigenIsComputed = false;
            matrixEigenVectors = null;
            matrixEigenValues = null;
        }

        //get tensor of the row and column entered
        public Tensor getTensor(int row, int column)
        {
            return matrixTensors[row, column];
        }

        //add tensor in the row and column entered
        public void setTensor(int row, int column,Tensor tensor)
        {
            matrixTensors[row, column] = tensor;
        }

        //basic function to create, compute and get a bitmap
        public void generateGridTensorField(Bitmap bmp,Graphics g,float theta)
        {            
            fillMatrixBasisField(theta, 1);
            computeTensorsEigenDecomposition();
            //exportEigenVectorsImage(bmp,g);//we have isolate this function so that we can have controll about visualizing 2 views
        }

        //fill the matrix with the theta and value L entered
        public void fillMatrixBasisField(float theta,int valueL)
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
                    t *= valueL;
                    t.theta = theta;
                    matrixTensors[i, j] = t;
                }
            }
            matrixFieldIsFilled = true;
        }

        //function that calculates the eigenVector and eigenValue of each point in the matrix and returns the number of degenerate points of the matrix
        public int computeTensorsEigenDecomposition()
        {
            if (!matrixFieldIsFilled)
            {
                throw new NotImplementedException("Erreur dans le calcul des points degénérés");
            }
            if (matrixEigenVectors == null || matrixEigenValues == null)
            {
                matrixEigenVectors = new EigenVector[Noise.size, Noise.size];
                matrixEigenValues = new EigenValue[Noise.size, Noise.size];
            }
            int numberOfDegeneratePoints = 0;
            for (int i = 0; i < Noise.size; i++)
            {
                for (int j = 0; j < Noise.size; j++)
                {
                    EigenVector v = new EigenVector { X = Math.Cos(matrixTensors[i, j].theta), Y = Math.Sin(matrixTensors[i, j].theta), Z = Math.Cos(matrixTensors[i, j].theta + Math.PI / 2), W = Math.Sin(matrixTensors[i, j].theta + Math.PI / 2) };
                    matrixEigenVectors[i, j] = v;
                    matrixEigenValues[i, j] = new EigenValue { X = 1, Y = -1 };
                    if (isDegenerate(matrixEigenVectors[i, j]))
                    {
                        numberOfDegeneratePoints++;
                    }
                }
            }
            matrixEigenIsComputed = true;
            return numberOfDegeneratePoints;
        }   

        //return back a bitmap with the final tensor field
        public void exportEigenVectorsImage(Bitmap bmp,Graphics g)
        {
            if (g == null)
            {
                throw new NotImplementedException("graphics are not initialed");
            }
            if (bmp == null)
            {
                bmp = new Bitmap(Noise.size, Noise.size, g);
            }
            if (!matrixFieldIsFilled)
            {
                throw new NotImplementedException("Matrix is not filled");
            }

            PointF origin = new PointF(0.5f, 0.5f);
            int scaleI = Convert.ToInt32(Noise.size / numberOfTensorsToDisplay);//size between tensors to display
            int scaleJ = scaleI;
            for(int i=5; i<Noise.size ; i=i+scaleI)
            {                
                for(int j=5; j<Noise.size; j=j+scaleJ)
                {
                    //Eigen Major Vector
                    PointF temp=new PointF(j, i);
                    PointF initPosLineMajor =new PointF(origin.X+temp.X,origin.Y+temp.Y);                        
                    EigenVector eigenMajorVector = matrixEigenVectors[i,j];
                    eigenMajorVector.X = eigenMajorVector.X * 0.5 / 2.0f * scaleI * 0.8;
                    eigenMajorVector.Y = eigenMajorVector.Y * 0.5 / 2.0f * scaleJ * 0.8;                    
                    PointF finalPosLineMajor =new PointF(initPosLineMajor.X+(float)eigenMajorVector.X,initPosLineMajor.Y+(float)eigenMajorVector.Y);
                    temp.X=initPosLineMajor.X;
                    temp.Y=initPosLineMajor.Y;
                    initPosLineMajor =new PointF(temp.X-(float) eigenMajorVector.X,temp.Y-(float)eigenMajorVector.Y);
                    g.DrawLine(new Pen(Color.Red),initPosLineMajor,finalPosLineMajor);

                    temp = new PointF(j , i);
                    PointF initPosLineMinor = new PointF(origin.X + temp.X, origin.Y + temp.Y);
                    //Eigen Minor Vector
                    EigenVector eigenMinorVector = matrixEigenVectors[i, j];
                    eigenMajorVector.X = eigenMinorVector.Z * 0.5 / 2.0f * scaleI * 0.8;
                    eigenMajorVector.Y = eigenMinorVector.W * 0.5 / 2.0f * scaleJ * 0.8;
                    PointF finalPosLineMinor = new PointF(initPosLineMinor.X + (float)eigenMinorVector.Z, initPosLineMinor.Y + (float)eigenMinorVector.W);
                    temp.X = initPosLineMinor.X;
                    temp.Y = initPosLineMinor.Y;
                    initPosLineMinor = new PointF(temp.X - (float)eigenMajorVector.X, temp.Y - (float)eigenMajorVector.Y);
                    g.DrawLine(new Pen(Color.Black), initPosLineMinor, finalPosLineMinor);
                }
            }         
        }    

        //function that says if a eigenVector is defenerate or not
        public static bool isDegenerate(EigenVector ev)
        {
            const double epsilon=1e-5;
            return ev.X < epsilon && ev.Y < epsilon && ev.Z < epsilon && ev.W < epsilon;
        }
    }
}
