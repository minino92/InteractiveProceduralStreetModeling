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
        List<FieldTensor> listPositionFieldTensor;
        private bool dragTensorField = false;
        private bool visualizaChoice = false;
        public IPSM()
        {
            InitializeComponent();
            comboBox1.SelectedIndex = 0;
            listPositionFieldTensor = new List<FieldTensor>();
            noise = new Noise();
            pictureZone.Size = new System.Drawing.Size(Noise.size, Noise.size);
            Size = new System.Drawing.Size(Noise.size + 250, Noise.size+50);
            numberTensorFields.Value = 16;
        }

        private void IPSM_Load(object sender, EventArgs e)
        {

        }

        private void IPSM_Paint(object sender, PaintEventArgs e)
        {
            //pictureZone.Size = new System.Drawing.Size(1000, 1000);
         
            th = new Thread(() =>
            {
                var g = pictureZone.CreateGraphics();
                g.Clear(Color.White);
                using (var bmp = new Bitmap(Noise.size, Noise.size, g)) 
                {
                    TensorField tf = new TensorField(Noise.size);
                    tf.NumberOfTensorsToDisplay = (int)numberTensorFields.Value;
                    tf.generateGridTensorField(bmp, g, (float)Math.PI / 6);
                    g.DrawImage(bmp, new PointF(0, 0));                   
                    if (!visualizaChoice)
                    {
                        tf.exportEigenVectorsImage(bmp, g);
                    }
                    else
                    {
                        StreetGraph sg = new StreetGraph(new PointF(0, Noise.size), new PointF(Noise.size, 0), tf, 30f);
                        sg.computeMajorHyperstreamlines(bmp, g);
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

        private void MouseClick(object sender, MouseEventArgs e)
        {
            if (dragTensorField)
            {
                dragTensorField = false;
                listPositionFieldTensor[listPositionFieldTensor.Count-1].finalPosition = e.Location;
                Invalidate();
            }
            else
            {
                if (e.Location.X < Noise.size && e.Location.Y < Noise.size)
                {
                    FieldTensor field = new FieldTensor();
                    field.position = e.Location;
                    listPositionFieldTensor.Add(field);
                    //log.Text = e.Location.ToString();
                    dragTensorField = true;
                    Invalidate();
                }
            }
        }

        private void ChangeNumberTensorFieldToDisplay(object sender, EventArgs e)
        {
            //Invalidate();
        }

        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            if (comboBox1.SelectedIndex == 0) visualizaChoice = false;
            else visualizaChoice = true;
            Invalidate();
        }        
    }
}
