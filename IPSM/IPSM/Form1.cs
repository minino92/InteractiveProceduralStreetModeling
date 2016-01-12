﻿using System;
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
        public IPSM()
        {
            InitializeComponent();
            listPositionFieldTensor = new List<FieldTensor>();
            noise = new Noise();
            pictureZone.Size = new System.Drawing.Size(Noise.size, Noise.size);
            Size = new System.Drawing.Size(Noise.size + 300, Noise.size+100);
        }

        private void IPSM_Load(object sender, EventArgs e)
        {

        }

        private void IPSM_Paint(object sender, PaintEventArgs e)
        {
            pictureZone.Size = new System.Drawing.Size(1000, 1000);
         
            th = new Thread(() =>
            {
                var g = pictureZone.CreateGraphics();
                using (var bmp = new Bitmap(Noise.size, Noise.size, g)) 
                {
                    TensorField tf = new TensorField(Noise.size);
                    tf.generateGridTensorField(bmp,g,(float)Math.PI*8/6);
                    g.DrawImage(bmp, new PointF(0, 0));
                    //using (var redPen = new Pen(Color.Red)) 
                    //{
                    //    for (int i = 0; i < Noise.size; i++ )
                    //    {
                    //        for (int j = 0; j < Noise.size; j++ )
                    //        {                                
                    //            bmp.SetPixel(i, j, Color.FromArgb(noise.noiseTable[i,j],noise.noiseTable[i,j],noise.noiseTable[i,j]));  
                    //        }
                    //    }
                    //    g.DrawImage(bmp, new PointF(0, 0));

                    //    foreach (FieldTensor fieldTensor in listPositionFieldTensor)
                    //    {
                    //        g.DrawRectangle(new Pen(Color.Red, 5), new Rectangle(fieldTensor.position, new System.Drawing.Size(7, 7)));
                    //        if (fieldTensor.finalPosition != new System.Drawing.Point())
                    //        {
                    //            g.DrawLine(new Pen(Color.LightGreen, 3), fieldTensor.position, fieldTensor.finalPosition);
                    //        }
                    //    }
                        
                    //}
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
    }
}
