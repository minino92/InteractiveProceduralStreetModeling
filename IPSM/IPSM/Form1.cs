using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace IPSM
{
    public partial class IPSM : Form
    {
        /// <summary>
        /// Bitmap is used to store the street graph when it is finished.
        /// </summary>
        Bitmap map;        
        public IPSM()
        {
            InitializeComponent();            
        }

        private void IPSM_Load(object sender, EventArgs e)
        {

        }

        private void IPSM_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = pictureZone.CreateGraphics();
            g.PageUnit = GraphicsUnit.Pixel;
            map = new Bitmap(pictureZone.Width, pictureZone.Height, g);
            g.DrawString("Draw something here", new Font("Arial Black", 40), new SolidBrush(Color.Black), new PointF(5, 10));
            g.DrawLine(new Pen(Color.Red), new PointF(0,250), new PointF(pictureZone.Width,250));
            g.Dispose();
        }        
    }
}
