using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Graph
{
    public partial class DrawPanel : Form
    {
        public delegate void DrawOnPictureBox();
        List<Tuple<string, string>> labels = new List<Tuple<string, string>>();
        public DrawOnPictureBox drawPic;
        public Bitmap bmp;
        public DrawPanel(string name,List<Tuple<string, string>> labelsR)
        {            
            InitializeComponent();
            labels = labelsR;
            this.Text = name;
        }

        public void CreateBMP()
        {
            if(pictureBox1.Width>0 && pictureBox1.Height>0)
                bmp = new Bitmap(pictureBox1.Width, pictureBox1.Height);
        }
        private void pictureBox1_Paint(object sender, PaintEventArgs e)
        {
            drawPic();           
            e.Graphics.DrawImage(bmp,0,0);
        }

        private void pictureBox1_Resize(object sender, EventArgs e)
        {
            CreateBMP();
            pictureBox1.Refresh();
            panel1.Refresh();
        }

        private void panel1_Paint(object sender, PaintEventArgs e)
        {

            Brush b = new SolidBrush(Color.Black);
            if (labels != null)
            {
                float stepY = bmp.Height / (labels.Count);
                float xx = stepY / 2;
                if (xx < 1)
                    xx = 1;
                Font f = new Font("arial", xx);
                float startPos = stepY / 4;
                if (labels != null)
                    for (int i = 0; i < labels.Count; i++)
                    {
                        e.Graphics.DrawString(labels[i].Item1, f, b, 10, startPos + i * stepY);
                        if (labels[i].Item2.Length > 0)
                            e.Graphics.DrawString(labels[i].Item2, f, b, bmp.Width - 20, startPos + i * stepY);
                    }
            }
        }
    }
}
