using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Graph
{
    public partial class Resolution : Form
    {
        int width, height,lineThickness;
        bool legend;
        Color linesColor;
        public int WidthR
        {
            get { return width;}
        }
        public int HeightR
        {
            get { return height; }
        }
        public bool ShowLegend
        {
            get { return legend;}
        }
        public int LineThickness
        {
            get { return lineThickness;}
        }
        public Color LinesColor
        {
            get { return linesColor; }
        }
        public Resolution(int Width,int Height,Color back)
        {
            InitializeComponent();
            width = Width;
            height = Height;
            pictureBox1.BackColor = back;
            linesColor = back;
            numericUpDown1.Value = width;
            numericUpDown2.Value = height;

        }

        private void button1_Click(object sender, EventArgs e)
        {
            width = (int)numericUpDown1.Value;
            height = (int)numericUpDown2.Value;
            legend = checkBox1.Checked;
            lineThickness = (int)numericUpDown3.Value;            
            this.Close();
            
        }

        private void button2_Click(object sender, EventArgs e)
        {
            this.Close();
        }

        private void button3_Click(object sender, EventArgs e)
        {
            DialogResult res=colorDialog1.ShowDialog();
            if(res==DialogResult.OK)
            {
                linesColor = colorDialog1.Color;
                pictureBox1.BackColor = linesColor;
            }
        }
    }
}
