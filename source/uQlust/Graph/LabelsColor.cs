using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace ClusterV.graph
{
    public partial class LabelsColor : Form
    {
        public LabelsColor(Dictionary <string,int[]> labelColors)
        {
            SolidBrush drawBrush = new System.Drawing.SolidBrush(System.Drawing.Color.Black);
            InitializeComponent();
            int x = 10, y = 25;
            foreach (var item in labelColors)
            {
                drawBrush = new System.Drawing.SolidBrush(Color.FromArgb(item.Value[0], item.Value[1], item.Value[2]));
               // e.Graphics.FillRectangle(drawBrush, x, y, 15, 10);
                drawBrush = new System.Drawing.SolidBrush(System.Drawing.Color.Black);
                //e.Graphics.DrawString(item.Key, drawFont, drawBrush, x + 20, y);
                y += 25;
                if (y > this.Size.Height)
                {
                    x += 150;
                    y = 25;
                }

            }

        }
    }
}
