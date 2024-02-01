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
    public partial class GeneSelection : Form
    {
        public HashSet<string> selection = new HashSet<string>();
        public GeneSelection(List<string> genes)
        {
            InitializeComponent();
            int posX = 10, posY = 10;
            for(int i=0;i<genes.Count;i++)
            {
                CheckBox b = new CheckBox();
                b.Location = new Point(posX, posY);
                b.Text = genes[i];
                b.Width = 200;
                this.Controls.Add(b);

                posX += 200;
                if(posX+200>this.Width)
                {
                    posX = 10;
                    posY += 40;
                }
            }

        }

        private void button1_Click(object sender, EventArgs e)
        {
            for(int i=0;i<this.Controls.Count;i++)
            {
                if(this.Controls[i] is CheckBox)
                {
                    CheckBox b = (CheckBox)Controls[i];
                    if(b.Checked)
                        selection.Add(b.Text);
                }
            }
            DialogResult = DialogResult.OK;
            this.Close();
        }
    }
}
