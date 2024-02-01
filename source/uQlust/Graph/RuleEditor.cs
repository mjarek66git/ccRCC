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
    public partial class RuleEditor : Form
    {
        public List<int> upRegulated = new List<int>();
        public List<int> downRegulated = new List<int>();
        public double threshold;
        public RuleEditor(List<string> labels,List<int> upReg,List<int> downReg,double thresh)
        {
            InitializeComponent();
            foreach (var item in labels)
            {
                listUp.Items.Add(item);
                listDown.Items.Add(item);
            }
            if (upReg != null)
                foreach (var item in upReg)
                    listUp.SelectedIndices.Add(item);

            if (downReg != null)
                foreach (var item in downReg)
                    listDown.SelectedIndices.Add(item);

            textBox1.Text = thresh.ToString();

        }

        private void button1_Click(object sender, EventArgs e)
        {
            this.DialogResult = DialogResult.Cancel;
            this.Close();
        }

        private void button2_Click(object sender, EventArgs e)
        {
            foreach(var item in listUp.SelectedIndices)
                upRegulated.Add((int)item);

            foreach (var item in listDown.SelectedIndices)
                downRegulated.Add((int)item);

            if (Double.TryParse(textBox1.Text, out threshold))
            {
                this.DialogResult = DialogResult.OK;
                this.Close();
            }
            else
                MessageBox.Show("Incorrect value for threshold. Must be (0,1].");

        }

        private void listUp_SelectedIndexChanged(object sender, EventArgs e)
        {
            foreach (var item in listUp.SelectedIndices)
                if (listDown.SelectedItems.Contains(item))
                    listDown.SelectedItems.Remove(item);
        }

        private void listDown_SelectedIndexChanged(object sender, EventArgs e)
        {
            foreach (var item in listDown.SelectedIndices)
                if (listUp.SelectedItems.Contains(item))
                    listUp.SelectedItems.Remove(item);
        }
    }
}
