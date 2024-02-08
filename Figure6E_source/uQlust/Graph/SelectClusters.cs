using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace Graph
{
    public partial class SelectClusters : Form
    {
        ClusterV.HClusterNode hNode;
        string name;
        string measure;
        public SelectionOpt sel;

        public SelectClusters(ClusterV.HClusterNode hNode,string name,string measure,SelectionOpt sel)
        {
            this.sel = sel;
            InitializeComponent();

            numericUpDown1.Value = sel.start;
            numericUpDown2.Value = sel.stop;

            if (sel.Name != "Not Defined")
            {
                radioButton1.Checked = sel.range;
                radioButton2.Checked = !sel.range;
            }
            this.hNode = hNode;
            this.name = name;
            this.measure = measure;
        }

        private void Enable(bool flag)
        {
            label2.Enabled = flag;           
            numericUpDown1.Enabled = flag;
            numericUpDown2.Enabled = flag;
            label1.Enabled = flag;
            button1.Enabled = !flag;
            label3.Enabled = !flag;
            label4.Enabled = !flag;

        }

        private void radioButton1_CheckedChanged(object sender, EventArgs e)
        {
            Enable(radioButton1.Checked);
        }

        private void radioButton2_CheckedChanged(object sender, EventArgs e)
        {
            Enable(!radioButton2.Checked);
        }

        private void button1_Click(object sender, EventArgs e)
        {
            visHierar win = new visHierar(hNode, name,measure);
            win.ShowCloseButton();
            win.ShowDialog();
            sel.start = Convert.ToInt32(numericUpDown1.Value);
            sel.stop = Convert.ToInt32(numericUpDown2.Value);
            sel.range = radioButton1.Checked;
            if (win.listNodes != null)
            {
                if (sel.clusters != null)
                    sel.clusters.Clear();
                else
                    sel.clusters = new List<List<string>>();
                label4.Text = win.listNodes.Count.ToString();
                foreach (var item in win.listNodes)
                {
                    List<string> cl = new List<string>(item.Key.setStruct);
                    sel.clusters.Add(cl);
                }
            }
            else
                label4.Text = "0";
           
           
        }

        private void button2_Click(object sender, EventArgs e)
        {
            this.DialogResult = DialogResult.OK;
            this.Close();
        }
    }
}
