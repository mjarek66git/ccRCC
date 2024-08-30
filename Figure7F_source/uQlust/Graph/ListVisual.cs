using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO;
using phiClustCore;
using phiClustCore.Interface;


namespace Graph
{
    public delegate void ClosingForm(string item);

    public partial class ListVisual : Form,IVisual
    {
        List<List<string>> clusters;
        string selectedItem = "";
        TextInput input = null;
        Dictionary<string, string> labels = null;
        public ClosingForm closeForm=null;
        ClusterOutput output;
        public ListVisual(ClusterOutput output,string item,Dictionary <string,string> labels)
        {
            InitializeComponent();
            this.output = output;
            this.clusters = output.clusters.list;
            for (int i = 1; i <= clusters.Count; i++)
            {
                if (clusters[i - 1].Count>1 && output.clusters.consistency!=null)
                    listBox1.Items.Add(String.Format("{0,12} {1,7} {2,8}", "Cluster_"+i ,clusters[i - 1].Count,output.clusters.consistency[i-1].ToString("0.00")));
                else
                    listBox1.Items.Add(String.Format("{0,12} {1,7} ","Cluster_"+ i, clusters[i - 1].Count));
            }
            this.Text = item;
            this.labels = labels;
            richTextBox1.SelectAll();
            if(listBox1.Items.Count>0)
                listBox1.SelectedIndex=0;
            //richTextBox1.SelectionColor = Color.Black;
         
        }
        public override string ToString()
        {
            return "List";
        }
        public void ToFront()
        {
            this.BringToFront();
        }
        private void listBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            richTextBox1.Text="";
            if (listBox1.SelectedIndex == -1)
                return;
           
            int size = 0;
            for (int i = 0; i < clusters[listBox1.SelectedIndex].Count; i++)
                size += clusters[listBox1.SelectedIndex][i].Length;

            StringBuilder st = new StringBuilder(size);
            int remLine=-1;
            for (int i = 0; i < clusters[listBox1.SelectedIndex].Count; i++)
            {
                string line = clusters[listBox1.SelectedIndex][i];

                if (labels != null)
                    if (labels.ContainsKey(line))
                        line += "\t" + labels[clusters[listBox1.SelectedIndex][i]];
                    
                st.AppendLine(line);
                if (selectedItem.Length > 0)
                    if (clusters[listBox1.SelectedIndex][i].Equals(selectedItem))
                        remLine = i;
            }
            richTextBox1.Text=st.ToString();
            if (remLine >= 0)
            {
                int startIndex = richTextBox1.GetFirstCharIndexFromLine(remLine);
                richTextBox1.Select(startIndex, selectedItem.Length);
                richTextBox1.ScrollToCaret();
                //richTextBox1.Select(1,10);
                richTextBox1.SelectionColor = System.Drawing.Color.White;
                richTextBox1.SelectionBackColor = System.Drawing.Color.Blue;
            }

        }

        private void BakerVisual_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (closeForm != null)
                closeForm(this.Text);
        }
        protected override bool ProcessCmdKey(ref Message msg, Keys keyData)
        {
            switch (keyData)
            {
                case Keys.Control | Keys.F:
                    {
                        if (input == null)
                        {
                            input = new TextInput("NEXT");
                            input.textBox = richTextBox1;
                            input.Show();
                            
                        }
                        input.BringToFront();
                        
                        return true;
                    }
            }
            return base.ProcessCmdKey(ref msg, keyData);
        }

        private void textBox1_TextChanged(object sender, EventArgs e)
        {
        }

        private void textBox1_Enter(object sender, EventArgs e)
        {

        }

        private void textBox1_KeyDown(object sender, KeyEventArgs e)
        {
            if (e.KeyCode == Keys.Enter)
            {
                for (int i = 0; i < clusters.Count; i++)
                {
                    for (int j = 0; j < clusters[i].Count; j++)
                        if (clusters[i][j].Equals(textBox1.Text))
                        {
                            selectedItem = textBox1.Text;
                            listBox1.SelectedIndex = i;                            
                            return;
                        }
                }
                selectedItem = "";
            }
        }

        private void button1_Click(object sender, EventArgs e)
        {
            DialogResult res = saveFileDialog1.ShowDialog();

            if(res==DialogResult.OK)
            {
                output.SaveTxt(saveFileDialog1.FileName);
            }
        }
        void DrawClusterData(Dictionary<string, string[]> data,int classNum)
        {
            double[] accCluster = new double[clusters.Count];
            for (int i = 0; i < clusters.Count; i++)
            {
                Dictionary<string, double> res = new Dictionary<string, double>();
                for (int j = 0; j < clusters[i].Count; j++)
                {
                    if (data.ContainsKey(clusters[i][j]))
                        if (res.ContainsKey(data[clusters[i][j]][classNum]))
                            res[data[clusters[i][j]][classNum]]++;
                        else
                            res.Add(data[clusters[i][j]][classNum], 1.0);
                }

                List<KeyValuePair<string, double>> xx = new List<KeyValuePair<string, double>>();
                xx = res.OrderBy(key => key.Value).ToList();
                double w= 0;
                foreach (var item in xx)
                {
                    w += item.Value;
                    // do something with item.Key and item.Value
                }
                accCluster[i] = xx[0].Value / w;
            }
        }
        private void button2_Click(object sender, EventArgs e)
        {
            DialogResult res = openFileDialog1.ShowDialog();

            if (res == DialogResult.OK)
            {
                StreamReader r = new StreamReader(openFileDialog1.FileName);
                Dictionary<string, string[]> data = new Dictionary<string, string[]>();
                string line = r.ReadLine();
                while(line!=null)
                {
                    string[] aux = line.Split(' ');
                    string[] classLabes = new string[aux.Length - 1];
                    for (int i = 1; i < aux.Length; i++)
                        classLabes[i - 1] = aux[i];
                    data.Add(aux[0], classLabes);
                    line = r.ReadLine();
                }
                r.Close();
                DrawClusterData(data, 0);
                DrawPanel pn = new DrawPanel("Leave  ",null);
                pn.Height = 500;
                pn.Width = 400;                

                pn.CreateBMP();
                //pn.drawPic = delegate { double res = DrawNodeProfiles(pn.bmp, clickNode.setProfiles, clickNode.consistency); pn.Text = res.ToString(); };
                pn.Show();
                pn.pictureBox1.Invalidate();

            }


        }
    }
}
