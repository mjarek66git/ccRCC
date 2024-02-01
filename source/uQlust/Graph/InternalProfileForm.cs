using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using phiClustCore;
using phiClustCore.Profiles;


namespace Graph
{
    public partial class InternalProfileForm : Form
    {
        public profileNode localNode = new profileNode();
        public InternalProfileForm(profileNode node,filterOPT filter,bool flag=true)
        {
            InitializeComponent();
            Settings set = new Settings();
            set.Load();
            if (set.mode == INPUTMODE.USER_DEFINED)
                button5.Visible = false;
            if(filter==filterOPT.SIMILARITY)
            {
                button4.Visible = false;
                button3.Visible = true;
            }
            else
            {
                button4.Visible = true;
                button3.Visible = false;

            }
            if (node == null)
                return;
           
            if(node.internalName==UserDefinedProfile.GetName())
                dataGridView1.Columns["Column1"].ReadOnly = false;
            else
                dataGridView1.Columns["Column1"].ReadOnly = flag;
            localNode.profProgram = node.profProgram;
            localNode.internalName = node.internalName;
            textBox1.Text = node.profName;
            if (node.GetNumberofStates()> 0)
            {
                dataGridView1.Rows.Add(node.GetNumberofStates());
                int i = 0;
                foreach (var item in node.states)
                {
                    dataGridView1.Rows[i].Cells[0].Value = item.Key;
                    dataGridView1.Rows[i++].Cells[1].Value = item.Value;
                }


                foreach (var item in node.profWeights.Keys)
                {
                    foreach (var item2 in node.profWeights[item].Keys)
                    {
                        dataGridView2.Rows.Add(1);
                        int index = dataGridView2.Rows.Count - 2;
                        dataGridView2.Rows[index].Cells[0].Value = item;
                        dataGridView2.Rows[index].Cells[1].Value = item2;
                        dataGridView2.Rows[index].Cells[2].Value = node.profWeights[item][item2].ToString();
                    }
                }
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            localNode.profName = textBox1.Text;
            for (int i = 0; i < dataGridView1.Rows.Count;i++)
                if (dataGridView1.Rows[i].Cells[0].Value != null && dataGridView1.Rows[i].Cells[1].Value != null)
                    localNode.AddStateItem((string)dataGridView1.Rows[i].Cells[0].Value, (string)dataGridView1.Rows[i].Cells[1].Value);

            localNode.profWeights.Clear();
            for (int i = 0; i < dataGridView2.Rows.Count; i++)
            {
                if (dataGridView2.Rows[i].Cells[0].Value != null && dataGridView2.Rows[i].Cells[1].Value != null)
                {
                    string item1 = (string)dataGridView2.Rows[i].Cells[0].Value;
                    string item2 = (string)dataGridView2.Rows[i].Cells[1].Value;
                    if (localNode.profWeights.ContainsKey(item1))
                    {
                        if (localNode.profWeights[item1].ContainsKey(item2))
                            localNode.profWeights[item1].Remove(item2);
                        localNode.profWeights[item1].Add(item2, Convert.ToDouble(dataGridView2.Rows[i].Cells[2].Value.ToString()));
                    }
                    else
                    {
                        SerializableDictionary<string, double> ww = new SerializableDictionary<string, double>();
                        ww.Add(item2, Convert.ToDouble(dataGridView2.Rows[i].Cells[2].Value.ToString()));
                        localNode.profWeights.Add(item1, ww);
                    }
                }
            }
            this.DialogResult = DialogResult.OK;
            this.Close();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            this.DialogResult = DialogResult.Cancel;

            this.Close();
        }

        private void GenerateWeights(bool similarity)
        {
            Dictionary<string, int> states = new Dictionary<string, int>();
            foreach(DataGridViewRow item in dataGridView1.Rows)
            {
                if (item.Cells[1].Value != null && !states.ContainsKey((string)item.Cells[1].Value))
                    states.Add(item.Cells[1].Value.ToString(),0);
            }
            dataGridView2.Rows.Clear();
            if (similarity)
            {
                foreach (var item in states)
                {
                    dataGridView2.Rows.Add(1);
                    int num = dataGridView2.Rows.Count - 2;
                    dataGridView2.Rows[num].Cells[0].Value = item.Key;
                    dataGridView2.Rows[num].Cells[1].Value = item.Key;
                    dataGridView2.Rows[num].Cells[2].Value = 1;
                }
            }
            else
            {
                foreach (var item in states.Keys)
                {
                    foreach (var item2 in states.Keys)
                    {
                        if (item!= item2)
                        {
                            dataGridView2.Rows.Add(1);
                            int num = dataGridView2.Rows.Count - 2;
                            dataGridView2.Rows[num].Cells[0].Value = item;
                            dataGridView2.Rows[num].Cells[1].Value = item2;
                            dataGridView2.Rows[num].Cells[2].Value = 1;
                        }
                    }
                }

            }

        }

        private void button3_Click(object sender, EventArgs e)
        {
            GenerateWeights(true);
        }

        private void button4_Click(object sender, EventArgs e)
        {
            GenerateWeights(false);
        }

        private void button5_Click(object sender, EventArgs e)
        {
            DialogResult res = openFileDialog1.ShowDialog();
            if(res==DialogResult.OK)
            {
                localNode.profName = textBox1.Text;
                ProfileTree t = ProfileAutomatic.AnalyseProfileFile(openFileDialog1.FileName,SIMDIST.SIMILARITY);

                if(t.masterNode.ContainsKey(localNode.profName))
                {
                    dataGridView1.Rows.Clear();
                    foreach(var item in t.masterNode[localNode.profName].states.Keys)
                    {
                        dataGridView1.Rows.Add(1);
                        dataGridView1.Rows[dataGridView1.Rows.Count - 2].Cells[0].Value = item;
                        dataGridView1.Rows[dataGridView1.Rows.Count - 2].Cells[1].Value = item;
                    }
                }

            }
           
        }
    }
}
