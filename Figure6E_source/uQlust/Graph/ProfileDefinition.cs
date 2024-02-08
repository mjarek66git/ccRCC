using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using phiClustCore;

namespace Graph
{
    public partial class ProfileDefinitionForm : Form
    {
        public profileNode profile;

        public ProfileDefinitionForm()
        {
            InitializeComponent();
            profile = new profileNode();
        }
        public ProfileDefinitionForm(profileNode profile)
        {
            InitializeComponent();
            this.profile = profile;
            textBox1.Text = profile.profName;
            textBox2.Text = profile.profProgram;
            textBox3.Text = profile.OutFileName;
            if (profile.progParameters == null)
                profile.progParameters = "input_file";
            textBox4.Text = profile.progParameters;
            checkBox1.Checked = profile.removeOutFile;
            foreach (var item in profile.profWeights.Keys)
            {
                foreach (var item2 in profile.profWeights[item].Keys)
                {
                    dataGridView1.Rows.Add();
                    dataGridView1.Rows[dataGridView1.Rows.Count-2].Cells[0].Value = item;
                    dataGridView1.Rows[dataGridView1.Rows.Count-2].Cells[1].Value = item2;
                    dataGridView1.Rows[dataGridView1.Rows.Count-2].Cells[2].Value = profile.profWeights[item][item2].ToString();
                }
            }

        }
        private void OkBtn_Click(object sender, EventArgs e)
        {

            if (textBox3.Text.Length == 0)
            {
                this.DialogResult = DialogResult.None;
                MessageBox.Show("Output filename must be provided");
                return;
            }
            if (textBox1.Text.Length == 0)
            {
                this.DialogResult = DialogResult.None;
                MessageBox.Show("Name of the profile must be provided");
                return;

            }

            profile.profName = textBox1.Text;
            profile.profProgram = textBox2.Text;
            profile.OutFileName = textBox3.Text;
            profile.removeOutFile = checkBox1.Checked;
            profile.progParameters =textBox4.Text;

            this.DialogResult = DialogResult.OK;

            profile.profWeights.Clear();
            for (int i = 0; i < dataGridView1.Rows.Count; i++)
            {
                if (dataGridView1.Rows[i].Cells[0].Value != null && dataGridView1.Rows[i].Cells[0].Value != null && dataGridView1.Rows[i].Cells[1].Value != null)
                {
                    string item1 = (string)dataGridView1.Rows[i].Cells[0].Value;
                    string item2 = (string)dataGridView1.Rows[i].Cells[1].Value;
                    if (profile.profWeights.ContainsKey(item1))
                    {
                        if (profile.profWeights[item1].ContainsKey(item2))
                            profile.profWeights[item1].Remove(item2);
                        profile.profWeights[item1].Add(item2, Convert.ToDouble((string)dataGridView1.Rows[i].Cells[2].Value));
                    }
                    else
                    {
                        SerializableDictionary<string, double> ww = new SerializableDictionary<string, double>();
                        ww.Add(item2, Convert.ToDouble((string)dataGridView1.Rows[i].Cells[2].Value));
                        profile.profWeights.Add(item1, ww);
                    }
                }

            }

            this.Close();
        }

        private string GetFileName(string currentString)
        {
            DialogResult res;
            res = openFileDialog1.ShowDialog();
            if (res == DialogResult.OK)
                return openFileDialog1.FileName;
            
            return currentString;
        }

        private void button1_Click(object sender, EventArgs e)
        {
            textBox2.Text = GetFileName(textBox2.Text);
        }

        private void button2_Click(object sender, EventArgs e)
        {
            textBox3.Text = GetFileName(textBox3.Text);
        }

        private void checkBox1_CheckedChanged(object sender, EventArgs e)
        {
            if (checkBox1.Checked)
            {
                button2.Visible = true;
           }
            else
            {
                button2.Visible = false;
            }
        }

        private void dataGridView1_KeyPress(object sender, KeyPressEventArgs e)
        {
        }

        private void ProfileDefinitionForm_Load(object sender, EventArgs e)
        {

        }        
    }
}
