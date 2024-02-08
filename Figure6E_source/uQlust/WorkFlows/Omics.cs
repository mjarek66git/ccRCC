using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;
using System.Windows.Forms;
using phiClustCore;
using phiClustCore.Profiles;

namespace WorkFlows
{
    public enum OMICS_CHOOSE
    {
        HNN,
        GUIDED_HASH,
        HEATMAP,
        CLUSTERING,
        NONE
    }
    public partial class Omics : Form
    {
        Form parent;
        Form c = null;
        OmicsInput om = new OmicsInput();
        bool previous = false;
        public string processName = null;
        OMICS_CHOOSE nextWindow;
        int numValue;
        static int counter = 0;
        List<int> sampleLabelsPos = new List<int>();
        public Omics()
        {
            parent = null;
            InitializeComponent();
        }
        public Omics(Form parent,OMICS_CHOOSE nextWindow)
        {
            this.parent = parent;
            this.nextWindow = nextWindow;
            InitializeComponent();
            foreach (var item in Enum.GetValues(typeof(CodingAlg)))            
                comboBox1.Items.Add(item);
            comboBox1.SelectedIndex = 0;
            numValue = (int)numericUpDown3.Value;
            counter++;
        }

        private void button4_Click(object sender, EventArgs e)
        {
            previous = true;
            parent.Show();
            this.Close();
        }

        private void Genome_FormClosed(object sender, FormClosedEventArgs e)
        {
            if (!previous)
                parent.Close();    
        }
        private string GetProcessName()
        {
            return Path.GetFileNameWithoutExtension(textBox1.Text);
        }
        public void SetOptions()
        {
            numericUpDown2.Value=om.numCol;
            numericUpDown1.Value=om.numRow;
            checkBox4.Checked= om.uLabelGene;
            geneLabelPosition.Value = Convert.ToInt32(om.labelGeneStartString);
            textBox2.Text=om.labelSampleStartString ;
            checkBox5.Checked=om.uLabelSample;
            numericUpDown3.Value=om.numStates;
            checkBox1.Checked=om.transpose;
            comboBox1.SelectedItem=om.coding;
            radioButton3.Checked=om.genePosition;
            checkBox2.Checked=om.zScore;
            checkBox3.Checked=om.quantile;
            textBox3.Text=om.fileSelectedGenes;
            numericUpDown4.Value=om.selectGenes;

        }
        public void SaveOptions()
        {

            om.numCol =(int) numericUpDown2.Value;
            om.numRow = (int)numericUpDown1.Value;
            om.uLabelGene = checkBox4.Checked;
            om.labelGeneStartString = geneLabelPosition.Value.ToString();
            om.labelSampleStartString = textBox2.Text;
            om.uLabelSample = checkBox5.Checked;
            om.numStates = (int)numericUpDown3.Value;
            om.transpose = checkBox1.Checked;
            om.coding = (CodingAlg)comboBox1.SelectedItem;
            om.processName = GetProcessName();
            om.genePosition = radioButton3.Checked;
            om.zScore = checkBox2.Checked;
            om.quantile = checkBox3.Checked;
            if(textBox3.Text.Length>0)
                om.fileSelectedGenes = textBox3.Text;
            if (checkBox6.Checked)
                om.selectGenes = (int)numericUpDown4.Value;
        }
        void CreateWindow(Settings set)
        {
        }
        public virtual void button1_Click(object sender, EventArgs e)
        {
            if(textBox1.Text.Length==0 || textBox2.Text.Length==0)
            {
                this.DialogResult = DialogResult.Ignore;
                return;
            }


            SaveOptions();
            Settings set = new Settings();
            set.Load();
            set.mode = INPUTMODE.OMICS;
            CreateWindow(set);
            c.Show();
            this.Hide();
            counter++;
        }

        private void radioButton3_CheckedChanged(object sender, EventArgs e)
        {
            if (radioButton3.Checked)
            {
                label_samples.Text = label_genes.Text;
                label_genes.Text = "Labels in column";                
            }
            else
            {
                label_samples.Text = label_genes.Text;
                label_genes.Text = "Labels in row";
            }
        }

        private void checkBox4_CheckedChanged(object sender, EventArgs e)
        {
            if (checkBox4.Checked)
            {
                geneLabelPosition.Enabled = true;
                label_genes.Enabled = true;
            }
            else
            {
                geneLabelPosition.Enabled = false;
                label_genes.Enabled = false;
            }
        }

        private void checkBox5_CheckedChanged(object sender, EventArgs e)
        {
            if(checkBox5.Checked)
            {
                label_samples.Enabled = true;
                textBox2.Enabled = true;
            }
            else
            {
                label_samples.Enabled = false;
                textBox2.Enabled = false;
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            DialogResult res = openFileDialog1.ShowDialog();

            if (res == DialogResult.OK)
            {
                textBox1.Text = openFileDialog1.FileName;
                if(Path.GetExtension(textBox1.Text).Contains("gct"))
                {
                    StreamReader r= new StreamReader(textBox1.Text);
                    string line = r.ReadLine();
                    if(line.Contains("#1.3"))
                    {
                        line = r.ReadLine();
                        string[] aux = line.Split('\t');
                        if (aux.Length == 4)
                        {
                            int metaSamples = Convert.ToInt32(aux[3]);
                            numericUpDown2.Value = Convert.ToInt32(aux[2]) + 2;
                            numericUpDown1.Value = metaSamples + 2;
                            textBox2.Text = "";
                            for (int i = 0; i < metaSamples - 1; i++)
                                textBox2.Text += (i + 2) + ";";
                            textBox2.Text += (metaSamples + 1);
                        }
                    }
                    r.Close();
                }
            }

        }

        private void button3_Click(object sender, EventArgs e)
        {
            DialogResult res = openFileDialog1.ShowDialog();

            if (res == DialogResult.OK)
                textBox3.Text = openFileDialog1.FileName;

        }

        private void checkBox6_CheckedChanged(object sender, EventArgs e)
        {
            label6.Enabled = !checkBox6.Checked;
            textBox3.Enabled = !checkBox6.Checked;
            button3.Enabled = !checkBox6.Checked;
            numericUpDown4.Enabled = checkBox6.Checked;
        }

        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            if ((CodingAlg)comboBox1.SelectedItem == CodingAlg.Z_SCORE)
                if (numericUpDown3.Value < 3)
                    numericUpDown3.Value = 3;
                else
                    if (((int)numericUpDown3.Value) % 2 == 0)
                        numericUpDown3.Value = numericUpDown3.Value + 1;
        }

        private void numericUpDown3_ValueChanged(object sender, EventArgs e)
        {
            if ((CodingAlg)comboBox1.SelectedItem == CodingAlg.Z_SCORE)
                if (numericUpDown3.Value < 3)
                    numericUpDown3.Value = 3;
                else                    
                    if (((int)numericUpDown3.Value) % 2 == 0)
                        if(numValue<numericUpDown3.Value)
                            numericUpDown3.Value = numericUpDown3.Value + 1;
                        else
                            numericUpDown3.Value = numericUpDown3.Value - 1;

            numValue = (int)numericUpDown3.Value;

        }
    }
}
