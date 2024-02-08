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
    public partial class SuperGenesForm : Form
    {
        List<string> genes;
        public string fileName;
        public Dictionary<string, HashSet<string>> superGenes = new Dictionary<string, HashSet<string>>();
        Dictionary<string, int> genePos = new Dictionary<string, int>();
        public SuperGenesForm(List <string> genes,Dictionary<string,HashSet<string>> superGenes,bool selection=false)
        {
            InitializeComponent();
            this.genes = genes;
            if (superGenes != null)
            {
                this.superGenes = superGenes;
                SuperGenesToList();
                saveB.Visible = true;
            }
            if(selection)
            {
                this.superGenes = null;
                saveB.Visible = false;
                button1.Visible = false;
                button6.Visible = false;
                button2.Visible = false;
                GenesToList();

            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            AddSuperGeneForm super = new AddSuperGeneForm(genes);
            DialogResult res = super.ShowDialog();
            if(res==DialogResult.OK)
            {
                if(super.superGene.Count>0)
                {
                    superGenes.Add(super.name, super.superGene);
                    SuperGenesToList();
                    saveB.Visible = true;
                }
            }
        }
        public List<int> GetSelectedGenes()
        {
            List<int> selected = new List<int>();
            foreach (var item in listBox1.Items)
                selected.Add(genePos[item.ToString()]);
            return selected;
        }
        private void button3_Click(object sender, EventArgs e)
        {
            List<string> toRemove = new List<string>();

            if (listBox1.SelectedIndex >= 0)
            {
                foreach (var item in listBox1.SelectedItems)
                {
                    string it = (string)item;
                    toRemove.Add(it);
                    //listBox1.Items.Remove(item);

                                        
                    if (superGenes?.Count == 0)
                        saveB.Visible = false;
                }
                foreach (var item in toRemove)
                {
                    listBox1.Items.Remove(item);
                    string[] aux = item.Split(':');
                    superGenes?.Remove(aux[0]);
                }
            }

        }
        void SaveSuperGenes(string fileName)
        {
            StreamWriter sw = new StreamWriter(fileName);

            foreach(var item in superGenes)
            {
                sw.WriteLine(">" + item.Key);
                for(int i=0;i<item.Value.Count;i++)
                {
                    sw.Write(item.Value.ElementAt(i) + " ");
                }
                sw.WriteLine(item.Value.ElementAt(item.Value.Count - 1));
            }

            sw.Close();

        }
        Dictionary<string, HashSet<string>> ReadSuperGenes(string fileName)
        {
            StreamReader sr = new StreamReader(fileName);

            Dictionary<string, int> hashGeneNames = new Dictionary<string, int>();

            for (int i = 0; i < genes.Count; i++)
                if (hashGeneNames.ContainsKey(genes[i]))
                    Console.WriteLine(genes[i]);
                else
                    hashGeneNames.Add(genes[i], i);

            Dictionary<int, int> check = new Dictionary<int, int>();

            string superGeneName = "";
            string line = sr.ReadLine();
            while (line != null)
            {
                if (line.StartsWith(">"))
                    superGeneName = line.Substring(1);
                else
                {
                    string[] aux = line.Split(' ');
                    List<string> l = aux.ToList();
                    HashSet<string> genesP = new HashSet<string>();
                    foreach (var item in l)
                        if (hashGeneNames.ContainsKey(item))
                            genesP.Add(item);
                       /* else
                            throw new Exception("Uknown gene: " + item);*/
                    if(genesP.Count>0 && !superGenes.ContainsKey(superGeneName))
                        superGenes.Add(superGeneName, genesP);
                }

                line = sr.ReadLine();
            }
            sr.Close();

            return superGenes;
        }
        void GenesToList()
        {
            listBox1.Items.Clear();
            genePos.Clear();
            for(int i=0;i<genes.Count;i++)
            {
                listBox1.Items.Add(genes[i]);
                genePos.Add(genes[i], i);
            }
        }
        void SuperGenesToList()
        {
            listBox1.Items.Clear();
            foreach (var item in superGenes)
            {
                string gene = item.Key + ": ";
                foreach (var it in item.Value)
                    gene += it + " ";

                listBox1.Items.Add(gene);
            }

        }
        private void button1_Click(object sender, EventArgs e)
        {
            DialogResult res=openFileDialog1.ShowDialog();

            if(res==DialogResult.OK)
            {
                ReadSuperGenes(openFileDialog1.FileName);
                fileName = openFileDialog1.FileName;
                SuperGenesToList();
                if (superGenes.Count > 0)
                    saveB.Visible = true;
            }
        }

        private void button5_Click(object sender, EventArgs e)
        {
            DialogResult = DialogResult.Cancel;
            this.Close();
        }

        private void button4_Click(object sender, EventArgs e)
        {
            DialogResult = DialogResult.OK;
            Close();
        }

        private void button6_Click(object sender, EventArgs e)
        {
            DialogResult res=saveFileDialog1.ShowDialog();
            if(res==DialogResult.OK)
            {
                SaveSuperGenes(saveFileDialog1.FileName);
            }
        }

        private void button6_Click_1(object sender, EventArgs e)
        {
            superGenes.Clear();
            listBox1.Items.Clear();
        }

        private void button7_Click(object sender, EventArgs e)
        {
            if (listBox1.SelectedIndex >= 0)
            {
                HashSet<string> xx = new HashSet<string>();
                HashSet<string> yy = new HashSet<string>();
                foreach (var item in listBox1.SelectedItems)
                {
                    xx.Add((string)item);
                    string[] aux = ((string)item).Split(':');
                    yy.Add(aux[0]);
                }

                listBox1.Items.Clear();
                foreach(var item in xx)
                    listBox1.Items.Add(item);

                if (superGenes != null)
                {
                    Dictionary<string, HashSet<string>> selected = new Dictionary<string, HashSet<string>>();

                    foreach (var item in superGenes)
                    {
                        if (yy.Contains(item.Key))
                            selected.Add(item.Key, item.Value);
                    }
                    superGenes = selected;
                }
            }

        }
    }
}
