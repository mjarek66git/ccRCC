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
    public partial class AddSuperGeneForm : Form
    {
        public HashSet<string> superGene = new HashSet<string>();
        BindingList<string> currentList = new BindingList<string>();
        public string name = "";
        Dictionary<string, int> genePos = new Dictionary<string, int>();
        Dictionary<string, BindingList<string>> sublists = new Dictionary<string, BindingList<string>>();
        BindingList<string> genes = new BindingList<string>();
        public AddSuperGeneForm(List<string> genes,string name=null,HashSet<string> superGene=null)
        {
            //this.genes = genes;
            InitializeComponent();
            if (name != null)
            {
                this.name = name;
                textBox1.Text = name;
            }
            foreach (var item in genes)
            {
                currentList.Add(item);
                for (int i = 1; i <= 3; i++)
                {
                    if (item.Length >= i)
                    {
                        string s = item.Substring(0, i);
                        if (!sublists.ContainsKey(s))
                            sublists.Add(s, new BindingList<string>());
                        sublists[s].Add(item);
                    }
                }
                this.genes.Add(item);
            }

            
            genesList.DataSource = currentList;

            for (int i = 0; i < genes.Count; i++)
                genePos.Add(genes[i], i);
            if (superGene != null)
            {
                this.superGene = superGene;
                foreach (var item in superGene)
                    superGeneList.Items.Add(item);
            }

        }

        private void button1_Click(object sender, EventArgs e)
        {
            if(genesList.SelectedIndex>=0)
            {
                superGeneList.Items.Add(genesList.SelectedItem);
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            if(superGeneList.SelectedIndex>=0)
            {
                superGeneList.Items.Remove(superGeneList.SelectedItem);
            }
        }

        private void button3_Click(object sender, EventArgs e)
        {
            DialogResult = DialogResult.Cancel;
            this.Close();
        }

        private void button4_Click(object sender, EventArgs e)
        {
            if (textBox1.Text.Length == 0)
            {
                MessageBox.Show("You have to provide super gene name");
                return;
            }
            name = textBox1.Text;
            foreach (var item in superGeneList.Items)
                superGene.Add((string)item);
            DialogResult = DialogResult.OK;
            this.Close();
        }

        private void textBox2_KeyPress(object sender, KeyPressEventArgs e)
        {
            
            //genesList.Refresh();
        }

        private void textBox2_KeyDown(object sender, KeyEventArgs e)
        {
            string lan = textBox2.Text;


            if (e.KeyCode == Keys.Back)
            {
                if (lan.Length > 0)
                    lan = lan.Remove(lan.Length - 1, 1);
            }
            else
                if (char.IsLetterOrDigit((char)e.KeyCode))
            {
                string key = e.KeyCode.ToString();
                if (key.Length > 1) key = key.Replace("NumPad", "").Replace("D", "");
                lan += key;
            }
            if (lan.Length == 0)
            {
                genesList.DataSource = genes;
                //currentList = genes;
                return;
            }
            if (lan.Length == 1 || lan.Length == 2 || lan.Length == 3)
            {
                currentList = sublists[lan];
                genesList.DataSource = currentList;
                //genesList.Refresh();
                return;
            }
            else
            {
                currentList=new BindingList<string>();
                genesList.DataSource = currentList;
                string subS = lan.Substring(0, 1);
                
                BindingList<string> list = sublists[subS];
                for (int i = 0; i < list.Count; i++)
                {
                    if (list[i].StartsWith(lan))
                        currentList.Add(list[i]);

                }
            }
        }
    }
}
