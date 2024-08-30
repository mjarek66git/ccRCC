using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;
using phiClustCore.Interface;

namespace Graph
{
    public partial class ResultsForm : Form,IShowResults
    {
        DataTable org;
        DataTable t;
        public ResultsForm()
        {
            InitializeComponent();
        }

        private void button2_Click(object sender, EventArgs e)
        {
            this.Close();
        }
        private void Show(KeyValuePair<string, DataTable> t)
        {
            this.Text = t.Key;
            this.t=t.Value.Copy();
            this.org = t.Value.Copy();
            dataGridView1.RowHeadersVisible = false;
            dataGridView1.DataSource = t.Value;
            for (int i = 0; i < dataGridView1.Columns.Count; i++)
                dataGridView1.Columns[i].AutoSizeMode = DataGridViewAutoSizeColumnMode.Fill;
            this.Show();

        }
        public void Show(List<KeyValuePair<string,DataTable>> t)
        {
            if (t == null)
                return;
            Show(t[0]);

            for (int i = 1; i < t.Count; i++)
            {
                ResultsForm f = new ResultsForm();
                f.Show(t[i]);
            }
            

        }
        public void ShowException(Exception ex)
        {
            if (ex != null)
                MessageBox.Show("Exception: " + ex.Message);
        }

        private void textBox1_KeyUp(object sender, KeyEventArgs e)
        {
            if (textBox1.Text.Length > 0)
            {
                dataGridView1.DataSource = null;
                t.Clear(); 
               
                for (int j = 0; j < org.Rows.Count; j++)
                    for (int i = 0; i < org.Rows[j].ItemArray.Length; i++)
                        if (org.Rows[j].ItemArray[i].ToString().StartsWith(textBox1.Text))
                        {
                            t.ImportRow(org.Rows[j]);
                            break;
                        }

            }
            else
                t = org.Copy();
            dataGridView1.DataSource = t;
            for (int i = 0; i < dataGridView1.Columns.Count; i++)
                dataGridView1.Columns[i].AutoSizeMode = DataGridViewAutoSizeColumnMode.Fill;
            dataGridView1.Update();

        }

        private void button1_Click(object sender, EventArgs e)
        {
            DialogResult res;
            saveFileDialog1.Filter = "out|*.out";
            res=saveFileDialog1.ShowDialog();
            if(res==DialogResult.OK)
            {
                List<int> columnSize = new List<int>();
                for (int i = 0; i < dataGridView1.Columns.Count; i++)                
                   columnSize.Add(dataGridView1.Columns[i].HeaderText.Length);
                
                for(int i=0;i<dataGridView1.Rows.Count;i++)
                {
                    for (int j = 0; j < dataGridView1.Rows[i].Cells.Count; j++)
                        if (dataGridView1.Rows[i].Cells[j].Value!=null)
                            if (columnSize[j] < dataGridView1.Rows[i].Cells[j].Value.ToString().Length)
                                columnSize[j] = dataGridView1.Rows[i].Cells[j].Value.ToString().Length;
                }

                using(StreamWriter wr = new StreamWriter(saveFileDialog1.FileName))
                {
                    for (int i = 0; i < dataGridView1.Columns.Count; i++)
                    {
                        wr.Write("{0,"+columnSize[i]+"}\t", dataGridView1.Columns[i].HeaderText);
                    }
                    wr.WriteLine();
                    for(int i=0;i<dataGridView1.Rows.Count;i++)
                    {
                        for (int j = 0; j < dataGridView1.Rows[i].Cells.Count; j++)
                            wr.Write("{0," + columnSize[j] + "}\t", dataGridView1.Rows[i].Cells[j].Value);
                        wr.WriteLine();
                    }
                }
            }
        }
    }
}
