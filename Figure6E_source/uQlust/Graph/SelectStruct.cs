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
    public partial class SelectStruct : Form
    {
        List<string> structures;
        public string selectedStruct=null;
        public SelectStruct(List <string> structures)
        {
            InitializeComponent();
            this.structures = structures;
            if (structures.Count > 0)
            {
                listBox1.BeginUpdate();
                foreach (var item in structures)
                    listBox1.Items.Add(item);

                listBox1.EndUpdate();
            }
            
        }

        private void textBox1_KeyUp(object sender, KeyEventArgs e)
        {
            listBox1.BeginUpdate();
            listBox1.Items.Clear();
            if (textBox1.Text.Length > 0)
            {
                
                foreach (var item in structures)
                    if(item.StartsWith(textBox1.Text))                    
                        listBox1.Items.Add(item);
                
            }
            else
                foreach (var item in structures)
                    listBox1.Items.Add(item);

            listBox1.EndUpdate();
        }

        private void saveBtn_Click(object sender, EventArgs e)
        {
            if (listBox1.Items.Count == 1)
            {
                selectedStruct = (string)listBox1.Items[0];
                this.DialogResult = DialogResult.OK;
                this.Close();
            }
            else
            {
                if (listBox1.SelectedItem != null)
                {
                    selectedStruct = (string)listBox1.SelectedItem;
                    this.DialogResult = DialogResult.OK;
                    this.Close();
             
                }
            }
        }

        private void listBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            textBox1.Text = (string)listBox1.SelectedItem;
        }

    }
}
