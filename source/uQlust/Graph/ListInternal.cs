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
    public partial class ListInternal : Form
    {
        public string selectedProfile = "";
        public ListInternal(List <profileNode> internalP)
        {
            InitializeComponent();
            foreach (var item in internalP)                
                listBox1.Items.Add(item.profName);            

        }

        private void button2_Click(object sender, EventArgs e)
        {
            selectedProfile = (string)listBox1.SelectedItem;
            this.DialogResult = DialogResult.OK;
            this.Close();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            this.DialogResult = DialogResult.Cancel;
            this.Close();
        }
    }
}
