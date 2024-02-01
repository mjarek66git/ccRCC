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
    public partial class GetProcessName : Form
    {
        public string name;
        public GetProcessName(string name)
        {
            InitializeComponent();
            this.name = name;
        }

        private void OkBtn_Click(object sender, EventArgs e)
        {
            name = textBox1.Text;
        }
    }
}
