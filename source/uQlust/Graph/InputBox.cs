using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace phiClustCore.Graph
{
    public partial class InputBox : Form
    {
        public string stringValue = "";
        public InputBox(string message)
        {
            InitializeComponent();
            label2.Text = message;
        }

        private void button2_Click(object sender, EventArgs e)
        {
            DialogResult = DialogResult.Cancel;
            this.Close();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            stringValue = textBox1.Text;
            DialogResult = DialogResult.OK;
            this.Close();
        }
    }
}
