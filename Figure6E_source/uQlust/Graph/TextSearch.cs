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
    public delegate void RunSearch(string x,bool backward);
    public partial class TextInput : Form
    {
        public RunSearch run=null;
        public RichTextBox textBox=null;
        int currentPosition = 0;
        public TextInput(string buttonStr)
        {
            InitializeComponent();
            OKBtn.Text = buttonStr;
        }

        private void CancelBtn_Click(object sender, EventArgs e)
        {
            this.Close();
        }

        private void OKBtn_Click(object sender, EventArgs e)
        {
            FindAndHighlight(textBox1.Text, checkBox1.Checked);
        }

        private void textBox1_KeyDown(object sender, KeyEventArgs e)
        {
            if (e.KeyCode == Keys.Enter)
            {
                FindAndHighlight(textBox1.Text, checkBox1.Checked);
            }
        }
        private void FindAndHighlight(string str, bool backward)
        {
            int index;

            if (backward)
                index = textBox.Find(str, 0,currentPosition, RichTextBoxFinds.Reverse);
            else
                index = textBox.Find(str, currentPosition,textBox.Text.Length, RichTextBoxFinds.MatchCase);
            
            if (index >= 0)
            {
                textBox.Select(index, str.Length);
                if (backward)
                    currentPosition = index;
                else
                    currentPosition = index + str.Length;
                textBox.Focus(); ;
            }
        }

    }
}
