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
using phiClustCore;

namespace WorkFlows
{
    public partial class UserOrOmics : Form
    {
        Form parent;
        OmicsInput om = new OmicsInput();
        bool previous = false;
        OMICS_CHOOSE prevWindow;
        public UserOrOmics(Form parent,OMICS_CHOOSE prevWindow)
        {
            this.parent = parent;
            this.prevWindow = prevWindow;
            InitializeComponent();         
        }

        private void button5_Click(object sender, EventArgs e)
        {
            Omics g = new Omics(this, prevWindow);
            g.Show();
            this.Hide();
        }
        string GetProcessName(object o,INPUTMODE mode)
        {
            return "WorkFlow_" + mode.ToString() + "_" + o.ToString();
        }
        private void button3_Click(object sender, EventArgs e)
        {
        }

        private void button4_Click(object sender, EventArgs e)
        {
            previous = true;
            parent.Show();
            this.Close();
        }

        private void UserOrOmics_FormClosed(object sender, FormClosedEventArgs e)
        {
            if (!previous)
                parent.Close();
        }
    }
}
