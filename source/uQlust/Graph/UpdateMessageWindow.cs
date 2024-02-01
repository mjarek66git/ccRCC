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
    public partial class UpdateMessageWindow : Form, MessageUpdate
    {
        public UpdateMessageWindow()
        {
            InitializeComponent();            
        }
        public void UpdateMessage(string message)
        {
            if (message != null)
            {
                this.Invoke((MethodInvoker)delegate
                {
                    label1.Text = message; // runs on UI thread
                });
            }
            else
            {
                this.Invoke((MethodInvoker)delegate
                {
                    this.Close();
                });
            }
        }
        public void ActivateUpdateting()
        {
            this.Show();
        }
        public void CloseUpdateting()
        {
            this.Close();
        }
    }
}
