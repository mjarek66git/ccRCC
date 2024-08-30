using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using phiClustCore.Interface;


namespace Graph
{
    public partial class Progress : Form
    {
        IProgressBar progress=null;
        IShowResults show=null;
        public Progress(IProgressBar progress,IShowResults res)
        {
            InitializeComponent();
            this.progress = progress;
            this.show = res;
        }
        public void Start()
        {
            timer1.Start();            
        }

        private void timer1_Tick(object sender, EventArgs e)
        {
            if (progress == null)
                return;
            Exception ex = progress.GetException();
            if(progress.GetException()!=null)
            {
                show.ShowException(ex);
                this.Close();
               
            }
            progressBar1.Value =(int)(progress.ProgressUpdate()*100);
            if (progressBar1.Value == progressBar1.Maximum)
            {
                if(show!=null)
                    show.Show(progress.GetResults());
                Close();
            }
            
        }
    }
}
