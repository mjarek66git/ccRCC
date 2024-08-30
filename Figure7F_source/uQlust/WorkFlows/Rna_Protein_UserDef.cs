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
using phiClustCore.Profiles;
using Graph;

namespace WorkFlows
{
    public partial class Rna_Protein_UserDef : Form
    {
        bool previous = false;    
        Form parent;
        static public ResultWindow results = new ResultWindow();
        public Settings set = new Settings();
        public Rna_Protein_UserDef(Form parent=null)
        {
            InitializeComponent();
            set.Load();
            if(!Directory.Exists(set.profilesDir))
            {
                MessageBox.Show("Directory for storing profiles does not exist!\n You need to provide the new one.");
                FormSettings s = new FormSettings(set);
                s.hideCancel = false;
                s.ShowDialog();

            }
            this.parent = parent;
            InternalProfilesManager.RemoveProfilesFile();
            List<InternalProfileBase> profiles = InternalProfilesManager.InitProfiles();
        }     

        private void button1_Click(object sender, EventArgs e)
        {
            switch(((Button)sender).Text)
            {                
                case "User defined":
                    set.mode = INPUTMODE.USER_DEFINED;
                    break;
                case "Omics":
                case "HeatMap Omics":
                    set.mode = INPUTMODE.OMICS;
                    break;

            }
            //ClusteringChoose cluster = new ClusteringChoose(set,data,this);
            //cluster.Show();
            this.Hide();

        }

        private void Rna_Protein_UserDef_FormClosed(object sender, FormClosedEventArgs e)
        {
            results.Close();
            if (!previous && parent!=null)
                parent.Close();
            

        }

        private void button5_Click(object sender, EventArgs e)
        {
            UserOrOmics g = new UserOrOmics(this,OMICS_CHOOSE.NONE);
            g.Show();
            this.Hide();
        }

        private void button6_Click(object sender, EventArgs e)
        {
            //UserOrOmics h = new UserOrOmics(this, OMICS_CHOOSE.HEATMAP);
            Omics h = new Omics(this, OMICS_CHOOSE.HEATMAP);
            h.Show();
            this.Hide();
        }

        private void button1_Click_1(object sender, EventArgs e)
        {
            UserOrOmics g = new UserOrOmics(this,OMICS_CHOOSE.HNN);
            g.Show();
            this.Hide();

        }

        private void button2_Click(object sender, EventArgs e)
        {
            UserOrOmics g = new UserOrOmics(this,OMICS_CHOOSE.GUIDED_HASH);
            g.Show();
            this.Hide();
        }

        private void button3_Click(object sender, EventArgs e)
        {
            FormSettings s = new FormSettings(set);
            s.ShowDialog();
            set.Load();
        }


    }
}
