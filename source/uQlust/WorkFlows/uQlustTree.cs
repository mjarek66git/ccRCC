using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using phiClustCore;
using phiClustCore.Profiles;

namespace WorkFlows
{
    public partial class uQlustTreeSimple : Form, IclusterType
    {
        protected Options opt = new Options();
        protected ResultWindow results;
        Settings set;
        public string processName;
        bool previous = false;
        static int counter = 0;
        CommonDialog dialog;
        Form parent;
        string dataFileName = "";
        OmicsDataSet data;
        ProfileTree tree = new ProfileTree();

        public uQlustTreeSimple(Form parent, OmicsDataSet data, Settings set,ResultWindow results)
        {
            InitializeComponent();
            this.data = data;

            distanceControl2.hideReference = true;
            distanceControl2.hideSetup = true;
            this.parent = parent;
            this.Location = parent.Location;
            dialog = folderBrowserDialog1;
            if (set.mode == INPUTMODE.USER_DEFINED || set.mode==INPUTMODE.OMICS)
            {
                dialog = openFileDialog1;
            }
          
            this.set = set; 
           
            numericUpDown1.Value = opt.hash.refPoints;
            this.results = results;
        }
        public override string ToString()
        {
            return "uQlustTree";
        }
        public INPUTMODE GetInputType()
        {
            return set.mode;
        }
        void SetProfileOptions()
        {
            relevantC.Value = opt.hash.relClusters;
            //distanceControl2.distDef = opt.hierarchical.distance;
            if (opt.hash.combine)
                radioButton1.Checked = true;
            else            
                Hash.Checked = true;
            
            switch (set.mode)
            {
                case INPUTMODE.USER_DEFINED:
                case INPUTMODE.OMICS:
                    //distanceControl2.distDef = DistanceMeasures.COSINE;
                    break;

            }
            if (opt.hash.profileName != null)
            {
                tree.LoadProfiles(opt.hash.profileName);
            }
        }
        public void SetProfileName(string name)
        {
            opt.ReadOptionFile(name);
            SetProfileOptions();
        }
        private void button2_Click(object sender, EventArgs e)
        {
            opt.dataDir.Clear();
            opt.profileFiles.Clear();

            opt.hash.relClusters = (int)relevantC.Value;
            opt.hash.perData = 90;
            if (radioButton1.Checked)
                opt.hash.combine = true;
            else
                opt.hash.combine = false;
            opt.clusterAlgorithm.Clear();
            /*if (!radioHTree.Checked)
                opt.clusterAlgorithm.Add(ClusterAlgorithm.uQlustTree);
            else
                opt.clusterAlgorithm.Add(ClusterAlgorithm.HTree);*/
           // opt.hierarchical.distance = distanceControl2.distDef;
           
           // opt.hierarchical.microCluster = radioHTree.Checked;
            opt.hash.refPoints = (int)numericUpDown1.Value;


            opt.hierarchical.distance = distanceControl2.distDef;

            if (opt.hierarchical.distance == DistanceMeasures.HAMMING || opt.hierarchical.distance==DistanceMeasures.COSINE)
                opt.hierarchical.reference1DjuryH = true;
            else
                opt.hierarchical.reference1DjuryH = false;
            if (Hash.Checked)
            {
                opt.hash.hashCluster = true;
                opt.hash.selectionMethod = COL_SELECTION.ENTROPY;
            }
            results.Show();
            results.Focus();
            results.BringToFront();
            set.Save();
            var joined = data.ApplyFilters(data.filters);

            results.Run(processName+"_"+counter++,joined, opt);
        }

        private void button4_Click(object sender, EventArgs e)
        {
            previous = true;
            parent.Show();
            this.Close();

        }

        private void uQlustTree_FormClosed(object sender, FormClosedEventArgs e)
        {
            if (!previous)
                parent.Close();     
        }


        private void distanceControl1_Load(object sender, EventArgs e)
        {

        }

        private void radioButton1_CheckedChanged(object sender, EventArgs e)
        {
            if(radioButton1.Checked)
            {
                label4.Visible = true;
                numericUpDown1.Visible = true;
            }
            else
            {
                label4.Visible = false;
                numericUpDown1.Visible = false;
            }
        }

        private void radioHTree_CheckedChanged(object sender, EventArgs e)
        {
           // distanceControl2.Visible = !radioHTree.Checked;

        }
    }
}
