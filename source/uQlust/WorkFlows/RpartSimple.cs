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
    public partial class RpartSimple : Form, IclusterType
    {
        protected Options opt=new Options();
        protected ResultWindow results;
        public string processName;
        public static int counter = 0;
        public bool previous = false;
        CommonDialog dialog;
        protected Settings set;
        public Form parent;
        ProfileTree tree = new ProfileTree();
        public string dataFileName;
        OmicsDataSet data;
        public RpartSimple()
        {
            InitializeComponent();
        }
        public INPUTMODE GetInputType()
        {
            return set.mode;
        }
        public override string ToString()
        {
            return "Rpart";
        }
        public void HideRmsdLike()
        {

        }
        public RpartSimple(Form parent,OmicsDataSet data, Settings set, ResultWindow results, string fileName = null,string dataFileName=null)
        {
            InitializeComponent();
            this.data = data;
            this.parent = parent;
            dialog = openFileDialog1;
            this.Location = parent.Location;
            this.set = set;
            this.dataFileName = dataFileName;
            if (fileName != null)
            {
                opt.ReadOptionFile(fileName);                
                SetProfileOptions();
            }

            this.results = results;
        }
        public void ShowLabels()
        {
            label4.Visible = true;
            label5.Visible = true;
        }
        public virtual void SetProfileOptions()
        {
            relevantC.Value = opt.hash.relClusters;
            percentData.Value = opt.hash.perData;
            refPoints.Value = opt.hash.refPoints;
            opt.clusterAlgorithm.Clear();
            opt.clusterAlgorithm.Add(ClusterAlgorithm.HashCluster);
        }
        public void SetProfileName(string name)
        {
        }
        
        public virtual void GetData()
        {
            opt.dataDir.Clear();
            opt.profileFiles.Clear();
        
            opt.hash.relClusters = (int)relevantC.Value;
            opt.hash.perData = (int)percentData.Value;
            opt.hash.refPoints = (int)refPoints.Value;
            opt.hash.selReference = selectReference.Checked;
            set.Save();            
        }
        public virtual void button2_Click(object sender, EventArgs e)
        {
            GetData();            
            results.Show();
            results.Focus();
            results.BringToFront();
            results.Run(processName + "_" + counter++,data, opt);

        }

        void button4_Click(object sender, EventArgs e)
        {
            previous = true;
            parent.Show();
            this.Close();
        }        
        private void RpartSimple_FormClosed(object sender, FormClosedEventArgs e)
        {
                if(!previous)
                    parent.Close();           
        }
      
    }
}
