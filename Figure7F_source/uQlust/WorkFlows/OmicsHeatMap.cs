using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using Graph;
using phiClustCore;
using phiClustCore.Profiles;

namespace WorkFlows
{
    public partial class  OmicsHeatMap:Form 
    {
        Form parent;
        bool previous = false;
        public Options opt = new Options();
        public string processName = null;
        static int counter=0;
        ResultWindow results;
        string dataFileName = "";
        ProfileTree tree = new ProfileTree();
        OmicsDataSet data;
        public OmicsHeatMap(Form parent,OmicsDataSet data, ResultWindow results)
        {
            this.data = data;
            this.parent = parent;
            this.results = results;
            InitializeComponent();
            SetProfileOptions();
            label1.Text = "Num. of genes: " + data.geneLabels.Count;
        }


        void SetProfileOptions()
        {

            opt.clusterAlgorithm.Clear();
            opt.clusterAlgorithm.Add(ClusterAlgorithm.OmicsHeatMap);
            relevantC.Value = opt.hash.relClusters;
            
            numericUpDown1.Value = opt.hash.refPoints;
        }


        private void button1_Click_1(object sender, EventArgs e)
        {

        }

        private void button1_Click(object sender, EventArgs e)
        {
            OmicsDataSet locData=null;
            if (data.selectedGenes == null)
                locData = data;
            else
                locData = data.SelectColumns(data,data.selectedGenes.ToArray());
            Settings set = new Settings();
            set.Load();
            set.mode = INPUTMODE.OMICS;
            opt.dataDir.Clear();
            opt.profileFiles.Clear();
            opt.hash.relClusters = (int)relevantC.Value;
            opt.hash.reqClusters = (int)numericUpDown4.Value;
            opt.hash.useConsensusStates = consensus.Checked;
            opt.hash.perData = 90;

            opt.hash.combine = radioButton1.Checked;
            opt.hash.hashCluster = !radioButton1.Checked;
            if (fastHamming.Checked)
                opt.hierarchical.distance = DistanceMeasures.FASTHAMMING;
            if(radioEucl.Checked)
                opt.hierarchical.distance = DistanceMeasures.EUCLIDIAN;
            if(radioCosine.Checked)
                opt.hierarchical.distance = DistanceMeasures.COSINE;
            if (radioPearson.Checked)
                opt.hierarchical.distance = DistanceMeasures.PEARSON;
            if (radio1DJury.Checked)
                opt.hierarchical.distance = DistanceMeasures.HAMMING;
            
            if(radioButton1.Checked)
                opt.hierarchical.microCluster = ClusterAlgorithm.HashCluster;
            if (radioButton2.Checked)
                opt.hierarchical.microCluster = ClusterAlgorithm.HKmeans;
            if(radioButton3.Checked)
                opt.hierarchical.microCluster = ClusterAlgorithm.FastJuryBased;
            /*if(radioButton1.Checked)
                opt.hierarchical.microCluster = ClusterAlgorithm.;*/            
            opt.hierarchical.dummyProfileName = textBox1.Text;
            opt.hierarchical.reference1DjuryH = true;
            opt.clusterAlgorithm.Clear();
            opt.clusterAlgorithm.Add(ClusterAlgorithm.OmicsHeatMap);
            opt.hash.refPoints = (int)numericUpDown1.Value;
            results.Show();
            results.Focus();
            results.BringToFront();
            set.Save();
            counter++;
            GetProcessName prName = new GetProcessName("Set process name");
            DialogResult res=prName.ShowDialog();
            if (res == DialogResult.OK)
                processName = prName.name;

            opt.omics.processName = processName + "-" + counter + ".genprof";
            opt.omics.heatmap = false;            
            results.Run(opt.omics.processName,locData, opt);
        }

        private void button4_Click(object sender, EventArgs e)
        {
            previous = true;
            parent.Show();
            this.Hide();
            //this.Close();
        }

        private void OmicsHeatMap_FormClosed(object sender, FormClosedEventArgs e)
        {
            if (parent != null && !previous)
                parent.Close();
        }



        private void radioButton2_CheckedChanged(object sender, EventArgs e)
        {

        }

        private void button2_Click(object sender, EventArgs e)
        {
            DialogResult res=openFileDialog1.ShowDialog();
            if(res==DialogResult.OK)
            {
                textBox1.Text = openFileDialog1.FileName;
                    
            }

        }

        private void button3_Click(object sender, EventArgs e)
        {
            SuperGenesForm super = new SuperGenesForm(data.geneLabels,null,true);
            DialogResult res=super.ShowDialog();
            if(res==DialogResult.OK)
            {
                data.selectedGenes = super.GetSelectedGenes();
                label1.Text = "Num. of genes: " + data.selectedGenes.Count;
            }
        }
    }
}
