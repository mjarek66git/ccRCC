using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

using OxyPlot;
using OxyPlot.Annotations;
using OxyPlot.Series;
using phiClustCore;

namespace Graph
{
    public partial class PlotForm : Form
    {
        static Dictionary<string, int> refGenes = new Dictionary<string, int>();
        OxyColor[] col = new OxyColor[10];
        OmicsDataSet omData;
        PlotModel mModel=new PlotModel();
//        double[,] data;
        List<Tuple<int, double>>[] points=null;
        double[] val;
        double[] refData;
        
        bool binomial = false;
        List<OmicsDataSet> data;
        Dictionary<string, int>[] freqSamples = new Dictionary<string, int>[2];

        
        public PlotForm(List<OmicsDataSet> data)
        {
            this.data = data;
            InitializeComponent();
            col[0] = OxyColors.Olive;
            col[1] = OxyColors.Red;
            col[2] = OxyColors.Yellow;
            col[3] = OxyColors.YellowGreen;
            col[4] = OxyColors.SkyBlue;
            col[5] = OxyColors.Azure;
            col[6] = OxyColors.Black;
            col[7] = OxyColors.BlueViolet;
            col[8] = OxyColors.Brown;
            col[9] = OxyColors.Coral;

            freqSamples[0] = new Dictionary<string, int>();
            freqSamples[1] = new Dictionary<string, int>();

            foreach (var item in data)            
                comboBox2.Items.Add(item.Name);
            comboBox2.SelectedIndex = 0;

            if (comboBox1.Items.Count==0)
            {
                foreach (var item in omData.geneLabels)
                    comboBox1.Items.Add(item);
                comboBox1.SelectedIndex = 0;

            }
        }
        void PlotHistogRank()
        {
            double[] hist = new double[omData.data.columns];
            for (int i = 0; i < omData.data.columns; i++)
            {
                double sum = 0;
                for (int j = 0; j < omData.data.rows; j++)
                {
                    sum += omData.data[j, i];
                }
                hist[i] = omData.data.rows;
            }
            int[] index = new int[hist.Length];

            Array.Sort(index, hist);
            Array.Reverse(hist);
            Array.Reverse(index);
            

            LineSeries l = new LineSeries();
            for(int i=0;i<hist.Length;i++)
                l.Points.Add(new DataPoint(i, hist[index[i]]));


            mModel.Series.Add(l);

            this.plotView1.Model = mModel;
            this.plotView1.Model.InvalidatePlot(true);
        }
        void PlotHistog()
        {
            mModel.Series.Clear();
            mModel.Annotations.Clear();
            mModel.ResetAllAxes();
            double[] dataSet;
            int[] index;
            if (!checkBox1.Checked)
            {
                dataSet = new double[omData.data.rows];
                index = new int[dataSet.Length];
            for (int i = 0; i < dataSet.Length; i++)
            {
                dataSet[i] = omData.data[i, comboBox1.SelectedIndex];
                index[i] = i;
            }
            }
            else
            {
                dataSet = new double[omData.data.columns];
                index = new int[dataSet.Length];
                for (int i = 0; i < dataSet.Length; i++)
                {                    
                    dataSet[i] = omData.data[comboBox1.SelectedIndex,i];
                    index[i] = i;
                }
            }


            Array.Sort(dataSet,index);
            Dictionary<double, int> freq = new Dictionary<double, int>() ;
            Dictionary<double, double> val = new Dictionary<double, double>();            
            if (omData.codes == null || omData.codes.Count== 0)
            {
                double num = (dataSet[dataSet.Length - 1] - dataSet[0])/((double)numericUpDown1.Value);
                if (num == 0)
                {
                    if (plotView1.Model != null)
                    {
                        plotView1.Model.Series.Clear();
                        this.plotView1.Model.InvalidatePlot(true);
                    }
                    return;
                }
                freq = new Dictionary<double, int>();
                double tmp = 0;
                int k = 0;
                for (int i = 0; i < dataSet.Length; i++)
                {
                    if (dataSet[i] < (dataSet[0] + num * (k + 1)))
                    {
                        if (!freq.ContainsKey(k))
                            freq.Add(k, 0);
                        freq[k]++;
                        tmp += dataSet[i];
                    }
                    else
                    {
                        if (freq.ContainsKey(k))
                        {
                            tmp /= freq[k];
                            val.Add(k, tmp);
                        }
                        tmp = 0;
                        k++;
                        i--;
                    }
                }
                if (freq.Count > val.Count)
                {
                    tmp /= freq[k];
                    val.Add(k, tmp);
                }
            }
            else
            {
                for (int i = 0, k = 0; i < dataSet.Length; i++)
                {
                    if (!freq.ContainsKey(dataSet[i]))
                        freq.Add(dataSet[i], 0);
                    freq[dataSet[i]]++;
                }
                foreach (var item in freq)
                    val.Add(item.Key, item.Key);
            }    
            
            LineSeries l = new LineSeries();
            foreach(var item in freq)
                if(item.Key>0)
                    l.Points.Add(new DataPoint(val[item.Key], item.Value));
            
            mModel.Series.Add(l);

            if(binomial)
            {
                

                OmicsDataSet oData= data[comboBox2.Items.Count-1];
                double maxV = double.MinValue;
                for (int i = 0; i < dataSet.Length; i++)
                    if (maxV < dataSet[i])
                        maxV = dataSet[i];

                for (int i = 0; i < dataSet.Length; i++)                    
                        dataSet[i]/=maxV;

                var res=Descritize.NegativeBinomialFit(dataSet, true, maxV);
                for (int i = 0; i < dataSet.Length; i++)
                    dataSet[i] *= maxV;

                int countM, countP, countZ;
                countM = countP = countZ = 0;
                for (int i = 0; i < dataSet.Length; i++)
                {
                    if (dataSet[i] < oData.intervals[comboBox1.SelectedIndex][1].min)
                    {
                        countM++;
                        string key = omData.sampleLabels[index[i]];
                        if (!freqSamples[0].ContainsKey(key))                        
                            freqSamples[0].Add(key, 0);
                        freqSamples[0][key]++;
                    }
                    if (dataSet[i] >= oData.intervals[comboBox1.SelectedIndex][1].min && dataSet[i] <= oData.intervals[comboBox1.SelectedIndex][1].max)
                        countZ++;
                    if (dataSet[i] > oData.intervals[comboBox1.SelectedIndex][1].max)
                    {
                        countP++;
                        string key = omData.sampleLabels[index[i]];
                        if (!freqSamples[1].ContainsKey(key))
                            freqSamples[1].Add(key, 0);
                        freqSamples[1][key]++;

                    }
                }
                var ordered0 = freqSamples[0].OrderBy(e => e.Value).Select(e => new { frequency = e.Value, name = e.Key }).ToList();

                var ordered1 = freqSamples[1].OrderBy(e => e.Value).Select(e => new { frequency = e.Value, name = e.Key }).ToList();

                int count0=0;
                for (int i = 0; i < ordered0.Count; i++)
                    if (ordered0[i].frequency > 40)
                        count0++;

                int count1 = 0;
                for (int i = 0; i < ordered1.Count; i++)
                    if (ordered1[i].frequency > 40)
                        count1++;


                LineSeries l2 = new LineSeries();
                double maxS = double.MinValue;
                for (int i = 0; i < res.Item3.GetLength(0); i++)
                {
                    if (res.Item4[i] > 0)
                        l2.Points.Add(new DataPoint(res.Item4[i], res.Item3[i]));
                    if (maxS < res.Item3[i])
                        maxS = res.Item3[i];
                }
                var textAnnotationZ = new TextAnnotation
                {
                    Text = "CounterZ="+countZ.ToString(),
                    TextPosition = new DataPoint(dataSet[dataSet.Length/2], maxS)
                };
                mModel.Annotations.Add(textAnnotationZ);
                var textAnnotationP = new TextAnnotation
                {
                    Text = "CounterP=" + countP.ToString(),
                    TextPosition = new DataPoint(dataSet[dataSet.Length / 2], maxS - 600)
                };
                mModel.Annotations.Add(textAnnotationP);
                var textAnnotationM = new TextAnnotation
                {
                    Text ="CounterM="+countM.ToString(),
                    TextPosition = new DataPoint(dataSet[dataSet.Length / 2], maxS-1200)
                };
                mModel.Annotations.Add(textAnnotationM);

                LineSeries l3 = new LineSeries();
                l3.Points.Add(new DataPoint(res.Item1, 0));
                l3.Points.Add(new DataPoint(res.Item1, maxS));

                    LineSeries l4 = new LineSeries { Color = OxyColor.FromRgb(0, 0, 1) };
                    l4.Points.Add(new DataPoint(oData.intervals[comboBox1.SelectedIndex][1].min, 0));
                    l4.Points.Add(new DataPoint(oData.intervals[comboBox1.SelectedIndex][1].min, maxS));
                    mModel.Series.Add(l4);
                    l4 = new LineSeries { Color = OxyColor.FromRgb(0, 0, 1) };
                    l4.Points.Add(new DataPoint(oData.intervals[comboBox1.SelectedIndex][1].max, 0));
                    l4.Points.Add(new DataPoint(oData.intervals[comboBox1.SelectedIndex][1].max, maxS));
                    mModel.Series.Add(l4);

                


                mModel.Series.Add(l2);
                mModel.Series.Add(l3);
            }

            this.plotView1.Model = mModel;
            this.plotView1.Model.InvalidatePlot(true);
            //this.Refresh();
        }
        void ResconBinomial()
        {

        }

        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            PlotHistog();
        }

        private void comboBox2_SelectedIndexChanged(object sender, EventArgs e)
        {
            omData=data[comboBox2.SelectedIndex];
            if (((string)comboBox2.SelectedItem).Contains("super") && !((string)comboBox2.SelectedItem).Contains("Desc"))
            {
                binomial = true;
            }
            else
                binomial = false;
            {
                int remSelected = comboBox1.SelectedIndex;
                comboBox1.Items.Clear();
                if (!checkBox1.Checked)
                {
                foreach (var item in omData.geneLabels)
                    comboBox1.Items.Add(item);
                }
                else
                    foreach (var item in omData.sampleLabels)
                        comboBox1.Items.Add(item);

                if (remSelected != -1 && remSelected<comboBox1.Items.Count)
                    comboBox1.SelectedIndex = remSelected;
                else
                    comboBox1.SelectedIndex = 0;

                PlotHistog();
            }

        }

        private void numericUpDown1_ValueChanged(object sender, EventArgs e)
        {
            PlotHistog();
        }

        private void checkBox1_CheckedChanged(object sender, EventArgs e)
        {
            comboBox1.Items.Clear();
            if(checkBox1.Checked)
            {
                foreach (var item in omData.sampleLabels)
                    comboBox1.Items.Add(item);
                
            }
            else
                foreach (var item in omData.geneLabels)
                    comboBox1.Items.Add(item);

            comboBox1.SelectedIndex = 0;
            PlotHistog();
        }
    }
}
