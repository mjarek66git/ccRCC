using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.IO;
using System.Diagnostics;
using System.Data;
using System.Drawing;
//using System.Data;
using phiClustCore;
using phiClustCore.Interface;
using phiClustCore.Distance;
using phiClustCore.Profiles;

namespace phiClustCore
{
    public delegate void UpdateJob(string item,bool errorFlag=false,bool finishAll=true);
    public delegate void StartJob(string name,string al,string dirName,string measure);
    public delegate void ErrorMessage(string item);
    public delegate void UpdateMessage(string m);
    public delegate void ErrorJob();

    class ThreadParam
    {
        public string name;
        public int num;
        public int start;
        public int stop;
        public OmicsDataSet data;
        //public string dirName;
        //public DistanceMeasure distance;
    }
    public interface MessageUpdate
    {
        void UpdateMessage(string message);
        void ActivateUpdateting();
        void CloseUpdateting();
    }
    class OmicsData
    {
        public List<Color> colorMap;
        public List<int> index;
        public List<string> labels;
    };
    public class JobManager
    {
        Dictionary<string, Thread> runnigThreads = new Dictionary<string, Thread>();
        public Dictionary<string, ClusterOutput> clOutput = new Dictionary<string, ClusterOutput>();
        Dictionary<string, IProgressBar> progressDic = new Dictionary<string, IProgressBar>();
        public Options opt = new Options();
        string clType = "";
        string currentProcessName="";
        Thread startProg;
       

        public UpdateJob updateJob;
        public StartJob beginJob;
        public ErrorMessage message;
        public UpdateMessage upMessage;
        public ErrorJob errorJob;
        private Object thisLock = new Object();


        public Dictionary<string,double> ProgressUpdate()
        {
            Dictionary<string, double> res = new Dictionary<string, double>();          

            lock (thisLock)
            {
                if (progressDic.Count == 0)
                    return null;

                foreach (var item in progressDic)
                {
                    res.Add(item.Key, item.Value.ProgressUpdate());
                }
                foreach (var item in res)
                    if (item.Value == 1)
                        progressDic.Remove(item.Key);

                return res;
            }            
            
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }

        Dictionary<string,List<double>> MakeDictionary(OmicsDataSet data)
        {
            Dictionary<string, List<double>> res = new Dictionary<string, List<double>>();
            List<string> lab;
            if (data.data.rows == data.geneLabels.Count)
                lab = data.geneLabels;
            else
                lab = data.sampleLabels;

            for(int i=0;i<data.data.rows;i++)
            {
                res.Add(lab[i], new List<double>());
                for (int j = 0; j < data.data.columns; j++)
                    res[lab[i]].Add(data.data[i, j]);

            }
            return res;
        }
        private void RunHashCluster(OmicsDataSet data,string name)
        {
            DateTime cpuPart1 = DateTime.Now;
            HashCluster hk = null;

            Dictionary<string, List<double>> dat = MakeDictionary(data);

            hk = new HashCluster(dat,opt);
            
            progressDic.Add(name, hk);
            if (beginJob != null)
                beginJob(currentProcessName, hk.ToString(), name, "HAMMING");
            hk.InitHashCluster();


            DateTime cpuPart2 = DateTime.Now;



            ClusterOutput output;
            output = hk.RunHashCluster();
            output.data = data;
            UpdateOutput(name, output, "HAMMING", cpuPart1, cpuPart2, hk);

        }
        private void RunGuidedHashCluster(OmicsDataSet dataS,string name)
        {
            DateTime cpuPart1 = DateTime.Now;
            GuidedHashCluster hk = null;

            var data = MakeDictionary(dataS);

            hk = new GuidedHashCluster(data, opt);


            progressDic.Add(name, hk);
            if (beginJob != null)
                beginJob(currentProcessName, hk.ToString(), name, "HAMMING");
            hk.InitHashCluster();

            DateTime cpuPart2 = DateTime.Now;

            ClusterOutput output;
            //hk.ReadClassLabels(opt.hNNLabels);
            output = hk.RunHashCluster();
            output.data = dataS;
            UpdateOutput(name, output, "HAMMING", cpuPart1, cpuPart2, hk);

        }
        public void RunJuryBasedCluster(OmicsDataSet dataS, string name)
        {
            DateTime cpuPart1 = DateTime.Now;
            List<int> thresholds = new List<int>();
            List<int> thresholdsMin = new List<int>();
            thresholds.Add(200);
            thresholds.Add(300);
            thresholds.Add(500);
            thresholds.Add(800);
            thresholds.Add(1000);
            thresholdsMin.Add(4000);
            thresholdsMin.Add(5000);
            thresholdsMin.Add(6000);
            thresholdsMin.Add(7000);
            thresholdsMin.Add(9000);

            ClusterOutput cl = null;
            FastRPart rp=null;
            foreach (var item in thresholdsMin)
            {
                foreach(var it in dataS.filters)
                {
                    if(it is MinMaxAmountThreshold)
                    {
                        MinMaxAmountThreshold m = (MinMaxAmountThreshold)it;
                        m.thresholdMax = 31000;
                        m.thresholdMin = item;

                    }
                }
                var joined = dataS.ApplyFilters(dataS.filters);

                var data = MakeDictionary(joined);
                continue;
                rp = new FastRPart(data);
                var clust = rp.FastCombineKeysNew(100);

                List<List<string>> xx = new List<List<string>>();

                foreach(var it in clust)
                {
                    if (it.Count > 10)
                        xx.Add(it);
                }


                cl = new ClusterOutput();

                cl.data = dataS;
                cl.clusters = new clusterRes();
                cl.clusters.list = xx;
             
            }
            DateTime cpuPart2 = DateTime.Now;
            UpdateOutput(name, cl, "Fast", cpuPart1, cpuPart2,rp );
        }

        private void RunHashDendrogCombine(OmicsDataSet dataS,string name)
        {
            DateTime cpuPart1 = DateTime.Now;
            
            HashClusterDendrog hk = null;

            var data = MakeDictionary(dataS);

            hk = new HashClusterDendrog(data, opt);
            

            ClusterOutput output;
            if (beginJob != null)
                beginJob(currentProcessName, hk.ToString(),name,"NONE");
            progressDic.Add(name, hk);
            hk.InitHashCluster();
            
            DateTime cpuPart2 = DateTime.Now;
            output = hk.RunHashDendrogCombine();
            output.data = dataS;
            if (opt.hash.profileName.Contains("omics"))
            {
                /*string intvName = "generatedProfiles/OmicsIntervals_" + opt.omics.processName + ".dat";
                if (File.Exists(intvName))
                {
                    OmicsData res = ReadOmicsData(intvName);
                    output.profilesColor = res.colorMap;
                    output.aux2 = res.labels;
                    
                }*/
            }
            UpdateOutput(name, output, hk.UsedMeasure(), cpuPart1, cpuPart2, hk);

        }
        OmicsData ReadOmicsData(string fileName)
        {
            OmicsData res=new OmicsData();
            List<Color> colorMap;
            List<int> index = new List<int>();
            List<string> labels=new List<string>();
            StreamReader wr = new StreamReader(fileName);
            Dictionary<int, double[]> codingInterv = null;
            codingInterv = new Dictionary<int, double[]>();
            string line = wr.ReadLine();
            while (line != null)
            {
                string[] aux = line.Split(' ');
                if (aux[0].Contains("Label"))
                    labels.Add(aux[1]);
                if (aux.Length == 4 && aux[0].Contains("Code"))
                {
                    double[] tab = new double[2];
                    tab[0] = Convert.ToDouble(aux[2]);
                    tab[1] = Convert.ToDouble(aux[3]);
                    codingInterv.Add(Convert.ToInt32(aux[1]), tab);
                }
                line = wr.ReadLine();
            }
            wr.Close();

            colorMap = new List<Color>(codingInterv.Count);
            for (int i = 0; i < codingInterv.Count; i++)
            {
                colorMap.Add(Color.Black);
                index.Add(0);
            }
            List<int> hot = new List<int>();
            List<int> cool = new List<int>();
            foreach (var item in codingInterv.Keys)
            {
                if (codingInterv[item][0] < 0 && codingInterv[item][1] < 0)
                    cool.Add(item);
                else
                    if (codingInterv[item][0] >= 0 && codingInterv[item][1] >= 0)
                        hot.Add(item);
            }

            hot.Sort((x, y) => x.CompareTo(y));
            cool.Sort((x, y) => y.CompareTo(x));
            if (hot.Count > 1)
            {

                Color stepper = Color.FromArgb(
                                       (byte)((255 - 0) / (hot.Count)),
                                       (byte)((255 - 0) / (hot.Count)),
                                       (byte)((0 - 0) / (hot.Count)));
                for (int i = 0; i < hot.Count; i++)
                {
                    colorMap[hot[i]] = Color.FromArgb(
                                                0 + (stepper.R * (i + 1)),
                                                0 + (stepper.G * (i + 1)),
                                                0 + (stepper.B * (i + 1)));
                    index[hot[i]] = i + 1;
                }
            }
            else
                if (hot.Count == 1)
                {
                    colorMap[hot[0]] = Color.FromArgb(255, 255, 0);
                    index[hot[0]] = 1;
                }

            if (cool.Count > 1)
            {
                Color stepper = Color.FromArgb(
                                       (byte)((20 - 0) / (cool.Count)),
                                       (byte)((0 - 0) / (cool.Count)),
                                       (byte)((255 - 0) / (cool.Count)));
                for (int i = 0; i < cool.Count; i++)
                {
                    colorMap[cool[i]] = Color.FromArgb(
                                                0 + (stepper.R * (i + 1)),
                                                0 + (stepper.G * (i + 1)),
                                                0 + (stepper.B * (i + 1)));
                    index[cool[i]] = -(i+1);
                }
            }
            else
                if (cool.Count == 1)
                {
                    colorMap[cool[0]] = Color.FromArgb(20, 0, 255);
                    index[cool[0]] = -1;
                }



            res.colorMap = colorMap;
            res.labels = labels;
            res.index = index;
            return res;
        }

        Dictionary<string,string> ReadLabels(string fileName)
        {
            StreamReader lFile = new StreamReader(fileName);
            Dictionary<string, string> classLabels;
            classLabels = new Dictionary<string, string>();
            string line = lFile.ReadLine();
            while (line != null)
            {
                string[] aux = line.Split(' ');
                if (aux.Length == 2)
                {
                    classLabels.Add(aux[0], aux[1]);
                }
                line = lFile.ReadLine();
            }

            lFile.Close();
            return classLabels;
        }
        void RunHTree(OmicsDataSet dataS,string name)
        {
            DateTime cpuPart1 = DateTime.Now;
            HashCluster hCluster;

            var data = MakeDictionary(dataS);

            hCluster = new HashCluster(data,  opt);

            HTree h = new HTree(data, hCluster);
            beginJob(currentProcessName, h.ToString(), name, "HAMMING");
            progressDic.Add(name, h);
            hCluster.InitHashCluster();
   
            DateTime cpuPart2 = DateTime.Now;            
          

            ClusterOutput output = new ClusterOutput();
            output = h.RunHTree();
            output.data = dataS;
            UpdateOutput(name, output, "NONE", cpuPart1, cpuPart2, h);


        }
        ClusterOutput MakeDummyClusters(Dictionary<string,List<double>> dic)
        {
            ClusterOutput output = new ClusterOutput();
            
            output.clusters = new clusterRes();
            output.clusters.list = new List<List<string>>();

            foreach (var item in dic)
            {
                List<string> ll = new List<string>();
                ll.Add(item.Key);
                output.clusters.list.Add(ll);
            }

            return output;
        }

        Tuple<ClusterOutput,DistanceMeasure> MakeClusters(Dictionary<string,List<double>> data,string name,int num)
        {
            ClusterOutput output = null;
            if (num <= data.Count)
            {
                switch (opt.hierarchical.microCluster)
                {
                    case ClusterAlgorithm.HashCluster:
                        HashCluster hashUpper = null;                        
                        hashUpper = new HashCluster(data, opt);
                        hashUpper.InitHashCluster();
                        hashUpper.StartProgress = 0;
                        hashUpper.EndProgress = 0.25;
                        progressDic.Add(name, hashUpper);

                        output = hashUpper.RunHashCluster();
                        break;
                    case ClusterAlgorithm.FastJuryBased:
                     //   var dataOut = MakeDictionary(data);                        
                        FastRPart rp = new FastRPart(data);
                        var clust = rp.FastCombineKeysNew(num);
                        output=new ClusterOutput();
                        output.clusters = new clusterRes();
                        output.clusters.list = clust;


                        break;
                    case ClusterAlgorithm.HKmeans:
                        DistanceMeasure dist = null;
                        switch (opt.hierarchical.distance)
                        {
                            case DistanceMeasures.HAMMING:
                                dist = new JuryDistance(data, true);
                                break;


                            case DistanceMeasures.COSINE:
                                dist = new CosineDistance(data,true);
                                break;
                            case DistanceMeasures.EUCLIDIAN:
                                dist = new Euclidian(data,true);
                                break;

                            case DistanceMeasures.PEARSON:
                                dist = new Pearson(data, true);
                                break;
                            case DistanceMeasures.FASTHAMMING:
                                dist = new Fast3States(1, data);
                                //dist = new FastDiscreteDist(1,data);
                                break;
                        }

                        dist.InitMeasure();
                        kMeans km = new kMeans(dist);
                        List<string> lData = new List<string>(data.Keys);
                        output = km.kMeansLevel(opt.hash.relClusters, 500, lData);
                        return new Tuple<ClusterOutput, DistanceMeasure>(output, dist);                        
                }


            }
            else
                output = MakeDummyClusters(data);

            return new Tuple<ClusterOutput, DistanceMeasure>(output,null);


        }
        ClusterOutput LoadClusters(string fileName)
        {
            StreamReader f = new StreamReader(fileName);
            Dictionary<int, List<string>> clDic = new Dictionary<int, List<string>>();
            string line = f.ReadLine() ;
            line = f.ReadLine();
            while (line!=null)
            {
                string[] tmp = line.Split(new char[] { ' ', '\t' });
                if(tmp.Length==2)
                {
                    int num=0;
                    int.TryParse(tmp[1],out num);
                    if(!clDic.ContainsKey(num))                    
                        clDic.Add(num, new List<string>());
                    clDic[num].Add(tmp[0]);
                }
                line = f.ReadLine();
            }
            ClusterOutput res = new ClusterOutput();
            res.clusters = new clusterRes();
            res.clusters.list = new List<List<string>>();
            foreach (var item in clDic.Keys)
                res.clusters.list.Add(clDic[item]);

            f.Close();

            return res;
        }
        private void RunOmicsHeatmap(OmicsDataSet dataS,string name)
        {
            DateTime cpuPart1 = DateTime.Now;
            HashClusterDendrog upper = null;
            HashClusterDendrog left = null;
            
            
            ClusterOutput aux,output;
            var joined = dataS;
            if (dataS.prev!=null)
                joined = dataS.prev;//.ApplyFilters(dataS.filters);
 //           joined.Save("Disc");
            // joined = joined.SelectTotestData();
            var data = MakeDictionary(joined);
            //  opt.hierarchical.distance = DistanceMeasures.PEARSON;
            int numLeftClusters, numUpperClusters;
            opt.omics.heatmap = false;
            numUpperClusters = opt.hash.reqClusters;
            numLeftClusters = opt.hash.relClusters;
            OmicsDataSet dataTrans = joined.Transpose();

            //dataTrans.Save("Disc_rot");

            var dataT = MakeDictionary(dataTrans);


            //joined.Save("test.dat");
            if (beginJob != null)
                beginJob(currentProcessName, "HeatMap", name, "HAMMING");

           
            ClusterOutput outputUpper;
            DateTime cpuPart2 = DateTime.Now;
            opt.hash.relClusters = numUpperClusters;
            outputUpper =MakeClusters(dataT, name, numUpperClusters).Item1;
            
            progressDic.Remove(name);

            ClusterOutput outputLeft;
            opt.omics.heatmap = true;
//            profOm.Save(OmicsInput.fileName);
            opt.hash.relClusters = numLeftClusters;
            Tuple<ClusterOutput, DistanceMeasure> cc = null;
                cc = MakeClusters(data, name, numLeftClusters);
                outputLeft = cc.Item1;
                if (opt.hierarchical.dummyProfileName.Length > 0)
                {
                    var xx=StaticDic.GetClusters(data);
                    if (xx != null)
                        outputLeft.clusters.list = xx;
                }
                //left.MakeDummyDendrog

            Settings set = new Settings();
            set.Load();



            List<string> listUpper = joined.geneLabels;
            List<string> listLeft = joined.sampleLabels;            


            var refStructUpper = HashClusterDendrog.StructuresToDenrogram(dataT,outputUpper.clusters.list);
            var refStructLeft = HashClusterDendrog.StructuresToDenrogram(data,outputLeft.clusters.list);

            output = new ClusterOutput();
     
           left = new HashClusterDendrog(data,opt);
           upper = new HashClusterDendrog(dataT,opt);


            output.clusters = outputUpper.clusters;
            output.clusters2 = outputLeft.clusters;
            if(cc!=null)
                output.distM = cc.Item2;
           output.aux1 = listLeft;
           output.aux2 = listUpper;
           output.nodes = new List<HClusterNode>();
           progressDic.Remove(name);
           upper.StartProgress = 0.5;
           upper.EndProgress = 0.75;
           progressDic.Add(name, upper);
            if(dataS.prev!=null)
                output.data = dataS;
            else
                output.data=joined ;
          aux=upper.DendrogUsingMicroClusters(dataT,refStructUpper.Item2, refStructUpper.Item1);

           List<HClusterNode> upperLeafs = aux.hNode.GetLeaves();

            /*Dictionary<string, int> toCheckLeft = new Dictionary<string, int>();
            foreach(var item in hashUpper.selectedColumnsHash)
            {
                toCheckLeft.Add(listLeft[item], 0);
            }*/


            

           output.nodes.Add(aux.hNode);
           progressDic.Remove(name);
           left.StartProgress = 0.75;
           left.EndProgress = 1.0;
           progressDic.Add(name, left);

           aux = left.DendrogUsingMicroClusters(data, refStructLeft.Item2, refStructLeft.Item1);
            
            /*Dictionary<string, int> toCheckUpper = new Dictionary<string, int>();
            foreach (var item in hashLeft.selectedColumnsHash)
            {
                toCheckUpper.Add(listUpper[item], 0);
            }


            
            foreach (var item in upperLeafs)
            {
                foreach (var it in item.setStruct)
                    if (toCheckUpper.ContainsKey(it))
                    {
                        item.flagSign = true;
                        break;
                    }
            }

            List<HClusterNode> leftLeafs = aux.hNode.GetLeaves();
            foreach (var item in leftLeafs)
            {
                foreach (var it in item.setStruct)
                    if (toCheckLeft.ContainsKey(it))
                    {
                        item.flagSign = true;
                        break;
                    }
            }*/



            output.nodes.Add(aux.hNode);

           UpdateOutput(name, output, upper.UsedMeasure(), cpuPart1, cpuPart2, left);

        }
        private void RunOmicsHeatmap2(Dictionary<string,List<double>> data,string name)
        {
            DateTime cpuPart1 = DateTime.Now;
            HashClusterDendrog upper = null;
            HashClusterDendrog left = null;
            OmicsProfile profOm = new OmicsProfile();
            //  opt.hierarchical.distance = DistanceMeasures.PEARSON;
            profOm.Load(OmicsInput.fileName);

            upper = new HashClusterDendrog(data, opt);

            //profOm.LoadOmicsSettings();

            string nameUpper = name + "upper";
            if (beginJob != null)
                beginJob(currentProcessName, upper.ToString(), name, "HAMMING");

            upper.EndProgress = 0.5;
            progressDic.Add(name, upper);
            upper.InitHashCluster();
            DateTime cpuPart2 = DateTime.Now;

            ClusterOutput output, aux;
            output = new ClusterOutput();
            aux = upper.RunHashDendrogCombine();

            output.nodes = new List<HClusterNode>();
            output.nodes.Add(aux.hNode);
            //  UpdateOutput(nameUpper, dirName, alignmentFile, output, "NONE", cpuPart1, cpuPart2, upper);
            opt.hash.relClusters = opt.hash.reqClusters;
            left = new HashClusterDendrog(data, opt);

            //profOm.transpose = true;
            profOm.oInput.heatmap = true;
            profOm.Save(OmicsInput.fileName);
            string nameLeft = name + "left";
            //progressDic.Add(nameLeft, left);
            left.StartProgress = 0.5;
            left.EndProgress = 1.0;
            progressDic[name] = left;
            left.InitHashCluster();

            aux = left.RunHashDendrogCombine();
            output.nodes.Add(aux.hNode);
            output.aux1 = left.dataKeys;
            UpdateOutput(name, output, upper.UsedMeasure(), cpuPart1, cpuPart2, left);

        }


        private void RunHierarchicalCluster(OmicsDataSet dataS,string name)
        {
            DateTime cpuPart1 = DateTime.Now;
            DistanceMeasure distance = null;
            //distance.CalcDistMatrix(distance.structNames);
            // opt.hierarchical.atoms = PDB.PDBMODE.ALL_ATOMS;
            var data = MakeDictionary(dataS);

            distance = CreateMeasure(data,opt.hierarchical.distance, opt.hierarchical.reference1DjuryAglom);

            DebugClass.WriteMessage("Measure Created");          
            hierarchicalCluster hk = new hierarchicalCluster(distance, opt);
            if (beginJob != null)
                beginJob(currentProcessName,hk.ToString(), name, distance.ToString());
            clType = hk.ToString();
            ClusterOutput output;
            progressDic.Add(name, hk);
            distance.InitMeasure();
            DateTime cpuPart2 = DateTime.Now;
            output = hk.HierarchicalClustering(new List<string>(distance.structNames.Keys));
            output.data = dataS;
            UpdateOutput(name,output, distance.ToString(), cpuPart1, cpuPart2, hk);

        }
        private void RunHKMeans(OmicsDataSet dataS,string name)
        {
            DateTime cpuPart1 = DateTime.Now;
            ClusterOutput clustOut = null;
            DistanceMeasure distance = null;

            var data = MakeDictionary(dataS);

            distance = CreateMeasure(data,opt.hierarchical.distance, opt.hierarchical.reference1DjuryKmeans);
         
            kMeans km;

            km = new kMeans(distance,true);
            if (beginJob != null)
                beginJob(currentProcessName, km.ToString(), name, distance.ToString());

            progressDic.Add(name, km);
            DateTime cpuPart2 = DateTime.Now;
            distance.InitMeasure();

            

            clType = km.ToString();
            km.BMIndex = opt.hierarchical.indexDB;
            km.threshold = opt.hierarchical.numberOfStruct;
            km.maxRepeat = opt.hierarchical.repeatTime;
            km.maxK = opt.hierarchical.maxK;
            clustOut = km.HierarchicalKMeans();
            clustOut.data = dataS;
            UpdateOutput(name,clustOut, distance.ToString(), cpuPart1, cpuPart2, km);
        }
       
     
        private void Run1DJury(OmicsDataSet dataS,string name)
        {
            DateTime cpuPart1 = DateTime.Now;
            ClusterOutput output;


            jury1D ju=new jury1D();
            if (beginJob != null)
                beginJob(currentProcessName, ju.ToString(), name, "NONE");

            progressDic.Add(name,ju);

            var data = MakeDictionary(dataS);
            //DistanceMeasure distance = CreateMeasure();
                if (opt.other.alignGenerate)
                    opt.other.alignFileName = "";
                ju.PrepareJury(data);

                
            clType = ju.ToString();
            DateTime cpuPart2 = DateTime.Now;
            //jury1D ju = new jury1D(opt.weightHE,opt.weightC,(JuryDistance) distance);
            //output = ju.JuryOpt(new List<string>(ju.stateAlign.Keys));
            if (ju.alignKeys != null)
            {
              
                output = ju.JuryOptWeights(ju.alignKeys);
            }
            else
            {
                UpadateJobInfo(name, true, false);
                throw new Exception("Alignment is epmty! Check errors");
            }
            output.data = dataS;
            UpdateOutput(name, output,ju.ToString(), cpuPart1,cpuPart2, ju);
        }
      
        private void UpdateOutput(string name, ClusterOutput output, string distStr, DateTime cpuPart1, DateTime cpuPart2, object obj)
        {           
            output.clusterType = obj.ToString();
            output.measure = distStr.ToString();
            DateTime cc = DateTime.Now;
            TimeSpan preprocess=new TimeSpan();
            TimeSpan cluster=new TimeSpan();
            if(cpuPart1!=null && cpuPart2!=null)
                preprocess = cpuPart2.Subtract(cpuPart1);
            if(cpuPart2!=null)
                cluster = cc.Subtract(cpuPart2);

            output.time = "Prep="+String.Format("{0:F2}", preprocess.TotalMinutes);
            
            if(cpuPart2!=null)
                output.time += " Clust=" + String.Format("{0:F2}", cluster.TotalMinutes);
            output.name = name;
            output.peekMemory = Process.GetCurrentProcess().PeakWorkingSet64;
            Process.GetCurrentProcess().Refresh();
            progressDic.Remove(name);
            //Process.GetCurrentProcess().
            clOutput.Add(output.name, output);
            UpadateJobInfo(name, false,false);
        }
        public void UpadateJobInfo(string processName, bool errorFlag,bool finishAll)
        {
            if (updateJob != null)
                updateJob(processName, errorFlag,finishAll);

        }
        public void FinishThread(string processName,bool errorFlag)
        {
            lock (runnigThreads)
            {
                UpadateJobInfo(currentProcessName, errorFlag,true);
                runnigThreads.Remove(processName);
                if (progressDic.ContainsKey(currentProcessName))
                    progressDic.Remove(currentProcessName);
            }
        }
        public void RemoveJob(string jobName)
        {
            if (runnigThreads.ContainsKey(jobName))
            {
                lock (runnigThreads)
                {
                    runnigThreads[jobName].Abort();
                    runnigThreads.Remove(jobName);
                }
            }
        }
        private DistanceMeasure CreateMeasure(Dictionary<string,List<double>> data,DistanceMeasures measure,bool jury1d)
        {
            DistanceMeasure dist=null;
            switch(measure)
            {
                case DistanceMeasures.HAMMING:
                        dist = new JuryDistance(data, jury1d);
                    break;
                case DistanceMeasures.COSINE:
                        dist = new CosineDistance(data, jury1d);
                    break;

                case DistanceMeasures.PEARSON:
                    dist = new Pearson(data, jury1d);
                    break;


            }
            return dist;
        }
        string MakeName(object processParams,ClusterAlgorithm alg,int counter)
        {
            string currentProcessName = "";
            if (((ThreadParam)processParams).name != null && ((ThreadParam)processParams).name.Length > 0)
                currentProcessName = ((ThreadParam)processParams).name + ";" + counter;
            else
                currentProcessName = alg.ToString() + ";" + counter;

            return currentProcessName;
        }
        public void StartAll(object processParams)
        {
            ErrorBase.ClearErrors();
            string orgProcessName = ((ThreadParam)processParams).name;
            currentProcessName = ((ThreadParam)processParams).name;
            OmicsDataSet data= ((ThreadParam)processParams).data;
            int counter = 1;
            try
            {
                    foreach (var alg in opt.clusterAlgorithm)
                    {                        
                            currentProcessName = MakeName(processParams, alg, counter);
                            //if (beginJob != null)
                              //  beginJob(currentProcessName, alg.ToString("g"), item, opt.GetDistanceMeasure(alg));

                            switch (alg)
                            {
                                case ClusterAlgorithm.FastJuryBased:
                                    RunJuryBasedCluster(data, currentProcessName);
                                    break;
                                case ClusterAlgorithm.uQlustTree:
                                    RunHashDendrogCombine(data,currentProcessName);
                                    break;                                
                                case ClusterAlgorithm.HashCluster:
                                    RunHashCluster(data,currentProcessName);
                                    break;
                                case ClusterAlgorithm.GuidedHashCluster:
                                    RunGuidedHashCluster(data,currentProcessName);
                                    break;
                                case ClusterAlgorithm.HierarchicalCluster:
                                    RunHierarchicalCluster(data,currentProcessName);
                                    break;
                                case ClusterAlgorithm.HKmeans:
                                    RunHKMeans(data,currentProcessName);
                                    break;
                                case ClusterAlgorithm.Jury1D:
                                    Run1DJury(data,currentProcessName);
                                    break;
                                case ClusterAlgorithm.HTree:
                                    RunHTree(data,currentProcessName);
                                    break;
                                 case ClusterAlgorithm.OmicsHeatMap:
                                    RunOmicsHeatmap(data, currentProcessName);
                                    break;
                            
                        }
                    }

                FinishThread(orgProcessName, false);
            
            }
            catch (Exception ex)
            {
                FinishThread(orgProcessName, true);
                message(ex.Message);
            }
            
        }
        public void RunJob(string processName,OmicsDataSet data)
        {
            ThreadParam tparam=new ThreadParam();

             startProg = new Thread(StartAll);
             tparam.name = processName;
            tparam.data = data;
             startProg.Start(tparam);            
            lock (runnigThreads)
            {
                runnigThreads.Add(processName, startProg);
            }
        }
        public void WaitAllNotFinished()
        {
            while (runnigThreads.Count > 0)
            {
                Thread.Sleep(1000);
            }
        }
        public void SaveOutput(string fileName)
        {            
            StreamWriter w = new StreamWriter(fileName);
            int count=0;
            foreach (var item in clOutput.Keys)
            {
                string name=fileName+count++;
                w.WriteLine(name);
                GeneralFunctionality.SaveBinary(name, clOutput[item]);
            }
            w.Close();

        }
        public void LoadOutput(string fileName)
        {
            ClusterOutput outP;

            if (File.Exists(fileName))
            {
                clOutput.Clear();
                string line;
                StreamReader r = new StreamReader(fileName);
                while (!r.EndOfStream)
                {
                    line = r.ReadLine();
                    if (File.Exists(line))
                    {
                        outP = ClusterOutput.Load(line);
                        clOutput.Add(outP.name, outP);
                    }
                }
                r.Close();
            }
        }

        public void LoadExternal(string fileName, Func<string, string, ClusterOutput> Load)
        {
            ClusterOutput aux;
            if (File.Exists(fileName))
            {
                string line;
                StreamReader r = new StreamReader(fileName);
                line = r.ReadLine();
                string dirNameOrg = line;
                while (!r.EndOfStream)
                {
                    line = r.ReadLine();
                    if (File.Exists(line))
                    {
                        string nn=Path.GetFileName(line);
                        if(nn.Contains("list"))
                        {
                            string[] a = nn.Split('.');
                            nn = a[0];
                        }
                        if (nn.Contains("_"))
                        {
                            string[] aa = nn.Split('_');
                            nn = aa[0];
                        }

                        string dirName = dirNameOrg+Path.DirectorySeparatorChar+nn;
                        aux = Load(line,dirName);
                        clOutput.Add(aux.name, aux);
                    }
                }
                r.Close();
            }
        }

        public void LoadExternalF(string fileName)
        {
            LoadExternal(fileName, ClusterOutput.LoadExternal);
        }
        public void LoadExternalPleiades(string fileName)
        {
            LoadExternal(fileName,ClusterOutput.LoadExternalPleiades);
        }
        public void LoadExternalPconsD(string fileName)
        {
            LoadExternal(fileName, ClusterOutput.LoadExternalPconsD);
        }

    }

}
