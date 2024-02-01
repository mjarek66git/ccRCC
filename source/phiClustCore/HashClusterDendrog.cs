using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using phiClustCore;
using System.Collections.Concurrent;
using System.Threading;
using System.Threading.Tasks;
using System.Diagnostics;
using phiClustCore.Interface;
using System.Data;
using phiClustCore.Distance;
using System.Text.RegularExpressions;

namespace phiClustCore
{
    public class HashClusterDendrog:HashCluster,IProgressBar
    {
         DistanceMeasures dMeasure;
         DistanceMeasure dist=null;
         AglomerativeType linkageType;        
         bool jury1d;
         string refJuryProfile;
         HierarchicalCInput hier;
         hierarchicalCluster hk = null;
         public HashClusterDendrog(Dictionary<string,List<double>> data, Options opt)
             : base(data, opt)
        {
            this.opt = opt;
            this.dMeasure = opt.hierarchical.distance;
            this.linkageType = opt.hierarchical.linkageType;
            this.jury1d = opt.hierarchical.reference1DjuryH;
            this.refJuryProfile = opt.hierarchical.jury1DProfileH;            
            hier = opt.hierarchical;
        }
        public void InitHashClusterDendrog()
         {
             base.InitHashCluster();
         }
        public new double ProgressUpdate()
        {
            double progress = 0;
            progress=0.5 * (double)currentV / maxV;

            if (hk != null)
                progress += 0.5 * hk.ProgressUpdate() ;

            return StartProgress + (EndProgress - StartProgress) * progress;
        }
        public new Exception GetException()
        {
            return null;
        }
        public new List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }

        public string UsedMeasure()
         {
             if (dist != null)
                 return dist.ToString();

             return "NONE";
         }
         public override string ToString()
         {        
             return "uQlust:Tree";
         }

         public ClusterOutput RunHashDendrogCombine()
         {
             ClusterOutput output = DendrogUsingMeasures(dataKeys);
             return output;
         }
        public void SaveDistances(string fileName)
        {
            if (dist != null)
                dist.SaveDistancedToFile(fileName);
        }
        static public List<string> ClustersReferences(Dictionary<string,List<double>> data,List<List<string>> structures,Options opt)
        {
            List<string> refCluster = new List<string>();
           
            jury1D juryLocal = new jury1D();
            juryLocal.PrepareJury(data);



            foreach (var item in structures)
            {
                ClusterOutput output = juryLocal.JuryOptWeights(item);
                if (item.Count > 2)
                    refCluster.Add(output.juryLike[0].Key);
                else
                    refCluster.Add(item[0]);
            }

            return refCluster;
        }
        static public Tuple<Dictionary<string,List<string>>, Dictionary<string, List<double>>> StructuresToDenrogram(Dictionary<string,List<double>> data,List<List<string>> structures)
         {
            Dictionary<string, List<string>> translateToCluster = new Dictionary<string, List<string>>();
            Dictionary<string, List<double>> dicStruct = new Dictionary<string, List<double>>();
             jury1D juryLocal = new jury1D();
             juryLocal.PrepareJury(data);

            

            foreach (var item in structures)
             {
                ClusterOutput output = null;
                try
                {
                    output = juryLocal.JuryOptWeights(item);
                }
                catch(Exception ex)
                {
//                    Console.WriteLine("UPS");
                }
                if (output == null)
                    continue;
                string rep = output.juryLike[0].Key;
                if (dicStruct.ContainsKey(rep))
                    if (output.juryLike.Count > 1)
                        rep = output.juryLike[1].Key;
                    else
                        continue;
                dicStruct.Add(rep, data[rep]);
                translateToCluster.Add(rep, new List<string>());               
                translateToCluster[rep].Add(rep);
                if (item.Count > 2)
                 {                     
                     foreach (var i in item)
                         if (!i.Equals(rep))
                            translateToCluster[rep].Add(i);                     
                 }
             }
            

             return new Tuple<Dictionary<string, List<string>>, Dictionary<string, List<double>>>(translateToCluster,dicStruct);
         }
        public static ClusterOutput MakeDummyDendrog(List<List<string>> clusters)
        {
            ClusterOutput res = new ClusterOutput();

            res.hNode = new HClusterNode();
            res.hNode.joined = new List<HClusterNode>();
            res.hNode.setStruct = new List<string>();
            res.hNode.refStructure = clusters[0][0];
            res.hNode.levelDist = 1.0;
            res.hNode.color = System.Drawing.Color.Black;
            for (int i=0;i<clusters.Count;i++)
            {
                HClusterNode aux = new HClusterNode();
                aux.setStruct = clusters[i];
                res.hNode.joined.Add(aux);
                res.hNode.setStruct.AddRange(clusters[i]);
                aux.parent = res.hNode;
                aux.refStructure = aux.setStruct[0];
                aux.levelDist = 0.0;
                aux.color = System.Drawing.Color.Black;
            }

            return res;

        }
        public ClusterOutput DendrogUsingMicroClusters(Dictionary<string, List<double>> dataFull,Dictionary<string, List<double>> data,Dictionary<string,List<string>> clustered)
         {
             ClusterOutput outC = null;
             switch (dMeasure)
             {
                 case DistanceMeasures.HAMMING:
                     if (refJuryProfile == null || !jury1d)
                         throw new Exception("Sorry but for jury measure you have to define 1djury profile to find reference structure");
                     dist = new JuryDistance(data,true);
                     //dist.InitMeasure();
                     break;


                 case DistanceMeasures.COSINE:
                         dist = new CosineDistance(data, jury1d);
                     break;
                 case DistanceMeasures.EUCLIDIAN:
                         dist = new Euclidian(data, jury1d);
                     break;

                 case DistanceMeasures.PEARSON:
                         dist = new Pearson(data, jury1d);
                     break;
                case DistanceMeasures.FASTHAMMING:
                    //dist = new FastDiscreteDist(1,data);
                    dist = new Fast3States(1, data);
                    break;
             }

             // return new ClusterOutput();
             DebugClass.WriteMessage("Start hierarchical");
             //Console.WriteLine("Start hierarchical " + Process.GetCurrentProcess().PeakWorkingSet64);
             currentV = maxV;
            //hk = new hierarchicalCluster(dist, opt);
             hierarchicalFast hk;
            hk = new hierarchicalFast(dist, opt);
             dist.InitMeasure();

             //Now just add strctures to the leaves   
             List<string> keys = new List<string>(data.Keys);

            outC = hk.HierarchicalClustering(new List<string>(data.Keys),null);
            //outC = hk.HierarchicalClustering(new List<string>(data.Keys), null);
            DebugClass.WriteMessage("Stop hierarchical");
             List<HClusterNode> hLeaves = outC.hNode.GetLeaves();
             
             foreach (var item in hLeaves)
             {
                 if (clustered.ContainsKey(item.refStructure))
                 {
                     item.setStruct.Clear();
                     item.setStruct.AddRange(clustered[item.refStructure]);

                    foreach (var str in item.setStruct)
                    {
                        if (dataFull.ContainsKey(str))
                            item.setProfiles.Add(dataFull[str]);
                    }

                     item.consistency = CalcClusterConsistency(item.setStruct);

                 }
                 else
                {
                    item.setStruct.Clear();
                    item.setStruct.Add(item.refStructure);
                }
             }
             outC.hNode.RedoSetStructures();
             outC.runParameters = hier.GetVitalParameters();
             outC.runParameters += input.GetVitalParameters();
             return outC;

         }
        Dictionary<string,List<int>> GetDummyKeys(List<string> structures)
        {
            Dictionary<string, List<int>> res = new Dictionary<string, List<int>>();
            Random r = new Random();
            char[] key = new char[100];
            for (int i=0;i<structures.Count;i++)
            {
                
                for (int j = 0; j < key.Length; j++)
                    if (r.Next(2) > 0)
                        key[j] = '1';
                    else
                        key[j] = '0';
                List<int> l = new List<int>();
                l.Add(i);
                string keyS = new string(key);
                if (!res.ContainsKey(keyS))
                    res.Add(keyS, l);
                else
                    i--;
            }


            return res;
        }
        Dictionary<string, List<int>> GetIEKeys(List<string> structures)
        {
            Dictionary<string, List<int>> res = new Dictionary<string, List<int>>();
            Random r = new Random();
            char[] key = new char[data[structures[0]].Count];
            for (int i = 0; i < structures.Count; i++)
            {

                for (int j = 0; j < data[structures[i]].Count; j++)
                    if (data[structures[i]][j] > 0)
                        key[j] = '1';
                    else
                        key[j] = '0';
                string keyS = new string(key);
                if (!res.ContainsKey(keyS))
                {
                    List<int> l = new List<int>();
                    l.Add(i);
                    res.Add(keyS, l);
                }
                else
                    res[keyS].Add(i);
            }


            return res;
        }
        Dictionary<string, List<int>> GetNoChangeKeys(List<string> structures)
        {
            Dictionary<string, List<int>> res = new Dictionary<string, List<int>>();
            Random r = new Random();
            char[] key = new char[data[structures[0]].Count];
            for (int i = 0; i < structures.Count; i++)
            {

                for (int j = 0; j < data[structures[i]].Count; j++)
                    key[j] = (char)('0' + data[structures[i]][j]);
                string keyS = new string(key);
                if (!res.ContainsKey(keyS))
                {
                    List<int> l = new List<int>();
                    l.Add(i);
                    res.Add(keyS, l);
                }
                else
                    res[keyS].Add(i);
            }


            return res;
        }

            public ClusterOutput DendrogUsingMeasures(List<string> structures)
         {

             jury1D juryLocal = new jury1D();
             juryLocal.PrepareJury(data);
             
             ClusterOutput outC = null;
             Dictionary<string, List<int>> dic;
             //Console.WriteLine("Start after jury " + Process.GetCurrentProcess().PeakWorkingSet64);
             maxV = 100;
             currentV = 0;
            if(input.relClusters>structures.Count)
            {
                //dic = GetDummyKeys(structures);
                //dic = GetIEKeys(structures);
                dic = GetNoChangeKeys(structures);
                StreamWriter st = new StreamWriter("test");
                foreach(var item in data)
                {
                    st.WriteLine($">{item.Key}");
                    for (int i = 0; i < item.Value.Count; i++)
                        st.Write($" {item.Value[i]}");
                    st.WriteLine();
                }
                st.WriteLine();
                st.Close();
            }
            else
             dic = PrepareKeys(structures,false,true);

             
             currentV+=5;
            //DebugClass.DebugOn();
            //  input.relClusters = input.reqClusters;
            //  input.perData = 90;
            Dictionary<string, string> translateToCluster = null;
            List<string> structuresToDendrogram = null;
            List<string> structuresFullPath = null;
            if (dic.Count > input.relClusters)
            {
                if (!input.combine)
                    dic = HashEntropyCombine(dic, structures, input.relClusters);
                else
                    //dic = FastCombineKeys(dic, structures, false);
                    dic = FastCombineKeysNew(dic, structures, false);
            }
                //Console.WriteLine("Entropy ready after jury " + Process.GetCurrentProcess().PeakWorkingSet64);
                DebugClass.WriteMessage("Entropy ready");
                //Alternative way to start of UQclust Tree must be finished
                //input.relClusters = 10000;


                //dic = FastCombineKeys(dic, structures, true);
                DebugClass.WriteMessage("dic size" + dic.Count);

                //Console.WriteLine("Combine ready after jury " + Process.GetCurrentProcess().PeakWorkingSet64);
                DebugClass.WriteMessage("Combine Keys ready");
                translateToCluster = new Dictionary<string, string>(dic.Count);
                structuresToDendrogram = new List<string>(dic.Count);
                structuresFullPath = new List<string>(dic.Count);
                DebugClass.WriteMessage("Number of clusters: " + dic.Count);
                int cc = 0;
                List<string> order = new List<string>(dic.Keys);
                //             order.Sort((a,b)=>dic[b].Count.CompareTo(dic[a].Count));
                order.Sort(delegate (string a, string b)
                {

                    if (dic[b].Count == dic[a].Count)
                        for (int i = 0; i < a.Length; i++)
                            if (a[i] != b[i])
                                if (a[i] == '0')
                                    return -1;
                                else
                                    return 1;


                    return dic[b].Count.CompareTo(dic[a].Count);
                });
                foreach (var item in order)
                {
                    if (dic[item].Count > 2)
                    {
                        List<string> cluster = new List<string>(dic[item].Count);                    
                        foreach (var str in dic[item])
                            cluster.Add(structures[str]);


                        ClusterOutput output = juryLocal.JuryOptWeights(cluster);
                        int end = output.juryLike.Count;
                        structuresToDendrogram.Add(output.juryLike[0].Key);
                        translateToCluster.Add(output.juryLike[0].Key, item);
                    }
                    else
                    {
                        structuresToDendrogram.Add(structures[dic[item][0]]);
                        translateToCluster.Add(structures[dic[item][0]], item);
                    }
                    cc++;
                }

             currentV+=10;
             DebugClass.WriteMessage("Jury finished");
             Regex rgx = new Regex("_similarity|_distance|profiles|profile");

             switch (dMeasure)
             {
                 case DistanceMeasures.HAMMING:
                     if (refJuryProfile == null || !jury1d)
                         throw new Exception("Sorry but for jury measure you have to define 1djury profile to find reference structure");
                         
                         dist = new JuryDistance(data, true);
                         break;


                 case DistanceMeasures.COSINE:
                            dist = new CosineDistance(data,jury1d);
                     break;
                 case DistanceMeasures.EUCLIDIAN:
                        dist = new Euclidian(data,jury1d);
                     break;

                 case DistanceMeasures.PEARSON:
                        dist=new Pearson(data,jury1d);
                     break;
                case DistanceMeasures.FASTHAMMING:
                    dist = new Fast3States(1, data);
                    //dist = new FastDiscreteDist(1,data);
                    break;


             }

            // return new ClusterOutput();
             DebugClass.WriteMessage("Start hierarchical");
             //Console.WriteLine("Start hierarchical " + Process.GetCurrentProcess().PeakWorkingSet64);
             currentV = maxV;
             hk = new hierarchicalCluster(dist, opt);
             dist.InitMeasure();
            
             //Now just add strctures to the leaves             
             outC = hk.HierarchicalClustering(structuresToDendrogram);
             DebugClass.WriteMessage("Stop hierarchical");
             List<HClusterNode> hLeaves = outC.hNode.GetLeaves();
            HashSet<string> checkNames = new HashSet<string>();
             foreach(var item in hLeaves)
             {
                 if (translateToCluster.ContainsKey(item.setStruct[0]))
                 {
                    foreach (var str in dic[translateToCluster[item.setStruct[0]]])
                        if (item.setStruct[0] != structures[str])
                        {
                            item.setStruct.Add(structures[str]);
                        }

                             foreach (var str in item.setStruct)
                             {
                                 item.setProfiles.Add(data[str]);
                             }
                         
                     


                     item.consistency = CalcClusterConsistency(item.setStruct);

                 }
                 else
                     throw new Exception("Cannot add structure. Something is wrong");
             }
             outC.hNode.RedoSetStructures();
            outC.profiles = data;
            foreach(var item in hLeaves)
                foreach(var it in item.setStruct)
                    if (!checkNames.Contains(it))
                        checkNames.Add(it);
                    else
                        Console.WriteLine("UPPPS");

            outC.runParameters = hier.GetVitalParameters();
             outC.runParameters += input.GetVitalParameters();
             return outC;
         }
         

    }
}
