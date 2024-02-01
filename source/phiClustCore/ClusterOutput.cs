using System;
using System.Diagnostics;
using System.IO;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;
using System.Collections;
using System.Collections.Specialized;
using System.Collections.Generic;
using System.Xml.Serialization;
using System.Drawing;
using phiClustCore.Distance;

namespace phiClustCore
{
    public interface ICluster
    {
        ClusterOutput RunCluster();
    }
    [Serializable]
    public class clusterRes
    {
        public clusterRes()
        {
            list = new List<List<string>>();
        }
        public List<List<string>> list = null;
        public List<string> labels = null;
        public List<double> consistency = null;


    }
    [Serializable]
    public class ClusterOutput
    {
        public OmicsDataSet data;
        public List<KeyValuePair<string, double>> juryLike;
        public Dictionary<string, string> hNNRes;
        public List<string> aux1;
        public List<string> aux2;
        public List<int> auxInt;
        public HClusterNode hNode;
        public List<HClusterNode> nodes;
        public clusterRes clusters = null;
        public clusterRes clusters2 = null;
        public Dictionary<string, List<double>> profiles = null;
        public DistanceMeasure distM=null;


        public clusterRes ClustersUpper {get{return clusters;} }
        public clusterRes ClustersLeft{get{return clusters2;} }

        public string clusterType;
        public string measure;
        public string time;
        public string dirName;
        public string name;
        public long peekMemory;
        public string runParameters="";

        public void SaveTxt(string fileName)
        {
            StreamWriter stream = new StreamWriter(fileName);
            if (runParameters.Length > 0)
                stream.WriteLine("======" + runParameters + "=====");
            if (juryLike != null)
            {               
                foreach(var item in juryLike)                
                    stream.WriteLine(item.Key+" "+item.Value);                
            }
            if(clusters!=null)
            {
                Save(clusters.list,clusters.consistency,stream,false,clusters.labels,auxInt);
            }
            if (hNode != null)
            {
                //Save(hNode.GetClusters(10),stream,true);
                ExportToAtr(hNode,Path.GetFileNameWithoutExtension(fileName)+".atr");
            }
            if(hNNRes!=null)
            {
                SaveHnn(hNNRes,fileName+"_prediction");
            }
            stream.Close();
        }
        public void SaveGraph(string fileName,int Width,int Height)
        {
            if (juryLike != null)
                throw new Exception("1D jury does not have graphical output");
            if (nodes != null && nodes.Count > 1)
            {
                HeatMapDraw heatMap = new HeatMapDraw(new Bitmap(Width, Height), nodes[1], nodes[0], null, this);
                heatMap.SavePaint(fileName);
            }
            else
                if(hNode!=null)
                {
                Dictionary<string, string> dic = ClusterOutput.ReadLabelsFile(GetLabelFile());
                Bitmap bit = new Bitmap(Width,Height);
                DrawHierarchical h = new DrawHierarchical(hNode, measure, dic, bit, true);
                h.Save(fileName, bit);
                }
            else
                throw new Exception("Graphics can be saved only for heatmap!");
        }
        static public void Save(List<List<string>> list,List<double> clCons,StreamWriter stream,bool reference,List<string> labels=null,List<int> selectedFeatures=null)
        {
            if (selectedFeatures != null)
            {
                stream.Write("Indexes of selected Features: ");
                foreach(var item in selectedFeatures)
                    stream.Write(item + " ");
                stream.WriteLine();
            }
            for(int i=0;i<list.Count;i++)
            {
                string line = "===========CLUSTER " + (i+1) + " === Size " + list[i].Count;
                if (labels != null)
                    line += " ====label:"+labels[i];
                if (reference)
                    line += "====Reference model:" + list[i][0];
                if (clCons != null && clCons.Count == list.Count)
                    line += "====Consitsency:" + clCons[i];
                stream.WriteLine(line);
               
                foreach (var it in list[i])
                    stream.WriteLine(it);
                stream.WriteLine("=======================================");
            }
        }
        static public void SaveHnn(Dictionary<string,string>hnnRes,string fileName)
        {
            StreamWriter wr = new StreamWriter(fileName);
            foreach(var item in hnnRes)
            {
                wr.WriteLine(item.Key + " " + item.Value);
            }
            wr.Close();
        }
        static public void ExportToAtr(HClusterNode node,string fileName)
        {
            List<HClusterNode > leafs=node.GetLeaves();
            List<List<string>> clusters = new List<List<string>>();

            foreach(var cl in leafs)
            {
                clusters.Add(cl.setStruct);
            }
            
            string newName = Path.GetFileNameWithoutExtension(fileName)+"_leafs.dat";
            StreamWriter stream = new StreamWriter(newName);
            Save(clusters,null,stream,true);
            stream.Close();
            stream = new StreamWriter(fileName);
            Queue<HClusterNode> ww = new Queue<HClusterNode>();
            Dictionary<HClusterNode, string> nodeName = new Dictionary<HClusterNode, string>();
            Dictionary<string,int> leafNum=new Dictionary<string,int>();
            int counter=1;

            ww.Enqueue(node);
            while (ww.Count != 0)
            {
                HClusterNode aux = ww.Dequeue();
                if(aux.joined!=null)
                    foreach (var item in aux.joined)
                    {
                        ww.Enqueue(item);
                        item.parent = aux;
                    }

            }


            foreach (var item in leafs)
            {
                //ww.Enqueue(item);
                leafNum.Add(item.refStructure,counter++);
            }
            ww.Enqueue(node);
            int nodeCounter = 0;
            while(ww.Count!=0)
            {
                HClusterNode aux = ww.Dequeue();


                string parentNode;
                if (!nodeName.ContainsKey(aux))
                {
                    parentNode = "Node_" + (++nodeCounter);
                    nodeName.Add(aux, parentNode);
                }
                else
                    parentNode = nodeName[aux];
                if(aux.joined!=null)
                {
                List<HClusterNode> hList = aux.joined;
                    
                    for (int i = 0; i < hList.Count;i++ )
                    {
                        string name1;
                        if (nodeName.ContainsKey(hList[i]))
                            name1 = nodeName[hList[i]];
                        else
                        {

                            name1 = hList[i].refStructure + "_" + leafNum[hList[i].refStructure];
                            if (hList[i].joined != null)
                                name1 = "Node_" + (++nodeCounter); 
                            nodeName.Add(hList[i], name1);
                        }
                        for (int j = i+1; j < hList.Count; j++)
                        {
                            string  name2;
                            if (nodeName.ContainsKey(hList[j]))
                                name2 = nodeName[hList[j]];
                            else
                            {

                                name2 = hList[j].refStructure + "_" + leafNum[hList[j].refStructure];
                                if (hList[j].joined != null)
                                    name2 = "Node_" + (++nodeCounter);
                                nodeName.Add(hList[j], name2);
                            }

                            stream.WriteLine(parentNode+" "+name1+" " + name2 + " " + aux.realDist);
                        }
                        ww.Enqueue(aux.joined[i]);
                    }
                
                }
            

            }
            stream.Close();

        }
        static public Dictionary<string, string> ReadLabelsFile(string fileName)
        {
            Dictionary<string, string> vecLabels = new Dictionary<string, string>();

            if (fileName == null)
                return null;

            StreamReader str = new StreamReader(fileName);
            string line = str.ReadLine();
            while (line != null)
            {
                string[] aux = line.Split(' ');
                if (aux.Length >= 2)
                {
                    if (aux[0].Length > 0 && aux[1].Length > 0)
                    {
                        
                        if (aux[0].Contains(Path.DirectorySeparatorChar.ToString()))
                        {
                            string[] xx = aux[0].Split(Path.DirectorySeparatorChar);
                            aux[0] = xx[xx.Length - 1];
                        }
                        vecLabels[aux[0]] = aux[1];
                    }
                }
                line = str.ReadLine();

            }
            str.Close();

            return vecLabels;

        }

        public string GetLabelFile()
        {
            string labelFile;

            if(File.Exists(dirName))
            {
                labelFile =Path.GetDirectoryName(dirName)+Path.DirectorySeparatorChar+ Path.GetFileNameWithoutExtension(dirName);
                labelFile+="_label.dat";
                if (File.Exists(labelFile))
                    return labelFile;
                return null;
            }
            else
            {
                labelFile = dirName + "_label.dat";
                if (File.Exists(labelFile))
                    return labelFile;

                return null;
            }
        }
        static public ClusterOutput Load(string fileName)
        {
            ClusterOutput outObject;
            Stream stream = File.Open(fileName, FileMode.Open);
            BinaryFormatter bFormatter = new BinaryFormatter();
            outObject=(ClusterOutput) bFormatter.Deserialize(stream);
            return outObject;
        }
        static public ClusterOutput LoadExternal(string fileName,string dirName)
        {
            ClusterOutput aux = new ClusterOutput() ;
            string line;
            StreamReader r = new StreamReader(fileName);
            aux.clusters = new clusterRes();
            aux.dirName = dirName;
            while (!r.EndOfStream)
            {
                line = r.ReadLine();
                if(!line.Contains("#") && line.Contains(":"))
                {
                    line = line.Replace("  ", "");
                    line = line.Replace("  ", "");
                    line=line.Replace(" : ", ":");
                    line=line.TrimEnd(' ');
                    string[] tmp = line.Split(' ');
                    List<string> auxList = new List<string>();
                 
                    string[] tmp2 = tmp[1].Split(':');
                    string[] loc = tmp2[1].Split('/');
                    
                    //aux.dirName = "/home/raad/work/casp10/" + loc[loc.Length - 2];
                    aux.name = loc[loc.Length - 2];
                    auxList.Add(loc[loc.Length - 1]+".pdb");
                    for (int i = 4; i < tmp.Length; i++)
                    {
                        loc = tmp[i].Split('/');
                        if (!auxList.Contains(loc[loc.Length - 1]))
                            auxList.Add(loc[loc.Length - 1]+".pdb");
                    }

                    aux.clusters.list.Add(auxList);
                }

            }
            return aux;
        }
        static public ClusterOutput LoadExternalPconsD(string fileName,string dirName)
        {
            ClusterOutput aux = new ClusterOutput();
            string[] tmp = null; ;
            string line;
            StreamReader r = new StreamReader(fileName);
            aux.clusters = new clusterRes();
            List<KeyValuePair<string,double>> auxList = new List<KeyValuePair<string,double>>();
            line = r.ReadLine();
            while (!r.EndOfStream)
            {
                
                tmp = line.Split(' ');
                if (tmp.Length>0)
                {
                    KeyValuePair<string, double> auxK = new KeyValuePair<string, double>(tmp[0], Convert.ToDouble(tmp[1]));
                    auxList.Add(auxK);
            
                }
                line = r.ReadLine();
            }
            aux.dirName = dirName;
            
            tmp =fileName.Split('\\');
            aux.name = Path.GetFileName(tmp[tmp.Length - 1]);
            aux.clusters = null;
            aux.juryLike = auxList;
            return aux;
        }
        static public ClusterOutput LoadExternalPleiades(string fileName,string dirName)
        {
            ClusterOutput aux = new ClusterOutput();
            string line;
            StreamReader r = new StreamReader(fileName);
            aux.clusters = new clusterRes();
            List<KeyValuePair<int, int>> clustSize = new List<KeyValuePair<int, int>>();
            aux.dirName = dirName;
            List<string> auxList = null;
            int clSize = 0,index=0;
            while (!r.EndOfStream)
            {
                line = r.ReadLine();
                if (line.Contains("Cluster"))
                {
                    if (auxList != null)
                    {
                        aux.clusters.list.Add(auxList);
                        clustSize.Add(new KeyValuePair<int,int>(clSize,index));
                        index++;
                    }

                    string[] aa = line.Split(',');
                    aa = aa[0].Split(':');
                    clSize = Convert.ToInt32(aa[1]);
                    auxList = new List<string>();
                }
                else
                    if (line.Contains("/") && !line.Contains("#")) 
                    {
                        line = line.TrimStart(' ');
                        line = line.Replace("\t", " ");
                        line = line.Replace("  ", " ");
                        line = line.Replace("  ", " ");
                        string[] tmp = line.Split(' ');

                        string[] loc = tmp[1].Split('/');
                        string[] xx = Path.GetFileName(fileName).Split('_');
                        //aux.dirName = Path.GetDirectoryName(fileName)+Path.DirectorySeparatorChar+xx[0];
                        //aux.dirName = "/home/raad/work/casp10/" + loc[loc.Length - 2];
                        aux.name = loc[loc.Length - 2];
                        auxList.Add(loc[loc.Length - 1]);
                    }

            }
            if (auxList.Count > 0)
            {
                aux.clusters.list.Add(auxList);
                clustSize.Add(new KeyValuePair<int, int>(clSize, index));
            }

            clustSize.Sort(delegate(KeyValuePair<int, int> first, KeyValuePair<int, int> second) { return second.Key.CompareTo(first.Key); });

            ClusterOutput final = new ClusterOutput();
            final.clusters = new clusterRes();
            final.name = aux.name;
            final.dirName = aux.dirName;
            foreach (var item in clustSize)
                final.clusters.list.Add(aux.clusters.list[item.Value]);

            return final;
        }


        public List<List<string>> GetClusters(int clustNum=0)
        {
            List<List<string>> aux=new List<List<string>>();
            if (clusters != null)
                return clusters.list;

            if (juryLike != null)
            {
                List<string> lista = new List<string>(juryLike.Count);
                foreach (var item in juryLike)
                    lista.Add(item.Key);
                aux.Add(lista);

                return aux;

            }
            if (hNode != null)
                return hNode.GetClusters(clustNum);


            return null;
        }

    }

    public enum INPUTMODE
    {
        //PROTEIN,
        //RNA,
        USER_DEFINED,
        OMICS
    };
    [Serializable]
    public class Settings
    {
        [NonSerialized]
        public static string gFileName = "Settings.clu";
        public string dir = "";
      
        public int numberOfCores { get; set; }
        public string extension { get; set; }
        public string profilesDir { get; set; }
        public INPUTMODE mode { get; set; }
        public bool iOTroubles { get; set; }
        public Settings()
        {
            //mode = INPUTMODE.PROTEIN;
            mode = INPUTMODE.USER_DEFINED;
            numberOfCores = 1;
            extension = "*";
            iOTroubles = false;
            dir = Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location);
            string dirTmp = Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location) + Path.DirectorySeparatorChar + "generatedProfiles";
            if (Directory.Exists(dirTmp))
                profilesDir = dirTmp;
            else
                profilesDir = Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location);
            
            if (!gFileName.Contains(Path.DirectorySeparatorChar.ToString()))
            {
                gFileName = dir + Path.DirectorySeparatorChar + gFileName;
            }
        }
        public Settings(string fileName)
        {
            Load(gFileName);
        }

        public void Load()
        {
            Load(gFileName);
        }

        public void Load(string fileName)
        {
            if (!File.Exists(fileName))
            {
                mode = INPUTMODE.USER_DEFINED;
                //mode = INPUTMODE.PROTEIN;
                numberOfCores = 1;
                extension = "*.pdb";
                return;
            }
                //throw new Exception("File "+fileName+" cannot be found");
            XmlSerializer ser = new XmlSerializer(typeof(Settings));
             StreamReader reader = new StreamReader(fileName);
            Settings localSettings = (Settings)ser.Deserialize(reader);                  
            profilesDir = localSettings.profilesDir;
            extension = localSettings.extension;
            numberOfCores = localSettings.numberOfCores;
            mode = localSettings.mode;
            iOTroubles = localSettings.iOTroubles;
            reader.Close();
        }
        public void Save()
        {
            StreamWriter writer = new StreamWriter(gFileName);
            XmlSerializer ser = new XmlSerializer(this.GetType());
            ser.Serialize(writer, this);
            writer.Close();
        }

    }

}

