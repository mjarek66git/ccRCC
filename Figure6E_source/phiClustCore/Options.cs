using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Xml.Serialization;
using System.Reflection;
using System.Linq.Expressions;
using System.IO;
using phiClustCore;
using phiClustCore.Interface;
using System.ComponentModel.DataAnnotations;

namespace phiClustCore
{
    public enum ClusterAlgorithm
    {
        [Display(Name ="Kmeans")]
        HKmeans,
        HNN,
        HierarchicalCluster,
        Jury1D,
        HashCluster,
        GuidedHashCluster,
        uQlustTree,
        OmicsHeatMap,
        [Display(Name = "Fast Jury")]
        FastJuryBased,        
        HTree,
     
    }
    public enum DistanceMeasures
    {
        HAMMING,
        COSINE,
        EUCLIDIAN,
        PEARSON,
        FASTHAMMING,
        NONE
    }
    public enum AglomerativeType
    {
        SINGLE,
        AVERAGE,
        COMPLETE,
        REFERENCE
    }
    public class DescriptionAttribute : Attribute
    {
        private string description;

        public string Description { get { return description; } }

        public DescriptionAttribute(string description)
        {
            this.description = description;
        }
    }
    [Serializable]
    [Description("Opcje")]
    public class Options:IAutomaticProfiles
    {
        [NonSerialized]
        public static string defaultFileName = "default.options";

        //[DescriptionAttribute("Path to settings file")]
        //public string settingsFile;
        [Description("Path to data directory")]
        public List<string> dataDir = new List<string>();
        [Description("List of aligned profile files")]
        public List<string> profileFiles = new List<string>();
        [Description("Output file name")]
        public string outputFile;
        [Description("PDB file name, needed when dcd file is be used")]
        public string pdbFileName;
        [Description("Temporary directory, needed when dcd file is used")]
        public string tempDirectory;
        [Description("Path to profiles definition file, needed when 1djury is used")]
        public string profiles1DJuryFile;
        [Description("Path to profiles definition file, needed when 1djury measure is used")]
        public string profiles1DJuryMeasureFile;
        [Description("Cluster algorithm to be used (available: uQlustTree, HKmeans, HierarchicalCluster, FastHCluster, Kmeans, ThresholdCluster, Jury1D, Jury3D, Sift)")]
        public List<ClusterAlgorithm> clusterAlgorithm=new List<ClusterAlgorithm>();
        [Description("Name of distance file")]
        public string distanceFile = "";
        [Description("Save object to binary file")]
        public string binaryFile = "";

        //[Description("Distance measure to be used [MAXSUB, RMSD, JURY1D]")]
        public RankingCInput other=new RankingCInput();
        public HierarchicalCInput hierarchical=new HierarchicalCInput();
        public ThresholdCInput threshold=new ThresholdCInput();        
        public KmeansInput kmeans=new KmeansInput();        
        public HashCInput hash = new HashCInput();
        public OmicsInput omics = new OmicsInput();
        public HNNCInput hnn = new HNNCInput();
        [NonSerialized]
        private Dictionary<string, string> dicField=new Dictionary<string,string>();
        [NonSerialized]
        private Dictionary<string, MemberInfo> dicMem=new Dictionary<string,MemberInfo>();
        public DistanceMeasures GetDistanceMeasure(ClusterAlgorithm clAlgorithm)
        {
            switch (clAlgorithm)
            {
                case ClusterAlgorithm.HierarchicalCluster:
                case ClusterAlgorithm.HKmeans:
                            return hierarchical.distance;
                default:
                            return DistanceMeasures.NONE;

            }
        }
        public Options()
        {
            PrepareAllFields();
        }

        public void GenerateAutomaticProfiles(string fileName)
        {
            Type t = this.GetType();
            MemberInfo [] members = t.GetMembers();

            if (profileFiles == null || profileFiles.Count == 0)
                return;

            foreach (MemberInfo mem in members)
            {
                if (mem.MemberType.ToString() != "Field")
                    continue;
                if (typeof(IAutomaticProfiles).IsAssignableFrom(mem.ReflectedType.GetField(mem.Name).FieldType))
                {
                    IAutomaticProfiles o = (IAutomaticProfiles)mem.ReflectedType.GetField(mem.Name).GetValue(this);
                    o.GenerateAutomaticProfiles(profileFiles[0]);
                }
            }

        }
        private void PrepareAllFields()
        {
            Type t = this.GetType();
            MemberInfo[] members = t.GetMembers();
            dicField.Clear();
            dicMem.Clear();

            foreach (MemberInfo mem in members)
            {
                object[] attributes = mem.GetCustomAttributes(true);

                if (attributes.Length != 0 && mem.MemberType.ToString() == "Field")
                {
                    string key = "";
                    foreach (object attribute in attributes)
                    {

                        //Console.Write("  {0} ", attribute.ToString());
                        DescriptionAttribute da = attribute as DescriptionAttribute;
                        if (da != null)
                            key= da.Description;
                    }
                    object o = mem.ReflectedType.GetField(mem.Name).GetValue(this);
                    if (o != null)
                        dicField.Add(key, mem.ReflectedType.GetField(mem.Name).GetValue(this).ToString());                    
                    else                    
                        dicField.Add(key, "");
                    
                    dicMem.Add(key, mem);

                }
            }
        }
        public void SaveDefault()
        {
            SaveOptions(PrepareDefaultName());
        }
        public void SaveOptions(string fileName)
        {
            StreamWriter file = new StreamWriter(fileName);
            PrepareAllFields();
            foreach (var item in dicField.Keys)
            {
                
                if (dicField[item].Contains("Generic.List"))
                {
                    file.Write(item + "# ");
                        if (item.Contains("profile"))
                        {
                            foreach (var itemStr in profileFiles)
                                file.Write(itemStr + ";");
                        }
                        if(item.Contains("algorithm"))
                        {
                            foreach (var it in clusterAlgorithm)
                                file.Write(it + ";");
                        }
                        if (item.Contains("Path"))
                            foreach (var itemStr in dataDir)
                                file.Write(itemStr + ";");
                    file.WriteLine();
                }
                else
                    file.WriteLine(item + "# " + dicField[item]);
            }

            if (hierarchical != null)
            {
                file.WriteLine("########## Parameters hierarchical ##########");

                hierarchical.SaveOptions(file);
                file.WriteLine("########## End ##########");
            }
            if (threshold != null)
            {
                file.WriteLine("########## Parameters threshold ##########");
                threshold.SaveOptions(file);
                file.WriteLine("########## End ##########");
            }
            if (kmeans != null)
            {
                file.WriteLine("########## Parameters kmeans ##########");
                kmeans.SaveOptions(file);
                file.WriteLine("########## End ##########");
            }
            if (other != null)
            {
                file.WriteLine("########## Parameters other ##########");
                other.SaveOptions(file);
                file.WriteLine("########## End ##########");
            }
            if (hash != null)
            {
                file.WriteLine("########## Parameters Hash ##########");
                hash.SaveOptions(file);
                file.WriteLine("########## End ##########");
            }
            if(omics!=null)
            {
                file.WriteLine("########## Parameters Omics ##########");
                omics.SaveOptions(file);
                file.WriteLine("########## End ##########");

            }
            if(hnn!=null)
            {
                file.WriteLine("########## Parameters HNN ##########");
                hnn.SaveOptions(file);
                file.WriteLine("########## End ##########");
            }
            file.Close();
        }
        private string PrepareDefaultName()
        {
            //string dir = Environment.GetEnvironmentVariable("CLUSTER");
            string dir = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location);
            if (dir != null && !defaultFileName.Contains(Path.DirectorySeparatorChar))
                defaultFileName = dir + Path.DirectorySeparatorChar + defaultFileName;
            return defaultFileName;
        }
        public void ReadDefaultFile()
        {
            ReadOptionFile(PrepareDefaultName());
        }
        public void ReadOptionFile(string fileName)
        {
            if (File.Exists(fileName))
            {   
                string line;
                MemberInfo memB;
                StreamReader r = new StreamReader(fileName);
               // DebugClass.DebugOn();
                DebugClass.WriteMessage("Start reading " + fileName);
                while(!r.EndOfStream)
                {
                    line = r.ReadLine();
                    if (line.Contains("#####"))
                    {
                        if (line.Contains("Parameters hierarchical"))
                            hierarchical.ReadOptionFile(r);
                        else
                            if (line.Contains("Parameters threshold"))
                            threshold.ReadOptionFile(r);
                        else
                                if (line.Contains("Parameters other"))
                            other.ReadOptionFile(r);
                        else
                                    if (line.Contains("Parameters Hash"))
                            hash.ReadOptionFile(r);
                        else
                                        if (line.Contains("Parameters kmeans"))
                            kmeans.ReadOptionFile(r);
                        else
                                            if (line.Contains("Omics"))
                            omics.ReadOptionFile(r);
                        else
                                                if (line.Contains("HNN"))
                            hnn.ReadOptionFile(r);
                                            
                                        

                       
                    }

                    
                    line=line.Replace("# ", "#");
                    line=line.TrimEnd(new char []{'\r','\n'});
                    string[] strTab = line.Split('#');
                    if (strTab[1].Length == 0)
                        continue;

                    if (dicField.ContainsKey(strTab[0]))
                    {
                        if (dicField[strTab[0]].Contains("Generic.List"))
                        {
                            if(strTab[1].EndsWith(";"))
                                strTab[1] = strTab[1].Remove(strTab[1].Length - 1, 1);
                            string[] aux = strTab[1].Split(';');

                                if (strTab[0].Contains("prof"))
                                {
                                    profileFiles.Clear();
                                    foreach (var item in aux)
                                        profileFiles.Add(item);
                                
                                }
                                    if(strTab[0].Contains("algorithm"))
                                    {
                                        clusterAlgorithm.Clear();
                                        foreach (var item in aux)
                                        {
                                            clusterAlgorithm.Add((ClusterAlgorithm)Enum.Parse(typeof(ClusterAlgorithm), item));
                                        }                                        
                                    }
                                if(strTab[0].Contains("data"))
                                {
                                    dataDir.Clear();
                                    foreach (var item in aux)
                                        dataDir.Add(item);
                                }
                         
                                continue;
                            
                        }
                        memB = dicMem[strTab[0]];
                        string ww = memB.ReflectedType.GetField(memB.Name).FieldType.Name;
                        switch (ww)
                        {
                            case "string":
                            case "String":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, strTab[1]);
                                break;
                            case "Boolean":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Convert.ToBoolean(strTab[1]));
                                break;
                            case "Single":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Convert.ToSingle(strTab[1]));
                                break;
                            case "Int32":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Convert.ToInt32(strTab[1]));
                                break;
                            case "Initialization":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Enum.Parse(typeof(Initialization), strTab[1]));
                                break;
                            case "ClusterAlgorithm":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Enum.Parse(typeof(ClusterAlgorithm), strTab[1]));
                                break;
                            case "DistanceMeasures":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Enum.Parse(typeof(DistanceMeasures), strTab[1]));
                                break;
                            case "AglomerativeType":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Enum.Parse(typeof(AglomerativeType), strTab[1]));
                                break;
                            case "COL_SELECTION":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Enum.Parse(typeof(COL_SELECTION), strTab[1]));
                                break;
                        }

                        //SetValue(this, strTab[1]);
                    }
                    else
                        DebugClass.WriteMessage("Not recognized: " + strTab[0]);
                }
                
                r.Close();
            }
			else
				throw new Exception("Config file not found");

            DebugClass.WriteMessage("Reading finished");

        }

    }

}
