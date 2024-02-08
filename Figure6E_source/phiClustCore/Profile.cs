using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;
using System.Xml.Serialization;
using System.Text.RegularExpressions;
using System.Data;
using phiClustCore.Interface;
using phiClustCore.Profiles;

namespace phiClustCore
{

    public enum wOpertion
    {
        SUM,
        MULT
    };
    public enum profileType
    {
        INTERNAL,
        EXTERNAL
    };
    [Serializable]
    public class ProfileTree:IProgressBar
    {
        string defaultProfile="default.profiles";
        public string listFile;
        public string currentProfileName = "";
        Settings set = new Settings();       
        public SerializableDictionary<string, profileNode> masterNode = new SerializableDictionary<string, profileNode>();
        public Dictionary<string,Dictionary<string, protInfo>> profiles = null;
        public Dictionary<string, List<byte>> protCombineStates = null;
        Dictionary<string, byte> protCombineCoding = new Dictionary<string, byte>();
        public InternalProfilesManager mgr;
        IProgressBar progressObject = null;
        List<double> profilesProgress = null;
        int profilesNum = 0;

        double startProgress = 0;
        double endProgress = 1;

        public double StartProgress { set { startProgress = value; } get { return startProgress; } }
        public double EndProgress { set { endProgress = value; } get { return endProgress; } }

        public ProfileTree()//(string profileFile)
        {
            mgr = new InternalProfilesManager();
            progressObject = mgr;
            set.Load();
        }

        public double ProgressUpdate()
        {
            double sum = 0;
            if (profilesProgress == null)
                return 0;

            foreach (var item in profilesProgress)
            {
                sum += item;
            }

            sum /= profilesNum*2;

            if(progressObject!=null)
                sum=sum+progressObject.ProgressUpdate()/2;
           

            return StartProgress+(EndProgress - StartProgress) *sum;
        }
        public Exception GetException()
        {
            return null;
        }
        public List<INPUTMODE> GetModes()
        {
            Dictionary<INPUTMODE,int> dic=new Dictionary<INPUTMODE,int>();
            List<INPUTMODE> modes=new List<INPUTMODE>();

            foreach(var item in masterNode)
            {

                profileNode n = InternalProfilesManager.GetNode(item.Value.internalName);
                if (n == null)
                    return null;
                Type t = Type.GetType(InternalProfilesManager.internalList[n]);
                InternalProfileBase c = Activator.CreateInstance(t) as InternalProfileBase;
                foreach (var it in c.destination)
                    if (!dic.ContainsKey(it))
                        dic.Add(it, 1);
                    else
                        dic[it]++;
                            
            }
            foreach(var item in dic)
            {
                if (item.Value == masterNode.Keys.Count)
                    modes.Add(item.Key);
            }
            return modes;
        }
        public void ClearTree()
        {
            if(masterNode!=null)
                masterNode.Clear();
            if(profiles!=null)
                profiles.Clear();
            if(protCombineStates!=null)
                protCombineStates.Clear();
            if(protCombineCoding!=null)
                protCombineCoding.Clear();
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }

        public void AdddNode(string Path,profileNode newNode)
        {
            profileNode aux;

            if (Path==null || Path.Length == 0 || Path == "/")
            {
                masterNode.Add(newNode.profName,newNode);
                return;
            }

            aux = FindNode(Path);
            if (aux != null && !aux.childrens.ContainsKey(newNode.profName))                
                aux.childrens.Add(newNode.profName, newNode);
            else
                throw new Exception("New node cannot be added, wrong Path!");

        }
        public void RemoveOutFiles()
        {
            Queue<profileNode> nodeList = new Queue<profileNode>();
            profileNode aux;

            foreach (var item in masterNode.Values)
                    nodeList.Enqueue(item);

            while(nodeList.Count>0)
            {
                aux = nodeList.Dequeue();

                if (File.Exists(aux.OutFileName) && (!aux.OutFileName.Contains("dssp")))
                {
                    if(aux.removeOutFile)
                        File.Delete(aux.OutFileName);
                }

                foreach (var item in aux.childrens.Values)
                            nodeList.Enqueue(item);
            }
            
        

        }
        public Dictionary<byte, Dictionary<byte, double>> GenerateWeights(wOpertion weightsOperation)
        {
            List<profileNode> active = GetActiveLeaves();
            Dictionary<string,int> dic = new Dictionary<string,int>();
            List<List<string>> item1 = new List<List<string>>();
            List<string> item2;
            Dictionary<byte, Dictionary<byte, double>> weights = new Dictionary<byte, Dictionary<byte, double>>();
            List <profileNode> activeNodes=new List<profileNode>();

            for (int i = 0; i < active.Count; i++)
            {
                dic.Clear();
                foreach (var itemK in active[i].profWeights.Keys)
                {
                    if (!dic.ContainsKey(itemK))
                        dic.Add(itemK, 0);
                    foreach(var itemK1 in active[i].profWeights[itemK].Keys)
                        if (!dic.ContainsKey(itemK1))
                            dic.Add(itemK1, 0);
                }

                List<string> aux = new List<string>(dic.Keys);
                item1.Add(aux);                
            }
            item2 = new List<string>(item1[0]); ;
            for (int i = 1; i < item1.Count; i++)
            {
                List<string> newList = new List<string>();
                for (int k = 0; k < item2.Count; k++)
                {
                    string lan = item2[k];
                    for (int j = 0; j < item1[i].Count; j++)
                    {
                        lan = item2[k]+"-" + item1[i][j];
                        newList.Add(lan);
                    }
                }
                item2 = new List<string>(newList);
            }

            for (int i = 0; i < item2.Count; i++)
            {
                string[] aux1 = item2[i].Split('-');
                for (int j = 0; j < item2.Count; j++)
                {
                    string[] aux2 = item2[j].Split('-');
                    for (int k = 0; k < aux1.Length; k++)
                    {
                        double wValue=0;
                        switch (weightsOperation)
                        {
                            case wOpertion.MULT:
                                wValue = 1;
                                break;
                            case wOpertion.SUM:
                                wValue = 0;
                                break;
                        }
                        for (int n = 0; n < active.Count; n++)
                        {
                            if (active[n].profWeights.ContainsKey(aux1[n]))
                                switch (weightsOperation)
                                {
                                    case wOpertion.MULT:
                                        if (active[n].profWeights[aux1[n]].ContainsKey(aux2[n]))
                                            wValue *= active[n].profWeights[aux1[n]][aux2[n]];
                                        else
                                            wValue = 0;
                                        break;
                                    case wOpertion.SUM:
                                        if (active[n].profWeights[aux1[n]].ContainsKey(aux2[n]))
                                            wValue += active[n].profWeights[aux1[n]][aux2[n]];
                                        break;
                                }   

                        }
                        if (wValue > 0)
                        {
                            
                            string[] taux1 =new string [aux1.Length];
                            string[] taux2 =new string [aux2.Length];
                            for (int n = 0; n < aux1.Length; n++)
                            {
                                 taux1[n]= active[n].GetCodedState(aux1[n]).ToString();
                                 taux2[n] = active[n].GetCodedState(aux2[n]).ToString();
                                                        }
                            if (!protCombineCoding.ContainsKey(String.Join("", taux1)) || !protCombineCoding.ContainsKey(String.Join("", taux2)))
                                continue;
                            byte k1=protCombineCoding[String.Join("",taux1)];
                            byte k2=protCombineCoding[String.Join("",taux2)];
                            if (!weights.ContainsKey(k1))
                            {
                                Dictionary<byte, double> sAux = new SerializableDictionary<byte, double>();
                                sAux.Add(k2, wValue);
                                weights.Add(k1, sAux);
                            }
                            else
                            {
                                Dictionary<byte, double> sAux = weights[k1];
                                if (!sAux.ContainsKey(k2))
                                    sAux.Add(k2, wValue);
                            }
                        }
                    }
                }
            }
            return weights;
        }
        private string ParentPath(string path)
        {
            string[] tmp = path.Split('/');
            path = String.Join("/", tmp, 0, tmp.Length - 1);
            if (path.Length == 0)
                path = "/";
            return path;

        }

        public bool CheckIfAllInternal()
        {
            Queue<profileNode> nodeList = new Queue<profileNode>();
            profileNode aux = null;
            foreach (profileNode item in masterNode.Values)
            {
                if (item.GetNumberofStates() == 0)
                    return false;
                nodeList.Enqueue(item);
                
            }
            while (nodeList.Count > 0)
            {
                aux = nodeList.Dequeue();
                foreach (var item in aux.childrens.Values)
                {
                    if (item.GetNumberofStates() == 0)
                        return false;
                    nodeList.Enqueue(item);
                }
            }
            
            return true;
        }
		public void ActivateParents(string Path)
		{
			string currentPath=Path;
			profileNode node;
			
			while((currentPath=ParentPath (currentPath))!="/")
			{
				node=FindNode(currentPath);
				if(node!=null)
					node.active=true;
					
			}
		}
		public void ActivateProfiles(List <string> profiles)
		{
            Queue<KeyValuePair<profileNode,string>> nodeList = new Queue<KeyValuePair<profileNode,string>>();
			KeyValuePair<profileNode,string> aux,auxP;

			
            foreach (profileNode item in masterNode.Values)
			{
				aux=new KeyValuePair<profileNode, string>(item,"/"+item.profName);
				//aux.Key=item;
				//aux.Value="/"+item.profName;
            	nodeList.Enqueue(aux);
			}

            while(nodeList.Count>0)
            {
                
                aux = nodeList.Dequeue();
				
				if(profiles.Contains(aux.Key.profName))
				{
					aux.Key.active=true;	
					ActivateParents(aux.Value);
					foreach (var item in aux.Key.childrens.Values)
					{		
						
						auxP=new KeyValuePair<profileNode, string>(item,aux.Value+"/"+item.profName);
						//auxP.Key=item;
						//auxP.Value=aux.Value+"/"+item.profName;
            			nodeList.Enqueue(auxP);
						
					}
				}
			}
		}
        public void RemoveNode(string Path)
        {
            profileNode aux;
            Dictionary<string, profileNode> callAux;
            string parentPath;
            Regex exp = new Regex("/$");
            string[] tmp = Path.Split('/');

            parentPath = String.Join("/", tmp, 0, tmp.Length - 1);
            if (parentPath.Length == 0)
                parentPath = "/";
            aux = FindNode(parentPath);
            if (aux == null)
                callAux = masterNode;
            else
                callAux=aux.childrens;
            if(callAux.ContainsKey(tmp[tmp.Length - 1]))
                    callAux.Remove(tmp[tmp.Length - 1]);
            


        }
        public int GetProfileLength()
        {

            if (profiles == null || profiles.Count == 0)
                return 0;


            List<string> profName = new List<string>(profiles.Keys);
            List<string> protName = new List<string>(profiles[profName[0]].Keys);
            if (protName == null || protName.Count == 0)
                return 0;

            protInfo info = profiles[profName[0]][protName[0]];

            return info.alignment.Count;
        }
        public Dictionary<byte, double> StatesStat()
        {
            Dictionary<byte, double> statesStat = new Dictionary<byte, double>();
            int sum = 0;
            foreach (var item in protCombineStates)
            {
                foreach (var it in item.Value)
                {
                    if (statesStat.ContainsKey(it))
                        statesStat[it]++;
                    else
                        statesStat.Add(it, 0);

                    sum++;
                }
            }
            List<byte> keys = new List<byte>(statesStat.Keys);
            foreach (var item in keys)
                statesStat[item] /= sum;

            return statesStat;
        }
        public void AddItemsCombineStates(string key,List<string> states)
        {
            if (states == null || states.Count == 0)
                return;
            if (!protCombineStates.Keys.Contains(key))
                protCombineStates.Add(key, new List<byte>(states.Count));
            else
                return;

            foreach(var item in states)
            {
                if (item == "0")
                {
                    protCombineStates[key].Add(0);
                    if (!protCombineCoding.ContainsKey(item))                    
                        protCombineCoding.Add(item, 0);

                    continue;

                }
                if (!protCombineCoding.ContainsKey(item))
                {
                    byte num = 0;
                    if(!Byte.TryParse(item,out num))
                    {
                        for (byte i = 0; i < byte.MaxValue; i++)
                        {
                            bool test = false;
                            foreach (var itemV in protCombineCoding.Values)
                                if (itemV == i)
                                {
                                    test = true;
                                    break;
                                }
                            if(!test)
                            {
                                num = i;
                                break;
                            }
                        }
                                
                    }                    
                    protCombineCoding.Add(item, num);
                }
                protCombineStates[key].Add(protCombineCoding[item]);

            }
        }
        public Dictionary<string, List<byte>> CombineProfiles()
        {
            List<string> keys = new List<string>(profiles.Keys);
            protCombineStates = new Dictionary<string, List<byte>>();
            protCombineCoding = new Dictionary<string, byte>();
            foreach(var item in profiles[keys[0]])
            {
                protCombineStates.Add(item.Key, item.Value.profile);                
                foreach(var k in item.Value.profile)
                {
                    string v = k.ToString();
                    if (!protCombineCoding.ContainsKey(v))
                        protCombineCoding.Add(v, k);
                }
            }            
            return protCombineStates;
        }

        public void GetProfiles(string currentPath,Dictionary<string, profileNode> node ,Dictionary<string, profileNode> aux,bool leaves)
        {            

            foreach (var item in node)
            {
                string localPath = currentPath+item.Key+"/";
                if (leaves)
                {
                    if (item.Value.childrens.Count == 0)
                        aux.Add(localPath, item.Value);
                }
                else
                    aux.Add(localPath, item.Value);

                GetProfiles(localPath, item.Value.childrens, aux,leaves);               
            }
            
        }        
        public profileNode FindNode(string Path)
        {
            Dictionary<string, profileNode> aux = null;
            profileNode correct=null;
            int i=0;

            if (Path == null)
                return null;

            string[] tmp = Path.Split('/');
            if (masterNode.Count>0)
            {
                aux = masterNode;
                for(i=1;i<tmp.Length;i++)
                {
                    correct=null;

                    if (aux.Count == 0 || !aux.ContainsKey(tmp[i]))
                        return null;

                    correct = aux[tmp[i]];
                    aux = aux[tmp[i]].childrens;
                }
            }
            return correct;
        }
		public Dictionary<string,protInfo> ReadProfile(profileNode node)
		{
			string line;
            bool flag = false;
            protInfo profile;
            Dictionary<string,protInfo> profList = new Dictionary<string,protInfo>();
			string protName;			
			if(File.Exists(node.OutFileName))
			{
				
				StreamReader file=null;
				try
				{
					file=new StreamReader (node.OutFileName);
					line=file.ReadLine();
					while(line!=null)
					{
						if(line.Contains(">"))
						{
                            string[] profState;
                            profile = new protInfo();
                            flag = true;
							string []nameTmp;
							protName=line.Replace(">","");
                            
							nameTmp=protName.Split(new char [] {'/','\\'});
							protName=nameTmp[nameTmp.Length-1];
                            profile.sequence=file.ReadLine();
                            line = file.ReadLine();
                            profState = line.Split(' ');

                            profile.profile = new List<byte>(profState.Length);
                            foreach (var item in profState)
                                profile.profile.Add(node.GetCodedState(item));
							if(profile.sequence.Length>5 || profile.profile.Count!=0)
                            	profList.Add(protName, profile);
                            line=file.ReadLine();	
						}
                        else
						    line=file.ReadLine();	
					}
				}
				finally
				{
					if(file!=null)
						file.Close();
						
				}
                if (!flag)
                    throw new Exception("Error: In the profile file: " + node.OutFileName+" could not find line that begins with > ");
			}
            
            return profList;
		}

        public void MakeProfiles(Options opt)
        {
            int active = 0;
            List<profileNode> prof = GetActiveProfiles();
            if (prof == null)
                return;
            profilesNum = GetActiveProfiles().Count;
            profilesProgress=new List<double>(profilesNum);            

            active = GenerateActiveProfiles(opt,masterNode);
            foreach(var item in profiles)
            {
                foreach(var it in item.Value)
                    if (it.Value.profile.Count == 0)
                        throw new Exception("Failed to generate profiles!");
            }
            if (active == 0)
                throw new Exception("Error: there are no active profiles!");

        }
        //public void PrepareProfil
        public void PrepareProfiles(Options opt,string fullPath)
        {
            string[] files;
            
            profiles = new Dictionary<string, Dictionary<string, protInfo>>();

            RemoveOutFiles();

            if(set.extension==null)
                files = Directory.GetFiles(fullPath);
            else
                files = Directory.GetFiles(fullPath, set.extension);
            

            if (files.Length == 0 && !EmptyProfiles())
                throw new Exception("In selected directory " + fullPath + " there are no files with extesion " + set.extension);

            Random r = new Random();
            int val = r.Next(100000);

            listFile = "list_" + val + ".txt";

            StreamWriter outputFile = new StreamWriter(listFile);
            foreach (string i in files)
                outputFile.WriteLine(i);

            outputFile.Close();

            if (masterNode == null || masterNode.Values == null || masterNode.Values.Count == 0)
                throw new Exception("Error: Profiles have not been defined");

            MakeProfiles(opt);

            if (File.Exists(listFile))
                File.Delete(listFile);


        }
        public void PrepareProfiles(Options opt,List<string> fileNames)
        {
            profiles = new Dictionary<string, Dictionary<string, protInfo>>();
            RemoveOutFiles();

            Random r = new Random();
            int val = r.Next(100000);

            listFile = "list_" + val + ".txt";

            StreamWriter outputFile = new StreamWriter(listFile);
            foreach (string i in fileNames)
                outputFile.WriteLine(i);
            outputFile.Close();

            if (masterNode == null || masterNode.Values == null || masterNode.Values.Count == 0)
                throw new Exception("Error: Profiles have not been defined");


            MakeProfiles(opt);

            if (File.Exists(listFile))
                File.Delete(listFile);


        }
        public bool EmptyProfiles()
        {
            List<profileNode> active = GetActiveProfiles();

            if(active!=null)
            foreach (var item in active)
                if (item.profProgram.Length > 0)
                    return false;

            return true;
        }
        public void GenerateProfile(profileNode item)
        {
            Process myProc = new Process();

            myProc.StartInfo.FileName = item.profProgram;
            myProc.StartInfo.RedirectStandardOutput = false;
            myProc.StartInfo.RedirectStandardError = false;
            myProc.StartInfo.UseShellExecute = false;
            //myProc.StartInfo.CreateNoWindow = true ;
/*            Environment.SetEnvironmentVariable("INSTALLDIR", set.installDir);
            Environment.SetEnvironmentVariable("INPUTFILE", listFile);
            Environment.SetEnvironmentVariable("OUTPUTFILE", item.OutFileName);*/
            string param = item.progParameters;
            myProc.StartInfo.Arguments = param + " " + item.OutFileName;
           

            myProc.Start();
            myProc.WaitForExit();           

        }
        public List<profileNode> GetActiveLeaves()
        {
            return ActiveProfiles(true);
        }
        public List<profileNode> GetActiveProfiles()
        {
            return ActiveProfiles(false);
        }
        public string GetStringActiveProfiles()
        {
            string outStr = "";
            List<profileNode> active=ActiveProfiles(false);
            if (active != null)
                for (int i = 0; i < active.Count; i++)
                {
                    outStr += active[i].profName;
                    if (i < active.Count - 1)
                        outStr += ", ";
                }
            return outStr;
        }
        private List<profileNode> ActiveProfiles(bool leaves)
        {
            List<profileNode> active = new List<profileNode>();

            Queue<profileNode> nodeList = new Queue<profileNode>();
            profileNode aux;
            if (masterNode.Count == 0)
                return null;
            foreach (var item in masterNode.Values)
                if (item.active)
                {
                    if (leaves)
                    {
                        if (item.childrens.Count == 0)
                            active.Add(item);
                    }
                    else
                        active.Add(item);
                    nodeList.Enqueue(item);
                }

            while (nodeList.Count > 0)
            {
                aux = nodeList.Dequeue();

                foreach (var item in aux.childrens.Values)
                    if (item.active)
                    {
                        if (leaves)
                        {
                            if (item.childrens.Count == 0)
                                active.Add(item);
                        }
                        else
                            active.Add(item);
                        nodeList.Enqueue(item);
                    }
            }
            


            return active;
        }
        public int GenerateActiveProfiles(Options opt,Dictionary<string, profileNode> node)
        {
            int activeProfiles = 0;
            
            if (node == null)
                node = masterNode;

            if (profiles == null)
                profiles = new Dictionary<string, Dictionary<string, protInfo>>();

            DebugClass.WriteMessage("Profiling started");

            foreach (var item in node.Values)
            {
                if (item.active)
                {
                    Dictionary<string, protInfo> profString;
                    mgr.ResetProgress();
                    if (item.GetNumberofStates()== 0)
                    {
                        if (item.profProgram.Length > 0)
                        {
                            if (!profiles.ContainsKey(item.profName))
                            {
                                activeProfiles++;
                                if (File.Exists(item.OutFileName))
                                    File.Delete(item.OutFileName);
                                GenerateProfile(item);
                                FileInfo f = new FileInfo(item.OutFileName);
                                if (!File.Exists(item.OutFileName) || f == null || f.Length == 0)
                                    throw new Exception("Error: profile " + item.profName + " did not generate output " + item.OutFileName);
                            }
                        }
                    }
                    else
                    {
                        DebugClass.WriteMessage("Run profiling");                        
                            mgr.RunProfile(opt,item.internalName, listFile);
                        activeProfiles++;
                    }

                    if (item.childrens.Count > 0)
                        activeProfiles += GenerateActiveProfiles(opt,item.childrens);
                    else
                    {
                        if (!profiles.ContainsKey(item.profName))
                        {
                            if (item.GetNumberofStates()== 0)
                                profString = ReadProfile(item);
                            else
                            {
                                    profString = mgr.GetProfile(opt,item, listFile);
                            }
                            if (profString!=null && profString.Count > 0)
                                profiles.Add(item.profName, profString);
                            else
                                throw new Exception("For profile " + item.profName + " profiles has not been generated!");

                        }

                    }
                    profilesProgress.Add(mgr.ProgressUpdate());

                }


            }
            progressObject = null;
            return activeProfiles;
        }
        private void SetupCodingAllNodes(Dictionary<string, profileNode> node)
        {
            foreach (var item in node.Values)
            {
                if (item.active)
                {
                    item.SetupCoding();
                    SetupCodingAllNodes(item.childrens);
                }

            }
//
        }
        public Stack<string> GetAllParentsPaths(string Path)
        {
            string subPath;
            string[] tmp;
            Regex expr = new Regex("/$");
            Stack<string> paths = new Stack<string>();

            Path = expr.Replace(Path, "");
            tmp = Path.Split('/');
            if (tmp.Length > 1)
            {
                for (int i = tmp.Length - 1; i > 1; i--)
                {                    
                    subPath = String.Join("/", tmp, 0, i);
                    paths.Push(subPath);
                }
            }
            else
                paths.Push(Path);

            return paths;

        }
        void LoadDefaultProfile()
        {
            if (File.Exists(defaultProfile))
            {
                StreamReader stream = new StreamReader(defaultProfile);
                string profFile = stream.ReadLine();
                LoadProfiles(profFile);
            }
        }
        void SaveDefaultProfile(string profileName)
        {
            StreamWriter stream = new StreamWriter(defaultProfile);
            stream.WriteLine(profileName);
            stream.Close();                            
        }
        public void SaveProfiles(string fileName)
        {
            SaveTree(fileName);          
        }
        public void LoadProfiles(string profileName)
        {
            if (File.Exists(profileName))
            {
                LoadTree(profileName);
                if (masterNode == null)
                    DebugClass.WriteMessage("Ups no master node");
                currentProfileName = profileName;
                //SetupCodingAllNodes(masterNode);
            }
            else
                throw new Exception("Profile: " + profileName + " not exists. Cannot be loaded!");
        }
        public void SaveTree(string fileName)
        {
                Stream stream = File.Open(fileName, FileMode.Create);
                
                //BinaryFormatter bFormatter = new BinaryFormatter();
            //    XmlSerializer ser = new XmlSerializer(typeof(Dictionary<string, profileNode>));
                XmlSerializer ser =new XmlSerializer(typeof(SerializableDictionary<string, profileNode>));
                //bFormatter.Serialize(stream, masterNode);
                ser.Serialize(stream, masterNode);
                stream.Close();
        }
       
     
        public void LoadTree(string fileName)
        {
            if (!File.Exists(fileName))
                return;
            Stream stream = File.Open(fileName, FileMode.Open);
            //BinaryFormatter bFormatter = new BinaryFormatter();
            XmlSerializer ser = new XmlSerializer(typeof(SerializableDictionary<string, profileNode>));
            //SerializableDictionary<string, profileNode> localProfiles = (SerializableDictionary<string, profileNode>)bFormatter.Deserialize(stream);
            SerializableDictionary<string, profileNode> localProfiles = (SerializableDictionary<string, profileNode>)ser.Deserialize(stream);
            masterNode = localProfiles;
            stream.Close();
            
        }

    }
    [Serializable]
    [XmlRoot("dictionary")]
    public class SerializableDictionary<TKey, TValue>
        : Dictionary<TKey, TValue>, IXmlSerializable
    {
        public SerializableDictionary()
        : base()
        {
        }
        public SerializableDictionary(SerializationInfo info, StreamingContext context)
  : base(info, context)
        {
        }
        #region IXmlSerializable Members
        public System.Xml.Schema.XmlSchema GetSchema()
        {
            return null;
        }

        public void ReadXml(System.Xml.XmlReader reader)
        {
            XmlSerializer keySerializer = new XmlSerializer(typeof(TKey));
            XmlSerializer valueSerializer = new XmlSerializer(typeof(TValue));

            bool wasEmpty = reader.IsEmptyElement;
            reader.Read();

            if (wasEmpty)
                return;

            while (reader.NodeType != System.Xml.XmlNodeType.EndElement)
            {
                reader.ReadStartElement("item");

                reader.ReadStartElement("key");
                TKey key = (TKey)keySerializer.Deserialize(reader);
                reader.ReadEndElement();

                reader.ReadStartElement("value");
                TValue value = (TValue)valueSerializer.Deserialize(reader);
                reader.ReadEndElement();

                this.Add(key, value);

                reader.ReadEndElement();
                reader.MoveToContent();
            }
            reader.ReadEndElement();
        }

        public void WriteXml(System.Xml.XmlWriter writer)
        {
            XmlSerializer keySerializer = new XmlSerializer(typeof(TKey));
            XmlSerializer valueSerializer = new XmlSerializer(typeof(TValue));

            foreach (TKey key in this.Keys)
            {
                writer.WriteStartElement("item");

                writer.WriteStartElement("key");
                keySerializer.Serialize(writer, key);
                writer.WriteEndElement();

                        
                writer.WriteStartElement("value");
                TValue value = this[key];
                valueSerializer.Serialize(writer, value);
                writer.WriteEndElement();

                writer.WriteEndElement();
            }
        }
        #endregion
    }
    [Serializable]
    public struct RangeWeight
    {
        public int start;
        public int stop;
        public double Weight;
    }
    [Serializable]
    public struct protInfo
    {
        public string sequence;
        public List<byte> profile;
        public List<byte> alignment;

    }
    [Serializable]
    public class profileNode
    {
        public string internalName { get; set; }
        public string profName { get; set; }
        public string profProgram { get; set; }
        public string progParameters { get; set; }
        public bool removeOutFile { get; set; }
        public string OutFileName { get; set; }
        public bool active = true;
        public bool deleteAble = true;
        public List<RangeWeight> rangeWeights = null;
        public SerializableDictionary<string, SerializableDictionary<string, double>> profWeights=new SerializableDictionary<string,SerializableDictionary<string,double>>();
        public SerializableDictionary<string, string> states = new SerializableDictionary<string, string>();        
        public SerializableDictionary<string, byte> codingToByte = new SerializableDictionary<string, byte>();        
        public SerializableDictionary<byte,string> codingToString = new SerializableDictionary<byte,string>();
        [NonSerialized]
        bool[] index = new bool[byte.MaxValue+1];
        public profileNode()
        {
        }
        public profileNode(string internalName)
        {
            this.internalName = internalName;
        }       

        public profileNode(profileNode copy)
        {
            CopyNode(copy);
        }
        public int GetNumberofStates()
        {
            return states.Keys.Count;
        }
        public bool ContainsState(string stateItem)
        {
            if (states.ContainsKey(stateItem))
                return true;

            return false;
        }
        public byte GetCodedState(string stateItem,bool force=false)
        {
            if (!codingToByte.ContainsKey(stateItem))
                 AddStateItem(stateItem,stateItem,force);


            return codingToByte[stateItem];                            
        }
        public void SetupCoding()
        {
            foreach(var item in states.Values)
            {
                if (!codingToByte.ContainsKey(item))
                {
                    byte counter = (byte)(codingToByte.Keys.Count + 1);

                    codingToByte.Add(item, counter);
                    codingToString.Add(counter, item);

                }
            }
        }
        public void AddStateItem(string itemKey,string itemValue, bool force=false )
        {
            byte stateValue=0;
            if (!states.ContainsKey(itemKey))                
                states.Add(itemKey,itemValue);        

            if (codingToByte.ContainsKey(itemValue))
                return;
            int num;
            if (!Int32.TryParse(itemValue, out num))
            {
                for (int i = 1; i < index.Length; i++)
                {
                    if (!index[i])
                    {
                        stateValue = (byte)i;
                        break;
                    }
                }
            }
            else
                if (num > 255)
                    stateValue = 255;
                else
                    stateValue = (byte)num;

            //if(!index[stateValue])
                
            index[stateValue] = true;
            codingToByte.Add(itemValue, stateValue);
            if(!codingToString.ContainsKey(stateValue))
                codingToString.Add(stateValue, itemValue);
        }
        public void CopyNode(profileNode copy)
        {
            profName = copy.profName;
            profProgram = copy.profProgram;
            internalName = copy.internalName;
            progParameters = copy.progParameters;
            removeOutFile = copy.removeOutFile;
            OutFileName = copy.OutFileName;
            active = copy.active;
            profWeights.Clear();
            foreach (var item in copy.profWeights.Keys)
            {
                SerializableDictionary<string, double> aux = new SerializableDictionary<string, double>();
                foreach (var item2 in copy.profWeights[item].Keys)
                    aux.Add(item2, copy.profWeights[item][item2]);

                profWeights.Add(item, aux);
            }
            states.Clear();
            foreach (var item in copy.states.Keys)
                states.Add(item, copy.states[item]);

            childrens.Clear();
            foreach (var item in copy.childrens)
                childrens.Add(item.Key,item.Value);

        }

        public SerializableDictionary<string, profileNode> childrens = new SerializableDictionary<string, profileNode>();

    }

}
    