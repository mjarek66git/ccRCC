using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading;
using System.Reflection;
using System.Xml.Serialization;
using System.Threading.Tasks;
using phiClustCore.Interface;

namespace phiClustCore.Profiles
{
    class ThreadFiles
    {
        public string fileName;
        public string dirName;
        public int threadNumber;
        //public DistanceMeasure distance;
    }
    [Serializable]
    public abstract class InternalProfileBase:IProgressBar
    {
        protected Exception ex = null;
        protected int threadNumbers;
        protected static string SEQprofile = "SEQ profile ";        
        public List<INPUTMODE> destination=null;
        Settings set = new Settings();
        abstract public int Run(object processParams);
        
//        virtual public void RunThreads(string listFile);
        abstract public void AddInternalProfiles();

        abstract public void RemoveInternalProfiles();

        abstract public List<INPUTMODE> GetMode();
        abstract public Dictionary<string, protInfo> GetProfile(profileNode node, string listFile);
        protected string profileName = "";
        protected List<string>[] threadingList;

        //protected InternalProfilesManager manager = new InternalProfilesManager();

        protected int maxV = 1;
        protected int currentProgress = 0;
        protected Object thisLock = new Object();

        double startProgress = 0;
        double endProgress = 1;

        public double StartProgress { set { startProgress = value; } get { return startProgress; } }
        public double EndProgress { set { endProgress = value; } get { return endProgress; } }
        
        public InternalProfileBase()
        {
            set.Load();
           // threadNumbers = 1;//22 set.numberOfCores;
            threadNumbers = set.numberOfCores;
        }
        public void LoadProfiles()
        {
            InternalProfilesManager.LoadProfiles();
        }
        public double ProgressUpdate()
        {
            return StartProgress + (EndProgress - StartProgress) * (double)currentProgress / maxV;
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }
        public virtual void CheckIfAvailable()
        {
        }

        private string GetDirectory(string fileName)
        {
            string dirName = "";
            if (!File.Exists(fileName))
                return null;

            using (StreamReader rd = new StreamReader(fileName))
            {
                string line = rd.ReadLine();
                DebugClass.WriteMessage("line=" + line);
                rd.Close();
             
                try
                {
                    if (line.Contains(Path.DirectorySeparatorChar))
                    {
                        string nn = line;
                        if(line.Contains("|"))
                        {
                            string[] aux = line.Split('|');
                            nn = aux[0];
                                
                        }
                        dirName = Path.GetDirectoryName(nn);
                    }
                }
                catch (Exception )
                {
                    return null;
                }
            }

            return dirName;
        }
        public void ClearProfiles()
        {
            InternalProfilesManager.ClearProfiles();
        }
        protected List<string> GetFileList(string fileName)
        {
            List<string> fileList = new List<string>();
            if (!File.Exists(fileName))
                return null;

            using (StreamReader rd = new StreamReader(fileName))
            {
                string line = rd.ReadLine();

                while (line != null)
                {
                    string[] aux = line.Split('|');
                    string file = line;
                    if (aux.Length == 2)
                        file = aux[0];
                    if (File.Exists(file))
                        fileList.Add(line);
                    line = rd.ReadLine();
                }

                rd.Close();
            }
            return fileList;
        }

        protected string GetProfileFileName(string fileName)
        {
            string[] aux;

            string dirName = GetDirectory(fileName);

            DebugClass.WriteMessage(" name=" + dirName);
            if (dirName != null)
            {
                string bDir = "";
                aux = dirName.Split(Path.DirectorySeparatorChar);
                if (set.profilesDir != null && set.profilesDir.Length > 0)
                    bDir = set.profilesDir + Path.DirectorySeparatorChar;

                DebugClass.WriteMessage("bdir=" + bDir);
                return bDir + aux[aux.Length - 1] + "_" + profileName;
            }
            return null;
        }

        protected List<string> CheckFile(string fileName)
        {
            List<string> fileToCheck = new List<string>();
            List<string> checkList = new List<string>();
            if (File.Exists(GetProfileFileName(fileName)))
            {
                List<string> files = GetFileList(fileName);
                Dictionary<string, string> checking = new Dictionary<string, string>();
                string[] aux;
                foreach (var item in files)
                {
                    string newItem = item;
                    if (item.Contains("|"))
                    {
                        aux = item.Split('|');
                        if (aux.Length == 2)
                        {
                            newItem = aux[0];
                        }
                    }
                    aux = newItem.Split(Path.DirectorySeparatorChar);
                    if (!checking.ContainsKey(aux[aux.Length - 1]))
                        checking.Add(aux[aux.Length-1], item);
                }

                using (StreamReader rf = new StreamReader(GetProfileFileName(fileName)))
                {
                    if (rf == null)
                        return null;

                    string line = rf.ReadLine();
                    while (line != null)
                    {
                        if (line.Contains('>'))
                        {
                            line = line.Remove(0, 1);
                            if (line.Contains("|"))
                            {
                                string[] tmp = line.Split('|');
                                line = tmp[0];
                            }

                            //string tmp = GetDirectory(fileName);
                            //string name = tmp + Path.DirectorySeparatorChar + line;
                            if (checking.ContainsKey(line))
                                checkList.Add(line);
                        }
                        line = rf.ReadLine();
                    }
                    rf.Close();
                }
                foreach (var item in checking)
                    if (!checkList.Contains(item.Key))
                        fileToCheck.Add(item.Value);

            }
            else
                fileToCheck = GetFileList(fileName);


            return fileToCheck;
        }
        private void PrepareDataForThreading(List<string> files,int threadNumbers)
        {
            threadingList = new List<string>[threadNumbers];

            for(int i=0;i<threadNumbers;i++)
                threadingList[i]=new List<string>();

            for(int i=0;i<threadNumbers;i++)
            {
                for (int j = i * files.Count / threadNumbers; j < (i + 1) * files.Count / threadNumbers; j++)
                    threadingList[i].Add(files[j]);

            }
        }
        public virtual void RunThreads(string fileName)
        {           
            List<string> files = CheckFile(fileName);

            if (files.Count == 0)
                return;
            maxV = files.Count;
            //Wrapper wrapper = new Wrapper();
            //List<Thread> runnigThreads = new List<Thread>();
            if (set.iOTroubles)// || this is RNASSProfile)
                threadNumbers = 1;
            Task[] runnigTask = new Task[threadNumbers];
            Task startProg;

            PrepareDataForThreading(files, threadNumbers);

            for (int i = 0; i < threadNumbers; i++)
            {
                ThreadFiles tparam = new ThreadFiles();


                tparam.fileName = GetProfileFileName(fileName) + "_" + i;
                tparam.dirName = GetDirectory(fileName);
                tparam.threadNumber = i;
                startProg = new Task<int>(n => Run((object) n),tparam);

                startProg.Start();
                //startProg.Start(tparam);
                runnigTask[i]=startProg;

                //while (!startProg.IsAlive) ;
                //RunSubList(auxFiles, fileN);
            }
            for (int i = 0; i < threadNumbers; i++)
                runnigTask[i].Wait();
                //Task.WaitAll(runnigTask);

            currentProgress = maxV;
            if (ex != null)
                    throw ex;
                //runnigTask[i].Join();
            if(files.Count>0)
                JoinFiles(fileName);

        }

        public virtual void JoinFiles(string fileName)
        {
            StreamWriter wr;
            
                if (File.Exists(GetProfileFileName(fileName)))
                    wr = File.AppendText(GetProfileFileName(fileName));
                else
                    wr = new StreamWriter(GetProfileFileName(fileName));

                if (wr == null)
                    throw new Exception("Cannot open file: " + GetProfileFileName(fileName));

                for (int i = 0; i < threadNumbers; i++)
                {
                    string fileN = GetProfileFileName(fileName) + "_" + i;
                    using (StreamReader rr = new StreamReader(fileN))
                    {
                        if (rr == null)
                            throw new Exception("Cannot open file: " + fileN);

                        string line = rr.ReadLine();
                        while (line != null)
                        {
                            wr.WriteLine(line);
                            line = rr.ReadLine();
                        }
                        rr.Close();
                        File.Delete(fileN);
                    }
                }
                wr.Close();
            

        }

    }
    [Serializable]
    public class InternalProfilesManager: IProgressBar
    {
        static string fileName = "internal.profiles";
        string dirName;
        public static SerializableDictionary<profileNode, string> internalList = new SerializableDictionary<profileNode, string>();
        IProgressBar progressInfo = null;
        double progress = 0;
        double startProgress = 0;
        double endProgress = 1;

        public double StartProgress { set { startProgress = value; } get { return startProgress; } }
        public double EndProgress { set { endProgress = value; } get { return endProgress; } }

        public static List<InternalProfileBase> InitProfiles()
        {
            //Type[] types = GetLoadableTypes(Assembly.GetExecutingAssembly());
            //Type[] types = Assembly.GetExecutingAssembly().GetTypes();
            Type [] types;
            List<InternalProfileBase> myTypes = new List<InternalProfileBase>();
            if (Assembly.GetExecutingAssembly() == null) throw new ArgumentNullException("Assembly cannot be readed");
            try
            {
                types=Assembly.GetExecutingAssembly().GetTypes();
            }
            catch (ReflectionTypeLoadException e)
            {
                types=e.Types;
            }
          

            RemoveProfilesFile();
            foreach (Type t in  types)
            {
                if (t!=null && t.IsSubclassOf(typeof(InternalProfileBase)))
                {
                    try
                    {
                        InternalProfileBase c = Activator.CreateInstance(t) as InternalProfileBase;
                        c.CheckIfAvailable();
                        myTypes.Add(c);
                    }
                    catch(Exception ex)
                    {
                        Console.WriteLine(ex);
                    }
                }
            }
            return myTypes;
        }

        public double ProgressUpdate()
        {
            if (progressInfo == null)
                return StartProgress + (EndProgress - StartProgress) * progress;

            double res = progress + progressInfo.ProgressUpdate();
            return StartProgress + (EndProgress - StartProgress) * res;
        }
        public void ResetProgress()
        {
            progress = 0;
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }



        public InternalProfilesManager(string dirName=null)
        {
            this.dirName = dirName;
            LoadProfiles();            
            
        }
        public static void LoadProfiles()
        {
            if (!File.Exists(fileName))
                return;
            Stream stream = File.Open(fileName, FileMode.Open);
            XmlSerializer ser = new XmlSerializer(typeof(SerializableDictionary<profileNode,string>));
            SerializableDictionary<profileNode,string> localProfiles = (SerializableDictionary<profileNode,string>)ser.Deserialize(stream);
            internalList = localProfiles;
            stream.Close();

        }
        static void SaveProfiles()
        {
            Stream stream = File.Open(fileName, FileMode.Create);
            XmlSerializer ser = new XmlSerializer(typeof(SerializableDictionary<profileNode,string>));
            ser.Serialize(stream, internalList);
            stream.Close();
        }
        public static void RemoveProfilesFile()
        {            
            if (File.Exists(fileName))
                File.Delete(fileName);
        }    
    

        public Dictionary<string, protInfo> GetProfile(Options opt,profileNode node,string listFile)
        {
            profileNode n=GetNode(node.internalName);
               Dictionary<string, protInfo> res=null;
               if (n != null)
               {
                   Type t = Type.GetType(internalList[n]);
                   InternalProfileBase c = Activator.CreateInstance(t,opt) as InternalProfileBase;
                   progressInfo = c;
                   res = c.GetProfile(node, listFile);
                   progressInfo = null;
                   progress++;
                   return res;
               }
            return null;
        }
        public static profileNode GetNode(string name)
        {
            foreach (var item in internalList.Keys)
                if (item.internalName == name)
                    return item;
           
            return null;
        }
        public static void RemoveNodeFromList(string profName)
        {
            List<profileNode> toRemove=new List<profileNode>();
            foreach(var item in internalList)
            {
                if (item.Key.profName == profName)
                    toRemove.Add(item.Key);
            }
            foreach (var item in toRemove)
                internalList.Remove(item);

            SaveProfiles();

        }
        public static void AddNodeToList(profileNode node,string type)
        {
            profileNode nn = GetNode(node.internalName);
           
            if (nn==null)
            {
                if (!internalList.ContainsKey(node))
                {
                    internalList.Add(node, type);
                    SaveProfiles();
                }
            }

        }
        public bool CheckAccessibility(profileNode profile,INPUTMODE type)
        {
            Type t = Type.GetType(internalList[profile]);
            InternalProfileBase c = Activator.CreateInstance(t) as InternalProfileBase;

            if (c.destination.Contains(type))
                return true;

            return false;
        }
        public static void ClearProfiles()
        {
            internalList.Clear();
        }
        public void RunProfile(Options opt,string name,string listFile)
        {
            profileNode node=GetNode(name);

            if (node != null)
            {
                Type t = Type.GetType(internalList[node]);

                InternalProfileBase c = Activator.CreateInstance(t, opt) as InternalProfileBase;
                c.CheckIfAvailable();
                progressInfo = c;
                c.RunThreads(listFile);
                progressInfo = null;
                progress++;
            }

        }
    }
}
