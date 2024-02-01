using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.IO;

namespace phiClustCore
{
    class RootHash:HNN
    {
        protected int binSizeG = 0;
        protected ManualResetEvent[] resetEvents;
        protected Dictionary<string, List<int>>[] hashDataBase = null;
        protected Dictionary<string, int[]>[] hashDataBaseTab = null;
        protected List<string> dataBaseKeys = null;
        protected Dictionary<string, string> results = new Dictionary<string, string>();
        protected int[][] dist;
        public int resSize = 20;
        //public StreamWriter wr ;
        public RootHash(int binSizeG, HashCluster hk, ClusterOutput outp, HNNCInput opt):base(hk, outp, opt)
        {
            this.binSizeG = binSizeG;
            //wr = new StreamWriter("log.txt");
            // number of bits to store the universe
        }
        protected List<string> GetKeyHashes(string item)
        {
            List<string> aux = new List<string>();
            for (int i = 0, counter = 0; i < item.Length; i += binSizeG, counter++)
            {
                string keyHash = "";
                if (i + binSizeG > item.Length)
                {
                    keyHash = item.Substring(i, item.Length - i);
                }
                else
                    keyHash = item.Substring(i, binSizeG);
                aux.Add(keyHash);
            }

            return aux;
        }

        protected string [] GetKeyHashesTab(string item)
        {
            string[] aux = new string [item.Length/binSizeG];
            for (int i = 0, counter = 0; i < item.Length; i += binSizeG, counter++)
            {
                string keyHash = "";
                if (i + binSizeG > item.Length)
                {
                    keyHash = item.Substring(i, item.Length - i);
                }
                else
                    keyHash = item.Substring(i, binSizeG);
                aux[i]=keyHash;
            }

            return aux;
        }



        public virtual void CreateBase(Dictionary<string, string> dataBase)
        {

            dataBaseKeys = new List<string>(dataBase.Keys);

            double x = dataBase[dataBaseKeys[0]].Length;
            x /= binSizeG;
            int tabSize = (int)x;

            if (x != Math.Round(x))
                tabSize++;

            hashDataBase = new Dictionary<string, List<int>>[tabSize];
            for (int i = 0; i < hashDataBase.Length; i++)
                hashDataBase[i] = new Dictionary<string, List<int>>();


            for (int j = 0; j < dataBaseKeys.Count; j++)
            {                
                List<string> aux = GetKeyHashes(dataBase[dataBaseKeys[j]]);
                for (int i = 0; i < aux.Count; i++)
                {
                    Dictionary<string, List<int>> locHash = hashDataBase[i];
                    if (locHash.ContainsKey(aux[i]))
                        locHash[aux[i]].Add(j);
                    else
                    {
                        List<int> xaux = new List<int>();
                        xaux.Add(j);
                        locHash.Add(aux[i], xaux);
                    }
                }
            }

            hashDataBaseTab = new Dictionary<string, int[]>[hashDataBase.Length];
            for(int i=0;i<hashDataBase.Length;i++)
            {
                hashDataBaseTab[i] = new Dictionary<string, int[]>();
                foreach(var item in hashDataBase[i])
                {
                    hashDataBaseTab[i].Add(item.Key, item.Value.ToArray());
                }
            }
        }
        public override void Preprocessing()
        {
            caseBase = new Dictionary<string, string>();
            foreach (var item in hk.dicFinal)
                for (int i = 0; i < item.Value.Count; i++)
                    caseBase.Add(hk.structNames[item.Value[i]], item.Key);
                    //caseBase.Add(hk.structNames[item.Value[0]], item.Key);




            dist = new int[set.numberOfCores][];

            CreateBase(caseBase);
            testData = caseBase;

            for (int i = 0; i < dist.Length; i++)
                dist[i] = new int[caseBase.Count];



        }
        public virtual string CalcDist(int threadNum, List<string> keys, int[] index, int num)
        {
            int[] locDist = dist[threadNum];
            
            for (int i = 0; i < keys.Count; i++)
            {
                int[] baseHash = null;
                if (hashDataBaseTab[i].TryGetValue(keys[i], out baseHash))
                //if (auxBase[i].TryGetValue(keys[i], out baseHash))
                {
                    for(int j=0;j<baseHash.Length;j++)
                    //foreach (var inx in baseHash)
                        locDist[baseHash[j]]++;
                }

            }
            //Array.Sort(locDist, index);
            //Array.Reverse(locDist);
            //Array.Reverse(index);
            string w = "";
            for (int i = 0; i < resSize; i++)
                w += dataBaseKeys[index[i]] + ":";
            w += dataBaseKeys[index[resSize]];

            return w;

        }
        public virtual void CalcHashDist(object o)
        {
            ThreadParam pp = (ThreadParam)o;
            int threadNum = pp.num;
            int start = pp.start;
            int stop = pp.stop;
            int[] index = new int[dist[threadNum].Length];
            int[] locDist = dist[threadNum];
            for (int n = start; n < stop; n++)
            {
                string testItem = testList[n];

                if (!testData.ContainsKey(testItem))
                    continue;

                List<string> keys = GetKeyHashes(testData[testItem]);
                for (int i = 0; i < locDist.Length; i++)
                {
                    locDist[i] = 0;
                    index[i] = i;
                }

                int num = 0;
                /*for (int i = 0; i < hashDataBase.Length; i++)
                {
                    if (hashDataBase[i].ContainsKey(testData[testItem][i].ToString()))
                        num++;
                }*/

                string w=CalcDist(threadNum, keys, index, num);
                lock (results)
                {
                    results.Add(testItem, w);
                }


            }
            resetEvents[threadNum].Set();
        }
        void Retrival(List<string> testList)
        {

            this.testList = testList;

            resetEvents = new ManualResetEvent[dist.Length];
            for (int n = 0; n < resetEvents.Length; n++)
            {
                ThreadParam pp = new ThreadParam();
                pp.num = n;
                pp.start = (int)(n * testList.Count / resetEvents.Length);
                pp.stop = (int)((n + 1) * testList.Count / resetEvents.Length);
                resetEvents[n] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(new WaitCallback(CalcHashDist), (object)pp);
            }

            for (int n = 0; n < resetEvents.Length; n++)
                resetEvents[n].WaitOne();
            //wr.Close();
        }
        public override Dictionary<string, string> HNNTest(List<string> testList)
        {           
            Retrival(testList);
            return results;

        }

    }
}
