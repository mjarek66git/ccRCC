using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore.Distance
{
    public class Fast3States : HammingBase
    {
        protected int binSizeG = 0;
        [NonSerialized]
        protected Dictionary<double, List<int>>[] hashDataBase = null;
        [NonSerialized]
        protected Dictionary<double, int[]>[] hashDataBaseTab = null;
        public string[] dataBaseKeys = null;
        [NonSerialized]
        protected Dictionary<string, string> results = new Dictionary<string, string>();
        [NonSerialized]
        protected int[][] dist;
        protected int[] minStateDataBase = null;
        protected int[] maxStateDataBase = null;
        [NonSerialized]
        protected Dictionary<double, KeyValuePair<int, List<int>>>[] locBase = null;
        [NonSerialized]
        protected Dictionary<double, KeyValuePair<int, int[]>>[] locBaseTab = null;
        protected int[,][] oTab = null;
        List<int>[] indexList;

        public int resSize = 20;
        string baseFileName = "";

        public Fast3States(int binSizeG, Dictionary<string, List<double>> data) : base(data, true)
        {
            this.binSizeG = binSizeG;
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
        public virtual int[] GetKeyHashesTabInt(string item)
        {
            int[] aux = new int[(int)Math.Ceiling(((double)item.Length) / binSizeG)];

            for (int i = 0, counter = 0; i < item.Length; i += binSizeG, counter++)
            {
                aux[counter] = Convert.ToInt32(item[i] - '0');
            }

            return aux;
        }
        public override void InitMeasure()
        {
            base.InitMeasure();
            CreateBase();
        }
        protected virtual string[] GetKeyHashesTab(string item)
        {
            string[] aux = new string[(int)Math.Ceiling(((double)item.Length) / binSizeG)];
            for (int i = 0, counter = 0; i < item.Length; i += binSizeG, counter++)
            {
                string keyHash = "";
                if (i + binSizeG > item.Length)
                {
                    keyHash = item.Substring(i, item.Length - i);
                }
                else
                    keyHash = item.Substring(i, binSizeG);
                aux[counter] = keyHash;
            }

            return aux;
        }



        public void CreateBase()
        {

            dataBaseKeys = data.Keys.ToArray();

            double x = data[dataBaseKeys[0]].Count;
            x /= binSizeG;
            int tabSize = (int)x;

            if (x != Math.Round(x))
                tabSize++;

            hashDataBase = new Dictionary<double, List<int>>[tabSize];
            for (int i = 0; i < hashDataBase.Length; i++)
                hashDataBase[i] = new Dictionary<double, List<int>>();


            for (int j = 0; j < dataBaseKeys.Length; j++)
            {
                List<double> aux = data[dataBaseKeys[j]];
                for (int i = 0; i < aux.Count; i++)
                {
                    Dictionary<double, List<int>> locHash = hashDataBase[i];
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

            hashDataBaseTab = new Dictionary<double, int[]>[hashDataBase.Length];
            for (int i = 0; i < hashDataBase.Length; i++)
            {
                hashDataBaseTab[i] = new Dictionary<double, int[]>();
                foreach (var item in hashDataBase[i])
                {
                    hashDataBaseTab[i].Add(item.Key, item.Value.ToArray());
                }
            }
            minStateDataBase = new int[dataBaseKeys.Length];
            maxStateDataBase = new int[dataBaseKeys.Length];
            for (int j = 0; j < minStateDataBase.Length; j++)
            {
                minStateDataBase[j] = 0;
                maxStateDataBase[j] = 0;
            }

            for (int j = 0; j < hashDataBase.Length; j++)
            {
                double maxKey = double.MinValue;
                double minKey = double.MinValue;
                int max = int.MinValue;
                int min = int.MaxValue;
                Dictionary<double, List<int>> locHash = hashDataBase[j];
                if (hashDataBase[j].Keys.Count > 1)
                {
                    foreach (var item in hashDataBase[j].Keys)
                    {
                        if (locHash[item].Count > max)
                        {
                            max = locHash[item].Count;
                            maxKey = item;
                        }
                        else
                        if (locHash[item].Count < min)
                        {
                            min = locHash[item].Count;
                            minKey = item;
                        }

                    }

                   // foreach (var item in locHash[minKey])
                     //   minStateDataBase[item]++;

                    locHash.Remove(maxKey);

                    foreach(var item in locHash.Keys)
                        foreach (var it in locHash[item])
                            minStateDataBase[it]++;

                }
                else
                    locHash.Clear();
            }
            // binSizeG = 1;


            locBase = new Dictionary<double, KeyValuePair<int, List<int>>>[hashDataBase.Length / binSizeG + 1];

            for (int i = 0; i < locBase.Length; i++)
                locBase[i] = new Dictionary<double, KeyValuePair<int, List<int>>>();

            dataBaseKeys = data.Keys.ToArray();
            //List<string> triples = GenerateTriplets();
            for (int j = 0; j < dataBaseKeys.Length; j++)
            {
                List<double> aux = data[dataBaseKeys[j]];

                for (int i = 0; i < aux.Count; i++)
                {
                    if (locBase[i].ContainsKey(aux[i]))
                        locBase[i][aux[i]].Value.Add(j);
                    else
                    {
                        List<int> xaux = new List<int>();
                        xaux.Add(j);
                        locBase[i].Add(aux[i], new KeyValuePair<int, List<int>>(1, xaux));
                    }
                }

            }

            locBaseTab = new Dictionary<double, KeyValuePair<int, int[]>>[locBase.Length];            
            for (int i = 0; i < locBase.Length; i++)
            {
                locBaseTab[i] = new Dictionary<double, KeyValuePair<int, int[]>>();
                foreach (var item in locBase[i])
                {
                    KeyValuePair<int, List<int>> aux = item.Value;
                    KeyValuePair<int, int[]> newAux = new KeyValuePair<int, int[]>(aux.Key, aux.Value.ToArray());
                    locBaseTab[i].Add(item.Key, newAux);
                    //oTab[i, (int)item.Key] = locBaseTab[i][item.Key].Value;
                }
            }

        }
        public override int[] GetDistance(string refStructure, List<string> structures)
        {
            List<double> prof1 = data[refStructure];

            int[] locDist = new int[structures.Count];
            int[] cDist = new int[structures.Count];
            Dictionary<double, List<int>> b = null;
            List<int> aux;
            int counter = 0;
            for (int i = 0; i < prof1.Count; i++)
            {
                b = hashDataBase[i];
                
                if (b.ContainsKey(prof1[i]))
                {
                    List<double> keys = new List<double>(b.Keys);
                    keys.Remove(prof1[i]);
                    aux = b[(int)prof1[i]];

                    for (int j = 0; j < aux.Count; j++)
                    {
                       // locDist[aux[j]]++;
                        cDist[aux[j]]++;
                    }

                    if (keys.Count > 0)
                    {
                        aux = b[keys[0]];
                        for (int j = 0; j < aux.Count; j++)
                        {
                            locDist[aux[j]]--;
                            //cDist[aux[j]]++;
                        }
                    }
                    counter++;
                }

            }
            for (int i = 0; i < locDist.Length; i++)
                //locDist[i] = (int)((1 - ((double)locDist[i]) / counter)*100);
                locDist[i] += counter + minStateDataBase[i] - 2 * cDist[i];

            return locDist;
        }
        public override int[][] GetDistance(List<string> refStructure, List<string> structures)
        {
            int[][] dist = new int[refStructure.Count][];
            rotStruct.Clear();
            for (int j = 0; j < refStructure.Count; j++)
            {
                dist[j] = new int[structures.Count];

                //foreach (string item in structures)
                dist[j] = GetDistance(refStructure[j], structures);

            }

            return dist;

        }
        protected override void CalcMatrix(object o)
        {
            int num = (int)o;
            foreach (var item in indexList[num])
            {
                int[] distC = GetDistance(structures[item], structures);

                for (int i = item + 1; i < structures.Count; i++)
                {
                    if (distC[i] < 0 || distC[i] >= int.MaxValue)
                        distC[i] = -1;
                    int index = GetIndex(item, i);
                    distanceMatrix[index] = distC[i];
                    Interlocked.Increment(ref currentV);
                }
            }
            if (resetEvents != null)
                resetEvents[num].Set();
        }
        public override void CalcDistMatrix(List<string> str)
        {
            currentV = 0;

            this.structures = new List<string>(data.Keys);

            hashIndex = new Dictionary<string, int>(structures.Count);
            for (int i = 0; i < structures.Count; i++)
                hashIndex.Add(structures[i], i);

            Settings set = new Settings();
            set.Load();
            int threadNumbers = set.numberOfCores;

            int part = structures.Count * (structures.Count + 1) / (2 * threadNumbers) + 1;

            distanceMatrix = new int[structures.Count * (structures.Count + 1) / 2];
            indexList = new List<int>[threadNumbers];

            for (int n = 0; n < threadNumbers; n++)
                indexList[n] = new List<int>(part);

            double step = structures.Count / threadNumbers;
            for (int i = 0; i < threadNumbers; i++)
                for (int j = (int)(i * step); j < step * i + step; j++)
                {
                    indexList[i].Add(j);
                }

            for (int n = 0; n < threadNumbers; n++)
                maxV += indexList[n].Count;

            resetEvents = new ManualResetEvent[threadNumbers];
            for (int n = 0; n < threadNumbers; n++)
            {
                int p;
                p = n;
                resetEvents[n] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(new WaitCallback(CalcMatrix), (object)p);

            }
            for (int j = 0; j < threadNumbers; j++)
                resetEvents[j].WaitOne();
            DebugClass.WriteMessage("Distance finished");
            currentV = maxV;
        }
        protected void FastSort(int[] dist, int[] index)
        {
            //int[] aux = new int[dist.Length];
            int[] aux = new int[101];
            int i = 0;
            for (i = 0; i < dist.Length; i++)
                aux[dist[i]]++;

            i = 0;
            int num = 0, border = 0;
            while (border < resSize && num < aux.Length)
                border += aux[num++];

            int[] sorted = new int[border];
            int counter = 0;
            for (i = 0; i < dist.Length; i++)
            {
                if (dist[i] < num)
                {
                    sorted[counter] = dist[i];
                    index[counter++] = i;
                }
            }

            Array.Sort(sorted, index);


        }

    }

}
