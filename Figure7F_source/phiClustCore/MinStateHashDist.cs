using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore
{
    class MinStateHash:RootHash
    {
        protected int[] minStateDataBase = null;
        protected Dictionary<string, KeyValuePair<int, List<int>>>[] locBase = null;
        protected Dictionary<string, KeyValuePair<int, int[]>>[] locBaseTab = null;

        public MinStateHash(int binSizeG, HashCluster hk, ClusterOutput outp, HNNCInput opt):base(binSizeG,hk, outp, opt)
        {

        }
/*        public List<string> GenerateTriplets()
        {
            List<string> tr = new List<string>();
            string alphabet = "ARDNCEQGHILKMFPSTWYV";
            for(int i=0;i<alphabet.Length;i++)
                for(int j=0;j<alphabet.Length;j++)
                    for(int s=0;s<alphabet.Length;s++)
                    {
                        string w = alphabet[i].ToString() + alphabet[j] + alphabet[s];
                        tr.Add(w);
                    }

            locBase = new Dictionary<string, KeyValuePair<int, List<int>>>[hashDataBase.Length / binSizeG + 1];

            for (int i = 0; i < locBase.Length; i++)
                locBase[i] = new Dictionary<string, KeyValuePair<int, List<int>>>();

            for (int i=0;i<locBase.Length;i++)
            {
                foreach(var item in tr)
                {
                    int num = 0;
                    for (int n = 0; n < item.Length; n++)
                        if (hashDataBase[i * item.Length + n].ContainsKey(item[n].ToString()))
                            num++;

                    if (num == 0)
                        continue;

                    if (locBase[i].ContainsKey(item))
                        locBase[i][item].Value.Add(j);
                    else
                    {
                        List<int> xaux = new List<int>();
                        xaux.Add(j);
                        locBase[i].Add(item, new KeyValuePair<int, List<int>>(num, xaux));
                    }
                }
            }

            return tr;
        }*/
        public override void CreateBase(Dictionary<string, string> dataBase)
        {
            int remBin = binSizeG;
            binSizeG = 1;

            base.CreateBase(dataBase);
            binSizeG = remBin;
            minStateDataBase = new int[dataBaseKeys.Count];
            for (int j = 0; j < minStateDataBase.Length; j++)
                minStateDataBase[j] = 0;

            for (int j = 0; j < hashDataBase.Length; j++)
            {
                string remKey = "";
                int min = int.MaxValue;
                Dictionary<string, List<int>> locHash = hashDataBase[j];
                foreach (var item in hashDataBase[j].Keys)
                {
                    if (locHash[item].Count < min)
                    {
                        min = locHash[item].Count;
                        remKey = item;
                    }
                }

                foreach (var item in locHash[remKey])
                    minStateDataBase[item]++;

                List<int> aux = hashDataBase[j][remKey];
                locHash.Clear();
                locHash.Add(remKey, aux);
            }
           // binSizeG = 1;


            locBase = new Dictionary<string, KeyValuePair<int, List<int>>>[hashDataBase.Length / binSizeG + 1];

            for (int i = 0; i < locBase.Length; i++)
                locBase[i] = new Dictionary<string, KeyValuePair<int, List<int>>>();

            dataBaseKeys = new List<string>(dataBase.Keys);
            //List<string> triples = GenerateTriplets();
            for (int j = 0; j < dataBaseKeys.Count; j++)
            {
                List<string> aux = GetKeyHashes(dataBase[dataBaseKeys[j]]);

                for (int i = 0; i < aux.Count; i++)
                {
                    int num = 0;
                    for (int n = 0; n < aux[i].Length; n++)
                        if (hashDataBase[i * aux[i].Length + n].ContainsKey(aux[i][n].ToString()))
                            num++;

                    if (num == 0)
                        continue;

                    if (locBase[i].ContainsKey(aux[i]))
                        locBase[i][aux[i]].Value.Add(j);
                    else
                    {
                        List<int> xaux = new List<int>();
                        xaux.Add(j);
                        locBase[i].Add(aux[i], new KeyValuePair<int, List<int>>(num, xaux));
                    }
                }

            }

            locBaseTab =new Dictionary<string, KeyValuePair<int, int[]>>[locBase.Length];
            for(int i=0;i<locBase.Length;i++)
            {
                locBaseTab[i] = new Dictionary<string, KeyValuePair<int, int[]>>();
                foreach(var item in locBase[i])
                {
                    KeyValuePair<int, List<int>> aux = item.Value;                    
                    KeyValuePair<int, int[]> newAux = new KeyValuePair<int, int[]>(aux.Key,aux.Value.ToArray());
                    locBaseTab[i].Add(item.Key, newAux);
                }
            }
        }
        public override string CalcDist(int threadNum, List<string> keys, int[] index, int num)
        {
            KeyValuePair<int, int[]> aux;
            int[] locDist = dist[threadNum];
            Dictionary<string, KeyValuePair<int, int[]>> b = null;
            string k = "";
            int i = 0;
            //foreach(var item in keys)
            for (i = 0; i < keys.Count; i++)
            {
                k = keys[i];
                b = locBaseTab[i];
                if (b.ContainsKey(k))
                {
                    aux = b[k];
                    int x = aux.Key;
                    int[] common = aux.Value;
                    for (int j = 0; j < common.Length; j++)
        //            foreach(var item in common)
                        locDist[common[j]] += x;
                }

            }
            for (i = 0; i < locDist.Length; i++)
                locDist[i] = minStateDataBase[index[i]] - 2 * locDist[i];

            //Array.Sort(dist[threadNum], index);
            string w = "";
            for (i = 0; i < resSize; i++)
                w += dataBaseKeys[index[i]] + ":";
            w += dataBaseKeys[index[resSize]];

            return w;
        }

        public string CalcDistN(int threadNum, List<string> keys, int[] index, int num)
        {            
            KeyValuePair<int, List<int>> aux;
            int []locDist = dist[threadNum];
            Dictionary<string, KeyValuePair<int, List<int>>> b = null;
            string k = "";
            int i = 0;
            //foreach(var item in keys)
            for (i = 0; i < keys.Count; i++)
            {
                k = keys[i];
                b = locBase[i];
                if (b.ContainsKey(k))
                {
                    aux = b[k];
                    int x = aux.Key;
                    List<int> common = aux.Value;
                    for(int j=0;j<common.Count;j++)
                        locDist[common[j]] += x;
                }

            }
            for (i = 0; i < locDist.Length; i++)
                locDist[i] = minStateDataBase[index[i]] - 2 * locDist[i];

            //Array.Sort(dist[threadNum], index);
             string w = "";
            for (i = 0; i < resSize; i++)
                w += dataBaseKeys[index[i]] + ":";
            w += dataBaseKeys[index[resSize]];

            return w;
        }


    }
}
