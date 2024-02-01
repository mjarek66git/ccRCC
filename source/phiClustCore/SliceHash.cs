using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore
{
    class SliceHash : MinStateHash
    {
        int rad = 10;
        Dictionary<int, Dictionary<string, KeyValuePair<int, List<int>>>[]> sliceDic = null;
        SliceHash(int binSizeG, HashCluster hk, ClusterOutput outp, HNNCInput opt) : base(binSizeG, hk, outp, opt)
        {

        }
        public override void CreateBase(Dictionary<string, string> dataBase)
        {
            base.CreateBase(dataBase);
            
            Dictionary<int, List<int>> sliceBase = new Dictionary<int, List<int>>();

            for (int j = 0; j < dataBaseKeys.Count; j++)
            {
                int count = 0;
                string xx = dataBase[dataBaseKeys[j]];
                for (int i = 0; i < xx.Length; i++)
                    if (hashDataBase[i].ContainsKey(xx[i].ToString()))
                        count++;
                if (sliceBase.ContainsKey(count))
                    sliceBase[count].Add(j);
                else
                {
                    List<int> inx = new List<int>();
                    inx.Add(j);
                    sliceBase.Add(count, inx);
                }
            }

            sliceDic = new Dictionary<int, Dictionary<string, KeyValuePair<int, List<int>>>[]>();


            foreach (var item in sliceBase)
            {
                Dictionary<string, KeyValuePair<int, List<int>>>[] aux = new Dictionary<string, KeyValuePair<int, List<int>>>[locBase.Length];
                for (int i = 0; i < locBase.Length; i++)
                    aux[i] = new Dictionary<string, KeyValuePair<int, List<int>>>();

                foreach (var sel in item.Value)
                {
                    List<string> auxW = GetKeyHashes(dataBase[dataBaseKeys[sel]]);
                    for (int i = 0; i < auxW.Count; i++)
                    {
                        int num = 0;
                        for (int n = 0; n < auxW[i].Length; n++)
                            if (hashDataBase[i * auxW[i].Length + n].ContainsKey(auxW[i][n].ToString()))
                                num++;

                        if (num == 0)
                            continue;

                        if (aux[i].ContainsKey(auxW[i]))
                            aux[i][auxW[i]].Value.Add(sel);
                        else
                        {
                            List<int> xaux = new List<int>();
                            xaux.Add(sel);
                            aux[i].Add(auxW[i], new KeyValuePair<int, List<int>>(num, xaux));
                        }
                    }

                }
                sliceDic.Add(item.Key, aux);

            }
        }

        public override string CalcDist(int threadNum, List<string> keys, int[] index, int num)
        {
            int[] locDist = dist[threadNum];
            List<int> baseHash = null;
            KeyValuePair<int, List<int>> aux;

            for (int z = Math.Max(0, num - rad); z < Math.Min(num + rad, sliceDic.Count); z++)
            {
                if (!sliceDic.ContainsKey(z))
                    continue;
                Dictionary<string, KeyValuePair<int, List<int>>>[] auxBase = sliceDic[z];

                for (int i = 0; i < keys.Count; i++)
                {
                    if (auxBase[i].TryGetValue(keys[i], out aux))
                    {
                        baseHash = aux.Value;
                        int val = aux.Key;
                        foreach (var inx in baseHash)
                            locDist[inx] += val;
                    }

                }
            }
            int vCounter = 0;
            Queue<int> q = new Queue<int>();

            for (int i = 0; i < dist[threadNum].Length; i++)
            {
                if (dist[threadNum][i] == 0)
                    q.Enqueue(i);
                else
                {
                    if (q.Count > 0)
                    {
                        int idx = q.Dequeue();
                        index[idx] = i;
                        dist[threadNum][idx] = dist[threadNum][i];
                        dist[threadNum][i] = 0;
                        q.Enqueue(i);
                    }
                    vCounter++;

                }
            }
            for (int i = 0; i < dist[threadNum].Length; i++)
                dist[threadNum][i] = minStateDataBase[index[i]] - 2 * dist[threadNum][i];

            Array.Sort(dist[threadNum], index,0,vCounter);
            string w = "";
            for (int i = 0; i < resSize; i++)
                w += dataBaseKeys[index[i]] + ":";
            w += dataBaseKeys[index[resSize]];

            return w;            
        }
    }
}
