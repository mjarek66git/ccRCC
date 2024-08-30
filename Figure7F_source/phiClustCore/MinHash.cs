using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace phiClustCore
{
    public class MinHash:HNN
    {
        Dictionary<string, Dictionary<string, int>> clusters = null;// Only for testing cath4.2
        // Constructor passed universe size and number of hash functions
        public MinHash(int universeSize, int numHashFunctions, HashCluster hk, ClusterOutput outp, HNNCInput opt):base(hk, outp, opt)
        {
            this.numHashFunctions = numHashFunctions;
            // number of bits to store the universe
        }

        private int numHashFunctions;

        // Returns number of hash functions defined for this instance
        public int NumHashFunctions
        {
            get { return numHashFunctions; }
        }

        public delegate uint Hash(int toHash);
        private Hash[] hashFunctions;

        // Public access to hash functions
        public Hash[] HashFunctions
        {
            get { return hashFunctions; }
        }

        // Generates the Universal Random Hash functions
        // http://en.wikipedia.org/wiki/Universal_hashing
        private void GenerateHashFunctions(int u)
        {
            hashFunctions = new Hash[numHashFunctions];

            // will get the same hash functions each time since the same random number seed is used
            Random r = new Random(10);
            for (int i = 0; i < numHashFunctions; i++)
            {
                uint a = 0;
                // parameter a is an odd positive
                while (a % 1 == 1 || a <= 0)
                    a = (uint)r.Next();
                uint b = 0;
                int maxb = 1 << u;
                // parameter b must be greater than zero and less than universe size
                while (b <= 0) b = (uint)r.Next(maxb); hashFunctions[i] = x => QHash(x, a, b, u);
            }
        }

        // Returns the number of bits needed to store the universe
        public int BitsForUniverse(int universeSize)
        {
            return (int)Math.Truncate(Math.Log((double)universeSize, 2.0)) + 1;
        }

        // Universal hash function with two parameters a and b, and universe size in bits
        private static uint QHash(int x, uint a, uint b, int u)
        {
            double aa = a * (uint)x + b;
            uint res = a * (uint)x + b;
            if (res != aa)
                Console.Write("ups");
            return (a * (uint)x + b) % 1753;// 179;//>> (32 - u);
        }

        // Returns the list of min hashes for the given set of word Ids
        public uint[] GetMinHash(List<int> wordIds)
        {
            uint[] minHashes = new uint[numHashFunctions];
            for (int h = 0; h < numHashFunctions; h++)
            {
                minHashes[h] = int.MaxValue;
            }
            foreach (int id in wordIds)
            {
                for (int h = 0; h < numHashFunctions; h++)
                {
                    uint hash = hashFunctions[h](id);
                    
                    minHashes[h] = Math.Min(minHashes[h], hash);
                }
            }
            return minHashes;
        }
        public double Similarity(IEnumerable<int> A, IEnumerable<int> B)
        {
            var CommonNumbers = from a in A.AsEnumerable<int>()
                                join b in B.AsEnumerable<int>() on a equals b
                                select a;
            double JaccardIndex = 2 * (((double)CommonNumbers.Count()) /
                                   ((double)(A.Count() + B.Count())));

            return JaccardIndex;
        }
        public double Similarity(uint [] A, uint []B)
        {
            double count = 0;
            for(int i=0;i<A.Length;i++)
                if(A[i]==B[i])
                {
                    count++;
                }

            return 2*count/(A.Length+B.Length);
        }

        Dictionary<string,List<int>> CreateBase(Dictionary<string,string> dataBase)
        {
            Dictionary<string,List<int>> minBase = new Dictionary<string, List<int>>();
            List<string> keys = new List<string>(dataBase.Keys);
            foreach (var item in keys)
            {
                List<int> tmp = new List<int>();
                for (int j = 0; j < item.Length; j++)
                    if (item[j] == '1')
                        tmp.Add(j);

                minBase.Add(dataBase[item],tmp);
            }
            return minBase;
        }
        public static Dictionary<string, Dictionary<string, int>> ReadClusters(List<string> data)
        {
            Dictionary<string, Dictionary<string, int>> res = new Dictionary<string, Dictionary<string, int>>();
            StreamReader st = new StreamReader("C:\\Projects\\cath4.2-domain-list.txt");
            Dictionary<string, List<string>> xx = new Dictionary<string, List<string>>();
            Dictionary<string, int> dataDic = new Dictionary<string, int>();

            foreach (var item in data)
                dataDic.Add(item, 0);

            string line = "";
            string clusterStr = "";
            line = st.ReadLine();
            while (line != null)
            {
                line = line.Replace('\t', ' ');
                line = Regex.Replace(line, @"\s+", " ");
                string[] aux = line.Split(' ');
                if (dataDic.ContainsKey(aux[0]))
                {
                    clusterStr = aux[1] + "." + aux[2] + "." + aux[3] + "." + aux[4];
                    if (xx.ContainsKey(clusterStr))
                        xx[clusterStr].Add(aux[0]);
                    else
                    {
                        List<string> qq = new List<string>();
                        qq.Add(aux[0]);
                        xx.Add(clusterStr, qq);
                    }
                }
                line = st.ReadLine();
            }

            st.Close();
            foreach (var item in xx)
            {
                Dictionary<string, int> w = new Dictionary<string, int>();
                for (int i = 0; i < item.Value.Count; i++)
                    w.Add(item.Value[i], 0);

                for (int i = 0; i < item.Value.Count; i++)
                    res.Add(item.Value[i], w);
            }
            return res;
        }

        //All data must be coded by binary vector
        public override Dictionary<string, string> HNNTest(List<string> testList)
        {
            Dictionary<string,string> aux=BasesForTest(testList);
            List<string> kk = new List<string>(aux.Keys);
            int u = BitsForUniverse(kk[0].Length+1);
            GenerateHashFunctions(u);
            clusters=ReadClusters(testList);

            Dictionary<string, string> res = new Dictionary<string, string>();
            if (aux.Count>0)
            {
                Dictionary<string,List<int>> minBase = CreateBase(caseBase);
                Dictionary<string, List<int>> retriveBase = CreateBase(aux);
                StringBuilder final = new StringBuilder();
                StreamWriter wr = new StreamWriter("xxx");

                List<string> keysminBase = new List<string>(minBase.Keys);
                Dictionary<string, int> keyPosition = new Dictionary<string, int>();

                int t = 0;
                foreach (var r in keysminBase)
                    keyPosition.Add(r, t++);


                double[] dist = new double[minBase.Count];

                uint[][] minHashDataBase = new uint[keysminBase.Count][];
                for (int j = 0; j < keysminBase.Count; j++)
                {
                    minHashDataBase[j] = GetMinHash(minBase[keysminBase[j]]);
                }

                    int[] index = new int[minBase.Count];
                foreach(var item in retriveBase.Keys)
                {
                    final.Clear();
                    uint[] st1 = GetMinHash(retriveBase[item]);
                    for (int j = 0; j < keysminBase.Count; j++)
                    {
                    
                        dist[j] = Similarity(st1, minHashDataBase[j]);

                        if (keysminBase[j] == item && dist[j]<1.0)
                            Console.Write("jjjd");

                    }
                    for (int j = 0; j < index.Length; j++)
                        index[j] = j;
                    Array.Sort(dist, index);
                    Array.Reverse(dist);
                    Array.Reverse(index);
//                    string str1 = item;
//                    string str2 = keysminBase[index[0]];
                    if (dist[0] < 1)
                        Console.Write("UPS");
                    if (clusters != null)
                    {
                        Dictionary<int, int> ix = new Dictionary<int, int>();

                        for (int g = 0; g < index.Length; g++)
                            ix.Add(index[g], g);
                        Dictionary<string, int> elements = clusters[item];
                        List<int> pos = new List<int>();
                        wr.WriteLine(">" +item);
                        foreach (var e in elements.Keys)
                        {
                            //int vv = locKeys.FindIndex(z => z.Contains(e));
                            int vv = -1;
                            if (keyPosition.ContainsKey(e))
                                vv = keyPosition[e];
                            if (vv == -1 || e == item)
                                continue;

                            if (vv != -1)
                            {
                                if (ix[vv] != 0)                                    
                                    pos.Add(ix[vv]);
                            }
                        }
                        pos.Sort();
                        double il = 1.0, norm = 1.0;
                        double p = 1.0 / pos.Count;
                        int co = 1;
                        foreach (var ss in pos)
                        {
                            wr.Write(ss + " ");
                            il *= ((double)ss) / co;

                            co++;
                        }
                        double avr = Math.Pow(il, p);
                        norm = Math.Pow(norm, p);
                        avr = norm / avr;
                        wr.WriteLine();
                        foreach (var ss in pos)
                            wr.Write(keysminBase[index[ss]] + " ");
                        wr.WriteLine();
                        foreach (var ss in pos)
                            wr.Write(dist[ss] + " ");
                        wr.WriteLine();
                        if (pos.Count > 0)
                            wr.WriteLine("Avr=" + avr);

                    }
                    else
                    {
                        Dictionary<string, int> c = new Dictionary<string, int>();
                        for (int j = 0; j < 10; j++)
                            if (!c.ContainsKey(keysminBase[index[j]]))
                                c.Add(keysminBase[index[j]], 1);
                        List<string> xx = new List<string>(c.Keys);

                        for (int j = 0; j < xx.Count - 1; j++)
                            final.Append(xx[j] + ":");
                        final.Append(xx[xx.Count - 1]);


                        if (final.Length == 0)
                            res.Add(item, "NOT CLASSIFIED");
                        else
                            res.Add(item, final.ToString());
                    }
                }



            }
           return res;
        }

    }
}
