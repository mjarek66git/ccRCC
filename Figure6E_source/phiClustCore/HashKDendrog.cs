using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore
{
    class HashKDendrog : HashCluster
    {
        public HashKDendrog(Dictionary<string,List<double>> data, Options input) : base(data, input) { }


        public override ClusterOutput Cluster(List<string> _structNames, bool dendrog = false)
        {

            structNames = new List<string>(_structNames.Count);
            ClusterOutput output = new ClusterOutput();

            PrepareClustering(_structNames);
            //output = PrepareClustersJuryLike(dicC, structNames);

            output = MakeDendrogram(dicC, this.input.relClusters, this.input.perData / 100.0);
            // output.clusterConsisten = CalcClustersConsistency(output.hNode.GetLeaves());

            currentV = maxV;
            return output;
        }

        private ClusterOutput MakeDendrogram(Dictionary<string, List<int>> dic, int k, double prec)
        {
            if (columns == null)
                return null;
            DebugClass.WriteMessage("Select");
            bool[] columnAvoid = new bool[columns.Length];
            for (int i = 0; i < columns.Length; i++)
                columnAvoid[i] = false;

            int[] indexes = new int[columnAvoid.Length];

            for (int j = 0; j < indexes.Length; j++)
                indexes[j] = j;

            double[] entropy = CalcEntropy(columns);

            Array.Sort(entropy, indexes);

            Dictionary<string, List<int>> hashClusters = new Dictionary<string, List<int>>();
            int position = 0;
            int positionLeft = 0;
            int positionRight = indexes.Length;

            do
            {
                position = (positionRight + positionLeft) / 2;

                hashClusters.Clear();
                foreach (var item in dic.Keys)
                {
                    StringBuilder keyB = new StringBuilder();
                    string key = "";
                    for (int i = 0; i < item.Length; i++)
                    {
                        if (columnAvoid[i])
                            continue;

                        keyB.Append(item[i]);
                    }
                    key = keyB.ToString();
                    if (!hashClusters.ContainsKey(key))
                        hashClusters.Add(key, new List<int>());

                    hashClusters[key].AddRange(dic[item]);
                }
                if (hashClusters.Count == dic.Count)
                    positionRight = position;
                else
                    positionLeft = position;

            }
            while (positionRight - positionLeft > 1);

            Dictionary<string, HClusterNode> hNodes = new Dictionary<string, HClusterNode>();
            foreach (var item in hashClusters)
            {
                HClusterNode aux = new HClusterNode();
                aux.setStruct = new List<string>();
                foreach (var index in item.Value)
                    aux.setStruct.Add(structNames[index]);

                hNodes.Add(item.Key, aux);
            }
            Dictionary<string, HClusterNode> prevLevel = hNodes;
            hNodes = new Dictionary<string, HClusterNode>();
            do
            {
                foreach (var item in prevLevel)
                {
                    string key = item.Key.Substring(0, item.Key.Length - 2);
                    if (!hNodes.ContainsKey(key))
                    {
                        HClusterNode aux = new HClusterNode();
                        aux.setStruct.AddRange(item.Value.setStruct);
                        aux.joined.Add(item.Value);
                        hNodes.Add(key, aux);
                    }
                    else
                    {
                        hNodes[key].setStruct.AddRange(item.Value.setStruct);
                        hNodes[key].joined.Add(item.Value);
                    }

                }
                if (hNodes.Count < prevLevel.Count)
                    prevLevel = hNodes;
            }
            while (hNodes.Count != 2);
            HClusterNode auxRoot = null;
            if (hNodes.Count > 0)
            {
                auxRoot = new HClusterNode();
                foreach (var item in hNodes)
                {
                    auxRoot.setStruct.AddRange(item.Value.setStruct);
                    auxRoot.joined.Add(item.Value);
                }
            }
            else
            {
                List<string> xx = new List<string>(hNodes.Keys);
                auxRoot = hNodes[xx[0]];
            }

            ClusterOutput outp = new ClusterOutput();

            outp.hNode = auxRoot;

            return outp;

        }
    }
}
