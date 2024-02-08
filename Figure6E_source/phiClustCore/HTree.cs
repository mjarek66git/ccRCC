using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using phiClustCore.Interface;
using System.Data;

namespace phiClustCore
{
    class HTree:IProgressBar
    {
        HashCluster hCluster;
        Options opt;
        jury1D juryLocal;
        Dictionary<string, List<double>> data;
        protected int maxV = 100;
        protected int currentV = 0;
        double startProgress = 0.5;
        double endProgress = 1;

        public double StartProgress { set { startProgress = value; } get { return startProgress; } }
        public double EndProgress { set { endProgress = value; } get { return endProgress; } }
        public HTree(Dictionary<string,List<double>> data, HashCluster hCl)
        {
            hCluster = hCl;
            this.data = data;
            opt = hCl.opt;
            juryLocal = new jury1D();

        }
        public override string ToString()
        {
            return "phiClust:HTree";
        }
        public double ProgressUpdate()
        {
            double res = 0;
            double progress = hCluster.ProgressUpdate();

            if (progress == 1)
                res = 0.05 + progress * 0.7;
            else
                return 0.5 * progress;
            return StartProgress + (EndProgress - StartProgress)  * ((double)currentV / maxV);
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }

        HClusterNode JoinNodes(List<HClusterNode> nodes)
        {
            ClusterOutput output;
            HClusterNode aux = new HClusterNode();
            aux.setStruct = new List<string>();            
            aux.stateFreq = nodes[0].stateFreq;
            if (nodes.Count > 1)
            {
                foreach (var it in nodes)
                {
                    aux.setStruct.AddRange(it.setStruct);
                    if (!it.Equals(nodes[0]))
                    {
                        for (int i = 0; i < aux.stateFreq.Length; i++)
                            foreach (var st in it.stateFreq[i])
                                if (aux.stateFreq[i].ContainsKey(st.Key))
                                    aux.stateFreq[i][st.Key] += st.Value;
                                else
                                    aux.stateFreq[i].Add(st.Key, 1);
                    }
                    it.parent = aux;
                }
                output = juryLocal.JuryOptWeights(aux.setStruct, aux.stateFreq);
                aux.consistency = hCluster.CalcClusterConsistency(aux.setStruct);
                aux.refStructure = output.juryLike[0].Key;
            }
            else
            {
                aux.setStruct = nodes[0].setStruct;
                aux.stateFreq = nodes[0].stateFreq;
                aux.refStructure = nodes[0].refStructure;
                aux.consistency = nodes[0].consistency;
            }            
            aux.joined = nodes;
            aux.parent = null;
            
            return aux;
        }

        public ClusterOutput RunHTree()
        {
            HClusterNode root = null;
            ClusterOutput outClust=hCluster.RunHashCluster();
            Dictionary<string, List<int>> clusters = hCluster.dicFinal;            
            
            juryLocal.PrepareJury();
            ClusterOutput output;
            List<HClusterNode> groundLevel = new List<HClusterNode>();

            foreach (var item in clusters)
            {
                HClusterNode aux = new HClusterNode();
                aux.parent = null;
                aux.joined = null;
                aux.setStruct = new List<string>(item.Value.Count+1);
                foreach (var index in item.Value)
                    aux.setStruct.Add(hCluster.structNames[index]);

                output = juryLocal.JuryOptWeights(aux.setStruct);
                aux.stateFreq = juryLocal.columns;
                aux.refStructure = output.juryLike[0].Key;                     
                
                aux.realDist = 0;
                aux.levelDist = 0;
                aux.dirName = item.Key;
                aux.consistency = hCluster.CalcClusterConsistency(aux.setStruct);
                groundLevel.Add(aux);

            }
            currentV = 20;            
            int size = groundLevel[0].dirName.Length;
            double step = 80.0 / size;
            for (int i = 1; i < size; i++)
            {

                Dictionary<string, List<HClusterNode>> join = new Dictionary<string, List<HClusterNode>>();
                foreach (var item in groundLevel)
                {
                    string key = item.dirName;
                    key = key.Substring(0, key.Length -1);
                    if (join.ContainsKey(key))
                        join[key].Add(item);
                    else
                    {
                        List<HClusterNode> aux = new List<HClusterNode>();
                        aux.Add(item);
                        join.Add(key, aux);
                    }
                }
                groundLevel = new List<HClusterNode>();
                foreach (var item in join)
                {
                    HClusterNode aux = JoinNodes(item.Value);
                    aux.levelDist = aux.levelNum = i;
                    aux.realDist = i;
                    aux.dirName = item.Key;
                    groundLevel.Add(aux);                
                }
                currentV += (int)step;
            }
            if (groundLevel.Count == 1)
                root = groundLevel[0];
            else
            {
                root = JoinNodes(groundLevel);
                root.levelDist = size;
                root.realDist = size;
            }
            currentV = 100;

            ClusterOutput outHTree = new ClusterOutput();
            outHTree.hNode = root;
            return outHTree;
        }
    }
}
