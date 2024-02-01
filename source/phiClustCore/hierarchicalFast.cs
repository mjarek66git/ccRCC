using phiClustCore.Distance;
using phiClustCore.Interface;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore
{
    class hierarchicalFast : IProgressBar
    {
        DistanceMeasure dMeasure;
        AglomerativeType linkageType;
        Options opt;
        public string mustRefStructure = null;
        int min;
        int currentV = 0;
        int maxV = 1;
        int progressRead = 0;
        HierarchicalCInput hierOpt;

        double startProgress = 0;
        double endProgress = 1;

        public double StartProgress { set { startProgress = value; } get { return startProgress; } }
        public double EndProgress { set { endProgress = value; } get { return endProgress; } }


        public hierarchicalFast(DistanceMeasure dMeasure, Options opt)
        {
            this.dMeasure = dMeasure;
            this.linkageType = opt.hierarchical.linkageType;
            progressRead = 0;
            hierOpt = opt.hierarchical;
            this.opt = opt;

        }
        public double ProgressUpdate()
        {
            double sumProgress = 0;
            double progress = dMeasure.ProgressUpdate();

            if (progressRead == 1)
                sumProgress = 0.05 + progress * 0.7;
            else
                sumProgress = 0.05 * progress;

            return StartProgress + (EndProgress - StartProgress) * (sumProgress + 0.25 * ((double)currentV / maxV));

        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }


        public override string ToString()
        {
            return "Agglomerative " + linkageType.ToString();
        }
        private int LMinimalDist(List<HClusterNode> levelNodes)
        {

            int min = Int32.MaxValue;
            for (int i = 0; i < levelNodes.Count; i++)
            {
                HClusterNode refStruct = levelNodes[i];
                for (int j = i + 1; j < levelNodes.Count; j++)
                {
                    int dist = dMeasure.FindMinimalDistance(refStruct, levelNodes[j], linkageType).Key;
                    if (min > dist)
                    {
                        min = dist;
                    }
                }
            }
            return min;
        }
        public ClusterOutput HierarchicalClustering(List<string> structures, List<Dictionary<double, int>[]> freq = null)
        {
            List<List<HClusterNode>> level = new List<List<HClusterNode>>();
            Dictionary<string, HClusterNode> dicH = new Dictionary<string, HClusterNode>();
            Dictionary<string, HClusterNode> dicOrg = new Dictionary<string, HClusterNode>();
            List<HClusterNode> rowNodes = new List<HClusterNode>();
            ClusterOutput outCl = new ClusterOutput();
            int levelCount = 0;
            bool end = false;
            HClusterNode node;
            if (structures.Count <= 1)
            {
                outCl.hNode = new HClusterNode();
                outCl.hNode.setStruct = structures;
                outCl.hNode.refStructure = structures[0];
                outCl.hNode.levelDist = 0;
                outCl.hNode.joined = null;
                return outCl;
            }

            progressRead = 1;
            dMeasure.CalcDistMatrix(structures);
            List<HClusterNode> levelNodes = new List<HClusterNode>();
            jury1D jury = new jury1D();
            jury.PrepareJury(((HammingBase)dMeasure).data);
            for (int i = 0; i < structures.Count; i++)
            {
                node = new HClusterNode();
                node.refStructure = structures[i];
                dicH.Add(structures[i], node);
                dicOrg.Add(structures[i], node);
                node.joined = null;
                node.setStruct.Add(structures[i]);
                node.levelNum = levelCount;
                if (freq != null && freq[i] != null)
                    node.stateFreq = freq[i];
                else
                {
                    //node.stateFreq = new Dictionary<byte, int>[];
                    jury.JuryOptWeights(node.setStruct);
                    node.stateFreq = jury.columns;

                }
                node.levelDist = dMeasure.maxSimilarity;
                node.realDist = dMeasure.GetRealValue(node.levelDist);
                levelNodes.Add(node);
            }
            List<HClusterNode> cc = new List<HClusterNode>();

            List<Tuple<int, HClusterNode, HClusterNode>> distKK = dMeasure.FindMinimalDistanceReference(cc);

            maxV = levelNodes.Count + 1;
            level.Add(levelNodes);
            HashSet<string> availableStruct = new HashSet<string>(structures);
            List<Tuple<int, HClusterNode, HClusterNode>> distK = dMeasure.FindMinimalDistanceReference(levelNodes);
            HashSet<string> tabuRefernce = new HashSet<string>();
            HashSet<string> levelRef = new HashSet<string>();
            int n = 0;
            while (n<distK.Count)
            {
                levelNodes = new List<HClusterNode>();
                int min = distK[n].Item1;
                int k = n;

                while (k+1<distK.Count-1 && distK[k+1].Item1 == min)
                    k++;
                levelRef.Clear();

                //allIndx=new 
                for(int i=n;i<=k;i++)
                {
         
                        HashSet<HClusterNode> indexes = new HashSet<HClusterNode>();
                    if (tabuRefernce.Contains(distK[i].Item2.refStructure) || tabuRefernce.Contains(distK[i].Item3.refStructure))
                        continue;

                    if (levelRef.Contains(distK[i].Item2.refStructure) || levelRef.Contains(distK[i].Item3.refStructure))
                        continue;

                    indexes.Clear();
                    indexes.Add(dicH[distK[i].Item2.refStructure]);
                    indexes.Add(dicH[distK[i].Item3.refStructure]);

                    for (int z = i + 1; z <= k; z++)
                    {
                        var itt = distK[z];

                        if (indexes.Contains(dicH[distK[z].Item2.refStructure]))
                        {                            
                            if (!levelRef.Contains(distK[z].Item3.refStructure))
                                if(!tabuRefernce.Contains(distK[z].Item3.refStructure))
                                    indexes.Add(dicH[distK[z].Item3.refStructure]);
                        }
                        if (indexes.Contains(dicH[distK[z].Item3.refStructure]))
                        {
                            if (!levelRef.Contains(distK[z].Item2.refStructure))
                                if(!tabuRefernce.Contains(distK[z].Item2.refStructure))
                                indexes.Add(dicH[distK[z].Item2.refStructure]);
                        }
                    }

                    foreach(var item in indexes)
                        levelRef.Add(item.refStructure);                    

                    if (indexes.Count > 0)
                    {
                        HClusterNode nodeTmp = new HClusterNode();
                        nodeTmp.levelDist = min;
                        nodeTmp.realDist = dMeasure.GetRealValue(min);
                        nodeTmp.joined = new List<HClusterNode>();
                        HashSet<string> test = new HashSet<string>();
                        foreach (var item in indexes)
                        {
                            nodeTmp.setStruct.AddRange(item.setStruct);
                            nodeTmp.joined.Add(item);
                            //item.parent = nodeTmp;                           
                            item.parent = nodeTmp;
                            //allIndx.Add(item.refStructure);
                        }
                        var orderList = dMeasure.GetReferenceList(nodeTmp.setStruct);
                        for (int s = 0; s < orderList.Count; s++)
                            if (!tabuRefernce.Contains(orderList[s].Key))
                                nodeTmp.refStructure = orderList[s].Key;

                        foreach (var item in indexes)
                        {
                            foreach(var str in item.setStruct)
                                if (str != nodeTmp.refStructure)
                                    tabuRefernce.Add(str);
                        }
                            //nodeTmp.refStructure = orderList[0].Key;                      
                            //allIndx.Remove(nodeTmp.refStructure);
                        dicH[nodeTmp.refStructure] = nodeTmp;
                        //nodeTmp.joined = new List<HClusterNode>(indexes);
                        levelNodes.Add(nodeTmp);
                        //if(levelRef.Contains(nodeTmp.refStructure))
                        //    Console.WriteLine();
                        //levelRef.Add(nodeTmp.refStructure);
                        //allIndx.Add(nodeTmp.refStructure);
                    }
                   
                }
                n = k + 1;
                if (levelNodes.Count > 0)
                    level.Add(levelNodes);

                
            }


            outCl.hNode = level[level.Count - 1][0];
            outCl.hNode.levelNum = 0;

         

            if (level[level.Count - 1][0].setStruct.Count != structures.Count)
                throw new Exception("Wrong hierarchical clustering!");
            //At the end level num must be set properly
            Queue<HClusterNode> qq = new Queue<HClusterNode>();
            HClusterNode h;
            for (int i = 0; i < level.Count; i++)
                for (int j = 0; j < level[i].Count; j++)
                    level[i][j].fNode = true;

            for (int i = 0; i < level.Count; i++)
                    for (int j = 0; j < level[i].Count; j++)
                    if (level[i][j].fNode)
                    {
                        level[i][j].levelDist = Math.Abs(level[i][j].levelDist - dMeasure.maxSimilarity);
                        level[i][j].realDist = dMeasure.GetRealValue(level[i][j].levelDist);
                        level[i][j].fNode = false;
                    }




            qq.Enqueue(level[level.Count - 1][0]);
            while (qq.Count != 0)
            {
                h = qq.Dequeue();

                if (h.joined != null)
                    foreach (var item in h.joined)
                    {
                        item.levelNum = h.levelNum + 1;
                        qq.Enqueue(item);

                    }
            }

            outCl.clusters = null;

            outCl.juryLike = null;
            currentV = maxV;
            outCl.runParameters = hierOpt.GetVitalParameters();

            return outCl;
        }

    }

}

