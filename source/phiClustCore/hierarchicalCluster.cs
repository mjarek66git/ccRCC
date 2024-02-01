using System;
using System.Collections;
using System.Collections.Specialized;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using phiClustCore;
using System.Data;
using phiClustCore.Interface;
using phiClustCore.Distance;


namespace phiClustCore
{
    [Serializable]
	public class HClusterNode
	{        
		public List <string> setStruct=new List<string>();
        public List<List<double>> setProfiles = new List<List<double>>();
		public List <HClusterNode> joined=null;        
        [NonSerialized]
        public HClusterNode parent = null;
		public string refStructure;
        public string dirName;
        public Dictionary<double, int>[] stateFreq = null;
		public bool fNode=false;
        public bool flagSign = false;
        public bool visible = true;
		public int num;
		public int iterNum=0;
		public int levelNum=0;
		public int counter=0;
        public double consistency = 0;
        public Color color=Color.Green;
        public double levelDist;
        public double realDist;
        [NonSerialized]
        public GraphNode gNode;        
        int kMin;

        public HClusterNode()
        {

        }

        public HClusterNode(HClusterNode node)
        {
            if(node.setStruct!=null)
                setStruct = new List<string>(node.setStruct);
            if(node.joined!=null)
                joined = new List<HClusterNode>(node.joined);
            parent = node.parent;
            refStructure = node.refStructure;
            dirName = node.dirName;
            fNode = node.fNode;
            num = node.num;
            iterNum = node.iterNum;
            levelNum = node.levelNum;
            counter = node.counter;
            color = node.color;
            levelDist = node.levelDist;
            realDist = node.realDist;
            gNode = new GraphNode(node.gNode);
            kMin = node.kMin;
        }


        public void ClearColors(Color color)
        {
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;

            st.Push(this);

            while (st.Count != 0)
            {
                current = st.Pop();
                current.color = color;
                if (current.joined != null)
                {                    
                    foreach (var item in current.joined)
                        st.Push(item);
                }

            }
        }
        public List<HClusterNode> GetAllClusters()
        {
            List<HClusterNode> res = new List<HClusterNode>();
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;

            st.Push(this);

            while (st.Count != 0)
            {
                current = st.Pop();
                res.Add(current);
                if (current.joined != null)
                {
                    foreach (var item in current.joined)
                        st.Push(item);
                }

            }

            return res;
        }
        public void ColorNode(string structMark,int color)
        {
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;

            st.Push(this);

            while (st.Count != 0)
            {
                current = st.Pop();
                if (current.joined != null)
                {
                    if (current.setStruct.Contains(structMark))
                        current.color = Color.Red;
                    else
                        current.color = Color.Green;
                    foreach (var item in current.joined)
                        st.Push(item);
                }

            }

        }
        public Dictionary<HClusterNode,Color> MarkNodesAccordingToRules(List<int> upRegulated,List<int> downRegulated,double threshold,
                                                                        Color color,Dictionary<string,List<double>> profiles)
        {
            Dictionary<HClusterNode, Color> returnList = new Dictionary<HClusterNode, Color>();
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;

            if (upRegulated.Count == 0 && downRegulated.Count == 0)
                return null;

            st.Push(this);

            while (st.Count != 0)
            {
                current = st.Pop();
                double counterUp = 0;
                double counterDown = 0;
                bool test = true;

                foreach (var index in upRegulated)
                {
                    if (test)
                    {
                        counterUp = 0;
                        foreach (var item in current.setStruct)
                        {
                            if (profiles[item][index] > 0)
                                counterUp++;
                        }
                        counterUp /= current.setStruct.Count;
                        if (counterUp < threshold)
                        {
                            test = false;
                            break;
                        }
                    }

                }

                foreach (var index in downRegulated)
                {
                    if (test)
                    {
                        counterDown = 0;
                        foreach (var item in current.setStruct)
                            if (profiles[item][index]< 0)
                                counterDown++;
                        counterDown /= current.setStruct.Count;
                        if (counterDown < threshold)
                        {
                            test = false;
                            break;
                        }

                    }

                }
                if (test)
                    returnList.Add(current,color);

                if (current.joined != null)
                    foreach (var item in current.joined)
                        st.Push(item);

            }
            return returnList;
        }

        public Dictionary<HClusterNode,System.Drawing.Color> MarkNodes(List<string> toMark,System.Drawing.Color color)
        {
            Dictionary<HClusterNode, System.Drawing.Color> returnList = new Dictionary<HClusterNode, System.Drawing.Color>();
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;

            st.Push(this);

            while (st.Count != 0)
            {
                current = st.Pop();
                if (current.joined == null || current.joined.Count == 0)
                {
                    foreach(var item in toMark)
                        if(current.setStruct.Contains(item))
                        {
                            returnList.Add(current,color);
                            break;
                        }
                }
                else
                    if (current.joined != null)
                        foreach (var item in current.joined)
                            st.Push(item);

            }
            return returnList;
        }
        public int GetMaxCLusterSize()
        {
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;
            int kMax;
            kMax = -(int)this.levelDist;
            
            st.Push(this);
            while (st.Count != 0)
            {
                current = st.Pop();
                if (current.setStruct.Count > kMax)
                    kMax = current.setStruct.Count;

                if (current.joined != null)
                    foreach (var item in current.joined)
                        st.Push(item);
            }
            return kMax;
        }
        public int SearchMaxDist()
        {
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;
            int kMax;
            kMax = -(int)this.levelDist;
            kMin = (int)this.levelDist;

            st.Push(this);
            while (st.Count != 0)
            {
                current = st.Pop();
                if (current.levelDist > kMax)
                    kMax = (int)current.levelDist;

                if(current.levelDist<kMin)
                    kMin = (int)current.levelDist;

                if (current.joined != null)
                    foreach (var item in current.joined)
                        st.Push(item);
            }
            return kMax;
        }
        public Dictionary<double, List<HClusterNode>> GetClustersByLevels()
        {
            Dictionary<double, List<HClusterNode>> dic = new Dictionary<double, List<HClusterNode>>();

            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;
            st.Push(this);
            while (st.Count != 0)
            {
                current = st.Pop();
                if (!dic.ContainsKey(current.levelDist))
                    dic.Add(current.levelDist, new List<HClusterNode>());
                                    
                dic[current.levelDist].Add(current);

                if (current.joined != null)
                    foreach (var item in current.joined)
                        st.Push(item);
            }



            return dic;
        }
        public List<List<string>> GetClusters(int k)
        {
            List<List<string>> results = new List<List<string>>();
            Dictionary<HClusterNode,System.Drawing.Color> currentRes;
            Dictionary<HClusterNode, System.Drawing.Color> bestRes;
            Dictionary<HClusterNode, System.Drawing.Color> newRes;
            bool NotFinished = true;
            int min, max;
            max=SearchMaxDist();
            min = kMin;
            //max = kMax;

            bool direction = true;

            if (levelDist < joined[0].levelDist)
                direction = false;

            int threshold = min+(max - min) / 2;
            currentRes = CutDendrog(threshold);
            bestRes=currentRes;            
            if(currentRes.Count!=k)
                while (NotFinished)
                {
                    if (direction)
                    {
                        if (currentRes.Count < k)
                            max = threshold;
                        else
                            min = threshold;
                    }
                    else
                        if (currentRes.Count > k)
                            max = threshold;
                        else
                            min = threshold;

                    threshold = min+(max - min) / 2;
                    newRes = CutDendrog(threshold);
                    if (max-min>1 && newRes.Count!=k)             
                        NotFinished = true;
                    else
                        NotFinished = false;

                    if (Math.Abs(newRes.Count - k) < Math.Abs(bestRes.Count-k))                                                    
                            bestRes=newRes;
                        
                    currentRes = newRes;

                }

            foreach (var item in bestRes)            
                results.Add(item.Key.setStruct);

            return results;
        }
        public HClusterNode CutDendrogToLeaves(int leavesNumber)
        {
            HClusterNode rootNode= new HClusterNode(this);
            Stack<HClusterNode> org = new Stack<HClusterNode>();
            Stack<HClusterNode> copy = new Stack<HClusterNode>();
            HClusterNode current = null;
            HClusterNode xxx = null;
            Dictionary<HClusterNode, Color> dic=null;
            int max=SearchMaxDist();
            int min=kMin;
            int dist;

            while (max - min > 1)
            {
                dist = (max + min) / 2;
                dic = CutDendrog(dist);
                if (dic.Keys.Count > leavesNumber)
                    min = dist;
                else
                    max = dist;
            }

            org.Push(this);
            copy.Push(rootNode);
            while (org.Count != 0)
            {
                current = org.Pop();
                xxx = copy.Pop();
                xxx.joined = null;

                if (current.joined != null && !dic.ContainsKey(current))
                {
                    xxx.joined = new List<HClusterNode>();
                    foreach (var item in current.joined)
                    {
                        org.Push(item);
                        HClusterNode q = new HClusterNode(item);
                        if (q.joined != null)
                            q.joined.Clear();

                        if (xxx.joined == null)
                            xxx.joined = new List<HClusterNode>();
                        xxx.joined.Add(q);
                        copy.Push(q);
                    }
                }
            }            

            return rootNode;

        }
        public List<HClusterNode> FinPath(HClusterNode node)
        {
            List<HClusterNode> path = new List<HClusterNode>();
            HashSet<HClusterNode> toRoot = new HashSet<HClusterNode>();
            HClusterNode current = this;
            HClusterNode common = null;
            toRoot.Add(current);
            while (current.parent != null)
            {
                current = current.parent;
                toRoot.Add(current);
            }

            current = node;
            while (current.parent != null)
            {
                current = current.parent;
                if (toRoot.Contains(current))
                {
                    common = current;
                    break;
                }
            }
            if (common == null)
                return null;

            current = this;
            path.Add(current);
            while (current != common)
            {
                current = current.parent;
                path.Add(current);
            }
            List<HClusterNode> aux = new List<HClusterNode>();
            current = node;
            aux.Add(current);
            while (current != common)
            {
                current = current.parent;
                aux.Add(current);
            }
            for (int i = aux.Count - 2; i >= 0; i--)
                path.Add(aux[i]);
            return path;
        }
        public void MakeInvisible()
        {
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            st.Push(this);
            while (st.Count != 0)
            {
                var current = st.Pop();
                current.visible = false;
                if (current.joined != null)
                    foreach (var item in current.joined)
                        st.Push(item);
            }

        }
        public Dictionary<HClusterNode,System.Drawing.Color> CutDendrog(int distThreshold)
        {
            Dictionary<HClusterNode, System.Drawing.Color> returnList = new Dictionary<HClusterNode,System.Drawing.Color>();
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;

            st.Push(this);
            while (st.Count != 0)
            {
                current = st.Pop();
                if (current.levelDist <= distThreshold)
                {
                    returnList.Add(current, System.Drawing.Color.Red);
                    if (current.joined != null)
                        foreach (var item in current.joined)
                            if(item.levelDist==current.levelDist)
                                returnList.Add(item, System.Drawing.Color.Red);

                }
                else
                    if (current.joined != null)
                        foreach (var item in current.joined)
                            st.Push(item);

            }
            return returnList;
        }

        public void SaveNode(StreamWriter file,int num)
        {
            if (file != null)
            {
                file.WriteLine("================Cluster Num=" + num + "=======================");
                if (refStructure != null && refStructure.Length != 0)
                    file.WriteLine(refStructure);
                foreach (var item in setStruct)
                    if(item!=refStructure)
                        file.WriteLine(item);
            }
        }
        private void MakeStructures()
        {
            
            int counter=0;
            if (joined != null)
            {
                setStruct.Clear();
                foreach (var item in joined)
                    counter += item.setStruct.Count;
                setStruct = new List<string>(counter);

                foreach (var item in joined)
                    foreach (var str in item.setStruct)
                        setStruct.Add(str);
            }
        }
        public void RedoSetStructures()
        {
            if(joined==null)
                return;
            foreach (var item in joined)            
                   item.RedoSetStructures();                
                            
            this.MakeStructures();

        }

        public List<HClusterNode> GetLeaves(bool visibilityFlag=true)
        {
            List<HClusterNode> listH = new List<HClusterNode>();
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;
            st.Push(this);
            while (st.Count != 0)
            {
                current = st.Pop();

                if (current.joined != null)
                    foreach (var item in current.joined)
                        st.Push(item);
                else
                    if (visibilityFlag)
                {
                    if (current.visible)
                        listH.Add(current);
                }
                else
                    listH.Add(current);
            }

            return listH;

        }
        public bool IsVisible(HClusterNode node)
        {
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;

            if (this == node)
                return true;


            st.Push(this);
            while (st.Count != 0)
            {
                current = st.Pop();
                if (current.joined != null)
                    foreach (var item in current.joined)
                    {
                        if (item == node)
                            return true;
                        st.Push(item);
                    }
            }
            return false;

        }
	}

	public class hierarchicalCluster:IProgressBar
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


		public hierarchicalCluster (DistanceMeasure dMeasure,Options opt)
		{
			this.dMeasure=dMeasure;
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
            return "Agglomerative "+linkageType.ToString();
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
        /*private List<List<HClusterNode>> LevelMinimalDist(List<HClusterNode> levelNodes)
        {
            List<List<HClusterNode>> res = new List<List<HClusterNode>>();
            Dictionary<HClusterNode, int> test = new Dictionary<HClusterNode, int>();

            foreach (var item in levelNodes)
                if (!test.ContainsKey(item))
                    test.Add(item, 1);

            min = LMinimalDist(levelNodes);
            for (int i = 0; i < levelNodes.Count; i++)
            {                
                List<HClusterNode> lista = new List<HClusterNode>();
                for (int j = i + 1; j < levelNodes.Count; j++)
                {
                    int dist = dMeasure.FindMinimalDistance(levelNodes[i], levelNodes[j], linkageType).Key;
                    if (min == dist)
                    {
                        if (!levelNodes[j].fNode && !levelNodes[i].fNode)
                        {
                            lista.Add(levelNodes[i]);
                            lista.Add(levelNodes[j]);
                            levelNodes[i].fNode = true;
                            levelNodes[j].fNode = true;
                        }
                    }
                }
                if (lista.Count >= 2)
                    res.Add(lista);
            }

            foreach (var item in res)
                foreach (var it in item)
                    it.fNode = true;

            return res;
        }*/

        private List<List<HClusterNode>> LevelMinimalDist(List<HClusterNode> levelNodes)
        {
            List<List<HClusterNode>> res = new List<List<HClusterNode>>();
            byte[] e = new byte[levelNodes.Count];
            Dictionary<HClusterNode, int> test = new Dictionary<HClusterNode, int>();

            foreach (var item in levelNodes)
                if (!test.ContainsKey(item))
                    test.Add(item, 1);

            min = Int32.MaxValue;
           
            for (int i = 0; i < levelNodes.Count; i++)
            {
                List<HClusterNode> lista = new List<HClusterNode>();
              
                    for (int j = i + 1; j < levelNodes.Count; j++)
                    {
                        //if (levelNodes[j].fNode)
                        //    continue;
                        int dist = dMeasure.FindMinimalDistance(levelNodes[i], levelNodes[j], linkageType).Key;
                        if (min > dist)
                        {
                            min = dist;
                           foreach (var item in lista)
                                item.fNode = false;
                            lista.Clear();
                            foreach (var item in res)
                                foreach (var it in item)
                                    it.fNode = false;
                            res.Clear();

                            lista.Add(levelNodes[i]);
                            lista.Add(levelNodes[j]);
                            levelNodes[i].fNode = true;
                            levelNodes[j].fNode = true;
                          

                        }
                        else
                            if (min == dist)
                            {
                               if (!levelNodes[j].fNode && !levelNodes[i].fNode)
                                {
                                    lista.Add(levelNodes[i]);
                                    lista.Add(levelNodes[j]);
                                    levelNodes[i].fNode = true;
                                    levelNodes[j].fNode = true;
                                }


                            }
                    }
                if (lista.Count >= 2)
                    res.Add(lista);
            }

            foreach (var item in res)
                foreach (var it in item)
                    it.fNode = true;

            return res;
        }

        public ClusterOutput HierarchicalClusteringNNN(List<string> structures, List<Dictionary<double, int>[]> freq = null)
        {
            List<HClusterNode> levelNodes = new List<HClusterNode>();
            ClusterOutput outCl = new ClusterOutput();
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
            int levelCount = 0;
            dMeasure.CalcDistMatrix(structures);
            HClusterNode node;
            jury1D jury = new jury1D();
            jury.PrepareJury(((HammingBase)dMeasure).data);
            for (int i = 0; i < structures.Count; i++)
            {
                node = new HClusterNode();
                node.refStructure = structures[i];
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
            SortedDictionary<int, List<Tuple<int, int>>> sDic = new SortedDictionary<int, List<Tuple<int, int>>>();
            int dist = 0;
            for (int i = 0, k = 0; i < levelNodes.Count; i++)
            {
                for (int j = i + 1; j < levelNodes.Count; j++, k++)
                {
                    dist = dMeasure.GetDistance(dMeasure.hashIndex[levelNodes[i].refStructure], dMeasure.hashIndex[levelNodes[j].refStructure]);

                    if(!sDic.ContainsKey(dist))
                        sDic.Add(dist, new List<Tuple<int, int>>());
                    sDic[dist].Add(new Tuple<int, int>(i, j));
                }
            }
            while(sDic.Count > 0)
            {
                List<int> s = new List<int>(sDic.Keys);
                List<Tuple<int, int>> t = sDic[s[0]];
                for(int i=0;i<t.Count;i++)
                {
                    HashSet<int> idx = new HashSet<int>();
                    idx.Add(t[i].Item1);
                    idx.Add(t[i].Item2);
                    for(int j=i+1;j<t.Count;j++)
                    {
                        if (idx.Contains(t[j].Item1))
                            idx.Add(t[j].Item2);
                        else
                            if(idx.Contains(t[j].Item2))
                                idx.Add(t[j].Item1);
                    }

                }
            }



            return outCl;
        }
        public ClusterOutput HierarchicalClustering(List <string> structures,List<Dictionary<double,int>[]> freq=null)
		{
			List<List<HClusterNode>> level =new List<List<HClusterNode>>();
			List<HClusterNode> levelNodes=new List<HClusterNode>();		
			List<HClusterNode> rowNodes=new List<HClusterNode>();
            ClusterOutput outCl = new ClusterOutput();
            int levelCount = 0;
			bool end=false;
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

            jury1D jury = new jury1D();
            jury.PrepareJury(((HammingBase)dMeasure).data);
            for (int i=0;i<structures.Count;i++)
			{
				node=new HClusterNode();
				node.refStructure=structures[i];
				node.joined=null;
				node.setStruct.Add(structures[i]);
                node.levelNum = levelCount;
                if (freq != null && freq[i]!=null)
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
            maxV = levelNodes.Count+1;
			level.Add(levelNodes);
           
			while(!end)
			{                
				levelNodes=new List<HClusterNode>();
                List<List<HClusterNode>> rowList = LevelMinimalDist(level[level.Count - 1]);
                //List<List<HClusterNode>> rowList = LevelMinimalDistFast(level[level.Count - 1]);
                if (rowList.Count > 0)
                {
                    for(int s=0;s<rowList.Count;s++)
                    //foreach (var item in rowList)
                    {
						node=new HClusterNode();
                        
						node.joined=rowList[s];
                        foreach (var item in rowList[s])
                            item.parent = node;
                        node.levelDist = min;
                        node.realDist = dMeasure.GetRealValue(min);
                        node.levelNum = level.Count;
                        
						for(int m=0;m<rowList[s].Count;m++)
						{
							node.setStruct.AddRange(rowList[s][m].setStruct);
                            rowList[s][m].fNode = true;

							
						}
                  
                        List<KeyValuePair<string, double>> orderList = null;
                        if (rowList[s][0].stateFreq != null)
                        {
                            node.stateFreq = new Dictionary<double, int>[rowList[s][0].stateFreq.Length];
                            for (int m = 0; m < node.stateFreq.Length; m++)
                                node.stateFreq[m] = new Dictionary<double, int>();
                            //node.stateFreq=rowList[s][0].stateFreq;
                            for (int m = 0; m < rowList[s].Count; m++)
                            {
                                for (int n = 0; n < rowList[s][m].stateFreq.Length; n++)
                                    foreach (var state in rowList[s][m].stateFreq[n])
                                        if (node.stateFreq[n].ContainsKey(state.Key))
                                            node.stateFreq[n][state.Key] += state.Value;
                                        else
                                            node.stateFreq[n].Add(state.Key, state.Value);
                            }
                            ClusterOutput outJury = jury.JuryOptWeights(node.setStruct, node.stateFreq);
                            orderList = outJury.juryLike;

                            //jury1D jury=new jury1D(node.stateFreq,((HammingBase)dMeasure).al)
                        }
                        else
                            orderList = dMeasure.GetReferenceList(node.setStruct);
                        

                        node.refStructure = orderList[0].Key;
                        //node.refStructure = node.setStruct[0];

                       Dictionary<string, byte> refList = new Dictionary<string, byte>();
                        foreach (var itemJoined in node.joined)
                            if(refList.ContainsKey(itemJoined.refStructure))
                                refList.Add(itemJoined.refStructure,0);

                        node.refStructure = null;
                        if (mustRefStructure != null)
                            if (refList.ContainsKey(mustRefStructure))
                                node.refStructure = mustRefStructure;

                        if (node.refStructure == null)
                        {
                            if (orderList != null)
                            {
                                List<KeyValuePair<string, int>> listK = new List<KeyValuePair<string, int>>();
                                node.refStructure = node.joined[0].refStructure;
                                foreach (var it in orderList)
                                {
                                    if (refList.ContainsKey(it.Key))
                                    {
                                        node.refStructure = it.Key;
                                        break;
                                    }
                                }
                            }
      //                      node.refStructure = dMeasure.GetReferenceStructure(node.setStruct, refList);                           
                        }


                        levelNodes.Add(node);
					}
					
					
				}
				if(levelNodes.Count>0)
				{
					level.Add(levelNodes);
					for(int i=0;i<level[level.Count-2].Count;i++)
					{
						if(!level[level.Count-2][i].fNode)
							level[level.Count-1].Add(level[level.Count-2][i]);
					}
                    currentV = maxV - levelNodes.Count;

				}				
				if(level[level.Count-1].Count==1)
					end=true;
				
			}
            

            outCl.hNode = level[level.Count - 1][0];
            outCl.hNode.levelNum = 0;

            


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

