using System;
using System.Collections.Generic;
using System.Data;
using System.Threading;
using phiClustCore.Interface;
using phiClustCore;
using System.IO;
using System.Linq;

namespace phiClustCore.Distance
{
    
    public abstract class DistanceMeasure:IProgressBar
    {
        protected int[] distanceMatrix = null;
        public Dictionary<string, int> hashIndex;// = new Dictionary<string, int>();
        [NonSerialized]
        protected ManualResetEvent[] resetEvents = null;
        protected Dictionary<string, Dictionary<string, float [,]>> rotStruct = new Dictionary<string, Dictionary<string, float [,]>>();
        public Dictionary<string,int> structNames;
        protected List<string> structures;
        public bool order; //false if best are biggest values otherwise like in regular distance true.
        public double maxSimilarity;
        private const int ErrorValue = int.MaxValue;
        public int errorValue { get { return ErrorValue; } }
        public string dirName;
        public string alignFile;
        protected List<KeyValuePair<int, int>>[] indexList;
        protected jury1D jury = null;
        protected string juryProfile = null;

        protected int maxV=1;
        protected int currentV;

        double startProgress = 0;
        double endProgress = 1;

        public double StartProgress { set { startProgress = value; } get { return startProgress; } }
        public double EndProgress { set { endProgress = value; } get { return endProgress; } }


        public DistanceMeasure()
        {

        }        
        public virtual void InitMeasure(Dictionary<string,List<double>> data, bool referenceFlag)
        {
            if (referenceFlag)
            {                
                jury = new jury1D();
                jury.PrepareJury(data);
            }
        }
        public abstract void InitMeasure();
        public virtual int GetDistance(string m, string n){ return 0;}
        public virtual double GetRealValue(double v){return v/100; }
        public virtual bool SimilarityThreshold(float threshold, float dist){return true;}


        public virtual double ProgressUpdate()
        {
            return StartProgress + (EndProgress - StartProgress) * (double)currentV / maxV;
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }

        protected int FindIndex(int matIndex)
        {
            return (2*hashIndex.Count - matIndex + 1) * matIndex / 2;
        }
        public virtual int GetDistance(float[,] refStructure, string modelStructure, out int[] index) {index = null; return 0; }
        //public virtual boo
        public virtual string GetReferenceStructure(List<string> structures)
        {
            ClusterOutput clustOut;

            clustOut = jury.JuryOptWeights(structures);
            if (clustOut == null)
                return null;
            //clustOut = jury.JuryOpt(structures);

            return clustOut.juryLike[0].Key;

        }

        public string GetReferenceStructure(List<string> structures, List<string> refStruct)
        {
            List<KeyValuePair<string, double>> refList = null;
            refList = GetReferenceList(structures);
           // return refList[0].Key;

            if (refStruct != null)
            {
                List<KeyValuePair<string, int>> listK = new List<KeyValuePair<string, int>>();
                foreach (var item in refList)
                {

                   // int remIndex = 0;
                    for (int j = 0; j < refStruct.Count; j++)
                        if (refStruct[j] == item.Key)
                        {
                            return item.Key;
                            //remIndex = j;
                            //break;
                        }
                    /*KeyValuePair<string,int> aux = new KeyValuePair<string, int>(item, remIndex);
                    listK.Add(aux);
                    listK.Sort((firstPair, nextPair) =>
                    {
                        return firstPair.Value.CompareTo(nextPair.Value);
                    });
                    return listK[0].Key;*/
                }

            }
            //return refList[0].Key;
            return refList[0].Key;
        }

        public virtual List<string> GetReferenceStructureList(List<string> structures)
        {
            List<string> refStruct = new List<string>();
            ClusterOutput op;
            if (jury != null)
            {
                //op = jury.JuryOpt(structures);
                op = jury.JuryOptWeights(structures);
                if (op == null)
                    return null;
                if (op.juryLike[0].Value - op.juryLike[op.juryLike.Count - 1].Value > 0.1)
                {
                    refStruct.Add(op.juryLike[0].Key);
                    refStruct.Add(op.juryLike[op.juryLike.Count - 1].Key);
                    return refStruct;
                }
            }
            return null;
        }
        protected virtual void CalcMatrix(object o)
        {
            int num = (int)o;
            foreach (var item in indexList[num])
            {
                int dd = FindIndex(item.Key);

                int distC = GetDistance(structures[item.Key], structures[item.Value]);
                if (distC < 0 || distC >= int.MaxValue)
                    distC = -1;

                distanceMatrix[dd + item.Value - item.Key] = distC;
                Interlocked.Increment(ref currentV);
            }
            resetEvents[num].Set();
        }


        public virtual void CalcDistMatrix(List<string> structures) 
        {
            currentV = 0;

            hashIndex = new Dictionary<string, int>(structures.Count);
            for (int i = 0; i < structures.Count; i++)               
                hashIndex.Add(structures[i], i);
            
            Settings set = new Settings();
            set.Load();
            int threadNumbers = set.numberOfCores;
            this.structures = new List<string>(structures);

            int part = structures.Count * (structures.Count + 1) / (2 * threadNumbers) + 1;          

            distanceMatrix = new int[structures.Count * (structures.Count + 1) / 2];
            indexList = new List<KeyValuePair<int, int>>[threadNumbers];
          
            for (int n = 0; n < threadNumbers; n++)            
                indexList[n] = new List<KeyValuePair<int, int>>(part);

            int count = 0;
            foreach (var dim1 in structures)
            {
                foreach (var dim2 in structures)
                {
                    if (hashIndex[dim1] <= hashIndex[dim2])
                    {
                        indexList[count].Add(new KeyValuePair<int, int>(hashIndex[dim1], hashIndex[dim2]));
                        if (part <= indexList[count].Count)
                            count++;
                    }
                }
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
        public void SaveDistancedToFile(string fileName)
        {
            StreamWriter w = new StreamWriter(fileName);
            foreach (var dim1 in structures)
            {
                foreach (var dim2 in structures)
                {
                    if (hashIndex[dim1] <= hashIndex[dim2])
                    {
                        w.WriteLine(dim1 + " " + dim2 + " " + GetDistance(hashIndex[dim1], hashIndex[dim2]));
                    }
                }
            }

            w.Close();
        }
        public virtual List<KeyValuePair<string, double>> GetReferenceList(List<string> structures)
        {

            ClusterOutput clustOut;
            if(structures.Count<=2)
            {
                List<KeyValuePair<string, double>> cc = new List<KeyValuePair<string, double>>();
                foreach(var item in structures)
                    cc.Add(new KeyValuePair<string, double>(item, 0));                
                return cc;
            }
            clustOut = jury.JuryOptWeights(structures);
            if (clustOut == null)
                return null;
            return clustOut.juryLike;

        }

        public virtual float[,] GetReferenceStructure(string refStruct,List <string> structures)
        {
            float[,] refPosition = new float[rotStruct[refStruct][structures[0]].GetLength(0), 3];
            int[] count = new int[rotStruct[refStruct][structures[0]].GetLength(0)];
            for (int i = 0; i < structures.Count;i++ )
            {
                float[,] aux = rotStruct[refStruct][structures[i]];
                for (int n = 0; n <aux.GetLength(0); n++)
                {
                    if (aux[n, 0] != float.MaxValue)
                    {                        
                        for (int j = 0; j < 3; j++)                        
                            refPosition[n, j] += aux[n, j];
            
                        count[n]++;
                    }
                }
            }
            for (int i = 0; i < refPosition.GetLength(0); i++)
                for (int j = 0; j < refPosition.GetLength(1); j++)
                    if (count[i] > 0 && refPosition[i,j]!=float.MaxValue)
                        refPosition[i, j] /= count[i];

            return refPosition;
        }
      public virtual int[][] GetDistance(List<string> refStructure, List<string> structures)
      {
          int[][] dist = new int[refStructure.Count][];
          rotStruct.Clear();
          for (int j = 0; j < refStructure.Count; j++)
          {
              dist[j] = new int[structures.Count];

              //foreach (string item in structures)
                  dist[j]= GetDistance(refStructure[j], structures);

          }

          return dist;

      }
       /* public virtual int[][] GetDistance(List<string> refStructure, List<string> structures)
        {
            int[][] dist = new int[refStructure.Count][];
            int i = 0;
			if(opt!=null)
				rotStruct.Clear();
            for (int j = 0; j < refStructure.Count; j++)
            {
                dist[j] = new int[structures.Count];
                i = 0;
				
                foreach (string item in structures)				
                    dist[j][i++] = GetDistance(refStructure[j], item);
									
            }

            return dist;

        }*/
        //private void int[] 
        /*public virtual int[] GetDistance(string refStructure, List<string> structures)
        {
            int[] dist = new int[structures.Count];
            int i = 0;
            foreach (string item in this.structures)
                dist[i++] = GetDistance(refStructure, item);

            return dist;

        }*/

        private void CalcMatrixList(object o)
        {
            int num = (int)o;

            for (int i = 0; i < indexList[num].Count; i++)
                distanceMatrix[indexList[num][i].Key] = GetDistance(structures[structures.Count - 1], structures[indexList[num][i].Key]);

            resetEvents[num].Set();
        }
        public virtual int[] GetDistance(string refStructure, List<string> structures)
        {
            int[] dist = new int[structures.Count];
            distanceMatrix = new int[structures.Count];
            Settings set = new Settings();
            set.Load();
            int threadNumbers = set.numberOfCores;
            this.structures = new List<string>(structures);
            double part = (double)this.structNames.Count / threadNumbers;
            indexList = new List<KeyValuePair<int, int>>[threadNumbers];

            for (int j = 0; j < indexList.Length; j++)
                indexList[j] = new List<KeyValuePair<int, int>>();

            int count = 0;
            for (int j = 0; j < structures.Count; j++)
            {
                if (indexList[count].Count < part*(j+1)-part*j)
                {
                    KeyValuePair<int, int> aux = new KeyValuePair<int, int>(j, 0);
                    indexList[count].Add(aux);
                }
                else
                    count++;
            }
            this.structures.Add(refStructure);
            resetEvents = new ManualResetEvent[threadNumbers];
            for (int n = 0; n < threadNumbers; n++)
            {
                int k = n;
                resetEvents[n] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(new WaitCallback(CalcMatrixList), (object)k);
                
            }
            for (int j = 0; j < threadNumbers; j++)
                resetEvents[j].WaitOne();


            for (int j = 0; j < this.structures.Count - 1; j++)
                dist[j] = distanceMatrix[j];

            return dist;

        }
		public virtual int[] GetDistance(float [,] refStructure, List<string> structures)
        {
            int[] dist = new int[structures.Count];
            int i = 0;
            foreach (string item in structures)
            {
                int[] index;
                dist[i++] = GetDistance(refStructure, item,out index);
            }

            return dist;

        }
        public int GetIndex(string i,string j)
        {
            return GetIndex(hashIndex[i], hashIndex[j]);
        }
		public int GetIndex(int i,int j)
        {
            int m;

            if (i < j)
            {
                return FindIndex(i) + j - i;
            }
            else
            {
                return FindIndex(j) + i - j;
            }

        }
		
        public int GetDistance(int i, int j)
        {
            if (distanceMatrix == null)
                return -1;

            return distanceMatrix[GetIndex(i, j)];
        }

        public int GetMaxDistance()
        {
            if (distanceMatrix == null)
                return -1;

            int max = distanceMatrix[0];

            for (int i = 1; i < distanceMatrix.GetLength(0); i++)
                    if (max < distanceMatrix[i])
                        max = distanceMatrix[i];

            return max;
        }
        public List<Tuple<int, HClusterNode,HClusterNode>> FindMinimalDistanceReference(List<HClusterNode> listNodes)
        {
            int dist = 0;
            Tuple<int, HClusterNode,HClusterNode>[] distArray = new Tuple<int, HClusterNode,HClusterNode>[listNodes.Count*(listNodes.Count+1)/2-listNodes.Count];
            for (int i = 0,k=0; i < listNodes.Count; i++)
            {
                for (int j = i + 1; j < listNodes.Count; j++,k++)
                {
                    dist = GetDistance(hashIndex[listNodes[i].refStructure], hashIndex[listNodes[j].refStructure]);
                    distArray[k] = new Tuple<int, HClusterNode, HClusterNode>(dist, listNodes[i], listNodes[j]);
                }
            }

            return distArray.OrderBy(x => x.Item1).ToList();
        }

        public KeyValuePair<int,int> FindMinimalDistance(HClusterNode clust1, HClusterNode clust2, AglomerativeType linkageType)
        {
            int dist=0;
            int index=0;
            int min = 0;
            switch (linkageType)
            {
                case AglomerativeType.SINGLE:
                    min = Int32.MaxValue;
                    foreach (var item1 in clust1.setStruct)
                    {
                        int v = hashIndex[item1];
                        foreach (var item2 in clust2.setStruct)
                        {
                            dist = GetDistance(v, hashIndex[item2]);
                            if (dist < min)
                            {
                                min = dist;
                                index = GetIndex(v, hashIndex[item2]);
                            }
                        }
                    }
                    break;
                case AglomerativeType.AVERAGE:
                    dist = GetDistance(hashIndex[clust1.refStructure], hashIndex[clust2.refStructure]);

                    index = GetIndex(hashIndex[clust1.refStructure], hashIndex[clust2.refStructure]);
                    break;
                case AglomerativeType.COMPLETE:
                    min = Int32.MinValue;
                    foreach (var item1 in clust1.setStruct)
                    {
                        int v = hashIndex[item1];
                        foreach (var item2 in clust2.setStruct)
                        {
                            dist = GetDistance(v, hashIndex[item2]);
                            if (dist > min)
                            {
                                min = dist;
                                index = GetIndex(v, hashIndex[item2]);
                            }
                        }
                    }
                    break;

            }
            return new KeyValuePair<int, int>(dist, index) ;
        }

    }

}
