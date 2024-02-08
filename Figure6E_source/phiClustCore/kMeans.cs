using System;
using System.Collections;
using System.Collections.Specialized;
using System.Collections.Generic;
using phiClustCore;
using phiClustCore.Interface;
using phiClustCore.Distance;
using System.Data;

namespace phiClustCore
{
    public enum Initialization
    {
        RANDOM,
        Jury1D
    };
	public class kMeans:IProgressBar
	{

        int loopCounter=0;
		int avrIterNumber=0;
        public int maxDist { get; set; }
        public int maxK { get; set; }
        public int maxRepeat {get;set;}
        public int threshold { get; set; }
        public double BMIndex{get;set;}
        public Initialization initialization { get; set; }
		int all=0;
        string clusterName = "K-Means";
		DistanceMeasure dMeasure;
        int currentV = 0,hcurrentV=0;
        int maxV = 1,hmaxV=1;
        bool hierarchical = false;

        double startProgress = 0;
        double endProgress = 1;

        public double StartProgress { set { startProgress = value; } get { return startProgress; } }
        public double EndProgress { set { endProgress = value; } get { return endProgress; } }

        public kMeans(DistanceMeasure dMeasure,Initialization initialization,bool hierarchical=false)
		{
            this.dMeasure = dMeasure;
            this.initialization = initialization;
            this.hierarchical = hierarchical;
		}

        public kMeans(DistanceMeasure dMeasure,bool hierarchical=false)
        {
            this.dMeasure = dMeasure;
//            this.initialization = initialization;
            this.hierarchical = hierarchical;
        }

        public override string ToString()
        {
            return clusterName;
        }
        public double ProgressUpdate()
        {
            double res = dMeasure.ProgressUpdate();

            if(!hierarchical)
                return StartProgress + (EndProgress - StartProgress) * (res * 0.5 + 0.5 * (double)currentV / maxV);            


            return res* 0.5+0.5*(double)hcurrentV/hmaxV;
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }

		public ClusterOutput HierarchicalKMeans()
		{
			HClusterNode node;
            maxDist = 0;
            currentV = 0;
            List<string> availStruct = new List<string>(dMeasure.structNames.Keys);
            hmaxV =availStruct.Count;
			node=MakeNodes(availStruct,0);
            hcurrentV = hmaxV;
            node.levelDist = maxDist;
            node.realDist = dMeasure.GetRealValue(maxDist);
            AddDistance(node);

            ClusterOutput outClust = new ClusterOutput();
            outClust.hNode = node;

            clusterName = "H-Kmeans";

            return outClust;
		
        }
		public float CalculateDaviesBouldinIndex(List<List <string>>  clusters)
		{
            //ClusterOutput clustOut;
            List<string> refStructues = new List<string>(clusters.Count);
			int []dist;
			float []avr=new float[clusters.Count];
			float measure=0;			
			for(int i=0;i<clusters.Count;i++)
			{
				float sum=0;
				if(clusters[i].Count==0)
					continue;
				//clustOut=jury.JuryOpt(clusters[i]);
                string refStr = dMeasure.GetReferenceStructure(clusters[i]);
				dist=dMeasure.GetDistance(refStr,clusters[i]);
				foreach(var item in dist)
					sum+=item;
				
				avr[i]=sum/dist.Length;

                refStructues.Add(refStr);									
			}
			for(int i=0;i<refStructues.Count;i++)
			{
				float max=0;
				for(int j=0;j<refStructues.Count;j++)
				{
					int cDist;
					float v;
					if(i==j)
						continue;
					cDist=dMeasure.GetDistance(refStructues[i],refStructues[j]);
					v=((float)(avr[i]+avr[j]))/cDist;
					if(v>max)
						max=v;					
				}
				measure+=max;
			}
					
			return measure/refStructues.Count;
			
			
		}
		
		public Dictionary <string,List<string>> CalckMeans(List <string> refStruct,List <string> allStruct)
		{
			int [][] refDist=new int [refStruct.Count][];
			Dictionary <string,List<string>> clusters=new Dictionary<string, List<string>>();

            refDist = dMeasure.GetDistance(refStruct, allStruct);
			//for(int i=0;i<refStruct.Count;i++)		
//				refDist[i]=dMeasure.GetDistance(refStruct[i],allStruct);
						
			for(int j=0;j<allStruct.Count;j++)
			{
				int min=refDist[0][j];
				int remIndex=0;
				for(int i=0;i<refStruct.Count;i++)
					if(refDist[i][j]<min)
					{
						remIndex=i;
						min=refDist[i][j];
					}
				clusters[refStruct[remIndex]].Add(allStruct[j]);
			}
			
			return clusters;			
		}
        private List<string> SimpleRefStruct(int k, List<string> allStruct)
        {
            List<string> refStruct = new List<string>(k);
            Random rand=new Random();

            for (int i = 0; i < k;)
            {
                string rStruct=allStruct[rand.Next(allStruct.Count-1)];
                if (!refStruct.Contains(rStruct))
                {
                    refStruct.Add(rStruct);
                    i++;
                }
            }
            return refStruct;
        }
		private List<string> RefRandom(int k, List<string> allStruct)	
		{
			List <string> refStruct=new List<string>(k);
			Random rand=new Random(1);
			for(int i=0;i<k;i++)
			{
				do
				{
					string s=allStruct[rand.Next(allStruct.Count)];
					if(!refStruct.Contains(s))
					{
						refStruct.Add(s);					
						break;
					}
						
				}
				while(true);
			}
						
			return refStruct;
		}
        /*private List<string> RefStructFind(int k, List<string> allStruct)
		{
            ClusterOutput clustOut;
			List <string> refStruct=new List<string>();
			Random rand=new Random();
			int [][] refDist=new int[k][];
            clustOut = jury.JuryOpt(allStruct);
            int step = clustOut.juryLike.Count / k;

			
			for(int i=0;i<k;i++)
            	//refStruct.Add(clustOut.juryLike[rand.Next(0,clustOut.juryLike.Count/4)].Key);
				refStruct.Add(clustOut.juryLike[i].Key);
			
			return refStruct;

   			int [] sum=new int[allStruct.Count];
			
			
			for(int n=1;n<k;n++)
			{
				for(int i=0;i<sum.Length;i++)
					sum[i]=0;
				
				for(int i=0;i<refStruct.Count;i++)
				{
					refDist[i]=dMeasure.GetDistance(refStruct[i],allStruct);												
					
					for(int m=0;m<sum.Length;m++)
						sum[m]+=refDist[i][m];
				}
				List <KeyValuePair<string,int>> li=new List<KeyValuePair<string, int>>();
				KeyValuePair <string, int> v;

				for(int m=0;m<sum.Length;m++)
				{
					v=new KeyValuePair<string, int>(allStruct[m],sum[m]);
					li.Add(v);
				}
					
				li.Sort((firstPair,nextPair) =>
    			{
        			return nextPair.Value.CompareTo(firstPair.Value);
    			});
								
				
				refStruct.Add(li[rand.Next(0,sum.Length/3)].Key);				
			}
			
			return refStruct;
			
		}*/
        public ClusterOutput kMeansL(int k,int maxIter, List<string> allStruct)
        {
            ClusterOutput clustOut,remClust=new ClusterOutput();
            float cost,remCost=100;
            for (int r = 0; r < 1; r++)
            {
                clustOut = kMeansLevel(k, maxIter,allStruct);
                cost = CalculateDaviesBouldinIndex(clustOut.clusters.list);
                if (remCost > cost)
                {
                     remCost = cost;
                     remClust = clustOut;
                }                
            }

            return remClust;
        }
        private double STDDev(int[] dist)
        {
            double stDev = 0,avr=0;

            for (int i = 0; i < dist.Length; i++)
                avr += dist[i];
            avr /= dist.Length;
            for (int i = 0; i < dist.Length; i++)
                stDev += (dist[i] - avr) * (dist[i] - avr);

            stDev = Math.Sqrt(stDev / dist.Length);

            return stDev;
        }
        public List<string> RefStructFind2(int k, List<string> allStruct)
		{
			List <string> refStruct=new List<string>(k);			
			int [] refDist;
       
			List<string> actualList=new List<string>(allStruct);
							
			for(int j=0;j<k;j++)
			{
				int [] indexTab;
                refStruct.Add(dMeasure.GetReferenceStructure(actualList));
				if(j<k-1)
				{
					refDist=dMeasure.GetDistance(refStruct[refStruct.Count-1],actualList);
					indexTab=new int [refDist.Length];
					for(int i=0;i<indexTab.Length;i++)
						indexTab[i]=i;
				    
					Array.Sort(refDist,indexTab);
                    if (!dMeasure.order)
                        Array.Reverse(indexTab);
                    
					List<string> newList=new List<string>();
					for(int i=allStruct.Count/k;i<actualList.Count;i++)					
						newList.Add(actualList[indexTab[i]]);
					


					actualList=newList;
				}
			}
			
			
			return refStruct;
			
		}
        public List<string> RefStructFind3(List<string> allStruct)
        {
            List<string> refStruct = new List<string>();
            List <KeyValuePair<string,double>> refS = dMeasure.GetReferenceList(allStruct);

            refStruct.Add(refS[0].Key);
            refStruct.Add(refS[refS.Count-1].Key);
            
            return refStruct;

        }
        public ClusterOutput kMeansRun(int iterNum, List<string> allStruct, List<string> refStruct)
        {
            ClusterOutput clustOut;
            int currentBest = allStruct.Count;
            int bestCounter = 0 ;
            bool end = false;
            int index = 0, changeCounter=1000;
            int[][] refDist = new int[refStruct.Count][];

            List<List<string>> clusters = new List<List<string>>(refStruct.Count);
            int[] status = new int[allStruct.Count];
            Random rand = new Random();

            if (refStruct == null || refStruct.Count<=1)
                return null;

            for (int j = 0; j < refStruct.Count; j++)            
                clusters.Add(new List<string>());
            


            //refStruct = SimpleRefStruct(k, allStruct);
            for (int i = 0; i < status.Length; i++)
                status[i] = -1;
            loopCounter = 0;
            while (!end)
            {
                for (int i = 0; i < refStruct.Count; i++)
                {
                    clusters[i].Add(refStruct[i]);
                }
                refDist = dMeasure.GetDistance(refStruct, allStruct);
                //for(int i=0;i<refStruct.Count;i++)
                //					refDist[i]=dMeasure.GetDistance(refStruct[i],allStruct);
                changeCounter = 0;
                for (int i = 0; i < allStruct.Count; i++)
                {
                    int min = refDist[0][i];
                    index = 0;
                    for (int j = 1; j < refStruct.Count; j++)
                    {
                        if (refDist[j][i] < min)
                        {
                            min = refDist[j][i];
                            index = j;
                        }
                    }
                    if (allStruct[i].Contains(clusters[index][0]))
                        continue;

                    clusters[index].Add(allStruct[i]);

                    if (status[i] != index)
                    {
                        status[i] = index;
                        changeCounter++;
                    }
                  
                }
                if (currentBest <= changeCounter)
                    bestCounter++;
                else
                {
                    bestCounter = 0;
                    currentBest = changeCounter;
                }
                if (changeCounter <= allStruct.Count / 100 || loopCounter > iterNum || bestCounter>=5)
                    end = true;
                else
                {
                    loopCounter++;
                    for (int i = 0; i < refStruct.Count; i++)
                    {
                        if (clusters[i].Count > 0)
                            refStruct[i] = dMeasure.GetReferenceStructure(clusters[i]);
                        else
                            refStruct[i] = allStruct[rand.Next(0, allStruct.Count)];
                        clusters[i].Clear();                
                    }
                    currentV++;
                }

            }
            List<List<string>> finalClusters = new List<List<string>>();
            for (int i = 0; i < clusters.Count; i++)
            {
                List<string> aux = new List<string>();
                aux.Add(refStruct[i]);
                foreach (var item in clusters[i])
                    if (item != refStruct[i])
                        aux.Add(item);
                if(aux.Count>=1)
                    finalClusters.Add(aux);
            }

            clustOut = new ClusterOutput();
            clustOut.clusters = new clusterRes();
            clustOut.clusters.list = finalClusters;
           //     clustOut.clusters = clusters;
            currentV = maxV;
            return clustOut;

        }
		public ClusterOutput kMeansLevel(int k, int maxIter, List <string> allStruct)
		{
			List <string> refStruct=new List<string>();
						
            if(initialization==Initialization.Jury1D)
                if(k==2)
                    refStruct = RefStructFind3(allStruct);
                else
                    refStruct = RefStructFind2(k, allStruct);
            else
			    refStruct=RefRandom(k,allStruct);

            maxV = maxIter+1;

            return kMeansRun(maxIter,allStruct, refStruct);

            
		}


		private KeyValuePair<float,List<List<string>>> FindClustersNumber(List<string> structures)
		{
            ClusterOutput clustOut;
			List<List <string>>   remClust=null;
			float cost=0,remCost=0;
			
			remCost=float.MaxValue;
			for(int k=2;k<=maxK;k++)
			{
				maxRepeat=1;
				for(int r=0;r<maxRepeat;r++)
				{
					bool test=false;
					clustOut=kMeansLevel(k,30,structures);
					for(int n=0;n<clustOut.clusters.list.Count;n++)
                        if (clustOut.clusters.list[n].Count < 10)
							test=true;

					if(!test)
					{
                        cost = CalculateDaviesBouldinIndex(clustOut.clusters.list);
						if(remCost>cost)
						{
							remCost=cost;
							remClust=clustOut.clusters.list;
						}	
					}
				}
			}
			return new KeyValuePair<float,List<List<string>>> (remCost,remClust);			
		}

        private void AddDistance(HClusterNode node)
        {
            if (node.joined == null)
                return;
            for (int i = 0; i < node.joined.Count; i++)
            {
                node.joined[i].levelDist = node.levelDist - 1;
                node.joined[i].realDist = dMeasure.GetRealValue(node.joined[i].levelDist);
                if (node.joined[i].joined != null && node.joined[i].joined.Count > 0)
                    AddDistance(node.joined[i]);
            }
        }
        private HClusterNode MakeNodes(List<string> sNames,int levelNum)
		{			
			HClusterNode node=new HClusterNode();
			HClusterNode nodeInside;
			
			List<List <string>> clusters=null;

			
			if(sNames.Count==0)
				return null;

            node.setStruct = sNames;
            node.iterNum = loopCounter;
            node.num = ++levelNum;
            if (node.num > maxDist)
                maxDist = node.num;

			KeyValuePair<float,List<List <string>>> cl=FindClustersNumber(sNames);
			if(sNames.Count==dMeasure.structNames.Count)
				clusters=cl.Value;
			else
				if(cl.Key<BMIndex)
					clusters=cl.Value;
				else
					clusters=null;
            if(clusters==null)
                hcurrentV+=sNames.Count;

			if(clusters!=null && clusters.Count>0)
			{
				for(int i=0;i<clusters.Count;i++)
				{
                    if (clusters[i].Count > threshold)
                    {
                        nodeInside = MakeNodes(clusters[i], levelNum);
                        if (nodeInside != null)
                        {
                            if (node.joined == null)
                                node.joined = new List<HClusterNode>();
                            node.joined.Add(nodeInside);
                        }
                    }
                    else
                        hcurrentV += clusters[i].Count;
							
				}
				avrIterNumber+=levelNum;
				all++;
			}

                
			return node;
		}		
		
	}
	
	
}

