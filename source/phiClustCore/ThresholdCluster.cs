using System;

using System.Collections;
using System.Collections.Generic;
using phiClustCore;
using phiClustCore.Interface;
using System.Data;
using phiClustCore.Distance;

namespace phiClustCore
{
	public class ThresholdCluster:IProgressBar
	{
		public DistanceMeasure dMeasure;
		bool [] pointMark;
		float threshold=5;
        int minCluster = 30;
        int progressRead = 0;
        int currentV = 0;
        int maxV = 1;

        double startProgress = 0;
        double endProgress = 1;

        public double StartProgress { set { startProgress = value; } get { return startProgress; } }
        public double EndProgress { set { endProgress = value; } get { return endProgress; } }

        public ThresholdCluster(DistanceMeasure dMeasure,float threshold,int minCluster)
        {
            this.dMeasure = dMeasure;
            this.threshold = threshold*100;
            this.minCluster = minCluster;
        

            //dist = cBase.CalcRmsd(cBase.structNames);
        }
        public double ProgressUpdate()
        {
            double sumProgress = 0;
            double progress = dMeasure.ProgressUpdate();

            if (progressRead == 1)
                sumProgress = 0.05 + progress * 0.7;
            else
                sumProgress = 0.05 * progress;

            return StartProgress+(EndProgress - StartProgress) * (sumProgress + 0.25 * ((double)currentV / maxV));
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }


        private List<string> CreateCluster(int index)
		{
			List <string> items=new List<string>();
            List<string> ii = new List<string>(dMeasure.structNames.Keys);
			for(int j=0;j<pointMark.Length;j++)
				if(!pointMark[j] && dMeasure.GetDistance(index,j)<threshold)
				{
					items.Add(ii[j]);
                    pointMark[j] = true;
				}			
			
			return items;
		}
		public ClusterOutput OrgClustering()
		{
			
			int []count;
			int []index;
			bool end;
            ClusterOutput output = new ClusterOutput();            
			List <List<string>> clusters=new List<List<string>>();
			List <string> items;

            pointMark = new bool[dMeasure.structNames.Count];
            for (int i = 0; i < pointMark.Length; i++)
                    pointMark[i] = false;


            progressRead = 1;
            dMeasure.CalcDistMatrix(new List <string>(dMeasure.structNames.Keys));

            maxV = dMeasure.structNames.Count;
			count=new int[dMeasure.structNames.Count];
            index = new int[dMeasure.structNames.Count];
			end=false;
			while(!end)
			{
			
				for(int i=0;i<pointMark.Length;i++)
				{
                    count[i] = 0;
                    if (pointMark[i])
                        continue;
					index[i]=i;					
					for(int j=0;j<pointMark.Length;j++)
						if(!pointMark[j] && dMeasure.GetDistance(i,j)<threshold)
							count[i]++;
				}
			
				Array.Sort<int>(index, (a,b) => count[b].CompareTo(count[a]));
                if (count[index[0]] < minCluster)
                {
                    end = true;
                    break;
                }
				items=CreateCluster(index[0]);
				if(items.Count>minCluster)
					clusters.Add(items);				
				else
					end=true;

                currentV += items.Count;
			}
            output.clusters = new clusterRes();
            output.clusters.list = clusters;
            currentV = maxV;
			return output;			
		}
		
/*		public ClusterOutput JuryBaker()
		{
			ClusterOutput juryRes;
			int index=0;
			List <List<string>> clusters=new List<List<string>>();			
			List <string> items;
			bool end;
			jury1D jury=new jury1D(dMeasure);				
			juryRes=jury.JuryOpt(dMeasure.structNames);
			
			end=false;
			while(!end)
			{
				for(int i=0;i<juryRes.juryLike.Count;i++)
				{
					if(juryRes.juryLike[i].Value>-1)
					{
						index=dMeasure.hashIndex[juryRes.juryLike[i].Key];
						juryRes.juryLike.Remove(juryRes.juryLike[i]);						
						break;
					}
				}
				items=CreateCluster(index);
				if(items.Count>0)
				{
					clusters.Add(items);				
					foreach(string name in items)
						for(int i=0;i<juryRes.juryLike.Count;i++)
							if(juryRes.juryLike[i].Key==name)
							{
								juryRes.juryLike.Remove(juryRes.juryLike[i]);
								break;
							}
				}
				else
					end=true;				
			}

            ClusterOutput output = new ClusterOutput();

            output.clusters = clusters;

			return output;
			
		}*/
			
	}
}

