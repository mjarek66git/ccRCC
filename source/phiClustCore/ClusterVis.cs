using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace phiClustCore
{
    public class ClusterVis
    {
        public ClusterOutput output;

        private StreamWriter wr;

        public ClusterVis() { }

        public ClusterVis(ClusterOutput output)
        {
            this.output = output;
        }

        public int SearchKmax(HClusterNode hNode)
        {
            int kMax = hNode.levelNum;
            if (hNode.joined == null)
                return kMax;
            foreach (var item in hNode.joined)
            {
                if (SearchKmax(item) > kMax)
                    kMax = item.levelNum;
            }
            return kMax;
        }
        public void SaveClusters(string fileName)
        {
            SaveClusters(output.clusters.list, fileName);
        }       
        public void SaveClusters(List<List<string>> clust,string fileName)
        {
            if (clust != null)
            {
                wr = new StreamWriter(fileName);
                int num = 1;
                if (clust != null)
                    foreach (var item in clust)
                    {
                        wr.WriteLine("========== Cluster Num=" + num + " ==========");
                        foreach (var cl in item)
                            wr.WriteLine(cl);
                        num++;
                    }
                wr.Close();
            }


        }
        public void SaveJuryLike(string fileName)
        {
            if (output.juryLike != null)
            {
                wr = new StreamWriter(fileName);
                wr.WriteLine("========== Structures order ==========");
                if(output.juryLike!=null)
                    foreach (var item in output.juryLike)
                        wr.WriteLine(item.Key + " " + item.Value);

                wr.Close();
            }
        }
		public void SaveLeafs(string fileName)
		{
            if (output.hNode != null)
            {
                wr = new StreamWriter(fileName);
				Queue<HClusterNode> queue = new Queue<HClusterNode>();
				HClusterNode hNode;
				queue.Enqueue(output.hNode);
                
                while (queue.Count != 0)
                {
                	hNode = queue.Dequeue();
	                if(hNode.setStruct!=null)// && hNode.joined==null)
					{
						wr.WriteLine("========== CLUSTER ==========");					
//						wr.WriteLine(hNode.refStructure);
                	    foreach (var item in hNode.setStruct)
        	        	    wr.WriteLine(item);
					}

                	if(hNode.joined!=null)
            	    	foreach (var item in hNode.joined)
                    		queue.Enqueue(item);
                
				}
				
				wr.Close();
			}			
		}
        public void SaveHierarchical(string fileName)
        {
            if (output.hNode != null)
            {
                wr = new StreamWriter(fileName);
                int num = 0;
                Queue<HClusterNode> queue = new Queue<HClusterNode>();
                HClusterNode hNode;
                queue.Enqueue(output.hNode);
                num = output.hNode.levelNum;
                wr.WriteLine("========== LEVEL " + num + " ==========");
                while (queue.Count != 0)
                {
                    hNode = queue.Dequeue();
                    if (num != hNode.levelNum)
                    {
                        num = hNode.levelNum;
                        wr.WriteLine("========== LEVEL " + num + " ==========");                        
                    }

                    wr.WriteLine("========== CLUSTER ==========");
                    if(hNode.setStruct!=null)
					{
						wr.WriteLine("Reference: "+hNode.refStructure);
                        foreach (var item in hNode.setStruct)
                            wr.WriteLine(item);
					}

                    if(hNode.joined!=null)
                        foreach (var item in hNode.joined)
                            queue.Enqueue(item);
                }
                wr.Close();
            }
        }
        public virtual void SCluster(string fileName)
        {
            if (output.hNode!=null)
            {                
        //        SaveHierarchical(fileName+"_Hnode.cl");
                List<List<string>> clust=output.hNode.GetClusters(10);
                SaveClusters(clust, fileName);
				//SaveLeafs(fileName+"_leaves.cl");
            }
            if (output.juryLike != null)
            {
                SaveJuryLike(fileName);   
            }
            if (output.clusters!=null)
            {
                SaveClusters(fileName);
            }
        }
    }
}
