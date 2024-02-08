using phiClustCore.Distance;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace phiClustCore
{
	public struct Edge
	{
		public int Source;
		public int Destination;
		public int Weight;
	}

	public struct GraphMST
	{
		public int VerticesCount;
		//public int EdgesCount;
		public Edge[] edge;
	}

	public struct Subset
	{
		public int Parent;
		public int Rank;
	}

	public class MST
    {
		List<List<string>> clusters = null;
		DistanceMeasure dist = null;
		int Find(Subset[] subsets, int i)
		{
			if (subsets[i].Parent != i)
				subsets[i].Parent = Find(subsets, subsets[i].Parent);

			return subsets[i].Parent;
		}

		void Union(Subset[] subsets, int x, int y)
		{
			int xroot = Find(subsets, x);
			int yroot = Find(subsets, y);

			if (subsets[xroot].Rank < subsets[yroot].Rank)
				subsets[xroot].Parent = yroot;
			else if (subsets[xroot].Rank > subsets[yroot].Rank)
				subsets[yroot].Parent = xroot;
			else
			{
				subsets[yroot].Parent = xroot;
				++subsets[xroot].Rank;
			}
		}
		public GraphMST CreateGraph(DistanceMeasure dist, List<HClusterNode> clusters)
        {
			List<List<string>> cl = new List<List<string>>();
			foreach(var item in clusters)
            {
				cl.Add(item.setStruct);
            }
			return CreateGraph(dist, cl);
        }
		
			public GraphMST CreateGraph(DistanceMeasure dist,List<List<string>> clusters)
        {
			this.clusters = clusters;
			this.dist = dist;
			GraphMST g = new GraphMST();

			g.VerticesCount = clusters.Count;
			g.edge = new Edge[g.VerticesCount * (g.VerticesCount - 1) / 2];

			for(int i=0,n=0;i<clusters.Count;i++)
            {
				for(int j=i+1;j<clusters.Count;j++,n++)
                {
					int di = dist.GetDistance(clusters[i][0], clusters[j][0]);
					g.edge[n].Weight = di;
					g.edge[n].Source = i;
					g.edge[n].Destination = j;
                }
            }

			return g;
        }
		void Print(Edge[] result, int e)
		{
			for (int i = 0; i < e; ++i)
				Console.WriteLine("{0} -- {1} == {2}", result[i].Source, result[i].Destination, result[i].Weight);
		}
		public ClusterOutput Closest(HClusterNode start)
        {
			ClusterOutput res = new ClusterOutput();
			List<string> aux = clusters[0];

			int bestStart = -1;
			for (int i = 0; i < clusters.Count; i++)
			{
				if (start.setStruct[0] == clusters[i][0])
					bestStart = i;
			}

			clusters[0] = clusters[bestStart];
			clusters[bestStart] = aux;
			res.clusters = new clusterRes();
			res.clusters.list = new List<List<string>>();
			res.clusters.list.Add(clusters[0]);
			for (int i = 0; i < clusters.Count; i++)
			{
				int min = int.MaxValue;
				int index = -1;
				for (int j = i + 1; j < clusters.Count; j++)
				{
					int di = dist.GetDistance(clusters[i][0], clusters[j][0]);
					if(di<min)
                    {
						min = di;
						index = j;
                    }
				}
				
				if (i + 1 < clusters.Count)
				{
					res.clusters.list.Add(clusters[index]);
					aux = clusters[i + 1];
					clusters[i + 1] = clusters[index];
					clusters[index] = aux;
				}
			}
			return res;
		}
		public List<List<string>> Kruskal(GraphMST graph,HClusterNode startNode)
		{
			int verticesCount = graph.VerticesCount;
			Edge[] result = new Edge[verticesCount];
			int i = 0;
			int e = 0;

			Array.Sort(graph.edge, delegate (Edge a, Edge b)
			{
				return a.Weight.CompareTo(b.Weight);
			});

			Subset[] subsets = new Subset[verticesCount];

			for (int v = 0; v < verticesCount; ++v)
			{
				subsets[v].Parent = v;
				subsets[v].Rank = 0;
			}

			while (e < verticesCount - 1)
			{
				Edge nextEdge = graph.edge[i++];
				int x = Find(subsets, nextEdge.Source);
				int y = Find(subsets, nextEdge.Destination);

				if (x != y)
				{
					result[e++] = nextEdge;
					Union(subsets, x, y);
				}
			}
			Dictionary<string, int> edgesDic = new Dictionary<string, int>();
			Dictionary<int, HashSet<int>> sourceDest = new Dictionary<int, HashSet<int>>();
			HashSet<int> used = new HashSet<int>();
			for(i=0;i<result.Length;i++)
            {
				if (!sourceDest.ContainsKey(result[i].Source))
					sourceDest.Add(result[i].Source, new HashSet<int>());
				if (!sourceDest.ContainsKey(result[i].Destination))
					sourceDest.Add(result[i].Destination, new HashSet<int>());
				sourceDest[result[i].Source].Add(result[i].Destination);
				sourceDest[result[i].Destination].Add(result[i].Source);
				string s = result[i].Source + "-" + result[i].Destination;
				if (!edgesDic.ContainsKey(s))
					edgesDic.Add(s, result[i].Weight);
				s= result[i].Destination + "-" + result[i].Source;
				if (!edgesDic.ContainsKey(s))
					edgesDic.Add(s, result[i].Weight);

			}
			int bestStart = -1;
			for(i=0;i<clusters.Count;i++)
            {
				if (startNode.setStruct[0] == clusters[i][0])
					bestStart = i;
            }
			List<List<string>> res = new List<List<string>>();
			Queue<int> q = new Queue<int>();
			List<string> dicKeys = new List<string>();
			if (sourceDest.ContainsKey(bestStart))
			{
				res.Add(clusters[bestStart]);
				used.Add(bestStart);
				foreach (var item in sourceDest[bestStart])
				{
					dicKeys.Add(bestStart + "-" + item);					
				}
				dicKeys.Sort(delegate  (string a, string b) { if (edgesDic[a] > edgesDic[b]) return 0; return 1; });
				foreach (var item in dicKeys)
				{
					string[] aux = item.Split('-');
					q.Enqueue(Convert.ToInt32(aux[1]));
				}
				sourceDest.Remove(bestStart);
			}
			else
				throw new Exception("Ups no best start");
			while(q.Count>0)
            {
				int w = q.Dequeue();
				if (sourceDest.ContainsKey(w))
				{
					dicKeys.Clear();					
					foreach(var item in sourceDest[w])
						if(!used.Contains(item))
							dicKeys.Add(w + "-" + item);

					if (dicKeys.Count > 0)
					{
						dicKeys.Sort(delegate (string a, string b) { if (edgesDic[a] > edgesDic[b]) return 0; return 1; });
						foreach (var item in dicKeys)
						{
							string[] aux = item.Split('-');
							q.Enqueue(Convert.ToInt32(aux[1]));
						}
					}
					sourceDest.Remove(w);
				}
				if (!used.Contains(w))
				{
					res.Add(clusters[w]);
					//if(used.Contains(w))
					used.Add(w);
				}
				
			}
			

			return res;
		}

	}
}
