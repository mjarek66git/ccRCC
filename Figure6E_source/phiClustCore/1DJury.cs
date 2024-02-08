using System;
using System.Collections;
using System.IO;
using System.Collections.Specialized;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using phiClustCore;
using System.Data;
using System.Threading;
using phiClustCore.Interface;
namespace phiClustCore
{
    [Serializable]
	public class jury1D:IProgressBar
	{
        Dictionary<string, List<double>> data;
                List<KeyValuePair<string, double>> res1Djury;
        public List<string> alignKeys
        {
            get
            {
                if (data != null && data.Keys != null)
                    return new List<string>(data.Keys);
                else
                    return null;
            }
        }
       
        public Dictionary<double, int>[] columns = null;
        List<string> allStructures;
        [NonSerialized]
        ManualResetEvent[] resetEvents = null;
        long maxV, currentV;
        Settings set = new Settings();
        int threadNumbers;

        double startProgress = 0;
        double endProgress = 1;

        public double StartProgress { set { startProgress = value; } get { return startProgress; } }
        public double EndProgress { set { endProgress = value; } get { return endProgress; } }

        public jury1D()
        {
            maxV = 1;
            currentV = 0;
            set.Load();
            threadNumbers = set.numberOfCores;
        }
        public jury1D(Dictionary<double, int>[] freq1,Dictionary<string,List<double>> data)
        {
            columns = freq1;          
            this.data = data;
            set.Load();
            threadNumbers = set.numberOfCores;

        }
        public double ProgressUpdate()
        {
            return StartProgress+(EndProgress - StartProgress) * ((double)currentV / maxV);
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }


        public void PrepareJury(Dictionary<string,List<double>> data)
        {
            this.data = data;

        }
        public void PrepareJury(List<KeyValuePair<string,string>> profiles,string profName, string profileFile)
        {
        }
        //all profiles in profileFile should be aligned
        public void PrepareJury()
        {
        }
        public override string ToString()
        {
            return "1DJury";
        }
        public Dictionary<double, List<int>>[] MakeColumnsLists(List<string> structNames)
        {
            double locState=0;
            Dictionary<double, List<int>>[] columns = null;

            if (structNames.Count == 0)
                return null;

            columns = new Dictionary<double, List<int>>[data[structNames[0]].Count];
            for (int i = 0; i < columns.Length; i++)
            {

                
                Dictionary<double, List<int>> dicCol = new Dictionary<double, List<int>>();
                for (int j = 0; j < structNames.Count; j++)
                {
                    if (data.ContainsKey(structNames[j]) && data[structNames[j]].Count > 0)
                        locState = data[structNames[j]][i];
                    else
                        continue;
                    if (locState == 0)
                        continue;

                    if (!dicCol.ContainsKey(locState))
                    {
                        List<int> lista = new List<int>();
                        lista.Add(j);
                        dicCol.Add(locState, lista);
                    }
                    else
                        dicCol[locState].Add(j);
                }

                columns[i] = dicCol;
            }
            
            return columns;
        }
        private void ThreadingMakeColumns(object o)
        {
            double locState;
            ThreadParam pp = (ThreadParam)o;
            int threadNum = pp.num;
            int start = pp.start;
            int stop = pp.stop;

            for (int i = start; i < stop; i++)
            {
                lock (columns[i])
                {
                   // int counter = 0;
                    foreach (var name in allStructures)
                    {
                        locState = data[name][i];
                        if (!columns[i].ContainsKey(locState))
                        {
                            columns[i].Add(locState, 1);
                            continue;
                        }
                        
                        //columns[i].AddOrUpdate(locState,0, (key, value) => value + 1);
                        columns[i][locState]++;//, 0, (key, value) => value + 1);
                    }

                }                
                Interlocked.Increment(ref currentV);
            }
            resetEvents[threadNum].Set();
        }
        public Dictionary<double, int>[] MakeColumns(List<string> structNames)
        {           

            if (structNames.Count == 0)
                return null;

            maxV = data[structNames[0]].Count + data.Count;

            columns = new Dictionary<double, int>[data[structNames[0]].Count];
            allStructures = structNames;

            resetEvents = new ManualResetEvent[threadNumbers];

            for (int i = 0; i < columns.Length; i++)
            {
                columns[i] = new Dictionary<double, int>();
           //     foreach (var item in weights.Keys)
             //       columns[i].Add(item,0);
            }
            for (int n = 0; n < threadNumbers; n++)
            {
                ThreadParam pp = new ThreadParam();
                pp.num = n;
                pp.start = (int)(n * columns.Length / Convert.ToDouble(threadNumbers));
                pp.stop=(int)( (n + 1) * columns.Length / Convert.ToDouble(threadNumbers));
                //int p = n;
                //int start = (int)(n * columns.Length /Convert.ToDouble(threadNumbers));
                //int stop =(int)( (n + 1) * columns.Length / Convert.ToDouble(threadNumbers));
                resetEvents[n] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(new WaitCallback(ThreadingMakeColumns), (object) pp);
                

            }
            for (int n = 0; n < threadNumbers; n++)
                resetEvents[n].WaitOne();
            return columns;
        }
        public List<double> GetStructureStates(string structName)
        {

            if(data.ContainsKey(structName))
                return data[structName];

            return null;
        }
        private void ThreadingScoreCalc(object o)
        {
            double score;
            double locState;

            //object[] array = o as object[];
            //int threadNum = (int)array[0];
            //int start = (int)array[1];
            //int stop = (int)array[2];
            ThreadParam k = (ThreadParam)o;
            int threadNum = k.num;
            int start = k.start;
            int stop = k.stop;
            for (int j = start; j < stop;j++ )
            {
                score = 0;
                string name = allStructures[j];
                List<double> listStates = data[name];
                for (int i = 0; i < listStates.Count; i++)
                {
                    score += columns[i][listStates[i]];
                }
                score /= allStructures.Count * listStates.Count;

                KeyValuePair<string, double> v = new KeyValuePair<string, double>(name, score);
                lock (res1Djury)
                {
                    res1Djury.Add(v);
                }
                Interlocked.Increment(ref currentV);
            }
            resetEvents[threadNum].Set();
            
        }
        public ClusterOutput JuryOptWeights(List<string> structNames,Dictionary<double,int>[]locColumns=null)
        {
            res1Djury = new List<KeyValuePair<string, double>>(structNames.Count);
            Dictionary<double, int>[] columns = null;

            List<string> aux = new List<string>(structNames);
            foreach (var item in structNames)
            {
                if (!data.ContainsKey(item) || data[item].Count == 0)
                    aux.Remove(item);
            }
            if (locColumns == null)
                columns = MakeColumns(aux);
            else
            {
                columns = locColumns;
                this.columns = locColumns;
            }

            
            currentV++;
            if (columns == null)
                return null;
            
            allStructures = new List<string>(aux);
                        
            resetEvents = new ManualResetEvent[threadNumbers];

            for (int n = 0; n < threadNumbers; n++)
            {
                ThreadParam pp = new ThreadParam();
                pp.num = n;
                pp.start = (int)(n * allStructures.Count / Convert.ToDouble(threadNumbers));
                pp.stop = (int)((n + 1) * allStructures.Count / Convert.ToDouble(threadNumbers));
                resetEvents[n] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(new WaitCallback(ThreadingScoreCalc), (object)pp);

            }
            for (int n = 0; n < threadNumbers; n++)
                resetEvents[n].WaitOne();
            
            currentV++;
            res1Djury.Sort((firstPair, nextPair) =>
            {
                return nextPair.Value.CompareTo(firstPair.Value);
            });

            ClusterOutput juryRes = new ClusterOutput();
            juryRes.runParameters = "Profile: ";
            juryRes.juryLike = res1Djury;
            currentV = maxV;
            
            return juryRes;
			

        }
        public ClusterOutput ConsensusJury(List<string> structNames)
        {
            List<KeyValuePair<string, double>> distCons = new List<KeyValuePair<string, double>>();
            double locState;
            List<Dictionary<double, int>> cons = new List<Dictionary<double,int>>(structNames.Count);
            List<double> finalCons = new List<double>();

            foreach (string name in structNames)
            {
                for (int i = 0; i < data[name].Count; i++)
                     cons.Add(new Dictionary<double, int>());
                break;
            }

            foreach (string name in structNames)
            {
                if (!data.ContainsKey(name))
                    continue;
                for (int i = 0; i < data[name].Count; i++)
                {
                    locState = data[name][i];
                    if (cons[i].ContainsKey(locState))
                        cons[i][locState]++;
                    else
                        cons[i].Add(locState, 1);
                }

            }
            foreach (var item in cons)
            {
                var items = from pair in item
                            orderby pair.Value descending
                            select pair;
                // Display results.
                foreach (KeyValuePair<double, int> pair in items)
                {
                    finalCons.Add(pair.Key);
                    break;
                }

            }
            foreach (string name in structNames)
            {
                if (!data.ContainsKey(name))
                    continue;
                float dist = 0;
                for (int i = 0; i < data[name].Count; i++)
                {
                    if (data[name][i] == finalCons[i])
                        dist++;
                }
                KeyValuePair<string, double> v = new KeyValuePair<string, double>(name, dist);
                distCons.Add(v);
            }
            distCons.Sort((firstPair, nextPair) =>
            {
                return firstPair.Value.CompareTo(nextPair.Value);
            });

            ClusterOutput juryRes = new ClusterOutput();

            juryRes.juryLike = distCons;
            return juryRes;

        }
/*        [Obsolete("Instead of this method use JuryOptWeights, this method has fixed weights values")]
		public ClusterOutput JuryOpt(List <string> structNames)
		{
			string locState;
			List<KeyValuePair<string,double>> res1Djury=new List<KeyValuePair<string,double>>();
			Dictionary <string, float> score= new Dictionary<string, float>();
            Dictionary<string, int>[] columns = null;


            columns = MakeColumns(structNames);

            if (columns == null)
                return null;

			foreach (string name in structNames)
				score[name]=0;
			
			
			foreach (string name in structNames)
			{
                for (int i = 0; i < stateAlign[name].Count; i++)
				{
                    locState = stateAlign[name][i];
					if(locState=="-")
						continue;
					if(locState.Contains("H") || locState.Contains("E"))
						score[name]+=weightHE*columns[i][locState];		
                    else
                        score[name] += weightC * columns[i][locState];		
					
					if(locState.Contains("H") || locState.Contains("E"))
					{
						//weight=0.5f;
						locState=locState.Replace("H","C");
						locState=locState.Replace("E","C");
						if(columns[i].ContainsKey(locState))
							score[name]+=weightC*columns[i][locState];														
					}
                    if (stateAlign[name][i].Contains("C"))
				    {
                        locState = stateAlign[name][i];					
						locState=locState.Replace("C","H");
						if(columns[i].ContainsKey(locState))
							if(score.ContainsKey(name))
								score[name]+=weightC*columns[i][locState];		
							else 
								score[name]=weightC*columns[i][locState];

                        locState = stateAlign[name][i];					
						locState=locState.Replace("C","E");
						if(columns[i].ContainsKey(locState))
							if(score.ContainsKey(name))
								score[name]+=weightC*columns[i][locState];		
							else 
								score[name]=weightC*columns[i][locState];
							
					}
				}
                score[name] /= structNames.Count * stateAlign[name].Count;
				
				KeyValuePair<string,double> v=new KeyValuePair<string, double>(name,score[name]);				
				res1Djury.Add(v);
			}			
			res1Djury.Sort((firstPair,nextPair) =>
    		{
        		return nextPair.Value.CompareTo(firstPair.Value);
    		});

			//Console.WriteLine("Next");
			//foreach (KeyValuePair<string,double> s in res1Djury)
			//	Console.WriteLine("Result "+s.Value+" VALUE="+s.Key);
            ClusterOutput juryRes = new ClusterOutput();

            juryRes.juryLike = res1Djury;
			return juryRes;
			
		}	*/
							
	}
	
}

