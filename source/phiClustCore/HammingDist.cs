using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using phiClustCore.Interface;
using System.Data;
using System.IO;

namespace phiClustCore
{
    class HammingDist:RootHash,IProgressBar
    {
        double startProgress = 0;
        double endProgress = 1;
        protected int maxV = 100;
        protected int currentV = 0;
        Dictionary<string, List<double>> baseDouble = new Dictionary<string, List<double>>();
        public HammingDist(HashCluster hk, ClusterOutput outp, HNNCInput opt):base(1,hk, outp, opt)
        {
        }

        public double StartProgress { set { startProgress = value; } get { return startProgress; } }
        public double EndProgress { set { endProgress = value; } get { return endProgress; } }
        public double ProgressUpdate()
        {
            double res = 0;
            return StartProgress + (EndProgress - StartProgress) * (res * 0.45 + 0.55 * (double)currentV / maxV);
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }
        public override void Preprocessing()
        {
            base.Preprocessing(); 
            foreach(var item in caseBase.Keys)
            {
                List<double> s = new List<double>();
                for (int i = 0; i < caseBase[item].Length; i++)
                {
                    int w = caseBase[item][i];
                    s.Add(Convert.ToDouble(w));
                }

                baseDouble.Add(item, s);
            }

        }
        public override void CalcHashDist(object o)
        {
            ThreadParam pp = (ThreadParam)o;
            int threadNum = pp.num;
            int start = pp.start;
            int stop = pp.stop;
            int[] index = new int[dist[threadNum].Length];
                Dictionary<string, string> locDic = new Dictionary<string, string>(caseBase);
            int[] locDist = dist[threadNum];
            for (int n = start; n < stop; n++)
            {
                string testItem = testList[n];
                try
                {

                    List<double> key;
                    if (!baseDouble.ContainsKey(testItem))
                        continue;
                    key = baseDouble[testItem];
                    string w = CalcDist(key, threadNum, index, locDic,locDist);

                    lock (results)
                    {

                        results.Add(testItem, w);
                        //wr.WriteLine("size=" + results.Count);
                    }
                }
                catch (Exception ex)
                {

                    StreamWriter qq = new StreamWriter("exception.dat");
                    qq.WriteLine("exception=" + ex.Message+" key=",testItem);
                    qq.Close();
                    //wr.Close();
                }
            }
            
            
            resetEvents[threadNum].Set();
            
        }
        public string CalcDist(List<double> key,int threadNum,int []index,Dictionary<string,string> locDic,int []locDist)
        {                    
            for (int i = 0; i < locDist.Length; i++)
                locDist[i] = 0;
            for (int j = 0; j < locDic.Count; j++)
            {
                List<double> baseKeys = baseDouble[dataBaseKeys[j]];
                double res = 0;
                for (int i = 0; i < key.Count; i++)
                {
                    /*    if (baseKeys[i] != key[i])
                            count++;                        */
                    double r = baseKeys[i]- key[i];
                    res += r;
                }
                //res = key.Count - res;
                locDist[j] = Convert.ToInt32(res);
                //locDist[j]=count;
                index[j] = j;
            }
//            Array.Sort(locDist, index);
//            Array.Reverse(dist[threadNum]);
//            Array.Reverse(index);
            string w = "";
            for (int i = 0; i < resSize; i++)
                w += dataBaseKeys[index[i]] + ":";
            w += dataBaseKeys[index[resSize]];

            return w;

        }

    }
}
