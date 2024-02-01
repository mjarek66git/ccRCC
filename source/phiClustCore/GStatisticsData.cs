using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore
{
    class GStatisticsData
    {
        int numFeatures;
        int numData;
        List<string> keys;

        public GStatisticsData(List <string> keys)
        {
            this.keys = keys;
            this.numFeatures = keys[0].Length;
            this.numData = keys.Count;
        }
        private double [] PrepareStatistics()
        {
            double[] prob = new double[numFeatures];

            foreach(var item in keys)
            {
                for (int i = 0; i < item.Length; i++)
                    if (item[i] == '1')
                        prob[i]++;                    
            }
            for (int i = 0; i < prob.Length; i++)
                prob[i] /= keys.Count;

                return prob;
        }
        public Dictionary<string,string> GenerateData()
        {
            Dictionary<string, string> dic = new Dictionary<string, string>();
            Random r = new Random();

            double[] probabilities = PrepareStatistics();

            for (int i = 0; i < numData; i++)
            {
                StringBuilder key = new StringBuilder() ;
                for (int j = 0; j < numFeatures; j++)
                {
                    double x = r.NextDouble();
                    if (x<=probabilities[j])
                        key.Append('1');
                    else
                        key.Append('0');
                }
                string name = "test" + i;
                dic.Add(name,key.ToString());
            }
            return dic;
        }
    }
}
