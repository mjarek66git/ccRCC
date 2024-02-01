



























using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace phiClustCore
{
    class GuidedHashCluster:HashCluster
    {
        Dictionary<string, string> classification = new Dictionary<string, string>();

        public GuidedHashCluster(Dictionary<string,List<double>> data,Options input):base(data,input)
        {
        }
        public void ReadClassLabels(string fileName)
        {
            StreamReader file = new StreamReader(fileName);
            string line = file.ReadLine();

            while(line!=null)
            {
                string []aux=line.Split(' ');
                if(aux.Length==2)
                {
                    classification.Add(aux[0], aux[1]);
                }
                line = file.ReadLine();
            }

            file.Close();
        }
        public override double[] CalcEntropy(Dictionary<double, int>[] locColumns)
        {
            double[] entropy = base.CalcEntropy(locColumns);
            Dictionary<string, List<string>> classLabels = new Dictionary<string,List<string>>();
            
            foreach(var item in classification)
            {
                if (!classLabels.ContainsKey(item.Value))
                {

                    List<string> xx = new List<string>();
                    xx.Add(item.Key);
                    classLabels.Add(item.Value, xx);
                }
                else
                    classLabels[item.Value].Add(item.Key);

            }
            Dictionary<string,double []> classDepColumns = new Dictionary<string,double[]>();
            foreach(var item in classLabels.Keys)
            {
                int counter=0;
                for (int i = 0; i < structNames.Count; i++)
                    if (classification.ContainsKey(structNames[i]))
                        counter++;
                Dictionary<double, int>[] aux = MakeColumnsLists(classLabels[item], this.data);
                counter = 0;

                double[] classEntropy = base.CalcEntropy(aux);

                classDepColumns.Add(item, classEntropy);

            }
            /*for (int i = 0; i < entropy.Length; i++)
            {
                double sum = 0;
                foreach (var item in classDepColumns.Keys)
                    sum += classDepColumns[item][i] / classLabels.Count;
                entropy[i] += sum / locColumns[i].Count;

            }*/

    
            return entropy;
        }
    }
}
