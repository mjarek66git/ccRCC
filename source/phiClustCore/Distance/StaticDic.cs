using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore.Distance
{
    public static class StaticDic
    {
        public static Dictionary<string, string> Id = new Dictionary<string, string>();
        

        public static List<List<string>> GetClusters(Dictionary<string,List<double>> refData)
        {
            List<List<string>> res = new List<List<string>>();
            if (Id.Count == 0)
                return null;
            var stage = LoadStage(@"C:\projects\bioinfo\Cluster_patient\newData\stage_id");
            Dictionary<int, List<string>> dic = new Dictionary<int, List<string>>();
            foreach(var item in Id)
            {
                int k = Convert.ToInt32(item.Value);
                if (!dic.ContainsKey(k))
                    dic.Add(k, new List<string>());

                if(refData.ContainsKey(item.Key))
                    dic[k].Add(item.Key);
            }
            for(int i=0;i<dic.Count;i++)
            {
                res.Add(dic[i]);
            }
            string[] order = { "Stage1", "Stage3","Metastatic" };
            for(int i=0;i<res.Count;i++)
            {
                Dictionary<string, List<string>> toOrder = new Dictionary<string, List<string>>();
                for(int j=0;j<res[i].Count;j++)
                {
                    if (stage.ContainsKey(res[i][j]))
                    { 
                        if (!toOrder.ContainsKey(stage[res[i][j]]))
                            toOrder.Add(stage[res[i][j]], new List<string>());
                        toOrder[stage[res[i][j]]].Add(res[i][j]);
                    }
                }
                res[i] = new List<string>();
                foreach(var item in order)
                {
                    if(toOrder.ContainsKey(item))
                        res[i].AddRange(toOrder[item]);
                }

            }

            return res;
        }
        
        public static void LoadId(string fileName)
        {
            StreamReader st = new StreamReader(fileName);
            string line = st.ReadLine();
            bool dataFlag = false;
            while(line!=null)
            {
                if (dataFlag)
                {
                    string[] aux = line.Split(' ');
                    Id.Add(aux[0], aux[1]);
                }
                if (line.Contains("[DATA]"))
                    dataFlag = true;

                line = st.ReadLine();
            }
            st.Close();
        }
        public static Dictionary<string, string> LoadStage(string fileName)
        {
            Dictionary<string, string> stage = new Dictionary<string, string>();
            StreamReader st = new StreamReader(fileName);
            string line = st.ReadLine();
            bool dataFlag = false;
            while (line != null)
            {
                if (dataFlag)
                {
                    string[] aux = line.Split(' ');
                    stage.Add(aux[0], aux[1]);
                }
                if (line.Contains("[DATA]"))
                    dataFlag = true;

                line = st.ReadLine();
            }
            st.Close();
            return stage;
        }


    }

}
