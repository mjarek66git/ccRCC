using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace phiClustCore.Profiles
{
    public enum SIMDIST
    {
        DISTANCE,
        SIMILARITY
    };


    public class ProfileAutomatic
    {
        public static string similarityProfileName="automatic_similarity.profile";
        public static string distanceProfileName = "automatic_distance.profile";
        public static SerializableDictionary<string, SerializableDictionary<string, double>> GenerateWeights(List<string> states, SIMDIST weightFlag)
        {
            SerializableDictionary<string, SerializableDictionary<string, double>> weights = new SerializableDictionary<string, SerializableDictionary<string, double>>();
            if (weightFlag == SIMDIST.SIMILARITY)
            {
                foreach (var item in states)
                {
                    weights.Add(item, new SerializableDictionary<string, double>());
                    weights[item].Add(item, 1.0);
                }
            }
            else
                foreach (var item in states)
                {
                    weights.Add(item, new SerializableDictionary<string, double>());
                    foreach (var item1 in states)
                    {
                        if (item != item1)                                                    
                            weights[item].Add(item1, 1.0);                        
                    }
                }


            return weights;
        }
        public static ProfileTree GenerateProfile(Dictionary<string,List<byte>> profile, SIMDIST similarityFlag, string internalName = null)
        {
            ProfileTree t = new ProfileTree();
            Dictionary<string, Dictionary<string, int>> dic = new Dictionary<string, Dictionary<string, int>>();
            string profName = "Dummy";
            dic.Add(profName, new Dictionary<string, int>());
            foreach(var item in profile)
            {
                for (int i = 0; i < item.Value.Count; i++)
                {
                    if (!dic[profName].ContainsKey(item.Value[i].ToString()))
                        dic[profName].Add(item.Value[i].ToString(), 0);
                }
            }

            foreach (var item in dic)
            {
                profileNode node = new profileNode();
                node.active = true;
                if (internalName != null)
                    node.internalName = internalName;
                else
                    node.internalName = "User defined profile";
                node.profName = item.Key;
                foreach (var itemK in item.Value)
                    node.AddStateItem(itemK.Key, itemK.Key);

                node.profWeights = GenerateWeights(new List<string>(item.Value.Keys), similarityFlag);
                t.AdddNode("/", node);
            }            
            return t;
        }
        public static ProfileTree AnalyseProfileFile(string fileName, SIMDIST similarityFlag,string internalName=null)
        {
            ProfileTree t = new ProfileTree();

            StreamReader wr;
            if (fileName == null || !File.Exists(fileName))
                throw new Exception("File:" + fileName + " not exists");

            wr = new StreamReader(fileName);
            string line = wr.ReadLine();
            
            Dictionary<string, Dictionary<string, int>> dic = new Dictionary<string, Dictionary<string, int>>();
            while (line != null)
            {
                if (line.Contains(">"))
                {
                    line = wr.ReadLine();
                    while (line != null && line[0] != '>')
                    {
                        if (line.Contains("profile") && !line.Contains("SEQ"))
                        {
                            string[] tmp = line.Split(new string[] { " profile " }, StringSplitOptions.None);
                            if (!dic.ContainsKey(tmp[0]))
                                dic.Add(tmp[0], new Dictionary<string, int>());
                            string[] aux;
                            if (tmp[1].Contains(" "))
                                aux = tmp[1].Split(' ');
                            else
                            {
                                aux = new string[tmp[1].Length];
                                for (int i = 0; i < tmp[1].Length; i++)
                                    aux[i] = tmp[1][i].ToString();
                            }
                            foreach (var item in aux)
                                if (item != "-" && item!="")
                                    if (!dic[tmp[0]].ContainsKey(item))
                                           dic[tmp[0]].Add(item, 0);

                        }
                        line = wr.ReadLine();
                    }
                }
                else
                    line = wr.ReadLine();
            }
            wr.Close();

            if (dic.Keys.Count == 0)
                throw new Exception("File " + fileName + " is not valid Profile file!");

            foreach (var item in dic)
            {
                profileNode node = new profileNode();
                node.active = true;
                if (internalName != null)
                    node.internalName = internalName;
                else
                    node.internalName = "User defined profile";
                node.profName = item.Key;
                foreach (var itemK in item.Value)
                    node.AddStateItem(itemK.Key, itemK.Key);

                node.profWeights = GenerateWeights(new List<string>(item.Value.Keys), similarityFlag);
                t.AdddNode("/", node);
            }

          //  t.pr


            return t;
        }
    }
}
