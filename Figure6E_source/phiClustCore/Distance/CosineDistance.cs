using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore.Distance
{
    class CosineDistance : JuryDistance
    {
        public CosineDistance(Dictionary<string,List<double>> data, bool flag):
                base(data,flag)
        {

        }
        public override string ToString()
        {
            return "Cosine";
        }
        public override List<KeyValuePair<string, double>> GetReferenceList(List<string> structures)
        {
            if (jury != null)
                //return jury.JuryOpt(structures).juryLike[0].Key;
                return jury.JuryOptWeights(structures).juryLike;
            //return jury.ConsensusJury(structures).juryLike;

            List<KeyValuePair<string, double>> refList = new List<KeyValuePair<string, double>>();
            double[] refPos = new double[data[structures[0]].Count];
            for(int i=0;i<structures.Count;i++)
            {
                List<double>mod1 = data[structures[i]];
                for (int j = 0; j < mod1.Count; j++)
                    refPos[j] += mod1[j];
            }
            for (int j = 0; j < refPos.Length; j++)
                refPos[j] /= structures.Count;

            double dl2 = 0;
            for (int j = 0; j < refPos.Length; j++)
                dl2 += refPos[j] * refPos[j];

            for (int i = 0; i < structures.Count; i++)
            {
                double dist = 0;
                List<double> mod1 = data[structures[i]];
                //for (int j = 0; j < mod1.Count; j++)
                  //  dist += (mod1[j] - refPos[j]) * (mod1[j] - refPos[j]);
                double il = 0, dl1 = 0;
                for (int j = 0; j < mod1.Count; j++)
                {
                    // dist += (mod1[j] - mod2[j]) * (mod1[j] - mod2[j]);
                    il += mod1[j] * refPos[j];
                    dl1 += mod1[j] * mod1[j];
                    
                }
                dist = (int)((1.0 - (double)(il / (Math.Sqrt(dl1) * Math.Sqrt(dl2)))) * 100);

                KeyValuePair<string, double> aux = new KeyValuePair<string, double>(structures[i], dist);
                refList.Add(aux);
            }
            refList.Sort((nextPair, firstPair) =>
            {
                return nextPair.Value.CompareTo(firstPair.Value);
            });
            return refList;

        }
        public override string GetReferenceStructure(List<string> structures)
        {
            List<KeyValuePair<string, double>> refList = null;
            refList = GetReferenceList(structures);

            return refList[0].Key;
        }
        public override  int GetDistance(string refStructure, string modelStructure)
        {
            int dist = 0;
            if (!data.ContainsKey(refStructure))
                throw new Exception("Structure: " + refStructure + " does not exists in the available list of structures");

            if (!data.ContainsKey(modelStructure))
                throw new Exception("Structure: " + modelStructure + " does not exists in the available list of structures");

            List<double> mod1 = data[refStructure];
            List<double> mod2 = data[modelStructure];
            double il=0, dl1=0, dl2=0;
            for (int j = 0; j < data[refStructure].Count; j++)
            {
               // dist += (mod1[j] - mod2[j]) * (mod1[j] - mod2[j]);
                il += mod1[j] * mod2[j];
                dl1 += mod1[j] * mod1[j];
                dl2 += mod2[j] * mod2[j];
            }
            dist = (int)((1.0 - (double)(il /( Math.Sqrt(dl1) * Math.Sqrt(dl2))))*100);
            return dist;
        }
    }
}
