using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore.Distance
{
    class Euclidian: JuryDistance

    {
    
        public Euclidian(Dictionary<string,List<double>> data, bool flag):
                base(data,flag)
        {

        }
        public override string ToString()
        {
            return "Euclidean";
        }
        public override List<KeyValuePair<string, double>> GetReferenceList(List<string> structures)
        {
            if (jury != null)
                //return jury.JuryOpt(structures).juryLike[0].Key;
                return jury.JuryOptWeights(structures).juryLike;
            //return jury.ConsensusJury(structures).juryLike;

            List<KeyValuePair<string, double>> refList = new List<KeyValuePair<string, double>>();
            double[] refPos = new double[data[structures[0]].Count];
            for (int i = 0; i < structures.Count; i++)
            {
                List<double> mod1 = data[structures[i]];
                for (int j = 0; j < mod1.Count; j++)
                    refPos[j] += mod1[j];
            }
            for (int j = 0; j < refPos.Length; j++)
                refPos[j] /= structures.Count;

            for (int i = 0; i < structures.Count; i++)
            {
                double dist = 0;
                List<double> mod1 = data[structures[i]];
                //for (int j = 0; j < mod1.Count; j++)
                //  dist += (mod1[j] - refPos[j]) * (mod1[j] - refPos[j]);
                for (int j = 0; j < mod1.Count; j++)
                {
                    // dist += (mod1[j] - mod2[j]) * (mod1[j] - mod2[j]);
                    dist += (mod1[j] - refPos[j]) * (mod1[j] - refPos[j]);

                }

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
        public override int GetDistance(string refStructure, string modelStructure)
        {
            double dist = 0;
            if (!data.ContainsKey(refStructure))
                throw new Exception("Structure: " + refStructure + " does not exists in the available list of structures");

            if (!data.ContainsKey(modelStructure))
                throw new Exception("Structure: " + modelStructure + " does not exists in the available list of structures");

            List<double> mod1 = data[refStructure];
            List<double> mod2 = data[modelStructure];
            for (int j = 0; j < data[refStructure].Count; j++)
            {
                 dist += (mod1[j] - mod2[j]) * (mod1[j] - mod2[j]);
            }
            return (int)(dist*100);
        }

    }
}
