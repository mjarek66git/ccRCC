using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore.Distance
{
    class Pearson : JuryDistance
    {
        public Pearson(Dictionary<string,List<double>> data, bool flag):
                base(data,flag)
        {
            //            if(StaticDic.Id.Count==0)
            StaticDic.Id.Clear();
            StaticDic.LoadId(@"C:\projects\bioinfo\Cluster_patient\newData\clusters_id53");
            //StaticDic.LoadId(@"C:\projects\bioinfo\Cluster_patient\newData\stage_id");
        }
        public override List<KeyValuePair<string, double>> GetReferenceList(List<string> structures)
        {
          //  if (jury != null)
                //return jury.JuryOpt(structures).juryLike[0].Key;
            //    return jury.JuryOptWeights(structures).juryLike;
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

            double refMod = 0;
            for (int j = 0; j < refPos.Length; j++)            
                refMod += refPos[j];
            refMod /= refPos.Length;

                for (int i = 0; i < structures.Count; i++)
            {
                double dist = 0;
                List<double> mod1 = data[structures[i]];
                double Sxx = 0;
                double Sxy = 0;
                double Syy = 0;

                double avrMod = 0;
              

                for (int j = 0; j < refPos.Length; j++)                
                    avrMod += mod1[j];

                avrMod /= mod1.Count;

                for (int j = 0; j < mod1.Count; j++)
                {
                    Sxx+=mod1[j]*mod1[j];
                    Syy+=refPos[j]*refPos[j];
                    Sxy+=mod1[j]*refPos[j];
                }
             //   Sxx -= mod1.Count * avrMod * avrMod;
                //Syy-= mod1.Count * avr * avr;
                //Sxy-=mod1.Count*avr*avrMod;
                dist = (Sxy - refPos.Length* refMod *avrMod) / (Math.Sqrt((Sxx - mod1.Count * avrMod * avrMod) * (Syy - refPos.Length * refMod * refMod)));
                dist = 1 - dist;

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
            double avrMod1=0,avrMod2=0;
            for(int j=0;j<mod1.Count;j++)
            {
                avrMod1 += mod1[j];
                avrMod2 += mod2[j];
            }

            avrMod1 /= mod1.Count;
            avrMod2 /= mod2.Count;


            double Sxx = 0;
            double Sxy = 0;
            double Syy = 0;
            for (int j = 0; j < mod1.Count; j++)
            {
                Sxx += mod1[j] * mod1[j];
                Syy += mod2[j] * mod2[j];
                Sxy += mod1[j] * mod2[j];
            }
            /*Sxx -= mod1.Count * avrMod1 * avrMod1;
            Syy -= mod2.Count * avrMod2 * avrMod2;
            Sxy -= mod1.Count * avrMod2 * avrMod1;*/
            double res1 = Sxx - mod1.Count * avrMod1 * avrMod1;
            double res2 = Syy - mod2.Count * avrMod2 * avrMod2;
            double vv = 0;
            if (res1 > 0 && res2 > 0)
                vv = (Sxy - mod1.Count * avrMod1 * avrMod2) / (Math.Sqrt((Sxx - mod1.Count * avrMod1 * avrMod1) * (Syy - mod2.Count * avrMod2 * avrMod2)));
            else
                vv = 0;
            dist = (1.0-vv)*100;
            /*if (StaticDic.Id.ContainsKey(refStructure) && StaticDic.Id.ContainsKey(modelStructure))
                if (StaticDic.Id[refStructure].Equals(StaticDic.Id[modelStructure]))
                    dist -= 10;
                else
                    dist += 10;*/
            if (dist > 100)
                dist = 100;
            if (dist < 0)
                dist = 0;

            return (int)dist;
        }
        public override string ToString()
        {
            return "Pearson";
        }
    }
}
