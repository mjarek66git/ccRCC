using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using phiClustCore;

namespace phiClustCore.Distance

{
    public class JuryDistance: HammingBase//,IDistance
    {

        public JuryDistance(Dictionary<string,List<double>> data, bool flag) 
            :base(data,flag)
        {

        }
        public void WriteStates()
        {
            StreamWriter wr = new System.IO.StreamWriter("symbols.txt");

            foreach(var item in lStates)
            {
                foreach (var iStates in states.Keys)
                    if (item.ContainsKey(iStates))
                        wr.Write(item[iStates] + " ");
                    else
                        wr.Write("0 ");
                wr.WriteLine();
            }

        }
      

        public override string ToString()
        {
            return "HAMMING";
        }
        public override bool SimilarityThreshold(float threshold, float dist)
        {
            if (dist < threshold)
                return true;
            return false;
        }


        public override int GetDistance(string refStructure, string modelStructure)
        {
            double lDist = 0;
            Dictionary <string,int> dicStates=new Dictionary<string,int>();

            if (!data.ContainsKey(refStructure) || !data.ContainsKey(modelStructure))
                return errorValue;

            List<double> refL = data[refStructure];
            List<double> modL = data[modelStructure];

         
            for (int j = 0; j < refL.Count; j++)
            {
                if (refL[j] !=modL[j])//= 0 || modL[j] == 0)
                {
                    lDist += 1;
                    continue;
                }

            }
            return (int)(lDist);

        }
         
    }
}
