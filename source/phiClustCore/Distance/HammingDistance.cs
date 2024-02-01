using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using phiClustCore;

namespace phiClustCore.Distance
{
    public class HammingDistance:HammingBase
    {
        public HammingDistance(Dictionary<string,List<double>>data, bool flag):
            base(data,flag)
        {
        }
        
        public override int GetDistance(string refStructure, string modelStructure)
        {
            int dist = 0;
            if(!data.ContainsKey(refStructure))
                    throw new Exception("Structure: "+refStructure+" does not exists in the available list of structures");

            if(!data.ContainsKey(modelStructure))
                    throw new Exception("Structure: "+modelStructure+" does not exists in the available list of structures");

            List<double> mod1 = data[refStructure];
            List<double> mod2 = data[modelStructure];
            for (int j = 0; j < data[refStructure].Count; j++)
            {
                if(mod1[j]!=mod2[j])// || mod1[j]==0)
                    dist++;
            }

            return dist;
        }
        public override string ToString()
        {
            return "HAMMING";
        }

    }

}
