using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using phiClustCore.Distance;

namespace phiClustCore
{
    class HammingConsensus : HammingBase
    {
        List<double> consensusStates = new List<double>();
        Dictionary<string, string> consensus = new Dictionary<string, string>();
        public Dictionary<string,int> distanceOrdered = new Dictionary<string,int>();
       
        public HammingConsensus(Dictionary<string,List<double>>data, bool flag)
            :base(data,flag)
        {
        }

        public override string ToString()
        {
            return "HAMMING";
        }
        private string TransformToConsensusStates(string structName)
        {
            double locState;
            string consensusRepr = "";

            if (!data.ContainsKey(structName))
                return null;

            for (int i = 0; i < data[structName].Count; i++)
            {
                locState = data[structName][i];
                if (locState == 0)
                    continue;

                if (consensusStates[i] == locState)
                    consensusRepr += "0";
                else
                    consensusRepr += "1";

            }
            return consensusRepr;
        }
        private int DistanceToConsensus(string structName)
        {
            int dist = 0;
            for (int i = 0; i < consensus[structName].Length; i++)
                if (consensus[structName][i] != '0')
                    dist++;

            return dist;
        }
        private void CalcAllDistances(List<string> structNames)
        {
            consensus.Clear();
            distanceOrdered.Clear();
            foreach (var item in structNames)
            {
                consensus.Add(item, TransformToConsensusStates(item));
                int dist = DistanceToConsensus(item);
                distanceOrdered.Add(item, dist);
            }

        }
        public void ToConsensusStates(List<string> structNames, string newConsensusStates)
        {
            List<string> states = new List<string>();
            if (!data.ContainsKey(newConsensusStates))
            {
                states=null;
                return;
            }

            consensusStates=data[newConsensusStates];

            CalcAllDistances(structNames);
        }
        public void ToConsensusStates(List<string> structNames, List<double> newConsensusStates)
        {
            consensusStates = newConsensusStates;
            CalcAllDistances(structNames);
        }
    }
}
