using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using phiClustCore;
using phiClustCore.Interface;
using System.Data;
using phiClustCore.Distance;

namespace phiClustCore
{
    class Jury3D:IProgressBar
    {
        DistanceMeasure dMeasure;
        int currentV, maxV;
        int progressRead = 0;

        double startProgress = 0;
        double endProgress = 1;

        public double StartProgress { set { startProgress = value; } get { return startProgress; } }
        public double EndProgress { set { endProgress = value; } get { return endProgress; } }

        public Jury3D(DistanceMeasure dMeasure)
        {
            this.dMeasure = dMeasure;
            maxV = 1;
            currentV = 0;
            
        }
        public override string ToString()
        {
            return "3DJury";
        }
        public double ProgressUpdate()
        {
            double sumProgress = 0;
           double progress = dMeasure.ProgressUpdate();

           if (progressRead == 1)
               sumProgress = 0.05 + progress * 0.7;
           else
               sumProgress = 0.05 * progress;

           return StartProgress+(EndProgress - StartProgress) * (sumProgress + 0.25 * ((double)currentV / maxV));
        }
        public Exception GetException()
        {
            return null;
        }
        public List<KeyValuePair<string, DataTable>> GetResults()
        {
            return null;
        }
        public ClusterOutput Run3DJury()
        {
            ClusterOutput output = new ClusterOutput();
            List<KeyValuePair<string, double>> li = new List<KeyValuePair<string, double>>();

            long[] distTab = new long[dMeasure.structNames.Count];

            progressRead = 1;

            dMeasure.CalcDistMatrix(new List <string>(dMeasure.structNames.Keys));

            maxV = dMeasure.structNames.Count + 1 ;

            for(int i=0;i<dMeasure.structNames.Count;i++)
            {
                long sum=0;
                for(int j=0;j<dMeasure.structNames.Count;j++)
                {
                    sum += dMeasure.GetDistance(i, j);
                }
                distTab[i] = sum;
                currentV++;
            }

            KeyValuePair<string, double> v;

            List<string> structKeys = new List<string>(dMeasure.structNames.Keys);
            for(int m=0;m<structKeys.Count;m++)
            {
                v = new KeyValuePair<string, double>(structKeys[m], (double)(distTab[m] / (100.0 * dMeasure.structNames.Count)));
                li.Add(v);
            }
            if (dMeasure.order == false)
            {
                li.Sort((firstPair, nextPair) =>
                {
                    return nextPair.Value.CompareTo(firstPair.Value);
                });
            }
            else
                li.Sort((firstPair, nextPair) =>
                {
                    return firstPair.Value.CompareTo(nextPair.Value);
                });

            output.juryLike=li;

            currentV = maxV;
            output.runParameters = "Distance measure: " + this.dMeasure;
            return output;
        }

    }
}
