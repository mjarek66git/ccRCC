using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using phiClustCore;

namespace phiClustCore.Distance
{
    public abstract class HammingBase:DistanceMeasure
    {
        public Dictionary<string, List<double>> data;
        Settings dirSettings = new Settings();

        protected Dictionary<string, int> states = new Dictionary<string, int>();
        protected List<Dictionary<string, int>> lStates = new List<Dictionary<string, int>>();

        Dictionary<double, List<int>>[] columns = null;

        protected bool flag;
        public  HammingBase(bool flag)
        {
            this.flag = flag;
        }
        public HammingBase(Dictionary<string,List<double>> data, bool flag)
        {
            this.data = data;
            this.flag = flag;
        }

        public override void InitMeasure()            
        {
            order = true;
            InitHamming();
        }
        public override double ProgressUpdate()
        {
            return base.ProgressUpdate();
        }
        private void InitHamming()
        {
            base.InitMeasure(data, true);
           structNames = new Dictionary<string,int>();
            foreach (string item in data.Keys)
            {
                string[] strTab = item.Split(Path.DirectorySeparatorChar);
                structNames.Add(strTab[strTab.Length - 1],1);
            }

            order = true;
        }
        
        public override int[][] GetDistance(List<string> refStructure, List<string> structures)
        {
            return base.GetDistance(refStructure, structures);
        }
      

        protected Dictionary<double, List<int>>[] MakeColumnsLists(List<string> structNames)
        {
            double locState = 0;

            if (structNames.Count == 0)
                return null;

            columns = new Dictionary<double, List<int>>[data[structNames[0]].Count];
           
            for (int i = 0; i < columns.Length; i++)
            {

                columns[i] = new Dictionary<double, List<int>>();
                try
                {
                    for (int j = 0; j < structNames.Count; j++)
                    {
                        if (data.ContainsKey(structNames[j]) && data[structNames[j]].Count > 0)
                            locState = data[structNames[j]][i];
                        else
                            continue;
                        if (locState == 0)
                            continue;


                        if (!columns[i].ContainsKey(locState))
                        {
                            List<int> lista = new List<int>();
                            lista.Add(j);
                            columns[i].Add(locState, lista);
                        }
                        else
                            columns[i][locState].Add(j);
                    }
                }
                catch (Exception ex)
                {
                    Console.WriteLine("Ups HammingBase :" + ex.Message);
                }
            }

            return columns;
        }

    }
}
