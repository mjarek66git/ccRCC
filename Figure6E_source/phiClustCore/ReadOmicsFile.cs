using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore
{
    public class OmicsDataSetSetup
    {
        public char[] separators;

        public bool columnFlag;

        public int rowPos;
        public int colPos;


        public int startCol;
        public int labelRow;
      
    }

    public class ReadOmicsFile
    {
        Dictionary<string, float[]> locData = new Dictionary<string, float[]>(3000);
        public List<string> labRow = new List<string>(1000);
        public List<string> labCol = new List<string>(1000);

        string[] keyTab;
        float[][] tabNum;

        OmicsDataSetSetup setup;

        public ReadOmicsFile(OmicsDataSetSetup setup)
        {
            this.setup = setup;
        }

        void ParseLines(int start, int end, string[] tabString)
        {

            string[] aux;
            for (int l = start; l < end; l++)
            {
                string line = tabString[l].Replace("\"", "");
                if (l == setup.labelRow)
                {
                    aux = line.Split(setup.separators);
                    lock (labRow)
                    {
                        for (int i = setup.colPos; i < aux.Length; i++)
                            labRow.Add(aux[i]);
                    }
                }

                if (l >= setup.rowPos)
                {
                    try
                    {


                        aux = line.Split(setup.separators);
                        float[] d = new float[aux.Length - setup.colPos];
                        NumberFormatInfo provider = new NumberFormatInfo();
                        provider.NumberDecimalSeparator = ".";
                        for (int i = setup.colPos, k = 0; i < aux.Length; i++, k++)
                            d[k] = Convert.ToSingle(aux[i], provider);

                        lock (tabNum)
                        {
                            tabNum[l] = d;
                        }
                        lock (keyTab)
                        {
                            keyTab[l] = aux[setup.startCol];
                        }
                    }
                    catch(Exception ex)
                    {
                        Console.WriteLine("Ups: " + ex.Message);
                    }

                }
            }
        }


        public OmicsDataSet GetOmicsFile(string fileName)
        {
            string[] tabString = new string[20000];
            StreamReader sr = new StreamReader(fileName, Encoding.UTF8, true, 1024);
            var data = new OmicsDataSet(Path.GetFileName(fileName));
            //locData.Clear();
//            Dictionary<string, float[]> locData = new Dictionary<string, float[]>(3000);

            labRow = new List<string>(1000);
            labCol = new List<string>(1000);
            string line = sr.ReadLine();
            string[] aux = null;
            int count = 0;
            tabString[0] = line;
            int x = 1;
            while (!sr.EndOfStream)
            {
                tabString[x] = sr.ReadLine();
                x++;
            }
            int cores = 5;
            Task[] tabTask = new Task[cores];

            keyTab = new string[x];
            tabNum = new float[x][];

            float step = x / ((float)cores);
            for (int i = 0; i < cores; i++)
            {
                int s = (int)step * i;
                int k = (int)step * (i + 1);

                tabTask[i] = Task.Run(() => ParseLines(s, k, tabString));
            }
            Task.WaitAll(tabTask);

            sr.Close();
            List<string> locKeys = new List<string>(locData.Keys);
            int rem = -1;
            for (int i = 0; i < keyTab.Length; i++)
                if (keyTab[i] != null)
                {
                    locKeys.Add(keyTab[i]);
                    labCol.Add(keyTab[i]);
                    rem = i;
                }

            if (setup.columnFlag)
            {
                data.geneLabels = labCol;
                data.sampleLabels = labRow;
            }
            else
            {
                data.geneLabels = labRow;
                data.sampleLabels = labCol;
            }


            //Debug.WriteLine("samples="+data.sampleLabels.Count);
            data.data = new SpareArray(data.sampleLabels.Count, data.geneLabels.Count, false);
            if (labCol.Contains(data.geneLabels[0]))
            {
                for (int i = 0; i < tabNum[rem].Length; i++)
                {
                    int counter = 0;
                    for (int j = 0; j < tabNum.Length; j++)
                    {
                        if (tabNum[j] == null)
                            continue;
                        data.data[i, counter++] = tabNum[j][i];
                    }
                }

            }
            else
            {
                for (int j = 0; j < locKeys.Count; j++)
                {
                    if (tabNum[j] != null)
                        for (int i = 0; i < tabNum[rem].Length; i++)
                            if (locData[locKeys[j]][i] > 0)
                                data.data[j, i] = tabNum[j][i];
                }

            }


            return data;
        }


        static public OmicsDataSet GetOmicsSpareFile(string fileName)
        {
            string[] tabString = new string[20000];
            StreamReader sr = new StreamReader(fileName, Encoding.UTF8, true, 1024);
            OmicsDataSet data = new OmicsDataSet(Path.GetFileName(fileName));
            //locData.Clear();

            string line = sr.ReadLine();
            string[] aux = null;

            tabString[0] = line;
            int x = 1;
            while (!sr.EndOfStream)
            {
                tabString[x] = sr.ReadLine();
                x++;
            }
            sr.Close();
            data.sampleLabels = new List<string>(1000);
            aux = tabString[0].Split('\t');
            for (int i = 1; i < aux.Length; i++)
                data.sampleLabels.Add(aux[i]);


            SpareArray sp = new SpareArray(data.sampleLabels.Count);

            data.geneLabels = new List<string>(20000);
            for (int i = 1; i < x; i++)
            {
                aux = tabString[i].Split('\t');
                data.geneLabels.Add(aux[0]);

                for (int j = 1; j < aux.Length; j++)
                {
                    if (aux[j].Contains(","))
                    {
                        string[] aux2 = aux[j].Split(',');
                        int index = int.Parse(aux2[0]);
                        sp[index, i - 1] = int.Parse(aux2[1]);
                    }
                }
            }

            data.data = sp;
            data.data.columns = data.geneLabels.Count;

            return data;
        }


    }
}
