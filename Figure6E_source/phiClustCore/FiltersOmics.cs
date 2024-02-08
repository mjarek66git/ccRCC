using phiClustCore.Distance;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
//using Accord.Statistics.Distributions.Univariate;

namespace phiClustCore
{
    public class Interval
    {
        public double min;
        public double max;
        public int code;

    }
    public class FileN
    {
        public string fileName;
    }
    public abstract class FilterOmics : ICloneable
    {
        public static List<OmicsDataSet> memoryFilteredData = new List<OmicsDataSet>();
        public static bool remData = false;
        protected string name = "";
        protected bool parameters = false;
        protected short cores = 5;

        protected readonly object lockObj = new object();
        public virtual string Name { get { return name; } }
        public bool Parameters { get { return parameters; } }
        public virtual Dictionary<string, Type> GetParameters() { return null; }
        public virtual void SetParameters(Dictionary<string, string> x) { }
        public virtual Dictionary<string, string> GetParametersValue() { return null; }
        public abstract OmicsDataSet ApplyFilter(OmicsDataSet data);

        public virtual void Clear()
        {

        }

        public object Clone()
        {
            var clone = this.MemberwiseClone();
            return clone;
        }

    }
    public class QuantileColumn : FilterOmics
    {
        public QuantileColumn()
        {
            name = ToString();
        }
        public override string ToString()
        {
            return "Quantile column normalization";
        }
        void QuantileThread(int[,] rankData, int[,] rank, SpareArray data, float[] avr, int start, int end)
        {
            Dictionary<double, int> valToRank = new Dictionary<double, int>(500);
            for (int i = start; i < end; i++)
            {
                valToRank.Clear();
                var colTuple = data.GetColumnValues(i);
                float[] columnValues = colTuple.Item1;
                int[] copyData = colTuple.Item2;

                Dictionary<double, int> diffValues = new Dictionary<double, int>();
                foreach (var item in columnValues)
                {
                    if (!diffValues.ContainsKey(item))
                        diffValues.Add(item, 0);

                    diffValues[item]++;

                }
                var val = diffValues.Keys.ToArray();

                int[] index = new int[diffValues.Count];
                for (int j = 0; j < index.Length; j++)
                    index[j] = j;

                Array.Sort<int>(index, (a, b) => val[a].CompareTo(val[b]));
                Array.Reverse(index);
                //index.Reverse();


                int k = 1;
                foreach (var item in index)
                {
                    valToRank.Add(val[item], k);
                    k += diffValues[val[item]];
                }

                for (int j = 0; j < copyData.Length; j++)
                    if (data[copyData[j], i] != 0)
                        lock (lockObj)
                        {
                            rank[copyData[j], i] = valToRank[data[copyData[j], i]];
                        }


                for (int j = 0; j < copyData.Length; j++)
                    if (data[copyData[j], i] != 0)
                        lock (lockObj)
                        {
                            rankData[copyData[j], i] = valToRank[data[copyData[j], i]]++;
                        }

                for (int j = 0; j < avr.Length; j++)
                {
                    if (rankData[j, i] > 0)
                    {
                        lock (lockObj)
                        {

                            avr[rankData[j, i] - 1] += data[j, i];
                            //count[j]++;
                        }
                    }
                }
            }

        }
        public SpareArray QuantileTransform(SpareArray data)
        {
            int[] copyData;
            SpareArray rankData;
            Dictionary<double, int> dic = new Dictionary<double, int>();
            float[] avr;

            //data.Print();

            avr = new float[data.rows];
            rankData = new SpareArray(data);
            int[,] rankAll = new int[data.rows, data.columns];
            int[,] rank = new int[data.rows, data.columns];

            Task[] t = new Task[cores];
            float step = data.columns / ((float)cores);

            for (int i = 0; i < cores; i++)
            {
                int s = (int)(step * i);
                int e = (int)(step * (i + 1));
                t[i] = Task.Run(() => QuantileThread(rankAll, rank, data, avr, s, e));
            }
            Task.WaitAll(t);

            for (int i = 0; i < avr.Length; i++)
                avr[i] /= data.columns;


            float[] vTab = new float[rankData.rows];
            int[] countI = new int[vTab.Length];
            for (int j = 0; j < rankData.columns; j++)
                {
                for (int i = 0; i < vTab.Length; i++)
                    {
                    vTab[i] = 0;
                    countI[i] = 0;
                    }
                for (int i = 0; i < rank.GetLength(0); i++)
                    if (rank[i, j] > 0)
                    {
                        vTab[rank[i, j] - 1] += avr[rankAll[i, j] - 1];
                        countI[rank[i, j] - 1]++;
                }

                for (int i = 0; i < rankData.rows; i++)
                    if (rank[i, j] > 0)
                        rankData[i, j] = vTab[rank[i, j] - 1] / countI[rank[i, j] - 1];

            }

            //rankData.Print();
            // rankData.Save("AfterQuantile.dat");
            return rankData;
        }

        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet(dataS.Name + "_" + ToString());

            res.data = QuantileTransform(dataS.data);
            res.geneLabels = dataS.geneLabels;
            res.sampleLabels = dataS.sampleLabels;

            if (remData)
                memoryFilteredData.Add(res);

            return res;
        }

    }
    public class QuantileRow : QuantileColumn
    {
        public QuantileRow()
        {
            name = ToString();
        }
        public override string ToString()
        {
            return "Quantile row normalization";
        }

        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet(dataS.Name + "_" + ToString());
            SpareArray tempData = dataS.data.Transpose();
            SpareArray d = QuantileTransform(tempData);
            res.data = d.Transpose();
            res.geneLabels = dataS.geneLabels;
            res.sampleLabels = dataS.sampleLabels;

            if (remData)
                memoryFilteredData.Add(res);
            //res.Save("Quantile");
            return res;
        }

    }
    public class BulkFilter : FilterOmics
    {
        public BulkFilter()
        {
            name = ToString();
        }
        public override string ToString()
        {
            return "Bulk normalization";
        }

        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet(dataS.Name + "_" + ToString());
            float[] hist = new float[dataS.data.columns];
            for (int i = 0; i < dataS.data.columns; i++)
            {
                float sum = 0;
                for (int j = 0; j < dataS.data.rows; j++)
                {
                    sum += dataS.data[j, i];
                }
                hist[i] = sum / dataS.data.rows;
            }
            int[] index = new int[hist.Length];

            res.data = new SpareArray(dataS.data.rows, dataS.data.columns, dataS.data.spareFlag);
            for (int i = 0; i < index.Length; i++)
                index[i] = i;

            Array.Sort(hist, index);
            Array.Reverse(index);
            Array.Reverse(hist);

            float[] tmp = new float[dataS.data.columns];
            int[] tmpIndex = new int[tmp.Length];
            for (int i = 0; i < dataS.data.rows; i++)
            {
                for (int j = 0; j < dataS.data.columns; j++)
                    tmpIndex[j] = j;
                for (int j = 0; j < dataS.data.columns; j++)
                    tmp[j] = dataS.data[i, j];

                Array.Sort(tmp, tmpIndex);
                Array.Reverse(tmpIndex);
                Array.Reverse(tmp);
                for (int j = 0; j < dataS.data.columns; j++)
                    if (tmp[j] > 0)
                        res.data[i, tmpIndex[j]] = hist[j];

            }

            res.geneLabels = dataS.geneLabels;
            res.sampleLabels = dataS.sampleLabels;

            if (remData)
                memoryFilteredData.Add(res);

            return res;
        }

    }

    public class HallmarkGenes : FilterOmics
    {
        public HallmarkGenes()
        {
            name = ToString();
        }
        public override string ToString()
        {
            return "Hallmark Genes";
        }
        HashSet<string> LoadHallmarkGenes()
        {
            HashSet<string> geneNames = new HashSet<string>();
            StreamReader st = new StreamReader(@"C:\projects\bioinfo\Cluster_patient\HALLMARK_test.txt");
            string line = st.ReadLine();
            while (line != null)
            {
                if (!line.Contains(">"))
                {
                    string[] aux = line.Split(' ');
                    foreach (var item in aux)
                        geneNames.Add(item);
                }
                line = st.ReadLine();
            }
            st.Close();

            return geneNames;
        }
        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            HashSet<string> validGenes = LoadHallmarkGenes();
            OmicsDataSet res = new OmicsDataSet(dataS.Name + "_" + ToString());

            res.data = new SpareArray(dataS.data.rows, validGenes.Count, dataS.data.spareFlag);

            int[] index = new int[validGenes.Count];
            res.geneLabels = new List<string>();
            for (int i = 0, k = 0; i < dataS.geneLabels.Count; i++)
            {
                if (validGenes.Contains(dataS.geneLabels[i]))
                {
                    index[k++] = i;
                    res.geneLabels.Add(dataS.geneLabels[i]);
                }
            }
            for (int i = 0; i < dataS.data.rows; i++)
                for (int j = 0; j < index.Length; j++)
                    res.data[i, j] = dataS.data[i, index[j]];

            res.sampleLabels = dataS.sampleLabels;

            if (remData)
                memoryFilteredData.Add(res);

            //memoryFilteredData.Add(res);

            return res;
        }

    }


    public class MinMaxAmountThreshold : FilterOmics
    {
        public double thresholdMin = 0;
        public double thresholdMax = 0;

        public MinMaxAmountThreshold()
        {
            name = ToString();
            parameters = true;
        }
        public override Dictionary<string, Type> GetParameters()
        {
            Dictionary<string, Type> res = new Dictionary<string, Type>();
            res.Add("Threshold Min", typeof(double));
            res.Add("Threshold Max", typeof(double));
            return res;
        }
        public override Dictionary<string, string> GetParametersValue()
        {
            Dictionary<string, string> v = new Dictionary<string, string>();
            v.Add("Threshold Min", thresholdMin.ToString());
            v.Add("Threshold Max", thresholdMax.ToString());

            return v;
        }

        public override void SetParameters(Dictionary<string, string> x)
        {
            thresholdMin = Convert.ToDouble(x["Threshold Min"]);
            thresholdMax = Convert.ToDouble(x["Threshold Max"]);
            name = ToString() + " " + thresholdMin + " " + thresholdMax;
        }
        public override string ToString()
        {
            return "Min Max amount threshold";
        }

        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet("MinMaxAmountThreshold");
            HashSet<int> toRemove = new HashSet<int>();
            /*            for (int i = 0; i < dataS.data.columns; i++)
                        {
                            int countW = 0;
                            for (int j = 0; j < dataS.data.rows; j++)
                            {
                                if (dataS.data[j, i] > 0)
                                    countW++;
                            }
                            if (countW < thresholdMin || countW>thresholdMax)
                                toRemove.Add(i);
                        }*/

            //            SpareArray data = new SpareArray(dataS.data.rows, dataS.data.columns - toRemove.Count,false);

            int[] freq = dataS.data.CountColumns();

            res.data = dataS.data.FastCopyColumnThreshold(freq, (int)thresholdMin, (int)thresholdMax);
            res.geneLabels = new List<string>();
            for (int i = 0; i < freq.Length; i++)
                if (freq[i] >= thresholdMin && freq[i] <= thresholdMax)
                    res.geneLabels.Add(dataS.geneLabels[i]);



            res.sampleLabels = dataS.sampleLabels;

            if (remData)
                memoryFilteredData.Add(res);

            return res;
        }

    }
    public class TopDevGenes : FilterOmics
                        {
        int threshold = 0;

        public TopDevGenes()
        {
            name = ToString();
            parameters = true;
        }
        public override Dictionary<string, Type> GetParameters()
        {
            Dictionary<string, Type> res = new Dictionary<string, Type>();
            res.Add("Threshold", typeof(int));
            return res;
        }
        public override void SetParameters(Dictionary<string, string> x)
        {
            threshold = Convert.ToInt32(x["Threshold"]);
            name = ToString() + " " + threshold;
        }
        public override Dictionary<string, string> GetParametersValue()
        {
            Dictionary<string, string> v = new Dictionary<string, string>();
            v.Add("Threshold", threshold.ToString());

            return v;
                        }
        public override string ToString()
        {
            return "Top best genes";
        }

        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet("TopBestGenes");
            SpareArray data = dataS.data;
            //List<int> omitGenes = new List<int>();

            var d = data.ColumnDataNumber();

            var omitGenes = Enumerable.Range(0, d.Length).Where(n => d[n] < 0.01 * data.rows).ToList();


            var colDevOrg = data.GetColDeviation();
            double[] colDev = new double[data.columns - omitGenes.Count()];
            if (colDev.Length == 0)
                throw new Exception("To many genes rejected");



            int[] index = new int[colDev.Length];
            for (int i = 0, k = 0; i < data.columns; i++)
                if (!omitGenes.Contains(i))
                        {
                    colDev[k] = colDevOrg[i];
                    index[k++] = i;

                }

            Array.Sort(colDev, index);
            Array.Reverse(index);
            if (index.Length < threshold)
                threshold = index.Length;

            int[] finalIndex = new int[threshold];

            Array.Copy(index, finalIndex, threshold);

            res.geneLabels = new List<string>();
            for (int i = 0; i < finalIndex.Length; i++)
                res.geneLabels.Add(dataS.geneLabels[finalIndex[i]]);

            res.data = data.CopyColumns(finalIndex);

            res.sampleLabels = dataS.sampleLabels;
            Interval[] intv = new Interval[2];

            intv[0] = new Interval();
            intv[0].code = 0;
            intv[1] = new Interval();
            intv[1].code = 1;
            res.intervals.Add(intv);

            if (remData)
                memoryFilteredData.Add(res);

            //res.Save("topGenes");

            return res;
        }

    }

    public class MaxActivationThreshold : FilterOmics
    {
        float threshold = 0;

        public MaxActivationThreshold()
        {
            name = ToString();
            parameters = true;
        }
        public override Dictionary<string, Type> GetParameters()
        {
            Dictionary<string, Type> res = new Dictionary<string, Type>();
            res.Add("Threshold", typeof(double));
            return res;
        }
        public override void SetParameters(Dictionary<string, string> x)
        {
            threshold = Convert.ToSingle(x["Threshold"]);
            name = ToString() + " " + threshold;
        }
        public override Dictionary<string, string> GetParametersValue()
        {
            Dictionary<string, string> v = new Dictionary<string, string>();
            v.Add("Threshold", threshold.ToString());

            return v;
        }
        public override string ToString()
        {
            return "Max value threshold";
        }

        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet("MaxValueThreshold");
            SpareArray data = dataS.data;
            res.data = new SpareArray(data.rows, data.columns, false);

            for (int j = 0; j < data.rows; j++)
            {
                for (int i = 0; i < data.columns; i++)
                {
                    res.data[j, i] = data[j, i];
                    if (res.data[j, i] >= threshold)
                    {
                        res.data[j, i] = threshold;
                    }
                }
            }



            res.geneLabels = dataS.geneLabels;
            res.sampleLabels = dataS.sampleLabels;
            Interval[] intv = new Interval[2];

            intv[0] = new Interval();
            intv[0].code = 0;
            intv[1] = new Interval();
            intv[1].code = 1;
            res.intervals.Add(intv);

            if (remData)
                memoryFilteredData.Add(res);

            return res;
        }

    }
    public class MinCountThresholdColumn : FilterOmics
    {
        float threshold = 0;

        public MinCountThresholdColumn()
        {
            name = ToString();
            parameters = true;
        }
        public override Dictionary<string, Type> GetParameters()
        {
            Dictionary<string, Type> res = new Dictionary<string, Type>();
            res.Add("Threshold", typeof(double));
            return res;
        }
        public override void SetParameters(Dictionary<string, string> x)
        {
            threshold = Convert.ToSingle(x["Threshold"]);
            name = ToString() + " " + threshold;
        }
        public override Dictionary<string, string> GetParametersValue()
        {
            Dictionary<string, string> v = new Dictionary<string, string>();
            v.Add("Threshold", threshold.ToString());

            return v;
        }
        public override string ToString()
        {
            return "Min count threshold Column";
        }
        private int CountTh(int startR, int endR, SpareArray data, bool[] index)
        {
            int count = 0;
            for (int j = startR; j < endR; j++)
            {
                double sum = 0;
                for (int i = 0; i < data.rows; i++)
                    if (data[i, j] > 0)
                        sum++;

                if (sum < threshold)
                {
                    count++;
                }
                else
                {
                    lock (lockObj)
                    {
                        index[j] = true;
                    }
                }
            }
            return count;
        }
        private void CopyData(int startR, int endR, SpareArray data, SpareArray outData, bool[] index, int[] co)
        {
            for (int j = startR; j < endR; j++)
            {
                if (index[j])
                {
                    int k;
                    k = co[j];

                    for (int i = 0; i < data.rows; i++)
                        outData[i, k] = data[i, j];

                }
            }
        }

        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet("MinCountThresholdColumn");
            SpareArray data = dataS.data;

            /*bool[] index = new bool[data.columns];
            int[] co = new int[data.columns];
            //res.geneLabels = new List<string>();
            //int count = 0;
            //for (int j = 0; j < data.columns; j++)
            //{
            //    double sum = 0;
            //    for (int i = 0; i < data.rows; i++)
            //        if (data[i, j] > 0)
            //            sum++;

            //    if (sum < threshold)
            //    {
            //        count++;
            //    }
            //    else
            //    {
            //        index[j] = true;
            //        res.geneLabels.Add(dataS.geneLabels[j]);
            //    }
            //}

            Task<int>[] t = new Task<int>[cores];

            float step = data.columns / ((float)cores);

            for (int i = 0; i < cores; i++)
            {
                int s = (int)(step * i);
                int e = (int)(step * (i + 1));
                t[i] = Task.Run(() => CountTh(s, e, data,index));
            }
            Task.WaitAll(t);
            int count = 0;
            for (int i = 0; i < cores; i++)
                count += t[i].Result;

            res.data = new SpareArray(data.rows, data.columns-count, data.spareFlag);
            for (int i = 0,counter=0; i < index.Length; i++)
                if (index[i])
                    co[i] = counter++;

            Task[] task = new Task[cores];
            for (int i = 0; i < cores; i++)
            {
                int s = (int)(step * i);
                int e = (int)(step * (i + 1));
                task[i] = Task.Run(() => CopyData(s, e, data, res.data,index,co));
            }
            Task.WaitAll(task);


                        //for (int j = 0; j < data.columns; j++)
                        //{
                        //    if (index[j])
                        //    {
                        //        for (int i = 0; i < data.rows; i++)
                        //            res.data[i,count] = data[i,j];
                        //        count++;
                        //        res.geneLabels.Add(dataS.geneLabels[j]);
                        //    }
                        //}
            
            */
            int[] freq = dataS.data.CountColumns();

            res.data = dataS.data.FastCopyColumnThreshold(freq, (int)threshold, int.MaxValue);
            res.geneLabels = new List<string>();
            for (int i = 0; i < freq.Length; i++)
                if (freq[i] >= threshold)
                    res.geneLabels.Add(dataS.geneLabels[i]);
            res.sampleLabels = dataS.sampleLabels;

            Interval[] intv = new Interval[2];

            intv[0] = new Interval();
            intv[0].code = 0;
            intv[1] = new Interval();
            intv[1].code = 1;
            res.intervals.Add(intv);
            if (remData)
                memoryFilteredData.Add(res);

            return res;
        }

    }

    public class MinCountThresholdRow : FilterOmics
    {
        float threshold = 0;

        public MinCountThresholdRow()
        {
            name = ToString();
            parameters = true;
        }
        public override Dictionary<string, Type> GetParameters()
        {
            Dictionary<string, Type> res = new Dictionary<string, Type>();
            res.Add("Threshold", typeof(double));
            return res;
        }
        public override void SetParameters(Dictionary<string, string> x)
        {
            threshold = Convert.ToSingle(x["Threshold"]);
            name = ToString() + " " + threshold;
        }
        public override Dictionary<string, string> GetParametersValue()
        {
            Dictionary<string, string> v = new Dictionary<string, string>();
            v.Add("Threshold", threshold.ToString());

            return v;
        }
        public override string ToString()
        {
            return "Min count threshold Row";
        }
        private int CountTh(int startR, int endR, SpareArray data)
        {
            int count = 0;
            for (int j = startR; j < endR; j++)
            {
                double sum = 0;
                for (int i = 0; i < data.columns; i++)
                    if (data[j, i] > 0)
                        sum++;

                if (sum < threshold)
                {
                    count++;
                }
            }
            return count;
        }


        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet("MinCountThresholdRow");
            SpareArray data = dataS.data;

            List<int> index = new List<int>();

            int count = 0;
            res.sampleLabels = new List<string>();
            for (int j = 0; j < data.rows; j++)
            {
                if (data.RowsDataNumber(j) < threshold)
                    count++;
                else
                {
                    index.Add(j);
                    res.sampleLabels.Add(dataS.sampleLabels[j]);
                }
            }
            if (count == 0)
                return dataS;
            res.data = data.CopyRows(index);
            count = 0;


            res.geneLabels = dataS.geneLabels;
            // res.sampleLabels = dataS.sampleLabels;
            Interval[] intv = new Interval[2];

            intv[0] = new Interval();
            intv[0].code = 0;
            intv[1] = new Interval();
            intv[1].code = 1;
            res.intervals.Add(intv);
            if (remData)
                memoryFilteredData.Add(res);

            return res;
        }

    }
    public class LoadSuperGenes : FilterOmics
    {
        public List<FileN> file = new List<FileN>();
        public override string Name {
            get
            {
                return ToString() + ":" + GetAllFiles();
            }
        }
        public Dictionary<string, HashSet<string>> superGenes = new Dictionary<string, HashSet<string>>();

        public LoadSuperGenes()
        {
            parameters = true;
        }
        public override Dictionary<string, Type> GetParameters()
        {
            Dictionary<string, Type> res = new Dictionary<string, Type>();
            for (int i = 0; i < file.Count; i++)
                res.Add("File name_" + i, typeof(string));
            return res;
        }
        string GetAllFiles()
        {
            string n = "";
            if (file.Count > 0)
                foreach (var item in file)
                    n += item.fileName + " ";

            return n;

        }
        public override void SetParameters(Dictionary<string, string> x)
        {
            foreach (var item in x.Keys)
            {
                FileN aux = new FileN();
                aux.fileName = x[item];
                file.Add(aux);
            }
            name = ToString() + " " + GetAllFiles();
        }
        public override Dictionary<string, string> GetParametersValue()
        {
            Dictionary<string, string> v = new Dictionary<string, string>();
            for (int i = 0; i < file.Count; i++)
                v.Add("File name_" + i, file[i].fileName);

            return v;
        }
        public override string ToString()
        {
            return "Super gene file";
        }
        Dictionary<string, HashSet<string>> ReadSuperGenes()
        {
            foreach (var itemFile in file)
            {
                StreamReader sr = new StreamReader(itemFile.fileName);
                string superGeneName = "";
                string line = sr.ReadLine();
                while (line != null)
                {
                    if (line.StartsWith(">"))
                        superGeneName = line.Substring(1);
                    else
                    {
                        string[] aux = line.Split(' ');
                        List<string> l = aux.ToList();
                        HashSet<string> genesHash = new HashSet<string>();
                        foreach (var item in l)
                            genesHash.Add(item);
                        /* else
                             throw new Exception("Uknown gene: " + item);*/
                        superGenes.Add(superGeneName, genesHash);
                    }

                    line = sr.ReadLine();
                }
                sr.Close();
            }
            return superGenes;
        }

        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            /*StreamWriter st = new StreamWriter("SuperIndex");
            foreach (var item in superGenes)
            {
                st.Write(item.Key + ": ");
                for (int i = 0; i < item.Value.Count; i++)
                    st.Write(item.Value[i] + " ");
                st.WriteLine();
            }
            st.Close();*/

            if (superGenes == null || superGenes.Count==0)
            {
                superGenes = ReadSuperGenes();
            }

            var res = dataS.CreateSuperGenesData(superGenes);

            Interval[] intv = new Interval[2];

            intv[0] = new Interval();
            intv[0].code = 0;
            intv[1] = new Interval();
            intv[1].code = 1;
            res.intervals.Add(intv);
            if (remData)
                memoryFilteredData.Add(res);

            //res.Save("SuperGenes");

            return res;
        }

    }
    public class SelectStage : FilterOmics
    {
        public string stageType = "";
        public SelectStage()
        {
            name = ToString();
            parameters = true;
        }
        public override string ToString()
        {
            return "Select samples based on stage";
        }
        public override Dictionary<string, Type> GetParameters()
        {
            Dictionary<string, Type> res = new Dictionary<string, Type>();
            res.Add("Set stage name", typeof(string));
            return res;
        }
        public override void SetParameters(Dictionary<string, string> x)
        {
            stageType = x["Set stage name"];
        }

        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet("Stage");
            SpareArray data = dataS.data;
            int count = 0;
            if (StaticDic.Id.Count == 0)
                throw new Exception("Stages not defined");
            for (int i = 0; i < dataS.sampleLabels.Count; i++)
                if (StaticDic.Id.ContainsKey(dataS.sampleLabels[i]))
                    if(StaticDic.Id[dataS.sampleLabels[i]].Equals(stageType))
                        count++;
            int[] index = new int[count];
            res.sampleLabels = new List<string>(count);
            count = 0;
            for (int i = 0; i < dataS.sampleLabels.Count; i++)
                if (StaticDic.Id.ContainsKey(dataS.sampleLabels[i]))
                    if (StaticDic.Id[dataS.sampleLabels[i]].Equals(stageType))
                    {
                        index[count] = i;
                        res.sampleLabels.Add(dataS.sampleLabels[i]);
                    }
            res.data = new SpareArray(index.Length, data.columns, false);
            for (int j = 0; j < index.Length; j++)
            {
                for (int i = 0; i < data.columns; i++)
                    res.data[j, i] = dataS.data[j, i];
            }

            res.geneLabels = dataS.geneLabels;
            if (remData)
                memoryFilteredData.Add(res);

            return res;
        }

    }

    public class ShiftDataToNonNegative : FilterOmics
    {

        public ShiftDataToNonNegative()
        {
            name = ToString();
        }
        public override string ToString()
        {
            return "Shift data to non negative values";
        }

        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet("ShiftNonNegative");
            SpareArray data = dataS.data;
            res.data = new SpareArray(data.rows, data.columns, false);

            double minValue = double.MaxValue;
            for (int j = 0; j < data.rows; j++)
            {             
                for (int i = 0; i < data.columns; i++)
                    if (data[j, i] <= minValue)
                        minValue = data[j, i];
            }
            for (int j = 0; j < data.rows; j++)
            {

                if (minValue < 0)
                {
                    for (int i = 0; i < data.columns; i++)
                        res.data[j, i] = (float)(data[j, i] + Math.Abs(minValue));
                }
                else
                    for (int i = 0; i < data.columns; i++)
                        res.data[j, i] = data[j, i];

            }
            
            res.geneLabels = dataS.geneLabels;
            res.sampleLabels = dataS.sampleLabels;
            if (remData)
                memoryFilteredData.Add(res);
            //res.Save("shifted");
            return res;
        }

    }

    public class Descritize : FilterOmics
    {
        [Description("Descret states based on column")]
        public bool codingColumn { get; set; } = true;
        CodingAlg coding;
        [Description("Number of descret states")]
        public int NumStates { get; set; }
        [Description("Coding algorithm")]
        public CodingAlg Coding { get { return coding; } set { coding = value; } }

        List<Interval[]> intervals = new List<Interval[]>();


        public Descritize()
        {
            name = ToString();
            parameters = true;
        }

        public override void Clear()
        {
            intervals.Clear();
        }

        public override Dictionary<string, Type> GetParameters()
        {
            Dictionary<string, Type> res = new Dictionary<string, Type>();

            res.Add("Number of states", typeof(int));

            res.Add("Coding algorithm", typeof(CodingAlg));
            res.Add("Based on columns", typeof(bool));
            return res;
        }
        public override void SetParameters(Dictionary<string, string> x)
        {
            NumStates = Convert.ToInt32(x["Number of states"]);
            Enum.TryParse<CodingAlg>(x["Coding algorithm"], out coding);
            codingColumn = Convert.ToBoolean(x["Based on columns"]);
            name = "Descritize " + NumStates + " " + Coding + " columns" + " " + codingColumn;
        }
        public override Dictionary<string, string> GetParametersValue()
        {
            Dictionary<string, string> v = new Dictionary<string, string>();
            v.Add("Number of states", NumStates.ToString());
            v.Add("Coding algorithm", coding.ToString());
            v.Add("Coding based on columns", codingColumn.ToString());
            return v;
        }

        public override string ToString()
        {
            return "Descritize";
        }
        Tuple<double, double> BinomialFit(double[] dat)
        {
            int[] frequency = new int[1000];
            int[] v = new int[dat.Length];
            for (int i = 0; i < dat.Length; i++)
            {
                v[i] = (int)(dat[i] * 1000);
                frequency[v[i]]++;
            }


            double mean;
            double sumF = 0;
            double sumFx = 0;
            int n = 0;
            for (int i = 0; i < frequency.Length; i++)
            {
                if (frequency[i] > 0)
                {
                    sumF += frequency[i];
                    sumFx = frequency[i] * i;
                    n++;
                }
            }
            mean = sumFx / sumF;
            double p = mean / n;
            double variance = mean * (1 - p);

            return new Tuple<double, double>(mean, variance);
        }

        public static Tuple<double, double, double[], double[]> NegativeBinomialFit(double[] dat, bool realValue, double maxV)
        {

            const float range = 10000;
            int[] frequency = new int[(int)range + 1];
            int[] v = new int[dat.Length];

            double max = int.MinValue;
            for (int i = 0; i < dat.Length; i++)
            {
                if (double.IsNaN(dat[i]))
                    Console.WriteLine();
                v[i] = (int)(dat[i] * range);
                
            }

            //frequency = new int[max+1];

            for (int i = 0; i < dat.Length; i++)
            {
                if(v[i]>=0 && v[i]<frequency.Length)
                frequency[v[i]]++;
            }
            int[] newFreq = new int[50];
            double[] realValues = new double[newFreq.Length];
            int step = 0;
            double s = frequency.Length / 49.0;
            for (double i = 0; i < frequency.Length; i += s, step++)
            {
                int sum = 0;
                double realV = 0;
                int count = 0;
                for (int j = (int)i; j < s * step + s; j++)
                {
                    sum += frequency[j];
                    realV += (j * maxV) / range;
                    count++;
                }
                newFreq[step] = sum;
                realValues[step] = realV / count;
            }
            frequency = newFreq;

            double mean;
            double sumF = 0;
            double sumFx = 0;
            double sumFx2 = 0;
            int n = 0;
            for (int i = 0; i < frequency.Length; i++)
            {
                if (frequency[i] > 0)
                {

                    double val;
                    if (realValue)
                        val = (i * maxV) / range;
                    else
                        val = i;
                    sumF += frequency[i];
                    sumFx += frequency[i] * i;
                    sumFx2 += frequency[i] * i * i;
                    n++;
                }
            }
            mean = sumFx / sumF;

            double deviation = sumFx2 / sumF - mean * mean;
            double p = mean / deviation;
            double q = 1 - p;
            double k = mean * p / q;
            double[] retriveFrequency = new double[frequency.Length];
            double px = Math.Pow(p, k);
            retriveFrequency[0] = px * sumF;
            for (int i = 1; i < frequency.Length; i++)
            {
                px = (i + k - 1) / (i) * q * px;
                retriveFrequency[i] = px * sumF;
            }
            double xx = realValues[(int)Math.Ceiling(mean)];
            double yy = realValues[(int)Math.Ceiling(Math.Sqrt(deviation))];
            return new Tuple<double, double, double[], double[]>(xx, yy, retriveFrequency, realValues);
        }
        public static Interval[] GetIntervals(double avr, double stdev, int NumStates)
        {
            Interval[] intervals = new Interval[NumStates];
            for (int n = 0; n < NumStates; n++)
                intervals[n] = new Interval();

            int s = 0;
            double st;
            if (NumStates == 3)
                st = stdev / (NumStates - 2);
            else
                st = stdev * 2 / (NumStates - 2);
            st /= 2;



            intervals[s].min = double.MinValue;
            intervals[s++].max = avr - st * (NumStates - 3) - st;
            for (int k = 0; k < NumStates - 2; k++, s++)
            {
                intervals[s].min = intervals[s - 1].max;
                intervals[s].max = intervals[s].min + 2 * st;
            }
            intervals[s].min = intervals[s - 1].max;
            intervals[s].max = double.MaxValue;


            return intervals;
        }
        SpareArray ZScoreCoding(SpareArray data, OmicsDataSet om)
        {

            SpareArray newData = new SpareArray(data.rows, data.columns, false);
            double[] col = new double[newData.rows];
            for (int i = 0; i < data.columns; i++)
            {

                double maxV = double.MinValue;
                for (int n = 0; n < data.rows; n++)
                {
                    col[n] = data[n, i];
                    if (maxV < col[n])
                        maxV = col[n];
                }
                //NormalDistribution b = new NormalDistribution();
                //BinomialDistribution b = new BinomialDistribution();  

                // NegativeBinomialDistribution b = new NegativeBinomialDistribution(1;

                double avr;
                double stdev;

                if(maxV==double.MinValue)
                {
                    for (int n = 0; n < col.Length; n++)                    
                                newData[n, i] = -1;
                                
                    continue;                            

                }
                //var h=BinomialFit(col);
                //Zróbmy normalizację col
                if (maxV != 0)
                {
                for (int n = 0; n < col.Length; n++)
                    col[n] /= maxV;

                var h = NegativeBinomialFit(col, true, maxV);

                for (int n = 0; n < col.Length; n++)
                    col[n] *= maxV;

                avr = h.Item1;
                stdev = h.Item2;
                }
                else
                {
                    avr = 0;
                    stdev = 0;
                }
                var intv = GetIntervals(avr, stdev, NumStates);

                intervals.Add(intv);

                AssignCodeToIntervals(avr, intv);
                om.AddCodes(intv);
                for (int n = 0; n < col.Length; n++)
                    for (int k = 0; k < intv.Length; k++)
                        if (col[n] > intv[k].min && col[n] <= intv[k].max)
                        {
                            newData[n, i] = intv[k].code;
                            break;
                        }
            }

            return newData;
        }
        void ApplyFilterThread(int s, int e)
        {

        }
        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet(dataS.Name + "_" + ToString());
            Dictionary<double, int> hashValues = new Dictionary<double, int>();
            SpareArray outData = null;
            SpareArray data = dataS.data;
            double max = Double.MinValue;
            Interval[] intv = null;
            if (Coding == CodingAlg.Z_SCORE)
                outData = ZScoreCoding(data, res);
            else
            {
                outData = new SpareArray(data.rows, data.columns, false);
                if (codingColumn)
                {

                    for (int i = 0; i < data.columns; i++)
                    {
                        hashValues.Clear();
                        for (int j = 0; j < data.rows; j++)
                        {

                            if (!hashValues.Keys.Contains(data[j, i]))
                                hashValues.Add(data[j, i], 1);
                            else
                                hashValues[data[j, i]]++;
                        }
                        intv = SetupIntervals(hashValues, NumStates, Coding);
                        var codedColumn = IntervalCodigPerGene(data, i, intv);
                        for (int k = 0; k < codedColumn.Length; k++)
                            outData[k, i] = codedColumn[k];
                        intervals.Add(intv);

                    }
                }
                else
                {
                    for (int i = 0; i < data.rows; i++)
                    {
                        hashValues.Clear();
                        for (int j = 0; j < data.columns; j++)
                        {

                            if (!hashValues.Keys.Contains(data[i, j]))
                                hashValues.Add(data[i, j], 1);
                            else
                                hashValues[data[i, j]]++;
                        }
                        intv = SetupIntervals(hashValues, NumStates, Coding);
                        var codedRow = IntervalCodigPerSample(data, i, intv);
                        for (int k = 0; k < codedRow.Length; k++)
                            outData[i, k] = codedRow[k];
                        intervals.Add(intv);


                    }
                }

            }
            Dictionary<double, int> freq = new Dictionary<double, int>();
            HashSet<int> toRemoveRows = new HashSet<int>();
            for (int i = 0; i < outData.rows; i++)
            {
                freq.Clear();
                for (int j = 0; j < outData.columns; j++)
                    if (freq.ContainsKey(outData[i, j]))
                        freq[outData[i, j]]++;
                    else
                        freq.Add(outData[i, j], 1);

                var maxValueKey = freq.Aggregate((x, y) => x.Value > y.Value ? x : y).Key;
                double xx = ((double)freq[maxValueKey]) / outData.columns;
                if (xx > 0.9)
                    toRemoveRows.Add(i);
            }
            Dictionary<string, double> dicTMP = new Dictionary<string, double>();
            HashSet<int> toRemoveColumns= new HashSet<int>();
            for (int i = 0; i < outData.columns; i++)
            {
                freq.Clear();
                for (int j = 0; j < outData.rows; j++)
                    if (freq.ContainsKey(outData[j,i]))
                        freq[outData[j,i]]++;
                    else
                        freq.Add(outData[j,i], 1);

                var maxValueKey = freq.Aggregate((x, y) => x.Value > y.Value ? x : y).Key;
                double xx = ((double)freq[maxValueKey]) / outData.rows;
                //dicTMP.Add(dataS.geneLabels[i], xx);
                if ((xx > 0.5 && maxValueKey==-1) || xx>0.9)
                {
                    toRemoveColumns.Add(i);
                }
            }


            if (toRemoveRows.Count > 0 || toRemoveColumns.Count>0)
            {
                SpareArray filteredData = new SpareArray(data.rows - toRemoveRows.Count, data.columns-toRemoveColumns.Count, false);
                res.sampleLabels = new List<string>();
                int k = 0;
                for (int i = 0; i < outData.rows; i++)
                {
                    if (!toRemoveRows.Contains(i))
                    {
                        int s = 0;
                        for (int j = 0; j < outData.columns; j++)
                            if(!toRemoveColumns.Contains(j))
                                filteredData[k, s++] = outData[i, j];

                        res.sampleLabels.Add(dataS.sampleLabels[i]);
                        k++;
                    }
                }
                res.data = filteredData;
                
                if(toRemoveColumns.Count>0)
                {
                    res.geneLabels = new List<string>();
                    for (int s = 0; s < dataS.geneLabels.Count; s++)
                        if(!toRemoveColumns.Contains(s))
                            res.geneLabels.Add(dataS.geneLabels[s]);
                }
                else
                    res.geneLabels = dataS.geneLabels;
            }
            else
            {
                res.data = outData;
                res.geneLabels = dataS.geneLabels;
                res.sampleLabels = dataS.sampleLabels;
            }
            res.intervals = intervals;
            res.AddCodes(intervals[0]);

            if (remData)
                memoryFilteredData.Add(res);
            //res.Save("Discr");
            //res.prev = dataS;
            return res;
        }
        public static void AssignCodeToIntervals(double avr, Interval[] intervals)
        {
            if (intervals.Length == 3)
            {
                intervals[0].code = -1;
                intervals[1].code = 0;
                intervals[2].code = 1;
            }
            else
                for (int i = 0; i < intervals.Length; i++)
                    intervals[i].code = i;


        }
        Interval[] SetupIntervals(Dictionary<double, int> dataValues, int numStates, CodingAlg coding)
        {
            if (numStates < 3)
                numStates = 3;
            Interval[] intervals = null;
            double max, min;
            max = double.MinValue;
            min = double.MaxValue;
            int n = 0;
            double step = 0;
            int i = 0;
            List<double> dlist = new List<double>(dataValues.Keys);
            dlist.Sort();
            min = dlist[0];
            max = dlist[dlist.Count - 1];


            switch (coding)
            {
                case CodingAlg.EQUAL_DIST:
                    intervals = new Interval[numStates];
                    step = (max - min) / numStates;
                    for (i = 0; i < numStates; i++)
                    {
                        intervals[i] = new Interval();
                        intervals[i].min = dlist[0] + i * step;
                        intervals[i].max = dlist[0] + (i + 1) * step;
                    }
                    intervals[0].min = double.MinValue;
                    intervals[intervals.Length - 1].max = double.MaxValue;
                    break;

                case CodingAlg.PERCENTILE:
                    int counter = 0;
                    double end = 0;

                    foreach (var item in dataValues.Values)
                        counter += item;
                    int amount = numStates * counter / 100;
                    counter = 0;
                    int k = 0;
                    intervals = new Interval[3];
                    for (int s = 0; s < intervals.Length; s++)
                        intervals[s] = new Interval();
                    double begin = dlist[0];
                    for (int s = 0; s < dlist.Count; s++)
                    {
                        counter += dataValues[dlist[s]];
                        if (counter >= amount)
                        {
                            int c = s + 1;
                            if (c < dlist.Count)
                                intervals[0].max = (dlist[s] + dlist[c])/2;
                            else
                                intervals[0].max = dlist[s];
                            intervals[0].min = double.MinValue;
                            break;
                        }
                    }
                    counter = 0;
                    for (int s = dlist.Count - 1; s > 0; s--)
                    {
                        counter += dataValues[dlist[s]];
                        if (counter > amount)
                        {
                            intervals[1].min = intervals[0].max;
                            int c = s-1;
                            if (c >=0)
                                intervals[1].max = (dlist[s] +dlist[c])/2;
                            else
                                intervals[1].max = dlist[s];
                            intervals[2].min = intervals[1].max;
                            intervals[2].max = double.MaxValue;
                            break;
                        }
                    }

                    break;


            }

            AssignCodeToIntervals(0, intervals);

            return intervals;

        }
        int[] IntervalCodigPerGene(SpareArray data, int num, Interval[] intv)
        {
            int[] codedData = new int[data.rows];
            for (int j = 0; j < data.rows; j++)
            {
                for (int k = 0; k < intv.GetLength(0); k++)
                {
                    if (data[j, num] >= intv[k].min && data[j, num] < intv[k].max)
                    {
                        codedData[j] = intv[k].code;
                        break;
                    }
                }


            }

            return codedData;
        }
        int[] IntervalCodigPerSample(SpareArray data, int num, Interval[] intv)
        {
            int[] codedData = new int[data.columns];
            for (int j = 0; j < data.columns; j++)
            {
                for (int k = 0; k < intv.GetLength(0); k++)
                {
                    if (data[num, j] >= intv[k].min && data[num, j] < intv[k].max)
                    {
                        codedData[j] = intv[k].code;
                        break;
                    }
                }


            }

            return codedData;
        }

    }
    public class GlobalZScore : FilterOmics
    {
        double[] dev;
        public double[] columns;

        public GlobalZScore()
        {
            name = ToString();
        }

        public override string ToString()
        {
            return "Global Z-score";
        }

        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            float sumAll = 0;
            OmicsDataSet res = new OmicsDataSet(dataS.Name + "_" + ToString());
            SpareArray data = (SpareArray)this.Clone();

            columns = new double[data.columns];

            for (int i = 0; i < data.rows; i++)
                for (int j = 0; j < data.columns; j++)
                    sumAll += data[i, j];

            for (int i = 0; i < data.columns; i++)
                for (int j = 0; j < data.rows; j++)
                {
                    data[j, i] /= sumAll;
                    columns[i] += data[j, i];
                }
            int[] index = new int[columns.Length];

            for (int i = 0; i < index.Length; i++)
            {
                index[i] = i;
                columns[i] /= data.rows;
            }


            /* double []locColumns =(double []) columns.Clone();
             Array.Sort(locColumns, index);

             StreamWriter cc = new StreamWriter("ccc");
             for (int i = 0; i < dataS.geneLabels.Count; i++)
                 cc.WriteLine(dataS.geneLabels[index[i]] + " " + locColumns[i]);
             cc.Close();*/

            //for (int i = 0; i < columns.Length; i++)
            //    columns[i] /= sumAll;
            res.geneLabels = dataS.geneLabels;
            res.sampleLabels = dataS.sampleLabels;
            res.data = data;
            if (remData)
                memoryFilteredData.Add(res);
            return res;
        }
    }
    public class ZScoreColumn : FilterOmics
    {
        double[] dev;
        double[] avr;

        public ZScoreColumn()
        {
            name = ToString();
        }

        public override string ToString()
        {
            return "Z-score column";
        }
        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet(dataS.Name + "_" + ToString());
            SpareArray data = (SpareArray)dataS.data.Clone();

            double[] xx = new double[data.rows];
            double[] weight = new double[xx.Length];
            dev = new double[data.columns];
            avr = new double[data.columns];

            for (int i = 0; i < data.columns; i++)
            {
                for (int j = 0; j < data.rows; j++)
                    xx[j] = data[j, i];

                var h = Descritize.NegativeBinomialFit(xx, true, 1);
                //BinomialDistribution b = new BinomialDistribution();
                //b.Fit(xx);
                //b.Fit(xx, weight);
                dev[i] = h.Item2;
                avr[i] = h.Item1;
                //b.Variance
            }

            for (int i = 0; i < data.columns; i++)
            {
                for (int j = 0; j < data.rows; j++)
                    if (data[j, i] != double.NaN && dev[i] > 0)
                        data[j, i] = (float)((data[j, i] - avr[i]) / dev[i]);
            }

            res.data = data;
            res.geneLabels = dataS.geneLabels;
            res.sampleLabels = dataS.sampleLabels;
            if (remData)
                memoryFilteredData.Add(res);
            return res;
        }

    }

    public class ZScoreRow : FilterOmics
    {

        float[] dev;
        float[] avr;

        public ZScoreRow()
        {
            name = ToString();
        }
        public override string ToString()
        {
            return "Z-score row";
        }

        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            float sumX = 0;
            float sumX2 = 0;
            OmicsDataSet res = new OmicsDataSet(dataS.Name + "_" + ToString());
            SpareArray data = (SpareArray)dataS.data.Clone();
            dev = new float[data.rows];
            avr = new float[data.rows];
            Dictionary<double, int> freqV = new Dictionary<double, int>();
            int counterP = 0;
            for (int i = 0; i < data.rows; i++)
            {
                sumX = 0;
                sumX2 = 0;
                freqV.Clear();
                int counter = 0;
                for (int j = 0; j < data.columns; j++)
                    if (data[i, j] != double.NaN)
                    {
                        sumX += data[i, j];
                        sumX2 += data[i, j] * data[i, j];
                        counter++;

                        if (!freqV.ContainsKey(data[i, j]))
                            freqV.Add(data[i, j], 0);

                        freqV[data[i, j]]++;
                    }

                int xx = int.MinValue;
                foreach(var item in freqV)
                {
                    if (item.Value > xx)
                        xx = item.Value;
                }
                if (xx > data.columns / 2)
                    counterP++;

                if (counter > 0)
                {
                    sumX /= counter;
                    sumX2 /= counter;
                    avr[i] = sumX;
                    dev[i] = (float)Math.Sqrt(sumX2 - avr[i] * avr[i]);
                }
            }


            for (int i = 0; i < data.rows; i++)
            {

                for (int j = 0; j < data.columns; j++)
                    if (data[i, j] != double.NaN  && dev[i] > 0)
                        data[i, j] = (data[i, j] - avr[i]) / dev[i];
            }
            res.data = data;
            res.geneLabels = dataS.geneLabels;
            res.sampleLabels = dataS.sampleLabels;
            if (remData)
                memoryFilteredData.Add(res);

// res.Save("NormalizationZScore");

            return res;
        }

    }
    public class RowNormalizationLog : FilterOmics
    {

        public RowNormalizationLog()
        {
            name = ToString();
        }
        public override string ToString()
        {
            return "Row normalization log";
        }

        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet(dataS.Name + "_" + ToString());

            SpareArray data = new SpareArray(dataS.data.rows, dataS.data.columns, dataS.data.spareFlag);
            res.data = dataS.data.NormalizeRowsLog();
            res.geneLabels = dataS.geneLabels;
            res.sampleLabels = dataS.sampleLabels;
            if (remData)
                memoryFilteredData.Add(res);
            //res.Save("NormalizationLog");
            return res;
        }

    }

    public class RowNormalization : FilterOmics
    {


        public RowNormalization()
        {
            name = ToString();
        }
        public override string ToString()
        {
            return "Row normalization - normalized counts";
        }
        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet(dataS.Name + "_" + ToString());

            res.data = dataS.data.NormalizeRows();
            res.geneLabels = dataS.geneLabels;
            res.sampleLabels = dataS.sampleLabels;
            if (remData)
                memoryFilteredData.Add(res);


            //res.Save("Normalization");

            return res;
        }

    }
    public class RowNormalizationGeneSelected : FilterOmics
    {
        string geneName = "";
        public RowNormalizationGeneSelected()
        {
            name = ToString();
            parameters = true;
        }
        public override Dictionary<string, Type> GetParameters()
        {
            Dictionary<string, Type> res = new Dictionary<string, Type>();

            res.Add("Gene name", typeof(string));

            return res;
        }
        public override void SetParameters(Dictionary<string, string> x)
        {
            geneName = x["Gene name"];
            name = "Gene name " + geneName;
        }
        public override Dictionary<string, string> GetParametersValue()
        {
            Dictionary<string, string> v = new Dictionary<string, string>();
            v.Add("Gene name", geneName);
            return v;
        }


        public override string ToString()
        {
            return "Row normalization based on selected genes- normalized counts";
        }
        public override OmicsDataSet ApplyFilter(OmicsDataSet dataS)
        {
            OmicsDataSet res = new OmicsDataSet(dataS.Name + "_" + ToString());


            HashSet<int> index = new HashSet<int>();
            for (int i = 0; i < dataS.geneLabels.Count; i++)
                if (dataS.geneLabels[i].StartsWith(geneName))
                    index.Add(i);


            res.data = dataS.data.NormalizeRows(index);
            res.geneLabels = dataS.geneLabels;
            res.sampleLabels = dataS.sampleLabels;
            if (remData)
                memoryFilteredData.Add(res);


            //res.Save("Normalization");

            return res;
        }
    }
 }