using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore
{
    public class SpareArray : ICloneable
    {
        float[,] data = null;
        Dictionary<int, float>[] dicData = null;
        public bool spareFlag { get; } = false;
        public int rows { get; }
        public int columns { get; set; }

        readonly object lockObj = new object();

        public SpareArray(SpareArray source)
        {
            columns = source.columns;
            rows = source.rows;
            spareFlag = source.spareFlag;
            if (spareFlag)
            {
                dicData = new Dictionary<int, float>[source.dicData.Length];
            for (int i = 0; i < dicData.Length; i++)
                    dicData[i] = new Dictionary<int, float>(source.dicData[i].Count);
            }
            else
                data = new float[rows, columns];


        }

        public SpareArray(int row)
        {
            rows = row;
            spareFlag = true;
            dicData = new Dictionary<int, float>[rows];
            for (int i = 0; i < rows; i++)
                dicData[i] = new Dictionary<int, float>(300);
        }
        public void PrepareRow(int row, int count)
        {
            dicData[row] = new Dictionary<int, float>(count);
        }
        public SpareArray(int row, int col, bool spare)
        {
            int dicCount = 300;
            rows = row;
            columns = col;

            if (col < 300)
                dicCount = col;

            spareFlag = spare;
            if (!spare)
                data = new float[row, col];
            else
            {
                dicData = new Dictionary<int, float>[rows];
                for (int i = 0; i < dicData.Length; i++)
                    dicData[i] = new Dictionary<int, float>(dicCount);

            }

        }
        public SpareArray(float[,] data, bool flag = false)
        {
            spareFlag = flag;
            rows = data.GetLength(0);
            columns = data.GetLength(1);
            if (spareFlag)
            {
                dicData = new Dictionary<int, float>[rows];
                for (int i = 0; i < dicData.Length; i++)
                    dicData[i] = new Dictionary<int, float>(columns);

                for (int i = 0; i < data.GetLength(0); i++)
                    for (int j = 0; j < data.GetLength(1); j++)
                        this[i, j] = data[i, j];
            }
            else
                this.data = data;

        }
        public int RowsDataNumber(int row)
        {
            int count = 0;
            if (spareFlag)
                return dicData[row].Count;
            else
            {
                for (int i = 0; i < columns; i++)
                    if (data[row, i] > 0)
                        count++;
            }

            return count;
        }
        public Tuple<float[], int[], int> GetColumnValues(int col)
        {
            float[] res;
            int[] index;
            int nonZeroCount = 0;
            if (spareFlag)
            {
                //res = new double[ColumnDataNumber(col)];
                res = new float[rows];
                index = new int[res.Length];

                for (int i = 0, j = 0; i < dicData.Length; i++)
                    if (dicData[i].ContainsKey(col))
                    {
                        res[j] = dicData[i][col];
                        index[j++] = i;
                        nonZeroCount++;
                    }
                    else
                    {
                        res[j] = 0;
                        index[j++] = i;
                    }



                return new Tuple<float[], int[], int>(res, index, nonZeroCount);
            }
            res = new float[rows];
            index = new int[res.Length];
            for (int i = 0, k = 0; i < data.GetLength(0); i++)
            {
                //if (data[i, col] != 0)
                {
                    res[k] = data[i, col];
                    index[k++] = i;
                    if (data[i, col] != 0)
                        nonZeroCount++;
                }
            }



            return new Tuple<float[], int[], int>(res, index, nonZeroCount);

        }
        public double[] GetColDeviation()
        {
            double[] avr = new double[columns];
            double[] dev = new double[columns];
            if (spareFlag)
            {
                int[] counts = new int[columns];
                for (int i = 0; i < dicData.Length; i++)
                    foreach (var item in dicData[i].Keys)
                    {
                        avr[item] += dicData[i][item];
                        counts[item]++;
                    }

                for (int i = 0; i < avr.Length; i++)
                    if (counts[i] > 0)
                        avr[i] /= counts[i];
                for (int i = 0; i < dicData.Length; i++)
                    foreach (var item in dicData[i].Keys)
                        dev[item] += Math.Pow(avr[item] - dicData[i][item], 2);

                for (int i = 0; i < dev.Length; i++)
                    //dev[i] /= counts[i];
                    dev[i] /= rows;

            }
            else
            {
                for (int j = 0; j < data.GetLength(1); j++)
                    for (int i = 0; i < data.GetLength(0); i++)
                        avr[j] += data[i, j];
                for (int i = 0; i < avr.Length; i++)
                    avr[i] /= data.GetLength(0);

                for (int j = 0; j < data.GetLength(1); j++)
                    for (int i = 0; i < data.GetLength(0); i++)
                        dev[j] += Math.Pow(avr[j] - data[i, j], 2);

                for (int i = 0; i < dev.Length; i++)
                    dev[i] /= rows;
            }

            return dev;
        }
        public double[] GetRowValues(int row)
        {
            double[] res;

            if (spareFlag)
            {
                res = new double[RowsDataNumber(row)];
                int i = 0;
                foreach (var item in dicData[row])
                    res[i++] = item.Value;
                return res;
            }            
            res = new double[columns];
            for (int i = 0; i < data.GetLength(1); i++)
                res[i] = data[row, i];

            return res;
        }
        public int[] ColumnDataNumber()
        {
            int[] counter = new int[columns];

            if (spareFlag)
            {
                for (int i = 0; i < dicData.Length; i++)
                {
                    foreach (var item in dicData[i].Keys)
                        counter[item]++;
                }

            }
            else
            {
                for (int j = 0; j < data.GetLength(1); j++)
                    for (int i = 0; i < data.GetLength(0); i++)
                        if (data[i, j] != 0)
                            counter[j]++;

            }

                return counter;
            }
        public SpareArray CopyRows(List<int> index)
        {
            SpareArray res = new SpareArray(index.Count, columns, spareFlag);
            if (spareFlag)
            {
                int k = 0;
                foreach (var item in index)
                    res.dicData[k++] = new Dictionary<int, float>(dicData[item]);
            }
            else
                for (int i = 0; i < index.Count; i++)
                {
                    for (int j = 0; j < data.GetLength(1); j++)
                        res.data[i, j] = data[index[i], j];
                   
                }
            return res;
        }
        public SpareArray CopyColumns(int[] index)
        {
            SpareArray res = new SpareArray(rows, index.Length, spareFlag);

            if (spareFlag)
            {
                Dictionary<int, int> h = new Dictionary<int, int>(index.Length);
                for (int i = 0; i < index.Length; i++)
                    h.Add(index[i], i);
                for (int i = 0; i < dicData.Length; i++)
                {

                    int count = 0;
                    foreach (var item in dicData[i])
                    {
                        if (h.ContainsKey(item.Key))
                            count++;
        }
                    res.dicData[i] = new Dictionary<int, float>(count);
                    foreach (var item in dicData[i])
                    {
                        if (h.ContainsKey(item.Key))
                            res.dicData[i].Add(h[item.Key], item.Value);

                    }
                }
            }
            else
            {
                for (int i = 0; i < data.GetLength(0); i++)
                {
                    for (int j = 0; j < index.Length; j++)
                        res.data[i, j] = data[i, index[j]];
                }
            }
            return res;
        }

        public SpareArray Transpose()
        {
            SpareArray trans;
            trans = new SpareArray(columns, rows, spareFlag);

            if (spareFlag)
            {
                for (int i = 0; i < dicData.Length; i++)
                {
                    foreach (var it in dicData[i])
                    {
                        trans.dicData[it.Key].Add(i, it.Value);
                    }
                }

            }
            else
            {
                trans.data = new float[data.GetLength(1), data.GetLength(0)];
                for (int i = 0; i < data.GetLength(0); i++)
                    for (int j = 0; j < data.GetLength(1); j++)
                        trans.data[j, i] = data[i, j];
            }


            return trans;
        }

        public int[] CountColumns()
        {
            int[] freq = new int[columns];

            if (spareFlag)
            {
                for (int i = 0; i < rows; i++)
                    foreach (var item in dicData[i])
                    {
                        freq[item.Key]++;
                    }
            }
            else
            {
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < columns; j++)
                        if (data[i, j] != 0)
                            freq[j]++;
                //throw new Exception("Not implemenetd yet");
            }
            return freq;
        }
        public SpareArray FastCopyColumnThreshold(int[] freq, int thresholdLow, int thresholdHeight)
        {
            SpareArray res = null;
            HashSet<int> validCoumns = new HashSet<int>();
            for (int i = 0; i < freq.Length; i++)
            {
                if (freq[i] >= thresholdLow && freq[i] <= thresholdHeight)
                    validCoumns.Add(i);
            }
            res = new SpareArray(rows, validCoumns.Count, spareFlag);

            if (spareFlag)
            {                
                res.dicData = new Dictionary<int, float>[dicData.Length];
                for (int i = 0; i < dicData.Length; i++)
                {
                    res.dicData[i] = new Dictionary<int, float>();
                    foreach (var item in dicData[i])
                    {
                        if (validCoumns.Contains(item.Key))
                            res.dicData[i].Add(item.Key, item.Value);
                    }

                }
            }
            else
            {
                for(int i=0;i<rows;i++)
                    for(int j=0;j<validCoumns.Count;j++)
                    {
                        res.data[i, j] = data[i, j];
                    }
            }
            return res;
        }

        void NormalizeRows(int start, int end, SpareArray res,HashSet<int> index)
        {
            if (spareFlag)
            {
                for (int i = start; i < end; i++)
                {
                    float sum = 0;
                    if (index == null)
                    {
                        foreach (var item in dicData[i])
                            sum += item.Value;
                    }
                    else
                        foreach (var item in dicData[i])
                            if(index.Contains(item.Key))
                                sum += item.Value;



                    foreach (var item in dicData[i])
                        res[i, item.Key] = item.Value / sum;
                }
            }
            else
            {
                for (int i = start; i < end; i++)
                {
                    float sum = 0;
                    if (index == null)
                    {
                        for (int j = 0; j < columns; j++)
                            sum += this[i, j];
                    }
                    else
                        foreach(var item in index)
                            sum += this[i, item];
                    // lock (res.values) lock (res.valuePos)

                    for (int j = 0; j < columns; j++)
                        res[i, j] = this[i, j] / sum;

                }
            }

        }
        void NormalizeRowsLog(int start, int end, SpareArray res)
        {
            if (spareFlag)
            {
                for (int i = start; i < end; i++)
                {
                    res.dicData[i] = new Dictionary<int, float>(dicData[i].Count);
                    double sum = 0;

                    foreach (var item in dicData[i])
                        res[i, item.Key] = (float)Math.Log(item.Value * 10000 + 1);
                }
            }
            else
            {
                for (int i = start; i < end; i++)
                {
                    for (int j = 0; j < columns; j++)
                        res[i, j] = (float)Math.Log(this[i, j] * 10000 + 1);
                }
            }

        }
        public SpareArray NormalizeRowsLog()
        {
            SpareArray res = null;
            res = new SpareArray(this);

            int cores = 5;
            Task[] taskTab = new Task[cores];
            float step = rows / ((float)cores);
            for (int i = 0; i < taskTab.Length; i++)
            {
                int s = (int)(i * step);
                int e = (int)((i + 1) * step);

                taskTab[i] = Task.Run(() => NormalizeRowsLog(s, e, res));
            }
            Task.WaitAll(taskTab);
            return res;
        }


        public SpareArray NormalizeRows(HashSet<int> index=null)
        {
            SpareArray res = null;
            res = new SpareArray(this);


            int cores = 5;
            Task[] taskTab = new Task[cores];
            float step = rows / ((float)cores);
            for (int i = 0; i < taskTab.Length; i++)
            {
                int s = (int)(i * step);
                int e = (int)((i + 1) * step);

                taskTab[i] = Task.Run(() => NormalizeRows(s, e, res,index));
            }
            Task.WaitAll(taskTab);
            return res;
        }
        public float this[int row, int col]
        {
            get
            {
                if (spareFlag)
                {
                    int cu = col;
                    if (dicData[row].ContainsKey(cu))
                        return dicData[row][cu];
                }
                else
                    return data[row, col];
                return 0;
            }
            set
            {

                if (spareFlag)
                {

                    if (value == 0)
                    {
                        lock (lockObj)
                        {
                            if (dicData[row].ContainsKey(col))
                                dicData[row].Remove(col);
                        }
                        return;
                    }
                    lock (lockObj)
                    {
                        /*if (!valuePos.ContainsKey(value))
                        {
                                values.Add(value);
                                valuePos.Add(value, (values.Count - 1));
                        }*/
                        dicData[row][col] = value;
                    }

                }
                else
                    data[row, col] = value;
            }
        }

        public float[,] GetArray()
        {
            if (spareFlag)
            {
                float[,] dat = new float[rows, columns];

                for (int i = 0; i < dat.GetLength(0); i++)
                    for (int j = 0; j < dat.GetLength(1); j++)
                        dat[i, j] = 0;

                for (int i = 0; i < dicData.Length; i++)
                    foreach (var it in dicData[i])
                        dat[i, it.Key] = it.Value;

                return dat;
            }
            return data;
        }
        public object Clone()
        {
            SpareArray clone = new SpareArray(rows, columns, spareFlag);
            if (data != null)
            {
                clone.data = new float[data.GetLength(0), data.GetLength(1)];
                for (int i = 0; i < data.GetLength(0); i++)
                    for (int j = 0; j < data.GetLength(1); j++)
                        clone.data[i, j] = data[i, j];

                //clone.data = (float[,])data.Clone();                

            }
            if (spareFlag)
            {
                //clone.valuePos = new Dictionary<float, int>();
                //foreach (var item in valuePos)
                //clone.valuePos.Add(item.Key, item.Value);

                //clone.values = new List<float>(values);

                for (int i = 0; i < dicData.Length; i++)
                    foreach (var item in dicData[i])
                    {
                        clone.dicData[i].Add(item.Key, item.Value);
                    }
                //clone.dicData =(Dictionary<int,float>[]) dicData.Clone();
            }


            return clone;
        }
        public void Add(SpareArray sp, int row)
        {
            if (spareFlag)
            {
                /*HashSet<double> vExist = new HashSet<double>();
                foreach (var item in values)
                    vExist.Add(item);

                foreach (var item in sp.values)
                    if (!vExist.Contains(item))
                    {
                        vExist.Add(item);
                        values.Add(item);
                        valuePos.Add(item, values.Count - 1);
                    }*/
                for (int i = 0; i < sp.dicData.Length; i++, row++)
                    foreach (var item in sp.dicData[i])
                    {
                        dicData[row].Add(item.Key, item.Value);
                    }
            }
            else
            {
                for (int i = row; i < sp.rows; i++)
                {
                    for (int j = 0; j < sp.columns; j++)
                        data[i, j] = sp.data[i - row, j];
                }
            }
        }
        public void Print()
        {
            if (spareFlag)
            {
                for (int i = 0; i < 3; i++)
                {
                    System.Diagnostics.Debug.Write(i + ": ");
                    foreach (var item in dicData[i])
                        System.Diagnostics.Debug.Write(" " + item.Key + ":" + item.Value + " ");
                    System.Diagnostics.Debug.WriteLine("");
                }
            }
            else

                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < columns; j++)
                    {
                        System.Diagnostics.Debug.Write(" " + this[i, j]);
                    }
                    System.Diagnostics.Debug.WriteLine("");
                }
        }
        public void Save(string fileName)
        {
            StreamWriter wr = new StreamWriter(fileName);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                    wr.Write(" " + this[i, j]);
                wr.WriteLine();
            }
            wr.Close();
        }
    }

    public class OmicsDataSet
    {

        //        HashSet<string> testData = new HashSet<string>() { "CAGATCAAGGTGCACA-1_10", "CTCGGAGAGTCATCCA-1_10", "TACTTACGTCTCTTAT-1_10", "GGAGCAATCGCGTAGC-1_30", "CCAATCCCAGGGTTAG-1_0", "TGCCCATTCCTACAGA-1_2", "CATCAGAGTGAGGCTA-1_2", "GCTGCAGAGGACAGCT-1_10", "TCACAAGTCCCAACGG-1_10", "GCCAAATGTATGAATG-1_2", "CCGGTAGAGGGTVGAT-1_6", "CTACACCTCGTCTGAA-1_10", "TCCACACTCATGTGGT-1_10", "TACCTTAAGGTGCTTT-1_10" };
        public string Name { get; set; }
        public SpareArray data = null;
        public List<string> geneLabels = null;
        public List<string> sampleLabels = null;
        public List<Interval[]> intervals = new List<Interval[]>();
        public HashSet<int> codes = new HashSet<int>();

        public List<int> selectedGenes = null;

        public List<string> metaData;

        public BindingList<FilterOmics> filters = new BindingList<FilterOmics>();

        public OmicsDataSet prev = null;


        public OmicsDataSet()
        {
        }
        public OmicsDataSet(OmicsDataSet data)
        {
            this.data = (SpareArray)data.data.Clone();
            if (data.geneLabels != null)
                this.geneLabels = new List<string>(data.geneLabels);
            if (data.sampleLabels != null)
                this.sampleLabels = new List<string>(data.sampleLabels);
            if (data.intervals != null)
            {

                this.intervals = new List<Interval[]>();
                foreach (var item in data.intervals)
                {
                    intervals.Add((Interval[])item.Clone());
                }
            }
            if (data.metaData != null)
                this.metaData = new List<string>(data.metaData);
            if (data.filters != null)
                this.filters = new BindingList<FilterOmics>(data.filters);
            if (data.codes != null)
                this.codes = new HashSet<int>(data.codes);
            this.Name = data.Name;

        }
        public OmicsDataSet(string fileName)
        {
            Name = fileName;
        }
        public void ClearFilters()
        {
            foreach (var item in filters)
                item.Clear();
        }
        public override string ToString()
        {
            return Name;
        }
        public int []GetSampleIndexes(List<string> samples)
        {
            HashSet<string> s = new HashSet<string>(samples);
            int[] res = new int[samples.Count];
            for(int i=0,n=0;i<sampleLabels.Count;i++)
            {
                if (s.Contains(sampleLabels[i]))
                    res[n++] = i;
            }
            return res;
        }
        public void AddCodes(Interval[] intV)
        {
            foreach (var item in intV)
                codes.Add(item.code);
        }
        public OmicsDataSet Transpose()
        {
            SpareArray transData = data.Transpose();

            OmicsDataSet res = new OmicsDataSet("Transpose");
            res.data = transData;
            res.geneLabels = geneLabels;
            res.sampleLabels = sampleLabels;

            return res;
        }

        public void Save(string fileName)
        {
            StreamWriter st = new StreamWriter(fileName);

            for (int i = 0; i < data.rows; i++)
            {
                if(data.rows==sampleLabels.Count)
                    st.WriteLine(">" + sampleLabels[i]);
                else
                    st.WriteLine(">" + geneLabels[i]);
                for (int j = 0; j < data.columns; j++)
                    st.Write(" " + data[i, j]);
                st.WriteLine();
            }

            st.Close();
        }
        public void Load(string fileName, char[] separators, string geneLocation, int genePos, int samplePos, int row, int column)
        {
            Name = Path.GetFileName(fileName);
            StreamReader sr = new StreamReader(fileName);
            Dictionary<string, List<int>> locData = new Dictionary<string, List<int>>();
            List<string> labRow = new List<string>();
            List<string> labCol = new List<string>();
            string line = sr.ReadLine();
            string[] aux = null;
            int count = 0;
            int startCol = samplePos;
            int labelRow = genePos;
            geneLabels = labRow;
            sampleLabels = labCol;

            if (geneLocation == "Column")
            {
                startCol = genePos;
                labelRow = samplePos;
                geneLabels = labCol;
                sampleLabels = labRow;

            }

            while (line != null)
            {
                line = line.Replace("\"", "");
                if (count == (int)labelRow)
                {
                    aux = line.Split(separators);
                    for (int i = column; i < aux.Length; i++)
                        labRow.Add(aux[i]);
                }
                if (count >= row)
                {
                    aux = line.Split(separators);
                    List<int> d = new List<int>();
                    for (int i = column; i < aux.Length; i++)
                    {
                        d.Add(Convert.ToInt32(aux[i]));
                    }
                    locData.Add(aux[startCol], d);
                    labCol.Add(aux[startCol]);
                }

                count++;
                line = sr.ReadLine();
            }
            sr.Close();


            List<string> locKeys = new List<string>(locData.Keys);

            if (locData.ContainsKey(geneLabels[0]))
            {
                data = new SpareArray(locData[locKeys[0]].Count, locData.Count, true);

                for (int i = 0; i < locData[locKeys[0]].Count; i++)
                    for (int j = 0; j < locKeys.Count; j++)
                    {
                        data[i, j] = locData[locKeys[j]][i];
                    }

            }
            else
            {
                data = new SpareArray(locData.Count, locData[locKeys[0]].Count, true);

                for (int j = 0; j < locKeys.Count; j++)
                {
                    for (int i = 0; i < locData[locKeys[j]].Count; i++)
                        data[j, i] = locData[locKeys[j]][i];
                }

            }


        }
        public static OmicsDataSet JoinOmicsData(List<OmicsDataSet> omicsData)
        {
            OmicsDataSet res = new OmicsDataSet("joined data");

            int count = 0;
            OmicsDataSet aux = null;
            for (int i = 0; i < omicsData.Count; i++)
            {

                count += omicsData[i].sampleLabels.Count;
            }
            res.data = new SpareArray(count, omicsData[0].geneLabels.Count, omicsData[0].data.spareFlag);
            res.geneLabels = new List<string>(omicsData[0].geneLabels);
            res.sampleLabels = new List<string>(count);
            int n = 0;


            int row = 0;
            for (int i = 0; i < omicsData.Count; i++)
            {
                aux = omicsData[i];
                res.sampleLabels.AddRange(aux.sampleLabels);
                res.data.Add(aux.data, row);
                row += aux.data.rows;
            }
            if (omicsData[0].codes != null)
            {
                res.codes = omicsData[0].codes;
                res.intervals = omicsData[0].intervals;
            }

            return res;
        }
        public OmicsDataSet RemoveGenes(List<string> genList)
        {
            OmicsDataSet removedData = new OmicsDataSet("Removed Genes");
            HashSet<int> ind = new HashSet<int>();
            for (int j = 0; j < genList.Count; j++)
            {
                for (int i = 0; i < geneLabels.Count; i++)
                {
                    if (geneLabels[i] == genList[j])
                    {
                        ind.Add(i);
                        break;
                    }
                }
            }
            List<string> newGenes = new List<string>();

            for (int i = 0; i < geneLabels.Count; i++)
                if (!ind.Contains(i))
                    newGenes.Add(geneLabels[i]);

            removedData.data = new SpareArray(data.rows, newGenes.Count, this.data.spareFlag);

            for (int i = 0; i < data.rows; i++)
                for (int j = 0, k = 0; j < data.columns; j++)
                    if (!ind.Contains(j))
                        removedData.data[i, k++] = data[i, j];

            removedData.sampleLabels = sampleLabels;
            removedData.filters = filters;

            return removedData;

        }
        public OmicsDataSet CreateSuperGenesData(Dictionary<string, HashSet<string>> superGenes)
        {
            OmicsDataSet superData = new OmicsDataSet(Name + " superGenes");

            superData.data = new SpareArray(data.rows, superGenes.Count, false);

            List<string> superGenesKeys = new List<string>(superGenes.Keys);
            superData.geneLabels = superGenesKeys;
            superData.sampleLabels = sampleLabels;
            superData.filters = filters;

            Dictionary<string, int> genePos = new Dictionary<string, int>();

            for (int i = 0; i < geneLabels.Count; i++)
                genePos.Add(geneLabels[i], i);

            for (int i = 0; i < data.rows; i++)
            {
                for (int n = 0; n < superGenesKeys.Count; n++)
                {
                    List<int> position = new List<int>();
                    foreach (var item in superGenes[superGenesKeys[n]])
                        if (genePos.ContainsKey(item))
                            position.Add(genePos[item]);

                    if (position.Count > 0)
                    {
                        float sum = 0;
                        for (int j = 0; j < position.Count; j++)
                            sum += data[i, position[j]];

                        superData.data[i, n] = sum / position.Count;
                    }
                    else
                        superData.data[i, n] = 0;
                }
            }

            return superData;
        }

        public Dictionary<double, int> ColFreq()
        {
            Dictionary<double, int> freq = new Dictionary<double, int>();

            for (int j = 0; j < data.columns; j++)
            {
                double sum = 0;
                for (int i = 0; i < data.rows; i++)
                    sum += data[i, j];

                if (!freq.ContainsKey(sum))
                    freq.Add(sum, 0);

                freq[sum]++;
            }
            return freq;
        }
        public static void CheckRanking(OmicsDataSet data, int order, string fileName)
        {
            Dictionary<string, int>[] freq = new Dictionary<string, int>[order];
            float[] d = new float[data.data.columns];
            int[] index = new int[d.Length];
            for (int i = 0; i < freq.Length; i++)
                freq[i] = new Dictionary<string, int>();

            for (int i = 0; i < data.data.rows; i++)
            {
                for (int j = 0; j < data.data.columns; j++)
                    index[j] = j;

                for (int j = 0; j < data.data.columns; j++)
                {
                    d[j] = data.data[i, j];
                }
                Array.Sort(index, (x, y) => d[y].CompareTo(d[x]));
                for (int j = 0; j < order; j++)
                {
                    if (!freq[j].ContainsKey(data.geneLabels[index[j]]))
                        freq[j].Add(data.geneLabels[index[j]], 0);
                    freq[j][data.geneLabels[index[j]]]++;
                }
            }
            StreamWriter st = new StreamWriter(fileName);
            for (int i = 0; i < freq.Length; i++)
            {
                var l = freq[i].ToList();
                l.Sort((p1, p2) => p2.Value.CompareTo(p1.Value));
                st.WriteLine("Position " + i + ":");
                foreach (var item in l)
                {
                    if (item.Value > 100)
                        st.WriteLine("\t" + item.Key + " " + item.Value);
                }
            }
            st.Close();
        }

        public OmicsDataSet ApplyFilters(BindingList<FilterOmics> filters)
        {
            OmicsDataSet current = this;

            foreach (var item in filters)
            {
                item.Clear();
                current = item.ApplyFilter(current);
                if (item.Name.Contains("Super"))
                {
                    //CheckRanking(current, current.geneLabels.Count, "superOrder");
                }
            }

            return current;
        }
        public void MakeStatistics()
        {
            double[] colStat = new double[data.rows];
            for (int j = 0; j < data.columns; j++)
            {
                double sum = 0;
                for (int i = 0; i < data.rows; i++)
                    sum += data[i, j];
                colStat[j] = sum;
            }
            double[] rowStat = new double[data.rows];
            for (int j = 0; j < data.rows; j++)
            {
                double sum = 0;
                for (int i = 0; i < data.columns; i++)
                    sum += data[j, i];
                rowStat[j] = sum;
            }


        }
        public OmicsDataSet SelectColumns(OmicsDataSet inputData, int[] index)
        {
            OmicsDataSet dataSet = new OmicsDataSet("selected columns");
            dataSet.data = new SpareArray(inputData.data.rows, index.Length, this.data.spareFlag);
            for (int i = 0; i < inputData.data.rows; i++)
            {
                for (int j = 0; j < index.Length; j++)
                    dataSet.data[i, j] = inputData.data[i, index[j]];

            }
            List<string> lab;
            if (geneLabels.Count != dataSet.data.rows)
                lab = geneLabels;
            else
                lab = sampleLabels;

            List<string> labels = new List<string>();
            for (int i = 0; i < index.Length; i++)
                labels.Add(lab[index[i]]);

            if (geneLabels.Count != dataSet.data.rows)
            {
                dataSet.geneLabels = labels;
                dataSet.sampleLabels = sampleLabels;
            }
            else
            {
                dataSet.geneLabels = geneLabels;
                dataSet.sampleLabels = labels;
            }
            dataSet.codes = codes;
            dataSet.intervals = intervals;


            return dataSet;

        }

    }
}
