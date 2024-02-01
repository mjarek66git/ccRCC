using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
//using System.Threading.Tasks;
using System.IO;
using System.Reflection;
using System.Globalization;
using System.Text.RegularExpressions;
//using DocumentFormat.OpenXml;
//using DocumentFormat.OpenXml.Packaging;
//using DocumentFormat.OpenXml.Spreadsheet;
using phiClustCore.Interface;
using System.Runtime.Serialization.Formatters.Binary;

//using phiClustCore;

namespace phiClustCore.Profiles
{
    [Serializable]
    public class OmicsProfile : UserDefinedProfile
    {
        public static string omicsSettings = "omicsSettings.dat";
        public OmicsInput oInput=new OmicsInput();

        public List<int> labelGeneStart = null;
        public List<int> labelSampleStart = null;       
//        profileNode localNode = null;
//        profileNode localNodeDist = null;
        List<string> labels = new List<string>();
        public List<string> []labelGenes =null;
        public List<string>[] labelSamples = null;
        public string[] labId = null;
        List<string> labelsData=new List<string>();
        static int profSize = 0;

        double []dev = null;
      
        double []avr=null;

        private static string ProfileName ="Omics profile";
        Settings set=new Settings();
        public OmicsProfile()
        {
            InitProfile();            
        }
        public OmicsProfile(Options opt):this()
        {
            InitProfile();
            //opt.SaveOptions("ggg");
           
            oInput = opt.omics;
            labelSampleStart = GetLabelsPositions(oInput.labelSampleStartString);
            labelGeneStart = GetLabelsPositions(oInput.labelGeneStartString);

            labelGenes = new List<string>[labelGeneStart.Count];
            for (int i = 0; i < labelGenes.Length; i++)
                labelGenes[i] = new List<string>();

            labelSamples = new List<string>[labelSampleStart.Count];
            for (int i = 0; i < labelSamples.Length; i++)
                labelSamples[i] = new List<string>();

            labId = new string[labelSamples.Length];
        }
        void InitProfile()
        {
            ProfileName = "Omics profile";
            profileName = ProfileName;
            AddInternalProfiles();
            destination = new List<INPUTMODE>();
            destination.Add(INPUTMODE.OMICS);
            maxV = 100;
            currentProgress = 0;
            set.Load();

        }

        List<int> GetLabelsPositions(string positions)
        {
            List<int> res = null;

            if (positions.Length > 0)
            {
                res = new List<int>();
                if (positions.Contains(";"))
                {
                    string[] aux = positions.Split(';');
                    foreach (var item in aux)
                    {
                        try
                        {
                            int x = Convert.ToInt32(item);
                            res.Add(x);
                        }
                        catch { }

                    }
                }
                else
                    try
                    {
                        int x = Convert.ToInt32(positions);
                        res.Add(x);
                    }
                    catch { }



            }
            return res;
        }
            public void Load(string fileName)
        {
            StreamReader file = new StreamReader(fileName);

            oInput.ReadOptionFile(file);
            labelSampleStart = GetLabelsPositions(oInput.labelSampleStartString);
            labelGeneStart = GetLabelsPositions(oInput.labelGeneStartString);

            labelGenes = new List<string>[labelGeneStart.Count];
            for (int i = 0; i < labelGenes.Length; i++)
                labelGenes[i] = new List<string>();

            labelSamples = new List<string>[labelSampleStart.Count];
            for (int i = 0; i < labelSamples.Length; i++)
                labelSamples[i] = new List<string>();

            labId = new string[labelSamples.Length];

            file.Close();
        }
        public void Save(string fileName)
        {
            StreamWriter file = new StreamWriter(fileName);
            if (labelSampleStart != null && labelSampleStart.Count > 0)
            {
                oInput.labelSampleStartString = "";
                int i = 0;
                for (i = 0; i < labelSampleStart.Count - 1; i++)
                    oInput.labelSampleStartString += labelSampleStart[i] + ";";
                oInput.labelSampleStartString += labelSampleStart[i];
            }
            if (labelGeneStart != null && labelGeneStart.Count > 0)
            {
                oInput.labelGeneStartString = "";
                int i = 0;
                for (i = 0; i < labelGeneStart.Count - 1; i++)
                    oInput.labelGeneStartString += labelGeneStart[i] + ";";
                oInput.labelGeneStartString += labelGeneStart[i];
            }

            oInput.SaveOptions(file);

            file.Close();
        }

        public override void AddInternalProfiles()
        {
            profileNode node = new profileNode();

            node.profName = ProfileName;
            node.internalName = ProfileName;
            InternalProfilesManager.AddNodeToList(node, this.GetType().FullName);
        }
        public void CombineTrainigTest(string fileName, string fileNameTrain, string fileNameTest)
        {
            System.Globalization.CultureInfo customCulture = (System.Globalization.CultureInfo)System.Threading.Thread.CurrentThread.CurrentCulture.Clone();
            customCulture.NumberFormat.NumberDecimalSeparator = ".";
            System.Threading.Thread.CurrentThread.CurrentCulture = customCulture;

            customCulture.NumberFormat.NumberDecimalSeparator = ".";
            StreamWriter wF = new StreamWriter(fileName);
            List<List<double>> dat = ReadOmicsFile(fileNameTest);
            List<string> auxRowTest = null;
            List<string> auxColTest = null;
            int i = 0;
            if (oInput.genePosition)
            {
                auxRowTest =new List<string>( labelGenes[0]);
                auxColTest = new List <string>(labelSamples[0]);
            }
            else
            {
                auxRowTest = new List<string>(labelSamples[0]);
                auxColTest = new List<string>(labelGenes[0]);
            }

            List<List<double>> trainDat = ReadOmicsFile(fileNameTrain);

            List<string> auxRow = null;
            List<string> auxCol = null;
            if (oInput.genePosition)
            {
                auxRow = labelGenes[0];
                auxCol = labelSamples[0];
            }
            else
            {
                auxRow = labelSamples[0];
                auxCol = labelGenes[0];
            }

            List<List<double>>  localDat=JoinData(trainDat, dat);

            wF.WriteLine();
            wF.Write("EMPTY ");
            for (i = 0; i < auxCol.Count; i++)
                wF.Write(auxCol[i] + " ");
            for (i = 0; i < auxColTest.Count - 1; i++)
                wF.Write(auxColTest[i] + " ");
            wF.WriteLine(auxColTest[auxColTest.Count - 1]);



            if (trainDat.Count == localDat.Count)
                for (i = 0; i < trainDat.Count; i++)
                {
                    wF.Write(auxRow[i] + " ");
                    for (int j = 0; j < localDat[i].Count; j++)
                        wF.Write(localDat[i][j] + " ");
                    wF.WriteLine();
                }
            else
            {
                for (i = 0; i < trainDat.Count; i++)
                {
                    wF.Write(auxRow[i] + " ");
                    for (int j = 0; j < localDat[i].Count; j++)
                        wF.Write(localDat[i][j] + " ");
                    wF.WriteLine();
                }
                for (int k = 0; k < dat.Count; k++)
                {
                    wF.Write(auxRowTest[k] + " ");
                    for (int j = 0; j < dat[k].Count; j++)
                        wF.Write(dat[k][j] + " ");
                    wF.WriteLine();
                }
            }

            wF.Close();
        }
        public static List<KeyValuePair< string,List<byte>>> ReadOmicsProfile(string fileName)
        {
            List<KeyValuePair<string, List<byte>>> dic = new List<KeyValuePair< string,List<byte>>>();
            Settings set = new Settings();
            set.Load();
            string key="";
            string value;

            StreamReader r = new StreamReader(set.profilesDir+Path.DirectorySeparatorChar+fileName);

            string line = r.ReadLine();
            while(line!=null)
            {
                if(line.Contains(">"))                
                    key = line.Remove(0, 1);
                else
                {
                    if(line.Contains(ProfileName))
                    {
                        value = line.Remove(0, ProfileName.Length+1);
                        string[] aux = value.Split(' ');
                        List<byte> v = new List<byte>(aux.Length);
                        foreach (var item in aux)
                            v.Add(Convert.ToByte(item, CultureInfo.InvariantCulture));
                        KeyValuePair<string, List<byte>> xx = new KeyValuePair<string, List<byte>>(key, v);
                        dic.Add(xx);
                    }
                }
                line = r.ReadLine();
            }
            r.Close();
            return dic;
        }
        public override void RunThreads(string fileName)
        {
            ThreadFiles ff = new ThreadFiles();

            //StreamReader str = new StreamReader(fileName);
            ff.fileName = fileName ;     
            Run((object)ff);
        }
        string GetProfileName(string fileName)
        {
            return set.profilesDir + Path.DirectorySeparatorChar +  oInput.processName;
        }
        public double [,] QuantileNorm(double[,] data)
        {
            int[][] copyData;
            double[,] rankData;
            Dictionary<double, int> dic = new Dictionary<double, int>();
            double[] avr;

            avr = new double[data.GetLength(0)];
            rankData = new double[data.GetLength(0), data.GetLength(1)];
            copyData = new int[data.GetLength(1)][];
            for (int i = 0; i < data.GetLength(1); i++)
                copyData[i] = new int[data.GetLength(0)];

            for (int i = 0; i < data.GetLength(1); i++)
                for (int j = 0; j < data.GetLength(0); j++)
                    copyData[i][j] = j;

            for (int i = 0; i < data.GetLength(1); i++)
            {
                Array.Sort<int>(copyData[i], (a, b) => data[a, i].CompareTo(data[b, i]));
                dic.Clear();
                for (int j = 0; j < data.GetLength(0); j++)
                {
                    if (data[copyData[i][j], i] != double.NaN)
                    {
                        if (!dic.ContainsKey(data[copyData[i][j], i]))
                            dic.Add(data[copyData[i][j], i], j + 1);
                        rankData[copyData[i][j], i] = dic[data[copyData[i][j], i]];
                    }
                    else
                        rankData[copyData[i][j], i] = double.NaN;
                }
            }


            for (int i = 0; i < avr.Length; i++)
            {
                double sum = 0;
                for (int j = 0; j < data.GetLength(1); j++)
                    if (data[copyData[j][i], j]!=double.NaN)
                        sum += data[copyData[j][i], j];
                avr[i] = sum / data.GetLength(1);
            }

            for (int i = 0; i < data.GetLength(0); i++)
                for (int j = 0; j < data.GetLength(1); j++)
                    if (rankData[i, j]!=double.NaN)
                        rankData[i, j] = avr[(int)rankData[i, j] - 1];


            return rankData;
        }
        void AddLabelsToGeneORSample(List <string> labels,int k)
        {
            
            if (!oInput.genePosition)
                labelGenes[k]=labels;
            else
                labelSamples[k]=labels;

        }
        List<int>  GetLabelPosition()
        {
            if (!oInput.genePosition)
                return labelGeneStart;
            else
                return labelSampleStart;

        }
        int ReadLabels(StreamReader r,string line,char delimiter=' ')
        {
            int i=0;
            int countPos=0;
            List<int> rowPositions;
            List<string> tmpAux;
            char locDelimiter = ' ';
            rowPositions = GetLabelPosition();
            rowPositions.Sort(delegate (int x,int y) { return x.CompareTo(y); });
            for (i=1;i<=rowPositions[rowPositions.Count-1]; i++)
            {
                line = r.ReadLine();
                if (!line.Contains(delimiter))
                    delimiter = locDelimiter;
                line = ProcessCSVLine(line, delimiter);
                if (i==rowPositions[countPos])
                {
                    tmpAux = new List<string>();
                    string[] aux = line.Split(delimiter);
                    labId[countPos]=aux[0];
                    for (int j = oInput.numCol, n = 0; j < aux.Length; j++, n++)
                        tmpAux.Add(aux[j]);

                    AddLabelsToGeneORSample(tmpAux,countPos);
                    countPos++;
                }

            }
            List<string>[] tmp = null;
            if (!oInput.genePosition)
                tmp=labelGenes;
            else
                tmp=labelSamples;

            bool generic = false;
            Dictionary<string, int> xx = new Dictionary<string, int>();
            for(int j=0;j<tmp[0].Count;j++)
            {
                if (!xx.ContainsKey(tmp[0][j]))
                    xx.Add(tmp[0][j], 1);
                else
                {
                    generic = true;
                    break;
                }
            }
            if(generic)
            {
                List<string>[] yy = new List<string>[tmp.Length + 1];
                List<string> zz = new List<string>();
                for (int j = 0; j < tmp[0].Count; j++)
                    zz.Add("Num_" + j);

                yy[0] = zz;
                for (int j = 1; j < yy.Length; j++)
                    yy[j] = tmp[j - 1];

                if (!oInput.genePosition)
                    labelGenes=yy;
                else
                    labelSamples=yy;

            }

            return i;
        }
        public List<string> ReadClassLabels(string fileName,bool column,int pos)
        {
            List<string> labels = new List<string>();
            string line="";
            char tab = '\t';
            StreamReader r = new StreamReader(fileName);
            string[] aux;

            if (r == null)
                throw new Exception("Cannot open file: " + fileName);

            if(!column)
            {
                for (int i=0; i < pos; i++)
                    line = r.ReadLine();
                line = line.Replace(tab.ToString(), " ");
                aux = line.Split(' ');
                for (int i = oInput.numCol; i < aux.Length; i++)
                    labels.Add(aux[i]);

            }
            else
            {
                for (int i=0; i <=oInput.numRow; i++)
                    line = r.ReadLine();
              
                while (line != null)
                {
                    line = line.Replace(tab.ToString(), " ");
                    aux = line.Split(' ');
                    labels.Add(aux[pos - 1]);
                    line = r.ReadLine();
                }
            }
            r.Close();

            return labels;
        }
       /* public List<List<double>> ReadOmicsExcelFile(string fileName)
        {
            List<List<double>> localData = new List<List<double>>();

            for (int i = 0; i < labelGenes.Length;i++ )
                labelGenes[i].Clear();
            for (int i = 0; i < labelSamples.Length; i++)
                labelSamples[i].Clear();
            List<int> labelPosition=GetLabelPosition();
            StringBuilder text = new StringBuilder();
            List<string> labels = new List<string>();
            int remProgress = currentProgress;
            int numLabel = 0;
            using (SpreadsheetDocument spr = SpreadsheetDocument.Open(fileName, false))
            {
                WorkbookPart wp = spr.WorkbookPart;
                WorksheetPart wsp = wp.WorksheetParts.First();

                
              //  IEnumerable<SheetData> sheetData = wsp.Worksheet.Elements<SheetData>();
//                IEnumerable<Row> row = sheetData.First().Elements<Row>();
                
                int maxRows = 0;
                 OpenXmlReader reader = OpenXmlReader.Create(wsp);
                 while (reader.Read())
                 {
                     if (reader.ElementType == typeof(Row) && reader.IsStartElement)
                     {
                         do
                         {
                             maxRows++;
                         } while (reader.ReadNextSibling()); // Skip to the next row
                         break;
                     }
                 }
                 reader = OpenXmlReader.Create(wsp);
                 double step = 20.0 /maxRows;
                 //currentProgress += 5;
                Cell c=null;
                int rem=0;
                int rowCounter = 0;
                
                while (reader.Read())
                {
                    if (reader.ElementType == typeof(Row) && reader.IsStartElement)
                    {
                        rowCounter++;
                        if (rowCounter < oInput.numRow - 1 && rowCounter != labelPosition[0])                                                 
                            continue;
               
                        if((int)(rowCounter*step)>rem)
                        {
                            if(currentProgress<20)
                                currentProgress++;
                            rem = (int)(rowCounter * step);
                        }
         
                        text.Clear();
                        reader.ReadFirstChild();                        
                        if (reader.ElementType == typeof(Cell))
                        {
                            int cellCounter = 0;
                            do
                            {
                                if (text.Length > 0)
                                    text.Append(";");
                                c = (Cell)reader.LoadCurrentElement();
                                for (numLabel = 0; numLabel < labelPosition.Count;numLabel++)
                                    if (rowCounter == labelPosition[numLabel])
                                    {
                                        if (cellCounter >= oInput.numCol)
                                            labels.Add(GetCellValue(c, wp));
                                        cellCounter++;
                                    }
                                    else
                                        text.Append(GetCellValue(c, wp));
                            }
                            while (reader.ReadNextSibling());
                            if(rowCounter>=oInput.numRow)
                                ProcessRow(text.ToString(), localData,';');
                        }                        
                    }
                }                
                AddLabelsToGeneORSample(labels,numLabel);

            }
            currentProgress += 20-(currentProgress-remProgress);
            GenerateDefaultLabels(localData.Count, localData[0].Count);
            
            return localData;
        }*/

        /*static string GetCellValue(Cell c, WorkbookPart workbookPart)
        {
            string cellValue = string.Empty;
            if (c.DataType != null && c.DataType == CellValues.SharedString)
            {
                SharedStringItem ssi =
                    workbookPart.SharedStringTablePart.SharedStringTable
                        .Elements<SharedStringItem>()
                        .ElementAt(int.Parse(c.CellValue.InnerText));
                if (ssi.Text != null)
                {
                    cellValue = ssi.Text.Text;
                }
            }
            else
            {
                if (c.CellValue != null)
                {
                    cellValue = c.CellValue.InnerText;
                }
            }
            return cellValue;
        }
        

        static string GetCellValue(SpreadsheetDocument document, Cell cell)
        {
            SharedStringTablePart stringTablePart = document.WorkbookPart.SharedStringTablePart;
            string value = cell.CellValue.InnerXml;

            if (cell.DataType != null && cell.DataType.Value == CellValues.SharedString)
            {
                return stringTablePart.SharedStringTable.ChildElements[Int32.Parse(value)].InnerText;
            }
            else
            {
                return value;
            }
        }*/
        private void GenerateDefaultLabels(int dataCount,int rowCount)
        {
            int i;
            if (!oInput.genePosition)
            {
                if (labelSamples[0] != null && labelSamples[0].Count == 0)
                {
                    labelSamples[0] = new List<string>();
                    for (i = 0; i < dataCount; i++)
                        labelSamples[0].Add("Sample_" + (i + 1));
                }
                if (labelGenes[0] == null || labelGenes[0].Count == 0)
                {
                    labelGenes[0] = new List<string>();
                    for (i = 0; i < rowCount; i++)
                        labelGenes[0].Add("Gene_" + (i + 1));
                }
            }
            else
            {
                if (labelGenes[0] == null || labelGenes[0].Count == 0)
                {
                    labelGenes[0] = new List<string>();
                    for (i = 0; i < dataCount; i++)
                        labelGenes[0].Add("Gene_" + (i + 1));
                }
                if (labelSamples[0] != null && labelSamples[0].Count == 0)
                {
                    labelSamples[0] = new List<string>();
                    for (i = 0; i < rowCount; i++)
                        labelSamples[0].Add("Sample_" + (i + 1));

                }
            }

        }
        private void ProcessRow(string line,List<List<double>> data,char delimiter=' ')
        {
            List<double> row = new List<double>();
            if(delimiter==' ')
                line = Regex.Replace(line, @"\s+", " ");
            line = line.TrimEnd();
            string[] aux = line.Split(delimiter);
            if (aux.Length == 0)
                return;
            try
            {
                for (int i = oInput.numCol; i < aux.Length; i++)
                {
                    double tmp = 0;
                    tmp = Convert.ToDouble(aux[i], CultureInfo.InvariantCulture);
                    row.Add(tmp);
                }
                if (oInput.genePosition && oInput.uLabelGene)
                {
                    int count=0;

                    foreach (var item in labelGeneStart)
                        labelGenes[count++].Add(aux[item - 1]);
                }
                else
                    if (!oInput.genePosition && oInput.uLabelSample)
                    {
                        int count = 0;
                        foreach (var item in labelSampleStart)
                            labelSamples[count++].Add(aux[item - 1]);

                    }
                if (row.Count > 0)
                    data.Add(row);

            }
            catch (FormatException)
            {

            }

        }
        public List<List<double>> JoinData(List<List<double>> data1, List<List<double>> data2)
        {
            List<List<double>> localData = new List<List<double>>();


            if (data1.Count == data2.Count)
            {


                for (int i = 0; i < data1.Count; i++)
                {
                    List<double> rowData = new List<double>();
                    foreach (var item in data1[i])
                        rowData.Add(item);
                    foreach (var item in data2[i])
                        rowData.Add(item);

                    localData.Add(rowData);
                }
            }
            else
            {
                if (data1[0].Count != data2[0].Count)
                    throw new Exception("nes! Cannot be combined");

                for (int i = 0; i < data1.Count; i++)
                {
                    List<double> rowData = new List<double>();
                    foreach (var item in data1[i])
                        rowData.Add(item);
                    localData.Add(rowData);
                }

                for (int i = 0; i < data2.Count; i++)
                {
                    List<double> rowData = new List<double>();
                    foreach (var item in data2[i])
                        rowData.Add(item);

                    localData.Add(rowData);
                }

            }


            return localData;

        }
        string ProcessCSVLine(string line,char delimiter)
        {
            StringBuilder newLine = new StringBuilder();
            Boolean q = false;
            for(int i=0;i<line.Length;i++)
            {
                if (line[i] == '"')
                {
                    q = !q;
                    continue;

                }
                if (q && line[i] == delimiter)
                    newLine.Append('-');
                else
                    newLine.Append(line[i]);                

            }

            return newLine.ToString();
        }
        public List<List<double>> ReadOmicsFile(string fileName)
        {
            List<List<double>> localData = new List<List<double>>();

            if (!File.Exists(fileName))
                throw new Exception("File :" + fileName + " cannot be open");


            //if (Path.GetExtension(fileName).Contains("xlsx"))
              //  return ReadOmicsExcelFile(fileName);

            StreamReader r = new StreamReader(fileName);
            char delimiter = '\t';
            if (Path.GetExtension(fileName).Contains("csv"))
                delimiter = ',';

            if (r == null)
                throw new Exception("Cannot open file: " + fileName);
            int i=0;
            
            string line= r.ReadLine();
            List<List<double>> data = new List<List<double>>();
            if(labelGenes!=null)
                for (int  v= 0; v < labelGenes.Length;v++ )
                    labelGenes[v].Clear();
            if(labelSamples!=null)
                for (int v = 0; v < labelSamples.Length; v++)
                    labelSamples[v].Clear();

            if(oInput.genePosition && oInput.uLabelSample || !oInput.genePosition && oInput.uLabelGene)
                i = ReadLabels(r,line,delimiter);

            for (; i < oInput.numRow-1; i++)
                line = r.ReadLine();


           // if(line.Contains())
            currentProgress += 5;
            string remLine = line;
            while (line != null)
            {
                if (line.Length > 5)
                {
                    if (!line.Contains(delimiter))
                        delimiter = ' ';
                    if (line.Contains("\""))
                        line = ProcessCSVLine(line, delimiter);
                        //line = Regex.Replace(line, "\"", "");
                    ProcessRow(line, localData, delimiter);
                }
                remLine = line;
                line = r.ReadLine();
            }
            r.Close();
            currentProgress += 15;
            remLine += "ll";
            GenerateDefaultLabels(localData.Count, localData[0].Count);
            
            return localData;
        }
        double [,] ApplySuperGenes(double [,] data,int geneIndex,string fileName)
        {            
            Dictionary<string, List<int>> superGenes = ReadSuperGenes(fileName,geneIndex);
            bool[] indexNonSuper = new bool[data.GetLength(0)];

            for (int i = 0; i < indexNonSuper.Length; i++)
                indexNonSuper[i] = true;

            int countAll = 0;
            foreach (var item in superGenes)
            {
                countAll += item.Value.Count;
                for (int i = 0; i < item.Value.Count; i++)
                    indexNonSuper[item.Value[i]] = false;
            }
            int count = 0;
            for (int i = 0; i < indexNonSuper.Length; i++)
                if (indexNonSuper[i])
                    count++;

            int finalGeneNum = superGenes.Count;// count + superGenes.Count;

            double[,] newData=new double[finalGeneNum,data.GetLength(1)];

            for (int i=0;i<newData.GetLength(1);i++)
            {
                int superIndex = 0;
                foreach(var item in superGenes)
                {
                    double sum = 0;
                    for (int j = 0; j < item.Value.Count; j++)
                        sum += data[item.Value[j], i];
                    newData[superIndex++, i] = sum;// (int)(sum/item.Value.Count);
                }
                /*for (int n = 0; n < indexNonSuper.Length; n++)
                    if (indexNonSuper[n])
                            newData[superIndex++, i] = data[n, i];*/
                
            }

            //Zscore 
            var zscore=ZscoreCodingRow(newData);

            double[,] dataZ = zscore.Item1;
            double[] avr = zscore.Item2;
            double[] std = zscore.Item3;
            //Coding
            for(int i=0;i<dataZ.GetLength(1);i++)
            {
                for(int j=0;j<dataZ.GetLength(0);j++)
                {
                    if (dataZ[j, i] >=- 1 && dataZ[j, i] < 1)
                        dataZ[j, i] = 2;
                    else
                    if (dataZ[j, i] <-1)
                        dataZ[j, i] = 1;
                    else
                        dataZ[j, i] = 3;

                }
            }
            newData = dataZ;
            List<string> newLabels = new List<string>();

            foreach(var item in superGenes)
                newLabels.Add(item.Key);

            /*for(int i=0;i<labelGenes[geneIndex].Count;i++)
                if (indexNonSuper[i])
                    newLabels.Add(labelGenes[geneIndex][i]);*/

            labelGenes[geneIndex] = newLabels;

            return newData;
        }

        Dictionary<string,List<int>> ReadSuperGenes(string fileName,int geneIndex)
        {
            Dictionary<string, List<int>> superGenes = new Dictionary<string, List<int>>();
            StreamReader sr = new StreamReader(fileName);

            Dictionary<string, int> hashGeneNames = new Dictionary<string, int>();

            for(int i=0;i<labelGenes[geneIndex].Count;i++)
                hashGeneNames.Add(labelGenes[geneIndex][i], i);

            Dictionary<int, int> check = new Dictionary<int, int>();

            string superGeneName = "";
            string line = sr.ReadLine();
            while(line!=null)
            {
                if(line.StartsWith(">"))                
                    superGeneName=line.Substring(1);
                else
                {
                    string[] aux = line.Split(' ');
                    List<string> l = aux.ToList();
                    List<int> genesPositions = new List<int>();
                    foreach (var item in l)
                        if (hashGeneNames.ContainsKey(item))
                            genesPositions.Add(hashGeneNames[item]);
                        else
                            throw new Exception("Uknown gene: " + item);
                    superGenes.Add(superGeneName, genesPositions);
                }

                line = sr.ReadLine();
            }
            sr.Close();

            return superGenes;
        }
        static double [,] TransposeData(double [,] data)
        {
            double [,]dataFinal = new double[data.GetLength(1), data.GetLength(0)];

            for (int i = 0; i < data.GetLength(0); i++)
                for (int j = 0; j < data.GetLength(1); j++)
                    dataFinal[j, i] = data[i, j];


                    return dataFinal;
        }
        void Check(double [,] data)
        {
            /*Console.WriteLine("Wiersze");
            for(int i=0;i<data.GetLength(0);i++)
            {
                double sum=0;                
                for (int j = 0; j < data.GetLength(1); j++)
                    sum += data[i, j];

                Console.WriteLine("sum="+sum);
            }*/
            Console.WriteLine("Columny");
            List<double> median = new List<double>();
            for (int j = 0; j < data.GetLength(1); j++)
            {
                double sum = 0;
                median.Clear();
                for (int i = 0; i < data.GetLength(0); i++)
                    median.Add(data[i, j]);
                    //sum += data[i, j];
                median.Sort((x, y) => x.CompareTo(y));
                sum = median[median.Count / 2];
                Console.WriteLine("sum=" + sum);
            }
        }
        double [,] SelectNMostDiff(double [,]data,int n)
        {
            double [,] resData=new double[n,data.GetLength(1)];
            double sumX;
            double sumX2;

            dev = new double[data.GetLength(0)];

            for (int i = 0; i < data.GetLength(0); i++)
            {
                sumX = 0;
                sumX2 = 0;
                for (int j = 0; j < data.GetLength(1); j++)
                    if (data[i, j] != double.NaN)
                    {
                        sumX += data[i, j];
                        sumX2 += data[i, j] * data[i, j];
                    }

                sumX /= data.GetLength(1);
                sumX2 /= data.GetLength(1);
                dev[i] = Math.Sqrt(sumX2 - avr[i] * avr[i]);

            }
        


            return resData;
        }
        double[,] QuantileCoding(double[,] dataFinal)
        {
            if (!oInput.genePosition)
                dataFinal = TransposeData(dataFinal);
            
            if (oInput.zScore)
                dataFinal = StandardData(dataFinal, false, oInput.selectGenes);
            currentProgress += 20;

            dataFinal = QuantileNorm(dataFinal);
            double[,] outData = dataFinal;
            if (oInput.coding!=CodingAlg.NONE)
                outData = IntervalCoding(dataFinal);

            return outData;
        }
        double[,] ZScoreCoding(double[,] dataFinal)
        {
            if (!oInput.genePosition)
                dataFinal = TransposeData(dataFinal);
           
            dataFinal = StandardData(dataFinal, false, oInput.selectGenes);

            double[,] cc = dataFinal;// TransposeData(dataFinal);
            currentProgress += 20;

            double[,] outData= dataFinal;
           // double[,] outData = IntervalCoding(dataFinal);
           if (oInput.coding!=CodingAlg.NONE)
                outData = IntervalCodingPerSample(dataFinal);


            //outData = TransposeData(outData);

            return outData;
        }

        public override int Run(object processParams)
        {
            string fileName = ((ThreadFiles)(processParams)).fileName;
            StreamReader r = new StreamReader(fileName);
            int i;


            if (oInput.heatmap)
                return 0;
            //Load(OmicsInput.fileName);
            List<List<double>> data = ReadOmicsFile(fileName);

            dev = new double[data.Count];
            avr = new double[data.Count];


            double [,]dataFinal;
            dataFinal=new double [data.Count,data[0].Count];
            for (i = 0; i < data.Count; i++)
                for (int j = 0; j < data[i].Count; j++)
                        dataFinal[i, j] = data[i][j];
            data = null;
            if (oInput.genePosition)
            {
                List<string>[] temp;
                if (labelGenes[0].Count > 0)
                    temp = labelGenes;
                else
                    temp = labelSamples;

                labelsData = new List<string>();
                string s = "";
                for (int n = 0; n < temp[0].Count; n++)
                {
                    s = "";
                    for (int k = 0; k < temp.Length - 1; k++)
                        s += temp[k][n] + ";";
                    s += temp[temp.Length - 1][n];
                    labelsData.Add(s);
                }

            }
            double[,] outData;

            if (oInput.fileSelectedGenes.Length > 0)
            {
                dataFinal = ApplySuperGenes(dataFinal, 0, oInput.fileSelectedGenes);
                HashSet<int> codes = new HashSet<int>();
                for (int n = 0; n < dataFinal.GetLength(0); n++)
                    for (int j = 0; j < dataFinal.GetLength(1); j++)
                        if (!codes.Contains((int)dataFinal[n, j]))
                            codes.Add((int)dataFinal[n, j]);
                SaveLabIdAndIntervals(codes);
            }

            if (oInput.quantile)
                outData = QuantileCoding(dataFinal);
            else
                if (oInput.zScore)
                outData = ZScoreCoding(dataFinal);
            else
                 if (oInput.coding == CodingAlg.MAX_VALUE)
                outData = MaxValueThreshold(dataFinal, oInput.numStates);
            else
                if (oInput.coding != CodingAlg.NONE)
                    outData = IntervalCoding(dataFinal);
            else           
                outData = dataFinal;
                         


            StreamWriter wr;
            string profFile = GetProfileName(fileName);

            wr = new StreamWriter(profFile);
            
            int l=0;
            for (i = 0; i < outData.GetLength(0); i++)
            {
                wr.WriteLine(">" + labelsData[i]);
                wr.Write(ProfileName + " ");
                for (l = 0; l < outData.GetLength(1) - 1; l++)
                    wr.Write((int)outData[i, l] + " ");
                wr.WriteLine((int)outData[i, l]);
            }
            wr.Close();
            currentProgress += 20;
            List<string>[] tx;
            if (!oInput.genePosition)            
                tx = labelGenes;
            else
                tx=labelSamples;
                

                labelsData = new List<string>();
                string st = "";


                for (int n = 0; n < tx[0].Count; n++)
                {
                    st = "";
                    for (int k = 0; k < tx.Length - 1; k++)
                        st += tx[k][n] + ";";
                    st += tx[tx.Length - 1][n];
                    labelsData.Add(st);
                }
            
/*            if (!genePosition)
                    labelsData = labelGenes[0];
            else
                    labelsData = labelSamples[0];*/

            wr.Close();
            string profFileTransp = profFile + "_transpose";
            wr = new StreamWriter(profFileTransp);
            for (i = 0; i < outData.GetLength(1); i++)
            {
                wr.WriteLine(">" + labelsData[i]);
                wr.Write(ProfileName + " ");
                for (l = 0; l < outData.GetLength(0) - 1; l++)
                    wr.Write((int)outData[l,i] + " ");
                wr.WriteLine((int)outData[l,i]);
            }
            wr.Close();
            currentProgress += 20;
            ProfileTree td=ProfileAutomatic.AnalyseProfileFile(profFile, SIMDIST.DISTANCE, ProfileName);
            List<string> keys=new List<string>(td.masterNode.Keys);

            ProfileTree ts = ProfileAutomatic.AnalyseProfileFile(profFile, SIMDIST.SIMILARITY, ProfileName);
                string locPath = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location) + Path.DirectorySeparatorChar + "profiles" + Path.DirectorySeparatorChar;
                td.SaveProfiles(locPath + Path.GetFileNameWithoutExtension(fileName) + "_distance.profile");
                ts.SaveProfiles(locPath + Path.GetFileNameWithoutExtension(fileName) + ".profiles");
                profSize = td.masterNode[keys[0]].codingToByte.Count;

                List<string> m = new List<string>(ts.masterNode.Keys);
                //localNode = ts.masterNode[m[0]];

                currentProgress = maxV;            

            return 0;

        }
        void PrepareStandardParams(double [,]data)
        {
            double sumX = 0;
            double sumX2 = 0;

            dev = new double[data.GetLength(0)];
            avr = new double[data.GetLength(0)];

            for (int i = 0; i < data.GetLength(0); i++)
            {
                sumX = 0;
                sumX2 = 0;
                int counter = 0;
                for (int j = 0; j < data.GetLength(1); j++)
                    if (data[i, j] != double.NaN)
                    {
                        sumX += data[i, j];
                        sumX2 += data[i, j] * data[i, j];
                        counter++;
                    }

                if (counter > 0)
                {
                    sumX /= counter;
                    sumX2 /= counter;
                    avr[i] = sumX;
                    dev[i] = Math.Sqrt(sumX2 - avr[i] * avr[i]);
                }
            }

        }
        double [,] StandardData(double [,]data,bool externalDataFlag=false,int nData=0)
        {
            if(!externalDataFlag)
                PrepareStandardParams(data);

            if (nData > data.GetLength(0))
                nData = 0;

            for (int i = 0; i < data.GetLength(0); i++)
            {

                    for (int j = 0; j < data.GetLength(1); j++)
                        if (data[i, j] != double.NaN && data[i,j]!=0 && dev[i] > 0 )
                            data[i, j] = (data[i, j] - avr[i]) / dev[i];
            }

            if(nData==0)            
                return data;           
            else
            {
                int[] index = new int[dev.Length];
                for (int i = 0; i < dev.Length; i++)
                    index[i] = i;

                Array.Sort(dev, index);
                Array.Reverse(index);
                double[,] locData = new double[nData, data.GetLength(1)];
                for(int i=0;i<nData;i++)
                {
                    for(int j=0;j<data.GetLength(1);j++)                    
                        locData[i, j] = data[index[i], j];                    
                }
                return locData;
            }
            
        }

        static double[,] LoadIntervals(string fileName)
        {
            double[,] intervals=null;
            StreamReader r = new StreamReader(fileName);
            string line = r.ReadLine();
            int num=0;
            while(line!=null)
            {
                if (line.Contains("Code"))
                    num++;
                line = r.ReadLine();
            }

            intervals = new double[num, 2];
            r.BaseStream.Seek(0, System.IO.SeekOrigin.Begin);
            num=0;
            line = r.ReadLine();
            while (line != null)
            {
                if (line.Contains("Code"))
                {
                    string[] aux = line.Split(' ');
                    if (aux.Length == 4)
                    {
                        intervals[num, 0] = Convert.ToDouble(aux[2]);
                        intervals[num++, 0] = Convert.ToDouble(aux[3]);
                    }
                }
                line = r.ReadLine();
            }            
            r.Close();

            return intervals;
        }
        static int[] IntervalCodigPerSample(double[,] data,int num, double[,] intervals)
        {
            int[] codedData = new int[data.GetLength(1)];
            for (int j = 0; j < data.GetLength(1); j++)
            {
                int code = intervals.GetLength(0);

                if (data[num,j] < intervals[0, 0])
                {
                    codedData[j] = 1;
                    continue;
                }
                if (data[num,j] > intervals[intervals.GetLength(0) - 1, 1])
                {                   
                    codedData[j]= intervals.GetLength(0) - 1;
                    continue;
                }
                for (int k = 0; k < intervals.GetLength(0); k++)
                {
                    if (data[num,j] == double.NaN)
                    {
                        code = 0;
                        break;
                    }
                    if (data[num,j] >= intervals[k, 0] && data[num,j] < intervals[k, 1])
                    {
                        code = k + 1;
                        break;
                    }
                }
       
                codedData[j] = code;
            }

            return codedData;
        }
        static double [,] IntervalCodig(double[,] data,double [,] intervals)
        {
            double[,] newData = new double[data.GetLength(0), data.GetLength(1)];
            int[] codedRow = new int[newData.GetLength(1)];
            for (int i = 0; i < data.GetLength(0); i++)
            {

                codedRow=IntervalCodigPerSample(data,i, intervals);
                /*for (int j = 0; j < data.GetLength(1); j++)
                {
                    int code = intervals.GetLength(0);

                    if (data[i, j] < intervals[0, 0])
                    {
                        code = 1;
                        continue;
                    }
                    if(data[i,j]>intervals[intervals.GetLength(0)-1,1])
                    {
                        code = intervals.Length - 1;
                        continue;
                    }
                    for (int k = 0; k < intervals.GetLength(0); k++)
                    {                        
                        if (data[i, j] == double.NaN)
                        {
                            code = 0;
                            break;
                        }
                        if (data[i, j] >= intervals[k, 0] && data[i, j] < intervals[k, 1])
                        {
                            code = k + 1;
                            break;
                        }
                    }
                    codedRow[j] = code;
                }*/
                for (int n = 0; n < newData.GetLength(1); n++)
                    newData[i, n] = codedRow[n];

            }

            return newData;
        }
        
        double [,] IntervalCodingPerSample(double [,]data)
        {
            double[,] newData = data;
            Dictionary<double, int> hashValues = new Dictionary<double, int>();
            double[] colValues = new double[data.GetLength(1)];
            double[,] intervals=null;
            double[,] outData = new double[data.GetLength(0), data.GetLength(1)];
            if (oInput.coding == CodingAlg.Z_SCORE)
                newData = ZscoreCodingColumn(data).Item1;

            for (int j = 0; j < newData.GetLength(0); j++)
            {
                hashValues.Clear();
                for (int i = 0; i < newData.GetLength(1); i++)
                {
                    colValues[i] = data[j, i];

                    if (!hashValues.Keys.Contains(data[j, i]))
                        hashValues.Add(data[j, i], 1);
                    else
                        hashValues[data[j, i]]++;
                }
                int[] coded=null;
                //Array.Sort(colValues);
                try
                {

                    intervals = SetupIntervals(hashValues);
                }
                catch (Exception ex)
                {
                    Console.Write("UPS:"+ex.Message);

                }
                coded = IntervalCodigPerSample(data, j, intervals);

                for (int i = 0; i < coded.Length; i++)
                    outData[j, i] = coded[i];
            }

            SaveLabIdAndIntervals(intervals);

            return outData;

        }
        void SaveLabIdAndIntervals(object inter)
        {
            StreamWriter ir;
            ir = new StreamWriter("generatedProfiles/labels.dat");
            ir.WriteLine("[GENES]");
            foreach (var item in labelGenes[0])
                ir.WriteLine(item+" ");
            ir.WriteLine();
            ir.WriteLine("[SAMPLES]");
            foreach (var item in labelSamples[0])
                ir.WriteLine(item + " ");

            ir.Close();
        }
        double [,] MaxValueThreshold(double [,]data,double threshold)
        {
            double[,] newData = data;

            for (int j = 0; j < newData.GetLength(0); j++)
            {
                for (int i = 0; i < newData.GetLength(1); i++)
                {
                    newData[j,i] = data[j,i];
                    if(newData[j,i]>threshold)
                    {
                        newData[j,i] = threshold;
                    }
                }
            }

            return newData;
        }
        double [,] IntervalCoding(double [,] data)
        {
            double[,] newData = data;
            Dictionary<double, int> hashValues = new Dictionary<double, int>();
            double[,] intervals;
            double[,] outData;
            if (oInput.coding == CodingAlg.Z_SCORE)
                newData=ZscoreCodingColumn(data).Item1;
            
            for (int j = 0; j < newData.GetLength(0); j++)
            {
                for (int i = 0; i < newData.GetLength(1)/10; i++)
                {
                    if (!hashValues.Keys.Contains(data[j, i]))
                        hashValues.Add(data[j, i], 1);
                    else
                        hashValues[data[j, i]]++;
                }

                //Array.Sort(colValues);
            }

            intervals = SetupIntervals(hashValues);
                       
            SaveLabIdAndIntervals(intervals);
            outData = IntervalCodig(data, intervals);

            return outData;
        }
        double[, ] SetupIntervals(Dictionary <double,int> dataValues)
        {
            if (oInput.numStates < 3)
                oInput.numStates = 3;
            double [,] intervals=new double [oInput.numStates,2];
            double max, min;
            max = double.MinValue;
            min = double.MaxValue;
//            double[] values = new double[dataValues.Keys.Count];
            double step = 0;
            int i = 0;
            List<double> dlist = new List<double>(dataValues.Keys);
            dlist.Sort();
            min = dlist[0];
            max = dlist[dlist.Count - 1];
            for (int q = 0; q < intervals.GetLength(0); q++)
            {
                intervals[q, 0] = double.MaxValue/2;
                intervals[q, 1] = double.MaxValue/2;
            }


            switch (oInput.coding)
            {
                case CodingAlg.Z_SCORE:
                    int s = 0;                   
                    double st = 3.0 / oInput.numStates;
                    intervals[s, 0] = -100000;
                    intervals[s++, 1] = -st *(int)(oInput.numStates/2);                  
                    for(i=(oInput.numStates-2)/2;i>=1;i--,s++)
                    {
                        intervals[s, 0] = -st * (i + 1);
                        intervals[s, 1] = -st * i;
                    }
                    intervals[s, 0] = -st;
                    intervals[s++, 1] = st;

                    for (i=1;i<=(oInput.numStates-2)/2;i++,s++)
                    {
                        intervals[s, 0] = st * i;
                        intervals[s, 1] = st * (i+1);
                    }
                    intervals[s, 0] = st*(int)(oInput.numStates/2);
                    intervals[s, 1] = 100000;

                    break;
                case CodingAlg.EQUAL_DIST:
                
                    step = (max - min) / oInput.numStates;
                    for (i = 0; i < oInput.numStates; i++)
                    {
                        intervals[i, 0] = dlist[0] + i * step;
                        intervals[i, 1] = dlist[0] + (i + 1) * step;
                    }
                    break;

                case CodingAlg.PERCENTILE:
                    int counter = 0;
                    double end = 0;


                    foreach (var item in dataValues.Values)
                        counter += item;
                    int amount = counter / oInput.numStates;
                    counter = 0;
                    int k = 0;
                    double begin = dlist[0];
                    foreach(var item in dlist)
                    {
                        counter += dataValues[item];
                        if(counter>=amount)
                        {
                            end = item;
                            counter = 0;
                            intervals[k, 0] = begin;
                            if (begin == end)
                                end = end + 0.01 * end;

                            intervals[k, 1] = end;                            
                            begin = intervals[k++, 1];

                        }
                        if(k>=intervals.GetLength(0))
                        {
                            intervals[k - 1, 1] = double.MaxValue/2;
                            break;
                        }
                    }
                    if (counter>0 && counter < amount)
                    {
                        if (k < intervals.GetLength(0))
                        {
                            intervals[k, 0] = begin;
                            end = double.MaxValue/2;
                            intervals[k, 1] = end;
                        }
                        else
                            intervals[k-1, 1] = double.MaxValue/2;
                    }

                    break;


            }
            return intervals;

        }
        static Tuple<double[,], double[], double[]> ZscoreCodingRow(double[,] data)
        {
            double[,] newData = new double[data.GetLength(0), data.GetLength(1)];
            double[] colValues = new double[data.GetLength(1)];
            double[] stdDev = new double[data.GetLength(0)];
            double[] avr = new double[data.GetLength(0)];
            for (int j = 0; j < data.GetLength(0); j++)
            {
                for (int i = 0; i < data.GetLength(1); i++)
                    colValues[i] = data[j, i];

                double sumX2 = 0, sumX = 0;

                foreach (var item in colValues)
                {
                    if (item != double.NaN)
                    {
                        sumX2 += item * item;
                        sumX += item;
                    }
                }
                avr[j] = sumX / colValues.Length;

                stdDev[j] = Math.Sqrt(sumX2 / colValues.Length - avr[j] * avr[j]);

                for (int i = 0; i < data.GetLength(1); i++)
                    if (data[j, i] != double.NaN)
                        newData[j, i] = (data[j, i] - avr[j]) / stdDev[j];

            }
            return new Tuple<double[,], double[], double[]>(newData, avr, stdDev);
        }

        static Tuple<double [,],double[],double[]> ZscoreCodingColumn(double [,] data)
        {
            double[,] newData = new double[data.GetLength(0), data.GetLength(1)];
            double []colValues = new double [data.GetLength(0)];
            double[] stdDev = new double[data.GetLength(1)];
            double[] avr = new double[data.GetLength(1)];
            for (int j = 0; j < data.GetLength(1); j++)
            {
                for (int i = 0; i < data.GetLength(0); i++)                
                    colValues[i]=data[i, j];                

                double sumX2=0, sumX=0;

                foreach(var item in colValues)
                {
                    if (item != double.NaN)
                    {
                        sumX2 += item * item;
                        sumX += item;
                    }
                }
                avr[j] = sumX / colValues.Length;

                stdDev[j] = Math.Sqrt(sumX2 / colValues.Length - avr[j]  * avr[j]);
                
                for(int i=0;i<data.GetLength(0);i++)                
                    if(data[i,j]!=double.NaN)
                        newData[i,j]= (data[i, j] - avr[j]) / stdDev[j];
                               
            }
            return new Tuple<double[,],double[],double[]>(newData,avr,stdDev);
        }
        public static List<string> GetOrderedProfiles(string fileName)
        {
            List<string> order = new List<string>();
            StreamReader st = new StreamReader(fileName);
            string line = st.ReadLine();

            while(line!=null)
            {
                if(line.Contains(">"))
                {
                    string name = line.Replace(">", "");
                    order.Add(name);
                }
                line = st.ReadLine();
            }

            return order;
        }
        public override Dictionary<string, protInfo> GetProfile(profileNode node, string fileNameProf)
        {
            Dictionary<string, protInfo> data;
            profileNode localNode;
            string fileName = GetProfileName(fileNameProf);

            if (oInput.heatmap || oInput.transpose)
                fileName = fileName + "_transpose";

            ProfileTree ts = new ProfileTree();
            string locPath = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location) + Path.DirectorySeparatorChar + "profiles" + Path.DirectorySeparatorChar;
            ts.LoadProfiles(locPath + Path.GetFileNameWithoutExtension(fileNameProf) + ".profiles");
            List<string> k = new List<string>(ts.masterNode.Keys);
            localNode = ts.masterNode[k[0]];
            if (node.profWeights.Count == localNode.profWeights.Count)
            {
                data = base.GetProfile(localNode, fileName);
            }
            else
            {
                ts.LoadProfiles(locPath + Path.GetFileNameWithoutExtension(fileNameProf) + "_distance.profile");
                k = new List<string>(ts.masterNode.Keys);
                localNode = ts.masterNode[k[0]];
                data = base.GetProfile(localNode, fileName);
            }
         
            return data;
        }
    }
}
