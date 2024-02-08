using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore
{
    class SetupException:Exception
    {
        public SetupException() { }
        public SetupException(string msg) : base(msg) { }
        public SetupException(string msg, Exception inner) : base(msg, inner) { }
    }
    public class LoadSetup
    {
        List<OmicsDataSet> data = new List<OmicsDataSet>();
        Dictionary<string, FilterOmics> availableFilters = new Dictionary<string, FilterOmics>();
        BindingList<FilterOmics> filters = new BindingList<FilterOmics>();
        List<OmicsDataSet> omicsFiltered = new List<OmicsDataSet>();
        public OmicsDataSet joined = null;
        public int numUpperCl = 0;
        public int numLeftCl = 0;
        public int resolutionWidth = 0;
        public int resolutionHeight = 0;
        public string clusterFileName = "";
        public string outputFileName = "";
        public bool simplifyUD = false;
        public bool horizontalLines = false;
        public string selectedGenesFile = "";
        public string rowAnnotationFile = "";
        public DistanceMeasures dist = DistanceMeasures.PEARSON;
        public LoadSetup()
        {
            Type[] xx = Assembly.GetAssembly(typeof(FilterOmics)).GetTypes();
            IEnumerable<Type> subclasses = xx.Where(t => t.IsSubclassOf(typeof(FilterOmics)));
            foreach (var item in subclasses)
            {
                FilterOmics aux = (FilterOmics)Activator.CreateInstance(item);
                availableFilters.Add(aux.ToString(), aux);
            }

        }

        public void ReadSetupFile(string fileName)
        {
            StreamReader sr = new StreamReader(fileName);
            string line = sr.ReadLine();
            try
            {
                while (line != null)
                {
                    switch (line)
                    {
                        case "[FILES]":
                            line = ReadFiles(sr);
                            break;
                        case "[FILTERS FOR JOINED DATA]":
                            line = ReadFilters(sr);
                            break;
                        case "[CLUSTER]":
                            line = ReadClusterInfo(sr);
                            break;
                        case "[OUTPUT]":
                            line = ReadOutputOptions(sr);
                            break;


                    }
                }
                sr.Close();
            }
            catch(SetupException s)
            {
                Console.WriteLine("Something wrong in Setup file in line: " + s);
                sr.Close();
                Environment.Exit(0);
            }
            
            if (data.Count > 1)
                joined = OmicsDataSet.JoinOmicsData(data);
            else
                joined = data[0];

            joined=joined.ApplyFilters(filters);
/*            OmicsDataSet aux = joined;
            foreach (var item in filters)
            {
                aux = item.ApplyFilter(aux);
            }*/
        }
        string ReadOutputOptions(StreamReader sr)
        {
            string line = sr.ReadLine();
            string[] aux = null;
            try
            {
                while (line != null && !line.Contains("["))
                {
                    switch (line)
                    {
                        case var _ when line.Contains("File"):
                            aux = line.Split('|');
                            this.outputFileName = aux[1];
                            break;
                        case var _ when line.Contains("Simpl"):
                            aux = line.Split('|');
                            if (aux[1].Contains("true"))
                                this.simplifyUD = true;
                            else
                                this.simplifyUD = false;
                            break;
                        case var _ when line.Contains("annotation"):
                            aux = line.Split('|');
                            this.rowAnnotationFile = aux[1];
                            break;
                        case var _ when line.Contains("Selected"):
                            aux = line.Split('|');
                            this.selectedGenesFile = aux[1];
                            break;
                        case var _ when line.Contains("Horizontal"):
                            aux = line.Split('|');
                            if (aux[1].Contains("true"))
                                this.horizontalLines = true;
                            else
                                this.horizontalLines = false;
                            break;
                        case var _ when line.Contains("Resolution"):
                            aux = line.Split('|');
                            if (aux[1].Contains("x"))
                            {
                                string[] tmp = aux[1].Split('x');
                                resolutionWidth = Convert.ToInt32(tmp[0]);
                                resolutionHeight = Convert.ToInt32(tmp[1]);
                            }
                            break;
                        default:
                            throw new SetupException("Unrecognized option");
                            

                    }
                    line = sr.ReadLine();
                }
            }
            catch(Exception ex)
            {
                throw new SetupException(ex+" "+line);
            }
            return line;
        }
        string ReadClusterInfo(StreamReader sr)
        {
            string line = sr.ReadLine();
            string[] aux = null;
            try
            {
                while (line != null && !line.Contains("["))
                {
                    switch (line)
                    {
                        case var _ when line.Contains("upper clusters"):
                            aux = line.Split('|');
                            numUpperCl = Convert.ToInt32(aux[1]);
                            break;
                        case var _ when line.Contains("Number of left clusters"):
                            aux = line.Split('|');
                            numLeftCl = Convert.ToInt32(aux[1]);
                            break;
                        case var _ when line.Contains("File"):
                            aux = line.Split('|');
                            this.clusterFileName = aux[1];
                            break;
                        case var _ when line.Contains("Distance"):
                            aux = line.Split('|');
                            Enum.TryParse<DistanceMeasures>(aux[1], out dist);
                            break;
                        default:
                            throw new SetupException("Unrecognized option");

                    }

                    line = sr.ReadLine();
                }
            }
            catch(Exception ex)
            {
                throw new SetupException(ex.Message+" "+line);
            }
            return line;
        } 
        string ReadFilters(StreamReader sr)
        {
            string line = sr.ReadLine();
            string[] auxS = null;
            try
            {


                while (line != null && !line.Contains("["))
                {
                    auxS = null;
                    string filterName = line;
                    if (line.Contains("|"))
                    {
                        auxS = line.Split('|');
                        filterName = auxS[0];
                    }
                    if (availableFilters.ContainsKey(filterName))
                    {
                        FilterOmics aux = (FilterOmics)availableFilters[filterName].Clone();
                        if (auxS != null)
                        {
                            auxS[1] = auxS[1].Trim();
                            string[] param = auxS[1].Split(' ');
                            Dictionary<string, string> paramDic = new Dictionary<string, string>();
                            foreach (var item in param)
                            {
                                string itemT = item.Trim();
                                string[] tmp = itemT.Split('=');
                                if (tmp[0].Contains("_"))
                                    tmp[0] = tmp[0].Replace('_', ' ');
                                paramDic.Add(tmp[0], tmp[1]);
                            }
                            aux.SetParameters(paramDic);
                        }

                        filters.Add(aux);
                    }
                    else
                        throw new Exception("Filter: " + line + " not avaiable");
                    line = sr.ReadLine();
                }
            }
            catch
            {
                throw new SetupException(line);
            }
            return line;
        }
        string ReadFiles(StreamReader sr)
        {

            OmicsDataSetSetup setup = new OmicsDataSetSetup();
            setup.colPos = 1;
            setup.columnFlag = true;
            setup.labelRow = 0;
            setup.rowPos = 1;
            setup.startCol = 0;
            setup.separators = new char[] { '\t' };

            string line = sr.ReadLine();
            while(line!=null && !line.Contains("["))
            {

                ReadOmicsFile file = new ReadOmicsFile(setup);

                OmicsDataSet aux = file.GetOmicsFile(line);
                data.Add(aux);
                line = sr.ReadLine();
            }
            return line;
        }
    }
}
