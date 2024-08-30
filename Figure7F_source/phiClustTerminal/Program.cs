using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Reflection;
using System.Linq.Expressions;
using System.Timers;
using phiClustCore;
using phiClustCore.Interface;
using phiClustCore.Distance;
using phiClustCore.Profiles;



namespace singleCellAnalisysTerminal
{
    class SingleCellAnalysis
    {
        JobManager manager = new JobManager();
        LoadSetup lSet;
        public bool errors = false;
        public bool times = false;
        public bool progress = false;
        public bool graphics = false;
        public string setupFile = "setup.dat";
        public void ErrorMessage(string message)
        {
            Console.WriteLine(message);
        }


        private void UpdateProgress(object sender, EventArgs e)
        {

            Dictionary<string, double> res = manager.ProgressUpdate();

            if (res == null)
            {
                //                TimeInterval.Stop();
                Console.Write("\r                                                             ");
                return;
            }
            Console.Write("\r                                                                              ");
            foreach (var item in res)
                Console.Write("\rProgress " + item.Key + " " + (item.Value * 100).ToString("0.00") + "%");

        }


        public void Run()
        {

            Options opt = new Options();
            ClusterVis clusterOut = new ClusterVis();

            lSet = new LoadSetup();

            lSet.ReadSetupFile(setupFile);


            Dictionary<string, ClusterOutput> clOut = new Dictionary<string, ClusterOutput>();
            Random r = new Random();
            
            opt.clusterAlgorithm.Add(ClusterAlgorithm.OmicsHeatMap);
            opt.omics.processName = "Batch_" + r.Next(1000);
            opt.hierarchical.dummyProfileName = lSet.clusterFileName;
            opt.hash.reqClusters = lSet.numUpperCl;
            opt.hash.relClusters = lSet.numLeftCl;
            opt.hierarchical.distance = lSet.dist;
            manager.opt = opt;
            manager.message = ErrorMessage;
            if (progress)
            {
                TimeIntervalTerminal.InitTimer(UpdateProgress);
                TimeIntervalTerminal.Start();
            }
            manager.RunJob(opt.omics.processName, lSet.joined);
            manager.WaitAllNotFinished();
            UpdateProgress(null, null);
            HeatmapDrawCore hCore = new HeatmapDrawCore(manager.clOutput[opt.omics.processName + ";1"], lSet.resolutionWidth, lSet.resolutionHeight);
            hCore.draw.PrepareDataForHeatMap();
            if (lSet.selectedGenesFile.Length > 0)
                hCore.ReadNoteGenes(lSet.selectedGenesFile);
            if (lSet.rowAnnotationFile.Length > 0)
                hCore.ReadLeftAnnotation(lSet.rowAnnotationFile);
            if (lSet.simplifyUD)
                hCore.LongLeavesClick();
            if (lSet.horizontalLines)
                hCore.horizontalLines = true;
            hCore.Save(lSet.outputFileName, lSet.resolutionWidth, lSet.resolutionHeight);
            if (progress)
                TimeIntervalTerminal.Stop();
            Console.Write("\r                                                                     ");
            clOut = manager.clOutput;

            if (times)
            {
                foreach (var item in clOut)
                    Console.WriteLine(item.Value.dirName + " " + item.Value.measure + " " + item.Value.time);
            }
            if (errors)
            {
                foreach (var item in ErrorBase.GetErrors())
                    Console.WriteLine(item);
            }
            Console.WriteLine();
        }

    }

    class Program
    {        

        static void Main(string[] args)
        {
            SingleCellAnalysis single = new SingleCellAnalysis();            
            if (args.Length == 0)
            {
                Console.WriteLine("Following argument is required:");
                Console.WriteLine("-f configuration_file");
                Console.WriteLine("Following options may be specified");
                Console.WriteLine("-e \n\t show all errors");
                Console.WriteLine("-n \n\t number of cores to be used");
                Console.WriteLine("-t \n\tshow time information");
                Console.WriteLine("-sg fileName resWidthxresHeigth \n\tSave results to png file (if possible)");
                Console.WriteLine("-p \n\tShow progres bar");
                return;
            }
            Settings set = new Settings();
            set.Load();
            if (set.profilesDir == null || set.profilesDir.Length == 0)
            {
                set.profilesDir = "generatedProfiles";
                //    set.Save();
            }
            //set.mode = INPUTMODE.PROTEIN;
            set.mode = INPUTMODE.USER_DEFINED;
            for (int i = 0; i < args.Length; i++)
            {
                switch (args[i])
                {
                    case "-f":
                        if (i + 1 >= args.Length)
                        {
                            Console.WriteLine("After -f option you have to provide configuration file");
                            return;
                        }
                        if (!File.Exists(args[i + 1]))
                        {
                            Console.WriteLine("File " + args[i + 1] + " does not exist");
                            return;
                        }
                        single.setupFile = args[i + 1];
                        i++;
                        break;
                    case "-n":
                        if (args.Length > i)
                        {
                            int num;
                            try
                            {
                                num = Convert.ToInt32(args[++i]);
                                set.numberOfCores = num;
                            }
                            catch (Exception ex)
                            {
                                Console.WriteLine("Wrong definition of number of cores: " + ex.Message);
                                return;
                            }
                        }
                        else
                            Console.WriteLine("Number of cores has been not provided");
                        break;
                    case "-e":
                        single.errors = true;
                        break;
                    case "-t":
                        single.times = true;
                        break;
                    case "-p":
                        single.progress = true;
                        break;
                    default:
                        if (args[i].Contains("-"))
                            Console.WriteLine("Unknown option " + args[i]);
                        break;

                }
            }
            set.Save();
            Console.WriteLine("Configuration file " + single.setupFile);
            single.Run();
        }
    }
}

