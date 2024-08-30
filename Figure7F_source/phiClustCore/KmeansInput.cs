using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using phiClustCore.Interface;
using phiClustCore.Profiles;

namespace phiClustCore
{
    [Serializable]
    public class KmeansInput:BaseCInput,IAutomaticProfiles
    {
        [Description("K-means algorithm: K in K-means")]
        public int maxK=10;
        [Description("K-means algorithm: How To initialize k-means (RANDOM,JURY1D)")]
        public Initialization kMeans_init=Initialization.Jury1D;
        [Description("K-means algorithm: name of the file with alignment needed when 1djury for initialization is used")]
        public string alignFileName;
        [Description("If you want to use automaticly generated alignment set this variable to true")]
        public bool alignGenerate = true;
        [Description("Distance measures used for clustering")]
        public DistanceMeasures kDistance = DistanceMeasures.HAMMING;   
        [Description("Path to the profile for weighted hamming distance [KMEANS]")]
        public string hammingProfile="profiles/SS-SA_distance.profile";
        [Description("Path to the profile for 1djury for initialization[KMEANS]")]
        public string jury1DProfile="profiles/SS-SA.profiles";
        [Description("Max number of iterations [KMEANS]")]
        public int maxIter=300;
        [Description("Use 1DJury to find reference vectors")]
        public bool reference1Djury =false;

        public void GenerateAutomaticProfiles(string fileName)
        {
            ProfileTree t = ProfileAutomatic.AnalyseProfileFile(fileName, SIMDIST.DISTANCE);
            string profileName = ProfileAutomatic.distanceProfileName;
            t.SaveProfiles(profileName);
            hammingProfile = profileName;
            t = ProfileAutomatic.AnalyseProfileFile(fileName, SIMDIST.SIMILARITY);
            profileName = ProfileAutomatic.similarityProfileName;
            t.SaveProfiles(profileName);
            jury1DProfile = profileName;
        }

        public override string GetVitalParameters()
        {
            string outLine = "";
            outLine = "Algorithm: Kmeans==";
            outLine += " K: " + maxK + "== Initialization: " + kMeans_init+"== MaxIterations: "+maxIter;
            if (kMeans_init == Initialization.Jury1D)
                outLine += "== profile: " + jury1DProfile;
            outLine += "== Distance measure: " + kDistance;
            if (kDistance == DistanceMeasures.COSINE || kDistance == DistanceMeasures.HAMMING)
                outLine += "== Profile: " + hammingProfile;

            return outLine;
        }
    }
}
