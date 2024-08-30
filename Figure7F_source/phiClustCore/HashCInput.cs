using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using phiClustCore.Interface;
using phiClustCore.Profiles;

namespace phiClustCore
{
    public enum COL_SELECTION
    {
        ENTROPY,
        META_COL
    };
    [Serializable]
    public class HashCInput:BaseCInput,IAutomaticProfiles
    {
        [Description("Hash clustering: Do you want to use regularization")]
        public bool regular=false;
        [Description("Regularization: size of the window ")]
        public int wSize=7;
        [Description("Regularization: Threshold distance below the threshold will lead to consensus state ")]
        public int regThreshold = 1;
        [Description("Combine clusters based on the hamming distance between them")]
        public bool combine = true;
        [Description("Use consensus states if true, otherwise orginal coding will be used")]
        public bool useConsensusStates = true;
        [Description("Use 1Djury to find consesnsus states ")]
        public bool jury=true;
        [Description("Profile name")]
        public string profileName="profiles/SS-SA.profiles";
        [Description("Profile name for regularization")]
        public string profileNameReg="profiles/SS-SA.profiles";
        [Description("Find clusters by automatic columns selection")]
        public bool hashCluster = false;
        [Description("Columns selection method")]
        public COL_SELECTION selectionMethod=COL_SELECTION.ENTROPY;
        [Description("Number of relevant clusters")]
        public int relClusters=100;
        [Description("Percent of data that must be in the relevant clusters")]
        public int perData=90;
        [Description("Number of requested clusters")]
        public int reqClusters = 1000;
        [Description("Number of requested reference Points")]
        public int refPoints = 5;
        [Description("Select refernce structure in the cluster")]
        public bool selReference=true;
        [Description("Use org data without projection")]
        public bool projection = true;

        public void GenerateAutomaticProfiles(string fileName)
        {
            ProfileTree t = ProfileAutomatic.AnalyseProfileFile(fileName, SIMDIST.SIMILARITY);
            string profileName;

            profileName = ProfileAutomatic.similarityProfileName;
            this.profileName = profileName;
            this.profileNameReg = profileName;
            t.SaveProfiles(profileName);
        }
        public override string GetVitalParameters()
        {
            string outLine="";
            if(combine)
                outLine="Algorithm: Rprop==";
            else
                outLine="Algorithm: Hash ==Column selection: "+selectionMethod;
            outLine += "== Profile: " + profileName+"== Number of relevant clusters: "+relClusters+"== Percent of the data: "+perData;
            if(regular)
                outLine+="== Regularization:  ON ==Window size: "+wSize+"== Threshold: "+regThreshold+"]";

            return outLine;

        }

    }
}
