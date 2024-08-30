using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore
{
    public enum CodingAlg
    {
        PERCENTILE,
        Z_SCORE,
        EQUAL_DIST,
        MAX_VALUE,
        NONE
    };
    [Serializable]
    public class OmicsInput :BaseCInput
    {
        public static string fileName = "Test.opt";
        [Description("Read data from column")]
        public int numCol = 10;
        [Description("Read data from row")]
        public int numRow = 18;
        [Description("Gene Position Rows")]
        public bool genePosition;
        [Description("Number of selected top best genes")]
        public int selectGenes = 0;
        [Description("Z-score")]
        public bool zScore;
        [Description("Quantile")] 
        public bool quantile;
        [Description("Number of states")]
        public int numStates = 6;
        [Description("Label gene positions")]
        public string labelGeneStartString = null;
        [Description("Label sample positions")]
        public string labelSampleStartString = null;
        [Description("Use gene labels")]
        public bool uLabelGene = false;
        [Description("Use sample labels")]
        public bool uLabelSample = false;
        [Description("OutputName")]
        public string processName;
        [Description("Heatmap (do not touch this!)")]
        public bool heatmap = false;
        [Description("File Selected Genes")]
        public string fileSelectedGenes = "";
        [Description("Transposition")]
        public bool transpose = false;
        [Description("Coding Algorithm")]        
        public CodingAlg coding = CodingAlg.EQUAL_DIST;
    }
}
