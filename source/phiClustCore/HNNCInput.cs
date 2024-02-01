using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore
{
    [Serializable]
    public class HNNCInput : BaseCInput
    {
        [Description("Clustering method (Rpart or Hash")]
        public bool Rpart = false;
        [Description("Labels for items in train file")]
        public string labelsFile = "";
        [Description("Test file")]
        public string testFile= "";
    }
}
