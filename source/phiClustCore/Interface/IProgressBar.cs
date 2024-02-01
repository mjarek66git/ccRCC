using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Data;

namespace phiClustCore.Interface
{
    public interface IProgressBar
    {
        double ProgressUpdate();
        List<KeyValuePair<string,DataTable>> GetResults();
        Exception GetException();
        double StartProgress { get; set; }
        double EndProgress { get; set; }
    }
}
