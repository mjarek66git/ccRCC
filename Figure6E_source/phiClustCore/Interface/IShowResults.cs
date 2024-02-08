using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Data;

namespace phiClustCore.Interface
{
    public interface IShowResults
    {
        void Show(List<KeyValuePair<string,DataTable>> T);
        void ShowException(Exception ex);
    }
}
