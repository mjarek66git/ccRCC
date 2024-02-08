using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using phiClustCore;

namespace WorkFlows
{
    public interface IclusterType
    {
        void SetProfileName(string name);
       void Show();
       INPUTMODE GetInputType();
       string ToString();
    }
}
