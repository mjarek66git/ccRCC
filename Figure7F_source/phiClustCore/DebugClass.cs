using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace phiClustCore
{
    public class DebugClass
    {
        static bool DEBUG = false;
        public static void WriteMessage(string message)
        {
            if(DEBUG)
                Console.WriteLine(message);
        }
        public static void DebugOn()
        {
            DEBUG = true;
        }
        public static void DebugOff()
        {
            DEBUG = false;
        }
    }
}
