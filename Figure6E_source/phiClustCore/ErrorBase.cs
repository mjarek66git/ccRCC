using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Reflection;

namespace phiClustCore
{
    public static class ErrorBase
    {
        private static List<string> errors = new List<string>();
        public static void ClearErrors()
        {
            errors.Clear();
        }
        public static void AddErrors(string error)
        {
            errors.Add(error);
        }
        public static List<string> GetErrors()
        {
            return errors;
        }
    }
}
