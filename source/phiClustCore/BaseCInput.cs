using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Reflection;
using System.Linq.Expressions;
using System.IO;

namespace phiClustCore
{
    [Serializable]
    public class BaseCInput
    {
        [NonSerialized]
        private Dictionary<string, string> dicField = new Dictionary<string, string>();
        [NonSerialized]
        private Dictionary<string, MemberInfo> dicMem = new Dictionary<string, MemberInfo>();

        //[Description("Use 1DJury to find reference vectors")]
        //public bool reference1Djury;
        public string alignmentFileName;

        public bool dummyProfile = false;
        public string dummyProfileName = "";

        public virtual string GetVitalParameters()
        {
            return "Parameters protocol not defined for ths clusterization!";
        }
        private void PrepareAllFields()
        {
            Type t = this.GetType();
            MemberInfo[] members = t.GetMembers();
            dicField.Clear();
            dicMem.Clear();

            foreach (MemberInfo mem in members)
            {
                object[] attributes = mem.GetCustomAttributes(true);

                if (attributes.Length != 0 && mem.MemberType.ToString() == "Field")
                {
                    string key = "";
                    foreach (object attribute in attributes)
                    {

                        //Console.Write("  {0} ", attribute.ToString());
                        DescriptionAttribute da = attribute as DescriptionAttribute;
                        if (da != null)
                            key = da.Description;
                    }
                    object o = mem.ReflectedType.GetField(mem.Name).GetValue(this);
                    if (o != null)
                        dicField.Add(key, mem.ReflectedType.GetField(mem.Name).GetValue(this).ToString());
                    else
                        dicField.Add(key, "");

                    dicMem.Add(key, mem);

                }
            }
        }
        public void SaveOptions(StreamWriter fileStream)
        {
            PrepareAllFields();
            foreach (var item in dicField.Keys)
                fileStream.WriteLine(item + "# " + dicField[item]);

        }
        public void ReadOptionFile(StreamReader fileRead)
        {
                PrepareAllFields();
                string line="";
                MemberInfo memB;
                while (!fileRead.EndOfStream && !line.Contains("End"))
                {
                    line = fileRead.ReadLine();
                    line = line.Replace("# ", "#");
                    string[] strTab = line.Split('#');

                    if (dicField.ContainsKey(strTab[0]))
                    {
                        memB = dicMem[strTab[0]];
                        string ww = memB.ReflectedType.GetField(memB.Name).FieldType.Name;
                        switch (ww)
                        {
                            case "string":
                            case "String":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, strTab[1]);
                                break;
                            case "Boolean":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Convert.ToBoolean(strTab[1]));
                                break;
                            case "Single":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Convert.ToSingle(strTab[1]));
                                break;
                            case "Int32":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Convert.ToInt32(strTab[1]));
                                break;
                            case "Initialization":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Enum.Parse(typeof(Initialization), strTab[1]));
                                break; 
                            case "ClusterAlgorithm":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Enum.Parse(typeof(ClusterAlgorithm), strTab[1]));
                                break;
                            case "DistanceMeasures":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Enum.Parse(typeof(DistanceMeasures), strTab[1]));
                                break;
                            case "AglomerativeType":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Enum.Parse(typeof(AglomerativeType), strTab[1]));
                                break;
                            case "CodingAlg":
                                memB.ReflectedType.GetField(memB.Name).SetValue(this, Enum.Parse(typeof(CodingAlg), strTab[1]));
                                break;

                    }

                        //SetValue(this, strTab[1]);
                    }
                    //else
                        //DebugMode.WriteMessage("Not recognized: " + strTab[0]);
                }

        }

    }
}
