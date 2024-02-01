using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Text.RegularExpressions;

namespace phiClustCore.Profiles
{
    class UserDefWeighted:UserDefinedProfile
    {
        public override Dictionary<string, protInfo> GetProfile(profileNode node, string fileName)
        {
            Dictionary<string, protInfo> dic = new Dictionary<string, protInfo>();
            StreamReader wr;
            DebugClass.WriteMessage("profile" + fileName);
            wr = new StreamReader(fileName);
            
            protInfo info;
            string line = wr.ReadLine();
            string name = "";
            string seq = "";
            List<string> profile = new List<string>();
            List<byte> newProfile = new List<byte>();
            node.rangeWeights=new List<RangeWeight>();
            bool testRange=false;
            int start=0;
            int end=0;
            while (line != null)
            {
                if (line.Contains(">"))
                {
                    if (name.Length > 0)
                    {
                        info = new protInfo();
                        info.sequence = null;
                        info.profile = newProfile;
                        if (dic.ContainsKey(name))
                            throw new Exception("The nameof profile must be unique, name: " + name + " already exists in " + fileName);
                        dic.Add(name, info);
                    }

                    name = line.Replace(">", "");
                    line = wr.ReadLine();
                }
                if (line.Contains(node.profName+" "))
                {
                    profile.Clear();

                    string cLine=line.Remove(0, (node.profName+" profile ").Length);
                    do
                    {
                        cLine = Regex.Replace(cLine, @"\s+", " ");
                        cLine = cLine.Trim();
                        string[] aux;
                        if (cLine.Contains(' '))
                            aux = cLine.Split(' ');
                        else
                        {
                            char[] charArray = cLine.ToCharArray();
                            aux = cLine.Select(x => x.ToString()).ToArray();                        
                        }

                        if(!testRange)
                        {
                            RangeWeight rw=new RangeWeight();
                            rw.start=start;
                            end+=aux.Length-1;
                            rw.stop=end;
                            start+=aux.Length;
                            rw.Weight=0;
                            node.rangeWeights.Add(rw);
                        }
    
                        foreach (var item in aux)
                            profile.Add(item);

                    
                        newProfile = new List<byte>();
                        for (int i = 0; i < profile.Count; i++)
                            if (node.ContainsState(profile[i].ToString()))
                                newProfile.Add(node.GetCodedState(node.states[profile[i].ToString()]));
                            else
                                if(profile[i].ToString()!="-" && profile[i].ToString()!="X")
                                    throw new Exception("Unknow state " + profile[i].ToString() + " in "+ node.profName + " profile!");
                                else
                                    newProfile.Add(0);
                       line = wr.ReadLine();
                    }
                    while(line.Contains(">"));
                    testRange=true;
                }
                else
                    line = wr.ReadLine();
            }
            info = new protInfo();
            info.sequence = seq;
            info.profile = newProfile;
            dic.Add(name, info);
            DebugClass.WriteMessage("number of profiles " + dic.Keys.Count);

            wr.Close();            

            return dic;
        

        }
    }
}
