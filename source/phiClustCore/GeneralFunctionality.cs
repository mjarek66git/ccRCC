using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using phiClustCore.Interface;
using System.Runtime.Serialization.Formatters.Binary;
using System.Threading.Tasks;

namespace phiClustCore
{
    public class GeneralFunctionality
    {
        static public void SaveBinary(string fileName,object o)
        {
            Stream stream = File.Open(fileName, FileMode.Create);
            BinaryFormatter bFormatter = new BinaryFormatter();
            bFormatter.Serialize(stream, o);
            stream.Close();

        }
        static public ISerialize LoadBinary(string fileName)
        {
            ISerialize outObject;
            Stream stream = File.Open(fileName, FileMode.Open);
            BinaryFormatter bFormatter = new BinaryFormatter();
            outObject = (ISerialize)bFormatter.Deserialize(stream);

            return outObject;
        }

    }
}
