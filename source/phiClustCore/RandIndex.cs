using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Data;
using phiClustCore.Interface;

namespace phiClustCore
{
    public class RandIndex
    {
        int[] dataEx;        
      
        DataTable resTable=null;
        public DataTable resT
        {
            get { return resTable; }
        }
        public RandIndex()
        {
            resTable = new DataTable();
            resTable.Columns.Add("Considered Clusters", typeof(string));
            resTable.Columns.Add("Clusters", typeof(string));
            resTable.Columns.Add("Rand Index", typeof(string));
            resTable.Columns.Add("Cluster Index", typeof(string));
           
        }
        public long  Size(List<List<string>> _out1, List<List<string>> _out2)
        {
            long  count1=0,count2=0;
            foreach (var item in _out1)
                count1 += item.Count;

            foreach (var item in _out2)
                count2 += item.Count;

            return  count1* (count2 + 1) / 2;
        }
        public void ClusterDistance(List<List<string>> _out1, List<List<string>> _out2,string name,ref long currentV,int consideredClusters)
        {
            int count1 = 0, count2 = 0;
            int maxx = 0;
            long remCurrent = currentV;
            Dictionary<string, int> allData = new Dictionary<string, int>();
            List<List<string>> out1 = new List<List<string>>();
            List<List<string>> out2 = new List<List<string>>();

            double P1 = 0, P2 = 0;
            double ED, EA;
            int i = 0;

         

            foreach (var item in _out1)
                foreach (var itemW in item)
                    allData.Add(itemW, 1);

            foreach (var item in _out2)
                foreach (var itemW in item)
                    if (allData.ContainsKey(itemW))
                        allData[itemW]++;
                    else
                        allData.Add(itemW, 1);

            foreach (var item in _out1)
            {
                List<string> aux = new List<string>();
                foreach (var itemW in item)
                    if (allData[itemW] == 2)
                        aux.Add(itemW);

                if (aux.Count > 0)
                    out1.Add(aux);
            }
            foreach (var item in _out2)
            {
                List<string> aux = new List<string>();
                foreach (var itemW in item)
                    if (allData[itemW] == 2)
                        aux.Add(itemW);

                if (aux.Count > 0)
                    out2.Add(aux);
            }
            if (out1.Count == 0 || out2.Count==0)
                throw new Exception("It looks that clusterization has been made on different data");

            allData.Clear();
            int counter = 0;
            foreach (var item in out1)
                foreach (var itemW in item)
                    allData.Add(itemW, counter++);

            foreach (var item in out1)
                count1 += item.Count;

            foreach (var item in out2)
                count2 += item.Count;

            maxx = count1 * (count2 + 1) / 2;
            dataEx = new int[count1 * (count2 + 1) / 2];

            for (i = 0; i < dataEx.GetLength(0); i++)
                dataEx[i] = 0;

            foreach (var item in out1)
                P1 += item.Count * (item.Count - 1);

            P1 = P1 / (count1 * (count1 - 1));

            foreach (var item in out2)
                P2 += item.Count * (item.Count - 1);

            P2 = P2 / (count2 * (count2 - 1));

            double pairs = count1 * (count1 - 1) / 2;

            ED = pairs * (P1 * (1 - P2) + P2 * (1 - P1));
            EA = pairs * P1 * P2;


            int a =0;//number of pairs that are in the same set in X and Y
            int b=0; //number of pairs that are in different set in X and Y
            int c=0;//number of pairs that are in the same set in X but in different sets in Y
            int d=0; //number of pairs that are in the same set in Y but in different sets in X

            for (int l = 0; l < out1.Count; l++)
            {
                for (int n = 0; n < out1[l].Count; n++)
                    for (int k = n + 1; k < out1[l].Count; k++,currentV++)
                        dataEx[GetIndex(allData.Count, allData[out1[l][n]], allData[out1[l][k]])] = 1;


            }

            for (int l = 0; l < out2.Count; l++)
            {
                for (int n = 0; n < out2[l].Count; n++)
                    for (int k = n + 1; k < out2[l].Count; k++)
                        if (dataEx[GetIndex(allData.Count, allData[out2[l][n]], allData[out2[l][k]])] == 0)
                        {
                            dataEx[GetIndex(allData.Count, allData[out2[l][n]], allData[out2[l][k]])] = -1;
                            d++;
                            currentV++;
                        }
                        else
                            a++; 
            }
            for (int l = 0; l < out2.Count; l++)
            {
                for (int s = l+1; s < out2.Count; s++)
                {
                    for (int n = 0; n < out2[l].Count; n++)
                        for (int k = 0; k < out2[s].Count; k++)
                            if (dataEx[GetIndex(allData.Count, allData[out2[l][n]], allData[out2[s][k]])] == 1)
                            {
                                c++;
                            }
                            else
                            {
                                if (dataEx[GetIndex(allData.Count, allData[out2[l][n]], allData[out2[s][k]])] != -1)
                                {
                                    b++;
                                    currentV++;
                                }
                            }
                }
            }

            

//            result.clusterDist = (c+d) / ED * EA /a;
//            result.randIndex = (a + b) /(float) (a + b + c + d);
            if(consideredClusters>0)
                resTable.Rows.Add(consideredClusters.ToString(),name, String.Format("{0:0.####}", (c + d) / ED * EA / a), String.Format("{0:0.####}", (a + b) / (float)(pairs)));
            else
                resTable.Rows.Add("All",name, String.Format("{0:0.####}", (c + d) / ED * EA / a), String.Format("{0:0.####}", (a + b) / (float)(pairs)));

//            return result;
            currentV=maxx+remCurrent;

        }
        int GetIndex(int n,int row, int col)
        {
            int index=0;
            int step = n;
            int x=row;
            if (row > col)
            {
                row = col;
                col = x;
            }

            index = row * (step + step - row + 1) / 2;
            index +=col-row;
            return index;
        }
    }
}
