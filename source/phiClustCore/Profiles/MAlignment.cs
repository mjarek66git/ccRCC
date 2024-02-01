using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace phiClustCore
{
    public class MAlignment: IDisposable
    {
        Dictionary<string, string> seqProt;
        Dictionary<string, Dictionary<string, int>> costMatrix = new Dictionary<string, Dictionary<string, int>>();
		cell[,] scoreTab;
		cell[,] scoreTab1;
		cell[,] scoreTab2;
        List<char> result1; 
        List<char> result2; 

        int size1, size2;
        int gapPenalty = 8;
		public struct alignSeq
		{
			public string seq1;
			public string seq2;
		};
        void IDisposable.Dispose()
        {
            GC.SuppressFinalize(this);
        }

        public MAlignment(Dictionary<string, string> seqProt)
        {
            this.seqProt = seqProt;            
        }
        public MAlignment(int size1)
        {
            result1=new List<char>(size1*size1);
            result2=new List<char>(size1*size1);

            size1++;
            scoreTab = new cell[size1, size1];
            scoreTab1 = new cell[size1, size1];
            scoreTab2 = new cell[size1, size1];

            for (int i = 0; i < size1; i++)
                for (int j = 0; j < size1; j++)
                {
                    scoreTab[i, j] = new cell();
                    scoreTab1[i, j] = new cell();
                    scoreTab2[i, j] = new cell();
                }
        }
		

        void MakeAlignment(string baseSeq)
        {
            foreach (var item in seqProt.Keys)
            {
                if (item != baseSeq)
                {
                    Align(baseSeq, item);
                }
            }
        }
        alignSeq TraceBack(string str1,string str2)
        {
            result1.Clear();
            result2.Clear();
			//string result1="";
			//string result2="";
			alignSeq alignRes;
			int val1,val2,val3;
			val1=scoreTab[size1 - 1, size2 - 1].value;							
			val2=scoreTab1[size1 - 1, size2 - 1].value;
			val3=scoreTab2[size1 - 1, size2 - 1].value;
			
			cell aux=scoreTab[size1 - 1, size2 - 1];
			if(val1>val2 && val1>val3)										
            	aux = scoreTab[size1 - 1, size2 - 1];							
			else
				if(val2>=val1 && val2>=val3)										
    	        	aux = scoreTab1[size1 - 1, size2 - 1];							
				else
					aux = scoreTab2[size1 - 1, size2 - 1];							
				
            do
            {        					
                if(aux.direction==2)
				{
                    result1.Add(str1[aux.ai-1]);
					result2.Add(str2[aux.aj-1]);
				}
				else
					if(aux.direction==0)				
					{
						result2.Add('-');
						result1.Add(str1[aux.ai-1]);
					}
					else
					{
						result1.Add('-');
						result2.Add(str2[aux.aj-1]);					
					}
				aux=aux.tab[aux.i, aux.j];				
                
            }
            while (aux.ai != 0 && aux.aj != 0);


            while (aux.ai > 0)
            {
                result1.Add(str1[--aux.ai]);
                result2.Add('-');
            }

            while (aux.aj > 0)
            {
                result2.Add(str1[--aux.aj]);
                result1.Add('-');
            }

            alignRes.seq1 = ReverseString(result1);
            alignRes.seq2 = ReverseString(result2);
			
            return alignRes;

        }
		string ReverseString(List <char>item)
		{
            char[] aux = new char[item.Count];
			string newString="";
			
			for(int i=item.Count-1;i>=0;i--)
				aux[item.Count-i-1]+=item[i];

            newString = new string(aux);

			return newString;
		}
        public alignSeq Align(string strReference, string str2)
        {
            size1 = strReference.Length + 1;
            size2 = str2.Length + 1;
            //scoreTab = new cell[strReference.Length+1, str2.Length+1];			
			//scoreTab1 = new cell[strReference.Length+1, str2.Length+1];
			//scoreTab2=new cell[strReference.Length+1, str2.Length+1];

            for (int i = 0; i < size1; i++)
            {
                //aux = new cell();
                scoreTab[i, 0].value = i * (-gapPenalty);
                scoreTab[i, 0].i = 0;
				if(i>0)
                    scoreTab[i, 0].i = i - 1;
                scoreTab[i, 0].j = 0;
                scoreTab[i, 0].ai = i; scoreTab[i, 0].aj = 0;
                scoreTab[i, 0].tab = scoreTab;
				//scoreTab[i, 0] = aux;  				
				//aux=new cell();
                scoreTab1[i, 0].value = -gapPenalty - i * gapPenalty / 2;
                scoreTab1[i, 0].i = 0;
				if(i>0)
                    scoreTab1[i, 0].i = i - 1;
                scoreTab1[i, 0].j = 0;
                scoreTab1[i, 0].ai = i; scoreTab1[i, 0].aj = 0;
                scoreTab1[i, 0].tab = scoreTab1;								
				//scoreTab1[i,0]=aux;				
				//aux=new cell();
                scoreTab2[i, 0].value = -10000000;
                scoreTab2[i, 0].i = 0;
				if(i>0)
                    scoreTab2[i, 0].i = i - 1;
                scoreTab2[i, 0].j = 0;
                scoreTab2[i, 0].ai = i; scoreTab2[i, 0].aj = 0;
                scoreTab2[i, 0].tab = scoreTab2;				
				//scoreTab2[i,0]=aux;
				
            }
            for (int i = 0; i < size2; i++)
            {
                //aux = new cell();
                scoreTab[0, i].value = i * (-gapPenalty);
                scoreTab[0, i].j = 0;
				if(i>0)
                    scoreTab[0, i].j = i - 1;
                scoreTab[0, i].i = 0;
                scoreTab[0, i].ai = 0; scoreTab[0, i].aj = i;
                //scoreTab[0,i] = aux;
                //aux = new cell();
                scoreTab1[0, i].value = -10000000;
                scoreTab1[0, i].j = 0;
				if(i>0)
                    scoreTab1[0, i].j = i - 1;
                scoreTab1[0, i].i = 0;
                scoreTab1[0, i].ai = 0; scoreTab1[0, i].aj = i;
                scoreTab1[0, i].tab = scoreTab1;								
				//scoreTab1[0,i]=aux;
                //aux = new cell();
                scoreTab2[0, i].value = -gapPenalty - i * gapPenalty / 2;
                scoreTab2[0, i].j = 0;
				if(i>0)
                    scoreTab2[0, i].j = i - 1;
                scoreTab2[0, i].i = 0;
                scoreTab2[0, i].ai = 0; scoreTab2[0, i].aj = i;
                scoreTab2[0, i].tab = scoreTab2;												
				//scoreTab2[0,i]=aux;
				
            }
            int v1, v2, v3;
			int cost;

            for(int i=1;i<size1;i++)
			{
				char ss=strReference[i-1];
                for (int j = 1; j < size2; j++)
                {
                    //aux = new cell();
                    scoreTab[i, j].ai = i; scoreTab[i, j].aj = j;
					
					if(ss==str2[j-1] || ss=='X' || str2[j-1]=='X')						
						cost=10;
					else
						cost=1;
				
					v1 = scoreTab[i - 1, j - 1].value + cost;
					v2 = scoreTab1[i-1, j-1].value +cost;				
					v3 = scoreTab2[i-1,j-1].value+cost;
                    scoreTab[i, j].i = i - 1;
                    scoreTab[i, j].j = j - 1;
                    scoreTab[i, j].direction = 2;				
                    if (v1 >= v2 && v1>=v3)
                    {
                        scoreTab[i, j].value = v1;
                        scoreTab[i, j].tab = scoreTab;
    
                    }
					else
						if(v2 >= v1 && v2>=v3)
                		{
                            scoreTab[i, j].value = v2;
                            scoreTab[i, j].tab = scoreTab1;
                    	}
						else
						{
                            scoreTab[i, j].value = v3;
                            scoreTab[i, j].tab = scoreTab2;
						}
//					scoreTab[i, j] = aux;
					
                    //aux = new cell();
                    scoreTab1[i, j].ai = i; scoreTab1[i, j].aj = j;

					v1=scoreTab[i-1,j].value-gapPenalty;
					v2=scoreTab1[i-1,j].value-gapPenalty/2;
					v3=scoreTab2[i-1,j].value-gapPenalty;
                    scoreTab1[i, j].i = i - 1;
                    scoreTab1[i, j].j = j;
                    scoreTab1[i, j].direction = 0;					
					if(v1>=v2 && v1>=v3)
					{
                        scoreTab1[i, j].value = v1;
                        scoreTab1[i, j].tab = scoreTab;
					}
					else
						if(v2>=v1 && v2>=v3)
						{
                            scoreTab1[i, j].value = v2;
                            scoreTab1[i, j].tab = scoreTab1;
						}
						else
						{
                            scoreTab1[i, j].value = v3;
                            scoreTab1[i, j].tab = scoreTab2;
						}
					
					//scoreTab1[i,j]=aux;
						
                    //aux = new cell();
                    scoreTab2[i, j].ai = i; scoreTab2[i, j].aj = j;

					v1=scoreTab[i,j-1].value-gapPenalty;
					v2=scoreTab1[i,j-1].value-gapPenalty;
					v3=scoreTab2[i,j-1].value-gapPenalty/2;
                    scoreTab2[i, j].i = i;
                    scoreTab2[i, j].j = j - 1;
                    scoreTab2[i, j].direction = 1;					
					if(v1>=v2 && v1>=v3)
					{
                        scoreTab2[i, j].value = v1;
                        scoreTab2[i, j].tab = scoreTab;
					}
					else
						if(v2>=v1 && v2>=v3)
						{
                            scoreTab2[i, j].value = v2;
                            scoreTab2[i, j].tab = scoreTab1;
						}
						else
						{
                            scoreTab2[i, j].value = v3;
                            scoreTab2[i, j].tab = scoreTab2;
						}
					
					//scoreTab2[i,j]=aux;										                    
                }
			}
//			PrintScore ();
            return TraceBack(strReference, str2);
        }
		void PrintScore()
		{
            for(int i=0;i<size1;i++)
			{
                for (int j = 0; j < size2; j++)                
				Console.Write(scoreTab[i,j].value+" ");				
				
				Console.WriteLine();
			}			
		}
        void ReadCostMatrix(string fileName)
        {
            StreamReader st = new StreamReader(fileName);
            Dictionary <string,int> dic;
            String line;
            while ((line = st.ReadLine()) != null)
            {
                string[] strA = line.Split(' ');
                if (costMatrix.ContainsKey(strA[0]))
                    costMatrix[strA[0]].Add(strA[1], Convert.ToInt16(strA[2]));
                else
                {
                    dic = new Dictionary<string, int>();
                    dic.Add(strA[1], Convert.ToInt16(strA[2]));
                    costMatrix.Add(strA[0], dic);
                }                
            }

            st.Close();
        }
    }
    class cell
    {
        public int i, j;
		public int ai,aj;
        public int value;
        public int direction;
		public cell  [,] tab;
		         
    }
}
