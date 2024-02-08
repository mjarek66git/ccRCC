using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.IO;
using System.Windows.Forms;
using System.Drawing.Drawing2D;
using phiClustCore;
using phiClustCore.Interface;

namespace Graph
{
    struct ColorStr
    {
        public int[] rgb;
        public int size;
        public int colorPosition;

        public ColorStr(int size)
        {
            rgb = new int [3];
            this.size = size;
            colorPosition = 0;
        }

        public void IncrementSize()
        {
            size++;
        }
        public void SetColor(int  []r)
        {
            for (int i = 0; i <3;i++ )
                rgb[i] = r[i];
        }
    }
    public partial class VisHierarCircle : Form, SavePic,IVisual
    {
        HClusterNode hnode;
        HClusterNode currentNode = null;
        char labelBreak = ':';
        public ClosingForm closeForm;
        Dictionary<string, GraphicsPath> dPaths = new Dictionary<string, GraphicsPath>();
        List <int []> colorMap=new List<int[]>();
        Stack<HClusterNode> memoryClicks = new Stack<HClusterNode>();
        Dictionary<Region,HClusterNode> allRegions = new Dictionary<Region,HClusterNode>();
        Dictionary<Region, string> colorRegions = new Dictionary<Region, string>();
        Dictionary<string, string> labels = new Dictionary<string, string>();
        Dictionary<string, int> labelColor = new Dictionary<string, int>();

        Bitmap buffer;
        double distStep=1;
        int stepBegin = 20;
       
        double angleStep;
        string winName;      
        int all;
        //double x1, y1, x2, y2,x3,y3,x4,y4;
        double maxDist;
        Random r = new Random();
        Point middleScreen;
        string measureName;
        ClusterOutput outp;
        public VisHierarCircle(ClusterOutput outp, string name, string measureName)
        {
            this.outp = outp;

            if (outp.hNode == null)
            {
                MessageBox.Show("There is no hierarchical clustering");
                return;
            }
            this.measureName=measureName;
            this.hnode = outp.hNode;
            currentNode = hnode;
            winName = name;
            this.Text = name;
            InitializeComponent();
            middleScreen = new Point(panel1.Size.Width / 2, panel1.Size.Height / 2);
            StartVis(hnode,panel1.Size.Height);
            int []tab=new int [4];

            tab[0] = 0; tab[1] = 85; tab[2] = 170; tab[3] = 255;

            for(int i=0;i<tab.Length;i++)
                for(int j=0;j<tab.Length;j++)
                    for (int n = 0; n < tab.Length; n++)
                    {
                        int[] aux = new int[3];
                        aux[0] = tab[i];
                        aux[1] = tab[j];
                        aux[2] = tab[n];

                        colorMap.Add(aux);
                    }

            FindAllLabels();
            if (labels.Count == 0)
            {
                string labelFileName = hnode.dirName+"_label.dat";
                if (File.Exists(labelFileName))
                    ReadLabels(labelFileName);
                else
                    ReadLabels();
            }
        }
        public override string ToString()
        {
            return "Circle";
        }
        public void ToFront()
        {
            this.BringToFront();
        }
        public void SavePicture(string fileName,Bitmap bmp)
        {

            if (!fileName.Contains(".png"))
                fileName += ".png";
            bmp.Save(fileName, System.Drawing.Imaging.ImageFormat.Png);           
        }
        void DefineColorsForLabels()
        {
            if (labels == null || labels.Count==0)
                return;

            foreach(var item in labels.Values)
                if(!labelColor.ContainsKey(item))                
                    labelColor.Add(item, 0);

            int step = (colorMap.Count-1) / labelColor.Count;
            int counter = 0;
            List<string> cList = new List<string>(labelColor.Keys);
            for (int i = 0; i < cList.Count;i++ )
            {

                labelColor[cList[i]] = step * counter + 1;

                counter++;
            }
            
        }
        void FindAllLabels()
        {
            labels.Clear();
            int counter=0;
            Dictionary<string, int> colorDic = new Dictionary<string, int>();
            foreach (var item in hnode.setStruct)
            {
                if (item.Contains(labelBreak))
                {
                    string[] aux = item.Split(labelBreak);
                    if(!colorDic.ContainsKey(aux[1]))
                        colorDic.Add(aux[1],counter++);
                }
            }

            if (counter > colorMap.Count)
                throw new Exception("Number of labels: " + counter + " is to big. The limit is: " + colorMap.Count);

            foreach (var item in hnode.setStruct)
            {
                if (item.Contains(labelBreak))
                {
                    string[] aux = item.Split(labelBreak);
                    if (!labels.ContainsKey(aux[0]))
                            labels.Add(aux[0],aux[1]);
                }
            }
            DefineColorsForLabels();
                
        }
        void StartVis(HClusterNode hN,int Height)
        {
            all = hN.setStruct.Count;

            maxDist = hN.realDist;
            angleStep = 360.0 / all;
            distStep = (Height / 2 - stepBegin) / maxDist;
        }
        void DrawPaintClusters(Graphics gr,HClusterNode hN,double startAngle)
        {
            if (hN == null || hN.joined==null)
                return;


            if (IsClusterClean(hN))
            {
                Dictionary<string, ColorStr> dicLab = GetSizeofLabels(hN);
                List<ColorStr> colorOrder = new List<ColorStr>();
                foreach (var it in dicLab)
                    colorOrder.Add(it.Value);
                double sweepAngle = angleStep * hN.setStruct.Count;
                DrawCluster(hN, colorOrder, stepBegin, gr, 2*stepBegin, startAngle, sweepAngle);
                return;
            }
            foreach (var item in hN.joined)
            {
              
                double sweepAngle = angleStep*item.setStruct.Count;
                Dictionary<string, ColorStr> dicLab = GetSizeofLabels(item);
                List<string> properOrder = GetLabelOrder(item);

                List<ColorStr> colorOrder = new List<ColorStr>();
                foreach (var it in properOrder)
                    colorOrder.Add(dicLab[it]);

                if(Math.Abs(item.realDist-hnode.realDist)>0.01 )
                    DrawCluster(item,colorOrder,(maxDist - item.realDist) * distStep + stepBegin, gr, (maxDist - hN.realDist) * distStep + stepBegin, startAngle, sweepAngle);
                if(!IsClusterClean(item))
                    DrawPaintClusters(gr, item, startAngle);
                startAngle += sweepAngle;
                //break;
            }


          
        }
        void DrawFrameClusters(Graphics gr, HClusterNode hN, double startAngle)
        {
            if (hN == null || hN.joined == null)
                return;

            if (IsClusterClean(hN))
            {
                Dictionary<string, ColorStr> dicLab = GetSizeofLabels(hN);
                double sweepAngle = angleStep * hN.setStruct.Count;
                DrawFrame(hN, stepBegin, gr, 2 * stepBegin, startAngle, sweepAngle);
                return;
            }

            foreach (var item in hN.joined)
            {

                double sweepAngle = angleStep * item.setStruct.Count;
                if (Math.Abs(item.realDist - hnode.realDist) > 0.01)
                    DrawFrame(item, (maxDist - item.realDist) * distStep + stepBegin, gr, (maxDist - hN.realDist) * distStep + stepBegin, startAngle, sweepAngle);
                if (!IsClusterClean(item))
                    DrawFrameClusters(gr, item, startAngle);
                startAngle += sweepAngle;
                //break;
            }



        }
        public void DrawFrame(HClusterNode node, double distEnd, Graphics gr, double pDist, double angleStart, double sweepAngle)
        {
            GraphicsPath gPath;

            gPath = PreparePath(distEnd, pDist, angleStart, sweepAngle);
            if (gPath != null)
            {
                Pen p = new Pen(Color.Black, 2);
                gr.DrawPath(p, gPath);
                allRegions.Add(new Region(gPath), node);
            }

        }
        public bool IsClusterClean(HClusterNode node)
        {
            Dictionary<string, int> freq = new Dictionary<string, int>();

            foreach (var item in node.setStruct)
            {
                string colorC;
                if (labels.ContainsKey(item))
                {
                    colorC = labels[item];
                    if (!freq.ContainsKey(colorC))
                        freq.Add(colorC, 0);
                    freq[colorC]++;
                }
            }

            if (freq.Keys.Count <= 1)
                return true;

            return false;
        }
        void ChangeColorMap(Dictionary<string,int[]> dicLabels)
        {

        }
        List <string> GetLabelOrder(HClusterNode hNode)
        {
            List<string> orderList = new List<string>();

            if (hNode.joined != null)
            {
                foreach (var item in hNode.joined)
                {
                    List<string> aux = GetLabelOrder(item);
                    foreach (var auxItem in aux)
                        if (!orderList.Contains(auxItem))
                            orderList.Add(auxItem);
                }
            }
            else
            {
                foreach (var item in hNode.setStruct)
                    if (!orderList.Contains(labels[item]))
                        orderList.Add(labels[item]);

                orderList.Sort((x, y) => x.CompareTo(y));
            }
            return orderList;
        }
        Dictionary<string,ColorStr> GetSizeofLabels(HClusterNode hNode)
        {
            Dictionary<string, ColorStr> dic = new Dictionary<string, ColorStr>();
            foreach (var item in hNode.setStruct)
            {
               /* if (labels == null)
                {
                    
                    string[] aux = item.Split(labelBreak);
                    if (!dic.ContainsKey(aux[0]))
                    {
                        ColorStr str = new ColorStr(1);
                        dic.Add(aux[0], str);
                    }
                    else
                    {
                        ColorStr copyStr = dic[aux[0]];
                        copyStr.IncrementSize();
                        dic[aux[0]]=copyStr;
                    }

                }
                else
                {*/
                    if (labels.ContainsKey(item))
                    {
                        string colorToString = labels[item];
                        if (!dic.ContainsKey(colorToString))
                        {
                            ColorStr str = new ColorStr(1);
                            dic.Add(colorToString, str);
                        }

                        else
                        {
                            ColorStr copyStr = dic[colorToString];
                            copyStr.IncrementSize();
                            dic[colorToString]=copyStr;
                        }
                    }
                //}
            }
            if (dic.Keys.Count > colorMap.Count)
                throw new Exception("Too many labels");
            
            foreach (var item in dic)
                item.Value.SetColor(colorMap[labelColor[item.Key]]);
            return dic;
        }
        GraphicsPath PreparePath(double distEnd,double pDist, double angleStart, double sweepAngle)
        {
                GraphicsPath gPath = new GraphicsPath();
                GraphicsPath auxPath1 = new GraphicsPath();
                GraphicsPath auxPath2 = new GraphicsPath();
                int pDistInt = (int)(pDist);
                int distEndInt = (int)(distEnd);
                int angleStartInt = (int)Math.Round(angleStart);
                int sweepAngleInt = (int)Math.Round(sweepAngle);
              //  if (sweepAngleInt >= 1)
                {
                    int midX=middleScreen.X - pDistInt;
                    int midY = middleScreen.Y - pDistInt;
                    int midXE = middleScreen.X - distEndInt;
                    int midYE = middleScreen.Y - distEndInt;
                    auxPath1.AddArc(midX, midY, (int)(2 * pDist), (int)(2 * pDist), (float)angleStart, (float)sweepAngle);
                    //auxPath2.AddArc(midXE, midYE, (int)(2 * distEnd), (int)(2 * distEnd), (int)Math.Round(angleStart + sweepAngle), -sweepAngleInt);
                    auxPath2.AddArc(midXE, midYE, (int)(2 * distEnd), (int)(2 * distEnd), (float)(angleStart + sweepAngle), -(float)sweepAngle);
                    if (auxPath1.PathPoints == null || auxPath2.PathPoints == null||auxPath1.PathPoints.Length <= 1 || auxPath2.PathPoints.Length <= 1)
                        return null;

                    gPath.AddArc(midX, midY, (int)(2 * pDist), (int)(2 * pDist), (float)angleStart, (float)sweepAngle);
                    gPath.AddLine(auxPath1.PathPoints[auxPath1.PathPoints.Length - 1].X, auxPath1.PathPoints[auxPath1.PathPoints.Length - 1].Y, auxPath2.PathPoints[0].X, auxPath2.PathPoints[0].Y);
                    //gPath.AddArc(midXE, midYE, (int)(2 * distEnd), (int)(2 * distEnd), (int)Math.Round(angleStart + sweepAngle), -sweepAngleInt);
                    gPath.AddArc(midXE, midYE, (int)(2 * distEnd), (int)(2 * distEnd),(float)( angleStart + sweepAngle), -(float)sweepAngle);
                    gPath.AddLine(auxPath2.PathPoints[auxPath2.PathPoints.Length - 1].X, auxPath2.PathPoints[auxPath2.PathPoints.Length - 1].Y, auxPath1.PathPoints[0].X, auxPath1.PathPoints[0].Y);
                    
                }
              // else
              //      return null;

            return gPath;

        }
        void DrawCluster(HClusterNode node, List<ColorStr> colorList,double distEnd,Graphics gr,double pDist,double angleStart,double sweepAngle)
        {

            double localAngleStart = angleStart;
            GraphicsPath gPath;

            int counter = 0;
            foreach (var item in colorList)
            {
                double localSweepAngle = item.size * angleStep;
                counter++;

              //  if (counter ==3)//|| counter==1)
                {
                    Color currentColor;

                    gPath = PreparePath(distEnd, pDist, localAngleStart, localSweepAngle);

                    if (gPath == null)
                        continue;
                    currentColor = Color.FromArgb(item.rgb[0], item.rgb[1], item.rgb[2]);
                    SolidBrush myBrush = new SolidBrush(currentColor);

                    gr.FillPath(myBrush, gPath);                        
                    
                    //localAngleStart += localSweepAngle;
                }
                localAngleStart += localSweepAngle;
                
            }
/*            gPath = PreparePath(distEnd, pDist, angleStart, sweepAngle);

            if (gPath != null)
            {
                Pen p = new Pen(Color.Black, 1);
                gr.DrawPath(p, gPath);
                allRegions.Add(new Region(gPath),node);
            }*/
        }
        public void DrawOnBuffer(Bitmap buff, bool drawLegend=true)
        {
            if (currentNode == null)
                return;
            middleScreen = new Point(buff.Width / 2, buff.Height / 2);
            Font drawFont = new System.Drawing.Font("Arial", 8);
            SolidBrush drawBrush = new System.Drawing.SolidBrush(System.Drawing.Color.Black);
            float frac = (float)(currentNode.setStruct.Count) / hnode.setStruct.Count;
            frac *= 100;
            string ss = frac.ToString("F1");
            Dictionary<string, ColorStr> dicLabels = GetSizeofLabels(currentNode);         
            Graphics g = Graphics.FromImage(buff);
            //Graphics g = e.Graphics;

            
            

            int x = 10, y = 25;
            colorRegions.Clear();
            if(drawLegend)
                foreach (var item in dicLabels)
                {
                    drawBrush = new System.Drawing.SolidBrush(Color.FromArgb(item.Value.rgb[0], item.Value.rgb[1], item.Value.rgb[2]));
                    g.FillRectangle(drawBrush, x, y, 15, 10);
                    drawBrush = new System.Drawing.SolidBrush(System.Drawing.Color.Black);
                    g.DrawString(item.Key, drawFont, drawBrush, x + 20, y);

                    Rectangle rect = new Rectangle(x, y, 15, 10);
                    colorRegions.Add(new Region(rect), item.Key);

                    y += 25;
                    if (y > buff.Height)
                    {
                        x += 150;
                        y = 25;
                    }

                }
            
            StartVis(currentNode,buff.Height);

            drawFont = new System.Drawing.Font("Arial", 11);
            g.DrawString(ss, drawFont, drawBrush, middleScreen.X - 20, middleScreen.Y - 10);
            DrawPaintClusters(g, currentNode, 0);
            DrawFrameClusters(g, currentNode, 0);
           
        }

        private void panel1_Paint(object sender, PaintEventArgs e)
        {

            /*if (currentNode == null)
                return;
            Font drawFont = new System.Drawing.Font("Arial", 8);
            SolidBrush drawBrush = new System.Drawing.SolidBrush(System.Drawing.Color.Black);
            float frac=(float)(currentNode.setStruct.Count)/hnode.setStruct.Count;
            frac *= 100;
            string ss=frac.ToString("F1");
            Dictionary <string,ColorStr> dicLabels=GetSizeofLabels(currentNode);
            buffer = new Bitmap(panel1.Width, panel1.Height);
            Graphics g = Graphics.FromImage(buffer);
            //Graphics g = e.Graphics;

            e.Graphics.Clear(this.BackColor);

            int x=10,y=25;
            colorRegions.Clear();
            foreach (var item in dicLabels)
            {
                drawBrush = new System.Drawing.SolidBrush(Color.FromArgb(item.Value.rgb[0], item.Value.rgb[1], item.Value.rgb[2]));
                g.FillRectangle(drawBrush, x, y, 15, 10);
                drawBrush = new System.Drawing.SolidBrush(System.Drawing.Color.Black);
                g.DrawString(item.Key, drawFont, drawBrush, x + 20, y);

                Rectangle rect = new Rectangle(x, y, 15, 10);
                colorRegions.Add(new Region(rect), item.Key);

                y += 25;
                if (y > this.Size.Height)
                {
                    x += 150;
                    y = 25;
                }
                
            }

            StartVis(currentNode,buffer.Height);
           
            drawFont = new System.Drawing.Font("Arial", 11);
            g.DrawString(ss, drawFont, drawBrush, middleScreen.X - 20, middleScreen.Y-10);
            DrawPaintClusters(g, currentNode, 0);
            DrawFrameClusters(g, currentNode, 0);*/
            buffer = new Bitmap(panel1.Width, panel1.Height);
            DrawOnBuffer(buffer);
            e.Graphics.DrawImage(buffer, 0, 0);
        }

        private void panel1_Resize(object sender, EventArgs e)
        {
            distStep = (this.Size.Height / 2 - stepBegin) / maxDist;
            middleScreen = new Point(this.Size.Width / 2, this.Size.Height / 2);
            panel1.Refresh();
            this.Invalidate();
        }

        private void VisHierarCircle_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (closeForm != null)
                closeForm(winName);
        }
        HClusterNode CheckRegion(int mouseX,int mouseY)
        {
            foreach (var item in allRegions)
            {
                if (item.Key.IsVisible(mouseX, mouseY))
                    return item.Value;
            }

            return null;
        }
        string CheckColorRegion(int mouseX, int mouseY)
        {
            foreach (var item in colorRegions)
            {
                if (item.Key.IsVisible(mouseX, mouseY))
                    return item.Value;
            }

            return null;
        }
        private void panel1_MouseClick(object sender, MouseEventArgs e)
        {
//            thr += 2;
//            CC = 0;
            if (e.Button == MouseButtons.Left)
            {
                HClusterNode remNode = currentNode;
                currentNode = CheckRegion(e.X, e.Y);
                if (currentNode != null)
                {
                    allRegions.Clear();
                    memoryClicks.Push(remNode);
                    this.Invalidate();
                    this.Refresh();
                }
                else
                    currentNode = remNode;

                string lab=CheckColorRegion(e.X,e.Y);
                if (lab != null)
                {
                    DialogResult res;
                    colorDialog1 = new ColorDialog();

                    res=colorDialog1.ShowDialog();
                    if(res==DialogResult.OK)
                    {
                        colorMap[labelColor[lab]][0] = colorDialog1.Color.R;
                        colorMap[labelColor[lab]][1] = colorDialog1.Color.G;
                        colorMap[labelColor[lab]][2] = colorDialog1.Color.B;

                        this.Invalidate();
                        this.Refresh();
                    }
                }
                toolStripButton2.Enabled = true;
            }
            else
            {
                HClusterNode remNode = CheckRegion(e.X, e.Y);

                if (remNode != null)
                {
                    visHierar dendrog = new visHierar(remNode, outp,"", measureName, null);
                    dendrog.Show();
                }
            }
        }
        private void ReadLabels(string fileName=null)
        {
            DialogResult res=DialogResult.OK;
            if (fileName == null)
            {
                openFileDialog1 = new OpenFileDialog();
                openFileDialog1.Title = "File name with the name of labels for each struture";
                res = openFileDialog1.ShowDialog();
                fileName = openFileDialog1.FileName;
            }
            if (res == DialogResult.OK)
            {
                labels = new Dictionary<string, string>();
                Dictionary<string, int> colorDic = new Dictionary<string, int>();
                StreamReader st = new StreamReader(fileName);
                string line;
                while ((line = st.ReadLine()) != null)
                {
                    string[] aux = line.Split(' ');

                    if (aux.Length == 2)
                    {
                        if(aux[0].Contains(Path.DirectorySeparatorChar))
                            aux[0]=Path.GetFileName(aux[0]);

                        if (hnode.setStruct.Contains(aux[0]))                        
                            labels.Add(aux[0], aux[1]);
                           
                        
                    }

                }

                st.Close();
              
                DefineColorsForLabels();

            }
        }
        private void toolStripButton1_Click(object sender, EventArgs e)
        {
            ReadLabels();
        }

        private void toolStripButton2_Click(object sender, EventArgs e)
        {
            if (memoryClicks.Count != 0)
            {
                currentNode = memoryClicks.Pop();
                allRegions.Clear();
                this.Invalidate();
                this.Refresh();
            }
            else
                toolStripButton2.Enabled = false;
        }

        private void toolStripButton3_Click(object sender, EventArgs e)
        {
            ClusterOutput output = new ClusterOutput();
            output.hNode = hnode;
            ClusterVis wrCluster = new ClusterVis(output);
            DialogResult res = saveFileDialog1.ShowDialog();
            if (res == DialogResult.OK && saveFileDialog1.FileName.Length > 0)
            {
                Resolution resForm = new Resolution(buffer.Width, buffer.Height, Color.Black);
                res = resForm.ShowDialog();
                if (res == DialogResult.OK)
                {
                    Bitmap bmp = new Bitmap(resForm.WidthR, resForm.HeightR);
                    DrawOnBuffer(bmp, false);
                    //DrawOnBuffer(bmp, resForm.ShowLegend, resForm.LineThickness, resForm.LinesColor);
                    SavePicture(saveFileDialog1.FileName, bmp);
                    //PrepareGraphNodes(buffer);
                }
                //this.SavePicture(saveFileDialog1.FileName, buffer);
            }
        }
    }
}
