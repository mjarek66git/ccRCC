using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.IO;
using phiClustCore;
using phiClustCore.Interface;

namespace Graph
{
    interface SavePic
    {
        void SavePicture(string fileName,Bitmap bmp);
    }
    public partial class visHierar : Form,IVisual,SavePic
    {        
        public Dictionary< HClusterNode,Color> listNodes = null;
        List<HClusterNode> sizeNodes = new List<HClusterNode>();
//        Dictionary<string, Color> classColor = null;
        Dictionary<string, string> vecColor = null;
        List <Color> colorMap=new List<Color>();
        Dictionary<double,Color> profilesColorMap = null;
        ClusterOutput outp;
        Bitmap buffer;
        bool linemode = false;
        bool markflag = false;
        bool viewType = true;
        public ClosingForm closeForm;
        Dictionary<string, List<double>> profiles;
        List<string> geneLabels = new List<string>();
        List<string> sampleLabels = new List<string>();

        List<Tuple<string, string>> labels = new List<Tuple<string, string>>();

        List<int> upReg = null;
        List<int> downReg = null;
        double threshold = 0.7;

        bool clearAll = false;
        string winName;
        int mposX, mposY;
        DrawHierarchical drawH;
        HClusterNode hnode;
        public override string ToString()
        {
            return "Dendrogram";
        }
        public void ToFront()
        {
            this.BringToFront();
        }
        visHierar(string name, string measureName, Dictionary<string, string> labels)
        {
            InitializeComponent();            
            buffer = new Bitmap(pictureBox1.Width, pictureBox1.Height);
            
            winName = name;
            this.Text = name;
            vecColor = labels;            

        }
        public visHierar(ClusterOutput outp, string name, string measureName, Dictionary<string, string> labels):this(name,measureName,labels)
        {
            this.outp = outp;
            hnode = outp.hNode;
            drawH = new DrawHierarchical(hnode, measureName, labels, buffer, true);
            vScrollBar1.Maximum = hnode.setStruct.Count;
            drawH.horizontalView = true;
            profiles = outp.profiles;
            CreateProfilesColorMap();
            InitVisHier();
            geneLabels = outp.data.geneLabels;
            sampleLabels = outp.data.sampleLabels;
        }
       
        public visHierar(HClusterNode hnode,ClusterOutput outp,string name,string measureName,Dictionary<string,string> labels):this(outp,name,measureName,labels)
        {
            this.hnode = hnode;
        }
        void InitVisHier()
        {
            int[] tab = new int[4];

            tab[0] = 0; tab[1] = 85; tab[2] = 170; tab[3] = 255;

            for (int i = 0; i < tab.Length; i++)
                for (int j = 0; j < tab.Length; j++)
                    for (int n = 0; n < tab.Length; n++)
                        colorMap.Add(Color.FromArgb(tab[i], tab[j], tab[n]));


            drawH.SetColors();
            

        }
        private void ClearBuffer()
        {
            Graphics g = Graphics.FromImage(buffer);
            g.Clear(pictureBox1.BackColor);
        }
        private void Form2_MouseClick(object sender, MouseEventArgs e)
        {
        
        }
        private void Form2_ResizeEnd(object sender, EventArgs e)
        {
            clearAll = true;
            if (pictureBox1.Width == 0 || pictureBox1.Height == 0)
                return;
            buffer = new Bitmap(pictureBox1.Width, pictureBox1.Height);
            drawH.PrepareGraphNodes(buffer);
            drawH.DrawOnBuffer(buffer,true,1,Color.Empty);
            pictureBox1.Refresh();
            this.Invalidate();
            drawH.maxGraphicsY = pictureBox1.Height-drawH.posStart-30;
            drawH.maxGraphicsX = pictureBox1.Width - drawH.posStart - 30;
        }
        public void SaveToFile(string fileName,bool ShowLegend,int LineThickness,Color lineColor,int resWidth=800,int resHeight=600)
        {
                Bitmap bmp = new Bitmap(resWidth, resHeight);
                drawH.PrepareGraphNodes(bmp);

                drawH.DrawOnBuffer(bmp,ShowLegend, LineThickness, lineColor);
                SavePicture(fileName, bmp);
                drawH.PrepareGraphNodes(buffer);

        }
        private void toolStripButton1_Click(object sender, EventArgs e)
        {
            ClusterOutput output = new ClusterOutput();
            output.hNode = drawH.hnode;
            ClusterVis wrCluster=new ClusterVis(output);
            saveFileDialog1.DefaultExt = "png";
            saveFileDialog1.Filter = "Png files|*.png";
            DialogResult res=saveFileDialog1.ShowDialog();
            if (res == DialogResult.OK && saveFileDialog1.FileName.Length > 0)
            {
                Resolution resForm = new Resolution(buffer.Width, buffer.Height,Color.Black);
                res=resForm.ShowDialog();
                if (res == DialogResult.OK)
                {
                    Bitmap bmp = new Bitmap(resForm.WidthR, resForm.HeightR);
                    drawH.PrepareGraphNodes(bmp);

                    drawH.DrawOnBuffer(bmp,resForm.ShowLegend,resForm.LineThickness,resForm.LinesColor);
                    SavePicture(saveFileDialog1.FileName, bmp);
                    drawH.PrepareGraphNodes(buffer);
                }
                //       this.SavePicture(saveFileDialog1.FileName);          
            }
        }

        private void toolStripButton2_Click(object sender, EventArgs e)
        {
            if (drawH.hnode != drawH.rootNode)
            {
                drawH.rootNode = drawH.hnode;
                drawH.PrepareGraphNodes(buffer);
                ClearBuffer();
                pictureBox1.Invalidate();
            }
        }

        private void Form2_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (closeForm != null)
                closeForm(winName);
        }
        public void SavePicture(string fileName,Bitmap bmp)
        {
            if (!fileName.Contains(".png"))
                fileName += ".png";
            bmp.Save(fileName, ImageFormat.Png);
        }
        private void pictureBox1_Paint(object sender, PaintEventArgs e)
        {

            if (buffer != null && drawH.currentRootNode == drawH.rootNode && buffer.Height == pictureBox1.Height && buffer.Width == pictureBox1.Width)
            {
                e.Graphics.Clear(this.BackColor);
                e.Graphics.DrawImage(buffer, 0, 0);
                if (linemode)
                {

                    Pen p = new Pen(Color.Brown);
                    Font drawFont = new Font("Arial", 10);
                    SolidBrush drawBrush = new SolidBrush(Color.Black);
                    StringFormat drawFormat = new StringFormat();


                    if (listNodes != null)
                        e.Graphics.DrawString(listNodes.Count.ToString(), drawFont, drawBrush, 50, mposY);
                    e.Graphics.DrawLine(p, 0, mposY, buffer.Width, mposY);
                    return;
                }
                if (listNodes != null)
                {
                    SolidBrush brush;
                    foreach (var item in listNodes)
                    {
                        if (item.Key.visible)
                        {
                            brush = new SolidBrush(item.Value);
                            if (drawH.rootNode.IsVisible(item.Key))
                                e.Graphics.FillEllipse(brush, item.Key.gNode.x, item.Key.gNode.y - 7, 7, 7);
                        }
                    }
                }
                if(sizeNodes.Count>0)
                {
                    SolidBrush brush;
                    int maxV = 0;
                    foreach (var item in sizeNodes)
                        if (item.setStruct.Count > maxV)
                            maxV = item.setStruct.Count;
                    foreach (var item in sizeNodes)
                    {
                       
                        brush = new SolidBrush(Color.Blue);
                        if (drawH.rootNode.IsVisible(item))
                        {
                            double res =( (double)item.setStruct.Count) / maxV;
                            res = 360 * res;
                            e.Graphics.FillPie(brush, item.gNode.x-25, item.gNode.y - 7, 15, 15,0,(float)res);
                        }
                    }

                }

            }
            else
            {
              
                if (drawH.hnode != null)
                {
                    if (clearAll)
                    {
                        //pictureBox1.Image = null;
                    }
/*                    e.Graphics.SmoothingMode = SmoothingMode.AntiAlias;
                    e.Graphics.TextRenderingHint = System.Drawing.Text.TextRenderingHint.AntiAlias;*/
                    if (buffer == null || buffer.Width != pictureBox1.Width || buffer.Height != pictureBox1.Height)
                    {
                        buffer = new Bitmap(pictureBox1.Width, pictureBox1.Height);
                    }
                    drawH.currentRootNode = drawH.rootNode;

                    drawH.DrawOnBuffer(buffer,true,1,Color.Empty);
                    e.Graphics.DrawImage(buffer, 0, 0);
                   // buffer.Save("proba.png", ImageFormat.Png);
                    //maxHDist = minHDist = rootNode;
                    //DrawGraph(rootNode, e.Graphics);
                    //DrawDistanceAx(rootNode, e.Graphics);

                    if (listNodes != null)
                    {
                        SolidBrush brush;
                        foreach (var item in listNodes)
                        {
                            brush= new SolidBrush(item.Value);
                            if (drawH.rootNode.IsVisible(item.Key))
                                e.Graphics.FillEllipse(brush, item.Key.gNode.x, item.Key.gNode.y - 7, 7, 7);
                        }
                    }
                    if (sizeNodes.Count > 0)
                    {
                        SolidBrush brush;
                        int maxV = 0;
                        foreach (var item in sizeNodes)
                            if (item.setStruct.Count > maxV)
                                maxV = item.setStruct.Count;
                        foreach (var item in sizeNodes)
                        {

                            brush = new SolidBrush(Color.Blue);
                            if (drawH.rootNode.IsVisible(item))
                            {
                                double res = ((double)item.setStruct.Count) / maxV;
                                res = 360 * res;
                                e.Graphics.FillPie(brush, item.gNode.x-25, item.gNode.y - 7, 15, 15, 0, (float)res);
                            }
                        }

                    }

                }
            }

        }

        private void pictureBox1_MouseClick(object sender, MouseEventArgs e)
        {
            HClusterNode clickNode;

            if (linemode)
            {
                linemode = false;
                if (listNodes.Count > 0)
                {
                    SaveMarkedClusters.Enabled = true;
                    orderVis.Enabled = true;
                }
                pictureBox1.Invalidate();
                return;
            }
            if (markflag)
            {
                clickNode = drawH.CheckClick(drawH.rootNode, e.X, e.Y);
                if (clickNode != null)
                {
                    if (e.Button == MouseButtons.Right)
                    {
                        if (listNodes.ContainsKey(clickNode))
                            listNodes.Remove(clickNode);
                    }
                    else
                    {
                        if (listNodes == null)
                            listNodes = new Dictionary<HClusterNode, Color>();
                        if (!listNodes.ContainsKey(clickNode))
                            listNodes.Add(clickNode,Color.Red);
                    }
                    if (listNodes.Count > 0)
                    {
                        SaveMarkedClusters.Enabled = true;
                        orderVis.Enabled = true;
                    }
                    pictureBox1.Invalidate();
                }


                return;
            }
            clickNode = drawH.CheckClick(drawH.rootNode, e.X, e.Y);

            if (clickNode != null)
            {
                if (e.Button == MouseButtons.Right)
                {
                    drawH.rootNode = clickNode;
                    drawH.PrepareGraphNodes(buffer);
                    ClearBuffer();
                    pictureBox1.Invalidate();

                    //this.Invalidate();
                }
                else
                {                    
                    FormText info = new FormText(clickNode);
                    info.Show();
                    if (labels.Count == 0)
                        foreach (var item in geneLabels)
                            labels.Add(new Tuple<string, string>(item, ""));
                        DrawPanel pn = new DrawPanel("Leave profile: " + clickNode.refStructure,labels);
                        pn.CreateBMP();
                        pn.drawPic = delegate { double res = DrawNodeProfiles(pn.bmp, clickNode.setStruct, clickNode.consistency); pn.Text = res.ToString(); };
                        pn.Show();
                        pn.pictureBox1.Invalidate();
                    
                    
                }
            }
            string lab = drawH.CheckColorRegion(e.X, e.Y);
            if (lab != null)
            {
                DialogResult res;
                colorDialog1 = new ColorDialog();
                int []colorTab=new int [drawH.classColor.Count];
                int i=0;
                foreach (var item in drawH.classColor)
                    colorTab[i++] = ColorTranslator.ToOle(item.Value);
                colorDialog1.CustomColors = colorTab;
   
                colorDialog1.Color = drawH.classColor[lab];
                
                res = colorDialog1.ShowDialog();
                if (res == DialogResult.OK)
                {
                    drawH.classColor[lab] = Color.FromArgb(colorDialog1.Color.R, colorDialog1.Color.G, colorDialog1.Color.B);
                    buffer = null;
                    pictureBox1.Invalidate();
                    pictureBox1.Refresh();
                }
            }

        }
        void CreateProfilesColorMap()
        {
            if (profiles == null)
                return;

            if (outp.data.intervals?.Count > 0)
            {
                List<int> counterM = new List<int>();
                List<int> counterP = new List<int>();
                int[] codes = new int[outp.data.intervals[0].Length];
                profilesColorMap = new Dictionary<double, Color>();
                for (int i = 0;i < outp.data.intervals[0].Length;i++)
                {
                    codes[i] = outp.data.intervals[0][i].code;
                    if (codes[i] < 0)
                        counterM.Add(codes[i]);
                    if (codes[i] == 0)
                        profilesColorMap.Add(0, Color.FromArgb(0, 0, 0));
                    if (codes[i] > 0)
                        counterP.Add(codes[i]);
                }
                
                if(counterM.Count>0)
                {
                    int[] aux = counterM.ToArray();
                    Array.Sort(aux);
                    aux.Reverse();

                    double st = 255.0 / aux.Length;
                    for(int i=0;i<aux.Length;i++)                    
                        profilesColorMap.Add(aux[i], Color.FromArgb(0,0,(int)(st*(i+1))));                    
                }
                if (counterP.Count > 0)
                {
                    int[] aux = counterP.ToArray();
                    Array.Sort(aux);                    

                    double st = 255.0 / aux.Length;
                    for (int i = 0; i < aux.Length; i++)
                        profilesColorMap.Add(aux[i], Color.FromArgb(0, (int)(st*(i+1)),0));
                }                    
            }
            else
                throw new Exception("Descrete state are not defined");
        }
        double DrawNodeProfiles(Bitmap bmp, List<string> profNames,double profConsistency)
        {
            Graphics g = Graphics.FromImage(bmp);
            int width = bmp.Width;
            int height = bmp.Height;
            float xStep = (float)width;
            float yStep = (float)height;

            if (profilesColorMap == null)
                return 0;
            //if(upperLeaves.Count>1)
            xStep = (float)(width-20) / (profNames.Count);

            //if(profiles.Count>1)
            yStep = (float)height / (profiles[profNames[0]].Count);
            double[] avrProf = new double[profiles[profNames[0]].Count];
            float currentX = 0;
            float currentY = 0;
            SolidBrush b = new SolidBrush(Color.Black);
            for (int i = 0; i < profiles[profNames[0]].Count; i++)
            {
                currentY = i * yStep;

                for (int j = 0; j < profNames.Count; j++)
                {
                    double ind = profiles[profNames[j]][i];

                    Color c = profilesColorMap[ind];
                    b.Color = c;
                    currentX = j * xStep;
                    g.FillRectangle(b, currentX, currentY, xStep, yStep);

                    avrProf[i] += profiles[profNames[j]][i];
                }

            }
            currentX =currentX+xStep;
            //average
            for (int j = 0; j < avrProf.Length; j++)
            {
                avrProf[j] =Math.Round(avrProf[j]/ profNames.Count);                
               

                Color c = profilesColorMap[(int)avrProf[j]];
                b.Color = c;
                currentY = j * yStep;
                g.FillRectangle(b, width-20, currentY, 20, yStep);                
            }
            b.Color = Color.Red;
            Pen p = new Pen(b);
            g.DrawRectangle(p, width-20, 0, 20, yStep * avrProf.Length);


            return profConsistency;
        }
        private void toolStripButton3_Click(object sender, EventArgs e)
        {
            if (!linemode)
            {
                linemode = true;
                /*buffer = new Bitmap(pictureBox1.Width, pictureBox1.Height);

                Graphics g = Graphics.FromImage(buffer); ;

                g.SmoothingMode = SmoothingMode.AntiAlias;
                g.TextRenderingHint = System.Drawing.Text.TextRenderingHint.AntiAlias;

                DrawGraph(rootNode, g);
                DrawDistanceAx(rootNode, g);*/
            }
        }

        private void pictureBox1_MouseMove(object sender, MouseEventArgs e)
        {
            if (linemode)
            {
                mposX = mposY;
                mposY = e.Y;
                listNodes=drawH.rootNode.CutDendrog(-(int)((e.Y-drawH.maxGraphicsY - drawH.posStart)/drawH.distanceStepY));
                pictureBox1.Invalidate();
            }
        }

        private void toolStripButton4_Click(object sender, EventArgs e)
        {
            markflag = !markflag;
            
            if(markflag)
                toolStripButton4.Image = ((System.Drawing.Image)phiClustCore.Properties.Resources.Flag2);
            else
                toolStripButton4.Image = ((System.Drawing.Image)phiClustCore.Properties.Resources.Flag);
        }

        private void Form2_KeyPress(object sender, KeyPressEventArgs e)
        {
            if (e.KeyChar == (char)27)
            {
                markflag = false ;
                linemode = false;
                if(listNodes!=null)
                    listNodes.Clear();
                SaveMarkedClusters.Enabled = false;
                pictureBox1.Invalidate();
            }
            
        }

        private void toolStripButton6_Click(object sender, EventArgs e)
        {
            if (listNodes != null)
            {
                listNodes.Clear();
                SaveMarkedClusters.Enabled = false;
                pictureBox1.Invalidate();
            }
        }

        private void toolStripButton5_Click(object sender, EventArgs e)
        {
            saveFileDialog1.Filter = "txt files|*.txt";
            saveFileDialog1.DefaultExt = "txt";
            if (saveFileDialog1.ShowDialog() == DialogResult.OK)
            {
                if (saveFileDialog1.FileName.Length > 0)
                {
                    string name = saveFileDialog1.FileName;
                    if (!name.Contains(".txt"))
                        name += ".txt";
                    StreamWriter file = new StreamWriter(name);

                    int i = 0;
                    foreach (var item in listNodes)
                    {
                        item.Key.SaveNode(file, i++);
                    }

                    file.Close();
                }
            }

        }
        public void ShowCloseButton()
        {
            toolStripButton5.Visible = true;
        }

        private void toolStripButton5_Click_1(object sender, EventArgs e)
        {
            this.DialogResult = DialogResult.OK;
            this.Close();  
        }


        private void toolStripButton8_Click(object sender, EventArgs e)
        {
            DialogResult res;
            SelectStruct s = new SelectStruct(drawH.hnode.setStruct);

            res=s.ShowDialog();
            if(res==DialogResult.OK)
            {
                drawH.hnode.ColorNode(s.selectedStruct, 2);
                buffer = null;
                pictureBox1.Invalidate();
                pictureBox1.Refresh(); 
            }
        }

        private void orderVis_Click(object sender, EventArgs e)
        {
            List<List<string>> clusterList=new List<List<string>>();
            List<int[]> colorMap = new List<int[]>();
            int[] tab = new int[4];
            tab[0] = 0; tab[1] = 65; tab[2] = 150; tab[3] = 235;

            for (int i = 0; i < tab.Length; i++)
                for (int j = 0; j < tab.Length; j++)
                    for (int n = 0; n < tab.Length; n++)
                    {
                        int[] aux = new int[3];
                        aux[0] = tab[i];
                        aux[1] = tab[j];
                        aux[2] = tab[n];

                        colorMap.Add(aux);
                    }
            int counter=0;
            List<Color> clusterColor = new List<Color>();
            List<HClusterNode> hk=new List<HClusterNode>(listNodes.Keys);
             VisOrder visFrame =null;
             if (hk.Count < colorMap.Count)
             {
                 for (int i = 0; i < hk.Count; i++)
                 {
                     clusterList.Add(hk[i].setStruct);
                     if (listNodes.Count < colorMap.Count)
                     {
                         Color r = Color.FromArgb(colorMap[counter][0], colorMap[counter][1], colorMap[counter][2]);
                         listNodes[hk[i]] = r;
                         clusterColor.Add(r);
                     }
                     counter += colorMap.Count / (hk.Count+10);
                 }
                 visFrame = new VisOrder(clusterList, null, clusterColor);
             }
             else
             {
                 for (int i = 0; i < hk.Count; i++)
                     clusterList.Add(hk[i].setStruct);
                 visFrame = new VisOrder(clusterList, null, null);
             }
            visFrame.Show();
            this.Invalidate();
            this.Refresh();
        }
/*        private void SetColors()
        {            
                if (vecColor == null)
                    return;
                classColor = new Dictionary<string, Color>();
                foreach (var item in vecColor)
                {
                    if (!classColor.ContainsKey(item.Value))
                        classColor.Add(item.Value, Color.Azure);
                }

                double step = (colorMap.Count - 2) / classColor.Keys.Count;

                List<string> ll = new List<string>(classColor.Keys);

                int count = 2;
                for (int i = 0; i < ll.Count; i++)
                {
                    classColor[ll[i]] = colorMap[count];
                    count = (int)((i + 1) * step);
                }
                buffer = null;
                pictureBox1.Invalidate();
                pictureBox1.Refresh();            

            
        }*/
        private void toolStripButton9_Click(object sender, EventArgs e)
        {
            DialogResult res= openFileDialog1.ShowDialog();

            if(res==DialogResult.OK)
            {
                drawH.ColorAndLabels(openFileDialog1.FileName);
                buffer = null;
                pictureBox1.Invalidate();
                pictureBox1.Refresh(); 
            }

        }

        private void toolStripButton7_Click(object sender, EventArgs e)
        {
            RuleEditor editor = new RuleEditor(geneLabels,upReg,downReg,threshold);
            DialogResult res=editor.ShowDialog();
           
            if(res==DialogResult.OK)
            {
                labels.Clear();
                listNodes = hnode.MarkNodesAccordingToRules(editor.upRegulated, editor.downRegulated, editor.threshold,Color.Red,profiles);
                if (listNodes?.Count > 0)
                {                    
                    for(int i=0;i<geneLabels.Count;i++)
                    {
                        if (editor.upRegulated.Contains(i))
                            labels.Add(new Tuple<string, string>(geneLabels[i], "U"));
                        else
                        if (editor.downRegulated.Contains(i))
                            labels.Add(new Tuple<string, string>(geneLabels[i], "D"));
                        else
                            labels.Add(new Tuple<string, string>(geneLabels[i], ""));

                    }
                    upReg = editor.upRegulated;
                    downReg = editor.downRegulated;
                    threshold = editor.threshold;
                    pictureBox1.Invalidate();
                }
            }
        }

        private void toolStripButton11_Click(object sender, EventArgs e)
        {
            sizeNodes = hnode.GetAllClusters();
            pictureBox1.Invalidate();
        }

        private void vScrollBar1_ValueChanged(object sender, EventArgs e)
        {
            drawH.sizeThreshold = vScrollBar1.Value;
            //            drawH.PrepareGraphNodes(buffer);            
            textBox1.Text = vScrollBar1.Value.ToString();
            drawH.DrawOnBuffer(buffer, true, 1, Color.Empty);

            pictureBox1.Invalidate();
            pictureBox1.Refresh();
        }

        private void toolStripButton10_Click(object sender, EventArgs e)
        {
            viewType = !viewType;
            drawH.PrepareGraphNodes(buffer);
            //ClearBuffer();
            buffer = null;
            pictureBox1.Invalidate();
           
        }



    }
}
