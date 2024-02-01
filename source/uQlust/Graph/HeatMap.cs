using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;
using phiClustCore;
using phiClustCore.Profiles;
using phiClustCore.Interface;
using System.Drawing.Imaging;
using phiClustCore.Distance;

namespace Graph
{
    public partial class HeatMap : Form, IVisual
    {
        HeatmapDrawCore hCore = new HeatmapDrawCore();

        int x_pos, y_pos;
        int mouseY, mouseX;
        bool cutHorizontalLine = false;
        bool cutVerticalLine = false;
        bool cutBranch = false;
        bool swap = false;
        DistanceMeasure distM = null;

        bool fullImage = false;
        HClusterNode leftRoot=null;
        HClusterNode nodeForKruskal = null;
        HClusterNode remleftRoot=null;
        HClusterNode remAuxLeft = null;
        HClusterNode colorNodeUpper = null;
        HClusterNode colorNodeLeft = null;

        //Bitmap heatmap=null;
        public void ToFront()
        {
            this.BringToFront();
        }

        public HeatMap(HClusterNode upperNode, HClusterNode leftNode, Dictionary<string, string> labels, ClusterOutput outp)
        {
            leftRoot=leftNode;
            hCore.upperLeaves = upperNode.GetLeaves();
            hCore.leftLeaves = leftRoot.GetLeaves();
            upperNode.ClearColors(Color.Black);
            leftRoot.ClearColors(Color.Black);
           
            InitializeComponent();
            this.Text = outp.name;
            //draw=new HeatMapDraw(new Bitmap(tableLayoutPanel1.Width, tableLayoutPanel1.Height), upperNode, leftNode, labels, outp);
            // draw = new HeatMapDraw(new Bitmap(pictureBox2.Width,pictureBox2.Height),new Bitmap(pictureBox3.Width,pictureBox3.Height),
            //     new Bitmap(pictureBox1.Width,pictureBox1.Height),new Bitmap(pictureBox4.Width,pictureBox4.Height),upperNode, leftNode, labels, outp);
            hCore.draw = new HeatMapDraw(new Bitmap(pictureBoxUpper.Width, pictureBoxUpper.Height), new Bitmap(pictureBoxLeft.Width, pictureBoxLeft.Height),
                new Bitmap(pictureBoxHeatMap.Width, pictureBoxHeatMap.Height), new Bitmap(pictureBox4.Width, pictureBox4.Height), upperNode, leftRoot, labels, outp);

            label1.Text = "Gene clusters: " + upperNode.GetLeaves().Count;
            label2.Text = "Sample clusters: " + leftRoot.GetLeaves().Count;
            distM = outp.distM;
            int max = 0;
            foreach (var item in leftRoot.GetLeaves())
                if (item.setStruct.Count > max)
                    max = item.setStruct.Count;

            label3.Text = "0";
            vScrollBar1.Minimum = 0;
            vScrollBar1.Maximum = max;
        }
        public void PrepareDataForHeatMap()
        {
            hCore.draw.PrepareDataForHeatMap();

            this.Name = "HeatMap " + hCore.draw.outp.dirName;

            comboBox1.Items.Add(".....................");
            for (int i = 1; i < hCore.draw.outp.aux2.Count; i++)
                comboBox1.Items.Add(hCore.draw.outp.aux2[i]);

            if (comboBox1.Items.Count > 0)
                comboBox1.SelectedIndex = 0;
            else
                comboBox1.Visible = false;


            pictureBoxUpper.Refresh();
        }
        void PrepareBarRegions(Graphics g)
        {
            int x = 50, y = 0;
            bool regionTest = false;

            if (hCore.draw.upper.labColor == null || hCore.draw.upper.labColor.Keys.Count == 0)
                return;

            if (hCore.regionBarColor.Count > 0)
                regionTest = true;

            List<string> labKeys = new List<string>(hCore.draw.upper.labColor.Keys);
            Font drawFont = new System.Drawing.Font("Arial", 8);
            foreach (var item in labKeys)
            {
                if (!regionTest)
                {
                    Region reg = new Region(new Rectangle(x, y, 15, 10));
                    hCore.regionBarColor.Add(reg, new KeyValuePair<bool, string>(false, item));
                    SizeF textSize = g.MeasureString(item, drawFont);
                    x += 25 + (int)textSize.Width;
                    if (x > this.Width)
                    {
                        y += 150;
                        x = 25;
                    }
                }
            }
            foreach (var regItem in hCore.regionBarColor)
            {
                if (regItem.Value.Key)
                    hCore.draw.upper.labColor[regItem.Value.Value] = Color.Empty;
            }

        }
        private void pictureBoxUpper_Paint(object sender, PaintEventArgs e)
        {
            Bitmap bit = hCore.PaintUpper(hCore.draw.upperBitMap);
            if (cutHorizontalLine)
            {
                Pen p = new Pen(Color.Red);
                Graphics g = Graphics.FromImage(bit);
                g.DrawLine(p, 0, mouseY, bit.Width, mouseY);
            }

            e.Graphics.DrawImage(bit, 0, 0);
        }
        private void pictureBox3_Paint(object sender, PaintEventArgs e)
        {
            e.Graphics.Clear(this.BackColor);
            Bitmap bit = hCore.PaintLeft(hCore.draw.leftBitMap);
            if (cutVerticalLine)
            {
                Pen p = new Pen(Color.Red);
                Graphics g = Graphics.FromImage(bit);
                g.DrawLine(p, mouseX, 0, mouseX, bit.Height);
            }

            e.Graphics.DrawImage(bit, 0, 0);

            //draw.leftBitMap.Save("test.jpg");
        }

        private void pictureBoxHeatMap_Paint(object sender, PaintEventArgs e)
        {
            Bitmap bit = hCore.PaintHeatMap(hCore.draw.heatMap);
            e.Graphics.DrawImage(bit, x_pos, y_pos);
            
        }
        void LeftClusterRearange()
        {
            MST ms = new MST();
            GraphMST gr = ms.CreateGraph(distM, hCore.leftLeaves);

            if (nodeForKruskal != null)
            {
                List<List<string>> kruskalOut = ms.Kruskal(gr, nodeForKruskal);

                remleftRoot = leftRoot;
                leftRoot = HashClusterDendrog.MakeDummyDendrog(kruskalOut).hNode;
                hCore.draw.leftNode = leftRoot;
                hCore.draw.left.rootNode = leftRoot;
                remAuxLeft = hCore.draw.auxLeft;
                hCore.draw.auxLeft = leftRoot;
                hCore.leftLeaves = leftRoot.GetLeaves();
                hCore.draw.left.PrepareGraphNodes(hCore.draw.leftBitMap, hCore.verticalCuttDistance);
                pictureBoxLeft.Refresh();
                pictureBoxHeatMap.Refresh();
            }

        }
        void LeftBackToHistogram()
        {
            if (remleftRoot != null)
            {
                leftRoot = remleftRoot;
                hCore.draw.leftNode = remleftRoot;
                hCore.draw.left.rootNode = remleftRoot;
                hCore.draw.auxLeft = remAuxLeft;
                hCore.leftLeaves = leftRoot.GetLeaves();
                hCore.draw.left.PrepareGraphNodes(hCore.draw.leftBitMap, hCore.verticalCuttDistance);

                pictureBoxLeft.Refresh();
                pictureBoxHeatMap.Refresh();
            }

        }
        private void HeatMap_ResizeEnd(object sender, EventArgs e)
        {
            hCore.ResizeHeatMap(pictureBoxUpper.Size, pictureBoxLeft.Size, pictureBox6.Size);
            if (!fullImage)
            {
                hCore.draw.upperBitMap = new Bitmap(pictureBoxUpper.Width, pictureBoxUpper.Height);
                hCore.draw.leftBitMap = new Bitmap(pictureBoxLeft.Width, pictureBoxLeft.Height);
                hCore.draw.heatMap = new Bitmap(pictureBoxUpper.Width, pictureBoxLeft.Height);
            }
            hCore.draw.left.PrepareGraphNodes(hCore.draw.leftBitMap, hCore.verticalCuttDistance);
            hCore.draw.upper.PrepareGraphNodes(hCore.draw.upperBitMap, hCore.horizontalCuttDistance);
            pictureBoxHeatMap.Refresh();
            pictureBoxUpper.Refresh();
            pictureBoxLeft.Refresh();
            hCore.PaintLeftDescription(pictureBox6.Width);
            pictureBox6.Refresh();
            pictureBox7.Refresh();

            //this.Invalidate();
        }
        private void pictureBox4_Paint(object sender, PaintEventArgs e)
        {
            Bitmap b = new Bitmap(pictureBox4.Width, pictureBox4.Height);
            hCore.PaintLegend(b);
            e.Graphics.DrawImage(b,0,0);
            //test.Paint();             
        }

        private void tableLayoutPanel1_Paint_1(object sender, PaintEventArgs e)
        {

        }

        private void pictureBoxUpper_MouseMove(object sender, MouseEventArgs e)
        {
            bool test = false;
            mouseY = e.Y;

            if (cutHorizontalLine)
                pictureBoxUpper.Refresh();
            else
            {

                if (colorNodeUpper != null)
                {
                    hCore.draw.upper.ChangeColors(colorNodeUpper, Color.Black);
                    test = true;
                }

                colorNodeUpper = hCore.draw.upper.FindClosestNode(e.X, e.Y);



                if (colorNodeUpper != null)
                {
                    hCore.draw.upper.ChangeColors(colorNodeUpper, Color.Red);
                    test = true;
                }
                if (test)
                {
                    Graphics g = Graphics.FromImage(hCore.draw.upperBitMap);
                    g.Clear(pictureBoxUpper.BackColor);

                    if (colorNodeUpper != null)
                    {
                        float v = ((float)colorNodeUpper.setStruct.Count) / hCore.draw.upperNode.setStruct.Count * 360;
                        g.FillPie(new SolidBrush(Color.Black), new Rectangle(e.X - 20, e.Y, 15, 15), v, 360 - v);
                        g.FillPie(new SolidBrush(Color.Red), new Rectangle(e.X - 20, e.Y, 15, 15), 0, v);
                    }
                    pictureBoxUpper.Refresh();

                }
            }
        }
        private void pictureBoxUpper_MouseClick(object sender, MouseEventArgs e)
        {
            switch (e.Button)
            {

                case MouseButtons.Left:

                    if (hCore.regionBarColor.Count > 0)
                    {
                        List<Region> l = new List<Region>(hCore.regionBarColor.Keys);
                        bool testFlag = false;
                        foreach (var item in l)
                            if (item.IsVisible(e.X, e.Y))
                            {
                                hCore.regionBarColor[item] = new KeyValuePair<bool, string>(!hCore.regionBarColor[item].Key, hCore.regionBarColor[item].Value);
                                testFlag = true;
                            }
                        if (testFlag)
                        {

                            pictureBoxUpper.Refresh();
                            return;
                        }
                    }
                    HClusterNode nodeC = colorNodeUpper;//upper.CheckClick(upper.rootNode,e.X,e.Y);
                    if (nodeC != null && nodeC.joined == null && !cutBranch)
                    {
                        TextBoxView rr = new TextBoxView(nodeC.setStruct);

                        rr.Show();

                        DrawPanel pn = new DrawPanel("Leave profile: " + nodeC.refStructure, null);
                        int tmp = pn.Height;
                        pn.Height = pn.Width;
                        pn.Width = tmp;
                        pn.CreateBMP();
                        List<HClusterNode> leftLeaves = hCore.draw.auxLeft.GetLeaves();
                        leftLeaves = leftLeaves.OrderByDescending(o => o.gNode.y).Reverse().ToList();
                        List<string> refList = new List<string>();
                        foreach (var item in leftLeaves)
                            refList.Add(item.refStructure);
                        pn.drawPic = delegate { hCore.draw.DrawHeatMapNode(pn.bmp, refList, nodeC.setStruct); pn.Text = nodeC.consistency.ToString(); };
                        pn.Show();
                    }
                    else
                    {
                        if (cutBranch)
                        {
                            cutBranch = false;
                            if (colorNodeUpper != null)
                            {
                                colorNodeUpper.MakeInvisible();
                                hCore.draw.upper.PrepareGraphNodes(hCore.draw.upperBitMap, hCore.horizontalCuttDistance);
                                pictureBoxUpper.Refresh();

                            }
                        }
                        else
                            if (colorNodeUpper != null && colorNodeUpper.joined != null)
                                hCore.draw.auxUpper = colorNodeUpper;
                    }
                    if (cutHorizontalLine)
                    {
                        hCore.draw.upper.MakeAllVisible();
                        cutHorizontalLine = !cutHorizontalLine;
                        double distScale = hCore.draw.upper.maxRealDist - hCore.draw.upper.minRealDist;
                        distScale /= 100;
                        hCore.horizontalCuttDistance = hCore.draw.upper.minRealDist + distScale * ((double)(hCore.draw.upper.maxGraphicsY - e.Y)) / hCore.draw.upper.maxGraphicsY;
                        hCore.draw.upper.MakeGraphInvisible(hCore.draw.upper.rootNode,hCore.horizontalCuttDistance);
                    }
                    label1.Text = "Gene clusters: " + hCore.draw.upper.rootNode.GetLeaves().Count;
                    break;
                case MouseButtons.Right:
                    if (!swap)
                    {
                        hCore.Swap();
                    } 
                    else
                    {
                        if (colorNodeUpper != null && colorNodeUpper.joined != null)
                        {
                            var aux = colorNodeUpper.joined[0];
                            colorNodeUpper.joined[0] = colorNodeUpper.joined[colorNodeUpper.joined.Count - 1];
                            colorNodeUpper.joined[colorNodeUpper.joined.Count - 1] = aux;
                        }
                    }
                    break;
            }
            if (hCore.draw.auxUpper != null)
            {
                if (colorNodeUpper != null)
                    hCore.draw.upper.ChangeColors(colorNodeUpper, Color.Black);
                colorNodeUpper = null;
                hCore.draw.upper.rootNode = hCore.draw.auxUpper;
                hCore.draw.upper.PrepareGraphNodes(hCore.draw.upperBitMap, hCore.horizontalCuttDistance);
                Graphics g = Graphics.FromImage(hCore.draw.upperBitMap);
                g.Clear(pictureBoxUpper.BackColor);
                pictureBoxUpper.Refresh();
                pictureBoxHeatMap.Refresh();
            }
        }

        private void pictureBox3_MouseMove(object sender, MouseEventArgs e)
        {
            bool test = false;
            mouseX = e.X;
            if (cutVerticalLine)
                pictureBoxLeft.Refresh();
            else
            {
                if (colorNodeLeft != null)
                {
                    test = true;
                    hCore.draw.left.ChangeColors(colorNodeLeft, Color.Black);
                }

                colorNodeLeft = hCore.draw.left.FindClosestNode(e.X, e.Y);


                if (colorNodeLeft != null)
                {
                    hCore.draw.left.ChangeColors(colorNodeLeft, Color.Red);
                    test = true;
                }
                if (test)
                {
                    Graphics g = Graphics.FromImage(hCore.draw.leftBitMap);
                    g.Clear(pictureBoxLeft.BackColor);

                    if (colorNodeLeft != null)
                    {
                        float v = ((float)colorNodeLeft.setStruct.Count) / hCore.draw.leftNode.setStruct.Count * 360;
                        g.FillPie(new SolidBrush(Color.Black), new Rectangle(e.X - 20, e.Y, 15, 15), v, 360 - v);
                        g.FillPie(new SolidBrush(Color.Red), new Rectangle(e.X - 20, e.Y, 15, 15), 0, v);
                    }
                    //g.DrawArc(new Pen(Color.DarkGreen),new Rectangle(e.X,e.Y,10,10), 0.0, v);
                    pictureBoxLeft.Refresh();
                }
            }

        }

        private void pictureBox3_MouseClick(object sender, MouseEventArgs e)
        {
            switch (e.Button)
            {

                case MouseButtons.Left:
                    HClusterNode nodeC = colorNodeLeft;//  left.CheckClick(left.rootNode,e.X,e.Y);

                    if (nodeC != null && nodeC.joined == null && !cutBranch)
                    {                        
                        nodeForKruskal = nodeC;
                        return;
                    }
                    else                    
                     if (cutBranch)
                    {
                        cutBranch = false;
                        if (colorNodeLeft != null)
                        {
                            colorNodeLeft.MakeInvisible();
                            hCore.draw.left.PrepareGraphNodes(hCore.draw.leftBitMap,hCore.verticalCuttDistance);
                            pictureBoxLeft.Refresh();

                        }
                    }
                    else
                        if (colorNodeLeft != null && colorNodeLeft.joined != null)
                            hCore.draw.auxLeft = colorNodeLeft;

                    if (cutVerticalLine)
                    {
                        cutVerticalLine = !cutVerticalLine;
                        double distScale = hCore.draw.left.maxRealDist - hCore.draw.left.minRealDist;
                        distScale /= 100;
                        hCore.verticalCuttDistance = hCore.draw.left.minRealDist + distScale * ((double)(hCore.draw.left.maxGraphicsX - e.X)) / hCore.draw.left.maxGraphicsX;
                    }
                    label2.Text = "Sample clusters: " + hCore.draw.left.rootNode.GetLeaves().Count;                    
                    break;
                case MouseButtons.Right:
                    if(!swap)
                        hCore.draw.auxLeft = hCore.draw.leftNode;
                    else
                    {
                        if (colorNodeLeft == null)
                            colorNodeLeft=hCore.draw.leftNode;
                        if(colorNodeLeft!=null && colorNodeLeft.joined!=null)
                        {
                            var aux = colorNodeLeft.joined[0];
                            colorNodeLeft.joined[0] = colorNodeLeft.joined[colorNodeLeft.joined.Count - 1];
                            colorNodeLeft.joined[colorNodeLeft.joined.Count - 1] = aux;
                        }
                        
                    }
                    break;
            }
            if (hCore.draw.auxLeft != null)
            {
                if (colorNodeLeft != null)
                    hCore.draw.left.ChangeColors(colorNodeLeft, Color.Black);
                colorNodeLeft = null;
                hCore.draw.left.rootNode = hCore.draw.auxLeft;
                hCore.leftLeaves = hCore.draw.left.rootNode.GetLeaves();
                hCore.draw.left.PrepareGraphNodes(hCore.draw.leftBitMap, hCore.verticalCuttDistance);
                Graphics g = Graphics.FromImage(hCore.draw.leftBitMap);
                g.Clear(pictureBoxLeft.BackColor);
                pictureBoxLeft.Refresh();
                pictureBoxHeatMap.Refresh();
                hCore.PaintLeftDescription(pictureBox6.Width);
                pictureBox6.Refresh();
                pictureBox7.Refresh();

            }
        }

        private void labels_Click(object sender, EventArgs e)
        {


        }

        private void toolStripLabel1_Click(object sender, EventArgs e)
        {
            hCore.showLabels = !hCore.showLabels;
            hCore.draw.upper.showLabels =hCore.showLabels;
            numericUpDown1.Visible = hCore.showLabels;
            numericUpDown1.Value = hCore.draw.upper.labelSize;
            hCore.showSelectedLabels = false;
            //draw.left.showLabels = showLabels;
            pictureBoxUpper.Refresh();
        }

        private void numericUpDown1_ValueChanged(object sender, EventArgs e)
        {
            hCore.draw.left.labelSize = (int)numericUpDown1.Value;
            hCore.draw.upper.labelSize = (int)numericUpDown1.Value;

            Graphics g = Graphics.FromImage(hCore.draw.upperBitMap);
            g.Clear(pictureBoxUpper.BackColor);
            g = Graphics.FromImage(hCore.draw.leftBitMap);
            g.Clear(pictureBoxLeft.BackColor);

            pictureBoxLeft.Refresh();
            pictureBoxUpper.Refresh();

        }

        private void pictureBoxUpper_MouseLeave(object sender, EventArgs e)
        {
            if (colorNodeUpper != null)
            {
                hCore.draw.upper.ChangeColors(colorNodeUpper, Color.Black);
                colorNodeUpper = null;
                Graphics g = Graphics.FromImage(hCore.draw.upperBitMap);
                g.Clear(pictureBoxUpper.BackColor);
                pictureBoxUpper.Refresh();
            }
        }

        private void toolStripButton1_Click(object sender, EventArgs e)
        {
            DialogResult res = saveFileDialog1.ShowDialog();
            if (res == DialogResult.OK)
            {
                string fileName = saveFileDialog1.FileName;
                List<List<string>> clusters = new List<List<string>>();
                List<HClusterNode> nodes = hCore.draw.leftNode.GetLeaves();
                List<double> cons = new List<double>();
                foreach (var item in nodes)
                {
                    clusters.Add(item.setStruct);
                    cons.Add(item.consistency);
                }

                StreamWriter wS = new StreamWriter(fileName + "_gene_microclusters.dat");
                ClusterOutput.Save(clusters, cons, wS, true);
                wS.Close();

                clusters.Clear();
                cons = new List<double>();
                nodes = hCore.draw.upperNode.GetLeaves();
                foreach (var item in nodes)
                {
                    clusters.Add(item.setStruct);
                    cons.Add(item.consistency);
                }

                wS = new StreamWriter(fileName + "_sample_microclusters.dat");
                ClusterOutput.Save(clusters, cons, wS, true);
                wS.Close();


            }
        }

        private void pictureBox1_Click(object sender, EventArgs e)
        {

        }

        private void pictureBox3_MouseLeave(object sender, EventArgs e)
        {
            if (colorNodeLeft != null)
            {
                hCore.draw.left.ChangeColors(colorNodeLeft, Color.Black);
                colorNodeLeft = null;
                Graphics g = Graphics.FromImage(hCore.draw.leftBitMap);
                g.Clear(pictureBoxLeft.BackColor);

                //g.DrawArc(new Pen(Color.DarkGreen),new Rectangle(e.X,e.Y,10,10), 0.0, v);
                pictureBoxLeft.Refresh();
            }
        }

        private void comboBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            if (comboBox1.SelectedItem == null)
                return;
            if (comboBox1.SelectedItem.ToString().Contains("..."))
            {
                hCore.draw.upper.rootNode.ClearColors(Color.Black);
            }
            else
            {
                List<HClusterNode> leaveas = hCore.draw.upperNode.GetLeaves();
                int index = 1;
                hCore.regionBarColor.Clear();
                foreach (var item in leaveas)
                {
                    if (item.refStructure.Contains("G2M"))
                        Console.WriteLine();
                    if (item.refStructure.Contains(comboBox1.SelectedItem.ToString()))
                    {
                        hCore.draw.upper.ChangeColors(item, Color.Red);
                        break;
                    }
                }
            }
            pictureBoxUpper.Refresh();
        }

        private void panel2_Scroll(object sender, ScrollEventArgs e)
        {
            x_pos = -e.NewValue;
            pictureBoxHeatMap.Invalidate();
        }

        private void panel3_Scroll(object sender, ScrollEventArgs e)
        {

            if (e.ScrollOrientation == ScrollOrientation.VerticalScroll)
            {
                y_pos = -e.NewValue;
                pictureBoxHeatMap.Invalidate();
            }
            pictureBoxLeft.Invalidate();

        }

        private void toolStrip1_ItemClicked(object sender, ToolStripItemClickedEventArgs e)
        {
            cutVerticalLine = !cutVerticalLine;
        }

        private void toolStripButton2_Click(object sender, EventArgs e)
        {
            cutVerticalLine = !cutVerticalLine;
        }

        private void vScrollBar1_ValueChanged(object sender, EventArgs e)
        {
            hCore.draw.left.clusterSize = vScrollBar1.Value;
            label3.Text = vScrollBar1.Value.ToString();
            pictureBoxLeft.Invalidate();
        }

        private void interactiveLabelColumnToolStripMenuItem_Click(object sender, EventArgs e)
        {
            interactiveLabelColumn.Checked = !interactiveLabelColumn.Checked;
            if (interactiveLabelColumn.Checked)
                interactiveLabelRow.Checked = false;

        }

        private void interactiveLabelRowToolStripMenuItem_Click(object sender, EventArgs e)
        {
            interactiveLabelRow.Checked = !interactiveLabelRow.Checked;
            if (interactiveLabelRow.Checked)
                interactiveLabelColumn.Checked = false;
        }
        string GetLabel(int x, int y, bool row, List<HClusterNode> nodes)
        {
            int pos = y;
            if (!row)
                pos = x;

            foreach (var item in nodes)
            {
                if (pos >= item.gNode.areaLeft && pos < item.gNode.areaRight)
                {
                    return item.refStructure;
                }
            }
            return "";
        }

        private void pictureBoxHeatMap_MouseMove(object sender, MouseEventArgs e)
        {
        }

        private void button1_Click(object sender, EventArgs e)
        {
            splitContainer1.SplitterDistance += 20;
            hCore.draw.leftBitMap = new Bitmap(pictureBoxLeft.Width, pictureBoxLeft.Height);
            hCore.draw.heatMap = new Bitmap(pictureBoxUpper.Width, pictureBoxLeft.Height);
            hCore.draw.left.PrepareGraphNodes(hCore.draw.leftBitMap);

            pictureBoxLeft.Invalidate();
            pictureBoxUpper.Invalidate();

        }


        Tuple<Dictionary<string, Tuple<Color, string>>, List<ColorRect>> CalcLeftBarsClusterAvr(Dictionary<string, Tuple<Color, string>> colorMap, int pos = 0, int width = 20)
        {
            List<ColorRect> r = new List<ColorRect>();
            foreach (var item in hCore.leftLeaves)
            {
                float startPos = item.gNode.areaLeft;
                float step = ((float)(item.gNode.areaRight - item.gNode.areaLeft)) / item.setStruct.Count;
                List<Color> temp = new List<Color>();
                for (int i = 0; i < item.setStruct.Count; i++)
                    if (colorMap.ContainsKey(item.setStruct[i]))
                        temp.Add(colorMap[item.setStruct[i]].Item1);

                Dictionary<Color, int> freq = new Dictionary<Color, int>();
                foreach (var it in temp)
                {
                    if (!freq.ContainsKey(it))
                        freq.Add(it, 0);

                    freq[it]++;
                }
                List<KeyValuePair<Color, int>> mappings = freq.ToList();
                if (mappings.Count > 0)
                {
                    mappings.Sort((x, y) => x.Value.CompareTo(y.Value));
                    ColorRect z = new ColorRect();
                    z.sY = item.gNode.areaLeft;
                    z.sX = 0 + pos;
                    z.width = width;
                    z.height = item.gNode.areaRight - item.gNode.areaLeft;
                    z.c = mappings[mappings.Count - 1].Key;
                    r.Add(z);
                }
            }
            return new Tuple<Dictionary<string, Tuple<Color, string>>, List<ColorRect>>(colorMap, r);
        }

        Tuple<Dictionary<string, Tuple<Color, string>>, List<ColorRect>> CalcLeftBars(Dictionary<string, Tuple<Color, string>> colorMap, int pos = 0, int width = 20)
        {
            List<ColorRect> r = new List<ColorRect>();
            foreach (var item in hCore.leftLeaves)
            {
                float startPos = item.gNode.areaLeft;
                float step = ((float)(item.gNode.areaRight - item.gNode.areaLeft)) / item.setStruct.Count;
                List<ColorRect> temp = new List<ColorRect>();
                for (int i = 0; i < item.setStruct.Count; i++)
                    if (colorMap.ContainsKey(item.setStruct[i]))
                    {
                        ColorRect x = new ColorRect();
                        x.sX = 0 + pos;
                        x.sY = startPos + step * i;
                        x.width = width;
                        x.height = step;

                        x.c = colorMap[item.setStruct[i]].Item1;
                        x.text = colorMap[item.setStruct[i]].Item2;
                        //r.Add(x);
                        temp.Add(x);
                    }
                //  if (step < 1)
                {
                    HashSet<int> ss = new HashSet<int>();
                    foreach (var it in temp)
                    {
                        ss.Add((int)it.sY);
                    }

                    Dictionary<int, Dictionary<Color, int>> freq = new Dictionary<int, Dictionary<Color, int>>();
                    foreach (var it in temp)
                    {
                        if (!freq.ContainsKey((int)Math.Floor(it.sY)))
                            freq.Add((int)Math.Floor(it.sY), new Dictionary<Color, int>());

                        if (!freq[(int)Math.Floor(it.sY)].ContainsKey(it.c))
                            freq[(int)Math.Floor(it.sY)].Add(it.c, 1);
                        else
                            freq[(int)Math.Floor(it.sY)][it.c]++;
                    }
                    foreach (var it in freq.Keys)
                    {
                        
                        List<KeyValuePair<Color, int>> mappings = freq[it].ToList();
                        KeyValuePair<Color, int> gkeyPair = mappings.OrderByDescending(MaxGuid => MaxGuid.Value).First();
                        double colorBrihtness = 0;
                        foreach(var am in mappings)                        
                            colorBrihtness += am.Value;
                        colorBrihtness =gkeyPair.Value/colorBrihtness;
                        if (mappings.Count > 0)
                        {
                            mappings.Sort((x, y) => x.Value.CompareTo(y.Value));
                            for (int i = 0; i < temp.Count; i++)
                                if ((int)Math.Floor(temp[i].sY) == it)
                                {
                                    Color col= mappings[mappings.Count - 1].Key;
                                    int re = col.R - System.Convert.ToInt32(colorBrihtness / (double)100 * col.R);
                                    int gr = col.G - System.Convert.ToInt32(colorBrihtness / (double)100 * col.G);
                                    int bl = col.B - System.Convert.ToInt32(colorBrihtness / (double)100 * col.B);
                                    //temp[i].c = mappings[mappings.Count - 1].Key;
                                    temp[i].c = Color.FromArgb(255,re, gr, bl);
                                    temp[i].height = 1;
                                    if (step > 1)
                                        temp[i].height = step;

                                    temp[i].sY = (int)Math.Floor(temp[i].sY);
                                    r.Add(temp[i]);
                                    break;
                                }
                        }
                    }

                }
                //else
                //    r.Add(temp[0]);


            }
            return new Tuple<Dictionary<string, Tuple<Color, string>>, List<ColorRect>>(colorMap, r);
        }
        private void noteGenesToolStripMenuItem_Click(object sender, EventArgs e)
        {

            DialogResult res = openFileDialog1.ShowDialog();
            if (res == DialogResult.OK)
            {
                hCore.ReadNoteGenes(openFileDialog1.FileName);
                pictureBoxUpper.Refresh();
            }
        }

        private void rowAnnotationClick(object sender, EventArgs e)
        {
            DialogResult res = openFileDialog1.ShowDialog();
            if (res == DialogResult.OK)
            {
                checkBox1.Visible = true;
                var annotInfo = hCore.ReadLeftAnnotation(openFileDialog1.FileName);
                
                LocalButton b = new LocalButton();

                b.Location = new Point(1, annotInfo.Item1);
                b.Text = "X";
                b.Size = new Size(10, 10);
                hCore.PaintLeftDescription(pictureBox6.Width);
                pictureBox6.Refresh();
                pictureBox7.Refresh();
                b.tupleName = annotInfo.Item2;
                b.Click += (s, ev) =>
                {
                    hCore.leftDescription.Remove(((LocalButton)s).tupleName);
                    pictureBox7.Controls.Remove((LocalButton)s);
                    pictureBox6.Invalidate();
                    pictureBox7.Invalidate();
                };
                pictureBox7.Controls.Add(b);
            }
        }
        private void pictureBox6_Paint(object sender, PaintEventArgs e)
        {
            Bitmap annotation = new Bitmap(pictureBox6.Width, pictureBox6.Height);
            Bitmap bit =hCore.PaintLeftAnnotation(annotation);
            e.Graphics.DrawImage(bit, 0, 0);
        }

        private void pictureBox6_MouseMove(object sender, MouseEventArgs e)
        {
            /*if (leftDescription != null)
            {
                mouseX = e.X + 5;
                mouseY = e.Y;

                foreach (var item in leftDescription)
                {
                    for (int j = 0; j < item.Value.Item2.Count; j++)
                    {
                        ColorRect r = item.Value.Item2[j];
                        if (r.text != null)
                            if (e.X >= r.sX && e.X <= r.sX + r.width)
                                if (e.Y >= r.sY && e.Y <= r.sY + r.height)
                                    labelNum = r.text.ToString();
                    }
                }
                pictureBox6.Invalidate();
            }*/
        }

        private void selectGenesToolStripMenuItem_Click(object sender, EventArgs e)
        {

            List<string> genes = new List<string>();
            foreach (var item in hCore.upperLeaves)
                genes.Add(item.refStructure);
            GeneSelection sel = new GeneSelection(genes);
            DialogResult res = sel.ShowDialog();
            if (res == DialogResult.OK)
            {
                hCore.geneSelection = sel.selection;
                pictureBoxHeatMap.Invalidate();
            }
        }
        Dictionary<string, Color> GetDistinctColors(Tuple<Dictionary<string, Tuple<Color, string>>, List<ColorRect>> item)
        {
            Dictionary<string, Color> x = new Dictionary<string, Color>();

            foreach (var it in item.Item1)
            {
                var y = it.Value;
                if (!x.ContainsKey(y.Item2.ToString()))
                {
                    x.Add(y.Item2.ToString(), y.Item1);
                }
            }
            return x;
        }
        private void pictureBox7_Paint(object sender, PaintEventArgs e)
        {
            
            Bitmap bit= new Bitmap(pictureBox7.Width, pictureBox7.Height);
            Graphics g = Graphics.FromImage(bit);
            g.Clear(Color.White);
            hCore.PaintRightAnnotation(bit);
            e.Graphics.DrawImage(bit, 0, 0);
        }

        private void checkBox1_CheckedChanged(object sender, EventArgs e)
        {
            if (checkBox1.Checked)
                hCore.avrRowAnnotations = true;
            else
                hCore.avrRowAnnotations = false;

            hCore.PaintLeftDescription(pictureBox6.Width);
            pictureBox6.Refresh();
            pictureBox7.Refresh();
        }

        private void histogramsToolStripMenuItem_Click(object sender, EventArgs e)
        {
            PlotForm pHist = new PlotForm(FilterOmics.memoryFilteredData);
            pHist.Show();
        }

        private void toolStripButton3_Click(object sender, EventArgs e)
        {
            cutHorizontalLine = !cutHorizontalLine;
        }


        private void HeatMap_KeyDown(object sender, KeyEventArgs e)
        {
            if (e.KeyCode == Keys.Escape)
            {
                cutHorizontalLine = false;
                cutVerticalLine = false;
            }
        }

        private void pictureBoxUpper_PreviewKeyDown(object sender, PreviewKeyDownEventArgs e)
        {
            if (e.KeyCode == Keys.Escape)
                cutHorizontalLine = false;
        }
        class LocalButton : Button
        {
            public string tupleName = "";

            private void InitializeComponent()
            {
                this.SuspendLayout();
                this.ResumeLayout(false);

            }
        }

        private void ClearDistance_Click(object sender, EventArgs e)
        {
            cutHorizontalLine = false;
            cutVerticalLine = false;
            hCore.ClearDistance();
            pictureBoxLeft.Refresh();
            pictureBoxUpper.Refresh();
            pictureBoxHeatMap.Refresh();
        }

        private void toolStripButton4_Click(object sender, EventArgs e)
        {
            cutBranch = true;
        }

        private void toolStripButton5_Click(object sender, EventArgs e)
        {           
            hCore.upperLeaves = hCore.draw.auxUpper.GetLeaves();
            double threshold = Convert.ToDouble(textBox1.Text);
            foreach (var up in hCore.upperLeaves)                
            {
                up.visible = true;
                bool test = false;
                foreach (var item in hCore.leftLeaves)
                {
                    if (item.setStruct.Count > vScrollBar1.Value)
                    {
                        double res = hCore.draw.KLIndex(item, up);
                        if (res > threshold)
                        {
                            test = true;
                            break;
                        }
                    }
                }
                    if (!test)
                        up.visible = false;
                
            }
            
            hCore.draw.upper.PrepareGraphNodes(hCore.draw.upperBitMap, hCore.horizontalCuttDistance);
            pictureBoxUpper.Refresh();
            pictureBoxHeatMap.Refresh();
        }

        private void kruskalToolStripMenuItem_Click(object sender, EventArgs e)
        {
            LeftClusterRearange();
            hCore.PaintLeftDescription(pictureBox6.Width);
            pictureBox6.Refresh();
            backToHistogramToolStripMenuItem.Enabled = true;
            kruskalToolStripMenuItem.Enabled = false;
            minClusterDistanceToolStripMenuItem.Enabled = true;
        }

        private void backToHistogramToolStripMenuItem_Click(object sender, EventArgs e)
        {
            LeftBackToHistogram();
            hCore.PaintLeftDescription(pictureBox6.Width);
            pictureBox6.Refresh();
            backToHistogramToolStripMenuItem.Enabled = false;
            kruskalToolStripMenuItem.Enabled = true;
            minClusterDistanceToolStripMenuItem.Enabled = true;

        }

        private void minClusterDistanceToolStripMenuItem_Click(object sender, EventArgs e)
        {
            MST ms = new MST();
            GraphMST gr = ms.CreateGraph(distM, hCore.leftLeaves);

            if (nodeForKruskal != null)
            {
                List<List<string>> kruskalOut = ms.Closest(nodeForKruskal).GetClusters();

                remleftRoot = leftRoot;
                leftRoot = HashClusterDendrog.MakeDummyDendrog(kruskalOut).hNode;
                hCore.draw.leftNode = leftRoot;
                hCore.draw.left.rootNode = leftRoot;
                remAuxLeft = hCore.draw.auxLeft;
                hCore.draw.auxLeft = leftRoot;
                hCore.leftLeaves = leftRoot.GetLeaves();
                hCore.draw.left.PrepareGraphNodes(hCore.draw.leftBitMap, hCore.verticalCuttDistance);
                pictureBoxLeft.Refresh();
                pictureBoxHeatMap.Refresh();

                hCore.PaintLeftDescription(pictureBox6.Width);
                pictureBox6.Refresh();
                backToHistogramToolStripMenuItem.Enabled = true;
                minClusterDistanceToolStripMenuItem.Enabled = false;
                kruskalToolStripMenuItem.Enabled = true;

            }

        }

        private void saveToolStripMenuItem_Click(object sender, EventArgs e)
        {
            saveFileDialog1.DefaultExt = "png";
            saveFileDialog1.Filter = "Png files|*.png";
            DialogResult res = saveFileDialog1.ShowDialog();
            if (res == DialogResult.OK && saveFileDialog1.FileName.Length > 0)
            {
                Resolution resForm = new Resolution(this.Width, this.Height, Color.Black);
                res = resForm.ShowDialog();
                if (res == DialogResult.OK)
                {

                    hCore.Save(saveFileDialog1.FileName,resForm.WidthR,resForm.HeightR);
                     
                    

                }
                //       this.SavePicture(saveFileDialog1.FileName);          
            }

        }

        private void LongLeaves_Click(object sender, EventArgs e)
        {
            hCore.LongLeavesClick();
            pictureBoxUpper.Refresh();
            pictureBoxHeatMap.Refresh();
        }

        private void toolStripButton7_Click(object sender, EventArgs e)
        {
            swap = !swap;
        }

        private void toolStripButton8_Click(object sender, EventArgs e)
        {
            hCore.showSelectedLabels = !hCore.showSelectedLabels;
            pictureBoxUpper.Refresh();
            pictureBoxLeft.Refresh();
        }

        private void toolStripButton9_Click(object sender, EventArgs e)
        {
            hCore.horizontalLines = !hCore.horizontalLines;
            pictureBoxHeatMap.Refresh();
        }
    }
}
