using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.Drawing;
using System.IO;
//using phiClustCore;


namespace phiClustCore
{
    public class DrawHierarchical
    {
        public HClusterNode hnode = null;
        public HClusterNode rootNode = null;
        public HClusterNode currentRootNode = null;
        public Dictionary<string, Color> classColor = null;
        Dictionary<string, string> vecColor = null;
        Dictionary<Region, string> colorRegions = new Dictionary<Region, string>();
        List<Color> colorMap = new List<Color>();
        Stack<HClusterNode> nodesStack = new Stack<HClusterNode>();
        public Dictionary<string, Color> labColor = null;
        public int maxGraphicsY;// = 2000;
        public int maxGraphicsX;
        public bool viewType = false;
        public bool horizontalView = true;
        public double distanceStepY;
        public double distanceStepX;
        public bool showLabels = false;
        public int labelSize = 5;
        public int sizeThreshold = 0;
        public int clusterSize = 0;
        double circleStep;
        public int posStart = 20;
        List<string> geneLabels = new List<string>();
        int maxDist;
        public int currentLabelIndex = -1;
        HClusterNode maxHDist;
        HClusterNode minHDist;
        public int minRealDist;
        public int maxRealDist;
        string measureName;
        public DrawHierarchical(HClusterNode hnode, string measureName, Dictionary<string, string> labels, Bitmap buffer, bool horizontalView)
        {
            this.hnode = this.rootNode = hnode;
            this.horizontalView = horizontalView;
            this.measureName = measureName;
            PrepareGraphNodes(buffer);

            colorMap = PrepareColorMap();

            vecColor = labels;
            SetColors();

        }
        public HClusterNode FindClosestNode(int x, int y)
        {
            HClusterNode returnNode = null;
            int remDist = int.MaxValue;
            int dist = int.MaxValue;
            Queue<HClusterNode> q = new Queue<HClusterNode>();
            q.Enqueue(rootNode);
            while (q.Count != 0)
            {
                HClusterNode node = q.Dequeue();
                if (node.visible)
                {
                    if (node.gNode != null)
                        dist = (node.gNode.x - x) * (node.gNode.x - x) + (node.gNode.y - y) * (node.gNode.y - y);

                    if (dist < remDist && dist < 100)
                    {
                        remDist = dist;
                        returnNode = node;
                    }

                }
                if (node.joined == null)
                    continue;
                foreach (var item in node.joined)
                    q.Enqueue(item);
            }

            return returnNode;
        }
        public int GetChildrenNum(HClusterNode hNode)
        {
            int counter = 0;
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;
            st.Push(hNode);
            while (st.Count != 0)
            {
                current = st.Pop();
                if (current.joined != null)
                {
                    foreach (var item in current.joined)
                        st.Push(item);
                }
                else
                    counter++;

            }
            return counter;
        }
        static public List<Color> PrepareColorMap()
        {
            List<Color> col = new List<Color>();
            int[] tab = new int[4];

            tab[0] = 0; tab[1] = 85; tab[2] = 170; tab[3] = 255;

            for (int i = 0; i < tab.Length; i++)
                for (int j = 0; j < tab.Length; j++)
                    for (int n = 0; n < tab.Length; n++)
                        col.Add(Color.FromArgb(tab[i], tab[j], tab[n]));

            return col;

        }
        public void SearchKmax(HClusterNode hNode)
        {
            Stack<HClusterNode> st = new Stack<HClusterNode>();
            HClusterNode current = null;
            maxDist = (int)hNode.levelDist;
            minRealDist = (int)hNode.levelDist;
            maxRealDist = (int)hNode.levelDist;
            st.Push(hNode);
            while (st.Count != 0)
            {
                current = st.Pop();
                if (current.levelDist > maxDist)
                    maxDist = (int)current.levelDist;
                if (current.levelDist > maxRealDist)
                    maxRealDist = (int)current.levelDist;

                if (current.levelDist < minRealDist)
                    minRealDist = (int)current.levelDist;

                if (current.joined != null)
                    foreach (var item in current.joined)
                        st.Push(item);
            }
        }


        public void PrepareGraphNodes(Bitmap bmp, double distanceThreshold = double.MaxValue)
        {

            circleStep = 40.0 / rootNode.setStruct.Count;
            SearchKmax(rootNode);
            rootNode.gNode = new GraphNode();
            rootNode.gNode.areaLeft = 20;

            maxGraphicsY = bmp.Size.Height - 20;// - posStart - 30;// -3 * posStart;
            maxGraphicsX = bmp.Size.Width - 20;// - posStart - 30;
            if (maxDist == 0)
                throw new Exception("Dendrogram cannot be build. Wrong distances!");
            distanceStepY = ((float)maxGraphicsY) / maxDist;
            distanceStepX = ((float)maxGraphicsX) / maxDist;

            if (horizontalView)
            {
                rootNode.gNode.areaRight = bmp.Size.Width - 20;
                rootNode.gNode.x = (rootNode.gNode.areaLeft + rootNode.gNode.areaRight) / 2;
                rootNode.gNode.y = maxGraphicsY - (int)(distanceStepY * rootNode.levelDist);// + posStart;
            }
            else
            {
                rootNode.gNode.areaRight = bmp.Size.Height - 20;
                rootNode.gNode.y = maxGraphicsY - (rootNode.gNode.areaLeft + rootNode.gNode.areaRight) / 2;// +posStart+50;
                rootNode.gNode.x = maxGraphicsX - (int)(distanceStepX * rootNode.levelDist);// + posStart;
            }
            rootNode.gNode.refName = rootNode.refStructure;

            //Console.WriteLine("node=" + rootNode.gNode.y + " ldist=" + rootNode.levelDist + " step=" + distanceStepY + " aa=" + distanceStepY * rootNode.levelDist);
            //maxX = hnode.gNode.areaRight;

            // FillGraphNodes(rootNode);
            MakeGraphNodes(rootNode, distanceThreshold);
            RecalcPositions(rootNode);
            if (horizontalView)
                RecalcPositionsN(3, bmp.Width - 3);
            else
                RecalcPositionsN(3, bmp.Height - 3);
        }
        void GetTreeStack(Stack<HClusterNode> stackNodes)
        {
            if (stackNodes.Count > 0)
            {
                HClusterNode node = stackNodes.Peek();
                if (node.joined != null)
                {
                    foreach (var item in node.joined)
                    {
                        if (item.joined != null)
                        {
                            stackNodes.Push(item);
                            GetTreeStack(stackNodes);
                        }
                    }
                }
            }
        }
        void RecalcPositionsN(int start, int stop)
        {
            List<HClusterNode> nodes = rootNode.GetLeaves();
            HashSet<HClusterNode> hNodes = new HashSet<HClusterNode>();
            int counter = 0;
            for (int i = 0; i < nodes.Count; i++)
                if (nodes[i].visible)
                    counter += nodes[i].setStruct.Count;


            //float step = ((float)(stop - start) )/ nodes.Count;
            float step = ((float)(stop - start)) / (counter + 2);
            if (horizontalView)
            {

                int pos = 1;
                for (int i = 0; i < nodes.Count; i++)
                {
                    if (nodes[i].visible)
                    {
                        nodes[i].gNode.areaLeft = (int)(pos * step);
                        pos += nodes[i].setStruct.Count;
                        nodes[i].gNode.areaRight = (int)(pos * step);
                        //nodes[i].gNode.x = (int)((step * pos) + start);
                        nodes[i].gNode.x = nodes[i].gNode.areaLeft + (int)(nodes[i].gNode.areaRight - nodes[i].gNode.areaLeft) / 2;
                    }
                    if (nodes[i].gNode.x < 0)
                        Console.WriteLine();
                }
            }
            else
            {
                int pos = 1;
                for (int i = 0; i < nodes.Count; i++)
                {
                    if (nodes[i].visible)
                    {
                        nodes[i].gNode.areaLeft = (int)(pos * step);
                        pos += nodes[i].setStruct.Count;
                        nodes[i].gNode.areaRight = (int)(pos * step);
                        nodes[i].gNode.y = nodes[i].gNode.areaLeft + (int)(nodes[i].gNode.areaRight - nodes[i].gNode.areaLeft) / 2;
                    }
                    if (nodes[i].gNode.y < 0)
                        Console.WriteLine();
                }
                
            }
            
            Stack<HClusterNode> sNodes = new Stack<HClusterNode>();
            sNodes.Push(rootNode);

            GetTreeStack(sNodes);
            while (sNodes.Count > 0)
            {
                var item = sNodes.Pop();

                float sum = 0;
                int count = 0;
                foreach (var it in item.joined)
                {
                    if (it.visible)
                    {
                        if (horizontalView)
                            sum += it.gNode.x;
                        else
                            sum += it.gNode.y;
                        count++;
                    }

                }
                if(count>0)
                    sum /= count;
                if (count == 0 && item.visible)
                {
                    item.visible = false;
                    continue;
                }
                if (horizontalView)
                    item.gNode.x = (int)sum;
                else
                    item.gNode.y = (int)sum;

                if (item.gNode.y < 0 || item.gNode.x<0)
                    Console.WriteLine();
            }
        }
            
        void RecalcPositions(int start,int stop)
        {
            List<HClusterNode> nodes=rootNode.GetLeaves();
            HashSet<HClusterNode> hNodes = new HashSet<HClusterNode>();
            List<HClusterNode> toRemove = new List<HClusterNode>();
            List<HClusterNode> toAdd = new List<HClusterNode>();
            int counter = 0;
            for (int i = 0; i < nodes.Count; i++)
                if(nodes[i].visible)
                    counter += nodes[i].setStruct.Count;


            //float step = ((float)(stop - start) )/ nodes.Count;
            float step = ((float)(stop - start)) / (counter+2);
            if (horizontalView)
            {

                int pos = 1;
                for (int i = 0; i < nodes.Count; i++)
                {
                    if (nodes[i].visible)
                    {
                        nodes[i].gNode.areaLeft = (int)(pos * step);
                        pos += nodes[i].setStruct.Count;
                        nodes[i].gNode.areaRight = (int)(pos * step);
                        //nodes[i].gNode.x = (int)((step * pos) + start);
                        nodes[i].gNode.x = nodes[i].gNode.areaLeft + (int)(nodes[i].gNode.areaRight - nodes[i].gNode.areaLeft) / 2;
                    }

                }
            }
            else
            {
                int pos = 1;
                for (int i = 0; i < nodes.Count; i++)
                {
                    if (nodes[i].visible)
                    {
                        nodes[i].gNode.areaLeft = (int)(pos * step);
                        pos += nodes[i].setStruct.Count;
                        nodes[i].gNode.areaRight = (int)(pos * step);
                        nodes[i].gNode.y = nodes[i].gNode.areaLeft + (int)(nodes[i].gNode.areaRight - nodes[i].gNode.areaLeft) / 2;
                    }
                }

            }
            foreach (var item in nodes)
            {
                if(item.visible)
                    hNodes.Add(item);
            }
            while (hNodes.Count > 1)
            {
                toAdd.Clear();
                toRemove.Clear();
                foreach (var item in hNodes)
                {
                    if (!item.visible)
                        continue;
                    bool test = true;
                    if (item.parent != null)
                    {
                        foreach (var it in item.parent.joined)
                            if (!hNodes.Contains(it))
                            {
                                test = false;
                                break;
                            }
                        if (test)
                        {
                            float sum = 0;
                            foreach (var it in item.parent.joined)
                            {
                                toRemove.Add(it);
                                if (horizontalView)
                                    sum += it.gNode.x;
                                else
                                    sum += it.gNode.y;

                            }
                            sum /= item.parent.joined.Count;
                            if (horizontalView)
                                item.parent.gNode.x = (int)sum;
                            else
                                item.parent.gNode.y = (int)sum;
                            toAdd.Add(item.parent);
                        }
                    }
                }
                foreach (var item in toAdd)
                        hNodes.Add(item);
                foreach (var item in toRemove)
                    hNodes.Remove(item);
            }
        }

        void RecalcPositions(HClusterNode node)
        {
            if (node.joined != null && node.joined.Count != 0)
            {
                int joinedCount = 0;
                foreach (var item in node.joined)
                    if (item.visible)
                    {
                        RecalcPositions(item);
                        joinedCount++;
                    }
                if (joinedCount == 0)
                    return;
                if (horizontalView)
                {
                    node.gNode.x = 0;
                    foreach (var item in node.joined)
                    {
                        if(item.visible)
                            node.gNode.x += item.gNode.x;
                    }
                    node.gNode.x /= joinedCount;
                }
                else
                {
                    node.gNode.y = 0;
                    foreach (var item in node.joined)
                    {
                        if(item.visible)
                            node.gNode.y += item.gNode.y;
                    }
                    node.gNode.y /= joinedCount;
                    
                }
                if (node.gNode.x < 0 || node.gNode.y < 0)
                    Console.WriteLine();
            }


        }
        private void DrawClassColor(Graphics g, Size sizeBuff)
        {
            if (classColor != null)
            {
                colorRegions.Clear();
                SolidBrush drawBrush;
                Font drawFont = new System.Drawing.Font("Arial", 8);
                g.PageUnit = GraphicsUnit.Pixel;
                int y = 20;
                int x = 0;
                foreach (var item in classColor)
                {
                    SizeF textSize = g.MeasureString(item.Key, drawFont);
                    if (x < textSize.Width)
                        x = (int)textSize.Width;
                }
                x = sizeBuff.Width - x - 20;
                foreach (var item in classColor)
                {

                    drawBrush = new System.Drawing.SolidBrush(item.Value);
                    g.FillRectangle(drawBrush, x, y, 15, 10);
                    drawBrush = new System.Drawing.SolidBrush(System.Drawing.Color.Black);
                    g.DrawString(item.Key, drawFont, drawBrush, x + 20, y);

                    Rectangle rect = new Rectangle(x, y, 15, 10);
                    colorRegions.Add(new Region(rect), item.Key);

                    y += 25;
                    if (y > sizeBuff.Height)
                    {
                        x += 150;
                        y = 25;
                    }

                }
            }
        }
        public void DrawOnBuffer(Bitmap bmp, bool legend, int lineThickness, Color linesColor)
        {
            Graphics g;
            
            g = Graphics.FromImage(bmp);
            g.SmoothingMode = SmoothingMode.AntiAlias;
            g.TextRenderingHint = System.Drawing.Text.TextRenderingHint.AntiAlias;
         //   maxGraphicsY = bmp.Size.Height;// -2 * posStart;
         //   maxGraphicsX = bmp.Size.Width;

            maxHDist = minHDist = rootNode;
            PrepareDrawGraph(rootNode);
            DrawStack(g, bmp.Size.Height-20,bmp.Size.Width, lineThickness, legend, linesColor);
            //DrawGraph(rootNode, g,bmp.Size.Height);
            DrawDistanceAx(rootNode, g, lineThickness);
            if (legend)
                DrawClassColor(g, bmp.Size);
        }
        void DrawDistanceAx(HClusterNode localNode, Graphics gr, int lineThickness)
        {
            Pen p;
            p = new Pen(Color.Black);
            p.Width = lineThickness;
            int max,min;
            if(horizontalView)
            {
                max = maxHDist.gNode.y;
                min = minHDist.gNode.y;
            }
            else
            {
                max = maxHDist.gNode.x;
                min = minHDist.gNode.x;
            }

            int fontSize = (max - min) / 90 + 5;
            System.Drawing.Font drawFont = new System.Drawing.Font("Arial", fontSize);

            if (horizontalView)
            {
                gr.DrawLine(p, 5, max, 5, min);
                gr.DrawLine(p, 5, min, 10, min);
                gr.DrawLine(p, 5, min + max / 2, 10, min + max / 2);
            }
            else
            {
                gr.DrawLine(p,max, 5, min,5);
                gr.DrawLine(p, min,5,  min,10);
                gr.DrawLine(p, min + max / 2,5, min + max / 2,10);

            }

            System.Drawing.SolidBrush drawBrush = new System.Drawing.SolidBrush(System.Drawing.Color.Black);
            System.Drawing.StringFormat drawFormat = new System.Drawing.StringFormat();

            if (horizontalView)
            {
                if (minHDist.realDist < 1)
                    gr.DrawString(minHDist.realDist.ToString("0.00"), drawFont, drawBrush, 10, min - drawFont.SizeInPoints / 2);
                else
                    gr.DrawString(minHDist.realDist.ToString("0.0"), drawFont, drawBrush, 10, min - drawFont.SizeInPoints / 2);
            }
            else
            {
                if (minHDist.realDist < 1)
                    gr.DrawString(minHDist.realDist.ToString("0.00"), drawFont, drawBrush, min - drawFont.SizeInPoints / 2,10);
                else
                    gr.DrawString(minHDist.realDist.ToString("0.0"), drawFont, drawBrush, min - drawFont.SizeInPoints / 2,10);

            }
            double val;
            val = minHDist.realDist - maxHDist.realDist;
            if (minHDist.realDist < maxHDist.realDist)
                val = maxHDist.realDist - minHDist.realDist;
            if (horizontalView)
            {
                if (val < 1)
                    gr.DrawString((val / 2).ToString("0.00"), drawFont, drawBrush, 10, min + max / 2 - drawFont.SizeInPoints / 2);
                else
                    gr.DrawString((val / 2).ToString("0.0"), drawFont, drawBrush, 10, min + max / 2 - drawFont.SizeInPoints / 2);
            }
            else
                if (val < 1)
                    gr.DrawString((val / 2).ToString("0.00"), drawFont, drawBrush, min + max / 2 - drawFont.SizeInPoints / 2,10);
                else
                    gr.DrawString((val / 2).ToString("0.0"), drawFont, drawBrush,  min + max / 2 - drawFont.SizeInPoints / 2,10);

            gr.DrawString(measureName, drawFont, drawBrush, 5, 0);
        }
        
        void DrawStack(Graphics gr, int SizeY, int SizeX,int lineThick, bool legend, Color linesColor)
        {
            HClusterNode aux;
            while (nodesStack.Count > 0)
            {
                aux = nodesStack.Pop();
                if (aux.setStruct.Count < sizeThreshold)
                {
                    aux.visible = false;
                    continue;
                }
                if (!aux.visible)
                    continue;
                //aux.visible = true;
                if (legend)
                    aux.gNode.DrawNode(gr, lineThick);
                if (nodesStack.Count == 0)
                    continue;
                Pen p;
                if (linesColor == Color.Empty)
                {
                    if(aux.color!=Color.Black)
                        p = new Pen(aux.color,2);
                    else
                        p = new Pen(aux.color, lineThick);
                }
                else
                {
                    p = new Pen(linesColor,1);
                    // p = new Pen(Color.Black);
                    p.Width = lineThick;
                }
                if (aux.parent != null)
                {
                    if (horizontalView)
                    {
                        int x = aux.gNode.x;
                        if (aux.parent.visible && aux.visible && aux.gNode.x>=0)
                        {
                            gr.DrawLine(p, x, aux.gNode.y, x, aux.parent.gNode.y - lineThick / 2);
                            int y = aux.parent.gNode.y;
                            if(aux.parent.gNode.x>=0)
                                gr.DrawLine(p, aux.gNode.x, y, aux.parent.gNode.x, y);
                        }
                    }
                    else
                    {
                        int y = aux.gNode.y;
                        if (aux.parent.visible && y>=0 && aux.gNode.x>=0 && aux.parent.gNode.y>=0)
                        {
                            gr.DrawLine(p, aux.gNode.x, y, aux.parent.gNode.x - lineThick / 2, y);
                            int x = aux.parent.gNode.x;
                            gr.DrawLine(p, x, aux.gNode.y, x, aux.parent.gNode.y);
                        }
                    }
                    //Console.WriteLine("x=" + x + " y=" + aux.gNode.y + " x=" + x + " y=" + (aux.parent.gNode.y - lineThick / 2));
                }
                if (aux.joined == null)
                {
                    if (vecColor != null)
                        if (vecColor.ContainsKey(aux.refStructure))
                        {
                            Pen r = new Pen(classColor[vecColor[aux.refStructure]]);
                            r.Width = lineThick;
                            gr.DrawLine(r, aux.gNode.x, aux.gNode.y, aux.gNode.x, SizeY);


                        }
                    if (aux.setStruct.Count > 1 && aux.visible)
                    {
                        Rectangle rect;
                        if (horizontalView)
                        {
                                rect = new Rectangle(aux.gNode.x - 2, aux.gNode.y, 5, 5);
                                gr.DrawLine(p, aux.gNode.areaLeft, aux.gNode.y, aux.gNode.areaRight, aux.gNode.y);
                                gr.DrawLine(p, aux.gNode.areaLeft, aux.gNode.y, aux.gNode.areaLeft, aux.gNode.y + 5);
                                gr.DrawLine(p, aux.gNode.areaRight, aux.gNode.y, aux.gNode.areaRight, aux.gNode.y + 5);
                        }
                        else
                        {

                                rect = new Rectangle(aux.gNode.x, aux.gNode.y - 2, 5, 5);
                                gr.DrawLine(p, aux.gNode.x, aux.gNode.areaLeft, aux.gNode.x, aux.gNode.areaRight);
                                gr.DrawLine(p, aux.gNode.x, aux.gNode.areaLeft, aux.gNode.x + 5, aux.gNode.areaLeft);
                                gr.DrawLine(p, aux.gNode.x, aux.gNode.areaRight, aux.gNode.x + 5, aux.gNode.areaRight);
                        }
                        Color c;
                        if (clusterSize == 0)
                        {
                            if (aux.consistency < 0.3)
                                c = Color.LightBlue;
                            else
                                if (aux.consistency < 0.5)
                                c = Color.Blue;
                            else
                                    if (aux.consistency < 0.7)
                                c = Color.MediumBlue;
                            else
                                c = Color.DarkBlue;

                            gr.FillEllipse(new SolidBrush(c), rect);
                        }
                        else
                            if(aux.setStruct.Count>clusterSize)
                            {
                                c = Color.Blue;
                                gr.FillEllipse(new SolidBrush(c), rect);
                            }
                    }
                    if (aux.flagSign)
                    {
                        Pen pl = new Pen(Color.Black);
                        GraphNode gn = aux.gNode;
                        if (horizontalView)
                            gr.DrawLine(pl, gn.x, gn.y+5, gn.x, gn.y + 10);
                        else
                            gr.DrawLine(pl, gn.x+5, gn.y , gn.x + 10, gn.y);
                    }

                    string label = aux.refStructure;
                    if (label.Contains(";"))
                    {
                        string[] strAux = label.Split(';');
                        label = strAux[0];
                    }
                    if(currentLabelIndex>0)
                    {
                        if (aux.refStructure.Contains(";"))
                        {
                            string[] xx = aux.refStructure.Split(';');
                            if (labColor[xx[currentLabelIndex]] != Color.Empty)
                            {
                                p = new Pen(labColor[xx[currentLabelIndex]]);
                                gr.DrawLine(p, aux.gNode.x, aux.gNode.y + 6, aux.gNode.x, aux.gNode.y + 16);
                            }
                        }
                    }
                    if(showLabels)
                    {
                       
                        Font drawFont = new System.Drawing.Font("Arial", labelSize);
                        SolidBrush drawBrush;
                        drawBrush = new System.Drawing.SolidBrush(System.Drawing.Color.Black);
                        SizeF ss = gr.MeasureString(label, drawFont);
                        if (horizontalView)
                        {
                                gr.TranslateTransform(aux.gNode.x, SizeY - ss.Width);
                                gr.RotateTransform(90);
                                gr.DrawString(label, drawFont, drawBrush, 0, 0);
                                gr.RotateTransform(-90);
                                gr.ResetTransform();
                        
                        }
                        else
                        {                                                       
//                            gr.DrawString(label, drawFont, drawBrush, aux.gNode.x-ss.Width, aux.gNode.y-6);
                            gr.DrawString(label, drawFont, drawBrush, SizeX - ss.Width, aux.gNode.y - 6);
                        }

                    }
                }


            }
        }

        public void PrepareDrawGraph(HClusterNode startNode)
        {
            HClusterNode aux;
            Stack<HClusterNode> localStack = new Stack<HClusterNode>();

            localStack.Push(startNode);
            nodesStack.Push(startNode);
           
            while (localStack.Count > 0)
            {
                aux = localStack.Pop();

                if (aux.joined != null && aux.joined.Count > 0)
                {
                    for (int i = 0; i < aux.joined.Count; i++)
                    {
                        if (aux.joined[i].visible)
                        {
                            if (horizontalView)
                            {
                                if (aux.joined[i].gNode.y < minHDist.gNode.y)
                                    minHDist = aux.joined[i];
                                if (aux.joined[i].gNode.y > maxHDist.gNode.y)
                                    maxHDist = aux.joined[i];
                            }
                            else
                            {
                                if (aux.joined[i].gNode.x < minHDist.gNode.x)
                                    minHDist = aux.joined[i];
                                if (aux.joined[i].gNode.x > maxHDist.gNode.x)
                                    maxHDist = aux.joined[i];

                            }
                        }
                        aux.joined[i].parent = aux;                        
                        localStack.Push(aux.joined[i]);
                        nodesStack.Push(aux.joined[i]);
                    }
                }

            }
        }
        public HClusterNode CheckClick(HClusterNode localNode, int mouseX, int mouseY)
        {
            HClusterNode clickNode = null;
            if (localNode.gNode.MouseClick(mouseX, mouseY))
                return localNode;
            if (localNode.joined == null)
                return null;
            for (int i = 0; i < localNode.joined.Count; i++)
            {
                clickNode = CheckClick(localNode.joined[i], mouseX, mouseY);
                if (clickNode != null && clickNode.visible)
                    return clickNode;
            }
            return null;

        }

        public string CheckColorRegion(int mouseX, int mouseY)
        {
            foreach (var item in colorRegions)
            {
                if (item.Key.IsVisible(mouseX, mouseY))
                    return item.Value;
            }

            return null;
        }
        public void MakeAllVisible()
        {
            Queue<HClusterNode> q = new Queue<HClusterNode>();
            List<HClusterNode> leafs = new List<HClusterNode>();
            q.Enqueue(rootNode);
            while (q.Count != 0)
            {
                HClusterNode node = q.Dequeue();
                if (node.joined == null)
                    continue;

                foreach (var item in node.joined)
                        item.visible = true;

                for (int i = 0; i < node.joined.Count; i++)
                        q.Enqueue(node.joined[i]);
            }

        }
        public void MakeGraphInvisible(HClusterNode localNode,double distanceThreshold)
        {
            Queue<HClusterNode> q = new Queue<HClusterNode>();
            List<HClusterNode> leafs = new List<HClusterNode>();
            q.Enqueue(localNode);
            while (q.Count != 0)
            {
                HClusterNode node = q.Dequeue();
                if (node.joined == null)
                    continue;

                foreach (var item in node.joined)
                    if (item.realDist > distanceThreshold)
                        item.visible = false;
                    else
                        if (node.realDist > distanceThreshold && item.joined == null)
                            item.visible = false;
                        else
                            item.visible = true;

                for (int i = 0; i < node.joined.Count; i++)
                        q.Enqueue(node.joined[i]);
            }

        }
        private void MakeGraphNodes(HClusterNode localNode,double distanceThreshold)
        {
            Queue<HClusterNode> q = new Queue<HClusterNode>();

            //MakeAllVisible();
            //MakeGraphInvisible(localNode, distanceThreshold);
            
            q.Enqueue(localNode);
            while (q.Count != 0)
            {
                HClusterNode node = q.Dequeue();
                if (node.joined == null)
                    continue;
              //  List<int> rangeTab = new List<int>();
                double sum = 0;
                int joinedNum = 0;

              

                    foreach (var item in node.joined)
                {
                    // sum += item.setStruct.Count;
                    if (item.visible)
                    {
                        sum += item.setStruct.Count;
                        //sum += GetChildrenNum(item);
                        joinedNum++;
                    }
                }

                int diffArea = node.gNode.areaRight - node.gNode.areaLeft;
                List<HClusterNode> auxNodes = new List<HClusterNode>();
                IEnumerable<HClusterNode> joinedSorted = node.joined.OrderBy(x=>x.setStruct.Count);
                foreach(var nodeJ in joinedSorted)
                {
                    GraphNode graph;

                    if (!nodeJ.visible)
                        continue;

                    

                    int range=0;
                    if (viewType)
                        range=diffArea / joinedNum;
                    else
                        //rangeTab.Add(diffArea *  node.joined[i].setStruct.Count / sum);
                        //range=(int)(diffArea * GetChildrenNum(node.joined[i]) / sum);
                        range = (int)(diffArea * nodeJ.setStruct.Count / sum);


                    nodeJ.gNode = new GraphNode();
                    graph = nodeJ.gNode;
                    graph.refName = nodeJ.refStructure;
                    if (auxNodes.Count == 0)
                        graph.areaLeft = node.gNode.areaLeft;
                    else
                        graph.areaLeft = auxNodes[auxNodes.Count-1].gNode.areaRight;

                    graph.areaRight = graph.areaLeft + range;

                    auxNodes.Add(nodeJ);

                    //if (localNode.joined[i].levelDist >0)
                    if (horizontalView)
                    {
                        graph.y = maxGraphicsY - (int)(distanceStepY * nodeJ.levelDist);// + posStart;
                        graph.x = (graph.areaRight + graph.areaLeft) / 2;
                    }
                    else
                    {
                        graph.y = maxGraphicsY - (graph.areaRight + graph.areaLeft) / 2;// +posStart+50;
                        graph.x = maxGraphicsX - (int)(distanceStepX * nodeJ.levelDist);// + posStart;
                    }
                    if (graph.y < 0 || graph.x < 0)
                        Console.WriteLine();

                    if (node.joined.Count > 0/* && node.joined[i].visible*/)
                        q.Enqueue(nodeJ);
                }
            }
        }
        public void ChangeColors(HClusterNode node,Color color)
        {
            Queue<HClusterNode> q = new Queue<HClusterNode>();
            q.Enqueue(node);
            while (q.Count != 0)
            {
                HClusterNode aux = q.Dequeue();
                aux.color = color;
                if (aux.joined == null)
                    continue;
                foreach (var item in aux.joined)
                    q.Enqueue(item);
            }
        }
        public void SetColors()
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
            //buffer = null;
            
        }
        public void ColorAndLabels(string fileName)
        {

            vecColor = ClusterOutput.ReadLabelsFile(fileName);
            SetColors();
        }
        public void Save(string fileName,Bitmap buffer)
        {
            this.DrawOnBuffer(buffer, true, 1, Color.Black);
            buffer.Save(fileName);
        }
    }
}
