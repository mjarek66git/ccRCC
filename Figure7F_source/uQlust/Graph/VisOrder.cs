using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Drawing.Drawing2D;
using System.Text.RegularExpressions;
using System.IO;
using phiClustCore.Interface;

namespace Graph
{
    public partial class VisOrder : Form,IVisual
    {
        enum LeftRight
        {
            LEFT,
            RIGHT,
            BOTH
        };
        public ClosingForm closeForm;
      
        int allCount = 0;
        double step;
        string winName;    
        Point middleScreen;
        List<Queue<PointInfo>> qp;
        List<KeyValuePair<int, int>> connection = null;
        Dictionary<string, int> orderNames = null;
        int r;
        class PointInfo
        {
            public LeftRight leftRight;
            public double angl;
            public double sweep;
            public double orgSweep;
            public PointF refPoint;
            public Point point;
        };

        List<List<string>> clusters;
        List<Color> color = null;
        Dictionary<Region, int> allRegions = new Dictionary<Region, int>();
        int activeCluster = -1;

        

        public VisOrder(List<List<string>> clusters,string winName,List<Color> color)
        {
            if (clusters == null)
            {
                MessageBox.Show("No clusterisation defined!");
                return;
            }
            InitializeComponent();
            
            this.color = color;
            this.winName = winName;
            this.clusters = clusters;
            foreach (var item in clusters)
                allCount += item.Count;
            step = 360.0 / allCount;
            middleScreen = new Point(panel1.Size.Width / 2, panel1.Size.Height / 2);
            qp = new List<Queue<PointInfo>>(clusters.Count);
            for (int i = 0; i < clusters.Count; i++)
                qp.Add(new Queue<PointInfo>());


            CheckOrder();
       
        }
        public override string ToString()
        {
            return "Order";
        }
        public void ToFront()
        {
            this.BringToFront();
        }
        void CheckOrder()
        {
            List<KeyValuePair<int, int>> data = new List<KeyValuePair<int, int>>();
            for(int i=0;i<clusters.Count;i++)
            {
                foreach (var item in clusters[i])
                {
                    if (orderNames == null)
                    {
                        string[] dd = Regex.Split(item, @"\D+");
                        int vv = 1000;
                        int res = 1;
                        foreach (var str in dd)
                        {
                            if (str.Length > 0)
                                res *= vv;
                        }
                        vv = res;
                        res = 0;
                        foreach (var str in dd)
                        {
                            if (str.Length > 0)
                            {
                                res += Convert.ToInt32(str) * vv;
                                vv /= 1000;
                            }
                        }
                        data.Add(new KeyValuePair<int, int>(res, i));
                    }
                    else
                        if(orderNames.ContainsKey(item))
                        {
                            data.Add(new KeyValuePair<int, int>(orderNames[item], i));
                        }
                        else
                        {
                            MessageBox.Show("There is no order for: " + item);
                        }
                }
            }
            data.Sort((firstPair, nextPair) =>
            {
                return nextPair.Key.CompareTo(firstPair.Key);
            });
            connection = new List<KeyValuePair<int, int>>();
            for(int i=0;i<data.Count-1;i++)
            {
                if (data[i].Value != data[i + 1].Value)
                    connection.Add(new KeyValuePair<int, int>(data[i].Value, data[i + 1].Value));
            }
        }
        void EnqueNode(PointInfo info, int i, LeftRight type)
        {
            PointInfo aux;
            aux = new PointInfo();
            aux.leftRight = type;
            aux.angl = info.angl;
            aux.refPoint = info.refPoint;
            aux.orgSweep = info.orgSweep;
            if (type == LeftRight.LEFT)
                aux.sweep = info.sweep - 0.1;
            else
                aux.sweep = info.sweep + 0.1;
            if (aux.sweep > 2 * aux.orgSweep)
                aux.sweep -= 0.1;
            if (aux.sweep < 0)
                aux.sweep += 0.1;
            aux.point = CalcPoint(aux.refPoint, aux.angl, aux.sweep, r - 20);
            qp[i].Enqueue(aux);
        }
        Point CalcPoint(PointF refP, double angle, double sweep, int r)
        {
            double tg = Math.Tan(Math.PI / 180 * (angle + sweep / 2));
            double xx = r / Math.Sqrt(1 + tg * tg);
            PointF rem = new PointF();
            rem.X = (refP.X - middleScreen.X);
            rem.Y = (refP.Y - +middleScreen.Y);
            double skalarny = rem.X * xx + rem.Y * xx * tg;
            skalarny /= Math.Sqrt(rem.X * rem.X + rem.Y * rem.Y) * Math.Sqrt(xx * xx + xx * xx * tg * tg);
            int skalar = Convert.ToInt32(100 * Math.Cos(Math.PI / 180 * sweep / 2));
            int cos = Convert.ToInt32(100 * skalarny);
            if (cos != skalar)
                xx = -xx;

            Point res = new Point();
            res.X = Convert.ToInt32(xx + middleScreen.X);
            res.Y = Convert.ToInt32(xx * tg + middleScreen.Y);

            return res;
        }

        Point GetQPoint(int i)
        {
            PointInfo info1;
            info1 = qp[i].Dequeue();
            if (info1.leftRight == LeftRight.BOTH)
            {
                EnqueNode(info1, i, LeftRight.LEFT);
                EnqueNode(info1, i, LeftRight.RIGHT);
            }
            else
                if (info1.leftRight == LeftRight.LEFT)
                    EnqueNode(info1, i, LeftRight.LEFT);

                else
                    EnqueNode(info1, i, LeftRight.RIGHT);

            return info1.point;

        }

        void DrawLines(Graphics g, Point p1, Point p2)
        {
            Pen p = new Pen(Color.Red, 1);
            Point[] all = { p1, middleScreen, p2 };
            g.DrawCurve(p, all);


        }

        private void panel1_Paint(object sender, PaintEventArgs e)
        {
            List<Point> points = new List<Point>();
            Graphics g = e.Graphics;
            g.Clear(this.BackColor);
            double angl = 0;
            r = Math.Min(panel1.Size.Height, panel1.Size.Width);
            r /= 2;
            r -= 5;
            Rectangle rec = new Rectangle(0, 0, r, r);
            Rectangle rec1 = new Rectangle(10, 10, r - 20, r - 20);
            GraphicsPath gPath;
            allRegions.Clear();
            for (int i = 0; i < clusters.Count;i++ )
            {
                double sweep = step * clusters[i].Count;
                gPath = PreparePath(r, r - 20, angl, sweep);
                qp[i].Clear();
                if (gPath != null)
                {
                    allRegions.Add(new Region(gPath), i);
                    Pen p;
                    Brush aBrush = (Brush)Brushes.Green;
                    p = new Pen(Color.Black, 2);
                    if (color == null)
                    {
                        if (activeCluster != i)
                            aBrush = (Brush)Brushes.BlueViolet;
                        else
                            aBrush = (Brush)Brushes.Yellow;
                    }
                    else
                    {
                        aBrush = new SolidBrush(color[i]);
                    }
                    g.FillPath(aBrush, gPath);
                    g.DrawPath(p, gPath);
                    //foreach (var po in gPath.PathData.Points) 
                    Point mid = CalcPoint(gPath.PathPoints[0], angl, sweep, r - 20);
                    g.FillRectangle(aBrush, mid.X, mid.Y, 4, 4);
                    points.Add(mid);
                    PointInfo aux = new PointInfo();
                    aux.leftRight = LeftRight.BOTH;
                    aux.angl = angl;
                    aux.sweep = sweep;
                    aux.orgSweep = aux.sweep;
                    aux.point = mid;
                    aux.refPoint = gPath.PathPoints[0];
                    qp[i].Enqueue(aux);

                }

                angl += sweep;
            }
            Point point1;
            Point point2;
            foreach(var item in connection)
            {
                if (activeCluster != -1 )                    
                {
                    if (activeCluster == item.Key || activeCluster == item.Value)
                    {
                        point1 = GetQPoint(item.Key);
                        point2 = GetQPoint(item.Value);
                        DrawLines(g, point1, point2);
                    }
                }
                else
                {
                    point1 = GetQPoint(item.Key);
                    point2 = GetQPoint(item.Value);
                    DrawLines(g, point1, point2);

                }

                
            }
        }
        GraphicsPath PreparePath(double distEnd, double pDist, double angleStart, double sweepAngle)
        {
            GraphicsPath gPath = new GraphicsPath();
            GraphicsPath auxPath1 = new GraphicsPath();
            GraphicsPath auxPath2 = new GraphicsPath();
            int pDistInt = (int)(pDist);
            int distEndInt = (int)(distEnd);
            int angleStartInt = (int)Math.Round(angleStart);
            int sweepAngleInt = (int)Math.Round(sweepAngle);
            int midX = middleScreen.X - pDistInt;
            int midY = middleScreen.Y - pDistInt;
            int midXE = middleScreen.X - distEndInt;
            int midYE = middleScreen.Y - distEndInt;
            auxPath1.AddArc(midX, midY, (int)(2 * pDist), (int)(2 * pDist), (float)angleStart, (float)sweepAngle);
            auxPath2.AddArc(midXE, midYE, (int)(2 * distEnd), (int)(2 * distEnd), (float)(angleStart + sweepAngle), -(float)sweepAngle);
            if (auxPath1.PathPoints == null || auxPath2.PathPoints == null || auxPath1.PathPoints.Length <= 1 || auxPath2.PathPoints.Length <= 1)
                return null;

            gPath.AddArc(midX, midY, (int)(2 * pDist), (int)(2 * pDist), (float)angleStart, (float)sweepAngle);
            gPath.AddLine(auxPath1.PathPoints[auxPath1.PathPoints.Length - 1].X, auxPath1.PathPoints[auxPath1.PathPoints.Length - 1].Y, auxPath2.PathPoints[0].X, auxPath2.PathPoints[0].Y);
            gPath.AddArc(midXE, midYE, (int)(2 * distEnd), (int)(2 * distEnd), (float)(angleStart + sweepAngle), -(float)sweepAngle);
            gPath.AddLine(auxPath2.PathPoints[auxPath2.PathPoints.Length - 1].X, auxPath2.PathPoints[auxPath2.PathPoints.Length - 1].Y, auxPath1.PathPoints[0].X, auxPath1.PathPoints[0].Y);


            return gPath;

        }

        private void VisOrder_Resize(object sender, EventArgs e)
        {
            middleScreen = new Point(panel1.Size.Width / 2, panel1.Size.Height / 2);
            this.Invalidate();
        }
        int CheckRegion(int mouseX, int mouseY)
        {
            foreach (var item in allRegions)
            {
                if (item.Key.IsVisible(mouseX, mouseY))
                    return item.Value;
            }

            return -1;
        }

        private void panel1_MouseClick(object sender, MouseEventArgs e)
        {
            if (e.Button == MouseButtons.Left)
            {
                activeCluster = CheckRegion(e.X, e.Y);
                this.Invalidate();
                this.Refresh();
            }
            else
                if (e.Button == MouseButtons.Right)
                {
                    activeCluster = -1;
                    this.Invalidate();
                    this.Refresh();
                }
        }

        private void VisOrder_FormClosed(object sender, FormClosedEventArgs e)
        {
            if (closeForm != null)
                closeForm(winName);
        }
        private void ReadOrder()
        {
            DialogResult res;

            openFileDialog1 = new OpenFileDialog();
            openFileDialog1.Title = "File name with the name of labels for each struture";
            res = openFileDialog1.ShowDialog();
            if (res == DialogResult.OK)
            {
                orderNames = new Dictionary<string, int>();
                StreamReader st = new StreamReader(openFileDialog1.FileName);
                string line;
                while ((line = st.ReadLine()) != null)
                {
                    string[] aux = line.Split(' ');

                    if (aux.Length == 2)
                    {
                        if (aux[0].Contains(Path.DirectorySeparatorChar))
                            aux[0] = Path.GetFileName(aux[0]);

                        orderNames.Add(aux[0], Convert.ToInt32(aux[1]));                          
                    }

                }

                st.Close();

                CheckOrder();
                this.Invalidate();
                this.Refresh();
            }
        }
        private void toolStripButton1_Click(object sender, EventArgs e)
        {
            ReadOrder();
        }

        private void toolStripButton2_Click(object sender, EventArgs e)
        {
            DialogResult res;
            saveFileDialog1 = new SaveFileDialog();
            saveFileDialog1.Filter = "Png files (*.png)|*.png";
            saveFileDialog1.DefaultExt = "*.png";
            res=saveFileDialog1.ShowDialog();

            if(res==DialogResult.OK)
                using (var bmp = new Bitmap(panel1.Width, panel1.Height))
                {
                    panel1.DrawToBitmap(bmp, new Rectangle(0, 0, bmp.Width, bmp.Height));
                    bmp.MakeTransparent(panel1.BackColor);
                    bmp.Save(saveFileDialog1.FileName);
                }
        }


    }
}
