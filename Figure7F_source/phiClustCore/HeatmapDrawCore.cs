using phiClustCore;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace phiClustCore
{
    public class ColorRect
    {
        public float sX;
        public float sY;
        public float width;
        public float height;
        public Color c;
        public string text;
    };
    public class ColorDist
    {
        public static Color[] colorTable;
        static string[] strColor = new string[]{ "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
        "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
        "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
        "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
        "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
        "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
        "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
        "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",

        "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
        "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
        "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
        "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
        "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
        "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
        "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
        "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58",

        "#7A7BFF", "#D68E01", "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
        "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", "#0AA3F7", "#E98176",
        "#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", "#02684E", "#962B75", "#8D8546", "#9695C5",
        "#E773CE", "#D86A78", "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
        "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", "#643127", "#513A01",
        "#6B94AA", "#51A058", "#A45B02", "#1D1702", "#E20027", "#E7AB63", "#4C6001", "#9C6966",
        "#64547B", "#97979E", "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
        "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", "#C6DC99", "#203B3C",

        "#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", "#CCAA35", "#374527", "#8BB400", "#797868",
        "#C6005A", "#3B000A", "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
        "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", "#71BB8C", "#00B433",
        "#789EC9", "#6D80BA", "#953F00", "#5EFF03", "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F",
        "#003109", "#0060CD", "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
        "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", "#EA8B66", "#958A9F",
        "#BDC9D2", "#9FA064", "#BE4700", "#658188", "#83A485", "#453C23", "#47675D", "#3A3F00",
        "#061203", "#DFFB71", "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",

        "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", "#866097", "#365D25",
        "#252F99", "#00CCFF", "#674E60", "#FC009C", "#92896B"};
        static ColorDist()
        {
            colorTable = new Color[strColor.Length];
            for (int i = 0; i < colorTable.Length; i++)
                colorTable[i] = ColorTranslator.FromHtml(strColor[i]);
        }
    };
    public class HeatmapDrawCore
    {
        public Dictionary<Region, KeyValuePair<bool, string>> regionBarColor { get; set; } = new Dictionary<Region, KeyValuePair<bool, string>>();

        public List<HClusterNode> upperLeaves { get; set; }
        public List<HClusterNode> leftLeaves { get; set; }
        Random rand = new Random();
        public Dictionary<string, Tuple<Dictionary<string, Tuple<Color, string>>, List<ColorRect>>> leftDescription { get; set; }
//        HClusterNode colorNodeUpper = null;
//        HClusterNode colorNodeLeft = null;

        public bool showLabels { get; set; } = false;
        public bool showSelectedLabels { get; set; } = false;
        Font drawFont = new Font("Arial", 10);
        SolidBrush brush = new SolidBrush(Color.Red);
        public HeatMapDraw draw { get; set; } = null;
        Dictionary<string, Tuple<string, int, int>> selectedGeneNames = new Dictionary<string, Tuple<string, int, int>>();
        public HashSet<string> geneSelection { get; set; }
        
        //bool cutVerticalLine = false;
        //int mouseX = 0;
        //int mouseY = 0;
        //int labelSize = 10;

        //bool interactiveLabelColumn = false;
        //bool interactiveLabelRow = false;
        public bool avrRowAnnotations { get; set; } = false;
        public bool horizontalLines { get; set; } = false;

        public double verticalCuttDistance { get; set; }
        public double horizontalCuttDistance { get; set; }

        public HeatmapDrawCore()
        {

        }
        public HeatmapDrawCore(ClusterOutput output,int resX,int resY)
        {
            this.upperLeaves = output.nodes[0].GetLeaves();
            this.leftLeaves = output.nodes[1].GetLeaves();
            draw = new HeatMapDraw(new Bitmap(resX, resY),output.nodes[0],output.nodes[1], null, output);
        }
        void FinSelectedGenesPositions()
        {
            upperLeaves = draw.auxUpper.GetLeaves();
            List<string> keys = new List<string>(selectedGeneNames.Keys);
            foreach (var item in keys)
            {
                string name = selectedGeneNames[item].Item1;
                selectedGeneNames[item] = new Tuple<string, int, int>(name, -1, -1);
            }

            foreach (var item in upperLeaves)
            {

                if (selectedGeneNames.ContainsKey(item.refStructure))
                {
                    string name = selectedGeneNames[item.refStructure].Item1;
                    selectedGeneNames[item.refStructure] = new Tuple<string, int, int>(name, item.gNode.x, item.gNode.y);
                }
            }
        }
        public void ReadNoteGenes(string fileName)
        {
            selectedGeneNames = new Dictionary<string, Tuple<string, int, int>>();
            //selectedGeneNames.Clear();
            StreamReader wr = new StreamReader(fileName);
            string line = wr.ReadLine();
            while (line != null)
            {
                string[] aux = line.Split(';');
                selectedGeneNames.Add(aux[0], new Tuple<string, int, int>(aux[1], -1, -1));
                line = wr.ReadLine();
            }
            wr.Close();

            showSelectedLabels = true;

        }
        void PrepareBarRegions(Graphics g)
        {
            int x = 50, y = 0;
            bool regionTest = false;

            if (draw.upper.labColor == null || draw.upper.labColor.Keys.Count == 0)
                return;

            if (regionBarColor.Count > 0)
                regionTest = true;

            List<string> labKeys = new List<string>(draw.upper.labColor.Keys);
            Font drawFont = new System.Drawing.Font("Arial", 8);
            foreach (var item in labKeys)
            {
                if (!regionTest)
                {
                    Region reg = new Region(new Rectangle(x, y, 15, 10));
                    regionBarColor.Add(reg, new KeyValuePair<bool, string>(false, item));
                    SizeF textSize = g.MeasureString(item, drawFont);
                    x += 25 + (int)textSize.Width;
                    if (x > g.VisibleClipBounds.Width)
                    {
                        y += 150;
                        x = 25;
                    }
                }
            }
            foreach (var regItem in regionBarColor)
            {
                if (regItem.Value.Key)
                    draw.upper.labColor[regItem.Value.Value] = Color.Empty;
            }

        }

        public Bitmap PaintUpper(Bitmap bit)
        {
            Graphics g = Graphics.FromImage(bit);
            Brush drawBrush = new System.Drawing.SolidBrush(Color.Black);
            g.Clear(Color.White);
            //selectedGeneNames = null;
            if (selectedGeneNames != null && showSelectedLabels)
            {

                FinSelectedGenesPositions();
                drawFont = new System.Drawing.Font("Arial", 22, FontStyle.Bold);
                Pen p = new Pen(Color.Red, 3);
                foreach (var item in selectedGeneNames)
                {
                    if (item.Value.Item2 >= 0)
                    {
                        SizeF ss = g.MeasureString(item.Value.Item1, drawFont);
                        g.TranslateTransform(item.Value.Item2, bit.Height);
                        g.RotateTransform(-45);

                        g.DrawString(item.Value.Item1, drawFont, drawBrush, 10, -15);
                        g.ResetTransform();
                        g.DrawEllipse(p, item.Value.Item2, bit.Height - 5, 3, 3);

                    }
                }

                return bit;
            }
            if (showLabels == true)
            {

                upperLeaves = draw.auxUpper.GetLeaves();
                drawFont = new System.Drawing.Font("Arial", draw.upper.labelSize);
                foreach (var item in upperLeaves)
                {
                    SizeF ss = g.MeasureString(item.refStructure, drawFont);
                    g.TranslateTransform(item.gNode.x, bit.Height);
                    g.RotateTransform(-45);

                    g.DrawString(item.refStructure, drawFont, drawBrush, 5, -5);
                    g.ResetTransform();
                }
            }

            if (!showLabels && !showSelectedLabels)
            {                
                PrepareBarRegions(g);
                //e.Graphics.Clear(pictureBox2.BackColor);

                draw.upper.DrawOnBuffer(bit, false, 3, Color.Empty);

                //                e.Graphics.DrawImage(draw.upperBitMap, 0, 0);

                //pictureBox1.Image = draw.upperBitMap;
                if (draw.upper.labColor != null && draw.upper.labColor.Count > 0)
                {
                    Font drawFont = new System.Drawing.Font("Arial", 8);
                    foreach (var item in regionBarColor)
                    {
                        RectangleF rec = item.Key.GetBounds(g);
                        drawBrush = new System.Drawing.SolidBrush(draw.upper.labColor[item.Value.Value]);
                        g.FillRectangle(drawBrush, rec.X, rec.Y, rec.Width, rec.Height);

                        SizeF textSize = g.MeasureString(item.Value.Value, drawFont);
                        drawBrush = new System.Drawing.SolidBrush(Color.Black);

                        g.DrawString(item.Value.Value, drawFont, drawBrush, rec.X + 20, rec.Y);
                        if (item.Value.Key)
                        {
                            Pen p = new Pen(Color.Black);
                            g.DrawLine(p, rec.X, rec.Y, rec.X + rec.Width, rec.Y + rec.Height);
                            g.DrawLine(p, rec.X + rec.Width, rec.Y, rec.X, rec.Y + rec.Height);
                        }
                    }

                }


            }


            return bit;
        }
        public Bitmap PaintLeft(Bitmap bit)
        {
            draw.left.posStart = 5;
            Graphics g = Graphics.FromImage(bit);
            g.Clear(Color.White);
            if (!showSelectedLabels)
                draw.left.DrawOnBuffer(bit, false, 3, Color.Empty);
            else
            {
                if (leftDescription != null && leftDescription.Count > 0)
                {
                    drawFont = new System.Drawing.Font("Arial", 20, FontStyle.Bold);
                    Pen p = new Pen(Color.Red, 3);
                    List<string> keys = new List<string>(leftDescription.Keys);
                    Brush drawBrush = new SolidBrush(Color.Black);
                    var format = new StringFormat() { Alignment = StringAlignment.Far, LineAlignment = StringAlignment.Center };

                    foreach (var item in leftLeaves)
                    {
                        if (item.setStruct.Count < 200)
                            continue;
                        if (leftDescription[keys[0]].Item1.ContainsKey(item.refStructure))
                        {
                            SizeF ss = g.MeasureString(leftDescription[keys[0]].Item1[item.refStructure].Item2, drawFont);
                            var rect = new RectangleF(0, item.gNode.y - 10, bit.Width, 22);
                            g.DrawString(leftDescription[keys[0]].Item1[item.refStructure].Item2, drawFont, drawBrush, rect, format);
                            p.Color = leftDescription[keys[0]].Item1[item.refStructure].Item1;
                            //g.DrawEllipse(p, item.gNode.x+2, item.gNode.y, 3, 3);

                        }
                    }
                }
            }


            return bit;
        }
        public Bitmap PaintHeatMap(Bitmap bit)
        {
            Graphics g = Graphics.FromImage(bit);
            g.Clear(Color.White);
            draw.DrawHeatMapN(bit);


            if (horizontalLines)
            {
                SolidBrush b = new SolidBrush(Color.White);
                Pen p = new Pen(b, 2);
                int max = int.MinValue;
                int min = int.MaxValue;
                foreach (var item in upperLeaves)
                {
                    if (item.gNode.x < min)
                        min = item.gNode.areaLeft;
                    if (item.gNode.x > max)
                        max = item.gNode.areaRight;
                }


                foreach (var item in leftLeaves)
                {
                    g.DrawLine(p, min, item.gNode.areaLeft, max, item.gNode.areaLeft);
                }
            }
            if (geneSelection != null && geneSelection.Count > 0)
            {
                Color newColor = Color.FromArgb(50, 255, 0, 0);
                Brush br = new SolidBrush(newColor);
                foreach (var item in upperLeaves)
                {
                    if (geneSelection.Contains(item.refStructure))
                    {
                        Rectangle r = new Rectangle(item.gNode.areaLeft, 0, item.gNode.areaRight - item.gNode.areaLeft, bit.Height);
                        g.FillRectangle(br, r);
                    }
                }
            }
            return bit;
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


        public void PaintLeftDescription(int widthB)
        {
            if (leftDescription != null && leftDescription.Count > 0)
            {
                int width = widthB / leftDescription.Count;
                Tuple<Dictionary<string, Tuple<Color, string>>, List<ColorRect>> element;
                List<string> keys = new List<string>(leftDescription.Keys);
                for (int i = 0; i < keys.Count; i++)
                {
                    int pos = i * width;
                    Dictionary<string, Tuple<Color, string>> x = leftDescription[keys[i]].Item1;
                    if (avrRowAnnotations)
                        element = CalcLeftBarsClusterAvr(x, pos, width);
                    else
                        element = CalcLeftBars(x, pos, width);

                    leftDescription[keys[i]] = element;

                }
            }

        }
        Tuple<Dictionary<string, Tuple<Color, string>>, List<ColorRect>> CalcLeftBarsClusterAvr(Dictionary<string, Tuple<Color, string>> colorMap, int pos = 0, int width = 20)
        {
            List<ColorRect> r = new List<ColorRect>();
            foreach (var item in leftLeaves)
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
            foreach (var item in leftLeaves)
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
                        foreach (var am in mappings)
                            colorBrihtness += am.Value;
                        colorBrihtness = gkeyPair.Value / colorBrihtness;
                        if (mappings.Count > 0)
                        {
                            mappings.Sort((x, y) => x.Value.CompareTo(y.Value));
                            for (int i = 0; i < temp.Count; i++)
                                if ((int)Math.Floor(temp[i].sY) == it)
                                {
                                    Color col = mappings[mappings.Count - 1].Key;
                                    int re = col.R - System.Convert.ToInt32(colorBrihtness / (double)100 * col.R);
                                    int gr = col.G - System.Convert.ToInt32(colorBrihtness / (double)100 * col.G);
                                    int bl = col.B - System.Convert.ToInt32(colorBrihtness / (double)100 * col.B);
                                    //temp[i].c = mappings[mappings.Count - 1].Key;
                                    temp[i].c = Color.FromArgb(255, re, gr, bl);
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
        Dictionary<string, Tuple<Color, string>> ReadClusterFile(string fileName)
        {
            StreamReader r = new StreamReader(fileName);
            Dictionary<string, Tuple<Color, string>> clDic = new Dictionary<string, Tuple<Color, string>>();
            Dictionary<string, string> labels = new Dictionary<string, string>();
            HashSet<string> distinc = new HashSet<string>();
            Dictionary<string, Color> dicColor = new Dictionary<string, Color>();
            string line = r.ReadLine();
            bool colorDef = false;
            bool dataDef = false;
            while (line != null)
            {
                string[] tmp = line.Split(new char[] { ' ', '\t' });
                if (colorDef)
                {
                    if (tmp.Length == 2)
                        dicColor.Add(tmp[0], ColorTranslator.FromHtml(tmp[1]));
                }
                if (dataDef)
                {
                    if (tmp.Length == 2)
                        labels.Add(tmp[0], tmp[1]);
                }
                if (line.Contains("[COLOR]"))
                {
                    colorDef = true;
                    dataDef = false;
                }
                if (line.Contains("[DATA]"))
                {
                    colorDef = false;
                    dataDef = true;
                }

                line = r.ReadLine();
            }
            r.Close();

            foreach (var item in labels)
                clDic.Add(item.Key, new Tuple<Color, string>(dicColor[item.Value], item.Value));

            return clDic;
        }

        public Tuple<int,string> ReadLeftAnnotation(string fileName)
        {
            if (leftDescription == null)
                 leftDescription = new Dictionary<string, Tuple<Dictionary<string, Tuple<Color, string>>, List<ColorRect>>>();
            var colorMap = ReadClusterFile(fileName);
            int pos = 1;
            if (leftDescription.Count > 0)
                pos = 35;

            foreach (var item in leftDescription)
            {
                pos += GetDistinctColors(item.Value).Count * 20;
            }
            string NameT = "Name" + rand.Next();
            leftDescription.Add(NameT, new Tuple<Dictionary<string, Tuple<Color, string>>, List<ColorRect>>(colorMap, null));
            return new Tuple<int, string>(pos,NameT);
        }
        public Bitmap PaintLeftAnnotation(Bitmap bit)
        {
            if (leftDescription != null)
            {
                Graphics g = Graphics.FromImage(bit);
                g.Clear(Color.White);
                Pen p = new Pen(Color.White, 2);
                foreach (var item in leftDescription)
                {
                    float lastX = 0;
                    float lastY = 0;
                    foreach (var it in item.Value.Item2)
                    {
                        Brush b = new SolidBrush(it.c);
                        g.FillRectangle(b, it.sX, it.sY, it.width, it.height);
                        lastX = it.sX;
                        lastY = it.sY + it.height;
                    }
                    if (leftDescription.Count > 1)
                    {
                        g.DrawLine(p, lastX, 0, lastX, lastY);
                    }
                }
            }
            return bit;

        }
        public Bitmap PaintRightAnnotation(Bitmap bit)
        {
            Dictionary<string, Color> x = new Dictionary<string, Color>();
            if (leftDescription != null)
            {
                Graphics g = Graphics.FromImage(bit);
                g.Clear(Color.White);
                int posy = 15;
                foreach (var item in leftDescription)
                {
                    x = GetDistinctColors(item.Value);

                    int posx = 0;
                    Brush stringB = new SolidBrush(Color.Black);
                    List<string> keys = new List<string>(x.Keys);

                    if (showSelectedLabels)
                    {
                        if (keys.Count != 3)
                            continue;
                    }


                    try
                    {
                        if (keys.Count == 3)
                        {
                            keys = new List<string>();
                            keys.Add("Stage1");
                            keys.Add("Stage3");
                            keys.Add("Metastatic");
                        }
                        else
                            keys.Sort((a, b) => Convert.ToInt32(a).CompareTo(Convert.ToInt32(b)));
                    }
                    catch (Exception ex)
                    {
                        keys.Sort();
                    }
                    foreach (var it in keys)
                    {
                        Brush b = new SolidBrush(x[it]);

                        g.FillRectangle(b, posx, posy, 10, 10);
                        g.DrawString(it, drawFont, stringB, posx + 10, posy - 3);
                        posy += 20;
                    }
                    //posx += 21;
                    posy += 40;
                }
            }
            return bit;

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
        public Bitmap PaintLegend(Bitmap bit)
        {
            Graphics g = Graphics.FromImage(bit);
            g.Clear(Color.White);
            SolidBrush b = new SolidBrush(Color.Black);
            int xPos, yPos;
            xPos = 5; yPos = 5;
            System.Drawing.Font drawFont = new System.Drawing.Font("Arial", 10);
            System.Drawing.SolidBrush drawBrush = new System.Drawing.SolidBrush(System.Drawing.Color.Black);
            System.Drawing.StringFormat drawFormat = new System.Drawing.StringFormat();

            List<double> ordered = new List<double>(draw.profilesColorMap.Keys);
            ordered.Sort();
            foreach (var item in ordered)
            {
                b.Color = draw.profilesColorMap[item];
                g.FillRectangle(b, xPos, yPos, 15, 10);
                //e.Graphics.DrawString(item.ToString(), drawFont, drawBrush, xPos+25,yPos-3);
                g.DrawString(item.ToString(), drawFont, drawBrush, xPos + 25, yPos - 3);
                yPos += 25;
                if (yPos > bit.Height)
                {
                    yPos = 5;
                    xPos += 40;
                }
            }
            return bit;
        }

        public void Save(string fileName, int widthP, int heightP)
        {
            double lW = 0.11;
            double aLW = 0.03;
            double aRW = 0.06;
            double wHM = 0.8;
            double hU = 0.2;
            double hHM = 0.8;
            Bitmap bmp = new Bitmap(widthP, heightP);
            byte[] resultingBytes;
            using (var stream = new MemoryStream())
            using (var graphicsFromSizedImage = Graphics.FromImage(bmp))
            using (var metafile = new Metafile(stream, graphicsFromSizedImage.GetHdc()))
            {
                Bitmap upperBitMap = new Bitmap((int)(bmp.Width*wHM), (int)(bmp.Height*hU));
                Bitmap leftBitMap = new Bitmap((int)(bmp.Width*lW), (int)(bmp.Height * hHM));
                Bitmap heatMap = new Bitmap(upperBitMap.Width, leftBitMap.Height);
                Bitmap leftAnnotation = new Bitmap((int)(bmp.Width * aLW), leftBitMap.Height);
                Bitmap rightAnnotation = new Bitmap((int)(bmp.Width * aRW), leftBitMap.Height);
                Bitmap legend = new Bitmap(leftBitMap.Width, upperBitMap.Height);
                draw.left.PrepareGraphNodes(leftBitMap, verticalCuttDistance);
                draw.upper.PrepareGraphNodes(upperBitMap, horizontalCuttDistance);

                PaintLeft(leftBitMap);
                PaintUpper(upperBitMap);
                PaintHeatMap(heatMap);
                PaintLeftDescription(leftAnnotation.Width);
                PaintLeftAnnotation(leftAnnotation);
                PaintRightAnnotation(rightAnnotation);
                PaintLegend(legend);
                Graphics g;

                using (fileName.Contains(".emf") ? g = Graphics.FromImage(metafile) : g = Graphics.FromImage(bmp))
                {
                    g.Clear(Color.White);

                    g.DrawImage(upperBitMap, leftBitMap.Width+leftAnnotation.Width, 0);
                    g.DrawImage(leftBitMap, 0, upperBitMap.Height);
                    g.DrawImage(leftAnnotation,leftBitMap.Width, upperBitMap.Height);
                    g.DrawImage(heatMap, leftBitMap.Width + leftAnnotation.Width,upperBitMap.Height);
                    //g.DrawImage(rightAnnotation, width + annotWidth + heatMap.Width, height);

                }

                if (fileName.Contains(".emf"))
                {
                    graphicsFromSizedImage.ReleaseHdc();
                    resultingBytes = stream.ToArray();
                    File.WriteAllBytes(fileName, resultingBytes);
                }
                else
                    bmp.Save(fileName, ImageFormat.Jpeg);
            }
            draw.left.PrepareGraphNodes(draw.leftBitMap, verticalCuttDistance);
            draw.upper.PrepareGraphNodes(draw.upperBitMap, horizontalCuttDistance);

        }
        public void ResizeHeatMap(Size bUpper,Size bLeft,Size b)
        {
            draw.upperBitMap = new Bitmap(bUpper.Width, bUpper.Height);
            draw.leftBitMap = new Bitmap(bLeft.Width, bLeft.Height);
            draw.heatMap = new Bitmap(bUpper.Width, bLeft.Height);
            draw.left.PrepareGraphNodes(draw.leftBitMap, verticalCuttDistance);
            draw.upper.PrepareGraphNodes(draw.upperBitMap, horizontalCuttDistance);
            PaintLeftDescription(b.Width);
        }
        public void ClearDistance()
        {
            horizontalCuttDistance = double.MaxValue;
            verticalCuttDistance = double.MaxValue;
            draw.auxUpper.ClearColors(Color.Black);
            draw.auxLeft.ClearColors(Color.Black);
            draw.upper.MakeAllVisible();
            draw.left.MakeAllVisible();
            draw.upper.PrepareGraphNodes(draw.upperBitMap, horizontalCuttDistance);
            draw.left.PrepareGraphNodes(draw.leftBitMap, verticalCuttDistance);
        }
        public void Swap()
        {
            upperLeaves = draw.upperNode.GetLeaves();
            draw.auxUpper = draw.upperNode;

        }
        public void LongLeavesClick()
        {
            draw.upper.MakeAllVisible();
            upperLeaves = draw.auxUpper.GetLeaves();
            for (int i = 0; i < upperLeaves.Count; i++)
            {
                if (upperLeaves[i].parent.realDist > 0.5)
                    upperLeaves[i].visible = false;
                if (upperLeaves[i].parent.parent != null)
                {
                    int l = upperLeaves[i].parent.GetLeaves().Count;
                    if (l < 3 && upperLeaves[i].parent.parent.realDist > 0.5 && upperLeaves[i].parent.realDist > 0.2)
                        upperLeaves[i].visible = false;
                }
            }
            draw.upper.PrepareGraphNodes(draw.upperBitMap, horizontalCuttDistance);
        }

    }
}
