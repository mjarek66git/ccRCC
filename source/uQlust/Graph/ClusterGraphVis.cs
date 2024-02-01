using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using Graph;
using System.Drawing;
using phiClustCore;
using phiClustCore.Interface;
namespace Graph
{

    public class ClusterGraphVis : ClusterVis
    {
        string Name;
        Random r = new Random();
        int randomV;
        IVisual active = null;
        Dictionary<string,ClusterOutput> lOut;
        public static List<string> hNodeOptions = new List<string>{"Dendrogram", "Sunburst chart"};
        public static List<string> clusterOptions = new List<string> { "Text List",  "Order Visual"};
        public ClusterGraphVis() { randomV = r.Next(); }
        public ClusterGraphVis(ClusterOutput output, string name, Dictionary<string, ClusterOutput> lOut = null) : base(output) { this.lOut = lOut; this.Name = name; randomV = r.Next(); }
        public ClosingForm Closing=null;

        public static List<string> GetVisOptions(ClusterOutput output)
        {
            if (output.hNode != null)
                return hNodeOptions;
            if (output.clusters != null || output.juryLike!=null)
                return clusterOptions;
               
            return null;
        }
        public override string ToString()
        {
            if (active != null)
                return Name + " " +randomV;       

            return "";
        }
        public void ActivateWindow()
        {
            if(active!=null)
                active.ToFront();           

        }
        public void CloseWindow()
        {
            if(active!=null)
                active.Close();          

        }
        public void SClusters(string item,string measureName,string option,bool graphics=false)
        {
            Dictionary<string, string> dic = ClusterOutput.ReadLabelsFile(output.GetLabelFile());
                if (output.hNode != null)
                {
                    // win = new visHierar(output.hNode,item,measureName);
                    if (option == null)
                        return;
                    if (option == "Dendrogram")
                        if (active == null || !(active is visHierar))
                        {
                            visHierar winH;
                            winH = new visHierar(output, item, measureName, dic);
                            winH.closeForm = Closing;
                            active = winH;
                            if (graphics)
                                winH.SaveToFile("result.png", true, 2, Color.Black, 1200, 800);
                            else
                                winH.Show();
                        }
                    if (option == "Sunburst chart")
                        if (active == null || !(active is VisHierarCircle))
                        {

                            VisHierarCircle winC;
                            winC = new VisHierarCircle(output, item, measureName);
                            winC.closeForm = Closing;
                            active = winC;
                            winC.Show();
                        }




                }
            if (output.clusters != null)
            {

                if (option == "Order Visual")
                {
                    if (active == null || !(active is VisOrder))
                    {
                        VisOrder visOrder;
                        visOrder = new VisOrder(output.clusters.list, item, null);
                        visOrder.closeForm = Closing;
                        active = visOrder;
                        visOrder.Show();
                    }
                }
                else
                {
                    if (output.nodes!=null && output.nodes.Count!=2 &&(active == null || !(active is ListVisual)))
                    {
                        ListVisual visBaker;
                        visBaker = new ListVisual(output, item, dic);
                        visBaker.closeForm = Closing;
                        active = visBaker;
                        visBaker.Show();
                    }

                }
            }

            if (output.nodes != null)
            {
                HeatMap heatRes = new HeatMap(output.nodes[0], output.nodes[1], null, output);
                heatRes.PrepareDataForHeatMap();
                active = heatRes;
                heatRes.Show();

            }
            else
            if (output.juryLike != null || output.hNNRes != null)
                {
                    if (active == null || !(active is FormText))
                    {

                        FormText showRes;
                        if (output.juryLike != null)
                            showRes = new FormText(output.juryLike, item);
                        else
                        {
                            showRes = new FormText(output.hNNRes, item);
                            showRes.Text = "Prediction";
                        }
                        showRes.closeForm = Closing;
                        active = showRes;
                        showRes.Show();
                    }
                }
              
            //}

        }

    }
}
