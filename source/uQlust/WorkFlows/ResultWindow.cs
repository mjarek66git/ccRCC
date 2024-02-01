using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.IO;
using System.Reflection;
using System.Threading.Tasks;
using System.Windows.Forms;
using Graph;
using phiClustCore;

namespace WorkFlows
{
    public partial class ResultWindow : Form
    {
        JobManager manager = new JobManager();
        public ResultWindow()
        {
            InitializeComponent();
            //System.Windows.Forms.DataGridViewCellStyle dataGridViewCellStyle1 = new System.Windows.Forms.DataGridViewCellStyle();
            manager.beginJob = AddItemToList;
            manager.updateJob = UpdateListBox;
            manager.message = ErrorMessage;
            TimeInterval.InitTimer(UpdateProgress);
        }

        void ErrorMessage(string message)
        {
            MessageBox.Show(message);
        }

        public void Run(string processName,OmicsDataSet data, Options opt)
        {
            manager.opt = opt;

            try
            {
                TimeInterval.Start();
                manager.RunJob(processName,data);
            }
            catch (Exception ex)
            {
                TimeInterval.Stop();
                MessageBox.Show("Exception: " + ex.Message);
            }
        }
        private void UpdateProgress(object sender, EventArgs e)
        {

            //    int xx=(int)(manager.ProgressUpdate()*100);
            //dataGridView1.Rows[0].Cells[6].Value = xx;


            Dictionary<string, double> res = manager.ProgressUpdate();

            if (res == null)
            {
                // TimeInterval.Stop();
                return;
            }

            for (int i = 0; i < dataGridView1.Rows.Count; i++)
                if (res.ContainsKey((string)dataGridView1.Rows[i].Cells[0].Value))
                    dataGridView1.Rows[i].Cells[6].Value = (int)(res[(string)dataGridView1.Rows[i].Cells[0].Value] * 100);

        }
        public delegate void AddListBoxItem(string item, string clType, string dirName, string measure);        
        public void AddItemToList(string item, string clType, string dirName, string measure)
        {
            int num;
            if (dataGridView1.InvokeRequired)
                dataGridView1.Invoke(new AddListBoxItem(AddItemToList), new object[] { item, clType, dirName, measure });
            else
            {

                dataGridView1.Rows.Add(1);
                num = dataGridView1.Rows.Count - 1;
                dataGridView1.Rows[num].DefaultCellStyle.ForeColor = Color.Red;
                dataGridView1.Rows[num].Cells[0].Value = item;
                dataGridView1.Rows[num].Cells[1].Value = clType;
                dataGridView1.Rows[num].Cells[2].Value = measure;
                dataGridView1.Rows[num].Cells[3].Value = dirName;
            }

        }
        private void UpdateGridView()
        {
            dataGridView1.Rows.Clear();
            foreach (var item in manager.clOutput.Keys)
            {
                int num;
                dataGridView1.Rows.Add(1);
                num = dataGridView1.Rows.Count - 1;
                dataGridView1.Rows[num].DefaultCellStyle.ForeColor = Color.Green;
                dataGridView1.Rows[num].Cells[0].Value = manager.clOutput[item].name;
                dataGridView1.Rows[num].Cells[1].Value = manager.clOutput[item].clusterType;
                dataGridView1.Rows[num].Cells[2].Value = manager.clOutput[item].measure;
                dataGridView1.Rows[num].Cells[4].Value = manager.clOutput[item].time;
                dataGridView1.Rows[num].Cells[3].Value = manager.clOutput[item].dirName;
            }
        }
        private void gotoProcessToolStripMenuItem_Click(object sender, EventArgs e)
        {
            menuClick();
        }
        private void closeAllToolStripMenuItem_Click(object sender, EventArgs e)
        {

        }
        private void endProcessToolStripMenuItem_Click(object sender, EventArgs e)
        {
            manager.RemoveJob((string)dataGridView1.Rows[dataGridView1.SelectedRows[0].Index].Cells[0].Value);

            if (manager.clOutput.ContainsKey((string)dataGridView1.Rows[dataGridView1.SelectedRows[0].Index].Cells[0].Value))
                manager.clOutput.Remove((string)dataGridView1.Rows[dataGridView1.SelectedRows[0].Index].Cells[0].Value);
            dataGridView1.Rows.RemoveAt(dataGridView1.SelectedRows[0].Index);
        }

        private void saveCurrentinTextModeToolStripMenuItem_Click(object sender, EventArgs e)
        {
            saveToolStripMenuItem1_Click(sender, e);
        }
        private void saveToolStripMenuItem1_Click(object sender, EventArgs e)
        {
            DialogResult res;
            //  saveFileDialog1 = new OpenFileDialog();
            saveFileDialog1.InitialDirectory =Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location)+Path.DirectorySeparatorChar+ "results";
            saveFileDialog1.FileName = "result_" + (string)dataGridView1.Rows[dataGridView1.SelectedRows[0].Index].Cells[0].Value+".out";

            res = saveFileDialog1.ShowDialog();
            if (res == DialogResult.OK)
            {
                ClusterOutput sel = GetSelectedOutput();
                if (sel != null)
                {
                    if (sel.hNode != null)
                    {
                       // MessageBox.Show("For hierarchical clustering only 10 top clusters are saved");
                    }
                    sel.SaveTxt(saveFileDialog1.FileName);
                }
                else
                    MessageBox.Show("Something went wrong. Cannot save results!");


            }
        }
        private void contextMenuStrip1_Opening(object sender, CancelEventArgs e)
        {
            if (dataGridView1.Rows.Count == 0)
            {

                for (int i = 0; i < contextMenuStrip1.Items.Count; i++)
                    if (contextMenuStrip1.Items[i].Text.Contains("Save") || contextMenuStrip1.Items[i].Text.Contains("Show") ||
                        contextMenuStrip1.Items[i].Text.Contains("Close"))
                        contextMenuStrip1.Items[i].Enabled = false;
                return;
            }
            else
                for (int i = 0; i < contextMenuStrip1.Items.Count; i++)
                    contextMenuStrip1.Items[i].Enabled = true;

            if (dataGridView1.Rows[dataGridView1.SelectedRows[0].Index].Cells[4].Value == null)
            {
                for (int i = 0; i < contextMenuStrip1.Items.Count; i++)
                    if (contextMenuStrip1.Items[i].Text.Contains("Save current") || contextMenuStrip1.Items[i].Text.Contains("Show"))
                        contextMenuStrip1.Items[i].Enabled = false;
                gotoProcessToolStripMenuItem.Enabled = false;
                finishProcessToolStripMenuItem.Text = "End Process";
            }
            else
            {
                finishProcessToolStripMenuItem.Text = "Close result";
                gotoProcessToolStripMenuItem.Enabled = true;
            }
            if (dataGridView1.Rows[dataGridView1.SelectedRows[0].Index].DefaultCellStyle.ForeColor == Color.Green)
            {
                if (manager.clOutput.ContainsKey((string)dataGridView1.Rows[dataGridView1.SelectedRows[0].Index].Cells[0].Value))
                {
                    List<string> options = ClusterGraphVis.GetVisOptions(manager.clOutput[(string)dataGridView1.Rows[dataGridView1.SelectedRows[0].Index].Cells[0].Value]);
                    gotoProcessToolStripMenuItem.DropDownItems.Clear();
                    if (options != null)
                    {
                        gotoProcessToolStripMenuItem.Click -= gotoProcessToolStripMenuItem_Click;
                        foreach (var item in options)
                        {
                            ToolStripItem it = gotoProcessToolStripMenuItem.DropDownItems.Add(item);
                            it.Click += new EventHandler(SubMenuClick);

                        }
                    }
                    else
                    {
                        gotoProcessToolStripMenuItem.Click -= gotoProcessToolStripMenuItem_Click;
                        gotoProcessToolStripMenuItem.Click += gotoProcessToolStripMenuItem_Click;
                    }


                }
            }
        }

        private ClusterOutput GetSelectedOutput()
        {

            DataGridViewCellCollection cell = dataGridView1.Rows[dataGridView1.SelectedRows[0].Index].Cells;
            string name = (string)cell[0].Value;
            if (cell[cell.Count - 1].Value != null)
                name = (string)cell[cell.Count - 1].Value;

            if (manager.clOutput.ContainsKey(name))
                return manager.clOutput[name];

            return null;
        }
        private void SubMenuClick(object sender, EventArgs e)
        {
            ToolStripItem item = sender as ToolStripItem;
            try
            {
                DataGridViewCellCollection cell = dataGridView1.Rows[dataGridView1.SelectedRows[0].Index].Cells;
                string name = (string)cell[0].Value;
                ClusterGraphVis vis = new ClusterGraphVis(manager.clOutput[name], name, manager.clOutput);
                vis.SClusters(name, cell[2].Value.ToString(),item.Text);
                // cell[cell.Count - 1].Value = vis.ToString();

            }
            catch (Exception ex)
            {
                MessageBox.Show("Error: " + ex.Message);
            }

        }
        private void finishProcessToolStripMenuItem_Click(object sender, EventArgs e)
        {

        }

        private void menuClick()
        {
            if (dataGridView1.Rows[dataGridView1.SelectedRows[0].Index].DefaultCellStyle.ForeColor == Color.Green)
            {
                DataGridViewCellCollection cell = dataGridView1.Rows[dataGridView1.SelectedRows[0].Index].Cells;
                string name = (string)cell[0].Value;
                //if(cell[cell.Count-1].Value!=null)
                // name = (string)cell[cell.Count-1].Value;

                if (manager.clOutput.ContainsKey(name))
                {
                    ClusterGraphVis vis = new ClusterGraphVis(manager.clOutput[name], name);
                    vis.SClusters(name, dataGridView1.Rows[dataGridView1.SelectedRows[0].Index].Cells[2].Value.ToString(), null);

                }
            }
        }
        public delegate void RemoveListBoxItem(string item, bool errorFlag, bool finishAll);

        private void UpdateListBox(string item, bool errorFlag, bool finishAll)
        {

            if (dataGridView1.InvokeRequired)
                dataGridView1.Invoke(new RemoveListBoxItem(UpdateListBox), new object[] { item, errorFlag, finishAll });
            else
            {
                if (finishAll)
                    TimeInterval.Stop();
                if (manager.clOutput.Count != 0)
                {
                    UpdateProgress(null, null);
                    TimeSpan span = DateTime.Now.Subtract(DateTime.Now);
                    if (manager.clOutput.ContainsKey(item))
                    {
                        int remIndex = 0;
                        //                        vis = new ClusterGraphVis(manager.clOutput[item]);
                        for (int i = 0; i < dataGridView1.Rows.Count; i++)
                            if ((string)dataGridView1.Rows[i].Cells[0].Value == item)
                            {
                                dataGridView1.Rows[i].DefaultCellStyle.ForeColor = Color.Green;
                                dataGridView1.Rows[i].Cells[4].Value = manager.clOutput[item].time;
                                dataGridView1.Rows[i].Cells[1].Value = manager.clOutput[item].clusterType;
                                dataGridView1.Rows[i].Cells[2].Value = manager.clOutput[item].measure;
                                dataGridView1.Rows[i].Cells[5].Value = manager.clOutput[item].peekMemory;
                                dataGridView1.Rows[i].Cells[6].Value = 100;
                                remIndex = i;
                                button1.Enabled = true;
                                button2.Enabled = true;
                            }
                        dataGridView1.Rows[remIndex].Selected = true;

                    }
                }
                else
                {
                    for (int i = 0; i < dataGridView1.Rows.Count; i++)
                    {
                        dataGridView1.Rows[i].DefaultCellStyle.ForeColor = Color.Red;
                        dataGridView1.Rows[i].Cells[4].Value = "Error";
                        dataGridView1.Rows[i].Cells[5].Value = "Error";
                    }

                }
            }
        }

        private void button1_Click(object sender, EventArgs e)
        {
            DataGridViewCellCollection cell = dataGridView1.Rows[dataGridView1.SelectedRows[0].Index].Cells;
            string name = (string)cell[0].Value;
            ClusterGraphVis vis = new ClusterGraphVis(manager.clOutput[name], name, manager.clOutput);

            vis.SClusters(name, cell[2].Value.ToString(), ""); ;

        }

        private void button2_Click(object sender, EventArgs e)
        {
            saveCurrentinTextModeToolStripMenuItem_Click(sender, e);
        }

        private void saveAllbinaryToolStripMenuItem_Click(object sender, EventArgs e)
        {
            saveFileDialog1.Filter = "Results file (*.cres)|*.cres";
            DialogResult res = saveFileDialog1.ShowDialog();
            if (res == DialogResult.OK)
                if (saveFileDialog1.FileName.Length > 0)
                    manager.SaveOutput(saveFileDialog1.FileName + ".cres");            
        }

        private void loadbinaryToolStripMenuItem_Click(object sender, EventArgs e)
        {
            openFileDialog1.Filter = "Results (.cres)|*.cres|All files (*.*)|*.*";
            openFileDialog1.FilterIndex = 1;
            DialogResult res = openFileDialog1.ShowDialog();
            if (res == DialogResult.OK)
                if (openFileDialog1.FileName.Length > 0)
                {
                    manager.LoadOutput(openFileDialog1.FileName);
                    UpdateGridView();
                }
        }



    }
}
