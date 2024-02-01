using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using phiClustCore;
using WorkFlows;

namespace Graph
{
    public partial class AddOmicsFile : Form
    {
        protected Dictionary<string, List<int>> superGenes=null;
        protected BindingList<FilterOmics> jFilters = new BindingList<FilterOmics>();
        OmicsDataSet joined = null;
        public AddOmicsFile()
        {
            InitializeComponent();
            filtersForJoined.DataSource = jFilters;
            filtersForJoined.DisplayMember = "Name";
        }

        private void FileButton_Click(object sender, EventArgs e)
        {
            DialogResult res=openFileDialog1.ShowDialog();
            if (res == DialogResult.OK)
            {
                if (checkBox1.Checked)
                {
                    groupBox1.Enabled = true;
                    groupBox2.Enabled = true;

                    omicsFiles.SelectedIndex = omicsFiles.Items.Count - 1;
                    button9.Enabled = true;


                    foreach (var item in openFileDialog1.FileNames)
                        omicsFiles.Items.Add(ReadOmicsFile.GetOmicsSpareFile(item));
                }
                else
                {
                    using (DataForm file = new DataForm(openFileDialog1.FileName))
                    {
                        res = file.ShowDialog();
                        if (res == DialogResult.OK)
                        {

                            groupBox1.Enabled = true;
                            groupBox2.Enabled = true;

                            omicsFiles.SelectedIndex = omicsFiles.Items.Count - 1;
                            button9.Enabled = true;

                            ReadOmicsFile omicsFile = new ReadOmicsFile(file.setup);
                            foreach (var item in openFileDialog1.FileNames)
                            {
                                //omicsFiles.Items.Add(file.GetOmicsFile(item));
                                omicsFiles.Items.Add(omicsFile.GetOmicsFile(item));
                            }

                        }
                    }
                }

            }
        }

        private void listBox1_SelectedIndexChanged(object sender, EventArgs e)
        {
            if(omicsFiles.Items.Count==0)
            {
                button1.Visible = false;
                button2.Visible = false;
            }
            else
            {
                button1.Visible = true;
                button2.Visible = true;
                if (omicsFiles.SelectedItem != null)
                {
                    filtersList.DataSource = ((OmicsDataSet)omicsFiles.SelectedItem).filters;
                    filtersList.DisplayMember = "Name";
                }
            }
        }

        private void filtersList_SelectedIndexChanged(object sender, EventArgs e)
        {
            if (filtersList.Items.Count == 0)
            {
                buttonUp.Enabled = false;
                buttonDown.Enabled = false;
                buttonRemoveFilter.Enabled = false;
            }
            else
            {
                buttonUp.Enabled = true;
                buttonDown.Enabled = true;
                buttonRemoveFilter.Enabled = true;
            }
        }

        private void buttonUp_Click(object sender, EventArgs e)
        {            
            int s = filtersList.SelectedIndex;
            if(s>0)
            {
                foreach (var item in omicsFiles.SelectedItems)
                {
                    OmicsDataSet aux = (OmicsDataSet)omicsFiles.SelectedItem;
                    FilterOmics o = aux.filters[s - 1];
                    aux.filters[s - 1] = aux.filters[s];
                    aux.filters[s] = o;
                    filtersList.SelectedIndex = s - 1;
                }
            }
        }

        private void buttonDown_Click(object sender, EventArgs e)
        {            
            int s = filtersList.SelectedIndex;
            if (s<filtersList.Items.Count-1)
            {
                foreach (var item in omicsFiles.SelectedItems)
                {
                    OmicsDataSet aux = (OmicsDataSet)omicsFiles.SelectedItem;
                    FilterOmics o = aux.filters[s + 1];
                    aux.filters[s + 1] = aux.filters[s];
                    aux.filters[s] = o;
                    filtersList.SelectedIndex = s + 1;
                }
            }
        }

        private void buttonRemoveFilter_Click(object sender, EventArgs e)
        {
            //            OmicsDataSet aux = (OmicsDataSet)omicsFiles.SelectedItem;
            if (filtersList.SelectedIndex > -1)
            {
                FilterOmics om = (FilterOmics)filtersList.SelectedItem;
                foreach (var item in omicsFiles.SelectedItems)
                {
                    OmicsDataSet aux = (OmicsDataSet)item;
                    aux.filters.Remove(om);
                }
            }
            
            
        }

        private void button2_Click(object sender, EventArgs e)
        {
            omicsFiles.Items.Remove(omicsFiles.SelectedItem);
            if(omicsFiles.Items.Count==0)
            {                                
                groupBox1.Enabled = false;
                groupBox2.Enabled = false;                
                button9.Enabled = false;
            }
        }

        private void button3_Click(object sender, EventArgs e)
        {
            
        }

        private void button1_Click(object sender, EventArgs e)
        {
            if (omicsFiles.SelectedItem == null)
            {
                MessageBox.Show("You must select file to which filter will be set");
                return;
            }
                FilterForm filters = new FilterForm((OmicsDataSet)omicsFiles.SelectedItem);
            DialogResult res = filters.ShowDialog();
            if (res == DialogResult.OK)
            {
                foreach (var item in omicsFiles.SelectedItems)
                {
                    OmicsDataSet aux = (OmicsDataSet)item;
                    aux.filters.Add(filters.selectedItem);
                }
                filtersList.Update();
                filtersList.Refresh();
            }
        }


        private void AddOmicsFile_Load(object sender, EventArgs e)
        {

        }

        private void listBox1_SelectedIndexChanged_1(object sender, EventArgs e)
        {
            if (filtersForJoined.SelectedIndex > -1)
                button8.Enabled = true;
        }

        private void joinedFilters_Click(object sender, EventArgs e)
        {
            FilterForm filters = new FilterForm((OmicsDataSet)omicsFiles.Items[0]);
            DialogResult res = filters.ShowDialog();
            if (res == DialogResult.OK)
            {             
                jFilters.Add(filters.selectedItem);
                filtersForJoined.Update();
                filtersForJoined.Refresh();
            }

        }

        private void button5_Click(object sender, EventArgs e)
        {            
            int s = filtersForJoined.SelectedIndex;
            if (s > 0)
            {
                FilterOmics o = jFilters[s - 1];
                jFilters[s - 1] = jFilters[s];
                jFilters[s] = o;
                filtersForJoined.SelectedIndex = s - 1;
            }

        }

        private void button6_Click(object sender, EventArgs e)
        {
            int s =filtersForJoined.SelectedIndex;
            if (s < filtersForJoined.Items.Count - 1)
            {
                FilterOmics o = jFilters[s + 1];
                jFilters[s + 1] = jFilters[s];
                jFilters[s] = o;
                filtersForJoined.SelectedIndex = s + 1;
            }

        }

        public virtual void button9_Click(object sender, EventArgs e)
        {
            List<OmicsDataSet> omicsFiltered = new List<OmicsDataSet>();
            FilterOmics.memoryFilteredData.Clear();
            if (omicsFiles.Items.Count > 0)
            {
                foreach (var item in omicsFiles.Items)
                {
                    OmicsDataSet aux = (OmicsDataSet)item;
                    aux.ClearFilters();
                    if (aux.filters.Count > 0)
                    {
                        var res = aux.ApplyFilters(aux.filters);
                    omicsFiltered.Add(res);
                }
                    else
                        omicsFiltered.Add(aux);
                }
                joined = OmicsDataSet.JoinOmicsData(omicsFiltered);
                //DataForm.CheckRanking(joined, 100);
                int s = joined.geneLabels.Count;
                if (s > 100)
                    s = 100;
                OmicsDataSet.CheckRanking(joined,s, "order");
                if (FilterOmics.remData)
                FilterOmics.memoryFilteredData.Add(joined);
                joined.filters = jFilters;
                Settings set;
                set = new Settings();
                set.Load();
                Form c = null;
                
                joined = joined.ApplyFilters(joined.filters);
                
                c = new ClusteringChoose(joined, set, this);
                c.Show();
                this.Hide();
            }
        }

        private void button10_Click(object sender, EventArgs e)
        {
            Close();
        }

        private void button8_Click(object sender, EventArgs e)
        {
            if (filtersForJoined.SelectedIndex > -1)
                jFilters.Remove((FilterOmics)filtersForJoined.SelectedItem);                
        }

        private void filtersList_DoubleClick(object sender, EventArgs e)
        {
            if (filtersList.SelectedIndex > -1)
            {
                FilterOmics selected = (FilterOmics)filtersList.SelectedItem;
                if (selected.Parameters)
                {
                    if(selected is LoadSuperGenes)
                    {
                        OmicsDataSet data = ((OmicsDataSet)omicsFiles.SelectedItem);
                        LoadSuperGenes cc = (LoadSuperGenes)selected;
                        SuperGenesForm super = new SuperGenesForm(data.geneLabels,cc.superGenes);
                        DialogResult res=super.ShowDialog();
                        if(res==DialogResult.OK)
                        {
                            FileN aux = new FileN();
                            aux.fileName = super.fileName;
                            cc.file.Add(aux);
                            cc.superGenes = super.superGenes;
                        }
                        
                    }
                    else
                        FilterForm.EditFilter(selected);
                    ((OmicsDataSet)omicsFiles.SelectedItem).filters[filtersList.SelectedIndex] = selected;
                    filtersList.Update();
                    filtersList.Refresh();
                }
                
            }
        }

        private void filtersForJoined_MouseDoubleClick(object sender, MouseEventArgs e)
        {
            if(filtersForJoined.SelectedIndex>-1)
            {
                FilterOmics selected = (FilterOmics)filtersForJoined.SelectedItem;
                if (selected.Parameters)
                {
                    FilterForm.EditFilter(selected);
                    jFilters[filtersForJoined.SelectedIndex] = selected;
                    filtersForJoined.Update();
                    filtersForJoined.Refresh();
                }

            }
        }
    }
}
    