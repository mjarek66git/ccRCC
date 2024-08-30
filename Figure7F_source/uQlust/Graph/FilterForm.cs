using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using phiClustCore;
using phiClustCore.Distance;

namespace Graph
{
    public partial class FilterForm : Form
    {
        public FilterOmics selectedItem;
        OmicsDataSet data;

        public FilterForm(OmicsDataSet data)
        {
            this.data = data;
            InitializeComponent();
            Type[] xx = Assembly.GetAssembly(typeof(FilterOmics)).GetTypes();
            IEnumerable<Type> subclasses = xx.Where(t => t.IsSubclassOf(typeof(FilterOmics)));
            foreach (var item in subclasses)
                listBox1.Items.Add((FilterOmics)Activator.CreateInstance(item));

        }

        public static void EditFilter(FilterOmics selectedItem)
        {
            if (selectedItem != null && selectedItem.Parameters)
            {
                Dictionary<string, Control> controls = new Dictionary<string, Control>();
                Dictionary<string, string> paramsValue = new Dictionary<string, string>();
                var types = selectedItem.GetParameters();
                Form form = new Form();
                int x = 10, y = 10;
                foreach (var item in types)
                {
                    var label = new Label();
                    label.Location = new Point(x, y);
                    label.Text = item.Key;
                    form.Controls.Add(label);
                    if (item.Value == typeof(int) || item.Value == typeof(double) || item.Value == typeof(string))
                    {
                        TextBox t = new TextBox();
                        controls.Add(item.Key, t);
                        t.Location = new Point(x + 130, y);
                        form.Controls.Add(t);
                    }
                    if (item.Value == typeof(CodingAlg))
                    {

                        ComboBox c = new ComboBox();
                        c.DataSource = Enum.GetValues(item.Value);
                        controls.Add(item.Key, c);
                        c.Location = new Point(x + 130, y);
                        form.Controls.Add(c);
                    }
                    if (item.Value == typeof(bool))
                    {
                        CheckBox b = new CheckBox();
                        controls.Add(item.Key, b);
                        b.Location = new Point(x + 130, y);
                        form.Controls.Add(b);
                    }

                    y += 30;
                }
                Button bOK = new Button();
                bOK.Text = "OK";
                bOK.Location = new Point(10, form.Height - 70);
                bOK.Click += (send, ex) =>
                {
                    foreach (var item in controls)
                    {
                        string val = "";
                        if (item.Value is TextBox)
                        {
                            if (item.Value.Text.Length == 0)
                                return;
                            val = item.Value.Text;
                        }
                        if (item.Value is ComboBox)
                        {
                            val = ((ComboBox)item.Value).SelectedItem.ToString();
                        }
                        if (item.Value is CheckBox)
                        {
                            val = ((CheckBox)item.Value).Checked.ToString();
                        }
                        paramsValue.Add(item.Key, val);
                    }
                    selectedItem.SetParameters(paramsValue);
                    form.DialogResult = DialogResult.OK;
                    form.Close();
                };
                Button bCancel = new Button();
                bCancel.Text = "Cancel";
                bCancel.Location = new Point(form.Width - 100, form.Height - 70);
                bCancel.Click += (send, ex) => { form.DialogResult = DialogResult.Cancel; form.Close(); };

                form.Controls.Add(bOK);
                form.Controls.Add(bCancel);
                form.ShowDialog();
            }
        }



        private void listBox1_MouseDoubleClick(object sender, MouseEventArgs e)
        {
            selectedItem = (FilterOmics)listBox1.SelectedItem;
            if (data == null)
                return;
            if (!(selectedItem is LoadSuperGenes) && !(selectedItem is SelectStage))
                EditFilter(selectedItem);
            else
            {
                if (selectedItem is LoadSuperGenes)
                {
                    SuperGenesForm super = new SuperGenesForm(data.geneLabels, null);
                    DialogResult res = super.ShowDialog();
                    if (res == DialogResult.OK)
                    {
                        LoadSuperGenes superF = (LoadSuperGenes)selectedItem;
                        FileN aux = new FileN();
                        aux.fileName = super.fileName;
                        superF.file.Add(aux);
                        superF.superGenes = super.superGenes;
                    }
                }
                if (selectedItem is SelectStage)
                {
                    OpenFileDialog fileStage = new OpenFileDialog();
                    DialogResult res = fileStage.ShowDialog();
                    if (res == DialogResult.OK)
                    {
                        StaticDic.Id.Clear();
                        StaticDic.LoadId(fileStage.FileName);
                        SelectStage stageFilter = (SelectStage)selectedItem;
                        HashSet<string> stages = new HashSet<string>();
                        foreach (var item in StaticDic.Id)
                            stages.Add(item.Value);

                        EditFilter(selectedItem);
                    }
                }
            }
                DialogResult = DialogResult.OK;
                this.Close();
            
        }
    }
}
