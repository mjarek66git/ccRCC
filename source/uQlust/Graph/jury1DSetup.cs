using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Data;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using phiClustCore;
//using ;
using Graph;

namespace Graph
{
    public delegate void ProfileNameChanged(string newName);
    public partial class jury1DSetup : UserControl
    {
        private string fileName;
        
        public ProfileNameChanged profileNameChanged=null;

        private string ProfileName=null;
        public List<INPUTMODE> inputMode = null;

        public jury1DSetup()
        {
            InitializeComponent();
        }
        
        public string profileName
        {
            set
            {
                ProfileName = value;
                label2.Text = ProfileName;
                UpdateLabels();
            }
            get
            {
                return ProfileName;
            }
        }
        public bool alignGenerate;
      
       
        private void button6_Click(object sender, EventArgs e)
        {
          
           ProfileForm profDef;

            profDef = new ProfileForm(profileName, "",filterOPT.SIMILARITY);

            DialogResult res = profDef.ShowDialog();
            if (res == DialogResult.OK)
            {
                try
                {

                    profileName = profDef.fileName;
                    if (profileNameChanged != null)
                        profileNameChanged(profileName);

                    profDef.SaveProfiles(profileName);
                   
                    fileName = profDef.alignFileName;
                    if (fileName != null && fileName.Length > 0)
                        alignGenerate = true;
                }
                catch
                {
                    MessageBox.Show("Profiles could not be saved!");
                }

                UpdateLabels();
            }

        }
        public void UpdateLabels()
        {
            if(ProfileName!=null && ProfileName.Length>0)
                label18.Text = LabelActiveProfiles(ProfileName);

        }
        private string LabelActiveProfiles(string fileName)
        {
            string outStr = "";
            ProfileTree tree = new ProfileTree();
            try
            {
                tree.LoadProfiles(fileName);
            }
            catch(Exception)
            {
                label2.Text = "";
                return null;
            }
            inputMode = tree.GetModes();
            List<profileNode> active = tree.GetActiveProfiles();
            outStr = "Active profiles: ";
            if(active!=null)
                for (int i = 0; i < active.Count; i++)
                {
                    outStr += active[i].profName;
                    if (i < active.Count - 1)
                        outStr += ", ";
                }

            return outStr;
        }

    }
}
