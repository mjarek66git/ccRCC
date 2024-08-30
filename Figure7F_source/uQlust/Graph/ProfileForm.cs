using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Text.RegularExpressions;
using phiClustCore;
using phiClustCore.Profiles;

namespace Graph
{
    public enum filterOPT
    {
        SIMILARITY,
        DISTANCE
    };
    public partial class ProfileForm : Form
    {        
        ProfileTree treeProfiles;
        public string fileName="defualt.profile";
        public string alignFileName;
        private filterOPT filter;
        public ProfileForm(string fileName,string alignFile,filterOPT filter)
        {
            InitializeComponent();
            if (fileName == null || fileName.Length == 0)
                this.Text = "Profile not defined";
            else
            {
                this.Text = fileName;
                this.fileName = fileName;
            }

            this.filter = filter;
            //profilesView.DrawMode = TreeViewDrawMode.OwnerDrawText;
            
            treeProfiles = new ProfileTree();
            if(fileName!=null && fileName.Length>0)
                if(File.Exists(fileName))
                    treeProfiles.LoadProfiles(fileName);
           
            ShowTree();
            //profilesView.
        }
        public void SaveProfiles(string fileName)
        {
            treeProfiles.SaveProfiles(fileName);
        }
        private void ShowTree()
        {
            profilesView.Nodes.Clear();
            profilesView.Nodes.Add("Profiles");
            if (treeProfiles.masterNode.Count > 0)
            {
                TreeNode tNode;
                TreeNodeCollection nodes = null;
                Dictionary<string, profileNode> aux = new Dictionary<string, profileNode>();
                treeProfiles.GetProfiles("/", treeProfiles.masterNode, aux, false);
                profilesView.PathSeparator = "/";
                profilesView.BeginUpdate();
                foreach (var item in aux)
                {
                    string w = "Profiles" + item.Key;
                    string[] tmp = w.Split('/');
                    string currentPath = "";
                    nodes = profilesView.Nodes;
                    foreach (var s in tmp)
                    {
                        if (s.Length == 0)
                            continue;
                        if (s != "Profiles")
                            currentPath += "/" + s;
                        //              else
                        //                   continue;
                        if (nodes != null)
                        {
                            bool flag = false;
                            foreach (TreeNode nItem in nodes)
                            {
                                if (nItem.Text == s)
                                {
                                    flag = true;
                                    nodes = nItem.Nodes;
                                    break;
                                }
                            }
                            if (!flag)
                            {
                                tNode = nodes.Add(s);
                                tNode.Tag = treeProfiles.FindNode(currentPath);
                                //tNode.Text.
                            }
                        }
                        else
                        {
                            tNode = nodes.Add(s);
                            tNode.Tag = treeProfiles.FindNode(currentPath);
                        }
                    }

                }
                profilesView.EndUpdate();
                profilesView.Nodes[0].Expand();
            }

        }
        private void GetNodes(TreeNode node, ProfileTree tree)
        {
            profileNode nodeP;
            string parentPath = node.FullPath;
            Regex reg = new Regex("^Profiles");
            parentPath = reg.Replace(parentPath, "");
            foreach (TreeNode item in node.Nodes)
            {                
                string Path=item.FullPath;
                Path=reg.Replace(Path,"");
                nodeP=treeProfiles.FindNode(Path);
                tree.AdddNode(parentPath, new profileNode(nodeP));
                if (item.Nodes != null)
                    GetNodes(item, tree);

            }
        }
        private void EditProfile(string Path, bool edit,profileType type)
        {
            ProfileDefinitionForm profDef = null ;
            InternalProfileForm internalDef = null;
            profileNode node = null;
            profileNode newNode = null;
            DialogResult res=DialogResult.Cancel;

            if (edit)
            {
                    node = treeProfiles.FindNode(Path);
                    if (node.internalName!=null)
                    {
                        internalDef = new InternalProfileForm(node,filter);
                        res = internalDef.ShowDialog();
                        if (res == DialogResult.OK)
                        {
                            newNode = internalDef.localNode;
                            toolSave.Enabled = true;
                            toolStripButton6.Enabled = true;
                        }
                    }
                    else
                    {
                        profDef = new ProfileDefinitionForm(node);                        
                        res = profDef.ShowDialog();
                        if (res == DialogResult.OK)
                        {
                            newNode = profDef.profile;
                            toolSave.Enabled = true;
                            toolStripButton6.Enabled = true;
                        }
                    }
            }
            else
            {
                if (Path == null)
                    Path = "Profiles";
                if (type == profileType.EXTERNAL)
                {
                    profDef = new ProfileDefinitionForm();
                    res = profDef.ShowDialog();
                    if (res == DialogResult.OK)
                    {
                        node = profDef.profile;
                        node.internalName = null;
                        toolSave.Enabled = true;
                        toolStripButton6.Enabled = true;
                    }
                }
                else
                {
                    Settings set=new Settings();
                    set.Load();
                    InternalProfilesManager manager = new InternalProfilesManager();
                    List<profileNode> validProfiles = new List<profileNode>();
                    foreach (var item in InternalProfilesManager.internalList.Keys)
                        if (manager.CheckAccessibility(item, set.mode))
                            validProfiles.Add(item);

                    ListInternal intr = new ListInternal(validProfiles);

                    res = intr.ShowDialog();

                    if (res == DialogResult.OK)
                    {
                        InternalProfileForm intForm;
                        profileNode localNode = InternalProfilesManager.GetNode(intr.selectedProfile);
                        if(localNode.profName.Contains("Load"))
                            intForm = new InternalProfileForm(localNode,filter,false);
                        else
                            intForm = new InternalProfileForm(localNode,filter);

                        res = intForm.ShowDialog();
                        if (res == DialogResult.OK)
                        {
                            node = intForm.localNode;
                            node.profProgram = intr.selectedProfile;
                            newNode = node;
                            toolSave.Enabled = true;
                            toolStripButton6.Enabled = true;

                        }
                    }
                }
            }

            if (res == DialogResult.OK)
            {
                if (node != null && File.Exists(node.OutFileName))
                    File.Delete(node.OutFileName);
                //profChanged = true;

                if (!edit)
                {
                    //profChanged = true;
                    Regex exp = new Regex("^/");

                    TreeNode nodeT;//= new TreeNode(profDef.profile.profName);

                    string nPath = "";
                    if (Path != null)
                        nPath = exp.Replace(Path, "");
                    if (profilesView.SelectedNode == null)
                    {
                        nodeT = profilesView.Nodes[0].Nodes.Add(node.profName);
                        profilesView.Nodes[0].Expand();
                    }
                    else
                    {
                        nodeT = profilesView.SelectedNode.Nodes.Add(node.profName);
                        profilesView.SelectedNode.Expand();
                    }
                    nodeT.Tag = newNode;
                    Regex exp2 = new Regex("^Profiles");
                    Path = exp2.Replace(Path, "");
                    treeProfiles.AdddNode(Path, node);
                }
                else
                {
                    profilesView.SelectedNode.Text = node.profName;
                    treeProfiles.RemoveNode(Path);
                    Regex exp2 = new Regex("^Profiles");
                    Path = exp2.Replace(Path, "");
                    treeProfiles.AdddNode(ParentPath(Path),newNode);
                    node.CopyNode(newNode);
                    ShowTree();
                }

            }
        }
        private string PathRemoveRoot(string Path)
        {
            string []tmp=Path.Split('/');

            return "/"+String.Join("/",tmp,1,tmp.Length-1);

        }
        private string ParentPath(string Path)
        {
            string[] tmp = Path.Split('/');
            Path = String.Join("/", tmp, 0, tmp.Length - 1);
            if (Path.Length == 0)
                Path = "/";
            return Path;

        }
        private void MakeViewCopy(ProfileTree tr,TreeNode node)
        {
            foreach (TreeNode item in node.Nodes)
            {
                string Path = item.FullPath;
                Path = PathRemoveRoot(Path);
                Path = ParentPath(Path);
                tr.AdddNode(Path, new profileNode((profileNode)item.Tag));

                if (item.Nodes.Count > 0)
                    MakeViewCopy(tr,item);


            }
        }

        private void toolStripButton1_Click(object sender, EventArgs e)
        {
            if (profilesView.SelectedNode == null)
                EditProfile(null, false,profileType.EXTERNAL);
            else
            {
                string Path = PathRemoveRoot(profilesView.SelectedNode.FullPath);
                EditProfile(Path,false,profileType.EXTERNAL);
            }
        }

        private void toolStripButton3_Click(object sender, EventArgs e)
        {
            if (profilesView.SelectedNode!=null)
            {
                string Path = profilesView.SelectedNode.FullPath;
                Regex exp2 = new Regex("^Profiles");
                profilesView.Nodes.Remove(profilesView.SelectedNode);
                Path = exp2.Replace(Path, "");
                treeProfiles.RemoveNode(Path);
            }
        }

        private void toolStripButton2_Click(object sender, EventArgs e)
        {
            if (profilesView.SelectedNode!=null)
            {
                string path = profilesView.SelectedNode.FullPath;
                path=path.Replace("\\","/");
                EditProfile(path,true,profileType.EXTERNAL);
            }
        }
        private void saveBtn_Click(object sender, EventArgs e)
        {
            bool test = true;
            if (fileName.Length > 0)
                treeProfiles.SaveProfiles(fileName);
            else
                test=SaveAs();
            if (test)
            {
                this.DialogResult = DialogResult.OK;
                this.Close();
            }
        }


        private void cancelBtn_Click(object sender, EventArgs e)
        {
            //profChanged = false;
            DialogResult = DialogResult.Cancel;
            this.Close();
        }
        private void ActiveNode()
        {
            if (profilesView.SelectedNode != null)
            {
                if (profilesView.SelectedNode.Tag != null)
                {
                    if (((profileNode)profilesView.SelectedNode.Tag).childrens.Count != 0)
                        return;

                     ((profileNode)profilesView.SelectedNode.Tag).active = !((profileNode)profilesView.SelectedNode.Tag).active;
                     TreeNode node = profilesView.SelectedNode;
                     while (node.Text != "Profiles")
                     {
                        node = node.Parent;
                        if (node.Tag != null)
                        {
                            if (((profileNode)node.Tag).active)
                            {
                                bool flag = false;
                                foreach (var item in ((profileNode)node.Tag).childrens)
                                    if (item.Value.active)
                                    {
                                        flag = true;
                                        break;
                                    }
                                    if (!flag)
                                        ((profileNode)node.Tag).active = false;
                             }
                             else
                                ((profileNode)node.Tag).active = !((profileNode)node.Tag).active;
                         }
                      }
                      profilesView.Invalidate();
                }
            }
        
        }
        private void profilesView_DrawNode(object sender, DrawTreeNodeEventArgs e)
        {
            Color nodeColor=Color.Black;

            if (e.Node.Tag == null)
                nodeColor = Color.Black;
            else
                if(((profileNode) (e.Node.Tag)).active)
                    nodeColor = Color.Green;
                else
                    nodeColor = Color.Red;


            if ((e.State & TreeNodeStates.Selected) != 0)
                nodeColor = SystemColors.HighlightText;

            TextRenderer.DrawText(e.Graphics,
                                  e.Node.Text,
                                  e.Node.NodeFont,
                                  e.Bounds,
                                  nodeColor,
                                  Color.Empty,
                                  TextFormatFlags.VerticalCenter);
        }

        private void ActiveButton_Click(object sender, EventArgs e)
        {
            ActiveNode();
        }


        private void profilesView_AfterSelect(object sender, TreeViewEventArgs e)
        {
            if (profilesView.Nodes.Count > 0)
            {
                if (profilesView.SelectedNode != null)
                {
                    if (profilesView.SelectedNode.Text == "Profiles")
                    {
                        ActiveButton.Enabled = false;
                        toolStripButton3.Enabled = false;
                        toolStripButton2.Enabled = false;

                    }
                    else
                    {
                        toolStripButton2.Enabled = true;
                        toolStripButton3.Enabled = true;
                        ActiveButton.Enabled = true;

                    }
                }
            }

        }
        private bool SaveAs()
        {
            DialogResult res;
            string dir = Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location);

            if (filter == filterOPT.SIMILARITY)
                saveFileDialog1.Filter = "Profiles file (*.profiles)|*.profiles";
            else
                saveFileDialog1.Filter = "Profiles file (*distance.profile)|*distance.profile";

            if (Directory.Exists(dir + Path.DirectorySeparatorChar + "profiles"))
                dir = dir + Path.DirectorySeparatorChar + "profiles";

            saveFileDialog1.InitialDirectory = dir;
            res = saveFileDialog1.ShowDialog();
            if (res == DialogResult.OK)
            {

                string profileNameFile = saveFileDialog1.FileName;
                if (filter == filterOPT.SIMILARITY)
                {
                    //if (!profileNameFile.EndsWith(".profiles"))
                      //  profileNameFile = profileNameFile.Replace(".profiles", "_jury.profiles");
                }
                else
                    if (!profileNameFile.EndsWith("distance.profile"))
                        profileNameFile = profileNameFile.Replace(".profile", "_distance.profile");

                this.Text = profileNameFile;
                treeProfiles.SaveProfiles(profileNameFile);
                fileName = profileNameFile;
                return true;
            }
            return false;
        }
        private void toolSave_Click(object sender, EventArgs e)
        {
            SaveAs();
        }

        private void toolLoad_Click(object sender, EventArgs e)
        {
            DialogResult res;
            if(filter==filterOPT.SIMILARITY)
                openFileDialog1.Filter = "Profiles file (*.profiles)|*.profiles";
            else
                openFileDialog1.Filter = "Profiles file (*distance.profile)|*distance.profile";
            string dir = Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location);

            if (Directory.Exists(dir + Path.DirectorySeparatorChar + "profiles"))
                dir = dir + Path.DirectorySeparatorChar + "profiles";
            openFileDialog1.InitialDirectory = dir;
            res = openFileDialog1.ShowDialog();
            if (res == DialogResult.OK)
            {
                fileName = openFileDialog1.FileName;
                treeProfiles.LoadProfiles(fileName);
                this.Text = fileName;
                ShowTree();
            }

        }


        private void toolStripButton6_Click(object sender, EventArgs e)
        {
            if (this.Text.Contains("Not defined"))
                toolSave_Click(sender, e);
            else
                treeProfiles.SaveProfiles(this.Text);                
        }



        private void toolStripButton4_Click(object sender, EventArgs e)
        {
            if (profilesView.SelectedNode == null)
                EditProfile(null, false, profileType.INTERNAL);
            else
            {
                string Path = PathRemoveRoot(profilesView.SelectedNode.FullPath);
                EditProfile(Path, false, profileType.INTERNAL);
            }

        }

        private void toolStripButton7_Click(object sender, EventArgs e)
        {
            profilesView.Nodes.Clear();
            treeProfiles.ClearTree();
            profilesView.Nodes.Add("Profiles");
            toolStripButton2.Enabled = false;
            toolStripButton3.Enabled = false;
            ActiveButton.Enabled = false;
            toolSave.Enabled = false;
            toolStripButton6.Enabled = false;
            fileName = "";
            this.Text = "Profile not defined";
           
        }

        private void button1_Click(object sender, EventArgs e)
        {

        }

    }
    
}
