namespace WorkFlows
{
    partial class ResultWindow
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.components = new System.ComponentModel.Container();
            this.dataGridView1 = new System.Windows.Forms.DataGridView();
            this.process = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.Cluster = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.Measure = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.Column1 = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.time = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.Column2 = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.Progress = new Graph.ProgressColumn();
            this.Column3 = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.contextMenuStrip1 = new System.Windows.Forms.ContextMenuStrip(this.components);
            this.gotoProcessToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.saveCurrentinTextModeToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.finishProcessToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.closeAllToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.progressColumn1 = new Graph.ProgressColumn();
            this.saveFileDialog1 = new System.Windows.Forms.SaveFileDialog();
            this.button1 = new System.Windows.Forms.Button();
            this.button2 = new System.Windows.Forms.Button();
            this.saveAllbinaryToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.loadbinaryToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.openFileDialog1 = new System.Windows.Forms.OpenFileDialog();
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).BeginInit();
            this.contextMenuStrip1.SuspendLayout();
            this.SuspendLayout();
            // 
            // dataGridView1
            // 
            this.dataGridView1.AllowUserToAddRows = false;
            this.dataGridView1.AllowUserToDeleteRows = false;
            this.dataGridView1.AllowUserToResizeRows = false;
            this.dataGridView1.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.dataGridView1.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.dataGridView1.Columns.AddRange(new System.Windows.Forms.DataGridViewColumn[] {
            this.process,
            this.Cluster,
            this.Measure,
            this.Column1,
            this.time,
            this.Column2,
            this.Progress,
            this.Column3});
            this.dataGridView1.ContextMenuStrip = this.contextMenuStrip1;
            this.dataGridView1.EditMode = System.Windows.Forms.DataGridViewEditMode.EditProgrammatically;
            this.dataGridView1.Location = new System.Drawing.Point(0, 0);
            this.dataGridView1.MultiSelect = false;
            this.dataGridView1.Name = "dataGridView1";
            this.dataGridView1.ReadOnly = true;
            this.dataGridView1.RowHeadersVisible = false;
            this.dataGridView1.SelectionMode = System.Windows.Forms.DataGridViewSelectionMode.FullRowSelect;
            this.dataGridView1.Size = new System.Drawing.Size(871, 169);
            this.dataGridView1.TabIndex = 10;
            // 
            // process
            // 
            this.process.AutoSizeMode = System.Windows.Forms.DataGridViewAutoSizeColumnMode.Fill;
            this.process.FillWeight = 160F;
            this.process.HeaderText = "Name";
            this.process.Name = "process";
            this.process.ReadOnly = true;
            // 
            // Cluster
            // 
            this.Cluster.AutoSizeMode = System.Windows.Forms.DataGridViewAutoSizeColumnMode.Fill;
            this.Cluster.FillWeight = 80F;
            this.Cluster.HeaderText = "Cluster";
            this.Cluster.Name = "Cluster";
            this.Cluster.ReadOnly = true;
            // 
            // Measure
            // 
            this.Measure.AutoSizeMode = System.Windows.Forms.DataGridViewAutoSizeColumnMode.Fill;
            this.Measure.FillWeight = 80F;
            this.Measure.HeaderText = "Measure";
            this.Measure.Name = "Measure";
            this.Measure.ReadOnly = true;
            // 
            // Column1
            // 
            this.Column1.HeaderText = "Dir Name";
            this.Column1.Name = "Column1";
            this.Column1.ReadOnly = true;
            // 
            // time
            // 
            this.time.AutoSizeMode = System.Windows.Forms.DataGridViewAutoSizeColumnMode.Fill;
            this.time.HeaderText = "time";
            this.time.Name = "time";
            this.time.ReadOnly = true;
            // 
            // Column2
            // 
            this.Column2.HeaderText = "Peek Memory";
            this.Column2.Name = "Column2";
            this.Column2.ReadOnly = true;
            // 
            // Progress
            // 
            this.Progress.AutoSizeMode = System.Windows.Forms.DataGridViewAutoSizeColumnMode.Fill;
            this.Progress.HeaderText = "Progress [%]";
            this.Progress.Name = "Progress";
            this.Progress.ProgressBarColor = System.Drawing.Color.Gray;
            this.Progress.ReadOnly = true;
            // 
            // Column3
            // 
            this.Column3.HeaderText = "Vis";
            this.Column3.Name = "Column3";
            this.Column3.ReadOnly = true;
            this.Column3.Visible = false;
            // 
            // contextMenuStrip1
            // 
            this.contextMenuStrip1.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.gotoProcessToolStripMenuItem,
            this.saveCurrentinTextModeToolStripMenuItem,
            this.saveAllbinaryToolStripMenuItem,
            this.loadbinaryToolStripMenuItem,
            this.finishProcessToolStripMenuItem,
            this.closeAllToolStripMenuItem});
            this.contextMenuStrip1.Name = "contextMenuStrip1";
            this.contextMenuStrip1.Size = new System.Drawing.Size(217, 136);
            this.contextMenuStrip1.Opening += new System.ComponentModel.CancelEventHandler(this.contextMenuStrip1_Opening);
            // 
            // gotoProcessToolStripMenuItem
            // 
            this.gotoProcessToolStripMenuItem.Name = "gotoProcessToolStripMenuItem";
            this.gotoProcessToolStripMenuItem.Size = new System.Drawing.Size(216, 22);
            this.gotoProcessToolStripMenuItem.Text = "Show results";
            this.gotoProcessToolStripMenuItem.Click += new System.EventHandler(this.gotoProcessToolStripMenuItem_Click);
            // 
            // saveCurrentinTextModeToolStripMenuItem
            // 
            this.saveCurrentinTextModeToolStripMenuItem.Name = "saveCurrentinTextModeToolStripMenuItem";
            this.saveCurrentinTextModeToolStripMenuItem.Size = new System.Drawing.Size(216, 22);
            this.saveCurrentinTextModeToolStripMenuItem.Text = "Save current (in text mode)";
            this.saveCurrentinTextModeToolStripMenuItem.Click += new System.EventHandler(this.saveCurrentinTextModeToolStripMenuItem_Click);
            // 
            // finishProcessToolStripMenuItem
            // 
            this.finishProcessToolStripMenuItem.Name = "finishProcessToolStripMenuItem";
            this.finishProcessToolStripMenuItem.Size = new System.Drawing.Size(216, 22);
            this.finishProcessToolStripMenuItem.Text = "Close current";
            this.finishProcessToolStripMenuItem.Click += new System.EventHandler(this.finishProcessToolStripMenuItem_Click);
            // 
            // closeAllToolStripMenuItem
            // 
            this.closeAllToolStripMenuItem.Name = "closeAllToolStripMenuItem";
            this.closeAllToolStripMenuItem.Size = new System.Drawing.Size(216, 22);
            this.closeAllToolStripMenuItem.Text = "Close All";
            this.closeAllToolStripMenuItem.Click += new System.EventHandler(this.closeAllToolStripMenuItem_Click);
            // 
            // progressColumn1
            // 
            this.progressColumn1.AutoSizeMode = System.Windows.Forms.DataGridViewAutoSizeColumnMode.Fill;
            this.progressColumn1.HeaderText = "Progress [%]";
            this.progressColumn1.Name = "progressColumn1";
            this.progressColumn1.ProgressBarColor = System.Drawing.Color.Gray;
            // 
            // button1
            // 
            this.button1.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Right)));
            this.button1.Enabled = false;
            this.button1.Location = new System.Drawing.Point(784, 175);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(75, 23);
            this.button1.TabIndex = 11;
            this.button1.Text = "Show";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // button2
            // 
            this.button2.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Right)));
            this.button2.Enabled = false;
            this.button2.Location = new System.Drawing.Point(674, 175);
            this.button2.Name = "button2";
            this.button2.Size = new System.Drawing.Size(75, 23);
            this.button2.TabIndex = 12;
            this.button2.Text = "Save";
            this.button2.UseVisualStyleBackColor = true;
            this.button2.Click += new System.EventHandler(this.button2_Click);
            // 
            // saveAllbinaryToolStripMenuItem
            // 
            this.saveAllbinaryToolStripMenuItem.Name = "saveAllbinaryToolStripMenuItem";
            this.saveAllbinaryToolStripMenuItem.Size = new System.Drawing.Size(216, 22);
            this.saveAllbinaryToolStripMenuItem.Text = "Save all (binary)";
            this.saveAllbinaryToolStripMenuItem.Click += new System.EventHandler(this.saveAllbinaryToolStripMenuItem_Click);
            // 
            // loadbinaryToolStripMenuItem
            // 
            this.loadbinaryToolStripMenuItem.Name = "loadbinaryToolStripMenuItem";
            this.loadbinaryToolStripMenuItem.Size = new System.Drawing.Size(216, 22);
            this.loadbinaryToolStripMenuItem.Text = "Load (binary)";
            this.loadbinaryToolStripMenuItem.Click += new System.EventHandler(this.loadbinaryToolStripMenuItem_Click);
            // 
            // openFileDialog1
            // 
            this.openFileDialog1.FileName = "openFileDialog1";
            // 
            // ResultWindow
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(871, 203);
            this.ControlBox = false;
            this.Controls.Add(this.button2);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.dataGridView1);
            this.Name = "ResultWindow";
            this.Text = "Results window";
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).EndInit();
            this.contextMenuStrip1.ResumeLayout(false);
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.DataGridView dataGridView1;
        private System.Windows.Forms.ContextMenuStrip contextMenuStrip1;
        private System.Windows.Forms.ToolStripMenuItem gotoProcessToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem saveCurrentinTextModeToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem finishProcessToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem closeAllToolStripMenuItem;
        private Graph.ProgressColumn progressColumn1;
        private System.Windows.Forms.SaveFileDialog saveFileDialog1;
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.Button button2;
        private System.Windows.Forms.DataGridViewTextBoxColumn process;
        private System.Windows.Forms.DataGridViewTextBoxColumn Cluster;
        private System.Windows.Forms.DataGridViewTextBoxColumn Measure;
        private System.Windows.Forms.DataGridViewTextBoxColumn Column1;
        private System.Windows.Forms.DataGridViewTextBoxColumn time;
        private System.Windows.Forms.DataGridViewTextBoxColumn Column2;
        private Graph.ProgressColumn Progress;
        private System.Windows.Forms.DataGridViewTextBoxColumn Column3;
        private System.Windows.Forms.ToolStripMenuItem saveAllbinaryToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem loadbinaryToolStripMenuItem;
        private System.Windows.Forms.OpenFileDialog openFileDialog1;
    }
}