namespace WorkFlows
{
    partial class RpartSimple
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
            System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(RpartSimple));
            this.folderBrowserDialog1 = new System.Windows.Forms.FolderBrowserDialog();
            this.relevantC = new System.Windows.Forms.NumericUpDown();
            this.percentData = new System.Windows.Forms.NumericUpDown();
            this.label7 = new System.Windows.Forms.Label();
            this.label6 = new System.Windows.Forms.Label();
            this.button2 = new System.Windows.Forms.Button();
            this.label5 = new System.Windows.Forms.Label();
            this.label4 = new System.Windows.Forms.Label();
            this.openFileDialog1 = new System.Windows.Forms.OpenFileDialog();
            this.refPoints = new System.Windows.Forms.NumericUpDown();
            this.label10 = new System.Windows.Forms.Label();
            this.selectReference = new System.Windows.Forms.CheckBox();
            this.button4 = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.relevantC)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.percentData)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.refPoints)).BeginInit();
            this.SuspendLayout();
            // 
            // relevantC
            // 
            this.relevantC.Location = new System.Drawing.Point(207, 34);
            this.relevantC.Maximum = new decimal(new int[] {
            100000000,
            0,
            0,
            0});
            this.relevantC.Name = "relevantC";
            this.relevantC.Size = new System.Drawing.Size(120, 20);
            this.relevantC.TabIndex = 64;
            this.relevantC.Value = new decimal(new int[] {
            5,
            0,
            0,
            0});
            // 
            // percentData
            // 
            this.percentData.Location = new System.Drawing.Point(207, 69);
            this.percentData.Name = "percentData";
            this.percentData.Size = new System.Drawing.Size(120, 20);
            this.percentData.TabIndex = 63;
            this.percentData.Value = new decimal(new int[] {
            30,
            0,
            0,
            0});
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(9, 71);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(186, 13);
            this.label7.TabIndex = 62;
            this.label7.Text = "Percent of data in relevant clusters (F)";
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(9, 34);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(152, 13);
            this.label6.TabIndex = 61;
            this.label6.Text = "Number of relevant clusters (K)";
            // 
            // button2
            // 
            this.button2.Anchor = System.Windows.Forms.AnchorStyles.Bottom;
            this.button2.Location = new System.Drawing.Point(569, 167);
            this.button2.Name = "button2";
            this.button2.Size = new System.Drawing.Size(53, 36);
            this.button2.TabIndex = 65;
            this.button2.Text = "RUN";
            this.button2.UseVisualStyleBackColor = true;
            this.button2.Click += new System.EventHandler(this.button2_Click);
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(201, 96);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(65, 13);
            this.label5.TabIndex = 69;
            this.label5.Text = "Entropy filter";
            this.label5.Visible = false;
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(9, 96);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(92, 13);
            this.label4.TabIndex = 68;
            this.label4.Text = "Key reduction by: ";
            this.label4.Visible = false;
            // 
            // openFileDialog1
            // 
            this.openFileDialog1.FileName = "openFileDialog1";
            // 
            // refPoints
            // 
            this.refPoints.Location = new System.Drawing.Point(511, 34);
            this.refPoints.Name = "refPoints";
            this.refPoints.Size = new System.Drawing.Size(120, 20);
            this.refPoints.TabIndex = 73;
            this.refPoints.Value = new decimal(new int[] {
            20,
            0,
            0,
            0});
            // 
            // label10
            // 
            this.label10.AutoSize = true;
            this.label10.Location = new System.Drawing.Point(370, 36);
            this.label10.Name = "label10";
            this.label10.Size = new System.Drawing.Size(135, 13);
            this.label10.TabIndex = 74;
            this.label10.Text = "Number of referance points";
            // 
            // selectReference
            // 
            this.selectReference.AutoSize = true;
            this.selectReference.Checked = true;
            this.selectReference.CheckState = System.Windows.Forms.CheckState.Checked;
            this.selectReference.Location = new System.Drawing.Point(13, 124);
            this.selectReference.Name = "selectReference";
            this.selectReference.Size = new System.Drawing.Size(170, 17);
            this.selectReference.TabIndex = 76;
            this.selectReference.Text = "Find reference for each cluster";
            this.selectReference.UseVisualStyleBackColor = true;
            // 
            // button4
            // 
            this.button4.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.button4.Image = ((System.Drawing.Image)(resources.GetObject("button4.Image")));
            this.button4.Location = new System.Drawing.Point(12, 167);
            this.button4.Name = "button4";
            this.button4.Size = new System.Drawing.Size(53, 36);
            this.button4.TabIndex = 4;
            this.button4.UseVisualStyleBackColor = true;
            this.button4.Click += new System.EventHandler(this.button4_Click);
            // 
            // RpartSimple
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(646, 215);
            this.Controls.Add(this.selectReference);
            this.Controls.Add(this.label10);
            this.Controls.Add(this.refPoints);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.button2);
            this.Controls.Add(this.relevantC);
            this.Controls.Add(this.percentData);
            this.Controls.Add(this.label7);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.button4);
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedSingle;
            this.Icon = ((System.Drawing.Icon)(resources.GetObject("$this.Icon")));
            this.MaximizeBox = false;
            this.MinimizeBox = false;
            this.Name = "RpartSimple";
            this.Text = "RpartSimple";
            this.FormClosed += new System.Windows.Forms.FormClosedEventHandler(this.RpartSimple_FormClosed);
            ((System.ComponentModel.ISupportInitialize)(this.relevantC)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.percentData)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.refPoints)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        public System.Windows.Forms.Button button4;
        private System.Windows.Forms.FolderBrowserDialog folderBrowserDialog1;
        public System.Windows.Forms.NumericUpDown relevantC;
        public System.Windows.Forms.NumericUpDown percentData;
        public System.Windows.Forms.Label label7;
        public System.Windows.Forms.Label label6;
        public System.Windows.Forms.Button button2;
        public System.Windows.Forms.Label label5;
        public System.Windows.Forms.Label label4;
        public System.Windows.Forms.OpenFileDialog openFileDialog1;
        public System.Windows.Forms.NumericUpDown refPoints;
        public  System.Windows.Forms.Label label10;
        private System.Windows.Forms.CheckBox selectReference;
    }
}