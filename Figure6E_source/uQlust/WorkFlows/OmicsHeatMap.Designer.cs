namespace WorkFlows
{
    partial class OmicsHeatMap
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
            System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(OmicsHeatMap));
            this.radioButton1 = new System.Windows.Forms.RadioButton();
            this.relevantC = new System.Windows.Forms.NumericUpDown();
            this.label11 = new System.Windows.Forms.Label();
            this.label12 = new System.Windows.Forms.Label();
            this.numericUpDown4 = new System.Windows.Forms.NumericUpDown();
            this.button4 = new System.Windows.Forms.Button();
            this.button1 = new System.Windows.Forms.Button();
            this.openFileDialog1 = new System.Windows.Forms.OpenFileDialog();
            this.consensus = new System.Windows.Forms.CheckBox();
            this.label4 = new System.Windows.Forms.Label();
            this.numericUpDown1 = new System.Windows.Forms.NumericUpDown();
            this.groupBox3 = new System.Windows.Forms.GroupBox();
            this.fastHamming = new System.Windows.Forms.RadioButton();
            this.radioEucl = new System.Windows.Forms.RadioButton();
            this.radioCosine = new System.Windows.Forms.RadioButton();
            this.label2 = new System.Windows.Forms.Label();
            this.radioPearson = new System.Windows.Forms.RadioButton();
            this.radio1DJury = new System.Windows.Forms.RadioButton();
            this.radioButton3 = new System.Windows.Forms.RadioButton();
            this.radioButton2 = new System.Windows.Forms.RadioButton();
            this.button2 = new System.Windows.Forms.Button();
            this.textBox1 = new System.Windows.Forms.TextBox();
            this.button3 = new System.Windows.Forms.Button();
            this.label1 = new System.Windows.Forms.Label();
            ((System.ComponentModel.ISupportInitialize)(this.relevantC)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.numericUpDown4)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.numericUpDown1)).BeginInit();
            this.groupBox3.SuspendLayout();
            this.SuspendLayout();
            // 
            // radioButton1
            // 
            this.radioButton1.AutoSize = true;
            this.radioButton1.Checked = true;
            this.radioButton1.Location = new System.Drawing.Point(12, 51);
            this.radioButton1.Name = "radioButton1";
            this.radioButton1.Size = new System.Drawing.Size(51, 17);
            this.radioButton1.TabIndex = 83;
            this.radioButton1.TabStop = true;
            this.radioButton1.Text = "Rpart";
            this.radioButton1.UseVisualStyleBackColor = true;
            // 
            // relevantC
            // 
            this.relevantC.Location = new System.Drawing.Point(199, 23);
            this.relevantC.Maximum = new decimal(new int[] {
            100000000,
            0,
            0,
            0});
            this.relevantC.Name = "relevantC";
            this.relevantC.Size = new System.Drawing.Size(120, 20);
            this.relevantC.TabIndex = 79;
            this.relevantC.Value = new decimal(new int[] {
            1000,
            0,
            0,
            0});
            // 
            // label11
            // 
            this.label11.AutoSize = true;
            this.label11.Location = new System.Drawing.Point(4, 25);
            this.label11.Name = "label11";
            this.label11.Size = new System.Drawing.Size(161, 13);
            this.label11.TabIndex = 78;
            this.label11.Text = "Number of required rows clusters";
            // 
            // label12
            // 
            this.label12.AutoSize = true;
            this.label12.Location = new System.Drawing.Point(359, 25);
            this.label12.Name = "label12";
            this.label12.Size = new System.Drawing.Size(178, 13);
            this.label12.TabIndex = 86;
            this.label12.Text = "Number of required columns clusters";
            // 
            // numericUpDown4
            // 
            this.numericUpDown4.Location = new System.Drawing.Point(543, 23);
            this.numericUpDown4.Maximum = new decimal(new int[] {
            100000000,
            0,
            0,
            0});
            this.numericUpDown4.Name = "numericUpDown4";
            this.numericUpDown4.Size = new System.Drawing.Size(120, 20);
            this.numericUpDown4.TabIndex = 87;
            this.numericUpDown4.Value = new decimal(new int[] {
            1000,
            0,
            0,
            0});
            // 
            // button4
            // 
            this.button4.Image = ((System.Drawing.Image)(resources.GetObject("button4.Image")));
            this.button4.Location = new System.Drawing.Point(7, 171);
            this.button4.Name = "button4";
            this.button4.Size = new System.Drawing.Size(53, 36);
            this.button4.TabIndex = 88;
            this.button4.UseVisualStyleBackColor = true;
            this.button4.Click += new System.EventHandler(this.button4_Click);
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(610, 171);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(53, 36);
            this.button1.TabIndex = 89;
            this.button1.Text = "RUN";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // openFileDialog1
            // 
            this.openFileDialog1.FileName = "openFileDialog1";
            // 
            // consensus
            // 
            this.consensus.AutoSize = true;
            this.consensus.Checked = true;
            this.consensus.CheckState = System.Windows.Forms.CheckState.Checked;
            this.consensus.Location = new System.Drawing.Point(228, 50);
            this.consensus.Name = "consensus";
            this.consensus.Size = new System.Drawing.Size(143, 17);
            this.consensus.TabIndex = 91;
            this.consensus.Text = "Use consenus projection";
            this.consensus.UseVisualStyleBackColor = true;
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(389, 50);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(135, 13);
            this.label4.TabIndex = 93;
            this.label4.Text = "Number of reference points";
            this.label4.Visible = false;
            // 
            // numericUpDown1
            // 
            this.numericUpDown1.Location = new System.Drawing.Point(543, 48);
            this.numericUpDown1.Name = "numericUpDown1";
            this.numericUpDown1.Size = new System.Drawing.Size(120, 20);
            this.numericUpDown1.TabIndex = 92;
            this.numericUpDown1.Value = new decimal(new int[] {
            20,
            0,
            0,
            0});
            this.numericUpDown1.Visible = false;
            // 
            // groupBox3
            // 
            this.groupBox3.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.groupBox3.Controls.Add(this.fastHamming);
            this.groupBox3.Controls.Add(this.radioEucl);
            this.groupBox3.Controls.Add(this.radioCosine);
            this.groupBox3.Controls.Add(this.label2);
            this.groupBox3.Controls.Add(this.radioPearson);
            this.groupBox3.Controls.Add(this.radio1DJury);
            this.groupBox3.Location = new System.Drawing.Point(7, 111);
            this.groupBox3.Name = "groupBox3";
            this.groupBox3.Size = new System.Drawing.Size(541, 54);
            this.groupBox3.TabIndex = 95;
            this.groupBox3.TabStop = false;
            this.groupBox3.Text = "Distance measures";
            // 
            // fastHamming
            // 
            this.fastHamming.AutoSize = true;
            this.fastHamming.Location = new System.Drawing.Point(409, 19);
            this.fastHamming.Name = "fastHamming";
            this.fastHamming.Size = new System.Drawing.Size(92, 17);
            this.fastHamming.TabIndex = 26;
            this.fastHamming.TabStop = true;
            this.fastHamming.Text = "Fast Hamming";
            this.fastHamming.UseVisualStyleBackColor = true;
            // 
            // radioEucl
            // 
            this.radioEucl.AutoSize = true;
            this.radioEucl.Location = new System.Drawing.Point(302, 19);
            this.radioEucl.Name = "radioEucl";
            this.radioEucl.Size = new System.Drawing.Size(72, 17);
            this.radioEucl.TabIndex = 25;
            this.radioEucl.TabStop = true;
            this.radioEucl.Text = "Euclidean";
            this.radioEucl.UseVisualStyleBackColor = true;
            // 
            // radioCosine
            // 
            this.radioCosine.AutoSize = true;
            this.radioCosine.Location = new System.Drawing.Point(207, 19);
            this.radioCosine.Name = "radioCosine";
            this.radioCosine.Size = new System.Drawing.Size(57, 17);
            this.radioCosine.TabIndex = 24;
            this.radioCosine.TabStop = true;
            this.radioCosine.Text = "Cosine";
            this.radioCosine.UseVisualStyleBackColor = true;
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(73, 67);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(0, 13);
            this.label2.TabIndex = 23;
            // 
            // radioPearson
            // 
            this.radioPearson.AutoSize = true;
            this.radioPearson.Location = new System.Drawing.Point(114, 19);
            this.radioPearson.Name = "radioPearson";
            this.radioPearson.Size = new System.Drawing.Size(64, 17);
            this.radioPearson.TabIndex = 2;
            this.radioPearson.TabStop = true;
            this.radioPearson.Text = "Pearson";
            this.radioPearson.UseVisualStyleBackColor = true;
            // 
            // radio1DJury
            // 
            this.radio1DJury.AutoSize = true;
            this.radio1DJury.Checked = true;
            this.radio1DJury.Location = new System.Drawing.Point(10, 19);
            this.radio1DJury.Name = "radio1DJury";
            this.radio1DJury.Size = new System.Drawing.Size(69, 17);
            this.radio1DJury.TabIndex = 1;
            this.radio1DJury.TabStop = true;
            this.radio1DJury.Text = "Hamming";
            this.radio1DJury.UseVisualStyleBackColor = true;
            // 
            // radioButton3
            // 
            this.radioButton3.AutoSize = true;
            this.radioButton3.Location = new System.Drawing.Point(83, 51);
            this.radioButton3.Name = "radioButton3";
            this.radioButton3.Size = new System.Drawing.Size(45, 17);
            this.radioButton3.TabIndex = 96;
            this.radioButton3.TabStop = true;
            this.radioButton3.Text = "Fast";
            this.radioButton3.UseVisualStyleBackColor = true;
            // 
            // radioButton2
            // 
            this.radioButton2.AutoSize = true;
            this.radioButton2.Location = new System.Drawing.Point(137, 51);
            this.radioButton2.Name = "radioButton2";
            this.radioButton2.Size = new System.Drawing.Size(63, 17);
            this.radioButton2.TabIndex = 97;
            this.radioButton2.TabStop = true;
            this.radioButton2.Text = "kMeans";
            this.radioButton2.UseVisualStyleBackColor = true;
            // 
            // button2
            // 
            this.button2.Location = new System.Drawing.Point(306, 83);
            this.button2.Name = "button2";
            this.button2.Size = new System.Drawing.Size(75, 23);
            this.button2.TabIndex = 98;
            this.button2.Text = "Cluster File";
            this.button2.UseVisualStyleBackColor = true;
            this.button2.Click += new System.EventHandler(this.button2_Click);
            // 
            // textBox1
            // 
            this.textBox1.Location = new System.Drawing.Point(17, 85);
            this.textBox1.Name = "textBox1";
            this.textBox1.Size = new System.Drawing.Size(266, 20);
            this.textBox1.TabIndex = 99;
            // 
            // button3
            // 
            this.button3.Location = new System.Drawing.Point(421, 85);
            this.button3.Name = "button3";
            this.button3.Size = new System.Drawing.Size(87, 23);
            this.button3.TabIndex = 100;
            this.button3.Text = "GeneSelection";
            this.button3.UseVisualStyleBackColor = true;
            this.button3.Click += new System.EventHandler(this.button3_Click);
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(528, 92);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(35, 13);
            this.label1.TabIndex = 101;
            this.label1.Text = "label1";
            // 
            // OmicsHeatMap
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(670, 218);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.button3);
            this.Controls.Add(this.textBox1);
            this.Controls.Add(this.button2);
            this.Controls.Add(this.radioButton2);
            this.Controls.Add(this.radioButton3);
            this.Controls.Add(this.groupBox3);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.numericUpDown1);
            this.Controls.Add(this.consensus);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.button4);
            this.Controls.Add(this.numericUpDown4);
            this.Controls.Add(this.label12);
            this.Controls.Add(this.radioButton1);
            this.Controls.Add(this.relevantC);
            this.Controls.Add(this.label11);
            this.Name = "OmicsHeatMap";
            this.Text = "HeatMap";
            this.FormClosed += new System.Windows.Forms.FormClosedEventHandler(this.OmicsHeatMap_FormClosed);
            ((System.ComponentModel.ISupportInitialize)(this.relevantC)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.numericUpDown4)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.numericUpDown1)).EndInit();
            this.groupBox3.ResumeLayout(false);
            this.groupBox3.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion
        private System.Windows.Forms.RadioButton radioButton1;
        private System.Windows.Forms.NumericUpDown relevantC;
        private System.Windows.Forms.Label label11;
        private System.Windows.Forms.Label label12;
        private System.Windows.Forms.NumericUpDown numericUpDown4;
        public System.Windows.Forms.Button button4;
        public System.Windows.Forms.Button button1;
        private System.Windows.Forms.OpenFileDialog openFileDialog1;
        private System.Windows.Forms.CheckBox consensus;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.NumericUpDown numericUpDown1;
        private System.Windows.Forms.GroupBox groupBox3;
        private System.Windows.Forms.RadioButton fastHamming;
        private System.Windows.Forms.RadioButton radioEucl;
        private System.Windows.Forms.RadioButton radioCosine;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.RadioButton radioPearson;
        private System.Windows.Forms.RadioButton radio1DJury;
        private System.Windows.Forms.RadioButton radioButton3;
        private System.Windows.Forms.RadioButton radioButton2;
        private System.Windows.Forms.Button button2;
        private System.Windows.Forms.TextBox textBox1;
        private System.Windows.Forms.Button button3;
        private System.Windows.Forms.Label label1;
    }
}