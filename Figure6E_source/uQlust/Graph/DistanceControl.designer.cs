namespace Graph
{
    partial class DistanceControl
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

        #region Component Designer generated code

        /// <summary> 
        /// Required method for Designer support - do not modify 
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.groupBox3 = new System.Windows.Forms.GroupBox();
            this.radioEucl = new System.Windows.Forms.RadioButton();
            this.radioCosine = new System.Windows.Forms.RadioButton();
            this.label2 = new System.Windows.Forms.Label();
            this.radioPearson = new System.Windows.Forms.RadioButton();
            this.radio1DJury = new System.Windows.Forms.RadioButton();
            this.referenceBox = new System.Windows.Forms.CheckBox();
            this.fastHamming = new System.Windows.Forms.RadioButton();
            this.groupBox3.SuspendLayout();
            this.SuspendLayout();
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
            this.groupBox3.Location = new System.Drawing.Point(7, 6);
            this.groupBox3.Name = "groupBox3";
            this.groupBox3.Size = new System.Drawing.Size(541, 54);
            this.groupBox3.TabIndex = 0;
            this.groupBox3.TabStop = false;
            this.groupBox3.Text = "Distance measures";
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
            this.radioCosine.CheckedChanged += new System.EventHandler(this.radio1DJury_CheckedChanged);
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
            this.radioPearson.CheckedChanged += new System.EventHandler(this.radioRmsd_CheckedChanged);
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
            this.radio1DJury.CheckedChanged += new System.EventHandler(this.radio1DJury_CheckedChanged);
            // 
            // referenceBox
            // 
            this.referenceBox.AutoSize = true;
            this.referenceBox.Location = new System.Drawing.Point(7, 83);
            this.referenceBox.Name = "referenceBox";
            this.referenceBox.Size = new System.Drawing.Size(115, 30);
            this.referenceBox.TabIndex = 7;
            this.referenceBox.Text = "Use 1DJury to find\r\nreference structure";
            this.referenceBox.UseVisualStyleBackColor = true;
            this.referenceBox.CheckedChanged += new System.EventHandler(this.referenceBox_CheckedChanged);
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
            // DistanceControl
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.Controls.Add(this.referenceBox);
            this.Controls.Add(this.groupBox3);
            this.Name = "DistanceControl";
            this.Size = new System.Drawing.Size(557, 141);
            this.groupBox3.ResumeLayout(false);
            this.groupBox3.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.GroupBox groupBox3;
        private System.Windows.Forms.RadioButton radioPearson;
        private System.Windows.Forms.RadioButton radio1DJury;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.CheckBox referenceBox;
        private System.Windows.Forms.RadioButton radioCosine;
        private System.Windows.Forms.RadioButton radioEucl;
        private System.Windows.Forms.RadioButton fastHamming;
    }
}
