
namespace Graph
{
    partial class DataForm
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
            this.dataGridView1 = new System.Windows.Forms.DataGridView();
            this.Column1 = new System.Windows.Forms.DataGridViewTextBoxColumn();
            this.Tab = new System.Windows.Forms.CheckBox();
            this.Comma = new System.Windows.Forms.CheckBox();
            this.Semicolon = new System.Windows.Forms.CheckBox();
            this.label1 = new System.Windows.Forms.Label();
            this.label2 = new System.Windows.Forms.Label();
            this.Row = new System.Windows.Forms.NumericUpDown();
            this.label3 = new System.Windows.Forms.Label();
            this.label4 = new System.Windows.Forms.Label();
            this.Column = new System.Windows.Forms.NumericUpDown();
            this.Space = new System.Windows.Forms.CheckBox();
            this.label5 = new System.Windows.Forms.Label();
            this.label6 = new System.Windows.Forms.Label();
            this.label7 = new System.Windows.Forms.Label();
            this.genePos = new System.Windows.Forms.NumericUpDown();
            this.samplePos = new System.Windows.Forms.NumericUpDown();
            this.Import = new System.Windows.Forms.Button();
            this.button1 = new System.Windows.Forms.Button();
            this.geneLabels = new System.Windows.Forms.ComboBox();
            this.sampleLabels = new System.Windows.Forms.ComboBox();
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.Row)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.Column)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.genePos)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.samplePos)).BeginInit();
            this.SuspendLayout();
            // 
            // dataGridView1
            // 
            this.dataGridView1.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.dataGridView1.Columns.AddRange(new System.Windows.Forms.DataGridViewColumn[] {
            this.Column1});
            this.dataGridView1.Location = new System.Drawing.Point(2, 316);
            this.dataGridView1.Name = "dataGridView1";
            this.dataGridView1.Size = new System.Drawing.Size(507, 132);
            this.dataGridView1.TabIndex = 0;
            // 
            // Column1
            // 
            this.Column1.AutoSizeMode = System.Windows.Forms.DataGridViewAutoSizeColumnMode.Fill;
            this.Column1.HeaderText = "Column1";
            this.Column1.MinimumWidth = 100;
            this.Column1.Name = "Column1";
            // 
            // Tab
            // 
            this.Tab.AutoSize = true;
            this.Tab.Checked = true;
            this.Tab.CheckState = System.Windows.Forms.CheckState.Checked;
            this.Tab.Location = new System.Drawing.Point(21, 62);
            this.Tab.Name = "Tab";
            this.Tab.Size = new System.Drawing.Size(45, 17);
            this.Tab.TabIndex = 3;
            this.Tab.Text = "Tab";
            this.Tab.UseVisualStyleBackColor = true;
            this.Tab.CheckedChanged += new System.EventHandler(this.Tab_CheckedChanged);
            // 
            // Comma
            // 
            this.Comma.AutoSize = true;
            this.Comma.Location = new System.Drawing.Point(86, 62);
            this.Comma.Name = "Comma";
            this.Comma.Size = new System.Drawing.Size(61, 17);
            this.Comma.TabIndex = 4;
            this.Comma.Text = "Comma";
            this.Comma.UseVisualStyleBackColor = true;
            // 
            // Semicolon
            // 
            this.Semicolon.AutoSize = true;
            this.Semicolon.Location = new System.Drawing.Point(174, 62);
            this.Semicolon.Name = "Semicolon";
            this.Semicolon.Size = new System.Drawing.Size(75, 17);
            this.Semicolon.TabIndex = 5;
            this.Semicolon.Text = "Semicolon";
            this.Semicolon.UseVisualStyleBackColor = true;
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Font = new System.Drawing.Font("Microsoft Sans Serif", 12F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(238)));
            this.label1.Location = new System.Drawing.Point(18, 23);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(153, 20);
            this.label1.TabIndex = 6;
            this.label1.Text = "Separator options";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Font = new System.Drawing.Font("Microsoft Sans Serif", 12F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(238)));
            this.label2.Location = new System.Drawing.Point(19, 93);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(151, 20);
            this.label2.TabIndex = 7;
            this.label2.Text = "Import Data From";
            // 
            // Row
            // 
            this.Row.Location = new System.Drawing.Point(86, 127);
            this.Row.Name = "Row";
            this.Row.Size = new System.Drawing.Size(76, 20);
            this.Row.TabIndex = 8;
            this.Row.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.Row.ValueChanged += new System.EventHandler(this.Row_ValueChanged);
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(18, 134);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(29, 13);
            this.label3.TabIndex = 9;
            this.label3.Text = "Row";
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(273, 134);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(42, 13);
            this.label4.TabIndex = 10;
            this.label4.Text = "Column";
            // 
            // Column
            // 
            this.Column.Location = new System.Drawing.Point(347, 132);
            this.Column.Name = "Column";
            this.Column.Size = new System.Drawing.Size(76, 20);
            this.Column.TabIndex = 11;
            this.Column.Value = new decimal(new int[] {
            1,
            0,
            0,
            0});
            this.Column.ValueChanged += new System.EventHandler(this.Column_ValueChanged);
            // 
            // Space
            // 
            this.Space.AutoSize = true;
            this.Space.Location = new System.Drawing.Point(276, 62);
            this.Space.Name = "Space";
            this.Space.Size = new System.Drawing.Size(57, 17);
            this.Space.TabIndex = 12;
            this.Space.Text = "Space";
            this.Space.UseVisualStyleBackColor = true;
            this.Space.CheckedChanged += new System.EventHandler(this.Space_CheckedChanged);
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Font = new System.Drawing.Font("Microsoft Sans Serif", 12F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(238)));
            this.label5.Location = new System.Drawing.Point(17, 162);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(202, 20);
            this.label5.TabIndex = 13;
            this.label5.Text = "Gene and sample labels";
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(18, 192);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(31, 13);
            this.label6.TabIndex = 16;
            this.label6.Text = "gene";
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(18, 223);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(40, 13);
            this.label7.TabIndex = 17;
            this.label7.Text = "sample";
            // 
            // genePos
            // 
            this.genePos.Location = new System.Drawing.Point(236, 190);
            this.genePos.Name = "genePos";
            this.genePos.Size = new System.Drawing.Size(63, 20);
            this.genePos.TabIndex = 18;
            this.genePos.ValueChanged += new System.EventHandler(this.genePos_ValueChanged);
            // 
            // samplePos
            // 
            this.samplePos.Location = new System.Drawing.Point(236, 221);
            this.samplePos.Name = "samplePos";
            this.samplePos.Size = new System.Drawing.Size(63, 20);
            this.samplePos.TabIndex = 19;
            this.samplePos.ValueChanged += new System.EventHandler(this.samplePos_ValueChanged);
            // 
            // Import
            // 
            this.Import.Font = new System.Drawing.Font("Microsoft Sans Serif", 9F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(238)));
            this.Import.Location = new System.Drawing.Point(23, 467);
            this.Import.Name = "Import";
            this.Import.Size = new System.Drawing.Size(75, 23);
            this.Import.TabIndex = 20;
            this.Import.Text = "Import";
            this.Import.UseVisualStyleBackColor = true;
            this.Import.Click += new System.EventHandler(this.Import_Click);
            // 
            // button1
            // 
            this.button1.Font = new System.Drawing.Font("Microsoft Sans Serif", 9F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(238)));
            this.button1.Location = new System.Drawing.Point(447, 467);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(75, 23);
            this.button1.TabIndex = 21;
            this.button1.Text = "Cancel";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // geneLabels
            // 
            this.geneLabels.FormattingEnabled = true;
            this.geneLabels.Items.AddRange(new object[] {
            "Row",
            "Column"});
            this.geneLabels.Location = new System.Drawing.Point(77, 189);
            this.geneLabels.Name = "geneLabels";
            this.geneLabels.Size = new System.Drawing.Size(121, 21);
            this.geneLabels.TabIndex = 22;
            this.geneLabels.SelectedIndexChanged += new System.EventHandler(this.geneLabels_SelectedIndexChanged);
            // 
            // sampleLabels
            // 
            this.sampleLabels.FormattingEnabled = true;
            this.sampleLabels.Items.AddRange(new object[] {
            "Row",
            "Column"});
            this.sampleLabels.Location = new System.Drawing.Point(77, 220);
            this.sampleLabels.Name = "sampleLabels";
            this.sampleLabels.Size = new System.Drawing.Size(121, 21);
            this.sampleLabels.TabIndex = 23;
            this.sampleLabels.SelectedIndexChanged += new System.EventHandler(this.sampleLabels_SelectedIndexChanged);
            // 
            // DataForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(556, 502);
            this.Controls.Add(this.sampleLabels);
            this.Controls.Add(this.geneLabels);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.Import);
            this.Controls.Add(this.samplePos);
            this.Controls.Add(this.genePos);
            this.Controls.Add(this.label7);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.Space);
            this.Controls.Add(this.Column);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.Row);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.Semicolon);
            this.Controls.Add(this.Comma);
            this.Controls.Add(this.Tab);
            this.Controls.Add(this.dataGridView1);
            this.Name = "DataForm";
            this.Text = "DataForm";
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.Row)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.Column)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.genePos)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.samplePos)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.DataGridView dataGridView1;
        private System.Windows.Forms.CheckBox Tab;
        private System.Windows.Forms.CheckBox Comma;
        private System.Windows.Forms.CheckBox Semicolon;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.NumericUpDown Row;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.NumericUpDown Column;
        private System.Windows.Forms.CheckBox Space;
        private System.Windows.Forms.DataGridViewTextBoxColumn Column1;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.NumericUpDown genePos;
        private System.Windows.Forms.NumericUpDown samplePos;
        private System.Windows.Forms.Button Import;
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.ComboBox geneLabels;
        private System.Windows.Forms.ComboBox sampleLabels;
    }
}