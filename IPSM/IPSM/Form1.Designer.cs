namespace IPSM
{
    partial class IPSM
    {
        /// <summary>
        /// Variable nécessaire au concepteur.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Nettoyage des ressources utilisées.
        /// </summary>
        /// <param name="disposing">true si les ressources managées doivent être supprimées ; sinon, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Code généré par le Concepteur Windows Form

        /// <summary>
        /// Méthode requise pour la prise en charge du concepteur - ne modifiez pas
        /// le contenu de cette méthode avec l'éditeur de code.
        /// </summary>
        private void InitializeComponent()
        {
            this.pictureZone = new System.Windows.Forms.PictureBox();
            this.log = new System.Windows.Forms.Label();
            this.groupBox1 = new System.Windows.Forms.GroupBox();
            this.textBox1 = new System.Windows.Forms.TextBox();
            this.label2 = new System.Windows.Forms.Label();
            this.label1 = new System.Windows.Forms.Label();
            this.comboBox1 = new System.Windows.Forms.ComboBox();
            this.numberTensorFields = new System.Windows.Forms.NumericUpDown();
            ((System.ComponentModel.ISupportInitialize)(this.pictureZone)).BeginInit();
            this.groupBox1.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.numberTensorFields)).BeginInit();
            this.SuspendLayout();
            // 
            // pictureZone
            // 
            this.pictureZone.Location = new System.Drawing.Point(10, 10);
            this.pictureZone.Name = "pictureZone";
            this.pictureZone.Size = new System.Drawing.Size(52, 47);
            this.pictureZone.SizeMode = System.Windows.Forms.PictureBoxSizeMode.CenterImage;
            this.pictureZone.TabIndex = 1;
            this.pictureZone.TabStop = false;
            this.pictureZone.MouseClick += new System.Windows.Forms.MouseEventHandler(this.MouseClick);
            // 
            // log
            // 
            this.log.AutoSize = true;
            this.log.Location = new System.Drawing.Point(4, 27);
            this.log.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.log.Name = "log";
            this.log.Size = new System.Drawing.Size(90, 13);
            this.log.TabIndex = 2;
            this.log.Text = "Tensor Fields row";
            // 
            // groupBox1
            // 
            this.groupBox1.Controls.Add(this.textBox1);
            this.groupBox1.Controls.Add(this.label2);
            this.groupBox1.Controls.Add(this.label1);
            this.groupBox1.Controls.Add(this.comboBox1);
            this.groupBox1.Controls.Add(this.numberTensorFields);
            this.groupBox1.Controls.Add(this.log);
            this.groupBox1.Dock = System.Windows.Forms.DockStyle.Right;
            this.groupBox1.Location = new System.Drawing.Point(1236, 0);
            this.groupBox1.Margin = new System.Windows.Forms.Padding(2);
            this.groupBox1.Name = "groupBox1";
            this.groupBox1.Padding = new System.Windows.Forms.Padding(2);
            this.groupBox1.Size = new System.Drawing.Size(200, 849);
            this.groupBox1.TabIndex = 3;
            this.groupBox1.TabStop = false;
            this.groupBox1.Text = "groupBox1";
            // 
            // textBox1
            // 
            this.textBox1.Location = new System.Drawing.Point(126, 63);
            this.textBox1.Name = "textBox1";
            this.textBox1.Size = new System.Drawing.Size(42, 20);
            this.textBox1.TabIndex = 7;
            this.textBox1.Text = "90";
            this.textBox1.TextChanged += new System.EventHandler(this.textBox1_TextChanged);
            this.textBox1.KeyDown += new System.Windows.Forms.KeyEventHandler(this.textBox1_KeyDown);
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(4, 63);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(117, 13);
            this.label2.TabIndex = 6;
            this.label2.Text = "Tensor orientation[en °]";
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(7, 107);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(68, 13);
            this.label1.TabIndex = 5;
            this.label1.Text = "Choose view";
            // 
            // comboBox1
            // 
            this.comboBox1.FormattingEnabled = true;
            this.comboBox1.Items.AddRange(new object[] {
            "Tensor field",
            "Stream line",
            "Flow visualization"});
            this.comboBox1.Location = new System.Drawing.Point(90, 107);
            this.comboBox1.Name = "comboBox1";
            this.comboBox1.Size = new System.Drawing.Size(98, 21);
            this.comboBox1.TabIndex = 4;
            this.comboBox1.SelectedIndexChanged += new System.EventHandler(this.comboBox1_SelectedIndexChanged);
            // 
            // numberTensorFields
            // 
            this.numberTensorFields.Location = new System.Drawing.Point(126, 25);
            this.numberTensorFields.Name = "numberTensorFields";
            this.numberTensorFields.Size = new System.Drawing.Size(42, 20);
            this.numberTensorFields.TabIndex = 3;
            this.numberTensorFields.ValueChanged += new System.EventHandler(this.ChangeNumberTensorFieldToDisplay);
            // 
            // IPSM
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.BackColor = System.Drawing.Color.White;
            this.ClientSize = new System.Drawing.Size(1436, 849);
            this.Controls.Add(this.groupBox1);
            this.Controls.Add(this.pictureZone);
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedSingle;
            this.Name = "IPSM";
            this.StartPosition = System.Windows.Forms.FormStartPosition.Manual;
            this.Text = "IPSM";
            this.FormClosing += new System.Windows.Forms.FormClosingEventHandler(this.IPSM_FormClosing);
            this.Load += new System.EventHandler(this.IPSM_Load);
            this.Paint += new System.Windows.Forms.PaintEventHandler(this.IPSM_Paint);
            ((System.ComponentModel.ISupportInitialize)(this.pictureZone)).EndInit();
            this.groupBox1.ResumeLayout(false);
            this.groupBox1.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.numberTensorFields)).EndInit();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.PictureBox pictureZone;
        private System.Windows.Forms.Label log;
        private System.Windows.Forms.GroupBox groupBox1;
        private System.Windows.Forms.NumericUpDown numberTensorFields;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.ComboBox comboBox1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.TextBox textBox1;


    }
}

