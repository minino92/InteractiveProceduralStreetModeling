﻿namespace IPSM
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
            this.button1 = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.pictureZone)).BeginInit();
            this.groupBox1.SuspendLayout();
            this.SuspendLayout();
            // 
            // pictureZone
            // 
            this.pictureZone.Location = new System.Drawing.Point(10, 11);
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
            this.log.Location = new System.Drawing.Point(65, 29);
            this.log.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.log.Name = "log";
            this.log.Size = new System.Drawing.Size(35, 13);
            this.log.TabIndex = 2;
            this.log.Text = "label1";
            // 
            // groupBox1
            // 
            this.groupBox1.Controls.Add(this.button1);
            this.groupBox1.Controls.Add(this.log);
            this.groupBox1.Dock = System.Windows.Forms.DockStyle.Right;
            this.groupBox1.Location = new System.Drawing.Point(1248, 0);
            this.groupBox1.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.groupBox1.Name = "groupBox1";
            this.groupBox1.Padding = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.groupBox1.Size = new System.Drawing.Size(188, 849);
            this.groupBox1.TabIndex = 3;
            this.groupBox1.TabStop = false;
            this.groupBox1.Text = "groupBox1";
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(5, 63);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(75, 23);
            this.button1.TabIndex = 3;
            this.button1.Text = "Show Street";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
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
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.PictureBox pictureZone;
        private System.Windows.Forms.Label log;
        private System.Windows.Forms.GroupBox groupBox1;
        private System.Windows.Forms.Button button1;


    }
}

