namespace emailer
{
    partial class Form1
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
            this.label1 = new System.Windows.Forms.Label();
            this.openFileDialog1 = new System.Windows.Forms.OpenFileDialog();
            this.button1 = new System.Windows.Forms.Button();
            this.label2 = new System.Windows.Forms.Label();
            this.export_label = new System.Windows.Forms.Label();
            this.load_btn = new System.Windows.Forms.Button();
            this.dataGridView1 = new System.Windows.Forms.DataGridView();
            this.create_email_btn = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).BeginInit();
            this.SuspendLayout();
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(23, 27);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(104, 13);
            this.label1.TabIndex = 0;
            this.label1.Text = "New Emailing Test 1";
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(45, 140);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(183, 23);
            this.button1.TabIndex = 3;
            this.button1.Text = "Select Excel Export";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(45, 78);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(111, 13);
            this.label2.TabIndex = 5;
            this.label2.Text = "Selected Excel Export";
            // 
            // export_label
            // 
            this.export_label.AutoSize = true;
            this.export_label.Location = new System.Drawing.Point(48, 111);
            this.export_label.Name = "export_label";
            this.export_label.Size = new System.Drawing.Size(0, 13);
            this.export_label.TabIndex = 6;
            // 
            // load_btn
            // 
            this.load_btn.Location = new System.Drawing.Point(48, 184);
            this.load_btn.Name = "load_btn";
            this.load_btn.Size = new System.Drawing.Size(180, 23);
            this.load_btn.TabIndex = 7;
            this.load_btn.Text = "Load Template";
            this.load_btn.UseVisualStyleBackColor = true;
            this.load_btn.Click += new System.EventHandler(this.load_btn_Click);
            // 
            // dataGridView1
            // 
            this.dataGridView1.ColumnHeadersHeightSizeMode = System.Windows.Forms.DataGridViewColumnHeadersHeightSizeMode.AutoSize;
            this.dataGridView1.Location = new System.Drawing.Point(45, 328);
            this.dataGridView1.Name = "dataGridView1";
            this.dataGridView1.Size = new System.Drawing.Size(670, 150);
            this.dataGridView1.TabIndex = 8;
            // 
            // create_email_btn
            // 
            this.create_email_btn.Location = new System.Drawing.Point(45, 237);
            this.create_email_btn.Name = "create_email_btn";
            this.create_email_btn.Size = new System.Drawing.Size(183, 23);
            this.create_email_btn.TabIndex = 9;
            this.create_email_btn.Text = "Create Emails";
            this.create_email_btn.UseVisualStyleBackColor = true;
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(982, 541);
            this.Controls.Add(this.create_email_btn);
            this.Controls.Add(this.dataGridView1);
            this.Controls.Add(this.load_btn);
            this.Controls.Add(this.export_label);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.label1);
            this.Name = "Form1";
            this.Text = "Form1";
            ((System.ComponentModel.ISupportInitialize)(this.dataGridView1)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.OpenFileDialog openFileDialog1;
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label export_label;
        private System.Windows.Forms.Button load_btn;
        private System.Windows.Forms.DataGridView dataGridView1;
        private System.Windows.Forms.Button create_email_btn;
    }
}

