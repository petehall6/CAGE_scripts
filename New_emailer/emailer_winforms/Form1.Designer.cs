namespace emailer_winforms
{
    partial class Form1
    {
        /// <summary>
        ///  Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        ///  Clean up any resources being used.
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
        ///  Required method for Designer support - do not modify
        ///  the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            FilePicker = new Button();
            excel_titlebx = new TextBox();
            load_btn = new Button();
            csv_dt_txtbx = new TextBox();
            SuspendLayout();
            // 
            // FilePicker
            // 
            FilePicker.Location = new Point(154, 175);
            FilePicker.Margin = new Padding(3, 4, 3, 4);
            FilePicker.Name = "FilePicker";
            FilePicker.Size = new Size(218, 31);
            FilePicker.TabIndex = 0;
            FilePicker.Text = "Choose Excel Template";
            FilePicker.UseVisualStyleBackColor = true;
            FilePicker.Click += button1_Click;
            // 
            // excel_titlebx
            // 
            excel_titlebx.AllowDrop = true;
            excel_titlebx.Location = new Point(398, 175);
            excel_titlebx.Margin = new Padding(3, 4, 3, 4);
            excel_titlebx.Name = "excel_titlebx";
            excel_titlebx.Size = new Size(217, 27);
            excel_titlebx.TabIndex = 1;
            // 
            // load_btn
            // 
            load_btn.Location = new Point(402, 223);
            load_btn.Margin = new Padding(3, 4, 3, 4);
            load_btn.Name = "load_btn";
            load_btn.Size = new Size(213, 31);
            load_btn.TabIndex = 2;
            load_btn.Text = "Load Template";
            load_btn.UseVisualStyleBackColor = true;
            load_btn.Click += load_btn_Click_1;
            // 
            // csv_dt_txtbx
            // 
            csv_dt_txtbx.Location = new Point(719, 180);
            csv_dt_txtbx.Multiline = true;
            csv_dt_txtbx.Name = "csv_dt_txtbx";
            csv_dt_txtbx.Size = new Size(567, 422);
            csv_dt_txtbx.TabIndex = 3;
            // 
            // Form1
            // 
            AutoScaleDimensions = new SizeF(8F, 20F);
            AutoScaleMode = AutoScaleMode.Font;
            ClientSize = new Size(1327, 664);
            Controls.Add(csv_dt_txtbx);
            Controls.Add(load_btn);
            Controls.Add(excel_titlebx);
            Controls.Add(FilePicker);
            Margin = new Padding(3, 4, 3, 4);
            Name = "Form1";
            Text = "Form1";
            ResumeLayout(false);
            PerformLayout();
        }

        #endregion
        private Button FilePicker;
        private TextBox excel_titlebx;
        private Button load_btn;
        private TextBox csv_dt_txtbx;
    }
}