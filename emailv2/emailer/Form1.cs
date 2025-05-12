using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.IO;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using Excel = Microsoft.Office.Interop.Excel;
using Outlook = Microsoft.Office.Interop.Outlook;
using ExcelDataReader;

namespace emailer
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        DataSet ds;
        private void button1_Click(object sender, EventArgs e)
        {


            OpenFileDialog ofd = new OpenFileDialog();
            ofd.Title = "Select SRM Template";
            ofd.InitialDirectory = "C:\\Users\\USER_NAME\\Downloads".Replace("USER_NAME", Environment.UserName);
            ofd.Filter = "Excel Files|*.xls;*.xlsx;*.xlsm";

            ofd.ShowDialog();

            if (ofd.FileName != "")
            {
                export_label.Text = ofd.FileName;
            }



        }

        private void load_btn_Click(object sender, EventArgs e)
        {
            //Need columns
            //SRM Order #, PI, Requested By, Project Scope, Species, Project Objective, Target Gene Name
            //Project Number

            FileStream fileStream = File.Open(export_label.Text, FileMode.Open, FileAccess.Read);
            IExcelDataReader reader = ExcelReaderFactory.CreateBinaryReader(fileStream);
            ds = reader.AsDataSet(new ExcelDataSetConfiguration()
            {
                ConfigureDataTable = (_) => new ExcelDataTableConfiguration()
                {
                    UseHeaderRow = true
                }
            });


            dataGridView1.DataSource = ds.Tables[0];



        }
    }
}
