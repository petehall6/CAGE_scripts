using System.Data;
using System.Windows.Forms.VisualStyles;

namespace emailer_winforms
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }


        private void button1_Click(object sender, EventArgs e)
        {

            int size = -1;
            OpenFileDialog openExcelDialog = new OpenFileDialog();
            openExcelDialog.InitialDirectory = Environment.GetFolderPath(Environment.SpecialFolder.Desktop);
            openExcelDialog.Filter = "(*.csv)|*.csv";
            DialogResult result = openExcelDialog.ShowDialog();


            if (result == DialogResult.OK)
            {


                string csv_name = openExcelDialog.FileName;
                try
                {
                    excel_titlebx.Text = csv_name;
                }
                catch (IOException)
                {
                }

            }
            Console.WriteLine(size);
            Console.WriteLine(result);
        }

        private void load_btn_Click_1(object sender, EventArgs e)
        {
            string csv_file_path = excel_titlebx.Text;
            char csvDelimiter = '\n';

            DataTable csv_dt = new DataTable();
            using (StreamReader stream_reader = new StreamReader(csv_file_path))
            {
                string[] headers = stream_reader.ReadLine().Split(csvDelimiter);
                foreach (string header in headers)
                {
                    try
                    {

                        csv_dt.Columns.Add(header);

                    }
                    catch { }
                }
                while (!stream_reader.EndOfStream)
                {

                    string[] rows = stream_reader.ReadLine().Split(csvDelimiter);
                    DataRow datarow = csv_dt.NewRow();
                    for (int i = 0; i < headers.Length; i++)
                    {
                        datarow[i] = rows[i];

                    }
                    csv_dt.Rows.Add(datarow);
                }

            }

            csv_dt_txtbx.Text = csv_dt.ToString();
        }

    }
}