##convert csv from wide to long.  Sqlite3 column limit = 2000.
import pandas as pd
import sqlite3 as sql
import csv

data_store = "copyNumSql.csv"
long_csv = "copy_number_long.csv"
copy_db = "copy_num.db"

#create sql db from extra wide csv by melting wide into long csv and then adding into generated sql db
def create_db(data_store):
    #open large csv via chunking
    print("Chunking CSV.  Could take a while\n")
    chunk = pd.read_csv(data_store, chunksize=1000)
    df = pd.concat(chunk)
    df = df.rename(columns={"cell_line_display_name": "cell_name"})
    #will convert to int
    df = df.round()
    print("Chunk complete\n")
    #print(len(columnNames))
    #melt table and convert from wide to long
    #columnNames = list(df.columns)
    melt = pd.melt(df, id_vars =['cell_name']).reset_index()
    #TODO is needed?
    melt = melt.drop(['index'], axis=1)
    melt = melt.rename(columns={"variable" : "gene", "value" : "copy_number" })
    #TODO check to see if this is necessary for primary key generation
    melt.index = melt.index+1

    print("\n\n")
    print(melt.info(verbose=False))
    print("\n\n")
    #melt.drop(index=melt.index[0], axis=0, inplace=True)
    print(melt.head(5))
    
    input("Press enter to melt: ")

    print("\ncreating long form csv")
    melt.to_csv("copy_number_long.csv")
    print("long form csv created.  Job's done!\n")

    #create SQLite3 database and upload copy number data

def create_table(copy_db):
    print("Connecting to db")
    conn = sql.connect(copy_db)
    cur = conn.cursor()

    create_fresh_table = '''DROP TABLE IF EXISTS copy_num'''
    cur.execute(create_fresh_table)
    print("Creating table")
    create_copy_num_table = '''CREATE TABLE copy_num(cell_id INTEGER AUTO_INCREMENT, cell_name TEXT NOT NULL, gene TEXT NOT NULL, num_of_copies INTEGER NOT NULL);'''
    cur.execute(create_copy_num_table)

    print("Table created")
    conn.commit()
    conn.close()

#upload csv into table
def uploadToSql(long_csv, copy_db):
    print("Reading csv\n")
    csvFile = open(long_csv)

    long_data = csv.reader(csvFile)
    print("Inserting data\n")
    insert_data = '''INSERT INTO copy_num (cell_id, cell_name, gene, num_of_copies) VALUES (?, ?, ?, ?)'''

    print("Connecting to db\n")
    conn = sql.connect(copy_db)
    cur = conn.cursor()
    print("connection established\n")
    cur.executemany(insert_data, long_data)
    print("Upload Complete\n")
    #table_info = '''PRAGMA table_info copy_num'''
    #cur.execute(table_info)
    
    conn.commit()
    conn.close()


def test_db(copy_db):
    print("Connecting to database. \n")
    conn = sql.connect(copy_db)
    cur = conn.cursor()
    print("Connection made \n")

    
    #EXPANDED TEST
    cell_choice = str(input("Enter cell line: ")).upper()
    gene_choice = str(input("Enter gene: ")).upper()
    #test_select = f''' SELECT num_of_copies FROM copy_num WHERE cell_name LIKE '{cell_choice}' AND gene LIKE '{gene_choice}' '''
    test_select = f''' SELECT cell_name, gene, num_of_copies FROM copy_num WHERE cell_name LIKE '%{cell_choice}%' AND gene LIKE '%{gene_choice}%' LIMIT 5'''
    
    cur.execute(test_select)
    query = cur.fetchall() 
    query_hits = str(len(query))
    print("Number of hits: "+query_hits)
    
    #TODO fix this
    for row in query:
        print(row)
    
    conn.commit()
    conn.close()

#add small interface so I dont keep creating and inserting into same table
while True:
    go = input("Choose task: \n 1) Create DB: \n 2) Create Table\n 3) Upload to DB\n 4) Test DB \n or quit\n\n").lower()
    
    if go == '1':
        create_db(data_store)
    if go == '2':
        create_table(copy_db)
    if go == '3':
        uploadToSql(long_csv, copy_db)
    if go == '4':
        test_db(copy_db)
    if go == 'quit' or go =='q':
        print("Goodbye")
        exit()











