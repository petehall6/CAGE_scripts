import csv
from io import StringIO

#creating
db = 'cellosaurus_offline_db_short.txt'

with open(db) as data_set, open('output.tsv', 'w', newline="") as output:
    csv_output = csv.writer(output, delimiter='\t')
    csv_output.writerow(['Line','ACCN','Amelogenin','CSF1PO','D13S317',
                         'D162539','D18S51','D21S11','D3S1358',
                         'D5S818','D7S820','D8S1179',
                         'FGA','Penta D','Penta E',
                         'TH01','TPOX','vWA'])
    
    for line in data_set:
        if line.startswith('ID'):
            type = 'Line'
            id = line.split('  ')
            id=str(id[1]).replace('\n','').replace('CVCL_',"")
        if line.startswith('AC'):
            type = 'ACCN'
            id = line.split('  ')
        csv_output.writerow([id,type])