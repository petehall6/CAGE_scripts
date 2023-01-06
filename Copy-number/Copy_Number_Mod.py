
import csv
import openpyxl as pyx
import sys
maxInt = sys.maxsize

wb = pyx.Workbook()
ws = wb.active
while True:
    
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)



with open('all_lines_test.csv') as file:
    reader = csv.reader(file, delimiter=",")
    for row in reader:
        ws.append(row)

wb.save('test.xlsx')


