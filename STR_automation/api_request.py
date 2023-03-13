import requests 
from urllib3.exceptions import InsecureRequestWarning
from urllib3 import disable_warnings
import os

cello_dom = "https://www.cellosaurus.org/str-search/api/"

disable_warnings(InsecureRequestWarning)
response = requests.get(cello_dom, verify=False)
print(response.status_code)

###Get Cellosaurus accession code





""" cell_line = "CVCL_K562"


req = requests.get (f"https://api.cellosaurus.org/cell-line/{cell_line}/?fields=ac,str&format=txt",verify=False)

print(req.text)

with open('cellosaurus_str.txt', 'w') as f:
    f.write(req.text) """