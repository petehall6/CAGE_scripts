from urllib import request
from urllib.request import urlopen
from urllib import response
from urllib3.exceptions import InsecureRequestWarning
from urllib3 import disable_warnings
from urllib.parse import urlencode
import os
#TODO see if urllib 3 can do all the urllib calls
cello_dom = "https://www.cellosaurus.org/str-search/api/"

disable_warnings(InsecureRequestWarning)
#response = request.Reget(cello_dom, verify=False)
status_code = urlopen(cello_dom)
print(status_code.getcode())

###Get Cellosaurus accession code

#urlparm = '?' + urlencode

#req = request.Request()

""" cell_line = "CVCL_K562"


req = requests.get (f"https://api.cellosaurus.org/cell-line/{cell_line}/?fields=ac,str&format=txt",verify=False)

print(req.text)

with open('cellosaurus_str.txt', 'w') as f:
    f.write(req.text) """