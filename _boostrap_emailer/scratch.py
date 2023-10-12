from html.parser import HTMLParser
import os
from bs4 import BeautifulSoup as bs
import re

sig_htm = "Email signature (phall3@stjude.org).htm"
with open(os.path.join(os.getcwd(), sig_htm), "r") as f:
    html = f.read()
    
    hlen = len(html)
    
    #print(html)

#print(html)


body_pattern = "(?:<body)(.|\n)*?<\/html>"

#returns match object.  matched group is accessed by .group() because life
text = re.search(body_pattern, html)

sig = text.group()

print(text.group())
