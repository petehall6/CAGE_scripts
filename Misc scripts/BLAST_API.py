import os
import pandas as pd
import requests
import lxml
from bs4 import BeautifulSoup




blast_url ="https://blast.ncbi.nlm.nih.gov/Blast.cgi"

source = requests.get(blast_url)

soup = BeautifulSoup(source, 'lxml')

