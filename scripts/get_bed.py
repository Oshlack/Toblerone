import requests
from bs4 import BeautifulSoup


# add ref to origin
# add option for switching between refernces automatically

url = 'http://genome-euro.ucsc.edu/cgi-bin/hgTables?'    
#url = 'http://genome-euro.ucsc.edu/cgi-bin/hgTables?hgsid=201790284_dkwVYFu7V6ISmTzFGlXzo23aUhXk'    
session = requests.Session()

params = {
#    'hgsid': '201790284_dkwVYFu7V6ISmTzFGlXzo23aUhXk',
    'jsh_pageVertPos': '0',
    'clade': 'mammal',
    'org': 'Human',
    'db': 'hg38',
#    'org': 'Mouse',
#    'db': 'mm10',
    'hgta_group': 'genes',
    'hgta_track': 'knownGene',
    'hgta_table': 'knownGene',
    'hgta_regionType': 'genome',
    #'position': 'chr9:21802635-21865969',
    'hgta_outputType': 'bed',
    'boolshad.sendToGalaxy': '0',
    'boolshad.sendToGreat': '0',
    'boolshad.sendToGenomeSpace': '0',
    'hgta_outFileName': '',
    'hgta_compressType': 'none',
    'hgta_doTopSubmit': 'get output'
}


response = session.post(url, data=params)
#print(response.content)


page = BeautifulSoup(response.text)
#print(page)
hidden_tags = page.find_all("input", type="HIDDEN")
hgsid="EMPTY"
for tag in hidden_tags:
  #print(tag.get("name"))
  if tag['name'] == "hgsid":
     hgsid = tag['value']
#print(hgsid)
paramsbed = {
    'hgsid': hgsid, 
    'fbQual': 'whole',
    'hgta_doGetBed': 'get BED'
}


response2 = session.get(url, data=paramsbed)
print(response2.text)
#print(response2)
