"""
Copyright: (c) 2025 Tatiana Lenskaia. All rights reserved.
E-mail: t.lenskaia@icloud.com
"""

# Retrieving sequence from NCBI in GenBank format


from Bio import Entrez
from Bio import SeqIO

# Specify email

Entrez.email = ""


def GetListFromFile(fInName):
    t = []
    fIn = open(fInName, "r")
    for line in fIn:
        line = line.strip()
        if line != "":
            if line not in t:
                t.append(line)
    #print(len(t))
    fIn.close()
    return t


def FetchFullGB(acc):
    # Search for the GenBank record by accession number
    handle = Entrez.efetch(db="nucleotide", id = acc, rettype="gbwithparts", retmode="text")
    record = handle.read()
    handle.close()
    return record



path = ""
fListName = "PAT_ids.txt"


t_list = cm.GetListFromFile(path+fListName)

print(len(t_list))

fLog = open(path+"fLog.txt","w")
fLog.close()

c = 0
for g_id in t_list:
    c += 1
    rec = FetchFullGB(g_id)
    #print(rec)
    print(c, g_id)
    
    fLog = open(path+"fLog.txt","a")
    fLog.write(g_id+"\n")
    fLog.close()
    
    fOut = open(path+g_id+".gb","w")
    fOut.write(rec)
    fOut.close()
    







