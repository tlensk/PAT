#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Tatiana Lenskaia. E-mail: t.lenskaia@utoronto.ca

"""

# Calculate GC-content
def CalcGC(seq):
    gc = ""

    n = len(seq)
    if n < 0:
        return gc
    
    seq = seq.lower()
    
    c_a = seq.count("a")
    c_t = seq.count("t")
    c_c = seq.count("c")
    c_g = seq.count("g")
    
    n_atcg = c_a + c_t + c_c + c_g
    
    gc = round((c_c+c_g)/n_atcg, 4)
    oth = n - n_atcg
    
    return [gc, oth]




# Check the data quality of a GenBank file

def SeqGB_sizecheck_singleGB(fInName):

    
    t_rec= []
    
    fIn = open(fInName, "r")
    lines = fIn.readlines()
    fIn.close()
    
    
    
    for line in lines:
        #print(line)
        if line != "\n":
            t_rec.append(line)
        if line == "//\n":
            break

    
    item = t_rec
    

    t_lc = item[0].strip().split()
    #print(t_lc)
    g_acc = t_lc[1]
    g_sz = int(t_lc[2])
    g_base = t_lc[4]
    g_type = t_lc[5]
    
    
    
    j = -1
    j_v = -1
    j_s = -1
    j_d = -1
    
    
    for ln in item:
        if ln[0:7] == "VERSION":
            j_v = item.index(ln)
        
        if ln[0:6] == "SOURCE":
            j_s = item.index(ln)
            
        if ln[0:10] == "DEFINITION":
            j_d = item.index(ln)
        
        if ln[0:6] == "ORIGIN":
            j = item.index(ln)
            #print(item[j-1:j+2])
            break
    t_seq = item[j+1:-1]
    #print(t_seq[0]+t_seq[-1])
    
    g_ver = item[j_v].strip().split()[-1]
    g_source = " ".join(item[j_s].strip().split()[1:])
    g_def = " ".join(item[j_d].strip().split()[1:])
    

    if g_ver == "VERSION":
        g_id = g_acc
    else:
        g_id = g_ver
    
    


    seq = ""
    for l in t_seq:
        t_l = l.strip().split()
        st = "".join(t_l[1:])
        #print(l)
        #print(st)
        seq = seq+st
    #print(seq[0:10])
    n_seq = len(seq)
    if n_seq != g_sz:
        fl = 1
    else:
        fl = 0
    res = CalcGC(seq)
    #t_stat = [g_id, g_source, g_sz,g_base, g_type,n_seq,fl,res[0], res[1]]
    t_stat = [g_id, g_def, g_sz,g_base, g_type,n_seq,fl,res[0], res[1]]
    

    return t_stat











path = ""
fListName = "PAT_ids.txt"

t_list = cm.GetListFromFile(path+fListName)

print(len(t_list))


sep = "\t"

fOutName = path + fListName.split(".")[0]+"_sizecheck.txt"
fOut = open(fOutName, "w") 
print("id","version","definition","rec_len","gbase","gtype","download_len","len_mismatch", "gc","Not_ACGT", sep = sep, file = fOut)
fOut.close()


c = 0

for g_id in t_list:
    
    
    fInName = path+g_id+".gb"
    
    c += 1
    
    t = SeqGB_sizecheck_singleGB(fInName)
    
    fOut = open(fOutName, "a")
    
    fOut.write(g_id+sep+sep.join([str(x) for x in t])+"\n")
    #print(g_id, g_sour, g_sz,g_base, g_type,n_seq,fl,res[0], res[1], sep = sep, file = fOut)
    fOut.close()
    
    print(c, g_id)
    