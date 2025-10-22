#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 20:15:54 2024

@author: tatianalenskaia
"""

import GB_lib as gbl

path = "/Users/tatianalenskaia/_Current_projects/TT_clustering/compare/"



fInName = "TT_prs.fasta"

d_rec = gbl.GetFata(path+fInName)[0]

print(len(d_rec))



def MaskAlign(seq, st, fn):
    n = len(seq)
    if (st < 0) or (fn >= n):
        return ""
    s = ""
    for i in range(n):
        if (i < st) or (i > fn):
            s = s+"."
        else:
            s = s+seq[i]
    return s



fInName = "TT_Clusters.txt"



fIn = open(path+fInName, "r")

lines = fIn.readlines()
lines = lines[1:]

d = {}
sep = "\t"
for line in lines:
    
    line = line.strip()
    t_line = line.split(sep)
    
    cl = int(t_line[0])
    cl_size = int(t_line[1])
    
    cl_desc = t_line[2]
    
    t_cl = cl_desc[1:-1].split(",")
    #print(len(t_cl))
    
    #print(t_line[0:2])
    
    #print(cl, cl_size, t_cl)
    
    if len(t_cl) != cl_size:
        print("size does not match!")
        break
    
    t_prs = []
    '''
    for it in t_cl:
        print(it)
    '''
    
    #fOutName = str(cl_size)+"_Cluster"+str(cl)+".fasta"
    
    #fOut = open(path+fOutName, "w")
    
    
    for it in t_cl:
        it = it.strip()
        #print("1",it)
        it = it.split(":")[0][1:-1]
        #print("2",it)
        
        t_it = it.split(" ", 1)
        p_id = t_it[0]
        it = t_it[1]
        t_it = it.split("[")
        phs = t_it[1][:-1]
        pr = t_it[0].strip()
        #print(p_id, pr, phs, sep=",")
        
        sq = ""
        if p_id in d_rec:
            sq = d_rec[p_id]
        else:
            print("no sequence!")
        
        if p_id not in d:
            d[p_id] = []
        d[p_id] = d[p_id] + [cl,pr,phs, sq]
        
        
        #fOut.write(">"+p_id+" "+pr+" ["+ phs+"]"+"\n"+sq+"\n")
    
#fOut.close()
    
fIn.close()

print(len(d))
#print(d)




path = "/Users/tatianalenskaia/HMM/All_HMMs/HMM_db/"
fInName = "A_TT_parse.csv"



fIn = open(path+fInName, "r")

lines = fIn.readlines()
lines = lines[1:]
fIn.close()


d_hmm = {}

sep = ","
for line in lines:
    
    line = line.strip()
    t_line = line.split(sep)
    print(t_line)
    
    t_id = t_line[0].split("|")

    g_id = t_id[0]
    p_id = t_id[1]
    
    pr_len = int(t_line[1])
    
    hmm_nm = t_line[2]
    hmm_ac = t_line[3]
    hmm_len = int(t_line[4])

    hmm_id = hmm_nm
    if hmm_ac != "-":
        hmm_id = hmm_ac

    e = float(t_line[5])
    
    pos_st = int(t_line[6])
    pos_fn = int(t_line[7])
    
    alg = ""
    if p_id in d_rec:
        print(p_id, d_rec[p_id][-1])
        seq = d_rec[p_id][-1]
        alg = MaskAlign(seq, pos_st-1, pos_fn-1)
        
    
    
    
    if hmm_id not in d_hmm:
        d_hmm[hmm_id] = {}
    dd = d_hmm[hmm_id]
    
    if p_id not in dd:
        dd[p_id] = [e, pos_st, pos_fn,alg]
    d_hmm[hmm_id] = dd




#print(d_hmm)
    

for hmm_id in d_hmm:
    print(hmm_id, len(d_hmm[hmm_id]))


hmm_id = "PF06841.15"


dd = d_hmm[hmm_id]



fOut = open(path+"AllDetails_"+hmm_id+".csv","w")
for pr in dd:
    #print(pr)
    
    tt = dd[pr]
    
    print(tt)
    
    
    if pr in d:
        #print(pr, d[pr][0])
        
        
        
        print(pr,d[pr][0], tt[0], len(d[pr][-1]), tt[1], tt[2], tt[3], sep = ",", file=fOut)
        
fOut.close()



"""
# Previous parse HMM output, only top hits with the smallest E-value

fInName = "TT_PAT_hits.txt"



fIn = open(path+fInName, "r")

lines = fIn.readlines()
lines = lines[1:]
fIn.close()


d_hmm = {}

sep = "\t"
for line in lines:
    
    line = line.strip()
    t_line = line.split(sep)
    #print(len(t_line))
    
    g_id = t_line[1]
    p_id = t_line[2]
    
    n_hits = int(t_line[3])
    hits = t_line[-1]
    
    
    #print(hits)
    
    t_hits = hits[3:-2].split("[")

    
    
    n = len(t_hits)
    
    if n != n_hits:
        print("mismatch!", n, n_hits)



    hit = t_hits[0]

    print(hit)    
    
    t_hit = hit.split(",")
    
    hmm_id = t_hit[2].strip()[1:-1]
    #hmm_ac = t_hit[3].strip()[1:-1]
    #hmm_nm = t_hit[4].strip()[1:-1]
    
    e = float(t_hit[5])
    
    pos_st = int(t_hit[6])
    pos_fn = int(t_hit[7])
    alg = t_hit[8][:-1]
    
    print(alg)
    
    #print(hmm_id, hmm_ac, hmm_nm, e, pos_st, pos_fn)
    
    
    
    if hmm_id not in d_hmm:
        d_hmm[hmm_id] = {}
    dd = d_hmm[hmm_id]
    
    if p_id not in dd:
        dd[p_id] = [e, pos_st, pos_fn,alg]
    d_hmm[hmm_id] = dd
    
    
    
#print(d_hmm)
    

for hmm_id in d_hmm:
    print(hmm_id, len(d_hmm[hmm_id]))


hmm_id = "PF06841.15"


dd = d_hmm[hmm_id]



fOut = open(path+"Details_"+hmm_id+".csv","w")
for pr in dd:
    #print(pr)
    
    tt = dd[pr]
    
    print(tt)
    
    
    if pr in d:
        #print(pr, d[pr][0])
        
        
        
        print(pr,d[pr][0], tt[0], len(d[pr][-1]), tt[1], tt[2], tt[3], sep = ",", file=fOut)
        
fOut.close()


"""





