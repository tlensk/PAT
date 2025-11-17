"""
@author: Tatiana Lenskaia
All rights reserved, 2025
E-mail: lensk010@gmail.com
"""




import core_methods as cm
import GB_lib as gbl
import PATgb_lib as pgbl


import os
import sys


# E-value cut-off

e_thr = 0.001



# Add path to PAT_categories_current.txt if not in the same directory as the main script
path = ""
fCatName = "PAT_categories_current.txt"
res_cat = cm.GetDictFromFile(path+fCatName, sep = "\t", header = 1, unique_col=3)
d_cat = res_cat[1]


# Add path to a phage genome if not in the same directory as the main script
path = ""
label = "PhG_"





fLog1 = open(path+"CDS_log.txt","w")
print("order","g_id","#rec","#unusual_start_codon","#CDS","#pseudogenes", file = fLog1, sep="\t")
cc = 0



    

# Step 1: CDS extraction



fInName = label+".gb"
print(fInName)


cc += 1
t = pgbl.SeqDB_joinCDS_update_stat_old_GI(fInName,path)

print(cc, fInName, "\t".join([str(x) for x in t]), file = fLog1, sep = "\t")
print(cc,fInName,t)
fLog1.close()

print("==============")
stage_res = input("Step1. CDS extraction is completed.")
print("==============")

#===============================================

# Step 2: CDS annotation


fInName_cds = label.rsplit(".",1)[0]+"_cds.faa"
print(fInName_cds)

os.system("sh test.sh "+path+" "+fInName_cds+" " + label)

fInName = label+".out"

print(fInName)

print("==============")
stage_res = input("Step2. CDS annotation is completed.")
print("==============")
#=====================

# Step 3: CDS annotation analysis

c = 0
l = 22



d_hits = {}
fIn = open(path+fInName, "r")
for line in fIn:
    if line[0] != "#":
        c += 1
        t_line = []
        ln = line.strip()
        for i in range(l):
            t = ln.split(" ",1)
            t_line = t_line + [t[0]]
            ln = t[1].strip()
        t_line = t_line + [ln]
        #print(t_line)
        
        name = t_line[0]+" "+t_line[-1]
        hmm_nm = t_line[3]
        hmm_id = t_line[4]
        e = float(t_line[6])
        pos_st = int(t_line[-6])
        pos_fn = int(t_line[-5])
        aln_ln = pos_fn-pos_st+1
        
        
        if e < e_thr:          
            if name not in d_hits:
                d_hits[name] = []
            d_hits[name] += [[e, hmm_nm, hmm_id, pos_st, pos_fn,aln_ln]]
print("d_hits",name, len(d_hits[name]))
        
        

fIn.close()
print("Number of entries in hmm.out:",c)

print("Number of proteins with at least one hmm hit:",len(d_hits))

d = {}
for it in d_hits:
    n = len(d_hits[it])
    if n not in d:
        d[n] = []
    d[n] += [it]
    


# Analysis of CDS in genomes

res_fa = gbl.GetFasta(path+fInName_cds)

d_fa = res_fa[0]
t_ids = res_fa[-1]
print("Number of proteins:", len(d_fa))

d_gens = {}
for it in d_fa:
    g_id = it[1:].split("|")[0]
    # print(g_id)
    if g_id not in d_gens:
        d_gens[g_id] = []
    d_gens[g_id] += [it]


fOutName = label+"_PAThits.txt"
fOutName1 = label+"_PATlabels.txt"
sp = "\t"
fOut = open(path+fOutName, "w")
fOut1 = open(path+fOutName1, "w")

print("g_id", "pid","e", "HMM_nm", "HMM_id", "aln_st","aln_fn", "aln_len", "All_hits",sep = sp, file = fOut)
print("g_id", "g_nm","PATlabels", sep = sp, file = fOut1)


for g_id in d_gens:
    #g_id = "NC_019516.2"
    cg = d_gens[g_id]
    #print(cg)
    st = ""
    for item in cg:
        #print(item)
        item = item.strip()[1:]
        lb = "*"
        if item in d_hits:
            #print(item, d_hits[item])
            tt = d_hits[item]
            # if e values are the same pick the longest alignmed part
            tt = sorted(tt, key=lambda x: x[0])
            print(g_id + sp+ item +sp + sp.join([str(x) for x in tt[0]]), tt, sep = sp, file = fOut)
            #fOut.write(g_id + sp+ item +sp + sp.join([str(x) for x in tt[0]])+sp+",".join(tt)+"\n")
            lb = tt[0][1].split("-")[0]
        st = st+"'"+lb+"'"+","
    st = "["+st[:-1]+"]"
    #print(g_id, d_phs[g_id][1], d_phs[g_id][2], round(float(d_phs[g_id][-2])*100, 2),  len(d_gens[g_id]), st, sep = sp, file=fOut1)
    print(g_id, label, st, sep = sp, file=fOut1)
fOut.close()
fOut1.close()


print("==============")
stage_res = input("Step3. CDS annotation analysis is completed.")
print("==============")

#=====================

# Step 4. GenBank update

# Update .gb with the hits



res2 = cm.GetDictFromFile(path+fOutName, sep = "\t", header=1, unique_col=1)

d_pat = res2[1]
t_pat = res2[0]

t_ind = []

d_conv = {}
d_loc = {}

for pid in t_pat:
    tt = pid.split("|")
    
    g_id = tt[0].split(".")[0]
    pr_id = tt[1]
    pr_annot = tt[2]
    pr_pos = tt[3]
    pr_len = tt[4]
    
    
    if g_id not in d_conv:
        d_conv[g_id] = {}
    if pr_id not in d_conv[g_id]:
        d_conv[g_id][pr_id] = pid
    
    
    st_loc = g_id+"|"+pr_pos
   
    if st_loc not in d_loc:
        d_loc[st_loc] = pid
    else:
        print("d_loc duplicate")
     
    




t_rec = gbl.GetGBRecords(path+label+".gb")



t_cds = []


fOut = open(path+label+"_updated.gb","w")

for rec in t_rec:
    n_rec = len(rec)
    g_id_cur = rec[0].strip().split()[1]
    #print(g_id_cur)
    
    if g_id_cur in d_conv:
        d_cds = d_conv[g_id_cur]
        #print(g_id_cur, d_cds)
    
    annot = "###"
    pr_id = "###"
    
    for ii in range(n_rec):
        ln = rec[ii]
        
        if ln.strip()[0:3] == "CDS":
            t_cds.append(ii)
            
            st_loc =g_id_cur+ "|"+ln.strip().split()[-1]
            
            if st_loc in d_loc:
                #print(d_loc[st_loc], "!!!")
                
                annot = d_loc[st_loc].split("|")[2]
                
                hmm_id = d_pat[d_loc[st_loc]][3]
                
                pr_id = d_loc[st_loc].split("|")[1]
                
                
                new_annot = d_cat[hmm_id][2]+" ### "+d_cat[hmm_id][3]
                
                #print(annot)
                #print(new_annot)
            else:
                annot = "###"
                pr_id = "###"
            
            fOut.write(ln)
            
        elif ln.strip()[0:9] == "/product=":
            if annot in ln:
                fOut.write(ln.split("=")[0]+'="'+new_annot+'"'+"\n")
            else:
                if annot != "###":
                    print("multi-line annotation!!!")
                    break
                else:
                    fOut.write(ln)
        elif ln.strip()[0:12] == "/protein_id=":
            if (pr_id != "###") and (pr_id not in ln):
                print("check updated protein id:", pr_id)
                break
            else:
                fOut.write(ln)
        else:
            fOut.write(ln)

fOut.close()

#print("==============")
stage_res = input("Step 4. GenBank update is completed.")
print("==============")

#===============================

# Step 5. GenBank update quality control.

# Check re-aanotation

fInName1 = label+ ".gb"
fInName2 = label+"_updated.gb"


fIn1 = open(path+fInName1, "r")
lines1 = fIn1.readlines()
fIn1.close()


fIn2 = open(path+fInName2, "r")
lines2 = fIn2.readlines()
fIn2.close()



print(len(lines1), len(lines2))

n = len(lines1)

c = 0
k = 0
for i in range(n):
    if lines1[i] == lines2[i]:
        c += 1
    else:
        k += 1
        print(k)
        print(lines1[i].strip())
        print(lines2[i].strip())
        print()
print(c, n-c, k)
    


print("==============")
stage_res = input("Step 5. GenBank update quality check is completed.")
print("==============")

    


fLog1.close()
