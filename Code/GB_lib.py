"""
@author: Tatiana Lenskaia
All rights reserved, 2025
E-mail: lensk010@gmail.com
"""



import os
import core_methods as cm

import gzip




def GetLoc(cds_loc):

    tt = []
    
    if cds_loc[0:16] == "complement(join(":
        cds_loc = cds_loc[16:-2]
        #print(cds_loc)
        tt = cds_loc.split(",")
        
        ttt = []
        
        for loc in tt:
            ttt.append([loc,"-"])
            
    else:
                
        if cds_loc[0:4] == "join":
            tt = cds_loc[5:-1].split(",")
        else:
            tt = [cds_loc]
    
    
        ttt = []    
        
        for loc in tt:
            loc = loc.strip()
            if loc[0:10] == "complement":
                ttt.append([loc[11:-1],"-"])
            else:
                ttt.append([loc, "+"])
            
    t_loc = []   
            
    for loc in ttt:
        
        strand = ""
        pos_st = -1
        pos_fn = -1
        
        if (">" in loc[0]) or ("<" in loc[0]):
            
            t_pos = loc[0].split("..")
            pos_st = t_pos[0]
            pos_fn = t_pos[-1]
            
            if pos_st[0] == "<":
                strand = "<"+loc[-1]
                pos_st = int(pos_st[1:])
            else:
                loc[-1] = "("+loc[-1]
                pos_st = int(pos_st)
                
                
            if pos_fn[0] == ">":
                strand += ">"
                pos_fn = int(pos_fn[1:])
            else:
                strand += ")"
                pos_fn = int(pos_fn)
            
            
        else:
            t_pos = loc[0].split("..")
            pos_st = int(t_pos[0])
            pos_fn = int(t_pos[-1])
            strand = "("+loc[-1]+")"
        
        if (strand != "") and (pos_st > 0) and (pos_fn > 0):
            t_loc.append([pos_st, pos_fn, strand])
        else:
            return []
                
    
    return t_loc








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



def GetFasta(fInName):
    fIn = open(fInName,"r")
    d_rec = {}
    t_rec = []
    t_ids = []
    for line in fIn:
        #print(line)
        if line[0] == ">":
            if t_rec != []:
                
                
                header = t_rec[0]
                
                sq_id = header
                #sq_id = header.split()[0][1:]
                
                t_ids.append(sq_id)
                
                
                sq = ""
                for it in t_rec[1:]:
                    it = it.strip()
                    sq = sq+it
                    
    
                if sq_id not in d_rec:
                    d_rec[sq_id] = sq
                #print(sq_id,len(sq), sq)
                
                t_rec = []
               
        t_rec.append(line)
    
    
    if t_rec != []:
        
        
        header = t_rec[0]
        
        sq_id = header
        #sq_id = header.split()[0][1:]
        
        t_ids.append(sq_id)
        
        
        sq = ""
        for it in t_rec[1:]:
            it = it.strip()
            sq = sq+it
            
    
        if sq_id not in d_rec:
            d_rec[sq_id] = sq
        #print(sq_id,len(sq), sq)

    return [d_rec, t_ids]


def GetFasta_header(fInName):
    fIn = open(fInName,"r")
    d_rec = {}
    t_rec = []
    t_ids = []
    for line in fIn:
        #print(line)
        if line[0] == ">":
            if t_rec != []:
                
                
                header = t_rec[0]
                
                
                sq_id = header.split()[0][1:]
                
                t_ids.append(sq_id)
                
                
                sq = ""
                for it in t_rec[1:]:
                    it = it.strip()
                    sq = sq+it
                    
    
                if sq_id not in d_rec:
                    d_rec[sq_id] = [header, sq]
                #print(sq_id,len(sq), sq)
                
                t_rec = []
               
        t_rec.append(line)
    
    
    if t_rec != []:
        
        
        header = t_rec[0]
        
        
        sq_id = header.split()[0][1:]
        
        t_ids.append(sq_id)
        
        
        sq = ""
        for it in t_rec[1:]:
            it = it.strip()
            sq = sq+it
            
    
        if sq_id not in d_rec:
            d_rec[sq_id] = [header, sq]
        #print(sq_id,len(sq), sq)

    return [d_rec, t_ids]





def GetFasta_lines(fInName):
    fIn = open(fInName,"r")
    d_rec = {}
    t_rec = []
    t_ids = []
    for line in fIn:
        #print(line)
        if line[0] == ">":
            if t_rec != []:
                
                
                header = t_rec[0]
                
                
                sq_id = header.split()[0][1:]
                
                t_ids.append(sq_id)
                

                if sq_id not in d_rec:
                    d_rec[sq_id] = t_rec
                #print(sq_id,len(sq), sq)
                
                t_rec = []
               
        t_rec.append(line)
    
    
    if t_rec != []:
        
        
        header = t_rec[0]
        
        
        sq_id = header.split()[0][1:]
        
        t_ids.append(sq_id)
        
    
            
    
        if sq_id not in d_rec:
            d_rec[sq_id] = t_rec
        #print(sq_id,len(sq), sq)

    return [d_rec, t_ids]










def CheckGB_complete(path, fInName, label):
    
    
    fIn = open(path+fInName,"r")
    lines = fIn.readlines()
    fIn.close()
    
    n = len(lines)

    
    
    
    # Count LOCUS entries
    fOut = open(path+label+"_LOCUS.txt", "w")
    ct = 0
    for i in range(n):
        line = lines[i]
        
        if line[0:5] == "LOCUS":
        #if line == "//\n":
            ct += 1   
            print(ct, i, line.strip(), sep = "\t", file = fOut)
    fOut.close()





    # Count // entries
    fOut = open(path+label+"_end.txt", "w")
    ct = 0
    for i in range(n):
        line = lines[i]
        

        if line == "//\n":
            ct += 1
            print(ct, i, line.strip(), sep = "\t", file = fOut)
            
            '''
            st = ""
            if lines[min(i+1,n-1)] == "\n":
                st = lines[min(i+2,n-1)]
            else:
                st = lines[min(i+1,n-1)]
                
            print(ct, i, st.strip(), sep = "\t", file = fOut)
            '''
    fOut.close()
    
    
    
    
    # Count ORIGIN entries
    
    fOut = open(path+label+"_ORIGIN", "w")
    ct = 0
    for i in range(n):
        line = lines[i]
        
        if "ORIGIN" in line:
            ct += 1
            
            print(ct, i, line.strip(), sep = "\t", file = fOut)
    fOut.close()
    
    return



def GetGBRecords(fInName):

    tt = []
    
    t_record = []
    
    fIn = open(fInName, "r")
    lines = fIn.readlines()
    fIn.close()
    
    
    
    for line in lines:
        #print(line)
        if line != "\n":
            tt.append(line)
        if line == "//\n":
            if tt != []:
                #print(len(tt))
                t_record.append(tt)
                tt = []
    
        
    
    return t_record

def GetGBFeatures(fInName):
    

    t_record = GetGBRecords(fInName)
    
    t_fts = []
    
    for item in t_record:
        t_lc = item[0].strip().split()
        #print(t_lc)
        g_acc = t_lc[1]
        g_sz = int(t_lc[2])
        g_base = t_lc[4]
        g_type = t_lc[5]
        #print(g_acc, g_sz, g_base, g_type)
        
        
        j = -1
        j_f = -1
        j_v = -1
        j_s = -1
        for ln in item:
            
    
            if ln[0:7] == "VERSION":
                j_v = item.index(ln)
            
            if ln[0:6] == "SOURCE":
                j_s = item.index(ln)
    
            if ln[0:8] == "FEATURES":
                j_f = item.index(ln)
            
            
            if "ORIGIN" in ln:
                j = item.index(ln)
                #print(item[j-1:j+2])
                break
    
    
        g_ver = item[j_v].strip().split()[-1]
        g_source = " ".join(item[j_s].strip().split()[1:])
        
        
        if g_ver == "VERSION":
            g_id = g_acc
        else:
            g_id = g_ver
    
    
        t = item[j_f+1:j]
        
        t_features = []
        tt = []
        for line in t:
            if line[5] != " ":
                #print(line.strip())
                if tt != []:
                    t_features.append(tt)
                    # The bug is caught (prev: tt = [line])
                    tt = []
            tt.append(line)
        
        
        if tt != []:
            t_features.append(tt)
            
        if t_features != []:
            t_fts.append(t_features)

    return t_fts

        

def GetField(t,i,s='"\n'):
    if i <0:
        return []
    j = -1
    #print("----------")
    #print(t,i, j)
    #print()
    for it in t[i:]:
        #print(it)
        if s in it:
            j = t.index(it)
            break
    #print(j)    
    return t[i:j+1]


def ExtractField(t):
    st = ""
    
    if t != []:
        n = len(t)
        
        #print(t)
        #print(n)
    
        
        if n == 1:
            ln0 = t[0]
            t_ln0 = ln0.strip().split('=')
            st = t_ln0[-1][1:-1]
        
        else:
            l1 = t[0]
            l1 = l1.strip()
            #print(l1)
            st = l1.split("=")[-1][1:]
            
            #print(st)
            
            for ln in t[1:]:
                ln = ln.strip()
                #print(ln)
                t_ln = ln.strip().split()
                if len(t_ln) > 1:
                    st = st+" "+ln
                else:
                    st = st+ln
            
            st = st[:-1]
             
    return st





def SeqGB_fasta(fInName):

    
    t_record = GetGBRecords(fInName)
    t_faa = []

    sep = "\t"
    sp = "|"
    
    #fOutName = fInName.split(".")[0]+"_sizecheck.txt"
    #fOut = open(fOutName, "w") 
    #print("id","source","rec_size","gbase","gtype","download_size","size_mismatch", sep = sep, file = fOut)
    
    c_mm = 0
    for item in t_record:
        t_lc = item[0].strip().split()
        #print(t_lc)
        g_acc = t_lc[1]
        g_sz = int(t_lc[2])
        g_base = t_lc[4]
        g_type = t_lc[5]
        
        
        
        j = -1
        j_v = -1
        j_s = -1
        for ln in item:
            if ln[0:7] == "VERSION":
                j_v = item.index(ln)
            
            if ln[0:6] == "SOURCE":
                j_s = item.index(ln)
            
            if "ORIGIN" in ln:
                j = item.index(ln)
                #print(item[j-1:j+2])
                break
        t_seq = item[j+1:-1]
        #print(t_seq[0]+t_seq[-1])
        
        g_ver = item[j_v].strip().split()[-1]
        g_source = " ".join(item[j_s].strip().split()[1:])
        

        if g_ver == "VERSION":
            g_id = g_acc
        else:
            g_id = g_ver
        
        
        
        #print(g_source)
        
        #print(g_acc, g_ver)
    
        seq = ""
        for l in t_seq:
            t_l = l.strip().split()
            st = "".join(t_l[1:])
            #print(l)
            #print(st)
            seq = seq+st
        n_seq = len(seq)
        if n_seq != g_sz:
            fl = 1
            c_mm += 1
        else:
            fl = 0
        #print(g_sz, len(seq))  
        
        if fl == 1:
            print(g_id, "size_mismatch!")
            return
        header = ">"+g_id+sp+g_source +sp + str(g_sz) +sp+ g_base +sp+g_type+"\n"
        t_faa.append([header, seq.upper()])
        
        
        #print(g_id, g_source, g_sz,g_base, g_type,n_seq,fl, sep = sep, file = fOut)
        #print(item[0].strip().split()[1],item[j].strip(),item[-1].strip())
        #break   
        
    
    #fOut.close()
    
    #print("Size mismatch:", c_mm)
    
    return t_faa




def SeqGB_sizecheck_list(fInName,gbpath=""):

    
    t_list = cm.GetListFromFile(fInName)

    t_record = []    

    for it in t_list:
        t_rec = GetGBRecords(gbpath+it)

        for rec in t_rec:
            t_record.append(rec)


    #t_record = GetGBRecords(fInName)


    sep = "\t"
    
    fOutName = fInName.split(".")[0]+"_sizecheck.txt"
    fOut = open(fOutName, "w") 
    print("id","source","rec_len","gbase","gtype","download_len","len_mismatch", "gc","Not_ACGT", sep = sep, file = fOut)
    
    c_mm = 0
    for item in t_record:
        t_lc = item[0].strip().split()
        #print(t_lc)
        g_acc = t_lc[1]
        g_sz = int(t_lc[2])
        g_base = t_lc[4]
        g_type = t_lc[5]
        
        
        
        j = -1
        j_v = -1
        j_s = -1
        for ln in item:
            if ln[0:7] == "VERSION":
                j_v = item.index(ln)
            
            if ln[0:6] == "SOURCE":
                j_s = item.index(ln)
            
            if ln[0:6] == "ORIGIN":
                j = item.index(ln)
                #print(item[j-1:j+2])
                break
        t_seq = item[j+1:-1]
        #print(t_seq[0]+t_seq[-1])
        
        g_ver = item[j_v].strip().split()[-1]
        g_source = " ".join(item[j_s].strip().split()[1:])
        

        if g_ver == "VERSION":
            g_id = g_acc
        else:
            g_id = g_ver
        
        
        
        #print(g_source)
        
        #print(g_acc, g_ver)
    
        seq = ""
        for l in t_seq:
            t_l = l.strip().split()
            st = "".join(t_l[1:])
            #print(l)
            #print(st)
            seq = seq+st
        n_seq = len(seq)
        if n_seq != g_sz:
            fl = 1
            c_mm += 1
        else:
            fl = 0
        #print(g_sz, len(seq)) 
        res = CalcGC(seq)
        print(g_id, g_source, g_sz,g_base, g_type,n_seq,fl,res[0], res[1], sep = sep, file = fOut)
        #print(item[0].strip().split()[1],item[j].strip(),item[-1].strip())
        #break   
        
    
    fOut.close()
    
    print("Size mismatch:", c_mm)
    
    return fOutName


















def SeqGB_sizecheck(fInName):

    
    t_record = GetGBRecords(fInName)


    sep = "\t"
    
    fOutName = fInName.split(".")[0]+"_sizecheck.txt"
    fOut = open(fOutName, "w") 
    print("id","source","rec_len","gbase","gtype","download_len","len_mismatch", "gc","Not_ACGT", sep = sep, file = fOut)
    
    c_mm = 0
    for item in t_record:
        t_lc = item[0].strip().split()
        #print(t_lc)
        g_acc = t_lc[1]
        g_sz = int(t_lc[2])
        g_base = t_lc[4]
        g_type = t_lc[5]
        
        
        
        j = -1
        j_v = -1
        j_s = -1
        for ln in item:
            if ln[0:7] == "VERSION":
                j_v = item.index(ln)
            
            if ln[0:6] == "SOURCE":
                j_s = item.index(ln)
            
            if ln[0:6] == "ORIGIN":
                j = item.index(ln)
                #print(item[j-1:j+2])
                break
        t_seq = item[j+1:-1]
        #print(t_seq[0]+t_seq[-1])
        
        g_ver = item[j_v].strip().split()[-1]
        g_source = " ".join(item[j_s].strip().split()[1:])
        

        if g_ver == "VERSION":
            g_id = g_acc
        else:
            g_id = g_ver
        
        
        
        #print(g_source)
        
        #print(g_acc, g_ver)
    
        seq = ""
        for l in t_seq:
            t_l = l.strip().split()
            st = "".join(t_l[1:])
            #print(l)
            #print(st)
            seq = seq+st
        n_seq = len(seq)
        if n_seq != g_sz:
            fl = 1
            c_mm += 1
        else:
            fl = 0
        #print(g_sz, len(seq)) 
        res = CalcGC(seq)
        print(g_id, g_source, g_sz,g_base, g_type,n_seq,fl,res[0], res[1], sep = sep, file = fOut)
        #print(item[0].strip().split()[1],item[j].strip(),item[-1].strip())
        #break   
        
    
    fOut.close()
    
    print("Size mismatch:", c_mm)
    
    return fOutName


def SeqGB_cutter(fInName, path = ""):
    # Updated to properly handle both new gb and old gb with GI (Feb 9, 2025)
    t_record = GetGBRecords(path+fInName)
    
    n_rec = len(t_record)
    
    subdir = str(n_rec)+"_gb/"
    
    new_dir = path+subdir
    

    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
        print(new_dir)


    fLogName = path+subdir+"flist.txt"
    
    fLog = open(fLogName,"w")
    
    
    for item in t_record:
        t_lc = item[0].strip().split()
        #print(t_lc)
        g_acc = t_lc[1]
        
    
        j_v = -1
        for ln in item:
            if ln[0:7] == "VERSION":
                j_v = item.index(ln)
                break
    
        # updated line to handle both new gb and old gb with GI
        g_ver = item[j_v].strip().split()[1]
        #g_ver = item[j_v].strip().split()[-1]
     
        
        fLog.write(g_ver+"\n")
        
        fOut = open(path+subdir+g_ver+".gb","w")
        for it in item:
            fOut.write(it)
        fOut.close()
        #break
    fLog.close()
    
    return fLogName



# Cutter for 2013 .gb format when the version field contained GI as follows:
# VERSION     NC_020859.1  GI:472341954
def SeqGB_cutter_GI_old(fInName, path = ""):
    
    t_record = GetGBRecords(path+fInName)
    
    n_rec = len(t_record)
    
    subdir = str(n_rec)+"_gb/"
    
    new_dir = path+subdir
    

    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
        print(new_dir)


    fLogName = path+subdir+"flist.txt"
    
    fLog = open(fLogName,"w")
    
    
    for item in t_record:
        t_lc = item[0].strip().split()
        #print(t_lc)
        g_acc = t_lc[1]
        
    
        j_v = -1
        for ln in item:
            if ln[0:7] == "VERSION":
                j_v = item.index(ln)
                break
    
        # Corrected!!!
        g_ver = item[j_v].strip().split()[1]
        # g_ver = item[j_v].strip().split()[-1]
     
        
        fLog.write(g_ver+"\n")
        
        fOut = open(path+subdir+g_ver+".gb","w")
        for it in item:
            fOut.write(it)
        fOut.close()
        #break
    fLog.close()
    
    return fLogName



def SeqDB_joinCDS_update(fInName, path, fListName = ""):
    
    
    t_list = []
    if fListName != "":      
        t_list = cm.GetListFromFile(path+fListName)
        print(len(t_list))
    
    
    #label = "log"
    #CheckGB_complete(path, fInName, label)
    t_record = GetGBRecords(path+fInName)
    
    print("Number of gb records:", len(t_record))
    
    sep = "\t"
    
    fOut = open(path+fInName.split(".")[0]+"_features.txt","w")
    print("g_acc","genes","cds","features", sep = sep, file = fOut)
    
    

    
    fOut3 = open(path+fInName.split(".")[0]+"_cds.txt","w")
    print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut3)
    
    
    
    fOut4 = open(path+fInName.rsplit(".")[0]+"_cds.faa","w")
    
    fLogName = path + "MissingTranslationField_1.txt"
    fLog = open(fLogName,"w")
    fLog.close()
    
    n_cds = 0
    
    n_psg = 0
    
    for item in t_record:
        t_lc = item[0].strip().split()
        print(t_lc)
        g_acc = t_lc[1]
        g_sz = int(t_lc[2])
        g_base = t_lc[4]
        g_type = t_lc[5]
        print(g_acc, g_sz, g_base, g_type)
        
        
        j = -1
        j_f = -1
        j_v = -1
        j_s = -1
        for ln in item:
            
    
            if ln[0:7] == "VERSION":
                j_v = item.index(ln)
            
            if ln[0:6] == "SOURCE":
                j_s = item.index(ln)
    
            if ln[0:8] == "FEATURES":
                j_f = item.index(ln)
            
            
            if "ORIGIN" in ln:
                j = item.index(ln)
                #print(item[j-1:j+2])
                break
    
    
        g_ver = item[j_v].strip().split()[-1]
        g_source = " ".join(item[j_s].strip().split()[1:])
        
        
        if g_ver == "VERSION":
            g_id = g_acc
        else:
            g_id = g_ver
    
    
        t = item[j_f+1:j]
        
        t_features = []
        tt = []
        for line in t:
            if line[5] != " ":
                #print(line.strip())
                if tt != []:
                    t_features.append(tt)
                    # The bug is caught (prev: tt = [line])
                    tt = []
            tt.append(line)
        
        
        if tt != []:
            t_features.append(tt)
            
        
        #print(len(t_features))  
        #print(t_features[-1])
        
        
        ct_cds = 0
        ct_gene = 0
        
        ct = 0
        
        t_cds_index = []
        t_cds = []
        
        
        for it in t_features:
            ln = it[0]
            #print(ln.strip())
            ct += 1
    
            if ln.strip()[0:4] == "gene":
                ct_gene += 1
                #print(ln.strip())
    
            if ln.strip()[0:3] == "CDS":
                ct_cds += 1
                
                t_cds_index.append(t_features.index(it))
                t_cds.append(it)
                
                # Check if the last field in cds feature is always translation
                # Not always, exception: NC_009799.3.gb, CDS             9713..14527 /db_xref="GeneID:5580369"
                i_l = -1
                fl = 0
                for l in it:
                    if l.strip()[0] == "/":
                        i_l = it.index(l)
                        
                        #if "/pseudo" in l:
                        #    fl = 1
                        
                '''        
                if (it[i_l].strip()[0:13] != "/translation="):
                    print(g_acc, ln.strip(), it[i_l].strip())
                '''
                #print(ln.strip())
    
        '''
        print("Number of genes:", ct_gene)
        print("Number of CDS:", ct_cds)
        
        print((ct-1) / 2)
        '''
        print(g_id, ct_gene, ct_cds, ct, sep = sep, file = fOut)
        
        res = [g_id, ct_gene, ct_cds, ct]
        
        if len(t_cds) == 0:
            print(g_id, "The number of CDS equals 0 !!!!!!!")
        
        
        # Processing of CDS list   
     
        #print(len(t_cds_index))    
        #print(len(t_cds))
    
    
    
    
        #fOut1 = open(path+subdir+g_ver+"_cds.txt","w")
        #print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut1)
    
    
        #fOut2 = open(path+subdir+g_ver+"_cds.faa","w")
        #print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut1)
    
    
    
        j = 0
    
        n_cds += len(t_cds)
        
        pr_i = 0
        
        ct_missed = 0
        
        c_psg = 0
        
        k = 0
    
        for it in t_cds:
            #print(t_cds.index(it)+1, it[0].strip())
            
            
            j_pr = -1
            j_pid = -1
            j_tr = -1
            j_psg = -1
            
            
            
            cds_loc = it[0].strip().split()[-1]
            
           
            #print(cds_loc)  
    
    
            for ln in it:
    
    
                if ln.strip()[0:9] == "/product=":
                    j_pr = it.index(ln)
                    #print(ln.strip())
                if ln.strip()[0:12] == "/protein_id=":
                    j_pid = it.index(ln)
                    #print(ln.strip())
                if ln.strip()[0:13] == "/translation=":
                    j_tr = it.index(ln)
                    #print(ln.strip())
                if ln.strip()[0:7] == "/pseudo":
                    j_psg = it.index(ln)
                    #print(ln.strip())
                
                
                
            t_pr = GetField(it,j_pr)
            st_pr = ExtractField(t_pr)
            
            
            t_pid = GetField(it,j_pid)
            st_pid = ExtractField(t_pid)
            
            
            t_tr = GetField(it,j_tr)
            st_tr = ExtractField(t_tr)
            
            
            if st_pid == "":
                pr_i += 1
                st_pid = g_id+"_"+"CDS"+str(pr_i)
            

            
            '''   
            if cds_loc == "9713..14527":
                print(j_tr)
                print(t_tr)
                print(it[0], it[1])
                break
            '''
                          
            
            #print(st_pr, st_pid, st_tr)
            
            
            print(g_id, st_pid, st_pr, cds_loc, st_tr, len(st_tr), sep = sep, file=fOut3)
            
            
            sp = "|"
            header = ">"+g_id+sp+st_pid +sp + st_pr +sp+ cds_loc +sp+str(len(st_tr))+ " aa"
            
            if st_tr != "":
                if (t_list == []) or ((t_list != []) and (st_pid in t_list)):
                    
                    if st_tr[0] == "-":
                        print(st_tr[0:3], st_pid)
                        st_tr = "M"+st_tr[1:]
                        k += 1
                    
                    
                    
                    fOut4.write(header+"\n")
                    fOut4.write(st_tr+"\n")
            else:
                if j_psg == -1:          
                    print("missing translation!", g_id, st_pid, header )
                    fLog = open(fLogName,"a")
                    fLog.write("111 "+header[1:]+"\n")
                    fLog.close()
                    ct_missed += 1
                else:
                    c_psg += 1
            
            
            
            '''
            if len(t_tr) > 1:
                #print(t_record.index(item), item[0], t_pr)
                print(t_tr)
             ''' 
    
    
            
            
            #print(t_pr)
            #break
    
        '''
        for i in t_cds_index:
            print(t_features[i][0].strip())
    
        '''
        #fOut1.close()
        #fOut2.close()
        #break
        
        n_psg += c_psg
        print(g_acc, ", CDS with unusual start aa", k)
        
        
        
     
    print("Total number of cds: ", n_cds)
    print("Total number of pseudogenes: ", n_psg)
    
    
 
    
    fOut3.close()
    fOut4.close()
    
    fOut.close()
    return res+[ct_missed]






def SeqDB_joinCDS_update_stat(fInName, path, fListName = ""):
    # return stat as a dictionary
    
    
    t_stat = []
    
    
    t_list = []
    if fListName != "":      
        t_list = cm.GetListFromFile(path+fListName)
        print(len(t_list))

    t_record = GetGBRecords(path+fInName)
    
    print("Number of gb records:", len(t_record))
    
    t_stat.append(len(t_record))
    
    sep = "\t"
    
    #fOut = open(path+fInName.split(".")[0]+"_features.txt","w")
    #print("g_acc","genes","cds","features", sep = sep, file = fOut)
    
    

    
    #fOut3 = open(path+fInName.split(".")[0]+"_cds.txt","w")
    #print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut3)
    
    
    
    fOut4 = open(path+fInName.rsplit(".",1)[0]+"_cds.faa","w")
    
    #fLogName = path + "MissingTranslationField_1.txt"
    #fLog = open(fLogName,"w")
    #fLog.close()
    
    n_cds = 0
    
    n_psg = 0
    
    for item in t_record:
        t_lc = item[0].strip().split()
        print(t_lc)
        g_acc = t_lc[1]
        g_sz = int(t_lc[2])
        g_base = t_lc[4]
        g_type = t_lc[5]
        print(g_acc, g_sz, g_base, g_type)
        
        
        j = -1
        j_f = -1
        j_v = -1
        j_s = -1
        for ln in item:
            
    
            if ln[0:7] == "VERSION":
                j_v = item.index(ln)
            
            if ln[0:6] == "SOURCE":
                j_s = item.index(ln)
    
            if ln[0:8] == "FEATURES":
                j_f = item.index(ln)
            
            
            if "ORIGIN" in ln:
                j = item.index(ln)
                #print(item[j-1:j+2])
                break
    
    
        g_ver = item[j_v].strip().split()[-1]
        g_source = " ".join(item[j_s].strip().split()[1:])
        
        
        if g_ver == "VERSION":
            g_id = g_acc
        else:
            g_id = g_ver
    
    
        t = item[j_f+1:j]
        
        t_features = []
        tt = []
        for line in t:
            if line[5] != " ":
                #print(line.strip())
                if tt != []:
                    t_features.append(tt)
                    # The bug is caught (prev: tt = [line])
                    tt = []
            tt.append(line)
        
        
        if tt != []:
            t_features.append(tt)
            
        

        
        ct_cds = 0
        ct_gene = 0
        
        ct = 0
        
        t_cds_index = []
        t_cds = []
        
        
        for it in t_features:
            ln = it[0]
            #print(ln.strip())
            ct += 1
    
            if ln.strip()[0:4] == "gene":
                ct_gene += 1
                #print(ln.strip())
    
            if ln.strip()[0:3] == "CDS":
                ct_cds += 1
                
                t_cds_index.append(t_features.index(it))
                t_cds.append(it)
                
                # Check if the last field in cds feature is always translation
                # Not always, exception: NC_009799.3.gb, CDS             9713..14527 /db_xref="GeneID:5580369"
                i_l = -1
                fl = 0
                for l in it:
                    if l.strip()[0] == "/":
                        i_l = it.index(l)
                        
                        #if "/pseudo" in l:
                        #    fl = 1
                        
                '''        
                if (it[i_l].strip()[0:13] != "/translation="):
                    print(g_acc, ln.strip(), it[i_l].strip())
                '''
                #print(ln.strip())
    
        '''
        print("Number of genes:", ct_gene)
        print("Number of CDS:", ct_cds)
        
        print((ct-1) / 2)
        '''
        #print(g_id, ct_gene, ct_cds, ct, sep = sep, file = fOut)
        
        res = [g_id, ct_gene, ct_cds, ct]
        
        if len(t_cds) == 0:
            print(g_id, "The number of CDS equals 0 !!!!!!!")
        
        
        # Processing of CDS list   
     
        #print(len(t_cds_index))    
        #print(len(t_cds))
    
    
    
    
        #fOut1 = open(path+subdir+g_ver+"_cds.txt","w")
        #print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut1)
    
    
        #fOut2 = open(path+subdir+g_ver+"_cds.faa","w")
        #print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut1)
    
    
    
        j = 0
    
        n_cds += len(t_cds)
        
        pr_i = 0
        
        ct_missed = 0
        
        c_psg = 0
        
        k = 0
    
        for it in t_cds:
            #print(t_cds.index(it)+1, it[0].strip())
            
            
            j_pr = -1
            j_pid = -1
            j_tr = -1
            j_psg = -1
            
            
            
            cds_loc = it[0].strip().split()[-1]
            
           
            #print(cds_loc)  
    
    
            for ln in it:
    
    
                if ln.strip()[0:9] == "/product=":
                    j_pr = it.index(ln)
                    #print(ln.strip())
                if ln.strip()[0:12] == "/protein_id=":
                    j_pid = it.index(ln)
                    #print(ln.strip())
                if ln.strip()[0:13] == "/translation=":
                    j_tr = it.index(ln)
                    #print(ln.strip())
                if ln.strip()[0:7] == "/pseudo":
                    j_psg = it.index(ln)
                    #print(ln.strip())
                
                
                
            t_pr = GetField(it,j_pr)
            st_pr = ExtractField(t_pr)
            
            
            t_pid = GetField(it,j_pid)
            st_pid = ExtractField(t_pid)
            
            
            t_tr = GetField(it,j_tr)
            st_tr = ExtractField(t_tr)
            
            
            if st_pid == "":
                pr_i += 1
                st_pid = g_id+"_"+"CDS"+str(pr_i)
            


            
            #print(g_id, st_pid, st_pr, cds_loc, st_tr, len(st_tr), sep = sep, file=fOut3)
            
            
            sp = "|"
            header = ">"+g_id+sp+st_pid +sp + st_pr +sp+ cds_loc +sp+str(len(st_tr))+ " aa"
            
            if st_tr != "":
                if (t_list == []) or ((t_list != []) and (st_pid in t_list)):
                    
                    if st_tr[0] == "-":
                        print(st_tr[0:3], st_pid)
                        st_tr = "M"+st_tr[1:]
                        k += 1
                    
                    
                    
                    fOut4.write(header+"\n")
                    fOut4.write(st_tr+"\n")
            else:
                if j_psg == -1:          
                    print("missing translation!", g_id, st_pid, header )
                    #fLog = open(fLogName,"a")
                    #fLog.write("111 "+header[1:]+"\n")
                    #fLog.close()
                    ct_missed += 1
                else:
                    c_psg += 1
            
            
            
            '''
            if len(t_tr) > 1:
                #print(t_record.index(item), item[0], t_pr)
                print(t_tr)
             ''' 
    
    
            
            
            #print(t_pr)
            #break
    
        '''
        for i in t_cds_index:
            print(t_features[i][0].strip())
    
        '''
        #fOut1.close()
        #fOut2.close()
        #break
        
        n_psg += c_psg
        print(g_acc, ", CDS with unusual start aa", k)
        t_stat.append(k)
        
        
        
     
    print("Total number of cds: ", n_cds)
    print("Total number of pseudogenes: ", n_psg)
    
    t_stat.append(n_cds)
    t_stat.append(n_psg)
    
    #fOut3.close()
    fOut4.close()
    
    #fOut.close()
    #return res+[ct_missed]
    return t_stat












def SeqDB_joinCDS_update_stat_new(fInName, path, fListName = ""):
    # return stat as a dictionary
    
    
    t_stat = []
    
    
    t_list = []
    if fListName != "":      
        t_list = cm.GetListFromFile(path+fListName)
        print(len(t_list))

    t_record = GetGBRecords(path+fInName)
    
    print("Number of gb records:", len(t_record))
    
    t_stat.append(len(t_record))
    
    sep = "\t"
    
    #fOut = open(path+fInName.split(".")[0]+"_features.txt","w")
    #print("g_acc","genes","cds","features", sep = sep, file = fOut)
    
    

    
    #fOut3 = open(path+fInName.split(".")[0]+"_cds.txt","w")
    #print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut3)
    
    
    
    fOut4 = open(path+fInName.rsplit(".",1)[0]+"_cds.faa","w")
    
    #fLogName = path + "MissingTranslationField_1.txt"
    #fLog = open(fLogName,"w")
    #fLog.close()
    
    n_cds = 0
    
    n_psg = 0
    
    for item in t_record:
        t_lc = item[0].strip().split()
        print(t_lc)
        g_acc = t_lc[1]
        g_sz = int(t_lc[2])
        g_base = t_lc[4]
        g_type = t_lc[5]
        print(g_acc, g_sz, g_base, g_type)
        
        
        j = -1
        j_f = -1
        j_v = -1
        j_s = -1
        j_l = -1
        
        for ln in item:

            if ln[0:5] == "LOCUS":
                j_l = item.index(ln)            
    
            if ln[0:7] == "VERSION":
                j_v = item.index(ln)
            
            if ln[0:6] == "SOURCE":
                j_s = item.index(ln)
    
            if ln[0:8] == "FEATURES":
                j_f = item.index(ln)
            
            
            if "ORIGIN" in ln:
                j = item.index(ln)
                #print(item[j-1:j+2])
                break
    
    
        g_ver = item[j_v].strip().split()[-1]
        g_source = "_".join(item[j_s].strip().split()[1:])

        
        
        if g_ver == "VERSION":
            g_id = g_acc
        else:
            g_id = g_ver
    
    
        t = item[j_f+1:j]
        
        t_features = []
        tt = []
        for line in t:
            if line[5] != " ":
                #print(line.strip())
                if tt != []:
                    t_features.append(tt)
                    # The bug is caught (prev: tt = [line])
                    tt = []
            tt.append(line)
        
        
        if tt != []:
            t_features.append(tt)
            
        

        
        ct_cds = 0
        ct_gene = 0
        
        ct = 0
        
        t_cds_index = []
        t_cds = []
        
        
        for it in t_features:
            ln = it[0]
            #print(ln.strip())
            ct += 1
    
            if ln.strip()[0:4] == "gene":
                ct_gene += 1
                #print(ln.strip())
    
            if ln.strip()[0:3] == "CDS":
                ct_cds += 1
                
                t_cds_index.append(t_features.index(it))
                t_cds.append(it)
                
                # Check if the last field in cds feature is always translation
                # Not always, exception: NC_009799.3.gb, CDS             9713..14527 /db_xref="GeneID:5580369"
                i_l = -1
                fl = 0
                for l in it:
                    if l.strip()[0] == "/":
                        i_l = it.index(l)
                        
                        #if "/pseudo" in l:
                        #    fl = 1
                        
                '''        
                if (it[i_l].strip()[0:13] != "/translation="):
                    print(g_acc, ln.strip(), it[i_l].strip())
                '''
                #print(ln.strip())
    
        '''
        print("Number of genes:", ct_gene)
        print("Number of CDS:", ct_cds)
        
        print((ct-1) / 2)
        '''
        #print(g_id, ct_gene, ct_cds, ct, sep = sep, file = fOut)
        
        res = [g_id, ct_gene, ct_cds, ct]
        
        if len(t_cds) == 0:
            print(g_id, "The number of CDS equals 0 !!!!!!!")
        
        
        # Processing of CDS list   
     
        #print(len(t_cds_index))    
        #print(len(t_cds))
    
    
    
    
        #fOut1 = open(path+subdir+g_ver+"_cds.txt","w")
        #print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut1)
    
    
        #fOut2 = open(path+subdir+g_ver+"_cds.faa","w")
        #print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut1)
    
    
    
        j = 0
    
        n_cds += len(t_cds)
        
        pr_i = 0
        
        ct_missed = 0
        
        c_psg = 0
        
        k = 0
    
    
        n_cds = len(t_cds)
        
        for jj in range(n_cds):
            
            it = t_cds[jj]
            
            #print(t_cds.index(it)+1, it[0].strip())
            
            
            j_pr = -1
            j_pid = -1
            j_tr = -1
            j_psg = -1
            
            
            
            cds_loc = it[0].strip().split()[-1]
            
           
            #print(cds_loc)  
    
    
            for ln in it:
    
    
                if ln.strip()[0:9] == "/product=":
                    j_pr = it.index(ln)
                    #print(ln.strip())
                if ln.strip()[0:12] == "/protein_id=":
                    j_pid = it.index(ln)
                    #print(ln.strip())
                if ln.strip()[0:13] == "/translation=":
                    j_tr = it.index(ln)
                    #print(ln.strip())
                if ln.strip()[0:7] == "/pseudo":
                    j_psg = it.index(ln)
                    #print(ln.strip())
                
                
                
            t_pr = GetField(it,j_pr)
            st_pr = ExtractField(t_pr)
            
            
            t_pid = GetField(it,j_pid)
            st_pid = ExtractField(t_pid)
            
            
            t_tr = GetField(it,j_tr)
            st_tr = ExtractField(t_tr)
            
            
            if st_pid == "":
                pr_i += 1
                st_pid = g_id+"_"+"CDS"+str(pr_i)
            


            
            #print(g_id, st_pid, st_pr, cds_loc, st_tr, len(st_tr), sep = sep, file=fOut3)
            
            
            sp = "|"
            #header = ">"+g_id+sp+st_pid +sp + st_pr +sp+ cds_loc +sp+str(len(st_tr))+ " aa"
            
            header = ">"+g_id+sp+st_pid+"_"+str(jj+1) +sp + "["+g_source+"]"+"_".join(st_pr.split()) +sp+ cds_loc +sp+str(len(st_tr))+ " aa"
            
            
            if st_tr != "":
                if (t_list == []) or ((t_list != []) and (st_pid in t_list)):
                    
                    if st_tr[0] == "-":
                        print(st_tr[0:3], st_pid)
                        st_tr = "M"+st_tr[1:]
                        k += 1
                    
                    
                    
                    fOut4.write(header+"\n")
                    fOut4.write(st_tr+"\n")
            else:
                if j_psg == -1:          
                    print("missing translation!", g_id, st_pid, header )
                    #fLog = open(fLogName,"a")
                    #fLog.write("111 "+header[1:]+"\n")
                    #fLog.close()
                    ct_missed += 1
                else:
                    c_psg += 1
            
            
            
            '''
            if len(t_tr) > 1:
                #print(t_record.index(item), item[0], t_pr)
                print(t_tr)
             ''' 
    
    
            
            
            #print(t_pr)
            #break
    
        '''
        for i in t_cds_index:
            print(t_features[i][0].strip())
    
        '''
        #fOut1.close()
        #fOut2.close()
        #break
        
        n_psg += c_psg
        print(g_acc, ", CDS with unusual start aa", k)
        t_stat.append(k)
        
        
        
     
    print("Total number of cds: ", n_cds)
    print("Total number of pseudogenes: ", n_psg)
    
    t_stat.append(n_cds)
    t_stat.append(n_psg)
    
    #fOut3.close()
    fOut4.close()
    
    #fOut.close()
    #return res+[ct_missed]
    return t_stat














#---------------- gz



def GetGBRecords_gz(fInName):

    tt = []
    
    t_record = []
    
    
    with gzip.open(fInName, "rt") as handle:
        for line in handle:
            if line != "\n":
                tt.append(line)
            if line == "//\n":
                if tt != []:
                    #print(len(tt))
                    t_record.append(tt)
                    tt = []

    return t_record




def SeqDB_joinCDS_gz(fInName, path, outpath = ""):
    #label = "log"
    #CheckGB_complete(path, fInName, label)
    t_record = GetGBRecords_gz(path+fInName)
    
    print("Number of gb records:", len(t_record))
    
    sep = "\t"
    
    prefix = fInName.rsplit(".",2)[0]
    
    
    fOut = open(path+prefix+"_features.txt","w")
    print("g_acc","seq","genes","cds","features-cds", sep = sep, file = fOut)
    
    

    
    fOut3 = open(path+prefix+"_cds.txt","w")
    print("Genome_id", "Sequence_id","Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut3)
    
    
    
    fOut4 = open(outpath+prefix+"_cds.faa","w")
    
    fLogName = path + "MissingTranslationField.txt"
    fLog = open(fLogName,"w")
    fLog.close()
    
    n_cds = 0
    
    for item in t_record:
        t_lc = item[0].strip().split()
        #print(t_lc)
        g_acc = t_lc[1]
        g_sz = int(t_lc[2])
        g_base = t_lc[4]
        g_type = t_lc[5]
        
        
        j = -1
        j_f = -1
        j_v = -1
        j_s = -1
        for ln in item:
            
    
            if ln[0:7] == "VERSION":
                j_v = item.index(ln)
            
            if ln[0:6] == "SOURCE":
                j_s = item.index(ln)
    
            if ln[0:8] == "FEATURES":
                j_f = item.index(ln)
            
            
            if "ORIGIN" in ln:
                j = item.index(ln)
                #print(item[j-1:j+2])
                break
    
    
        g_ver = item[j_v].strip().split()[-1]
        g_source = " ".join(item[j_s].strip().split()[1:])
        
        
        if g_ver == "VERSION":
            g_id = g_acc
        else:
            g_id = g_ver
    
    
        t = item[j_f+1:j]
        
        t_features = []
        tt = []
        for line in t:
            if line[5] != " ":
                #print(line.strip())
                if tt != []:
                    t_features.append(tt)
                    # The bug is caught (prev: tt = [line])
                    tt = []
            tt.append(line)
        
        
        if tt != []:
            t_features.append(tt)
            
        
        #print(len(t_features))  
        #print(t_features[-1])
        
        
        ct_cds = 0
        ct_gene = 0
        
        ct = 0
        
        t_cds_index = []
        t_cds = []
        
        
        for it in t_features:
            ln = it[0]
            #print(ln.strip())
            ct += 1
    
            if ln.strip()[0:4] == "gene":
                ct_gene += 1
                #print(ln.strip())
    
            if ln.strip()[0:3] == "CDS":
                ct_cds += 1
                
                t_cds_index.append(t_features.index(it))
                t_cds.append(it)
                
                # Check if the last field in cds feature is always translation
                # Not always, exception: NC_009799.3.gb, CDS             9713..14527 /db_xref="GeneID:5580369"
                i_l = -1
                for l in it:
                    if l.strip()[0] == "/":
                        i_l = it.index(l)
                        
                if it[i_l].strip()[0:13] != "/translation=":
                    print(g_acc, ln.strip(), it[i_l].strip())
                #print(ln.strip())
    
        '''
        print("Number of genes:", ct_gene)
        print("Number of CDS:", ct_cds)
        
        print((ct-1) / 2)
        '''
        print(prefix, g_id, ct_gene, ct_cds, ct-ct_cds, sep = sep, file = fOut)
        
     
        
        if len(t_cds) == 0:
            print(prefix, g_id, "The number of CDS equals 0 !!!!!!!")
        
        
        # Processing of CDS list   
     
        #print(len(t_cds_index))    
        #print(len(t_cds))
    
    
    
    
        #fOut1 = open(path+subdir+g_ver+"_cds.txt","w")
        #print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut1)
    
    
        #fOut2 = open(path+subdir+g_ver+"_cds.faa","w")
        #print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut1)
    
    
    
        j = 0
    
        n_cds += len(t_cds)
        
        pr_i = 0
    
        for it in t_cds:
            #print(t_cds.index(it)+1, it[0].strip())
            
            
            j_pr = -1
            j_pid = -1
            j_tr = -1
            
            cds_loc = it[0].strip().split()[-1]
            
           
            #print(cds_loc)  
    
    
            for ln in it:
    
    
                if ln.strip()[0:9] == "/product=":
                    j_pr = it.index(ln)
                    #print(ln.strip())
                if ln.strip()[0:12] == "/protein_id=":
                    j_pid = it.index(ln)
                    #print(ln.strip())
                if ln.strip()[0:13] == "/translation=":
                    j_tr = it.index(ln)
                    #print(ln.strip())
                    
    
                
                
                
            t_pr = GetField(it,j_pr)
            st_pr = ExtractField(t_pr)
            
            
            t_pid = GetField(it,j_pid)
            st_pid = ExtractField(t_pid)
            
            
            t_tr = GetField(it,j_tr)
            st_tr = ExtractField(t_tr)
            
            
            if st_pid == "":
                pr_i += 1
                st_pid = g_id+"_"+"CDS"+str(pr_i)
            
            
            '''   
            if cds_loc == "9713..14527":
                print(j_tr)
                print(t_tr)
                print(it[0], it[1])
                break
            '''
                          
            
            #print(st_pr, st_pid, st_tr)
            
            
            print(prefix,g_id, st_pid, st_pr, cds_loc, st_tr, len(st_tr), sep = sep, file=fOut3)
            
            
            sp = "|"
            header = ">"+prefix+sep+g_id+sp+st_pid +sp + st_pr +sp+ cds_loc +sp+str(len(st_tr))+ " aa"
            
            if st_tr != "":
                fOut4.write(header+"\n")
                fOut4.write(st_tr+"\n")
            else:
                print("missing translation!", g_id, st_pid, header )
                fLog = open(fLogName,"a")
                fLog.write("111 "+header[1:]+"\n")
                fLog.close()
            
            
            
            '''
            if len(t_tr) > 1:
                #print(t_record.index(item), item[0], t_pr)
                print(t_tr)
             ''' 
    
    
            
            
            #print(t_pr)
            #break
    
        '''
        for i in t_cds_index:
            print(t_features[i][0].strip())
    
        '''
        #fOut1.close()
        #fOut2.close()
        #break
        
     
    print("Total number of cds: ", n_cds)
    
 
    
    fOut3.close()
    fOut4.close()
    
    fOut.close()
    return


#-----------------------------------------------



def SeqDB_update(fInName, path):

    t_record = GetGBRecords(path+fInName)
    
    print("Number of gb records:", len(t_record))
    
    sep = "\t"
    

    n_rec = len(t_record)
    c_rec = 0
    
    if n_rec == 1:
    
        for item in t_record:
            c_rec += 1
            
            n_item = len(item)
            t_lc = item[0].strip().split()
            #print(t_lc)
            g_acc = t_lc[1]
            g_sz = int(t_lc[2])
            g_base = t_lc[4]
            g_type = t_lc[5]
            
            
            j = -1
            j_f = -1
            j_v = -1
            j_s = -1
            for ln in item:
                
        
                if ln[0:7] == "VERSION":
                    j_v = item.index(ln)
                
                if ln[0:6] == "SOURCE":
                    j_s = item.index(ln)
        
                if ln[0:8] == "FEATURES":
                    j_f = item.index(ln)
                
                
                if "ORIGIN" in ln:
                    j = item.index(ln)
                    #print(item[j-1:j+2])
                    break
        
        
            g_ver = item[j_v].strip().split()[-1]
            g_source = " ".join(item[j_s].strip().split()[1:])
        
            # extract FEATURES section
            t = item[j_f+1:j]
            
            #print(t)
            
            if g_ver == "VERSION":
                g_id = g_acc
            else:
                g_id = g_ver
        
            
            
            
            
            
            # separates individual FEATURES
            t_features = []
            tt = []
            for line in t:
                if line[5] != " ":
                    #print(line.strip())
                    if tt != []:
                        t_features.append(tt)
                        # The bug is caught (prev: tt = [line])
                        tt = []
                tt.append(line)
            
            
            if tt != []:
                t_features.append(tt)
                
            
            print(len(t_features))  
            #print(t_features[-1])
            
            
            
            
            # extract CDS
            ct_cds = 0
            ct_gene = 0
            
            ct = 0
            
            t_cds_index = []
            t_cds = []
            
    
            t_cds_pos = []
            
            t_new_cds = []
            
            for it in t_features:
                ln = it[0]
                #print(ln.strip())
                ct += 1
        
                if ln.strip()[0:4] == "gene":
                    ct_gene += 1
                    #print(ln.strip())
        
                if ln.strip()[0:3] == "CDS":
                    ct_cds += 1
                    
                    
                    
                    t_cds_index.append(t_features.index(it))
                    t_cds_pos.append([item.index(it[0]), item.index(it[-1])])
                    t_cds.append(it)
                    
                 
                    k_product = -1
                    k_proteinid = -1
                    
                    for item1 in it:
                        if "/product=" in item1:
                            k_product = it.index(item1)
                        if "/protein_id=" in item1:
                            k_proteinid = it.index(item1)
                    print(k_product, k_proteinid)
                    k = -1
                    if k_proteinid < 0:
                        #k = -1 
                        for itt in it[k_product+1:]:
                            if itt.strip()[0] == "/":
                                k = it.index(itt)
                                break
                    print(k_product, k)
                    kk = it[k_product].find("/")
                    st_proteinid = it[k_product][:kk] +'/protein_id="'+g_id+'_CDS'+str(ct_cds)+'"\n'
                    new_it = it[:k] + [st_proteinid] + it[k:]
                    
                    
                    t_new_cds.append(new_it)
                    
                    
                    #k_insert = item.index(it[k])
                        
                    #print(k_insert)
                    
                    #print(t_cds_pos)
                                         
                    
                    
                    
                    
                    # Check if the last field in cds feature is always translation
                    # Not always, exception: NC_009799.3.gb, CDS             9713..14527 /db_xref="GeneID:5580369"
                    i_l = -1
                    for l in it:
                        if l.strip()[0] == "/":
                            i_l = it.index(l)
                            
                    if it[i_l].strip()[0:13] != "/translation=":
                        print(g_acc, ln.strip(), it[i_l].strip())
                    #print(ln.strip())
                    print(it)
            
                    print(new_it)
        
        
                    #break
                
            
            c_cds = 0
            #fOut = open(path+"gbrecord"+str(c_rec)+".gb","w")
            
            fOutName = path+fInName.split(".")[0]+"_updated.gb"
            fOut = open(fOutName, "w")
            
            ii = 0
            while ii < n_item:
                if (ii <= j_f) or (ii >= j):
                    fOut.write(item[ii])
                else:
                    
                    if item[ii].strip()[0:3] == "CDS":
                        cds_pos = t_cds_pos[c_cds]
                        new_cds = t_new_cds[c_cds]
                        for itt in new_cds:
                            fOut.write(itt)
                        c_cds += 1
                        ii = cds_pos[1]
                    else:
                        fOut.write(item[ii])
                ii += 1
                        
                    
            
            
            fOut.close()
            
            print(c_cds)
            

    return fOutName








#-----------------------------------------------




