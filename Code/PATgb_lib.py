"""
@author: Tatiana Lenskaia
All rights reserved, 2025
E-mail: lensk010@gmail.com
"""


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


def SeqDB_joinCDS_update_stat_old_GI(fInName, path, fListName = ""):
    # return stat as a dictionary
    
    
    t_stat = []
    
    
    t_list = []
    if fListName != "":      
        t_list = GetListFromFile(path+fListName)
        print(len(t_list))

    t_record = GetGBRecords(path+fInName)
    
    #print("Number of gb records:", len(t_record))
    
    t_stat.append(len(t_record))
    
    sep = "\t"
    
    #fOut = open(path+fInName.split(".")[0]+"_features.txt","w")
    #print("g_acc","genes","cds","features", sep = sep, file = fOut)
    
    

    
    #fOut3 = open(path+fInName.split(".")[0]+"_cds.txt","w")
    #print("Genome_id", "Protein_id", "GB_annot", "Location", "Sequence", "Seq_len",sep=sep, file=fOut3)
    
    
    
    fOut4 = open(path+fInName.split(".")[0]+"_cds.faa","w")
    
    #fLogName = path + "MissingTranslationField_1.txt"
    #fLog = open(fLogName,"w")
    #fLog.close()
    
    n_cds = 0
    
    n_psg = 0
    
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
    
    
        g_ver = item[j_v].strip().split()[1]
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
        #print(g_acc, ", CDS with unusual start aa", k)
        t_stat.append(k)
        
        
        
     
    #print("Total number of cds: ", n_cds)
    #print("Total number of pseudogenes: ", n_psg)
    
    t_stat.append(n_cds)
    t_stat.append(n_psg)
    
    #fOut3.close()
    fOut4.close()
    
    #fOut.close()
    #return res+[ct_missed]
    return t_stat
