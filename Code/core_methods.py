#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Tatiana Lenskaia, email: lensk010@umn.edu
All rights reserved, 2021
"""

import random
import copy as cp






def ConcatFiles(fListName, fOutName, path = ""):
    
    fOut = open(path+fOutName, "w")
    
    
    fList = open(fListName,"r")
    for line in fList:
        line = line.strip()
        fIn = open(path+line, "r")
        for ln in fIn:
            fOut.write(ln)
        fIn.close()
        #break
    fList.close()
    
    fOut.close()
    
    return fOutName








def GetRandomSeq_unif(n_len, base = "ACGT"):
    seq = ""
    num_letters = len(base)
    for i in range(n_len):
        num = random.randint(0,num_letters-1)
        seq = seq+base[num]
    return seq

def Rev_seq(seq):
    rev_seq = seq[::-1]
    return rev_seq

def Rev_cmp(st):
    #st = st.lower()
    cmp_st = st.translate(str.maketrans("ACGT","TGCA"))
    rev_cmp_st = cmp_st[::-1]
    return rev_cmp_st;


def Cmp(st):
    #st = st.lower()
    cmp_st = st.translate(str.maketrans("ACGT","TGCA"))
    return cmp_st;



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


def CompareLists(fInName1, fInName2):
    t1 = GetListFromFile(fInName1)
    t2 = GetListFromFile(fInName2)

    n1 = len(t1)
    n2 = len(t2)
    
    
    print(n1, fInName1)
    print(n2, fInName2)
    
    t = []
    t1_rem = cp.copy(t1)
    t2_rem = cp.copy(t2)
    

    if n1 < n2:
        for it in t1:
            
            if it in t2:
                t2_rem.remove(it)
                t.append(it)
                t1_rem.remove(it)
    
    else:
        for it in t2:
            if it in t1:
    
                t1_rem.remove(it)
                t.append(it)
                t2_rem.remove(it)    
     
    n1_r = len(t1_rem)
    n2_r = len(t2_rem)
    n = len(t)

    
    fOutName = fInName1.rsplit("/",1)[0]+"/compare_lists.txt"
    fOut = open(fOutName,"w")
    
    print(n1, fInName1, file = fOut)
    
    for it in t1_rem:
        fOut.write(it+"\n")
        
    print(file = fOut)
    print(n2, fInName2, file = fOut)
    
    for it in t2_rem:
        fOut.write(it+"\n")

 
    print(file = fOut)
    print("Shared elements between two lists:", n, file = fOut)
    
    for it in t:
        fOut.write(it+"\n")   
        
    fOut.close()
    
    return [fOutName, n1_r, n2_r, n]



def GetDictFromFile(fInName, sep, header, unique_col = 0):
    fIn = open(fInName,"r")
    lines = fIn.readlines()
    
    header_line = ""
    
    if header == 1:
        header_line = lines[0]
        lines = lines[1:]
    
    
    d = {}
    t = []
    
    n = len(lines[0].split(sep))
    
    for line in lines:
        line = line.strip()
        if line != "":
            t_line = line.split(sep)
            if len(t_line) != n:
                print("Check format!", line, t_line)
                if t_line[0] not in d:
                    t.append(t_line[0])
                    d[t_line[0]] = t_line[1:]
                
            else:
                
                col_id = t_line[unique_col]
                tt = t_line[0:unique_col]+t_line[unique_col+1:]
                
                if col_id not in d:
                    t.append(col_id)
                    d[col_id] = tt
                else:
                    print("First colum has non-unique values!")
    fIn.close()
    return [t,d, header_line]





#Updated: May 4, 2019
# Needs update: base on the file extention (fasta or gb)    
def GetText(finName):
    '''Extracts text from a single fasta file'''
    fin = open(finName, 'r')
    text = ''
    for line in fin:
        line = line.strip()
        if (line != "") and (line[0] != '>'):
            text = text + line
    fin.close()
    return text


def CreateDictLocD_upper_count(text, m, gtp = "l"):
    if (m <= 0) or (m > len(text)):
        print("(n = "+str(m)+") is not a valid window size for this genome!!!");
        return {}
	
    d_g = dict()
    nn = len(text);
    gtype = gtp.lower();
    gtype = gtype[0];
	
    if gtype == "c":
        text = text + text[0:(m-1)];
        lastpos = nn;
    elif gtype == "l":
        lastpos = nn-m+1;
    else:
        print("Is this genome linear or circular?");
        return d_g;
		
    for ii in range (lastpos):
        bl = text[ii:(ii+m)]; 
        bl = bl.upper();
        if bl in d_g :
            d_g[bl] = d_g[bl] + 1
        else:
            d_g[bl] = 1
    return d_g;	



# Create a dictionary of unique strings of length mm with frequencies and locations for texttext
#Linear genome or circular genome!!!
def CreateDictLoc(text, mm, gtp = "l"):
	if (mm< 0) or (mm > len(text)):
		#print "N is bigger than genome size!!!";
		return {};
	
	
	d_g = dict()
	nn = len(text);
	gtype = gtp.lower();
	gtype = gtype[0];
	
	if gtype == "c":
		text = text + text[0:(mm-1)];
		lastpos = nn;
	elif gtype == "l":
			lastpos = nn-mm+1;
			#print lastpos;
	else:
		#print "Is this genome linear or circular?";
		return d_g;

		
	for ii in range (lastpos):
		bl = text[ii:(ii+mm)]; 
		#bl = bl.lower();
		#print ii, bl;
		if bl in d_g :
			tt = d_g[bl];
			#tt[0] = tt[0] + 1;
			tt.append(ii);
			d_g[bl] = tt;
		else:
			tt = list();
			#tt.append(1);
			tt.append(ii);
			d_g[bl] = tt;
			
	return d_g;	




def FindIntersection(d_g11, d_g22):
	
	t_ints = list()
	nn1 = len(d_g11);
	nn2 = len(d_g22);
	
	if nn1 <= nn2:
		d_first = d_g11;
		d_last = d_g22;
	else:
		d_first = d_g22;
		d_last = d_g11;
	
	for s in d_first:
		if s in d_last:
			t_ints.append(s);
	return t_ints;
