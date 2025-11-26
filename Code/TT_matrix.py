"""
Copyright: (c) 2025 Tatiana Lenskaia. All rights reserved.
E-mail: t.lenskaia@icloud.com
"""

import core_methods as cm

path = "/Users/tatianalenskaia/_PAT/TT/TT_compare/"
fInName = "PAT_TThits_36HMMs.txt"

fIn = open(path+fInName, "r")
lines = fIn.readlines()
fIn.close()



fInName = "PAT_TT_HMMs.txt"
t_hmm = cm.GetListFromFile(path+fInName)
n_hmm = len(t_hmm)

print(n_hmm)



lines = lines[1:]

thr = 0.0001

#d_hmm = {}


M = [[0 for i in range(n_hmm)] for j in range(n_hmm)]

print(M)

for line in lines:
    t_line = line.strip().split("\t")
    #print(t_line)
    g_id = t_line[3]
    pid = t_line[4]
    hmm = t_line[5]
    e = float(t_line[8])
    pos_st = int(t_line[-2])
    pos_fn = int(t_line[-1])


    



    
    
    '''
    
    if hmm not in d_hmm:
        d_hmm[hmm] = 0
    d_hmm[hmm] += 1
    '''
'''
for hmm in t_hmm:
    print(hmm, d_hmm[hmm])
'''
