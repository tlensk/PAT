"""
Copyright: (c) 2025 Tatiana Lenskaia. All rights reserved.
E-mail: t.lenskaia@icloud.com
"""

import core_methods as cm


def HMMinfo(path,fname):
    res = cm.GetDictFromFile(path+fname, "\t", 1, unique_col = 1)
    
    # d_hmm contains hmm description: key = hmm_id, value = array of two elements:hmm_nm and hmm_desc.
    # For example: 'PT001001': ['ST-Lam', 'Small Terminase subunit; TerS']
    d_hmm = res[1]
    
    #print(len(d_hmm))
    #print(d_hmm)
    
    
    
    
    # d_gr contains hmm ids grouped by category id: key = cat_id, value = array of hmm_ids in this category.
    # For example, 'PT033': ['PT033001', 'PT033002', 'PT033003', 'PT033004', 'PT033005']
    d_gr = {}
    t = sorted(list(d_hmm.keys()))
    
    for hmm in t:
        cat = hmm[0:5]
        #print(cat)
        if cat not in d_gr:
            d_gr[cat] = []
        d_gr[cat].append(hmm)
    
    
    #print(d_gr)
    
    
    # d_cat contains info about categories: key = cat_id, value = dictionary of cat_nm with counts in a given category.
    # For example, 'PT011': {'HT1a': 19, 'HT1b': 3}
    d_cat = {}
    for gr in d_gr:
        #print(gr, len(d_gr[gr]))
        for hmm in d_gr[gr]:
            cat_id = hmm[0:5]
            hmm_nm = d_hmm[hmm][0]
            cat_nm = hmm_nm.split("-")[0]
            #print(hmm_nm)
            if cat_id not in d_cat:
                d_cat[cat_id] = {}
            dd = d_cat[cat_id]
            if cat_nm not in dd:
                dd[cat_nm] = 1
            else:
                dd[cat_nm] += 1
                
            d_cat[cat_id] = dd
            
    #print(d_cat)
      
    
    # d_catnm2id is for cat_nm to cat_id conversion: key = cat_nm, value=cat_id.
    # For example, 'ST': 'PT001'
    d_catnm2id = {}
    for cat in d_cat:
        #print(cat, d_cat[cat])
        
        for item in d_cat[cat]:
            #print(item)
            if item not in d_catnm2id:
                d_catnm2id[item] = cat
    
    #print(d_catnm2id)
    return [d_hmm, d_gr, d_cat, d_catnm2id]

"""

res = HMMinfo(path, fname)
#print(len(res))
'''
res[0]
    # d_hmm contains hmm description: key = hmm_id, value = array of two elements:hmm_nm and hmm_desc.
    # For example: 'PT001001': ['ST-Lam', 'Small Terminase subunit; TerS']

res[1]
    # d_gr contains hmm ids grouped by category id: key = cat_id, value = array of hmm_ids in this category.
    # For example, 'PT033': ['PT033001', 'PT033002', 'PT033003', 'PT033004', 'PT033005']

res[2]
    # d_cat contains info about categories: key = cat_id, value = dictionary of cat_nm with counts in a given category.
    # For example, 'PT011': {'HT1a': 19, 'HT1b': 3}

res[3]
    # d_catnm2id is for cat_nm to cat_id conversion: key = cat_nm, value=cat_id.
    # For example, 'ST': 'PT001'
'''
for it in res:
    print(len(it))
"""
