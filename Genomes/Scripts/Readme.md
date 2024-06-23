# Scripts

This folder contains a set of custom scripts for obtaining and processing genomic sequences.
Below, one can find brief instructions for getting the PAT genomic sequences.

## 1. Download

**PAT_ids.txt** contains a list of NCBI IDs for the PAT genomic sequences.
There are multiple ways of obtaining sequences from NCBI.
Here we describe, one of the possible solutions using our custom script with Biopython.
We assume that Python 3 and Biopython are already installed on the local machine where genomes are to be downloaded.

Save and **PATgb_fetch.py** in the same directory where PAT_ids.txt file is saved on the local machine.
Let this directory be the current working directory by opening a terminal (or a command line window) and navigating to this directory.
Next, run the following command:

```
python3 PATgb_fetch.py
```

The genomic sequences in GenBank format will be uploaded to the current working directory. 
Please make sure that there is enough space on the local machine to accomodate the sequences (**about 18 Gb**).  
