#!/usr/bin/env python
# coding: utf-8

# # Download [RNA base triple database](http://rna.bgsu.edu/triples/triples.php)

# In[1]:


import os, sys, shutil
import pathlib
import glob as glob
import numpy as np
import json
import wget
from itertools import product
from zipfile import ZipFile
from pdbfixer import PDBFixer
import warnings
from openbabel import openbabel
import time


# In[2]:


#warnings.filterwarnings("ignore")
#sys.stderr = sys.__stderr__


# In[3]:


url = "http://rna.bgsu.edu/triples/zip"
release_version = "v1.4"


# In[4]:


base_path = os.path.dirname(os.path.abspath("__file__")).strip('notebooks')
output_path = os.path.join(base_path, "pdb", "triplebase")


# In[5]:


if os.path.isdir(output_path):
    print(">remove directory: {}".format(output_path))
    shutil.rmtree(output_path)
    
pathlib.Path(output_path).mkdir(parents=True, exist_ok=True) 


# In[6]:


d = "GCAU"
arr = list(product(d, repeat=3))


# In[7]:


for a in arr:
    seq = ''.join(a)    
    _output_path = os.path.join(output_path, seq)
    _url = os.path.join(url, release_version, seq + ".zip")
    
    print('{}.zip'.format(seq))
    wget.download(_url, out=output_path, bar=None)     
    shutil.unpack_archive('{}.zip'.format(_output_path), _output_path)


# In[8]:


# delete zip file
for a in arr:
    seq = ''.join(a)
    os.remove(os.path.join(output_path, seq + ".zip"))


# ### load pdb with openbabel and resave as pdb

# In[9]:


files = glob.glob(output_path + "/*/*.pdb")


# In[10]:


files


# In[11]:


for f in files:
    print(os.path.basename(f))
    
    f_org = f + ".org"
    shutil.move(f, f_org)

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "pdb")
    mol = openbabel.OBMol()    
    print("n_atoms: {}, n_bonds: {}, n_residues: {}".format(mol.NumAtoms(), mol.NumBonds(), mol.NumResidues()))

    obConversion.ReadFile(mol, f_org)
    #mol.AddHydrogens()
    mol.DeleteHydrogens()    
    obConversion.WriteFile(mol, f)

    time.sleep(1)


# In[12]:


for f in files:
    basename = os.path.basename(f)

    try:
        # check converted pdb with PDBFixer
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("default")
            fixer = PDBFixer(filename=f)    

        # raise warning if duplicate residue exists
        if len(w) != 0:
            print("{}: {}".format(basename, w[0]))  
    except:
        print("{}: Invalid".format(basename))


# In[ ]:





# In[ ]:




