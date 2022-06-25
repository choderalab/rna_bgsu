#!/usr/bin/env python
# coding: utf-8

# ## Minimize structures
# 
# Add missing hydrogen using pdbfixer. 5' base will be automatically detected so there is no need to change the first residue name.
# 
# **Notes**
# Alternatively, missing hydrogens could be added by pdb4amber. Note that the first residue needs to be modify (e.g. A --> A5).
# >import pdb4amber  
# >pdb4amber.run(arg_pdbin='input.pdb', arg_pdbout='pdb4amber.pdb', arg_add_missing_atoms=True)
# 
# **References**
# - [OpenMM-Tricks-and-Recipes](https://github-wiki-see.page/m/ParmEd/ParmEd/wiki/OpenMM-Tricks-and-Recipes)
# - [OPENMM_TUTORIAL](https://gpantel.github.io/assets/PDF/OpenMM_Tutorial.pdf)  

import os, sys, shutil
import pathlib
import glob as glob
import numpy as np
import re
import warnings
import mdtraj as md
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from openmm.app import PDBFile
from pdbfixer import PDBFixer




def test(f):
    fixer = PDBFixer(filename=f)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.findMissingAtoms()
    
    if fixer.missingAtoms:
        print("{}: missing atoms".format(f))
        shutil.move(f, f + ".warning")
        status = "failed"
    else:
        status = "success"

    return status




if __name__ == "__main__":
    base_path = os.path.dirname(os.path.abspath("__file__")).strip('scripts')
    output_path = os.path.join(base_path, "minimized")

    # motif
    _path = os.path.join(base_path, "pdb", "motif", "cluster", "triplebase")
    motif_files = glob.glob(_path + "/*/centroid/rep*.pdb")
    
    # triplebase
    _path = os.path.join(base_path, "pdb", "triplebase")
    triplebase_files = glob.glob(_path + "/*.pdb" )
        
    # basepairs
    _path = os.path.join(base_path, "pdb", "bpcatalog")
    bpcatalog_files = glob.glob(_path + "/*.pdb" )

    
    files = motif_files + triplebase_files + bpcatalog_files
    #files = triplebase_files + bpcatalog_files
    #files = motif_files
    print(">{} files found".format(len(files)))
    
    
    # check
    for f in files:
        status = test(f)
    
