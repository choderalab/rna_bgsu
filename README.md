# RNA structure information from [RNA BGSU](https://www.bgsu.edu/research/rna.html)


Download json/csv files
------
[Internal loop motifs](http://rna.bgsu.edu/rna3dhub/motifs/release/il/3.57)
- InternalLoopMotifAtlasRelease3.57.json

[Hairpin loop motifs](http://rna.bgsu.edu/rna3dhub/motifs/release/hl/3.57)
- HairpinLoopMotifAtlasRelease3.57.json

[Representative RNA structures](http://rna.bgsu.edu/rna3dhub/nrlist/release/3.233)
- nrlist_3.233_2.5A.csv



Extract structure data using [RNA BGSU APIs](https://www.bgsu.edu/research/rna/APIs.html)
------

- `download_motifs.ipynb`  
    - Download internal/hairpin loop motifs
    - Note that these are not properly formatted in CIF format where html/css codes are also included in the file.
    
- `download_junctions_from_nrlist.ipynb`  
    - Download junction loops 
    - Since junction loops are not provided as part of the [RNA 3D Motif Atlas](http://rna.bgsu.edu/rna3dhub/motifs) at the time of analysis, junction loops were extracted from [Representative RNA structures (version 3.233, resolution<2.5)](http://rna.bgsu.edu/rna3dhub/nrlist/release/3.233).

- `download_basepair.ipynb`  
    - Download base pairs from [BPCatalog](http://ndbserver.rutgers.edu/ndbmodule/services/BPCatalog/bpCatalog.html)  

- `download_triplebase.ipynb`  
    - Download base triples from [RNA Base Triple Database](http://rna.bgsu.edu/triples/triples.php)
    - Note that files are downloaded as pdb and most structures are modeled. Experimental structures are tagged with "exemplar".
    - Inconsistent residue numbers were roughly tested and fixed when possible



Check downloaded cif/pdb files
------

- `check_motif.ipynb`  
    - Check for missing atoms, inconsistent residue numbers, residue names. Following files failed to pass the test. Filename added with "warning" or "missing_atoms" and ignored for further process. 
        - HL_44730.2.cif, HL_51090.1.cif, HL_01181.4.cif, HL_70505.1.cif, HL_35188.1.cif, HL_49873.1.cif, HL_19239.1.cif, HL_48810.1.cif, IL_46464.1.cif, IL_51971.1.cif, IL_77341.1.cif, IL_39900.1.cif, IL_16160.1.cif, IL_50112.1.cif, IL_57188.3.cif, IL_89836.1.cif, IL_22427.1.cif, IL_80617.1.cif, IL_64414.1.cif, J3_7RQB_030.cif, J3_7RQ8_031.cif, J3_7RQ8_003.cif, J3_7RQB_003.cif

- `check_base.ipynb`  
    - Check for inconsistency. Same residue numbers were found for BP_cWH_GG.cif and BP_cWH_UU.cif. Residue numbers fixed manually.



Split loops into three consecutive bases and perform clustering for each base set
------

- `split_motif.ipynb`  
    - Hairpins, internal, and junction loops were split into three connected bases (triple bases).
    - Atoms `P`, `OP1`, and `OP2` were deleted from 5' base.

- `cluster_motifs.ipynb`
    - Each triple base set was clustered with bottom-up (agglomerative) hierarchical clustering with average linking.
    - Euclidean distances between set of pre-defined atom pairs were used as input features.
    - Distance threhold of "1" was used for clustering. This was defined so that the number of cluster centroids from each triple base set were ~4000 structures. The maximum RMSD measured from the centroid structure for each cluster was also monitored to insure structures are well-clustered. Maximum RMSD from all clusters are below ~2.2 Angstroms.



Edit 5' base residues and add TER from base pair catalogs and base triple database structures
------

- `edit_base.ipynb`  
    - Atoms `P`, `OP1`, and `OP2` were deleted from 5' base.
    - TER was added to insure bases were disconnected



Add hydrogen atoms and minimize hydrogens
------

- `script/check_files.py`
    - All files were checked if they could be loaded with PDBFixer
    - Warning found for two BP files. TER was inserted in the wrong position. Manually fixed.

- `script/openmm_implicit_minimizer.py`
    - Added hydrogen prior to minimization
    - Implicit solvent applied (GBN2)
    - Heavy atom restraint applied (5 kcal/mol/A^2)
    - 10 step NVT to check nothing fuzzy with the input structure
        -  Errors raised for some Triple base structures downloaded from [RNA Base Triple Database](http://rna.bgsu.edu/triples/triples.php)
    - It turns out that some triple base structures (e.g. Triple_tWH_cSS_CAG.pdb/cif) have overlapping atoms... These errors were raised for modeled structures. No problems detected for experimental structures (i.e. "exemplar").   
    - Triple base structures were moved to directory `minimized/_obsolete` 

- `script/convert_pdb2sdf.sh`
    - convert minimized pdb to sdf format using schrodinger software (pdbconvert: pdb -> mae / canvasConvert: mae -> sdf)
