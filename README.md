# RNA structure information from 


Extract json/csv file
------
[Internal loop motifs](http://rna.bgsu.edu/rna3dhub/motifs/release/il/3.57)
- InternalLoopMotifAtlasRelease3.57.json

[Hairpin loop motifs](http://rna.bgsu.edu/rna3dhub/motifs/release/hl/3.57)
- HairpinLoopMotifAtlasRelease3.57.json

[Representative RNA structures](http://rna.bgsu.edu/rna3dhub/nrlist/release/3.230)
- nrlist_3.230_4.0A.csv



Extract coordinate data (cif) from json/csv data
------

- Internal/Hairpin loop motifs  
    
    1. Extract cif data using representative RNA stuctures
        >wget http://rna.bgsu.edu/rna3dhub/rest/getCoordinates?coord=[motif id]
    2. Modify and split cif file into subsequent base sets (i.e. 2 connected bases)
    3. Categorize connected bases into groups (UG, CA, etc.)  
    4. Calculate internal distance and perform clustering for each group
    5. Select cluster centroid for each group

- Interaction pairs (Specify set of residues) 

    1. Extract cif data for base-pairs, base-stacking, base-phosphate, base-ribose interacions using 
        >wget http://rna.bgsu.edu/rna3dhub/pdb/2QBG/interactions/fr3d/basepairs/csv  
        >wget http://rna.bgsu.edu/rna3dhub/rest/getCoordinates?coord=2QBG|1|A|G|69,2QBG|1|A|G|107
    2. Modify and split cif file into base sets and interaction type (UC-s55, GU-s53)
    4. Calculate internal distance and perform clustering for each group
    5. Select cluster centroid for each group