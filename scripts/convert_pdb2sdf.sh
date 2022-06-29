#!/bin/bash


cd ../minimized

ls *.pdb > files

while read line
do
    basename=`echo $line | sed -e 's/.pdb//g'`
    echo "$basename"

    pdbconvert -ipdb ${line} -omae ${basename}.mae
    canvasConvert -imae ${basename}.mae -osd ${basename}.sdf
    #mol2convert -imae ${basename}.mae -omol2 ${basename}.mol2
done < files

rm files
