#!/bin/bash

ls ../minimized/*.pdb > list

while read line
do
    filename=`echo ${line} | tr '/' ' ' | awk '{ print $3 }'`
    basename=`echo ${filename} | sed -e 's/.pdb//g'`

    echo $basename
    pdbconvert -ipdb ${line} -omae ../minimized/${basename}.mae
    canvasConvert -imae ../minimized/${basename}.mae -osd ../minimized/${basename}.sdf
done < list

rm list
rm ../minimized/*.mae