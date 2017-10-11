#!/bin/bash

# pass the name of the folder Something.mda and the name of the directory
#  that the folder lives in, e.g. NoiseRetino

directory=$2
temp=_MountainSort
saveDirectory=$directory$temp

folderName=$1

fullpath=/home/byron/CloudStation/ByronExp/$directory/$folderName

cd $fullpath

mp-daemon-start byron

var=0

for f in raw*.mda
   do
      var=$((var+1))
      datasetName="channel$var"
      prv-create $fullpath/$f /home/byron/CloudStation/ByronExp/$saveDirectory/datasets/$datasetName/raw.mda.prv
done

var=0
# kron-view results ms3 $datasetAbbrev
for f in event_times*.mda
   do
      var=$((var+1))
      cd /home/byron/CloudStation/ByronExp/$saveDirectory/
      datasetAbbrev="ch$var"
      kron-run ms3 --prescribed_event_times=$fullpath/$f $datasetAbbrev
     
  
      datasetName="channel$var"
      cd /home/byron/CloudStation/ByronExp/$saveDirectory/datasets/$datasetName
      rm raw.mda.prv
done

cd /tmp/mountainlab/mpdaemon/completed_processes
rm *.json

cd /tmp/mountainlab/tmp_long_term
rm *.tmp

cd /tmp/mountainlab/tmp_short_term
rm *

cd /home/byron/CloudStation/ByronExp/$saveDirectory/

temp=-mda.mat
fileName=${folderName:0:(-4)}$temp

touch fileName.txt
echo $fileName >> fileName.txt 

matlab -r -nodisplay "firingsmda2mat;exit"

rm fileName.txt

cd ..