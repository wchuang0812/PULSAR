#!/bin/bash

#Observation ID
obsID='1145367872'
#ncpus
ncpus=1
#Your Presto route
prestor='$PRESTO'
#Data file formate (.fits, .fil)
filetype='.fits'
#Choose to use normal search (zmax=0) or accelerate search (zmax=?int)
zmax=200
#Length of time used in rfifind (float in seconds)
rfitime=12.0
#Would you like to remove the 0-dispersion signal? (birds=1, Yes; birds=0, No)
birds=0
#Number of subbands used in de-dispersion
numsubbands=512
#Zap channels (zc=0 if no gaps among channels; zc=1, sign zap channels)
zc=1
#Fold all candidates (fcand=1) or not (fcand=0)
fcand=1
#Use single pulse search (singlep=1) or not (singlep=0)
singlep=1

source presto.sh
export OMP_NUM_THREADS=${ncpus}
cd /ssd/Pulsar
if [ ! -d "$obsID" ]; then
    mkdir $obsID
fi
cd $obsID
rm -rf *


echo ${prestor} > parameters.txt
echo ${obsID} >> parameters.txt
echo ${filetype} >> parameters.txt
echo ${zmax} >> parameters.txt
echo ${rfitime} >> parameters.txt
echo ${birds} >> parameters.txt
echo ${numsubbands} >> parameters.txt
echo ${zc} >> parameters.txt
echo ${fcand} >> parameters.txt
echo ${singlep} >> parameters.txt
echo ${ncpus} >> parameters.txt
echo "/ssd/Pulsar/${obsID}/fits/" >> parameters.txt

cp -r /ssd/Pulsar/Pulsar_script/autosearch_pulsar3.2opt.py .
cp -r /ssd/Pulsar/Pulsar_script/single_pulse_search_opt.py .
cp -r /ssd/Pulsar/Pulsar_script/dedispersion.py .
cp -r /ssd/Pulsar/Pulsar_script/accelsearch.py .
cp -r /ssd/Pulsar/Pulsar_script/prepfold.py .
cp -r /o9000/MWA/Pulsar/ICS_data/${obsID}/fits .

mkdir command.txt
python autosearch_pulsar3.2opt.py 
