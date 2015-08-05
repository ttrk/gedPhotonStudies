#!/bin/sh

progName="gedPhoton";
#inputFile="/afs/cern.ch/user/k/katatar/CMSSW_7_5_0/src/HiForest.root";
#inputFile="~/CMSSW_7_5_0/src/HiForest.root";
inputFile="/mnt/hadoop/cms/store/user/luck/GEDPhotonCommissioning/MBHydjet_502_photonCommish_partialMerge/0.root";

outputFile="gedPhoton_MBHydjet_502_photonCommish_partialMerge.root";

g++ $progName.C $(root-config --cflags --libs) -Wall -O2 -o $progName.exe || exit 1

./$progName.exe $inputFile $outputFile || exit 1
