#!/bin/sh

progName="gedPhoton_2012bins";
#inputFile="/afs/cern.ch/user/k/katatar/CMSSW_7_5_0/src/HiForest.root";
#inputFile="~/CMSSW_7_5_0/src/HiForest.root";
inputFile="/mnt/hadoop/cms/store/user/luck/GEDPhotonCommissioning/MBHydjet_502_photonCommish_partialMerge/0.root";
outputFile="gedPhoton_MBHydjet_502_photonCommish_partialMerge_2012bins.root";

inputFile1="/mnt/hadoop/cms/store/user/tatar/GEDPhotonCommissioning/merged_AllQCDPhoton30_standard_forest_2nd.root";
inputFile2="/mnt/hadoop/cms/store/user/tatar/GEDPhotonCommissioning/merged_EmEnrichedDijet30_standard_forest.root";

outputFile1="gedPhoton_AllQCDPhoton30_standard_forest_2nd_2012bins_v2.root";
outputFile2="gedPhoton_EmEnrichedDijet30_standard_forest_2012bins_v2.root";

g++ $progName.C $(root-config --cflags --libs) -Wall -O2 -o $progName.exe || exit 1

#./$progName.exe $inputFile $outputFile || exit 1
./$progName.exe $inputFile1 $outputFile1 || exit 1
./$progName.exe $inputFile2 $outputFile2 || exit 1
