#!/bin/sh

progName="gedPhoton";
inputFile="/afs/cern.ch/user/k/katatar/CMSSW_7_5_0/src/HiForest.root";
# outputFile="test_GammaJetAnalyzer_v1507.root";

g++ $progName.C $(root-config --cflags --libs) -Wall -O2 -o $progName.exe || exit 1

./$progName.exe $inputFile || exit 1
