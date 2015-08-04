#!/bin/sh

progName="gedPhoton";
#inputFile="/afs/cern.ch/user/k/katatar/CMSSW_7_5_0/src/HiForest.root";
inputFile="~/CMSSW_7_5_0/src/HiForest.root";

g++ $progName.C $(root-config --cflags --libs) -Wall -O2 -o $progName.exe || exit 1

./$progName.exe $inputFile || exit 1
