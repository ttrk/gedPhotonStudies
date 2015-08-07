/*
 * macro to study different photon Reconstruction algorithms
 * */

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TList.h>
#include <TString.h>

#include <iostream>     // std::cout
#include <ctime>        // std::clock()
#include <algorithm>    // std::find()
#include <iomanip>      // std::setprecision()
#include <vector>

#include "/net/hisrv0001/home/tatar/code/HIUtils/histoUtil.h"
//#include "../../HeavyIonUtils/HIUtils/histoUtil.h"


TList* drawSame(TList* histos1, TList* histos2);

void plotOnTop(const char* inputFileName, const char* dirName)
{
    TFile* inputFile = new TFile(inputFileName, "READ");
    std::cout << "input file for Histograms : " << inputFile->GetName() << std::endl;

    TList* histos = getListOfALLHistograms(inputFile);


    TList* histos_old = new TList();    // histograms for old RECO
    TList* histos_ged = new TList();    // histograms for old RECO
    TH1D* h;
    TString h_name;
    TIter* iter = new TIter(histos);
    while ((h=(TH1D*)iter->Next())) {

        h_name=h->GetName();
        if(!h_name.EndsWith("_GED"))
        {
            histos_old->Add(h);
        }
        else {
            histos_ged->Add(h);
        }
    }

    TList* canvases = drawSame(histos_old, histos_ged);

    saveAllCanvasesToPicture(canvases, "pdf", dirName);
}

TList* drawSame(TList* histos1, TList* histos2)
{
    TList* canvasList = new TList();

    TH1D* h1;
    TH1D* h2;
    TCanvas* c;
    TIter* iter1 = new TIter(histos1);
    TIter* iter2 = new TIter(histos2);
    while((h1=(TH1D*)iter1->Next())) {
        h2=(TH1D*)iter2->Next();

        h1->SetMarkerStyle(21);
        h2->SetMarkerStyle(20);

        h1->SetMarkerColor(kBlack);
        h2->SetMarkerColor(kRed);

        h2->SetMarkerSize(h1->GetMarkerSize()*0.8);

        c=new TCanvas(h1->GetName());

        h1->Draw("e");
        h2->Draw("e SAME");

        canvasList->Add(c);
    }

    return canvasList;
}

int main(int argc, char** argv)
{
    plotOnTop("~/Desktop/gedPhoton_merged_AllQCDPhoton30_standard_forest_2nd.root","~/Desktop/AllQCDPhoton30_standard_forest_2nd");
//    plotOnTop("~/Desktop/gedPhoton_merged_EmEnrichedDijet30_standard_forest.root","~/Desktop/EmEnrichedDijet30_standard_forest");
    //plotOnTop("~/Desktop/gedPhoton_MBHydjet_502_photonCommish_partialMerge.root","~/Desktop/MBHydjet_502_photonCommish_partialMerge");

    return 0;
}
