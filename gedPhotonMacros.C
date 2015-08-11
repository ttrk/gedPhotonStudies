/*
 * macro to study different photon Reconstruction algorithms
 * */

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TString.h>

#include <iostream>     // std::cout
#include <ctime>        // std::clock()
#include <algorithm>    // std::find()
#include <iomanip>      // std::setprecision()
#include <vector>

//#include "/net/hisrv0001/home/tatar/code/HIUtils/histoUtil.h"
#include "../../HeavyIonUtils/HIUtils/histoUtil.h"


TList* drawSame(TList* histos1, TList* histos2);
TList* draw2D(TList* histos);

void plotOnTop(const char* inputFileName, const char* dirName)
{
    TFile* inputFile = new TFile(inputFileName, "READ");
    std::cout << "input file for Histograms : " << inputFile->GetName() << std::endl;

    const char* suffix_dependence[2] = {"_cent", "_GENpT"};
    const char* suffix_pt[2] = {"", "_pt60"};
    const char* suffix_eta[4] = {"_barrel", "_endcap", "_endcap1", "_endcap2"};

    std::vector<std::string> histoNameList;

    // energyScale
    for (int i=0; i<2; ++i){
        for(int j=0; j<4; ++j){
            histoNameList.push_back(Form("energyScale%s%s", suffix_dependence[i], suffix_eta[j]));
            histoNameList.push_back(Form("widthEnergyScale%s%s", suffix_dependence[i], suffix_eta[j]));
        }
    }

    // energy scale distribution
    for(int j=0; j<4; ++j){
        histoNameList.push_back(Form("energyScale%s", suffix_eta[j]));
    }

    // position res.
    for(int j=0; j<4; ++j){
        histoNameList.push_back(Form("deltaPhi_cent%s", suffix_eta[j]));
        histoNameList.push_back(Form("widthDeltaPhi_cent%s", suffix_eta[j]));
        histoNameList.push_back(Form("deltaEta_cent%s", suffix_eta[j]));
        histoNameList.push_back(Form("widthDeltaEta_cent%s", suffix_eta[j]));
    }

    // deltaR
    const char* suffix_pt2[3] = {"", "_pt40", "_pt60"};
    for (int i=0; i<3; ++i){
        for(int j=0; j<4; ++j){
            histoNameList.push_back(Form("deltaR_cent%s%s", suffix_eta[j], suffix_pt2[i]));
            histoNameList.push_back(Form("widthDeltaR_cent%s%s", suffix_eta[j], suffix_pt2[i]));
        }
    }

    // isolation
    for (int i=0; i<2; ++i){
        for(int j=0; j<4; ++j){
            histoNameList.push_back(Form("trkIso%s%s", suffix_dependence[i], suffix_eta[j]));
            histoNameList.push_back(Form("widthTrkIso%s%s", suffix_dependence[i], suffix_eta[j]));
            histoNameList.push_back(Form("calIso%s%s", suffix_dependence[i], suffix_eta[j]));
            histoNameList.push_back(Form("widthCalIso%s%s", suffix_dependence[i], suffix_eta[j]));
        }
    }

    const char* suffix_cent[5] = {"", "_hiBin20", "_hiBin60", "_hiBin100", "_hiBin200"};
    const char* suffix_recoType[3] = {"", "_matched", "_fake"};

    // isolation (raw distributions)
    for (int i=0; i<5; ++i){
        for (int j=0; j<4; ++j){
            for (int k=0; k<3; ++k){
                histoNameList.push_back(Form("calIso_recoPhotons%s%s%s", suffix_recoType[k], suffix_eta[j], suffix_cent[i]));
                histoNameList.push_back(Form("trkIso_recoPhotons%s%s%s", suffix_recoType[k], suffix_eta[j], suffix_cent[i]));
            }
        }
    }

    // matching efficiency , fake ratio
    const char* suffix_dependence2[3] = {"_cent", "_GENpT", "_eta"};
    const char* suffix_dependence3[3] = {"_cent", "_RECOpT", "_eta"};
    for (int i=0; i<4; ++i){
        for (int j=0; j<2; ++j){
            for (int k=0; k<3; ++k){
                histoNameList.push_back(Form("matchRatio%s%s%s", suffix_dependence2[k], suffix_eta[i], suffix_pt[j]));
                histoNameList.push_back(Form("fakeRatio%s%s%s", suffix_dependence3[k], suffix_eta[i], suffix_pt[j]));
            }
        }
    }

    // shower shape
    for (int k=0; k<2; ++k){
        for (int j=0; j<4; ++j){
            for (int i=0; i<3; ++i){
                histoNameList.push_back(Form("phoR9_recoPhotons%s%s%s", suffix_recoType[i], suffix_eta[j], suffix_pt[k]));
                histoNameList.push_back(Form("phoHoverE_recoPhotons%s%s%s", suffix_recoType[i], suffix_eta[j], suffix_pt[k]));
                histoNameList.push_back(Form("phoSigmaIEtaIEta_recoPhotons%s%s%s", suffix_recoType[i], suffix_eta[j], suffix_pt[k]));
            }
        }
        histoNameList.push_back(Form("phoR9_recoPhotons%s", suffix_pt[k]));
        histoNameList.push_back(Form("phoHoverE_recoPhotons%s", suffix_pt[k]));
        histoNameList.push_back(Form("phoSigmaIEtaIEta_recoPhotons%s", suffix_pt[k]));
    }

    for (int j=0; j<4; ++j){
        for (int i=0; i<3; ++i){
            histoNameList.push_back(Form("phoEt_recoPhotons%s%s", suffix_recoType[i], suffix_eta[j]));
        }
    }

    histoNameList.push_back("phoEt_recoPhotons");

    const char* suffix_genType[3] = {"", "_matched", "_missing"};
    for (int j=0; j<4; ++j){
        for (int i=0; i<3; ++i){
            histoNameList.push_back(Form("genPhotons%s%s", suffix_genType[i], suffix_eta[j]));
        }
    }

    // 2D histograms
    // isolation correlation
    std::vector<std::string> histo2DNameList;
    for (int i=0; i<4; ++i){
        histo2DNameList.push_back(Form("corr_trkIso%s", suffix_eta[i]));
        histo2DNameList.push_back(Form("corr_calIso%s", suffix_eta[i]));
    }

    // add the GED version of all these histograms
    int len = histoNameList.size();
    for (int i = 0 ; i<len ; ++i){
        std::string histoName = histoNameList.at(i);
        histoNameList.push_back(Form("%s_GED",histoName.data()));
    }

    int len2D = histo2DNameList.size();
    for (int i = 0 ; i<len2D ; ++i){
        std::string histoName = histo2DNameList.at(i);
        histo2DNameList.push_back(Form("%s_GED",histoName.data()));
    }

    TList* histos = getListOfGIVENHistograms(inputFile, histoNameList);

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
    saveAllCanvasesToPicture(canvases, "gif", dirName);

    // TH2D
    TList* histos2D = getListOfGIVENHistograms(inputFile, histo2DNameList);
    TList* canvases2D = draw2D(histos2D);
    saveAllCanvasesToPicture(canvases2D, "pdf", dirName);
    saveAllCanvasesToPicture(canvases2D, "gif", dirName);
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
        c->Clear();

        h1->Draw("e");
        h2->Draw("e SAME");

        // set Y-axis ranges
        // default Y-axis range is that of h1
        double max1 = h1->GetMaximum();
        double max2 = h2->GetMaximum();
        double min1 = h1->GetMinimum();
        double min2 = h2->GetMinimum();

        if (max2 > max1)
            h1->SetMaximum(max2+TMath::Abs(max2)*0.2);

        if (min2 < min1)
            h1->SetMinimum(min2-TMath::Abs(min2)*0.2);

        // special cases
        TString histoName = h1->GetName();
        if(histoName.Contains("calIso_recoPhotons_fake") || histoName.Contains("trkIso_recoPhotons_matched")
                                                         || histoName.Contains("trkIso_recoPhotons_fake")
//                                                         || histoName.Contains("phoSigmaIEtaIEta_recoPhotons_matched")
                                                         || histoName.Contains("phoEt_recoPhotons")
                                                         || (histoName.Contains("genPhotons") && !histoName.Contains("missing"))
                                                         || histoName.Contains("phoHoverE"))
        {
            c->SetLogy();
        }
        if(histoName.Contains("_cent") || histoName.Contains("_GENpT")
                                       || histoName.Contains("_eta")
                                       || histoName.Contains("_RECOpT"))
        {
            h1->SetStats(0);
            h2->SetStats(0);
        }

        canvasList->Add(c);
    }

    return canvasList;
}

TList* draw2D(TList* histos)
{
    TList* canvasList = new TList();

    TH2D* h;
    TCanvas* c;
    TIter* iter = new TIter(histos);
    while((h=(TH2D*)iter->Next())) {

        c=new TCanvas(h->GetName());

        h->Draw("colz");
        c->SetLogz();

        canvasList->Add(c);
    }

    return canvasList;
}

int main(int argc, char** argv)
{
    plotOnTop("~/Desktop/gedPhotonResults/gedPhoton_AllQCDPhoton30_standard_forest_2nd_v2.root","~/Desktop/gedPhotonResults/gedPhoton_AllQCDPhoton30_standard_forest_2nd");
    plotOnTop("~/Desktop/gedPhotonResults/gedPhoton_AllQCDPhoton30_standard_forest_2nd_2012bins_v2.root","~/Desktop/gedPhotonResults/gedPhoton_AllQCDPhoton30_standard_forest_2nd_2012bins");

    plotOnTop("~/Desktop/gedPhotonResults/gedPhoton_EmEnrichedDijet30_standard_forest_v2.root","~/Desktop/gedPhotonResults/gedPhoton_EmEnrichedDijet30_standard_forest");
    plotOnTop("~/Desktop/gedPhotonResults/gedPhoton_EmEnrichedDijet30_standard_forest_2012bins_v2.root","~/Desktop/gedPhotonResults/gedPhoton_EmEnrichedDijet30_standard_forest_2012bins");

    return 0;
}
