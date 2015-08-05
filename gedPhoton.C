/*
 * macro to study different photon Reconstruction algorithms
 * */

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1D.h>
#include <TMath.h>

#include <iostream>     // std::cout
#include <ctime>        // std::clock()
#include <algorithm>    // std::find()
#include <iomanip>      // std::setprecision()
#include <vector>

Double_t getDR( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
Double_t getDPHI( Double_t phi1, Double_t phi2);
Double_t getDETA(Double_t eta1, Double_t eta2);
void gedPhoton(const char* hiForestfileName = "HiForest.root", const char* outputFileName = "gedPhoton.root");
std::vector<TH1*> gedPhotonAnalyzer(TFile* inputFile, const char* treePath);

const int MAXGENPARTICLES = 50000;  // number of gen particles can be large
const int MAXPHOTONS = 500;
const int PDG_PHOTON = 22;
const double cutdeltaR = 0.2; // 0.3    // cut for matching gen and reco. particles
const double cutptGEN  = 15;
const double cutptRECO = 15;

const int numHistos = 31;

void gedPhoton(const char* hiForestfileName, const char* outputFileName)
{
    TFile* inputFile = new TFile(hiForestfileName, "READ");
    std::cout << "input HiForest : " << inputFile->GetName() << std::endl;

    std::vector<TH1*> histos =gedPhotonAnalyzer(inputFile, "ggHiNtuplizer/EventTree");
    std::vector<TH1*> histosGED =gedPhotonAnalyzer(inputFile, "ggHiNtuplizerGED/EventTree");

    TFile* outputFile=new TFile(outputFileName, "RECREATE");
    outputFile->cd();

    // save histograms for old RECO
    for(int i=0; i<numHistos; ++i)
    {
        TH1D* h = (TH1D*)histos.at(i);
        h->Write();
    }
    // save histograms for GED RECO
    for(int i=0; i<numHistos; ++i)
    {
        TH1D* h_GED = (TH1D*)histosGED.at(i);
        h_GED->SetName(Form("%s_GED",h_GED->GetName()));
        h_GED->Write();
    }

    outputFile->Close();
    inputFile->Close();
}

std::vector<TH1*> gedPhotonAnalyzer(TFile* inputFile, const char* treePath)
{
    TTree* hiEvtAnalyzerTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
    TTree* ggHiNtuplizerTree = (TTree*)inputFile->Get(treePath);

    Int_t hiBin;
    hiEvtAnalyzerTree->SetBranchAddress("hiBin", &hiBin);

    // GEN particles
    Int_t nMC;
    std::vector<float>* mcPt=0;
    std::vector<float>* mcEta=0;
    std::vector<float>* mcPhi=0;
    std::vector<int>*   mcPID=0;
    std::vector<float>* mcCalIsoDR04=0;
    std::vector<float>* mcTrkIsoDR04=0;

    ggHiNtuplizerTree->SetBranchAddress("nMC",&nMC);
    ggHiNtuplizerTree->SetBranchAddress("mcPt",&mcPt);
    ggHiNtuplizerTree->SetBranchAddress("mcEta",&mcEta);
    ggHiNtuplizerTree->SetBranchAddress("mcPhi",&mcPhi);
    ggHiNtuplizerTree->SetBranchAddress("mcPID",&mcPID);
    ggHiNtuplizerTree->SetBranchAddress("mcCalIsoDR04",&mcCalIsoDR04);
    ggHiNtuplizerTree->SetBranchAddress("mcTrkIsoDR04",&mcTrkIsoDR04);

    // RECO photons
    Int_t nPho;
    std::vector<float>* phoEt=0;
    std::vector<float>* phoEta=0;
    std::vector<float>* phoPhi=0;
    std::vector<float>* pho_ecalClusterIsoR4=0;
    std::vector<float>* pho_hcalRechitIsoR4=0;
    std::vector<float>* pho_trackIsoR4PtCut20=0;

    ggHiNtuplizerTree->SetBranchAddress("nPho",&nPho);
    ggHiNtuplizerTree->SetBranchAddress("phoEt",&phoEt);
    ggHiNtuplizerTree->SetBranchAddress("phoEta",&phoEta);
    ggHiNtuplizerTree->SetBranchAddress("phoPhi",&phoPhi);
    ggHiNtuplizerTree->SetBranchAddress("pho_ecalClusterIsoR4",&pho_ecalClusterIsoR4);
    ggHiNtuplizerTree->SetBranchAddress("pho_hcalRechitIsoR4",&pho_hcalRechitIsoR4);
    ggHiNtuplizerTree->SetBranchAddress("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20);

    const int numGENptBins = 20;
    const double maxGENpt  = 200;
    double GENptBins_energyScale[numGENptBins];     // upper edges of energy scale histograms
    double energyScale[3][numGENptBins];            // energyScale[0][] = sum of energy scales
                                                    // energyScale[1][] = sum of square of energy scales
                                                    // energyScale[2][] = number of energy scales

    double GENptBins_Iso_ratio[numGENptBins];     // upper edges of isolation ratio histograms
    double trkIso_ratio_GENpT[3][numGENptBins];
    double calIso_ratio_GENpT[3][numGENptBins];

    for (int i=0; i<numGENptBins; ++i)
    {
        GENptBins_energyScale[i] = maxGENpt/numGENptBins*(i+1);
        GENptBins_Iso_ratio  [i] = maxGENpt/numGENptBins*(i+1);

        for (int j=0; j<3; ++j){
            energyScale[j][i]=0;
            trkIso_ratio_GENpT[j][i]=0;
            calIso_ratio_GENpT[j][i]=0;
        }
    }

    const int numCentBins = 22;
    const int maxCent     = 220;
    double CentBins_energyScale[numCentBins];     // upper edges of energy scale histograms
    double energyScale_cent[3][numCentBins];      // energyScale[0][] = sum of energy scales
                                                  // energyScale[1][] = sum of square of energy scales
                                                  // energyScale[2][] = number of energy scales

    double CentBins_pos_res[numCentBins];  // upper edges of position resolution histograms
    double deltaPhi_cent[3][numCentBins];
    double deltaEta_cent[3][numCentBins];
    double deltaR_cent[3][numCentBins];     // deltaR_cent[0][] = sum of abs(deltaR)
                                            // deltaR_cent[1][] = sum of square of abs(deltaR)
                                            // deltaR_cent[2][] = number of abs(deltaR)

    double CentBins_Iso_ratio[numCentBins];       // upper edges of isolation ratio histograms
    double trkIso_ratio_cent[3][numCentBins];
    double calIso_ratio_cent[3][numCentBins];           // calIso_ratio_cent[0][] = sum of isolation ratios
                                                        // calIso_ratio_cent[1][] = sum of square of isolation ratios
                                                        // calIso_ratio_cent[2][] = number of isolation ratios
    for (int i=0; i<numCentBins; ++i)
    {
        CentBins_energyScale[i] = (double)maxCent/numCentBins*(i+1);
        CentBins_pos_res    [i] = (double)maxCent/numCentBins*(i+1);
        CentBins_Iso_ratio  [i] = (double)maxCent/numCentBins*(i+1);

        // initialization
        for (int j=0; j<3; ++j){
            energyScale_cent [j][i]=0;

            deltaPhi_cent[j][i]=0;
            deltaEta_cent[j][i]=0;
            deltaR_cent  [j][i]=0;

            trkIso_ratio_cent[j][i]=0;
            calIso_ratio_cent[j][i]=0;
        }
    }

    // histograms for RECO photons
    TH1D* h[numHistos];
    h[0] = new TH1D("recoPhotons",         "RECO photons;p_{T} (GeV)",200,0,200);
    h[1] = new TH1D("recoPhotons_matched", "matched RECO photons (matched to GEN);p_{T} (GeV)",200,0,200);
    h[2] = new TH1D("recoPhotons_fake",    "fake RECO photons (not matched to GEN);p_{T} (GeV)",200,0,200);
    h[3] = new TH1D("genPhotons",          "GEN photons;p_{T} (GeV)",200,0,200);
    h[4] = new TH1D("genPhotons_matched",  "matched GEN photons (matched to RECO);p_{T} (GeV)",200,0,200);
    h[5] = new TH1D("genPhotons_missing",  "missing GEN photons (not matched to RECO);p_{T} (GeV)",200,0,200);
    h[6] = new TH1D("energyScale_GENpT",     "energy scale;GEN p_{T} (GeV);<RECO p_{T} / GEN p_{T}>",                numGENptBins,0,maxGENpt);
    h[7] = new TH1D("widthEnergyScale_GENpT","width of energy scale;GEN p_{T} (GeV);#sigma(RECO p_{T} / GEN p_{T})", numGENptBins,0,maxGENpt);
    h[8] = new TH1D("energyScale_cent",      "energy scale;hiBin;<RECO p_{T} / GEN p_{T}>",numCentBins,0,maxCent);
    h[9] = new TH1D("widthEnergyScale_cent", "width of energy scale;hiBin;#sigma(RECO p_{T} / GEN p_{T})",numCentBins,0,maxCent);
    h[10] = new TH1D("deltaPhi",       "#phi^{RECO} - #phi^{GEN};#Delta#phi",                    100,-cutdeltaR, cutdeltaR);
    h[11] = new TH1D("deltaEta",       "#eta^{RECO} - #eta^{GEN};#Delta#eta",                    100,-cutdeltaR, cutdeltaR);
    h[12] = new TH1D("deltaR",        "#DeltaR = #sqrt{#Delta#eta^{2}+#Delta#phi^{2}};#DeltaR",  100,0, cutdeltaR);
    h[13] = new TH1D("deltaPhi_cent",      "#Delta#phi = #phi^{RECO} - #phi^{GEN};hiBin;<|#Delta#phi|>",     numCentBins, 0, maxCent);
    h[14] = new TH1D("deltaEta_cent",      "#Delta#eta = #eta^{RECO} - #eta^{GEN};hiBin;<|#Delta#eta|>",     numCentBins, 0, maxCent);
    h[15] = new TH1D("deltaR_cent",        "#DeltaR = #sqrt{#Delta#eta^{2}+#Delta#phi^{2}};hiBin;<#DeltaR>", numCentBins, 0, maxCent);
    h[16] = new TH1D("widthDeltaPhi_cent", "#sigma(|#Delta#phi|);hiBin;#sigma(|#Delta#phi|)",   numCentBins, 0, maxCent);
    h[17] = new TH1D("widthDeltaEta_cent", "#sigma(|#Delta#eta|);hiBin;#sigma(|#Delta#eta|)",   numCentBins, 0, maxCent);
    h[18] = new TH1D("widthDeltaR_cent",   "#sigma(#DeltaR)   ;hiBin;#sigma(#DeltaR)",          numCentBins, 0, maxCent);
    h[19] = new TH1D("trkIso",            "#DeltaE_{track}^{ISO} = E_{track}^{ISO}(RECO) - E_{track}^{ISO}(GEN);#DeltaE_{track}^{ISO}", 200,-100, 100);
    h[20] = new TH1D("calIso",            "#DeltaE_{calo}^{ISO} = E_{calo}^{ISO}(RECO) - E_{calo}^{ISO}(GEN);#DeltaE_{calo}^{ISO}",     200,-100, 100);
    h[21] = new TH1D("trkIso_ratio",      "E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN);E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN)", 100, -5, 5);
    h[22] = new TH1D("calIso_ratio",      "E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN);E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN)",     100, -5, 5);
    h[23] = new TH1D("trkIso_ratio_GENpT",      ";GEN p_{T} (GeV);<E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN)>", numGENptBins,0,maxGENpt);
    h[24] = new TH1D("calIso_ratio_GENpT",      ";GEN p_{T} (GeV);<E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN)>",   numGENptBins,0,maxGENpt);
    h[25] = new TH1D("widthTrkIso_ratio_GENpT", ";GEN p_{T} (GeV);#sigma(E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN))", numGENptBins,0,maxGENpt);
    h[26] = new TH1D("widthCalIso_ratio_GENpT", ";GEN p_{T} (GeV);#sigma(E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN))",   numGENptBins,0,maxGENpt);
    h[27] = new TH1D("trkIso_ratio_cent",       ";hiBin;<E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN)>", numCentBins,0,maxCent);
    h[28] = new TH1D("calIso_ratio_cent",       ";hiBin;<E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN)>",   numCentBins,0,maxCent);
    h[29] = new TH1D("widthTrkIso_ratio_cent",  ";hiBin;#sigma(E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN))", numCentBins,0,maxCent);
    h[30] = new TH1D("widthCalIso_ratio_cent",  ";hiBin;#sigma(E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN))",   numCentBins,0,maxCent);

    std::cout << "entering event loop" << std::endl;
    Long64_t entries = ggHiNtuplizerTree->GetEntries();
    std::cout << "number of entries = " << entries << std::endl;

    std::clock_t    start_loop, end_loop;
    start_loop = std::clock();
    for(Long64_t jj = 0; jj < entries; ++jj)
    {
        if (jj % 100000 == 0)  {
          std::cout << "current entry = " <<jj<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)jj/entries*100<<" %"<<std::endl;
        }

        hiEvtAnalyzerTree->GetEntry(jj);
        ggHiNtuplizerTree->GetEntry(jj);

        bool isFake[nPho];
        int  matched[nPho];           // index of the matched GEN photon in  to this RECO photon
        double deltaPhi[nPho];        // deltaR between the RECO photon and the matched GEN photon
        double deltaEta[nPho];        // deltaR between the RECO photon and the matched GEN photon
        double deltaR[nPho];          // deltaR between the RECO photon and the matched GEN photon
        double deltaTrkIso[nPho];
        double deltaCalIso[nPho];
        double ratioTrkIso[nPho];
        double ratioCalIso[nPho];

        // find GEN photons that match to RECO photons
        for (int i=0; i < nPho; ++i)
        {
            bool passedRECOselection = (phoEt->at(i) > cutptRECO) ;
            if(!passedRECOselection) continue;

            // default values
            isFake[i]  = true;
            matched[i] = -1;
            deltaPhi[i]  = 999;
            deltaEta[i]  = 999;
            deltaR[i]    = 999;
            deltaTrkIso[i]    = 999;
            deltaCalIso[i]    = 999;
            ratioTrkIso[i]    = 999;
            ratioCalIso[i]    = 999;
            double deltaRMin  = 999;

            bool passedGENselection;      // selections for GEN photon
            bool passedDR;
            double deltaRtmp;
            for (int j=0; j<nMC; ++j)
            {
                deltaRtmp = getDR(phoEta->at(i), phoPhi->at(i), mcEta->at(j), mcPhi->at(j));
                passedGENselection = (mcPID->at(j) == PDG_PHOTON) && (mcPt->at(j) > cutptGEN);
                passedDR           = (deltaRtmp < cutdeltaR);

                if (passedGENselection && passedDR)
                {
                    // matched GEN photon is the one closest to the RECO photon
                    if (deltaRtmp < deltaRMin)
                    {
                        deltaRMin  = deltaRtmp;
                        matched[i] = j;
                    }
                }
            }
            isFake[i] = (matched[i] == -1);     // if no matched GEN photon, then this RECO photon is fake.

            h[0]->Fill(phoEt->at(i));         // all RECO photons
            if (isFake[i] == false) {       // RECO photons matched to GEN photon

                int j = matched[i];         // index of matched GEN photon
                deltaPhi[i]  = getDPHI(phoPhi->at(i),mcPhi->at(j));
                deltaEta[i]  = phoEta->at(i) - mcEta->at(j);
                deltaR  [i]  = deltaRMin;
                deltaTrkIso[i] = pho_trackIsoR4PtCut20->at(i) - mcTrkIsoDR04->at(j);
                deltaCalIso[i] = (pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)) - mcCalIsoDR04->at(j);
                ratioTrkIso[i] = pho_trackIsoR4PtCut20->at(i) / mcTrkIsoDR04->at(j);
                ratioCalIso[i] = (pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)) / mcCalIsoDR04->at(j);

                // energy scale as fnc. of GEN pT
                for(int r=0; r < numGENptBins; ++r){
                    if(mcPt->at(j)<=GENptBins_energyScale[r])
                    {
                        double val = phoEt->at(i)/mcPt->at(j);
                        energyScale[0][r] += val;
                        energyScale[1][r] += val*val;
                        energyScale[2][r] += 1;
                        break;
                    }
                }

                // energy scale as fnc. of centrality
                for(int r=0; r < numCentBins; ++r){
                    if(hiBin <= CentBins_energyScale[r])
                    {
                        double val = phoEt->at(i)/mcPt->at(j);
                        energyScale_cent[0][r] += val;
                        energyScale_cent[1][r] += val*val;
                        energyScale_cent[2][r] += 1;
                        break;
                    }
                }

                // position resolution as fnc. of centrality
                for(int r=0; r < numCentBins; ++r){
                    if(hiBin <= CentBins_pos_res[r])
                    {
                        double val = TMath::Abs(deltaPhi[i]);
                        deltaPhi_cent[0][r] += val;
                        deltaPhi_cent[1][r] += val*val;
                        deltaPhi_cent[2][r] += 1;

                        val = TMath::Abs(deltaEta[i]);
                        deltaEta_cent[0][r] += val;
                        deltaEta_cent[1][r] += val*val;
                        deltaEta_cent[2][r] += 1;

                        deltaR_cent  [0][r] += deltaR[i];
                        deltaR_cent  [1][r] += deltaR[i]*deltaR[i];
                        deltaR_cent  [2][r] += 1;
                        break;
                    }
                }

                // isolation ratios as fnc. of GEN pT
                for(int r=0; r < numGENptBins; ++r){
                    if(mcPt->at(j)<=GENptBins_Iso_ratio[r])
                    {
                        trkIso_ratio_GENpT[0][r] += ratioTrkIso[i];
                        trkIso_ratio_GENpT[1][r] += ratioTrkIso[i]*ratioTrkIso[i];
                        trkIso_ratio_GENpT[2][r] += 1;

                        calIso_ratio_GENpT[0][r] += ratioCalIso[i];
                        calIso_ratio_GENpT[1][r] += ratioCalIso[i]*ratioCalIso[i];
                        calIso_ratio_GENpT[2][r] += 1;
                        break;
                    }
                }

                // isolation ratios as fnc. of centrality
                for(int r=0; r < numCentBins; ++r){
                    if(hiBin <= CentBins_Iso_ratio[r])
                    {
                        trkIso_ratio_cent[0][r] += ratioTrkIso[i];
                        trkIso_ratio_cent[1][r] += ratioTrkIso[i]*ratioTrkIso[i];
                        trkIso_ratio_cent[2][r] += 1;

                        calIso_ratio_cent[0][r] += ratioCalIso[i];
                        calIso_ratio_cent[1][r] += ratioCalIso[i]*ratioCalIso[i];
                        calIso_ratio_cent[2][r] += 1;
                        break;
                    }
                }

                h[1] ->Fill(phoEt->at(i));
                h[10]->Fill(deltaPhi[i]);
                h[11]->Fill(deltaEta[i]);
                h[12]->Fill(deltaR[i]);
                h[19]->Fill(deltaTrkIso[i]);
                h[20]->Fill(deltaCalIso[i]);
                h[21]->Fill(ratioTrkIso[i]);
                h[22]->Fill(ratioCalIso[i]);
            }
            else {
                // fake RECO photons
                h[2]->Fill(phoEt->at(i));
            }
        }

        // find RECO photons that match to GEN photons
        bool isMiss[nMC];
        for(int i = 0; i < nMC; ++i){

            // default values
            isMiss[i]=true;

            bool passedGENselection;      // selections for GEN photon
            passedGENselection = (mcPID->at(i) == PDG_PHOTON) && (mcPt->at(i) > cutptGEN);

            if(passedGENselection)
            {
                // check if that GEN photon was matched to a RECO photon in the previous loop over RECO photons
                int* indexRECO = std::find(matched,matched+nPho, i);
                bool found = (indexRECO < matched+nPho);
                isMiss[i] = !found;

                // fill histograms
                h[3]->Fill(mcPt->at(i));
                if (isMiss[i] == false)
                    h[4]->Fill(mcPt->at(i));
                else
                    h[5]->Fill(mcPt->at(i));
            }
        }
    }

    // fill remaining histograms
    // histograms as fnc. of GEN pT
    for(int r=0; r < numGENptBins; ++r){
        int n;
        double mean;
        double meanOfSquares;

        // energy scale as fnc. of GEN pT
        n = energyScale[2][r];
        if(n>0)
        {
            mean          = energyScale[0][r]/n;
            meanOfSquares = energyScale[1][r]/n;
            h[6]->SetBinContent(r+1,mean);
            h[7]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
        }

        // isolation ratios as fnc. of GEN pT
        n = trkIso_ratio_GENpT[2][r];
        if (n>0) {
            mean          = trkIso_ratio_GENpT[0][r]/n;
            meanOfSquares = trkIso_ratio_GENpT[1][r]/n;
            h[23]->SetBinContent(r+1,mean);
            h[25]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
        }

        n = calIso_ratio_GENpT[2][r];
        if (n>0) {
            mean          = calIso_ratio_GENpT[0][r]/n;
            meanOfSquares = calIso_ratio_GENpT[1][r]/n;
            h[24]->SetBinContent(r+1,mean);
            h[26]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
        }
    }

    // histograms as fnc. of centrality
    for(int r=0; r < numCentBins; ++r){
        int n;
        double mean;
        double meanOfSquares;

        // energy scale as fnc. of centrality
        n = energyScale_cent[2][r];
        if (n>0) {
            mean          = energyScale_cent[0][r]/n;
            meanOfSquares = energyScale_cent[1][r]/n;
            h[8]->SetBinContent(r+1,mean);
            h[9]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
        }

        // position resolution as fnc. of centrality
        n = deltaPhi_cent[2][r];
        if (n>0) {
            mean          = deltaPhi_cent[0][r]/n;
            meanOfSquares = deltaPhi_cent[1][r]/n;
            h[13]->SetBinContent(r+1,mean);
            h[16]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
        }

        n = deltaEta_cent[2][r];
        if (n>0) {
            mean          = deltaEta_cent[0][r]/n;
            meanOfSquares = deltaEta_cent[1][r]/n;
            h[14]->SetBinContent(r+1,mean);
            h[17]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
        }

        n = deltaR_cent[2][r];
        if (n>0) {
            mean          = deltaR_cent[0][r]/n;
            meanOfSquares = deltaR_cent[1][r]/n;
            h[15]->SetBinContent(r+1,mean);
            h[18]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
        }
        // isolation ratios as fnc. of centrality
        n = trkIso_ratio_cent[2][r];
        if (n>0) {
            mean          = trkIso_ratio_cent[0][r]/n;
            meanOfSquares = trkIso_ratio_cent[1][r]/n;
            h[27]->SetBinContent(r+1,mean);
            h[29]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
        }

        n = calIso_ratio_cent[2][r];
        if (n>0) {
            mean          = calIso_ratio_cent[0][r]/n;
            meanOfSquares = calIso_ratio_cent[1][r]/n;
            h[28]->SetBinContent(r+1,mean);
            h[30]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
        }
    }

    end_loop = std::clock();
    std::cout.precision(6);      // get back to default precision
    std::cout << "LOOP finished in             : " << (end_loop - start_loop) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "exited event loop" << std::endl;

    std::vector<TH1*> histos;
    for (int i=0; i<numHistos; ++i)
    {
        histos.push_back(h[i]);
    }
    return histos;
}

int main(int argc, char** argv)
{
    if(argc == 1)
    {
        gedPhoton();
        return 0;
    }
    else if(argc == 2)
    {
        gedPhoton(argv[1]);
        return 0;
    }
    else if(argc == 3)
    {
        gedPhoton(argv[1], argv[2]);
        return 0;
    }
    else
    {
        std::cout<<"wrong input"<<std::endl;
    }
}

/*
 * copied from https://github.com/CmsHI/HiForestAnalysis/blob/master/commonUtility.h
 */
Double_t getDR( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2){
  Double_t theDphi = getDPHI( phi1, phi2);
  Double_t theDeta = eta1 - eta2;
  return TMath::Sqrt ( theDphi*theDphi + theDeta*theDeta);
}

/*
 * copied from https://github.com/CmsHI/HiForestAnalysis/blob/master/commonUtility.h
 */
Double_t getDPHI( Double_t phi1, Double_t phi2) {
  Double_t dphi = phi1 - phi2;

  if ( dphi > 3.141592653589 )
    dphi = dphi - 2. * 3.141592653589;
  if ( dphi <= -3.141592653589 )
    dphi = dphi + 2. * 3.141592653589;

  if ( TMath::Abs(dphi) > 3.141592653589 ) {
    std::cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 " << std::endl;
  }

  return dphi;
}

Double_t getDETA(Double_t eta1, Double_t eta2){
    return eta1 - eta2;
}
