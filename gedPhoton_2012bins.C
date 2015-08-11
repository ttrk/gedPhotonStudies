/*
 * macro to study different photon Reconstruction algorithms
 * version of gedPhoton.C to match pT binning convention in http://arxiv.org/pdf/1201.3093v2.pdf
 * */

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include <iostream>     // std::cout
#include <ctime>        // std::clock()
#include <algorithm>    // std::find()
#include <iomanip>      // std::setprecision()
#include <vector>

const int MAXGENPARTICLES = 50000;  // number of gen particles can be large
const int MAXPHOTONS = 500;
const int PDG_PHOTON = 22;
const int cutmcStatus = 1;
const double cutdeltaR = 0.2; // 0.3    // cut for matching gen and reco. particles
const float cutptGEN  = 20;
const float cutptRECO = 0;				// switch off the cut for RECO photon pT, only cut on the GEN photon pt,
										// no RECO pT cut so that you keep the downward-biased RECO photons

const float cutetaBarrel = 1.4791;			// cut to separate photons into Barrel and Endcap photons
const float cutetaEndCap = 2;				// cut to separate photons in Endcap into 2.

const int cutmcMomPID_pi0 = 111;

const int numHistos = 70;
const int numHistos2D = 2;

Double_t getDR( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
Double_t getDPHI( Double_t phi1, Double_t phi2);
Double_t getDETA(Double_t eta1, Double_t eta2);
void gedPhoton(const char* hiForestfileName = "HiForest.root", const char* outputFileName = "gedPhoton.root");
std::vector<TH1*> gedPhotonAnalyzer(TFile* inputFile, const char* treePath,
                                        float eta_gt = 0, float eta_lt = cutetaBarrel,
                                        int mcMomPID_gt = -999999, int mcMomPID_lt = 999999,
                                        int hiBin_gt = -999999, int hiBin_lt = 999999,
                                        float ptGEN = cutptGEN);

void gedPhoton(const char* hiForestfileName, const char* outputFileName)
{
    TFile* inputFile = new TFile(hiForestfileName, "READ");
    std::cout << "input HiForest : " << inputFile->GetName() << std::endl;


    float eta_gt[5] = {     0,            0, cutetaBarrel, cutetaBarrel, cutetaEndCap};
    float eta_lt[5] = {999999, cutetaBarrel,       999999, cutetaEndCap,       999999};
    int mcMomPID_gt[5] = {-999999, 21, -999999,     22, 110};
    int mcMomPID_lt[5] = { 999999, 23,      22, 999999, 112};
    int hiBin_gt[4] = {-999999,  0, 20,     100};   // int hiBin_gt[5] = {-999999,  0, 20,  60,     100};
    int hiBin_lt[4] = { 999999, 20, 60, 9999999};   // int hiBin_lt[5] = { 999999, 20, 60, 100, 9999999};
    float ptGEN[5] = {20, 15, 30, 40, 60};

    std::vector<TH1*> histos[2][5][5][5][5];
    // [2][][][][] : OLD Reco, GED Reco
    // [][5][][][] : eta cuts : no cut, Barrel, Endcap, Endcap1, Endcap2
    // [][][5][][] : mcMomPID cuts : all photons, prompt photon, fragmentation photon, decay photon (bkg.), decay photon from pi0
    // [][][][5][] : hiBin cuts : all centralities, 0-10%, 10-30%, 30-50%, 50-100%
    // [][][][][5] : GEN pT cuts

    // get histograms
    for (int j=0; j<5; ++j)
    {
        for (int k=0; k<5; ++k)
        {
            for (int l=0; l<4; ++l)
            {
                for (int m=0; m<5; ++m)
                {
                    // do not run it for all combinations, othwerwise will get 10K many histograms
                    bool skipAnalyzer = (k>0 && l>0 && m>0);
                    if(skipAnalyzer) continue;

                    histos[0][j][k][l][m] = gedPhotonAnalyzer(inputFile, "ggHiNtuplizer/EventTree",
                            eta_gt[j],      eta_lt[j],
                            mcMomPID_gt[k], mcMomPID_lt[k],
                            hiBin_gt[l], hiBin_lt[l],
                            ptGEN[m]);

                    histos[1][j][k][l][m] = gedPhotonAnalyzer(inputFile, "ggHiNtuplizerGED/EventTree",
                            eta_gt[j],      eta_lt[j],
                            mcMomPID_gt[k], mcMomPID_lt[k],
                            hiBin_gt[l], hiBin_lt[l],
                            ptGEN[m]);
                }
            }
        }
    }

//    std::vector<TH1*> histosOLD = gedPhotonAnalyzer(inputFile, "ggHiNtuplizer/EventTree", 0, 999999);
//    std::vector<TH1*> histosOLD_Barrel = gedPhotonAnalyzer(inputFile, "ggHiNtuplizer/EventTree", 0, cutetaBarrel);
//    std::vector<TH1*> histosOLD_Endcap = gedPhotonAnalyzer(inputFile, "ggHiNtuplizer/EventTree", cutetaBarrel, 999999);
//    std::vector<TH1*> histosOLD_Endcap1 = gedPhotonAnalyzer(inputFile, "ggHiNtuplizer/EventTree", cutetaBarrel, cutetaEndCap);
//    std::vector<TH1*> histosOLD_Endcap2 = gedPhotonAnalyzer(inputFile, "ggHiNtuplizer/EventTree", cutetaEndCap, 999999);
//
//    std::vector<TH1*> histosGED = gedPhotonAnalyzer(inputFile, "ggHiNtuplizerGED/EventTree", 0, 9999);
//    std::vector<TH1*> histosGED_Barrel = gedPhotonAnalyzer(inputFile, "ggHiNtuplizerGED/EventTree", 0, cutetaBarrel);
//    std::vector<TH1*> histosGED_Endcap = gedPhotonAnalyzer(inputFile, "ggHiNtuplizerGED/EventTree", cutetaBarrel, 9999);
//    std::vector<TH1*> histosGED_Endcap1 = gedPhotonAnalyzer(inputFile, "ggHiNtuplizerGED/EventTree", cutetaBarrel, cutetaEndCap);
//    std::vector<TH1*> histosGED_Endcap2 = gedPhotonAnalyzer(inputFile, "ggHiNtuplizerGED/EventTree", cutetaEndCap, 9999);

    TFile* outputFile=new TFile(outputFileName, "RECREATE");
    outputFile->cd();

    const char* eta_hname[5] = {"","_barrel", "_endcap", "_endcap1", "_endcap2"};
    const char* eta_title[5] = {"",", |#eta| < 1.4791 (Barrel)", ", |#eta| >= 1.4791 (Endcap)", ", |#eta| < 2 (Endcap)", ", |#eta| >= 2 (Endcap)"};

    const char* mcMomPID_hname[5] = {"","_prompt", "_frag", "_decay", "_decaypi0"};
    const char* mcMomPID_title[5] = {"",", #gamma^{prompt} (mcMomPID==22)", ", #gamma^{frag} (abs(mcMomPID)<22)", ", #gamma^{decay} (abs(mcMomPID)>22)", ", #gamma^{decay}_{#pi^{0}} (mcMomPID==111)"};

    const char* hiBin_hname[4] = {"","_hiBin20", "_hiBin60", "_hiBin200"};//{"","_hiBin20", "_hiBin60", "_hiBin100", "_hiBin200"};
    const char* hiBin_title[4] = {"",", hiBin:0-20", ", hiBin:20-60", ", hiBin:60-200"};//{"",", hiBin:0-20", ", hiBin:20-60", ", hiBin:60-100", ", hiBin:100-200"};

    const char* ptGEN_hname[5] = {"","_pt15", "_pt30", "_pt40", "_pt60"};
    const char* ptGEN_title[5] = {"",", p_{T}^{#gamma}(GEN)>15", ", p_{T}^{#gamma}(GEN)>30", ", p_{T}^{#gamma}(GEN)>40", ", p_{T}^{#gamma}(GEN)>60"};
    // rename histograms and save them
    TH1* h;
    TH1* hGED;
    for (int j=0; j<5; ++j)
    {
        for (int k=0; k<5; ++k)
        {
            for (int l=0; l<4; ++l){

                for (int m=0; m<5; ++m)
                {
                    for(int i=0; i<numHistos+numHistos2D; ++i){
                        // do not run it for all combinations, othwerwise will get 10K many histograms
                        bool skipAnalyzer = (k>0 && l>0 && m>0);
                        if(skipAnalyzer) continue;

                        h    = (TH1*)histos[0][j][k][l][m].at(i);
                        hGED = (TH1*)histos[1][j][k][l][m].at(i);

                        h->SetName(Form("%s%s%s%s%s", h->GetName(),
                                eta_hname[j],
                                mcMomPID_hname[k],
                                hiBin_hname[l],
                                ptGEN_hname[m]));
                        hGED->SetName(Form("%s%s%s%s%s_GED", hGED->GetName(),
                                eta_hname[j],
                                mcMomPID_hname[k],
                                hiBin_hname[l],
                                ptGEN_hname[m]));

                        h->SetTitle(Form("%s%s%s%s%s", h->GetTitle(),
                                eta_title[j],
                                mcMomPID_title[k],
                                hiBin_title[l],
                                ptGEN_title[m]));

                        hGED->SetTitle(Form("%s%s%s%s%s", hGED->GetTitle(),
                                eta_title[j],
                                mcMomPID_title[k],
                                hiBin_title[l],
                                ptGEN_title[m]));

                        h->Write();
                        hGED->Write();
                    }
                }
            }
        }
    }

    outputFile->Close();
    inputFile->Close();
}

std::vector<TH1*> gedPhotonAnalyzer(TFile* inputFile, const char* treePath,
                                        float eta_gt, float eta_lt,
                                        int mcMomPID_gt, int mcMomPID_lt,
                                        int hiBin_gt, int hiBin_lt,
                                        float ptGEN)
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
    std::vector<int>*   mcMomPID=0;
    std::vector<int>*   mcStatus=0;
    std::vector<float>* mcCalIsoDR04=0;
    std::vector<float>* mcTrkIsoDR04=0;

    ggHiNtuplizerTree->SetBranchAddress("nMC",&nMC);
    ggHiNtuplizerTree->SetBranchAddress("nMC",&nMC);
    ggHiNtuplizerTree->SetBranchAddress("mcPt",&mcPt);
    ggHiNtuplizerTree->SetBranchAddress("mcEta",&mcEta);
    ggHiNtuplizerTree->SetBranchAddress("mcPhi",&mcPhi);
    ggHiNtuplizerTree->SetBranchAddress("mcPID",&mcPID);
    ggHiNtuplizerTree->SetBranchAddress("mcMomPID",&mcMomPID);
    ggHiNtuplizerTree->SetBranchAddress("mcStatus",&mcStatus);
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
    std::vector<float>* phoR9=0;
    std::vector<float>* phoHoverE=0;
    std::vector<float>* phoSigmaIEtaIEta=0;
//    std::vector<float>* phoSigmaIPhiIPhi=0;

    ggHiNtuplizerTree->SetBranchAddress("nPho",&nPho);
    ggHiNtuplizerTree->SetBranchAddress("phoEt",&phoEt);
    ggHiNtuplizerTree->SetBranchAddress("phoEta",&phoEta);
    ggHiNtuplizerTree->SetBranchAddress("phoPhi",&phoPhi);
    ggHiNtuplizerTree->SetBranchAddress("pho_ecalClusterIsoR4",&pho_ecalClusterIsoR4);
    ggHiNtuplizerTree->SetBranchAddress("pho_hcalRechitIsoR4",&pho_hcalRechitIsoR4);
    ggHiNtuplizerTree->SetBranchAddress("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20);
    ggHiNtuplizerTree->SetBranchAddress("phoR9",&phoR9);
    ggHiNtuplizerTree->SetBranchAddress("phoHoverE",&phoHoverE);
    ggHiNtuplizerTree->SetBranchAddress("phoSigmaIEtaIEta",&phoSigmaIEtaIEta);
//    ggHiNtuplizerTree->SetBranchAddress("phoSigmaIPhiIPhi_2012",&phoSigmaIPhiIPhi);

    // GEN pT dependent variables
    const int numGENptBins = 5;//    const int numGENptBins = 20;
    const double maxGENpt  = 80;//    const double maxGENpt  = 200;
    double GENptBins_energyScale[numGENptBins];     // upper edges of energy scale histograms
    double energyScale[3][numGENptBins];            // energyScale[0][] = sum of energy scales
                                                    // energyScale[1][] = sum of square of energy scales
                                                    // energyScale[2][] = number of energy scales

    double GENptBins_Iso_ratio[numGENptBins];     // upper edges of isolation ratio histograms
    double trkIso_ratio_GENpT[3][numGENptBins];
    double calIso_ratio_GENpT[3][numGENptBins];

    double GENptBins_Iso[numGENptBins];     // upper edges of isolation ratio histograms
    double trkIso_GENpT[3][numGENptBins];
    double calIso_GENpT[3][numGENptBins];

    double GENptBins_matchRatio[numGENptBins];     // upper edges of matching efficiency histograms
    double matchRatio_GENpT[3][numGENptBins];     // CHANGE of convention for array indices :
                                                     // matchRatio_GENpT[0][] = # of matched GEN photons
                                                     // matchRatio_GENpT[1][] = empty
                                                     // matchRatio_GENpT[2][] = # of all     GEN photons

    double RECOptBins_fakeRatio[numGENptBins];     // upper edges of fake rate histograms
    double fakeRatio_RECOpT[3][numGENptBins];      // CHANGE of convention for array indices :
                                                     // fakeRatio_GENpT[0][] = # of fake RECO photons
                                                     // fakeRatio_GENpT[1][] = empty
                                                     // fakeRatio_GENpT[2][] = # of all  RECO photons

    double pTBins_12013093v2[5]={25, 30, 40, 50, 80};
    for (int i=0; i<numGENptBins; ++i)
    {
        GENptBins_energyScale[i] = pTBins_12013093v2[i]; //maxGENpt/numGENptBins*(i+1);
        GENptBins_Iso_ratio  [i] = pTBins_12013093v2[i]; //maxGENpt/numGENptBins*(i+1);
        GENptBins_Iso        [i] = pTBins_12013093v2[i]; //maxGENpt/numGENptBins*(i+1);
        GENptBins_matchRatio [i] = pTBins_12013093v2[i]; //maxGENpt/numGENptBins*(i+1);
        RECOptBins_fakeRatio  [i] = pTBins_12013093v2[i]; //maxGENpt/numGENptBins*(i+1);

        for (int j=0; j<3; ++j){
            energyScale[j][i]=0;
            trkIso_ratio_GENpT[j][i]=0;
            calIso_ratio_GENpT[j][i]=0;
            trkIso_GENpT      [j][i]=0;
            calIso_GENpT      [j][i]=0;
            matchRatio_GENpT  [j][i]=0;
            fakeRatio_RECOpT   [j][i]=0;
        }
    }

    // centrality dependent variables
    const int numCentBins = 4; //11;
    const int maxCent     = 200; //220;
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

    double CentBins_Iso[numCentBins];       // upper edges of isolation difference histograms
    double trkIso_cent[3][numCentBins];
    double calIso_cent[3][numCentBins];           // calIso_cent[0][] = sum of isolation differences
                                                  // calIso_cent[1][] = sum of square of isolation differences
                                                  // calIso_cent[2][] = number of isolation differences

    double CentBins_matchRatio[numCentBins];     // upper edges of matching efficiency histograms
    double matchRatio_cent[3][numCentBins];

    double CentBins_fakeRatio[numCentBins];     // upper edges of fake rate histograms
    double fakeRatio_cent[3][numCentBins];

    double CentBins_12013093v2[4] = {20, 60, 100, 200};
    for (int i=0; i<numCentBins; ++i)
    {
        CentBins_energyScale[i] = CentBins_12013093v2[i];  //(double)maxCent/numCentBins*(i+1);
        CentBins_pos_res    [i] = CentBins_12013093v2[i];  //(double)maxCent/numCentBins*(i+1);
        CentBins_Iso_ratio  [i] = CentBins_12013093v2[i];  //(double)maxCent/numCentBins*(i+1);
        CentBins_Iso        [i] = CentBins_12013093v2[i];  //(double)maxCent/numCentBins*(i+1);
        CentBins_matchRatio [i] = CentBins_12013093v2[i];  //(double)maxCent/numCentBins*(i+1);
        CentBins_fakeRatio  [i] = CentBins_12013093v2[i];  //(double)maxCent/numCentBins*(i+1);

        // initialization
        for (int j=0; j<3; ++j){
            energyScale_cent [j][i]=0;

            deltaPhi_cent[j][i]=0;
            deltaEta_cent[j][i]=0;
            deltaR_cent  [j][i]=0;

            trkIso_ratio_cent[j][i]=0;
            calIso_ratio_cent[j][i]=0;

            trkIso_cent[j][i]=0;
            calIso_cent[j][i]=0;

            matchRatio_cent[j][i]=0;
            fakeRatio_cent [j][i]=0;
        }
    }

    // eta dependent variables
    const int numEtaBins = 20;
    double tmpMaxEta = eta_lt;
    if (eta_lt > 5)
        tmpMaxEta = 5;
    const double maxEta  = tmpMaxEta;
    double etaBins_matchRatio[numEtaBins];     // upper edges of matching efficiency histograms
    double matchRatio_eta[3][numEtaBins];

    double etaBins_fakeRatio[numEtaBins];     // upper edges of fake rate histograms
    double fakeRatio_eta[3][numEtaBins];

    for (int i=0; i<numEtaBins; ++i)
    {
        etaBins_matchRatio[i] = -maxEta + maxEta/numEtaBins*(i+1)*2;
        etaBins_fakeRatio [i] = -maxEta + maxEta/numEtaBins*(i+1)*2;

        for (int j=0; j<3; ++j){
            matchRatio_eta[j][i]=0;
            fakeRatio_eta [j][i]=0;
        }
    }


    const double GENptBins_arr[6] = {20, 25, 30, 40, 50, 80};
    const double CentBins_arr[5]  = {0, 20, 60, 100, 200};
    TH1::SetDefaultSumw2();
    // histograms for photons
    TH1D* h[numHistos];
    h[0] = new TH1D("phoEt_recoPhotons",         "RECO photons;p_{T} (GeV)",100,0,200);
    h[1] = new TH1D("phoEt_recoPhotons_matched", "matched RECO photons (matched to GEN);p_{T} (GeV)",100,0,200);
    h[2] = new TH1D("phoEt_recoPhotons_fake",    "fake RECO photons (not matched to GEN);p_{T} (GeV)",100,0,200);
    h[3] = new TH1D("genPhotons",          "GEN photons;p_{T} (GeV)",100,0,200);
    h[4] = new TH1D("genPhotons_matched",  "matched GEN photons (matched to RECO);p_{T} (GeV)",100,0,200);
    h[5] = new TH1D("genPhotons_missing",  "missing GEN photons (not matched to RECO);p_{T} (GeV)",100,0,200);
    h[6] = new TH1D("energyScale_GENpT",     "energy scale;GEN p_{T} (GeV);<RECO p_{T} / GEN p_{T}>",                numGENptBins, GENptBins_arr);
    h[7] = new TH1D("widthEnergyScale_GENpT","width of energy scale;GEN p_{T} (GeV);#sigma(RECO p_{T} / GEN p_{T})", numGENptBins, GENptBins_arr);
    h[8] = new TH1D("energyScale_cent",      "energy scale;hiBin;<RECO p_{T} / GEN p_{T}>",numCentBins, CentBins_arr);
    h[9] = new TH1D("widthEnergyScale_cent", "width of energy scale;hiBin;#sigma(RECO p_{T} / GEN p_{T})",numCentBins, CentBins_arr);
    h[10] = new TH1D("deltaPhi",       "#phi^{RECO} - #phi^{GEN};#Delta#phi",                    100,-cutdeltaR/2, cutdeltaR/2);
    h[11] = new TH1D("deltaEta",       "#eta^{RECO} - #eta^{GEN};#Delta#eta",                    100,-cutdeltaR/2, cutdeltaR/2);
    h[12] = new TH1D("deltaR",        "#DeltaR = #sqrt{#Delta#eta^{2}+#Delta#phi^{2}};#DeltaR",  100,0, cutdeltaR/2);
    h[13] = new TH1D("deltaPhi_cent",      "#Delta#phi = #phi^{RECO} - #phi^{GEN};hiBin;<|#Delta#phi|>",     numCentBins, CentBins_arr);
    h[14] = new TH1D("deltaEta_cent",      "#Delta#eta = #eta^{RECO} - #eta^{GEN};hiBin;<|#Delta#eta|>",     numCentBins, CentBins_arr);
    h[15] = new TH1D("deltaR_cent",        "#DeltaR = #sqrt{#Delta#eta^{2}+#Delta#phi^{2}};hiBin;<#DeltaR>", numCentBins, CentBins_arr);
    h[16] = new TH1D("widthDeltaPhi_cent", "#sigma(|#Delta#phi|);hiBin;#sigma(|#Delta#phi|)",   numCentBins, CentBins_arr);
    h[17] = new TH1D("widthDeltaEta_cent", "#sigma(|#Delta#eta|);hiBin;#sigma(|#Delta#eta|)",   numCentBins, CentBins_arr);
    h[18] = new TH1D("widthDeltaR_cent",   "#sigma(#DeltaR)   ;hiBin;#sigma(#DeltaR)",          numCentBins, CentBins_arr);
    h[19] = new TH1D("trkIso",            "#DeltaE_{track}^{ISO} = E_{track}^{ISO}(RECO) - E_{track}^{ISO}(GEN);#DeltaE_{track}^{ISO}", 100,-100, 100);
    h[20] = new TH1D("calIso",            "#DeltaE_{calo}^{ISO} = E_{calo}^{ISO}(RECO) - E_{calo}^{ISO}(GEN);#DeltaE_{calo}^{ISO}",     100,-100, 100);
    h[21] = new TH1D("trkIso_ratio",      "E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN);E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN)", 100, -5, 5);
    h[22] = new TH1D("calIso_ratio",      "E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN);E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN)",     100, -5, 5);
    h[23] = new TH1D("trkIso_ratio_GENpT",      "<E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN)>;GEN p_{T} (GeV);<E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN)>", numGENptBins, GENptBins_arr);
    h[24] = new TH1D("calIso_ratio_GENpT",      "<E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN)>;GEN p_{T} (GeV);<E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN)>",     numGENptBins, GENptBins_arr);
    h[25] = new TH1D("widthTrkIso_ratio_GENpT", "#sigma(E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN));GEN p_{T} (GeV);#sigma(E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN))", numGENptBins, GENptBins_arr);
    h[26] = new TH1D("widthCalIso_ratio_GENpT", "#sigma(E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN));GEN p_{T} (GeV);#sigma(E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN))",     numGENptBins, GENptBins_arr);
    h[27] = new TH1D("trkIso_ratio_cent",       "<E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN)>;hiBin;<E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN)>", numCentBins, CentBins_arr);
    h[28] = new TH1D("calIso_ratio_cent",       "<E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN)>;hiBin;<E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN)>",     numCentBins, CentBins_arr);
    h[29] = new TH1D("widthTrkIso_ratio_cent",  "#sigma(E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN));hiBin;#sigma(E_{track}^{ISO}(RECO) / E_{track}^{ISO}(GEN))", numCentBins, CentBins_arr);
    h[30] = new TH1D("widthCalIso_ratio_cent",  "#sigma(E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN));hiBin;#sigma(E_{calo}^{ISO}(RECO) / E_{calo}^{ISO}(GEN))",     numCentBins, CentBins_arr);

    h[31] = new TH1D("trkIso_GENpT",      "<#DeltaE_{track}^{ISO} = E_{track}^{ISO}(RECO) - E_{track}^{ISO}(GEN)>;GEN p_{T} (GeV);<#DeltaE_{track}^{ISO}>", numGENptBins, GENptBins_arr);
    h[32] = new TH1D("calIso_GENpT",      "<#DeltaE_{calo}^{ISO} = E_{calo}^{ISO}(RECO) - E_{calo}^{ISO}(GEN)>;GEN p_{T} (GeV);<#DeltaE_{calo}^{ISO}>",     numGENptBins, GENptBins_arr);
    h[33] = new TH1D("widthTrkIso_GENpT", "#sigma(#DeltaE_{track}^{ISO} = E_{track}^{ISO}(RECO) - E_{track}^{ISO}(GEN));GEN p_{T} (GeV);#sigma(#DeltaE_{track}^{ISO})", numGENptBins, GENptBins_arr);
    h[34] = new TH1D("widthCalIso_GENpT", "#sigma(#DeltaE_{calo}^{ISO} = E_{calo}^{ISO}(RECO) - E_{calo}^{ISO}(GEN));GEN p_{T} (GeV);#sigma(#DeltaE_{calo}^{ISO})",     numGENptBins, GENptBins_arr);
    h[35] = new TH1D("trkIso_cent",       "<#DeltaE_{track}^{ISO} = E_{track}^{ISO}(RECO) - E_{track}^{ISO}(GEN)>;hiBin;<#DeltaE_{track}^{ISO}>", numCentBins,0,maxCent);
    h[36] = new TH1D("calIso_cent",       "<#DeltaE_{calo}^{ISO} = E_{calo}^{ISO}(RECO) - E_{calo}^{ISO}(GEN)>;hiBin;<#DeltaE_{calo}^{ISO}>",     numCentBins,0,maxCent);
    h[37] = new TH1D("widthTrkIso_cent",  "#sigma(#DeltaE_{track}^{ISO} = E_{track}^{ISO}(RECO) - E_{track}^{ISO}(GEN));hiBin;#sigma(#DeltaE_{track}^{ISO})", numCentBins,0,maxCent);
    h[38] = new TH1D("widthCalIso_cent",  "#sigma(#DeltaE_{calo}^{ISO} = E_{calo}^{ISO}(RECO) - E_{calo}^{ISO}(GEN));hiBin;#sigma(#DeltaE_{calo}^{ISO})",     numCentBins,0,maxCent);

    h[39] = new TH1D("matchRatio_GENpT",  "matching efficiency = N_{#gamma}^{GEN}(matched) / N_{#gamma}^{GEN}(all);GEN p_{T} (GeV);N_{#gamma}^{GEN}(matched) / N_{#gamma}^{GEN}(all)", numGENptBins, GENptBins_arr);
    h[40] = new TH1D("fakeRatio_RECOpT",   "fake rate = N_{#gamma}^{RECO}(fake) / N_{#gamma}^{RECO}(all);RECO p_{T} (GeV);N_{#gamma}^{RECO}(fake) / N_{#gamma}^{RECO}(all)",             numGENptBins, GENptBins_arr);
    h[41] = new TH1D("matchRatio_cent",   "matching efficiency = N_{#gamma}^{GEN}(matched) / N_{#gamma}^{GEN}(all);hiBin;N_{#gamma}^{GEN}(matched) / N_{#gamma}^{GEN}(all)",  numCentBins, CentBins_arr);
    h[42] = new TH1D("fakeRatio_cent",    "fake rate = N_{#gamma}^{RECO}(fake) / N_{#gamma}^{RECO}(all);hiBin;N_{#gamma}^{RECO}(fake) / N_{#gamma}^{RECO}(all)",              numCentBins, CentBins_arr);
    h[43] = new TH1D("matchRatio_eta",   "matching efficiency = N_{#gamma}^{GEN}(matched) / N_{#gamma}^{GEN}(all);#eta^{GEN};N_{#gamma}^{GEN}(matched) / N_{#gamma}^{GEN}(all)",  numEtaBins,-maxEta,maxEta);
    h[44] = new TH1D("fakeRatio_eta",    "fake rate = N_{#gamma}^{RECO}(fake) / N_{#gamma}^{RECO}(all);#eta^{RECO};N_{#gamma}^{RECO}(fake) / N_{#gamma}^{RECO}(all)",             numEtaBins,-maxEta,maxEta);

    h[45] = new TH1D("trkIso_recoPhotons",         "track Iso - RECO photons;track Iso",100,-50,150);
    h[46] = new TH1D("trkIso_recoPhotons_matched", "track Iso - matched RECO photons (matched to GEN);track Iso" ,100,-50,150);
    h[47] = new TH1D("trkIso_recoPhotons_fake",    "track Iso - fake RECO photons (not matched to GEN);track Iso",100,-50,150);

//    h[48] = new TH1D("calIso_recoPhotons",         "calo Iso (pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4) - RECO photons;caloIso",100,-100,100);
//    h[49] = new TH1D("calIso_recoPhotons_matched", "calo Iso (pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4) - matched RECO photons (matched to GEN);caloIso" ,100,-100,100);
//    h[50] = new TH1D("calIso_recoPhotons_fake",    "calo Iso (pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4) - fake RECO photons (not matched to GEN);caloIso",100,-100,100);

    h[48] = new TH1D("calIso_recoPhotons",         "calo Iso - RECO photons;caloIso",100,-50,150);
    h[49] = new TH1D("calIso_recoPhotons_matched", "calo Iso - matched RECO photons (matched to GEN);caloIso" ,100,-50,150);
    h[50] = new TH1D("calIso_recoPhotons_fake",    "calo Iso - fake RECO photons (not matched to GEN);caloIso",100,-50,150);

    h[51] = new TH1D("phoR9_recoPhotons",         "phoR9 - RECO photons;phoR9",100,0,1);
    h[52] = new TH1D("phoR9_recoPhotons_matched", "phoR9 - matched RECO photons (matched to GEN);phoR9" ,100,0,1);
    h[53] = new TH1D("phoR9_recoPhotons_fake",    "phoR9 - fake RECO photons (not matched to GEN);phoR9",100,0,1);

    h[54] = new TH1D("phoHoverE_recoPhotons",         "phoHoverE - RECO photons;phoHoverE",100, 0, 2.5);
    h[55] = new TH1D("phoHoverE_recoPhotons_matched", "phoHoverE - matched RECO photons (matched to GEN);phoHoverE" ,100, 0, 2.5);
    h[56] = new TH1D("phoHoverE_recoPhotons_fake",    "phoHoverE - fake RECO photons (not matched to GEN);phoHoverE",100, 0, 2.5);

//    h[57] = new TH1D("phoSigmaIEtaIEta_recoPhotons",         "#sigma_{#eta#eta} (phoSigmaIEtaIEta) - RECO photons;phoSigmaIEtaIEta",100,0,0.1);
//    h[58] = new TH1D("phoSigmaIEtaIEta_recoPhotons_matched", "#sigma_{#eta#eta} (phoSigmaIEtaIEta) - matched RECO photons (matched to GEN);phoSigmaIEtaIEta" ,100,0,0.1);
//    h[59] = new TH1D("phoSigmaIEtaIEta_recoPhotons_fake",    "#sigma_{#eta#eta} (phoSigmaIEtaIEta) - fake RECO photons (not matched to GEN);phoSigmaIEtaIEta",100,0,0.1);

    h[57] = new TH1D("phoSigmaIEtaIEta_recoPhotons",         "#sigma_{#eta#eta} - RECO;phoSigmaIEtaIEta",100,0,0.1);
    h[58] = new TH1D("phoSigmaIEtaIEta_recoPhotons_matched", "#sigma_{#eta#eta} - matched RECO (matched to GEN);phoSigmaIEtaIEta" ,100,0,0.1);
    h[59] = new TH1D("phoSigmaIEtaIEta_recoPhotons_fake",    "#sigma_{#eta#eta} - fake RECO (not matched to GEN);phoSigmaIEtaIEta",100,0,0.1);

    h[60] = new TH1D("phoSigmaIPhiIPhi_recoPhotons",         "#sigma_{#phi#phi} (phoSigmaIPhiIPhi) - RECO photons;phoSigmaIPhiIPhi",100,-10,10);
    h[61] = new TH1D("phoSigmaIPhiIPhi_recoPhotons_matched", "#sigma_{#phi#phi} (phoSigmaIPhiIPhi) - matched RECO photons (matched to GEN);phoSigmaIPhiIPhi" ,100,-10,10);
    h[62] = new TH1D("phoSigmaIPhiIPhi_recoPhotons_fake",    "#sigma_{#phi#phi} (phoSigmaIPhiIPhi) - fake RECO photons (not matched to GEN);phoSigmaIPhiIPhi",100,-10,10);

    h[63] = new TH1D("energyScale", "energy scale - p_{T}^{#gamma}(RECO)/p_{T}^{#gamma}(GEN);RECO p_{T} / GEN p_{T}" ,100,0,3);

    h[64] = new TH1D("matchRatio_GENpT_v2",  "matching efficiency = N_{#gamma}^{GEN}(matched) / N_{#gamma}^{GEN}(all);GEN p_{T} (GeV);N_{#gamma}^{GEN}(matched) / N_{#gamma}^{GEN}(all)", 100,0,200);
    h[65] = new TH1D("fakeRatio_RECOpT_v2",   "fake rate = N_{#gamma}^{RECO}(fake) / N_{#gamma}^{RECO}(all);RECO p_{T} (GeV);N_{#gamma}^{RECO}(fake) / N_{#gamma}^{RECO}(all)",             100,0,200);
    h[66] = new TH1D("matchRatio_cent_v2",   "matching efficiency = N_{#gamma}^{GEN}(matched) / N_{#gamma}^{GEN}(all);hiBin;N_{#gamma}^{GEN}(matched) / N_{#gamma}^{GEN}(all)",  100,0,maxCent);
    h[67] = new TH1D("fakeRatio_cent_v2",    "fake rate = N_{#gamma}^{RECO}(fake) / N_{#gamma}^{RECO}(all);hiBin;N_{#gamma}^{RECO}(fake) / N_{#gamma}^{RECO}(all)",              100,0,maxCent);
    h[68] = new TH1D("matchRatio_eta_v2",   "matching efficiency = N_{#gamma}^{GEN}(matched) / N_{#gamma}^{GEN}(all);#eta^{GEN};N_{#gamma}^{GEN}(matched) / N_{#gamma}^{GEN}(all)",  100,-maxEta,maxEta);
    h[69] = new TH1D("fakeRatio_eta_v2",    "fake rate = N_{#gamma}^{RECO}(fake) / N_{#gamma}^{RECO}(all);#eta^{RECO};N_{#gamma}^{RECO}(fake) / N_{#gamma}^{RECO}(all)",             100,-maxEta,maxEta);

    // correlation plots
    TH2D* h2D[numHistos2D];
    h2D[0] = new TH2D("corr_trkIso", "E_{track}^{ISO}(RECO) vs. E_{track}^{ISO}(GEN);E_{track}^{ISO}(GEN);E_{track}^{ISO}(RECO)", 100, 0, 200, 125, -50, 200);
    h2D[1] = new TH2D("corr_calIso", "E_{calo}^{ISO}(RECO) vs. E_{calo}^{ISO}(GEN);E_{calo}^{ISO}(GEN);E_{calo}^{ISO}(RECO)", 100, 0, 200, 125, -50, 200);

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

        // event selection
        bool passedEvent = (hiBin >= hiBin_gt) && (hiBin < hiBin_lt);
        if(!passedEvent) continue;

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
            bool passedRECOselection = (phoEt->at(i) > cutptRECO); // && (TMath::Abs(phoEta->at(i)) >= cutEta_gt )
                                                                   // && (TMath::Abs(phoEta->at(i)) <  cutEta_lt );
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
                passedGENselection = (mcPID->at(j) == PDG_PHOTON) && (mcPt->at(j) > ptGEN) && (mcStatus->at(j) == cutmcStatus)
                                                                  && (TMath::Abs(mcEta->at(j)) >= eta_gt)
                                                                  && (TMath::Abs(mcEta->at(j)) <  eta_lt)
                                                                  && (TMath::Abs(mcMomPID->at(j)) > mcMomPID_gt)
                                                                  && (TMath::Abs(mcMomPID->at(j)) < mcMomPID_lt);

                deltaRtmp = getDR(phoEta->at(i), phoPhi->at(i), mcEta->at(j), mcPhi->at(j));
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

            bool passedRECOetaSelection = (TMath::Abs(phoEta->at(i)) >= eta_gt) && (TMath::Abs(phoEta->at(i)) <  eta_lt);
            if (passedRECOetaSelection)
            {
                h[0]->Fill(phoEt->at(i));         // all RECO photons
                h[45]->Fill(pho_trackIsoR4PtCut20->at(i));
                h[48]->Fill((pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)));
                h[51]->Fill(phoR9->at(i));
                h[54]->Fill(phoHoverE->at(i));
                h[57]->Fill(phoSigmaIEtaIEta->at(i));
                //            h[60]->Fill(phoSigmaIPhiIPhi->at(i));
            }

            if (isFake[i] == false) {       // RECO photons matched to GEN photon

                int j = matched[i];         // index of matched GEN photon
                deltaPhi[i]  = getDPHI(phoPhi->at(i),mcPhi->at(j));
                deltaEta[i]  = phoEta->at(i) - mcEta->at(j);
                deltaR  [i]  = deltaRMin;
                deltaTrkIso[i] = pho_trackIsoR4PtCut20->at(i) - mcTrkIsoDR04->at(j);
                deltaCalIso[i] = (pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)) - mcCalIsoDR04->at(j);
                ratioTrkIso[i] = 0;
                if(mcTrkIsoDR04->at(j)!=0)
                    ratioTrkIso[i] = pho_trackIsoR4PtCut20->at(i) / mcTrkIsoDR04->at(j);
                ratioCalIso[i] = 0;
                if(mcCalIsoDR04->at(j)!=0)
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

                // isolation differences as fnc. of GEN pT
                for(int r=0; r < numGENptBins; ++r){
                    if(mcPt->at(j)<=GENptBins_Iso[r])
                    {
                        trkIso_GENpT[0][r] += deltaTrkIso[i];
                        trkIso_GENpT[1][r] += deltaTrkIso[i]*deltaTrkIso[i];
                        trkIso_GENpT[2][r] += 1;

                        calIso_GENpT[0][r] += deltaCalIso[i];
                        calIso_GENpT[1][r] += deltaCalIso[i]*deltaCalIso[i];
                        calIso_GENpT[2][r] += 1;
                        break;
                    }
                }

                // isolation differences as fnc. of centrality
                for(int r=0; r < numCentBins; ++r){
                    if(hiBin <= CentBins_Iso[r])
                    {
                        trkIso_cent[0][r] += deltaTrkIso[i];
                        trkIso_cent[1][r] += deltaTrkIso[i]*deltaTrkIso[i];
                        trkIso_cent[2][r] += 1;

                        calIso_cent[0][r] += deltaCalIso[i];
                        calIso_cent[1][r] += deltaCalIso[i]*deltaCalIso[i];
                        calIso_cent[2][r] += 1;
                        break;
                    }
                }

                // fake ratio as fnc. of RECO pT
                for(int r=0; r < numGENptBins; ++r){
                    if(phoEt->at(i)<=RECOptBins_fakeRatio[r])
                    {
                        // matched RECO do not go into index 0 // fakeRatio_RECOpT[0][r] += 1;
                        fakeRatio_RECOpT[2][r] += 1;
                        break;
                    }
                }

                // fake ratio as fnc. of eta
                for(int r=0; r < numEtaBins; ++r){
                    if(phoEta->at(i) <= etaBins_fakeRatio[r])
                    {
                        // matched RECO do not go into index 0 // fakeRatio_eta[0][r] += 1;
                        fakeRatio_eta[2][r] += 1;
                        break;
                    }
                }

                // fake ratio as fnc. of centrality
                for(int r=0; r < numCentBins; ++r){
                    if(hiBin <= CentBins_fakeRatio[r])
                    {
                        // matched RECO do not go into index 0 // fakeRatio_cent[0][r] += 1;
                        fakeRatio_cent[2][r] += 1;
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

                h[46]->Fill(pho_trackIsoR4PtCut20->at(i));
                h[49]->Fill((pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)));
                h[52]->Fill(phoR9->at(i));
                h[55]->Fill(phoHoverE->at(i));
                h[58]->Fill(phoSigmaIEtaIEta->at(i));
//                h[61]->Fill(phoSigmaIPhiIPhi->at(i));
                h[63]->Fill(phoEt->at(i)/mcPt->at(j));

                h2D[0]->Fill(mcTrkIsoDR04->at(j), pho_trackIsoR4PtCut20->at(i));
                h2D[1]->Fill(mcCalIsoDR04->at(j), (pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)));
            }
            else {
                if(passedRECOetaSelection){
                    // fake RECO photons
                    h[2]->Fill(phoEt->at(i));
                    h[47]->Fill(pho_trackIsoR4PtCut20->at(i));
                    h[50]->Fill((pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)));
                    h[53]->Fill(phoR9->at(i));
                    h[56]->Fill(phoHoverE->at(i));
                    h[59]->Fill(phoSigmaIEtaIEta->at(i));
                    //                h[62]->Fill(phoSigmaIPhiIPhi->at(i));

                    // fake ratio as fnc. of RECO pT
                    for(int r=0; r < numGENptBins; ++r){
                        if(phoEt->at(i)<=RECOptBins_fakeRatio[r])
                        {
                            fakeRatio_RECOpT[0][r] += 1;
                            fakeRatio_RECOpT[2][r] += 1;
                            break;
                        }
                    }

                    // fake ratio as fnc. of eta
                    for(int r=0; r < numEtaBins; ++r){
                        if(phoEta->at(i) <= etaBins_fakeRatio[r])
                        {
                            fakeRatio_eta[0][r] += 1;
                            fakeRatio_eta[2][r] += 1;
                            break;
                        }
                    }

                    // fake ratio as fnc. of centrality
                    for(int r=0; r < numCentBins; ++r){
                        if(hiBin <= CentBins_fakeRatio[r])
                        {
                            fakeRatio_cent[0][r] += 1;
                            fakeRatio_cent[2][r] += 1;
                            break;
                        }
                    }
                }
            }
        }

        // find RECO photons that match to GEN photons
        bool isMiss[nMC];
        for(int i = 0; i < nMC; ++i){

            // default values
            isMiss[i]=true;

            bool passedGENselection;      // selections for GEN photon
            passedGENselection = (mcPID->at(i) == PDG_PHOTON) && (mcPt->at(i) > ptGEN) && (mcStatus->at(i) == cutmcStatus)
                                                              && (TMath::Abs(mcEta->at(i)) >= eta_gt)
                                                              && (TMath::Abs(mcEta->at(i)) <  eta_lt)
                                                              && (TMath::Abs(mcMomPID->at(i)) > mcMomPID_gt)
                                                              && (TMath::Abs(mcMomPID->at(i)) < mcMomPID_lt);

            if(passedGENselection)
            {
                // check if that GEN photon was matched to a RECO photon in the previous loop over RECO photons
                int* indexRECO = std::find(matched,matched+nPho, i);
                bool found = (indexRECO < matched+nPho);
                isMiss[i] = !found;

                // fill histograms
                h[3]->Fill(mcPt->at(i));
                if (isMiss[i] == false) {
                    h[4]->Fill(mcPt->at(i));

                    // matching efficiency as fnc. of GEN pT
                    for(int r=0; r < numGENptBins; ++r){
                        if(mcPt->at(i) <= GENptBins_matchRatio[r])
                        {
                            matchRatio_GENpT[0][r] += 1;
                            matchRatio_GENpT[2][r] += 1;
                            break;
                        }
                    }

                    // matching efficiency as fnc. of centrality
                    for(int r=0; r < numCentBins; ++r){
                        if(hiBin <= CentBins_matchRatio[r])
                        {
                            matchRatio_cent[0][r] += 1;
                            matchRatio_cent[2][r] += 1;
                            break;
                        }
                    }

                    // matching efficiency as fnc. of eta
                    for(int r=0; r < numEtaBins; ++r){
                        if(mcEta->at(i) <= etaBins_matchRatio[r])
                        {
                            matchRatio_eta[0][r] += 1;
                            matchRatio_eta[2][r] += 1;
                            break;
                        }
                    }
                }
                else {
                    h[5]->Fill(mcPt->at(i));

                    // matching efficiency as fnc. of GEN pT
                    for(int r=0; r < numGENptBins; ++r){
                        if(mcPt->at(i) <= GENptBins_matchRatio[r])
                        {
                            // non-matched GEN do not go into index 0 // matchRatio_GENpT[0][r] += 1;
                            matchRatio_GENpT[2][r] += 1;
                            break;
                        }
                    }

                    // matching efficiency as fnc. of centrality
                    for(int r=0; r < numCentBins; ++r){
                        if(hiBin <= CentBins_matchRatio[r])
                        {
                            // non-matched GEN do not go into index 0 // matchRatio_cent[0][r] += 1;
                            matchRatio_cent[2][r] += 1;
                            break;
                        }
                    }

                    // matching efficiency as fnc. of eta
                    for(int r=0; r < numEtaBins; ++r){
                        if(mcEta->at(i) <= etaBins_matchRatio[r])
                        {
                            // non-matched GEN do not go into index 0 // matchRatio_eta[0][r] += 1;
                            matchRatio_eta[2][r] += 1;
                            break;
                        }
                    }
                }
            }
        }
    } // exited event loop

    // fill remaining histograms
    // histograms as fnc. of GEN pT
    for(int r=0; r < numGENptBins; ++r){
        int n;
        double mean;
        double meanOfSquares;
        double ratio;

        // energy scale as fnc. of GEN pT
        n = energyScale[2][r];
        if(n>0)
        {
            mean          = energyScale[0][r]/n;
            meanOfSquares = energyScale[1][r]/n;
            h[6]->SetBinContent(r+1,mean);
            h[6]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[7]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[7]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean)/TMath::Sqrt(n));
        }

        // isolation ratios as fnc. of GEN pT
        n = trkIso_ratio_GENpT[2][r];
        if (n>0) {
            mean          = trkIso_ratio_GENpT[0][r]/n;
            meanOfSquares = trkIso_ratio_GENpT[1][r]/n;
            h[23]->SetBinContent(r+1,mean);
            h[23]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[25]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[25]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean)/TMath::Sqrt(n));
        }

        n = calIso_ratio_GENpT[2][r];
        if (n>0) {
            mean          = calIso_ratio_GENpT[0][r]/n;
            meanOfSquares = calIso_ratio_GENpT[1][r]/n;
            h[24]->SetBinContent(r+1,mean);
            h[24]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[26]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[26]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean)/TMath::Sqrt(n));
        }

        // isolation differences as fnc. of GEN pT
        n = trkIso_ratio_GENpT[2][r];
        if (n>0) {
            mean          = trkIso_GENpT[0][r]/n;
            meanOfSquares = trkIso_GENpT[1][r]/n;
            h[31]->SetBinContent(r+1,mean);
            h[31]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[33]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[33]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean)/TMath::Sqrt(n));
        }

        n = calIso_ratio_GENpT[2][r];
        if (n>0) {
            mean          = calIso_GENpT[0][r]/n;
            meanOfSquares = calIso_GENpT[1][r]/n;
            h[32]->SetBinContent(r+1,mean);
            h[32]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[34]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[34]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean)/TMath::Sqrt(n));
        }

        // matching efficiency as fnc. of GEN pT
        n = matchRatio_GENpT[2][r];
        if (n>0) {
            ratio = matchRatio_GENpT[0][r]/n;
            h[39]->SetBinContent(r+1,ratio);
            h[39]->SetBinError  (r+1,ratio/TMath::Sqrt(n));
        }

        // fake ratio as fnc. of RECO pT
        n = fakeRatio_RECOpT[2][r];
        if (n>0) {
            ratio = fakeRatio_RECOpT[0][r]/n;
            h[40]->SetBinContent(r+1,ratio);
            h[40]->SetBinError  (r+1,ratio/TMath::Sqrt(n));
        }
    }

    // histograms as fnc. of centrality
    for(int r=0; r < numCentBins; ++r){
        int n;
        double mean;
        double meanOfSquares;
        double ratio;

        // energy scale as fnc. of centrality
        n = energyScale_cent[2][r];
        if (n>0) {
            mean          = energyScale_cent[0][r]/n;
            meanOfSquares = energyScale_cent[1][r]/n;
            h[8]->SetBinContent(r+1,mean);
            h[8]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[9]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[9]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean)/TMath::Sqrt(n));
        }

        // position resolution as fnc. of centrality
        n = deltaPhi_cent[2][r];
        if (n>0) {
            mean          = deltaPhi_cent[0][r]/n;
            meanOfSquares = deltaPhi_cent[1][r]/n;
            h[13]->SetBinContent(r+1,mean);
            h[13]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[16]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[16]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean)/TMath::Sqrt(n));
        }

        n = deltaEta_cent[2][r];
        if (n>0) {
            mean          = deltaEta_cent[0][r]/n;
            meanOfSquares = deltaEta_cent[1][r]/n;
            h[14]->SetBinContent(r+1,mean);
            h[14]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[17]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[17]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean)/TMath::Sqrt(n));
        }

        n = deltaR_cent[2][r];
        if (n>0) {
            mean          = deltaR_cent[0][r]/n;
            meanOfSquares = deltaR_cent[1][r]/n;
            h[15]->SetBinContent(r+1,mean);
            h[15]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[18]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[18]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean)/TMath::Sqrt(n));
        }
        // isolation ratios as fnc. of centrality
        n = trkIso_ratio_cent[2][r];
        if (n>0) {
            mean          = trkIso_ratio_cent[0][r]/n;
            meanOfSquares = trkIso_ratio_cent[1][r]/n;
            h[27]->SetBinContent(r+1,mean);
            h[27]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[29]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[29]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean)/TMath::Sqrt(n));
        }

        n = calIso_ratio_cent[2][r];
        if (n>0) {
            mean          = calIso_ratio_cent[0][r]/n;
            meanOfSquares = calIso_ratio_cent[1][r]/n;
            h[28]->SetBinContent(r+1,mean);
            h[28]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[30]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[30]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean)/TMath::Sqrt(n));
        }

        // isolation differences as fnc. of centrality
        n = trkIso_cent[2][r];
        if (n>0) {
            mean          = trkIso_cent[0][r]/n;
            meanOfSquares = trkIso_cent[1][r]/n;
            h[35]->SetBinContent(r+1,mean);
            h[35]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[37]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[37]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean)/TMath::Sqrt(n));
        }

        n = calIso_cent[2][r];
        if (n>0) {
            mean          = calIso_cent[0][r]/n;
            meanOfSquares = calIso_cent[1][r]/n;
            h[36]->SetBinContent(r+1,mean);
            h[36]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[38]->SetBinContent(r+1,TMath::Sqrt(meanOfSquares-mean*mean));
            h[38]->SetBinError(r+1,TMath::Sqrt(meanOfSquares-mean*mean)/TMath::Sqrt(n));
        }

        // matching efficieny as fnc. of centrality
        n = matchRatio_cent[2][r];
        if (n>0) {
            ratio = matchRatio_cent[0][r]/n;
            h[41]->SetBinContent(r+1,ratio);
            h[41]->SetBinError  (r+1,ratio/TMath::Sqrt(n));
        }

        // fake ratio as fnc. of centrality
        n = fakeRatio_cent[2][r];
        if (n>0) {
            ratio = fakeRatio_cent[0][r]/n;
            h[42]->SetBinContent(r+1,ratio);
            h[42]->SetBinError  (r+1,ratio/TMath::Sqrt(n));
        }
    }

    // histograms as fnc. of eta
    for(int r=0; r < numEtaBins; ++r){
        int n;
        double ratio;

        // matching efficieny as fnc. of eta
        n = matchRatio_eta[2][r];
        if (n>0) {
            ratio = matchRatio_eta[0][r]/n;
            h[43]->SetBinContent(r+1,ratio);
            h[43]->SetBinError  (r+1,ratio/TMath::Sqrt(n));
        }

        // fake ratio as fnc. of eta
        n = fakeRatio_eta[2][r];
        if (n>0) {
            ratio = fakeRatio_cent[0][r]/n;
            h[44]->SetBinContent(r+1,ratio);
            h[44]->SetBinError  (r+1,ratio/TMath::Sqrt(n));
        }
    }

    end_loop = std::clock();
    std::cout.precision(6);      // get back to default precision
    std::cout << "LOOP finished in             : " << (end_loop - start_loop) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "exited event loop" << std::endl;

    std::vector<TH1*> histos;
    for (int i=0; i<numHistos; ++i) {
        histos.push_back(h[i]);
    }
    for (int i=0; i<numHistos2D; ++i) {
        histos.push_back(h2D[i]);
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
