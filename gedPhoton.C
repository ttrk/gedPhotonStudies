/*
 */

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include <iostream>
#include <ctime>
#include <iomanip>
#include <vector>

Double_t getDR( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
Double_t getDPHI( Double_t phi1, Double_t phi2);
void gedPhoton(const char* hiForestfileName = "HiForest.root", const char* outputFileName = "gedPhotonTest.root");

const int MAXGENPARTICLES = 50000;  // number of gen particles can be larges
const int MAXPHOTONS = 500;
const int PDG_PHOTON = 22;
const double cutdeltaR = 0.2; // 0.3    // cut for matching gen and reco. particles
const double cutptGEN = 15 ;

void gedPhoton(const char* hiForestfileName, const char* outputFileName)
{
    TFile* inputFile = new TFile(hiForestfileName, "READ");
    std::cout << "input HiForest : " << inputFile->GetName() << std::endl;

    TTree* HiGenParticleAnaTree = (TTree*)inputFile->Get("HiGenParticleAna/hi");
    TTree* ggHiNtuplizerTree = (TTree*)inputFile->Get("ggHiNtuplizer/EventTree");
    TTree* ggHiNtuplizerGEDTree = (TTree*)inputFile->Get("ggHiNtuplizerGED/EventTree");

    Int_t npart;
    Int_t mult;
    Float_t pt_gen[MAXGENPARTICLES];
    Float_t eta_gen[MAXGENPARTICLES];
    Float_t phi_gen[MAXGENPARTICLES];
    Int_t pdg_gen[MAXGENPARTICLES];

    HiGenParticleAnaTree->SetBranchAddress("npart",&npart);
    HiGenParticleAnaTree->SetBranchAddress("mult",&mult);
    HiGenParticleAnaTree->SetBranchAddress("pt",pt_gen);
    HiGenParticleAnaTree->SetBranchAddress("eta",eta_gen);
    HiGenParticleAnaTree->SetBranchAddress("phi",phi_gen);
    HiGenParticleAnaTree->SetBranchAddress("pdg",pdg_gen);

    Int_t nPho;
    std::vector<float>* phoEt=0;
    std::vector<float>* phoEta=0;
    std::vector<float>* phoPhi=0;

    ggHiNtuplizerTree->SetBranchAddress("nPho",&nPho);
    ggHiNtuplizerTree->SetBranchAddress("phoEt",&phoEt);
    ggHiNtuplizerTree->SetBranchAddress("phoEta",&phoEta);
    ggHiNtuplizerTree->SetBranchAddress("phoPhi",&phoPhi);

    Int_t nPho_GED;
    std::vector<float>* phoEt_GED=0;
    std::vector<float>* phoEta_GED=0;
    std::vector<float>* phoPhi_GED=0;

    ggHiNtuplizerGEDTree->SetBranchAddress("nPho",&nPho_GED);
    ggHiNtuplizerGEDTree->SetBranchAddress("phoEt",&phoEt_GED);
    ggHiNtuplizerGEDTree->SetBranchAddress("phoEta",&phoEta_GED);
    ggHiNtuplizerGEDTree->SetBranchAddress("phoPhi",&phoPhi_GED);


    TFile* outputFile=new TFile(outputFileName, "RECREATE");

    const int numGENptBins = 20;
    const double maxGENpt  = 200;
    double GENptBins_energyScale[numGENptBins];     // upper edges of energy scale histograms
    double energyScale[3][numGENptBins][2];            // energyScale[0][][] = sum of energy scales
                                                       // energyScale[1][][] = sum of square of energy scales
                                                       // energyScale[2][][] = number of energy scales
                                                       // energyScale[][][0] = RECO photon
                                                       // energyScale[][][1] = GED RECO photon
    for (int i=0; i<numGENptBins; ++i)
    {
        GENptBins_energyScale[i]= maxGENpt/numGENptBins*(i+1);
        for (int j=0; j<3; ++j){
            for (int k=0; k<2; ++k){
                energyScale[j][i][k]=0;
            }
        }
    }

    const int numCentBins = 22;
    const int maxCent  = 440;
    double CentBins_energyScale[numCentBins];     // upper edges of energy scale histograms
    double energyScale_cent[3][numCentBins][2];        // energyScale[0][][] = sum of energy scales
                                                       // energyScale[1][][] = sum of square of energy scales
                                                       // energyScale[2][][] = number of energy scales
                                                       // energyScale[][][0] = RECO photon
                                                       // energyScale[][][1] = GED RECO photon
    for (int i=0; i<numCentBins; ++i)
    {
        CentBins_energyScale[i]= (double)maxCent/numCentBins*(i+1);
        for (int j=0; j<3; ++j){
            for (int k=0; k<2; ++k){
                energyScale_cent[j][i][k]=0;
            }
        }
    }

    TH1D* h[10][2];
    // histograms for RECO photons
    h[0][0] = new TH1D("recoPhotons","RECO photons;p_{T} (GeV)",200,0,200);
    h[1][0] = new TH1D("fakePhotons","fake RECO photons (not matched to GEN);p_{T} (GeV)",200,0,200);
    h[2][0] = new TH1D("energyScale","energy scale;GEN p_{T} (GeV);<RECO p_{T} / GEN p_{T}>",numGENptBins,0,maxGENpt);
    h[3][0] = new TH1D("widthEnergyScale","width of energy scale;GEN p_{T} (GeV);#sigma(RECO p_{T} / GEN p_{T})",numGENptBins,0,maxGENpt);
    h[4][0] = new TH1D("energyScale","energy scale;N_{Part};<RECO p_{T} / GEN p_{T}>",numCentBins,0,maxCent);
    h[5][0] = new TH1D("widthEnergyScale","width of energy scale;N_{Part};#sigma(RECO p_{T} / GEN p_{T})",numCentBins,0,maxCent);
    h[6][0] = new TH1D("genPhotons","GEN photons;p_{T} (GeV)",200,0,200);
    h[7][0] = new TH1D("missingPhotons","missing GEN photons (not matched to RECO);p_{T} (GeV)",200,0,200);

    // histograms for GED RECO photons
    for (int i=0; i < 8; ++i)
    {
        h[i][1] = (TH1D*)h[i][0]->Clone(Form("%s_GED",h[i][0]->GetName()));
        h[i][1]->SetTitle(Form("GED - %s",h[i][0]->GetTitle()));
    }

    TH1D* recoPhotons[2];
    recoPhotons[0] = new TH1D("recoPhotons","RECO photons;p_{T} (GeV)",100,0,100);
    recoPhotons[1] = (TH1D*)recoPhotons[0]->Clone(Form("%s_GED",recoPhotons[0]->GetName()));

    TH1D* genPhotons;
    genPhotons = new TH1D("genPhotons","GEN photons;p_{T} (GeV)",100,0,100);

    TH1D* fakePhotons[2];
    fakePhotons[0] = new TH1D("fakePhotons","# of non-matched reco photons (fakes);Reco p_{T} (GeV)",100,0,100);
    fakePhotons[1] = (TH1D*)fakePhotons[0]->Clone(Form("%s_GED",fakePhotons[0]->GetName()));

    TH1D* missedPhotons[2];
    missedPhotons[0] = new TH1D("missedPhotons","# of non-matched gen photons (misses);Gen p_{T} (GeV)",100,0,100);
    missedPhotons[1] = (TH1D*)missedPhotons[0]->Clone(Form("%s_GED",missedPhotons[0]->GetName()));

    std::cout << "entering event loop" << std::endl;
    Long64_t entries = HiGenParticleAnaTree->GetEntries();
    std::cout << "number of entries = " << entries << std::endl;

    std::clock_t    start_loop, end_loop;
    start_loop = std::clock();
    for(Long64_t jj = 0; jj < entries; ++jj)
    {
        if (jj % 100000 == 0)  {
          std::cout << "current entry = " <<jj<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)jj/entries*100<<" %"<<std::endl;
        }

        HiGenParticleAnaTree->GetEntry(jj);
        ggHiNtuplizerTree->GetEntry(jj);
        ggHiNtuplizerGEDTree->GetEntry(jj);

        int nReco = nPho;
        if (nPho_GED > nPho)
            nReco = nPho_GED;

        int nPhotons [2] = {nPho, nPho_GED};
        // [][0]   = RECO photon
        // [][1]   = GED RECO photon
        bool isFake[nReco][2];
        int  matched[nReco][2];         // index of the matched GEN photon in  to this RECO photon
        double deltaR[nReco][2];        // deltaR between the RECO photon and the matched GEN photon

        // find GEN photons that match to RECO photons
        for (int k = 0; k<2; ++k)
        {
            int n = nPhotons[k];
            for (int i=0; i < n; ++i)
            {
                isFake[i][k] = true;
                matched[i][k] = -1;
                deltaR[i][k] = 999;

                bool passedGENselection;      // selections for GEN photon
                bool passedDR;                // selections for GEN photon
                double deltaRtmp;
                for (int j=0; j<mult; ++j)
                {
                    deltaRtmp = getDR(eta_gen[j],phi_gen[j], phoEta->at(i), phoPhi->at(i));
                    passedGENselection = (pdg_gen[j] == PDG_PHOTON) && (pt_gen[j] > cutptGEN);
                    passedDR           = (deltaRtmp < cutdeltaR);

                    if (passedGENselection && passedDR)
                    {
                        // matched GEN photon is the one closest to the RECO photon
                        if (deltaRtmp < deltaR[i])
                        {
                            matched[i][k] = j;
                            deltaR [i][k] = deltaRtmp;
//                            isMiss [j][k] = false;

                            // energy scale as fnc. of GEN pT
                            for(int r=0; r<numGENptBins; ++r){
                                if(pt_gen[j]<=GENptBins_energyScale[r])
                                {
                                    energyScale[0][r][0] += phoEt->at(i)/pt_gen[j];
                                    energyScale[1][r][0] += (phoEt->at(i)/pt_gen[j])*(phoEt->at(i)/pt_gen[j]);
                                    energyScale[2][r][0] += 1;
                                    break;
                                }
                            }

                            // energy scale as fnc. of centrality
                            for(int r=0; r<numCentBins; ++r){
                                if(npart<=CentBins_energyScale[r])
                                {
                                    energyScale_cent[0][r][0] += phoEt->at(i)/pt_gen[j];
                                    energyScale_cent[1][r][0] += (phoEt->at(i)/pt_gen[j])*(phoEt->at(i)/pt_gen[j]);
                                    energyScale_cent[2][r][0] += 1;
                                    break;
                                }
                            }
                        }
                    }
                }
                isFake[i][k] = (matched[i][k] == -1);     // if no matched GEN photon, then this RECO photon is fake.

                // fill histograms
                recoPhotons[k]->Fill(phoEt->at(i));
                if (isFake[i][k] == true)  fakePhotons[k]->Fill(phoEt->at(i));

                h[0][k]->Fill(phoEt->at(i));         // all RECO photons
                if (isFake[i] == false) {
                    // RECO photons matched to GEN photon

                }
                else {
                    // fake RECO photons
                    h[1][k]->Fill(phoEt->at(i));
                }
            }
        }

        // find RECO photons that match to GEN photons
        bool isMiss[mult][2];
        // [][0]   = RECO photon
        // [][1]   = GED RECO photon
        for (int k = 0; k<2; ++k){
            for(int i = 0; i < mult; ++i){

                isMiss[i][2]=true;

                bool passedGENselection;      // selections for GEN photon
                passedGENselection = (pdg_gen[j] == PDG_PHOTON) && (pt_gen[j] > cutptGEN);

                if(passedGENselection)
                {
                    // check if that GEN photon was matched to a RECO photon in the previous loop over RECO photons
                    int* indexRECO = std::find(matched[0][k],matched[0][k]+nReco, i);
                    bool found = (indexRECO < matched[0][k]+nReco);
                    isMiss[i][k] = !found;


                    // fill histograms
                    h[6][k]->Fill(pt_gen[i]);
                    if (isMiss[i][k] == true)  h[7][k]->Fill(pt_gen[i]);
                }
            }
        }
    }
    // fill remaining histograms
    // energy scale as fnc. of GEN pT
    for(int r=1; r<=numGENptBins; ++r){
        for(int k=0; k<2; ++k)
        {
            int n = energyScale[2][r][k];
            double mean = energyScale[0][r][k]/n;
            double meanOfSquares = energyScale[1][r][k]/n;
            h[2][k]->SetBinContent(r,mean);
            h[3][k]->SetBinContent(r,TMath::Sqrt(meanOfSquares-mean*mean));
        }
    }
    // energy scale as fnc. of centrality
    for(int r=1; r<=numCentBins; ++r){
        for(int k=0; k<2; ++k)
        {
            int n = energyScale_cent[2][r][k];
            double mean = energyScale_cent[0][r][k]/n;
            double meanOfSquares = energyScale_cent[1][r][k]/n;
            h[4][k]->SetBinContent(r,mean);
            h[5][k]->SetBinContent(r,TMath::Sqrt(meanOfSquares-mean*mean));
        }
    }


    end_loop = std::clock();
    std::cout.precision(6);      // get back to default precision
    std::cout << "LOOP finished in             : " << (end_loop - start_loop) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "exited event loop" << std::endl;

    outputFile->cd();
    outputFile->Write();    // write all histograms to the output file.

    outputFile->Close();
    inputFile->Close();
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
