#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <TSystem.h>
#include <TROOT.h>
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldInvert.h"
#include <TH1D.h>
#include <TTree.h>
#include <TChain.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPave.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TLine.h>
#include <TText.h>
#include <iostream>
#include <fstream>
#include <TClonesArray.h>
#include <cmath>
#include <string>
#include "TRandom.h"
#include "TString.h"
#include "TLeaf.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "Math/GenVector/PxPyPzM4D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TVector2.h"
#endif

using namespace std;

// Debug flag - set to true for verbose output
const bool DEBUG = false;

TChain* CreateChain(const string& inFile, const char* treeName) {
    TChain* chain = new TChain(treeName);
    if (inFile.find(".txt") != string::npos) {
        ifstream file(inFile);
        if (!file.is_open()) {
            cerr << "Error: Could not open file list " << inFile << endl;
            return chain;
        }
        string rootFile;
        while (getline(file, rootFile)) {
            if (!rootFile.empty()) {
                if (DEBUG) cout << "Adding: " << rootFile << endl;
                int added = chain->Add(rootFile.c_str());
                if (added == 0) {
                    cerr << "Error: Could not add file " << rootFile << endl;
                }
            }
        }
    } else {
        if (DEBUG) cout << "Adding single file: " << inFile << endl;
        int added = chain->Add(inFile.c_str());
        if (added == 0) {
            cerr << "Error: Could not add file " << inFile << endl;
        }
    }
    if (DEBUG) cout << "Chain has " << chain->GetEntries() << " entries" << endl;
    return chain;
}
// Calculate the angle between two 2D vectors defined by their (PX, PY) components
float GetAngle(float vPX1, float vPY1, float vPX2, float vPY2) {
    TVector2 DielectronVect(vPX1, vPY1);
    TVector2 ElectronVect(vPX2, vPY2);
    
    float DielectronMag = sqrt(DielectronVect.X()*DielectronVect.X() + DielectronVect.Y()*DielectronVect.Y());
    float ElectronMag = sqrt(ElectronVect.X()*ElectronVect.X() + ElectronVect.Y()*ElectronVect.Y());
    
    float cosAngle = (DielectronVect.X()*ElectronVect.X() + DielectronVect.Y()*ElectronVect.Y()) / (DielectronMag * ElectronMag);
    float sinAngle = (DielectronVect.X()*ElectronVect.Y() - DielectronVect.Y()*ElectronVect.X()) / (DielectronMag * ElectronMag);
    float tanAngle = atan2(sinAngle, cosAngle);
    
    return tanAngle;
}

// Configuration parameters
const float xrangeMax = TMath::Pi();
const float xrangeMin = -TMath::Pi();
const int xbinz = 18;
const float maxpT = 1.0;
const float minpT = 0.3;
const float MaxRapidity = 2.4;
const float eleMass = 0.000511;
const float minMass = 5.0;
const float maxMass = 40.0;

void ManuelElectronResponseMatrix() {
    string inputMC = "/afs/cern.ch/user/j/jmbagwu/RooUnfold/examples/lblfiles/cms/my_file.txt";
    string outfile = "PEE2SCenhanced_superchic2018_dielectrons_full.root";

    // Check output directory exists
    if (gSystem->AccessPathName(gSystem->DirName(outfile.c_str()))) {
        cerr << "Error: Output directory does not exist: " << gSystem->DirName(outfile.c_str()) << endl;
        return;
    }

    TChain *t1 = CreateChain(inputMC, "ggHiNtuplizer/EventTree");
    if (t1->GetEntries() == 0) {
        cerr << "Error: No entries found in the chain!" << endl;
        return;
    }

    TFile f(outfile.c_str(), "recreate");
    if (!f.IsOpen()) {
        cerr << "Error: Could not create output file " << outfile << endl;
        return;
    }

    // Create the SuperEventTree
    TTree* SuperEventTree = new TTree("SuperEventTree", "Tree with Event Variables");

    // --- Generator level electron branches
    vector<float> *gen_elePt = nullptr;
    vector<float> *gen_elePhi = nullptr;
    vector<float> *gen_eleEta = nullptr;
    vector<int>   *gen_elePID = nullptr;
    vector<float> *gen_trkvx = nullptr;
    vector<float> *gen_trkvy = nullptr;
    vector<float> *gen_trkvz = nullptr;
    int gen_nEle = 0;

    t1->SetBranchAddress("mcPt",   &gen_elePt);
    t1->SetBranchAddress("mcPhi",  &gen_elePhi);
    t1->SetBranchAddress("mcEta",  &gen_eleEta);
    t1->SetBranchAddress("nMC",    &gen_nEle);
    t1->SetBranchAddress("mcPID",  &gen_elePID);
    t1->SetBranchAddress("mcVtx_x", &gen_trkvx);
    t1->SetBranchAddress("mcVtx_y", &gen_trkvy);
    t1->SetBranchAddress("mcVtx_z", &gen_trkvz);

    // --- Reconstructed level electron branches
    vector<float> *reco_elePt = nullptr;
    vector<float> *reco_elePhi = nullptr;
    vector<float> *reco_eleEta = nullptr;
    vector<int>   *reco_eleCharge = nullptr;
    vector<float> *reco_eleDz = nullptr;
    vector<float> *reco_trkvx = nullptr;
    vector<float> *reco_trkvy = nullptr;
    vector<float> *reco_trkvz = nullptr;
    int reco_nEle = 0;
    int reco_nTrk = 0;

    t1->SetBranchAddress("elePt",    &reco_elePt);
    t1->SetBranchAddress("elePhi",   &reco_elePhi);
    t1->SetBranchAddress("eleEta",   &reco_eleEta);
    t1->SetBranchAddress("nEle",     &reco_nEle);
    t1->SetBranchAddress("eleCharge",&reco_eleCharge);
    t1->SetBranchAddress("eleDz",    &reco_eleDz);
    t1->SetBranchAddress("nTrk",     &reco_nTrk);
    t1->SetBranchAddress("trkvx",    &reco_trkvx);
    t1->SetBranchAddress("trkvy",    &reco_trkvy);
    t1->SetBranchAddress("trkvz",    &reco_trkvz);

    // --- Calorimeter branches
    vector<float> *CaloTower_e = nullptr;
    vector<float> *CaloTower_et = nullptr;
    vector<float> *CaloTower_eta = nullptr;
    vector<float> *CaloTower_hadE = nullptr;
    vector<float> *CaloTower_emE = nullptr;
    vector<float> *CaloTower_phi = nullptr;
    int nTower = 0;
    t1->SetBranchAddress("CaloTower_e",    &CaloTower_e);
    t1->SetBranchAddress("CaloTower_et",   &CaloTower_et);
    t1->SetBranchAddress("CaloTower_eta",  &CaloTower_eta);
    t1->SetBranchAddress("CaloTower_hadE", &CaloTower_hadE);
    t1->SetBranchAddress("CaloTower_emE",  &CaloTower_emE);
    t1->SetBranchAddress("CaloTower_phi",  &CaloTower_phi);
    t1->SetBranchAddress("nTower",         &nTower);

    // Event variables to store
    float maxHFp = 0;
    float maxHFm = 0;
    float maxTower_E = 0;
    float maxTower_eta = 0;
    float maxTower_phi = 0;
    SuperEventTree->Branch("maxHFp", &maxHFp);
    SuperEventTree->Branch("maxHFm", &maxHFm);
    SuperEventTree->Branch("maxTower_E", &maxTower_E);
    SuperEventTree->Branch("maxTower_eta", &maxTower_eta);
    SuperEventTree->Branch("maxTower_phi", &maxTower_phi);

    // Load delta phi histogram for electrons (data)
    TFile *dataFile = new TFile("/afs/cern.ch/user/j/jmbagwu/RooUnfold/examples/lblfiles/cms/DATA/RooUnfold/examples/delta_phi_histogram_electron.root", "READ");
    if (!dataFile || dataFile->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // 2. Get the original histogram (read-only)
    TH1F *hDeltaPhiData = (TH1F*)dataFile->Get("hDeltaPhi");
    if (!hDeltaPhiData) {
        std::cerr << "Histogram hDeltaPhi not found in file!" << std::endl;
        return;
    }

    // 3. Create working copies with clear naming
    TH1F *hDeltaPhiRaw = (TH1F*)hDeltaPhiData->Clone("hDeltaPhiRaw"); // Exact copy of original
    TH1F *hDeltaPhiPurityCorrected = (TH1F*)hDeltaPhiRaw->Clone("hDeltaPhiPurityCorrected");

    // Initialize response matrix
    RooUnfoldResponse response(xbinz, xrangeMin, xrangeMax);

    // Event counters
    int count_True = 0;
    int count_Fake = 0;
    int count_Miss = 0;
    int Event_count = 0;
    int GEventPass = 0;
    int REventPass = 0;

    // Histograms
    TH1F *hGenMass = new TH1F("hGenMass", "Generator Level Dielectron Mass", 200, 0, 50);
    TH1F *hGenMassCuts = new TH1F("hGenMassCuts", "Generator Level Dielectron Mass After Cuts", 200, 0, maxMass);
    TH1F *hRecoMass = new TH1F("hRecoMass", "Reconstructed Dielectron Mass", 200, 0, 50);
    TH1F *hRecoMassCuts = new TH1F("hRecoMassCuts", "Reconstructed Dielectron Mass After Cuts", 200, 0, maxMass);

    TH1F *hreco = new TH1F("hreco", "Reconstructed #Delta#phi", xbinz, xrangeMin, xrangeMax);
    TH1F *hrecoTrue = new TH1F("hrecoTrue", "True Reconstructed #Delta#phi", xbinz, xrangeMin, xrangeMax);
    TH1F *hrecoTrue_split1 = new TH1F("hrecoTrue_split1", "True Reco Split 1 (even events)", xbinz, xrangeMin, xrangeMax);
    TH1F *hrecoTrue_split2 = new TH1F("hrecoTrue_split2", "True Reco Split 2 (odd events)", xbinz, xrangeMin, xrangeMax);
    TH1F *hrecoFake = new TH1F("hrecoFake", "Fake Reconstructed #Delta#phi", xbinz, xrangeMin, xrangeMax);
    TH1F *hgen = new TH1F("hgen", "Generator Level #Delta#phi", xbinz, xrangeMin, xrangeMax);
    TH1F *hgenTrue = new TH1F("hgenTrue", "True Generator Level #Delta#phi", xbinz, xrangeMin, xrangeMax);
    TH1F *hgenMiss = new TH1F("hgenMiss", "Missed Generator Level #Delta#phi", xbinz, xrangeMin, xrangeMax);
    
    // Efficiency and purity histograms
    TH1F* hEfficiency = new TH1F("hEfficiency", "Efficiency;#Delta#phi;Efficiency", xbinz, xrangeMin, xrangeMax);
    TH1F* hPurity = new TH1F("hPurity", "Purity;#Delta#phi;Purity", xbinz, xrangeMin, xrangeMax);

    float RDeltaPhi = 0, GDeltaPhi = 0;
    int numberentry = t1->GetEntries();
    int split1_count = 0;
    int split2_count = 0;

    // Main event loop
    for (int j = 0; j < numberentry; j++) {
        t1->GetEntry(j);

        // Reset event flags
        int FakeEvt = 0;
        int TrueEvt = 0;
        int MissEvt = 0;
        int passGenCuts = 0;
        int passRecoCuts = 0;

        // Initialize indices to invalid values
        int iRp = -1, iRn = -1;

        // Declare reconstructed variables at the start of event processing
        TLorentzVector RLep1, RLep2, RDiele;

        // --- Generator-level dielectron construction
        int idx_ele_minus = -1;  // index of electron (e-)
        int idx_ele_plus = -1;   // index of positron (e+)

        for (size_t k = 0; k < gen_elePID->size(); ++k) {
            int pid = gen_elePID->at(k);

            // Skip any particle with |PID| >= 20 to filter out unwanted particles
            if (abs(pid) >= 20) continue;

            if (pid == -11 && idx_ele_minus == -1) {
                idx_ele_minus = k;  // Found electron
            }
            else if (pid == 11 && idx_ele_plus == -1) {
                idx_ele_plus = k;   // Found positron
            }
        }

        // If either electron or positron not found, skip this event
        if (idx_ele_minus == -1 || idx_ele_plus == -1) continue;
        // Construct gen-level 4-vectors
        TLorentzVector GLep1, GLep2, GDiele;
        GLep1.SetPtEtaPhiM(gen_elePt->at(idx_ele_minus), gen_eleEta->at(idx_ele_minus),
                           gen_elePhi->at(idx_ele_minus), eleMass);
        GLep2.SetPtEtaPhiM(gen_elePt->at(idx_ele_plus), gen_eleEta->at(idx_ele_plus),
                           gen_elePhi->at(idx_ele_plus), eleMass);
        GDiele = GLep1 + GLep2;

        // Gen-level cuts
        bool genTwoelectron = (gen_nEle >= 2);
        bool genElectroncharge = (gen_elePID->at(idx_ele_minus) + gen_elePID->at(idx_ele_plus) == 0);
        bool GRapidity = abs(GLep1.Eta()) <= MaxRapidity && abs(GLep2.Eta()) <= MaxRapidity;
        bool GMassCut = (GDiele.M() > minMass && GDiele.M() < maxMass);
        
        //generator-level pT cut
//        bool GpTcut = (GDiele.Pt() < maxpT);
        bool GpTcut = (GDiele.Pt() < maxpT && GDiele.Pt() > minpT);
        bool GVertexCut = (gen_trkvz->at(idx_ele_minus) <= 20);
        bool GenSingleElectronPtCut = GLep1.Pt() > 1 && GLep2.Pt() > 1;

        // Calculate gen-level delta phi
        float GenDPx = GDiele.Pt() * cos(GDiele.Phi());
        float GenDPy = GDiele.Pt() * sin(GDiele.Phi());
        float GenMPx = GLep1.Pt() * cos(GLep1.Phi());
        float GenMPy = GLep1.Pt() * sin(GLep1.Phi());
        GDeltaPhi = GetAngle(GenDPx, GenDPy, GenMPx, GenMPy);

        hGenMass->Fill(GDiele.M());

        // Apply gen-level selection
        if (GMassCut && GRapidity && GpTcut && genTwoelectron && genElectroncharge && GVertexCut && GenSingleElectronPtCut) {
            passGenCuts = 1;
            GEventPass++;
            hgen->Fill(GDeltaPhi);
            hGenMassCuts->Fill(GDiele.M());
        }

        // --- Reco-level selections ---
        bool Twoelectron = (reco_nEle == 2);
        bool Twotrk = (reco_nTrk == 2);

        if (Twoelectron && Twotrk) {
            
            // Check charge consistency
            if (reco_eleCharge->size() < 2) continue;

            bool Electroncharge = (reco_eleCharge->at(0) * reco_eleCharge->at(1) == -1);
            if (!Electroncharge) continue;

            // Determine charge indices
            if (reco_eleCharge->at(0) == 1 && reco_eleCharge->at(1) == -1) {
                iRp = 0;
                iRn = 1;
            } else if (reco_eleCharge->at(1) == 1 && reco_eleCharge->at(0) == -1) {
                iRp = 1;
                iRn = 0;
            } else {
                continue;
            }

            // Construct reco-level 4-vectors
            TLorentzVector RLep1, RLep2, RDiele;
            RLep1.SetPtEtaPhiM(reco_elePt->at(iRp), reco_eleEta->at(iRp),
                               reco_elePhi->at(iRp), eleMass);
            RLep2.SetPtEtaPhiM(reco_elePt->at(iRn), reco_eleEta->at(iRn),
                               reco_elePhi->at(iRn), eleMass);
            RDiele = RLep1 + RLep2;

            // Calculate reco-level delta phi
            float RecoDPx = RDiele.Pt() * cos(RDiele.Phi());
            float RecoDPy = RDiele.Pt() * sin(RDiele.Phi());
            float RecoMPx = RLep1.Pt() * cos(RLep1.Phi());
            float RecoMPy = RLep1.Pt() * sin(RLep1.Phi());
            RDeltaPhi = GetAngle(RecoDPx, RecoDPy, RecoMPx, RecoMPy);

            hRecoMass->Fill(RDiele.M());

            // Tower energy calculations
            maxHFp = 0;
            maxHFm = 0;
            for (int t = 0; t < nTower; t++) {
                float tower_eta = CaloTower_eta->at(t);
                float tower_E = CaloTower_e->at(t);

                if (tower_eta > 2.9 && tower_eta < 5.2)
                    maxHFp = TMath::Max(maxHFp, tower_E);
                if (tower_eta < -2.9 && tower_eta > -5.2)
                    maxHFm = TMath::Max(maxHFm, tower_E);
            }

            // Reco-level cuts
            bool RRapidity = abs(reco_eleEta->at(0)) <= MaxRapidity && abs(reco_eleEta->at(1)) <= MaxRapidity;
            bool RVertexCut = (reco_trkvz->at(0) <= 20);
            bool RMassCut = (RDiele.M() > minMass && RDiele.M() < maxMass);
            
            //reconstructed-level pT cut,
//            bool RpTCut = (RDiele.Pt() < maxpT);
            bool RpTCut = (RDiele.Pt() < maxpT && RDiele.Pt() > minpT);
            bool HFCuts = (maxHFm <= 6.0 && maxHFp <= 6.0);
            bool RecoSingleElectronPtCut = RLep1.Pt() > 2 && RLep2.Pt() > 2;

            // ✅ Apply reco-level selections and fill
            if (RMassCut && HFCuts && RVertexCut && RRapidity && RecoSingleElectronPtCut && RpTCut) {
                hreco->Fill(RDeltaPhi);
                hRecoMassCuts->Fill(RDiele.M());
                REventPass++;
                passRecoCuts = 1;
            }
        }

        // --- Event classification
        if (passGenCuts == 1 && passRecoCuts == 0) MissEvt = 1;
        if (passGenCuts == 0 && passRecoCuts == 1) FakeEvt = 1;
        if (passGenCuts == 1 && passRecoCuts == 1) TrueEvt = 1;

        // Debug output if enabled
        if (DEBUG && TrueEvt) {
            cout << "GDeltaPhi is " << GDeltaPhi << endl;
            cout << "RDeltaPhi is " << RDeltaPhi << endl;
            cout << "G lep 1 charge is " << gen_elePID->at(idx_ele_plus) << endl;
            if (iRn >= 0) {  // Check if reconstruction was successful
                cout << "R lep1 charge is " << reco_eleCharge->at(iRn) << endl;
                cout << "R lep1 is " << RLep2.Phi() << endl;
                cout << "RDielectron phi is " << RDiele.Phi() << endl;
            }
            cout << "G lep 1 is " << GLep1.Phi() << endl;
            cout << "GDielectron phi is " << GDiele.Phi() << endl;
        }

        // --- Fill histograms based on event classification
            if (MissEvt == 1) {
                count_Miss++;
                hgenMiss->Fill(GDeltaPhi);
            }
            if (FakeEvt == 1) {
                count_Fake++;
                hrecoFake->Fill(RDeltaPhi);
            }
            if (TrueEvt == 1) {
                count_True++;
                response.Fill(RDeltaPhi, GDeltaPhi);
                hrecoTrue->Fill(RDeltaPhi);
//                hgenTrue->Fill(GDeltaPhi);
                
                
//                if (gRandom->Rndm() < 0.5) {
//                    split1_count++;
//                    hrecoTrue_split1->Fill(RDeltaPhi);
//                } else {
//                    split2_count++;
//                    hrecoTrue_split2->Fill(RDeltaPhi);
//                }

                // --- Split sample for cross-checks
                if (j % 2 == 0) {
                    split1_count++;
                    hrecoTrue_split1->Fill(RDeltaPhi);
                } else {
                    split2_count++;
                    hrecoTrue_split2->Fill(RDeltaPhi);
                }
                
                hgenTrue->Fill(GDeltaPhi);
            }

            Event_count++;
            SuperEventTree->Fill();
        } // end loop over entries
    
    
    if (DEBUG) {
        cout << "Split1 entries: " << split1_count << endl;
        cout << "Split2 entries: " << split2_count << endl;
        cout << "True events: " << count_True << endl;
        cout << "Fake events: " << count_Fake << endl;
        cout << "Missed events: " << count_Miss << endl;
    }
    
    // Correct generator-level histogram for inefficiencies
    for (int iBinX = 1; iBinX <= hreco->GetNbinsX(); ++iBinX) {
        double genVal = hgen->GetBinContent(iBinX);
        double trueGen = hgenTrue->GetBinContent(iBinX);
        double missedGen = hgenMiss->GetBinContent(iBinX);
        double denomGen = trueGen + missedGen;

        if (denomGen > 0) {
            hgen->SetBinContent(iBinX, genVal * (trueGen / denomGen));
        }
    }
       
    // Apply purity corrections: purity = true / (true + fake)
    for (int iBinX = 1; iBinX <= hDeltaPhiPurityCorrected->GetNbinsX(); iBinX++) {
        double trueCounts = hrecoTrue->GetBinContent(iBinX);
        double fakeCounts = hrecoFake->GetBinContent(iBinX);
        double denom = trueCounts + fakeCounts;
        
        if (denom > 0) {
            double purity = trueCounts / denom;
            
            // Get current bin content and error before correction
            double oldContent = hDeltaPhiPurityCorrected->GetBinContent(iBinX);
            double oldError = hDeltaPhiPurityCorrected->GetBinError(iBinX);
            
            // Apply purity correction to content and propagate error (scaling)
            hDeltaPhiPurityCorrected->SetBinContent(iBinX, oldContent * purity);
            hDeltaPhiPurityCorrected->SetBinError(iBinX, oldError * purity);
        }
        else {
            // Handle empty bins - set to zero or keep original?
            hDeltaPhiPurityCorrected->SetBinContent(iBinX, 0);
            hDeltaPhiPurityCorrected->SetBinError(iBinX, 0);
        }
    }


    // Calculate and fill efficiency and purity histograms
    for (int iBinX = 1; iBinX <= hreco->GetNbinsX(); ++iBinX) {
        double trueGen = hgenTrue->GetBinContent(iBinX);
        double missedGen = hgenMiss->GetBinContent(iBinX);
        double trueReco = hrecoTrue->GetBinContent(iBinX);
        double fakeReco = hrecoFake->GetBinContent(iBinX);

        // Efficiency = trueReco / (trueGen + missedGen)
        double efficiency = 0.;
        double denomEff = trueGen + missedGen;
        if (denomEff > 0) efficiency = trueReco / denomEff;
        hEfficiency->SetBinContent(iBinX, efficiency);

        // Purity = trueReco / (trueReco + fakeReco)
        double purity = 0.;
        double denomPur = trueReco + fakeReco;
        if (denomPur > 0) purity = trueReco / denomPur;
        hPurity->SetBinContent(iBinX, purity);
    }

    // --- Split test check on hrecoTrue
    TCanvas* cSplitTest = new TCanvas("cSplitTest", "Split Test Check", 800, 600);
    
    // Style settings
    hrecoTrue_split1->SetLineColor(kBlue);
    hrecoTrue_split2->SetLineColor(kRed);
    hrecoTrue_split1->SetMarkerStyle(20);
    hrecoTrue_split2->SetMarkerStyle(21);
    
    // Normalize to same area for comparison
    
    if (hrecoTrue_split1->Integral() > 0) hrecoTrue_split1->Scale(hrecoTrue->Integral() / hrecoTrue_split1->Integral());
    if (hrecoTrue_split2->Integral() > 0) hrecoTrue_split2->Scale(hrecoTrue->Integral() / hrecoTrue_split2->Integral());
    
    // Draw the comparison
    hrecoTrue->SetLineColor(kBlack);
    hrecoTrue->SetLineWidth(2);
    hrecoTrue->Draw("HIST");
    hrecoTrue_split1->Draw("HIST SAME");
    hrecoTrue_split2->Draw("HIST SAME");
    
    // Add legend
    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(hrecoTrue, "2018 Superchic recoTrue", "l");
    leg->AddEntry(hrecoTrue_split1, "Split 1 (recoTrue)", "l");
    leg->AddEntry(hrecoTrue_split2, "Split 2 (recoTrue)", "l");
    leg->Draw();
    
    // Calculate chi2 between splits
    double chi2 = 0;
    int ndf = 0;
    for (int i=1; i<=hrecoTrue_split1->GetNbinsX(); i++) {
        double diff = hrecoTrue_split1->GetBinContent(i) - hrecoTrue_split2->GetBinContent(i);
        double err1 = hrecoTrue_split1->GetBinError(i);
        double err2 = hrecoTrue_split2->GetBinError(i);
        double err = sqrt(err1*err1 + err2*err2);
        if (err > 0) {
            chi2 += (diff*diff)/(err*err);
            ndf++;
        }
    }
    
    // Add chi2 information to the plot
    TPaveText* pt1 = new TPaveText(0.15, 0.75, 0.45, 0.85, "NDC");
    pt1->SetFillColor(0);
    pt1->SetBorderSize(1);
    pt1->AddText(Form("#chi^{2}/ndf = %.2f/%d", chi2, ndf));
    pt1->AddText(Form("P-value = %.3f", TMath::Prob(chi2, ndf)));
    pt1->Draw();
    
    cSplitTest->Update();

    // --- Unfolding with RooUnfoldBayes for data
    RooUnfoldBayes unfold_bayes_data(&response, hDeltaPhiPurityCorrected, 4);
    TH1* UnfoldData = unfold_bayes_data.Hunfold();
    UnfoldData->SetName("UnfoldData2018");
    UnfoldData->SetTitle("Unfold 2018 Data");
    UnfoldData->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e-}");

    // Forward fold the unfolded data to compare with original data
    TH1* hForwardFoldedData = (TH1*)response.ApplyToTruth(UnfoldData, "hForwardFoldedData");
    hForwardFoldedData->SetTitle("Forward Fold Unfolded 2018 Data");
    hForwardFoldedData->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e-}");

    // Create ratio of forward folded data to original data
    TH1* hRatioData = (TH1*)hForwardFoldedData->Clone("hRatioData");
    hRatioData->Divide(hDeltaPhiPurityCorrected);
    hRatioData->SetTitle("Ratio: Forward Fold Data / 2018 Data");
    hRatioData->GetYaxis()->SetTitle("Ratio");

    // --- Monte Carlo closure test
    RooUnfoldBayes unfold_bayes_mc(&response, hrecoTrue, 4);
    TH1* UnfoldMC = unfold_bayes_mc.Hunfold();
    UnfoldMC->SetName("UnfoldSCMC2018");
    UnfoldMC->SetTitle("Unfold 2018 Superchic MC True");
    UnfoldMC->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e-}");
    
    // Forward fold the unfolded MC
    TH1* hForwardFoldedMC = (TH1*)response.ApplyToTruth(UnfoldMC, "hForwardFoldedMC");
    hForwardFoldedMC->SetTitle("Forward Fold 2018 Superchic MC True");
    hForwardFoldedMC->GetXaxis()->SetTitle("#Delta#phi = #phi_{ee} - #phi_{e-}");
    
    // Create ratio of forward folded MC to original MC
    TH1* hRatioMC = (TH1*)hForwardFoldedMC->Clone("hRatioMC");
    hRatioMC->Divide(hrecoTrue);
    hRatioMC->SetTitle("Ratio: Forward Fold 2018 Superchic MC / True Reco 2018 Superchic MC");
    hRatioMC->GetYaxis()->SetTitle("Ratio");
    

    // Bottom-line test plotting
    TCanvas *cDataUnfold = new TCanvas("cDataUnfold", "Bottom Line Test - Data", 900, 700);

    if (!hDeltaPhiPurityCorrected || !hForwardFoldedData) {
        std::cerr << "One or both histograms are null!" << std::endl;
        return;
    }

    std::cout << "Corrected Integral: " << hDeltaPhiPurityCorrected->Integral() << std::endl;
    std::cout << "Forward Folded Integral: " << hForwardFoldedData->Integral() << std::endl;

    // Disable stats box only for this canvas
    hDeltaPhiPurityCorrected->SetStats(kFALSE);
    hForwardFoldedData->SetStats(kFALSE);


    hDeltaPhiPurityCorrected->SetLineColor(kBlack);
    hDeltaPhiPurityCorrected->SetMarkerStyle(20);
    hDeltaPhiPurityCorrected->SetLineWidth(2);
    hDeltaPhiPurityCorrected->SetTitle("2018 Data vs Forward Fold with Superchic 2018");
    hDeltaPhiPurityCorrected->GetXaxis()->SetTitle("#Delta#phi");

    double maxY = std::max(hDeltaPhiPurityCorrected->GetMaximum(), hForwardFoldedData->GetMaximum());
    if (maxY > 0)
        hDeltaPhiPurityCorrected->SetMaximum(1.2 * maxY);

    hDeltaPhiPurityCorrected->Draw("E");

    hForwardFoldedData->SetLineColor(kRed);
    hForwardFoldedData->SetMarkerStyle(21);
    hForwardFoldedData->SetLineWidth(2);
    hForwardFoldedData->Draw("E SAME");

    TLegend *legData = new TLegend(0.65, 0.7, 0.9, 0.85);
    legData->AddEntry(hDeltaPhiPurityCorrected, "2018 Data", "lep");
    legData->AddEntry(hForwardFoldedData, "Forward Fold", "lep");
    legData->Draw();

    cDataUnfold->Update();
    cDataUnfold->Draw();
    cDataUnfold->SaveAs("bottom_line_test.pdf"); // optional


    // Add chi2 information box on this canvas
    TPaveText* pt2 = new TPaveText(0.15, 0.75, 0.45, 0.85, "NDC");
    pt2->SetFillColor(0);
    pt2->SetBorderSize(1);
    pt2->AddText(Form("#chi^{2}/ndf = %.2f/%d", chi2, ndf));
    pt2->AddText(Form("P-value = %.3f", TMath::Prob(chi2, ndf)));
    pt2->Draw();

    cDataUnfold->Update();

    // Canvas for data ratio
    TCanvas *cDataRatio = new TCanvas("cDataRatio", "2018 Data Ratio ForwardFolded/Corrected", 900, 700);
    hRatioData->SetLineColor(kBlue);
    hRatioData->SetMarkerStyle(20);
    hRatioData->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatioData->Draw("E");
    cDataRatio->Update();

    // Canvas to compare MC before and after unfolding+folding (closure test)
    TCanvas *cMCUnfold = new TCanvas("cMCUnfold", "Bottom Line Test - MC Closure", 900, 700);
    hrecoTrue->SetLineColor(kBlack);
    hrecoTrue->SetMarkerStyle(20);
    hrecoTrue->SetTitle("2018 Superchic MC: recoTrue Split Test Check");
    hrecoTrue->GetXaxis()->SetTitle("#Delta#phi");
    hrecoTrue->Draw("E");

    hForwardFoldedMC->SetLineColor(kRed);
    hForwardFoldedMC->SetMarkerStyle(21);
    hForwardFoldedMC->Draw("E SAME");

    TLegend *legMC = new TLegend(0.65,0.7,0.9,0.85);
    legMC->AddEntry(hrecoTrue, "Reco True MC", "lep");
    legMC->AddEntry(hForwardFoldedMC, "Forward Folded Unfolded 2018 Superchic MC", "lep");
    legMC->Draw();

    cMCUnfold->Update();

    // Canvas for MC ratio
    TCanvas *cMCRatio = new TCanvas("cMCRatio", "MC Ratio ForwardFolded/RecoTrue", 900, 700);
    hRatioMC->SetLineColor(kBlue);
    hRatioMC->SetMarkerStyle(20);
    hRatioMC->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatioMC->Draw("E");
    cMCRatio->Update();

    // --- Write all objects to output file
    f.cd();
    
    // Response matrix
    auto* hResponse = response.HresponseNoOverflow();
    if (hResponse) {
        // First write the response matrix to file as before
        hResponse->SetTitle("");
        hResponse->GetXaxis()->SetTitle("Superchic 2018 Gen-Level #Delta#phi = #phi_{ee} - #phi_{e-}");
        hResponse->GetYaxis()->SetTitle("Superchic 2018 Reco-Level #Delta#phi = #phi_{ee} - #phi_{e-}");
        hResponse->Write("hResponseMatrix");
        response.Write("responseObject");

        // Then create and save a nice plot
        TCanvas *c = new TCanvas("cResponse", "Response Matrix", 800, 600);
        gStyle->SetOptStat();  // Disable stats box for 2D plots (optional)
        hResponse->Draw("COLZ");

        gPad->Update();
        TPaveStats *stat = (TPaveStats*)hResponse->FindObject("stats");
        if (stat) {
            stat->SetX1NDC(0.60);
            stat->SetX2NDC(0.90);
            stat->SetY1NDC(0.28);
            stat->SetY2NDC(0.42);
            stat->SetTextSize(0.03);
            gPad->Modified();
            gPad->Update();
        }

        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(62);
        latex.SetTextSize(0.04);
        latex.SetTextAlign(13);
        latex.DrawLatex(0.10, 0.93, "CMS");
        latex.SetTextFont(42);
        latex.SetTextSize(0.035);
        latex.DrawLatex(0.16, 0.93, "#it{work in progress}");
        latex.SetTextSize(0.038);
        latex.SetTextAlign(33);
        latex.DrawLatex(0.89, 0.935, "PbPb 2018 #sqrt{#it{s}_{NN}} = 5.02 TeV");
        latex.SetTextAlign(31);
        latex.SetTextSize(0.030);
        latex.DrawLatex(0.84, 0.24, "p_{T}_{ee} < 1.0 GeV");
//        latex.DrawLatex(0.84, 0.22, "p_{T}^{ee} < 1.0 GeV");     // First line
        latex.DrawLatex(0.84, 0.20, "p_{T}_{e} > 2.0 GeV");      // Second line
        latex.DrawLatex(0.84, 0.17, "|#eta_{e}| < 2.4");             // Third line
        latex.DrawLatex(0.84, 0.13, "M_{ee} > 5 GeV");           // Fourth line
//        latex.DrawLatex(0.84, 0.28, "p_{T}_{ee} < 1.0 GeV");
//        latex.DrawLatex(0.84, 0.24, "p_{T}_{e} > 2.0 GeV");
//        latex.DrawLatex(0.84, 0.20, "|#eta| < 2.4");
//        latex.DrawLatex(0.84, 0.16, "M_{ee} > 5 GeV");

        // Save the canvas
        c->Update();
        c->Draw();
        c->Print("PEE2SCSuperchicResponseMatrix.pdf");

        // Cleanup
        delete c;
    }

    // Ensure correct axis titles before writing
    hDeltaPhiPurityCorrected->GetXaxis()->SetTitle("#Delta#phi"); // After purity correction
    hDeltaPhiPurityCorrected->Write();
    hDeltaPhiRaw->Write();                     // Original data
    hForwardFoldedData->GetXaxis()->SetTitle("#Delta#phi");
    
    // Data unfolding results
    UnfoldData->Write();
    hForwardFoldedData->Write();
    hRatioData->Write();


    // MC unfolding results
    UnfoldMC->Write();
    hForwardFoldedMC->Write();
    hRatioMC->Write();

    // Split test results
    hrecoTrue->Write();
    hrecoTrue_split1->Write();
    hrecoTrue_split2->Write();
    cSplitTest->Write();

    // Other histograms
    hGenMass->Write();
    hGenMassCuts->Write();
    hRecoMass->Write();
    hRecoMassCuts->Write();
    hreco->Write();
    hrecoFake->Write();
    hgen->Write();
    hgenTrue->Write();
    hgenMiss->Write();
    hEfficiency->Write();
    hPurity->Write();

    // Event tree
    SuperEventTree->Write();
    
    f.Close();
}

#ifndef __CINT__
int main() { 
    ManuelElectronResponseMatrix(); 
    return 0; 
}
#endif

