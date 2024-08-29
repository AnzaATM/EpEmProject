#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "fastjet/PseudoJet.hh"
#include "download/header.h"
#include "fastjet/ClusterSequence.hh"
#endif

void Delphes()
{
    gSystem->Load("libDelphes");
    TString infile = "/home/anza-tshilidzi/Downloads/MG5_aMC_v3_5_5/Z95bb_mumu_back/Events/run_08/tag_1_delphes_events.root";
    //TString infile = "/home/anza-tshilidzi/Downloads/MG5_aMC_v3_5_5/Higgs95_Sig/Events/run_07/tag_1_delphes_events.root"; 
  
    // Create chain of root trees
    TChain chain("Delphes");
    chain.Add(infile);
    chain.Print();

    // Create object of class ExRootTreeReader
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();
    cout << numberOfEntries << endl;
    Double_t weight;

    // Declare variables
    double Ne = 0, Nb = 0, Nj = 0;
    double ThetaMuP, ThetaMuM;
    double PhiMuP, PhiMuM;
    double ThetaBjet1, ThetaBjet2;
    double PhiBjet1, PhiBjet2;
    double EnrMuP, EnrMuM; 
    double EnrBjet1, EnrBjet2;  
    double Mrecoil; 
    double InMmuPM;
    double InMbb;

    // Histograms and trees
    TH1 *histPTe = new TH1F("pte", "pte", 20, 0.0, 400.0);
    TFile *hfile = TFile::Open("/home/anza-tshilidzi/Downloads/MG5_aMC_v3_5_5/Z95bb_mumu_back/Events/run_08/Back_events300k2.root", "RECREATE");
    //TFile *hfile = TFile::Open("/home/anza-tshilidzi/Downloads/MG5_aMC_v3_5_5/Higgs95_Sig/Events/run_07/tag_1_delphes_events.root", "RECREATE");
    TTree *tree = new TTree("tree","tree");

    // Tree branches for DNN features
    tree->Branch("EnrMuP", &EnrMuP, "EnrMuP/D");
    tree->Branch("EnrMuM", &EnrMuM, "EnrMuM/D");
    tree->Branch("ThetaMuP", &ThetaMuP, "ThetaMuP/D");
    tree->Branch("ThetaMuM", &ThetaMuM, "ThetaMuM/D");
    tree->Branch("PhiMuP", &PhiMuP, "PhiMuP/D");
    tree->Branch("PhiMuM", &PhiMuM, "PhiMuM/D");
    tree->Branch("EnrBjet1", &EnrBjet1, "EnrBjet1/D");
    tree->Branch("EnrBjet2", &EnrBjet2, "EnrBjet2/D");
    tree->Branch("ThetaBjet1", &ThetaBjet1, "ThetaBjet1/D");
    tree->Branch("ThetaBjet2", &ThetaBjet2, "ThetaBjet2/D");
    tree->Branch("PhiBjet1", &PhiBjet1, "PhiBjet1/D");
    tree->Branch("PhiBjet2", &PhiBjet2, "PhiBjet2/D");
    tree->Branch("InMmuPM", &InMmuPM, "InMmuPM/D");
    tree->Branch("InMbb", &InMbb, "InMbb/D");
    tree->Branch("Mrecoil", &Mrecoil, "Mrecoil/D");
    tree->Branch("weight", &weight);

    // Get branches from the Delphes file
    TClonesArray *branchEvent = treeReader->UseBranch("Event");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchJet = treeReader->UseBranch("Jet");

    for (Int_t entry = 0; entry < numberOfEntries; ++entry) {
        
        treeReader->ReadEntry(entry);
        
        HepMCEvent *event = (HepMCEvent*) branchEvent->At(0);
        weight = event->Weight;

        TLorentzVector vecMuP, vecMuM, vecBJet1, vecBJet2;    
        std::vector<Muon*> muPlusVec;
        std::vector<Muon*> muMinusVec;
        std::vector<Jet*> bJetVec; 

        // Initial four-momentum of e+e- system 
        TLorentzVector vece1(0.0, 0.0, 125.0, 125.0); // (px, py, pz, E)
        TLorentzVector vece2(0.0, 0, -125.0, 125.0);
        Double_t SqrtS = 2 * vece1.E(); // COM energy of e+e- system

        //=============================================================
        // Muon branch
        //=============================================================
        // Separate muons into mu+ and mu-
        for(int i = 0; i < branchMuon->GetEntries(); ++i)
        {
            Muon *muon = (Muon*) branchMuon->At(i);
            if(muon->Charge == 1)
            {
                muPlusVec.push_back(muon);
            }
            else if(muon->Charge == -1)
            {
                muMinusVec.push_back(muon);
            }
        }
                                
        // Sort mu+ and mu- by PT
        std::sort(muPlusVec.begin(), muPlusVec.end(), [](Muon* a, Muon* b) { return a->PT > b->PT; });
        std::sort(muMinusVec.begin(), muMinusVec.end(), [](Muon* a, Muon* b) { return a->PT > b->PT; });                        
        
        if(muPlusVec.empty() || muMinusVec.empty()) {
            continue;
        }
        
        // Get the highest PT mu+ and mu-
        Muon *muonPlus = muPlusVec[0];
        Muon *muonMinus = muMinusVec[0];                        
        vecMuP = muonPlus->P4();
        vecMuM = muonMinus->P4();                        
        
        InMmuPM = (vecMuP + vecMuM).M(); // Invariant mass of final state mu+ and mu-
        EnrMuP = vecMuP.E(); // Energy of final mu+
        EnrMuM = vecMuM.E(); // Energy of final mu-
        ThetaMuP = vecMuP.Theta(); // Polar angle of mu+
        ThetaMuM = vecMuM.Theta(); // Polar angle of mu-
        PhiMuP = vecMuP.Phi(); // Azimuthal angle of mu+
        PhiMuM = vecMuM.Phi(); // Azimuthal angle of mu-
        Mrecoil = sqrt(pow(SqrtS, 2) + pow(InMmuPM, 2.0) - 2 * (EnrMuP + EnrMuM) * SqrtS);  // Recoil mass of mu+ mu-                      

        //=====================================================                       
        // Jet branch
        //=====================================================
        // Separate b-jets
        for(int i = 0; i < branchJet->GetEntries(); ++i)
        {
            Jet *jet = (Jet*) branchJet->At(i);
            if(jet->BTag)
            {
                bJetVec.push_back(jet);
            }
        }
        
        // Sort b-jets by PT
        std::sort(bJetVec.begin(), bJetVec.end(), [](Jet* a, Jet* b) { return a->PT > b->PT; });
        if(bJetVec.size() < 2) {       
            continue;
        }      
        
        // Get the two highest PT b-jets
        Jet *Bjet1 = bJetVec[0];
        Jet *Bjet2 = bJetVec[1];                   
        vecBJet1 = Bjet1->P4();
        vecBJet2 = Bjet2->P4();                    
        
        InMbb = (vecBJet1 + vecBJet2).M();  // Invariant mass of b1 and b2 jets                      
        EnrBjet1 = vecBJet1.E(); // Energy of b1 jet
        EnrBjet2 = vecBJet2.E(); // Energy of b2 jet                  
        ThetaBjet1 = vecBJet1.Theta(); // Polar angle of b1 jet
        ThetaBjet2 = vecBJet2.Theta(); // Polar angle of b2 jet
        PhiBjet1 = vecBJet1.Phi(); // Azimuthal angle of b1 jet
        PhiBjet2 = vecBJet2.Phi(); // Azimuthal angle of b2 jet                       
                                          
        if(Mrecoil < 120 && InMbb < 100 && EnrMuP > 5 && EnrMuM > 5 && EnrBjet1 > 5 && EnrBjet2 > 5)
        { 
            tree->Fill();
        }
    }

    // Write and close the tree
    tree->Write();
    hfile->Close();
}
