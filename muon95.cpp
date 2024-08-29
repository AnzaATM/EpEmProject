/*
Simple macro showing how to access branches from the delphes output root file for Snowmass studies,
loop over events, and plot simple quantities such as the leading electron and jet pt.

root -l examples/muon.C'("/home/pramod/Desktop/Higgs95/back/Events/run_01/tag_1_delphes_events.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY (libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm> 
#include "TLegend.h"
#include "TH1F.h"
using namespace std;
//#include <bits/stdc++.h>
#else
class ExRootTreeReader;
class ExRootResult;
#endif
//------------------------------------------------------------------------------

void muon95()//(const char *inputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  //chain.Add(inputFile);
  //chain.Add("/home/anza-tshilidzi/Downloads/MG5_aMC_v3_5_5/Higgs95_Sig/Events/run_05/tag_1_delphes_events.root"); 
  chain.Add("/home/anza-tshilidzi/Downloads/MG5_aMC_v3_5_5/Z95bb_mumu_back/Events/run_06/tag_1_delphes_events.root"); 
  
  

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");


  TH1D *hMreco = new TH1D("Recoil Mass", "Recoil mass of mu+ mu-", 100, 70.0, 120.0);
  TH1D *hInvMuMu = new TH1D("Invarinat mass of mu+ mu-", "Invarinat mass of mu+ mu-", 100, 70.0, 130.0);
  TH1D *hInvbjj = new TH1D("Invarinat mass of b-jets", "Invarinat mass of b-jets", 100, 0.0, 130.0);
  TH1D *hEMup = new TH1D("Energy of mu+", "Energy of final mu+", 100, 0.0, 120.0);
  TH1D *hEMum = new TH1D("Energy of mu-", "Energy of final mu-", 100, 0.0, 120.0);


// Int_t num;
// num=0;
  Double_t eventWeight;
  // Loop over all entries
  for(Long64_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
//    cout << numberOfEntries << endl;
//    std::cout << "Total number of muons in event " << entry << ": " << branchMuon->GetEntries() << std::endl;

    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    eventWeight = event->Weight;
//      cout << event->Weight << endl;
    
    TLorentzVector vecMuP, vecMuM, vecBJet1, vecBJet2;
    
    std::vector<Muon*> muPlusVec;
    std::vector<Muon*> muMinusVec;
    std::vector<Jet*> bJetVec;   
    
/*      std::cout << "Event " << entry << ":" << std::endl;
      for(int i = 0; i < branchMuon->GetEntries(); ++i)
      {
        Muon *muon = (Muon*) branchMuon->At(i);
        std::cout << "  Muon " << i << ": PT = " << muon->PT << ", Charge = " << muon->Charge << std::endl;
      }
*/
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
    
    // Separate b-jets
    for(int i = 0; i < branchJet->GetEntries(); ++i)
    {
      Jet *jet = (Jet*) branchJet->At(i);
      if(jet->BTag)
      {
        bJetVec.push_back(jet);
      }
    }
    
    // Sort mu+ and mu- by PT
    std::sort(muPlusVec.begin(), muPlusVec.end(), [](Muon* a, Muon* b) { return a->PT > b->PT; });
    std::sort(muMinusVec.begin(), muMinusVec.end(), [](Muon* a, Muon* b) { return a->PT > b->PT; });
    
// Sort b-jets by PT
    std::sort(bJetVec.begin(), bJetVec.end(), [](Jet* a, Jet* b) { return a->PT > b->PT; });


      if(muPlusVec.empty() || muMinusVec.empty() || bJetVec.size() < 2) {
//      std::cout << "No mu+ or mu- found in event " << entry << std::endl;
      continue;
    }
    // Get the highest PT mu+ and mu-
    Muon *muonPlus = muPlusVec[0];
    Muon *muonMinus = muMinusVec[0];
    
    // Get the two highest PT b-jets
    Jet *Bjet1 = bJetVec[0];
    Jet *Bjet2 = bJetVec[1];

//      std::cout << "  Highest PT mu+: PT = " << muonPlus->PT << ", Charge = " << muonPlus->Charge << std::endl;
//      std::cout << "  Highest PT mu-: PT = " << muonMinus->PT << ", Charge = " << muonMinus->Charge << std::endl;

    vecMuP = muonPlus->P4();
    vecMuM = muonMinus->P4();
    vecBJet1 = Bjet1->P4();
    vecBJet2 = Bjet2->P4();
    
    //initial four-momentum of e+e- system 
    TLorentzVector vece1(0.0, 0.0, 125.0, 125.0); // (px, py, pz, E)
    TLorentzVector vece2(0.0, 0, -125.0, 125.0); 
    
    Double_t InMmu12 = (vecMuP+vecMuM).M(); // Invarinat mass of final state mu+ and mu-
    Double_t InMbb = (vecBJet1 + vecBJet2).M();
    Double_t EmuP = vecMuP.E(); //Energy of final mu+
    Double_t EmuM = vecMuM.E(); //Energy of final mu-
    Double_t Eb1 = vecBJet1.E(); //Energy of b1 and b2 jets
    Double_t Eb2 = vecBJet2.E();
    Double_t SqrtS = 2*vece1.E(); //COM energy of e+e- system
    Double_t Mreco = sqrt(pow(SqrtS,2)+pow(InMmu12, 2.0)-2*(EmuP+EmuM)*SqrtS);
    if(Mreco < 120 & InMbb < 100 & EmuP > 5 & EmuM > 5 & Eb1 > 5 & Eb2 > 5)
    {
          hMreco->Fill(Mreco); 
          hInvMuMu->Fill(InMmu12);
          hInvbjj->Fill(InMbb);
          hEMup->Fill(EmuP);
          hEMum->Fill(EmuM);
//        num = num+1;
    }
    }   
//cout << num << endl;
  // Draw histogram
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  //hMreco->Draw();
  //hInvMuMu->Draw();
  //hInvbjj->Draw();
  //hEMup->Draw();
  hEMum->Draw();
  // Save histogram to a file
  /*c1->SaveAs("Recoils_mass_Sig.pdf");
  c1->SaveAs("InvMmumu_Sig.pdf");
  c1->SaveAs("InvMbjj_Sig.pdf");
  c1->SaveAs("EmuP_Sig.pdf");
  c1->SaveAs("EmuM_Sig.pdf");*/
  ///////////////////////////////////////////////
  //c1->SaveAs("Recoils_mass_Back.pdf");
  //c1->SaveAs("InvMmumu_Back.pdf");
  //c1->SaveAs("InvMbjj_Back.pdf");
  //c1->SaveAs("EmuP_Back.pdf");
  c1->SaveAs("EmuM_Back.pdf");
 
 
//Saving histo in a file  
        ofstream outfile;
        //outfile.open("/home/anza-tshilidzi/Downloads/MG5_aMC_v3_5_5/Higgs95_Sig/sig.tsv");
        outfile.open("/home/anza-tshilidzi/Downloads/MG5_aMC_v3_5_5/Z95bb_mumu_back/back.tsv");
  
        for(int i=0; i<100; i++)
        {
        std:outfile << std::fixed << std::setprecision(6) <<70.25+i*0.5<<"\t"<<eventWeight*(hMreco->GetBinContent(i))/numberOfEntries<<endl;   
        }
        outfile.close();
  
  
}


