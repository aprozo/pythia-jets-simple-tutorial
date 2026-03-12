#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector2.h"

#include <vector>
#include <iostream>
#include <cmath>

void anaTree(TString inputFileName="jets_ptHat_10_20_jetR0.4.root") {

  TFile inputFile(inputFileName, "READ");
  if (inputFile.IsZombie()) {
    std::cout << "Cannot open input file " << inputFileName << std::endl;
    return;
  }

  TTree *tree = (TTree*)inputFile.Get("events");
  if (!tree) {
    std::cout << "Cannot find tree 'events'" << std::endl;
    return;
  }

  TFile outFile("anaTrees.root", "RECREATE");

  // event / jet histograms
  TH1D hNJetsPerEvent("hNJetsPerEvent", ";N_{jets};Counts", 12, 0, 12);
  TH1D hJetPt("hJetPt", ";jet p_{T} (GeV/c);Counts", 100, 0, 100);
  TH1D hJetEta("hJetEta", ";jet #eta;Counts", 60, -3, 3);
  TH1D hJetPhi("hJetPhi", ";jet #phi;Counts", 64, -TMath::Pi(), TMath::Pi());

  // dijet histograms
  TH1D hAsymmetry("hAsymmetry", ";A_{J};Counts", 50, 0, 1);
  TH1D hDeltaPhi("hDeltaPhi", ";#Delta#phi;Counts", 64, 0, TMath::Pi());

  // one constituent histogram for dijet events
  TH1D hConstPtDijet("hConstPtDijet", ";constituent p_{T} (GeV/c);Counts", 100, 0, 20);

  // branches
  int nJets = 0;
  int nConstituents = 0;

  std::vector<double> *jetPt = 0;
  std::vector<double> *jetEta = 0;
  std::vector<double> *jetPhi = 0;

  std::vector<int>    *constJetIndex = 0;
  std::vector<double> *constPt = 0;

  tree->SetBranchAddress("nJets", &nJets);
  tree->SetBranchAddress("nConstituents", &nConstituents);

  tree->SetBranchAddress("jetPt",  &jetPt);
  tree->SetBranchAddress("jetEta", &jetEta);
  tree->SetBranchAddress("jetPhi", &jetPhi);

  tree->SetBranchAddress("constituentJetIndex", &constJetIndex);
  tree->SetBranchAddress("constituentPt",       &constPt);

  Long64_t nEntries = tree->GetEntries();

  for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {
    tree->GetEntry(iEvent);

    hNJetsPerEvent.Fill(nJets);

    // single-jet histograms
    for (int iJet = 0; iJet < nJets; ++iJet) {
      hJetPt.Fill(jetPt->at(iJet));
      hJetEta.Fill(jetEta->at(iJet));
      hJetPhi.Fill(jetPhi->at(iJet));
    }

    // dijet analysis - only look at events with exactly 2 jets
    if (nJets == 2) { 
      double pt1 = jetPt->at(0);
      double pt2 = jetPt->at(1);

      double AJ = (pt1 - pt2) / (pt1 + pt2);
      hAsymmetry.Fill(AJ);

      double dphi = std::fabs(TVector2::Phi_mpi_pi(jetPhi->at(0) - jetPhi->at(1)));
      hDeltaPhi.Fill(dphi);

      // fill one constituent histogram using all constituents from both jets
      for (int iConst = 0; iConst < nConstituents; ++iConst) {
        const int jetIndex = constJetIndex->at(iConst);
        if (jetIndex != 0 && jetIndex != 1) continue;

        hConstPtDijet.Fill(constPt->at(iConst));
      }
      
    }
  }

  outFile.Write();
  outFile.Close();
  inputFile.Close();
}