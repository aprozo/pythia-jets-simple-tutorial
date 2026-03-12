// main.cc
#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TString.h"
#include "TVector2.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace Pythia8;

// ---------------------------------------------------------
// One simple container for everything we want to save
// ---------------------------------------------------------
struct MyEvent {
  // counters
  int nJets;
  int nConstituents;

  // jet information
  std::vector<double> jetPt;
  std::vector<double> jetEta;
  std::vector<double> jetPhi;
  std::vector<int>    jetNConstituent;

  // constituent information
  // constituentJetIndex tells us which jet this constituent belongs to
  std::vector<int>    constituentJetIndex;
  std::vector<int>    constituentPdgId;
  std::vector<double> constituentCharge;
  std::vector<double> constituentPx;
  std::vector<double> constituentPy;
  std::vector<double> constituentPz;
  std::vector<double> constituentE;
  std::vector<double> constituentPt;
  std::vector<double> constituentEta;
  std::vector<double> constituentPhi;
  std::vector<double> constituentMass;

  void clear() {
    nJets = 0;
    nConstituents = 0;

    jetPt.clear();
    jetEta.clear();
    jetPhi.clear();
    jetNConstituent.clear();

    constituentJetIndex.clear();
    constituentPdgId.clear();
    constituentCharge.clear();
    constituentPx.clear();
    constituentPy.clear();
    constituentPz.clear();
    constituentE.clear();
    constituentPt.clear();
    constituentEta.clear();
    constituentPhi.clear();
    constituentMass.clear();
  }
};

// ---------------------------------------------------------
// Create ROOT branches
// ---------------------------------------------------------
void makeBranches(TTree* tree, MyEvent& data) {
  tree->Branch("nJets",           &data.nJets,           "nJets/I");
  tree->Branch("nConstituents",   &data.nConstituents,   "nConstituents/I");

  tree->Branch("jetPt",     &data.jetPt);
  tree->Branch("jetEta",    &data.jetEta);
  tree->Branch("jetPhi",    &data.jetPhi);
  tree->Branch("jetNConstituent", &data.jetNConstituent);

  tree->Branch("constituentJetIndex", &data.constituentJetIndex);
  tree->Branch("constituentPdgId",    &data.constituentPdgId);
  tree->Branch("constituentCharge",   &data.constituentCharge);
  tree->Branch("constituentPx",       &data.constituentPx);
  tree->Branch("constituentPy",       &data.constituentPy);
  tree->Branch("constituentPz",       &data.constituentPz);
  tree->Branch("constituentE",        &data.constituentE);
  tree->Branch("constituentPt",       &data.constituentPt);
  tree->Branch("constituentEta",      &data.constituentEta);
  tree->Branch("constituentPhi",      &data.constituentPhi);
  tree->Branch("constituentMass",     &data.constituentMass);
}

int main(int argc, char *argv[]) {

  // -----------------------------
  // Settings
  // -----------------------------
  const double jetRadius      = 0.4;
  const double jetPtMin       = 3.0;
  const double particlePtMin  = 0.15;
  const double particleEtaMax = 1.5;
  const double jetEtaMax      = particleEtaMax - jetRadius;

  // -----------------------------
  // Read command line
  // -----------------------------
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0]
              << " pTHatMin pTHatMax|inf [nEvents=100]\n";
    return 1;
  }

  const double ptHatMin = std::stod(argv[1]);

  double ptHatMax = -1.0;
  std::string s = argv[2];
  if (s != "inf" && s != "Inf" && s != "INF") {
    ptHatMax = std::stod(s);
  }

  if (ptHatMax > 0.0 && ptHatMax < ptHatMin) {
    std::cerr << "[error] pTHatMax < pTHatMin\n";
    return 1;
  }

  const int nEvents = (argc > 3) ? std::atoi(argv[3]) : 100;

  const TString outFile = Form("jets_ptHat_%.0f_%.0f_jetR%.1f.root",
                               ptHatMin,
                               (ptHatMax > 0.0 ? ptHatMax : 999.0),
                               jetRadius);

  // -----------------------------
  // Pythia setup
  // -----------------------------
  Pythia pythia8;
  pythia8.readString("Beams:idA = 2212");
  pythia8.readString("Beams:idB = 2212");
  pythia8.readString("Beams:eCM = 200.");
  pythia8.readString("HardQCD:all = on");

  pythia8.readString(Form("PhaseSpace:pTHatMin = %.1f", ptHatMin));
  if (ptHatMax > 0.0) {
    pythia8.readString(Form("PhaseSpace:pTHatMax = %.1f", ptHatMax));
  }

  if (!pythia8.init()) {
    std::cerr << "[error] PYTHIA init() failed.\n";
    return 2;
  }

  // -----------------------------
  // ROOT output
  // -----------------------------
  TFile *fout = new TFile(outFile, "RECREATE");
  TTree *tree = new TTree("events", "Jet events");

  MyEvent data;
  makeBranches(tree, data);

  // -----------------------------
  // Event loop
  // -----------------------------
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia8.next()) continue;

    data.clear();

    // Build particles for FastJet
    std::vector<fastjet::PseudoJet> particlesForJets;

    for (int iParticle = 0; iParticle < pythia8.event.size(); ++iParticle) {
      const Particle &particle = pythia8.event[iParticle];

      if (!particle.isFinal() || !particle.isVisible() || !particle.isCharged()) continue;
      if (std::abs(particle.eta()) > particleEtaMax) continue;
      if (particle.pT() < particlePtMin) continue;

      fastjet::PseudoJet fjParticle(
        particle.px(), particle.py(), particle.pz(), particle.e()
      );

      fjParticle.set_user_index(iParticle);

      particlesForJets.push_back(fjParticle);
    }

    // If no particles, skip event
    if (particlesForJets.empty())      continue;

    // -----------------------------
    // Run FastJet
    // -----------------------------
    fastjet::JetDefinition jetDefinition(fastjet::antikt_algorithm, jetRadius);
    fastjet::ClusterSequence clusterSequence(particlesForJets, jetDefinition);

    fastjet::Selector selectEta = fastjet::SelectorAbsEtaMax(jetEtaMax);
    fastjet::Selector selectPt  = fastjet::SelectorPtMin(jetPtMin);
    fastjet::Selector selectJets = selectPt && selectEta;

    std::vector<fastjet::PseudoJet> inclusiveJets = clusterSequence.inclusive_jets();
    std::vector<fastjet::PseudoJet> sortedJets    = fastjet::sorted_by_pt(inclusiveJets);
    std::vector<fastjet::PseudoJet> selectedJets  = selectJets(sortedJets);

    // -----------------------------
    // Loop over jets
    // -----------------------------
    for (size_t iJet = 0; iJet < selectedJets.size(); ++iJet) {
      const fastjet::PseudoJet &jet = selectedJets[iJet];

      // save jet info
      data.jetPt.push_back(jet.pt());
      data.jetEta.push_back(jet.eta());
      data.jetPhi.push_back(TVector2::Phi_mpi_pi(jet.phi()));
      data.jetNConstituent.push_back(0); // we will count constituents below

      // get constituents of this jet
      std::vector<fastjet::PseudoJet> constituents =
          fastjet::sorted_by_pt(jet.constituents());

      // loop over constituents of this jet
      for (size_t iConstituent = 0; iConstituent < constituents.size(); ++iConstituent) {
        const fastjet::PseudoJet &c = constituents[iConstituent];

        int pythiaIndex = c.user_index();
        if (pythiaIndex < 0 || pythiaIndex >= pythia8.event.size()) continue;

        const Particle &p = pythia8.event[pythiaIndex];

        data.constituentJetIndex.push_back((int)iJet);
        data.constituentPdgId.push_back(p.id());
        data.constituentCharge.push_back(p.charge());

        data.constituentPx.push_back(p.px());
        data.constituentPy.push_back(p.py());
        data.constituentPz.push_back(p.pz());
        data.constituentE.push_back(p.e());

        data.constituentPt.push_back(p.pT());
        data.constituentEta.push_back(p.eta());
        data.constituentPhi.push_back(TVector2::Phi_mpi_pi(p.phi()));
        data.constituentMass.push_back(p.m());

        data.jetNConstituent.back()++;
      }
    }

    data.nJets = (int)data.jetPt.size();
    data.nConstituents = (int)data.constituentPt.size();

    // -----------------------------
    // Print event content
    // -----------------------------
    std::cout << "=============================================\n";
    std::cout << "Event " << iEvent
              << "  nJets = " << data.nJets
              << "  nConstituents = " << data.nConstituents
              << "\n";

    for (int j = 0; j < data.nJets; ++j) {
      std::cout << "  Jet " << j
                << "  pt="   << data.jetPt[j]
                << "  eta="  << data.jetEta[j]
                << "  phi="  << data.jetPhi[j]
                << "  nConstituent=" << data.jetNConstituent[j]
                << "\n";
    }

    tree->Fill();
  }

  // -----------------------------
  // Save generator information
  // -----------------------------
  const double sigmaGen = pythia8.info.sigmaGen();
  const double sigmaErr = pythia8.info.sigmaErr();

  TH1D *stats = new TH1D("stats", "stats", 4, 0, 4);
  std::vector<TString> statNames;
  statNames.push_back("nEvents");
  statNames.push_back("sigmaGen_mb");
  statNames.push_back("ptHatMin");
  statNames.push_back("ptHatMax");

  for (size_t i = 0; i < statNames.size(); ++i) {
    stats->GetXaxis()->SetBinLabel(i + 1, statNames[i]);
  }

  stats->SetBinContent(1, nEvents);
  stats->SetBinContent(2, sigmaGen);
  stats->SetBinError(2, sigmaErr);
  stats->SetBinContent(3, ptHatMin);
  stats->SetBinContent(4, ptHatMax);

  std::cout << "[done] Wrote " << outFile << "\n"
            << "       sigmaGen = " << sigmaGen << " mb (± " << sigmaErr << ")\n"
            << "       ptHatMin = " << ptHatMin << "\n"
            << "       ptHatMax = " << ptHatMax << "\n";

  pythia8.stat();

  fout->Write();
  fout->Close();
  delete fout;

  return 0;
}