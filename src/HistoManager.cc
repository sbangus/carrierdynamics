//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file HistoManager.cc
/// \brief Implementation of the HistoManager class

#include "HistoManager.hh"

#include "DetectorConstruction.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);  // enable inactivation of histograms

  // Default values (to be reset via /analysis/h1/set command or auto-configured)
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  const G4String procCatNames[] = {"hadElastic", "nInelastic", "Compton scattering",
                                   "photoelectric absorption", "pair production"};

  // Create all histograms as inactivated. Run::Run() reconfigures binning and
  // activates the subset that is meaningful for the current geometry.
  for (G4int k = 0; k < kMaxHisto; k++) {
    G4String id = std::to_string(k);
    G4String title = "h" + id;  // generic placeholder

    if (k >= 1 && k <= kMaxAbsor) {
      title = "Edep in absorber " + std::to_string(k);
    }
    else if (k > kMaxAbsor && k <= 2 * kMaxAbsor) {
      title = "Edep longit. profile (MeV/event) in absorber "
              + std::to_string(k - kMaxAbsor);
    }
    else if (k == 2 * kMaxAbsor + 1) {
      title = "energy flow (MeV/event)";
    }
    else if (k == 2 * kMaxAbsor + 2) {
      title = "total energy deposited";
    }
    else if (k == 2 * kMaxAbsor + 3) {
      title = "total energy leakage";
    }
    else if (k == 2 * kMaxAbsor + 4) {
      title = "total energy released : Edep + Eleak";
    }
    else if (k >= kDepthHistoBase && k < kDepthHistoEnd) {
      G4int relIdx = k - kDepthHistoBase;
      G4int absNum = relIdx / kNbProcessCat + 1;
      G4int procCat = relIdx % kNbProcessCat;
      title = "Depth of " + procCatNames[procCat];
      if (procCat == kProcCatNInelastic) title += " (incl. nCapture)";
      title += " in absorber " + std::to_string(absNum);
    }
    else if (k >= kEdepByParticleBase && k < kCaptureGammaBase) {
      G4int relIdx = k - kEdepByParticleBase;
      G4int absNum = relIdx / kNbFixedParticles + 1;
      G4int partIdx = relIdx % kNbFixedParticles;
      title = "Edep by " + FixedParticleName(partIdx)
              + " in absorber " + std::to_string(absNum);
    }
    else if (k >= kCaptureGammaBase && k < kEhpBase) {
      G4int absNum = k - kCaptureGammaBase + 1;
      title = "Capture-gamma KE spectrum in absorber " + std::to_string(absNum);
    }
    else if (k >= kEhpBase && k < kSecDepthBase) {
      G4int absNum = k - kEhpBase + 1;
      title = "EHP per event (Edep / 6 eV) in absorber " + std::to_string(absNum);
    }
    else if (k >= kSecDepthBase && k < kSecESpecBase) {
      G4int absNum = k - kSecDepthBase + 1;
      title = "Secondary production depth in absorber " + std::to_string(absNum);
    }
    else if (k >= kSecESpecBase && k < kNeutronEnergyBase) {
      G4int partIdx = k - kSecESpecBase;
      title = "KE spectrum of secondary " + FixedParticleName(partIdx);
    }
    else if (k >= kNeutronEnergyBase && k < kInteractionFractionId) {
      G4int absNum = k - kNeutronEnergyBase + 1;
      title = "Neutron KE at interaction in absorber " + std::to_string(absNum);
    }
    else if (k == kInteractionFractionId) {
      title = "Fraction of incident neutrons that interact (per absorber)";
    }
    else if (k == kComparativeResponseId) {
      title = "Mean Edep per event vs. absorber index";
    }
    else if (k >= kFirstGenSecESpecBase &&
             k < kFirstGenSecESpecBase + kNbFixedParticles) {
      G4int partIdx = k - kFirstGenSecESpecBase;
      title = "KE at generation of first-generation secondary "
              + FixedParticleName(partIdx);
    }
    else if (k == kDetectorFrontKEAll) {
      title = "Neutron KE at detector front face (all neutrons, after environment)";
    }
    else if (k == kDetectorFrontKEPrimary) {
      title = "Neutron KE at detector front face (primary neutrons only)";
    }
    else if (k == kDetectorFrontGammaKE) {
      title = "Gamma KE spectrum at detector front face (just upstream of absorber 1)";
    }
    else if (k >= kEdepElectronGammaBase && k < kEdepElectronIonicBase) {
      G4int absNum = k - kEdepElectronGammaBase + 1;
      title = "Edep by e- (gamma-mediated) in absorber " + std::to_string(absNum);
    }
    else if (k >= kEdepElectronIonicBase && k < kModeratorCaptureGammaKE) {
      G4int absNum = k - kEdepElectronIonicBase + 1;
      title = "Edep by e- (ion/hadronic-mediated) in absorber " + std::to_string(absNum);
    }
    else if (k == kModeratorCaptureGammaKE) {
      title = "Capture-gamma KE at birth in moderator (0-3 MeV line spectrum)";
    }
    else if (k >= kModeratorSecESpecBase && k < kModeratorFirstGenSecESpecBase) {
      G4int partIdx = k - kModeratorSecESpecBase;
      title = "KE spectrum of moderator-born secondary "
              + FixedParticleName(partIdx);
    }
    else if (k >= kModeratorFirstGenSecESpecBase &&
             k < kModeratorFirstGenSecESpecBase + kNbFixedParticles) {
      G4int partIdx = k - kModeratorFirstGenSecESpecBase;
      title = "KE at generation of first-generation moderator secondary "
              + FixedParticleName(partIdx);
    }
    else if (k >= kTrackPathLengthBase && k < kTotalTrackLengthBase) {
      G4int relIdx = k - kTrackPathLengthBase;
      G4int absNum = relIdx / kNbFixedParticles + 1;
      G4int partIdx = relIdx % kNbFixedParticles;
      title = "Path length of " + FixedParticleName(partIdx)
              + " in absorber " + std::to_string(absNum);
    }
    else if (k >= kTotalTrackLengthBase && k < kPathElectronGammaBase) {
      G4int partIdx = k - kTotalTrackLengthBase;
      title = "Total track length of " + FixedParticleName(partIdx);
    }
    else if (k >= kPathElectronGammaBase && k < kPathElectronIonicBase) {
      title = "Path length of e- (gamma-mediated) in absorber "
              + std::to_string(k - kPathElectronGammaBase + 1);
    }
    else if (k >= kPathElectronIonicBase && k < kTotalPathElectronGamma) {
      title = "Path length of e- (ion/hadronic-mediated) in absorber "
              + std::to_string(k - kPathElectronIonicBase + 1);
    }
    else if (k == kTotalPathElectronGamma) {
      title = "Total track length of e- (gamma-mediated)";
    }
    else if (k == kTotalPathElectronIonic) {
      title = "Total track length of e- (ion/hadronic-mediated)";
    }
    else if (k == kSecESpecElectronGamma) {
      title = "KE spectrum of secondary e- (gamma-mediated)";
    }
    else if (k == kSecESpecElectronIonic) {
      title = "KE spectrum of secondary e- (ion/hadronic-mediated)";
    }
    else if (k == kFirstGenSecESpecElectronGamma) {
      title = "KE at generation of first-generation secondary e- (gamma-mediated)";
    }
    else if (k == kFirstGenSecESpecElectronIonic) {
      title = "KE at generation of first-generation secondary e- (ion/hadronic-mediated)";
    }
    else if (k >= kPrimaryPathLengthBase && k < kTotalPrimaryPathLength) {
      G4int absNum = k - kPrimaryPathLengthBase + 1;
      title = "Path length of primary in absorber " + std::to_string(absNum);
    }
    else if (k == kTotalPrimaryPathLength) {
      title = "Total track length of primary";
    }
    else if (k >= kPrimaryEdepBase && k < kPrimaryEdepBase + kMaxAbsor) {
      G4int absNum = k - kPrimaryEdepBase + 1;
      title = "Edep by primary in absorber " + std::to_string(absNum);
    }
    G4int ih = analysisManager->CreateH1(id, title, nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }
}
