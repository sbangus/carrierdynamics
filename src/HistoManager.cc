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

#include <string>

namespace {

// Decode a LET H1 ID (kMaxHisto .. kMaxHistoLet-1) into a descriptive title.
G4String LetH1Title(G4int k)
{
  auto perCat = [](const G4String& base, G4int relBase) {
    G4int rel = relBase;
    G4int a = rel / kNbDepositCategories + 1;
    G4int c = rel % kNbDepositCategories;
    return base + ", absorber " + std::to_string(a) + ", " + DepositCategoryName(c);
  };

  if (k >= kLetStepCountBase && k < kLetStepEWeightedBase)
    return perCat("LETdep step-count spectrum", k - kLetStepCountBase);
  if (k >= kLetStepEWeightedBase && k < kLetCalcEWeightedBase)
    return perCat("LETdep Eion-weighted spectrum", k - kLetStepEWeightedBase);
  if (k >= kLetCalcEWeightedBase && k < kTrackEionCatBase)
    return perCat("LETcalc Eion-weighted spectrum", k - kLetCalcEWeightedBase);
  if (k >= kTrackEionCatBase && k < kTrackLetCatBase)
    return perCat("Track Eion", k - kTrackEionCatBase);
  if (k >= kTrackLetCatBase && k < kTrackDepthSpanCatBase)
    return perCat("Track LET", k - kTrackLetCatBase);
  if (k >= kTrackDepthSpanCatBase && k < kEionEventAllBase)
    return perCat("Track depth span", k - kTrackDepthSpanCatBase);

  auto perAbs = [](const G4String& base, G4int rel) {
    return base + ", absorber " + std::to_string(rel + 1);
  };
  if (k >= kEionEventAllBase && k < kEionEventHitBase)
    return perAbs("Eion per event (all events, incl. zero)", k - kEionEventAllBase);
  if (k >= kEionEventHitBase && k < kNielEventBase)
    return perAbs("Eion per event (hit-only)", k - kEionEventHitBase);
  if (k >= kNielEventBase && k < kLetTEventBase)
    return perAbs("NIEL per event", k - kNielEventBase);
  if (k >= kLetTEventBase && k < kLetDdepEventBase)
    return perAbs("Event track-length-averaged LET", k - kLetTEventBase);
  if (k >= kLetDdepEventBase && k < kLetDcalcEventBase)
    return perAbs("Event Eion-weighted Ldep", k - kLetDdepEventBase);
  if (k >= kLetDcalcEventBase && k < kLetMaxEventBase)
    return perAbs("Event Eion-weighted Lcalc", k - kLetDcalcEventBase);
  if (k >= kLetMaxEventBase && k < kFracAbove1Base)
    return perAbs("Event maximum step LET", k - kLetMaxEventBase);
  if (k >= kFracAbove1Base && k < kFracAbove2Base)
    return perAbs("Eion fraction Lcalc>=10 keV/um", k - kFracAbove1Base);
  if (k >= kFracAbove2Base && k < kFracAbove3Base)
    return perAbs("Eion fraction Lcalc>=100 keV/um", k - kFracAbove2Base);
  if (k >= kFracAbove3Base && k < kEionDepthBase)
    return perAbs("Eion fraction Lcalc>=1000 keV/um", k - kFracAbove3Base);
  if (k >= kEionDepthBase && k < kEionNormDepthBase)
    return perAbs("dEion/dz absolute depth profile", k - kEionDepthBase);
  if (k >= kEionNormDepthBase && k < kTrackContainBase)
    return perAbs("dEion/du normalized depth profile", k - kEionNormDepthBase);
  if (k >= kTrackContainBase && k < kPrimaryLetStepCountBase)
    return perAbs("Track containment (0=unk,1=cont,2=front,3=back,4=lat,5=multi)",
                  k - kTrackContainBase);
  if (k >= kPrimaryLetStepCountBase && k < kPrimaryLetStepEWeightedBase)
    return perAbs("Primary LETdep step-count spectrum", k - kPrimaryLetStepCountBase);
  if (k >= kPrimaryLetStepEWeightedBase && k < kPrimaryLetCalcEWeightedBase)
    return perAbs("Primary LETdep Eion-weighted spectrum",
                  k - kPrimaryLetStepEWeightedBase);
  if (k >= kPrimaryLetCalcEWeightedBase && k < kPrimaryTrackEionBase)
    return perAbs("Primary LETcalc Eion-weighted spectrum",
                  k - kPrimaryLetCalcEWeightedBase);
  if (k >= kPrimaryTrackEionBase && k < kPrimaryTrackLetBase)
    return perAbs("Primary track Eion", k - kPrimaryTrackEionBase);
  if (k >= kPrimaryTrackLetBase && k < kPrimaryTrackDepthSpanBase)
    return perAbs("Primary track LET", k - kPrimaryTrackLetBase);
  if (k >= kPrimaryTrackDepthSpanBase && k < kMaxHistoLet)
    return perAbs("Primary track depth span", k - kPrimaryTrackDepthSpanBase);

  return "letH" + std::to_string(k);
}

// H2 family index (0..7) and absorber -> descriptive name / title.
const char* H2FamilyTag(G4int fam)
{
  static const char* tags[8] = {
      "evtEionVsLetCalc", "evtEionVsLetDep", "letCalcVsDepth", "letDepVsDepth",
      "letCalcVsKE",      "trkEionVsLetCalc", "trkEionVsLetDep", "depthVsTransverse"};
  return (fam >= 0 && fam < 8) ? tags[fam] : "h2";
}

}  // namespace

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
    else if (k >= kEdepByOriginBase && k < kEdepByOriginBase + kMaxAbsor * kNbOrigin) {
      G4int relIdx = k - kEdepByOriginBase;
      G4int depositAbs = relIdx / kNbOrigin + 1;
      G4int originIdx = relIdx % kNbOrigin;
      if (originIdx == 0) {
        title = "Edep in absorber " + std::to_string(depositAbs)
                + " from external/primary lineage";
      }
      else {
        title = "Edep in absorber " + std::to_string(depositAbs)
                + " from lineage originating in absorber " + std::to_string(originIdx);
      }
    }
    G4int ih = analysisManager->CreateH1(id, title, nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }

  // -------------------------------------------------------------------------
  // LET / track-structure H1 histograms (IDs kMaxHisto .. kMaxHistoLet-1).
  // Created inactive with placeholder binning; Run::Run() reconfigures binning
  // and activates the ones belonging to charge-active absorbers.
  // -------------------------------------------------------------------------
  for (G4int k = kMaxHisto; k < kMaxHistoLet; ++k) {
    G4int ih = analysisManager->CreateH1(std::to_string(k), LetH1Title(k),
                                         nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, false);
  }

  // -------------------------------------------------------------------------
  // LET H2 histograms (independent H2 ID space). Created in family-major,
  // absorber-minor order so the sequential H2 IDs match the H2XxxId() helpers.
  // Placeholder binning; Run::Run() reconfigures/activates active absorbers.
  // -------------------------------------------------------------------------
  for (G4int fam = 0; fam < 8; ++fam) {
    for (G4int a = 1; a <= kMaxAbsor; ++a) {
      const G4String name =
          G4String(H2FamilyTag(fam)) + "_abs" + std::to_string(a);
      G4int ih = analysisManager->CreateH2(name, name, 10, 0., 1., 10, 0., 1.);
      analysisManager->SetH2Activation(ih, false);
    }
  }

  // -------------------------------------------------------------------------
  // LET ROOT ntuples. Enable merging BEFORE creating them (and before the file
  // is opened in RunAction::BeginOfRunAction) so per-thread ntuples are merged
  // into a single set in the master file. Column order here MUST match the fill
  // order in EventAction / TrackingAction / SteppingAction.
  // -------------------------------------------------------------------------
  // Column-wise (the default) produces one named ROOT branch per column, which
  // is what ROOT/uproot analysis expects; row-wise packs everything into a
  // single opaque branch that downstream tools cannot split by column.
  analysisManager->SetNtupleMerging(true);

  // EventLET (id 0): one row per event per charge-active absorber.
  analysisManager->CreateNtuple("EventLET", "Event-level LET scoring");
  analysisManager->CreateNtupleIColumn("runID");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->CreateNtupleIColumn("absorberID");
  analysisManager->CreateNtupleIColumn("primaryPDG");
  analysisManager->CreateNtupleDColumn("primaryEnergy_MeV");
  analysisManager->CreateNtupleDColumn("trackWeight");
  analysisManager->CreateNtupleDColumn("edep_keV");
  analysisManager->CreateNtupleDColumn("eion_keV");
  analysisManager->CreateNtupleDColumn("niel_keV");
  analysisManager->CreateNtupleDColumn("chargedPath_um");
  analysisManager->CreateNtupleDColumn("letT_keV_per_um");
  analysisManager->CreateNtupleDColumn("letDdep_keV_per_um");
  analysisManager->CreateNtupleDColumn("letDcalc_keV_per_um");
  analysisManager->CreateNtupleDColumn("maxLetDep_keV_per_um");
  analysisManager->CreateNtupleDColumn("maxLetCalc_keV_per_um");
  analysisManager->CreateNtupleDColumn("fractionEionAbove10");
  analysisManager->CreateNtupleDColumn("fractionEionAbove100");
  analysisManager->CreateNtupleDColumn("fractionEionAbove1000");
  analysisManager->CreateNtupleDColumn("meanDepth_um");
  analysisManager->CreateNtupleDColumn("depthSigma_um");
  analysisManager->CreateNtupleIColumn("nDepositingSteps");
  analysisManager->CreateNtupleIColumn("nChargedTracks");
  analysisManager->CreateNtupleDColumn("initialPairs");
  analysisManager->FinishNtuple();

  // TrackLET (id 1): one row per depositing track per charge-active absorber.
  analysisManager->CreateNtuple("TrackLET", "Track-level LET scoring");
  analysisManager->CreateNtupleIColumn("runID");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->CreateNtupleIColumn("trackID");
  analysisManager->CreateNtupleIColumn("parentID");
  analysisManager->CreateNtupleIColumn("absorberID");
  analysisManager->CreateNtupleIColumn("PDG");
  analysisManager->CreateNtupleDColumn("charge");
  analysisManager->CreateNtupleIColumn("atomicZ");
  analysisManager->CreateNtupleIColumn("atomicA");
  analysisManager->CreateNtupleIColumn("depositCategory");
  analysisManager->CreateNtupleIColumn("reactionCategory");
  analysisManager->CreateNtupleIColumn("electronLineage");
  analysisManager->CreateNtupleIColumn("originAbsorber");
  analysisManager->CreateNtupleDColumn("vertexKE_MeV");
  analysisManager->CreateNtupleDColumn("entryKE_MeV");
  analysisManager->CreateNtupleDColumn("exitKE_MeV");
  analysisManager->CreateNtupleDColumn("edep_keV");
  analysisManager->CreateNtupleDColumn("eion_keV");
  analysisManager->CreateNtupleDColumn("niel_keV");
  analysisManager->CreateNtupleDColumn("path_um");
  analysisManager->CreateNtupleDColumn("trackLET_keV_per_um");
  analysisManager->CreateNtupleDColumn("letDdep_keV_per_um");
  analysisManager->CreateNtupleDColumn("letDcalc_keV_per_um");
  analysisManager->CreateNtupleDColumn("maxLetDep_keV_per_um");
  analysisManager->CreateNtupleDColumn("maxLetCalc_keV_per_um");
  analysisManager->CreateNtupleDColumn("startDepth_um");
  analysisManager->CreateNtupleDColumn("endDepth_um");
  analysisManager->CreateNtupleDColumn("minDepth_um");
  analysisManager->CreateNtupleDColumn("maxDepth_um");
  analysisManager->CreateNtupleDColumn("meanDepth_um");
  analysisManager->CreateNtupleIColumn("escapeCategory");
  analysisManager->CreateNtupleIColumn("stoppedInActiveLayer");
  analysisManager->CreateNtupleIColumn("nDepositingSteps");
  analysisManager->CreateNtupleDColumn("trackWeight");
  analysisManager->FinishNtuple();

  // StepLET (id 2): sampled/gated per-step diagnostic ntuple.
  analysisManager->CreateNtuple("StepLET", "Sampled step-level LET scoring");
  analysisManager->CreateNtupleIColumn("runID");
  analysisManager->CreateNtupleIColumn("eventID");
  analysisManager->CreateNtupleIColumn("trackID");
  analysisManager->CreateNtupleIColumn("parentID");
  analysisManager->CreateNtupleIColumn("absorberID");
  analysisManager->CreateNtupleIColumn("PDG");
  analysisManager->CreateNtupleIColumn("atomicZ");
  analysisManager->CreateNtupleIColumn("atomicA");
  analysisManager->CreateNtupleIColumn("depositCategory");
  analysisManager->CreateNtupleIColumn("reactionCategory");
  analysisManager->CreateNtupleDColumn("preKE_MeV");
  analysisManager->CreateNtupleDColumn("postKE_MeV");
  analysisManager->CreateNtupleDColumn("midpointKE_MeV");
  analysisManager->CreateNtupleDColumn("edep_keV");
  analysisManager->CreateNtupleDColumn("eion_keV");
  analysisManager->CreateNtupleDColumn("niel_keV");
  analysisManager->CreateNtupleDColumn("stepLength_um");
  analysisManager->CreateNtupleDColumn("letDep_keV_per_um");
  analysisManager->CreateNtupleDColumn("letCalc_keV_per_um");
  analysisManager->CreateNtupleDColumn("localX_um");
  analysisManager->CreateNtupleDColumn("localY_um");
  analysisManager->CreateNtupleDColumn("localZ_um");
  analysisManager->CreateNtupleDColumn("normalizedDepth");
  analysisManager->CreateNtupleDColumn("globalTime_ns");
  analysisManager->CreateNtupleIColumn("processSubType");
  analysisManager->CreateNtupleDColumn("weight");
  analysisManager->FinishNtuple();
}
