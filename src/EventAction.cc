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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"

#include "DetectorConstruction.hh"
#include "HistoManager.hh"
#include "Run.hh"

#include "G4Event.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

#include <algorithm>
#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* det) : fDetector(det) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  // initialize EnergyDeposit per event
  //
  for (G4int k = 0; k < kMaxAbsor; k++) {
    fEnergyDeposit[k] = fTrackLengthCh[k] = 0.0;
    fPrimaryEdep[k] = 0.0;
    fEdepByParticle[k].clear();
    fNeutronInteracted[k] = false;
    fEdepElectronGamma[k] = 0.;
    fEdepElectronIonic[k] = 0.;
    for (G4int o = 0; o < kNbOrigin; ++o) fEdepByOrigin[k][o] = 0.;
  }

  // initialize EnergyLeakage per event
  //
  fEnergyLeak[0] = fEnergyLeak[1] = 0.0;
  fMonitorCrossedTrackIds.clear();
  fElectronLineage.clear();
  fTrackParticleName.clear();
  fTrackOriginAbsor.clear();

  // Reset LET / track-structure scoring state.
  for (G4int k = 0; k < kMaxAbsor; k++) {
    fEventLet[k] = EventLetScore{};
  }
  fTrackLetInAbsor.clear();
  fTrackReaction.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SumEnergy(G4int k, G4double de, G4double dl)
{
  fEnergyDeposit[k] += de;
  fTrackLengthCh[k] += dl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SumPrimaryEnergy(G4int k, G4double de)
{
  if (k >= 1 && k <= kMaxAbsor && de > 0.) fPrimaryEdep[k - 1] += de;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SumEnergyByParticle(G4int k, const G4String& particleName, G4double de)
{
  fEdepByParticle[k][particleName] += de;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SumEnergyLeak(G4double eleak, G4int index)
{
  fEnergyLeak[index] += eleak;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::MarkNeutronInteraction(G4int k)
{
  if (k >= 0 && k < kMaxAbsor) fNeutronInteracted[k] = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool EventAction::HasMonitorCrossed(G4int trackId) const
{
  return (fMonitorCrossedTrackIds.find(trackId) != fMonitorCrossedTrackIds.end());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::MarkMonitorCrossed(G4int trackId)
{
  fMonitorCrossedTrackIds.insert(trackId);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::TagElectronLineage(G4int trackId, G4int lineage)
{
  fElectronLineage[trackId] = lineage;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int EventAction::GetElectronLineage(G4int trackId) const
{
  const auto it = fElectronLineage.find(trackId);
  if (it == fElectronLineage.end()) return kElectronLineageNone;
  return it->second;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::RegisterTrackParticle(G4int trackId, const G4String& particleName)
{
  fTrackParticleName[trackId] = particleName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String EventAction::GetTrackParticleName(G4int trackId) const
{
  const auto it = fTrackParticleName.find(trackId);
  if (it == fTrackParticleName.end()) return "";
  return it->second;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SetTrackOrigin(G4int trackId, G4int originAbsor)
{
  fTrackOriginAbsor[trackId] = originAbsor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int EventAction::GetTrackOrigin(G4int trackId) const
{
  const auto it = fTrackOriginAbsor.find(trackId);
  if (it == fTrackOriginAbsor.end()) return 0;  // external / primary lineage
  return it->second;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SumEdepByOrigin(G4int depositAbsor, G4int originAbsor, G4double de)
{
  if (depositAbsor < 1 || depositAbsor > kMaxAbsor) return;
  if (originAbsor < 0 || originAbsor >= kNbOrigin) return;
  if (de <= 0.) return;
  fEdepByOrigin[depositAbsor - 1][originAbsor] += de;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddTrackPathLength(G4int absorNum, G4double len)
{
  if (absorNum < 1 || absorNum > kMaxAbsor || len <= 0.) return;
  fTrackPathInAbsor[absorNum] += len;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SumElectronEdepByLineage(G4int absorNum, G4int lineage, G4double de)
{
  if (absorNum < 1 || absorNum > kMaxAbsor || de <= 0.) return;
  const G4int k = absorNum - 1;
  if (lineage == kElectronLineageGamma) {
    fEdepElectronGamma[k] += de;
  }
  else if (lineage == kElectronLineageIon) {
    fEdepElectronIonic[k] += de;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AccumulateLetStep(G4int absorNum, const LetStepData& s)
{
  if (absorNum < 1 || absorNum >= kMaxAbsor) return;
  EventLetScore& ev = fEventLet[absorNum];
  TrackLetScore& trk = fTrackLetInAbsor[absorNum];

  ev.eDep += s.edep;
  ev.eIon += s.eIon;
  ev.niel += s.niel;
  trk.eDep += s.edep;
  trk.eIon += s.eIon;
  trk.niel += s.niel;

  if (s.charged) {
    ev.chargedPath += s.dl;
    trk.chargedPath += s.dl;
  }

  // Track entry/exit kinematics and boundary crossings (containment).
  if (s.enteredAbsor && trk.entryKE < 0.) trk.entryKE = s.preKE;
  trk.exitKE = s.postKE;
  if (s.exitedFront) trk.crossedFront = true;
  if (s.exitedBack) trk.crossedBack = true;
  if (s.exitedLateral) trk.crossedLateral = true;

  if (s.eIon > 0.) {
    const G4int c = static_cast<G4int>(s.depCat);

    ev.sumEionLetDep += s.eIon * s.letDep;
    ev.sumEionLetCalc += s.eIon * s.letCalc;
    ev.maxLetDep = std::max(ev.maxLetDep, s.letDep);
    ev.maxLetCalc = std::max(ev.maxLetCalc, s.letCalc);
    ev.sumEionDepth += s.eIon * s.depth;
    ev.sumEionDepth2 += s.eIon * s.depth * s.depth;
    ev.eIonByCategory[c] += s.eIon;
    ++ev.nDepositingSteps;
    if (s.letCalc >= kLetThreshold1) ev.eIonAbove1 += s.eIon;
    if (s.letCalc >= kLetThreshold2) ev.eIonAbove2 += s.eIon;
    if (s.letCalc >= kLetThreshold3) ev.eIonAbove3 += s.eIon;

    trk.sumEionLetDep += s.eIon * s.letDep;
    trk.sumEionLetCalc += s.eIon * s.letCalc;
    trk.maxLetDep = std::max(trk.maxLetDep, s.letDep);
    trk.maxLetCalc = std::max(trk.maxLetCalc, s.letCalc);
    trk.sumEionDepth += s.eIon * s.depth;
    trk.minDepth = std::min(trk.minDepth, s.depth);
    trk.maxDepth = std::max(trk.maxDepth, s.depth);
    if (!trk.hasStart) {
      trk.startDepth = s.depth;
      trk.hasStart = true;
    }
    trk.endDepth = s.depth;
    ++trk.nDepositingSteps;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::CountCompletedChargedTrack(G4int absorNum)
{
  if (absorNum >= 1 && absorNum < kMaxAbsor) ++fEventLet[absorNum].nChargedTracks;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::SetTrackReaction(G4int trackId, ReactionCategory category)
{
  fTrackReaction[trackId] = category;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ReactionCategory EventAction::GetTrackReaction(G4int trackId) const
{
  const auto it = fTrackReaction.find(trackId);
  return (it == fTrackReaction.end()) ? ReactionCategory::None : it->second;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DepositCategory EventAction::ClassifyDeposit(const G4Track* track) const
{
  const G4ParticleDefinition* p = track->GetDefinition();
  const G4int pdg = p->GetPDGEncoding();

  if (pdg == 11) {
    const G4int lineage = GetElectronLineage(track->GetTrackID());
    return (lineage == kElectronLineageGamma) ? DepositCategory::ElectronGamma
                                              : DepositCategory::ElectronIon;
  }
  if (pdg == -11) return DepositCategory::Positron;

  const G4String& name = p->GetParticleName();
  if (name == "proton") return DepositCategory::Proton;
  if (name == "deuteron" || name == "triton" || name == "He3")
    return DepositCategory::LightIon;
  if (name == "alpha") return DepositCategory::Alpha;

  const G4double q = p->GetPDGCharge();
  const G4int Z = p->GetAtomicNumber();
  if (Z == 3) return DepositCategory::Lithium;
  if (Z == 6) return DepositCategory::Carbon;
  if (Z == 4 || Z == 5) return DepositCategory::BeB;
  if (q != 0. && Z > 0) return DepositCategory::OtherHeavyIon;
  if (q == 0.) return DepositCategory::NeutralLocal;
  return DepositCategory::Other;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
  // get Run
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();

  G4double EdepTot = 0.;
  const G4int nAbs = fDetector->GetNbOfAbsor();
  for (G4int k = 1; k <= nAbs; k++) {
    run->SumEdepPerAbsorber(k, fEnergyDeposit[k], fTrackLengthCh[k]);
    run->SumEdepByParticle(k, fEdepByParticle[k]);
    if (fEnergyDeposit[k] > 0.) analysis->FillH1(k, fEnergyDeposit[k]);
    EdepTot += fEnergyDeposit[k];

    // Per-particle Edep histograms (per event, per absorber). Only the fixed
    // particles tracked by FixedParticleIdx() get histogrammed; the rest
    // remain in the run's text summary via fEdepByParticle.
    for (const auto& entry : fEdepByParticle[k]) {
      G4int pi = FixedParticleIdx(entry.first);
      if (pi >= 0 && entry.second > 0.) {
        analysis->FillH1(EdepByParticleId(k, pi), entry.second);
      }
    }

    // Estimated electron-hole pairs per event in PNDI-T10 / PNDI-T10-B4C
    // using a fixed pair-creation energy of 6 eV. Filled only for PNDI*
    // absorbers; the histogram is inactive elsewhere so a stray fill is
    // harmless.
    G4String matName = fDetector->GetAbsorMaterial(k)->GetName();
    if (fEnergyDeposit[k] > 0. && matName.find("PNDI") != std::string::npos) {
      G4double ehp = fEnergyDeposit[k] / (6.0 * eV);
      analysis->FillH1(EhpId(k), ehp);
    }

    const G4int idx = k - 1;
    if (fEdepElectronGamma[idx] > 0.) {
      analysis->FillH1(EdepElectronGammaId(k), fEdepElectronGamma[idx]);
    }
    if (fEdepElectronIonic[idx] > 0.) {
      analysis->FillH1(EdepElectronIonicId(k), fEdepElectronIonic[idx]);
    }
    if (fPrimaryEdep[idx] > 0.) {
      analysis->FillH1(PrimaryEdepId(k), fPrimaryEdep[idx]);
    }

    // Cross-absorber provenance: Edep in absorber k split by origin absorber
    // (0 = external/primary, 1..nAbs = born in that absorber's lineage).
    for (G4int o = 0; o <= nAbs; ++o) {
      const G4double de = fEdepByOrigin[idx][o];
      if (de > 0.) {
        analysis->FillH1(EdepByOriginId(k, o), de);
        run->SumEdepByOrigin(k, o, de);
      }
    }
  }
  run->SumEnergies(EdepTot, fEnergyLeak[0], fEnergyLeak[1]);
  run->AddInteractedFlags(fNeutronInteracted);

  G4double EleakTot = fEnergyLeak[0] + fEnergyLeak[1];
  G4double ETot = EdepTot + EleakTot;
  G4int id = 2 * kMaxAbsor + 1;
  analysis->FillH1(++id, EdepTot);
  analysis->FillH1(++id, EleakTot);
  analysis->FillH1(++id, ETot);

  // ---- LET / track-structure: event-level histograms + EventLET ntuple ----
  // Write one EventLET row per charge-active absorber, including events with
  // zero deposition (the all-events Eion spectrum answers "energy deposited per
  // incident primary" and must include zeros). The units stored are explicit
  // human units (keV, um, keV/um, MeV) matching the ntuple column names.
  const G4double keV_um = keV / um;

  // Actual sampled primary for this event (first primary of the first vertex).
  G4int primaryPDG = 0;
  G4double primaryEnergy = 0.;
  G4double primaryWeight = 1.;
  if (evt->GetNumberOfPrimaryVertex() > 0) {
    const G4PrimaryVertex* pv = evt->GetPrimaryVertex(0);
    if (pv != nullptr && pv->GetNumberOfParticle() > 0) {
      const G4PrimaryParticle* pp = pv->GetPrimary(0);
      if (pp != nullptr) {
        primaryPDG = pp->GetPDGcode();
        primaryEnergy = pp->GetKineticEnergy();
        primaryWeight = pp->GetWeight();
      }
    }
  }

  const G4int runID =
      G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
  const G4int eventID = evt->GetEventID();

  for (G4int k = 1; k <= nAbs; ++k) {
    if (!fDetector->IsChargeActiveAbsorber(k)) continue;
    const EventLetScore& sc = fEventLet[k];

    // Run-level energy-partition totals (NIEL/Edep provenance + regression guard).
    run->AccumulateLetTotals(k, sc.eDep, sc.eIon, sc.niel, sc.eIonByCategory);

    // Derived event quantities (guarded divisions).
    const G4double letT =
        (sc.chargedPath > 0.) ? sc.eIon / sc.chargedPath : 0.;
    const G4double letDdep = (sc.eIon > 0.) ? sc.sumEionLetDep / sc.eIon : 0.;
    const G4double letDcalc = (sc.eIon > 0.) ? sc.sumEionLetCalc / sc.eIon : 0.;
    const G4double meanDepth = (sc.eIon > 0.) ? sc.sumEionDepth / sc.eIon : 0.;
    const G4double depthVar =
        (sc.eIon > 0.)
            ? std::max(0., sc.sumEionDepth2 / sc.eIon - meanDepth * meanDepth)
            : 0.;
    const G4double depthSigma = std::sqrt(depthVar);
    const G4double f10 = (sc.eIon > 0.) ? sc.eIonAbove1 / sc.eIon : 0.;
    const G4double f100 = (sc.eIon > 0.) ? sc.eIonAbove2 / sc.eIon : 0.;
    const G4double f1000 = (sc.eIon > 0.) ? sc.eIonAbove3 / sc.eIon : 0.;
    const G4double W = fDetector->GetPairCreationEnergy(k);
    const G4double initialPairs = (W > 0.) ? sc.eIon / W : 0.;

    // Event-level histograms.
    analysis->FillH1(EionEventAllId(k), sc.eIon / keV, primaryWeight);
    if (sc.eIon > 0.) {
      analysis->FillH1(EionEventHitId(k), sc.eIon / keV, primaryWeight);
      analysis->FillH1(NielEventId(k), sc.niel / keV, primaryWeight);
      analysis->FillH1(LetTEventId(k), letT / keV_um, primaryWeight);
      analysis->FillH1(LetDdepEventId(k), letDdep / keV_um, primaryWeight);
      analysis->FillH1(LetDcalcEventId(k), letDcalc / keV_um, primaryWeight);
      analysis->FillH1(LetMaxEventId(k), sc.maxLetCalc / keV_um, primaryWeight);
      analysis->FillH1(FracAbove1Id(k), f10, primaryWeight);
      analysis->FillH1(FracAbove2Id(k), f100, primaryWeight);
      analysis->FillH1(FracAbove3Id(k), f1000, primaryWeight);
      analysis->FillH2(H2EventEionVsLetCalcId(k), sc.eIon / MeV,
                       letDcalc / keV_um, primaryWeight);
      analysis->FillH2(H2EventEionVsLetDepId(k), sc.eIon / MeV,
                       letDdep / keV_um, primaryWeight);
    }

    // EventLET ntuple row (column order MUST match HistoManager booking).
    const G4int nt = kNtupleEventLet;
    G4int col = 0;
    analysis->FillNtupleIColumn(nt, col++, runID);
    analysis->FillNtupleIColumn(nt, col++, eventID);
    analysis->FillNtupleIColumn(nt, col++, k);
    analysis->FillNtupleIColumn(nt, col++, primaryPDG);
    analysis->FillNtupleDColumn(nt, col++, primaryEnergy / MeV);
    analysis->FillNtupleDColumn(nt, col++, primaryWeight);
    analysis->FillNtupleDColumn(nt, col++, sc.eDep / keV);
    analysis->FillNtupleDColumn(nt, col++, sc.eIon / keV);
    analysis->FillNtupleDColumn(nt, col++, sc.niel / keV);
    analysis->FillNtupleDColumn(nt, col++, sc.chargedPath / um);
    analysis->FillNtupleDColumn(nt, col++, letT / keV_um);
    analysis->FillNtupleDColumn(nt, col++, letDdep / keV_um);
    analysis->FillNtupleDColumn(nt, col++, letDcalc / keV_um);
    analysis->FillNtupleDColumn(nt, col++, sc.maxLetDep / keV_um);
    analysis->FillNtupleDColumn(nt, col++, sc.maxLetCalc / keV_um);
    analysis->FillNtupleDColumn(nt, col++, f10);
    analysis->FillNtupleDColumn(nt, col++, f100);
    analysis->FillNtupleDColumn(nt, col++, f1000);
    analysis->FillNtupleDColumn(nt, col++, meanDepth / um);
    analysis->FillNtupleDColumn(nt, col++, depthSigma / um);
    analysis->FillNtupleIColumn(nt, col++, sc.nDepositingSteps);
    analysis->FillNtupleIColumn(nt, col++, sc.nChargedTracks);
    analysis->FillNtupleDColumn(nt, col++, initialPairs);
    analysis->AddNtupleRow(nt);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
