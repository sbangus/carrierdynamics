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
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

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

void EventAction::EndOfEventAction(const G4Event*)
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
