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
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "Run.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4StepStatus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* det, EventAction* evt)
  : fDetector(det), fEventAct(evt)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  // get Run
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  // Energy flow initialisation for primary particle
  //
  if (track->GetTrackID() == 1) {
    G4int Idnow = 1;
    if (track->GetVolume() != fDetector->GetphysiWorld()) {
      // unique identificator of layer+absorber
      const G4VTouchable* touchable = track->GetTouchable();
      G4int absorNum = touchable->GetCopyNumber();
      G4int layerNum = touchable->GetReplicaNumber(1);
      Idnow = (fDetector->GetNbOfAbsor()) * layerNum + absorNum;
    }

    G4double Eflow = track->GetKineticEnergy();
    /// if (track->GetDefinition() == G4Positron::Positron()) {
    ///   Eflow += 2*electron_mass_c2;
    /// }

    // flux artefact, if primary vertex is inside the calorimeter
    for (G4int pl = 1; pl <= Idnow; ++pl) {
      run->SumEnergyFlow(pl, Eflow);
    }
  }

  // Reset the per-track path-length accumulator for the track about to be
  // transported (Geant4 finishes one track before starting the next).
  fEventAct->BeginTrackPath();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  // energy leakage
  G4StepStatus status = aTrack->GetStep()->GetPostStepPoint()->GetStepStatus();
  if (status == fWorldBoundary) {
    G4int parentID = aTrack->GetParentID();
    G4int index = 0;
    if (parentID > 0) index = 1;  // primary=0, secondaries=1
    G4double eleak = aTrack->GetKineticEnergy();
    fEventAct->SumEnergyLeak(eleak, index);
  }

  // Flush per-track path-length accumulation into per-particle histograms.
  // Secondaries only (the primary is the neutron source, not among the fixed
  // particle types). Lengths are raw internal units (mm); histograms carry the
  // "mm" unit so the axis reads in mm.
  if (aTrack->GetParentID() > 0) {
    const G4int pi = FixedParticleIdx(aTrack->GetDefinition()->GetParticleName());
    if (pi >= 0) {
      auto* analysis = G4AnalysisManager::Instance();
      for (const auto& entry : fEventAct->GetTrackPathMap()) {
        if (entry.second > 0.) {
          analysis->FillH1(TrackPathLengthId(entry.first, pi), entry.second);
        }
      }
      const G4double total = aTrack->GetTrackLength();
      if (total > 0.) analysis->FillH1(TotalTrackLengthId(pi), total);

      // e- path length split by ancestry, mirroring the e- Edep lineage split.
      if (pi == 5) {
        const G4int lineage = fEventAct->GetElectronLineage(aTrack->GetTrackID());
        if (lineage == kElectronLineageGamma || lineage == kElectronLineageIon) {
          for (const auto& entry : fEventAct->GetTrackPathMap()) {
            if (entry.second > 0.) {
              const G4int id = (lineage == kElectronLineageGamma)
                                 ? PathElectronGammaId(entry.first)
                                 : PathElectronIonicId(entry.first);
              analysis->FillH1(id, entry.second);
            }
          }
          if (total > 0.) {
            analysis->FillH1((lineage == kElectronLineageGamma)
                               ? kTotalPathElectronGamma
                               : kTotalPathElectronIonic,
                             total);
          }
        }
      }
    }
  }

  // Primary-particle path length (ParentID == 0): per-absorber accumulation and
  // total track length, particle-type agnostic (alpha, neutron, gamma, etc.).
  if (aTrack->GetParentID() == 0) {
    auto* analysis = G4AnalysisManager::Instance();
    for (const auto& entry : fEventAct->GetTrackPathMap()) {
      if (entry.second > 0.) {
        analysis->FillH1(PrimaryPathLengthId(entry.first), entry.second);
      }
    }
    const G4double total = aTrack->GetTrackLength();
    if (total > 0.) analysis->FillH1(kTotalPrimaryPathLength, total);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
