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
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "G4StepStatus.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

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
  fEventAct->BeginTrackLet();
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

  // ---- LET / track-structure: flush per-track TrackLET rows ----------------
  // One row per absorber in which this track deposited ionizing energy. Track
  // LET metrics are unweighted physical quantities; the statistical track weight
  // is stored as its own column. Track-category QA histograms are filled here.
  {
    auto* analysis = G4AnalysisManager::Instance();
    const std::map<G4int, TrackLetScore>& trackLet = fEventAct->GetTrackLetMap();
    if (!trackLet.empty()) {
      const G4ParticleDefinition* def = aTrack->GetDefinition();
      const G4int trackID = aTrack->GetTrackID();
      const G4int parentID = aTrack->GetParentID();
      const G4int pdg = def->GetPDGEncoding();
      const G4double charge = def->GetPDGCharge();
      const G4int atomicZ = def->GetAtomicNumber();
      const G4int atomicA = def->GetAtomicMass();
      const DepositCategory depCat = fEventAct->ClassifyDeposit(aTrack);
      const ReactionCategory rxnCat = fEventAct->GetTrackReaction(trackID);
      const G4int lineage = fEventAct->GetElectronLineage(trackID);
      const G4int originAbsor = fEventAct->GetTrackOrigin(trackID);
      const G4double vertexKE = aTrack->GetVertexKineticEnergy();
      const G4double trackWeight = aTrack->GetWeight();

      const G4Event* evt =
          G4EventManager::GetEventManager()->GetConstCurrentEvent();
      const G4int eventID = (evt != nullptr) ? evt->GetEventID() : -1;
      const G4int runID =
          G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
      const G4double keV_um = keV / um;

      for (const auto& [absorNum, score] : trackLet) {
        if (score.eIon <= 0.) continue;  // depositing tracks only

        const G4double trackLET =
            (score.chargedPath > 0.) ? score.eIon / score.chargedPath : 0.;
        const G4double letDdep = score.sumEionLetDep / score.eIon;
        const G4double letDcalc = score.sumEionLetCalc / score.eIon;
        const G4double meanDepth = score.sumEionDepth / score.eIon;
        const EscapeCategory escape = ResolveEscapeCategory(score);
        const G4int stoppedInLayer =
            (escape == EscapeCategory::Contained) ? 1 : 0;
        const G4double depthSpan =
            (score.maxDepth >= score.minDepth) ? (score.maxDepth - score.minDepth) : 0.;

        // Track-category QA histograms.
        const G4int c = static_cast<G4int>(depCat);
        analysis->FillH1(TrackEionCatId(absorNum, c), score.eIon / keV, trackWeight);
        analysis->FillH1(TrackLetCatId(absorNum, c), trackLET / keV_um, trackWeight);
        analysis->FillH1(TrackDepthSpanCatId(absorNum, c), depthSpan / um, trackWeight);
        // Containment: categorical histogram (bin = EscapeCategory value).
        analysis->FillH1(TrackContainmentId(absorNum),
                         static_cast<G4double>(static_cast<G4int>(escape)),
                         trackWeight);
        analysis->FillH2(H2TrackEionVsLetCalcId(absorNum),
                         score.eIon / keV, letDcalc / keV_um, trackWeight);
        analysis->FillH2(H2TrackEionVsLetDepId(absorNum),
                         score.eIon / keV, letDdep / keV_um, trackWeight);

        // Primary-only track LET (incident particle; ParentID == 0, charged).
        if (parentID == 0 && charge != 0.) {
          analysis->FillH1(PrimaryTrackEionId(absorNum), score.eIon / keV, trackWeight);
          analysis->FillH1(PrimaryTrackLetId(absorNum), trackLET / keV_um, trackWeight);
          analysis->FillH1(PrimaryTrackDepthSpanId(absorNum), depthSpan / um, trackWeight);
        }

        // TrackLET ntuple row (column order MUST match HistoManager booking).
        const G4int nt = kNtupleTrackLet;
        G4int col = 0;
        analysis->FillNtupleIColumn(nt, col++, runID);
        analysis->FillNtupleIColumn(nt, col++, eventID);
        analysis->FillNtupleIColumn(nt, col++, trackID);
        analysis->FillNtupleIColumn(nt, col++, parentID);
        analysis->FillNtupleIColumn(nt, col++, absorNum);
        analysis->FillNtupleIColumn(nt, col++, pdg);
        analysis->FillNtupleDColumn(nt, col++, charge);
        analysis->FillNtupleIColumn(nt, col++, atomicZ);
        analysis->FillNtupleIColumn(nt, col++, atomicA);
        analysis->FillNtupleIColumn(nt, col++, static_cast<G4int>(depCat));
        analysis->FillNtupleIColumn(nt, col++, static_cast<G4int>(rxnCat));
        analysis->FillNtupleIColumn(nt, col++, lineage);
        analysis->FillNtupleIColumn(nt, col++, originAbsor);
        analysis->FillNtupleDColumn(nt, col++, vertexKE / MeV);
        analysis->FillNtupleDColumn(nt, col++, score.entryKE / MeV);
        analysis->FillNtupleDColumn(nt, col++, score.exitKE / MeV);
        analysis->FillNtupleDColumn(nt, col++, score.eDep / keV);
        analysis->FillNtupleDColumn(nt, col++, score.eIon / keV);
        analysis->FillNtupleDColumn(nt, col++, score.niel / keV);
        analysis->FillNtupleDColumn(nt, col++, score.chargedPath / um);
        analysis->FillNtupleDColumn(nt, col++, trackLET / keV_um);
        analysis->FillNtupleDColumn(nt, col++, letDdep / keV_um);
        analysis->FillNtupleDColumn(nt, col++, letDcalc / keV_um);
        analysis->FillNtupleDColumn(nt, col++, score.maxLetDep / keV_um);
        analysis->FillNtupleDColumn(nt, col++, score.maxLetCalc / keV_um);
        analysis->FillNtupleDColumn(nt, col++, score.startDepth / um);
        analysis->FillNtupleDColumn(nt, col++, score.endDepth / um);
        analysis->FillNtupleDColumn(nt, col++, score.minDepth / um);
        analysis->FillNtupleDColumn(nt, col++, score.maxDepth / um);
        analysis->FillNtupleDColumn(nt, col++, meanDepth / um);
        analysis->FillNtupleIColumn(nt, col++, static_cast<G4int>(escape));
        analysis->FillNtupleIColumn(nt, col++, stoppedInLayer);
        analysis->FillNtupleIColumn(nt, col++, score.nDepositingSteps);
        analysis->FillNtupleDColumn(nt, col++, trackWeight);
        analysis->AddNtupleRow(nt);

        fEventAct->CountCompletedChargedTrack(absorNum);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
