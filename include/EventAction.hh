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
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "DetectorConstruction.hh"
#include "HistoManager.hh"  // kMaxAbsor, kNbOrigin, histogram ID helpers

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <map>
#include <set>

class G4Track;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction(DetectorConstruction*);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    void SumEnergy(G4int k, G4double de, G4double dl);
    void SumPrimaryEnergy(G4int k, G4double de);
    void SumEnergyByParticle(G4int k, const G4String& particleName, G4double de);
    void SumEnergyLeak(G4double eleak, G4int index);

    // Mark that a neutron underwent a non-Transportation interaction in
    // absorber k during this event. Idempotent within an event.
    void MarkNeutronInteraction(G4int k);

    // Track-once gating for neutron entries into the front-face monitor,
    // so the same neutron is not scored multiple times on back-and-forth
    // boundary crossings.
    G4bool HasMonitorCrossed(G4int trackId) const;
    void MarkMonitorCrossed(G4int trackId);

    // Ancestry tag for e- tracks (gamma-mediated vs ion/hadronic-mediated).
    void TagElectronLineage(G4int trackId, G4int lineage);
    G4int GetElectronLineage(G4int trackId) const;
    void SumElectronEdepByLineage(G4int absorNum, G4int lineage, G4double de);

    // Particle name at track birth (first step), for e- ancestry classification.
    void RegisterTrackParticle(G4int trackId, const G4String& particleName);
    G4String GetTrackParticleName(G4int trackId) const;

    // Origin-absorber tag for cross-absorber energy-deposit provenance. A
    // track's origin is the absorber where the earliest ancestor of its lineage
    // was created (0 = external/primary lineage, i.e. never born in an
    // absorber). Set once at the track's first step in SteppingAction.
    void SetTrackOrigin(G4int trackId, G4int originAbsor);
    G4int GetTrackOrigin(G4int trackId) const;
    void SumEdepByOrigin(G4int depositAbsor, G4int originAbsor, G4double de);

    // Per-track path-length accumulation: sum of step lengths inside each
    // absorber for the track currently being transported. Cleared per track in
    // TrackingAction::PreUserTrackingAction and flushed (into per-particle,
    // per-absorber histograms) in PostUserTrackingAction. Relies on Geant4
    // transporting one track to completion before starting the next.
    void BeginTrackPath() { fTrackPathInAbsor.clear(); }
    void AddTrackPathLength(G4int absorNum, G4double len);
    const std::map<G4int, G4double>& GetTrackPathMap() const { return fTrackPathInAbsor; }

    // ---- LET / track-structure scoring --------------------------------------
    // Accumulate one charged/neutral step's contribution into both the current
    // event's per-absorber EventLetScore and the current track's per-absorber
    // TrackLetScore. Called from SteppingAction for charge-active absorbers.
    void AccumulateLetStep(G4int absorNum, const LetStepData& step);

    // Per-track LET map lifecycle (mirrors BeginTrackPath / GetTrackPathMap).
    // Cleared per track in TrackingAction::PreUserTrackingAction and flushed to
    // TrackLET ntuple rows in PostUserTrackingAction.
    void BeginTrackLet() { fTrackLetInAbsor.clear(); }
    const std::map<G4int, TrackLetScore>& GetTrackLetMap() const { return fTrackLetInAbsor; }

    // Per-event EventLetScore access (for EndOfEventAction row writing).
    const EventLetScore& GetEventLet(G4int absorNum) const { return fEventLet[absorNum]; }
    void CountCompletedChargedTrack(G4int absorNum);

    // Reaction-channel tag per track (set at the interaction, inherited by
    // descendants on their first step; see SteppingAction).
    void SetTrackReaction(G4int trackId, ReactionCategory category);
    ReactionCategory GetTrackReaction(G4int trackId) const;

    // Classify a depositing track into a stable DepositCategory (uses PDG,
    // particle name, atomic number, and the electron gamma/ion lineage tag).
    // Shared by SteppingAction (per step) and TrackingAction (per track flush).
    DepositCategory ClassifyDeposit(const G4Track* track) const;

  private:
    DetectorConstruction* fDetector = nullptr;

    G4double fEnergyDeposit[kMaxAbsor];
    G4double fPrimaryEdep[kMaxAbsor] = {0.};
    G4double fTrackLengthCh[kMaxAbsor];
    G4double fEnergyLeak[2] = {0., 0.};
    std::map<G4String, G4double> fEdepByParticle[kMaxAbsor];
    G4bool   fNeutronInteracted[kMaxAbsor] = {false};
    std::set<G4int> fMonitorCrossedTrackIds;

    std::map<G4int, G4int> fElectronLineage;
    std::map<G4int, G4String> fTrackParticleName;
    std::map<G4int, G4double> fTrackPathInAbsor;
    G4double fEdepElectronGamma[kMaxAbsor] = {0.};
    G4double fEdepElectronIonic[kMaxAbsor] = {0.};

    // Cross-absorber Edep provenance (per event): fEdepByOrigin[d-1][o] is the
    // energy deposited in absorber d by tracks whose lineage originated in
    // absorber o (o = 0 -> external/primary). Origin tags are per track.
    std::map<G4int, G4int> fTrackOriginAbsor;
    G4double fEdepByOrigin[kMaxAbsor][kNbOrigin] = {{0.}};

    // LET / track-structure scoring state (1-based absorber indexing).
    // fEventLet[absNum] is reset each event; fTrackLetInAbsor is reset per track;
    // fTrackReaction maps trackID -> its lineage's reaction channel.
    EventLetScore fEventLet[kMaxAbsor];
    std::map<G4int, TrackLetScore> fTrackLetInAbsor;
    std::map<G4int, ReactionCategory> fTrackReaction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
