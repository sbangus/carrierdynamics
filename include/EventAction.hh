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

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <map>
#include <set>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction(DetectorConstruction*);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    void SumEnergy(G4int k, G4double de, G4double dl);
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

  private:
    DetectorConstruction* fDetector = nullptr;

    G4double fEnergyDeposit[kMaxAbsor];
    G4double fTrackLengthCh[kMaxAbsor];
    G4double fEnergyLeak[2] = {0., 0.};
    std::map<G4String, G4double> fEdepByParticle[kMaxAbsor];
    G4bool   fNeutronInteracted[kMaxAbsor] = {false};
    std::set<G4int> fMonitorCrossedTrackIds;

    std::map<G4int, G4int> fElectronLineage;
    std::map<G4int, G4String> fTrackParticleName;
    G4double fEdepElectronGamma[kMaxAbsor] = {0.};
    G4double fEdepElectronIonic[kMaxAbsor] = {0.};
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
