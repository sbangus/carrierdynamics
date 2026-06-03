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
/// \file Run.hh
/// \brief Definition of the Run class

#ifndef Run_h
#define Run_h 1

#include "DetectorConstruction.hh"

#include "G4Run.hh"
#include "G4VProcess.hh"
#include "globals.hh"

#include <map>

class DetectorConstruction;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(DetectorConstruction*);
    ~Run() override = default;

  public:
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);
    void SetPrimaryDescription(const G4String& desc);
    void CountProcesses(const G4VProcess* process);
    void CountSecondaryBirth(G4int absIdx, const G4String& particleName, G4int creatorCat);

    // Count one gamma crossing into the front-face monitor (detector face,
    // just upstream of absorber 1). Used to integrate the gamma fluence.
    void CountGammaAtDetectorFace() { ++fGammaAtDetectorFace; }

    void SumEdepPerAbsorber(G4int, G4double, G4double);
    void SumEdepByParticle(G4int kAbs, const std::map<G4String, G4double>& edepMap);
    void SumEnergies(G4double edeptot, G4double eleak0, G4double eleak1);
    void SumEnergyFlow(G4int plane, G4double Eflow);

    // Per-event helpers used by EventAction at end-of-event.
    void AddInteractedFlags(const G4bool flags[]);

    void Merge(const G4Run*) override;
    void EndOfRun();

    /// Called by RunAction (master) before EndOfRun(); enables CSV append using maps below.
    void SetBatchCsvExport(G4bool enable, const G4String& path, const G4String& tag);

  private:
    void AppendBatchSummaryCsv();
    void PrintSecondaryBirthSummary() const;

    DetectorConstruction* fDetector = nullptr;
    G4ParticleDefinition* fParticle = nullptr;
    G4double fEkin = 0.;
    G4String fPrimaryDesc = "";

    G4double fSumEAbs[kMaxAbsor], fSum2EAbs[kMaxAbsor];
    G4double fSumLAbs[kMaxAbsor], fSum2LAbs[kMaxAbsor];

    G4double fEdepTot = 0., fEdepTot2 = 0.;

    G4double fEnergyLeak[2] = {0., 0.};
    G4double fEleakTot = 0., fEleakTot2 = 0.;

    G4double fEtotal = 0., fEtotal2 = 0.;

    std::map<G4String, G4int> fProcCounter;
    std::map<G4String, G4long> fSecBirthByPart[kMaxAbsor];
    std::map<G4int, G4long> fSecBirthByCreatorCat[kMaxAbsor];
    std::vector<G4double> fEnergyFlow;
    std::map<G4String, G4double> fEdepByParticle[kMaxAbsor];

    // Number of events in which a primary neutron interacted in absorber k
    // (one entry per absorber slot). Used to compute interaction fraction.
    G4long fNeutronInteractedEvents[kMaxAbsor] = {0};

    // Total number of gammas that crossed into the front-face monitor slab
    // (detector face). Integrated over the run; divided by the monitor face
    // area in EndOfRun to report the gamma fluence.
    G4long fGammaAtDetectorFace = 0;

    // Batch CSV (Option B): set each run end by RunAction before EndOfRun().
    G4bool fBatchCsvEnable = false;
    G4String fBatchCsvPath;
    G4String fBatchCsvTag;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
