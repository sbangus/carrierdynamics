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
/// \file HistoManager.hh
/// \brief Definition of the HistoManager class

#ifndef HistoManager_h
#define HistoManager_h 1

#include "G4AnalysisManager.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

// ---------------------------------------------------------------------------
// Histogram ID layout
// ---------------------------------------------------------------------------
//
// Existing histograms (legacy IDs, do not renumber):
//   1 .. kMaxAbsor                     : Edep per absorber (per event)
//   kMaxAbsor+1 .. 2*kMaxAbsor         : Longitudinal Edep profile per absorber
//   2*kMaxAbsor+1                      : Energy flow
//   2*kMaxAbsor+2                      : Total energy deposited
//   2*kMaxAbsor+3                      : Total energy leakage
//   2*kMaxAbsor+4                      : Total energy released
//
// Depth of neutron interactions (elastic, inelastic+capture) per absorber:
//   kDepthHistoBase + (a-1)*kNbProcessCat + c
// where a = absorber index (1..kMaxAbsor), c = process category (0..kNbProcessCat-1)

const G4int kNbProcessCat = 2;  // 0: hadElastic, 1: nInelastic (incl. nCapture)
const G4int kDepthHistoBase = 2 * kMaxAbsor + 5;
const G4int kDepthHistoEnd  = kDepthHistoBase + kMaxAbsor * kNbProcessCat;

inline G4int DepthHistoId(G4int absIdx, G4int processCat) {
    return kDepthHistoBase + (absIdx - 1) * kNbProcessCat + processCat;
}

// ---------------------------------------------------------------------------
// Per-particle Edep per absorber: kEdepByParticleBase + (a-1)*kNbFixedParticles + p
const G4int kNbFixedParticles = 7;  // alpha, Li7, Li6, proton, gamma, e-, C12
const G4int kEdepByParticleBase = kDepthHistoEnd;

inline G4int EdepByParticleId(G4int absIdx, G4int partIdx) {
    return kEdepByParticleBase + (absIdx - 1) * kNbFixedParticles + partIdx;
}

// Capture gamma spectrum per absorber: kCaptureGammaBase + (a-1)
const G4int kCaptureGammaBase = kEdepByParticleBase + kMaxAbsor * kNbFixedParticles;

inline G4int CaptureGammaId(G4int absIdx) {
    return kCaptureGammaBase + (absIdx - 1);
}

// Estimated EHP per event per absorber: kEhpBase + (a-1)
const G4int kEhpBase = kCaptureGammaBase + kMaxAbsor;

inline G4int EhpId(G4int absIdx) {
    return kEhpBase + (absIdx - 1);
}

// Secondary production depth per absorber: kSecDepthBase + (a-1)
const G4int kSecDepthBase = kEhpBase + kMaxAbsor;

inline G4int SecDepthId(G4int absIdx) {
    return kSecDepthBase + (absIdx - 1);
}

// Energy spectrum of each secondary particle type (summed across absorbers):
//   kSecESpecBase + p
const G4int kSecESpecBase = kSecDepthBase + kMaxAbsor;

inline G4int SecESpecId(G4int partIdx) {
    return kSecESpecBase + partIdx;
}

// Neutron energy at each interaction per absorber (log binning):
//   kNeutronEnergyBase + (a-1)
const G4int kNeutronEnergyBase = kSecESpecBase + kNbFixedParticles;

inline G4int NeutronEnergyId(G4int absIdx) {
    return kNeutronEnergyBase + (absIdx - 1);
}

// Single-histogram summaries (one bin per absorber):
const G4int kInteractionFractionId = kNeutronEnergyBase + kMaxAbsor;
const G4int kComparativeResponseId = kInteractionFractionId + 1;

// KE at generation of first-generation secondaries (ParentID == 1):
//   kFirstGenSecESpecBase + p
const G4int kFirstGenSecESpecBase = kComparativeResponseId + 1;

inline G4int FirstGenSecESpecId(G4int partIdx) {
    return kFirstGenSecESpecBase + partIdx;
}

// Neutron kinetic energy at first entry into the front-face monitor (i.e.
// the spectrum of neutrons reaching the detector after any environmental
// interactions: bud-box wall scattering, air, etc.). Two histograms:
//   kDetectorFrontKEAll:     all neutron tracks (primaries + secondaries)
//   kDetectorFrontKEPrimary: primary neutrons only (ParentID == 0)
const G4int kDetectorFrontKEAll = kFirstGenSecESpecBase + kNbFixedParticles;
const G4int kDetectorFrontKEPrimary = kDetectorFrontKEAll + 1;

// Gamma kinetic-energy spectrum at the detector front face: every gamma that
// crosses into the front-face monitor slab (immediately upstream of absorber 1
// = the first silver layer), scored once per track. This histogram is the
// energy distribution of the gamma field reaching the detector face; its
// integrated entry count divided by the monitor face area gives the gamma
// fluence reported in Run::EndOfRun.
const G4int kDetectorFrontGammaKE = kDetectorFrontKEPrimary + 1;

// Per-absorber e- Edep split by ancestry (in addition to total e- in EdepByParticleId):
//   kEdepElectronGammaBase + (a-1)  gamma-mediated (compt/phot/conv chains, incl. deltas)
//   kEdepElectronIonicBase + (a-1)  ion/hadronic-mediated (recoils, hIoni, etc.)
const G4int kElectronLineageNone = 0;
const G4int kElectronLineageGamma = 1;
const G4int kElectronLineageIon = 2;

const G4int kEdepElectronGammaBase = kDetectorFrontGammaKE + 1;
const G4int kEdepElectronIonicBase = kEdepElectronGammaBase + kMaxAbsor;

inline G4int EdepElectronGammaId(G4int absIdx) {
    return kEdepElectronGammaBase + (absIdx - 1);
}

inline G4int EdepElectronIonicId(G4int absIdx) {
    return kEdepElectronIonicBase + (absIdx - 1);
}

const G4int kMaxHisto = kEdepElectronIonicBase + kMaxAbsor;

// Creator-process categories for secondary-birth run summary (Run::PrintSecondaryBirthSummary).
const G4int kNbSecCreatorCat = 5;  // hadronic, em_ion, em_photon, em_other, other

inline const G4String& SecCreatorCategoryName(G4int cat) {
    static const G4String names[kNbSecCreatorCat] =
        {"hadronic", "em_ion", "em_photon", "em_other", "other"};
    return names[cat];
}

inline G4int SecCreatorCategory(const G4String& procName) {
    if (procName == "nCapture" || procName == "neutronInelastic"
        || procName == "hadElastic" || procName == "protonInelastic"
        || procName == "ionInelastic" || procName == "alphaInelastic"
        || procName == "nuclElastic" || procName == "hFritiofCaptureAtRest")
    {
        return 0;
    }
    if (procName == "eIoni" || procName == "muIoni" || procName == "hIoni"
        || procName == "ionIoni")
    {
        return 1;
    }
    if (procName == "phot" || procName == "compt" || procName == "conv"
        || procName == "Rayleigh")
    {
        return 2;
    }
    if (procName == "eBrem" || procName == "muBrems" || procName == "hBrems"
        || procName == "ionBrems" || procName == "annihil" || procName == "ePairProd")
    {
        return 3;
    }
    return 4;
}

// Map a Geant4 particle name to a fixed-particle index.
// Returns -1 if the particle is not in the tracked list.
// Strips ion excitation suffix (e.g. "Li7[478.0]" -> "Li7") before matching.
inline G4int FixedParticleIdx(const G4String& name) {
    G4String stripped = name.substr(0, name.find('['));
    if (stripped == "alpha")  return 0;
    if (stripped == "Li7")    return 1;
    if (stripped == "Li6")    return 2;
    if (stripped == "proton") return 3;
    if (stripped == "gamma")  return 4;
    if (stripped == "e-")     return 5;
    if (stripped == "C12")   return 6;
    return -1;
}

inline const G4String& FixedParticleName(G4int idx) {
    static const G4String names[kNbFixedParticles] =
        {"alpha", "Li7", "Li6", "proton", "gamma", "e-", "C12"};
    return names[idx];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
  public:
    HistoManager();
    ~HistoManager() = default;

  private:
    void Book();
    G4String fFileName = "Hadr05";
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
