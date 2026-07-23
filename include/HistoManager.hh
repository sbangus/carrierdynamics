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
#include "LetScoring.hh"  // DepositCategory, kNbDepositCategories, thresholds

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
// Depth of selected interactions per absorber:
//   kDepthHistoBase + (a-1)*kNbProcessCat + c
// where a = absorber index (1..kMaxAbsor), c = process category (0..kNbProcessCat-1)
// Categories: 0 neutron hadElastic, 1 neutron inelastic (incl. nCapture),
//             2 gamma Compton (compt), 3 gamma photoelectric (phot),
//             4 gamma pair production (conv).

const G4int kNbProcessCat = 5;
const G4int kProcCatHadElastic = 0;
const G4int kProcCatNInelastic = 1;
const G4int kProcCatCompton = 2;
const G4int kProcCatPhotoelectric = 3;
const G4int kProcCatPairProd = 4;
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

// Capture-gamma line spectrum for gammas born in the moderator stack
// (nCapture / neutronInelastic): linear 0-3 MeV for H(n,gamma) at 2.223 MeV.
const G4int kModeratorCaptureGammaKE = kEdepElectronIonicBase + kMaxAbsor;

// KE at birth in moderator, by particle type (summed across moderator layers):
//   kModeratorSecESpecBase + p
const G4int kModeratorSecESpecBase = kModeratorCaptureGammaKE + 1;

inline G4int ModeratorSecESpecId(G4int partIdx) {
    return kModeratorSecESpecBase + partIdx;
}

// First-generation (ParentID == 1) moderator-born secondaries:
//   kModeratorFirstGenSecESpecBase + p
const G4int kModeratorFirstGenSecESpecBase =
    kModeratorSecESpecBase + kNbFixedParticles;

inline G4int ModeratorFirstGenSecESpecId(G4int partIdx) {
    return kModeratorFirstGenSecESpecBase + partIdx;
}

// Per-absorber, per-particle secondary track path length (sum of step lengths
// while the track is inside the absorber; observable A):
//   kTrackPathLengthBase + (a-1)*kNbFixedParticles + p
const G4int kTrackPathLengthBase =
    kModeratorFirstGenSecESpecBase + kNbFixedParticles;

inline G4int TrackPathLengthId(G4int absIdx, G4int partIdx) {
    return kTrackPathLengthBase + (absIdx - 1) * kNbFixedParticles + partIdx;
}

// Total track length per particle type, summed over all volumes (observable B):
//   kTotalTrackLengthBase + p
const G4int kTotalTrackLengthBase =
    kTrackPathLengthBase + kMaxAbsor * kNbFixedParticles;

inline G4int TotalTrackLengthId(G4int partIdx) {
    return kTotalTrackLengthBase + partIdx;
}

// Per-absorber e- path length split by ancestry (analog of EdepElectron*Id):
//   kPathElectronGammaBase + (a-1)  gamma-mediated (compt/phot/conv chains, incl. deltas)
//   kPathElectronIonicBase + (a-1)  ion/hadronic-mediated (recoils, hIoni, etc.)
const G4int kPathElectronGammaBase = kTotalTrackLengthBase + kNbFixedParticles;
const G4int kPathElectronIonicBase = kPathElectronGammaBase + kMaxAbsor;

inline G4int PathElectronGammaId(G4int absIdx) {
    return kPathElectronGammaBase + (absIdx - 1);
}

inline G4int PathElectronIonicId(G4int absIdx) {
    return kPathElectronIonicBase + (absIdx - 1);
}

// Total (all-volumes) e- track length split by the same ancestry buckets:
const G4int kTotalPathElectronGamma = kPathElectronIonicBase + kMaxAbsor;
const G4int kTotalPathElectronIonic = kTotalPathElectronGamma + 1;

// Secondary e- KE-at-birth spectrum split by ancestry (split of SecESpecId(e-),
// summed across absorbers): gamma-mediated (Compton/photo/pair) vs ion/hadronic.
const G4int kSecESpecElectronGamma = kTotalPathElectronIonic + 1;
const G4int kSecESpecElectronIonic = kSecESpecElectronGamma + 1;

// First-generation (ParentID == 1) secondary e- KE-at-birth split by ancestry
// (split of FirstGenSecESpecId(e-)):
const G4int kFirstGenSecESpecElectronGamma = kSecESpecElectronIonic + 1;
const G4int kFirstGenSecESpecElectronIonic = kFirstGenSecESpecElectronGamma + 1;

// Primary-particle observables (ParentID == 0; particle-type agnostic):
//   kPrimaryPathLengthBase + (a-1)  path length inside absorber a (per primary track)
//   kTotalPrimaryPathLength         full track length across all volumes (per primary track)
//   kPrimaryEdepBase + (a-1)        Edep in absorber a summed per event (primary steps only)
const G4int kPrimaryPathLengthBase = kFirstGenSecESpecElectronIonic + 1;
const G4int kTotalPrimaryPathLength = kPrimaryPathLengthBase + kMaxAbsor;

inline G4int PrimaryPathLengthId(G4int absIdx) {
    return kPrimaryPathLengthBase + (absIdx - 1);
}

const G4int kPrimaryEdepBase = kTotalPrimaryPathLength + 1;

inline G4int PrimaryEdepId(G4int absIdx) {
    return kPrimaryEdepBase + (absIdx - 1);
}

// Cross-absorber energy-deposit provenance: energy deposited in absorber
// `depositAbs` (1..kMaxAbsor) split by the "origin absorber" of the depositing
// track's lineage. The origin index is:
//   0             : external / primary lineage (the track descends only from
//                   the source particle; it was never born inside an absorber)
//   1..kMaxAbsor  : the absorber in which the earliest ancestor of the track
//                   was created (e.g. a gamma born in the silver-epoxy blob that
//                   Compton-scatters an e- which deposits in a downstream layer
//                   -> origin = the epoxy absorber index).
// Layout: kEdepByOriginBase + (depositAbs-1)*kNbOrigin + originIdx.
// The diagonal-like case originIdx == depositAbs is "born and deposited in the
// same absorber"; originIdx == 0 is "deposited directly by the primary lineage".
const G4int kNbOrigin = kMaxAbsor + 1;  // origin 0 (external) + 1..kMaxAbsor
const G4int kEdepByOriginBase = kPrimaryEdepBase + kMaxAbsor;

inline G4int EdepByOriginId(G4int depositAbs, G4int originIdx) {
    return kEdepByOriginBase + (depositAbs - 1) * kNbOrigin + originIdx;
}

const G4int kMaxHisto = kEdepByOriginBase + kMaxAbsor * kNbOrigin;

// ===========================================================================
// LET / track-structure H1 histogram IDs (new family; all IDs >= kMaxHisto).
// ---------------------------------------------------------------------------
// The legacy histograms above are retired (booked but deactivated in Run.cc).
// These new histograms carry the LET analysis. Two ID layouts are used:
//   * Per-(absorber, deposit-category) families: kMaxAbsor*kNbDepositCategories
//     slots each; indexed LetXxxId(a, c) with a = 1-based absorber number and
//     c = DepositCategory integer value (0..kNbDepositCategories-1).
//   * Per-absorber families: kMaxAbsor slots each; indexed XxxId(a).
// Everything is booked generically (over the full kMaxAbsor range) in
// HistoManager::Book(); only charge-active absorbers are configured and
// activated per-geometry in Run::Run().
// ===========================================================================

const G4int kLetCatStride = kMaxAbsor * kNbDepositCategories;

// Per-(absorber, category) step- and track-level LET spectra.
const G4int kLetStepCountBase      = kMaxHisto;                          // Ldep number spectrum
const G4int kLetStepEWeightedBase  = kLetStepCountBase     + kLetCatStride;  // Eion-weighted Ldep
const G4int kLetCalcEWeightedBase  = kLetStepEWeightedBase + kLetCatStride;  // Eion-weighted Lcalc
const G4int kTrackEionCatBase      = kLetCalcEWeightedBase + kLetCatStride;  // per-track Eion
const G4int kTrackLetCatBase       = kTrackEionCatBase     + kLetCatStride;  // per-track LET
const G4int kTrackDepthSpanCatBase = kTrackLetCatBase      + kLetCatStride;  // per-track depth span

inline G4int LetStepCountId(G4int a, G4int c) {
    return kLetStepCountBase + (a - 1) * kNbDepositCategories + c;
}
inline G4int LetStepEWeightedId(G4int a, G4int c) {
    return kLetStepEWeightedBase + (a - 1) * kNbDepositCategories + c;
}
inline G4int LetCalcEWeightedId(G4int a, G4int c) {
    return kLetCalcEWeightedBase + (a - 1) * kNbDepositCategories + c;
}
inline G4int TrackEionCatId(G4int a, G4int c) {
    return kTrackEionCatBase + (a - 1) * kNbDepositCategories + c;
}
inline G4int TrackLetCatId(G4int a, G4int c) {
    return kTrackLetCatBase + (a - 1) * kNbDepositCategories + c;
}
inline G4int TrackDepthSpanCatId(G4int a, G4int c) {
    return kTrackDepthSpanCatBase + (a - 1) * kNbDepositCategories + c;
}

// Convenience overloads taking a DepositCategory directly.
inline G4int LetStepCountId(G4int a, DepositCategory c) {
    return LetStepCountId(a, static_cast<G4int>(c));
}
inline G4int LetStepEWeightedId(G4int a, DepositCategory c) {
    return LetStepEWeightedId(a, static_cast<G4int>(c));
}
inline G4int LetCalcEWeightedId(G4int a, DepositCategory c) {
    return LetCalcEWeightedId(a, static_cast<G4int>(c));
}
inline G4int TrackEionCatId(G4int a, DepositCategory c) {
    return TrackEionCatId(a, static_cast<G4int>(c));
}
inline G4int TrackLetCatId(G4int a, DepositCategory c) {
    return TrackLetCatId(a, static_cast<G4int>(c));
}
inline G4int TrackDepthSpanCatId(G4int a, DepositCategory c) {
    return TrackDepthSpanCatId(a, static_cast<G4int>(c));
}

// Per-absorber event/track singleton families.
const G4int kEionEventAllBase  = kTrackDepthSpanCatBase + kLetCatStride;  // all events (incl. zero)
const G4int kEionEventHitBase  = kEionEventAllBase  + kMaxAbsor;          // hit-only Eion
const G4int kNielEventBase     = kEionEventHitBase  + kMaxAbsor;          // NIEL per event
const G4int kLetTEventBase     = kNielEventBase     + kMaxAbsor;          // track-length-avg LET
const G4int kLetDdepEventBase  = kLetTEventBase     + kMaxAbsor;          // Eion-weighted Ldep
const G4int kLetDcalcEventBase = kLetDdepEventBase  + kMaxAbsor;          // Eion-weighted Lcalc
const G4int kLetMaxEventBase   = kLetDcalcEventBase + kMaxAbsor;          // max step LET
const G4int kFracAbove1Base    = kLetMaxEventBase   + kMaxAbsor;          // Eion frac Lcalc>=thr1
const G4int kFracAbove2Base    = kFracAbove1Base    + kMaxAbsor;          // >=thr2
const G4int kFracAbove3Base    = kFracAbove2Base    + kMaxAbsor;          // >=thr3
const G4int kEionDepthBase     = kFracAbove3Base    + kMaxAbsor;          // dEion/dz (absolute)
const G4int kEionNormDepthBase = kEionDepthBase     + kMaxAbsor;          // dEion/du (normalized)
const G4int kTrackContainBase  = kEionNormDepthBase + kMaxAbsor;          // containment (categorical)

inline G4int EionEventAllId(G4int a)     { return kEionEventAllBase  + a - 1; }
inline G4int EionEventHitId(G4int a)     { return kEionEventHitBase  + a - 1; }
inline G4int NielEventId(G4int a)        { return kNielEventBase     + a - 1; }
inline G4int LetTEventId(G4int a)        { return kLetTEventBase     + a - 1; }
inline G4int LetDdepEventId(G4int a)     { return kLetDdepEventBase  + a - 1; }
inline G4int LetDcalcEventId(G4int a)    { return kLetDcalcEventBase + a - 1; }
inline G4int LetMaxEventId(G4int a)      { return kLetMaxEventBase   + a - 1; }
inline G4int FracAbove1Id(G4int a)       { return kFracAbove1Base    + a - 1; }
inline G4int FracAbove2Id(G4int a)       { return kFracAbove2Base    + a - 1; }
inline G4int FracAbove3Id(G4int a)       { return kFracAbove3Base    + a - 1; }
inline G4int EionDepthId(G4int a)        { return kEionDepthBase     + a - 1; }
inline G4int EionNormDepthId(G4int a)    { return kEionNormDepthBase + a - 1; }
inline G4int TrackContainmentId(G4int a) { return kTrackContainBase  + a - 1; }

// Upper bound of the LET H1 ID space (one past the last LET H1 ID).
const G4int kMaxHistoLet = kTrackContainBase + kMaxAbsor;

// ---------------------------------------------------------------------------
// LET H2 histogram IDs (independent H2 ID space, starting at 0). Per-absorber
// families of kMaxAbsor slots each. Booked in Book() (in this exact order so
// the helper IDs match g4tools' sequential H2 counter) and configured/activated
// for charge-active absorbers in Run::Run().
// ---------------------------------------------------------------------------
const G4int kH2EventEionVsLetCalcBase = 0;                                  // x=Eion, y=event Lcalc
const G4int kH2EventEionVsLetDepBase  = kH2EventEionVsLetCalcBase + kMaxAbsor;  // x=Eion, y=event Ldep
const G4int kH2LetCalcVsDepthBase     = kH2EventEionVsLetDepBase  + kMaxAbsor;  // x=depth, y=Lcalc
const G4int kH2LetDepVsDepthBase      = kH2LetCalcVsDepthBase     + kMaxAbsor;  // x=depth, y=Ldep
const G4int kH2LetCalcVsKEBase        = kH2LetDepVsDepthBase      + kMaxAbsor;  // x=KE, y=Lcalc
const G4int kH2TrackEionVsLetCalcBase = kH2LetCalcVsKEBase        + kMaxAbsor;  // x=track Eion, y=track LET
const G4int kH2TrackEionVsLetDepBase  = kH2TrackEionVsLetCalcBase + kMaxAbsor;  // x=track Eion, y=track Ldep
const G4int kH2DepthVsTransverseBase  = kH2TrackEionVsLetDepBase  + kMaxAbsor;  // x=depth, y=local y/z
const G4int kNbH2Let                  = kH2DepthVsTransverseBase  + kMaxAbsor;

inline G4int H2EventEionVsLetCalcId(G4int a) { return kH2EventEionVsLetCalcBase + a - 1; }
inline G4int H2EventEionVsLetDepId(G4int a)  { return kH2EventEionVsLetDepBase  + a - 1; }
inline G4int H2LetCalcVsDepthId(G4int a)     { return kH2LetCalcVsDepthBase     + a - 1; }
inline G4int H2LetDepVsDepthId(G4int a)      { return kH2LetDepVsDepthBase      + a - 1; }
inline G4int H2LetCalcVsKEId(G4int a)        { return kH2LetCalcVsKEBase        + a - 1; }
inline G4int H2TrackEionVsLetCalcId(G4int a) { return kH2TrackEionVsLetCalcBase + a - 1; }
inline G4int H2TrackEionVsLetDepId(G4int a)  { return kH2TrackEionVsLetDepBase  + a - 1; }
inline G4int H2DepthVsTransverseId(G4int a)  { return kH2DepthVsTransverseBase  + a - 1; }

// ---------------------------------------------------------------------------
// LET ntuple IDs (independent ntuple ID space, created in Book() in this order).
// EventLET and TrackLET are written for every production run; StepLET is a
// sampled/gated diagnostic ntuple.
// ---------------------------------------------------------------------------
const G4int kNtupleEventLet = 0;
const G4int kNtupleTrackLet = 1;
const G4int kNtupleStepLet  = 2;

// Human-readable names for the deposit categories (for H1 titles / QA).
inline const G4String& DepositCategoryName(G4int c) {
    static const G4String names[kNbDepositCategories] = {
        "e-(gamma)", "e-(ion)", "e+", "proton", "lightIon", "alpha",
        "Li", "C", "Be/B", "otherHeavyIon", "neutralLocal", "other"};
    return (c >= 0 && c < kNbDepositCategories) ? names[c] : names[kNbDepositCategories - 1];
}

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
