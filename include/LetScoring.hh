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
/// \file LetScoring.hh
/// \brief Shared data model for LET / track-structure scoring
//
// This header defines the small, self-contained scoring model shared by
// SteppingAction, EventAction, TrackingAction, and HistoManager. Collecting
// the calculation definitions here (rather than scattering more arrays and
// maps across EventAction) keeps the meaning of each accumulated quantity
// explicit and lets the fill sites stay thin.
//
// Terminology (see the LET review, appendix C):
//   Edep  - total energy deposited in a step (G4Step::GetTotalEnergyDeposit)
//   NIEL  - non-ionizing energy deposit reported for the step
//           (G4Step::GetNonIonizingEnergyDeposit). Under the active EM list
//           (G4EmStandardPhysics_option4) this is populated by G4NuclearStopping
//           for ions (alpha/He3/GenericIon -> Li, C, B recoils, capture
//           products) and is ~0 by construction for protons, e+/-, and gammas
//           (proton nuclear stopping is negligible above ~keV). So Eion ~ Edep
//           for light/leptonic deposits, while heavy low-energy recoils carry a
//           real non-ionizing fraction (e.g. ~18% for a 0.2 MeV carbon recoil).
//           Verify per run: Run::EndOfRun prints the NIEL/Edep partition.
//   Eion  - max(0, Edep - NIEL); proxy for energy available to electronic
//           excitation / ionization.
//   Ldep  - Eion / step length; a step- and cut-dependent local deposition
//           density proxy.
//   Lcalc - electronic stopping power from G4EmCalculator at the step-midpoint
//           kinetic energy; smoother and less tied to step segmentation.
//
// Absorber indexing convention (IMPORTANT):
//   Absorbers are 1-based, matching DetectorConstruction (fAbsorThickness[i],
//   fAbsorMaterial[i], GetAbsorberLogical(i) with i in 1..kMaxAbsor-1). Every
//   per-absorber array declared with kMaxAbsor entries here is indexed
//   DIRECTLY by the absorber number (fArr[absNum]); index 0 is reserved and
//   unused. Do not mix in (absNum-1) indexing.

#ifndef LetScoring_h
#define LetScoring_h 1

#include "globals.hh"

#include "G4SystemOfUnits.hh"

#include <cfloat>
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Default pair-creation energy W (estimated) for the active semiconductor.
// This is a proxy, not a tabulated constant like Si's 3.6 eV; it is
// overridable per absorber via DetectorConstruction / the messenger. Reported
// output derived from it (N0 = Eion / W) is an ESTIMATE of initial pairs, not
// collected charge.
constexpr G4double kDefaultPairCreationEnergy = 6.8 * eV;

// LET thresholds (on Lcalc) separating low- from high-LET deposition, used for
// the high-LET energy-fraction observables.
constexpr G4double kLetThreshold1 = 10.0 * keV / micrometer;
constexpr G4double kLetThreshold2 = 100.0 * keV / micrometer;
constexpr G4double kLetThreshold3 = 1000.0 * keV / micrometer;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Depositing-particle category. Stable integer values (persisted in ntuples),
// so append new categories at the end and never renumber. kNbDepositCategories
// must stay in sync with the last enumerator.
enum class DepositCategory : G4int
{
  ElectronGamma = 0,  // e- of gamma-mediated lineage (Compton/photo/pair)
  ElectronIon = 1,    // e- of ion/hadronic-mediated lineage (delta rays)
  Positron = 2,
  Proton = 3,
  LightIon = 4,       // deuteron, triton, He3
  Alpha = 5,
  Lithium = 6,        // Z == 3
  Carbon = 7,         // Z == 6
  BeB = 8,            // Z == 4 or 5 (Be, B)
  OtherHeavyIon = 9,  // any other charged ion (Z > 0)
  NeutralLocal = 10,  // neutral local deposition
  Other = 11
};

constexpr G4int kNbDepositCategories = 12;

// Reaction / interaction channel of the track's lineage. Stable integer values;
// append only.
enum class ReactionCategory : G4int
{
  None = 0,
  HElastic,             // neutron elastic on hydrogen (recoil proton)
  CElastic,             // neutron elastic on carbon (recoil carbon)
  OtherElastic,
  B10CaptureAlphaLi7,   // 10B(n,alpha)7Li
  C12NAlphaBe9,         // 12C(n,alpha)9Be
  C12Breakup3Alpha,     // 12C breakup into 3 alpha (+ n)
  OtherNeutronInelastic,
  GammaCompton,
  GammaPhotoelectric,
  GammaPairProduction,
  Other
};

// Track boundary-escape classification within an active absorber.
enum class EscapeCategory : G4int
{
  Unknown = 0,
  Contained,
  Front,
  Back,
  Lateral,
  Multiple
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Per-track LET accumulators for one absorber. One TrackLetScore is kept per
// absorber for the track currently being transported; flushed to a TrackLET
// ntuple row in TrackingAction::PostUserTrackingAction.
struct TrackLetScore
{
  G4double eDep = 0.;            // sum Edep
  G4double eIon = 0.;           // sum Eion
  G4double niel = 0.;           // sum NIEL
  G4double chargedPath = 0.;    // sum step length over charged steps
  G4double sumEionLetDep = 0.;  // sum(Eion * Ldep) -> energy-weighted Ldep
  G4double sumEionLetCalc = 0.; // sum(Eion * Lcalc) -> energy-weighted Lcalc
  G4double maxLetDep = 0.;      // max step Ldep
  G4double maxLetCalc = 0.;     // max step Lcalc
  G4double sumEionDepth = 0.;   // sum(Eion * depth) -> energy-weighted depth
  G4double minDepth = DBL_MAX;  // shallowest depositing depth
  G4double maxDepth = -DBL_MAX; // deepest depositing depth
  G4double startDepth = 0.;     // depth at first depositing step
  G4double endDepth = 0.;       // depth at last depositing step
  G4double entryKE = -1.;       // kinetic energy on entering the absorber
  G4double exitKE = -1.;        // kinetic energy on leaving / stopping
  G4int nDepositingSteps = 0;
  G4bool hasStart = false;
  G4bool crossedFront = false;  // exited the front face (+depth boundary)
  G4bool crossedBack = false;   // exited the back face
  G4bool crossedLateral = false;// exited a lateral (transverse) boundary
};

// Per-step scoring payload passed from SteppingAction to
// EventAction::AccumulateLetStep. Bundling the fields (rather than a long
// argument list) keeps the call site readable and the model easy to extend.
struct LetStepData
{
  G4double edep = 0.;
  G4double eIon = 0.;
  G4double niel = 0.;
  G4double dl = 0.;            // step length
  G4double letDep = 0.;       // Ldep = eIon / dl (charged steps only)
  G4double letCalc = 0.;      // Lcalc from G4EmCalculator at midpoint KE
  G4double depth = 0.;        // local depth of the step midpoint in the absorber
  G4double preKE = 0.;        // pre-step kinetic energy
  G4double postKE = 0.;       // post-step kinetic energy
  DepositCategory depCat = DepositCategory::Other;
  ReactionCategory rxnCat = ReactionCategory::None;
  G4bool charged = false;
  G4double weight = 1.;
  G4bool enteredAbsor = false;   // first step of this track inside the absorber
  G4bool exitedFront = false;    // step exits via the front (-depth) face
  G4bool exitedBack = false;     // step exits via the back (+depth) face
  G4bool exitedLateral = false;  // step exits via a transverse face
  G4bool stopped = false;        // track stops (KE -> 0) in this step
};

// Per-event LET accumulators for one absorber. One EventLetScore per absorber
// is reset at BeginOfEventAction and written as one EventLET ntuple row per
// active absorber at EndOfEventAction (including zero-deposition events).
struct EventLetScore
{
  G4double eDep = 0.;
  G4double eIon = 0.;
  G4double niel = 0.;
  G4double chargedPath = 0.;
  G4double sumEionLetDep = 0.;
  G4double sumEionLetCalc = 0.;
  G4double maxLetDep = 0.;
  G4double maxLetCalc = 0.;
  G4double eIonAbove1 = 0.;     // Eion on steps with Lcalc >= kLetThreshold1
  G4double eIonAbove2 = 0.;     // ... >= kLetThreshold2
  G4double eIonAbove3 = 0.;     // ... >= kLetThreshold3
  G4double sumEionDepth = 0.;   // sum(Eion * depth)
  G4double sumEionDepth2 = 0.;  // sum(Eion * depth^2) -> depth variance
  G4int nDepositingSteps = 0;
  G4int nChargedTracks = 0;
  std::map<G4int, G4double> eIonByCategory;  // key = DepositCategory value
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Resolve a track's boundary-escape classification from its crossed-face flags.
inline EscapeCategory ResolveEscapeCategory(const TrackLetScore& s)
{
  const G4int n = (s.crossedFront ? 1 : 0) + (s.crossedBack ? 1 : 0)
                  + (s.crossedLateral ? 1 : 0);
  if (n == 0) return EscapeCategory::Contained;
  if (n > 1) return EscapeCategory::Multiple;
  if (s.crossedFront) return EscapeCategory::Front;
  if (s.crossedBack) return EscapeCategory::Back;
  return EscapeCategory::Lateral;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
