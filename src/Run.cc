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
/// \file Run.cc
/// \brief Implementation of the Run class

#include "Run.hh"

#include "DetectorConstruction.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"
#include "G4UnitsTable.hh"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>

namespace {

// Per-particle log binning for secondary KE-at-birth spectra (absorber and
// moderator). Ranges tuned to typical DD / moderation / (n,gamma) physics.
void ConfigureSecondaryESpec(G4AnalysisManager* analysis, G4int ih, G4int partIdx)
{
  G4int nbins = 100;
  G4double vmin = 1.e-2 * MeV;
  G4double vmax = 2.e1 * MeV;
  switch (partIdx) {
    case 4:  // gamma
      nbins = 200;
      vmin = 1.e-3 * MeV;
      vmax = 5. * MeV;
      break;
    case 5:  // e-
      nbins = 150;
      vmin = 1.e-3 * MeV;
      vmax = 20. * MeV;
      break;
    case 3:  // proton
      nbins = 100;
      vmin = 50. * keV;
      vmax = 10. * MeV;
      break;
    case 0:
    case 1:
    case 2:  // alpha, Li7, Li6
      nbins = 150;
      vmin = 100. * keV;
      vmax = 5. * MeV;
      break;
    case 6:  // C12
      nbins = 100;
      vmin = 1. * MeV;
      vmax = 20. * MeV;
      break;
    default:
      break;
  }
  analysis->SetH1(ih, nbins, vmin, vmax, "MeV", "none", "log");
}

// Upper bound for per-absorber path-length histos: 10x layer thickness, floored
// at 10 nm and capped at 1 mm so thin layers are not drowned by empty axis.
G4double AbsorberTrackLengthVmax(G4double thickness)
{
  constexpr G4double kCap = 1. * mm;
  constexpr G4double kFloor = 10. * nm;
  return std::min(kCap, std::max(10. * thickness, kFloor));
}

// Per-particle log vmin/nbins for path-length histos (bounds in mm).
void TrackLengthParticleBinning(G4int partIdx, G4int& nbins, G4double& vmin)
{
  nbins = 120;
  vmin = 1.e-3 * mm;  // 1 um
  switch (partIdx) {
    case 0:
    case 1:
    case 2:
    case 6:  // alpha, Li7, Li6, C12 (heavy ions: short, sub-mm ranges)
      nbins = 150;
      vmin = 1.e-5 * mm;  // 10 nm
      break;
    case 3:  // proton
      nbins = 120;
      vmin = 1.e-4 * mm;  // 100 nm
      break;
    case 4:  // gamma
      nbins = 120;
      vmin = 1.e-2 * mm;  // 10 um
      break;
    case 5:  // e-
      nbins = 150;
      // 1 um default collides with AbsorberTrackLengthVmax on 100 nm layers
      // (vmax capped at 1 um), which leaves g4tools unable to stream the histo.
      vmin = 1.e-5 * mm;  // 10 nm
      break;
    default:
      break;
  }
}

// Observable A: path length accumulated inside one absorber (per track).
void ConfigureAbsorberTrackLength(G4AnalysisManager* analysis, G4int ih,
                                  G4int partIdx, G4double absorberThickness)
{
  G4int nbins = 120;
  G4double vmin = 1.e-3 * mm;
  TrackLengthParticleBinning(partIdx, nbins, vmin);

  G4double vmax = AbsorberTrackLengthVmax(absorberThickness);
  constexpr G4double kMinLogSpan = 1.e-2;
  if (vmin >= vmax) {
    vmin = vmax * kMinLogSpan;
  }

  analysis->SetH1(ih, nbins, vmin, vmax, "mm", "none", "log");
}

// Observable B: full track length across all volumes (per track).
void ConfigureTotalTrackLength(G4AnalysisManager* analysis, G4int ih, G4int partIdx)
{
  G4int nbins = 120;
  G4double vmin = 1.e-3 * mm;
  G4double vmax = 1.e2 * mm;  // 10 cm
  TrackLengthParticleBinning(partIdx, nbins, vmin);

  switch (partIdx) {
    case 0:
    case 1:
    case 2:
    case 6:
      vmax = 0.1 * mm;  // 100 um
      break;
    case 3:
      vmax = 1. * mm;
      break;
    case 4:
      vmax = 5.e2 * mm;  // 50 cm
      break;
    case 5:
      vmax = 50. * mm;
      break;
    default:
      break;
  }

  analysis->SetH1(ih, nbins, vmin, vmax, "mm", "none", "log");
}

// e- path length split by ancestry (gamma- vs ion/hadronic-mediated). Explicit
// log bounds in mm; the caller picks the range appropriate to the lineage.
void ConfigureElectronLineagePath(G4AnalysisManager* analysis, G4int ih,
                                  G4double vmin, G4double vmax, G4int nbins)
{
  if (vmin >= vmax) vmin = vmax * 1.e-3;
  analysis->SetH1(ih, nbins, vmin, vmax, "mm", "none", "log");
}

// Primary-particle path length inside one absorber (per track, any primary type).
void ConfigurePrimaryAbsorberPath(G4AnalysisManager* analysis, G4int ih,
                                  G4double absorberThickness)
{
  G4int nbins = 150;
  G4double vmin = 1.e-5 * mm;
  G4double vmax = AbsorberTrackLengthVmax(absorberThickness);
  constexpr G4double kMinLogSpan = 1.e-2;
  if (vmin >= vmax) vmin = vmax * kMinLogSpan;
  analysis->SetH1(ih, nbins, vmin, vmax, "mm", "none", "log");
}

// Primary-particle total track length across all volumes (per track).
void ConfigurePrimaryTotalPath(G4AnalysisManager* analysis, G4int ih)
{
  analysis->SetH1(ih, 150, 1.e-5 * mm, 1.e3 * mm, "mm", "none", "log");
}

}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det) : fDetector(det)
{
  // initialize energy deposited per absorber
  //
  for (G4int k = 0; k < kMaxAbsor; k++) {
    fSumEAbs[k] = fSum2EAbs[k] = fSumLAbs[k] = fSum2LAbs[k] = 0.;
  }

  // initialize total energy deposited
  //
  fEdepTot = fEdepTot2 = 0.;

  // initialize leakage
  //
  fEnergyLeak[0] = fEnergyLeak[1] = 0.;
  fEleakTot = fEleakTot2 = 0.;

  // initialize total energy released
  //
  fEtotal = fEtotal2 = 0.;

  // initialize Eflow
  //
  G4int nbPlanes = (fDetector->GetNbOfLayers()) * (fDetector->GetNbOfAbsor()) + 2;
  fEnergyFlow.resize(nbPlanes);
  for (G4int k = 0; k < nbPlanes; k++) {
    fEnergyFlow[k] = 0.;
  }

  // auto-configure and activate per-absorber histograms
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  const G4int nAbs = fDetector->GetNbOfAbsor();

  for (G4int k = 1; k <= nAbs; k++) {
    // Edep per absorber (per event)
    analysis->SetH1(k, 100, 0., 1.0, "MeV");
    analysis->SetH1Activation(k, true);

    // Longitudinal Edep profile per absorber (scaled to MeV/event in EndOfRun)
    {
      G4int ih = kMaxAbsor + k;  // IDs 21..40 when kMaxAbsor=20
      analysis->SetH1(ih, 100, 0., 100., "MeV");
      analysis->SetH1Activation(ih, true);
    }

    G4double thickness = fDetector->GetAbsorThickness(k);
    G4String matName = fDetector->GetAbsorMaterial(k)->GetName();

    // Neutron interaction depth histograms (elastic, inelastic+capture)
    for (G4int c = 0; c < kNbProcessCat; c++) {
      G4int ih = DepthHistoId(k, c);
      analysis->SetH1(ih, 100, 0., thickness / mm, "mm");
      analysis->SetH1Activation(ih, true);
    }

    // Per-particle Edep per absorber (MeV)
    for (G4int p = 0; p < kNbFixedParticles; p++) {
      G4int ih = EdepByParticleId(k, p);
      analysis->SetH1(ih, 100, 0., 2.5, "MeV");
      analysis->SetH1Activation(ih, true);
    }

    // Capture-gamma KE spectrum: dominant B-10(n,a)Li-7* line is at 478 keV,
    // so 0..2 MeV with 200 bins gives ~10 keV resolution.
    {
      G4int ih = CaptureGammaId(k);
      analysis->SetH1(ih, 200, 0., 2.0, "MeV");
      analysis->SetH1Activation(ih, true);
    }

    // EHP per event: only meaningful for the polymer (PNDI-T10*) layers,
    // where charge generation actually happens. Activated only for those.
    {
      G4int ih = EhpId(k);
      analysis->SetH1(ih, 200, 0., 200000., "none");  // counts of EHP / event
      G4bool isPNDI = (matName.find("PNDI") != std::string::npos);
      analysis->SetH1Activation(ih, isPNDI);
    }

    // Secondary production depth in absorber (mm)
    {
      G4int ih = SecDepthId(k);
      analysis->SetH1(ih, 100, 0., thickness / mm, "mm");
      analysis->SetH1Activation(ih, true);
    }

    // Neutron KE at interaction (log binning, 1 meV .. 14 MeV, labelled in eV).
    // Bounds are passed as physics quantities in G4 internal units so that
    // after G4AnalysisManager's per-fill unit conversion the displayed axis
    // reads in eV and matches the (un-divided) fill value from SteppingAction.
    {
      G4int ih = NeutronEnergyId(k);
      analysis->SetH1(ih, 100, 1.e-3 * eV, 1.4e7 * eV, "eV", "none", "log");
      analysis->SetH1Activation(ih, true);
    }
  }

  // Energy spectrum of each secondary type (one per particle, summed across
  // absorbers). Per-particle log ranges tuned in ConfigureSecondaryESpec().
  for (G4int p = 0; p < kNbFixedParticles; p++) {
    G4int ih = SecESpecId(p);
    ConfigureSecondaryESpec(analysis, ih, p);
    analysis->SetH1Activation(ih, true);
  }

  // KE at generation of first-generation secondaries only (ParentID == 1),
  // one histogram per fixed particle type, summed across absorbers.
  for (G4int p = 0; p < kNbFixedParticles; p++) {
    G4int ih = FirstGenSecESpecId(p);
    ConfigureSecondaryESpec(analysis, ih, p);
    analysis->SetH1Activation(ih, true);
  }

  // e- KE-at-birth spectrum split by ancestry, the KE analog of the e- Edep and
  // path-length lineage splits. The two lineages span different energy scales
  // (investigation, 30 um PNDI-T10-B4C, mixed n/gamma source):
  //   - gamma-mediated: ~30 eV .. ~2 MeV (Compton/photo e- off the 0.5/2 MeV
  //     gammas); median ~150 eV with an energetic tail.
  //   - ion/hadronic-mediated: delta rays off recoil ions; peak ~0.27 keV,
  //     essentially all below ~4 keV.
  // vmin is dropped to 10 eV so the sub-keV bulk (otherwise lost to underflow
  // with the 1 keV all-e- ConfigureSecondaryESpec floor) is resolved. The
  // first-generation ion-mediated histogram is expected to stay ~empty:
  // primaries are gammas/neutrons and neutrons do not directly produce e-.
  analysis->SetH1(kSecESpecElectronGamma, 150, 1.e-5 * MeV, 3. * MeV, "MeV", "none", "log");
  analysis->SetH1Activation(kSecESpecElectronGamma, true);
  analysis->SetH1(kSecESpecElectronIonic, 150, 1.e-5 * MeV, 1.e-2 * MeV, "MeV", "none", "log");
  analysis->SetH1Activation(kSecESpecElectronIonic, true);
  analysis->SetH1(kFirstGenSecESpecElectronGamma, 150, 1.e-5 * MeV, 3. * MeV, "MeV", "none", "log");
  analysis->SetH1Activation(kFirstGenSecESpecElectronGamma, true);
  analysis->SetH1(kFirstGenSecESpecElectronIonic, 150, 1.e-5 * MeV, 1.e-2 * MeV, "MeV", "none", "log");
  analysis->SetH1Activation(kFirstGenSecESpecElectronIonic, true);

  // Neutron spectrum entering the detector after environment interactions.
  // Log binning from 1 meV to 20 MeV, axis labelled in eV (matches the
  // existing per-absorber NeutronEnergyId convention so SteppingAction can
  // pass raw internal-unit kinetic energies).
  analysis->SetH1(kDetectorFrontKEAll, 100,
                  1.e-3 * eV, 2.e7 * eV, "eV", "none", "log");
  analysis->SetH1Activation(kDetectorFrontKEAll, true);
  analysis->SetH1(kDetectorFrontKEPrimary, 100,
                  1.e-3 * eV, 2.e7 * eV, "eV", "none", "log");
  analysis->SetH1Activation(kDetectorFrontKEPrimary, true);

  // Gamma energy spectrum at the detector front face (front-face monitor slab,
  // immediately upstream of absorber 1 = first silver layer). Log binning from
  // 1 keV to 20 MeV (covers capture/inelastic de-excitation lines), axis in
  // MeV; SteppingAction fills it with the raw internal-unit kinetic energy.
  analysis->SetH1(kDetectorFrontGammaKE, 100,
                  1. * keV, 2.e1 * MeV, "MeV", "none", "log");
  analysis->SetH1Activation(kDetectorFrontGammaKE, true);

  // Per-event e- Edep split by ancestry (gamma- vs ion/hadronic-mediated). Both
  // span a wide dynamic range and are dominated by low-energy deposits, so log
  // binning is used (vs. the linear all-e- histogram EdepByParticleId, id 150).
  // Investigation (30 um PNDI-T10-B4C, mixed n/gamma source):
  //   - gamma-mediated: peak ~0.5 keV, tail to ~0.15 MeV (Compton/photo e- that
  //     mostly traverse the thin layer depositing little).
  //   - ion/hadronic-mediated: peak ~23 keV (B-10(n,alpha)Li-7 delta-ray sum),
  //     tail to ~0.5 MeV.
  for (G4int k = 1; k <= nAbs; k++) {
    analysis->SetH1(EdepElectronGammaId(k), 100, 1.e-5 * MeV, 0.2 * MeV, "MeV", "none", "log");
    analysis->SetH1Activation(EdepElectronGammaId(k), true);
    analysis->SetH1(EdepElectronIonicId(k), 100, 1.e-5 * MeV, 0.5 * MeV, "MeV", "none", "log");
    analysis->SetH1Activation(EdepElectronIonicId(k), true);
  }

  // Moderator-born secondary spectra (activated only when a stack is placed).
  const G4bool hasModerators = fDetector->HasModeratorStack();
  analysis->SetH1(kModeratorCaptureGammaKE, 300, 0., 3.0, "MeV");
  analysis->SetH1Activation(kModeratorCaptureGammaKE, hasModerators);
  for (G4int p = 0; p < kNbFixedParticles; p++) {
    ConfigureSecondaryESpec(analysis, ModeratorSecESpecId(p), p);
    analysis->SetH1Activation(ModeratorSecESpecId(p), hasModerators);
    ConfigureSecondaryESpec(analysis, ModeratorFirstGenSecESpecId(p), p);
    analysis->SetH1Activation(ModeratorFirstGenSecESpecId(p), hasModerators);
  }

  // Secondary track path length per absorber, per particle type (observable A),
  // and total track length per particle across all volumes (observable B).
  for (G4int k = 1; k <= nAbs; k++) {
    const G4double thickness = fDetector->GetAbsorThickness(k);
    for (G4int p = 0; p < kNbFixedParticles; p++) {
      G4int ih = TrackPathLengthId(k, p);
      ConfigureAbsorberTrackLength(analysis, ih, p, thickness);
      analysis->SetH1Activation(ih, true);
    }
  }
  for (G4int p = 0; p < kNbFixedParticles; p++) {
    G4int ih = TotalTrackLengthId(p);
    ConfigureTotalTrackLength(analysis, ih, p);
    analysis->SetH1Activation(ih, true);
  }

  // e- path length split by ancestry (gamma- vs ion/hadronic-mediated), the
  // path-length analog of the e- Edep lineage split. The two lineages live on
  // very different length scales, so they get different ranges:
  //   - gamma-mediated: energetic Compton/photo/pair e- (um .. cm), scaled like
  //     the all-e- path-length histograms.
  //   - ion/hadronic-mediated: low-energy delta rays from heavily-ionizing
  //     recoil ions (alpha/Li7/proton). Ranges peak near ~6 nm, so the axis is
  //     fixed at the nm scale (1 nm .. 1 um per absorber, 1 nm .. 100 um total).
  // vmin is dropped to 1 nm for all four so the sub-micron delta-ray population
  // (otherwise lost to underflow) is resolved.
  const G4double kElectronPathVmin = 1.e-6 * mm;  // 1 nm
  for (G4int k = 1; k <= nAbs; k++) {
    const G4double thickness = fDetector->GetAbsorThickness(k);
    ConfigureElectronLineagePath(analysis, PathElectronGammaId(k),
                                 kElectronPathVmin, AbsorberTrackLengthVmax(thickness), 160);
    analysis->SetH1Activation(PathElectronGammaId(k), true);
    ConfigureElectronLineagePath(analysis, PathElectronIonicId(k),
                                 kElectronPathVmin, 1.e-3 * mm, 150);
    analysis->SetH1Activation(PathElectronIonicId(k), true);
  }
  ConfigureElectronLineagePath(analysis, kTotalPathElectronGamma,
                               kElectronPathVmin, 50. * mm, 160);
  analysis->SetH1Activation(kTotalPathElectronGamma, true);
  ConfigureElectronLineagePath(analysis, kTotalPathElectronIonic,
                               kElectronPathVmin, 1.e-1 * mm, 150);
  analysis->SetH1Activation(kTotalPathElectronIonic, true);

  // Primary-particle path length (per track) and Edep (per event, primary steps).
  for (G4int k = 1; k <= nAbs; k++) {
    ConfigurePrimaryAbsorberPath(analysis, PrimaryPathLengthId(k),
                                 fDetector->GetAbsorThickness(k));
    analysis->SetH1Activation(PrimaryPathLengthId(k), true);
    analysis->SetH1(PrimaryEdepId(k), 100, 0., 20., "MeV");
    analysis->SetH1Activation(PrimaryEdepId(k), true);
  }
  ConfigurePrimaryTotalPath(analysis, kTotalPrimaryPathLength);
  analysis->SetH1Activation(kTotalPrimaryPathLength, true);

  // Run-level scalar histograms: energy flow, total Edep, leakage, released.
  // IDs 41..44 (= 2*kMaxAbsor+1 .. 2*kMaxAbsor+4 when kMaxAbsor=20).
  {
    G4int base = 2 * kMaxAbsor;
    for (G4int k = 1; k <= 4; k++) {
      analysis->SetH1(base + k, 100, 0., 100., "MeV");
      analysis->SetH1Activation(base + k, true);
    }
  }

  // Per-absorber summary histograms (one bin per absorber, centered at k).
  // Configured to hold up to nAbs bins; bin contents are written in EndOfRun.
  if (nAbs > 0) {
    analysis->SetH1(kInteractionFractionId, nAbs, 0.5, nAbs + 0.5, "none");
    analysis->SetH1Activation(kInteractionFractionId, true);
    analysis->SetH1(kComparativeResponseId, nAbs, 0.5, nAbs + 0.5, "MeV");
    analysis->SetH1Activation(kComparativeResponseId, true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{
  fParticle = particle;
  fEkin = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimaryDescription(const G4String& desc)
{
  fPrimaryDesc = desc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountProcesses(const G4VProcess* process)
{
  if (process == nullptr) return;
  G4String procName = process->GetProcessName();
  std::map<G4String, G4int>::iterator it = fProcCounter.find(procName);
  if (it == fProcCounter.end()) {
    fProcCounter[procName] = 1;
  }
  else {
    fProcCounter[procName]++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::CountSecondaryBirth(G4int absIdx, const G4String& particleName, G4int creatorCat)
{
  if (absIdx < 1 || absIdx > kMaxAbsor) return;
  const G4int k = absIdx - 1;
  fSecBirthByPart[k][particleName]++;
  if (creatorCat >= 0 && creatorCat < kNbSecCreatorCat) {
    fSecBirthByCreatorCat[k][creatorCat]++;
  }
  else {
    fSecBirthByCreatorCat[k][4]++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumEdepPerAbsorber(G4int kAbs, G4double EAbs, G4double LAbs)
{
  // accumulate statistic with restriction
  //
  fSumEAbs[kAbs] += EAbs;
  fSum2EAbs[kAbs] += EAbs * EAbs;
  fSumLAbs[kAbs] += LAbs;
  fSum2LAbs[kAbs] += LAbs * LAbs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumEdepByParticle(G4int kAbs, const std::map<G4String, G4double>& edepMap)
{
  for (const auto& entry : edepMap) {
    fEdepByParticle[kAbs][entry.first] += entry.second;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumEnergies(G4double edeptot, G4double eleak0, G4double eleak1)
{
  fEdepTot += edeptot;
  fEdepTot2 += edeptot * edeptot;

  fEnergyLeak[0] += eleak0;
  fEnergyLeak[1] += eleak1;
  G4double eleaktot = eleak0 + eleak1;
  fEleakTot += eleaktot;
  fEleakTot2 += eleaktot * eleaktot;

  G4double etotal = edeptot + eleaktot;
  fEtotal += etotal;
  fEtotal2 += etotal * etotal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SumEnergyFlow(G4int plane, G4double Eflow)
{
  fEnergyFlow[plane] += Eflow;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddInteractedFlags(const G4bool flags[])
{
  // Each event sets flags[k] = true if a primary neutron underwent any
  // non-Transportation interaction in absorber k during this event.
  for (G4int k = 0; k < kMaxAbsor; k++) {
    if (flags[k]) ++fNeutronInteractedEvents[k];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin = localRun->fEkin;
  fPrimaryDesc = localRun->fPrimaryDesc;

  // accumulate sums
  //
  for (G4int k = 0; k < kMaxAbsor; k++) {
    fSumEAbs[k] += localRun->fSumEAbs[k];
    fSum2EAbs[k] += localRun->fSum2EAbs[k];
    fSumLAbs[k] += localRun->fSumLAbs[k];
    fSum2LAbs[k] += localRun->fSum2LAbs[k];
  }

  fEdepTot += localRun->fEdepTot;
  fEdepTot2 += localRun->fEdepTot2;

  fEnergyLeak[0] += localRun->fEnergyLeak[0];
  fEnergyLeak[1] += localRun->fEnergyLeak[1];

  fEleakTot += localRun->fEleakTot;
  fEleakTot2 += localRun->fEleakTot2;

  fEtotal += localRun->fEtotal;
  fEtotal2 += localRun->fEtotal2;

  G4int nbPlanes = (fDetector->GetNbOfLayers()) * (fDetector->GetNbOfAbsor()) + 2;
  for (G4int k = 0; k < nbPlanes; k++) {
    fEnergyFlow[k] += localRun->fEnergyFlow[k];
  }

  // map: per-particle edep
  for (G4int k = 0; k < kMaxAbsor; k++) {
    for (const auto& entry : localRun->fEdepByParticle[k]) {
      fEdepByParticle[k][entry.first] += entry.second;
    }
  }

  // per-absorber interaction-event counters
  for (G4int k = 0; k < kMaxAbsor; k++) {
    fNeutronInteractedEvents[k] += localRun->fNeutronInteractedEvents[k];
  }

  // gammas crossing the detector front face
  fGammaAtDetectorFace += localRun->fGammaAtDetectorFace;

  // map: processes count
  std::map<G4String, G4int>::const_iterator itp;
  for (itp = localRun->fProcCounter.begin(); itp != localRun->fProcCounter.end(); ++itp) {
    G4String procName = itp->first;
    G4int localCount = itp->second;
    if (fProcCounter.find(procName) == fProcCounter.end()) {
      fProcCounter[procName] = localCount;
    }
    else {
      fProcCounter[procName] += localCount;
    }
  }

  for (G4int k = 0; k < kMaxAbsor; k++) {
    for (const auto& entry : localRun->fSecBirthByPart[k]) {
      fSecBirthByPart[k][entry.first] += entry.second;
    }
    for (const auto& entry : localRun->fSecBirthByCreatorCat[k]) {
      fSecBirthByCreatorCat[k][entry.first] += entry.second;
    }
  }

  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{
  // run condition
  //
  if (!fPrimaryDesc.empty()) {
    G4cout << "\n ---> The run is " << numberOfEvent << " events, source: "
           << fPrimaryDesc << G4endl;
  }
  else if (fParticle) {
    G4cout << "\n ---> The run is " << numberOfEvent << " "
           << fParticle->GetParticleName() << " of "
           << G4BestUnit(fEkin, "Energy") << " through calorimeter" << G4endl;
  }

  // frequency of processes
  //
  G4cout << "\n Process calls frequency :" << G4endl;
  G4int index = 0;
  std::map<G4String, G4int>::iterator it;
  for (it = fProcCounter.begin(); it != fProcCounter.end(); it++) {
    G4String procName = it->first;
    G4int count = it->second;
    G4String space = " ";
    if (++index % 3 == 0) space = "\n";
    G4cout << " " << std::setw(22) << procName << "=" << std::setw(10) << count << space;
  }

  G4cout << G4endl;
  G4int nEvt = numberOfEvent;
  G4double norm = G4double(nEvt);
  if (norm > 0) norm = 1. / norm;
  G4double qnorm = std::sqrt(norm);

  // energy deposit per absorber
  //
  G4double beamEnergy = fEkin;
  G4double sqbeam = std::sqrt(beamEnergy / GeV);

  G4double MeanEAbs, MeanEAbs2, rmsEAbs, resolution, rmsres;
  G4double MeanLAbs, MeanLAbs2, rmsLAbs;

  std::ios::fmtflags mode = G4cout.flags();
  G4int prec = G4cout.precision(2);
  G4cout << "\n------------------------------------------------------------\n";
  G4cout << std::setw(16) << "material" << std::setw(22) << "Edep        rmsE" << std::setw(31)
         << "sqrt(E0(GeV))*rmsE/Edep" << std::setw(23) << "total tracklen \n \n";

  for (G4int k = 1; k <= fDetector->GetNbOfAbsor(); k++) {
    MeanEAbs = fSumEAbs[k] * norm;
    MeanEAbs2 = fSum2EAbs[k] * norm;
    rmsEAbs = std::sqrt(std::abs(MeanEAbs2 - MeanEAbs * MeanEAbs));

    resolution = 100. * sqbeam * rmsEAbs / MeanEAbs;
    rmsres = resolution * qnorm;

    MeanLAbs = fSumLAbs[k] * norm;
    MeanLAbs2 = fSum2LAbs[k] * norm;
    rmsLAbs = std::sqrt(std::abs(MeanLAbs2 - MeanLAbs * MeanLAbs));

    // print
    //
    G4cout << std::setw(2) << k << std::setw(14) << fDetector->GetAbsorMaterial(k)->GetName()
           << std::setprecision(5) << std::setw(10) << G4BestUnit(MeanEAbs, "Energy")
           << std::setprecision(4) << std::setw(8) << G4BestUnit(rmsEAbs, "Energy") << std::setw(10)
           << resolution << " +- " << std::setprecision(3) << std::setw(5) << rmsres << " %"
           << std::setprecision(4) << std::setw(12) << G4BestUnit(MeanLAbs, "Length") << " +- "
           << std::setprecision(3) << std::setw(5) << G4BestUnit(rmsLAbs, "Length") << G4endl;
  }

  // per-particle energy deposit breakdown
  //
  G4cout << "\n------------------------------------------------------------\n";
  G4cout << " Energy deposited by particle type (per event):\n";
  for (G4int k = 1; k <= fDetector->GetNbOfAbsor(); k++) {
    if (fEdepByParticle[k].empty()) continue;
    G4cout << "\n   Absorber " << k << " ("
           << fDetector->GetAbsorMaterial(k)->GetName() << "):\n";
    for (const auto& entry : fEdepByParticle[k]) {
      G4double meanEdep = entry.second * norm;
      G4cout << "     " << std::setw(16) << entry.first << " : "
             << std::setprecision(4) << G4BestUnit(meanEdep, "Energy") << "\n";
    }
  }

  G4cout << "\n Total energy deposited by particle type (sum over all events), "
            "by absorber layer:\n";
  for (G4int k = 1; k <= fDetector->GetNbOfAbsor(); k++) {
    if (fEdepByParticle[k].empty()) continue;
    G4cout << "\n   Absorber " << k << " ("
           << fDetector->GetAbsorMaterial(k)->GetName() << "):\n";
    for (const auto& entry : fEdepByParticle[k]) {
      const G4double totalEdep = entry.second;
      G4cout << "     " << std::setw(16) << entry.first << " : "
             << std::setprecision(4) << G4BestUnit(totalEdep, "Energy") << "\n";
    }
  }
  G4cout << "------------------------------------------------------------\n";

  PrintSecondaryBirthSummary();

  // total energy deposited
  //
  fEdepTot /= nEvt;
  fEdepTot2 /= nEvt;
  G4double rmsEdep = std::sqrt(std::abs(fEdepTot2 - fEdepTot * fEdepTot));

  G4cout << "\n Total energy deposited = " << std::setprecision(4) << G4BestUnit(fEdepTot, "Energy")
         << " +- " << G4BestUnit(rmsEdep, "Energy") << G4endl;

  // Energy leakage
  //
  fEnergyLeak[0] /= nEvt;
  fEnergyLeak[1] /= nEvt;
  fEleakTot /= nEvt;
  fEleakTot2 /= nEvt;
  G4double rmsEleak = std::sqrt(std::abs(fEleakTot2 - fEleakTot * fEleakTot));

  G4cout << " Leakage :  primary = " << G4BestUnit(fEnergyLeak[0], "Energy")
         << "   secondaries = " << G4BestUnit(fEnergyLeak[1], "Energy")
         << "  ---> total = " << G4BestUnit(fEleakTot, "Energy") << " +- "
         << G4BestUnit(rmsEleak, "Energy") << G4endl;

  // total energy released
  //
  fEtotal /= nEvt;
  fEtotal2 /= nEvt;
  G4double rmsEtotal = std::sqrt(std::abs(fEtotal2 - fEtotal * fEtotal));

  G4cout << " Total energy released :  Edep + Eleak = " << G4BestUnit(fEtotal, "Energy") << " +- "
         << G4BestUnit(rmsEtotal, "Energy") << G4endl;
  G4cout << "------------------------------------------------------------\n";

  // Gamma field at the detector front face (front-face monitor slab, placed
  // immediately upstream of absorber 1 = the first silver layer). The energy
  // distribution is stored in histogram kDetectorFrontGammaKE; here we report
  // the integrated flux (fluence) from the gamma crossing count and the
  // monitor face area.
  {
    const G4double areaCm2 = fDetector->GetCalorTransverseArea() / cm2;
    const G4double meanPerEvent =
        (nEvt > 0) ? G4double(fGammaAtDetectorFace) / G4double(nEvt) : 0.;

    G4cout << " Gamma flux at detector front face (just upstream of absorber 1):\n";
    G4cout << "   gammas crossing face    = " << fGammaAtDetectorFace << "\n";
    G4cout << "   per primary event       = " << std::setprecision(4)
           << meanPerEvent << "\n";
    G4cout << "   monitor face area       = " << std::setprecision(4)
           << areaCm2 << " cm2\n";
    if (areaCm2 > 0.) {
      G4cout << "   fluence (whole run)     = " << std::setprecision(4)
             << G4double(fGammaAtDetectorFace) / areaCm2 << " /cm2\n";
      G4cout << "   fluence per event       = " << std::setprecision(4)
             << meanPerEvent / areaCm2 << " /cm2/event\n";
    }
    G4cout << "   energy spectrum         : histogram id "
           << kDetectorFrontGammaKE << " (1 keV - 20 MeV, log) in Hadr05.root\n";
    G4cout << "------------------------------------------------------------\n";
  }

  if (fDetector->HasModeratorStack()) {
    G4cout << " Moderator-born secondary spectra (KE at birth in moderator stack):\n";
    G4cout << "   capture-gamma line hist : id " << kModeratorCaptureGammaKE
           << " (0-3 MeV linear, 300 bins)\n";
    G4cout << "   all secondary types     : ids " << kModeratorSecESpecBase << "-"
           << (kModeratorSecESpecBase + kNbFixedParticles - 1)
           << " (per-particle log ranges)\n";
    G4cout << "   first-generation only   : ids " << kModeratorFirstGenSecESpecBase
           << "-" << (kModeratorFirstGenSecESpecBase + kNbFixedParticles - 1)
           << " (per-particle log ranges)\n";
    G4cout << "------------------------------------------------------------\n";
  }

  // Energy flow
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  G4int Idmax = (fDetector->GetNbOfLayers()) * (fDetector->GetNbOfAbsor());
  for (G4int Id = 1; Id <= Idmax + 1; Id++) {
    analysis->FillH1(2 * kMaxAbsor + 1, (G4double)Id, fEnergyFlow[Id] / nEvt);
  }

  // normalize histograms
  //
  for (G4int ih = kMaxAbsor + 1; ih < 2 * kMaxAbsor + 1; ih++) {
    analysis->ScaleH1(ih, norm / MeV);
  }

  // Fill summary histograms (one bin per absorber). FillH1(id, x, w) deposits
  // weight w at coordinate x; with bin centers at 1..nAbs, this puts the
  // value into the right bin exactly once per absorber.
  for (G4int k = 1; k <= fDetector->GetNbOfAbsor(); k++) {
    G4double meanEdep = fSumEAbs[k] * norm;  // MeV / event
    analysis->FillH1(kComparativeResponseId, (G4double)k, meanEdep / MeV);

    G4double frac = (nEvt > 0)
        ? G4double(fNeutronInteractedEvents[k]) / G4double(nEvt)
        : 0.;
    analysis->FillH1(kInteractionFractionId, (G4double)k, frac);
  }

  AppendBatchSummaryCsv();

  // remove all contents in fProcCounter and fEdepByParticle
  fProcCounter.clear();
  for (G4int k = 0; k < kMaxAbsor; k++) {
    fEdepByParticle[k].clear();
    fSecBirthByPart[k].clear();
    fSecBirthByCreatorCat[k].clear();
    fNeutronInteractedEvents[k] = 0;
  }
  fGammaAtDetectorFace = 0;

  G4cout.setf(mode, std::ios::floatfield);
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetBatchCsvExport(G4bool enable, const G4String& path, const G4String& tag)
{
  fBatchCsvEnable = enable;
  fBatchCsvPath = path;
  fBatchCsvTag = tag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::PrintSecondaryBirthSummary() const
{
  G4cout << "\n Secondary track births (first step in absorber, ParentID > 0):\n";
  for (G4int k = 1; k <= fDetector->GetNbOfAbsor(); k++) {
    const G4int idx = k - 1;
    G4long total = 0;
    for (const auto& entry : fSecBirthByPart[idx]) {
      total += entry.second;
    }
    if (total == 0 && fSecBirthByCreatorCat[idx].empty()) continue;

    G4cout << "\n   Absorber " << k << " ("
           << fDetector->GetAbsorMaterial(k)->GetName() << "): " << total
           << " secondary tracks"
           << "  [histogram SecDepthId = " << SecDepthId(k) << "]\n";

    if (!fSecBirthByCreatorCat[idx].empty()) {
      G4cout << "     By creator category:\n";
      for (G4int c = 0; c < kNbSecCreatorCat; c++) {
        const auto it = fSecBirthByCreatorCat[idx].find(c);
        if (it == fSecBirthByCreatorCat[idx].end() || it->second == 0) continue;
        G4cout << "       " << std::setw(10) << SecCreatorCategoryName(c) << " : "
               << it->second << "\n";
      }
    }

    if (!fSecBirthByPart[idx].empty()) {
      G4cout << "     By particle (top contributors):\n";
      for (const auto& entry : fSecBirthByPart[idx]) {
        G4cout << "       " << std::setw(10) << entry.first << " : " << entry.second
               << "\n";
      }
    }
  }
  G4cout << "------------------------------------------------------------\n";
}

namespace
{
G4String CsvEscape(const G4String& s)
{
  G4String out = "\"";
  for (std::size_t i = 0; i < s.size(); ++i) {
    if (s[i] == '"') {
      out += "\"\"";
    }
    else {
      out += s[i];
    }
  }
  out += "\"";
  return out;
}
}  // namespace

void Run::AppendBatchSummaryCsv()
{
  if (!fBatchCsvEnable || fBatchCsvPath.empty()) return;
  if (!G4Threading::IsMasterThread()) return;

  const G4long nEvt = numberOfEvent;
  if (nEvt <= 0) return;

  static const char* kExpectedHdr =
    "run_id,n_events,batch_tag,absorber,material,particle,total_edep_MeV,"
    "gammas_crossing_face";

  std::ifstream probe(fBatchCsvPath.c_str());
  G4bool writeHeader = true;
  if (probe.good()) {
    std::string firstLine;
    std::getline(probe, firstLine);
    while (!firstLine.empty() && firstLine.back() == '\r') {
      firstLine.pop_back();
    }
    if (!firstLine.empty()) {
      if (firstLine == kExpectedHdr) {
        writeHeader = false;
      }
      else {
        G4cerr << "\n *** Run::AppendBatchSummaryCsv: file \"" << fBatchCsvPath
               << "\" exists but the first line is not this example's CSV header."
               << "\n     Delete/rename the file or use a new path. Skipping append."
               << G4endl;
        return;
      }
    }
  }

  std::ofstream out(fBatchCsvPath.c_str(), std::ios::app);
  if (!out.good()) {
    G4cerr << "\n *** Run::AppendBatchSummaryCsv: cannot open \"" << fBatchCsvPath
           << "\" for append." << G4endl;
    return;
  }

  out << std::scientific << std::setprecision(17);
  if (writeHeader) {
    out << "run_id,n_events,batch_tag,absorber,material,particle,total_edep_MeV,"
           "gammas_crossing_face\n";
  }

  const G4int runId = GetRunID();
  const G4double invMeV = 1. / MeV;

  // The gamma-crossing count is a single per-batch scalar; it is repeated on
  // every (absorber, particle) row of this batch so each batch's value can be
  // read directly or de-duplicated when pivoting the long-format table.
  for (G4int k = 1; k <= fDetector->GetNbOfAbsor(); ++k) {
    const G4String& matName = fDetector->GetAbsorMaterial(k)->GetName();
    for (const auto& entry : fEdepByParticle[k]) {
      const G4double totalMeV = entry.second * invMeV;
      out << runId << ',' << nEvt << ',' << CsvEscape(fBatchCsvTag) << ',' << k << ','
          << CsvEscape(matName) << ',' << CsvEscape(entry.first) << ',' << totalMeV << ','
          << fGammaAtDetectorFace << '\n';
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
