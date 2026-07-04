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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "Run.hh"

#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Neutron.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

namespace {

G4bool IsIonHadronicParentParticle(const G4String& name)
{
  const G4String stripped = name.substr(0, name.find('['));
  return stripped == "neutron" || stripped == "proton" || stripped == "alpha"
         || stripped == "Li7" || stripped == "Li6" || stripped == "C12"
         || stripped == "deuteron" || stripped == "triton" || stripped == "He3";
}

G4int ClassifyElectronLineage(const G4Track* track, EventAction* evt)
{
  if (track->GetParentID() <= 0) return kElectronLineageNone;

  const G4int parentId = track->GetParentID();
  const G4int parentLineage = evt->GetElectronLineage(parentId);
  if (parentLineage == kElectronLineageGamma) return kElectronLineageGamma;
  if (parentLineage == kElectronLineageIon) return kElectronLineageIon;

  const G4String parentName = evt->GetTrackParticleName(parentId);
  if (parentName == "gamma") return kElectronLineageGamma;
  if (IsIonHadronicParentParticle(parentName)) return kElectronLineageIon;

  const G4VProcess* creator = track->GetCreatorProcess();
  if (creator == nullptr) return kElectronLineageNone;

  const G4int cat = SecCreatorCategory(creator->GetProcessName());
  if (cat == 2) return kElectronLineageGamma;
  if (cat == 0 || cat == 1) return kElectronLineageIon;
  if (cat == 3) {
    if (parentName == "e-" || parentName == "gamma") return kElectronLineageGamma;
    if (IsIonHadronicParentParticle(parentName)) return kElectronLineageIon;
  }
  return kElectronLineageNone;
}

}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* evt)
  : fDetector(det), fEventAct(evt)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // count processes
  //
  const G4StepPoint* prePoint = aStep->GetPreStepPoint();
  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  const G4VProcess* process = endPoint->GetProcessDefinedStep();
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->CountProcesses(process);

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4Track* track = aStep->GetTrack();
  const G4ParticleDefinition* particle = track->GetDefinition();

  // Register secondary tracks and tag e- lineage on first step (before any
  // absorber Edep so the same step can be split into gamma vs ionic buckets).
  if (track->GetCurrentStepNumber() == 1 && track->GetParentID() > 0) {
    fEventAct->RegisterTrackParticle(track->GetTrackID(), particle->GetParticleName());
    if (particle == G4Electron::Definition()) {
      const G4int lineage = ClassifyElectronLineage(track, fEventAct);
      if (lineage != kElectronLineageNone) {
        fEventAct->TagElectronLineage(track->GetTrackID(), lineage);
      }
    }
  }

  // Front-face monitor (thin slab +Z of calorimeter): first entry from the
  // cavity / source side scores forward-incident spectrum.
  if (particle == G4Neutron::Definition()
      && endPoint->GetStepStatus() == fGeomBoundary)
  {
    G4VPhysicalVolume* postVol = endPoint->GetPhysicalVolume();
    if (postVol != nullptr
        && postVol->GetLogicalVolume() == fDetector->GetMonitorFrontLogical())
    {
      const G4double ke = endPoint->GetKineticEnergy();
      G4int trackId = aStep->GetTrack()->GetTrackID();
      if (!fEventAct->HasMonitorCrossed(trackId)) {
        analysisManager->FillH1(kDetectorFrontKEAll, ke);
        if (aStep->GetTrack()->GetParentID() == 0) {
          analysisManager->FillH1(kDetectorFrontKEPrimary, ke);
        }
        fEventAct->MarkMonitorCrossed(trackId);
      }
    }
  }

  // Gamma field at the detector face: every gamma crossing into the front-face
  // monitor slab (immediately upstream of absorber 1 = first silver layer) is
  // scored once per track. The KE fills the energy-distribution histogram and
  // the run-level crossing counter (integrated for the gamma fluence). Track
  // IDs are unique within an event, so the neutron monitor's track-once gating
  // set is reused without collision.
  if (particle == G4Gamma::Definition()
      && endPoint->GetStepStatus() == fGeomBoundary)
  {
    G4VPhysicalVolume* postVol = endPoint->GetPhysicalVolume();
    if (postVol != nullptr
        && postVol->GetLogicalVolume() == fDetector->GetMonitorFrontLogical())
    {
      G4int trackId = aStep->GetTrack()->GetTrackID();
      if (!fEventAct->HasMonitorCrossed(trackId)) {
        analysisManager->FillH1(kDetectorFrontGammaKE, endPoint->GetKineticEnergy());
        run->CountGammaAtDetectorFace();
        fEventAct->MarkMonitorCrossed(trackId);
      }
    }
  }

  // Identify whether this step's pre-step volume is one of the absorbers.
  // In the new geometry the world is air, the bud-box wall is aluminium,
  // the cavity is air, and the monitor is a thin air slab -- none of which
  // are absorbers. Pointer comparison against the detector's absorber
  // logicals is the only robust way to distinguish them.
  G4VPhysicalVolume* volume = prePoint->GetTouchableHandle()->GetVolume();
  const G4LogicalVolume* preLV =
      (volume != nullptr) ? volume->GetLogicalVolume() : nullptr;

  // Cross-absorber Edep provenance: on a secondary's first step, tag its origin
  // absorber. Origin = the absorber in which the earliest ancestor of the
  // lineage was created. Rule: inherit the parent's origin if it has one;
  // otherwise, if this track is itself born inside an absorber, its origin is
  // that absorber (else it stays 0 = external/primary lineage). preLV at the
  // first step is the volume containing the creation vertex.
  if (track->GetCurrentStepNumber() == 1 && track->GetParentID() > 0) {
    G4int origin = fEventAct->GetTrackOrigin(track->GetParentID());
    if (origin == 0 && preLV != nullptr) {
      for (G4int i = 1; i <= fDetector->GetNbOfAbsor(); ++i) {
        if (preLV == fDetector->GetAbsorberLogical(i)) { origin = i; break; }
      }
    }
    if (origin > 0) fEventAct->SetTrackOrigin(track->GetTrackID(), origin);
  }

  // Moderator-born secondaries: score at first step inside a moderator layer.
  if (track->GetCurrentStepNumber() == 1 && track->GetParentID() > 0 && preLV != nullptr
      && fDetector->IsModeratorLogical(preLV))
  {
    const G4VProcess* creator = track->GetCreatorProcess();
    G4double vtxKE = track->GetVertexKineticEnergy();

    if (particle == G4Gamma::Definition() && creator != nullptr) {
      const G4String& creatorName = creator->GetProcessName();
      if (creatorName == "nCapture" || creatorName == "neutronInelastic") {
        analysisManager->FillH1(kModeratorCaptureGammaKE, vtxKE);
      }
    }

    G4int pi = FixedParticleIdx(particle->GetParticleName());
    if (pi >= 0) {
      analysisManager->FillH1(ModeratorSecESpecId(pi), vtxKE);
      if (track->GetParentID() == 1) {
        analysisManager->FillH1(ModeratorFirstGenSecESpecId(pi), vtxKE);
      }
    }
  }

  G4int absorNum = -1;
  for (G4int i = 1; i <= fDetector->GetNbOfAbsor(); ++i) {
    if (preLV == fDetector->GetAbsorberLogical(i)) { absorNum = i; break; }
  }
  if (absorNum < 0) return;

  // here we are in an absorber. Locate it
  //
  G4int layerNum = prePoint->GetTouchableHandle()->GetCopyNumber(1);

  // Per-track path length in this absorber (all particle types, charged and
  // neutral). Accumulated here and flushed once per track in
  // TrackingAction::PostUserTrackingAction into per-particle, per-absorber
  // path-length histograms and the total-track-length histograms.
  fEventAct->AddTrackPathLength(absorNum, aStep->GetStepLength());

  // collect energy deposit (taking into account track weight)
  G4double edep = aStep->GetTotalEnergyDeposit();
  ////edep *= (aStep->GetTrack()->GetWeight());

  // collect step length of charged particles
  G4double stepl = 0.;
  if (particle->GetPDGCharge() != 0.) stepl = aStep->GetStepLength();

  // sum up per event
  fEventAct->SumEnergy(absorNum, edep, stepl);
  if (track->GetParentID() == 0 && edep > 0.) {
    fEventAct->SumPrimaryEnergy(absorNum, edep);
  }
  if (edep > 0.) {
    fEventAct->SumEnergyByParticle(absorNum, particle->GetParticleName(), edep);
    if (particle == G4Electron::Definition()) {
      const G4int lineage = fEventAct->GetElectronLineage(aStep->GetTrack()->GetTrackID());
      if (lineage != kElectronLineageNone) {
        fEventAct->SumElectronEdepByLineage(absorNum, lineage, edep);
      }
    }
    // Cross-absorber provenance: attribute this deposit to the origin absorber
    // of the depositing track's lineage (0 = external/primary).
    const G4int origin = fEventAct->GetTrackOrigin(track->GetTrackID());
    fEventAct->SumEdepByOrigin(absorNum, origin, edep);
  }

  // longitudinal profile of edep per absorber
  if (edep > 0.) {
    G4AnalysisManager::Instance()->FillH1(kMaxAbsor + absorNum, G4double(layerNum + 1), edep);
  }

  // energy flow
  //
  //  unique identificator of layer+absorber
  G4int Idnow = (fDetector->GetNbOfAbsor()) * layerNum + absorNum;
  G4int plane;
  //
  // leaving the absorber ?
  if (endPoint->GetStepStatus() == fGeomBoundary) {
    // The calorimeter is rotated in the global frame, so the global X axis
    // is no longer the stack direction. Use the local-frame X component of
    // the momentum direction (positive = deeper into the stack).
    G4ThreeVector dirGlobal = endPoint->GetMomentumDirection();
    G4ThreeVector dirLocal =
        prePoint->GetTouchableHandle()->GetHistory()
            ->GetTopTransform().TransformAxis(dirGlobal);
    G4double Eflow = endPoint->GetKineticEnergy();
    if (dirLocal.x() >= 0.)
      run->SumEnergyFlow(plane = Idnow + 1, Eflow);
    else
      run->SumEnergyFlow(plane = Idnow, -Eflow);
  }

  // neutron interaction depth + neutron-energy-at-interaction histograms,
  // and per-event interaction flag for the "fraction interacted" summary.
  //
  if (particle == G4Neutron::Definition() && process != nullptr) {
    G4String procName = process->GetProcessName();
    if (procName != "Transportation") {
      G4int category = -1;
      if (procName == "hadElastic") category = kProcCatHadElastic;
      else if (procName == "neutronInelastic" || procName == "nCapture")
        category = kProcCatNInelastic;

      if (category >= 0) {
        G4ThreeVector localPos = prePoint->GetTouchableHandle()->GetHistory()
            ->GetTopTransform().TransformPoint(endPoint->GetPosition());
        G4double depth = localPos.x() + 0.5 * fDetector->GetAbsorThickness(absorNum);
        analysisManager->FillH1(DepthHistoId(absorNum, category), depth / mm);
      }

      // Any non-Transportation interaction counts towards the per-absorber
      // "fraction interacted" summary and the neutron-KE-at-interaction plot.
      // KE is taken at the pre-step point (just before the interaction).
      // Pass raw internal-unit KE; G4AnalysisManager applies the "eV" unit
      // conversion for display (axis labelled in eV).
      analysisManager->FillH1(NeutronEnergyId(absorNum),
                              prePoint->GetKineticEnergy());

      fEventAct->MarkNeutronInteraction(absorNum);
    }
  }

  // Gamma interaction depth per absorber: Compton scattering (compt),
  // photoelectric absorption (phot), and pair production (conv). Scored at the
  // interaction point (post-step) expressed in the absorber's local frame, so
  // the depth runs 0 (front face) -> thickness (back face) along the stack.
  if (particle == G4Gamma::Definition() && process != nullptr) {
    const G4String& procName = process->GetProcessName();
    G4int category = -1;
    if (procName == "compt") category = kProcCatCompton;
    else if (procName == "phot") category = kProcCatPhotoelectric;
    else if (procName == "conv") category = kProcCatPairProd;

    if (category >= 0) {
      G4ThreeVector localPos = prePoint->GetTouchableHandle()->GetHistory()
          ->GetTopTransform().TransformPoint(endPoint->GetPosition());
      G4double depth = localPos.x() + 0.5 * fDetector->GetAbsorThickness(absorNum);
      analysisManager->FillH1(DepthHistoId(absorNum, category), depth / mm);
    }
  }

  // ----- Per-track first-step plots (capture gamma, secondary spectrum, depth)
  //
  // Many secondary-tracking quantities only need to be filled once per track,
  // at its very first step. Doing it here (rather than in a TrackingAction)
  // keeps the file count down at the cost of one cheap test per step.
  if (track->GetCurrentStepNumber() == 1 && track->GetParentID() > 0) {
    const G4VProcess* creator = track->GetCreatorProcess();
    G4double vtxKE = track->GetVertexKineticEnergy();

    // Capture-gamma KE spectrum: gammas born from a neutron-induced reaction,
    // scored in the absorber where they were produced. Includes both nCapture
    // and neutronInelastic creators because the 478 keV de-excitation line
    // from B-10(n,a)Li-7* is tagged as neutronInelastic in FTFP_BERT_HP.
    if (particle == G4Gamma::Definition() && creator != nullptr) {
      const G4String& creatorName = creator->GetProcessName();
      if (creatorName == "nCapture" || creatorName == "neutronInelastic") {
        analysisManager->FillH1(CaptureGammaId(absorNum), vtxKE);
      }
    }

    // Secondary production depth in absorber (any secondary, any process).
    {
      G4ThreeVector localPos = prePoint->GetTouchableHandle()->GetHistory()
          ->GetTopTransform().TransformPoint(prePoint->GetPosition());
      G4double depth = localPos.x() + 0.5 * fDetector->GetAbsorThickness(absorNum);
      analysisManager->FillH1(SecDepthId(absorNum), depth / mm);

      G4String partName = particle->GetParticleName();
      G4int creatorCat = 4;
      if (creator != nullptr) {
        creatorCat = SecCreatorCategory(creator->GetProcessName());
      }
      run->CountSecondaryBirth(absorNum, partName, creatorCat);
    }

    // KE spectrum of secondary, by particle type (summed across absorbers).
    G4int pi = FixedParticleIdx(particle->GetParticleName());
    if (pi >= 0) {
      analysisManager->FillH1(SecESpecId(pi), vtxKE);
      if (track->GetParentID() == 1) {
        analysisManager->FillH1(FirstGenSecESpecId(pi), vtxKE);
      }

      // e- KE-at-birth split by ancestry (gamma- vs ion/hadronic-mediated),
      // mirroring the e- Edep and path-length lineage splits. The lineage tag
      // is set earlier in this same first step.
      if (pi == 5) {
        const G4int lineage = fEventAct->GetElectronLineage(track->GetTrackID());
        if (lineage == kElectronLineageGamma || lineage == kElectronLineageIon) {
          analysisManager->FillH1(lineage == kElectronLineageGamma
                                    ? kSecESpecElectronGamma
                                    : kSecESpecElectronIonic,
                                  vtxKE);
          if (track->GetParentID() == 1) {
            analysisManager->FillH1(lineage == kElectronLineageGamma
                                      ? kFirstGenSecESpecElectronGamma
                                      : kFirstGenSecESpecElectronIonic,
                                    vtxKE);
          }
        }
      }
    }
  }

  ////  example of Birk attenuation
  /// G4double destep   = aStep->GetTotalEnergyDeposit();
  /// G4double response = BirksAttenuation(aStep);
  /// G4cout << " Destep: " << destep/keV << " keV"
  ///       << " response after Birks: " << response/keV << " keV" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SteppingAction::BirksAttenuation(const G4Step* aStep)
{
  // Example of Birk attenuation law in organic scintillators.
  // adapted from Geant3 PHYS337. See MIN 80 (1970) 239-244
  //
  const G4Material* material = aStep->GetTrack()->GetMaterial();
  G4double birk1 = material->GetIonisation()->GetBirksConstant();
  G4double destep = aStep->GetTotalEnergyDeposit();
  G4double stepl = aStep->GetStepLength();
  G4double charge = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
  //
  G4double response = destep;
  if (birk1 * destep * stepl * charge != 0.) {
    response = destep / (1. + birk1 * destep / stepl);
  }
  return response;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
