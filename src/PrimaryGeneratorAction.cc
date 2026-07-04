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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4RunManager.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SingleParticleSource.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det) : fDetector(det)
{
  fGPS = new G4GeneralParticleSource();

  // Default: DD-fusion neutron source at the configured source standoff on
  // the +Z exterior side of the bud box, emitting a 5 deg cone of mono
  // 2.45 MeV neutrons aimed at the detector. Macros override any/all of
  // these via /gps/... commands; nothing here is enforced per-event.
  G4SingleParticleSource* src = fGPS->GetCurrentSource();

  src->SetParticleDefinition(
    G4ParticleTable::GetParticleTable()->FindParticle("neutron"));

  src->GetEneDist()->SetEnergyDisType("Mono");
  src->GetEneDist()->SetMonoEnergy(2.45 * MeV);

  src->GetPosDist()->SetPosDisType("Point");
  src->GetPosDist()->SetCentreCoords(fDetector->GetDefaultSourcePosition());

  // 5 deg cone aimed along -Z (toward the box / detector). For GPS "iso"
  // angular distribution, theta is measured such that theta = 0 emits along
  // -AngRef3, with AngRef3 = AngRef1 x AngRef2 (see G4SPSAngDistribution::
  // GenerateIsotropicFlux). Choosing rot1 = +X, rot2 = +Y gives AngRef3 = +Z,
  // hence theta = 0 emits along -Z.
  G4SPSAngDistribution* ang = src->GetAngDist();
  ang->SetAngDistType("iso");
  ang->DefineAngRefAxes("angref1", G4ThreeVector(1., 0., 0.));
  ang->DefineAngRefAxes("angref2", G4ThreeVector(0., 1., 0.));
  ang->SetMinTheta(0. * deg);
  ang->SetMaxTheta(5. * deg);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGPS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fGPS->GeneratePrimaryVertex(anEvent);

  // Record the actual sampled primary (type + kinetic energy) into the current
  // (thread-local) run. This is the reliable source of the run's primary energy
  // for the EndOfRun summary; GPS::GetParticleEnergy() read in BeginOfRunAction
  // returns the GPS default because no vertex has been sampled yet.
  const G4PrimaryVertex* vertex = anEvent->GetPrimaryVertex(0);
  if (vertex != nullptr) {
    const G4PrimaryParticle* primary = vertex->GetPrimary(0);
    if (primary != nullptr) {
      auto* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
      if (run != nullptr) {
        run->AddPrimary(primary->GetParticleDefinition(),
                        primary->GetKineticEnergy());
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
