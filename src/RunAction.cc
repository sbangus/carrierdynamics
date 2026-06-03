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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "Run.hh"
#include "RunMessenger.hh"

#include "G4AnalysisManager.hh"
#include "G4GeneralParticleSource.hh"
#include "G4RunManager.hh"
#include "G4SingleParticleSource.hh"
#include "G4SystemOfUnits.hh"
#include "G4Timer.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim)
  : fDetector(det), fPrimary(prim)
{
  fHistoManager = new HistoManager();
  fRunMessenger = new RunMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fRunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetBatchSummaryPath(const G4String& path) { fBatchSummaryPath = path; }

void RunAction::SetBatchSummaryEnable(G4bool enable) { fBatchSummaryEnable = enable; }

void RunAction::SetBatchTag(const G4String& tag) { fBatchTag = tag; }

void RunAction::SetRotateRootPerRun(G4bool rotate) { fRotateRootPerRun = rotate; }

void RunAction::SetRootFileBase(const G4String& base) { fRootFileBase = base; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{
  fRun = new Run(fDetector);
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* run)
{
  // keep run condition
  if (fPrimary) {
    G4GeneralParticleSource* gps = fPrimary->GetGPS();
    G4int nSrc = gps->GetNumberofSource();

    if (nSrc == 1) {
      fRun->SetPrimary(gps->GetParticleDefinition(),
                        gps->GetParticleEnergy());
    }
    else {
      std::ostringstream oss;
      oss << nSrc << " sources:\n";
      G4double totalIntensity = 0.;
      std::vector<G4double> intensities(nSrc);
      for (G4int i = 0; i < nSrc; i++) {
        gps->SetCurrentSourceto(i);
        intensities[i] = gps->GetCurrentSourceIntensity();
        totalIntensity += intensities[i];
      }
      for (G4int i = 0; i < nSrc; i++) {
        gps->SetCurrentSourceto(i);
        G4SingleParticleSource* src = gps->GetCurrentSource();
        G4double pct = 100. * intensities[i] / totalIntensity;
        oss << "     source " << i << ": "
            << src->GetParticleDefinition()->GetParticleName() << " at "
            << G4BestUnit(src->GetEneDist()->GetMonoEnergy(), "Energy")
            << " (" << std::fixed << std::setprecision(0) << pct << "%)";
        if (i < nSrc - 1) oss << "\n";
      }
      fRun->SetPrimaryDescription(oss.str());
      fRun->SetPrimary(gps->GetParticleDefinition(),
                        gps->GetParticleEnergy());
    }
  }

  // histograms
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  if (analysis->IsActive()) {
    if (fRotateRootPerRun && run != nullptr) {
      G4String stem = fRootFileBase;
      if (stem.empty()) {
        stem = "Hadr05";
      }
      analysis->SetFileName(stem + std::to_string(run->GetRunID()));
    }
    analysis->OpenFile();
  }

  // save Rndm status and open the timer

  if (isMaster) {
    //    G4Random::showEngineStatus();
    fTimer = new G4Timer();
    fTimer->Start();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  // compute and print statistic
  if (isMaster) {
    fTimer->Stop();
    if (!((G4RunManager::GetRunManager()->GetRunManagerType() == G4RunManager::sequentialRM))) {
      G4cout << "\n"
             << "Total number of events:  " << fRun->GetNumberOfEvent() << G4endl;
      G4cout << "Master thread time:  " << *fTimer << G4endl;
    }
    delete fTimer;
    if (fRun != nullptr) {
      const G4bool csvOn = fBatchSummaryEnable && !fBatchSummaryPath.empty();
      fRun->SetBatchCsvExport(csvOn, fBatchSummaryPath, fBatchTag);
      fRun->EndOfRun();
    }
  }
  // save histograms
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  if (analysis->IsActive()) {
    analysis->Write();
    analysis->CloseFile();
  }

  // show Rndm status
  //  if (isMaster)  G4Random::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
