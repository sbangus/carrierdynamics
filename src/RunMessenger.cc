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
/// \file RunMessenger.cc
/// \brief Implementation of RunMessenger

#include "RunMessenger.hh"

#include "RunAction.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunMessenger::RunMessenger(RunAction* ra) : fRunAction(ra)
{
  fRunCmdDir = new G4UIdirectory("/testhadr/run/");
  fRunCmdDir->SetGuidance("Run-level commands (batch summaries, etc.)");

  fBatchSummaryEnableCmd =
    new G4UIcmdWithABool("/testhadr/run/setBatchSummaryEnable", this);
  fBatchSummaryEnableCmd->SetGuidance(
    "Enable/disable CSV batch summaries. Output uses setBatchSummaryFile path "
    "when enabled and path is non-empty; default is false.");
  fBatchSummaryEnableCmd->SetParameterName("enable", false);
  fBatchSummaryEnableCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fBatchSummaryEnableCmd->SetToBeBroadcasted(false);

  fBatchSummaryFileCmd =
    new G4UIcmdWithAString("/testhadr/run/setBatchSummaryFile", this);
  fBatchSummaryFileCmd->SetGuidance(
    "CSV append path (does not enable writing by itself). "
    "Enable with setBatchSummaryEnable; empty path disables.");
  fBatchSummaryFileCmd->SetParameterName("path", true);
  fBatchSummaryFileCmd->SetDefaultValue("");
  fBatchSummaryFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fBatchSummaryFileCmd->SetToBeBroadcasted(false);

  fBatchTagCmd = new G4UIcmdWithAString("/testhadr/run/setBatchTag", this);
  fBatchTagCmd->SetGuidance(
    "Label written in each CSV row for this batch (e.g. seed pair)."
    " Set before each /run/beamOn when using multiple batches.");
  fBatchTagCmd->SetParameterName("tag", true);
  fBatchTagCmd->SetDefaultValue("");
  fBatchTagCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fBatchTagCmd->SetToBeBroadcasted(false);

  fBatchCmdDir = new G4UIdirectory("/testhadr/batch/");
  fBatchCmdDir->SetGuidance(
    "Multi-/run/beamOn helpers (distinct ROOT output per run).");

  fRotateRootPerRunCmd =
    new G4UIcmdWithABool("/testhadr/batch/setRotateRootPerRun", this);
  fRotateRootPerRunCmd->SetGuidance(
    "If true, before each run set analysis file name to "
    "<setRootFileBase><runID> (extension from analysis format). "
    "Default false.");
  fRotateRootPerRunCmd->SetParameterName("rotate", false);
  fRotateRootPerRunCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fRotateRootPerRunCmd->SetToBeBroadcasted(false);

  fRootFileBaseCmd =
    new G4UIcmdWithAString("/testhadr/batch/setRootFileBase", this);
  fRootFileBaseCmd->SetGuidance(
    "Stem for per-run ROOT files when setRotateRootPerRun is true "
    "(e.g. hadr05_batch -> hadr05_batch0.root). If empty while rotation is on, "
    "falls back to Hadr05.");
  fRootFileBaseCmd->SetParameterName("stem", false);
  fRootFileBaseCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fRootFileBaseCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunMessenger::~RunMessenger()
{
  delete fRootFileBaseCmd;
  delete fRotateRootPerRunCmd;
  delete fBatchTagCmd;
  delete fBatchSummaryFileCmd;
  delete fBatchSummaryEnableCmd;
  delete fBatchCmdDir;
  delete fRunCmdDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fBatchSummaryEnableCmd) {
    fRunAction->SetBatchSummaryEnable(fBatchSummaryEnableCmd->GetNewBoolValue(newValue));
  }
  else if (command == fBatchSummaryFileCmd) {
    fRunAction->SetBatchSummaryPath(newValue);
  }
  else if (command == fBatchTagCmd) {
    fRunAction->SetBatchTag(newValue);
  }
  else if (command == fRotateRootPerRunCmd) {
    fRunAction->SetRotateRootPerRun(fRotateRootPerRunCmd->GetNewBoolValue(newValue));
  }
  else if (command == fRootFileBaseCmd) {
    fRunAction->SetRootFileBase(newValue);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
