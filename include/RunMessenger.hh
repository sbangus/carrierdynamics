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
/// \file RunMessenger.hh
/// \brief UI commands for optional per-run batch CSV summaries and ROOT rotation

#ifndef RunMessenger_h
#define RunMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class RunAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunMessenger : public G4UImessenger
{
  public:
    explicit RunMessenger(RunAction*);
    ~RunMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    RunAction* fRunAction = nullptr;

    G4UIdirectory* fRunCmdDir = nullptr;
    G4UIdirectory* fBatchCmdDir = nullptr;
    G4UIcmdWithABool* fBatchSummaryEnableCmd = nullptr;
    G4UIcmdWithAString* fBatchSummaryFileCmd = nullptr;
    G4UIcmdWithAString* fBatchTagCmd = nullptr;
    G4UIcmdWithABool* fRotateRootPerRunCmd = nullptr;
    G4UIcmdWithAString* fRootFileBaseCmd = nullptr;
};

#endif
