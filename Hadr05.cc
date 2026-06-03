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
/// \file Hadr05.cc
/// \brief Main program of the hadronic/Hadr05 example

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "HadronInelasticShieldingBIC.hh"

#include "G4DeexPrecoParameters.hh"
#include "G4EmParameters.hh"
#include "G4IonPhysicsPHP.hh"
#include "G4NuclearLevelData.hh"
#include "G4NuclideTable.hh"
#include "G4ParticleHPManager.hh"
#include "G4PhysListFactory.hh"
#include "G4ProductionCutsTable.hh"
#include "G4RunManagerFactory.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4SteppingVerbose.hh"
#include "G4SystemOfUnits.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

#include <cmath>
#include <cstdlib>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

namespace
{
  // The QGSP_BIC_AllHP proton / light-ion libraries (HadronInelasticShieldingBIC
  // proton override + G4IonPhysicsPHP) drive G4ParticleHPInelastic for protons,
  // deuterons, tritons, He3 and alphas. Those models read the G4TENDL data set,
  // located via G4PARTICLEHPDATA (umbrella) or the per-particle *HPDATA
  // variables. G4NDL (G4NEUTRONHPDATA) only covers neutrons. If the TENDL data
  // is not visible, Geant4 aborts deep inside ConstructProcess with a terse
  // "setenv G4PARTICLEHPDATA" message; warn early and clearly instead.
  void CheckParticleHPData()
  {
    const bool haveUmbrella = (std::getenv("G4PARTICLEHPDATA") != nullptr);
    const char* perParticle[] = {"G4PROTONHPDATA", "G4DEUTERONHPDATA", "G4TRITONHPDATA",
                                 "G4HE3HPDATA", "G4ALPHAHPDATA"};
    bool haveAllPerParticle = true;
    for (const char* v : perParticle) {
      if (std::getenv(v) == nullptr) haveAllPerParticle = false;
    }
    if (haveUmbrella || haveAllPerParticle) return;

    // NB: do NOT setenv() G4PARTICLEHPDATA here. Calling setenv() after Geant4
    // has started can reallocate the process environment and invalidate the
    // G4NEUTRONHPDATA pointer that ParticleHP cached at startup, which leads to
    // spurious "data does not exist in G4NEUTRONHPDATA" errors and intermittent
    // crashes in MT mode. The variable must be exported in the shell before the
    // program starts (see README) -- just warn here.
    G4cout
      << "\n*************************************************************************\n"
      << "*  WARNING: ParticleHP data for protons / light ions was not found.     *\n"
      << "*                                                                       *\n"
      << "*  This build uses the QGSP_BIC_AllHP proton and light-ion libraries,   *\n"
      << "*  which require the G4TENDL data set. Set ONE of:                      *\n"
      << "*     export G4PARTICLEHPDATA=/path/to/G4TENDL<version>                 *\n"
      << "*  (or the per-particle G4PROTONHPDATA / G4DEUTERONHPDATA /             *\n"
      << "*   G4TRITONHPDATA / G4HE3HPDATA / G4ALPHAHPDATA variables).            *\n"
      << "*                                                                       *\n"
      << "*  Without it Geant4 will abort during /run/initialize. G4NDL           *\n"
      << "*  (neutrons) alone is not sufficient for the proton/ion HP models.     *\n"
      << "*************************************************************************\n"
      << G4endl;
  }
}  // namespace

int main(int argc, char** argv)
{
  // detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Use SteppingVerbose with Unit
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  // Creating run manager
  auto runManager = G4RunManagerFactory::CreateRunManager();

  if (argc == 3) {
    G4int nThreads = G4UIcommand::ConvertToInt(argv[2]);
    runManager->SetNumberOfThreads(nThreads);
  }

  // Half-life threshold for the nuclide table — must be set BEFORE the
  // physics list is constructed so that the radioactive decay constructor
  // sees it. Kept here (was previously in the retired custom PhysicsList).
  const G4double meanLife = 1 * nanosecond;
  const G4double halfLife = meanLife * std::log(2.);
  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(halfLife);

  // Detector
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);

  // Physics list: Shielding_HPT_EMZ (Geant4 reference modular list)
  //
  //   - Shielding: HP neutrons below 20 MeV (elastic, inelastic, capture,
  //     fission), Bertini/FTFP above; tuned for shielding and neutron
  //     transport applications.
  //   - HPT: NeutronHP thermal scattering (S(alpha,beta)) for neutrons
  //     below ~4 eV. Use TS-named elements in materials (e.g. water_TS,
  //     polyethylene_TS in DetectorConstruction) so thermal kernels bind.
  //   - EMZ: G4EmStandardPhysics_option4 — highest-accuracy standard EM
  //     package in Geant4 (Livermore/Penelope low-energy models, advanced
  //     msc, full atomic deexcitation). This is the EM configuration used
  //     by the reference _EMZ lists; there is no higher-precision *standard*
  //     EM alternative in the factory.
  //
  // For neutrons specifically, Geant4 documents the HP/ParticleHP data
  // as the precision option; Shielding + HPT is the reference combination
  // for thermal + fast neutrons. LEND-based lists (e.g. ShieldingLEND_HP)
  // trade different evaluated data for some isotopes — not universally
  // "higher precision" — so Shielding_HPT_EMZ remains the requested default.
  G4PhysListFactory factory;
  G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("Shielding_HPT_EMZ");
  physicsList->SetVerboseLevel(1);

  // High-precision proton + light-ion transport in the absorbers, taken from
  // QGSP_BIC_AllHP while leaving the Shielding HP neutron treatment (and the
  // HPT thermal-scattering elastic) untouched:
  //
  //   - HadronInelasticShieldingBIC inherits G4HadronPhysicsShielding (so the
  //     neutron inelastic / capture / fission HP models are kept verbatim) and
  //     overrides only Proton() with the QGSP_BIC_AllHP proton chain: data-
  //     driven ParticleHP proton inelastic below 200 MeV + Binary cascade.
  //   - G4IonPhysicsPHP provides ParticleHP inelastic for d, t, He3, alpha
  //     below 200 MeV/n + Binary cascade, i.e. the light-ion library from
  //     QGSP_BIC_AllHP.
  //
  // ReplacePhysics swaps these in by physics type (hadron-inelastic, ions),
  // so only the proton and light-ion libraries change; everything else in
  // Shielding_HPT_EMZ (neutron HP, HPT thermal scattering, EMZ / option4 EM,
  // gamma/lepto-nuclear, decay, ...) is preserved.
  physicsList->ReplacePhysics(new HadronInelasticShieldingBIC(1));
  physicsList->ReplacePhysics(new G4IonPhysicsPHP(1));

  // Enforce the per-volume G4UserLimits max-step set on the absorber logical
  // volumes (DetectorConstruction::DefineRegionsAndCuts). The Shielding_HPT_EMZ
  // reference list does NOT register a step limiter, so without this the
  // 0.1 um absorber step caps would be silently ignored.
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());

  // The AllHP proton/light-ion models need the G4TENDL data set; warn clearly
  // if it is not reachable (otherwise Geant4 aborts cryptically at init).
  CheckParticleHPData();

  // Lower the production-cut table energy floor (default 990 eV) so the very
  // small absorber range cuts actually resolve to distinct, low energies.
  // Without this, the 0.1 um (absorber) and 0.01 um (metal gamma) range cuts
  // both clamp to the 990 eV floor and the tighter metal gamma cut has no
  // effect. 100 eV matches the EMZ / Livermore low-energy table limit.
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100 * eV, 10 * GeV);

  // Baseline range cuts for lab / air / water volumes without a dedicated
  // G4Region (see DetectorConstruction::DefineRegionsAndCuts() for thin-film
  // absorbers, monitor, box wall, moderator).
  physicsList->SetDefaultCutValue(10.0 * mm);
  runManager->SetUserInitialization(physicsList);

  // EM step-function tuning previously applied by the local
  // ElectromagneticPhysics helper — EmStandardPhysics_option4 (EMZ) does
  // not set these, so re-apply for ion / heavy-charged-particle stepping.
  G4EmParameters* emParam = G4EmParameters::Instance();
  emParam->SetStepFunction(0.2, 100 * um);
  emParam->SetStepFunctionMuHad(0.1, 10 * um);
  emParam->SetStepFunctionLightIons(0.1, 10 * um);
  emParam->SetStepFunctionIons(0.1, 1 * um);
  emParam->SetDeexcitationIgnoreCut(true);
  // Atomic relaxation: fluorescence + Auger + Auger-cascade + PIXE so the
  // K/L line emission from photoelectric absorption in the high-Z metal
  // contacts and oxide layers is produced (dominant low-energy gamma channel).
  emParam->SetFluo(true);
  emParam->SetAuger(true);
  emParam->SetAugerCascade(true);
  emParam->SetPixe(true);
  // Mott correction improves e- backscatter in the high-Z metal contacts.
  emParam->SetUseMottCorrection(true);

  // Nuclear de-excitation: emit correlated gamma cascades (proper line
  // multiplicities / angular correlations) following neutron radiative
  // capture and inelastic excitation. This directly improves the (n,gamma)
  // gamma-line structure that is the focus of the gamma study.
  G4DeexPrecoParameters* deex = G4NuclearLevelData::GetInstance()->GetParameters();
  deex->SetCorrelatedGamma(true);
  deex->SetStoreAllLevels(true);

  // set user action classes
  runManager->SetUserInitialization(new ActionInitialization(detector));

  // Replaced HP environmental variables with C++ calls
  G4ParticleHPManager::GetInstance()->SetSkipMissingIsotopes(false);
  G4ParticleHPManager::GetInstance()->SetDoNotAdjustFinalState(true);
  G4ParticleHPManager::GetInstance()->SetUseOnlyPhotoEvaporation(true);
  G4ParticleHPManager::GetInstance()->SetNeglectDoppler(false);
  G4ParticleHPManager::GetInstance()->SetProduceFissionFragments(false);
  G4ParticleHPManager::GetInstance()->SetUseWendtFissionModel(false);
  G4ParticleHPManager::GetInstance()->SetUseNRESP71Model(true);

  // initialize visualization
  G4VisManager* visManager = nullptr;

  // get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (ui) {
    // interactive mode
    visManager = new G4VisExecutive();
    visManager->Initialize();
    ui->SessionStart();
    delete ui;
  }
  else {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }

  // job termination
  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
