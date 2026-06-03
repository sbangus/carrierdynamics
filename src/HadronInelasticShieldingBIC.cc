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
/// \file HadronInelasticShieldingBIC.cc
/// \brief Implementation of the HadronInelasticShieldingBIC class

#include "HadronInelasticShieldingBIC.hh"

#include "G4BertiniProtonBuilder.hh"
#include "G4BinaryProtonBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4HadronicParameters.hh"
#include "G4ParticleHPInelastic.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4ParticleInelasticXS.hh"
#include "G4ProcessManager.hh"
#include "G4Proton.hh"
#include "G4QGSPProtonBuilder.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadronInelasticShieldingBIC::HadronInelasticShieldingBIC(G4int verbose)
  : G4HadronPhysicsShielding("hInelastic ShieldingBIC", verbose)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadronInelasticShieldingBIC::Proton()
{
  // Reproduces G4HadronPhysicsQGSP_BIC_AllHP::Proton(). The energy boundaries
  // are declared locally because this class derives from
  // G4HadronPhysicsShielding (FTFP_BERT branch), whose inherited proton members
  // carry FTFP_BERT semantics rather than the QGSP_BIC ladder reproduced here:
  //   ParticleHP : 0      -> 200 MeV   (data-driven inelastic)
  //   Binary     : 190 MeV -> 1.5 GeV
  //   Bertini    : 1.0 GeV -> maxFTFP_Cascade transition
  //   FTFP       : minFTFP_Cascade transition -> maxQGS_FTF transition
  //   QGSP       : minQGS_FTF transition -> maxQGS_FTF transition
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  const G4bool useFactorXS = param->ApplyFactorXS();

  const G4double maxHP_proton = 200.0 * MeV;
  const G4double minBIC_proton = 190.0 * MeV;
  const G4double maxBIC_proton = 1.5 * GeV;
  const G4double minBERT_proton = 1.0 * GeV;
  const G4double maxBERT_proton = param->GetMaxEnergyTransitionFTF_Cascade();
  const G4double minFTFP_proton = param->GetMinEnergyTransitionFTF_Cascade();
  const G4double maxFTFP_proton = param->GetMaxEnergyTransitionQGS_FTF();
  const G4double minQGSP_proton = param->GetMinEnergyTransitionQGS_FTF();
  const G4bool quasiElasticQGS = true;
  const G4bool quasiElasticFTF = false;

  G4ParticleDefinition* proton = G4Proton::Proton();
  auto* inel = new G4HadronInelasticProcess("protonInelastic", proton);
  proton->GetProcessManager()->AddDiscreteProcess(inel);

  G4QGSPProtonBuilder qgs(quasiElasticQGS);
  qgs.SetMinEnergy(minQGSP_proton);
  qgs.Build(inel);

  G4FTFPProtonBuilder ftf(quasiElasticFTF);
  ftf.SetMinEnergy(minFTFP_proton);
  ftf.SetMaxEnergy(maxFTFP_proton);
  ftf.Build(inel);

  if (maxBERT_proton > minBERT_proton) {
    G4BertiniProtonBuilder bert;
    bert.SetMinEnergy(minBERT_proton);
    bert.SetMaxEnergy(maxBERT_proton);
    bert.Build(inel);
  }

  if (maxBIC_proton > 0.0) {
    G4BinaryProtonBuilder bic;
    bic.SetMinEnergy(minBIC_proton);
    bic.SetMaxEnergy(maxBIC_proton);
    bic.Build(inel);
  }

  inel->AddDataSet(new G4ParticleInelasticXS(proton));

  inel->AddDataSet(new G4ParticleHPInelasticData(proton));
  auto* mod = new G4ParticleHPInelastic(proton, "ProtonHPInelastic");
  mod->SetMaxEnergy(maxHP_proton);
  inel->RegisterMe(mod);

  if (useFactorXS) {
    inel->MultiplyCrossSectionBy(param->XSFactorNucleonInelastic());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
