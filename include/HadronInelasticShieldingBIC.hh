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
/// \file HadronInelasticShieldingBIC.hh
/// \brief Definition of the HadronInelasticShieldingBIC class

#ifndef HadronInelasticShieldingBIC_h
#define HadronInelasticShieldingBIC_h 1

#include "G4HadronPhysicsShielding.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Hadron inelastic physics that keeps the full G4HadronPhysicsShielding
// neutron treatment (ParticleHP inelastic / radiative capture / fission with
// ENDF cross sections below 20 MeV) but replaces the proton builder with the
// QGSP_BIC_AllHP proton chain: data-driven ParticleHP proton inelastic below
// 200 MeV, Binary cascade up to 1.5 GeV, then Bertini / FTFP / QGSP at higher
// energies. This gives high-precision protons (the dominant charged secondary
// from fast-neutron recoils and (n,p) reactions in the organic absorbers)
// without disturbing the neutron transport, which stays on Shielding HP.
//
// Pions, kaons and other hadrons remain the inherited FTFP_BERT/Shielding
// builders (high energy only; irrelevant for a DD-neutron / gamma source).
class HadronInelasticShieldingBIC : public G4HadronPhysicsShielding
{
  public:
    explicit HadronInelasticShieldingBIC(G4int verbose = 1);
    ~HadronInelasticShieldingBIC() override = default;

    HadronInelasticShieldingBIC(HadronInelasticShieldingBIC&) = delete;
    HadronInelasticShieldingBIC& operator=(const HadronInelasticShieldingBIC& right) = delete;

  protected:
    // Mirror of G4HadronPhysicsQGSP_BIC_AllHP::Proton().
    void Proton() override;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
