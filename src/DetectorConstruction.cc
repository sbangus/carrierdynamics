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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "DetectorMessenger.hh"

#include <cmath>

#include "G4Box.hh"
#include "G4Ellipsoid.hh"
#include "G4GeometryManager.hh"
#include "G4Isotope.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Polyhedra.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4RunManager.hh"
#include "G4SolidStore.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"

#include <algorithm>
#include <iomanip>
#include <string>

namespace
{
  // Metallic contact / electrode layers, which receive an even tighter gamma
  // production cut (0.01 um) so photoelectric K/L lines from high-Z metals are
  // resolved. Covers both NIST (G4_*) and locally-defined metal materials.
  G4bool IsMetalAbsorberMaterial(const G4Material* m)
  {
    if (m == nullptr) return false;
    const G4String& n = m->GetName();
    return (n == "G4_Ag" || n == "G4_Au" || n == "G4_Cu" || n == "G4_Al" || n == "G4_Fe"
            || n == "G4_W" || n == "G4_Pb" || n == "G4_Ti" || n == "Gold" || n == "Copper"
            || n == "Aluminium" || n == "Iron" || n == "Tungsten" || n == "Lead"
            || n == "Titanium");
  }
}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  for (G4int i = 0; i < kMaxAbsor; ++i) {
    fAbsorMaterial[i] = nullptr;
    fAbsorThickness[i] = 0.0;
    fSolidAbsor[i] = nullptr;
    fLogicAbsor[i] = nullptr;
    fPhysiAbsor[i] = nullptr;
  }

  // default parameter values of the calorimeter
  fNbOfAbsor = 4;
  fAbsorThickness[1] = 40 * nm;
  fAbsorThickness[2] = 75 * micrometer;
  fAbsorThickness[3] = 100 * nm;
  fAbsorThickness[4] = 100 * micrometer;
  fNbOfLayers = 1;
  fCalorSizeYZ = 1 * mm;
  fCalorSizeY = 1 * mm;
  fCalorSizeZ = 1 * mm;
  ComputeCalorParameters();

  // Bud box (aluminum) external dimensions and wall thickness.
  // The detector lies parallel to a 103x53 face; the stack normal is along
  // global Z, and the front face is fDetectorStandoff away from the inside
  // of the +Z face.
  fBoxExternalX = 103. * mm;
  fBoxExternalY = 53. * mm;
  fBoxExternalZ = 25. * mm;
  fBoxWallThickness = 1.8 * mm;
  fDetectorStandoff = 10. * mm;

  // Surrounding lab/air world: large enough to comfortably contain external
  // sources at typical calibration standoffs without leakage artefacts at
  // the world boundary.
  fLabSizeX = 2000. * mm;
  fLabSizeY = 2000. * mm;
  fLabSizeZ = 2000. * mm;

  // Default source standoff from the box +Z exterior face. Used by
  // PrimaryGeneratorAction when no GPS position is given by macro.
  fSourceStandoff = 100. * mm;

  // Silver-epoxy contact blob (localized half-ellipsoid) defaults. Only used
  // when /testhadr/det/setSilverEpoxyBlob true; otherwise absorber 1 is a
  // full-area box like the other layers.
  fEpoxyBlobRadius = 2. * mm;
  fEpoxyBlobHeight = 1. * mm;

  // materials
  DefineMaterials();
  SetWorldMaterial("G4_AIR");
  SetCavityMaterial("G4_AIR");
  SetBoxMaterial("G4_Al");
  SetAbsorMaterial(1, "G4_Ag");
  SetAbsorMaterial(2, "PNDI-T10");
  SetAbsorMaterial(3, "ITO");
  SetAbsorMaterial(4, "QuartzGlass");

  // Default moderator material (legacy single-slab fallback when
  // setModeratorEnable true and no explicit layer stack is configured).
  fModeratorMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
  fModeratorUpstreamGap = 1. * mm;
  for (G4int i = 0; i <= kMaxModeratorLayers; ++i) {
    fModeratorLayerMat[i] = nullptr;
    fModeratorLayerThick[i] = 0.;
    fSolidModeratorLayer[i] = nullptr;
    fLogicModeratorLayer[i] = nullptr;
    fPhysiModeratorLayer[i] = nullptr;
  }

  // Polycarbonate encasement defaults (shine-dd-box); geometry off unless macro
  // enables it.
  fEncasementOuterX = 336.3 * mm;
  fEncasementOuterY = 336.3 * mm;
  fEncasementOuterZ = 183.39 * mm;
  fEncasementWallThickness = 4.6 * mm;

  // create commands for interactive definition of the calorimeter
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
  delete fCalorRotation;
  delete fEpoxyBlobRotation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // This function illustrates the possible ways to define materials using
  // G4 database on G4Elements
  G4NistManager* manager = G4NistManager::Instance();
  manager->SetVerbose(0);
  //
  // define Elements
  //
  G4Element* H = manager->FindOrBuildElement(1);
  G4Element* C = manager->FindOrBuildElement(6);
  G4Element* N = manager->FindOrBuildElement(7);
  G4Element* O = manager->FindOrBuildElement(8);
  G4Element* Si = manager->FindOrBuildElement(14);
  G4Element* S = manager->FindOrBuildElement(16);
  G4Element* In = manager->FindOrBuildElement(49);
  G4Element* Sn = manager->FindOrBuildElement(50);

  // Bound hydrogen for thermal neutron scattering kernels (HPT).
  // The element name "TS_H_of_Polyethylene" is recognized by
  // G4NeutronHPThermalScattering (Shielding_HPT_* lists), which then uses the
  // polyethylene S(alpha,beta) tables for neutrons below ~4 eV in any
  // material that contains this element. (FindOrBuildElement returns nullptr
  // for these special TS names, so we construct the element directly — same
  // pattern as Hadr04 / Hadr06 / Hadr05Adapted.)
  G4Element* H_TS_Poly =
    new G4Element("TS_H_of_Polyethylene", "H", 1., 1.0079 * g / mole);
  //
  // define an Element from isotopes, by relative abundance
  //
  G4int iz, n;  // iz=number of protons  in an isotope;
                //  n=number of nucleons in an isotope;
  G4int ncomponents;
  G4double z, a;
  G4double abundance;

  G4Isotope* U5 = new G4Isotope("U235", iz = 92, n = 235, a = 235.01 * g / mole);
  G4Isotope* U8 = new G4Isotope("U238", iz = 92, n = 238, a = 238.03 * g / mole);

  G4Element* U = new G4Element("enriched Uranium", "U", ncomponents = 2);
  U->AddIsotope(U5, abundance = 90. * perCent);
  U->AddIsotope(U8, abundance = 10. * perCent);

  G4Isotope* B10 = new G4Isotope("B10", iz = 5, n = 10, a = 10.0129 * g / mole);
  G4Isotope* B11 = new G4Isotope("B11", iz = 5, n = 11, a = 11.0093 * g / mole);

  G4Element* enrichedB = new G4Element("enriched Boron", "B_enr", ncomponents = 2);
  enrichedB->AddIsotope(B10, abundance = 96. * perCent);
  enrichedB->AddIsotope(B11, abundance = 4. * perCent);

  //
  // define simple materials
  //
  G4double density;

  new G4Material("liquidH2", z = 1., a = 1.008 * g / mole, density = 70.8 * mg / cm3);
  new G4Material("Aluminium", z = 13., a = 26.98 * g / mole, density = 2.700 * g / cm3);
  new G4Material("liquidArgon", z = 18, a = 39.948 * g / mole, density = 1.396 * g / cm3);
  new G4Material("Titanium", z = 22., a = 47.867 * g / mole, density = 4.54 * g / cm3);
  new G4Material("Iron", z = 26., a = 55.85 * g / mole, density = 7.870 * g / cm3);
  new G4Material("Copper", z = 29., a = 63.55 * g / mole, density = 8.960 * g / cm3);
  new G4Material("Tungsten", z = 74., a = 183.85 * g / mole, density = 19.30 * g / cm3);
  new G4Material("Gold", z = 79., a = 196.97 * g / mole, density = 19.32 * g / cm3);
  new G4Material("Lead", z = 82., a = 207.20 * g / mole, density = 11.35 * g / cm3);
  new G4Material("Uranium", z = 92., a = 238.03 * g / mole, density = 18.95 * g / cm3);

  //
  // define a material from elements.   case 1: chemical molecule
  //
  G4int natoms;

  G4Material* H2O = new G4Material("Water", density = 1.000 * g / cm3, ncomponents = 2);
  H2O->AddElement(H, natoms = 2);
  H2O->AddElement(O, natoms = 1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0 * eV);
  H2O->SetChemicalFormula("H_2O");

  G4Material* CH = new G4Material("Polystyrene", density = 1.032 * g / cm3, ncomponents = 2);
  CH->AddElement(C, natoms = 1);
  CH->AddElement(H, natoms = 1);

  G4Material* Sci = new G4Material("Scintillator", density = 1.032 * g / cm3, ncomponents = 2);
  Sci->AddElement(C, natoms = 9);
  Sci->AddElement(H, natoms = 10);

  G4Material* B4C = new G4Material("enrichedB4C", density = 2.52 * g / cm3, ncomponents = 2);
  B4C->AddElement(enrichedB, natoms = 4);
  B4C->AddElement(C, natoms = 1);

  Sci->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);

  // PNDI-T10: (C62H88N2O4S2)_0.9 . (C58H86N2O4S)_0.1
  // Effective composition by mass fraction from weighted atomic counts:
  //   0.9*(C62 H88 N2 O4 S2) + 0.1*(C58 H86 N2 O4 S)
  //   C: 0.9*62 + 0.1*58 = 61.6,  H: 0.9*88 + 0.1*86 = 87.8
  //   N: 0.9*2  + 0.1*2  = 2.0,   O: 0.9*4  + 0.1*4  = 4.0
  //   S: 0.9*2  + 0.1*1  = 1.9
  G4double totalMassPNDI = 61.6 * 12.011 + 87.8 * 1.008 + 2.0 * 14.007
                         + 4.0 * 15.999 + 1.9 * 32.065;
  G4Material* PNDI = new G4Material("PNDI-T10", density = 1.5 * g / cm3, ncomponents = 5);
  PNDI->AddElement(C, 61.6 * 12.011 / totalMassPNDI);
  // Documented proxy: PNDI-T10 has no dedicated S(alpha,beta) library, so we
  // wire its hydrogen as TS_H_of_Polyethylene. This is justified because both
  // are saturated aliphatic / aromatic-with-alkyl-side-chain hydrogenous
  // polymers; the polyethylene kernel is the closest available approximation
  // for the bound-hydrogen vibrational/rotational modes below ~4 eV.
  PNDI->AddElement(H_TS_Poly, 87.8 * 1.008 / totalMassPNDI);
  PNDI->AddElement(N, 2.0 * 14.007 / totalMassPNDI);
  PNDI->AddElement(O, 4.0 * 15.999 / totalMassPNDI);
  PNDI->AddElement(S, 1.9 * 32.065 / totalMassPNDI);

  // PNDI-T10 loaded with 10 wt% enriched B4C
  // 10 wt% means: for 100g of PNDI-T10, add 10g of B4C → total 110g
  // mass fractions: PNDI = 100/110, B4C = 10/110
  G4double fPNDI = 100. / 110.;
  G4double fB4C = 10. / 110.;
  G4double densityBlend = 1. / (fPNDI / 1.5 + fB4C / 2.52) * g / cm3;
  G4Material* PNDI_B4C =
    new G4Material("PNDI-T10-B4C", densityBlend, ncomponents = 2);
  PNDI_B4C->AddMaterial(PNDI, fPNDI);
  PNDI_B4C->AddMaterial(B4C, fB4C);

  // Indium Tin Oxide (ITO): 90% In2O3 + 10% SnO2 by mass
  G4double fIn2O3 = 0.90;
  G4double fSnO2 = 0.10;
  G4double massIn2O3 = 2 * 114.818 + 3 * 15.999;
  G4double massSnO2 = 118.710 + 2 * 15.999;
  G4double totalMassITO = fIn2O3 * massIn2O3 + fSnO2 * massSnO2;
  G4double wIn = fIn2O3 * 2 * 114.818 / totalMassITO;
  G4double wSn = fSnO2 * 118.710 / totalMassITO;
  G4double wO_ITO = (fIn2O3 * 3 * 15.999 + fSnO2 * 2 * 15.999) / totalMassITO;
  G4Material* ITO = new G4Material("ITO", density = 7.14 * g / cm3, ncomponents = 3);
  ITO->AddElement(In, wIn);
  ITO->AddElement(Sn, wSn);
  ITO->AddElement(O, wO_ITO);

  // Quartz glass (SiO2)
  G4Material* QuartzGlass = new G4Material("QuartzGlass", density = 2.20 * g / cm3, ncomponents = 2);
  QuartzGlass->AddElement(Si, natoms = 1);
  QuartzGlass->AddElement(O, natoms = 2);

  // Silver-loaded conductive epoxy (bisphenol-F diglycidyl ether-like resin
  // filled with Ag flakes). Modeled as a single homogeneous mixture: the flake
  // scale (~um) is far below the 662 keV gamma mean free path and the secondary
  // electron range, so a mass-fraction average is accurate. Composition and
  // density from the supplied datasheet.
  G4Element* Ag = manager->FindOrBuildElement(47);
  G4Material* silverEpoxy = new G4Material("SilverEpoxy", density = 2.34 * g / cm3, ncomponents = 4);
  silverEpoxy->AddElement(Ag, 0.5700);
  silverEpoxy->AddElement(C, 0.3140);
  silverEpoxy->AddElement(H, 0.0278);
  silverEpoxy->AddElement(O, 0.0881);

  // Polyethylene (C2H4)n — repeat unit C2H4, density for validation slabs.
  // Hydrogen is wired as TS_H_of_Polyethylene so the bound-atom
  // S(alpha,beta) thermal-scattering kernel (n < ~4 eV) is used in this
  // organic material, consistent with PNDI-T10 / PNDI-T10-B4C.
  G4Material* polyethylene = new G4Material("polyethylene", density = 0.93 * g / cm3, ncomponents = 2);
  polyethylene->AddElement(C, natoms = 2);
  polyethylene->AddElement(H_TS_Poly, natoms = 4);

  // ABS (acrylonitrile-butadiene-styrene) — bulk approximation C15H17N at
  // molded-plastic density (~1.04 g/cm3). Used by shine-dd*.mac moderator slab.
  {
    const G4double massC = 15.0 * 12.011 * g / mole;
    const G4double massH = 17.0 * 1.008 * g / mole;
    const G4double massN = 1.0 * 14.007 * g / mole;
    const G4double massTot = massC + massH + massN;
    G4Material* ABS = new G4Material("ABS", density = 1.04 * g / cm3, ncomponents = 3);
    ABS->AddElement(C, massC / massTot);
    ABS->AddElement(H, massH / massTot);
    ABS->AddElement(N, massN / massTot);
  }

  //
  // Thermal-scattering ("TS_*") material variants.
  //
  // To activate the G4ParticleHPThermalScattering model below ~4 eV the
  // moderator elements must carry the canonical TS names that map to the
  // bound-atom S(alpha,beta) data files in G4NDL (see
  // G4ParticleHPThermalScatteringNames.cc):
  //   - "TS_H_of_Water"        : H bound in light water
  //   - "TS_H_of_Polyethylene" : H bound in polyethylene
  //
  // The materials below are parallel to the standard G4_WATER / polyethylene
  // variants so existing macros keep working unchanged. Select them in a
  // macro for maximum thermal-neutron accuracy in moderators:
  //   /testhadr/det/setWorldMaterial     water_TS
  //   /testhadr/det/setModeratorMaterial water_TS
  //   /testhadr/det/setModeratorMaterial polyethylene_TS
  //
  G4Element* H_TS_water =
    new G4Element("TS_H_of_Water", "H", z = 1., a = 1.0079 * g / mole);
  G4Material* water_TS = new G4Material("water_TS", density = 1.000 * g / cm3,
                                        ncomponents = 2, kStateLiquid,
                                        293.15 * kelvin, 1 * atmosphere);
  water_TS->AddElement(H_TS_water, natoms = 2);
  water_TS->AddElement(O, natoms = 1);
  water_TS->GetIonisation()->SetMeanExcitationEnergy(78.0 * eV);
  water_TS->SetChemicalFormula("H_2O");

  G4Material* polyethylene_TS = new G4Material("polyethylene_TS",
                                               density = 0.94 * g / cm3,
                                               ncomponents = 2);
  polyethylene_TS->AddElement(C, natoms = 2);
  polyethylene_TS->AddElement(H_TS_Poly, natoms = 4);

  // Heavy water (D2O) with TS deuterium for thermal neutron scattering below ~4 eV.
  G4Isotope* H2 = new G4Isotope("H2", 1, 2);
  G4Element* D_TS = new G4Element("TS_D_of_Heavy_Water", "D", 1);
  D_TS->AddIsotope(H2, 100 * perCent);
  G4Material* HeavyWater = new G4Material("HeavyWater", 1.11 * g / cm3, ncomponents = 2,
                                          kStateLiquid, 293.15 * kelvin, 1 * atmosphere);
  HeavyWater->AddElement(D_TS, natoms = 2);
  HeavyWater->AddElement(O, natoms = 1);

  //
  // examples of gas in non STP conditions
  //
  G4double temperature, pressure;

  G4Material* CO2 =
    new G4Material("CarbonicGas", density = 27. * mg / cm3, ncomponents = 2, kStateGas,
                   temperature = 325. * kelvin, pressure = 50. * atmosphere);
  CO2->AddElement(C, natoms = 1);
  CO2->AddElement(O, natoms = 2);

  new G4Material("ArgonGas", z = 18, a = 39.948 * g / mole, density = 1.782 * mg / cm3, kStateGas,
                 273.15 * kelvin, 1 * atmosphere);
  //
  // example of vacuum
  //
  density = universe_mean_density;  // from PhysicalConstants.h
  pressure = 3.e-18 * pascal;
  temperature = 2.73 * kelvin;
  new G4Material("Galactic", z = 1., a = 1.008 * g / mole, density, kStateGas, temperature,
                 pressure);

  //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::MaterialWithSingleIsotope(G4String name, G4String symbol,
                                                            G4double density, G4int Z, G4int A)
{
  // define a material from an isotope
  //
  G4int ncomponents;
  G4double abundance, massfraction;

  G4Isotope* isotope = new G4Isotope(symbol, Z, A);

  G4Element* element = new G4Element(name, symbol, ncomponents = 1);
  element->AddIsotope(isotope, abundance = 100. * perCent);

  G4Material* material = new G4Material(name, density, ncomponents = 1);
  material->AddElement(element, massfraction = 100. * perCent);

  return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
  fLayerThickness = 0.;
  for (G4int iAbs = 1; iAbs <= fNbOfAbsor; iAbs++) {
    fLayerThickness += fAbsorThickness[iAbs];
  }
  fCalorThickness = fNbOfLayers * fLayerThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool DetectorConstruction::IsModeratorLogical(const G4LogicalVolume* lv) const
{
  if (lv == nullptr) return false;
  for (G4int i = 1; i <= kMaxModeratorLayers; ++i) {
    if (lv == fLogicModeratorLayer[i]) return true;
  }
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector DetectorConstruction::GetDefaultSourcePosition() const
{
  // Bud box centered at origin; default source sits on the +Z side at
  // sourceStandoff from the outermost envelope on +Z (PC shell if enabled).
  const G4double zOuterPlusZ = 0.5 * fBoxExternalZ
                             + (fEncasementEnable ? fEncasementWallThickness : 0.);
  const G4double zSrc = zOuterPlusZ + fSourceStandoff;
  return G4ThreeVector(0., 0., zSrc);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  if (fPhysiWorld) {
    return fPhysiWorld;
  }
  // complete the Calor parameters definition
  ComputeCalorParameters();

  //
  // Lab world (large air volume containing the bud box and any external
  // sources)
  //
  fSolidWorld = new G4Box("World", fLabSizeX / 2, fLabSizeY / 2, fLabSizeZ / 2);
  fLogicWorld = new G4LogicalVolume(fSolidWorld, fWorldMaterial, "World");
  fPhysiWorld = new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World",
                                  0, false, 0);

  //
  // Bud box (aluminum wall) and inner air cavity. If wall thickness is 0,
  // the box and cavity are skipped entirely; the calorimeter is then
  // placed directly inside the lab world (legacy bare-detector geometry).
  //
  G4LogicalVolume* calorMother = fLogicWorld;
  G4double frontZ = 0.;

  G4bool encActive = fEncasementEnable && (fBoxWallThickness > 0.);
  if (fEncasementEnable && fBoxWallThickness <= 0.) {
    G4cout << "\n ---> Warning from DetectorConstruction: encasement requested but "
              "bud box is disabled (wall thickness = 0). Encasement not built."
           << G4endl;
    encActive = false;
  }

  G4double budOffsetZ = 0.;
  if (encActive) {
    const G4double t = fEncasementWallThickness;
    const G4double ix = fEncasementOuterX - 2. * t;
    const G4double iy = fEncasementOuterY - 2. * t;
    const G4double iz = fEncasementOuterZ - 2. * t;
    if (ix <= fBoxExternalX || iy <= fBoxExternalY || iz <= fBoxExternalZ || ix <= 0.
        || iy <= 0. || iz <= 0.)
    {
      G4cout << "\n ---> Warning from DetectorConstruction: encasement inner cavity "
                "must exceed bud-box externals in X, Y, Z (and wall thickness). "
                "Encasement not built."
             << G4endl;
      encActive = false;
    }
    else {
      // Inner +Z face of the air cavity flush with Al bud box +Z exterior (z = +halfZbox).
      const G4double zInnerPlusZ = 0.5 * fBoxExternalZ;
      const G4double innerHz = 0.5 * iz;
      const G4double zAirCenter = zInnerPlusZ - innerHz;
      budOffsetZ = innerHz - 0.5 * fBoxExternalZ;

      G4Material* airEnc = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
      fSolidEncasementAir = new G4Box("EncasementAirCavity", ix / 2, iy / 2, innerHz);
      fLogicEncasementAir =
        new G4LogicalVolume(fSolidEncasementAir, airEnc, "EncasementAirCavity");
      fPhysiEncasementAir = new G4PVPlacement(0,
                                              G4ThreeVector(0., 0., zAirCenter),
                                              fLogicEncasementAir, "EncasementAirCavity",
                                              fLogicWorld, false, 0);

      G4Material* pcMat = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYCARBONATE");
      auto* solidPCOuter = new G4Box("PC_enc_outer",
                                     fEncasementOuterX / 2, fEncasementOuterY / 2,
                                     fEncasementOuterZ / 2);
      auto* solidPCInner = new G4Box("PC_enc_inner", ix / 2, iy / 2, innerHz);
      auto* solidPCShell = new G4SubtractionSolid(
        "PC_enc_shell", solidPCOuter, solidPCInner, nullptr, G4ThreeVector(0., 0., 0.));
      fLogicPCShell = new G4LogicalVolume(solidPCShell, pcMat, "PC_enc_shell");
      fPhysiPCShell = new G4PVPlacement(0,
                                        G4ThreeVector(0., 0., zAirCenter),
                                        fLogicPCShell, "PC_enc_shell",
                                        fLogicWorld, false, 0);
    }
  }

  if (fBoxWallThickness > 0.) {
    fSolidBox = new G4Box("BudBox",
                          fBoxExternalX / 2, fBoxExternalY / 2, fBoxExternalZ / 2);
    fLogicBox = new G4LogicalVolume(fSolidBox, fBoxMaterial, "BudBox");
    G4LogicalVolume* budMother = encActive ? fLogicEncasementAir : fLogicWorld;
    fPhysiBox = new G4PVPlacement(0,
                                  G4ThreeVector(0., 0., budOffsetZ),
                                  fLogicBox, "BudBox",
                                  budMother, false, 0);

    G4double cavX = fBoxExternalX - 2. * fBoxWallThickness;
    G4double cavY = fBoxExternalY - 2. * fBoxWallThickness;
    G4double cavZ = fBoxExternalZ - 2. * fBoxWallThickness;
    fSolidCavity = new G4Box("BudBoxCavity", cavX / 2, cavY / 2, cavZ / 2);
    fLogicCavity = new G4LogicalVolume(fSolidCavity, fCavityMaterial, "BudBoxCavity");
    fPhysiCavity = new G4PVPlacement(0, G4ThreeVector(), fLogicCavity, "BudBoxCavity",
                                     fLogicBox, false, 0);

    calorMother = fLogicCavity;
    frontZ = 0.5 * cavZ - fDetectorStandoff;
  }
  else {
    fSolidBox = nullptr;
    fLogicBox = nullptr;
    fPhysiBox = nullptr;
    fSolidCavity = nullptr;
    fLogicCavity = nullptr;
    fPhysiCavity = nullptr;
    // Legacy bare-detector geometry: place the detector front face at +Z
    // of origin. Standoff is interpreted as "from origin" so the source
    // (at z = fBoxExternalZ/2 + fSourceStandoff by default) still reaches
    // a sensible distance.
    frontZ = 0.5 * fCalorThickness;
  }

  //
  // Moderator stack (+Z of bud box or, if encased, +Z of polycarbonate shell):
  // layer 1 adjacent to the box +Z face; each subsequent layer stacks toward +Z
  // (source). Primary travels toward -Z; GPS sits at z = boxHalfZ + (PC wall
  // if encased) + sourceStandoff.
  //
  const G4bool wantModerators =
    (fNbModeratorLayers > 0) || (fModeratorEnable && fModeratorMaterial != nullptr);

  if (wantModerators && fBoxWallThickness <= 0.) {
    G4cout << "\n ---> Warning from DetectorConstruction: moderator requested but "
              "bud box is disabled (wall thickness = 0). Moderators not placed."
           << G4endl;
  }
  else if (wantModerators) {
    G4int nLayers = fNbModeratorLayers;
    G4double layerThick[kMaxModeratorLayers + 1];
    G4Material* layerMat[kMaxModeratorLayers + 1];

    if (nLayers > 0) {
      for (G4int i = 1; i <= nLayers; ++i) {
        layerThick[i] = fModeratorLayerThick[i];
        layerMat[i] = fModeratorLayerMat[i];
      }
    }
    else {
      nLayers = 1;
      layerThick[1] = std::max(0., fSourceStandoff - fModeratorUpstreamGap);
      layerMat[1] = fModeratorMaterial;
    }

    G4double stackSum = 0.;
    for (G4int i = 1; i <= nLayers; ++i) {
      stackSum += layerThick[i];
    }
    const G4double minStandoff = stackSum + fModeratorUpstreamGap;
    if (stackSum <= 0.) {
      G4cout << "\n ---> Warning from DetectorConstruction: moderator stack total "
                "thickness is <= 0. Moderators not placed."
             << G4endl;
    }
    else if (fSourceStandoff < minStandoff) {
      G4cout << "\n ---> Warning from DetectorConstruction: sourceStandoff ("
             << G4BestUnit(fSourceStandoff, "Length")
             << ") is less than moderator stack + gap ("
             << G4BestUnit(minStandoff, "Length")
             << "). Moderators not placed."
             << G4endl;
    }
    else {
      const G4double hx =
        (fModeratorFullX > 0. ? fModeratorFullX : fBoxExternalX) / 2.;
      const G4double hy =
        (fModeratorFullY > 0. ? fModeratorFullY : fBoxExternalY) / 2.;
      const G4double zModPlane = fBoxExternalZ / 2.
                               + (encActive ? fEncasementWallThickness : 0.);

      G4double zStackBase = zModPlane;
      for (G4int i = 1; i <= nLayers; ++i) {
        if (layerMat[i] == nullptr || layerThick[i] <= 0.) {
          G4cout << "\n ---> Warning from DetectorConstruction: moderator layer " << i
                 << " has invalid material or thickness. Layer skipped." << G4endl;
          continue;
        }
        const G4double hz = layerThick[i] / 2.;
        const G4double zCenter = zStackBase + hz;

        G4String solidName = "ModeratorSolid_";
        solidName += std::to_string(i);
        G4String logicName = "Moderator_";
        logicName += std::to_string(i);
        G4String physName = logicName;

        fSolidModeratorLayer[i] = new G4Box(solidName, hx, hy, hz);
        fLogicModeratorLayer[i] =
          new G4LogicalVolume(fSolidModeratorLayer[i], layerMat[i], logicName);
        fPhysiModeratorLayer[i] = new G4PVPlacement(
          0, G4ThreeVector(0., 0., zCenter), fLogicModeratorLayer[i], physName,
          fLogicWorld, false, i);

        zStackBase += layerThick[i];
      }
    }
  }

  //
  // Calorimeter: rotated -90 deg about Y so the (locally X-stacked) absorber
  // layout becomes Z-stacked in the global frame, with the first absorber
  // (local -X end) facing the +Z (source) side.
  //
  // With rotateY(-90), the local -X end (where absorber 1 sits) maps to global
  // +Z = front face, i.e. the side closest to the default source position, so
  // the stack order from the source inward is absorber 1, 2, 3, ...  (e.g.
  // silver -> polymer -> ITO -> glass). rotateY(+90) reverses this and puts the
  // last absorber at the front, which is NOT what the geometry intends.
  if (fCalorRotation == nullptr) {
    fCalorRotation = new G4RotationMatrix();
    fCalorRotation->rotateY(-90. * deg);
  }
  G4double calorCenterZ = frontZ - 0.5 * fCalorThickness;

  fSolidCalor = new G4Box("Calorimeter",
                          fCalorThickness / 2, fCalorSizeY / 2, fCalorSizeZ / 2);
  fLogicCalor = new G4LogicalVolume(fSolidCalor, fCavityMaterial, "Calorimeter");
  fPhysiCalor = new G4PVPlacement(fCalorRotation,
                                  G4ThreeVector(0., 0., calorCenterZ),
                                  fLogicCalor, "Calorimeter",
                                  calorMother, false, 0);

  //
  // Layers (filled with absorbers; lateral and longitudinal extent fully
  // tile the calorimeter, so cavity material here is only a placeholder
  // for any unfilled internal volume)
  //
  fSolidLayer = new G4Box("Layer", fLayerThickness / 2, fCalorSizeY / 2, fCalorSizeZ / 2);
  fLogicLayer = new G4LogicalVolume(fSolidLayer, fCavityMaterial, "Layer");
  if (fNbOfLayers > 1) {
    fPhysiLayer =
      new G4PVReplica("Layer", fLogicLayer, fLogicCalor, kXAxis, fNbOfLayers, fLayerThickness);
  }
  else {
    fPhysiLayer =
      new G4PVPlacement(0, G4ThreeVector(), fLogicLayer, "Layer", fLogicCalor, false, 0);
  }
  //
  // Absorbers
  //
  G4double xfront = -0.5 * fLayerThickness;
  fLogicEpoxyContainer = nullptr;
  for (G4int k = 1; k <= fNbOfAbsor; ++k) {
    // Full-area box slab occupying this absorber's slot in the stack. Always
    // built so the stacking geometry (thicknesses, downstream layer positions)
    // is unchanged whether or not the blob option is active.
    fSolidAbsor[k] = new G4Box("Absorber",
                               fAbsorThickness[k] / 2, fCalorSizeY / 2, fCalorSizeZ / 2);

    G4double xcenter = xfront + 0.5 * fAbsorThickness[k];
    xfront += fAbsorThickness[k];

    const G4bool blobHere = (fSilverEpoxyBlob && k == 1);
    if (blobHere) {
      // Localized silver-epoxy contact: the absorber-1 slot is filled with
      // cavity material (air) and a half-ellipsoid of absorber-1 material
      // (SilverEpoxy) is placed inside it, resting on the front face of the
      // next layer (the pure Ag film). fLogicAbsor[1] is the ellipsoid so all
      // scoring that keys on GetAbsorberLogical(1) targets the epoxy blob; the
      // surrounding air is inert filler that maps to no absorber.
      fLogicEpoxyContainer =
        new G4LogicalVolume(fSolidAbsor[k], fCavityMaterial, "EpoxyContainer");
      new G4PVPlacement(0, G4ThreeVector(xcenter, 0., 0.), fLogicEpoxyContainer,
                        "EpoxyContainer", fLogicLayer, false, k);

      // Half-ellipsoid: circular base of radius fEpoxyBlobRadius (in the layer
      // Y-Z plane) and semi-axis fEpoxyBlobHeight along the stack (layer X).
      // Ellipsoid local axis is Z; keep the top half (z in [0, height]) so the
      // flat base sits at z = 0, then rotate so local +Z -> layer -X (dome
      // pointing toward the +Z source side).
      G4double blobHeight = std::min(fEpoxyBlobHeight, fAbsorThickness[k]);
      G4Ellipsoid* solidBlob =
        new G4Ellipsoid("SilverEpoxyBlob", fEpoxyBlobRadius, fEpoxyBlobRadius,
                        blobHeight, 0., blobHeight);
      fLogicAbsor[k] = new G4LogicalVolume(solidBlob, fAbsorMaterial[k],
                                           fAbsorMaterial[k]->GetName());

      // G4PVPlacement takes a passive (frame) rotation, so rotateY(+90 deg)
      // maps the ellipsoid's local +Z (dome axis) onto the layer -X direction,
      // i.e. the dome points toward the +Z source side while the flat base
      // rests against the +X face (the Ag-film interface).
      delete fEpoxyBlobRotation;
      fEpoxyBlobRotation = new G4RotationMatrix();
      fEpoxyBlobRotation->rotateY(90. * deg);
      // Base sits at the container +X face (interface with the Ag film);
      // container-local +X face is at +fAbsorThickness[1]/2.
      fPhysiAbsor[k] =
        new G4PVPlacement(fEpoxyBlobRotation,
                          G4ThreeVector(0.5 * fAbsorThickness[k], 0., 0.),
                          fLogicAbsor[k], fAbsorMaterial[k]->GetName(),
                          fLogicEpoxyContainer, false, k, true);  // checkOverlaps
    }
    else {
      fLogicAbsor[k] = new G4LogicalVolume(fSolidAbsor[k], fAbsorMaterial[k],
                                           fAbsorMaterial[k]->GetName());
      fPhysiAbsor[k] = new G4PVPlacement(0, G4ThreeVector(xcenter, 0., 0.), fLogicAbsor[k],
                                         fAbsorMaterial[k]->GetName(), fLogicLayer, false,
                                         k);  // copy number
    }
  }

  //
  // Front-face monitor: thin slab placed immediately in front (+Z) of the
  // calorimeter front face, lateral extent matching the calorimeter, used
  // to score the neutron spectrum entering the detector after any
  // environment interactions. A 1 nm gap avoids coincident-face issues.
  //
  const G4double monitorThickness = 0.5 * micrometer;
  const G4double monitorGap = 1. * nm;
  G4double monitorCenterZ = frontZ + monitorGap + 0.5 * monitorThickness;

  // Lateral extent in global X-Y matches the rotated calorimeter front face
  // (local Z -> global X, local Y -> global Y after rotateY(+90 deg)).
  fSolidMonitor = new G4Box("MonitorFront",
                            fCalorSizeZ / 2, fCalorSizeY / 2, monitorThickness / 2);
  fLogicMonitor = new G4LogicalVolume(fSolidMonitor, fCavityMaterial, "MonitorFront");
  fPhysiMonitor = new G4PVPlacement(0,
                                    G4ThreeVector(0., 0., monitorCenterZ),
                                    fLogicMonitor, "MonitorFront",
                                    calorMother, false, 0);

  //
  // Cs-137 check source envelope: regular hexagonal prism (local axis Z, then
  // rotated +90 deg about Y so the axis lies along global +X). phiStart = 30 deg
  // (half of 360/6) orients a rectangular side face (10 mm x 7 mm) with its
  // outward normal along local +X; after rotateY(+90 deg) that face lies in the
  // detector plane (XY) with its bottom at the bud-box +Z exterior
  // (z = zOuterPlusZ). GPS: /gps/pos/confine Cs137Source.
  //
  if (fCs137SourceEnable) {
    constexpr G4double kHexSide = 7. * mm;
    constexpr G4double kHalfLen = 5. * mm;       // 1 cm total length along +/- X
    constexpr G4double kBoreRadius = 4. * mm;    // 0.4 cm hollow core
    const G4double apothem = kHexSide * std::cos(30. * deg);

    const G4double zOuterPlusZ = 0.5 * fBoxExternalZ
                               + (fEncasementEnable ? fEncasementWallThickness : 0.);
    // Side face at local +X (distance apothem) maps to global z = centreZ - apothem
    // after rotateY(+90 deg). Set that face on the box +Z exterior.
    const G4double centreZ = zOuterPlusZ + apothem;

    G4double zPlane[] = {-kHalfLen, kHalfLen};
    G4double rInner[] = {0., 0.};
    G4double rOuter[] = {apothem, apothem};

    // phiStart = 30 deg: side face normal along local +X (phiStart = 0 puts a
    // vertex on +X instead of the 10 mm x 7 mm face).
    auto* hex = new G4Polyhedra("Cs137Hex", 30. * deg, twopi, 6, 2, zPlane, rInner,
                                rOuter);
    auto* bore = new G4Tubs("Cs137Bore", 0., kBoreRadius, kHalfLen + 1. * nm, 0., twopi);
    fSolidCs137Source = new G4SubtractionSolid("Cs137SourceSolid", hex, bore);
    fLogicCs137Source =
      new G4LogicalVolume(fSolidCs137Source, fWorldMaterial, "Cs137Source");

    if (fCs137SourceRotation == nullptr) {
      fCs137SourceRotation = new G4RotationMatrix();
      fCs137SourceRotation->rotateY(90. * deg);
    }
    fPhysiCs137Source = new G4PVPlacement(fCs137SourceRotation,
                                          G4ThreeVector(0., 0., centreZ),
                                          fLogicCs137Source, "Cs137Source",
                                          fLogicWorld, false, 0);
  }

  DefineRegionsAndCuts();

  PrintCalorParameters();

  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineRegionsAndCuts()
{
  // Range-based production cuts (gamma, e-, e+, proton) control secondary
  // production thresholds; they are not tracking cuts.
  //
  // Active detector stack (all absorbers + front monitor, which sits at the
  // detector front interface): the tightest settings, since these are the
  // sensitive volumes / interfaces of interest.
  //   - all particles: 0.1 um production cut
  //   - metal contacts: gamma 0.01 um (to resolve high-Z photoelectric lines)
  //   - all absorbers: 0.1 um max-step user limit (enforced via
  //     G4StepLimiterPhysics in PhysicsList)
  constexpr G4double kAbsorberRange = 0.1 * micrometer;
  constexpr G4double kMetalGammaRange = 0.01 * micrometer;
  constexpr G4double kAbsorberMaxStep = 0.1 * micrometer;

  // Surrounding passive volumes keep coarser cuts (they only shape the
  // incident field; their internal secondaries are not scored).
  constexpr G4double kBoxWallRange = 5.0 * micrometer;
  constexpr G4double kModeratorRange = 5.0 * mm;

  auto attachRegion = [](const G4String& name, G4LogicalVolume* lv, G4ProductionCuts* cuts) {
    auto* reg = new G4Region(name);
    reg->AddRootLogicalVolume(lv);
    reg->SetProductionCuts(cuts);
  };

  // Thin front-face monitor slab sits at the detector front interface: give it
  // absorber-level cuts and max-step so the entering field is resolved at the
  // same precision as the stack.
  if (fLogicMonitor != nullptr) {
    auto* cuts = new G4ProductionCuts();
    cuts->SetProductionCut(kAbsorberRange, "gamma");
    cuts->SetProductionCut(kAbsorberRange, "e-");
    cuts->SetProductionCut(kAbsorberRange, "e+");
    cuts->SetProductionCut(kAbsorberRange, "proton");
    attachRegion("Reg_MonitorFront", fLogicMonitor, cuts);
    fLogicMonitor->SetUserLimits(new G4UserLimits(kAbsorberMaxStep));
  }

  if (fLogicModeratorLayer[1] != nullptr) {
    for (G4int i = 1; i <= kMaxModeratorLayers; ++i) {
      G4LogicalVolume* lv = fLogicModeratorLayer[i];
      if (lv == nullptr) continue;
      auto* cuts = new G4ProductionCuts();
      cuts->SetProductionCut(kModeratorRange);
      G4String regName = "Reg_Moderator_";
      regName += std::to_string(i);
      attachRegion(regName, lv, cuts);
    }
  }

  if (fLogicBox != nullptr) {
    auto* cuts = new G4ProductionCuts();
    cuts->SetProductionCut(kBoxWallRange, "gamma");
    cuts->SetProductionCut(kBoxWallRange, "e-");
    cuts->SetProductionCut(kBoxWallRange, "e+");
    cuts->SetProductionCut(0.5 * mm, "proton");
    attachRegion("Reg_BudBoxWall", fLogicBox, cuts);
  }

  if (fLogicPCShell != nullptr) {
    constexpr G4double kPCShellRange = 5.0 * micrometer;
    auto* cuts = new G4ProductionCuts();
    cuts->SetProductionCut(kPCShellRange, "gamma");
    cuts->SetProductionCut(kPCShellRange, "e-");
    cuts->SetProductionCut(kPCShellRange, "e+");
    cuts->SetProductionCut(0.2 * mm, "proton");
    attachRegion("Reg_PC_Encasement", fLogicPCShell, cuts);
  }

  // Every absorber: 0.1 um cuts (all particles) + 0.1 um max-step. Metal
  // contacts additionally get a 0.01 um gamma cut.
  for (G4int k = 1; k <= fNbOfAbsor; ++k) {
    G4LogicalVolume* lv = fLogicAbsor[k];
    if (lv == nullptr) continue;
    G4Material* mat = fAbsorMaterial[k];
    if (mat == nullptr) continue;

    const G4bool isMetal = IsMetalAbsorberMaterial(mat);
    G4String regName = "Reg_Absor_";
    regName += std::to_string(k).c_str();
    regName += (isMetal ? "_metal" : "_absorber");

    auto* cuts = new G4ProductionCuts();
    cuts->SetProductionCut(isMetal ? kMetalGammaRange : kAbsorberRange, "gamma");
    cuts->SetProductionCut(kAbsorberRange, "e-");
    cuts->SetProductionCut(kAbsorberRange, "e+");
    cuts->SetProductionCut(kAbsorberRange, "proton");
    attachRegion(regName, lv, cuts);

    // Upper bound on step length in the absorber; transport may still use
    // shorter steps when msc/ioni models require it (never coarsened).
    lv->SetUserLimits(new G4UserLimits(kAbsorberMaxStep));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintCalorParameters()
{
  G4int prec = 4, wid = prec + 2;
  G4int dfprec = G4cout.precision(prec);

  G4double totLength(0.), totRadl(0.), totNuclear(0.);

  G4cout << "\n-------------------------------------------------------------"
         << "\n ---> The calorimeter is " << fNbOfLayers << " layers of:";
  for (G4int i = 1; i <= fNbOfAbsor; ++i) {
    G4Material* material = fAbsorMaterial[i];
    G4double radl = material->GetRadlen();
    G4double nuclearl = material->GetNuclearInterLength();
    G4double sumThickness = fNbOfLayers * fAbsorThickness[i];
    G4double nbRadl = sumThickness / radl;
    G4double nbNuclearl = sumThickness / nuclearl;
    totLength += sumThickness;
    totRadl += nbRadl;
    totNuclear += nbNuclearl;
    G4cout << "\n   " << std::setw(12) << fAbsorMaterial[i]->GetName() << ": " << std::setw(wid)
           << G4BestUnit(fAbsorThickness[i], "Length") << "  --->  sum = " << std::setw(wid)
           << G4BestUnit(sumThickness, "Length") << " = " << std::setw(wid) << nbRadl << " Radl "
           << " = " << std::setw(wid) << nbNuclearl << " NuclearInteractionLength ";
  }
  G4cout << "\n\n                       total thickness = " << std::setw(wid)
         << G4BestUnit(totLength, "Length") << " = " << std::setw(wid) << totRadl << " Radl "
         << " = " << std::setw(wid) << totNuclear << " NuclearInteractionLength " << G4endl;

  G4cout << "                     transverse sizeY  = " << std::setw(wid)
         << G4BestUnit(fCalorSizeY, "Length") << G4endl;
  G4cout << "                     transverse sizeZ  = " << std::setw(wid)
         << G4BestUnit(fCalorSizeZ, "Length") << G4endl;
  if (fSilverEpoxyBlob) {
    G4cout << " Absorber 1 geometry : localized half-ellipsoid blob ("
           << fAbsorMaterial[1]->GetName() << ")   base radius = "
           << G4BestUnit(fEpoxyBlobRadius, "Length")
           << "   height = " << G4BestUnit(fEpoxyBlobHeight, "Length")
           << "   (surrounded by cavity material in the absorber-1 slot)"
           << G4endl;
  }
  else {
    G4cout << " Absorber 1 geometry : full-area box" << G4endl;
  }
  G4cout << "-------------------------------------------------------------\n";

  G4cout << "\n Bud box (external) : "
         << G4BestUnit(fBoxExternalX, "Length") << " x "
         << G4BestUnit(fBoxExternalY, "Length") << " x "
         << G4BestUnit(fBoxExternalZ, "Length")
         << "   wall = " << G4BestUnit(fBoxWallThickness, "Length");
  if (fBoxWallThickness <= 0.) G4cout << "  (box DISABLED)";
  G4cout << G4endl;
  G4cout << " Detector standoff (inside +Z face -> detector front) = "
         << G4BestUnit(fDetectorStandoff, "Length") << G4endl;
  G4cout << " Lab (world) size : "
         << G4BestUnit(fLabSizeX, "Length") << " x "
         << G4BestUnit(fLabSizeY, "Length") << " x "
         << G4BestUnit(fLabSizeZ, "Length") << G4endl;
  G4cout << " Default source standoff (from outer +Z of Al or PC shell) = "
         << G4BestUnit(fSourceStandoff, "Length") << G4endl;
  if (fPhysiPCShell != nullptr) {
    G4cout << " Polycarbonate encasement : enabled   outer "
           << G4BestUnit(fEncasementOuterX, "Length") << " x "
           << G4BestUnit(fEncasementOuterY, "Length") << " x "
           << G4BestUnit(fEncasementOuterZ, "Length")
           << "   wall = " << G4BestUnit(fEncasementWallThickness, "Length")
           << "   material = G4_POLYCARBONATE" << G4endl;
  }
  else {
    G4cout << " Polycarbonate encasement : disabled" << G4endl;
  }
  if (fPhysiCs137Source != nullptr) {
    const G4ThreeVector pos = fPhysiCs137Source->GetTranslation();
    G4cout << " Cs-137 check source : enabled   hex side = 7 mm   length = 1 cm"
           << " (axis || global +X)   bore radius = 4 mm"
           << "   flat face on box +Z at z = "
           << G4BestUnit(0.5 * fBoxExternalZ
                        + (fEncasementEnable ? fEncasementWallThickness : 0.),
                        "Length")
           << "   centre = (0, 0, "
           << G4BestUnit(pos.z(), "Length") << ")" << G4endl;
  }
  else {
    G4cout << " Cs-137 check source : disabled" << G4endl;
  }
  if (fPhysiModeratorLayer[1] != nullptr) {
    G4double stackSum = 0.;
    const G4double modX = (fModeratorFullX > 0. ? fModeratorFullX : fBoxExternalX);
    const G4double modY = (fModeratorFullY > 0. ? fModeratorFullY : fBoxExternalY);
    G4cout << " Moderator stack : enabled   upstream gap (GPS side) = "
           << G4BestUnit(fModeratorUpstreamGap, "Length") << G4endl;
    G4cout << "   transverse (full X x Y) = " << G4BestUnit(modX, "Length") << " x "
           << G4BestUnit(modY, "Length");
    if (fModeratorFullX <= 0. && fModeratorFullY <= 0.) {
      G4cout << "  (bud-box externals)";
    }
    G4cout << G4endl;
    for (G4int i = 1; i <= kMaxModeratorLayers; ++i) {
      if (fPhysiModeratorLayer[i] == nullptr) continue;
      G4double dz = fModeratorLayerThick[i];
      if (fNbModeratorLayers == 0 && i == 1) {
        dz = std::max(0., fSourceStandoff - fModeratorUpstreamGap);
      }
      stackSum += dz;
      G4cout << "   layer " << i << " : dz = " << G4BestUnit(dz, "Length")
             << "   material = "
             << fLogicModeratorLayer[i]->GetMaterial()->GetName() << G4endl;
    }
    G4cout << "   total stack height = " << G4BestUnit(stackSum, "Length")
           << G4endl;
  }
  else if (GetModeratorEnable()) {
    G4cout << " Moderator stack : requested but not placed (see warnings)."
           << G4endl;
  }
  else {
    G4cout << " Moderator stack : disabled" << G4endl;
  }
  G4cout << "-------------------------------------------------------------\n";

  G4cout << "\n" << fWorldMaterial << G4endl;
  G4cout << "\n" << fCavityMaterial << G4endl;
  if (fBoxWallThickness > 0.) G4cout << "\n" << fBoxMaterial << G4endl;
  for (G4int j = 1; j <= fNbOfAbsor; ++j) {
    G4cout << "\n" << fAbsorMaterial[j] << G4endl;
  }
  G4cout << "\n-------------------------------------------------------------\n";

  // restore default format
  G4cout.precision(dfprec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& material)
{
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) {
    fWorldMaterial = pttoMaterial;
    if (fLogicWorld) {
      fLogicWorld->SetMaterial(fWorldMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCavityMaterial(const G4String& material)
{
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) {
    fCavityMaterial = pttoMaterial;
    if (fLogicCavity) fLogicCavity->SetMaterial(fCavityMaterial);
    if (fLogicCalor) fLogicCalor->SetMaterial(fCavityMaterial);
    if (fLogicLayer) fLogicLayer->SetMaterial(fCavityMaterial);
    if (fLogicMonitor) fLogicMonitor->SetMaterial(fCavityMaterial);
    if (fLogicWorld) {
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetBoxMaterial(const G4String& material)
{
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) {
    fBoxMaterial = pttoMaterial;
    if (fLogicBox) {
      fLogicBox->SetMaterial(fBoxMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetBoxExternal(G4double xExt, G4double yExt, G4double zExt)
{
  if (xExt <= 0. || yExt <= 0. || zExt <= 0.) {
    G4cout << "\n ---> warning from SetBoxExternal: dimensions ("
           << xExt << "," << yExt << "," << zExt
           << ") must all be positive. Command refused" << G4endl;
    return;
  }
  fBoxExternalX = xExt;
  fBoxExternalY = yExt;
  fBoxExternalZ = zExt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetBoxWallThickness(G4double t)
{
  if (t < 0.) {
    G4cout << "\n ---> warning from SetBoxWallThickness: thickness " << t
           << " must be >= 0. Command refused" << G4endl;
    return;
  }
  fBoxWallThickness = t;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorStandoff(G4double d)
{
  if (d < 0.) {
    G4cout << "\n ---> warning from SetDetectorStandoff: standoff " << d
           << " must be >= 0. Command refused" << G4endl;
    return;
  }
  fDetectorStandoff = d;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetLabSize(G4double xExt, G4double yExt, G4double zExt)
{
  if (xExt <= 0. || yExt <= 0. || zExt <= 0.) {
    G4cout << "\n ---> warning from SetLabSize: dimensions ("
           << xExt << "," << yExt << "," << zExt
           << ") must all be positive. Command refused" << G4endl;
    return;
  }
  fLabSizeX = xExt;
  fLabSizeY = yExt;
  fLabSizeZ = zExt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSourceStandoff(G4double d)
{
  if (d < 0.) {
    G4cout << "\n ---> warning from SetSourceStandoff: standoff " << d
           << " must be >= 0. Command refused" << G4endl;
    return;
  }
  fSourceStandoff = d;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCs137SourceEnable(G4bool on) { fCs137SourceEnable = on; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSilverEpoxyBlob(G4bool on) { fSilverEpoxyBlob = on; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSilverEpoxyBlobSize(G4double radius, G4double height)
{
  if (radius > 0.) fEpoxyBlobRadius = radius;
  if (height > 0.) fEpoxyBlobHeight = height;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetModeratorEnable(G4bool on)
{
  fModeratorEnable = on;
}

void DetectorConstruction::SetNbOfModeratorLayers(G4int ival)
{
  if (ival < 0) {
    G4cout << "\n ---> warning from SetNbOfModeratorLayers: " << ival
           << " must be >= 0. Command refused" << G4endl;
    return;
  }
  if (ival > kMaxModeratorLayers) {
    G4cout << "\n ---> warning from SetNbOfModeratorLayers: " << ival
           << " exceeds kMaxModeratorLayers (" << kMaxModeratorLayers
           << "). Command refused" << G4endl;
    return;
  }
  fNbModeratorLayers = ival;
}

void DetectorConstruction::SetModeratorLayer(G4int i, const G4String& material,
                                             G4double thickness)
{
  if (i < 1 || i > kMaxModeratorLayers) {
    G4cout << "\n ---> warning from SetModeratorLayer: layer index " << i
           << " out of range [1, " << kMaxModeratorLayers << "]. Command refused"
           << G4endl;
    return;
  }
  if (thickness <= 0.) {
    G4cout << "\n ---> warning from SetModeratorLayer: thickness must be > 0. "
              "Command refused"
           << G4endl;
    return;
  }
  G4Material* m = G4Material::GetMaterial(material, false);
  if (m == nullptr) {
    m = G4NistManager::Instance()->FindOrBuildMaterial(material);
  }
  if (m == nullptr) {
    G4cout << "\n ---> warning from SetModeratorLayer: unknown material \""
           << material << "\". Command refused" << G4endl;
    return;
  }
  fModeratorLayerMat[i] = m;
  fModeratorLayerThick[i] = thickness;
  if (fLogicModeratorLayer[i] != nullptr) {
    fLogicModeratorLayer[i]->SetMaterial(m);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

void DetectorConstruction::SetModeratorUpstreamGap(G4double gap)
{
  if (gap < 0.) {
    G4cout << "\n ---> warning from SetModeratorUpstreamGap: gap " << gap
           << " must be >= 0. Command refused" << G4endl;
    return;
  }
  fModeratorUpstreamGap = gap;
}

void DetectorConstruction::SetModeratorMaterial(const G4String& material)
{
  G4Material* m = G4Material::GetMaterial(material, false);
  if (m == nullptr) {
    m = G4NistManager::Instance()->FindOrBuildMaterial(material);
  }
  if (m != nullptr) {
    fModeratorMaterial = m;
    if (fLogicModeratorLayer[1] != nullptr && fNbModeratorLayers == 0) {
      fLogicModeratorLayer[1]->SetMaterial(fModeratorMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
  else {
    G4cout << "\n ---> warning from SetModeratorMaterial: unknown material \""
           << material << "\". Command refused" << G4endl;
  }
}

void DetectorConstruction::SetModeratorTransverseXY(G4double fullX, G4double fullY)
{
  if (fullX <= 0. || fullY <= 0.) {
    fModeratorFullX = 0.;
    fModeratorFullY = 0.;
    return;
  }
  fModeratorFullX = fullX;
  fModeratorFullY = fullY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetEncasementEnable(G4bool on) { fEncasementEnable = on; }

void DetectorConstruction::SetEncasementOuter(G4double xExt, G4double yExt, G4double zExt)
{
  if (xExt <= 0. || yExt <= 0. || zExt <= 0.) {
    G4cout << "\n ---> warning from SetEncasementOuter: dimensions must be positive."
           << G4endl;
    return;
  }
  fEncasementOuterX = xExt;
  fEncasementOuterY = yExt;
  fEncasementOuterZ = zExt;
}

void DetectorConstruction::SetEncasementWallThickness(G4double t)
{
  if (t <= 0.) {
    G4cout << "\n ---> warning from SetEncasementWallThickness: thickness must be > 0."
           << G4endl;
    return;
  }
  fEncasementWallThickness = t;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfLayers(G4int ival)
{
  // set the number of Layers
  //
  if (ival < 1) {
    G4cout << "\n --->warning from SetfNbOfLayers: " << ival
           << " must be at least 1. Command refused" << G4endl;
    return;
  }
  fNbOfLayers = ival;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfAbsor(G4int ival)
{
  // set the number of Absorbers
  //
  if (ival < 1 || ival > (kMaxAbsor - 1)) {
    G4cout << "\n ---> warning from SetfNbOfAbsor: " << ival << " must be at least 1 and and most "
           << kMaxAbsor - 1 << ". Command refused" << G4endl;
    return;
  }
  fNbOfAbsor = ival;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorMaterial(G4int ival, const G4String& material)
{
  // search the material by its name
  //
  if (ival > fNbOfAbsor || ival <= 0) {
    G4cout << "\n --->warning from SetAbsorMaterial: absor number " << ival
           << " out of range. Command refused" << G4endl;
    return;
  }

  G4Material* pttoMaterial = G4Material::GetMaterial(material, false);
  if (!pttoMaterial) {
    pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(material);
  }
  if (pttoMaterial) {
    fAbsorMaterial[ival] = pttoMaterial;
    if (fLogicAbsor[ival]) {
      fLogicAbsor[ival]->SetMaterial(pttoMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorThickness(G4int ival, G4double val)
{
  // change Absorber thickness
  //
  if (ival > fNbOfAbsor || ival <= 0) {
    G4cout << "\n --->warning from SetAbsorThickness: absor number " << ival
           << " out of range. Command refused" << G4endl;
    return;
  }
  if (val <= DBL_MIN) {
    G4cout << "\n --->warning from SetAbsorThickness: thickness " << val
           << " out of range. Command refused" << G4endl;
    return;
  }
  fAbsorThickness[ival] = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCalorSizeYZ(G4double val)
{
  // Square transverse cross-section (sets both Y and Z).
  //
  if (val <= DBL_MIN) {
    G4cout << "\n --->warning from SetCalorSizeYZ: size " << val
           << " out of range. Command refused" << G4endl;
    return;
  }
  fCalorSizeYZ = val;
  fCalorSizeY = val;
  fCalorSizeZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCalorSizeY(G4double val)
{
  if (val <= DBL_MIN) {
    G4cout << "\n --->warning from SetCalorSizeY: size " << val
           << " out of range. Command refused" << G4endl;
    return;
  }
  fCalorSizeY = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCalorSizeZ(G4double val)
{
  if (val <= DBL_MIN) {
    G4cout << "\n --->warning from SetCalorSizeZ: size " << val
           << " out of range. Command refused" << G4endl;
    return;
  }
  fCalorSizeZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4AutoDelete.hh"
#include "G4GlobalMagFieldMessenger.hh"

void DetectorConstruction::ConstructSDandField()
{
  if (fFieldMessenger.Get() == nullptr) {
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue = G4ThreeVector();
    G4GlobalMagFieldMessenger* msg = new G4GlobalMagFieldMessenger(fieldValue);
    // msg->SetVerboseLevel(1);
    G4AutoDelete::Register(msg);
    fFieldMessenger.Put(msg);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
