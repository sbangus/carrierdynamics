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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det) : fDetector(Det)
{
  fTestemDir = new G4UIdirectory("/testhadr/");
  fTestemDir->SetGuidance("UI commands specific to this example");

  fDetDir = new G4UIdirectory("/testhadr/det/");
  fDetDir->SetGuidance("detector construction commands");

  fSizeYZCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/det/setSizeYZ", this);
  fSizeYZCmd->SetGuidance("Set tranverse size of the calorimeter");
  fSizeYZCmd->SetParameterName("Size", false);
  fSizeYZCmd->SetRange("Size>0.");
  fSizeYZCmd->SetUnitCategory("Length");
  fSizeYZCmd->AvailableForStates(G4State_PreInit);
  fSizeYZCmd->SetToBeBroadcasted(false);

  fNbLayersCmd = new G4UIcmdWithAnInteger("/testhadr/det/setNbOfLayers", this);
  fNbLayersCmd->SetGuidance("Set number of layers.");
  fNbLayersCmd->SetParameterName("NbLayers", false);
  fNbLayersCmd->SetRange("NbLayers>0");
  fNbLayersCmd->AvailableForStates(G4State_PreInit);
  fNbLayersCmd->SetToBeBroadcasted(false);

  fNbAbsorCmd = new G4UIcmdWithAnInteger("/testhadr/det/setNbOfAbsor", this);
  fNbAbsorCmd->SetGuidance("Set number of Absorbers.");
  fNbAbsorCmd->SetParameterName("NbAbsor", false);
  fNbAbsorCmd->SetRange("NbAbsor>0");
  fNbAbsorCmd->AvailableForStates(G4State_PreInit);
  fNbAbsorCmd->SetToBeBroadcasted(false);

  fAbsorCmd = new G4UIcommand("/testhadr/det/setAbsor", this);
  fAbsorCmd->SetGuidance("Set the absor nb, the material, the thickness.");
  fAbsorCmd->SetGuidance("  absor number : from 1 to NbOfAbsor");
  fAbsorCmd->SetGuidance("  material name");
  fAbsorCmd->SetGuidance("  thickness (with unit) : t>0.");
  //
  G4UIparameter* AbsNbPrm = new G4UIparameter("AbsorNb", 'i', false);
  AbsNbPrm->SetGuidance("absor number : from 1 to NbOfAbsor");
  AbsNbPrm->SetParameterRange("AbsorNb>0");
  fAbsorCmd->SetParameter(AbsNbPrm);
  //
  G4UIparameter* MatPrm = new G4UIparameter("material", 's', false);
  MatPrm->SetGuidance("material name");
  fAbsorCmd->SetParameter(MatPrm);
  //
  G4UIparameter* ThickPrm = new G4UIparameter("thickness", 'd', false);
  ThickPrm->SetGuidance("thickness of absorber");
  ThickPrm->SetParameterRange("thickness>0.");
  fAbsorCmd->SetParameter(ThickPrm);
  //
  G4UIparameter* unitPrm = new G4UIparameter("unit", 's', false);
  unitPrm->SetGuidance("unit of thickness");
  G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm"));
  unitPrm->SetParameterCandidates(unitList);
  fAbsorCmd->SetParameter(unitPrm);
  //
  fAbsorCmd->AvailableForStates(G4State_PreInit);
  fAbsorCmd->SetToBeBroadcasted(false);

  fIsotopeCmd = new G4UIcommand("/testhadr/det/setIsotopeMat", this);
  fIsotopeCmd->SetGuidance("Build and select a material with single isotope");
  fIsotopeCmd->SetGuidance("  symbol of isotope, Z, A, density of material");
  //
  G4UIparameter* symbPrm = new G4UIparameter("isotope", 's', false);
  symbPrm->SetGuidance("isotope symbol");
  fIsotopeCmd->SetParameter(symbPrm);
  //
  G4UIparameter* ZPrm = new G4UIparameter("Z", 'i', false);
  ZPrm->SetGuidance("Z");
  ZPrm->SetParameterRange("Z>0");
  fIsotopeCmd->SetParameter(ZPrm);
  //
  G4UIparameter* APrm = new G4UIparameter("A", 'i', false);
  APrm->SetGuidance("A");
  APrm->SetParameterRange("A>0");
  fIsotopeCmd->SetParameter(APrm);
  //
  G4UIparameter* densityPrm = new G4UIparameter("density", 'd', false);
  densityPrm->SetGuidance("density of material");
  densityPrm->SetParameterRange("density>0.");
  fIsotopeCmd->SetParameter(densityPrm);
  //
  G4UIparameter* unitPrm1 = new G4UIparameter("unit", 's', false);
  unitPrm1->SetGuidance("unit of density");
  G4String unitList1 = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("g/cm3"));
  unitPrm1->SetParameterCandidates(unitList1);
  fIsotopeCmd->SetParameter(unitPrm1);
  //
  fIsotopeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // ---- Bud-box / lab environment commands ---------------------------------
  // All geometry-mutating commands are PreInit-only and not broadcast in MT,
  // matching the workspace rule for DetectorMessenger.

  fBoxExternalCmd = new G4UIcmdWith3VectorAndUnit("/testhadr/det/setBoxExternal", this);
  fBoxExternalCmd->SetGuidance("Set bud-box external dimensions (X Y Z). "
                               "Wall thickness is set separately.");
  fBoxExternalCmd->SetParameterName("xExt", "yExt", "zExt", false);
  fBoxExternalCmd->SetUnitCategory("Length");
  fBoxExternalCmd->SetRange("xExt>0. && yExt>0. && zExt>0.");
  fBoxExternalCmd->AvailableForStates(G4State_PreInit);
  fBoxExternalCmd->SetToBeBroadcasted(false);

  fBoxWallCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/det/setBoxWallThickness", this);
  fBoxWallCmd->SetGuidance("Set bud-box wall thickness. "
                           "Setting it to 0 disables the box (legacy bare detector).");
  fBoxWallCmd->SetParameterName("Thickness", false);
  fBoxWallCmd->SetUnitCategory("Length");
  fBoxWallCmd->SetRange("Thickness>=0.");
  fBoxWallCmd->AvailableForStates(G4State_PreInit);
  fBoxWallCmd->SetToBeBroadcasted(false);

  fDetStandoffCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/det/setDetectorStandoff", this);
  fDetStandoffCmd->SetGuidance("Distance from inside of +Z box face to detector front face.");
  fDetStandoffCmd->SetParameterName("Standoff", false);
  fDetStandoffCmd->SetUnitCategory("Length");
  fDetStandoffCmd->SetRange("Standoff>=0.");
  fDetStandoffCmd->AvailableForStates(G4State_PreInit);
  fDetStandoffCmd->SetToBeBroadcasted(false);

  fBoxMatCmd = new G4UIcmdWithAString("/testhadr/det/setBoxMaterial", this);
  fBoxMatCmd->SetGuidance("Set bud-box wall material (any NIST or custom name).");
  fBoxMatCmd->SetParameterName("BoxMaterial", false);
  fBoxMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fBoxMatCmd->SetToBeBroadcasted(false);

  fCavityMatCmd = new G4UIcmdWithAString("/testhadr/det/setCavityMaterial", this);
  fCavityMatCmd->SetGuidance("Set cavity / calorimeter mother / monitor material.");
  fCavityMatCmd->SetParameterName("CavityMaterial", false);
  fCavityMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fCavityMatCmd->SetToBeBroadcasted(false);

  fWorldMatCmd = new G4UIcmdWithAString("/testhadr/det/setWorldMaterial", this);
  fWorldMatCmd->SetGuidance("Set lab/world material (the volume containing the bud box).");
  fWorldMatCmd->SetParameterName("WorldMaterial", false);
  fWorldMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fWorldMatCmd->SetToBeBroadcasted(false);

  fLabSizeCmd = new G4UIcmdWith3VectorAndUnit("/testhadr/det/setLabSize", this);
  fLabSizeCmd->SetGuidance("Set lab world dimensions (X Y Z). Should fit the bud box "
                           "and any external sources with margin.");
  fLabSizeCmd->SetParameterName("xLab", "yLab", "zLab", false);
  fLabSizeCmd->SetUnitCategory("Length");
  fLabSizeCmd->SetRange("xLab>0. && yLab>0. && zLab>0.");
  fLabSizeCmd->AvailableForStates(G4State_PreInit);
  fLabSizeCmd->SetToBeBroadcasted(false);

  fSourceStandoffCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/det/setSourceStandoff", this);
  fSourceStandoffCmd->SetGuidance("Default source standoff from the outer +Z face used for "
                                  "GPS defaults: Al bud box, or polycarbonate shell if "
                                  "encasement is enabled. Used by PrimaryGeneratorAction if "
                                  "the macro does not set a GPS position.");
  fSourceStandoffCmd->SetParameterName("Standoff", false);
  fSourceStandoffCmd->SetUnitCategory("Length");
  fSourceStandoffCmd->SetRange("Standoff>=0.");
  fSourceStandoffCmd->AvailableForStates(G4State_PreInit);
  fSourceStandoffCmd->SetToBeBroadcasted(false);

  // Moderator stack between bud-box +Z exterior and GPS plane (lab/world frame).
  fModeratorEnableCmd =
    new G4UIcmdWithABool("/testhadr/det/setModeratorEnable", this);
  fModeratorEnableCmd->SetGuidance(
    "Legacy single-slab mode: if true and setNbOfModeratorLayers is 0, place one "
    "moderator slab from the bud-box +Z outer face toward +Z with thickness = "
    "sourceStandoff - moderatorGap. Prefer setNbOfModeratorLayers + "
    "setModeratorLayer for explicit multi-layer stacks. Requires bud box "
    "(wall thickness > 0).");
  fModeratorEnableCmd->SetParameterName("Enable", false);
  fModeratorEnableCmd->AvailableForStates(G4State_PreInit);
  fModeratorEnableCmd->SetToBeBroadcasted(false);

  fNbModeratorLayersCmd =
    new G4UIcmdWithAnInteger("/testhadr/det/setNbOfModeratorLayers", this);
  fNbModeratorLayersCmd->SetGuidance(
    "Number of moderator layers in the stack (0 = none). Layer 1 sits adjacent to "
    "the bud-box +Z face; each higher index stacks toward +Z (source). Configure "
    "each layer with setModeratorLayer.");
  fNbModeratorLayersCmd->SetParameterName("NbLayers", false);
  fNbModeratorLayersCmd->SetRange("NbLayers>=0");
  fNbModeratorLayersCmd->AvailableForStates(G4State_PreInit);
  fNbModeratorLayersCmd->SetToBeBroadcasted(false);

  fModeratorLayerCmd = new G4UIcommand("/testhadr/det/setModeratorLayer", this);
  fModeratorLayerCmd->SetGuidance("Set moderator layer index, material, and thickness "
                                  "along global Z.");
  fModeratorLayerCmd->SetGuidance("  layer index : from 1 to NbOfModeratorLayers");
  fModeratorLayerCmd->SetGuidance("  material name");
  fModeratorLayerCmd->SetGuidance("  thickness (with unit) : t>0.");
  G4UIparameter* modLayerIdx = new G4UIparameter("LayerNb", 'i', false);
  modLayerIdx->SetGuidance("layer index : from 1 to NbOfModeratorLayers");
  modLayerIdx->SetParameterRange("LayerNb>0");
  fModeratorLayerCmd->SetParameter(modLayerIdx);
  G4UIparameter* modLayerMat = new G4UIparameter("material", 's', false);
  modLayerMat->SetGuidance("material name");
  fModeratorLayerCmd->SetParameter(modLayerMat);
  G4UIparameter* modLayerThick = new G4UIparameter("thickness", 'd', false);
  modLayerThick->SetGuidance("layer thickness along Z");
  modLayerThick->SetParameterRange("thickness>0.");
  fModeratorLayerCmd->SetParameter(modLayerThick);
  G4UIparameter* modLayerUnit = new G4UIparameter("unit", 's', false);
  modLayerUnit->SetGuidance("unit of thickness");
  modLayerUnit->SetParameterCandidates(G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm")));
  fModeratorLayerCmd->SetParameter(modLayerUnit);
  fModeratorLayerCmd->AvailableForStates(G4State_PreInit);
  fModeratorLayerCmd->SetToBeBroadcasted(false);

  fModeratorGapCmd =
    new G4UIcmdWithADoubleAndUnit("/testhadr/det/setModeratorGap", this);
  fModeratorGapCmd->SetGuidance(
    "Thickness of lab/world material (usually air) between GPS plane "
    "(z = boxHalfZ + sourceStandoff) and the upstream (+Z) face of the "
    "moderator stack. Requires sourceStandoff >= sum(layer thicknesses) + gap.");
  fModeratorGapCmd->SetParameterName("Gap", false);
  fModeratorGapCmd->SetUnitCategory("Length");
  fModeratorGapCmd->SetRange("Gap>=0.");
  fModeratorGapCmd->AvailableForStates(G4State_PreInit);
  fModeratorGapCmd->SetToBeBroadcasted(false);

  fModeratorMatCmd = new G4UIcmdWithAString("/testhadr/det/setModeratorMaterial", this);
  fModeratorMatCmd->SetGuidance(
    "Legacy single-slab moderator material (default G4_WATER). Used only when "
    "setModeratorEnable true and setNbOfModeratorLayers is 0.");
  fModeratorMatCmd->SetParameterName("Material", false);
  fModeratorMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fModeratorMatCmd->SetToBeBroadcasted(false);

  fModeratorXYCmd = new G4UIcommand("/testhadr/det/setModeratorTransverse", this);
  fModeratorXYCmd->SetGuidance(
    "Full moderator lateral sizes X Y. Use two zeros with unit to revert "
    "to bud-box external X x Y.");
  G4UIparameter* modXP =
    new G4UIparameter("fullX", 'd', false);
  modXP->SetGuidance("full width X");
  modXP->SetParameterRange("fullX>=0.");
  fModeratorXYCmd->SetParameter(modXP);
  G4UIparameter* modYP =
    new G4UIparameter("fullY", 'd', false);
  modYP->SetGuidance("full width Y");
  modYP->SetParameterRange("fullY>=0.");
  fModeratorXYCmd->SetParameter(modYP);
  G4UIparameter* modUnt =
    new G4UIparameter("unit", 's', false);
  modUnt->SetGuidance("length unit");
  modUnt->SetParameterCandidates(G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm")));
  fModeratorXYCmd->SetParameter(modUnt);
  fModeratorXYCmd->AvailableForStates(G4State_PreInit);
  fModeratorXYCmd->SetToBeBroadcasted(false);

  fEncasementEnableCmd = new G4UIcmdWithABool("/testhadr/det/setEncasementEnable", this);
  fEncasementEnableCmd->SetGuidance(
    "If true, build a hollow G4_POLYCARBONATE shell around the bud box. Inner cavity is "
    "air; Al +Z exterior is flush with the inner +Z face of the shell. Requires bud box "
    "(wall thickness > 0). Outer dimensions and wall are set separately.");
  fEncasementEnableCmd->SetParameterName("Enable", false);
  fEncasementEnableCmd->AvailableForStates(G4State_PreInit);
  fEncasementEnableCmd->SetToBeBroadcasted(false);

  fEncasementOuterCmd = new G4UIcmdWith3VectorAndUnit("/testhadr/det/setEncasementOuter", this);
  fEncasementOuterCmd->SetGuidance("Outer dimensions (X Y Z) of the polycarbonate encasement.");
  fEncasementOuterCmd->SetParameterName("xExt", "yExt", "zExt", false);
  fEncasementOuterCmd->SetUnitCategory("Length");
  fEncasementOuterCmd->SetRange("xExt>0. && yExt>0. && zExt>0.");
  fEncasementOuterCmd->AvailableForStates(G4State_PreInit);
  fEncasementOuterCmd->SetToBeBroadcasted(false);

  fEncasementWallCmd =
    new G4UIcmdWithADoubleAndUnit("/testhadr/det/setEncasementWallThickness", this);
  fEncasementWallCmd->SetGuidance("Uniform wall thickness of the hollow polycarbonate shell.");
  fEncasementWallCmd->SetParameterName("Thickness", false);
  fEncasementWallCmd->SetUnitCategory("Length");
  fEncasementWallCmd->SetRange("Thickness>0.");
  fEncasementWallCmd->AvailableForStates(G4State_PreInit);
  fEncasementWallCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fSizeYZCmd;
  delete fNbLayersCmd;
  delete fNbAbsorCmd;
  delete fAbsorCmd;
  delete fIsotopeCmd;
  delete fBoxExternalCmd;
  delete fBoxWallCmd;
  delete fDetStandoffCmd;
  delete fBoxMatCmd;
  delete fCavityMatCmd;
  delete fWorldMatCmd;
  delete fLabSizeCmd;
  delete fSourceStandoffCmd;
  delete fModeratorEnableCmd;
  delete fNbModeratorLayersCmd;
  delete fModeratorLayerCmd;
  delete fModeratorGapCmd;
  delete fModeratorMatCmd;
  delete fModeratorXYCmd;
  delete fEncasementEnableCmd;
  delete fEncasementOuterCmd;
  delete fEncasementWallCmd;
  delete fDetDir;
  delete fTestemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fSizeYZCmd) {
    fDetector->SetCalorSizeYZ(fSizeYZCmd->GetNewDoubleValue(newValue));
  }

  if (command == fNbLayersCmd) {
    fDetector->SetNbOfLayers(fNbLayersCmd->GetNewIntValue(newValue));
  }

  if (command == fNbAbsorCmd) {
    fDetector->SetNbOfAbsor(fNbAbsorCmd->GetNewIntValue(newValue));
  }

  if (command == fAbsorCmd) {
    G4int num;
    G4double tick;
    G4String unt, mat;
    std::istringstream is(newValue);
    is >> num >> mat >> tick >> unt;
    G4String material = mat;
    tick *= G4UIcommand::ValueOf(unt);
    fDetector->SetAbsorMaterial(num, material);
    fDetector->SetAbsorThickness(num, tick);
  }

  if (command == fIsotopeCmd) {
    G4int Z;
    G4int A;
    G4double dens;
    G4String name, unt;
    std::istringstream is(newValue);
    is >> name >> Z >> A >> dens >> unt;
    dens *= G4UIcommand::ValueOf(unt);
    fDetector->MaterialWithSingleIsotope(name, name, dens, Z, A);
  }

  if (command == fBoxExternalCmd) {
    G4ThreeVector v = fBoxExternalCmd->GetNew3VectorValue(newValue);
    fDetector->SetBoxExternal(v.x(), v.y(), v.z());
  }

  if (command == fBoxWallCmd) {
    fDetector->SetBoxWallThickness(fBoxWallCmd->GetNewDoubleValue(newValue));
  }

  if (command == fDetStandoffCmd) {
    fDetector->SetDetectorStandoff(fDetStandoffCmd->GetNewDoubleValue(newValue));
  }

  if (command == fBoxMatCmd) {
    fDetector->SetBoxMaterial(newValue);
  }

  if (command == fCavityMatCmd) {
    fDetector->SetCavityMaterial(newValue);
  }

  if (command == fWorldMatCmd) {
    fDetector->SetWorldMaterial(newValue);
  }

  if (command == fLabSizeCmd) {
    G4ThreeVector v = fLabSizeCmd->GetNew3VectorValue(newValue);
    fDetector->SetLabSize(v.x(), v.y(), v.z());
  }

  if (command == fSourceStandoffCmd) {
    fDetector->SetSourceStandoff(fSourceStandoffCmd->GetNewDoubleValue(newValue));
  }

  if (command == fModeratorEnableCmd) {
    fDetector->SetModeratorEnable(G4UIcmdWithABool::GetNewBoolValue(newValue.c_str()));
  }

  if (command == fNbModeratorLayersCmd) {
    fDetector->SetNbOfModeratorLayers(fNbModeratorLayersCmd->GetNewIntValue(newValue));
  }

  if (command == fModeratorLayerCmd) {
    G4int num;
    G4double tick;
    G4String unt, mat;
    std::istringstream is(newValue);
    is >> num >> mat >> tick >> unt;
    tick *= G4UIcommand::ValueOf(unt);
    fDetector->SetModeratorLayer(num, mat, tick);
  }

  if (command == fModeratorGapCmd) {
    fDetector->SetModeratorUpstreamGap(fModeratorGapCmd->GetNewDoubleValue(newValue));
  }

  if (command == fModeratorMatCmd) {
    fDetector->SetModeratorMaterial(newValue);
  }

  if (command == fModeratorXYCmd) {
    G4double x = 0.;
    G4double y = 0.;
    G4String unt;
    std::istringstream is(newValue);
    is >> x >> y >> unt;
    x *= G4UIcommand::ValueOf(unt);
    y *= G4UIcommand::ValueOf(unt);
    fDetector->SetModeratorTransverseXY(x, y);
  }

  if (command == fEncasementEnableCmd) {
    fDetector->SetEncasementEnable(G4UIcmdWithABool::GetNewBoolValue(newValue.c_str()));
  }

  if (command == fEncasementOuterCmd) {
    G4ThreeVector v = fEncasementOuterCmd->GetNew3VectorValue(newValue);
    fDetector->SetEncasementOuter(v.x(), v.y(), v.z());
  }

  if (command == fEncasementWallCmd) {
    fDetector->SetEncasementWallThickness(fEncasementWallCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
