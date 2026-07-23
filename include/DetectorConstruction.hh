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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4Cache.hh"
#include "G4RotationMatrix.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

class G4GlobalMagFieldMessenger;

const G4int kMaxAbsor = 20;  // 0 + 19
const G4int kMaxModeratorLayers = 10;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

  public:
    G4Material* MaterialWithSingleIsotope(G4String, G4String, G4double, G4int, G4int);

    void SetNbOfAbsor(G4int);
    void SetAbsorMaterial(G4int, const G4String&);
    void SetAbsorThickness(G4int, G4double);

    void SetWorldMaterial(const G4String&);
    void SetCavityMaterial(const G4String&);
    void SetBoxMaterial(const G4String&);
    void SetCalorSizeYZ(G4double);
    void SetCalorSizeY(G4double);
    void SetCalorSizeZ(G4double);
    void SetNbOfLayers(G4int);

    void SetBoxExternal(G4double xExt, G4double yExt, G4double zExt);
    void SetBoxWallThickness(G4double);
    void SetDetectorStandoff(G4double);
    void SetLabSize(G4double xExt, G4double yExt, G4double zExt);
    void SetSourceStandoff(G4double);

    void SetModeratorEnable(G4bool);
    void SetNbOfModeratorLayers(G4int);
    void SetModeratorLayer(G4int, const G4String&, G4double);
    void SetModeratorUpstreamGap(G4double);
    void SetModeratorMaterial(const G4String&);
    void SetModeratorTransverseXY(G4double fullX, G4double fullY);

    /// Optional hollow polycarbonate shell around the bud box (inner cavity = air).
    void SetEncasementEnable(G4bool);
    void SetEncasementOuter(G4double xExt, G4double yExt, G4double zExt);
    void SetEncasementWallThickness(G4double);

    /// Optional solid slab (e.g. concrete shielding) flush against the bud box
    /// -Z (back) exterior face, in the lab world. Same transverse X/Y as the
    /// box; thickness and material set separately. Requires the bud box and is
    /// incompatible with the encasement.
    void SetBackSlabEnable(G4bool);
    void SetBackSlabThickness(G4double);
    void SetBackSlabMaterial(const G4String&);
    /// Full transverse X/Y of the back slab. Pass 0 for either to revert that
    /// dimension to the bud-box external size (the default).
    void SetBackSlabTransverse(G4double fullX, G4double fullY);

    /// Cs-137 check-source envelope: hexagonal prism (axis || global +Z) with a
    /// central cylindrical bore. Used with /gps/pos/confine Cs137Source.
    void SetCs137SourceEnable(G4bool);

    /// Silver-epoxy contact modeled as a localized half-ellipsoid blob (base
    /// radius fEpoxyBlobRadius, height fEpoxyBlobHeight) sitting on the front
    /// face of absorber 1, instead of a full-area box. Absorber 1 must be
    /// SilverEpoxy and its thickness sets the front air region depth.
    void SetSilverEpoxyBlob(G4bool);
    void SetSilverEpoxyBlobSize(G4double radius, G4double height);

    G4bool GetModeratorEnable() const
    {
      return fModeratorEnable || fNbModeratorLayers > 0;
    }
    G4int GetNbOfModeratorLayers() const { return fNbModeratorLayers; }
    G4bool GetEncasementEnable() const { return fEncasementEnable; }

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

  public:
    void PrintCalorParameters();

    G4double GetWorldSizeX() const { return fLabSizeX; };
    G4double GetWorldSizeYZ() const { return fLabSizeY; };

    G4double GetLabSizeX() const { return fLabSizeX; };
    G4double GetLabSizeY() const { return fLabSizeY; };
    G4double GetLabSizeZ() const { return fLabSizeZ; };

    G4double GetBoxExternalX() const { return fBoxExternalX; };
    G4double GetBoxExternalY() const { return fBoxExternalY; };
    G4double GetBoxExternalZ() const { return fBoxExternalZ; };
    G4double GetBoxWallThickness() const { return fBoxWallThickness; };
    G4double GetDetectorStandoff() const { return fDetectorStandoff; };
    G4double GetSourceStandoff() const { return fSourceStandoff; };

    G4double GetCalorThickness() const { return fCalorThickness; };
    G4double GetCalorSizeYZ() const { return fCalorSizeYZ; };
    G4double GetCalorSizeY() const { return fCalorSizeY; };
    G4double GetCalorSizeZ() const { return fCalorSizeZ; };
    G4double GetCalorTransverseArea() const { return fCalorSizeY * fCalorSizeZ; };

    G4int GetNbOfLayers() const { return fNbOfLayers; };

    G4int GetNbOfAbsor() const { return fNbOfAbsor; };
    G4double GetAbsorThickness(G4int i) const { return fAbsorThickness[i]; };
    const G4Material* GetAbsorMaterial(G4int i) const { return fAbsorMaterial[i]; };
    const G4LogicalVolume* GetAbsorberLogical(G4int i) const { return fLogicAbsor[i]; };

    // Per-absorber charge-sensitive (semiconductor) scoring configuration.
    // Absorbers flagged charge-active receive the full LET / track-structure
    // scoring; W is the (estimated) pair-creation energy used only for the
    // initial-pair estimate N0 = Eion / W. Both are 1-based (see LetScoring.hh);
    // defaults are derived from the absorber material (PNDI* -> active, W =
    // kDefaultPairCreationEnergy) in Construct(), and explicit setter calls
    // (e.g. from the messenger) override those defaults.
    G4bool IsChargeActiveAbsorber(G4int i) const
    {
      return (i >= 1 && i < kMaxAbsor) ? fChargeActive[i] : false;
    }
    G4double GetPairCreationEnergy(G4int i) const
    {
      return (i >= 1 && i < kMaxAbsor) ? fPairCreationEnergy[i] : 0.;
    }
    void SetChargeActiveAbsorber(G4int i, G4bool active);
    void SetPairCreationEnergy(G4int i, G4double value);

    // Convergence controls: absorber production range cut and max step limit.
    // Applied in DefineRegionsAndCuts() (Construct-time), so set PreInit. The
    // 0.1 um baseline is the default; sweep coarse/baseline/fine for convergence.
    void SetAbsorberRangeCut(G4double v) { fAbsorberRangeCut = v; }
    void SetAbsorberMaxStep(G4double v) { fAbsorberMaxStep = v; }
    G4double GetAbsorberRangeCut() const { return fAbsorberRangeCut; }
    G4double GetAbsorberMaxStep() const { return fAbsorberMaxStep; }

    // Sampled StepLET ntuple gating (read by SteppingAction). Disabled by
    // default; when enabled, StepLET rows are written for events with
    // eventID < maxStepEvents (maxStepEvents <= 0 means all events).
    void SetWriteStepNtuple(G4bool v) { fWriteStepNtuple = v; }
    void SetMaxStepEvents(G4int n) { fMaxStepEvents = n; }
    G4bool GetWriteStepNtuple() const { return fWriteStepNtuple; }
    G4int GetMaxStepEvents() const { return fMaxStepEvents; }

    const G4VPhysicalVolume* GetphysiWorld() const { return fPhysiWorld; };
    const G4Material* GetWorldMaterial() const { return fWorldMaterial; };
    const G4Material* GetCavityMaterial() const { return fCavityMaterial; };
    const G4Material* GetBoxMaterial() const { return fBoxMaterial; };
    const G4VPhysicalVolume* GetAbsorber(G4int i) const { return fPhysiAbsor[i]; };
    const G4VPhysicalVolume* GetMonitorFront() const { return fPhysiMonitor; };
    const G4LogicalVolume* GetMonitorFrontLogical() const { return fLogicMonitor; };
    const G4VPhysicalVolume* GetModeratorPhysical() const
    {
      return fPhysiModeratorLayer[1];
    }
    G4bool HasModeratorStack() const { return fLogicModeratorLayer[1] != nullptr; }
    const G4LogicalVolume* GetModeratorLogical(G4int layer) const
    {
      return (layer >= 1 && layer <= kMaxModeratorLayers)
                 ? fLogicModeratorLayer[layer]
                 : nullptr;
    }
    G4bool IsModeratorLogical(const G4LogicalVolume* lv) const;

    // Position (global) where the source default sits: +Z at sourceStandoff from
    // the outer +Z face of the bud box, or from the polycarbonate shell if
    // encasement is enabled (see GetDefaultSourcePosition).
    G4ThreeVector GetDefaultSourcePosition() const;

  private:
    void DefineMaterials();
    void ComputeCalorParameters();
    /// Per-volume G4Region production cuts + G4UserLimits on polymer absorbers.
    void DefineRegionsAndCuts();

    G4int fNbOfAbsor = 0;
    G4Material* fAbsorMaterial[kMaxAbsor];
    G4double fAbsorThickness[kMaxAbsor];

    // Per-absorber charge-active flags and (estimated) pair-creation energies,
    // indexed by the 1-based absorber number (index 0 reserved/unused). The
    // *Set arrays record explicit user overrides so Construct() does not clobber
    // messenger-configured values with material-derived defaults.
    G4bool fChargeActive[kMaxAbsor] = {false};
    G4bool fChargeActiveSet[kMaxAbsor] = {false};
    G4double fPairCreationEnergy[kMaxAbsor] = {0.};
    G4bool fPairCreationEnergySet[kMaxAbsor] = {false};

    // Convergence controls (absorber production cut + max step) and StepLET
    // sampled-output gating. Defaults set in the constructor.
    G4double fAbsorberRangeCut = 0.;
    G4double fAbsorberMaxStep = 0.;
    G4bool fWriteStepNtuple = false;
    G4int fMaxStepEvents = 0;

    G4int fNbOfLayers = 0;
    G4double fLayerThickness = 0.;

    G4double fCalorSizeYZ = 0.;
    G4double fCalorSizeY = 0.;
    G4double fCalorSizeZ = 0.;
    G4double fCalorThickness = 0.;

    // Lab world (large air-filled volume containing the bud box)
    G4Material* fWorldMaterial = nullptr;
    G4double fLabSizeX = 0.;
    G4double fLabSizeY = 0.;
    G4double fLabSizeZ = 0.;

    // Aluminum bud box wall (external dimensions)
    G4Material* fBoxMaterial = nullptr;
    G4double fBoxExternalX = 0.;
    G4double fBoxExternalY = 0.;
    G4double fBoxExternalZ = 0.;
    G4double fBoxWallThickness = 0.;

    // Air cavity inside the bud box (filled with cavity material; the
    // calorimeter, layer, and front-face monitor all use this material).
    G4Material* fCavityMaterial = nullptr;

    // Distance from the inside of the box's +Z (front) face to the front
    // face of the detector stack, measured along global +Z.
    G4double fDetectorStandoff = 0.;

    // Default source standoff: distance from the +Z exterior of the outermost
    // solid on the source side (bud box, or polycarbonate encasement if enabled)
    // to the default GPS position, along global +Z.
    G4double fSourceStandoff = 0.;

    // Optional hollow outer shell (e.g. polycarbonate) around the Al bud box.
    // Inner cavity is air; Al +Z exterior is flush with the inner +Z face of
    // the shell. Disabled by default.
    G4bool fEncasementEnable = false;
    G4double fEncasementOuterX = 0.;
    G4double fEncasementOuterY = 0.;
    G4double fEncasementOuterZ = 0.;
    G4double fEncasementWallThickness = 0.;

    // Optional slab (e.g. concrete shielding) flush against the bud box -Z
    // (back) exterior face, in the lab world. Transverse X/Y match the box.
    G4bool fBackSlabEnable = false;
    G4double fBackSlabThickness = 0.;
    G4Material* fBackSlabMaterial = nullptr;
    // Full transverse sizes; 0 means "use the bud-box external size".
    G4double fBackSlabFullX = 0.;
    G4double fBackSlabFullY = 0.;

    // Optional moderator stack in the lab world between the bud box +Z exterior
    // face and the GPS plane. Layer 1 sits adjacent to the box; higher indices
    // stack toward +Z (source). Requires wall thickness > 0.
    //
    // Legacy single-slab mode: setModeratorEnable true with fNbModeratorLayers=0
    // uses one layer of thickness (sourceStandoff - moderatorGap).
    G4bool fModeratorEnable = false;
    G4int fNbModeratorLayers = 0;
    G4Material* fModeratorLayerMat[kMaxModeratorLayers + 1];
    G4double fModeratorLayerThick[kMaxModeratorLayers + 1];
    G4double fModeratorUpstreamGap = 0.;
    G4Material* fModeratorMaterial = nullptr;
    G4double fModeratorFullX = 0.;
    G4double fModeratorFullY = 0.;

    G4Box* fSolidWorld = nullptr;
    G4LogicalVolume* fLogicWorld = nullptr;
    G4VPhysicalVolume* fPhysiWorld = nullptr;

    G4Box* fSolidBox = nullptr;
    G4LogicalVolume* fLogicBox = nullptr;
    G4VPhysicalVolume* fPhysiBox = nullptr;

    G4Box* fSolidCavity = nullptr;
    G4LogicalVolume* fLogicCavity = nullptr;
    G4VPhysicalVolume* fPhysiCavity = nullptr;

    G4Box* fSolidCalor = nullptr;
    G4LogicalVolume* fLogicCalor = nullptr;
    G4VPhysicalVolume* fPhysiCalor = nullptr;

    G4Box* fSolidLayer = nullptr;
    G4LogicalVolume* fLogicLayer = nullptr;
    G4VPhysicalVolume* fPhysiLayer = nullptr;

    G4Box* fSolidAbsor[kMaxAbsor];
    G4LogicalVolume* fLogicAbsor[kMaxAbsor];
    G4VPhysicalVolume* fPhysiAbsor[kMaxAbsor];

    // Thin scoring slab placed immediately in front of (i.e. +Z of) the
    // calorimeter front face, used to score the neutron spectrum entering
    // the detector after any environment interactions.
    G4Box* fSolidMonitor = nullptr;
    G4LogicalVolume* fLogicMonitor = nullptr;
    G4VPhysicalVolume* fPhysiMonitor = nullptr;

    G4Box* fSolidModeratorLayer[kMaxModeratorLayers + 1];
    G4LogicalVolume* fLogicModeratorLayer[kMaxModeratorLayers + 1];
    G4VPhysicalVolume* fPhysiModeratorLayer[kMaxModeratorLayers + 1];

    // Polycarbonate (or other) encasement: inner air volume + hollow shell solid.
    G4Box* fSolidEncasementAir = nullptr;
    G4LogicalVolume* fLogicEncasementAir = nullptr;
    G4VPhysicalVolume* fPhysiEncasementAir = nullptr;
    G4LogicalVolume* fLogicPCShell = nullptr;
    G4VPhysicalVolume* fPhysiPCShell = nullptr;

    // Concrete (or other) back slab flush against the bud box -Z exterior.
    G4Box* fSolidBackSlab = nullptr;
    G4LogicalVolume* fLogicBackSlab = nullptr;
    G4VPhysicalVolume* fPhysiBackSlab = nullptr;

    // Cs-137 check source: hex shell (G4Polyhedra) minus central bore (G4Tubs).
    G4bool fCs137SourceEnable = false;
    G4VSolid* fSolidCs137Source = nullptr;
    G4LogicalVolume* fLogicCs137Source = nullptr;
    G4VPhysicalVolume* fPhysiCs137Source = nullptr;
    G4RotationMatrix* fCs137SourceRotation = nullptr;

    // Silver-epoxy contact as a localized half-ellipsoid on the front absorber.
    // When enabled, absorber 1's box region is filled with cavity material and a
    // half-ellipsoid of absorber 1's material (SilverEpoxy) is placed inside it;
    // fLogicAbsor[1] is set to the ellipsoid so scoring targets the epoxy.
    G4bool fSilverEpoxyBlob = false;
    G4double fEpoxyBlobRadius = 0.;  // set in constructor (default 2 mm)
    G4double fEpoxyBlobHeight = 0.;  // set in constructor (default 1 mm)
    G4LogicalVolume* fLogicEpoxyContainer = nullptr;
    G4RotationMatrix* fEpoxyBlobRotation = nullptr;

    // Rotation applied to the calorimeter so its (locally X-stacked) layout
    // becomes Z-stacked in the global frame, with the first absorber facing
    // the +Z (source) side. Owned by this class.
    G4RotationMatrix* fCalorRotation = nullptr;

    DetectorMessenger* fDetectorMessenger = nullptr;
    G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
