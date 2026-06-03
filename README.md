\page ExampleHadr05 Example Hadr05

 How to collect energy deposition in a sampling calorimeter.
 How to survey energy flow.
 Hadr05 is the hadronic equivalent of TestEm3.


## GEOMETRY DEFINITION
 
  The calorimeter is a box made of a given number of layers.
  A layer consists of a sequence of various absorbers (maximum MaxAbsor=9).
  The layer is replicated.
 
  Parameters defining the calorimeter :
    - the number of layers,
    - the number of absorbers within a layer,		
    - the material of the absorbers,
    - the thickness of the absorbers,
    - the transverse size of the calorimeter (the input face is a square). 
 
  In addition a transverse uniform magnetic field can be applied.
 
  The default geometry is constructed in DetectorConstruction class, but all
  of the above parameters can be modified interactively via the commands 
  defined in the DetectorMessenger class.

```
        |<----layer 0---------->|<----layer 1---------->|<----layer 2---------->|
        |           |           |                       |                       |
        ==========================================================================
        ||          |           ||          |           ||          |           ||
        ||          |           ||          |           ||          |           ||
        ||   abs 1  | abs 2     ||   abs 1  | abs 2     ||   abs 1  | abs 2     ||
        ||          |           ||          |           ||          |           ||
        ||          |           ||          |           ||          |           ||
 beam   ||          |           ||          |           ||          |           ||
======> ||          |           ||          |           ||          |           ||
        ||          |           ||          |           ||          |           ||
        ||          |           ||          |           ||          |           ||
        ||          |           ||          |           ||          |           ||
        ||          |           ||          |           ||          |           ||
        ||   cell 1 | cell 2    ||   cell 3 | cell 4    ||   cell 5 | cell 6    ||
        ==========================================================================
        ^           ^           ^           ^           ^           ^           ^
        pln1        pln2        pln3       pln4        pln5        pln6       pln7
```
 
  NB. The number of absorbers and the number of layers can be set to 1.
  In this case we have a unique homogeneous block of matter, which looks like 
  a bubble chamber rather than a calorimeter ...
  (see the macro emtutor.mac)
  
  A function, and its associated UI command, allows to build a material
  directly from a single isotope.
  
  To be identified by the ThermalScattering module, the elements composing a
  material must have a specific name (see G4ParticleHPThermalScatteringNames.cc)
  Examples of such materials are build in Hadr06/src/DetectorConstruction.cc
    
## PHYSICS LISTS

  The executable uses the Geant4 **reference** modular physics list
  `Shielding_HPT_EMZ`, instantiated in `Hadr05.cc` via
  `G4PhysListFactory::GetReferencePhysList("Shielding_HPT_EMZ")`.

  - **Shielding** — Reference hadronic configuration for shielding and
    neutron applications: ParticleHP (ENDF-based) for neutrons below 20 MeV
    (elastic, inelastic, radiative capture, fission where applicable),
    Bertini cascade and FTFP at higher energies, plus the usual ion,
    stopping, decay, gamma-nuclear, and radioactive-decay constructors bundled
    with the Shielding family.
  - **HPT** — Thermal neutron scattering on bound atoms: the
    `G4NeutronHPThermalScattering` process and S(alpha,beta) data from `G4NDL`
    are registered for neutrons below about 4 eV. Materials must use the
    canonical `TS_*` element names for the kernels to apply (see below).
  - **EMZ** — `G4EmStandardPhysics_option4`, the most accurate *standard*
    electromagnetic package in the reference-list factory (Livermore /
    Penelope low-energy models, advanced msc, atomic deexcitation). There is
    no higher-precision **standard** EM variant in `G4PhysListFactory`; lists
    such as `ShieldingLEND_HP` swap evaluated neutron data (LEND) for some
    isotopes rather than strictly increasing precision everywhere.

  **High-precision proton / light-ion transport (QGSP_BIC_AllHP libraries).**
  After the reference list is built, `Hadr05.cc` calls `ReplacePhysics()` twice
  to swap in only the proton and light-ion inelastic libraries from
  `QGSP_BIC_AllHP`, leaving the Shielding HP neutron treatment and the HPT
  thermal scattering untouched:

  - `HadronInelasticShieldingBIC` (in `src/`) inherits `G4HadronPhysicsShielding`
    — so the neutron HP inelastic / capture / fission models are kept verbatim —
    and overrides **only** `Proton()` with the `QGSP_BIC_AllHP` proton chain:
    data-driven `G4ParticleHPInelastic` below 200 MeV plus Binary cascade
    (190 MeV–1.5 GeV), Bertini / FTFP / QGSP above.
  - `G4IonPhysicsPHP` replaces the default ion inelastic with `ParticleHP`
    inelastic for d, t, He3, alpha below 200 MeV/n plus Binary cascade.

  > **Data requirement:** these `ParticleHP` models read the **G4TENDL** data
  > set (`G4TENDL1.4` for Geant4 11.4). It is an *optional* dataset and is **not**
  > auto-discovered through `GEANT4_DATA_DIR`, so the `G4PARTICLEHPDATA`
  > environment variable must point at it. This install's `bin/geant4.sh` now
  > exports `G4PARTICLEHPDATA=$GEANT4_DATA_DIR/G4TENDL1.4` automatically when the
  > directory exists, so simply sourcing `geant4.sh` is sufficient. `G4NDL`
  > (`G4NEUTRONHPDATA`) only covers neutrons.
  >
  > `Hadr05.cc::CheckParticleHPData()` prints a clear warning at startup if the
  > data is not visible (without it Geant4 aborts during `/run/initialize`).
  > Note: the variable must be set in the **shell** before launching — do not
  > `setenv()` it from inside the program, as a late environment reallocation
  > can invalidate Geant4's cached `G4NEUTRONHPDATA` pointer and cause
  > intermittent crashes in multithreaded mode.

  After the list is built, `G4EmParameters` step functions,
  `SetDeexcitationIgnoreCut(true)`, atomic relaxation (fluorescence + Auger +
  Auger-cascade + PIXE), and the Mott e- backscatter correction are applied in
  `main` (the EMZ constructor alone does not set the ion / mu-had stepping or
  these extras). Nuclear de-excitation uses
  `G4DeexPrecoParameters::SetCorrelatedGamma(true)` so (n,gamma) and inelastic
  gamma cascades carry realistic line multiplicities / correlations.

  Legacy sources under `src/PhysicsList.cc` and `HadronElasticPhysicsHP.*`
  are **not** linked from `main`; they remain in the tree for reference or
  local experiments. The UI command `/testhadr/phys/thermalScattering` from
  the old custom elastic physics is **not** available when using the
  reference list — thermal scattering follows the Shielding_HPT configuration.

  Several hadronic physics options are controlled by environment variables.
  To select them, see Hadr07.cc

## EM PROCESSES (Shielding_HPT_EMZ)

  The **EMZ** suffix selects `G4EmStandardPhysics_option4`, which registers
  the standard charged-particle EM processes including, for the relevant
  particle types: **G4hMultipleScattering**, **G4hIonisation**,
  **G4ionIonisation** (with `G4IonParametrisedLossModel` for `GenericIon`),
  **G4NuclearStopping** (e.g. on alpha / heavy ions), **G4eMultipleScattering**,
  **G4eIonisation**, **G4eBremsstrahlung**, **G4eplusAnnihilation**, and the
  low-energy gamma processes (Livermore photoelectric, Compton, conversion,
  Rayleigh). You can confirm at runtime with `/process/list` after
  `/run/initialize`.

## PRODUCTION CUTS AND STEP LIMITS (thin-film stack)

  In `Hadr05.cc`, `SetDefaultCutValue(10 mm)` sets a coarse baseline for the
  large lab / air / water volumes. **Region-specific** range cuts and step
  limits are attached in `DetectorConstruction::DefineRegionsAndCuts()`:

  - **All absorbers** (every layer in the stack, any material): **0.1 µm**
    production cut for gamma, e-, e+, and proton; **G4UserLimits** max step
    **0.1 µm** (never forces a coarser step than the physics would already use).
  - **Metal contacts** (`G4_Ag/Au/Cu/Al/Fe/W/Pb/Ti` and the locally-defined
    `Gold/Copper/Aluminium/Iron/Tungsten/Lead/Titanium`): same as above, but an
    even tighter **0.01 µm** gamma cut so high-Z photoelectric K/L lines are
    resolved.
  - **MonitorFront** slab (sits at the detector front interface): **0.1 µm**
    for gamma / e- / e+ / proton and **0.1 µm** max step.
  - **Bud box** aluminum shell: **5 µm** gamma / e± (passive environment).
  - **Polycarbonate encasement** (if enabled): **5 µm** gamma / e±.
  - **Moderator** layers (if placed): **5 mm** (environment class).

  Step limits require `G4StepLimiterPhysics`. The `Shielding_HPT_EMZ` reference
  list does **not** register one, so `Hadr05.cc` adds it explicitly via
  `RegisterPhysics(new G4StepLimiterPhysics())`; without it the absorber
  `G4UserLimits` max-step caps would be silently ignored.

## LAB WORLD AND MODERATOR STACK

  The **default lab/world material is `G4_AIR`** (`DetectorConstruction` constructor).
  Water immersion (`G4_WATER`) is opt-in per macro via `/testhadr/det/setWorldMaterial`.

  Moderators are optional slabs in the lab world between the bud-box **+Z exterior**
  and the GPS plane. Neutrons travel from the source (+Z) toward the detector (-Z).

  **Multi-layer stack** (preferred): layer 1 sits adjacent to the box; each higher
  index stacks toward +Z (source). Configure before `/run/initialize`:

```
/testhadr/det/setNbOfModeratorLayers 3
/testhadr/det/setModeratorLayer 1 G4_Pb 5 mm
/testhadr/det/setModeratorLayer 2 HeavyWater 85 mm
/testhadr/det/setModeratorLayer 3 polyethylene 95 mm
/testhadr/det/setModeratorGap 10 mm
/testhadr/det/setSourceStandoff 200 mm
```

  `sourceStandoff` must be **≥ sum(layer thicknesses) + moderatorGap**. See
  [`dd-moderator-stack.mac`](dd-moderator-stack.mac) for a full example.

  **Legacy single-slab mode** (backward compatible): `/testhadr/det/setModeratorEnable
  true` with `/testhadr/det/setModeratorMaterial` and no explicit layer count places
  one slab of thickness `sourceStandoff - moderatorGap`. See [`dd-water.mac`](dd-water.mac).

## THERMAL-SCATTERING-AWARE MATERIALS

  `G4ParticleHPThermalScattering` requires elements with the canonical
  `TS_*` names that map to the S(alpha,beta) tables in `G4NDL` (see
  `G4ParticleHPThermalScatteringNames.cc`). `DetectorConstruction.cc`
  defines two ready-to-use variants:

  - `water_TS`        — H2O with `TS_H_of_Water` (use instead of `G4_WATER`)
  - `polyethylene_TS` — (C2H4)n with `TS_H_of_Polyethylene`
  - `polyethylene`    — same TS hydrogen as `polyethylene_TS`
  - `HeavyWater`      — D2O with `TS_D_of_Heavy_Water` for thermal D scattering

  Macros can opt in via, e.g.,

```
/testhadr/det/setWorldMaterial      water_TS
/testhadr/det/setModeratorLayer 1 G4_Pb 5 mm
/testhadr/det/setModeratorLayer 2 HeavyWater 85 mm
/testhadr/det/setModeratorLayer 3 polyethylene 95 mm
/testhadr/det/setModeratorMaterial  water_TS
/testhadr/det/setModeratorMaterial  polyethylene_TS
```

  The standard `G4_WATER` / `polyethylene` materials remain available and
  unchanged; they will still be transported correctly, but without the
  bound-atom S(alpha,beta) treatment for thermal neutrons.
    
## AN EVENT : THE PRIMARY GENERATOR
 
  The primary kinematic consists of a single particle which hits the calorimeter
  perpendicular to the input face. The type of the particle and its energy are 
  set in the PrimaryGeneratorAction class, and can be changed via the 
  G4 build-in commands of G4ParticleGun class (see the macros provided with this 
  example).
 	
  In addition one can choose randomly the impact point of the incident particle.
  The corresponding interactive command is built in PrimaryGeneratorAction.
 	
  A RUN is a set of events.
  
  Hadr05 computes the energy deposited per absorber and the energy flow through
  the calorimeter.
 				
## VISUALIZATION
 
  The Visualization Manager is set in the main() (see Hadr05.cc).
  The initialisation of the drawing is done via the commands :
  /vis/... in the macro vis.mac. In interactive session:
```
  PreInit or Idle > /control/execute vis.mac
```
 	
  The default view is a longitudinal view of the calorimeter.
 	
## PHYSICS DEMO
 
  The particle's type and the physics processes which will be available
  in this example are set in PhysicsList class.
 	
  In addition a built-in interactive command (/process/inactivate processName)
  allows to activate/inactivate the processes one by one.
  Then one can well visualize the processes one by one, especially 
  in the bubble chamber setup with a transverse magnetic field.
 
## HOW TO START ?
 
  - Execute Hadr05 in 'batch' mode from macro files
```
% ./Hadr05  Cu-lAr.mac
```
 
  - Execute Hadr05 in 'interactive mode' with visualization
```
% ./Hadr05
....
Idle> type your commands. For instance:
Idle> /control/execute vis.mac
....
Idle> exit
```

  Macros provided in this example:
  - hadr05.in: macro used in Geant4 testing
  - Fe-Sci.mac, Cu-lAr.mac, Pb-lAr.mac, W-lAr.mac : names are self explanatory
  - Pb-lAr-em.mac : electromagnetic calorimeter
  - emtest.mac, emtutor.mac : to be run interactively
  - vis.mac: to activate visualization

## HISTOGRAMS
 
 Hadr05 can produce histograms :

```
  histo 1 : energy deposit in absorber 1
  histo 2 : energy deposit in absorber 2
  ...etc...........
    
  histo 11 : longitudinal profile of energy deposit in absorber 1 (MeV/event)
  histo 12 : longitudinal profile of energy deposit in absorber 2 (MeV/event)  
  ...etc...........  
  
  histo 21 : energy flow (MeV/event)
 
  histo 22 : total energy deposited 
  histo 23 : total energy leakage  
  histo 24 : total energy released : Edep + Eleak
        
  NB. Numbering scheme for histograms:
  layer     : from 1 to NbOfLayers (included)
  absorbers : from 1 to NbOfAbsor (included)
  planes    : from 1 to NbOfLayers*NbOfAbsor + 1 (included)
```
  
 One can control the binning of the histo with the command:
```
/analysis/h1/set   idAbsor  nbin  Emin  Emax  unit 
```
  etc.,  
  where unit is the desired energy unit for that histo (see Hadr05.in).
  
  One can control the name of the histograms file with the command:
```
/analysis/setFileName  name  (default hadr05)
```
   
  It is possible to choose the format of the histogram file : root (default),
  xml, csv, by using namespace in HistoManager.hh 
 
 It is also possible to print selected histograms on an ascii file:
```
/analysis/h1/setAscii id
```
 All selected histos will be written on a file name.ascii  (default hadr05)

## BATCH CSV SUMMARIES

  Multiple `/run/beamOn` calls in one session can append **per-absorber,
  per-particle** energy-deposit totals to a CSV file for offline averaging
  (e.g. independent batches with different RNG seeds).

```
/testhadr/run/setBatchSummaryEnable true           # default false; gate CSV writing
/testhadr/run/setBatchSummaryFile edep_batches.csv # target path (writing needs enable + non-empty path)
/testhadr/run/setBatchTag my_label                 # optional; set before each beamOn
/run/beamOn 100000
```

  Columns: `run_id`, `n_events`, `batch_tag`, `absorber`, `material`, `particle`,
  `total_edep_MeV`. Strings are CSV-quoted. The header is
  written when the file is new or empty. If an existing file's first line does not
  match this header, that run skips CSV append (rename/delete the file or pick a
  new path).

  For distinct ROOT output per `/run/beamOn` without repeating `/analysis/setFileName`:

```
/testhadr/batch/setRotateRootPerRun true
/testhadr/batch/setRootFileBase hadr05_batch   # writes hadr05_batch0.root, ...
```

  See `shine-dd-batch-example.mac` for two batches with distinct seeds, CSV output,
  and per-run ROOT files via those commands.

  For many batches (100), see `shine-dd-box-batched.mac`.


  