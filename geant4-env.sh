#!/bin/sh
# Source before running hadr05 in a fresh shell:
#   source /path/to/NeutronGammaComplete/geant4-env.sh
#
# geant4make.sh sets G4NEUTRONHPDATA, G4LEDATA, etc. G4PARTICLEHPDATA (G4TENDL)
# is required by this example's proton/light-ion HP physics and is included there
# after GEANT4_INSTALL_DATASETS_TENDL=ON was enabled at build time.

G4MAKE=/Users/samuelangus/geant4-v11.3.2/build/geant4make.sh
if [ ! -f "$G4MAKE" ]; then
  echo "geant4-env.sh: missing $G4MAKE" >&2
  return 1 2>/dev/null || exit 1
fi

# shellcheck disable=SC1090
. "$G4MAKE"

# Fallback if an older geant4make.sh lacks G4PARTICLEHPDATA.
if [ -z "${G4PARTICLEHPDATA-}" ]; then
  if [ -d /Users/samuelangus/geant4-v11.3.2/datasets/G4TENDL1.4 ]; then
    export G4PARTICLEHPDATA=/Users/samuelangus/geant4-v11.3.2/datasets/G4TENDL1.4
  elif [ -d /Users/samuelangus/geant4-v11.3.2/build/data/G4TENDL1.4 ]; then
    export G4PARTICLEHPDATA=/Users/samuelangus/geant4-v11.3.2/build/data/G4TENDL1.4
  else
    echo "geant4-env.sh: G4TENDL not found; run: cd geant4-v11.3.2/build && cmake -DGEANT4_INSTALL_DATASETS_TENDL=ON . && cmake --build . --target G4TENDL" >&2
    return 1 2>/dev/null || exit 1
  fi
fi
