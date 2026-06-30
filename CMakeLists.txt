#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(Hadr05)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(Hadr05 Hadr05.cc ${sources} ${headers})
target_link_libraries(Hadr05 ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build Hadr05. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(Hadr05_SCRIPTS
    Cu-lAr.mac Fe-Sci.mac  Pb-lAr-em.mac Pb-lAr.mac W-lAr.mac
    emtest.mac emtutor.mac hadr05.in plotHisto.C vis.mac
    fast-detector.mac thermal-detector.mac stacked-detector.mac test.mac test2.mac
    validation.mac
    dd-test.mac
    dd-water.mac
    dd-water-vis.mac
    dd-reactor.mac
    fission-reactor.mac
    fission-reactor-vis.mac
    dd-moderated.mac
    dd-moderated-vis.mac
    dd-reactor-vis.mac
    dd-test-vis.mac
    shine-dd.mac
    shine-dd-vis.mac
    shine-dd-box.mac
    shine-dd-box-vis.mac
    shine-dd-batch-example.mac
    shine-dd-box-batched.mac
    shine-box-spectrum.mac
    shine-box-spectrum-vis.mac
    shine-box-spectrum-single.mac
    shine-dd-box-sectrace.mac
    dd-moderated-vis.mac
    dd-moderator-stack.mac
    dd-moderator-stack-vis.mac
    dd-moderator-stack-batched.mac
    test3.mac
    rotation-test.mac
    rotation-test-vis.mac
    perpendicular-test.mac
    perpendicular-test-vis.mac
    mixedsource.mac
    mixedsource-vis.mac
    Cs-137.mac
    Cs-137-vis.mac
    PuBe.mac
    PuBe-vis.mac
    PuBeGamma.mac
    Am-241.mac
    dd-single-moderator-vis.mac
    dd-single-moderator.mac
    Am-241-test.mac
  )

foreach(_script ${Hadr05_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS Hadr05 DESTINATION bin)

