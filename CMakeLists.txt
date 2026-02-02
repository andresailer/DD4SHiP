cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

#---------------------------
set( PackageName DD4SHIP )
#---------------------------

project(${PackageName} LANGUAGES CXX)

# project version
SET( ${PackageName}_VERSION_MAJOR 0 )
SET( ${PackageName}_VERSION_MINOR 1 )
SET( ${PackageName}_VERSION_PATCH 0 )

SET( ${PackageName}_VERSION "${${PackageName}_VERSION_MAJOR}.${${PackageName}_VERSION_MINOR}" )
SET( ${PackageName}_SOVERSION "${${PackageName}_VERSION_MAJOR}.${${PackageName}_VERSION_MINOR}" )

set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${PROJECT_SOURCE_DIR}/cmake )
set(LIBRARY_OUTPUT_PATH    ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

include(GNUInstallDirs)
set(CMAKE_INSTALL_LIBDIR lib)
set(CMAKE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

#------------- set the default installation directory to be the source directory

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET( CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR} CACHE PATH
    "install prefix path  - overwrite with -D CMAKE_INSTALL_PREFIX = ..."
    FORCE )
  MESSAGE(STATUS "CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX} - overwrite with -D CMAKE_INSTALL_PREFIX" )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)



#---Options-------------------------------------------------------------------------

option(BUILD_TESTING "Enable and build tests" ON)
option(INSTALL_COMPACT_FILES "Copy compact files to install area" OFF)
option(INSTALL_BEAMPIPE_STL_FILES "Download CAD files for building the detailed beampipe" OFF)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

find_package(DD4hep REQUIRED COMPONENTS DDRec DDG4 DDParsers)

find_package ( ROOT REQUIRED COMPONENTS Geom GenVector)
message ( STATUS "ROOT_VERSION: ${ROOT_VERSION}" )

find_package( Geant4 REQUIRED ) 
file(GLOB sources
  ./src/*.cpp
  )

add_dd4hep_plugin(${PackageName} SHARED ${sources})

target_link_libraries(${PackageName} DD4hep::DDCore DD4hep::DDRec DD4hep::DDParsers ROOT::Core)

#Create this_package.sh file, and install
dd4hep_instantiate_package(${PackageName})

# Destination directories are hardcoded because GNUdirectories are not included
install(TARGETS ${PackageName}
  EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT shlib
)
