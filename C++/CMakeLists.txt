#=========================================================================
#Program:   Diffusion Weighted MRI Reconstruction
#Module:    $RCSfile: CMakeLists.txt,v $
#Language:  C++
#Author:    $Author: bjian $
#Date:      $Date: 2008/09/10 02:12:44 $
#Version:   $Revision: 1.26 $
#=========================================================================

CMAKE_MINIMUM_REQUIRED(VERSION 1.6)
PROJECT(DW-MRI)

INCLUDE_REGULAR_EXPRESSION(".*")

# Tell CMake's FIND_PACKAGE command where to look for our Find<pkg>
# modules.
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${DW-MRI_SOURCE_DIR})

IF(NOT WIN32)
  OPTION(BUILD_SHARED_LIBS "Build Application with shared libraries." OFF)
ENDIF(NOT WIN32)

# Don't build anything unless the version of CMake is high enough.
# The matching ELSE/ENDIF should be the last lines in the file.
IF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 1.5)



FIND_PACKAGE(VXL)
IF(VXL_FOUND)
  INCLUDE(${VXL_CMAKE_DIR}/UseVXL.cmake)
ENDIF(VXL_FOUND)


INCLUDE_DIRECTORIES( ${VXL_VNL_INCLUDE_DIR} )

IF(WIN32)
  ADD_DEFINITIONS(-DWIN32)
ENDIF(WIN32)

ADD_LIBRARY(port_ini STATIC port_ini.c)


OPTION(BUILD_INI_SHARED "Build the port_ini libraries shared." OFF)
OPTION(BUILD_REC_SHARED "Build the Reconstructor libraries shared." OFF)
OPTION(BUILD_API_SHARED "Build the Reconstructor API libraries shared." OFF)

IF(BUILD_INI_SHARED)
  SET(INI_LIBRARY_TYPE SHARED)
ELSE(BUILD_INI_SHARED)
  SET(INI_LIBRARY_TYPE STATIC)
ENDIF(BUILD_INI_SHARED)

IF(BUILD_REC_SHARED)
  SET(REC_LIBRARY_TYPE SHARED)
ELSE(BUILD_REC_SHARED)
  SET(REC_LIBRARY_TYPE STATIC)
ENDIF(BUILD_REC_SHARED)

IF(BUILD_API_SHARED)
  SET(API_LIBRARY_TYPE SHARED)
ELSE(BUILD_API_SHARED)
  SET(API_LIBRARY_TYPE STATIC)
ENDIF(BUILD_API_SHARED)


SET(RECONSTRUCTION_SRCS
    Reconstructor.cpp
    QBIReconstructor.cpp
    DOTReconstructor.cpp
    MoWReconstructor.cpp
    utilities.cpp
    spherical_harmonics.cpp
    special_functions.cpp
    nnls.cc
    myendian.c
    FltFile.cpp
)


IF(WIN32)
    ADD_DEFINITIONS(-DWIN32 -D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
ELSE(WIN32)  
    ADD_LIBRARY(port_ini ${INI_LIBRARY_TYPE} port_ini.c)  
ENDIF(WIN32)



# generate STATIC library
ADD_LIBRARY(reconstructor ${REC_LIBRARY_TYPE} ${RECONSTRUCTION_SRCS})
TARGET_LINK_LIBRARIES(reconstructor vnl_algo vnl vcl)
IF(UNIX)
    TARGET_LINK_LIBRARIES(reconstructor m port_ini)
ENDIF(UNIX)

#ADD_LIBRARY(reconstructor_api ${API_LIBRARY_TYPE} reconstructor_api.cpp)
#TARGET_LINK_LIBRARIES(reconstructor_api reconstructor)

ADD_EXECUTABLE(test_dot test_dot.cpp)
TARGET_LINK_LIBRARIES(test_dot reconstructor)

ADD_EXECUTABLE(test_mow test_mow.cpp)
TARGET_LINK_LIBRARIES(test_mow reconstructor)

#ADD_EXECUTABLE(test_interp test_interp.cpp)
#TARGET_LINK_LIBRARIES(test_interp reconstructor)

ADD_EXECUTABLE(test_qbi test_qbi.cpp )
TARGET_LINK_LIBRARIES(test_qbi reconstructor)

ADD_EXECUTABLE(reconstruction reconstruction.cpp )
TARGET_LINK_LIBRARIES(reconstruction reconstructor)

#ADD_EXECUTABLE(test_interpolation test_interpolation.cpp)
#TARGET_LINK_LIBRARIES(test_interpolation port_ini vnl_algo vnl vcl)

#ADD_EXECUTABLE(computeS0 computeS0.cpp)
#TARGET_LINK_LIBRARIES(computeS0 vnl_algo vnl vcl)

#ADD_EXECUTABLE(test_api test_api.c)
#TARGET_LINK_LIBRARIES(test_api reconstructor_api)

#ADD_EXECUTABLE(testNNLS testNNLS.cpp nnls.cc utilities.cpp)
#TARGET_LINK_LIBRARIES(testNNLS vnl_algo vnl vcl)

#IF(UNIX)
#  TARGET_LINK_LIBRARIES(flt2itk4d m pthread g2c )
#ENDIF(UNIX)

INSTALL_TARGETS(/bin port_ini test_qbi test_dot test_mow)

#make install
SET(CMAKE_INSTALL_PREFIX ${DW-MRI_SOURCE_DIR})

# If the version of CMake was too old, complain and build nothing.
# These should be the last lines in this file.
ELSE("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 1.5)
  MESSAGE(SEND_ERROR
          "This version of CMake is too old to build InsightApplications.  "
          "Please upgrade to at least CMake 1.6.")
ENDIF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 1.5)
