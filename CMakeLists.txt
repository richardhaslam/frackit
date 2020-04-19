# minimum required cmake version
cmake_minimum_required(VERSION 3.0)

# project name
project(frackit VERSION 1.0.0)

# project version
set(FRACKIT_VERSION ${frackit_VERSION} CACHE STRING "project version")
set(FRACKIT_VERSION_MAJOR ${frackit_VERSION_MAJOR} CACHE STRING "major project version")
set(FRACKIT_VERSION_MINOR ${frackit_VERSION_MINOR} CACHE STRING "minor project version")

# make sure our own modules are found
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules")

#include the macros
include(FrackitMacros)

# compiler definitions
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS} -std=c++0x -O3 -g -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "$GXX_RELEASE_OPTS $GXX_RELEASE_WARNING_OPTS -g -ggdb -Wall")
set(CMAKE_BUILD_TYPE Release)

# Find OpenCascade
find_path(OCC_INC "Standard_Version.hxx"
          PATH_SUFFIXES opencascade include/opencascade occt/)

if (OCC_INC)
    include_directories( ${OCC_INC} )

    # get OCC version
    file(STRINGS ${OCC_INC}/Standard_Version.hxx
         OCC_MAJOR REGEX "#define OCC_VERSION_MAJOR.*")
    file(STRINGS ${OCC_INC}/Standard_Version.hxx
         OCC_MINOR REGEX "#define OCC_VERSION_MINOR.*")
    file(STRINGS ${OCC_INC}/Standard_Version.hxx
         OCC_MAINT REGEX "#define OCC_VERSION_MAINTENANCE.*")
    if(OCC_MAJOR AND OCC_MINOR AND OCC_MAINT)
      string(REGEX MATCH "[0-9]+" OCC_MAJOR "${OCC_MAJOR}")
      string(REGEX MATCH "[0-9]+" OCC_MINOR "${OCC_MINOR}")
      string(REGEX MATCH "[0-9]+" OCC_MAINT "${OCC_MAINT}")
      set(OCC_VERSION "${OCC_MAJOR}.${OCC_MINOR}.${OCC_MAINT}")
      message(STATUS "Found OpenCASCADE version ${OCC_VERSION} in ${OCC_INC}")
    endif()

    # link required libraries
    set(OCC_LIBS_REQUIRED
        # ModelingAlgorithms
        TKBO TKTopAlgo TKBRep TKGeomAlgo TKPrim
        # ModelingData
        TKGeomBase TKG3d TKG2d
        # FoundationClasses
        TKMath TKernel)
    list(LENGTH OCC_LIBS_REQUIRED NUM_OCC_LIBS_REQUIRED)

    if(OCC_LIBS)
      message(STATUS "OCC libraries specified explicitly: " ${OCC_LIBS})
      list(LENGTH OCC_LIBS_REQUIRED NUM_OCC_LIBS)
    else()
      set(OCC_LIBS)
      foreach(OCC ${OCC_LIBS_REQUIRED})
        find_library(OCC_LIB ${OCC}
                     PATH_SUFFIXES lib ${OCC_SYS_NAME}/lib ${OCC_SYS_NAME}/vc8/lib)
        if(OCC_LIB)
	        list(APPEND OCC_LIBS ${OCC_LIB})
        else()
	        message(STATUS "OCC lib " ${OCC} " not Found")
        endif()
        unset(OCC_LIB CACHE)
      endforeach()
      list(LENGTH OCC_LIBS NUM_OCC_LIBS)
    endif()

    if(NUM_OCC_LIBS EQUAL NUM_OCC_LIBS_REQUIRED)
      set(HAVE_OCC TRUE)
      list(APPEND CONFIG_OPTIONS "OpenCASCADE")
      message(STATUS "Found " "OpenCASCADE")
      list(APPEND EXTERNAL_LIBRARIES ${OCC_LIBS})
      list(APPEND EXTERNAL_INCLUDES ${OCC_INC})
      if(HAVE_64BIT_SIZE_T)
        add_definitions(-D_OCC64)
      endif()
    endif()
else()
    message(FATAL_ERROR "OpenCascade sources not found")
endif()

# indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc_doxygen
                       COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       COMMENT "Generating API documentation with Doxygen"
                       VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

# include the sources
include_directories( ${CMAKE_SOURCE_DIR} )

# include subdirectories
add_subdirectory(appl)
add_subdirectory(cmake)
add_subdirectory(doc)
add_subdirectory(frackit)
add_subdirectory(test)
