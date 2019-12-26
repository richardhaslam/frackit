# minimum required cmake version
cmake_minimum_required(VERSION 2.8)

# project name
project(frackit)

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
find_path(OCC_INC "Standard_Version.hxx" HINTS ~/Software/opencascade-7.3.0/occ-build/include/opencascade)
          # PATH_SUFFIXES opencascade include/opencascade)

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
        find_library(OCC_LIB ${OCC} HINTS ~/Software/opencascade-7.3.0/occ-build/lin64/gcc/lib)
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

# include the sources
include_directories( ${CMAKE_SOURCE_DIR} )

# include subdirectories
add_subdirectory(cmake)
add_subdirectory(frackit)
add_subdirectory(test)
