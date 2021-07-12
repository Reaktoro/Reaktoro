# Include the ExternalProject module
include(ExternalProject)

# Set the list of compiler flags for the external projects
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(CFLAGS "${CFLAGS} -fPIC")
    set(CXXFLAGS "${CXXFLAGS} -fPIC")
else()
    set(CFLAGS "${CMAKE_C_FLAGS}")
    set(CXXFLAGS "${CMAKE_CXX_FLAGS}")
endif()

if (NOT ${REAKTORO_THIRDPARTY_EXTRA_BUILD_ARGS} STREQUAL "")
    separate_arguments(REAKTORO_THIRDPARTY_EXTRA_BUILD_ARGS NATIVE_COMMAND "${REAKTORO_THIRDPARTY_EXTRA_BUILD_ARGS}")
endif()

# Set the common cmake arguments to all external projects
set(REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS
    -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
    $<$<BOOL:CMAKE_GENERATOR_PLATFORM>:-DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}>
    -DCMAKE_INSTALL_PREFIX=${REAKTORO_THIRDPARTY_INSTALL_PREFIX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_C_COMPILER_LAUNCHER=${CMAKE_C_COMPILER_LAUNCHER}
    -DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER}
    -DCMAKE_C_FLAGS=${CFLAGS}
    -DCMAKE_CXX_FLAGS=${CXXFLAGS}
    -DCMAKE_INSTALL_INCLUDEDIR=include
    -DCMAKE_INSTALL_LIBDIR=lib
    -DCMAKE_INSTALL_BINDIR=bin
    -DBUILD_SHARED_LIBS=OFF
    -DCMAKE_CXX_STANDARD=14
    -DCMAKE_CXX_STANDARD_REQUIRED=ON
    # Disables the cmake warning message `Manually-specified
     # variables were not used by the project` because of
     # unused arguments above
    --no-warn-unused-cli
    # Allows the specification of extra args via command-line
    ${REAKTORO_THIRDPARTY_EXTRA_BUILD_ARGS}
    )

# Build and install the pugixml library
ExternalProject_Add(PUGIXML
    PREFIX build
    URL ${CMAKE_CURRENT_SOURCE_DIR}/pugixml/pugixml-1.9.tar.gz
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS})

# Build and install the miniz library
ExternalProject_Add(MINIZ
    PREFIX build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/miniz
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS})

    # Build and install the cvode library
if(NOT DEFINED REAKTORO_PATH_TO_CVODE)
    ExternalProject_Add(CVODE
        PREFIX build
        URL ${CMAKE_CURRENT_SOURCE_DIR}/cvode/cvode-2.8.2.tar.gz
        CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS}
            -DEXAMPLES_ENABLE=OFF
            -DEXAMPLES_INSTALL=OFF)
endif()

# Build and install the Phreeqc library
ExternalProject_Add(PHREEQC
    PREFIX build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/phreeqc
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS})

# Build and install the Gems library
ExternalProject_Add(GEMS
    PREFIX build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gems
    CMAKE_ARGS ${REAKTORO_THIRDPARTY_COMMON_INSTALL_ARGS})
