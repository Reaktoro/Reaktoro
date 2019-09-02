# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp *.hxx)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# The name of the shared and static libraries
set(REAKTORO_SHARED_LIB ${PROJECT_NAME}${SUFFIX_SHARED_LIBS})
set(REAKTORO_STATIC_LIB ${PROJECT_NAME}${SUFFIX_STATIC_LIBS})

# Set the list of names of the third-party targets and libraries
set(THIRDPARTY_TARGETS PUGIXML MINIZ CVODE PHREEQC GEMS)
set(THIRDPARTY_LIBS pugixml miniz sundials_cvode sundials_nvecserial phreeqc gems)

# Enable automatic creation of a module definition (.def) file for a SHARED library on Windows.
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# Add the library path of the external dependencies to link directories
# Note: cmake v3.13 is needed to use target_link_directories command!
link_directories(${REAKTORO_THIRDPARTY_LIBRARY_PATH})

# Build a Reaktoro C++ library
add_library(Reaktoro ${HEADER_FILES} ${SOURCE_FILES})

# Add an alias Reaktoro::Reaktoro to the target library Reaktoro
add_library(Reaktoro::Reaktoro ALIAS Reaktoro)

# Specify the external dependencies of Reaktoro to ensure proper build sequence
add_dependencies(Reaktoro ${THIRDPARTY_TARGETS})

# Add the include paths to the Reaktoro target
target_include_directories(Reaktoro
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE ${REAKTORO_THIRDPARTY_INCLUDE_PATH})

# Link Reaktoro library against external dependencies
target_link_libraries(Reaktoro
    PRIVATE ${THIRDPARTY_LIBS}
    PUBLIC Boost::boost)
# Set compilation features to be propagated to dependent codes.
target_compile_features(Reaktoro PUBLIC cxx_std_17)

# Install Reaktoro C++ library
install(TARGETS Reaktoro
    EXPORT ReaktoroTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT libraries)

# Install debug symbols
if(MSVC)
    install(
        FILES $<TARGET_PDB_FILE:Reaktoro>
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT libraries
        OPTIONAL
    )
endif()

# Install Reaktoro header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/Reaktoro
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT headers
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx")

# Install deps header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/Reaktoro/deps
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Reaktoro COMPONENT headers)
