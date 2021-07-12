# Recursively collect all header files from the current directory
file(GLOB_RECURSE HEADER_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.hpp *.hxx)

# Recursively collect all source files from the current directory
file(GLOB_RECURSE SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# The name of the shared and static libraries
set(REAKTORO_SHARED_LIB ${PROJECT_NAME}${SUFFIX_SHARED_LIBS})
set(REAKTORO_STATIC_LIB ${PROJECT_NAME}${SUFFIX_STATIC_LIBS})

# Set the list of names of the third-party targets and libraries
set(THIRDPARTY_TARGETS PUGIXML MINIZ PHREEQC GEMS)
set(THIRDPARTY_LIBS pugixml miniz phreeqc gems nlohmann_json::nlohmann_json)

if(DEFINED REAKTORO_PATH_TO_CVODE)
    list(APPEND THIRDPARTY_LIBS ${SUNDIALS_CVODE} ${SUNDIALS_NVECSERIAL})
else()
    list(APPEND THIRDPARTY_TARGETS CVODE)
    list(APPEND THIRDPARTY_LIBS sundials_cvode sundials_nvecserial)
endif()

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

if(REAKTORO_USE_OPENLIBM)
    configure_target_to_use_openlibm(Reaktoro)
    target_compile_definitions(Reaktoro PUBLIC REAKTORO_USE_OPENLIBM=1)
endif()

if(REAKTORO_PATH_TO_CVODE)
    target_include_directories(Reaktoro PRIVATE ${REAKTORO_PATH_TO_CVODE}/include)
    message(${REAKTORO_PATH_TO_CVODE}/include)
endif()

# Link Reaktoro library against ThermoFun if found
if(ThermoFun_FOUND)
    target_link_libraries(Reaktoro PUBLIC ThermoFun::ThermoFun)
endif()

# Set compilation features to be propagated to dependent codes.
target_compile_features(Reaktoro PUBLIC cxx_std_17)

# Set some target properties
set_target_properties(Reaktoro PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CXX_EXTENSIONS OFF)

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
