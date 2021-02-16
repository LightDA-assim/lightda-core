cmake_minimum_required(VERSION 3.10)

project(lightda)

enable_language(Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/")

set(Fortran_MODULE_DIRECTORY include/fortran_modules/${CMAKE_Fortran_COMPILER_ID}/${CMAKE_Fortran_COMPILER_VERSION})
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${Fortran_MODULE_DIRECTORY})

include(lightda_functions)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 9)
    # Compiler doesn't support overriding finalizers
  else()
    # Compiler supports overriding finalizers
    add_definitions(-DOVERRIDABLE_FINALIZERS)
  endif()
else()
  # Compiler supports overriding finalizers
  add_definitions(-DOVERRIDABLE_FINALIZERS)
endif()

find_package(system_mpi REQUIRED)
find_package(fortran_exceptions REQUIRED)

add_subdirectory(src)

get_target_property(lightda_INTERFACE_INCLUDE_DIRECTORIES lightda INTERFACE_INCLUDE_DIRECTORIES)

message("lightda_INTERFACE_INCLUDE_DIRECTORIES=${lightda_INTERFACE_INCLUDE_DIRECTORIES}")

get_target_property(lightda_INCLUDE_DIRECTORIES lightda INCLUDE_DIRECTORIES)

message("lightda_INCLUDE_DIRECTORIES=${lightda_INCLUDE_DIRECTORIES}")

configure_file(animate.py animate.py COPYONLY)
configure_file(animate_mpi.py animate_mpi.py COPYONLY)
configure_file(advect1d.py advect1d.py COPYONLY)
configure_file(enkf.py enkf.py COPYONLY)
configure_file(libsuffix.py libsuffix.py COPYONLY)

include(CMakePackageConfigHelpers)

# Add all targets to the build-tree export set
export(TARGETS lightda lightda_util lightda_random
  FILE "${PROJECT_BINARY_DIR}/lightdaTargets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE lightda)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/lightdaConfig.cmake"
  INSTALL_DESTINATION lib/cmake/lightda
  )

install(FILES
          "${CMAKE_CURRENT_BINARY_DIR}/lightdaConfig.cmake"
        DESTINATION lib/cmake/lightda
)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${Fortran_MODULE_DIRECTORY}/..)

install(EXPORT lightda
  FILE lightdaTargets.cmake
  DESTINATION lib/cmake/lightda
)

# -------------------
# Build documentation
# -------------------

set ( SKIP_DOC_GEN FALSE CACHE BOOL
  "Disable building the API documentation with FORD" )
if ( NOT SKIP_DOC_GEN )
  find_program(FORD ford)

  if(FORD)

    set(DOC_ROOT "${PROJECT_BINARY_DIR}/doc")
    set(DOC_DIR "${DOC_ROOT}/lightda")
    set(PROJ_DIR "${PROJECT_SOURCE_DIR}")
    set(FORD_PROJECT_FILE "${CMAKE_CURRENT_BINARY_DIR}/lightda.md")
    set(FORD_PROJECT_EXCLUDE_DIRS "${MPI_Fortran_INCLUDE_DIRS};${PROJECT_BINARY_DIR};${PROJ_DIR}/build;${PROJ_DIR}/tests;${PROJ_DIR}/sphinx")
    set(FORD_PROJECT_EXCLUDE_FILES "random.f90;random_integer.f90;tree.f90")
    set(MACRO_FLAG "")
    string(REPLACE ";" "\ninclude:" FORD_PROJECT_INCLUDES "${PROJ_DIR}/src")
    string(REPLACE ";" "\nexclude_dir:" FORD_PROJECT_EXCLUDE_DIRS "${FORD_PROJECT_EXCLUDE_DIRS}")
    string(REPLACE ";" "\nexclude:" FORD_PROJECT_EXCLUDE_FILES "${FORD_PROJECT_EXCLUDE_FILES}")

    configure_file(lightda.md.in lightda.md)

    # Pick the preprocessor to use based on the Fortran compiler
    if ( "${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel" )
      set ( FPP "fpp\n" )
    else ()
      set ( FPP "gfortran -E\n" ) # default to gfortran -E for gfortran and unsupported compilers
    endif ()
    file ( WRITE "${PROJECT_BINARY_DIR}/.PREPROCESSOR" "${FPP}" )

    # Dynamically generate the FORD outputs list
    message ( STATUS "Dynamically computing FORD output information..." )
    if ( NOT (DEFINED FORD_OUTPUTS_CACHED) )
      message ( STATUS "Running FORD to dynamically compute documentation outputs, this could take a while..." )
      execute_process ( COMMAND ${CMAKE_COMMAND} -E remove_directory ${DOC_DIR}
	COMMAND ${CMAKE_COMMAND} -E make_directory ${DOC_DIR}
	COMMAND "${FORD}" --debug -q ${MACRO_FLAG} -d "${PROJ_DIR}" -o "${DOC_DIR}" -p "${PAGE_DIR}" "${FORD_PROJECT_FILE}" OUTPUT_QUIET )
    else ()
      message ( STATUS "Re-using cached FORD outputs, rather than regenerating them" )
    endif()

    # Compile a list of FORD output files
    file ( GLOB_RECURSE FORD_OUTPUTS
      "${DOC_DIR}/*.*" )
    file ( GLOB_RECURSE FORD_CLEAN_OUTPUTS
      "${DOC_DIR}/*.*" )

    # Cache FORD outputs
    if ( (DEFINED FORD_OUTPUTS) AND ( NOT ( "${FORD_OUTPUTS}" STREQUAL "" ) ) )
      message ( STATUS "Caching FORD outputs" )
      set ( FORD_OUTPUTS_CACHED "${FORD_OUTPUTS}"
	CACHE STRING "variable containing FORD outputs to prevent rebuilding FORD docs" FORCE )
    endif ()
    message ( STATUS "Done dynamically computing FORD outputs." )

    foreach ( DOC_SRC_FILE ${PROJECT_BINARY_DIR}/lightda.md ${SYSTEM_MPI_SOURCES} ${UTIL_SOURCES} ${ASSIM_SOURCES} ${ASSIM_PARALLEL_SOURCES} ${ADVECT1D_ASSIMILATE_INTERFACES_SOURCES} ${ADVECT1D_ASSIMILATE_SOURCES} ${ADVECT1D_SOURCES} ${PAGES_FILES} )
      list ( APPEND FORD_DEPENDS "${DOC_SRC_FILE}" )
    endforeach ()

    # Re-build FORD output if needed
    add_custom_command ( OUTPUT ${FORD_OUTPUTS_CACHED}
      COMMAND "${FORD}" --debug ${MACRO_FLAG} -d "${PROJ_DIR}" -o "${DOC_DIR}" "${FORD_PROJECT_FILE}"
      MAIN_DEPENDENCY "${FORD_PROJECT_FILE}"
      DEPENDS ${FORD_DEPENDS}
      COMMENT "Building HTML documentation for ${PROJECT_NAME} using FORD" )

    add_custom_target ( documentation ALL
      DEPENDS ${FORD_OUTPUTS_CACHED} )

    # Install documentation
    set ( INSTALL_API_DOCUMENTATION TRUE
      CACHE BOOL "Install FORD generated documentation?" )
    if ( INSTALL_API_DOCUMENTATION )
      install ( DIRECTORY "${DOC_ROOT}/" TYPE DOC )
    else ()

      # FORD Not found
      message ( WARNING
	"FORD not found. Please set the CMake cache variable FORD to point to the installed FORD executable, and reconfigure or disable building the documentation. FORD can be installed from PYPI with `sudo pip install FORD` or from <https://github.com/cmacmackin/ford> If you do not wish to install FORD and build the JSON-Fortran documentation, then please set the CMake cache variable SKIP_DOC_GEN to TRUE." )

    endif ()
  endif ()

endif()

find_program(FPRETTIFY_EXECUTABLE NAMES "fprettify" DOC "Path to fprettify")

if(NOT FPRETTIFY_EXECUTABLE)
  message(STATUS "fprettify not found.")
else()
  set(DO_FPRETTIFY "${FPRETTIFY_EXECUTABLE}" "--diff" "-i" "2" "-l" "81")

  list(TRANSFORM ALL_SOURCES PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/)

  add_custom_target(
    fprettify
    COMMAND ${DO_FPRETTIFY}
    ${ALL_SOURCES}
    )

endif()

enable_testing()
add_subdirectory(tests)
