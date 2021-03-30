find_package(LAPACK REQUIRED)

set(UTIL_SOURCES util.F90 tree.f90 list.f90 darray.F90)
add_library(lightda_util SHARED ${UTIL_SOURCES})
target_include_directories(lightda_util INTERFACE
  $<INSTALL_INTERFACE:${Fortran_MODULE_DIRECTORY}>)
target_link_libraries(lightda_util PUBLIC fortran_exceptions)
target_link_libraries(lightda_util PUBLIC system_mpi)
install(TARGETS lightda_util EXPORT lightda
  LIBRARY DESTINATION lib ${CMAKE_INSTALL_LIBDIR})

set(RANDOM_SOURCES random.f90 random_integer.f90)
add_library(lightda_random SHARED ${RANDOM_SOURCES})
install(TARGETS lightda_random EXPORT lightda
  LIBRARY DESTINATION lib ${CMAKE_INSTALL_LIBDIR})

set(LIGHTDA_SOURCES dummy_assimilator.f90 base_assimilation_manager.F90 assimilation_filter.F90 assimilation_batch_manager.F90 assimilation_model_interface.F90 observations.F90 localization.f90 observation_manager.F90 forward_operator.F90 assimilation_manager.F90)
add_library(lightda SHARED ${LIGHTDA_SOURCES})
target_link_libraries(lightda PUBLIC lightda_random lightda_util)
target_include_directories(lightda INTERFACE
  $<INSTALL_INTERFACE:${Fortran_MODULE_DIRECTORY}>)
install(TARGETS lightda EXPORT lightda
  LIBRARY DESTINATION lib ${CMAKE_INSTALL_LIBDIR})

set(ALL_SOURCES ${UTIL_SOURCES};${RANDOM_SOURCES};${LIGHTDA_SOURCES})

# Prepend path to sources and propagate to parent scope
list_transform_prepend(ALL_SOURCES "src/")
set(ALL_SOURCES ${ALL_SOURCES} PARENT_SCOPE)
