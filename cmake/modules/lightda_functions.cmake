cmake_minimum_required(VERSION 3.10)

enable_language(Fortran)

find_package(HDF5 COMPONENTS Fortran REQUIRED)

function(configure_hdf5_target target)

  target_link_libraries(${target} PRIVATE ${HDF5_Fortran_LIBRARIES} ${HDF5_LIBRARIES})

  if(HDF5_Fortran_INCLUDE_DIRS)
    target_include_directories(${target} PRIVATE ${HDF5_Fortran_INCLUDE_DIRS})
  elseif(HDF5_Fortran_INCLUDE_DIR)
    target_include_directories(${target} PRIVATE ${HDF5_Fortran_INCLUDE_DIR})
  else()
    message(FATAL_ERROR "No HDF5 Fortran include directory found.")
  endif()

endfunction()
