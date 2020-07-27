cmake_minimum_required(VERSION 3.10)

function(configure_mpi_target target)
  target_compile_definitions(${target} PRIVATE ${MPI_Fortran_COMPILE_DEFINITIONS})
  target_compile_options(${target} PRIVATE ${MPI_Fortran_COMPILE_OPTIONS})
  target_include_directories(${target} PRIVATE ${MPI_Fortran_INCLUDE_DIRS})
  target_link_libraries(${target} PRIVATE ${MPI_Fortran_LINK_FLAGS})
  target_link_libraries(${target} PRIVATE ${MPI_Fortran_LIBRARIES})
endfunction()

function(configure_hdf5_target target)

  cmake_policy(SET CMP0079 NEW)

  target_link_libraries(${target} PRIVATE ${HDF5_Fortran_LIBRARIES} ${HDF5_LIBRARIES})

  if(HDF5_Fortran_INCLUDE_DIRS)
    target_include_directories(${target} PRIVATE ${HDF5_Fortran_INCLUDE_DIRS})
  elseif(HDF5_Fortran_INCLUDE_DIR)
    target_include_directories(${target} PRIVATE ${HDF5_Fortran_INCLUDE_DIR})
  else()
    message(FATAL_ERROR "No HDF5 Fortran include directory found.")
  endif()

endfunction()
