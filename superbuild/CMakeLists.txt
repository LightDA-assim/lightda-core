cmake_minimum_required(VERSION 3.10)

include(ExternalProject)

project(lenkf-rsm)

ExternalProject_Add(
  system_mpi
  GIT_REPOSITORY git@gitlab.hpc.mil:john.haiducek.ctr/system_mpi.git
  GIT_TAG main
  CMAKE_CACHE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/system_mpi)

ExternalProject_Add(
  fortran_exceptions
  GIT_REPOSITORY git@gitlab.hpc.mil:john.haiducek.ctr/fortran_exceptions.git
  GIT_TAG main
  CMAKE_CACHE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/fortran_exceptions)

ExternalProject_Add(lightda
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.."
  CMAKE_CACHE_ARGS
    -Dsystem_mpi_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/system_mpi/lib/cmake/system_mpi
    -Dfortran_exceptions_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/fortran_exceptions/lib/cmake/fortran_exceptions
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/lightda
  DEPENDS system_mpi fortran_exceptions)
