cmake_minimum_required(VERSION 3.10)

include(ExternalProject)

project(lenkf-rsm)

ExternalProject_Add(
  system_mpi
  GIT_REPOSITORY /Users/jhaiduce/Development/repositories/system_mpi.git
  CMAKE_CACHE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/system_mpi)

ExternalProject_Add(lenkf-rsm
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.."
  CMAKE_CACHE_ARGS
    -Dsystem_mpi_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/system_mpi/lib/cmake/system_mpi
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/lenkf-rsm
  DEPENDS system_mpi)
