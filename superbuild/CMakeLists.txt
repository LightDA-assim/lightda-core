cmake_minimum_required(VERSION 3.10)

include(ExternalProject)

project(lenkf-rsm)

find_package(PythonInterp)

set(
  system_mpi_GIT_URL
  "https://github.com/LightDA-assim/system_mpi.git" CACHE STRING
  "URL of system_mpi git repository")

set(
  fortran_exceptions_GIT_URL
  "https://github.com/LightDA-assim/fortran_exceptions.git" CACHE STRING
  "URL of fortran_exceptions git repository")

ExternalProject_Add(
  system_mpi
  GIT_REPOSITORY ${system_mpi_GIT_URL}
  GIT_TAG main
  CMAKE_CACHE_ARGS
        -DPYTHON_EXECUTABLE:PATH=${PYTHON_EXECUTABLE}
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/system_mpi
	-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE})

ExternalProject_Add(
  fortran_exceptions
  GIT_REPOSITORY ${fortran_exceptions_GIT_URL}
  GIT_TAG main
  CMAKE_CACHE_ARGS
        -DPYTHON_EXECUTABLE:PATH=${PYTHON_EXECUTABLE}
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/fortran_exceptions
	-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE})

ExternalProject_Add(lightda
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.."
  CMAKE_CACHE_ARGS
    -Dsystem_mpi_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/system_mpi/lib/cmake/system_mpi
    -Dfortran_exceptions_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/fortran_exceptions/lib/cmake/fortran_exceptions
    -DPYTHON_EXECUTABLE:PATH=${PYTHON_EXECUTABLE}
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/lightda
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  DEPENDS system_mpi fortran_exceptions)
