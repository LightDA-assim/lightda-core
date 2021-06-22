LightDA is a lightweight, extensible data assimilation library. It is designed with the following goals in mind:

- The data assimilation process be model-agnostic, and adding a new model should be as simple as possible
- Assimilation can be done in parallel if needed
- Parallel models are supported without making any assumptions about how the model's parallelism is implemented

This repository contains the LightDA core library. It provides the infrastructure for parallel data assimilation and interface definitions for integrating models, filters, observation types, forward operators, and localization schemes.

## Requirements

LightDA requires CMake 3.0.2 or later, a Fortran 2008 compiler (gcc 4.9.4 and later and ifort 2021 have been tested) and an MPI library. Building the LightDA documentation additionally requires Python 3.5 or later. The LightDA examples also require HDF5.

## Quick start

The easiest way to get started with LightDA is to build the [lightda-examples package](https://github.com/LightDA-assim/lightda-examples), which provides usage examples for LightDA.

Use the following commands to build lightda-examples and its dependencies (including the LightDA core library):
```bash
git clone https://github.com/LightDA-assim/lightda-examples.git
cd lightda-examples
mkdir build
cmake ../superbuild
make
```

## Compiling

The LightDA core library can be built by itself (without the examples) as follows:

```bash
git clone https://github.com/LightDA-assim/lightda.git
cd lightda
mkdir build
cd build
cmake ../superbuild
make
```

The above commands build [system_mpi](https://github.com/LightDA-assim/system_mpi) and [fortran_exceptions](https://github.com/LightDA-assim/fortran_exceptions) as part of the LightDA build. If you prefer to build these separately, this can be accomplished with by replacing the above ```cmake ../superbuild``` with ```cmake ..```. This requires building and installing system_mpi and fortran_exceptions first.

## Extending

LightDA makes extensive use of Fortran 2003 [abstract types](https://gist.github.com/n-s-k/de4af7ce6cc8f2c85e4b33cedb51fd88#file-oop_f2003_part_2-md) and [type-bound procedures](https://gist.github.com/n-s-k/de4af7ce6cc8f2c85e4b33cedb51fd88#file-oop_f2003_part_2-md) to define interfaces for models, observation types, forward operators, and assimilation algorithms. For those new to object-oriented programming in Fortran, Mark Leair's [tutorial on the subject](https://gist.github.com/n-s-k/522f2669979ed6d0582b8e80cf6c95fd) may prove helpful.

The abstract types and type-bound procedures defined in LightDA provide a standard interface between the assimilation system and models, enabling new models to be added without requiring changes to LightDA itself.

LightDA uses the `exceptions` module of the [fortran_exceptions](file://${fortran_exceptions_DOCDIR}/index.html) library for error handling. Most procedures defined in LightDA accept an optional `status` argument of type `error_status`. This optional argument provides a means to send error information back to the calling procedure, where the error can be handled if possible. If the `status` argument is not given in the call, any error will result in program termination. Details on how to use this facility are provided in the [fortran_exceptions documentation](file://${fortran_exceptions_DOCDIR}/index.html).
