
if(FPRETTIFY_EXECUTABLE)

  set(FPRETTIFY_CMD ${DO_FPRETTIFY})
  list(APPEND FPRETTIFY_CMD "${ALL_SOURCES}")
  string(REPLACE ";" "\;" FPRETTIFY_CMD "${FPRETTIFY_CMD}")
  add_test(NAME test_fprettify
    COMMAND ${CMAKE_COMMAND}
    -DCMD=${FPRETTIFY_CMD}
    -DFAIL_PATTERN=.
    -P ${CMAKE_CURRENT_SOURCE_DIR}/test_and_check_regex.cmake)
endif()

add_library(test_types dummy_model_forward_operator.F90 random_observations.F90 dummy_model_interfaces.F90)
target_link_libraries(test_types PRIVATE lightda)
configure_mpi_target(test_types)

add_executable(dummy_assimilator_tests dummy_assimilator_tests.F90)
target_link_libraries(dummy_assimilator_tests PRIVATE lightda)

add_library(model_interface_tests model_interface_tests.f90)
target_link_libraries(model_interface_tests PRIVATE lightda)

add_executable(dummy_model_interface_tests dummy_model_interface_tests.f90)
target_link_libraries(dummy_model_interface_tests PRIVATE model_interface_tests system_mpi test_types)

add_test(test_dummy_model_interface ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ./dummy_model_interface_tests)

add_test(test_dummy_assimilator dummy_assimilator_tests)

add_executable(batch_manager_tests batch_manager_tests.F90)
target_link_libraries(batch_manager_tests PRIVATE test_types lightda)

add_test(test_batch_manager ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ./batch_manager_tests)

add_executable(assimilation_manager_tests assimilation_manager_tests.F90)
target_link_libraries(assimilation_manager_tests PRIVATE test_types lightda)

add_test(test_assimilation_manager ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ./assimilation_manager_tests)

add_executable(util_tests util_tests.f90)
target_link_libraries(util_tests lightda_util)

add_test(test_util ./util_tests)

add_library(darray_tests darray_tests.F90)
target_link_libraries(darray_tests lightda_util)

add_executable(darray_test test_darray.F90)
target_link_libraries(darray_test PUBLIC darray_tests)
configure_mpi_target(darray_test)

add_test(test_darray ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ./darray_test)

add_custom_target(all_tests ALL
  DEPENDS darray_test util_tests assimilation_manager_tests batch_manager_tests
  dummy_model_interface_tests dummy_assimilator_tests)

add_custom_command(
     TARGET all_tests
     COMMENT "Run tests"
     POST_BUILD
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
     COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION> --output-on-failures
)
