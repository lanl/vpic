find_package( PythonInterp )

# we chose to pass full paths rather than copying files around
set(GRID_HEATING_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/gridHeatingTestElec.py")

set(MPI_NUMPROC 1)

set(test "gridHeatingTestElec")
add_executable(${test} ./${test}.cxx)
target_link_libraries(${test} vpic)
add_test(NAME ${test} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI_NUMPROC} ${MPIEXEC_PREFLAGS} ./${test})

set(PY_TEST_NAME "grid_heating_check")
set(TEST_GROUP_LABEL "${dim}d_heating")

add_test(NAME ${PY_TEST_NAME} COMMAND ${PYTHON_EXECUTABLE} ${GRID_HEATING_SCRIPT} ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(${PY_TEST_NAME} PROPERTIES DEPENDS ${test})
set_property(TEST ${test} PROPERTY FIXTURES_SETUP ${TEST_GROUP_LABEL})
set_property(TEST ${PY_TEST_NAME} PROPERTY FIXTURES_REQUIRED ${TEST_GROUP_LABEL})
