find_package( PythonInterp )

# we chose to pass full paths rather than copying files around
set(GRID_HEATING_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/gridHeatingTest.py")

LIST(APPEND dims_to_test 1 2)

foreach(dim IN LISTS dims_to_test)
    set(test "gridHeatingTest${dim}D")
    add_executable(${test} ./${test}.cxx)
    target_link_libraries(${test} vpic)
    add_test(NAME ${test} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ./${test})

    set(PY_TEST_NAME "grid_heating_check_${dim}d")
    set(TEST_GROUP_LABEL "${dim}d_heating")

    add_test(NAME ${PY_TEST_NAME} COMMAND ${PYTHON_EXECUTABLE} ${GRID_HEATING_SCRIPT} ${CMAKE_CURRENT_BINARY_DIR} ${dim})

    set_tests_properties(${PY_TEST_NAME} PROPERTIES DEPENDS ${test})
    set_property(TEST ${test} PROPERTY FIXTURES_SETUP ${TEST_GROUP_LABEL})
    set_property(TEST ${PY_TEST_NAME} PROPERTY FIXTURES_REQUIRED ${TEST_GROUP_LABEL})
endforeach(dim)
