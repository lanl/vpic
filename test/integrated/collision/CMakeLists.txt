set(test test_collision)
build_a_vpic(${test} ${CMAKE_CURRENT_SOURCE_DIR}/${test}.deck)
add_test(NAME ${test} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ./${test})

include(FindPythonInterp)

if (${PYTHONINTERP_FOUND})
    # for now, only check the correctness if we can find python
    add_test(
        NAME ${test}_script
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_col_jx.py
    )
    set_tests_properties(${test}_script PROPERTIES DEPENDS ${test})
else()
    message("Python not found. Collision test will not check answers against theory")
endif()
