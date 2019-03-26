LIST(APPEND TESTS weibel_driver)

list(APPEND gold_file "${CMAKE_CURRENT_SOURCE_DIR}/energies_gold")
add_definitions(-DGOLD_ENERGY_FILE=${gold_file})

foreach (test ${TESTS})
    add_executable(${test} ./${test}.cc)
    target_link_libraries(${test} vpic)
    add_test(NAME ${test} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ./${test})
endforeach(test)

set(test weibel_driver_non_vaccum)
add_executable(${test} ./weibel_driver.cc)
target_compile_definitions(${test} PRIVATE ADD_EXTRA_MATERIAL=1)
target_link_libraries(${test} vpic)
add_test(NAME ${test} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ./${test})
