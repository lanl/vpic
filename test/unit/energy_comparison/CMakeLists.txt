LIST(APPEND TESTS weibel_driver 3d_test)

list(APPEND gold_file "${CMAKE_CURRENT_SOURCE_DIR}/energies_gold")

foreach (test ${TESTS})
    add_executable(${test} ./${test}.cc)
    #add_definitions(-DGOLD_ENERGY_FILE=${gold_file})
    target_compile_definitions(${test} PRIVATE GOLD_ENERGY_FILE=${gold_file}.${test})
    target_link_libraries(${test} vpic)
    add_test(NAME ${test} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ./${test})
endforeach(test)

# TODO: this could be collapsed using if defs into the original weibel
set(test weibel_driver_non_vaccum)
add_executable(${test} ./weibel_driver.cc)
target_compile_definitions(${test} PRIVATE ADD_EXTRA_MATERIAL=1 GOLD_ENERGY_FILE=${gold_file}.weibel_driver)
target_link_libraries(${test} vpic)
add_test(NAME ${test} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ./${test})

set(test 3d_test_non_vaccum)
add_executable(${test} ./3d_test.cc)
target_compile_definitions(${test} PRIVATE ADD_EXTRA_MATERIAL=1 GOLD_ENERGY_FILE=${gold_file}.3d_test)
target_link_libraries(${test} vpic)
add_test(NAME ${test} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ./${test})

set(test 3d_test_non_vaccum_threaded)
add_executable(${test} ./3d_test.cc)
target_compile_definitions(${test} PRIVATE ADD_EXTRA_MATERIAL=1 GOLD_ENERGY_FILE=${gold_file}.3d_test_threaded)
target_link_libraries(${test} vpic)
add_test(NAME ${test} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ./${test} --tpp 8)
