set(target uncenter)
add_executable(${target} ./uncenter.cpp)
target_link_libraries(${target} vpic)
add_test(NAME ${target} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ./${target})
