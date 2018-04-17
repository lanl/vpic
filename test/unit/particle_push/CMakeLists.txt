add_executable(array_syntax ./array_syntax.cc)
target_link_libraries(array_syntax vpic)
add_test(NAME array_syntax COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ./array_syntax)
