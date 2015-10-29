# build the executable
execute_process(COMMAND
                ${BINDIR}/bin/vpic
                ${SRCDIR}/test/integrated/${NAME}.deck
		RESULT_VARIABLE HAD_ERROR)
if(HAD_ERROR)
  message(FATAL_ERROR "Test ${NAME} failed to build.")
endif(HAD_ERROR)

# run the executable
execute_process(COMMAND
                ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI_NUM_RANKS}
                ${NAME}.deck.${CMAKE_SYSTEM_NAME} ${ARGS}
                RESULT_VARIABLE HAD_ERROR)
if(HAD_ERROR)
  message(FATAL_ERROR "Test ${NAME} failed.")
endif(HAD_ERROR)

# compare the results
# ??
# this file could easily be extended to compare generated output against gold
# standard output.
