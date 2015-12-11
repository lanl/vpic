# turn on -O3
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# turn on avx
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
#    -mavx -qopt-report -qopt-report-phase=vec,loop")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
#    -mavx -qopt-report -qopt-report-phase=vec,loop")

# generate assembler
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -save-temps -fsource-asm")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -save-temps -fsource-asm")


# hacks for using allinea map on trinitite. load the allinea map module first.
#execute_process(COMMAND
#                make-profiler-libraries --platform=cray --lib-type=shared
#                RESULT_VARIABLE HAD_ERROR)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
#    -dynamic -L${CMAKE_CURRENT_BINARY_DIR} -lmap-sampler-pmpi -lmap-sampler \
#    -Wl,--eh-frame-hdr")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
#    -dynamic -L${CMAKE_CURRENT_BINARY_DIR} -lmap-sampler-pmpi -lmap-sampler \
#    -Wl,--eh-frame-hdr")
#set(VPIC_CXX_FLAGS "${VPIC_CXX_FLAGS} \
#    -dynamic -L${CMAKE_CURRENT_BINARY_DIR} -lmap-sampler-pmpi -lmap-sampler \
#    -Wl,--eh-frame-hdr")
