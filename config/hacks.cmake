# turn on -O3
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# turn on avx
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -mavx -qopt-report -qopt-report-phase=vec,loop")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -mavx -qopt-report -qopt-report-phase=vec,loop")

# generate assembler
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -save-temps -fsource-asm")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -save-temps -fsource-asm")
