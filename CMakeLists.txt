cmake_minimum_required(VERSION 2.8)

project(cs211-lab1)

include_directories(/act/opt/intel/mkl/include)
link_directories(/act/opt/intel/mkl/lib/intel64)
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath=/act/opt/intel/mkl/lib/intel64")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -W -DMKL_ILP64 -lmkl_avx2 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -m64")
add_executable(reg_reuse reg_reuse.c util.c mygemm.c)
add_executable(cache_part3 cache_part3.c util.c mygemm.c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -W -DMKL_ILP64 -lmkl_avx2 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -m64")
add_executable(cache_part4_o0 cache_part4.c util.c mygemm.c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O1 -W -DMKL_ILP64 -lmkl_avx2 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -m64")
add_executable(cache_part4_o1 cache_part4.c util.c mygemm.c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -W -DMKL_ILP64 -lmkl_avx2 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -m64")
add_executable(cache_part4_o2 cache_part4.c util.c mygemm.c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -W -DMKL_ILP64 -lmkl_avx2 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -m64")
add_executable(cache_part4_o3 cache_part4.c util.c mygemm.c)
