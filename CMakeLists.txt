project(ihpc)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpic++)

add_executable(ihpc main.cpp)
