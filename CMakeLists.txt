project(ihpc)

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpic++)

set(CMAKE_CXX_RELEASE_FLAGS "${CMAKE_CXX_RELEASE_FLAGS} -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(ihpc_seq SEQ_Poisson.c)
add_executable(ihpc MPI_Poisson.cpp)
add_executable(ihpc2 MPI_Poisson2.cpp)
add_executable(ihpc3 MPI_Fempois.cpp)
add_executable(GridDist GridDist.c)

target_link_libraries(GridDist m)
