cmake_minimum_required(VERSION 3.2)
project(CMSC502Project1)

# find header & source
file(GLOB_RECURSE SOURCE "sequential.cpp")
file(GLOB_RECURSE SOURCE "threaded.cpp")
file(GLOB_RECURSE SOURCE "mpi.cpp")

add_executable(${PROJECT_NAME}Sequential
        sequential.cpp
        ${HEADER}
        )

add_executable(${PROJECT_NAME}Threaded
        threaded.cpp
        ${HEADER}
        )

add_executable(${PROJECT_NAME}MPI
        mpi.cpp
        ${HEADER}
        )


find_package(MPI REQUIRED)

set(CMAKE_C_COMPILER /usr/local/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/local/bin/mpic++)

include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME}MPI ${MPI_LIBRARIES})


source_group("Header include" FILES ${HEADER})
source_group("Source src"     FILES ${SOURCE})
