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


source_group("Header include" FILES ${HEADER})
source_group("Source src"     FILES ${SOURCE})
