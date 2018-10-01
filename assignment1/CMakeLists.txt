cmake_minimum_required(VERSION 3.2)
project(CMSC502Project1)

# find header & source
file(GLOB_RECURSE SOURCE "main.cpp")

add_executable(${PROJECT_NAME}
        ${SOURCE}
        ${HEADER}
        )


source_group("Header include" FILES ${HEADER})
source_group("Source src"     FILES ${SOURCE})
