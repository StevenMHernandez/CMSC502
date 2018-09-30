cmake_minimum_required(VERSION 3.2)
project(DosomethinghereWhyNot)

# includes
#include_directories(
#        include
#)

# find header & source
#file(GLOB_RECURSE HEADER "include/*.h")
file(GLOB_RECURSE SOURCE "main.cpp")

add_executable(${PROJECT_NAME}
        ${SOURCE}
        ${HEADER}
        )


source_group("Header include" FILES ${HEADER})
source_group("Source src"     FILES ${SOURCE})

# link libraries
#target_link_libraries(${PROJECT_NAME} ${SOME_LIB})