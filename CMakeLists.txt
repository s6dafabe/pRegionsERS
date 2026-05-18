cmake_minimum_required(VERSION 3.5)

project(aggregation CXX)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(default_build_type "Release")

#cmake module for finding gurobi
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

#add gurobi 
find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})

# Source Files
set(SOURCE_FILES )
add_subdirectory(src)

add_executable(aggregation src/main.cpp ${SOURCE_FILES})

#gurobi libraries
find_package(OpenMP)
target_link_libraries(${CMAKE_PROJECT_NAME}  PUBLIC OpenMP::OpenMP_CXX optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY} ${GUROBI_LIBRARY})


