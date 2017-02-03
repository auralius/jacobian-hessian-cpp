cmake_minimum_required(VERSION 2.6)
project(jacobian-hessian-cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ../bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ../bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CMAKE_BUILD_TYPE Debug)

add_library (jacobian-hessian-cpp ./src/fx.cpp)

target_include_directories (jacobian-hessian-cpp PUBLIC ./src)

add_executable(demo1 ./src/main1.cpp)

install(TARGETS demo1 RUNTIME DESTINATION ../bin)

target_link_libraries (demo1 LINK_PUBLIC jacobian-hessian-cpp armadillo)   