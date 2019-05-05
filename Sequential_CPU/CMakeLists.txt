
cmake_minimum_required(VERSION 3.0)
project(barnesHutSimulation)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pthread")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB SRC_FILES *.cpp *.hpp)

add_executable(main.exe ${SRC_FILES})
target_link_libraries(main.exe)