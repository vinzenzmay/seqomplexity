cmake_minimum_required(VERSION 3.12)

project(sequence_complexity)

set(CMAKE_CXX_STANDARD 11)

# Add your source files
set(SOURCE_FILES
    src/sequence_complexity.cpp
)

add_executable(sequence_complexity ${SOURCE_FILES})
