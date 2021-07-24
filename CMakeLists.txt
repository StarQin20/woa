cmake_minimum_required(VERSION 2.8)
project(woa)
add_compile_options(-std=c++11)

include_directories(${PROJECT_SOURCE_DIR}/inc/)
add_executable(woa src/woa.cpp)
