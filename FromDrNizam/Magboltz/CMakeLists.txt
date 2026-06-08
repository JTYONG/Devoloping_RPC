cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(make_gas)

find_package(Garfield REQUIRED)

add_executable(make_gas make_gas.C)
target_link_libraries(make_gas Garfield::Garfield)
