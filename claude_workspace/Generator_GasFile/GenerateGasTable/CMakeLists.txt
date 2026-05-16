cmake_minimum_required(VERSION 3.16)
project(generate_gas LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Garfield REQUIRED)

add_executable(generate_gas generate_gas.cc)
target_link_libraries(generate_gas PRIVATE Garfield::Garfield)
