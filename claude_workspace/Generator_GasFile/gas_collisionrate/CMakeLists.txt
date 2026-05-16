# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(generate_gas)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(generate_gas generate_gas.cc)
target_link_libraries(generate_gas Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file )
  configure_file(${_file} ${_file} COPYONLY)
endforeach()
