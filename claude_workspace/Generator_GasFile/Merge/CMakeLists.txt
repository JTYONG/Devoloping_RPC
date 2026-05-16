# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(Merge)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(MergeGasFile MergeGasFile.cpp)
target_link_libraries(MergeGasFile Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file TIFRH_40000_40500.gas TIFRH_40500_41000.gas TIFRH_41000_41500.gas TIFRH_41500_42000.gas TIFRH_42000_42500.gas TIFRH_42500_43000.gas)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()
