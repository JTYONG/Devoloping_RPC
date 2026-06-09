# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(RPC)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(ParallelPlate ParallelPlate.C)
target_link_libraries(ParallelPlate Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file c2h2f4_95-5_ic4h10_4-2_sf6_0-3.gas)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()
