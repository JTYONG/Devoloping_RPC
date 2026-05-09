# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(RPC)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
find_package(OpenMP REQUIRED)
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(Debug_10kV Debug_10kV.C)
target_link_libraries(Debug_10kV Garfield::Garfield OpenMP::OpenMP_CXX)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file c2h2f4_ic4h10_sf6.gas TimingRPCGas.gas)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()
