# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(RPC)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(AvalMC_Nebem AvalMC_Nebem.C)
target_link_libraries(AvalMC_Nebem Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file c2h2f4_ic4h10_sf6.gas TimingRPCGas.gas)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()
