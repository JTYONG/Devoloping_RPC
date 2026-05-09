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
add_executable(AvalMC_addIons AvalMC_addIons.C)
target_link_libraries(AvalMC_addIons Garfield::Garfield OpenMP::OpenMP_CXX)

# ---Build executable------------------------------------------------------------
#add_executable(STDdebug STDdebug.C)
#target_link_libraries(STDdebug Garfield::Garfield OpenMP::OpenMP_CXX)


# ---Copy all files locally to the build directory-------------------------------
foreach(_file c2h2f4_ic4h10_sf6.gas TimingRPCGas.gas)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()

#file(COPY neBEMOut DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/neBEMOut)
	file(COPY neBEMOut DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  message(STATUS "Copied neBEMOut to build directory")
else()
  message(STATUS "neBEMOut not found, skipping copy")
endif()

if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/neBEMInp)
	file(COPY neBEMInp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  message(STATUS "Copied neBEMInp to build directory")
else()
  message(STATUS "neBEMInp not found, skipping copy")
endif()
