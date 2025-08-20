# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(Gem)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()


# ---Build executable------------------------------------------------------------
add_executable(gem_hole gem_hole.C)
target_link_libraries(gem_hole Garfield::Garfield ROOT::Tree ROOT::Physics)

add_executable(sim main.C)
target_link_libraries(sim Garfield::Garfield ROOT::Tree ROOT::Physics)

add_executable(quick quick.C)
target_link_libraries(quick Garfield::Garfield ROOT::Tree ROOT::Physics)

# add_executable(shower shower.C)
# target_link_libraries(shower Garfield::Garfield ROOT::Tree ROOT::Physics)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file ELIST.lis MPLIST.lis NLIST.lis PRNSOL.lis field.txt dielectrics.dat mesh.mphtxt field3.txt mesh3.mphtxt gas.gas ar80co220.gas)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()

