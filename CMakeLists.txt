cmake_minimum_required (VERSION 3.14)

project (parvmec C CXX Fortran)
set (CMAKE_CXX_STANDARD 11)

add_library (vmec)
add_dependencies (vmec stell)

target_link_libraries (vmec stell)
target_include_directories (vmec PUBLIC $<TARGET_PROPERTY:stell,BINARY_DIR>)

add_executable (xvmec)
add_dependencies (xvmec vmec)
target_link_libraries (xvmec vmec)

add_subdirectory (Sources)

################################################################################
#  Testing                                                                     #
################################################################################

#  Build test utilities.
add_executable (xwout_diff)
add_dependencies (xwout_diff stell)

target_link_libraries (xwout_diff stell)

add_subdirectory (Testing)
