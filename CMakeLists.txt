cmake_minimum_required (VERSION 3.14)

project (parvmec C CXX Fortran)

add_library (vmec)

target_link_libraries (vmec PUBLIC stell)
target_include_directories (vmec

                            PUBLIC

                            $<TARGET_PROPERTY:vmec,BINARY_DIR>
)

add_executable (xvmec)
target_link_libraries (xvmec PUBLIC vmec)

add_subdirectory (Sources)

install (TARGETS xvmec)

################################################################################
#  Testing                                                                     #
################################################################################

#  Build test utilities.
add_executable (xwout_diff)

target_link_libraries (xwout_diff PUBLIC stell)
target_compile_features (xwout_diff PRIVATE cxx_std_11)

add_subdirectory (Testing)
