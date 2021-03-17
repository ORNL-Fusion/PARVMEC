#  The python interpreter may not be in the same location on all systems.
find_package (PythonInterp)

set (TEST_SCRIPT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities)

add_subdirectory (tests)
add_subdirectory (test_utilities)
