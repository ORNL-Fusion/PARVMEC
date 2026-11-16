#  The python interpreter may not be in the same location on all systems.
find_package (Python3 COMPONENTS Interpreter REQUIRED)

set (TEST_SCRIPT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_utilities)

add_subdirectory (tests)
