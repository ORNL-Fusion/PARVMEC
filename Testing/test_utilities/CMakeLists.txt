target_sources(xwout_diff
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/wout_diff.cpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vmec_test_commandline_parser.hpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vmec_test_commandline_parser.cpp>
)
