target_sources(xdiff_wout
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/wout_diff.cpp>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vmec_test_commandline_parser.hpp>
)
