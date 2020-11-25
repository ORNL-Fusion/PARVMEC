target_sources(xdiff_wout
               PRIVATE
               $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/diff_wout.cpp>
)
