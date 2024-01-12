cmake_minimum_required (VERSION 3.4)
project (seqan3_tutorial CXX)
 
# add the Sharg Parser to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../sharg-parser/build_system")
 
find_package (sharg 1.0 REQUIRED)

# add_executable (seqomplexity seqomplexity.cpp)
# target_link_libraries (seqomplexity sharg::sharg)

# add_executable (GC_content GC_content.cpp)
# target_link_libraries (GC_content sharg::sharg)

add_executable (fast_sequence_complexity fast_sequence_complexity.cpp)
target_link_libraries (fast_sequence_complexity sharg::sharg)

add_executable (fast_GC_content fast_GC_content.cpp)
target_link_libraries (fast_GC_content sharg::sharg)