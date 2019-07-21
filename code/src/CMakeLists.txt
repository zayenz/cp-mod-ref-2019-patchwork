find_library(GECODE Gecode)
if (NOT GECODE)
    message(FATAL_ERROR "Gecode not found")
endif()

# add the binary tree to the search path for include files
# so that we will find config.h
include_directories("${PROJECT_BINARY_DIR}/src/")
include_directories ("${PROJECT_SOURCE_DIR}/patchwork")
include_directories("${GECODE}")
link_libraries()
add_subdirectory (patchwork)

set(GECODE_LIBS  ${GECODE})

add_executable(patchwork-test-placements test-placements.cpp)
target_link_libraries (patchwork-test-placements PatchworkLib ${GECODE_LIBS} )

add_executable(patchwork-test-game-play-statistics test-game-play-statistics.cpp)
target_link_libraries (patchwork-test-game-play-statistics PatchworkLib ${GECODE_LIBS} )
