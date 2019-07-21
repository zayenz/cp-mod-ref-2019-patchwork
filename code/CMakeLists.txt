cmake_minimum_required(VERSION 3.10)
project( Patchwork )

# C++ standard version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=leak -g")

# The version number.
set (Patchwork_VERSION_MAJOR 0)
set (Patchwork_VERSION_MINOR 1)
set (Patchwork_VERSION_PATCH 0)

add_subdirectory (src)


