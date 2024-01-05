cmake_minimum_required(VERSION 3.5)

project(dsla LANGUAGES CXX VERSION 1.0.0 DESCRIPTION "distributed sparse linear algebra")
add_library(dsla SHARED)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


target_include_directories(dsla PRIVATE inc)
add_subdirectory(utils)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -g -Werror -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -Wall -Werror -pedantic")

