cmake_minimum_required(VERSION 3.31)
project(op)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_INSTALL_PREFIX /Users/plewis/bin)

add_executable(op src/main.cpp)
include_directories(/Users/plewis/Documents/software/libraries/boost_1_86_0)
include_directories(/Users/plewis/Documents/software/libraries/ncl-2.1.18)

find_library(BOOST_PROGRAM_OPTIONS libboost_program_options.a /Users/plewis/lib/static)
target_link_libraries(op LINK_PUBLIC ${BOOST_PROGRAM_OPTIONS})

find_library(BOOST_FILESYSTEM libboost_filesystem.a /Users/plewis/lib/static)
target_link_libraries(op LINK_PUBLIC ${BOOST_FILESYSTEM})

find_library(BOOST_SYSTEM libboost_system.a /Users/plewis/lib/static)
target_link_libraries(op LINK_PUBLIC ${BOOST_SYSTEM})

find_library(NEXUS_CLASS_LIBRARY libncl.a /Users/plewis/lib/static)
target_link_libraries(op LINK_PUBLIC ${NEXUS_CLASS_LIBRARY})

add_subdirectory(Google_tests)
