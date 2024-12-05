cmake_minimum_required (VERSION 3.0)

project (WeFold)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -static -w -fno-omit-frame-pointer -Wall -Wextra -Werror" )

string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)
if( ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Og -O0 -frtti" )
    message("Building the debug mode ....... ")
else()
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3" )
    message("Building the release mode ....... ")
endif()

option(NOT_USE_PHMAP "use the parallel hash map" OFF)
if (NOT_USE_PHMAP)
    message("Building folder without using parallel hash map")
else()
    add_definitions(-DUSE_PHMAP)
endif()

if(EXISTS "/home/longxing/include/boost_1_67_0/rst.css")
    include_directories( "/home/longxing/include/boost_1_67_0/" )
    link_directories( "/home/longxing/include/boost_1_67_0/stage/lib" )
else()
    include_directories( "/home/caolongxingLab/caolongxing/include/boost_1_67_0/" )
    link_directories( "/home/caolongxingLab/caolongxing/include/boost_1_67_0/libs/" )
endif()

include_directories(external)
include_directories(src)

# need to compile?? interesting
add_subdirectory(external/gzip)
add_subdirectory(external/dssp)

add_subdirectory(src)
