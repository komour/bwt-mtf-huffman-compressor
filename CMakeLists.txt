cmake_minimum_required(VERSION 3.20)
project(BWT_MTF_Huffman)

set(CMAKE_CXX_STANDARD 23)
set (CMAKE_CXX_FLAGS "-O3")
#set (CMAKE_CXX_FLAGS "-DFULLPIPELINE")
#set (CMAKE_CXX_FLAGS "-DCOMPRESS")
#set (CMAKE_CXX_FLAGS "-DDECOMPRESS")

add_executable(BWT_MTF_Huffman main.cpp)

