cmake_minimum_required(VERSION 3.1)
project(cpp-utils)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


add_library(${PROJECT_NAME} INTERFACE)

find_package(Eigen3 REQUIRED)

target_include_directories(${PROJECT_NAME} INTERFACE include/)

