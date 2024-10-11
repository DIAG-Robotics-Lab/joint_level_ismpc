cmake_minimum_required(VERSION 2.8.12)

project(jointISMPC)

# should be unnecessary for plain CMake if bashrc correctly configured
list(APPEND CMAKE_PREFIX_PATH "/opt/openrobots/")

find_package(DART 6.1.1 REQUIRED COMPONENTS utils-urdf gui-osg)
find_package(pinocchio REQUIRED)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

#add_compile_options(-O3)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
add_compile_options(${pinocchio_CFLAGS_OTHER})

include_directories(include)
include_directories(${DART_INCLUDE_DIRS})
include_directories(${pinocchio_INCLUDE_DIRS})

link_directories(${PROJECT_NAME} PUBLIC ${pinocchio_LIBRARY_DIRS})

add_subdirectory(labrob_qpsolvers)
add_subdirectory(cpp-utils build/cpp-utils)


file(GLOB srcs "src/*.cpp")
add_executable(${PROJECT_NAME} ${srcs})

target_link_libraries(${PROJECT_NAME} labrob_qpsolvers ${DART_LIBRARIES} pinocchio::pinocchio cpp-utils)

target_compile_options(${PROJECT_NAME} PUBLIC ${pinocchio_CFLAGS_OTHER})

