cmake_minimum_required(VERSION 3.0.2)
project(labrob_qpsolvers)

## Compile as C++11, supported in ROS Kinetic and newer
add_compile_options(-std=c++11)

add_library(${PROJECT_NAME} INTERFACE)

# Differentiate between catkin and plain CMake build
if(CATKIN_DEVEL_PREFIX)

    find_package(catkin REQUIRED)
    ## System dependencies are found with CMake's conventions
    # find_package(Boost REQUIRED COMPONENTS system)
    find_package(blasfeo REQUIRED)
    find_package(hpipm REQUIRED)

    # Manually set INCLUDE_DIRS and LIBRARIES since not handled in blasfeo and hpipm
    set(blasfeo_INCLUDE_DIRS /opt/blasfeo/include)
    set(hpipm_INCLUDE_DIRS /opt/hpipm/include)
    set(blasfeo_LIBRARIES /opt/blasfeo/lib/libblasfeo.a)
    set(hpipm_LIBRARIES /opt/hpipm/lib/libhpipm.a)


    catkin_package(
            INCLUDE_DIRS include
            DEPENDS blasfeo hpipm
    )

    target_include_directories(${PROJECT_NAME} INTERFACE
            include
            ${catkin_INCLUDE_DIRS}
            ${blasfeo_INCLUDE_DIRS}
            ${hpipm_INCLUDE_DIRS}
            )

    target_link_libraries(${PROJECT_NAME} INTERFACE
            ${catkin_LIBRARIES}
            blasfeo
            hpipm
            )

else()

    # Manually set INCLUDE_DIRS and LIBRARIES since not handled in blasfeo and hpipm
    set(blasfeo_INCLUDE_DIRS /opt/blasfeo/include)
    set(hpipm_INCLUDE_DIRS /opt/hpipm/include)
    set(blasfeo_LIBRARIES_PATH /opt/blasfeo/lib)
    set(hpipm_LIBRARIES_PATH /opt/hpipm/lib)

    find_library(BLASFEO_LIB blasfeo HINTS ${blasfeo_LIBRARIES_PATH} REQUIRED)
    find_library(HPIPM_LIB hpipm HINTS ${hpipm_LIBRARIES_PATH} REQUIRED)

    target_include_directories(${PROJECT_NAME} INTERFACE
            include/
            ${hpipm_INCLUDE_DIRS}
            ${blasfeo_INCLUDE_DIRS}
            )

    target_link_libraries(${PROJECT_NAME} INTERFACE
            ${HPIPM_LIB}
            ${BLASFEO_LIB}
            )

endif()