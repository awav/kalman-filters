project(ekf)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++14 -g -Werror)

set(SOURCES
    src/FusionEKF.cpp
    src/kalman_filter.cpp
    src/main.cpp
    src/tools.cpp)

include_directories(src/)

add_executable(${PROJECT_NAME} ${SOURCES})

#set(SOURCES src/eigen_test.cpp)
#add_executable(eigen_test ${SOURCES})
