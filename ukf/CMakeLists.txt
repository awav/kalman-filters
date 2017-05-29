project(ukf)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++14 -g -Werror)

set(SOURCES
   src/ukf.cpp
   src/main.cpp
   src/tools.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})
