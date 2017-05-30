project(ekf)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++14)

set(CXX_FLAGS "-Wall -g -Werror")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  include_directories(/usr/local/include)
  include_directories(/usr/local/opt/openssl/include)
  link_directories(/usr/local/lib)
  link_directories(/usr/local/opt/openssl/lib)
  link_directories(/usr/local/Cellar/libuv/1.11.0/lib)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

set(SOURCES
    src/main.cpp
    src/kalman_filter.cpp
    src/FusionEKF.cpp
    src/tools.cpp)

include_directories(src/)

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} z ssl uv uWS)

#set(SOURCES src/eigen_test.cpp)
#add_executable(eigen_test ${SOURCES})
