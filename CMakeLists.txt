cmake_minimum_required(VERSION 2.8.3)
project(optimize)

set(CMAKE_CXX_FLAGS "-std=c++11")

find_package(Eigen3 REQUIRED)
include_directories(
  ${EIGEN3_INCLUDE_DIR}
)

find_package( Ceres REQUIRED )
include_directories( ${CERES_INCLUDE_DIRS} )

add_executable(optimize optimize.cpp)
target_link_libraries( optimize ${CERES_LIBRARIES})
