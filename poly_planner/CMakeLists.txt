cmake_minimum_required(VERSION 2.8.3)

project(poly_planner)

set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -fPIC")

find_package(Eigen3 REQUIRED)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  message_generation
  message_filters
  std_msgs
  geometry_msgs
  nav_msgs
  tf
  image_transport
  cv_bridge
  quadrotor_msgs
  visualization_msgs
  decomp_ros_utils
)

include_directories(
    ${catkin_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    include
)

catkin_package()

add_executable(${PROJECT_NAME}_node src/fsto_node.cpp src/fsto.cpp)

target_link_libraries(${PROJECT_NAME}_node
  ${catkin_LIBRARIES}
)
