cmake_minimum_required(VERSION 3.16)
project(ur_ik_solver)

set(CMAKE_BUILD_TYPE Debug)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(ament_cmake REQUIRED)
find_package(ik_solver   REQUIRED)

add_library(ur5_ik_solver SHARED
  src/ur_ik_solver/ur_ik_solver.cpp
  src/ur_ik_solver/ur_kinematics/ur_kin.cpp)
target_include_directories(ur5_ik_solver PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>")
ament_target_dependencies(ur5_ik_solver PUBLIC
  rclcpp
  ik_solver
  pluginlib)
target_compile_definitions(ur5_ik_solver PUBLIC UR5_PARAMS)

pluginlib_export_plugin_description_file(ik_solver plugin.xml)

install(TARGETS ur5_ik_solver
        EXPORT export_ur5_ik_solver
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin)

ament_export_targets(export_ur5_ik_solver)
ament_export_dependencies(ik_solver pluginlib)
ament_package()
