cmake_minimum_required(VERSION 3.16)
project(ur_ik_solver)

set(CMAKE_BUILD_TYPE Debug)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


set(ROS_VERSION $ENV{ROS_VERSION})
if(ROS_VERSION STREQUAL "")
  message(FATAL_ERROR "No supported version of ROS found!")
endif()
if(ROS_VERSION EQUAL 1)
  message(FATAL_ERROR "ROS 1 support missing")
endif()

find_package(Eigen3 REQUIRED)

# Define map {ur_model : ur_definition_macro}
list(APPEND UR_MODELS
  ur3 ur3e
  ur5 ur5e
  ur10 ur10e
  ur16e
  ur20
  ur30)
set(urmap_ur3 UR3_PARAMS)
set(urmap_ur3e UR3e_PARAMS)
set(urmap_ur5 UR5_PARAMS)
set(urmap_ur5e UR5e_PARAMS)
set(urmap_ur10 UR10_PARAMS)
set(urmap_ur10e UR10e_PARAMS)
set(urmap_ur16e UR16e_PARAMS)
set(urmap_ur20 UR20_PARAMS)
set(urmap_ur30 UR30_PARAMS)

foreach(ur_type IN LISTS UR_MODELS)
  message("ur_ik_solver: Building ROS 2 Version")
  find_package(ament_cmake REQUIRED)
  find_package(ik_solver   REQUIRED)
  find_package(pluginlib   REQUIRED)

  add_library(${ur_type}_ik_solver SHARED
  src/ur_ik_solver/ur_ik_solver.cpp
  src/ur_ik_solver/ur_kinematics/ur_kin.cpp)
  target_include_directories(${ur_type}_ik_solver PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>")
  ament_target_dependencies(${ur_type}_ik_solver PUBLIC rclcpp
                                                 ik_solver
                                                 pluginlib)
  target_compile_definitions(${ur_type}_ik_solver PUBLIC ${urmap_${ur_type}})
  install(TARGETS ${ur_type}_ik_solver
          EXPORT export_ur_ik_solver
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib
          RUNTIME DESTINATION bin)
endforeach()

pluginlib_export_plugin_description_file(ik_solver plugin2.xml)

ament_export_targets(export_ur_ik_solver)
ament_export_dependencies(ik_solver pluginlib)
ament_package()

install(DIRECTORY config launch
  DESTINATION share/${PROJECT_NAME}
)
