cmake_minimum_required(VERSION 3.0.2)
project(ur_ik_solver)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(catkin REQUIRED COMPONENTS
  eigen_conversions
  ik_solver_msgs
  ur_kinematics
  roscpp
  tf_conversions
  ik_solver
)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ur_ik_solver
  CATKIN_DEPENDS eigen_conversions ur_kinematics ik_solver_msgs roscpp tf_conversions ik_solver
#  DEPENDS system_lib
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
)

# ur_kinematics individual libraries
find_library(libur3_kin ur3_kin)
find_library(libur5_kin ur5_kin)
find_library(libur10_kin ur10_kin)
find_library(libur3e_kin ur3e_kin)
find_library(libur5e_kin ur5e_kin)
find_library(libur10e_kin ur10e_kin)

add_library(ur3_ik_solver
  src/${PROJECT_NAME}/ur_ik_solver.cpp
)
add_dependencies(ur3_ik_solver ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

add_library(ur5_ik_solver
  src/${PROJECT_NAME}/ur_ik_solver.cpp
)
add_dependencies(ur5_ik_solver ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

add_library(ur10_ik_solver
  src/${PROJECT_NAME}/ur_ik_solver.cpp
)
add_dependencies(ur10_ik_solver ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

add_library(ur3e_ik_solver
  src/${PROJECT_NAME}/ur_ik_solver.cpp
)
add_dependencies(ur3e_ik_solver ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

add_library(ur5e_ik_solver
  src/${PROJECT_NAME}/ur_ik_solver.cpp
)
add_dependencies(ur5e_ik_solver ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

add_library(ur10e_ik_solver
  src/${PROJECT_NAME}/ur_ik_solver.cpp
)
add_dependencies(ur10e_ik_solver ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})


 # add_executable(${PROJECT_NAME}_node src/node.cpp)
 #set_target_properties(${PROJECT_NAME}_node PROPERTIES OUTPUT_NAME ur_ik PREFIX "")
 #add_dependencies(${PROJECT_NAME}_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

target_link_libraries(ur3_ik_solver
    ${libur3_kin} # Order matters!! Place before catkin_LIBRARIES to specify which UR is selected
    ${catkin_LIBRARIES}
)

target_link_libraries(ur5_ik_solver
    ${libur5_kin} # Order matters!!
    ${catkin_LIBRARIES}
)

target_link_libraries(ur10_ik_solver
    ${libur10_kin} # Order matters!!
    ${catkin_LIBRARIES}
)

target_link_libraries(ur3e_ik_solver
    ${libur3e_kin} # Order matters!!
    ${catkin_LIBRARIES}    
)

target_link_libraries(ur5e_ik_solver
    ${libur5e_kin} # Order matters!!
    ${catkin_LIBRARIES}
)

target_link_libraries(ur10e_ik_solver
    ${libur10e_kin} # Order matters!!
    ${catkin_LIBRARIES}
)
