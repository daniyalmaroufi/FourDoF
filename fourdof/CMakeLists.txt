cmake_minimum_required(VERSION 3.0.2)
project(fourdof)

find_package(catkin REQUIRED COMPONENTS
  rospy
  std_msgs
  message_generation
)

catkin_python_setup()

add_message_files(
  FILES
  JointCommand.msg
  JointState.msg
)

add_service_files(
  FILES
  SetHome.srv
  SetTorque.srv
)

generate_messages(
  DEPENDENCIES
  std_msgs
)

catkin_package(
  CATKIN_DEPENDS rospy std_msgs message_runtime
)

catkin_install_python(PROGRAMS
  fourdof/rotational_controller_node.py
  fourdof/translational_controller_node.py
  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

install(DIRECTORY launch/
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}/launch
)

install(DIRECTORY config/
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}/config
)

install(DIRECTORY libs/dynamixel_easy_sdk/control_table/
  DESTINATION ${CATKIN_PACKAGE_PYTHON_DESTINATION}/dynamixel_easy_sdk/control_table
  FILES_MATCHING PATTERN "*.model"
)
