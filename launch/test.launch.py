from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument, ExecuteProcess
from launch.conditions import IfCondition
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node
from ament_index_python.packages import get_package_share_directory
import os

def generate_launch_description():
    # Arguments

    # Get config file path
    config_file = os.path.join(
        get_package_share_directory('ur_ik_solver'),
        'config',
        'params.yaml'
    )

    ik_node = Node(
        package='ik_solver',
        executable='ik_solver_node',
        name='ur_ik',
        output='screen',
        namespace='ur_ik_solver'
        # parameters=[config_file]
    )

    load_params = ExecuteProcess(cmd=['cnr_param_server', '-p', str(config_file)])



    return LaunchDescription([
        load_params,
        ik_node
    ])
