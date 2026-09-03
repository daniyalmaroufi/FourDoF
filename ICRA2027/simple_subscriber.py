#!/usr/bin/env python  
import rospy

import math
import tf2_ros
import geometry_msgs.msg
import turtlesim.srv
import os 
from transform import quaternion_matrix
import yaml
import numpy as np
import struct

if __name__ == '__main__':
    rospy.init_node('tool_tracking_listener')
    data = []
    tfBuffer = tf2_ros.Buffer()
    listener = tf2_ros.TransformListener(tfBuffer)

    rate = rospy.Rate(10.0)
    
    while not rospy.is_shutdown():
        try:
            # ans = input('Enter id to continue (q to quit)...')
            # if ans == 'q':
            #     print('Exiting because quit requested')
            #     break
            trans = tfBuffer.lookup_transform('tool3', 'tool2', rospy.Time())
            trans = trans.transform
            mat = quaternion_matrix([trans.rotation.x, trans.rotation.y, trans.rotation.z,trans.rotation.w])
            mat[:3, 3] = np.array([trans.translation.x, trans.translation.y,trans.translation.z])
            print(mat)



            # with open('/home/shuojue/handeye_data/pos'+ans+'.yaml', 'w') as f:
            #     data_str = '[\n'
            #     for m in mat:
            #         data_str +='['
            #         data_str += ', '.join([str(x) for x in m]) +'],\n'
            #     data_str +=']'
            #     data = {
            #     '__version__':{'serializer':1, 'data':1},
            #     'FloatMatrix':{'Data': data_str}
            #     }
            #     yaml.dump(data, f,sort_keys=False) 
            # trans1 = tfBuffer.lookup_transform('optical_tracker', 'tool2', rospy.Time())
        except (tf2_ros.LookupException, tf2_ros.ConnectivityException, tf2_ros.ExtrapolationException):
        
            continue
            
            # f.write('['+str(trans.transform.translation.x)+','+\
            #             str(trans.transform.translation.y)+','+
            #             str(trans.transform.translation.z)+'],\n')

        # rate.sleep()
