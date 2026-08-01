#!/usr/bin/env python3
import sys
import os

# Add libs to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../libs')))

from dynamixel_easy_sdk.connector import Connector
from dynamixel_easy_sdk.dynamixel_error import DxlRuntimeError

def scan_motors():
    # Settings from config or defaults
    port = '/dev/ttyUSB4'
    baud = 57600
    
    print(f"Scanning on {port} at {baud} baud...")
    
    try:
        connector = Connector(port, baud)
        print("Port opened successfully.")
        
        print("Performing Broadcast Ping...")
        ids = connector.broadcastPing()
        
        if not ids:
            print("No motors found.")
        else:
            print(f"Found {len(ids)} motors:")
            for mid in ids:
                try:
                    # Ping individual to get model number
                    model = connector.ping(mid)
                    print(f" - ID: {mid}, Model: {model}")
                except Exception as e:
                    print(f" - ID: {mid}, Error: {e}")
                    
        connector.closePort()
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == '__main__':
    scan_motors()
