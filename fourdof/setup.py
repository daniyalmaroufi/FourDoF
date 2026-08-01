from distutils.core import setup
from catkin_pkg.python_setup import generate_distutils_setup

d = generate_distutils_setup(
    packages=['fourdof', 'dynamixel_easy_sdk'],
    package_dir={
        '': '.', 
        'dynamixel_easy_sdk': 'libs/dynamixel_easy_sdk'
    }
)

setup(**d)
