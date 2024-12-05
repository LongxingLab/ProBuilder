#########################################################################
# File Name: compile.sh
# Author: longxing
# mail: longxing@uw.edu
# Created Time: 2021年09月28日 星期二 12时51分30秒
#########################################################################
#!/bin/bash

c++ -O3 -Wall -shared -std=c++17 -I../external/ -fPIC $(python3 -m pybind11 --includes) xbin.cpp -o xbin$(python3-config --extension-suffix)
