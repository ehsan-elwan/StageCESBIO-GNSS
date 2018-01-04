#!/usr/bin/python
# coding=utf-8
"""
Created on Tue, 09 Jul 2017 11:01:31

Short Description:
callback function to start listing on serial port from the web page

@author: ELWAN

"""
from Serialreader import *
print "Content-type:text/html\r\n\r\n"
print '<meta http-equiv="Content-type" content="text/html;charset=UTF-8">'

start_serial('/dev/ttyACM0', 9600, 1)
