#!/usr/bin/python
# coding=utf-8
"""
Created on Tue, 09 Jul 2017 11:01:31

Short Description:
callback function to close the serial port/ save log and history files of the ublox frame's
from the web page

@author: ELWAN

"""
import os

print "Content-type:text/html\r\n\r\n"
print '<meta http-equiv="Content-type" content="text/html;charset=UTF-8">'

path = os.path.dirname(os.path.realpath(__file__)) + "/"
print path
f = open(path + "offline", "w")
f.close()
