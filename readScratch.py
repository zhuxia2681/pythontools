#!/usr/bin/env python
#
# JDY 
#


import sys,os,glob,re,datetime
import numpy as np
import h5py
from  scipy.io import FortranFile as FF
import xml.etree.ElementTree as ET

# in future hdf5? 
g = FF(sys.argv[1],"r")
lines = FF.read_record(g,dtype=complex)

print(lines[0])

