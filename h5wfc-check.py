#!/usr/bin/env python
#
# JDY 
# 17JDY2020
#

import sys,os,glob,re
import numpy as np
import h5py

def PrintAllObjects(name):

  print(name)

  return

def PrintAllDataset(name, obj):

  print(name)

  return


pwd = os.getcwd()
hdf5path = pwd + "/" + sys.argv[1]


with h5py.File(hdf5path, 'r') as f:

  f.visititems(PrintAllDataset)

# wfns
  MI = f['MillerIndices']
  EVC = f['evc']
  IK = f.attrs['ik']
  XK = f.attrs['xk']
  ISPIN = f.attrs['ispin']
  GAMMA_ONLY = f.attrs['gamma_only']
  SCALE_FACTOR = f.attrs['scale_factor']
  NGW = f.attrs['ngw']
  IGWX = f.attrs['igwx']
  NPOL = f.attrs['npol']
  NBND = f.attrs['nbnd']


  print("--- MillerIndices ---")
  print(len(MI[...]))
  print(np.shape(MI))
  print(MI[...])
  print()
  print("--- evc ---")
  print(len(EVC[...]))
  print(np.shape(EVC))
  print((EVC[...]))
  print(np.argmax(EVC[:,:]))
  print()
  print("--- ik ---")
  print(IK[...])
  print("--- xk ---")
  print(XK[...])
  print("--- ispin ---")
  print(ISPIN[...])
  print("--- gamma_only ---")
  print(GAMMA_ONLY)
  print("--- scale_factor ---")
  print(SCALE_FACTOR)
  print("--- ngw ---")
  print(NGW[...])
  print(NGW)
  print("--- igwx ---")
  print(IGWX[...])
  print("--- npol ---")
  print(NPOL[...])
  print("--- nbnd ---")
  print(NBND[...])


f.close()
