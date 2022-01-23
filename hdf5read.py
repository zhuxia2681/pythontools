#!/usr/bin/env python
#
# JDY 
# 09JDY2020
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

  for k in f:
    print(k)
  ALAT = f['mf_header/crystal/alat']
  ADOT = f['mf_header/crystal/adot']
  APOS = f['mf_header/crystal/apos']
  ATYP = f['mf_header/crystal/atyp']
#  ATYP = f['']
  FFTGRID = f['mf_header/gspace/FFTgrid']
  NG = f['mf_header/gspace/ng']

# kpoints
  NGK = f['mf_header/kpoints/ngk']
  OCC= f['mf_header/kpoints/occ']
  KGRID = f['mf_header/kpoints/kgrid']
  RK = f['mf_header/kpoints/rk']

# symmetry
  CELL_SYMMETRY = f['mf_header/symmetry/cell_symmetry']
  MTRX = f['mf_header/symmetry/mtrx']
  NTRAN = f['mf_header/symmetry/ntran']
  TNP = f['mf_header/symmetry/tnp']

# wfns
  GVECS = f['wfns/gvecs']
  COEFFS = f['wfns/coeffs']


#version number
  VERSIONNUMBER = f['mf_header/versionnumber']

  print("--- adot (r-metric tensor) ---")
  print(len(ADOT))
  print(ADOT[:])
  print()
  print("--- apos ---")
  print(len(APOS))
  print(APOS[:])
  print()
  print("--- alat ---")
  print(ALAT[...])
  print()
  print("--- atyp ---")
  print(ATYP[...])
  print()
  print("--- FFTgrid ---")
  print(len(FFTGRID))
  print(FFTGRID[:])
  print()
  print("--- ng ---")
  print(NG[...])
  print()
  print("--- ngk ---")
  print(NGK[...])
  ngk = np.sum(NGK)
  print(ngk)
  print()
  print("--- cell_symmetry ---")
  print(CELL_SYMMETRY[...])
  print()
  print("--- s matrix ---")
  print(len(MTRX))
  print(MTRX[:])
  print()
  print("--- translation ---")
  print(NTRAN[...])
  print()
  print("--- tnp ---")
  print(len(TNP))
  print(TNP[...])
  print()
  print("--- version number ---")
  print(VERSIONNUMBER[...])
  print()
  print("--- gvecs ---")
  print(len(GVECS[...]))
  print(np.shape(GVECS))
  print()
  print("--- coeffs ---")
  print(len(COEFFS[...]))
  print(COEFFS[0,0,0,0]) # [:,1,:,:] is spin
  print(COEFFS[0,0,0,1]) # [:,:,:,1] is Re or Im
  print(COEFFS[0,1,12904,0])
  print(COEFFS[0,1,12904,1])
  print(np.shape(COEFFS))
#  print((COEFFS[0,0,0,:]))
  print()
  print("--- kgird ---")
  print(len(KGRID[...]))
  print(np.shape(KGRID))
  print((KGRID[...]))
  print("--- rk ---")
  print(np.shape(RK))
  print((RK[...]))


f.close()
