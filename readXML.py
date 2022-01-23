#!/usr/bin/env python
#
# JDY 
# 18JDY2020
#

import sys,os,glob,re
import numpy as np
import defusedxml.ElementTree as ET
from ase import Atoms


def get_xml_data(dir_name):

  global xml_nprocs,xml_ndiag,xml_nat, xml_alat,atom, \
  xml_species,xml_atom_index,xml_atoms

  tree = ET.parse(dir_name)
  root = tree.getroot()

# parallel info
  for i in root.findall("./parallel_info/nprocs"):
    xml_nprocs = i.text
  for i in root.findall("./parallel_info/ndiag"):
    xml_ndiag = i.text

# output atomic_structure
  for i in root.findall("./output/atomic_structure"):
    xml_nat = int(i.attrib['nat'])
    xml_alat= i.attrib['alat']
  n = 0
  xml_atoms = np.zeros((xml_nat,5)) # index, atomic number, pos x,y,z
  for i in root.findall("./output/atomic_structure/atomic_positions/atom"):
    tmp = i.text
    tmp = np.array(tmp.split()).astype(float)
    xml_species = Atoms(symbols=i.attrib['name']).numbers
    xml_atom_index = int(i.attrib['index'])
    xml_atoms[n] = [xml_atom_index,xml_species,tmp[0],tmp[1],tmp[2]]
    n += 1

  return

if __name__ == "__main__":

  get_xml_data(sys.argv[1])
  print(int(xml_atoms[0,0]))
