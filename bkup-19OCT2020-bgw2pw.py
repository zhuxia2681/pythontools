#!/usr/bin/env python3
#
# bgw2pw 
# original code has developped by D. Strubbe
# JDY 
# 16JDY2020
#

#
# - write_evc
#   read wfc from WFN which is BGW fromat and
#   write out to evec file in prefix.save/K000**. 
#
# - write_cd
#
# - check_inversion
#   check flavor selection is appropriate or not.
#

import datetime as dt
import sys,os,glob,re
import numpy as np
import h5py
import defusedxml.ElementTree as ET
import subprocess as sp

def get_xml_data(path):

  tree = ET.parse(path +'/data-file-schema.xml')
  root = tree.getroot()

  return root


# input file

prefix = '1H-MoS2' 
outdir = './work'
real_or_complex = 2
wfng_flag = False
# set sapo wfc hdf5 file
wfng_file = "./4.1/wfn.sapo.cplx.hdf5"
wfng_nband = 2538
rhog_flag = False
rhog_file = "RHO"
eig_flag = True
# set sapo directory where ener.dat files
# exist.
eig_dir = "./4.1"

#
restart_dir = outdir + "/" + prefix + ".save"

# constant
eps6 = 1.0E-6

# get parameter for header
with h5py.File(wfng_file, 'r') as h:
  tmp = h['mf_header/flavor']
  flavor = tmp[...]
  tmp = h['mf_header/kpoints/nspin']
  nspin = tmp[...]
  tmp = h['mf_header/kpoints/nrk']
  nrk = tmp[...]
  tmp = h['mf_header/kpoints/mnband']
  mnband = tmp[...]
  tmp = h['mf_header/kpoints/ngkmax']
  ngkmax = tmp[...]
  tmp = h['mf_header/kpoints/ecutwfc']
  ecutwfc = tmp[...]
  tmp = h['mf_header/kpoints/ngk']
  ngk = tmp[...]
  tmp = h['mf_header/kpoints/occ']
  occ = tmp[...]
  tmp = h['mf_header/gspace/ng']
  ng = tmp[...]
  tmp = h['mf_header/gspace/ecutrho']
  ecutrho = tmp[...]
  tmp = h['mf_header/gspace/FFTgrid']
  FFTgrid = tmp[...]
  tmp = h['mf_header/symmetry/cell_symmetry']
  cell_symmetry = tmp[...]
  tmp = h['mf_header/crystal/nat']
  na = tmp[...]
  tmp = h['mf_header/crystal/celvol']
  celvol = tmp[...]
  tmp = h['mf_header/crystal/alat']
  alat = tmp[...]
  tmp = h['mf_header/crystal/avec']
  avec = tmp[...]
  tmp = h['mf_header/crystal/adot']
  adot = tmp[...]
  tmp = h['mf_header/crystal/recvol']
  recvol = tmp[...]
  tmp = h['mf_header/crystal/blat']
  blat = tmp[...]
  tmp = h['mf_header/crystal/bvec']
  bvec = tmp[...]
  tmp = h['mf_header/crystal/bdot']
  bdot = tmp[...]
  tmp = h['mf_header/symmetry/ntran']
  nsym = tmp[...]
  tmp = h['mf_header/symmetry/mtrx']
  s = tmp[...]
if wfng_flag:
  with h5py.File(wfng_file, 'r') as h:
    print("read wfc from hdf5 of BGW")
    tmp = h['wfns/coeffs']
    coeffs = tmp[...]
    tmp = h['wfns/gvecs']
    gvecs = tmp[...]
h.close()

if eig_flag:
# backup data.xml
  cmd = "cp " + restart_dir + "/data-file-schema.xml " + restart_dir +"/bkup-data.xml"
  proc = sp.Popen(cmd,shell=True)
  proc.wait()

  eig = np.zeros((wfng_nband*nspin,nrk))
  files = glob.glob(eig_dir+"/ener*")
  if len(files) != nrk*nspin:
    print("system inconsistent : read ener.dat")
    sys.exit()

  eig = np.zeros((mnband*nspin,nrk))
  for i in range(nrk):
    fnameup = "ener_k"+str(i+1)+"_s1.dat"
    fnamedw = "ener_k"+str(i+1)+"_s2.dat"
    fup = open(eig_dir+"/"+fnameup)
    fdw = open(eig_dir+"/"+fnamedw)
    linesup = fup.readlines()    
    linesdw = fdw.readlines()
    print("ik = ",i+1)
    print("nbup = ",len(linesup)-1)
    print("nbdw = ",len(linesdw)-1)
    n = 0
    m = 1
    for k in range(mnband):
      linesup[m] = linesup[m].lstrip()
      linesup[m] = re.sub(r"\s+"," ",linesup[m])
      tmpup = linesup[m].split()
      eig[n,i] = float(tmpup[1])
      m += 1
      n += 1
    m = 1
    for k in range(mnband):
      linesdw[m] = linesdw[m].lstrip()
      linesdw[m] = re.sub(r"\s+"," ",linesdw[m])
      tmpdw = linesdw[m].split()
      eig[n,i] = float(tmpdw[1])
      m += 1
      n += 1

  with h5py.File(wfng_file, 'r') as h:
    tmp= h['mf_header/kpoints/occ']
    occ = tmp[...]
  h.close()
  
#  eig_root = put_xml_data(restart_dir)
  


def volume(alat, at1, at2, at3):

  omega = 0.0
  s = 1.0
  i = 0
  j = 1
  k = 2
  # where s is sign of a permutation
  for cnt in range(3):
    for iperm in range(3): 
      omega = omega + s * at1[i] * at2[j] * at3[k]
      l = i
      i = j
      j = k
      k = l
    i = 1
    j = 0
    k = 2
    s = -s

  omega = np.abs(omega) * alat ** 3

  return omega 

def check_inversion(real_or_complex, ntran, mtrx, nspin, warn, real_need_inv):

# warn must be set false when one want to suppress warnings

  mtrx = np.zeros((3,3,48))

  invflag = 0
  for isym in range(ntran):
    itest = 0
    for ii in range(3):
      for jj in range(3):
        if ii == jj:
          itest = itest + (mtrx[ii,jj,isym]+1)**2
        else:
          itest = itest + mtrx[ii,jj,isym]**2

  if itest == 0:
    invflag = invflag + 1
  if invflag > 1:
    print("check_inversion : \
           More than one inversion symmetry operation is present",invflag)
    sys.exit()  
  
  
  if real_or_complex == 2:
    if invflag != 0 and warn:
      print("WARNING : Inversion symmetry is present. The real version would be faster.")
  else:
    if invflag == 0:
      if real_need_inv :
        print("check_inversion : The real version cannot be used without inversion symmetry.")
      print("WARNING: Inversion symmetry is absent in symmetries used to reduce k-grid.")
      print("Be sure inversion is still a spatial symmetry, or you must use complex version instead.")
    if nspin == 2:
      print("check_inversion : Time-reversal symmetry is absent in spin-polarized calculation.\
 Complex version must be used.")

  return

def write_evc(input_file_name,real_or_complex,wfng_nband,output_dir_name ):

  print("check inversion")
  check_inversion(real_or_complex,nsym,s,nspin,True,True)
  print("check inversion done")

  root = get_xml_data(output_dir_name)
  for i in root.findall("./output/magnetization/lsda"):
    xml_lsda = i.text
  for i in root.findall("./output/magnetization/spinorbit"):
    xml_soc = i.text

  if flavor == 1 and real_or_complex == 1:
    stitle = "WFN-Real" 
  if flavor == 2 and real_or_complex == 2:
    stitle = "WFN-Complex" 
  else:
    print("error write_evc : file header error")
    sys.exit()  
    
  sdate = str(dt.date.today())
  stime = str(dt.datetime.now().time())

  if xml_lsda == "true" :
    xml_nspin = 2
  else:
    xml_nspin = 1
  if xml_soc == "true" :
    xml_nspin = 4
  if nspin != xml_nspin:
    print("error write_evc : ns")
    sys.exit()

  for i in root.findall("./output/basis_set/ngm"):
    xml_ng = int(i.text)
  if ng != xml_ng:
    print("error write_evc : ng")
    sys.exit()
  
  for i in root.findall("./output/atomic_structure"):
    xml_na = int(i.attrib["nat"])
  if na != xml_na:
    print("error write_evc : na")
    sys.exit()

  for i in root.findall("./output/band_structure/nks"):
    xml_nks = int(i.text)
#JDY 17OCT2020
# nks in the data.~~xml file is indepenedent 
# with nspin. 
  if nrk != xml_nks :
#  if nrk != xml_nks / nspin:
    print("error write_evc : nk")
    sys.exit()

  for i in root.findall("./output/basis_set/fft_grid"):
    xml_nr1 = int(i.attrib["nr1"])
    xml_nr2 = int(i.attrib["nr2"])
    xml_nr3 = int(i.attrib["nr3"])
  if FFTgrid[0] != xml_nr1 or FFTgrid[1] != xml_nr2 or FFTgrid[2] != xml_nr3:
    print("error write_evc : nr")
    sys.exit()

  for i in root.findall("./output/atomic_structure"):
    xml_alat = float(i.attrib['alat'])
  for i in root.findall("./output/atomic_structure/cell/a1"):
    tmp = i.text
    tmp = tmp.split()
    xml_a1 = np.array(tmp).astype(float)
  for i in root.findall("./output/atomic_structure/cell/a2"):
    tmp = i.text
    tmp = tmp.split()
    xml_a2 = np.array(tmp).astype(float)
  for i in root.findall("./output/atomic_structure/cell/a3"):
    tmp = i.text
    tmp = tmp.split()
    xml_a3 = np.array(tmp).astype(float)
    xml_a1 = xml_a1 / xml_alat
    xml_a2 = xml_a2 / xml_alat
    xml_a3 = xml_a3 / xml_alat
  xml_a = np.zeros((3,3))
  xml_a[0,:] = xml_a1 
  xml_a[1,:] = xml_a2 
  xml_a[2,:] = xml_a3
 
  omega = volume(xml_alat,xml_a1,xml_a2,xml_a3)
  if (celvol - omega) > eps6:
    print("error write_evc : unit cell volume")
    sys.exit()
  
  xdel = 0.0
  for j in range(3):
    for i in range(3):
      xdel = xdel + np.abs(alat* avec[i,j]-xml_alat*xml_a[i,j])
  if xdel > eps6:
    print("error write_evc : reciprocal lattice vectors")
    sys.exit()

#JDY 17OCT2020
  if nspin == 1:
    for i in root.findall("./output/band_structure/nbnd"):
      nb = int(i.text)
    nbgw = nb
    if wfng_nband > 0 and wfng_nband < nb:
      nb = wfng_nband      
  if nspin == 2:
    for i in root.findall("./output/band_structure/nbnd_up"):
      nb_up = int(i.text)
    for i in root.findall("./output/band_structure/nbnd_dw"):
      nb_dw = int(i.text)
    if nb_up != nb_dw:
      print("error write_evc : nb_up neq nb_dw")
      sys.exit()
    nbgw = nb_up
    if wfng_nband > 0 and wfng_nband < nb_up:
      nb = wfng_nband      
          

# JDY 16OCT2020
# int after divide or dvide after int?

  for i in root.findall("./parallel_info/nprocs"):
    nproc = int(i.text)
  if np.mod(ngkmax,nproc) == 0 :
    tmp = ngkmax / nproc
    ngkdist_1 = int(tmp)
  else:
    tmp = ngkmax / nproc + 1
    ngkdist_1 = int(tmp)
  ngkdist_g = ngkdist_1 * nproc
# end

  nk = nrk
  nb = mnband
  ns = nspin
  nr = [FFTgrid[i] for i in range(3)]
  al = alat 

#  ngk_g = np.zeros((nk))
#  k = np.zeros((3,nk))
#  en = np.zeros((nb,nk,ns))
#  oc = np.zeros((nb,nk,ns))
#  gvec = np.zeros((3,ng))
#  igk_buf = np.zeros((ngkdist_g))
#  igk_dist = np.zeros((ngkdist_1,nk))
#  gk_buf = np.zeros((3,ngkdist_g))
#  gk_dist = np.zeros((3,ngkdist_1))
#  if real_or_complex == 1 :
#    wfngr = np.zeros((ngkmax,ns))
#  else:
#    wfngc = np.zeros((ngkmax,ns))
#  wfng_buf = np.zeros((ngkdist_g,ns))
#  wfng_dist = np.zeros((ngkdist_1,nb,ns,nk))

#
# gvec write
# 

  print("write wfc.hdf5")
  if ns == 1:
    for i in range(nk):
      fname = output_dir_name + "wfc"+str(i+1)+".hdf5"
      print(fname)
      with h5py.File(fname,"w") as hh:

        hh.create_dataset('MillerIndices',data=gvecs[:ngk[i],:])
        hh.create_dataset('evc',data=coeffs[:,0,:ngk[i],ns])

      hh.close()

  return


def write_cd(input_file_name,real_or_complex,wfng_nband,output_dir_name):

  print("check inversion")
  check_inversion(real_or_complex,nsym,s,nspin,True,True)
  print("check inversion done")
   

  return

def write_eig(input_file_name,real_or_complex,output_dir_name):

  print("check inversion")
  check_inversion(real_or_complex,nsym,s,nspin,True,True)
  print("check inversion done")

  print("write eig and occ on out.xml")
  root = get_xml_data(output_dir_name)
  xml_ngk = np.zeros((nrk)).astype(int)
  n = 0
  for i in root.findall("./output/band_structure/ks_energies/npw"):
    xml_ngk[n] = int(i.text)
    n += 1
  if ngk.all() != xml_ngk.all():
    print("error write_eig : ngk")
    sys.exit()

  print(ngk)
  print(xml_ngk)

  for i in root.findall("./general_info/job"):
    print(i.text)
    i.text = "(-w-)"
    print(i.text)

  tree.write(encoding = 'utf-8',xml_declaration = True,\
         method = 'xml',short_empty_elements = False)


if __name__ == "__main__":

  pwd = os.getcwd()
  output_dir_name = pwd + "/" + outdir + "/" + prefix + ".save/"
  print(output_dir_name)  


  if wfng_flag  :
    write_evc(wfng_file,real_or_complex,wfng_nband,output_dir_name)

  if rhog_flag :
    write_cd(wfng_file,real_or_complex,wfng_nband,output_dir_name)

  if eig_flag :
    write_eig(eig_dir,real_or_complex,output_dir_name)

