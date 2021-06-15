#!/bin/env python3
# -*- coding: utf-8 -*-

# Author: Carlos A M de Melo Neto
# Contact: cammneto@gmail.com

import sys
import numpy as np
from scipy import constants as spc
from math import factorial as fct
from cmtools import *

np.set_printoptions(precision=4)

### carrier mobility input files ###
xyzfile = 'ethylene-0-dimer.xyz'                ###  Path of geometry file .xyz
log     =  'ethylene-neutral-displacement.log'  ###  Path of neutral to excited displacement Gaussian .log  file
log1    =  'ethylene-anion-displacement.log'    ###  Path of excited to neutral displacement Gaussian .log  file
### CATNIP input files ### For more details check CATNIP page
pun_file   = 'ethylene-0-dimer.pun'             ###  Path of Dimer pun file for CATNIP
pun_file_1 = 'ethylene-0-m1.pun'                ###  Path of 1st monomer pun file for CATNIP
pun_file_2 = 'ethylene-0-m2.pun'                ###  Path of 2nd monomer pun file for CATNIP

print('\n','Open Files','\n',log,'\n',log1,'\n',xyzfile,'\n',pun_file,'\n',pun_file_1,'\n',pun_file_2)

### Choose orbital for coupling calculation with CATNIP
orb_ty_1 = 'LUMO'                                           ###  1st orbital chosen for coupling (ex: LUMO, HOMO...)
orb_ty_2 = 'LUMO'                                           ###  2nd orbital chosen for coupling (ex: LUMO, HOMO...)

### Constants ###
c = spc.c*2*spc.pi*100
kb = spc.value(u'Boltzmann constant in eV/K')
hbar = spc.value(u'Planck constant over 2 pi in eV s')
bohrme = spc.value('Bohr radius')*(np.sqrt(spc.m_e))

### Defined Variables ###
T  = 298    ### Temperature in Kelvin
G0 = 0.00   ### Gibbs free energy
ls = 0.04   ### Type here your sample external reorganization energy in eV.

### Computed Variables
W  = freq(log)
dQ = shift(log)
Si=((dQ**2)*W)/(2*spc.hbar)
weff=np.sum(Si*W)/np.sum(Si)
W1 = freq(log1)
dQ1 = shift(log1)
Si1=((dQ1**2)*W1)/(2*spc.hbar)
n = len(W)
print('\n','Frequencies (\u03C9i) ','\n','\u03C9i =',W)
print('\n','Shift Vector (dQ) times sqrt(\u03BCi) ','\n','dQ*sqrt(\u03BC) =',dQ)
print('\n','Effective Frequency (\u03C9_eff)','\n','\u03C9_eff =', round(weff/c,4),'cm-\N{SUPERSCRIPT ONE}')
lv0=np.sum(Si*hbar*W)
lv1=np.sum(Si1*hbar*W1)
lv=lv0+lv1
print('\n',' Vibronic Internal Reorganization Energy (\u03BBv) ','\n','\u03BB_v =',round(lv,4),'eV')
J_eff=CATNIP(pun_file_1,orb_ty_1,pun_file_2,orb_ty_2,pun_file)### Calling CATNIP to compute transfer integral (J_eff) between orbitals defined in begining of this file.

### MLJ calculation
Had=float(J_eff[1])**2
SOMA = 0
C = spc.pi/(hbar*np.sqrt(spc.pi*kb*T*ls))
S = lv/(hbar*weff)
S1 = np.exp(-S)
for ni in range(len(W)):
    S2 = (S**ni)/fct(ni)
    S3n = (-G0 + ls + ni*hbar*weff)**2
    S3d = 4*ls*kb*T
    S3 = np.exp(-S3n/S3d)
    SOMA += S2*S3
K_mlj=C*Had*S1*SOMA
### Semi-Classical Marcus (SCM) calculation
lamb=lv+ls
Cm = 2*spc.pi/(hbar*np.sqrt(4*spc.pi*lamb*kb*T))
Smn=(lamb-G0)**2
Smd=(4*lamb*kb*T)
Sm=np.exp(-Smn/Smd)
K_scm=Cm*Had*Sm
### Array with transfer rates
ket=np.array([K_scm,K_mlj])
print('\n','SCM and MJL rates respectivelly','\n',ket,'s-\N{SUPERSCRIPT ONE}')
### Carrier Mobility
R = centerDistance_vec(xyzfile)
d = np.sqrt((R[0]**2)+(R[1]**2)+(R[2]**2))/10**8
print('\n','Site Distance =',round(d*10**8,4),'Angstrons')
mob=(spc.e*ket*(d**2))/(2*spc.k*T)
print('\n','SCM and MJL mobilities respectivelly','\n',mob,'cm\u00b2 V-\N{SUPERSCRIPT ONE} s-\N{SUPERSCRIPT ONE}', 'for T =', T, 'K')
