#!/bin/python
import sys
import numpy as np
import subprocess
from scipy import constants as spc
from math import factorial as fct
from mtools import centerDistance as R_vec
from vibration import *

np.set_printoptions(precision=4)
system  = 'ss-c3'
dim_ty  = 'q'
orb_ty1 = 'LUMO'
orb_ty2 = 'LUMO'
log = 'inputs/'+system+'/dis-'+system+'.log'
log1 = 'inputs/'+system+'/dis-'+system+'-ani.log'
xyzfile = 'inputs/'+system+'/'+system+'-'+dim_ty+'.xyz'
print('\n',system, 'dimer', dim_ty,'\n','Reading', log,'\n','Reading', log1,'\n','Reading', xyzfile)
c = spc.c*2*spc.pi*100
kb = spc.value(u'Boltzmann constant in eV/K')
hbar = spc.value(u'Planck constant over 2 pi in eV s')
bohrme = (spc.value('Bohr radius')*(np.sqrt(spc.m_e)))#5.2917721067e-11
T  = 296
G0 = 0
ls = np.array([0.001,0.0015,0.0022,0.0031,0.0034,0.0037,0.0043,0.0049,0.0051,0.0077,0.0078])
ls = np.mean(ls)
W  = freq(log)
dQ = shift(log)
Si=((dQ**2)*W)/(2*spc.hbar)
weff=np.sum(Si*W)/np.sum(Si)
W1 = freq(log1)
dQ1 = shift(log1)
Si1=((dQ1**2)*W1)/(2*spc.hbar)
print('\n','Number of Normal Modes (n) ','\n' ,'n =',n,'\n','Frequencies (\u03C9i) ','\n','\u03C9i =',W,'\n','Reduced Masses (\u03BCi) ','\n','\u03BCi =', mu)
print('\n','Shift Vector (dQ) times sqrt(\u03BCi) ','\n','dQ*sqrt(\u03BC) =',dQ,'\n',' Effective Frequency (\u03C9eff)','\n','\u03C9eff =', round(weff/c,4),'cm-\N{SUPERSCRIPT ONE}')
lv0=np.sum(Si*hbar*W)
lv1=np.sum(Si1*hbar*W1)
lv=lv0+lv1
print('\n',' Vibronic Internal Reorganization Energy (\u03BBv) ','\n','\u03BBv =',round(lv,4),'eV')
print('\n', 'Calculating electronic coupling...')
args = ("./calc_J", "-p_1", 'inputs/'+system+'/'+system+'-'+dim_ty+"m1.pun", '-orb_ty_1', orb_ty1, '-p_2', 'inputs/'+system+'/'+system+'-'+dim_ty+'m2.pun', '-orb_ty_2', orb_ty2, '-p_P', 'inputs/'+system+'/'+system+'-'+dim_ty+'.pun')
popen = subprocess.Popen(args, stdout=subprocess.PIPE)
popen.wait()
output = popen.stdout.read()
Jeff = output.split()[-3:]
outline = output.splitlines()
print('\n','Electronic Coupling (H_ad)','\n','H_ad =',(float(Jeff[1]))*1000, 'eV')
Had=float(Jeff[1])**2
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
Kmlj=C*Had*S1*SOMA
lamb=lv+ls
Cm = 2*spc.pi/(hbar*np.sqrt(4*spc.pi*lamb*kb*T))
Smn=(lamb-G0)**2
Smd=(4*lamb*kb*T)
Sm=np.exp(-Smn/Smd)
Kmarcus=Cm*Had*Sm
ket=np.array([Kmarcus,Kmlj])
print('\n','Marcus and MJL rates respectivelly','\n',ket,'s-\N{SUPERSCRIPT ONE}')
R = R_vec(xyzfile)
d = np.sqrt((R[0]**2)+(R[1]**2)+(R[2]**2))/10**8
print('\n','Site Distance =',round(d*10**8,4),'Angstrons')
mob=(spc.e*ket*(d**2))/(2*spc.k*T)
print('\n','Marcus and MJL mobilities respectivelly','\n',mob,'cm\u00b2 V-\N{SUPERSCRIPT ONE} s-\N{SUPERSCRIPT ONE}', 'for T =', T, 'K')
