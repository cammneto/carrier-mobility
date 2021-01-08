#!/bin/python
import numpy as np
from scipy import constants as spc


#np.set_printoptions(precision=4)
#system  = 'os-c3'
#log = 'inputs/'+system+'/dis-'+system+'.log'
#
#log1 = 'inputs/'+system+'/dis-'+system+'-ani'+'.log'
#
#print('\n','Reading', log)
#print('\n','Reading', log1)

##############################################################################################################################################################
c=spc.c*2*spc.pi*100
kb=spc.value(u'Boltzmann constant in eV/K')
hbar=spc.value(u'Planck constant over 2 pi in eV s')
bohrme = (spc.value('Bohr radius')*(np.sqrt(spc.m_e)))#5.2917721067e-11
##############################################################################################################################################################
def readFile(dispFile):
    infile = open(dispFile, "r")
    return infile

def freq(dispFile):
    f=readFile(dispFile)
    freqs=[]
    for line in f:
        if "Frequencies --" in line:
            freqs.extend([float(i) for i in line.split()[-3:]])
    W  = np.asarray(freqs[:int(len(freqs)/3)])*c
    return W

def shift(dispFile):
    f=readFile(dispFile)
    masses=[]
    shift=[]
    for line in f:
        if line.find("Red. masses --") != -1:
            masses.extend([float(i) for i in line.split()[-3:]])
        if line.find("Shift") != -1:
            for line in f:
                if line.find("1") != -1:
                    for line in f:
                        if line.find('0') != -1:
                            shift.extend(line.split()[-1:])
        mu = np.asarray(masses[:int(len(masses)/3)])
        shift = shift[:int(len(mu))]
    for i in range(len(shift)):
        shift[i] = float(shift[i].replace('D',"E"))
    dQ = np.asarray(shift)*bohrme
    return dQ


#print(shift(log1))
##############################################################################################################################################################
#W = freq(log)
#dQ = shift(log)
#Si=((dQ**2)*W)/(2*spc.hbar)
#weff=np.sum(Si*W)/np.sum(Si)
##print('\n',' Effective Frequency (\u03C9eff)','\n','\u03C9eff =', round(weff/c,4),'cm-\N{SUPERSCRIPT ONE}')
#dQ1 = shift(log1)
#W1 = freq(log1)
#Si1=((dQ1**2)*W1)/(2*spc.hbar)
#weff1=np.sum(Si1*W1)/np.sum(Si1)
#print('\n',' Effective Frequency (\u03C9eff)','\n','\u03C9eff =', round(weff1/c,4),'cm-\N{SUPERSCRIPT ONE}')
##############################################################################################################################################################
#li0=np.sum(Si*hbar*W)
#print(li0)
#li1=np.sum(Si1*hbar*W1)
#print(li1)
#li=li0+li1
#print('\n',' Vibronic Internal Reorganization Energy (\u03BBv) ','\n','\u03BBv =',round(li,4),'eV')
