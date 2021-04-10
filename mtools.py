#!/bin/env python3
# -*- coding: utf-8 -*-

# Author: Carlos A M de Melo Neto
# Contact: cammneto@gmail.com

import numpy as np
import pandas as pd
import subprocess
np.set_printoptions(precision=5)

def importXYZGeom(fileName):
    xyzGeom = pd.DataFrame(columns=['atomName', 'x', 'y', 'z'])
    with open(fileName,'r') as openFile:
        lines = openFile.readlines()[2:]
    for line in lines:
        line = line.split()
        if line == []:
            continue
        else:
            xyzGeom = xyzGeom.append({'atomName':line[0], 'x':float(line[1]), 'y':float(line[2]), 'z':float(line[3])} , ignore_index=True)
    return xyzGeom

def centerOFmass(mol):
    atomSym = mol['atomName'].tolist()
    weights = {'H':1,'C':12} #Atomic symbol and number
    totalMass = np.sum([weights[i] for i in atomSym])
    centerOFmassx = np.sum(np.array([weights[i] for i in atomSym])*np.array([mol[['x']][i] for i in mol[['x']]]))/totalMass
    centerOFmassy = np.sum(np.array([weights[i] for i in atomSym])*np.array([mol[['y']][i] for i in mol[['y']]]))/totalMass
    centerOFmassz = np.sum(np.array([weights[i] for i in atomSym])*np.array([mol[['z']][i] for i in mol[['z']]]))/totalMass
    centerOFmass = np.array([centerOFmassx,centerOFmassy,centerOFmassz])
    #print(centerOFmass)
    return centerOFmass

def centerDistance_vec(fileName):
    mol1 = importXYZGeom(fileName).loc[:int((len(importXYZGeom(fileName))-1)/2)]
    mol1.index = mol1.index + 1
    mol2 = importXYZGeom(fileName).loc[int((len(importXYZGeom(fileName)))/2):]
    mol2.index = mol2.index + 1
    R = centerOFmass(mol1) - centerOFmass(mol2)
    return R

def CATNIP(pun_file_1,orb_ty_1,pun_file_2,orb_ty_2,pun_file):
    print('\n', 'Calculating electronic coupling...')
    args = ("./calc_J", "-p_1", pun_file_1, '-orb_ty_1', orb_ty_1, '-p_2', pun_file_2, '-orb_ty_2', orb_ty_2, '-p_P', pun_file)
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    output = popen.stdout.read()
    outline = output.splitlines()
    #for i in outline:
    #    print(i)
    J_eff = output.split()[-3:]
    print('\n','Electronic Coupling (H_ad)','\n','H_ad =',(float(J_eff[1])), 'eV')
    return J_eff
