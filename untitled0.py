# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 14:14:39 2017

@author: dfs1
"""

import numpy as np
import CDSAXSfunctions as CD
import math as math

#QxQz map definition
Angle = 121
Pitch = 128
Wavelength = 0.7
Spacing = 2*np.pi/Pitch
QxI = np.arange(Spacing,1.26,Spacing)
DetectorAngle=np.zeros([len(QxI),1])
UQz=np.zeros([len(QxI),1])
i=0
R=math.radians(60)
while i < len(QxI):
    DetectorAngle[i] = math.asin(QxI[i]*Wavelength/(4*np.pi))
    UQz[i]=2*QxI[i]*math.sin(R-DetectorAngle[i]/2)
    i=i+1
    
for i in range(len(UQz)):
    if UQz[i,0] > 1.1:
        UQz[i,0]=1.1-(QxI[i]-0.65);
        
Qx = np.zeros([Angle,len(UQz)])
Qz = np.zeros([Angle,len(UQz)])

for i in range(len(UQz)):
    Qx[:,i]=QxI[i]
    Qz[:,i]=np.arange(-1*UQz[i,0],UQz[i,0]+0.0001,2*UQz[i,0]/(Angle-1))
            