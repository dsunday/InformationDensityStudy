# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 14:14:39 2017

@author: dfs1
"""

import numpy as np
import CDSAXSfunctions as CD
import math as math
import CDplot as CDp

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
    
    
#Sample Parameter Definition

Trapnumber = 3

DW = 1.3
I0 = 0.0025
Bk =1
SLD1 = 1;
TPAR=np.zeros([Trapnumber+1,2])
SLD=np.zeros([Trapnumber+1,1])
SPAR=np.zeros(3)
SPAR[0]=DW; SPAR[1]=I0; SPAR[2]=Bk;


TPAR[0,0]=25; TPAR[0,1]=15; SLD[0,0]=SLD1;
TPAR[1,0]=21; TPAR[1,1]=123;SLD[1,0]=SLD1;
TPAR[2,0]=20; TPAR[2,1]=12; SLD[2,0]=SLD1;
TPAR[3,0]=10; TPAR[3,1]=0; SLD[3,0]=SLD1;

Coord=CD.ID1CoordAssign(TPAR,SLD,Trapnumber,Pitch)
#CDp.plotID1(Coord,Trapnumber,Pitch)
(FITPAR,FITPARLB,FITPARUB)=CD.PBA_ID1(TPAR,SPAR,Trapnumber)

def SimInt_ID1(FITPAR):
    TPARs=np.zeros([Trapnumber+1,2])
    TPARs[:,0:2]=np.reshape(FITPAR[0:(Trapnumber+1)*2],(Trapnumber+1,2))
    SPAR=FITPAR[Trapnumber*2+2:Trapnumber*2+5]
    (Coord)= CD.ID1CoordAssign(TPAR,SLD,Trapnumber,Pitch)
    F1 = CD.FreeFormTrapezoid(Coord[:,:,0],Qx,Qz,Trapnumber) 
    M=np.power(np.exp(-1*(np.power(Qx,2)+np.power(Qz,2))*np.power(SPAR[0],2)),0.5)
    Formfactor=F1*M
    Formfactor=abs(Formfactor)
    SimInt = np.power(Formfactor,2)*SPAR[1]+SPAR[2]
    return (SimInt,Formfactor)



R = np.random.normal(0, 0.225, [len(Qx[:,0]),len(Qx[0,:])])
(DummyIntensity,Amplitude)=SimInt_ID1(FITPAR)
PreInt=abs(Amplitude)
PreInt=np.power(PreInt,2)


M=np.amax(PreInt)
I0=(20000)/M # Scales the intensity so that the max is 20000
SPAR[1]=I0
Intensity=PreInt*I0+Bk

(FITPAR,FITPARLB,FITPARUB)=CD.PBA_ID1(TPAR,SPAR,Trapnumber) #regenerates FITPAR with proper intensity scaling

N=(1/(np.power(Intensity,0.5)))*Intensity # Generates  noise
N=N*R
Intensity2=Intensity+R; # Applies Noise

C=CD.Misfit(Intensity,Intensity2)
Chi2=np.sum((C),0)

MCPAR=np.zeros([7])
MCPAR[0] = 12 # Chainnumber
MCPAR[1] = len(FITPAR)
MCPAR[2] = 10 #stepnumber
MCPAR[3] = 1 #randomchains
MCPAR[4] = 1 # Resampleinterval
MCPAR[5] = 100 # stepbase
MCPAR[6] = 100 # steplength





