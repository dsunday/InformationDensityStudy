# -*- coding: utf-8 -*-
"""
Functions to be used in analyzing CDSAXS data
"""
import numpy as np

from scipy.interpolate import interp1d
    
    
def ID7MCoordAssign(TPAR,SLD,Trapnumber,Pitch,X1):
    Coord=np.zeros([Trapnumber+1,5,6])
    for T in range (Trapnumber+1):
        if T==0:
            Coord[T,0,0]=0
            Coord[T,1,0]=TPAR[0,0]
            Coord[T,2,0]=TPAR[0,1]
            Coord[T,3,0]=0
            Coord[T,4,0]=SLD[0,0]
        else:
            Coord[T,0,0]=Coord[T-1,0,0]+0.5*(TPAR[T-1,0]-TPAR[T,0])
            Coord[T,1,0]=Coord[T,0,0]+TPAR[T,0]
            Coord[T,2,0]=TPAR[T,1]
            Coord[T,3,0]=0
            Coord[T,4,0]=SLD[T,0]
    
    Coord[:,:,1]=Coord[:,:,0]
    Coord[:,0:2,1]=Coord[:,0:2,1]+2*X1
    Coord[:,:,2]=Coord[:,:,1]
    Coord[:,0:2,2]=Coord[:,0:2,2]+X1
    Coord[:,:,3]=Coord[:,:,2]
    Coord[:,0:2,3]=Coord[:,0:2,3]+X1
    Coord[:,:,4]=Coord[:,:,3]
    Coord[:,0:2,4]=Coord[:,0:2,4]+X1
    Coord[:,:,5]=Coord[:,:,4]
    Coord[:,0:2,5]=Coord[:,0:2,5]+X1
    return (Coord)
         
         
           
                        
def FreeFormTrapezoid(Coord,Qx,Qz,Trapnumber):
    H1 = Coord[0,3]
    H2 = Coord[0,3]
    form=np.zeros([len(Qx[:,1]),len(Qx[1,:])])
    for i in range(int(Trapnumber)):
        H2 = H2+Coord[i,2];
        if i > 0:
            H1 = H1+Coord[i-1,2]
        x1 = Coord[i,0]
        x4 = Coord[i,1]
        x2 = Coord[i+1,0]
        x3 = Coord[i+1,1]
        if x2==x1:
            x2=x2-0.000001
        if x4==x3:
            x4=x4-0.000001
        SL = Coord[i,2]/(x2-x1)
        SR = -Coord[i,2]/(x4-x3)
        
        A1 = (np.exp(1j*Qx*((H1-SR*x4)/SR))/(Qx/SR+Qz))*(np.exp(-1j*H2*(Qx/SR+Qz))-np.exp(-1j*H1*(Qx/SR+Qz)))
        A2 = (np.exp(1j*Qx*((H1-SL*x1)/SL))/(Qx/SL+Qz))*(np.exp(-1j*H2*(Qx/SL+Qz))-np.exp(-1j*H1*(Qx/SL+Qz)))
        form=form+(1j/Qx)*(A1-A2)*Coord[i,4]
    return form
  

def Misfit(Exp,Sim):
    Chi2= abs(np.log(Exp)-np.log(Sim))
    #ms=np.zeros([len(Exp[:,1]),len(Exp[1,:]),2])
    #ms[:,:,0]=Sim
    #ms[:,:,1]=Exp
    #MS= np.nanmin(ms,2)
    #Chi2=np.power((D/MS),2)
    Chi2[np.isnan(Chi2)]=0
    return Chi2
    
def PBA_ID2(TPAR,SPAR,Trapnumber,X1):
     
    SPARLB=SPAR[0:4]*0.8
    SPARUB=SPAR[0:4]*1.2

    FITPAR=TPAR[:,0:2].ravel()
    FITPARLB=FITPAR*0.05
    FITPARUB=FITPAR*20
    FITPAR=np.append(FITPAR,SPAR)
       
    FITPARLB=np.append(FITPARLB,SPARLB)
    
    FITPARUB=np.append(FITPARUB,SPARUB)
    
    FITPAR=np.append(FITPAR,X1)
       
    FITPARLB=np.append(FITPARLB,(X1-5))
    
    FITPARUB=np.append(FITPARUB,(X1+5))
    
    return (FITPAR,FITPARLB,FITPARUB)
    
   


        
