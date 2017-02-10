# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 14:14:39 2017

@author: dfs1
"""

import numpy as np
import CDSAXSfunctions as CD
import math as math
import CDplot as CDp
from multiprocessing import Pool
import time
Param=np.loadtxt('SingleLineParameters.txt')

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
    
def MCMCInit_ID1(FITPAR,FITPARLB,FITPARUB,MCPAR):
    
    MCMCInit=np.zeros([int(MCPAR[0]),int(MCPAR[1])+1])
    for i in range(int(MCPAR[0])):
        if i <MCPAR[3]: #reversed from matlab code assigns all chains below randomnumber as random chains
            for c in range(int(MCPAR[1])):
                MCMCInit[i,c]=FITPARLB[c]+(FITPARUB[c]-FITPARLB[c])*np.random.random_sample()
                (SimInt,Amplitude)=SimInt_ID1(MCMCInit[i,:])
            C=np.sum(CD.Misfit(Intensity2,SimInt))
            print(C)
            MCMCInit[i,int(MCPAR[1])]=C
            
        else:
            MCMCInit[i,0:int(MCPAR[1])]=FITPAR
            (SimInt,Amplitude)=SimInt_ID1(MCMCInit[i,:])
            C=np.sum(CD.Misfit(Intensity2,SimInt))
            MCMCInit[i,int(MCPAR[1])]=C
            
           
    return MCMCInit


    
    
def MCMC_ID1(MCMC_List):
    
    MCMCInit=MCMC_List
    
    L = int(MCPAR[1])
    Stepnumber= int(MCPAR[2])
        
    SampledMatrix=np.zeros([Stepnumber,L+1]) 
    SampledMatrix[0,:]=MCMCInit
    Move = np.zeros([L+1])
    
    ChiPrior = MCMCInit[L]
    for step in np.arange(1,Stepnumber,1): 
        Temp = SampledMatrix[step-1,:].copy()
        for p in range(L-1):
            StepControl = MCPAR[5]+MCPAR[6]*np.random.random_sample()
            Move[p] = (FITPARUB[p]-FITPARLB[p])/StepControl*(np.random.random_sample()-0.5) # need out of bounds check
            Temp[p]=Temp[p]+Move[p]
            if Temp[p] < FITPARLB[p]:
                Temp[p]=FITPARLB[p]+(FITPARUB[p]-FITPARLB[p])/1000
            elif Temp[p] > FITPARUB[p]:
                Temp[p]=FITPARUB[p]-(FITPARUB[p]-FITPARLB[p])/1000
        (SimPost,AmpPost)=SimInt_ID1(Temp)
        ChiPost=np.sum(CD.Misfit(Intensity2,SimPost))
        if ChiPost < ChiPrior:
            SampledMatrix[step,0:L]=Temp[0:L]
            SampledMatrix[step,L]=ChiPost
            ChiPrior=ChiPost
            
        else:
            MoveProb = np.exp(-0.5*np.power(ChiPost-ChiPrior,2))
            if np.random.random_sample() < MoveProb:
                SampledMatrix[step,0:L]=Temp[0:L]
                SampledMatrix[step,L]=ChiPost
                ChiPrior=ChiPost
            else:
                SampledMatrix[step,:]=SampledMatrix[step-1,:]
    
    ReSampledMatrix=np.zeros([int(MCPAR[2])/int(MCPAR[4]),len(SampledMatrix[1,:])])

    c=-1
    for i in np.arange(0,len(SampledMatrix[:,1]),MCPAR[4]):
        c=c+1
        ReSampledMatrix[c,:]=SampledMatrix[i,:]
    return (ReSampledMatrix)

for SampleNumber in range(1):
    #QxQz map definition
    Angle = Param[SampleNumber,3]
    Pitch = Param[SampleNumber,2]
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
    # DW, I0 and Bk paramters are defined internally, structural pramters are defined externally
    DW = 1.3
    I0 = 0.0025
    Bk =1
    SLD1 = 1;
    TPAR=np.zeros([Trapnumber+1,2])
    SLD=np.zeros([Trapnumber+1,1])
    SPAR=np.zeros(3)
    SPAR[0]=DW; SPAR[1]=I0; SPAR[2]=Bk;
    
    W=Param[SampleNumber,0]
    H=Param[SampleNumber,1]
    WidthLoad='W'+str(int(W))+'.txt'
    HeightLoad='H'+str(int(H))+'.txt'
    Width=np.loadtxt(WidthLoad)
    Height=np.loadtxt(HeightLoad)
    TPAR[:,0]=Width
    TPAR[:,1]=Height
    SLD[0,0]=SLD1;
    SLD[1,0]=SLD1;
    SLD[2,0]=SLD1;
    SLD[3,0]=SLD1;

    Coord=CD.ID1CoordAssign(TPAR,SLD,Trapnumber,Pitch)
    #CDp.plotID1(Coord,Trapnumber,Pitch)
    (FITPAR,FITPARLB,FITPARUB)=CD.PBA_ID1(TPAR,SPAR,Trapnumber)

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
    Chi2=np.sum((C))

    MCPAR=np.zeros([7])
    MCPAR[0] = 1 # Chainnumber
    MCPAR[1] = len(FITPAR)
    MCPAR[2] = 10 #stepnumber
    MCPAR[3] = 0 #randomchains
    MCPAR[4] = 1 # Resampleinterval
    MCPAR[5] = 14 # stepbase
    MCPAR[6] = 14 # steplength 
  
    
    MCMCInitial=MCMCInit_ID1(FITPAR,FITPARLB,FITPARUB,MCPAR)
    Acceptprob=0;
    while Acceptprob < 0.3 or Acceptprob > 0.4:
        L = int(MCPAR[1])
        Stepnumber= int(MCPAR[2])
        
        SampledMatrix=np.zeros([Stepnumber,L+1]) 
        SampledMatrix[0,:]=MCMCInitial[0,:]
        Move = np.zeros([L+1])
    
        ChiPrior = MCMCInitial[0,L]
        for step in np.arange(1,Stepnumber,1): 
            Temp = SampledMatrix[step-1,:].copy()
            for p in range(L-1):
                StepControl = MCPAR[5]+MCPAR[6]*np.random.random_sample()
                Move[p] = (FITPARUB[p]-FITPARLB[p])/StepControl*(np.random.random_sample()-0.5) # need out of bounds check
                Temp[p]=Temp[p]+Move[p]
                if Temp[p] < FITPARLB[p]:
                    Temp[p]=FITPARLB[p]+(FITPARUB[p]-FITPARLB[p])/1000
                elif Temp[p] > FITPARUB[p]:
                    Temp[p]=FITPARUB[p]-(FITPARUB[p]-FITPARLB[p])/1000
            (SimPost,AmpPost)=SimInt_ID1(Temp)
            ChiPost=np.sum(CD.Misfit(Intensity2,SimPost))
            if ChiPost < ChiPrior:
                SampledMatrix[step,0:L]=Temp[0:L]
                SampledMatrix[step,L]=ChiPost
                ChiPrior=ChiPost
            
            else:
                MoveProb = np.exp(-0.5*np.power(ChiPost-ChiPrior,2))
                if np.random.random_sample() < MoveProb:
                    SampledMatrix[step,0:L]=Temp[0:L]
                    SampledMatrix[step,L]=ChiPost
                    ChiPrior=ChiPost
                else:
                    SampledMatrix[step,:]=SampledMatrix[step-1,:]
        AcceptanceNumber=0
        Acceptancetotal=len(SampledMatrix[:,1])

        for i in np.arange(1,len(SampledMatrix[:,1]),1):
            if SampledMatrix[i,0] != SampledMatrix[i-1,0]:
                AcceptanceNumber=AcceptanceNumber+1
        Acceptprob=AcceptanceNumber/Acceptancetotal
        
        if Acceptprob < 0.3:
            MCPAR[5]=MCPAR[5]+2
            MCPAR[6]=MCPAR[6]+2
        if Acceptprob > 0.4:
            MCPAR[5]=MCPAR[5]-2
            MCPAR[6]=MCPAR[6]-2
        
    start_time = time.perf_counter()
    
    
    MCMC_List=[0]*int(MCPAR[0])
    for i in range(int(MCPAR[0])):
        MCMC_List[i]=MCMCInitial[i,:]
    if __name__ =='__main__':  
        pool = Pool(processes=12)
              
        F=pool.map(MCMC_ID1,MCMC_List)
        F=tuple(F)
        Savename='P'+str(int(Pitch))+'_'+'W'+str(int(W))+'_'+'H'+str(int(H))+'_'+'A'+str(int(Angle))
        np.save(Savename,F) 
        end_time=time.perf_counter()   
        print(end_time-start_time)    
        
        ReSampledMatrix=F[0]
        Xi = np.zeros([102,2,len(ReSampledMatrix[:,0])])
        Yi = np.zeros([102,1,len(ReSampledMatrix[:,0])])
        PopWidth= np.zeros([len(ReSampledMatrix[:,0])])
        PopHeight=np.zeros([len(ReSampledMatrix[:,0])])
        
        for PopNumber in range(len(ReSampledMatrix[:,0])):
            
            
            TPARU=np.zeros([Trapnumber+1,2])
            TPARU[:,0:2]=np.reshape(ReSampledMatrix[PopNumber,0:(Trapnumber+1)*2],(Trapnumber+1,2))
           
            (CoordUnc)= CD.ID1CoordAssign(TPARU,SLD,Trapnumber,Pitch)
            
            TrapHeight=np.zeros([Trapnumber+1,len(ReSampledMatrix[:,0])])
            for i in np.arange(1,Trapnumber+1,1):
                TrapHeight[i,PopNumber]=TrapHeight[i-1,PopNumber]+TPARU[i-1,1]
            
            for LineNumber in range(1):
                c=-1
                EffTrapnumber=0
                #LeftSide
                X1L = CoordUnc[EffTrapnumber,0,LineNumber]
                X2L = CoordUnc[EffTrapnumber+1,0,LineNumber]
                X1R = CoordUnc[EffTrapnumber,1,LineNumber]
                X2R = CoordUnc[EffTrapnumber+1,1,LineNumber]
                Y1 = TrapHeight[EffTrapnumber,PopNumber]
                Y2 = TrapHeight[EffTrapnumber+1,PopNumber]
                
                for Disc in np.arange(0,TrapHeight[Trapnumber,PopNumber]+0.0000001,TrapHeight[Trapnumber,PopNumber]/101):
                    c=c+1
                    if Disc > Y2:
                        EffTrapnumber=EffTrapnumber+1
                        
                        X1L = CoordUnc[EffTrapnumber,0,LineNumber]
                        X2L = CoordUnc[EffTrapnumber+1,0,LineNumber]
                        X1R = CoordUnc[EffTrapnumber,1,LineNumber]
                        X2R = CoordUnc[EffTrapnumber+1,1,LineNumber]
                        Y1 = TrapHeight[EffTrapnumber,PopNumber]
                        Y2 = TrapHeight[EffTrapnumber+1,PopNumber]
                    ML = (Y2-Y1)/(X2L-X1L)
                    MR=  (Y2-Y1)/(X2R-X1R)
                    BL = Y1-ML*X1L
                    BR = Y1-MR*X1R
                    Xi[c,0,PopNumber]=(Disc-BL)/ML
                    Xi[c,1,PopNumber]=(Disc-BR)/MR
                    Yi[c,0,PopNumber]=Disc
            Xi[:,:,PopNumber]=Xi[:,:,PopNumber]-(Xi[0,1,PopNumber]-Xi[0,0,PopNumber])/2
            PopWidth[PopNumber]=Xi[50,1,PopNumber]-Xi[50,0,PopNumber]
            PopHeight[PopNumber]=Yi[101,0,PopNumber]
