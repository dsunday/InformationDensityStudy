# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 14:14:39 2017

@author: dfs1
"""

import numpy as np
import CDSAXSfunctions as CD
import math as math
from multiprocessing import Pool
import time
import pickle
Param=np.loadtxt('3Angle_W12P32H25.txt')

def SimInt_ID1(FITPAR):
    TPARs=np.zeros([Trapnumber+1,2])
    TPARs[:,0:2]=np.reshape(FITPAR[0:(Trapnumber+1)*2],(Trapnumber+1,2))
    SPAR=FITPAR[Trapnumber*2+2:Trapnumber*2+5]
    (Coord)= CD.ID1CoordAssign(TPARs,SLD,Trapnumber,Pitch)
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
        
    SampledMatrixI=np.zeros([Stepnumber,L+1]) 
    SampledMatrixI[0,:]=MCMCInit
    
    
    Move = np.zeros([L+1])
    
    ChiPrior = MCMCInit[L]
    for step in np.arange(1,Stepnumber,1): 
        Temp = SampledMatrixI[step-1,:].copy()
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
            SampledMatrixI[step,0:L]=Temp[0:L]
            SampledMatrixI[step,L]=ChiPost
            ChiPrior=ChiPost
            
        else:
            MoveProb = np.exp(-0.5*np.power(ChiPost-ChiPrior,2))
            if np.random.random_sample() < MoveProb:
                SampledMatrixI[step,0:L]=Temp[0:L]
                SampledMatrixI[step,L]=ChiPost
                ChiPrior=ChiPost
            else:
                SampledMatrixI[step,:]=SampledMatrixI[step-1,:]
    
    ReSampledMatrixI=np.zeros([int(MCPAR[2])/int(MCPAR[4]),len(SampledMatrixI[1,:])])

    c=-1
    for i in np.arange(0,len(SampledMatrixI[:,1]),MCPAR[4]):
        c=c+1
        ReSampledMatrixI[c,:]=SampledMatrixI[i,:]
    (UNCT_Param)=Uncertainty1T(ReSampledMatrixI)
    return (UNCT_Param) #ReSampledMatrixI
    
    
def Uncertainty1T(ReSampledMatrix):
        Xi = np.zeros([101,2,len(ReSampledMatrix[:,0])])
        Yi = np.zeros([101,1,len(ReSampledMatrix[:,0])])
        PopWidth= np.zeros([len(ReSampledMatrix[:,0])])
        PopHeight=np.zeros([len(ReSampledMatrix[:,0])])
        TrapHeight=np.zeros([Trapnumber+1,len(ReSampledMatrix[:,0])])
        for PopNumber in range(len(ReSampledMatrix[:,0])):
            
            
            TPARU=np.zeros([Trapnumber+1,2])
            TPARU[:,0:2]=np.reshape(ReSampledMatrix[PopNumber,0:(Trapnumber+1)*2],(Trapnumber+1,2))
           
            (CoordUnc)= CD.ID1CoordAssign(TPARU,SLD,Trapnumber,Pitch)
            
            
            for i in np.arange(1,Trapnumber+1,1):
                TrapHeight[i,PopNumber]=TrapHeight[i-1,PopNumber]+TPARU[i-1,1]
            
            for LineNumber in range(1):
                
                EffTrapnumber=0
                #LeftSide
                X1L = CoordUnc[EffTrapnumber,0,LineNumber]
                X2L = CoordUnc[EffTrapnumber+1,0,LineNumber]
                X1R = CoordUnc[EffTrapnumber,1,LineNumber]
                X2R = CoordUnc[EffTrapnumber+1,1,LineNumber]
                Y1 = TrapHeight[EffTrapnumber,PopNumber]
                Y2 = TrapHeight[EffTrapnumber+1,PopNumber]
                Disc=0
                for c in np.arange(0,101,1):
                    
                    if Disc > Y2 and Disc <TrapHeight[Trapnumber,PopNumber]:
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
                    Disc=Disc+TrapHeight[Trapnumber,PopNumber]/100
            Xi[:,:,PopNumber]=Xi[:,:,PopNumber]-(Xi[0,1,PopNumber]-Xi[0,0,PopNumber])/2
            PopWidth[PopNumber]=Xi[50,1,PopNumber]-Xi[50,0,PopNumber]
            PopHeight[PopNumber]=Yi[100,0,PopNumber]
            
        S=np.std(Xi,2)*1.96
        Center = np.average(Xi,2)
        Sy=np.std(Yi,2)*1.96
        YC=np.average(Yi,2)
        OuterEdge=Center
        YInner=YC-Sy
        YOuter=YC+Sy
        OuterEdge[:,0]=OuterEdge[:,0]-S[:,0]
        OuterEdge[:,1]=OuterEdge[:,1]+S[:,1]

        InnerEdge=Center
        InnerEdge[:,0]=InnerEdge[:,0]+S[:,0]
        InnerEdge[:,1]=InnerEdge[:,1]-S[:,1]

        LinePlot=np.zeros([2*101,2])
        InnerPlot=np.zeros([2*101,2])
        OuterPlot=np.zeros([2*101,2])
        
        LinePlot[0:101,0]=Center[:,0]
        LinePlot[101:202,0]=np.flipud(Center[:,1])
        LinePlot[0:101,1]=YC[:,0]
        LinePlot[101:202,1]=np.flipud(YC[:,0])
        
        InnerPlot[0:101,0]=InnerEdge[:,0]
        InnerPlot[101:202,0]=np.flipud(InnerEdge[:,1])
        InnerPlot[0:101,1]=YInner[:,0]
        InnerPlot[101:202,1]=np.flipud(YInner[:,0])
        
        
        OuterPlot[0:101,0]=OuterEdge[:,0]
        OuterPlot[101:202,0]=np.flipud(OuterEdge[:,1])
        OuterPlot[0:101,1]=YOuter[:,0]
        OuterPlot[101:202,1]=np.flipud(YOuter[:,0])
        
        vi=np.zeros([100,1])
        vo=np.zeros([100,1])
        vc=np.zeros([100,1])
        for l in np.arange(1,1+0.0001,2):
            for h in np.arange(1,101,1):
                vi[h-1,l/2]=0.5*(YInner[h,l-1]-YInner[h-1,l-1])*((InnerEdge[h-1,l]-InnerEdge[h-1,l-1])+(InnerEdge[h,l]-InnerEdge[h,l-1]))
                vo[h-1,l/2]=0.5*(YOuter[h,l-1]-YOuter[h-1,l-1])*((OuterEdge[h-1,l]-OuterEdge[h-1,l-1])+(OuterEdge[h,l]-OuterEdge[h,l-1]))
                vc[h-1,l/2]=0.5*(YC[h,l-1]-YC[h-1,l-1])*((Center[h-1,l]-Center[h-1,l-1])+(Center[h,l]-Center[h,l-1]))
        vd=vo-vi
        vt=np.sum(vd)
        vct=np.sum(vc)
        WidthAvg=np.average(PopWidth)
        HeightAvg=np.average(PopHeight)
        WidthStd=np.std(PopWidth)
        HeightStd=np.std(PopHeight)
        
        viLower=np.zeros([10,1])
        voLower=np.zeros([10,1])
        vcLower=np.zeros([10,1])
        for l in np.arange(1,1+0.0001,2):
            for h in np.arange(1,11,1):
                viLower[h-1,l/2]=0.5*(YInner[h,l-1]-YInner[h-1,l-1])*((InnerEdge[h-1,l]-InnerEdge[h-1,l-1])+(InnerEdge[h,l]-InnerEdge[h,l-1]))
                voLower[h-1,l/2]=0.5*(YOuter[h,l-1]-YOuter[h-1,l-1])*((OuterEdge[h-1,l]-OuterEdge[h-1,l-1])+(OuterEdge[h,l]-OuterEdge[h,l-1]))
                vcLower[h-1,l/2]=0.5*(YC[h,l-1]-YC[h-1,l-1])*((Center[h-1,l]-Center[h-1,l-1])+(Center[h,l]-Center[h,l-1]))
        vdLower=voLower-viLower
        vtLower=np.sum(vdLower)
        vctLower=np.sum(vcLower)
        
        viUpper=np.zeros([10,1])
        voUpper=np.zeros([10,1])
        vcUpper=np.zeros([10,1])
        for l in np.arange(1,1+0.0001,2):
            for h in np.arange(91,101,1):
                viUpper[h-91,l/2]=0.5*(YInner[h,l-1]-YInner[h-1,l-1])*((InnerEdge[h-1,l]-InnerEdge[h-1,l-1])+(InnerEdge[h,l]-InnerEdge[h,l-1]))
                voUpper[h-91,l/2]=0.5*(YOuter[h,l-1]-YOuter[h-1,l-1])*((OuterEdge[h-1,l]-OuterEdge[h-1,l-1])+(OuterEdge[h,l]-OuterEdge[h,l-1]))
                vcUpper[h-91,l/2]=0.5*(YC[h,l-1]-YC[h-1,l-1])*((Center[h-1,l]-Center[h-1,l-1])+(Center[h,l]-Center[h,l-1]))
        vdUpper=voUpper-viUpper
        vtUpper=np.sum(vdUpper)
        vctUpper=np.sum(vcUpper)
        RA=vt/vct
        RAL=vtLower/vctLower
        RAU=vtUpper/vctUpper
        return(vt,vct,RA,vtLower,vctLower,RAL,vtUpper,vctUpper,RAU,WidthAvg,HeightAvg,WidthStd,HeightStd,LinePlot,InnerPlot,OuterPlot)
        
BR=np.zeros([len(Param[:,0]),17])

for SampleNumber in range(len(Param[:,0])):   #len(Param[:,0])
    #QxQz map definition
    Angle = 121
    Pitch = Param[SampleNumber,2]
    SelectAngle1=Param[SampleNumber,3]
    SelectAngle2=Param[SampleNumber,4]
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
        
    Qxt = np.zeros([Angle,len(UQz)])
    Qzt = np.zeros([Angle,len(UQz)])

    for i in range(len(UQz)):
        Qxt[:,i]=QxI[i]
        Qzt[:,i]=np.arange(-1*UQz[i,0],UQz[i,0]+0.0001,2*UQz[i,0]/(Angle-1))
    
    Qx = np.zeros([5,len(UQz)])
    Qz = np.zeros([5,len(UQz)])
    
    Qx[0,:] = Qxt[60-SelectAngle1,:];
    Qz[0,:] = Qzt[60-SelectAngle1,:];
    Qx[1,:] = Qxt[60-SelectAngle2,:];
    Qz[1,:] = Qzt[60-SelectAngle2,:];
    Qx[2,:] = Qxt[60+SelectAngle1,:];
    Qz[2,:] = Qzt[60+SelectAngle1,:];
    Qx[3,:] = Qxt[60+SelectAngle2,:];
    Qz[3,:] = Qzt[60+SelectAngle2,:];
    Qx[4,:] = Qxt[60,:];
    Qz[4,:] = Qzt[60,:];
    #Sample Parameter Definition

    Trapnumber = 3
    # DW, I0 and Bk paramters are defined internally, structural pramters are defined externally
    DW = 1.3
    I0 = 0.01
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
    (Intensity,Amplitude)=SimInt_ID1(FITPAR)


    N=(1/(np.power(Intensity,0.5)))*Intensity # Generates  noise
    N=N*R
    Intensity2=Intensity+R; # Applies Noise

    C=CD.Misfit(Intensity,Intensity2)
    Chi2=np.sum((C))

    MCPAR=np.zeros([7])
    MCPAR[0] = 1 # Chainnumber
    MCPAR[1] = len(FITPAR)
    MCPAR[2] = 5000 #stepnumber
    MCPAR[3] = 0 #randomchains
    MCPAR[4] = 20 # Resampleinterval
    MCPAR[5] = 100 # stepbase
    MCPAR[6] = 100 # steplength 
  
    
    MCMCInitial=MCMCInit_ID1(FITPAR,FITPARLB,FITPARUB,MCPAR)
    Acceptprob=0;
    while Acceptprob < 0.3 or Acceptprob > 0.45:
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
        print(Acceptprob)
        if Acceptprob < 0.3:
            print(MCPAR[5])
            MCPAR[5]=MCPAR[5]+10
            MCPAR[6]=MCPAR[6]+10
        if Acceptprob > 0.45:
            print(MCPAR[5])
            MCPAR[5]=MCPAR[5]-10
            MCPAR[6]=MCPAR[6]-10
        
    start_time = time.perf_counter()
    MCPAR[0]=24
    MCPAR[2]=400000
    MCMCInitial=MCMCInit_ID1(FITPAR,FITPARLB,FITPARUB,MCPAR)
    MCMC_List=[0]*int(MCPAR[0])
    

    for i in range(int(MCPAR[0])):
        MCMC_List[i]=MCMCInitial[i,:]
    if __name__ =='__main__':  
        pool = Pool(processes=24)
              
        F=pool.map(MCMC_ID1,MCMC_List)
    
        Savename='P'+str(int(Pitch))+'_'+'W'+str(int(W))+'_'+'H'+str(int(H))+'_'+'A'+str(int(Angle))
        end_time=time.perf_counter()   
        print(end_time-start_time)    
        nonzerocount=0;      
        for i in range(int(MCPAR[0])):
            A=F[i]
            if A[0] != 0:
                nonzerocount=nonzerocount+1
        print(nonzerocount)
        UNCT=np.zeros([nonzerocount,13])
        nzc=-1        
        for i in range(int(MCPAR[0])):
            A=F[i]
            if A[0]!=0:
                nzc=nzc+1                
                UNCT[nzc,0:13]=A[0:13]
        UNCTAvg=np.average(UNCT,0)
        BR[SampleNumber,0]=W; BR[SampleNumber,1]=H; BR[SampleNumber,2]=Pitch; BR[SampleNumber,3]=Angle; BR[SampleNumber,4:17]=UNCTAvg;
        np.savetxt('3Angle_W12P32H25.csv',BR);
