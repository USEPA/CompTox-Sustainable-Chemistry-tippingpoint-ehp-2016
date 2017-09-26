import numpy as np
import math
import scipy.optimize as optz
import pandas as pd
import matplotlib.gridspec as gridspec
from util.misc import *
import scipy.optimize as optz

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.text as text
import matplotlib.font_manager as fm 
import pylab as pl
from matplotlib.patches import Ellipse, Circle, Rectangle
from matplotlib.patches import Ellipse, Circle, FancyArrowPatch

from collections import *
import matplotlib.cm as cm
import matplotlib as mpl
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import pdist

def calcCellRespSurface(ch_id,V,N=20,N2=200,use_t0=True,conc_deg=2,time_deg=2,
                        Feats=['p53','SK','OS','Mt','MM','MMP','MA','CCA','NS','CN','norm']):
    Vi = V.select(lambda x: x[0]==ch_id)
    AllZ = pd.DataFrame()
    
    for ft in Feats:
        Z1 = pd.DataFrame()
        for c in sorted(set([i[4] for i in Vi.index])):
            Vc = Vi[ft].select(lambda x:x[4]==c).reset_index()
            T0 = np.array(Vc.timeh,np.float)
            Y0 = np.array(Vc[ft],np.float)
            T1 = np.linspace(T0.min(),T0.max(),num=N)
            t_deg = time_deg
            if len(Vc.timeh.unique())==3:
                t_deg=2
            else:
                t_deg=1
            FY_splrp= interp.splrep(T0,Y0,s=1,k=t_deg)
            FY  = lambda x: interp.splev(x,FY_splrp,der=0)
            Y1  = FY(T1)
            Z1  = Z1.append(zip(T1,np.ones(N)*c,Y1,[ft]*N),ignore_index=True)

        Z1.columns=['timeh','lc','resp','ft']
        Z1=pd.DataFrame(Z1) 
        Z1['chem_name']=Vc.chem_name[0]
        Z1['ch_id']=Vc.ch_id[0]
        AllZ = AllZ.append(Z1)    

    FV= dict()
    DFV_DT=dict()
    DFV_DC=dict()
    D2FV_DCDT=dict()
    SPLRP=dict()

    for ft in Feats:
        Z2 = AllZ[(AllZ['ft']==ft)]
        #print# ft,Z2.shape
        SPLRP[ft]= interp.bisplrep(Z2.lc,Z2.timeh,Z2.resp,
                                 #w=np.array(Z1.w,dtype=np.float),
                                 kx=conc_deg,ky=time_deg,
                                 s=3)
                                 #quiet=1)
        FV[ft]=lambda c,t: interp.bisplev(c,t,SPLRP[ft])
        DFV_DT[ft] = lambda c,t: interp.bisplev(c,t,SPLRP[ft],dx=0,dy=1)
        DFV_DC[ft] = lambda c,t: interp.bisplev(c,t,SPLRP[ft],dx=1,dy=0)
        D2FV_DCDT[ft] = lambda c,t: interp.bisplev(c,t,SPLRP[ft],dx=1,dy=1)
        
    return dict(fx=FV,dx_dt=DFV_DT,dx_dc=DFV_DC,d2x_dcdt=D2FV_DCDT,splrp=SPLRP,allz=AllZ,chem=Vc.chem_name[0])


 
def calcTrajectoryTimeGrad(ch_id,V,use_times=[24,72],labs=['V_by_t']):
    Vc = V.select(lambda x: x[0]==ch_id)
    Z0 = list()
    for c in sorted(set([i[4] for i in Vc.index])):
        V_t = Vc.select(lambda x:x[4]==c and x[3] in use_times)
        D_t = [i[3] for i in V_t.index]
        
        if 0 in use_times:
            X1 = [0.0] + D_t
            Y1 = [0.0] + list(V_t.norm)
            
        else:
            X1 = D_t
            Y1 = list(V_t.norm)
        
        M1 = [c] + list(60*np.diff(Y1)/np.diff(X1))
        #M1 = 60*np.diff(Y1)/np.diff(X1)
        Z0.append(M1)
        #ax.plot(M1)
    Conc= set([i[4] for i in Vc.index])        
    M2 =pd.DataFrame(Z0,columns=['conc']+labs)
    M2['chem_name']=Vc.index[0][1]
    M2['ch_id']=ch_id
    # Add the 72 h norm and state
    Vct = Vc.select(lambda x: x[3]==use_times[-1])
    Vct = Vct.reset_index()
    M2['V_raw']=Vct['norm']
    #M2['S']=Vct['State']
    return M2

import random

def calcTrajectoryTimeGradSample(ch_id,V,use_times=[24,72],use_conc=9,labs=['V_by_t']):
    Vc = V.select(lambda x: x[0]==ch_id)
    Z0 = list()
    
    LC = list(set([i[4] for i in Vc.index]))
    LCi = sorted([LC[i] for i in random.sample(range(10),use_conc)])
    
    for c in LCi:
        V_t = Vc.select(lambda x:x[4]==c and x[3] in use_times)
        D_t = [i[3] for i in V_t.index]
        
        if 0 in use_times:
            X1 = [0.0] + D_t
            Y1 = [0.0] + list(V_t.norm)
            
        else:
            X1 = D_t
            Y1 = list(V_t.norm)
        
        M1 = [c] + list(60*np.diff(Y1)/np.diff(X1))
        #M1 = 60*np.diff(Y1)/np.diff(X1)
        Z0.append(M1)
        #ax.plot(M1)
    #return Z0
    M2 =pd.DataFrame(Z0,columns=['conc']+labs)
    M2['chem_name']=Vc.index[0][1]
    M2['ch_id']=ch_id
    # Add the 72 h norm and state
    Vct = Vc.select(lambda x: x[3]==use_times[-1] and x[4] in LCi)
    Vct = Vct.reset_index()
    M2['V_raw']=Vct['norm']
    #M2['S']=Vct['State']
    return M2

import numpy.linalg as LA
def residuals(X,Y,k=2,s=1):
    FY_splrp= interp.splrep(X,Y,s=s,k=k)
    FY  = lambda x: interp.splev(x,FY_splrp,der=0)
    return LA.norm(FY(X)-Y)
    
def calcTrajectoryProps(V_grad2):
    # Find the critical concentrations at which dV/dt ==0
    CC = optz.fsolve(DFM,cc0)

def calcTrajectoryDerivs(ch_id,V,k=3,cc0=-6,npnts=500,**kwargs):
    # This computes the time derivative of V 
    V_grad = calcTrajectoryTimeGrad(ch_id,V,**kwargs)
    # M contains V, V_by_t and S(tate)
    # Smooth gradient 
    FM_splrp= interp.splrep(V_grad.conc.tolist(),V_grad.V_by_t.tolist(),k=k,s=1)#,w=[3,3,3,2],s=10,k=2)
    FM  = lambda x: interp.splev(x,FM_splrp,der=0)
    # Smooth V
    FV_splrp= interp.splrep(V_grad.conc.tolist(),V_grad.V_raw.tolist(),k=k,s=1)#,w=[3,3,3,2],s=10,k=2)
    FV  = lambda x: interp.splev(x,FV_splrp,der=0)
    # Derivative of V with respect to concentration -- not really needed 
    DFV = lambda x: interp.splev(x,FV_splrp,der=1)
    # Derivatives of dV/dt with respect to concentration  
    DFM = lambda x: interp.splev(x,FM_splrp,der=1)
    D2FM = lambda x: interp.splev(x,FM_splrp,der=2)
    if k>2:
        D3FM = lambda x: interp.splev(x,FM_splrp,der=3)
    else:
        D3FM = None
        

    LC = np.linspace(V_grad.conc.min(),V_grad.conc.max(),npnts)
    lc1 = V_grad.conc.tolist()

    V_grad['V'] =FV(lc1)
    V_grad['dV_dc']=DFV(lc1)
    V_grad['dV_dt']=FM(lc1)
    V_grad['d2V_dtdc']=DFM(lc1)
    V_grad['d3V_dtdc2']=D2FM(lc1)
        
    #V_grad2 = pd.DataFrame(zip(LC,FV(LC),FM(LC),DFM(LC),D2FM(LC),D3FM(LC)),
    #                       columns=['conc','V','dV_dt','d2V_dtdc','d3V_dtdc2','d4V_dtdc3'])
    #Derivs = dict(dV_dc=DFV,dV_dt=FM,d2V_dtdc=DFM,d3V_dtdc2=D2FM,d4V_dtdc3=D3FM)

    V_grad2 = pd.DataFrame(zip(LC,FV(LC),DFV(LC),FM(LC),DFM(LC),D2FM(LC)),
                           columns=['conc','V','dV_dc','dV_dt','d2V_dtdc','d3V_dtdc2'])
    Funcs = dict(V=FV,dV_dc=DFV,dV_dt=FM,d2V_dtdc=DFM,d3V_dtdc2=D2FM)
    if k>2:
        V_grad['d4V_dtdc3']=D3FM(lc1)
        V_grad2['d4V_dtdc3']=D3FM(LC)
        Funcs['d4V_dtdc3']=D3FM
        
    return V_grad,V_grad2,Funcs

import scipy.interpolate as interp

def calcTrajectoryDerivsSample(ch_id,V,k=3,use_conc=9,cc0=-6,npnts=500,**kwargs):
    # This computes the time derivative of V 
    V_grad = calcTrajectoryTimeGradSample(ch_id,V,use_conc=use_conc,**kwargs)
    # M contains V, V_by_t and S(tate)
    # Smooth gradient 
    FM_splrp= interp.splrep(V_grad.conc.tolist(),V_grad.V_by_t.tolist(),k=k,s=1)#,w=[3,3,3,2],s=10,k=2)
    FM  = lambda x: interp.splev(x,FM_splrp,der=0)
    # Smooth V
    FV_splrp= interp.splrep(V_grad.conc.tolist(),V_grad.V_raw.tolist(),k=k,s=1)#,w=[3,3,3,2],s=10,k=2)
    FV  = lambda x: interp.splev(x,FV_splrp,der=0)
    # Derivative of V with respect to concentration -- not really needed 
    DFV = lambda x: interp.splev(x,FV_splrp,der=1)
    # Derivatives of dV/dt with respect to concentration  
    DFM = lambda x: interp.splev(x,FM_splrp,der=1)
    D2FM = lambda x: interp.splev(x,FM_splrp,der=2)
    if k>2:
        D3FM = lambda x: interp.splev(x,FM_splrp,der=3)
    else:
        D3FM = None
        
    LC = np.linspace(V_grad.conc.min(),V_grad.conc.max(),npnts)
    lc1 = V_grad.conc.tolist()

    V_grad['V'] =FV(lc1)
    V_grad['dV_dc']=DFV(lc1)
    V_grad['dV_dt']=FM(lc1)
    V_grad['d2V_dtdc']=DFM(lc1)
    V_grad['d3V_dtdc2']=D2FM(lc1)
        
    #V_grad2 = pd.DataFrame(zip(LC,FV(LC),FM(LC),DFM(LC),D2FM(LC),D3FM(LC)),
    #                       columns=['conc','V','dV_dt','d2V_dtdc','d3V_dtdc2','d4V_dtdc3'])
    #Derivs = dict(dV_dc=DFV,dV_dt=FM,d2V_dtdc=DFM,d3V_dtdc2=D2FM,d4V_dtdc3=D3FM)

    V_grad2 = pd.DataFrame(zip(LC,FV(LC),DFV(LC),FM(LC),DFM(LC),D2FM(LC)),
                           columns=['conc','V','dV_dc','dV_dt','d2V_dtdc','d3V_dtdc2'])
    Funcs = dict(V=FV,dV_dc=DFV,dV_dt=FM,d2V_dtdc=DFM,d3V_dtdc2=D2FM)
    if k>2:
        V_grad['d4V_dtdc3']=D3FM(lc1)
        V_grad2['d4V_dtdc3']=D3FM(LC)
        Funcs['d4V_dtdc3']=D3FM
        
    return V_grad,V_grad2,Funcs
    
def calcTrajCritical(Func,k=2,Vm=2):
    LC0 = np.linspace(-7,-3,num=200)
    I = np.where(np.abs(Func['d2V_dtdc'](LC0))<0.1)
    guess = np.unique(np.round(LC0[I],decimals=1))

    if len(guess)==0: 
        guess = [-8,-1]
        
    status = 1
    Roots = None
    Roots=optz.fsolve(Func['d2V_dtdc'],guess)
    Res = None
    #10**(Roots.x+6)

    if len(Roots)<0: return 0,None
    
    Res= pd.DataFrame(zip(Roots,10**(Roots+6),Func['V'](Roots),Func['d2V_dtdc'](Roots),np.round(Func['d3V_dtdc2'](Roots),decimals=2)),
                          #np.round(Func['d4V_dtdc3'](Roots),decimals=2)),
                      columns=['c_crit','c_uM','V','d2V_dtdc','d3V_dtdc2'])
    Res = np.round(Res,decimals=2).drop_duplicates()
    # conc where V=Vm
    
    FVm = lambda x:Func['V'](x)-Vm
    Cm=optz.fsolve(FVm,[-3])
    if abs(FVm(Cm[0]))<0.1:
        Res['c_max']=Cm[0]
        Res['c_max_uM']=10**(Cm[0]+6)
    else:
        Res['c_max']=None
        Res['c_max_uM']=None
        
    
    # Find the minimum concentration at which d3V_dtdc2 is positive 
    #X = Res[(Res['d3V_dtdc2']>0)]
    
    return 1,Res


def calcCritConc(CC_i):
    LC = sorted(CC_i.keys())
    c_max=LC[-1]
    d3V  =CC_i[c_max]
    Cls=dict(recov=0,norecov=0)
    Crit=dict(recov=[],norecov=[])
    
    if d3V>0:
       Cls['norecov']+=1.0
       Crit['norecov'].append(c_max)
    else:
       Cls['recov']+=1.0
       Crit['recov'].append(c_max)
    Cls['nroots']=len(CC_i)
    
    for k,v in Crit.iteritems(): 
        v = 10**(6+np.array(v))
        Cls[k+'_mn']=np.mean(v)
        Cls[k+'_sd']=np.std(v)
    return Cls

def calcCritConcSummary(CC):
    Cls=dict(rs_recov=0,rs_norecov=0)
    Crit=dict(rs_recov=[],rs_norecov=[])
    Roots=[]
    for CC_i in CC:
        LC = sorted(CC_i.keys())
        c_max=LC[-1]
        d3V  =CC_i[c_max]
        Roots.append(len(CC_i))
        if d3V>0:
           Cls['rs_norecov']+=1.0
           Crit['rs_norecov'].append(c_max)
        else:
           Cls['rs_recov']+=1.0
           Crit['rs_recov'].append(c_max)
    Cls['rs_recov'] /= len(CC)
    Cls['rs_norecov'] /= len(CC)
    Cls['rs_nroots']=int(np.median(Roots))
    for k,v in Crit.iteritems(): 
        v = 10**(6+np.array(v))
        Cls[k+'_mn']=np.mean(v)
        Cls[k+'_sd']=np.std(v)
    return Cls

FT=['p53Act','StressKinase','OxidativeStress','MicrotubuleCSK','MitoMass','MitoMembPot',
    'MitoticArrest','CellCycleArrest','NuclearSize','CellLoss']
LB=['p53 Act','Stress Kinase','Oxidative Stress','Microtubule','Mito Mass','Mito Memb Pot',
    'Mitotic Arrest','Cell Cycle Arrest','Nuclear Size','Cell Number']
FTLB = dict(zip(FT,LB))
