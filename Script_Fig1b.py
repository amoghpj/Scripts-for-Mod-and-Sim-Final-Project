#####################################################
############# Defining the two cell model ###########
from PyDSTool import *
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

plt.rcParams['figure.figsize']=12,8
T0=4e8
rt=0.5
TwoTargetCells=args(name='Tropism_model',
                    varspecs={'Td':'-beta*Td*V',
                              'Ts':'-rb*beta*Ts*V',
                              'Ed':'beta*Td*V-k*Ed',
                              'Es':'rb*beta*Ts*V-k*Es',
                              'Id':'k*Ed-delta*Id',
                              'Is':'k*Es-delta*Is',
                              'V':'p*Id+rp*p*Is-c*V'
                             },
                    pars={'k':0.23809523809523808,
                          'beta':1e-5/24.0,
                          'delta':0.3448275862068966,
                          'c':0.3448275862068966,
                          'rb':0,
                          'rp':0,
                          'p':0.00833333333
                         },
                    ics={'Td':(1-rt)*T0,
                         'Ts':rt*T0,
                         'Id':0.0,
                         'Is':0.0,
                         'Ed':0.0,
                         'Es':0.0,
                         'V':7.6e-2
                        }
                   )
TwoTargetCellsDS=Generator.Vode_ODEsystem(TwoTargetCells)

##########################################################################
###################### Generating Fig 1b from the paper ##################
rp=[4.0,40.0,400.0]

#rp=[10.0,100.0,1000.0]
# rp=[40]
rb=5e-3
f,(ax1,ax2,ax3)=plt.subplots(1,3)
# f,ax2=plt.subplots(1,1)
def normalize(X):
#     minX=m
    maxX=max(X)
    Y=[]
    for i in range(0,len(X)):
        Y.append((X[i])/(maxX))
    return Y

Td=[]
Ts=[]
time=[]
V=[]
for i in range(0,len(rp)):
    TwoTargetCellsDS.set(pars={'rp':rp[i],'rb':rb},tdata=[0,24*14])
    points=TwoTargetCellsDS.compute('twotarget').sample(dt=0.1)
    Td.append(points['Td'])
    Ts.append(points['Ts'])
    time.append(points['t']/24)
    V.append(points['V'])

ax1.set_ylabel('Viral Titer')
ax1.plot(time[0],V[0],'k',lw=2)
ax11=ax1.twinx()
ax11.plot(time[0],normalize(Td[0]),'b--')
ax11.plot(time[0],normalize(Ts[0]),'b')
ax11.set_ylim([0,1])
ax1.set_yscale('log')
ax1.set_ylim([1e0,1e9])
ax1.set_xlabel('Time (dpi)')

ax22=ax2.twinx()
ax22.plot(time[1],normalize(Td[1]),'b--')
ax22.plot(time[1],normalize(Ts[1]),'b')
ax22.set_ylim([0,1])
ax2.plot(time[1],V[1],'k',lw=2)
ax2.set_yscale('log')
ax2.set_ylim([1e0,1e9])
ax2.set_xlabel('Time (dpi)')

ax33=ax3.twinx()
ax33.plot(time[2],normalize(Td[2]),'b--')
ax33.plot(time[2],normalize(Ts[2]),'b')
ax33.set_ylim([0,1])
ax3.plot(time[2],V[2],'k',lw=2)
ax3.set_yscale('log')
ax3.set_ylim([1e0,1e10])
ax3.set_xlabel('Time (dpi)')
ax33.set_ylabel('Target Cell Abundance')
