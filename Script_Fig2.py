############### Headers ########################
from PyDSTool import *
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
plt.rcParams['figure.figsize']=10,6
T0=4e8

#######################################################################
###### Vary r_T such that r_t=[0.91,0.92,0.93,0.94]
rt=0.91
#######################################################################
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

########## Resolution of heatmap ################
num_points=100
#################################################
rbrange=np.logspace(-3,3,num_points)
rprange=np.logspace(-3,3,num_points)
max_V=[]
times_const_RB_RP=[]
times_const_RP=[]
for RB in tqdm(rbrange):
    times_const_RP=[]
    for RP in rprange:
        TwoTargetCellsDS.set(pars={'rp':RP,'rb':RB},tdata=[0,24*15]) 
        points=TwoTargetCellsDS.compute('twotarget'+str(RB+RP)).sample(dt=0.1)
        times_const_RP.append(points['t'][points['V']==max(points['V'])]/24.0)
    times_const_RB_RP.append(times_const_RP)

f1=open('times_at_max_viral_titer_'+str(num_points)+'_rt='+str(rt)+'_.txt','w')
for row in times_const_RB_RP:
    for value in row:
        f1.write('%s\t'%str(float(value)))
f1.close()

fname='times_at_max_viral_titer_100_rt='+str(rt)+'_'
D1=pd.read_csv(fname+'.txt',sep='\t',header=None)
A=list(D1.stack())
B=reshape(A,(100,100))
f,ax=plt.subplots()
xlabels=['1e-3','1e-2','1e-1','0','1e1','1e2','1e3']
ax.set_xticks(np.linspace(0,100,7),minor=False)
ax.set_xticklabels(xlabels)
ax.set_yticks(np.linspace(0,100,7),minor=False)
ax.set_yticklabels(xlabels)
ax.set_xlabel('$r_p$')
ax.set_ylabel('$r_\\beta$')
p=plt.pcolor(B,vmin=1.0, vmax=B.max())
cbar=f.colorbar(p)
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
plt.savefig(fname+'.pdf')

