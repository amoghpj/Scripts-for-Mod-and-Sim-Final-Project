###############################################################################################
####################### data used for fiting in this script are provided in separate .txt files
from PyDSTool import *
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
plt.rcParams['figure.figsize']=12,8

# Two cell model
T0=7e9
rt=0.9
MiceTwoTargetCells=args(name='Tropism_model',
                    varspecs={'Td':'-beta26*Td*V',
                              'Ts':'-beta23*Ts*V',
                              'Ed':'beta26*Td*V-k*Ed',
                              'Es':'beta23*Ts*V-k*Es',
                              'Id':'k*Ed-delta*Id',
                              'Is':'k*Es-delta*Is',
                              'V':'p26*Id+p23*Is-c*V'
                             },
                    pars={'k':0.23809523809523808,
                          'beta26':1e-5/24.0,
                          'beta23':1e-5/24.0,
                          'delta':0.3448275862068966,
                          'c':0.3448275862068966,
                          'p26':0.00833333333,
                          'p23':0.00833333333
                         },
                    ics={'Td':(1-rt)*T0,
                         'Ts':rt*T0,
                         'Id':0.0,
                         'Is':0.0,
                         'Ed':0.0,
                         'Es':0.0,
                         'V':7.6e-2
                        },
                        tdata=[0,24.0*10]
                   )
MiceTwoTargetCellsDS=Generator.Vode_ODEsystem(MiceTwoTargetCells)


# One cell model
T0=7.0e9
MiceOneTargetCell=args(name='Tropism_model',
                    varspecs={'Td':'-beta*Td*V',
                              'Ed':'beta*Td*V-k*Ed',
                              'Id':'k*Ed-delta*Id',
                              'V':'p*Id-c*V'
                             },
                    pars={'k':0.23809523809523808,
                          'beta':1e-5/24.0,
                          'delta':0.3448275862068966,
                          'c':0.3448275862068966,
                          'rb':0,
                          'rp':0,
                          'p':0.00833333333
                         },
                    ics={'Td':T0,                         
                         'Id':0.0,                        
                         'Ed':0.0,                         
                         'V':7.6e-2},
                       tdata=[0,24*10]
                       )
                       
                   
MiceOneTargetCellDS=Generator.Vode_ODEsystem(MiceOneTargetCell)

#k,delta,c,beta,p
SingleCellParameters=[[1.0/(4.4*1e5), 1.0/7.3, 1.0/1.2, (1.1*1e-3)/24.0, (7.7)/24.0],
                      [1.0/240.0, 1.0/130.0, 1.0/0.41, (6.9*1e-5)/24.0, (0.019)/24.0],
                      [1.0/33.0, 1.0/32.0, 1.0/34.0, (9.8*1e-5)/24.0, (1.2*1e-3)/24.0],
                      [1.0/16.0, 1e-15, 1.0/7.1, (5.9*1e-7)/24.0, (2.7*1e-3)/24.0]
                     ]
SingleCellIcs=[5.9*1e-3, 8.2, 65.0, 1.5*1e5]

# k,delta, c,beta23,beta26,p23,p26
TwoCellParameters=[[1.0/24.0, 1/0.74, 1.0/0.8, 2.4*1e-8/24.0, 8.2*1e-4/24.0, 77.0/24.0, 0.3*1e-2/24.0],
                   [1.0/9.2, 1.0/1.1, 1.0/1.3, 2.7*1e-7/24.0, 3.3*1e-5/24.0, 2.9/24.0, 3.9*1e-3/24.0],
                   [1.0/4.9,1.0/5.0, 1.0/5.0, 5.2*1e-8/24.0, 2.4*1e-5/24.0, 0.9/24.0, 3.3*1e-3/24.0],
                   [1.0/2.9, 1.0/3.5, 1.0/3.8, 4.1*1e-8/24.0, 9.4*1e-5/24.0, 2.4/24.0,2.2*1e-4/24.0]
                  ]
TwoCellIcs=[4.1*1e-3,1.4,19.0,4.8*1e3]


pointsOneCell=[]
pointsTwoCell=[]
for i in range(0,len(SingleCellParameters)):
    MiceOneTargetCellDS.set(ics={'V':SingleCellIcs[i]},
                            pars={'k':SingleCellParameters[i][0],
                                  'delta':SingleCellParameters[i][1],
                                  'c':SingleCellParameters[i][2],
                                  'beta':SingleCellParameters[i][3],
                                  'p':SingleCellParameters[i][4],
                                 }
                           )
    pts=MiceOneTargetCellDS.compute('params'+str(i)).sample(dt=0.1)
    times=pts['t']/24.0
    pointsOneCell.append(pts['V'])
    
for i in range(0,len(SingleCellParameters)):
    MiceTwoTargetCellsDS.set(ics={'V':TwoCellIcs[i]},
                            pars={'k':TwoCellParameters[i][0],
                                  'delta':TwoCellParameters[i][1],
                                  'c':TwoCellParameters[i][2],
                                  'beta26':TwoCellParameters[i][3],
                                  'beta23':TwoCellParameters[i][4],
                                  'p26':TwoCellParameters[i][5],
                                  'p23':TwoCellParameters[i][6]
                                 }
                           )
    pts=MiceTwoTargetCellsDS.compute('params'+str(i)).sample(dt=0.1)
    #times=pts['t']/24.0
    pointsTwoCell.append(pts['V'])


# This is data extracted from the plots directly. The processed values are stored in seprate .txt files
data=[[884.0,450.0,320.0,282.0,251.0],[896.0,395.0,467.0,369.0,346.0],[367.0,218.0,164.0,117.0,117.0],[378.0,148.0,155.0,142.0,192.0]]
zero_def=1018.0
max_def=80.0
abs_scale=[]
rowise=[]
for row in data:
    rowise=[]
    for item in row:
        rowise.append((zero_def-item)/133.0)#-max_def
    abs_scale.append(rowise)
abs_scale
scaled_data=[]
nrow=[]
for row in abs_scale:
    nrow=[]
    for item in row:
        nrow.append(10**item)
    scaled_data.append(nrow)
f,((ax1,ax2),(ax3,ax4))=plt.subplots(2,2)

ax1.plot(times,pointsOneCell[0],'r')
ax1.plot(times,pointsTwoCell[0],'b')
ax1.set_yscale('log')
ax1.plot(data_times,scaled_data[1],'ks')
ax1.set_ylim([1e0,1e8])
ax1.set_title('1991')
ax1.set_xlabel('Time (dpi)')
ax1.set_ylabel('Viral Titer')

ax2.plot(times,pointsOneCell[1],'r')
ax2.plot(times,pointsTwoCell[1],'b')
ax2.set_yscale('log')
ax2.plot(data_times,scaled_data[0],'ks')
ax2.set_ylim([1e0,1e8])
ax2.set_title('Thailand/83')
ax2.set_xlabel('Time (dpi)')
ax2.set_ylabel('Viral Titer')

ax3.plot(times,pointsOneCell[2],'r')
ax3.plot(times,pointsTwoCell[2],'b')
ax3.set_yscale('log')
ax3.plot(data_times,scaled_data[3],'ks')
ax3.set_ylim([1e0,1e8])
ax3.set_title('1918')
ax3.set_xlabel('Time (dpi)')
ax3.set_ylabel('Viral Titer')

ax4.plot(times,pointsOneCell[3],'r')
ax4.plot(times,pointsTwoCell[3],'b')
ax4.set_yscale('log')
ax4.plot(data_times,scaled_data[2],'ks')
ax4.set_ylim([1e0,1e8])
ax4.set_title('Thailand/16')
ax4.set_xlabel('Time (dpi)')
ax4.set_ylabel('Viral Titer')
    #ax[i].plot(data_fit[i])
plt.savefig('Dobrovolny_Fig4.pdf')
