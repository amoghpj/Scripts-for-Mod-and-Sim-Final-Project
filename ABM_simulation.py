import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
############### LOAD GRID
####### Add string of grid file name here
Grid_file_name='GRID_1clusters_0.5-2.txt'
Gridsize=80 ## Modify grid size according to you grid input

grid_vals=pd.read_csv(Grid_dimension,sep='\t',header=None)
del grid_vals[Gridsize] ####Modify this according to your grid size 
GRID=grid_vals.values.tolist()
Viral_Abstract=[]
Virus_tracker=[]
Td_counter=[]
Ed_counter=[]
Id_counter=[]    
Ts_counter=[]
Es_counter=[]
Is_counter=[]

#Default of 10 iterations

for iteration in tqdm(range(0,10)):
    Grid_dimension=Gsize
    # Counters
    
    V_levels=[]
    
    Td=[]
    Id=[]
    Ed=[]
    Ts=[]
    Is=[]
    Es=[]
    # Time 
    Total_time=300
    current_time=0

    # State variables
    Cell_grid=[]
    Virus_grid=[]
    Virus_grid_update=[]
    Time_series_grid=[]
    row_states=[]
    viral_states=[]
    prev_state=[]
    Initial_Cell_states=GRID
 
    Time_series_grid.append(Initial_Cell_states)
    
    
    # Transition rules:
    beta=1.0/24.0
    k=1.0/4.2
    delta=1.0/2.9
    c=1.0e-1/2.9
    p=7.6e-2
    r_beta=1e-3
    r_p=1e+9
    Initial_viral_concentration=p/2#p/2.0
    
    # Virus grid: Initialization
    for i in range(0,Grid_dimension):
#         row_states=[]
        viral_states=[]
        for j in range(0,Grid_dimension):
            #row_states.append('T')
            viral_states.append(Initial_viral_concentration)
        #Initial_Cell_states.append(row_states)
        Virus_grid.append(viral_states)
    
    V_levels.append(sum(sum(np.array(Virus_grid))))
    V_levels.append(sum(sum(np.array(Virus_grid))))

    for t in range(1,Total_time):
        Cell_grid=[]
        prev_state=[]
        prev_state=Time_series_grid[t-1]
        for i in range(0,Grid_dimension):
            row_states=[]
            for j in range(0,Grid_dimension):
                ############ Flip a coin #########
                rand=random.random()
                ##################################
                # Td-> Ed
                if prev_state[i][j]=='Td':
                    if rand>(1-beta*viral_influence(i,j,Virus_grid,Grid_dimension)):
                        row_states.append('Ed')
                    else:
                        row_states.append(prev_state[i][j])
                
                # Ed-> Id
                elif prev_state[i][j]=='Ed':
                    if rand>(1-k):
                        row_states.append('Id')
                    else:
                        row_states.append(prev_state[i][j])
                        
                # Id -> D
                elif prev_state[i][j]=='Id':
                    if rand>(1-delta):
                        row_states.append('D')
                    else:
                        row_states.append(prev_state[i][j])
                        
                # Ts -> Es
                elif prev_state[i][j]=='Ts':
                    if rand>(1-r_beta*beta*viral_influence(i,j,Virus_grid,Grid_dimension)):
                        row_states.append('Es')
                    else:
                        row_states.append(prev_state[i][j])
                        
                # Es -> Is
                elif prev_state[i][j]=='Es':
                    if rand>(1-k):
                        row_states.append('Is')
                    else:
                        row_states.append(prev_state[i][j])
                
                # Is -> D
                elif prev_state[i][j]=='Is':
                    if rand>(1-delta):
                        row_states.append('D')
                    else:
                        row_states.append(prev_state[i][j])
                
                # D -> D
                elif prev_state[i][j]=='D':
                    row_states.append('D')
                
            Cell_grid.append(row_states)
        Time_series_grid.append(Cell_grid)
        
        Virus_grid_update=[]
        
        for i in range(0,Grid_dimension):
            row_states=[]
            for j in range(0,Grid_dimension):
                if Cell_grid[i][j]=='Id':
                    row_states.append(Virus_grid[i][j]*(1-c)+10*p)
                elif Cell_grid[i][j]=='Is':
                    row_states.append(Virus_grid[i][j]*(1-c)+10*p*r_p)
                else:
                    row_states.append(Virus_grid[i][j]*(1-c))
            Virus_grid_update.append(row_states)
        Virus_grid=Virus_grid_update
        V_levels.append(sum(sum(np.array(Virus_grid))))
        Virus_tracker.append(Virus_grid)

    Viral_Abstract.append(V_levels)

    for i in range(0,len(Time_series_grid)):
        Td.append(celltype_counter('Td',Time_series_grid[i],Grid_dimension))
        Id.append(celltype_counter('Id',Time_series_grid[i],Grid_dimension))
        Ed.append(celltype_counter('Ed',Time_series_grid[i],Grid_dimension))
        Ts.append(celltype_counter('Ts',Time_series_grid[i],Grid_dimension))
        Is.append(celltype_counter('Is',Time_series_grid[i],Grid_dimension))
        Es.append(celltype_counter('Es',Time_series_grid[i],Grid_dimension))
    Td_counter.append(Td)
    Id_counter.append(Id)
    Is_counter.append(Is)
    Ts_counter.append(Ts)
    Ed_counter.append(Ed)
    Es_counter.append(Es)


f,((ax1,ax2),(ax3,ax4))=plt.subplots(2,2)
# for i in range(0,len(Viral_Abstract)):
def normalized(A,maxval):
    minA=min(A)
    norm=[]
    for i in range(0,len(A)):
        norm.append(float(A[i]-minA)/float(maxval-minA))
    return norm
for i in range(0,len(Viral_Abstract)):
    ax1.plot(Viral_Abstract[i])
#     ax2.plot(Td_counter[i])
#     ax2.plot(Ts_counter[i],'--')
    ax2.plot(normalized(Td_counter[i],3200))
    ax2.plot(normalized(Ts_counter[i],3200),'--',alpha=0.75)
    ax3.plot(Id_counter[i])
    ax3.plot(Is_counter[i],'--',alpha=0.75)
    ax4.plot(Ed_counter[i])
    ax4.plot(Es_counter[i],'--',alpha=0.75)

ax1.set_title('V titer')
# ax1.set_ylim([1e-1,1e3])
# ax1.set_yscale('log')
ax2.set_title('T')
# ax2.set_ylim([0,3500])
# ax2.set_yscale('log')
ax3.set_title('I')
ax3.set_yscale('log')
ax4.set_title('E')
# ax4.set_yscale('log')
plt.savefig('rbeta'+str(r_beta)+'_rp='+str(r_p)+'.pdf')
