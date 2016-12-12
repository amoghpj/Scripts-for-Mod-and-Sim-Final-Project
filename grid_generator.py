Gsize=80
Cluster_coordinates=[]
num_clusters=1
abundance=0.5
num_obj=abundance*Gsize**2
neighborhood=0
while neighborhood**2<int(num_obj/num_clusters):
    if neighborhood**2<int(num_obj/num_clusters):
        neighborhood=neighborhood+1
radius=int(neighborhood*0.3)
# Randomly section the 'ordered' objects
# Partitions=random.randint(0,num_obj,num_clusters)
# sorted_parts=sort(Partitions)
objs_in_clusters=[]
# for i in range(1,len(sorted_parts)):
#     objs_in_clusters.append(sorted_parts[i]-sorted_parts[i-1])
for i in range(0,num_clusters):
    objs_in_clusters.append(int(num_obj/num_clusters))

# Assign random positions to clusters on the grid
# Cluster_positions_values=sort(random.choice(range(0,Gsize**2),num_clusters+1))
x=0
y=0
Cluster_rand_c=set((x,y) for x in range(0,Gsize) for y in range(0,Gsize))
Choose_coords=random.choice(range(0,len(Cluster_rand_c)),num_clusters,replace=False)
for i in range(0,num_clusters):
    Cluster_coordinates.append(list(Cluster_rand_c)[Choose_coords[i]])
i=0
# x_positions=[]
# y_positions=[]
# for i in range(0,num_clusters):
    
#     random.shuffle(Cluster_coordinates)
#     x_positions.append(Cluster_coordinates[0])
#     random.shuffle(Cluster_coordinates)
#     y_positions.append(Cluster_coordinates[0])

# x_positions=random.choice(,num_clusters,replace=False)
# y_positions=random.choice(range(0,Gsize),num_clusters,replace=False)
# Convert Cluster values to positions
counter=0
cum_counter=0
Cluster_positions_values=[]
Cluster_positions_values=[[40,40]]#Cluster_coordinates
# for i in range(0,len(x_positions)):
#     Cluster_positions_values.append([x_positions[i],y_positions[i]])
# for i in range(0,Gsize):
#     for j in range(0,Gsize):
#         cum_counter=cum_counter+1
#         if cum_counter==Cluster_positions_values[counter] and counter<len(Cluster_positions_values)-1:
#             Cluster_positions.append([i,j])
#             counter=counter+1

# Generating the grid
GRID=[]
row=[]
for i in range(0,Gsize):
    row=[]
    for j in range(0,Gsize):
        row.append('Td')
    GRID.append(row)

def place_clusters(pos,num):
    cluster_x=pos[0]
    cluster_y=pos[1]
#     print(type(cluster_x))
#     print(cluster_x,cluster_y)
    neighborhood=0.0
    while neighborhood**2<num:
        if neighborhood**2<num:
            neighborhood=neighborhood+1
    print('Neighborhood=%i'%neighborhood)
    radius=int(neighborhood*0.7)
#     print(radius)
    object_positions=[]
#     xrand=[]
#     yrand=[]
#     X_RANGE=[]
#     Y_RANGE=[]%%!
    x=0
    y=0
    Superset_coordinates=set((x,y) for x in range(cluster_x-radius,cluster_x+radius-1) for y in range(cluster_y-radius,cluster_y+radius-1))
#     X_RANGE=list(range(x-radius,x+radius-1))
#     Y_RANGE=list(range(y-radius,y+radius-1))
    print(num)
    indices=random.choice(range(0,len(Superset_coordinates)),num,replace=False)
    for i in range(0,len(indices)):
        object_positions.append(list(Superset_coordinates)[indices[i]])
#     for i in range(0,num):
#         random.shuffle(X_RANGE)
#         xrand.append(X_RANGE[0])
#         random.shuffle(Y_RANGE)
#         yrand.append(Y_RANGE[0])
#     xrand=random.choice(range(x-radius,x+radius),num,replace=False)
#     yrand=random.choice(range(y-radius,y+radius),num,replace=False)
#     for i in range(0,num):
#         object_positions.append([xrand[i],yrand[i]])
    return object_positions
coordinates=[]
count=0


for i in range(0,num_clusters):
    coordinates.append(place_clusters(Cluster_positions_values[i],objs_in_clusters[i]))

for i in range(0,len(coordinates)):
    for j in range(0,len(coordinates[i])):
        GRID[coordinates[i][j][0]][coordinates[i][j][1]]='Ts'
f,ax=plt.subplots()
ax.pcolor(heatmap_generator(GRID,Gsize),vmin=0,vmax=6)
