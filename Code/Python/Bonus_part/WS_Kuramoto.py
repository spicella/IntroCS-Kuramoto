#!/usr/bin/env python
# coding: utf-8

# In[9]:


import pandas as pd
import numpy as np
import math
import random 
import matplotlib.pyplot as plt
plt.style.use("seaborn-darkgrid")
import os


# In[16]:


def coordinates_on_circle(n):
    """Returns x,y coordinates of points on an unit circle with spacing 2π/n"""
    if type(n)!=int:
        raise Exception("Wrong input: \n the argument must be an integer number of points on the unit circle!")
    x,y = [],[]
    step_angle = 2*math.pi/n
    for i in range(0,n):
        x.insert(len(x),math.cos(i*step_angle))
        y.insert(len(y),math.sin(i*step_angle))
    return x,y


def create_starting_graph(n,r):
    if type(r)!=int:
        raise Exception("Wrong input: \n r must be an integer number of edges between vertices")
    if r>n-1:
        raise Exception("Wrong input: \n r must not exceed n-1!")
    coords = coordinates_on_circle(n)
    
    #create adjacency_matrix as pandas df
    #Initialize adjacency matrix
    adj_mat_df = pd.DataFrame(np.zeros([n,n]),dtype='bool')
    #Make starting connections with pbc
    for i in range(0,n):
        #left
        if(i-r>=0):
            adj_mat_df.iloc[i][i-r:i] = True
        else:
            diff = r-i
            adj_mat_df.iloc[i][0:i] = True
            adj_mat_df.iloc[i][n-diff:n+1] = True
        #right
        if(i+r<n):
            adj_mat_df.iloc[i][i+1:i+r+1] = True #+1 to avoid self loop and up to sym value
        else:
            diff = i+r-n
            adj_mat_df.iloc[i][i+1:n+1] = True
            adj_mat_df.iloc[i][0:diff+1] = True

    return adj_mat_df


def create_strogatz(n, r ,p, place_labels=False ):
    """Plots the graph of the Strogatz model on a unit circle."""
    print("Creating WS graph for p=%.2f"%(p))
    #Procedure to create results folder automatically
    path = os.getcwd()
    results_dir = "/results_WS" 
    try:
        os.mkdir(path+results_dir)
    except OSError:
        print ("Creation of the directory %s failed" % results_dir)
    else:
        print ("Successfully created the directory %s " % results_dir)
    #names for file paths
    name_plot = "/plot_n%d_r%d_p%.3f.png" %(n,r,p)
    name_csv = "/data_n%d_r%d_p%.3f.csv" %(n,r,p)
    name_plot_rewired = "/plot_rewired_n%d_r%d_p%.3f.png" %(n,r,p)
    name_csv_rewired = "/data_rewired_n%d_r%d_p%.3f.csv" %(n,r,p)
    
    
    #check for errors
    if p>1 or p<0:
        raise Exception("Wrong input: \n p must be in [0,1]")
    coords = coordinates_on_circle(n)
    adj_mat = create_starting_graph(n,r)
 
    
    for i in range(0,n):
        for j in range(0,r):  #for each link to vertex i

            if (random.random()<p): #attempt a rewire
                #performing the rewire
                #    - Choose which of the connected edge to rewire->deleated_edge
                #    - Choose were to rewire it among the available positions->candidates
                #    - Perform the connection/deleate old connection/update mirror adjmat

                #choose which edge to remove: [+periodic boundary conditions]
                deleated_edge = i+1+j
                if deleated_edge>n-1:
                    deleated_edge = deleated_edge-n
                
                #chose available position:
                candidates = list(adj_mat[adj_mat.iloc[i] == False].index.tolist())
                candidates.remove(i) #take out self loop
                new_edge = random.choice(candidates)

                #print("candidates list = ",candidates)
                #print("new edge chosen = ",new_edge)
                
                #create new wire
                adj_mat.iloc[i][new_edge]=True
                adj_mat.iloc[new_edge][i]=True

                #deleate old wire
                adj_mat.iloc[i][deleated_edge]=False
                adj_mat.iloc[deleated_edge][i]=False
                


    return adj_mat

def write_extended(adj,p_ext):
    """Write the edge list and the associated degree for each vertex
    Starting from an adjacency matrix and the p value used (just for outputting the right name)
    col[0]-->from, col[1]-->to, col[2]-->freq"""
    result = []
    for i in range(0,n): #rows
        for j in range(0,n): #cols
            if(adj[i][j] == True):
                result.append([i,j,0])#third column is the degree of distribution for each vertex, useful for c computation  
    result = pd.DataFrame(result,index=None,columns=['from','to','freq']) #columns: [from,to]
    a = result.iloc[:,0:1]
    b = a.apply(pd.value_counts, sort=False)
    b['index'] = b.index
    for i in range(0,len(a)):
        result['freq'][i]=(b.iloc[a.iloc[i,0]][0])
    result.to_csv(path+results_dir+"/extended_adj_mat_p=%.2f.csv"%(p_ext),header=False, index=False)
    return result


path = os.getcwd()
results_dir = "/results_WS" 
n,r=2000,2 #insert same values as in C code!
p_list = [0,0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, .75, 0.8, .85, 0.9,.95,1]

for i in p_list:
    adj_mat = create_strogatz(n,r,i)
    write_extended(adj_mat,i)
