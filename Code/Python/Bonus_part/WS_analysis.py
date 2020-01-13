#!/usr/bin/env python
# coding: utf-8

# In[199]:


import os
from sympy import *
import pandas as pd
import numpy as np
import scipy.fftpack
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use("seaborn-paper")
from mpl_toolkits.axes_grid1 import make_axes_locatable



def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def find_nearest_value(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


# ## Canvas palette

# In[2]:


#Canvas for single plot
x = np.linspace(0,10,100)
y = np.sin(x)
plt.figure(figsize=[14,6])
plt.grid(True)
plt.title("Change-me!",fontsize=20)
plt.plot(x,y,label="testvalue")
plt.legend(fontsize=16)
plt.xlabel("XLABEL (unit)",fontsize=18)
plt.ylabel("YLABEL (unit)",fontsize=18)
plt.show()


# In[3]:


#Canvas for side by side
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14,6))
fig.suptitle("test",y=1.05,fontsize=20)

axes[0].grid(True)
axes[0].plot(x,y,label="testvalue")
axes[0].legend(fontsize=16)
axes[0].set_title("TESTTITLE",fontsize=18)
axes[0].set_xlabel("XLABEL (unit)",fontsize=18)
axes[0].set_ylabel("YLABEL (unit)",fontsize=18)
axes[0].legend(fontsize=16)
axes[0].tick_params(axis='both', which='major', labelsize=15)


axes[1].grid(True)
axes[1].plot(x,y,label="testvalue")
axes[1].legend(fontsize=16)
axes[1].set_title("TESTTITLE",fontsize=18)
axes[1].set_xlabel("XLABEL (unit)",fontsize=18)
axes[1].set_ylabel("YLABEL (unit)",fontsize=18)
axes[1].legend(fontsize=16)
axes[1].tick_params(axis='both', which='major', labelsize=15)

fig.tight_layout()
plt.show()


# In[4]:


#Canvas for side by side
fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(14,6))
fig.suptitle("test",y=1.05,fontsize=20)

axes[0,0].grid(True)
axes[0,0].plot(x,y,label="testvalue")
axes[0,0].legend(fontsize=16)
axes[0,0].set_title("TESTTITLE",fontsize=18)
axes[0,0].set_xlabel("XLABEL (unit)",fontsize=18)
axes[0,0].set_ylabel("YLABEL (unit)",fontsize=18)
axes[0,0].legend(fontsize=16)
axes[0,0].tick_params(axis='both', which='major', labelsize=15)


axes[0,1].grid(True)
axes[0,1].plot(x,y,label="testvalue")
axes[0,1].legend(fontsize=16)
axes[0,1].set_title("TESTTITLE",fontsize=18)
axes[0,1].set_xlabel("XLABEL (unit)",fontsize=18)
axes[0,1].set_ylabel("YLABEL (unit)",fontsize=18)
axes[0,1].legend(fontsize=16)
axes[0,1].tick_params(axis='both', which='major', labelsize=15)

fig.tight_layout()
plt.show()


# ## Read data

# In[8]:


#Folder and paths definitions
main_path  = os.getcwd()
datafolder_path = main_path+"/results"
results_dir = "/output_py" 
output_dir = main_path+results_dir
try:
    os.mkdir(output_dir)
except OSError:
    print ("Creation of the directory %s failed" % results_dir)
else:
    print ("Successfully created the directory %s " % results_dir)


# In[ ]:





# In[277]:


#Simulation parameters
N = 2000
T = 20000
n_runs = 20
dt = .01
freq = "gfreq"
gamma = 10
MF = "NOMF"

if(freq =="gfreq"):
    freq_plot="$\\vec{\\omega_{0}} = \\mathcal{N}(0,1)$"
else:
    freq_plot="$\\vec{\\omega_{0}} = U(-%.1f,%.1f)$"%(gamma,gamma)
if(MF =="MF"):
    MF_plot="MeanField"
else:
    MF_plot="non-meanField"


# In[278]:


#OutputFileNames
#S/N --> |r(t)|/sigma(r(t))
sn_name = "S_N"
#(Mod&Phase)(t)
modphase_name = "ModPhase_t"
#(Mod)(t)
mod_name = "Mod_t"
#Spectrum
spectrum_name = "Spectrum"
#r_inf
rinf_name = "r_inf"
#Configuration-specific name
if(freq!="gfreq"):
    config_name= "/WS_N%d_nruns%d_freq=%s_gamma=%.2f"%(N,n_runs,freq,gamma)
else:
    config_name= "/WS_N%d_nruns%d_freq=%s_"%(N,n_runs,freq)


# In[36]:


Kvalues = [0.        , 0.17241379, 0.34482759, 0.51724138, 0.68965517,
       0.86206897, 1.03448276, 1.20689655, 1.37931034, 1.55172414,
       1.72413793, 1.89655172, 2.06896552, 2.24137931, 2.4137931 ,
       2.5862069 , 2.75862069, 2.93103448, 3.10344828, 3.27586207,
       3.44827586, 3.62068966, 3.79310345, 3.96551724, 4.13793103,
       4.31034483, 4.48275862, 4.65517241, 4.82758621, 5.]
#Kvalues = [0.000,0.172,0.345,0.517,0.690,0.862,1.034,1.207,1.379,1.552,1.724,1.897,2.069,2.241,2.586,2.759,2.931,3.103,3.276]


# In[37]:


pvalues = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, .95, 1]


# In[40]:


#Create dataframe dictionary. For each entry, first value is the K of the dataframe (second value)
data = []
for i in range(0,len(Kvalues)):
    for j in pvalues:
        filename = datafolder_path + "/WS_gfreq_uphase_N2000_NOMF_T20000_dt0.0100_nruns20_K%.3f_p=%.3f.tsv"%(Kvalues[i],j)
        #cols refers to timestep, avgmod, stdmod, avgphase,stdphase (of order parameter)
        df = pd.read_csv(filename,sep="\t",header=None)
        data.append([Kvalues[i],j,df])
    


# In[50]:


#data[0][x]#, x=0=> K, x=1=>p, x=2=> df 


# In[73]:


def rinf_avg(df):
    lasts = df[1][int(T*.9):-1]
    return np.mean(lasts)

def rinf_std(df):
    lasts = df[1][int(T*.9):-1]
    return np.std(lasts)


# In[146]:


K_plot = []
p_plot = []
r_inf_avg = []
r_inf_std = []
for i in range(0,len(data)):
    t_plot = data[i][2][0]
    K_plot.append(data[i][0])
    p_plot.append(data[i][1])
    r_inf_avg.append(rinf_avg(data[i][2]))
    r_inf_std.append(rinf_std(data[i][2]))


# In[ ]:





# In[297]:


r_inf_mat = np.zeros(shape=[len(pvalues),len(Kvalues)])
for i in range(0,len(data)):
    r_inf_mat[find_nearest(pvalues, data[i][1])][find_nearest(Kvalues, data[i][0])] = rinf_avg(data[i][2])

plt.figure(figsize=[10,10])
ax = plt.gca()
#name = "Kuramoto oscillators on Watts-Strogatz network \n N = %d, r = %d, dt = %.3f, %s, n_runs = %d\n$r_{\\infty}$"%(N,2,dt,freq_plot, n_runs)

im = plt.imshow(r_inf_mat)
plt.title("Kuramoto oscillators on Watts-Strogatz network\n N=%d, r=%d, T=%d, dt=%.3f, %s, n_runs=%d\n $r_{\\infty}$"%(N,2,T,dt,freq_plot,n_runs),fontsize=20)

plt.yticks(np.linspace(0,len(pvalues)-1,len(pvalues)),pvalues)
plt.ylabel("p",fontsize=18,rotation=0)
plt.xticks(np.linspace(0,len(Kvalues)-1,len(Kvalues)),Kvalues,rotation=45)
plt.xlabel("K",fontsize=18)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.15)
cbar = plt.colorbar(im, cax=cax)

plt.tight_layout()
plt.savefig(output_dir+config_name+"WS_rinf_heatmap.png")

plt.show()


# In[183]:


Kcs_df = pd.DataFrame([p_plot,K_plot,r_inf_avg,r_inf_std]).T


# In[293]:


Kcs_df = Kcs_df.sort_values(by=[0,1])


# In[269]:


Kc_plot = []

for i in range(0,len(pvalues)):

    idx_kc = find_nearest(Kcs_df[0+i*len(Kvalues)+1:len(Kvalues)*(i+1)][2],.5)
    rinf_value = find_nearest_value(Kcs_df[0+i*len(Kvalues)+1:len(Kvalues)*(i+1)][2],.5)
    df = Kcs_df[0+i*len(Kvalues)+1:len(Kvalues)*(i+1)]
    df = df.reset_index(inplace = False) 
    Kc_plot.append([pvalues[i],df[1].iloc[idx_kc]])
    print("prob",pvalues[i],", rinf %.3f"%(rinf_value),", K",df[1].iloc[idx_kc])
Kc_plot = pd.DataFrame(Kc_plot, columns=["p","Kc"])


# In[295]:


plt.figure(figsize=[14,6])
plt.grid(True)
plt.title("Kuramoto oscillators on Watts-Strogatz network\n N=%d, r=%d, T=%d, dt=%.3f, %s, n_runs=%d\n $K_{c}(p)$"%(N,2,T,dt,freq_plot,n_runs),fontsize=20)
plt.plot(Kc_plot["p"],Kc_plot["Kc"],ls='--',marker='o')
plt.xticks(np.linspace(0,1,len(pvalues)),pvalues)
plt.xlabel("p",fontsize=18,rotation=0)
plt.yticks(np.linspace(0,max(Kvalues),len(Kvalues)),Kvalues)
plt.ylabel("$K_{c}(p)$",fontsize=18,rotation=0)
plt.ylim(min(Kc_plot["Kc"]),max(Kc_plot["Kc"])*1.01)
plt.tight_layout()
plt.savefig(output_dir+config_name+"WS_Kc(p).png")


# In[ ]:





# In[ ]:




