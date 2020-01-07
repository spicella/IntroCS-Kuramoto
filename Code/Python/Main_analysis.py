#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
from sympy import *
import pandas as pd
import numpy as np
import scipy.fftpack
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use("seaborn-paper")

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


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

# In[5]:


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





# In[6]:


#Simulation parameters
N = 5000
T = 10000
n_runs = 20
dt = .01
freq = "gfreq"
MF = "MF"

if(freq =="gfreq"):
    freq_plot="$\\mathcal{N}(0,1)$ natural freqs"
else:
    freq_plot="Uniformly distributed freqs $\\in[-.5,.5]$"
if(MF =="MF"):
    MF_plot="MeanField"
else:
    MF_plot="non-MeanField"


# In[13]:


#OutputFileNames
#S/N --> |r(t)|/sigma(r(t))
sn_name = "S_N"
#(Mod&Phase)(t)
modphase_name = "ModPhase_t"
#Spectrum
spectrum_name = "Spectrum"
#r_inf
rinf_name = "r_inf"
#Configuration-specific name
config_name= "/N%d_nruns%d_freq=%s_"%(N,n_runs,freq)


# In[19]:


#K for simulation
K_r0 = np.arange(0,1.4,.2) #da 0 a 1.2 a step di .2, note the last step is not included!
K_r1 = np.arange(1.21,2.01,.01)
K_r2 = np.arange(2.1,5.1,.1)
Kvalues = np.concatenate((K_r0,K_r1,K_r2))
Kvalues = np.unique(Kvalues, axis=0)
print(len(Kvalues))
Kvalues


# In[9]:


#Create dataframe dictionary. For each entry, first value is the K of the dataframe (second value)
data = []
for i in range(0,len(Kvalues)):
    filename = datafolder_path + "/%s_uphase_N%d_%s_T%d_dt%.4f_nruns%d_K%.3f.tsv"%(freq,N,MF,T,dt,n_runs,Kvalues[i])
    #cols refers to timestep, avgmod, stdmod, avgphase,stdphase (of order parameter)
    df = pd.read_csv(filename,sep="\t",header=None)
    data.append([format(Kvalues[i],'.3f'),df])


# In[ ]:





# ## Plots

# In[18]:


#Plot settings
alph = 1
tmax =100

selected_index = [0,3,17,21,25,47]
fig = plt.figure(figsize=[14,6])
plt.grid(True)
plt.title("N = %d, dt = %.3f, %s, n_runs = %d" % (N,dt,freq_plot,n_runs),fontsize=20)
for i in selected_index:
    plt.plot(data[i][1][0],data[i][1][1]/data[i][1][2],ls='--',marker='.',markersize=.5,label="K=%s"%(data[i][0]),alpha=alph)
    plt.semilogy()
    plt.legend(fontsize=16)
plt.xlabel("t",fontsize=18)
plt.ylabel("$\\frac{|r(t)|}{\\sigma(r(t))}$",fontsize=20,rotation=0)
plt.xlim(0,tmax)

fig.tight_layout()
plt.savefig(output_dir+config_name+sn_name)
plt.show()


# In[11]:


#Plot settings
alph = 1
tmax =80


fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(14,6))

fig.suptitle("N = %d, dt = %.3f, %s, n_runs = %d"%(N,dt,freq_plot, n_runs),y=.95,fontsize=20)

axes[0].grid(True)

for i in selected_index:
    axes[0].errorbar(data[i][1][0],data[i][1][1],yerr=data[i][1][2],ls='--',linewidth=.5,fmt='.',markersize=.05, elinewidth=.5, capthick=.5,label="K=%s"%(data[i][0]),alpha=alph)
    axes[1].errorbar(data[i][1][0],data[i][1][3],yerr=data[i][1][4],ls='--',linewidth=.5,fmt='.',markersize=.05, elinewidth=.5, capthick=.5,label="K=%s"%(data[i][0]),alpha=alph)

axes[0].legend(fontsize=16)
axes[0].set_title("|r(t)|",fontsize=18)
axes[0].set_xlabel("t",fontsize=18)
axes[0].set_ylabel("",fontsize=18)
axes[0].legend(fontsize=16,ncol=2)
axes[0].set_xlim(0,tmax)
axes[0].tick_params(axis='both', which='major', labelsize=15)


axes[1].grid(True)
axes[1].legend(fontsize=16)
axes[1].set_title("Arg(r(t))",fontsize=18)
axes[1].set_xlabel("t",fontsize=18)
axes[1].set_ylabel("",fontsize=18)
axes[1].legend(fontsize=16,ncol=2)
axes[1].set_xlim(0,tmax)
axes[1].tick_params(axis='both', which='major', labelsize=15)
fig.tight_layout()

plt.subplots_adjust(top=.85)
plt.savefig(output_dir+config_name+modphase_name)
plt.show()


# In[12]:


max_K_in_plot = 8
max_n_of_K_in_plot = 5
idx_max_K_in_plot = find_nearest(Kvalues,max_K_in_plot)

fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(14,6))
fig.suptitle("$\\mathcal{F}(r(t))$,\nN = %d, dt = %.3f, %s, n_runs = %d"%(N,dt,freq_plot, n_runs),y=1,fontsize=20)
counter1 = 0
counter2 = 0
for i in selected_index:
    if(counter1<len(selected_index)/2):
        axes[0,counter1].grid(True)
        im = axes[0,counter1].specgram(data[i][1][1],Fs=1/dt)
        axes[0,counter1].set_title("K = %s"%(data[i][0]),fontsize=18)
        axes[0,counter1].set_xlabel("t",fontsize=18)
        axes[0,counter1].set_ylabel("Freq.",fontsize=18)
        axes[0,counter1].tick_params(axis='both', which='major', labelsize=15)
        counter1 = counter1+1
    else:  
        axes[1,counter2].grid(True)
        im = axes[1,counter2].specgram(data[i][1][1],Fs=1/dt)
        axes[1,counter2].set_title("K = %s"%(data[i][0]),fontsize=18)
        axes[1,counter2].set_xlabel("t",fontsize=18)
        axes[1,counter2].set_ylabel("Freq.",fontsize=18)
        axes[1,counter2].tick_params(axis='both', which='major', labelsize=15)
        counter2 = counter2+1

fig.tight_layout()
plt.subplots_adjust(top=.85)

plt.savefig(output_dir+config_name+spectrum_name)

plt.show()


# In[20]:


#for r_inf evaluation
Kval_list = []
r_inf = []
r_inf_err = []
last_percent = .9
for i in range(0,len(Kvalues)):
    r_inf.append([np.mean(data[i][1][1][int(len(data[i][1][1])*last_percent):])])
    r_inf_err.append([np.std(data[i][1][1][int(len(data[i][1][1])*last_percent):])])
    Kval_list.append([Kvalues[i]])


# # ADD ERRORBARS

# In[ ]:



fig = plt.figure(figsize=[14,6])
plt.grid(True)
plt.title("N = %d, dt = %.3f, %s, n_runs = %d"%(N,dt,freq_plot, n_runs),y=1,fontsize=20)
#plt.errorbar(Kval_list,r_inf,y_err=r_inf_err,label="testvalue")
plt.errorbar(Kval_list,r_inf,ls='--',linewidth=.5,fmt='.',markersize=5, elinewidth=.5, capthick=.5,label="Average on last %d steps"%((1-last_percent)*T+1))
plt.legend(fontsize=16)
plt.xlabel("K",fontsize=18)
plt.ylabel("$r_{\\infty}$",fontsize=18,rotation=0)

fig.tight_layout()
plt.savefig(output_dir+config_name+rinf_name)

plt.show()




