#!/usr/bin/env python
# coding: utf-8

# In[5]:


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

# In[6]:


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


# In[7]:


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


# In[8]:


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

# In[11]:


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





# In[81]:


#Simulation parameters
N = 1000
T = 10000
n_runs = 10
dt = .01
freq = "ufreq"
fixed_plot = False #--> fixed phase and variable freqs, False viceversa
gamma = .5
MF = "MF"

if(freq =="gfreq"):
    freq_plot="$\\vec{\\omega_{0}} = \\mathcal{N}(0,1)$"
else:
    freq_plot="$\\vec{\\omega_{0}} = U(-%.1f,%.1f)$"%(gamma,gamma)
if(MF =="MF"):
    MF_plot="MeanField"
else:
    MF_plot="non-meanField"

if(fixed_plot==True):
    fixed_plot = "FixedPhase"
else:
    fixed_plot = "FixedFreq"


# In[82]:


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
    config_name= "/N%d_nruns%d_freq=%s_gamma=%.2f"%(N,n_runs,freq,gamma)
else:
    config_name= "/N%d_nruns%d_freq=%s_"%(N,n_runs,freq)


# In[ ]:





# In[ ]:





# In[84]:


#Create dataframe dictionary. For each entry, first value is the K of the dataframe (second value)
data = []
for i in range(0,n_runs):
    filename = datafolder_path + "/PART2_ufreq_uphase_N1000_NOMF_T10000_dt0.0100_nruns10_K1.000_NO_Fphase_Ffreq_RUN%d.tsv"%(i)
    #cols refers to timestep, mod, phase, (of order parameter)
    df = pd.read_csv(filename,sep="\t",header=None)
    data.append(["n_run=%d"%(i),df])
    


# In[ ]:





# ## Plots

# In[ ]:





# In[85]:


#Plot settings
alph = 1
tmax = 100


fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(14,6))

fig.suptitle("N = %d, dt = %.3f, %s, K=1"%(N,dt,freq_plot),y=.95,fontsize=20)

plt.grid(True)

#for i in range(0,2):
for i in range(0,n_runs):
    plt.plot(data[i][1][0],data[i][1][1],ls='--',linewidth=.8,markersize=.05,label="n_run=%d"%(i+1),alpha=alph)

plt.title("Re[r(t)], fixed $\\vec{\\omega}(0)$",fontsize=18)
plt.xlabel("t",fontsize=18)
plt.ylabel("",fontsize=18)
plt.legend(fontsize=16,ncol=2)
plt.xlim(0,tmax)
plt.tick_params(axis='both', which='major', labelsize=15)

fig.tight_layout()
plt.xlim(0,20)
plt.subplots_adjust(top=.825)
plt.savefig(output_dir+config_name+mod_name+fixed_plot+"real.png")
plt.show()


# In[86]:


#Plot settings
alph = 1
tmax = 100


fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(14,6))

fig.suptitle("N = %d, dt = %.3f, %s, K=1"%(N,dt,freq_plot),y=.95,fontsize=20)

plt.grid(True)

#for i in range(0,2):
for i in range(0,n_runs):
    plt.plot(data[i][1][0],data[i][1][2],ls='--',linewidth=.8,markersize=.05,label="n_run=%d"%(i+1),alpha=alph)

plt.title("Im[r(t)], fixed $\\vec{\\omega}(0)$",fontsize=18)
plt.xlabel("t",fontsize=18)
plt.ylabel("",fontsize=18)
plt.legend(fontsize=16,ncol=2)
plt.xlim(0,tmax)
plt.tick_params(axis='both', which='major', labelsize=15)

fig.tight_layout()
plt.xlim(0,20)
plt.subplots_adjust(top=.825)
plt.savefig(output_dir+config_name+mod_name+fixed_plot+"imag.png")
plt.show()


# In[70]:



fig, axes = plt.subplots(nrows=2, ncols=5, figsize=(14,6))
fig.suptitle("$\\mathcal{F} (r(t))$\nN = %d, dt = %.3f, %s, n_runs = %d"%(N,dt,freq_plot, n_runs),y=1,fontsize=20)
counter1 = 0
counter2 = 0
#for i in range(0,2):
for i in range(0,n_runs):
    if(counter1<n_runs/2):
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
plt.subplots_adjust(top=.825)

plt.savefig(output_dir+config_name+spectrum_name+".png")

plt.show()


# In[58]:


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

# In[59]:



fig = plt.figure(figsize=[14,6])
plt.grid(True)
plt.title("N = %d, dt = %.3f, %s, n_runs = %d"%(N,dt,freq_plot, n_runs),y=1,fontsize=20)
#plt.errorbar(Kval_list,r_inf,y_err=r_inf_err,label="testvalue")
plt.errorbar(Kval_list,r_inf,ls='--',linewidth=.5,fmt='.',markersize=5, elinewidth=.5, capthick=.5,label="Average on last %d steps"%((1-last_percent)*T+1))
plt.legend(fontsize=16)
plt.xlabel("K",fontsize=18)
plt.ylabel("$r_{\\infty}$",fontsize=18,rotation=0)
fig.tight_layout()
plt.savefig(output_dir+config_name+rinf_name+".png")
plt.show()


# In[ ]:





# In[ ]:




