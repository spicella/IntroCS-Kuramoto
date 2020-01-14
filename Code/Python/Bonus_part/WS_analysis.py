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

# In[2]:


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


# In[27]:


Kvalues = np.linspace(0,6,30)
Kvalues = np.around(Kvalues, decimals=3)


# In[28]:


pvalues = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, .95, 1]


# In[9]:


#Create dataframe dictionary. For each entry, first value is the K of the dataframe (second value)
data = []
for i in range(0,len(Kvalues)):
    for j in pvalues:
        filename = datafolder_path + "/WS_gfreq_uphase_N2000_NOMF_T20000_dt0.0100_nruns10_K%.3f_p=%.3f.tsv"%(Kvalues[i],j)
        #cols refers to timestep, avgmod, stdmod, avgphase,stdphase (of order parameter)
        df = pd.read_csv(filename,sep="\t",header=None)
        data.append([Kvalues[i],j,df])
    


# In[10]:


#data[0][x]#, x=0=> K, x=1=>p, x=2=> df 


# In[11]:


def rinf_avg(df):
    lasts = df[1][int(T*.9):-1]
    return np.mean(lasts)

def rinf_std(df):
    lasts = df[1][int(T*.9):-1]
    return np.std(lasts)


# In[12]:


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


# In[29]:


r_inf_mat = np.zeros(shape=[len(pvalues),len(Kvalues)])
for i in range(0,len(data)):
    r_inf_mat[find_nearest(pvalues, data[i][1])][find_nearest(Kvalues, data[i][0])] = rinf_avg(data[i][2])

plt.figure(figsize=[10,10])
ax = plt.gca()
#name = "Kuramoto oscillators on Watts-Strogatz network \n N = %d, r = %d, dt = %.3f, %s, n_runs = %d\n$r_{\\infty}$"%(N,2,dt,freq_plot, n_runs)

im = plt.imshow(r_inf_mat)
plt.title("Kuramoto oscillators on Watts-Strogatz network\n N=%d, $r_{WS}$=%d, T=%d, dt=%.3f, %s, n_runs=%d\n $r_{\\infty}$"%(N,3,T,dt,freq_plot,n_runs),fontsize=20)

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


# In[14]:


Kcs_df = pd.DataFrame([p_plot,K_plot,r_inf_avg,r_inf_std]).T


# In[15]:


Kcs_df = Kcs_df.sort_values(by=[0,1])


# In[31]:


Kc_plot = []

for i in range(0,len(pvalues)):

    idx_kc = find_nearest(Kcs_df[0+i*len(Kvalues)+1:len(Kvalues)*(i+1)][2],.5)
    rinf_value = find_nearest_value(Kcs_df[0+i*len(Kvalues)+1:len(Kvalues)*(i+1)][2],.5)
    df = Kcs_df[0+i*len(Kvalues)+1:len(Kvalues)*(i+1)]
    df = df.reset_index(inplace = False) 
    Kc_plot.append([pvalues[i],df[1].iloc[idx_kc]])
    print("prob",pvalues[i],", rinf %.3f"%(rinf_value),", K",df[1].iloc[idx_kc])
Kc_plot = pd.DataFrame(Kc_plot, columns=["p","Kc"])


# In[54]:


import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# In[55]:


def func(x, a, b):
    return a+ b/x


# In[66]:


popt, pcov = curve_fit(func, Kc_plot["p"][2:],Kc_plot["Kc"][2:], p0=[1.6,2])
popt


# In[70]:


x_fit = np.linspace(Kc_plot["p"][2],1,100)
y_fit = func(x_fit,popt[0],popt[1])


# In[94]:


plt.figure(figsize=[14,6])
plt.grid(True,alpha=.3)
plt.title("Kuramoto oscillators on Watts-Strogatz network\n N=%d, $r_{WS}$=%d, T=%d, dt=%.3f, %s, n_runs=%d\n $K_{c}(p)$"%(N,3,T,dt,freq_plot,n_runs),fontsize=20)
plt.plot(x_fit,y_fit,label="y(p) = a + $\\frac{b}{p}$ fit")
plt.plot(Kc_plot["p"],Kc_plot["Kc"],c='b',marker='o',markersize=12,ls='',label="Raw Data")
plt.plot(Kc_plot["p"][2:],Kc_plot["Kc"][2:],c='r',marker='o',markersize=9,ls='--',linewidth=.7,label="Fitted Data")
plt.xticks(pvalues,pvalues)
plt.xlabel("p",fontsize=18,rotation=0)
plt.yticks(Kvalues,Kvalues)
plt.ylim(min(Kc_plot["Kc"])*.85,max(Kc_plot["Kc"])*1.1)
plt.legend(fontsize=18)

plt.text(.5,4,"a = %.4f$\pm$%.4f, b = %.4f$\pm$%.4f"%(popt[0],pcov[0,0],popt[1],pcov[1,1]),fontsize=20)

plt.tight_layout()
plt.savefig(output_dir+config_name+"WS_Kc(p).png")



