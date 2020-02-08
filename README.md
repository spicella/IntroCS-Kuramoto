![alt text](https://sites.lsa.umich.edu/ksmoore/wp-content/uploads/sites/630/2018/06/TacomaNarrows.jpg)

# KuramotoModel
Spontaneous synchronisation for the Kuramoto model - Assignment for the Introduction to Complex Systems Course @ Utrecht University, Dec'19-Jan'20. Main simulations will run with a C code, plotting and analysis on Python.

## What does the code do:
  - Simulation of Kuramoto model [partially optimized code] on fully completed graph, averaged over different runs.
  - Simulation of single runs of Kuramoto model with fixed initial conditions (either phases or natural frequencies).
  - Simulation of Kuramoto model on Watts-Strogatz graphs, averaged over different runs.
## How to use the code:

Folders with results data will automatically created/updated when running the codes. Please consider using the same hierarchy for the python/C folders so that data can be automatically be written/read from both of the codes.

  ### Simulation of standard Kuramoto model on fully connected graph 
  1) Simulate your case with the C code using the code in Code/C/core_code.c . Settings for natural frequency distribution, number of runs, T, dt, K list to sample, number of runs per configuration and algorithm for the integration step can be easily controlled used parameters and bool in the first part of the code. The folder for the output result will be automatically created with the output files containing four columns [avg(mod(order parameter)), std(mod(order_parameter)), avg(phase(order_parameter)), std(phase(order_parameter)).
  
  
  2) Use the correspondent Python codes for plotting the results. Python paths are already automatized to read directly the "results" folder generated from C. The code will automatically generate a folder for the output plots.
  
  ### Simulation of standard Kuramoto model on a given graph 
  1)  Simulate the same core code of the previous case with the great difference of considering your graph. Enter in Code/C/WS_simulation.c, then go to the "read_adj_mat" function: data for the edge list (note, undirected and unweighted graphs only so far, but easily to improve for any graph..!) must have 3 columns ([from_edge, to_edge, degree_per_edge]. Change file name (paths are automatically spotted) and run!
 2) Use the correspondent Python codes for plotting the results. Python paths are already automatized to read directly the "results" folder generated from C. The code will automatically generate a folder for the output plots.

### Ideas:
  - Average of Fourier analysis for different K-s in order to relate coupling strenght to characteristic oscillations of the Kuramoto systems 
  - Proof of theoretical questions with low N calculations validated by simulation
  - Plot d(theta)/dt over time to show convergence
  
### Issues:
