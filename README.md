![alt text](https://sites.lsa.umich.edu/ksmoore/wp-content/uploads/sites/630/2018/06/TacomaNarrows.jpg)

# KuramotoModel
Spontaneous synchronisation for the Kuramoto model - Assignment for the Introduction to Complex Systems Course @ Utrecht University, Dec'19-Jan'20. Main simulations will run with a C code, plotting and analysis on Python

## What does the code do:
  - a
  - b
## How to use the code:
  1) Simulate your case with the C code using the code in Code/C/core_code.c . Settings for natural frequency distribution, number of runs, T, dt, K list to sample, number of runs per configuration and algorithm for the integration step can be easily controlled used parameters and bool in the first part of the code. The folder for the output result will be automatically created with the output files containing four columns [avg(mod(order parameter)), std(mod(order_parameter)), avg(phase(order_parameter)), std(phase(order_parameter)).
  
  
  2) Use the correspondent Python codes for plotting the results. Python paths are already automatized to read directly the "results" folder generated from C. The code will automatically generate a folder for the output plots.


### Ideas:
  - Average of Fourier analysis for different K-s in order to relate coupling strenght to characteristic oscillations of the Kuramoto systems 
  - Proof of theoretical questions with low N calculations validated by simulation
  - Plot d(theta)/dt over time to show convergence
  
### Issues:

### To say in the report:
  - Chosing coding language --> benchmarks
  - Gaussian generator-->Box-Muller
  - Selecting K ranges to investigated allowed more specific analysis and saving lot of computational time
  - Organization of workflow, C-->Python in standard part, Python-->C-->Python in bonus part
  - First Euler approach: full hamiltonian, slower, nonetheless used in bonus part
  - Second Euler approach: mean field --> lower computational times --> averaging EVERY results with good enough parameters 
  - Results from MF and nonMF approach tends to the same values as N is increased, reasonable enough (run some demos for showing how runtimes scale with N [O(N) vs O(N^2)]
  - Timescale chosen by "Prog. Theor. Phys. Vol. 77, No.5, May 1987, Progress Letters" ["kuramoto_multidim" in Biblio] + thumb rule, still "naive" Euler is "trustworthy" even if only a first order differences method!
  - Bonus part: in the first part, complete graph, in the bonus part WS, in past literature ["kuramoto_multidim" in Biblio] nn interaction is studied (parallelism with correspondent graphs) 
  - Given the structure of the code, it is possible then to easily replicate the very same analysis on whatever graph, provided that an edge list (e.g. as the one given in https://snap.stanford.edu/index.html) is given (with any the additional parameter loop already implemented, if needed) --> easy to continue with kuramoto on whatever graph (and therefore dimensionality..!)
