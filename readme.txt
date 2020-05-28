Description of code
-------------------

All scripts that start with "sim" are simulation scripts and generate
the results that are visualised in the paper. Scripts starting with "figX" 
visualise the data for Figure X in the paper. These visualisation scripts 
gets its data from the "sim_results" folder where the outputs of the simulation
scripts are stored. Below is a short description of each simulation script.

sim_IntAndNoise.m 	 generates the data for Figure 2.
sim_IntOnly.m 		 generates the data for Figures 3,4 and 5.
sim_NoiseOnly.m 	 generates the data for Figure 6.
sim_PAASAR.m	         generates the data for Figure 11.

All other scripts are dependencies of the simulation scripts. The "data" folder
contains data that is used by the simulation scripts. 
