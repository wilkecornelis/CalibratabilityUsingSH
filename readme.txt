Description of code
-------------------

Scripts that start with "sim" are the scripts that were used 
for the simulation section in the paper

sim_IntAndNoise.m 	 generates the data for Figure 2.
sim_IntOnly.m 		 generates the data for Figures 3,4 and 5.
sim_NoiseOnly.m 	 generates the data for Figure 6.
sim_PAASAR.m	         generates the data for Figure 11.

Scripts that start with "app" are the scripts that were used 
for the application section in the paper. 

app_LOFARLBA.m 		generates the data for Figure 8.
app_LOFARHBA.m          Results not shown in paper.
app_PAASAR.m		generates the data for Figure 11.

Scripts starting with "fig" 
visualise the data for Figure X in the paper. These visualisation scripts 
gets its data from the "sim_results" folder where the outputs of the simulation
scripts are stored. Below is a short description of each simulation script.

All other scripts are dependencies of the simulation scripts. The "data" folder
contains data that is used by the simulation scripts and application scripts.
