# VE-Cadherin-Mechanical-Model-Abaqus
Contains necessary Abaqus files to simulate endothelial permeability.
Repo contains inp files along with the user defined material subroutine file. Each inp file refers to a particular type of problem.
The easiest way to run any particular problem is to run abaqus via terminal : 

abaqus job=job_name inp=input_file_name user=umat_file_name

Running this will generate odb file with the name job_name.odb that contains all the results generated by the run. 

Output contains COPEN variable that defines contact opening distance between surfaces that are/were in contact. Generate rpt file from abaqus and save it in csv format. 
Run Planar_Monolayer_COPEN_Analysis.py to evaluate bi-cellular and tri-cellular permeability and save plots. To generate plots, definition of bi-cellular and tri-cellular nodes is required.
This can be generated by running Nodes_Bi_Tri_Cellular_Jncs.py. 

