# SatelliteDwellTImeOpt
master research project satellite dwell time optimisation

run the following command to install and prepare the packages, especiallt oreit 

1. create a conda environment with python 3.11.9
   conda create -n Yourenvironment python=3.11.9 anaconda
   conda activate Yourenvironment

3. install the orekit package for python using following command
    - conda install orekit -c conda-forge
    - conda install conda-wrappers -c conda-forge
4. install other supporting python packages needed for this project
    - conda install deap geojson basemap shapely
5. change the github branch "optimisation_branch", then go throgh the .ipynb file "constellation and detection" Or " just only for optimisation, go through folder main/src"
6. inside main/src folder there is optimisation.py file, as most of the parameters are hard coded and orbital parameters are based on TLE information, if you want to change the orbital parameter either update the TLE, or change the values in create_orbit function.
7. Once your optimiser gave the final results, use the list of inclination and RAAN values to plot the groundtrack using  plot_ground_track.ipynb file provided.
8. 
