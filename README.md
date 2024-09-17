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
5. go throgh the .ipynb file "constellation and detection"
