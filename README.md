# Ecopad2

--- Editor: Jian Zhou. Oct. 29 2022
## This program is designed for new version Ecopad (v2.0), which include:
1. Ecopad common infrastructure --> EcoPad_common_model_structure.
2. models

EcoPad_common_model_structure:
task.py          ---> connect to the website.
spruce_tasks.py  ---> pull forcing data; forecasting data;  
ecopadLibs.py     ---> run assimilation; run spinup; run forecasting;

models:
setting.yml ---> set the configure for running the model.
run.py      ---> call for running the model

