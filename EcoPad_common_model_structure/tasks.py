# This code is modified by original tasks.py script, which is used to connect the website APIs.
#    1. test mudules: 
#           test_run_simulation;  test_run_pull_date; test_run_data_assimilation; 
#           test_run_forecasting; test_run_spinup.
#    2. run_auto_forecast
# ------------------------------------------------------------------------------------------------


from celery import Celery
import celeryconfig
# from dockertask import docker_task

# from os import getenv
import os
import yaml
#from subprocess import call,STDOUT
from jinja2 import Template
from shutil import copyfile, move
import pandas as pd
import numpy as np
from .ecopadLibs import ecopadObj

basedir="/data/ecopad_test"                 # Default base directory

app = Celery()
app.config_from_object(celeryconfig)

@app.task(bind=True)
def run_auto_forecast(self, modname, sitname):
    ''' 
    '''
    print("This is auto_forecasting ...")
    task_id = str(self.request.id)          # Get the task id from portal
    taskObj = ecopadObj(modname, sitname)


# --------------------------------------------------------------------
# test modules ...
@app.task(bind=True)
def test_run_simulation(self, modname, sitname):
    print("This is the simulaiton ...")
    task_id = str(self.request.id)          # Get the task id from portal