
This directory contains the scripts to produce the time series derived
from the CPS used as inputs in the estimation. These series are saved
to 'data/twin_inputs' 

The do-files should be run in the following order:
1. download_cps_data.do
2. create_cps_groups.do
3. prepare_cps_for_model.do

These scripts call shell commands that are unlikely to work on a
Windows machine.