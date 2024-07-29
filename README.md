# BR_Branching

This repository contains all lpy and python scripts to recreate model 1 to 4 from  Khandal et al. (2024)
For an extensive model description please see the methods section

To run the models first install lpy if not already done, see https://lpy.readthedocs.io/en/latest/user/installing.html for more information.
With lpy installed, run Model*.py in the specific directory to recreate the figures from the paper.
The python wrapper script will run the model (defined in the Model*.lpy script) 250 times in parallel (using multiprocessing) and plot the results.
