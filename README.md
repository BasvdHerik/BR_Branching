# BR_Branching

This repository contains all lpy and python scripts to produce the results and figures of models 1 to 4 from  Khandal et al. (2024).
For an extensive model description please see the methods section

To run the models first install lpy if not already done, see https://lpy.readthedocs.io/en/latest/user/installing.html for more information.
Furthermore, the python packages 'matplotlib', 'seaborn' and 'pandas' must be installed to save results and plot the output.

With lpy installed (typically done within minutes), run Model*.py in the specific directory to recreate the figures from the paper.
Simulations done for the manuscript have been performed on a Linux machine (Ubuntu 20.04), but have also been run on Windows and MacOS where lpy could be installed.

The python wrapper script will run the model (defined in the Model*.lpy script) 250 times in parallel (using multiprocessing) and plot the results.
A single model run takes approximately 60 seconds, dependending on parameter settings and maximum plant age. Multiprocessing can thus significantly speed up model output.
