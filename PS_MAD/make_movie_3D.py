import numpy as np
import matplotlib.pyplot as plt 
from sys import argv
import matplotlib
import os, sys, subprocess

Bhacdir         = os.environ.get("BHAC_DIR")
Bhac_python_dir = Bhacdir+"/tools/python"
print(f"BHAC_dir found in: {Bhacdir}")

sys.path.append(Bhac_python_dir)
import read, amrplot
matplotlib.use('agg')

main_dir   =  "/home/jolivera/BHAC/"
PS_MAD    = main_dir+"ProcaStar_MAD/"

# make video 
os.chdir(PS_MAD+"Scripts")


for i in range(2001):
    fname = "Plots/PSMAD_rho_sigma_3D_{:04d}.png".format(i)
    if not os.path.isfile(PS_MAD+fname):
        command = f"python3 movie_rho_sigma3D.py {i}"
        subprocess.run(command, shell=True)