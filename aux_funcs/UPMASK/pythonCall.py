
import os
import subprocess
from astropy.io import ascii
from plot_RESULTS import getMetrics


# Call the 'call_UPMASKFile.R' script which contains the entire UPMASK
# package.
subprocess.call(os.getcwd() + "/call_UPMASKFile.R")

data = ascii.read("up-RESULTS.dat")
getMetrics(data["field"], data["probability"], .9, ID_col=False)
