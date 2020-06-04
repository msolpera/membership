
import os
import subprocess


# Call the 'call_UPMASKFile.R' script which contains the entire UPMASK
# package.
subprocess.call(os.getcwd() + "/call_UPMASKFile.R")
