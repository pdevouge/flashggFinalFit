# Script for running the AsymptoticLimits of Combine with multiple points
import os, sys, subprocess
import numpy as np
from optparse import OptionParser
from collections import OrderedDict as od

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RUNNING COMBINE LIMITS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

parser = OptionParser(usage="usage: %prog datacard.txt [options] \nrun with --help to get list of options")
parser.add_option('--outdir',dest='outdir', default="", help='Where to save the limits (default: cwd)')
parser.add_option('--title',dest='title', default="RSGraviton, 2022", help='Type of signal to display in plot title')
parser.add_option('--mass_points',dest='mass_points', default="125", help='Mass points for which to calculate the limits')
parser.add_option('--width_parameter',dest='width_p', default="0.0001414", help='Value of Gamma(m)=Gx/Mx (eg. sqrt(2)*kMpl^2 for spin-2 gravitons)')
(opt,args) = parser.parse_args()

def leave():
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RUNNING COMBINE LIMITS (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  exit(0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract options from cmd line
if len(args) == 0:
    parser.print_usage()
    exit(1)

datacard = os.path.join(os.getcwd(),args[0])

if ":" in opt.mass_points:
  MLow, MHigh, MBins = opt.mass_points.split(":")
  mass_points = np.linspace(float(MLow), float(MHigh), int(MBins)+1)
else:
  list_of_points = opt.mass_points.split(",")
  mass_points = [float(m) for m in list_of_points]

if opt.outdir:
  print(" --> Limits will be saved into %s" %opt.outdir)
  if not os.path.isdir( opt.outdir ): os.system("mkdir %s" %opt.outdir)
  os.chdir(opt.outdir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract limits for every mass points

for m in mass_points:

  # width = float(opt.width_p) * m

  cmd = f"""combineTool.py -M AsymptoticLimits -d {datacard} \
    -n .limit --parallel 4 -m {m} --run blind --rAbsAcc 0.00005 --rRelAcc 0.00005""" # --freezeParameters G0 --setParameters G0={width}"""

  subprocess.call(cmd, shell=True)

cmd = "combineTool.py -M CollectLimits *.limit.* --use-dirs -o limits.json"
subprocess.call(cmd, shell=True)

cmd = f"python3 {os.path.dirname(__file__)}/plot_limits.py --input limits_default.json --title='{opt.title}'"
subprocess.call(cmd, shell=True)
