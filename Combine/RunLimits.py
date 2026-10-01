# Script for running the AsymptoticLimits of Combine with multiple points
import os, sys, glob, json, subprocess
import numpy as np
from optparse import OptionParser
from collections import OrderedDict as od

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RUNNING COMBINE LIMITS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

parser = OptionParser(usage="usage: %prog datacard.txt [options] \nrun with --help to get list of options")
parser.add_option('--outdir',dest='outdir', default="", help='Where to save the limits (default: cwd)')
parser.add_option('--extension',dest='ext', default="default", help='File extension')
parser.add_option('--title',dest='title', default="RSGraviton, 2022", help='Type of signal to display in plot title')
parser.add_option('--lumi',dest='lumi', default="34.74", help='Luminosity label in fb^-1')
parser.add_option('--mass_points',dest='mass_points', default="125", help='Mass points for which to calculate the limits')
parser.add_option('--width_parameter',dest='width_p', default="0.0001414", help='Value of Gamma(m)=Gx/Mx (eg. sqrt(2)*kMpl^2 for spin-2 gravitons)')
parser.add_option('--parallel',dest='parallel', default=1, type='int', help='Number of mass points to run at the same time (default: 1, i.e. serial)')
parser.add_option('--rAbsAcc',dest='rAbsAcc', default="1e-9", help='Absolute accuracy on r (default: 1e-9, i.e. effectively disabled)')
parser.add_option('--rRelAcc',dest='rRelAcc', default="0.001", help='Relative accuracy on r (default: 0.001)')
# Passthrough for anything else combine needs eg for the interference:
#   --combine_opts "--freezeParameters delta --setParameters delta=0"
parser.add_option('--combine_opts',dest='combine_opts', default="", help='Extra options appended verbatim to the combine command line')
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
# Extract limits for every mass point, into files tagged with this extension
tag = "limit_%s" % opt.ext
grid = ",".join("%g" % m for m in mass_points)

cmd = f"""combineTool.py -M AsymptoticLimits -d {datacard} \
  -n .{tag} --parallel {opt.parallel} -m {grid} --run blind --rAbsAcc {opt.rAbsAcc} --rRelAcc {opt.rRelAcc} {opt.combine_opts}"""

subprocess.call(cmd, shell=True)

# Check if the grid is complete
missing = [m for m in mass_points if not os.path.exists("higgsCombine.%s.AsymptoticLimits.mH%g.root" % (tag, m))]
if missing:
  print(" --> [ERROR] %d of %d mass point(s) produced no limit file: %s" % (len(missing), len(mass_points), ", ".join("%g" % m for m in missing)))
  print(" --> Not collecting a partial scan.")
  sys.exit(1)

cmd = f"combineTool.py -M CollectLimits *.{tag}.* -o limits_{opt.ext}.json"
subprocess.call(cmd, shell=True)

# Catch it missing quantiles in the collected json.
quantiles = ["exp-2", "exp-1", "exp0", "exp+1", "exp+2"]
with open("limits_%s.json" % opt.ext) as f: collected = json.load(f)
incomplete = od()
for m in mass_points:
  entry = collected.get("%g" % m, collected.get("%.1f" % m, {}))
  absent = [q for q in quantiles if q not in entry]
  if absent: incomplete["%g" % m] = absent
if incomplete:
  print(" --> [ERROR] %d of %d mass point(s) have an incomplete set of quantiles:" % (len(incomplete), len(mass_points)))
  for m, absent in incomplete.items():
    print("       mH%s missing %s" % (m, ", ".join(absent)))
  print(" --> This normally means the AsymptoticLimits bisection could not resolve r.")
  print(" --> Retry with a smaller --rRelAcc (currently %s); not plotting a partial scan." % opt.rRelAcc)
  sys.exit(1)

cmd = f"python3 {os.path.dirname(__file__)}/plot_limits.py --input limits_{opt.ext}.json --title='{opt.title}' --lumi='{opt.lumi}' --output limits_{opt.ext}"
subprocess.call(cmd, shell=True)
