# Script for running combineTool.py -M Impacts (nuisance-parameter impact ranking)
import os, sys, subprocess
from optparse import OptionParser

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RUNNING COMBINE IMPACTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

parser = OptionParser(usage="usage: %prog datacard.txt [options] \nrun with --help to get list of options")
parser.add_option('--outdir',dest='outdir', default="", help='Where to save the impacts (default: cwd)')
parser.add_option('--extension',dest='ext', default="default", help='File extension')
parser.add_option('--mass',dest='mass', default="125", help='Mass point to run the impacts at')
parser.add_option('--parallel',dest='parallel', default=1, type='int', help='Number of nuisance-parameter fits to run at the same time (default: 1, i.e. serial)')
# Passthrough for anything else combine needs eg for the interference:
#   --combine_opts "-t -1 --setParameters r=1,delta=0 --freezeParameters delta"
parser.add_option('--combine_opts',dest='combine_opts', default="", help='Extra options appended verbatim to the combine command line')
(opt,args) = parser.parse_args()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract options from cmd line
if len(args) == 0:
    parser.print_usage()
    exit(1)

datacard = os.path.join(os.getcwd(),args[0])

if opt.outdir:
  print(" --> Impacts will be saved into %s" %opt.outdir)
  if not os.path.isdir( opt.outdir ): os.makedirs(opt.outdir)
  os.chdir(opt.outdir)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3-step combineTool.py workflow: doInitialFit -> doFits -> collect to json.
# See Combine/README.md#Impacts.
tag = opt.ext
initial_fit = "higgsCombine_initialFit_%s.MultiDimFit.mH%s.root" % (tag, opt.mass)

cmd = f"""combineTool.py -M Impacts -d {datacard} -m {opt.mass} \
  -n {tag} --doInitialFit --robustFit 1 {opt.combine_opts}"""
subprocess.call(cmd, shell=True)

if not os.path.exists(initial_fit):
  print(" --> [ERROR] Initial fit produced no %s. Not running the per-nuisance scans." % initial_fit)
  sys.exit(1)

cmd = f"""combineTool.py -M Impacts -d {datacard} -m {opt.mass} \
  -n {tag} --robustFit 1 --doFits --parallel {opt.parallel} {opt.combine_opts}"""
subprocess.call(cmd, shell=True)

json_out = "impacts_%s.json" % opt.ext
cmd = f"combineTool.py -M Impacts -d {datacard} -m {opt.mass} -n {tag} -o {json_out}"
subprocess.call(cmd, shell=True)

if not os.path.exists(json_out):
  print(" --> [ERROR] Impacts collection produced no %s. Some per-nuisance fit(s) likely failed "
        "-- check higgsCombine_paramFit_%s_*.MultiDimFit.mH%s.root in this directory." % (json_out, tag, opt.mass))
  sys.exit(1)

# ROOT files stay in this (Impacts) dir; the plot goes into a Plots dir shared with the Limits step
plots_dir = os.path.join("..", "Plots")
if not os.path.isdir(plots_dir): os.makedirs(plots_dir)
cmd = f"plotImpacts.py -i {json_out} -o {os.path.join(plots_dir, f'impacts_{opt.ext}')}"
subprocess.call(cmd, shell=True)
