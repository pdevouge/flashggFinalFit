# Script for running the entire pipeline for flashggFinalFit
import os, sys, yaml, subprocess
from optparse import OptionParser
from collections import OrderedDict as od

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RUNNING FINALFIT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

def get_options():
  parser = OptionParser()
  parser.add_option('--runOnly', dest='run_only', default='', help="Run only given steps (trees, signal, background, datacard, combine)")
  parser.add_option('--doSystematics', dest='do_syst', action='store_true', help="Run with systematics")
  return parser.parse_args()
(opt,args) = get_options()

def leave():
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RUNNING FINALFIT (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  exit(0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract options from config file
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EXTRACTING CONFIG FROM YAML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
if opt.do_syst: print("-----> 'doSystematics': FinalFit will run with systematics included.")
else: print("You have chosen to run FinalFit without systematics. If this was an error please use '--doSystematics'")
if (args_count := len(sys.argv)) > 1:
    config_file = sys.argv[1]
    if os.path.exists( config_file ):
        with open(config_file) as f:
            cfg = yaml.safe_load(f)
else:
  print("[ERROR] no config file was specified. Leaving...")
  leave()

year = cfg["common"]["year"]
ext = cfg["common"]["extension"]
input_dir = os.path.join(os.getcwd(),cfg["tree"]["input_dir"])
cat = cfg["common"]["cats"]
MLow, MHigh = cfg["common"]["binning"].split(",")
MBins = cfg["common"]["nbins"]

# Create config files for each step
# For signal
py_config = f"""_year = '{year}'

signalScriptCfg = {{

  # Setup
  'inputWSDir': '{input_dir}/signal/ws_{cfg["signal"]["procs"]}/',
  'procs': '{cfg["signal"]["procs"]}', # if auto: inferred automatically from filenames
  'cats': '{cat}', # if auto: inferred automatically from workspace
  'ext': '{ext}',
  'analysis': 'highMassAnalysis', # To specify which replacement dataset mapping (defined in ./python/replacementMap.py)
  'year': '%s' % _year, # Use 'combined' if merging all years: not recommended
  'width': '{cfg["signal"]["width"]}',
  'massPoints': '{cfg["signal"]["mass_points"]}',

  # Photon shape systematics
  'scales': 'Scale', # separate nuisance per year
  'scalesCorr': '', # correlated across years
  'scalesGlobal': '', # affect all processes equally, correlated across years
  'smears': 'Smearing', # separate nuisance per year

  # Job submission options
  'batch': 'local', # ['condor','SGE','IC','local']
  'queue': 'espresso',

}}
"""

with open("Signal/config_high_mass.py", "w") as f:
    f.write(py_config)

with open("Interference/config_high_mass.py", "w") as f:
    f.write(py_config)


# For background
py_config = f"""

backgroundScriptCfg = {{

  # Setup
  'inputWS': '{input_dir}/data/ws/{cfg["background"]["file"]}',  # Input RooWorkspace

  'cats': '{cfg["common"]["cats"]}',
  'catOffset': 0,         # No offset needed since this is a single file
  'ext': '{ext}', # Will name output directories like: outputs_fTest_2022inclusive/
  'year': '{year}',         # Shown in plots; adjust if merging multiple years

  # Job submission options
  'batch': 'local',       # You can change to 'condor' if you'd prefer distributed batch mode
  'queue': 'espresso'     # Ignored when batch is 'local'

}}
"""

with open("Background/config_high_mass.py", "w") as f:
    f.write(py_config)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
os.chdir("Trees2WS")
if len(opt.run_only) == 0 or "trees" in opt.run_only:
    # Run Tree2Workspace
    syst_opt = '--doSystematics' if opt.do_syst else ''

    cmd = f"python3 RunWSScripts.py --inputDir {input_dir}/signal/ --inputConfig config_high_mass.py \
        --year {year} --mode trees2ws --batch local --modeOpts \"--minMass {MLow} --maxMass {MHigh} {syst_opt} \""
    subprocess.call(cmd, shell=True)

    cmd = f"python3 RunWSScripts.py --inputDir {input_dir}/data/ --inputConfig config_high_mass.py \
        --year {year} --mode trees2ws_data --batch local --modeOpts \"--applyMassCut --massCutRange {MLow},{MHigh} \""
    subprocess.call(cmd, shell=True)

if len(opt.run_only) == 0 or "signal" in opt.run_only:
    # Run Signal fit
    os.chdir("../Signal")

    if opt.do_syst:
        cmd = f"""python3 RunSignalScripts.py --inputConfig config_high_mass.py --mode calcPhotonSyst \
            --modeOpts \" --nBins {MBins}  --minMass {MLow} --maxMass {MHigh}\""""
        subprocess.call(cmd, shell=True)

    syst_opt = '' if opt.do_syst else '--skipSystematics'

    cmd = f"""python3 RunSignalScripts.py --inputConfig config_high_mass.py --mode signalFit \
        --modeOpts \" --doPlots {syst_opt} --skipVertexScenarioSplit --skipBeamspotReweigh --nBins {MBins}  --minMass {MLow} --maxMass {MHigh} {cfg['signal']['options']} \""""
    subprocess.call(cmd, shell=True)

    cmd = f"""python3 RunPackager.py --cats {cat} --exts {ext} --year {year} \
        --massPoints {cfg["signal"]["mass_points"]} --batch local"""
    subprocess.call(cmd, shell=True)

if len(opt.run_only) == 0 or "interference" in opt.run_only:
    # Run Background fit
    os.chdir("../Interference")

    cmd = f"""python3 RunInterferenceScripts.py --inputConfig config_high_mass.py --mode computeIntf \
        --modeOpts \"  --minMass {MLow} --maxMass {MHigh} \""""
    subprocess.call(cmd, shell=True)

if len(opt.run_only) == 0 or "background" in opt.run_only:
    # Run Background fit
    os.chdir("../Background")

    cmd = f"python3 RunBackgroundScripts.py --inputConfig config_high_mass.py --mode fTestParallel"
    subprocess.call(cmd, shell=True)

if len(opt.run_only) == 0 or "datacard" in opt.run_only:
    # Run Datacard maker
    os.chdir("../Datacard")
    syst_opt = '--doSystematics' if opt.do_syst else ''

    cmd = f"""python3 RunYields.py --inputWSDirMap {year}={input_dir}/signal/ws_{cfg["signal"]["procs"]} {syst_opt}\
        --cats {cat} --procs {cfg["signal"]["procs"]} --ext {ext} --skipCOWCorr --batch local --mass 700 --width {cfg["signal"]["width"]}"""
    print("------>", cmd)
    subprocess.call(cmd, shell=True)

    cmd = f"""python3 makeDatacard.py --ext {ext} --years {year} --prune {syst_opt} \
        --skipCOWCorr --doMCStatUncertainty --saveDataFrame --output Datacard_{ext} --mass 700"""
    print("------>", cmd)
    subprocess.call(cmd, shell=True)

if len(opt.run_only) == 0 or "combine" in opt.run_only:
    # Move everything to Combine area
    os.chdir("../Combine")

    if not os.path.isdir("Models"):
        os.system("mkdir -p Models/signal")
        os.system("mkdir -p Models/background")
        os.system("mkdir -p Models/interference")
    else:
        print("Models directory already exists: models will not be copied into Combine dir. Leaving...")
        leave()

    os.system(f"cp ../Signal/outdir_packaged/CMS-HGG_sigfit_packaged_{cat}_{year}*.root Models/signal/")
    os.system(f"cp ../Background/outdir_{ext}/CMS-HGG_multipdf_{cat}_{year}*.root Models/background/")
    os.system(f"cp ../Interference/outdir_{ext}/computeIntf/output/CMS-HGG_intfm_{cat}_{year}*.root Models/interference/")
    os.system(f"cp ../Datacard/Datacard_{ext}.txt .")